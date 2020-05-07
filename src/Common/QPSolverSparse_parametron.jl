"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

using OSQP
using LinearAlgebra
using Parametron

include("Config.jl")
include("Quaternions.jl")

struct MPCControllerParams
	# discretization length
	dt::Float64

	# planning horizon length
	N::Int64

	# Parametron Variables:
	model::Model
	x::Array{Variable, 1}
	u::Array{Variable, 1}

	# Parameter Values (to be updated in place)
	x_ref_reshaped::Vector{Float64}
	u_ref::Vector{Float64}
	A_d::Array{Array{Float64, 2}, 1}
	B_d::Array{Array{Float64, 2}, 1}
	d_d::Array{Array{Float64, 1}, 1}
	Q::Array{Float64, 2}
	R::Array{Float64, 2}

	function MPCControllerParams(dt::Float64, N::Int64; n::Int64=12, m::Int64=12)
		# initialize model and variables
		model = Model(OSQP.Optimizer(verbose=false))
		x = [Variable(model) for _ = 1:((N+1)*n)]
		u = [Variable(model) for _ = 1:((N)*n)]

		# initialize quadratic cost parameters
		Q = Diagonal(repeat([1e4, 1e4, 5e4, 4e3, 4e3, 1e2, 1e4, 1e4, 1e2, 1, 1, 1e2], N+1))
		R = Diagonal(repeat([1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4], N))

		x_ref_reshaped = zeros((N+1)*n)
		u_ref = zeros((N)*n)

		A_d = [zeros(n,n) for i=1:N]
		B_d = [zeros(n,m) for i=1:N]
		d_d = [zeros(n) for i=1:N]

		Q_param = Parameter(()->Q, model)
		R_param = Parameter(()->R, model)

		x_ref_param = Parameter(()->x_ref_reshaped, model)
		u_ref_param = Parameter(()->u_ref, model)

		A_d_param = [Parameter(()->A_d[i], model) for i=1:N]
		B_d_param = [Parameter(()->B_d[i], model) for i=1:N]
		d_d_param = [Parameter(()->d_d[i], model) for i=1:N]

		# constants
		μ = 0.7
		min_vert_force = 1
		max_vert_force = 133

		# objective
		@objective(model, Minimize, transpose(x-x_ref_param)*Q_param*(x-x_ref_param) + transpose(u-u_ref_param)*R_param*(u-u_ref_param))

		# define all constraints
		@constraint(model, x[select12(0)] == x_ref_param[select12(0)])
		for i=0:N-1
			# Control constraints
			for j=1:4
				# convert absolute value constraint to linear inequality:
				@constraint(model, u[select12_3(i,j,1)] <= μ*u[select12_3(i,j,3)])
				@constraint(model, u[select12_3(i,j,1)] >= -μ*u[select12_3(i,j,3)])
				@constraint(model, u[select12_3(i,j,2)] <= μ*u[select12_3(i,j,3)])
				@constraint(model, u[select12_3(i,j,2)] >= -μ*u[select12_3(i,j,3)])

				@constraint(model, u[select12_3(i,j,3)] >= min_vert_force)
				@constraint(model, u[select12_3(i,j,3)] <= max_vert_force)
			end

			# Dynamics constraints
			@constraint(model, x[select12(i+1)] == A_d_param[i+1]*x[select12(i)] + B_d_param[i+1]*u[select12(i)] + d_d_param[i+1])
		end

		new(dt, N, model, x, u, x_ref_reshaped, u_ref, A_d, B_d, d_d, Q, R)
	end
end

function NonLinearContinuousDynamics(x, u, r, J)
	x_dot = zeros(12)

	x_dot[1:3] = x[7:9]
	x_dot[4:6] = 0.5*V()*L_q(ThreeParamToQuat(x[4:6]))*VecToQuat(x[10:12])
	x_dot[7:9] = 1/WOOFER_CONFIG.SPRUNG_MASS * sum(reshape(u, 3, 4), dims=2) + [0; 0; -9.81]

	torque_sum = zeros(3)
	for i=1:4
		torque_sum += SkewSymmetricMatrix(r[LegIndexToRange(i)])*QuatToRotMatrix(ThreeParamToQuat(x[4:6]))'*u[LegIndexToRange(i)]
	end
	x_dot[10:12] = inv(J)*(-SkewSymmetricMatrix(x[10:12])*J*x[10:12] + torque_sum)

	return x_dot
end

function LinearizedContinuousDynamicsA(x, J)
    x_ϕ = x[4:6]
    x_w = x[10:12]
    dvdv = Diagonal(ones(3))
    dphidphi = 0.5 * V() * R_q(ThreeParamToQuat(x_w)) * [-x_ϕ' / sqrt(abs(1 - x_ϕ' * x_ϕ)); Diagonal(ones(3))]
    dphidw = 0.5 * V() * L_q(ThreeParamToQuat(x_ϕ)) * V()'
    dwdw = -inv(J) * (SkewSymmetricMatrix(x_w) * J - SkewSymmetricMatrix(Vector(J * x_w)))

    r1 = [zeros(3, 6) dvdv zeros(3, 3)]
    r2 = [zeros(3, 3) dphidphi zeros(3, 3) dphidw]
    r3 = zeros(3, 12)
    r4 = [zeros(3,9) dwdw]
    R = [r1; r2; r3; r4]
    return R
end

function LinearizedContinuousDynamicsB(x, r, contacts, J)
	b1 = zeros(3, 12)
	b2 = zeros(3, 12)
	b3 = 1/WOOFER_CONFIG.SPRUNG_MASS*repeat(Diagonal(ones(3)), 1, 4)

	b4 = zeros(3,12)
	for i=1:4
		b4[:, LegIndexToRange(i)] = inv(J)*contacts[i]*SkewSymmetricMatrix(r[LegIndexToRange(i)])*QuatToRotMatrix(ThreeParamToQuat(x[4:6]))'
	end

	return [b1; b2; b3; b4]
end

function generateReferenceTrajectory!(x_ref::Array{T, 2}, x_curr::Vector{T}, x_des::Vector{T}, mpc_config::MPCControllerParams, integrate::Bool=false) where {T<:Number}
	# TODO: integrate the x,y,ψ position from the reference
	α = collect(range(0, 1, length=mpc_config.N+1))
	# interpolate everything but x, y position
	if integrate
		interp_indices = 3:12
	else
		interp_indices = 1:12
	end

	for i in 1:mpc_config.N+1
		x_ref[interp_indices, i] .= (1-α[i]) * x_curr[interp_indices] + α[i]*x_des[interp_indices]
	end

	if integrate
		# integrate x, y position
		integ_indices = 1:2
		x_ref[integ_indices, :] .= repeat(x_curr[integ_indices, 1], 1, mpc_config.N+1) + cumsum(x_ref[7:8, :]*mpc_config.dt, dims=2)
	end
end

function solveFootForces!(forces::Vector{T}, x_ref::Array{T, 2}, contacts::Array{Int,2}, foot_locs::Array{T,2}, mpc_config::MPCControllerParams, use_lqr::Bool=false) where {T<:Number}
	# x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
	# contacts: 4xN+1 matrix of foot contacts over the planning horizon
	# foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon
	N = mpc_config.N
	dt = mpc_config.dt
	n = 12
	m = 12

	mpc_config.x_ref_reshaped .= reshape(x_ref, (N+1)*n, 1)[:]

	for i=0:N-1
		A_c_i = LinearizedContinuousDynamicsA(x_ref[select12(i)], WOOFER_CONFIG.INERTIA)
		B_c_i = LinearizedContinuousDynamicsB(x_ref[select12(i)], foot_locs[:, i+1], contacts[:, i+1], WOOFER_CONFIG.INERTIA)
		d_c_i = NonLinearContinuousDynamics(x_ref[select12(i)], mpc_config.u_ref[select12(i)], foot_locs[:, i+1], WOOFER_CONFIG.INERTIA)

		# Euler Discretization
		mpc_config.A_d[i+1] .= Array(I + A_c_i*dt)
		mpc_config.B_d[i+1] .= Array(B_c_i*dt)
		mpc_config.d_d[i+1] .= Array(d_c_i*dt)

		if use_lqr
			mpc_config.u_ref[select12(i)] .= pinv(mpc_config.B_d[i+1])*((I - mpc_config.A_d[i+1])*mpc_config.x_ref_reshaped[select12(i)] - mpc_config.d_d[i+1])
		end
	end

	if use_lqr
		V = dare(mpc_config.A_d[N], mpc_config.B_d[N], mpc_config.Q[1:n, 1:n], mpc_config.R[1:m, 1:m])
		mpc_config.Q[(N*n+1):((N+1)*n), (N*n+1):((N+1)*n)] .= V
	else
		mpc_config.Q[(N*n+1):((N+1)*n), (N*n+1):((N+1)*n)] .= mpc_config.Q[1:n, 1:n]
	end

	solve!(mpc_config.model)

	forces .= value.(mpc_config.model, mpc_config.u)[1:12]
end

function select12(i)
	return (12*i+1):(12*i+12)
end

function select12_3(i,j,k)
	return 12*i+3*(j-1)+k
end
