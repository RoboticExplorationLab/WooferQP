"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

using ECOS
using OSQP
using LinearAlgebra
using Convex

include("Config.jl")
include("Quaternions.jl")

@with_kw struct MPCControllerParams
	# discretization length
	dt::Float64

	# planning horizon length
	N::Int64
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

#TODO fix StaticArray/DynamicArray mix
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

function generateReferenceTrajectory!(x_ref::Array{T, 2}, x_curr::Vector{T}, x_des::Vector{T}, mpc_config::MPCControllerParams) where {T<:Number}
	# TODO: integrate the x,y,ψ position from the reference

	x_diff = 1/mpc_config.N * (x_des - x_curr)
	x_ref[:,1] .= x_curr
	for i in 1:mpc_config.N
		if i<N
			x_ref[:, i+1] .= x_curr + i*x_diff
		else
			x_ref[:, i+1] .= x_des
		end
	end
end

function solveFootForces!(forces::Vector{T}, x_ref::Array{T, 2}, contacts::Array{Int,2}, foot_locs::Array{T,2}, mpc_config::MPCControllerParams, use_lqr::Bool=false) where {T<:Number}
	# x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
	# contacts: 4xN+1 matrix of foot contacts over the planning horizon
	# foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon

	n = 12
	m = 12

	N = mpc_config.N
	dt = mpc_config.dt

	x = Variable((N+1)*n)
	u = Variable((N)*n)

	if use_lqr
		u_ref_i = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.SPRUNG_MASS*9.81/4
		@warn "need to add cost to go matrix"
	else
		Q = Diagonal(repeat([1e3, 1e3, 5e4, 1e3, 1e3, 1e3, 1e2, 1e2, 1e2, 1, 1, 1e2], N+1))
		u_ref_i = zeros(12)
	end

	u_ref = repeat(u_ref_i, N)

	R = Diagonal(repeat([1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4], N))

	x_ref_reshaped = reshape(x_ref, (N+1)*n, 1)

	# problem = minimize(quadform(x-x_ref, Q) + quadform(u-u_ref, R))
	problem = minimize(quadform(x, Q))

	# constants
	μ = 0.7
	min_vert_force = 1
	max_vert_force = 133

	for i=0:N-1
		# Control constraints
		for j=1:4
			# use absolute value friction constraint:
			# problem.constraints += abs(u[select12_3(i,j,1)]) <= μ*u[select12_3(i,j,3)]
			# problem.constraints += abs(u[select12_3(i,j,2)]) <= μ*u[select12_3(i,j,3)]

			# convert absolute value constraint to linear inequality:
			problem.constraints += u[select12_3(i,j,1)] <= μ*u[select12_3(i,j,3)]
			problem.constraints += u[select12_3(i,j,1)] >= -μ*u[select12_3(i,j,3)]
			problem.constraints += u[select12_3(i,j,2)] <= μ*u[select12_3(i,j,3)]
			problem.constraints += u[select12_3(i,j,2)] >= -μ*u[select12_3(i,j,3)]

			problem.constraints += u[select12_3(i,j,3)] >= min_vert_force
			problem.constraints += u[select12_3(i,j,3)] <= max_vert_force
		end

		A_c_i = LinearizedContinuousDynamicsA(x_ref[select12(i)], WOOFER_CONFIG.INERTIA)
		B_c_i = LinearizedContinuousDynamicsB(x_ref[select12(i)], foot_locs[:, i+1], contacts[:, i+1], WOOFER_CONFIG.INERTIA)
		d_c_i = NonLinearContinuousDynamics(x_ref[select12(i)], u_ref[select12(i)], foot_locs[:, i+1], WOOFER_CONFIG.INERTIA)

		# Euler Discretization
		A_d_i = I + A_c_i*dt
		B_d_i = B_c_i*dt
		d_d_i = d_c_i*dt

		# Dynamics constraints
		problem.constraints += x[select12(i+1)] == A_d_i*x[select12(i)] + B_d_i*u[select12(i)] + d_d_i

		print(d_d_i)
	end

	solve!(problem, OSQP.Optimizer)

	return (u.value, x.value)
end

function select12(i)
	return (12*i+1):(12*i+12)
end

function select12_3(i,j,k)
	return 12*i+3*(j-1)+k
end
