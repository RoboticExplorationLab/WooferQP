"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

using ForwardDiff

function mass(contacts)
	return woofer.inertial.sprung_mass + sum([2*woofer.inertial.lower_link_mass*(1-contacts[i]) for i=1:4])
end

function NonLinearContinuousDynamics(x, u, r, contacts)
	J = woofer.inertial.body_inertia

	x_d_1_3 = x[7:9]
	x_d_4_6 = 0.5*V()*L_q(ThreeParamToQuat(x[4:6]))*VecToQuat(x[10:12])
	x_d_7_9 = 1/woofer.inertial.sprung_mass * sum(repeat(contacts', 3, 1) .* reshape(u, 3, 4), dims=2) + [0; 0; -9.81]

	torque_sum = zeros(3)
	for i=1:4
		torque_sum += contacts[i]*SkewSymmetricMatrix(r[LegIndexToRange(i)])*QuatToRotMatrix(ThreeParamToQuat(x[4:6]))'*u[LegIndexToRange(i)]
	end
	x_d_10_12 = inv(J)*(-SkewSymmetricMatrix(x[10:12])*J*x[10:12] + torque_sum)

	return [x_d_1_3..., x_d_4_6..., x_d_7_9..., x_d_10_12...]
end

function LinearizedContinuousDynamicsA(x, u, r, contacts)
	return ForwardDiff.jacobian((x_var)->NonLinearContinuousDynamics(x_var, u, r, contacts), x)
end

function LinearizedContinuousDynamicsB(x, u, r, contacts)
	return ForwardDiff.jacobian((u_var)->NonLinearContinuousDynamics(x, u_var, r, contacts), u)
end

function generateReferenceTrajectory!(x_curr::Vector{T}, param::ControllerParams) where {T<:Number}
	# TODO: integrate the x,y,ψ position from the reference
	α = collect(range(0, 1, length=param.N+1))
	# interpolate everything but x, y position
	if param.vel_ctrl
		interp_indices = 3:12
	else
		interp_indices = 1:12
	end

	for i in 1:param.N+1
		param.x_ref[interp_indices, i] .= (1-α[i]) * x_curr[interp_indices] + α[i]*param.x_des[interp_indices]
	end

	if param.vel_ctrl
		# integrate x, y position
		integ_indices = 1:2
		param.x_ref[integ_indices, :] .= repeat(x_curr[integ_indices, 1], 1, param.N+1) + cumsum(param.x_ref[7:8, :]*param.optimizer.dt, dims=2)
	end
end

function solveFootForces!(param::ControllerParams) where {T<:Number}
	# x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
	# contacts: 4xN+1 matrix of foot contacts over the planning horizon
	# foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon
	N = param.N
	dt = param.optimizer.dt
	n = 12
	m = 12

	param.optimizer.x_ref_reshaped .= reshape(param.x_ref, (N+1)*n, 1)[:]

	for i=0:N-1
		A_c_i = LinearizedContinuousDynamicsA(param.x_ref[select12(i)], param.optimizer.u_ref[select12(i)], param.foot_locs[:, i+1], param.contacts[:, i+1])
		B_c_i = LinearizedContinuousDynamicsB(param.x_ref[select12(i)], param.optimizer.u_ref[select12(i)], param.foot_locs[:, i+1], param.contacts[:, i+1])
		d_c_i = NonLinearContinuousDynamics(param.x_ref[select12(i)], param.optimizer.u_ref[select12(i)], param.foot_locs[:, i+1], param.contacts[:, i+1])

		# Euler Discretization
		param.optimizer.A_d[i+1] .= Array(I + A_c_i*dt)
		param.optimizer.B_d[i+1] .= Array(B_c_i*dt)
		param.optimizer.d_d[i+1] .= Array(d_c_i*dt + param.x_ref[select12(i)])

		if param.use_lqr
			param.optimizer.u_ref[select12(i)] .= pinv(param.optimizer.B_d[i+1])*((I - param.optimizer.A_d[i+1])*param.optimizer.x_ref_reshaped[select12(i)] - param.optimizer.d_d[i+1])
		end
	end

	if param.use_lqr
		V = dare(param.optimizer.A_d[N], param.optimizer.B_d[N], param.optimizer.Q[1:n, 1:n], param.optimizer.R[1:m, 1:m])
		param.optimizer.Q[(N*n+1):((N+1)*n), (N*n+1):((N+1)*n)] .= V
	end

	solve!(param.optimizer.model)

	param.forces .= value.(param.optimizer.model, param.optimizer.u)[select12(0)]
end

function select12(i)
	return (12*i+1):(12*i+12)
end

function select12_3(i,j,k)
	return 12*i+3*(j-1)+k
end
