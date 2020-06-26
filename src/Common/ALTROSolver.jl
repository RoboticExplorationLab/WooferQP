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

function TrajectoryOptimization.dynamics(model, x, u, t)
	k = searchsortedfirst(model.times, t)
	if model.times[k] != t
		k -= 1
	end

	# A = LinearizedContinuousDynamicsA(model.x_ref[k], model.u_ref[k], model.foot_locs[k], model.contacts[k])
	# B = LinearizedContinuousDynamicsB(model.x_ref[k], model.u_ref[k], model.foot_locs[k], model.contacts[k])
	# d = NonLinearContinuousDynamics(model.x_ref[k], model.u_ref[k], model.foot_locs[k], model.contacts[k])
	# x_dot = A*(x - model.x_ref[k]) + B*(u - model.u_ref[k]) + d 

	x_dot = NonLinearContinuousDynamics(x, u, model.foot_locs[k], model.contacts[k])

	return x_dot
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
	# contacts: 4xN matrix of foot contacts over the planning horizon
	# foot_locs: 12xN matrix of foot location in body frame over planning horizon
	opt = param.optimizer

	N = param.N
	dt = opt.dt
	n = 12
	m = 12

	opt.model.x_ref .= reshape(param.x_ref, (N+1)*n, 1)[:]
	opt.model.contacts .= param.contacts
	opt.model.foot_locs .= param.foot_locs

	opt.problem = Problem()

	solve!(param.optimizer.model)

	param.forces .= value.(param.optimizer.model, param.optimizer.u)[select12(0)]
end

function select12(i)
	return (12*i+1):(12*i+12)
end

function select12_3(i,j,k)
	return 12*i+3*(j-1)+k
end
