"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

function reference_trajectory!(x_curr::AbstractVector{T}, param::ControllerParams) where {T<:Number}
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

function foot_forces!(x_curr::AbstractVector{T}, param::ControllerParams) where {T<:Number}
	# x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
	# contacts: 4xN+1 matrix of foot contacts over the planning horizon
	# foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon
	opt = param.optimizer

	tf = param.N * opt.dt
	n = 12
	m = 12

	# TrajectoryOptimization.reset!(opt.constraints)

	# TODO: standardize between QP and ALTRO
	for i=1:(param.N+1)
		opt.model.x_ref[i] = SVector{12}(param.x_ref[:, i])
		opt.model.contacts[i] = SVector{4}(param.contacts[:, i])
		opt.model.foot_locs[i] = SVector{12}(param.foot_locs[:, i])
	end

	opt.solver.solver_uncon.x0 .= x_curr

	initial_states!(opt.problem, opt.X0)
	initial_controls!(opt.problem, opt.U0)
	@time solve!(opt.solver)

	# allocs = @allocated(solve!(opt.solver))
	# println("Allocations: ", allocs)

	opt.X0 = states(opt.solver)
	opt.U0 = controls(opt.solver)

	param.forces .= opt.U0[1]
end
