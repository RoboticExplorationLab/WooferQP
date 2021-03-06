include("Structs/OptimizerParams.jl")
include("Structs/SwingLegParams.jl")
include("Structs/GaitParams.jl")
include("Structs/ControllerParams.jl")
include("QPSolverSparse.jl")
include("SwingLegController.jl")
include("Gait.jl")
include("FootstepPlanner.jl")

function mpcControlWoofer!(torques::Vector{T}, x_est::Vector{T}, t::T, joint_pos::Vector{T}, joint_vel::Vector{T}, param::ControllerParams) where {T<:Number}
	# get current leg positions
	param.cur_foot_loc = ForwardKinematicsAll(joint_pos)

	# prev phase -> cur_phase check contacts to regenerate swing
	(param.cur_phase, param.cur_phase_time) = getPhase(t, param, return_time=true)
	param.active_feet = param.gait.contact_phases[:, param.cur_phase]
	coordinateExpander!(param.active_feet_12, param.active_feet)

	R = QuatToRotMatrix(ThreeParamToQuat(x_est[4:6]))' # inertial -> body
	v_b = R*x_est[7:9]
	ω = R*x_est[10:12]

	# swing leg
	for i in 1:4
		# calculate footstep and generate trajectory (stored in swing params) if needed
		if param.gait.contact_phases[i, param.prev_phase] == 1
			if param.gait.contact_phases[i, param.cur_phase] == 0
				param.next_foot_loc[LegIndexToRange(i)] = nextFootstepLocation(v_b, ω[3], param.cur_phase, i, param)

				# make sure MPC accounts for this next foot location
				param.planner_foot_loc[LegIndexToRange(i)] .= param.next_foot_loc[LegIndexToRange(i)]

				generateFootTrajectory(	-x_est[7:9],
										-x_est[7:9],
										t,
										t+param.gait.phase_times[param.cur_phase],
										i,
										param,
										regen_z=true)
				param.last_replan_t = t
			end
		end

		# actually calculate swing torques
		if param.gait.contact_phases[i, param.cur_phase] == 0
			# calculate current foot tip velocity
			J = LegJacobian(joint_pos[LegIndexToRange(i)], i)
			cur_foot_vel_i = J * joint_vel[LegIndexToRange(i)]

			if (t - param.last_replan_t) > param.replan_update
				param.next_foot_loc[LegIndexToRange(i)] = nextFootstepLocation(v_b, ω[3], param.cur_phase, i, param)

				# make sure MPC accounts for this next foot location
				param.planner_foot_loc[LegIndexToRange(i)] .= param.next_foot_loc[LegIndexToRange(i)]

				generateFootTrajectory( cur_foot_vel_i,
										-x_est[7:9],
										t,
										(t-param.cur_phase_time)+param.gait.phase_times[param.cur_phase],
										i,
										param,
										regen_z=false)
				param.last_replan_t = t
			end

			swing_torque_i = calcSwingTorques(cur_foot_vel_i, joint_pos[LegIndexToRange(i)], t, i, param)
			param.swing_torques[LegIndexToRange(i)] .= swing_torque_i
		end
	end
	param.prev_phase = param.cur_phase

	if (t - param.last_t) >= param.mpc_update
		# update MPC forces
		generateReferenceTrajectory!(x_est, param)
		constructFootHistory!(t, param)
		solveFootForces!(param)

		println("X Velocity: ", x_est[7])

		param.last_t = t
	end

	# needs to be negative so force is exerted by body on world
	param.mpc_torques = Force2Torque(-param.forces, joint_pos)

	torques .= param.active_feet_12 .* param.mpc_torques + (ones(12)-param.active_feet_12) .* param.swing_torques
end
