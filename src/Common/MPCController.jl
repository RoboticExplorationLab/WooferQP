include("Structs.jl")
include("QPSolverSparse.jl")
include("SwingLegController.jl")
include("Gait.jl")
include("FootstepPlanner.jl")

function mpcControlWoofer!(torques::Vector{T}, x_est::Vector{T}, t::T, joint_pos::Vector{T}, joint_vel::Vector{T}, param::ControllerParams) where {T<:Number}
	# get current leg positions
	param.cur_foot_loc = ForwardKinematicsAll(joint_pos)

	# prev phase -> cur_phase check contacts to regenerate swing
	param.cur_phase = getPhase(t, param)
	param.active_feet = param.gait.contact_phases[:, param.cur_phase]
	coordinateExpander!(param.active_feet_12, param.active_feet)

	# swing leg
	for i in 1:4
		# calculate footstep and generate trajectory (stored in swing params) if needed
		if param.gait.contact_phases[i, param.prev_phase] == 1
			if param.gait.contact_phases[i, param.cur_phase] == 0
				R = QuatToRotMatrix(ThreeParamToQuat(x_est[4:6]))'
				v_b = R*x_est[7:9]
				ω = R*x_est[10:12]

				next_foot_phase = nextPhase(param.cur_phase, param)
				param.next_foot_loc[LegIndexToRange(i)] = nextFootstepLocation(v_b, ω[3], next_foot_phase, i, param)

				# make sure MPC accounts for this next foot location
				param.planner_foot_loc[LegIndexToRange(i)] .= param.next_foot_loc[LegIndexToRange(i)]

				print("Next Leg $i Position: ")
				print(param.next_foot_loc[LegIndexToRange(i)])
				println()
				generateFootTrajectory(x_est[7:9], t, t+param.gait.phase_times[param.cur_phase], i, param)
			end
		end

		# actually calculate swing torques
		if param.gait.contact_phases[i, param.cur_phase] == 0
			# calculate current foot tip velocity
			J = LegJacobian(joint_pos[LegIndexToRange(i)], i)
			cur_foot_vel_i = J * joint_vel[LegIndexToRange(i)]

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

		param.last_t = t
	end

	# needs to be negative so force is exerted by body on world
	param.mpc_torques = Force2Torque(-param.forces, joint_pos)

	torques .= param.active_feet_12 .* param.mpc_torques + (ones(12)-param.active_feet_12) .* param.swing_torques
end
