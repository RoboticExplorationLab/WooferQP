function control!(
    torques::AbstractVector,
    q::AbstractVector,
    q̇::AbstractVector,
    t::AbstractFloat,
    param::ControllerParams,
) 
    rot = UnitQuaternion(q[4], q[5], q[6], q[7])
    mrp = MRP(rot)
    ω = rot \ q̇[SUnitRange(4, 6)]

    x_est = [
                q[SUnitRange(1, 3)]; 
                Rotations.params(mrp); 
                q̇[SUnitRange(1, 3)]; 
                ω
            ]

    # annoying way to get rid of knee joint measurements
    joint_pos = q[@SVector [8,9,11,13,14,16,18,19,21,23,24,26]]
    joint_vel = q̇[@SVector [7,8,10,12,13,15,17,18,20,22,23,25]]

    # get current leg positions
    param.cur_foot_loc = ForwardKinematicsAll(joint_pos)

    # prev phase -> cur_phase check contacts to regenerate swing
    param.cur_phase = get_phase(t, param)
    param.cur_phase_time = get_phase_time(t, param.cur_phase, param)
    
    param.active_feet = param.gait.contact_phases[param.cur_phase]
    coordinate_expander!(param.active_feet_12, param.active_feet)

    v_i = @SVector [x_est[7], x_est[8], x_est[9]]
    v_b = rot \ v_i # inertial -> body

    # swing leg
    for i = 1:4
        # calculate footstep and generate trajectory (stored in swing params) if needed
        if param.gait.contact_phases[param.prev_phase][i] == 1
            if param.gait.contact_phases[param.cur_phase][i] == 0
                param.next_foot_loc[i] =
                    footstep_location(v_b, ω[3], param.cur_phase, i, param)

                # make sure MPC accounts for this next foot location
                param.planner_foot_loc[i] = param.next_foot_loc[i]

                foot_trajectory(
                    -v_b,
                    -v_b,
                    t,
                    t + param.gait.phase_times[param.cur_phase],
                    i,
                    param,
                    regen_z = true,
                )
                param.last_replan_t = t
            end
        end

        # actually calculate swing torques
        if param.gait.contact_phases[param.cur_phase][i] == 0
            # calculate current foot tip velocity
            J = LegJacobian(joint_pos[SLegIndexToRange(i)], i)
            cur_foot_vel_i = J * joint_vel[SLegIndexToRange(i)]

            if (t - param.last_replan_t) > param.replan_update
                param.next_foot_loc[i] =
                    footstep_location(v_b, ω[3], param.cur_phase, i, param)

                # make sure MPC accounts for this next foot location
                param.planner_foot_loc[i] = param.next_foot_loc[i]

                foot_trajectory(
                    cur_foot_vel_i,
                    -v_b,
                    t,
                    (t - param.cur_phase_time) +
                    param.gait.phase_times[param.cur_phase],
                    i,
                    param,
                    regen_z = false,
                )
                param.last_replan_t = t
            end

            swing_torque_i = swing_torques(
                cur_foot_vel_i,
                joint_pos[SLegIndexToRange(i)],
                t,
                i,
                param,
            )
            param.swing_torques[LegIndexToRange(i)] .= swing_torque_i
        end
    end
    param.prev_phase = param.cur_phase

    if (t - param.last_t) >= param.mpc_update
        # update MPC forces
        reference_trajectory!(x_est, param)
        foot_history!(t, param)
        # @time foot_forces!(x_est, param)
        foot_forces!(x_est, param)

        param.last_t = t
    end

    # needs to be negative so force is exerted by body on world
    param.mpc_torques = Force2Torque(-param.forces, joint_pos)

    torques .=
        param.active_feet_12 .* param.mpc_torques +
        (ones(12) - param.active_feet_12) .* param.swing_torques
end
