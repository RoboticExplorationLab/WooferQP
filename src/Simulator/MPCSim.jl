using LinearAlgebra
using Parametron
using OSQP
using StaticArrays
using DataStructures

include("XMLParser.jl")
include("../Common/QuadrupedDynamics.jl")
include("../Common/MPCControl/MPCControl.jl")
include("../Common/Utilities.jl")
include("../Common/Quaternions.jl")

using .QuadrupedDynamics
import .MPCControl

function simulate()
	woofer = WooferConfig("../Common/Woofer.yaml")

    # Create the robot XML file
    ParseXML(woofer)
    # Get model from local directory
    s = loadmodel("woofer_out.xml", 1200, 900)

    # Pre-allocate memory for time-delayed torques, joint pos/vel, and mocap data
    lagged_control = zeros(12)
    lagged_joint = zeros(24)
    lagged_mocap = zeros(7)

    # Pre-allocate memory for actuator torques
    actuator_torques = zeros(12)

    last_actuator_timestamp = 0.0

    # Mujoco data and model
    d = s.d
    m = s.m

    # Initialize update rates
    ODRIVE_DT = 0.001
    SIM_DATA_DT = 0.001
    MOCAP_DT = 0.009
    BUFFER_UPDATE_DT = 0.001

    last_t = 0.0
    last_sim_data_t = 0.0
    last_mocap_t = 0.0
    last_buffer_update_t = 0.0

	# CircularBuffer for simulating control latency
	control_latency_frames = 1 # 1ms
    control_buffer = CircularBuffer{Vector{Float64}}(control_latency_frames)
    fill!(control_buffer, zeros(12))

	# CircularBuffer for simulating mocap latency
	mocap_latency_frames = 1 # 1ms
	init_mocap_measurement = [0, 0, 0.32, 1, zeros(3)...]
    mocap_buffer = CircularBuffer{Vector{Float64}}(mocap_latency_frames)
    fill!(mocap_buffer, init_mocap_measurement)

	# CircularBuffer for simulating measurement latency
	joint_latency_frames = 1 # 1ms
	init_joint_measurement = zeros(24)
    joint_buffer = CircularBuffer{Vector{Float64}}(joint_latency_frames)
    fill!(joint_buffer, init_joint_measurement)


	# Allocate everything once

	# desired position
	xy_vel = [0.5, 0.0]
	ω = 0.0
	x_des = [0.0, 0.0, 0.28, zeros(3)..., xy_vel..., zeros(3)..., ω]
	# true state information
	x = zeros(13)
	x_true = zeros(12)
	joint_pos = zeros(12)
	joint_vel = zeros(12)
	contacts = zeros(4)

	# Δt for MPC dynamics
	planning_dt = 0.05
	# number of planning horizon
	N = 15
	# frequency of force calculation
	mpc_update = 0.05

	# state penalty in QP
	# q = [1e4, 1e4, 5e4, 4e3, 4e3, 1e3, 1e4, 1e4, 1e2, 1e2, 1e2, 1e3]
	q = [0, 0, 5e4, 1e3, 1e5, 1e3, 1e4, 1e4, 1e2, 1e2, 1e2, 1e4]
	# q = [0, 0, 5e4, 1e5, 1e5, 1e2, 1e4, 1e4, 1e2, 1e4, 1e4, 1e3]
	# control penalty in QP
	r = [1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4]

	optimizer = MPCControl.OptimizerParams(planning_dt, N, q, r)

	nom_foot_loc = ForwardKinematicsAll(zeros(12))
	offset = [1 -1 1 -1]
	Δx = 0.0
	Δy = 0.025
	for i=1:4
		nom_foot_loc[LegIndexToRange(i)] += [Δx, Δy*offset[i], 0]
	end

	# gait = createStandingGait()
	# gait = createTrotGait(stance_time=0.15, swing_time=0.15)
	gait = MPCControl.createPronkGait(stance_time=0.2, flight_time=0.1)
	# gait = createPaceGait(stance_time=0.1, swing_time=0.15)
	# gait = createBoundGait(front_time=0.15, back_time=0.15, stance_time=0.05)
	swing = MPCControl.SwingLegParams(-0.20, 100, 1)

	use_lqr = false # use lqr in cost to go
	vel_ctrl = false # integrate positions, interpolate velocities
	param = MPCControl.ControllerParams(N, mpc_update, x_des, use_lqr, vel_ctrl, nom_foot_loc, optimizer, gait, swing)

    # Loop until the user closes the window
    WooferSim.alignscale(s)
    println("Press enter to begin simulation...")
    readline(stdin)
    println("Starting sim.")

    while !GLFW.WindowShouldClose(s.window)
        ### basically sim step so things don't have to be defined multiple times
        if s.paused
            if s.pert[].active > 0
                mjv_applyPerturbPose(m, d, s.pert, 1)  # move mocap and dynamic bodies
                mj_forward(m, d)
            end
        else
            # Slow motion factor: 10x
            # Regular time needs a factor of 0.5
            factor = s.slowmotion ? 10.0 : 0.5

            # Advance effective simulation time by 1/refreshrate
            startsimtm = d.d[].time
            starttm = time()
            refreshtm = 1.0 / (factor * s.refreshrate)
            updates = refreshtm / m.m[].opt.timestep

            steps = round(Int, round(s.framecount + updates) - s.framecount)
            s.framecount += updates

            for i = 1:steps
                # clear old perturbations, apply new
                d.xfrc_applied .= 0.0

                if s.pert[].select > 0
                    mjv_applyPerturbPose(m, d, s.pert, 0) # move mocap bodies only
                    mjv_applyPerturbForce(m, d, s.pert)
                end

                t = d.d[].time

                # Check if it's time to update the delay buffers for control, joint pos/vel, and mocap
                if (t - last_buffer_update_t) > BUFFER_UPDATE_DT
                    # Add lagged motion capture data
					push!(mocap_buffer, s.d.qpos[1:7])
                    lagged_mocap = popfirst!(mocap_buffer)

                    # Add lagged joint pos/vel data
					push!(joint_buffer, s.d.sensordata[7:30])
                    lagged_joint = popfirst!(joint_buffer)

                    # Add lagged torques
                    push!(control_buffer, actuator_torques)
                    lagged_control = popfirst!(control_buffer)

                    last_buffer_update_t = t
                end

                # Publish simulation data LCM message. Sim data does not include any latency effects
                if t - last_sim_data_t > SIM_DATA_DT
					# ground truth states
	               x[1:3] .= s.d.qpos[1:3]
	               x[4:7] .= s.d.qpos[4:7]
	               x[8:10] .= s.d.qvel[1:3]
	               x[11:13] .= s.d.qvel[4:6]

	               x_true[1:3] .= s.d.qpos[1:3]
	               x_true[4:6] .= s.d.qpos[5:7]
	               x_true[7:9] .= s.d.qvel[1:3]
	               x_true[10:12] .= QuatToRotMatrix(s.d.qpos[4:7])'*s.d.qvel[4:6]

	               joint_pos   .= s.d.sensordata[7:18]
	               joint_vel   .= s.d.sensordata[19:30]

				   MPCControl.control!(actuator_torques, x_true, t, joint_pos, joint_vel, param)
                   s.d.ctrl .= actuator_torques
                end

                # # Publish mocap data at rate determined by MOCAP_DT. Includes latency
                # if (t - last_mocap_t) > MOCAP_DT
				#
                # end

                s.d.ctrl .= lagged_control
                mj_step(s.m, s.d)

                # Break on reset
                (d.d[].time < startsimtm) && break
            end
        end

        render(s, s.window)
        GLFW.PollEvents()
    end

    GLFW.DestroyWindow(s.window)
end
