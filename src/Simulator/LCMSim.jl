using LinearAlgebra
using Rotations
using Parameters
using StaticArrays
using Plots
using LCMCore
using DataStructures

include("XMLParser.jl")
include("../Common/Dynamics.jl")
include("../Common/Config.jl")
include("../Common/LCMTypes.jl")
include("../Common/StateEstimator.jl")
include("../Common/LQRController.jl")
include("../Common/Rotations.jl")

function simulate()
    # Create the robot XML file
    ParseXML()
    # Get model from local directory
    s = loadmodel("woofer_out.xml", 1200, 900)

    # Instantiate Lightweight Communication & Marshalling socket
    lcm = LCM()

    # Start asynchronous handler. This works because the handler is operating on a file descriptor 
    # which uses processor condition variables and event notifications.
    @async while true
        handle(lcm)
    end

    # Pre-allocate memory for time-delayed torques, joint pos/vel, and mocap data
    lagged_control = zeros(12)
    lagged_joint = zeros(24)
    lagged_mocap = zeros(7)

    # Pre-allocate memory for actuator torques
    actuator_torques = zeros(12)

    last_actuator_timestamp = 0.0

    # LCM callback for when torque command is received
    function TorqueCallback(channel::String, received_command::TorqueCommand)
        actuator_torques = Vector{Float64}(received_command.torque)
        last_actuator_timestamp = received_command.timestamp
	end
	subscribe(lcm, "TORQUE_COMMAND", TorqueCallback, TorqueCommand)

    # Mujoco data and model
    d = s.d
    m = s.m

    # Initialize update rates 
    ODRIVE_DT = 0.001
    SIM_DATA_DT = 0.01
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

    # Initial LCM packet
	odrive_message = ODriveData()
    mocap_message = MocapData()
    sim_message = SimData()
    sim_message.gyro = @SVector zeros(3)
	sim_message.accelerometer = @SVector zeros(3)
    sim_message.run_loop = true

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

                # Publish fake LCM ODriveData message that includes latency
                if (t - last_t) > ODRIVE_DT
					odrive_message.timestamp = t - joint_latency_frames
					odrive_message.joint_position = SVector{12}(lagged_joint[1:12])
					odrive_message.joint_velocity = SVector{12}(lagged_joint[13:24])
                    publish(lcm, "ODRIVE_DATA", encode(odrive_message))

                    last_t = t
                end

                # Publish simulation data LCM message. Sim data does not include any latency effects
                if t - last_sim_data_t > SIM_DATA_DT
                    sim_message.position = SVector{3}(s.d.qpos[1:3])
                    sim_message.orientation = SVector{4}(s.d.qpos[4], s.d.qpos[5], s.d.qpos[6], s.d.qpos[7])
                    sim_message.linear_velocity = SVector{3}(s.d.qvel[1:3])
                    robot_quat_orientation = UnitQuaternion(s.d.qpos[4], s.d.qpos[5], s.d.qpos[6], s.d.qpos[7])
                    sim_message.angular_velocity = SVector{3}(rotmat(inv(robot_quat_orientation)) * SVector{3}(s.d.qvel[4:6]))
                    sim_message.timestamp = t

                    # sensor measurements include noise based on sensornoise flag in mujoco XML
                    # these also include latency as set by sim
                    sim_message.joint_position = SVector{12}(s.d.sensordata[7:18]...)
                    sim_message.joint_velocity = SVector{12}(s.d.sensordata[19:30]...)
                    publish(lcm, "SIM_DATA", encode(sim_message))

                    last_sim_data_t = t

                    # println(actuator_torques)
                end

                # Publish mocap data at rate determined by MOCAP_DT. Includes latency
                if (t - last_mocap_t) > MOCAP_DT
                    # Send motion capture lcm data
                    mocap_message.timestamp = t
                    mocap_message.position = SVector{3}(lagged_mocap[1:3])
                    mocap_message.orientation = SVector{4}(lagged_mocap[4:7])
                    publish(lcm, "MOCAP_DATA", encode(mocap_message))
                    
                    last_mocap_t = t
                    println("Sim time: ", t)
                    flush(stdout)
                end

                s.d.ctrl .= lagged_control
                mj_step(s.m, s.d)
                
                # Break on reset
                (d.d[].time < startsimtm) && break
            end
        end

        render(s, s.window)
        GLFW.PollEvents()
    end

	# Kill control/estimation loop
	sim_message.run_loop = false
	publish(lcm, "SIM_DATA", encode(sim_message))

    GLFW.DestroyWindow(s.window)
end
