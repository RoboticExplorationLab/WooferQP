module Woofer

using Parameters
using LinearAlgebra
using StaticArrays
using PyCall
using Revise
using DelimitedFiles

include("Config.jl")
include("Dynamics.jl")
include("LQRController.jl")
include("StateEstimator.jl")
include("../Hardware/HardwareInterface.jl")
include("../Hardware/MoCap.jl")
include("PID.jl")
include("Safety.jl")
include("WooferInterface.jl")


@with_kw struct DataLog
    n::Int
    timing::Matrix{Float64} = zeros(6, n)
    mocap_r::Matrix{Float64} = zeros(3, n)
    mocap_q::Matrix{Float64} = zeros(4, n)
    state_estimate::Matrix{Float64} = zeros(12, n)
    lqr_forces::Matrix{Float64} = zeros(12, n)
    actuator_torques::Matrix{Float64} = zeros(12, n)
end

function RunClosedLoopControl(wd, mocap_client, data_log; n, use_mocap, use_pid, max_torque, max_position)
    println("Running closed loop control for ", n, " steps.")
    for i = 1:n
        (propogation_time, read_mocap_time, mocap_time, lqr_time, pid_time, actuator_time) = ClosedLoopIteration(wd, mocap_client, max_torque=max_torque, max_position=max_position, use_mocap=use_mocap, use_pid=use_pid)
        
        data_log.timing[1, i] = propogation_time
        data_log.timing[2, i] = read_mocap_time
        data_log.timing[3, i] = mocap_time
        data_log.timing[4, i] = lqr_time
        data_log.timing[5, i] = pid_time
        data_log.timing[6, i] = actuator_time
        data_log.mocap_r[:, i] .= wd.mocap_r
        data_log.mocap_q[:, i] .= wd.mocap_q
        data_log.state_estimate[:, i] .= wd.est_params.x
        data_log.lqr_forces[:, i] .= wd.lqr_forces
        data_log.actuator_torques[:, i] .= wd.actuator_torques
    end
end

function run(; do_calibration, max_torque, max_position, use_mocap=true, calibrate_sequentially=false, closed_loop_iters=4000)
    """
    Constants - to be checked
    """

    dt = 0.001
    wd = WooferData(
        dt = dt,
        lqr_params = initLQRParams(
            dt = dt,
            x0 = SVector{12,Float64}([0.0, 0.0, 0.32, zeros(9)...]),
            rxy_penalty = 1e2,
            rz_penalty = 1e2,
            phi_penalty = 1e2,
            v_penalty = 1,
            w_penalty = 1,
            fxy_penalty = 1e-3,
            fz_penalty = 1e-3,
        ),
        est_params = StateEstimatorParams(
            dt = dt,
            x = SVector{12}([0.0, 0.0, 0.30, zeros(9)...]),
            P = 1e-3 * Diagonal(@SVector ones(12)),
            Q = 1e-1 * Diagonal(@SVector ones(12)),
            R_mocap = 1e-3 * Diagonal(@SVector ones(6)),
            R_joint = 1e-2 * Diagonal(@SVector ones(12)),
        ),
        odrive_hardware = CreateODriveHardware(
            serial_numbers = [
                "35619029856331",
                "35799416713291",
                "59752448340023",
                "35726434842440",
                "35799416778827",
                "35795121942603",
            ],
            odrive_motor_mapping = [(2, 3), (4, 1), (6, 5), (12, 11), (10, 7), (8, 9)],
            direction_corrections = [1, -1, 1, 1, -1, 1, -1, -1, 1, -1, -1, 1],
            counts_per_radian = 2000 / 2pi, # encoder counts per motor radian
            torque_constant = 0.082, # motor Nm per A, should include torque loss from belt drive
            speed_reduction = 4.0, # speed reduction of the belt drive
        ),
        local_ip = "192.168.0.9",
        server_ip = "192.168.0.2",
        body_id = 1,
    )
    data_log = DataLog(n=closed_loop_iters)

    mocap_client = nothing
    try
        if use_mocap
            mocap_client = pyimport("pymocap")
            mocap_client.startDataStream(local_ip = wd.local_ip,server_ip = wd.server_ip, body_id = wd.body_id)
        end

        # Calibrate odrives
        if do_calibration
            println("Press enter to calibrate.")
            readline(stdin)
            CalibrateODrives(wd.odrive_hardware, calibrate_sequentially)
        else
            println("Skipping calibration. Press enter to confirm and continue.")
            readline(stdin)
        end

        println("Waiting to run closed loop code with zero torque. Press enter to continue.")
        readline(stdin)

        SetClosedLoopCurrentControl(wd.odrive_hardware)
        ActuatorTask(wd)
        println("Initial actuator positions: ", wd.joint_pos)

        # Start PID loop and State Estimator
        StateEstimatorDynamicsPropogation(
            wd.dt,
            SVector{12}(wd.joint_pos),
            SVector{12}(wd.joint_vel),
            wd.est_params,
            wd.lqr_forces,
        )

        # Get mocap measurement
        if use_mocap
            wd.mocap_prev_timestamp = 0.0
            mocapdata = mocap_client.getData()
            wd.mocap_timestamp = mocapdata[8]
            if wd.mocap_prev_timestamp != wd.mocap_timestamp
                wd.mocap_start_r = mocapdata[1:3]
                wd.mocap_r = zeros(3)
                wd.mocap_q = mocapdata[4:7]
                StateEstimatorMocapUpdate(SVector{3}(wd.mocap_r), UnitQuaternion(wd.mocap_q...), wd.est_params)
            end
        end

        # Run LQR with torque = 0
        CalculateLQRTorques!(wd)

        # Run PID with torque = 0
        wd.actuator_torques = PID(
            SVector{12}(wd.joint_pos),
            SVector{12}(wd.joint_vel),
            torque_limit = 0.0,
            kp = 0.4,
            kd = 0.01,
        )

        ActuatorTask(wd)
        @show max_position
        RunSafetyChecks(wd, max_torque=0.1, max_position=max_position, max_vel=1.0)

        ClosedLoopIteration(wd, mocap_client, max_torque=0.0, max_position=max_position)

        println("Data from the first loop:")
        @show wd.mocap_r, wd.mocap_q
        @show wd.est_params.x
        println("Waiting to run closed loop controller with activated motors. Press enter to continue.")
        readline(stdin)

        @time closed_loop_time = @elapsed RunClosedLoopControl(wd, mocap_client, data_log, n=closed_loop_iters, use_mocap=use_mocap, use_pid=true, max_position=max_position, max_torque=max_torque)
        combined_data_log = [data_log.timing; data_log.mocap_r; data_log.mocap_q; data_log.state_estimate; data_log.lqr_forces; data_log.actuator_torques]
        open("data_log.txt", "w") do io
            writedlm(io, combined_data_log', ',')
        end

        @show closed_loop_time
    catch e
        @show e
    finally
        SetODrivesIdle(wd.odrive_hardware)
        CloseODriveHardware(wd.odrive_hardware)
        println("Succeeded in closing ODrive hardware.")
    end
end
end
