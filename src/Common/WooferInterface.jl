using StaticArrays
using Parameters

include("LQRController.jl")
include("StateEstimator.jl")
include("Safety.jl")
include("../Hardware/HardwareInterface.jl")

"""
Collection of parameters and cache variables
"""
@with_kw mutable struct WooferData
    joint_pos::Vector{Float64} = zeros(12)
    joint_vel::Vector{Float64} = zeros(12)

    mocap_r::Vector{Float64} = zeros(3)
    mocap_q::Vector{Float64} = [1, 0, 0, 0]
    mocap_timestamp::Float64 = 0.0
    mocap_prev_timestamp::Float64 = 0.0

    mocap_start_r::Vector{Float64} = zeros(3)

    dt::Float64
    lqr_params::LQRParams
    est_params::StateEstimatorParams

    odrive_hardware::ODriveHardware

    lqr_forces::SVector{12,Float64} = @SVector zeros(12)
    actuator_torques::SVector{12,Float64} = @SVector zeros(12)

    local_ip::String
    server_ip::String
    body_id::Int
end

"""
Sends torques to ODrives and then waits for the responses using the supplied WooferData.
Throws SafetyException is joint state exceeds the safety boundaries.
Blocking.
"""
function ActuatorTask(wd::WooferData)
    SetActuatorTorques(wd.actuator_torques, wd.odrive_hardware)
    WaitForEncoderMeasurements(wd.odrive_hardware)
    ReadActuatorStates!(wd.joint_pos, wd.joint_vel, wd.odrive_hardware)
end

"""
Run safety checks on the given WooferData. If safety is violated, throws a SafetyException.
"""
function RunSafetyChecks(wd::WooferData; max_torque::AbstractFloat, max_position::AbstractFloat, max_vel::AbstractFloat)
    # TODO get rid of this allocation
    ValidActuatorTorques(wd.actuator_torques, max_torque) # max 1Nm torque
    ValidActuatorPositions(wd.joint_pos, max_position) # max 1.0 rad
    ValidActuatorRates(wd.joint_vel, max_vel) # max 10 rad/s
end

"""
Uses LQR to find the stabilizing torques. Modifies WooferData in-place.
"""
function CalculateLQRTorques!(wd::WooferData)
    wd.lqr_forces = LQRBalance(wd.est_params.x, wd.lqr_params)
    q = reconstructq(wd.est_params.x[@SVector [4, 5, 6]])
    wd.actuator_torques = InertialLegForcesToActuatorTorques(
        wd.lqr_forces,
        q,
        SVector{12}(wd.joint_pos),
    )
end

"""
Go through one iteration of closed loop control
"""
function ClosedLoopIteration(wd::WooferData, mocap_client; max_torque = 2.5, max_position = 0.25, use_mocap = true, use_pid = true)
    propogation_time = @elapsed StateEstimatorDynamicsPropogation(
        wd.dt,
        SVector{12}(wd.joint_pos),
        SVector{12}(wd.joint_vel),
        wd.est_params,
        wd.lqr_forces,
    )

    mocap_time = 0.0
    read_mocap_time = 0.0

    if use_mocap
        read_mocap_time = @elapsed mocapdata = mocap_client.getData()
        wd.mocap_timestamp = mocapdata[8]
        if wd.mocap_prev_timestamp != wd.mocap_timestamp
            wd.mocap_r = mocapdata[1:3] - wd.mocap_start_r
            wd.mocap_q = mocapdata[4:7]            
            mocap_time = @elapsed StateEstimatorMocapUpdate(
                SVector{3}(wd.mocap_r),
                UnitQuaternion(wd.mocap_q...),
                wd.est_params,
            )
            wd.mocap_prev_timestamp = wd.mocap_timestamp
        end
    end

    lqr_time = @elapsed CalculateLQRTorques!(wd)
    if use_pid
        # TODO: remove this for when we do real robot stuff
        # Run PID with torque limit 2.0Nm, kp of 1.2
        # Overrides the lqr torques calculated above
        pid_time = @elapsed wd.actuator_torques = PID(
            SVector{12}(wd.joint_pos),
            SVector{12}(wd.joint_vel),
            torque_limit = max_torque,
            kp = 6.0,
            kd = 0.02,
        )
        @show wd.actuator_torques
    end

    RunSafetyChecks(wd, max_torque = max_torque, max_position = max_position, max_vel = 10.0)
    actuator_time = @elapsed ActuatorTask(wd)

    return (propogation_time, read_mocap_time, mocap_time, lqr_time, pid_time, actuator_time)
end
