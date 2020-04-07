mutable struct SimData <: LCMType
    position::SVector{3, Float64}
    orientation::SVector{4, Float64}
    linear_velocity::SVector{3, Float64}
    angular_velocity::SVector{3, Float64} # body frame

    timestamp::Float64
    run_loop::Bool

    # sensor measurements include noise based on sensornoise flag in mujoco XML
    # these also include latency as set by sim
    joint_position::SVector{12, Float64}
    joint_velocity::SVector{12, Float64}
    gyro::SVector{3, Float64}
    accelerometer::SVector{3, Float64}
end

@lcmtypesetup(SimData)

mutable struct ODriveData <: LCMType
    timestamp::Float64
    joint_position::SVector{12, Float64}
    joint_velocity::SVector{12, Float64}
end

@lcmtypesetup(ODriveData)

mutable struct MocapData <: LCMType
    timestamp::Float64
    position::SVector{3, Float64} # lagged, noisy version of position
    orientation::SVector{4, Float64}; # lagged, noisy version of orientation
end

@lcmtypesetup(MocapData)

mutable struct TorqueCommand <: LCMType
    timestamp::Float64
    torque::SVector{12, Float64}
end

@lcmtypesetup(TorqueCommand)
