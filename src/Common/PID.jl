"""
PID around zero positions -> maybe should be feedforward?
"""
function PID!(
    torques::Vector{T},
    q::Vector{T},
    qdot::Vector{T},
    torque_limit::T = 0.0,
    kp::T = 0.4,
    kd::T = 0.01,
) where {T<:Real}
    torques .= kp * (0 .- q) + kd * (0 .- qdot)
    clamp!(torques, -torque_limit, torque_limit)
end

function PID(q::SVector, qdot::SVector; torque_limit, kp, kd)
    unsaturated = kp * (0 .- q) + kd * (0 .- qdot)
    return clamp.(unsaturated, -torque_limit, torque_limit)
end


"""
PID around zero positions -> maybe should be feedforward?
"""
function PID!(
    torques::Vector{T},
    q::Vector{T},
    qdot::Vector{T},
    torque_limit::Vector{T},
    kp::T = 0.4,
    kd::T = 0.01,
) where {T<:Real}
    torques .= kp * (0 .- q) + kd * (0 .- qdot)
    torques .= clamp.(torques, -torque_limit, torque_limit)
end
