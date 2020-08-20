struct SafetyException <: Exception end

"""
Returns false if any of the torques are above the max_torque value. Max_torque can either be a scalar or a vector.
TODO: Verify whether the union is the fastest way to do this.
"""
function ValidActuatorTorques(torques::AbstractVector, max_torque)
    for i=1:length(torques)
        if torques[i] < -max_torque || torques[i] > max_torque
            println("FAILSAFE. Unsafe torque for actuators # ", i)
            println("All torques: ", torques)
            throw(SafetyException(torques, "Torques out of range."))
        end
    end
end

"""
Returns false if any of the actuator positions are above their specified maximum. 
max_positions is a vector.
"""
function ValidActuatorPositions(q::AbstractVector, max_position)
    for i=1:length(q)
        if q[i] < -max_position || q[i] > max_position
            println("FAILSAFE. Unsafe positoin for actuators # ", i)
            println("All positions: ", q)
            throw(SafetyException(q, "Positions out of range."))
        end
    end
end

"""
Returns false if any of the actuator angular rates are above a maximum magnitude. Otherwise returns true.
max_rates is a vector.
"""
function ValidActuatorRates(qdot::AbstractVector, max_rate)
    for i=1:length(qdot)
        if qdot[i] < -max_rate || qdot[i] > max_rate
            println("FAILSAFE. Unsafe ang velocities for actuators # ", i)
            println("All ang velocities: ", qdot)
            throw(SafetyException(qdot, "Ang velocities out of range."))
        end
    end
end