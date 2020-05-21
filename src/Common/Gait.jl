function getPhase(t::AbstractFloat, param::ControllerParams)
	phase_time = t % param.gait.phase_length

	for i in 1:param.gait.num_phases
		if phase_time < sum(param.gait.phase_times[1:i])
			return i
		end
	end
end

function nextPhase(phase::Integer, param::ControllerParams)
	if (phase == param.gait.num_phases)
		return 1
	else
		return phase+1
	end
end

# function nextFootTimeOnGround(i, phase::Integer, gait_params::GaitParams)
# """
# Finds the next length of time for which foot i will be on the ground
# """
# 	# TODO
# end

function coordinateExpander!(expanded::Vector, compact::Vector)
	expanded[1:3] .= compact[1]
	expanded[4:6] .= compact[2]
	expanded[7:9] .= compact[3]
	expanded[10:12] .= compact[4]
end
