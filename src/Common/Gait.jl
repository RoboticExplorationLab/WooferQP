struct GaitParams
	# add in a cyclic array for the phases?
	num_phases::Int64

	# 4xnum_phase array of contacts for each gait phase
	contact_phases::Array{Int64}
	phase_times::Vector{Float64}

	phase_length::Float64
	alpha::Float64
	beta::Float64
end

GaitParams(num_phases, contact_phases, phase_times) = GaitParams(	num_phases,
																	contact_phases,
																	phase_times,
																	sum(phase_times),
																	0.5,
																	0.0	)

function createTrotGait()
	num_phases = 4
	contact_phases = [	1 1 1 0;
						1 0 1 1;
						1 0 1 1;
						1 1 1 0	]
	phase_times = [0.6, 0.2, 0.6, 0.2]

	return GaitParams(num_phases, contact_phases, phase_times)
end

function createStandingGait()
	return GaitParams(1, [1;1;1;1], [1.0])
end

function getPhase(t::AbstractFloat, gait_params::GaitParams)
	phase_time = t % gait_params.phase_length

	for i in 1:gait_params.num_phases
		if phase_time < sum(gait_params.phase_times[1:i])
			return i
		end
	end
end

function nextPhase(phase::Integer, gait_params::GaitParams)
	if (phase == gait_params.num_phases)
		return 1
	else
		return phase+1
	end
end

function coordinateExpander!(expanded::Vector, compact::Vector)
	expanded[1:3] .= compact[1]
	expanded[4:6] .= compact[2]
	expanded[7:9] .= compact[3]
	expanded[10:12] .= compact[4]
end
