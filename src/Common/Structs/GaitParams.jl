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
																	0.5	)

function createTrotGait(;stance_time=0.6, swing_time=0.2)
	num_phases = 4
	contact_phases = [	1 1 1 0;
						1 0 1 1;
						1 0 1 1;
						1 1 1 0	]
	phase_times = [stance_time, swing_time, stance_time, swing_time]

	return GaitParams(num_phases, contact_phases, phase_times)
end

function createStandingGait()
	return GaitParams(1, [1;1;1;1], [1.0])
end

function createBoundGait(;stance_time=0.2, flight_time=0.1)
	num_phases = 2
	contact_phases = [	1 0;
						1 0;
						1 0;
						1 0	]
	phase_times = [stance_time, flight_time]

	return GaitParams(num_phases, contact_phases, phase_times)
end
