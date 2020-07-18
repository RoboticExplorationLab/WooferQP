struct GaitParams{T,S}
    # add in a cyclic array for the phases?
    num_phases::S

    # 4xnum_phase array of contacts for each gait phase
    contact_phases::Matrix{S}
    phase_times::Vector{T}

    phase_length::T
    alpha::T
    beta::T
end

function GaitParams(
    num_phases::S,
    contact_phases::Matrix{S},
    phase_times::Vector{T},
) where {T<:Number,S<:Integer}
    @assert num_phases == size(contact_phases, 2) == length(phase_times)
    GaitParams{T,S}(
        num_phases,
        contact_phases,
        phase_times,
        sum(phase_times),
        convert(T, 0.5),
        convert(T, 0.5),
    )
end

function createTrotGait(; stance_time = 0.6, swing_time = 0.2)
    num_phases = 4
    contact_phases = [
        1 1 1 0
        1 0 1 1
        1 0 1 1
        1 1 1 0
    ]
    phase_times = [stance_time, swing_time, stance_time, swing_time]

    return GaitParams(num_phases, contact_phases, phase_times)
end

function createStandingGait()
    return GaitParams(2, [1 1; 1 1; 1 1; 1 1], [1.0, 1.0])
end

function createPronkGait(; stance_time = 0.2, flight_time = 0.1)
    num_phases = 2
    contact_phases = [
        1 0
        1 0
        1 0
        1 0
    ]
    phase_times = [stance_time, flight_time]

    return GaitParams(num_phases, contact_phases, phase_times)
end

function createPaceGait(; stance_time = 0.6, swing_time = 0.2)
    num_phases = 4
    contact_phases = [
        1 1 1 0
        1 0 1 1
        1 1 1 0
        1 0 1 1
    ]
    phase_times = [stance_time, swing_time, stance_time, swing_time]

    return GaitParams(num_phases, contact_phases, phase_times)
end

function createBoundGait(; front_time = 0.2, back_time = 0.2, stance_time = 0.1)
    num_phases = 4
    contact_phases = [
        1 1 1 0
        1 1 1 0
        1 0 1 1
        1 0 1 1
    ]
    phase_times = [stance_time, front_time, stance_time, back_time]

    return GaitParams(num_phases, contact_phases, phase_times)
end

function createFlyingTrot(; stance_time=0.2, flight_time=0.1)
    num_phases = 4
    contact_phases = [
        1 0 0 0
        0 0 1 0
        0 0 1 0
        1 0 0 0
    ]
    phase_times = [stance_time, flight_time, stance_time, flight_time]

    return GaitParams(num_phases, contact_phases, phase_times)
end
