using Parameters
using LinearAlgebra
using StaticArrays

@with_kw struct WooferConfig
    FRAME_MASS::Float64 = 3.0
    MODULE_MASS::Float64 = 1.033
    THIGH_MASS::Float64 = 0.070
    SHIN_MASS::Float64 = 0.059
    TOTAL_LEG_MASS::Float64 = (THIGH_MASS + SHIN_MASS) * 2

    TOTAL_MASS::Float64 = FRAME_MASS + 4 * MODULE_MASS + 4 * TOTAL_LEG_MASS
    SPRUNG_MASS::Float64 = FRAME_MASS + 4 * MODULE_MASS + THIGH_MASS * 8

    # Robot joint limits
    MAX_JOINT_TORQUE::Float64 = 12
    MAX_LEG_FORCE::Float64 = 133
    REVOLUTE_RANGE::Float64 = 3
    PRISMATIC_RANGE::Float64 = 0.18

    # Robot geometry
    LEG_FB::Float64 = 0.23# front-back distance from center line to leg axis
    LEG_LR::Float64 = 0.109 # left-right distance from center line to leg plane
    ABDUCTION_OFFSET::Float64 = 0.064 # distance from abduction axis to leg
    FOOT_RADIUS::Float64 = 0.02

    HIP_LAYOUT::SMatrix{4, 3, Float64, 12} = @SMatrix  [LEG_FB -LEG_LR 0;
                        LEG_FB  LEG_LR 0;
                        -LEG_FB -LEG_LR 0;
                        -LEG_FB  LEG_LR 0]
    ABDUCTION_LAYOUT::SVector{4, Float64} = @SVector [-ABDUCTION_OFFSET, ABDUCTION_OFFSET, -ABDUCTION_OFFSET, ABDUCTION_OFFSET]

    L::Float64 = 0.66
    W::Float64 = 0.176
    T::Float64 = 0.092
    Ix::Float64 = SPRUNG_MASS / 12 * (W^2 + T^2)
    Iy::Float64 = SPRUNG_MASS / 4 * (L^2 + T^2)
    Iz::Float64 = SPRUNG_MASS / 4 * (L^2 + W^2)
    INERTIA::SMatrix{3, 3, Float64, 9} = Diagonal(@SVector [Ix, Iy, Iz])
    l0::Float64 = 0.18
    l1::Float64 = 0.32
end

WOOFER_CONFIG = WooferConfig()

@with_kw struct EnvironmentConfig
    MU::Float64 = 1.5 # coeff friction
    SIM_STEPS::Int = 10000 # simulation steps to take
    DT::Float64 = 0.0001 # timestep [s]
end

ENVIRONMENT_CONFIG = EnvironmentConfig()
