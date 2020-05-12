using ForwardDiff
using StaticArrays
include("Config.jl")
include("Utilities.jl")

#########################################################################################################
# Static Array Version
function ForwardKinematics(α::SVector{3}, i::Integer)
    base = ForwardHipRelativeKinematics(α, i)
    return base + (WOOFER_CONFIG::WooferConfig).HIP_LAYOUT[i, :]
end

function ForwardKinematicsAll(α::SVector{12})
    return SVector{12}([ForwardKinematics(α[SLegIndexToRange(1)], 1);
                        ForwardKinematics(α[SLegIndexToRange(2)], 2);
                        ForwardKinematics(α[SLegIndexToRange(3)], 3);
                        ForwardKinematics(α[SLegIndexToRange(4)], 4)])
end

function rotx(theta)
    @SMatrix [1 0 0;
              0 cos(theta) -sin(theta);
              0 sin(theta) cos(theta)]
end

function ForwardHipRelativeKinematics(α::SVector{3}, i::Integer)
    γ = 0.5 * (α[3] - α[2]) + 0.5π
    θ = - 0.5 * (α[2] + α[3])

    d = (WOOFER_CONFIG::WooferConfig).l0 * sin(γ)
    h1 = (WOOFER_CONFIG::WooferConfig).l0 * cos(γ)
    h2 = sqrt((WOOFER_CONFIG::WooferConfig).l1^2 - d^2)
    L = h1 + h2
    unrotated = SVector{3}(L * sin(θ), (WOOFER_CONFIG::WooferConfig).ABDUCTION_LAYOUT[i], -L * cos(θ))
    return rotx(α[1]) * unrotated
end

function LegJacobian(q::SVector{3}, i::Int)
    ForwardDiff.jacobian(q->ForwardKinematics(q, i), q)
end

function Force2Torque(f::SVector{12}, α::SVector{12})
    return [LegJacobian(α[SLegIndexToRange(1)], 1)' * f[SLegIndexToRange(1)];
            LegJacobian(α[SLegIndexToRange(2)], 2)' * f[SLegIndexToRange(2)];
            LegJacobian(α[SLegIndexToRange(3)], 3)' * f[SLegIndexToRange(3)];
            LegJacobian(α[SLegIndexToRange(4)], 4)' * f[SLegIndexToRange(4)]]
end

#########################################################################################################
# Dynamic Array Versions
function ForwardHipRelativeKinematics(α::Vector, i::Integer)
    γ = 0.5 * (α[3] - α[2]) + 0.5π
    θ = - 0.5 * (α[2] + α[3])

    d = (WOOFER_CONFIG::WooferConfig).l0 * sin(γ)
    h1 = (WOOFER_CONFIG::WooferConfig).l0 * cos(γ)
    h2 = sqrt((WOOFER_CONFIG::WooferConfig).l1^2 - d^2)
    L = h1 + h2
    unrotated = [L * sin(θ), (WOOFER_CONFIG::WooferConfig).ABDUCTION_LAYOUT[i], -L * cos(θ)]
    return Vector(rotx(α[1]) * unrotated)
end

function ForwardKinematics(α::Vector, i::Integer)
    return ForwardHipRelativeKinematics(α, i) + Vector((WOOFER_CONFIG::WooferConfig).HIP_LAYOUT[i, :])
end

function ForwardKinematicsAll(α::Vector)
    return [ForwardKinematics(α[LegIndexToRange(1)], 1);
            ForwardKinematics(α[LegIndexToRange(2)], 2);
            ForwardKinematics(α[LegIndexToRange(3)], 3);
            ForwardKinematics(α[LegIndexToRange(4)], 4)]
end

function Force2Torque(f::Vector, α::Vector)
    return [LegJacobian(α[LegIndexToRange(1)], 1)' * f[LegIndexToRange(1)];
            LegJacobian(α[LegIndexToRange(2)], 2)' * f[LegIndexToRange(2)];
            LegJacobian(α[LegIndexToRange(3)], 3)' * f[LegIndexToRange(3)];
            LegJacobian(α[LegIndexToRange(4)], 4)' * f[LegIndexToRange(4)]]
end

function LegJacobian(q::Vector, i::Int)
    return ForwardDiff.jacobian((q)->ForwardKinematics(q, i), q)
end

#########################################################################################################



"""
Fills A with skew symmetric version of a
"""
function SkewSymmetricMatrix!(A::Matrix{T}, a::Vector{T}) where {T<:Real}
    A[1, 1] = T(0)
    A[2, 2] = T(0)
    A[3, 3] = T(0)
    A[1, 2] = -a[3]
    A[1, 3] = a[2]
    A[2, 1] = a[3]
    A[2, 3] = -a[1]
    A[3, 1] = -a[2]
    A[3, 2] = a[1]
end
