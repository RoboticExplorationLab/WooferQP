using LinearAlgebra

function SkewSymmetricMatrix(q::Vector{T}) where {T<:Real}
    return [T(0) -q[3] q[2]; q[3] T(0) -q[1]; -q[2] q[1] T(0)]
end

function L_q(q::Vector{T}) where {T<:Real}
    return [q[1] -q[2:4]'; q[2:4] q[1] * I + SkewSymmetricMatrix(q[2:4])]::Matrix{T}
end

function R_q(q::Vector{T}) where {T<:Real}
    return [q[1] -q[2:4]'; q[2:4] q[1] * I - SkewSymmetricMatrix(q[2:4])]::Matrix{T}
end

function QuatToRotMatrix(q::Vector{T}) where {T<:Real}
    V = [0 0 0; 1 0 0; 0 1 0; 0 0 1]
    return V' * L_q(q) * R_q(QuatConjugate(q)) * V
end

function QuatConjugate(q::Vector{T}) where {T<:Real}
    return [q[1], -q[2], -q[3], -q[4]]
end

function RotateVectorByQuat!(v::Vector{T}, q::Vector{T}) where {T<:Real}
    v .= (L_q(q) * R_q(QuatConjugate(q)) * VecToQuat(v))[2:4]
end

function RotateVectorByQuat(v::Vector, q::Vector)
    return (L_q(q) * R_q(QuatConjugate(q)) * VecToQuat(v))[2:4]
end

function ThreeParamToQuat(ϕ::Vector{T}) where {T<:Real}
    scalar_squared = T(1.0) - dot(ϕ, ϕ)
    scalar = sign(scalar_squared) * sqrt(abs(scalar_squared))
    return [scalar; ϕ]
end

function QuatToThreeParam(ϕ::Vector{T}) where {T<:Real}
    return ϕ[2:4]
end

function VecToQuat(w::Vector{T}) where {T<:Real}
    return [T(0); w]
end

function V()
    return [    0 1 0 0;
                0 0 1 0;
                0 0 0 1    ]
end
