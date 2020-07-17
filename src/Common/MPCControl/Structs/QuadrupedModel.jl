abstract type DiscreteQuadruped <: TO.RobotDynamics.Explicit end

struct Quadruped{T,S} <: TO.RobotDynamics.AbstractModel
    dt::T
    x_ref::Vector{SVector{12,T}}
    u_ref::Vector{SVector{12,T}}
    foot_locs::Vector{SVector{12,T}}
    contacts::Vector{SVector{4,T}}
    A::Vector{SMatrix{12,12,T,144}}
    B::Vector{SMatrix{12,12,T,144}}
    d::Vector{SVector{12,T}}
    J::SMatrix{3,3,T,9}
    sprung_mass::T

    function Quadruped(dt::T, N::S) where {T<:Real,S<:Integer}
        x_ref = [@SVector zeros(T, 12) for i = 1:(N+1)]
        u_ref = [@SVector zeros(T, 12) for i = 1:(N+1)]
        foot_locs = [@SVector zeros(T, 12) for i = 1:(N+1)]
        contacts = [@SVector zeros(T, 4) for i = 1:(N+1)]
        A = [@SMatrix zeros(T, 12, 12) for _ = 1:N+1]
        B = [@SMatrix zeros(T, 12, 12) for _ = 1:N+1]
        d = [@SVector zeros(T, 12) for _ = 1:N+1]

        new{T,S}(
            dt,
            x_ref,
            u_ref,
            foot_locs,
            contacts,
            A,
            B,
            d,
            woofer.inertial.body_inertia,
            woofer.inertial.sprung_mass,
        )
    end
end

TO.RobotDynamics.state_dim(::Quadruped) = 12
TO.RobotDynamics.control_dim(::Quadruped) = 12

function TO.RobotDynamics.discrete_dynamics(
    ::Type{DiscreteQuadruped},
    model::Quadruped,
    x::StaticVector,
    u::StaticVector,
    t,
    dt
)
    k = trunc(Integer, t / model.dt) + 1

    x = model.A[k] * (x - model.x_ref[k]) + model.B[k] * (u - model.u_ref[k]) + model.d[k]

    return x
end

function TO.RobotDynamics.discrete_jacobian!(
    ::Type{DiscreteQuadruped},
    ∇f,
    model::Quadruped,
    z::TO.AbstractKnotPoint{T,N,M},
) where {T,N,M,Q<:TO.RobotDynamics.Explicit}
    ix, iu, idt = z._x, z._u, N + M + 1
    t = z.t
    k = trunc(Integer, t / model.dt) + 1
    ∇f .= [model.A[k] model.B[k]]
    return nothing
end
