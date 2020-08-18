struct Quadruped{T,S} <: RD.DiscreteLTV
    dt::T
    x_ref::Vector{SVector{12,T}}
    u_ref::Vector{SVector{12,T}}
    foot_locs::Vector{SVector{12,T}}
    contacts::Vector{SVector{4,T}}
    A::Vector{SMatrix{12,12,T,144}}
    B::Vector{SMatrix{12,12,T,144}}
    d::Vector{SVector{12,T}}
    times::Vector{T}
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

        tf = (N - 1) * dt
        times = collect(0.0:dt:(tf))

        new{T,S}(
            dt,
            x_ref,
            u_ref,
            foot_locs,
            contacts,
            A,
            B,
            d,
            times,
            woofer.inertial.body_inertia,
            woofer.inertial.sprung_mass,
        )
    end
end

RD.is_affine(model::Quadruped) = Val(true)

RD.state_dim(::Quadruped) = 12
RD.control_dim(::Quadruped) = 12

RD.get_A(model::Quadruped, k::Integer) = model.A[k]
RD.get_B(model::Quadruped, k::Integer) = model.B[k]
RD.get_d(model::Quadruped, k::Integer) = model.d[k]
RD.get_times(model::Quadruped) = model.times