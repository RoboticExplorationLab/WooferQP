struct Quadruped{T,S} <: TO.RobotDynamics.AbstractModel
    dt::T
    x_ref::Vector{SVector{12,T}}
    u_ref::Vector{SVector{12,T}}
    foot_locs::Vector{SVector{12,T}}
    contacts::Vector{SVector{4,T}}
    J::SMatrix{3,3,T,9}
    sprung_mass::T

    function Quadruped(dt::T, N::S) where {T<:Real,S<:Integer}
        x_ref = [@SVector zeros(T, 12) for i = 1:(N+1)]
        u_ref = [@SVector zeros(T, 12) for i = 1:(N+1)]
        foot_locs = [@SVector zeros(T, 12) for i = 1:(N+1)]
        contacts = [@SVector zeros(T, 4) for i = 1:(N+1)]

        new{T,S}(
            dt,
            x_ref,
            u_ref,
            foot_locs,
            contacts,
            woofer.inertial.body_inertia,
            woofer.inertial.sprung_mass,
        )
    end
end

TO.RobotDynamics.state_dim(::Quadruped) = 12
TO.RobotDynamics.control_dim(::Quadruped) = 12

function NonLinearContinuousDynamics(
    x::SVector,
    u::SVector,
    r::SVector,
    contacts::SVector,
    J::SMatrix,
    sprung_mass::AbstractFloat,
)
    rot = MRP(x[4], x[5], x[6])

    v = @SVector [x[7], x[8], x[9]]
    ω = @SVector [x[10], x[11], x[12]]

    x_d_1_3 = v
    x_d_4_6 = Rotations.kinematics(rot, ω)

    torque_sum = @SVector zeros(3)
    force_sum = @SVector [0, 0, -9.81]
    for i = 1:4
        force_sum += 1 / sprung_mass * contacts[i] * u[SLegIndexToRange(i)]
        torque_sum +=
            contacts[i] *
            Rotations.skew(r[SLegIndexToRange(i)]) *
            rot' *
            u[SLegIndexToRange(i)]
    end
    x_d_7_9 = force_sum
    x_d_10_12 = inv(J) * (-Rotations.skew(ω) * J * ω + torque_sum)

    return [x_d_1_3; x_d_4_6; x_d_7_9; x_d_10_12]
end

function LinearizedContinuousDynamicsA(
    x::SVector{n,T},
    u::SVector{m,T},
    r,
    contacts,
    J,
    sprung_mass,
)::SMatrix{n,n,T,n * n} where {T,n,m}
    return ForwardDiff.jacobian(
        (x_var) ->
            NonLinearContinuousDynamics(x_var, u, r, contacts, J, sprung_mass),
        x,
    )
end

function LinearizedContinuousDynamicsB(
    x::SVector{n,T},
    u::SVector{m,T},
    r,
    contacts,
    J,
    sprung_mass,
)::SMatrix{n,m,T,n * m} where {T,n,m}
    return ForwardDiff.jacobian(
        (u_var) ->
            NonLinearContinuousDynamics(x, u_var, r, contacts, J, sprung_mass),
        u,
    )
end

function TO.RobotDynamics.dynamics(model::Quadruped, x, u, t)
    k = trunc(Integer, t / model.dt) + 1

    A = LinearizedContinuousDynamicsA(
        model.x_ref[k],
        model.u_ref[k],
        model.foot_locs[k],
        model.contacts[k],
        model.J,
        model.sprung_mass,
    )
    B = LinearizedContinuousDynamicsB(
        model.x_ref[k],
        model.u_ref[k],
        model.foot_locs[k],
        model.contacts[k],
        model.J,
        model.sprung_mass,
    )
    d = NonLinearContinuousDynamics(
        model.x_ref[k],
        model.u_ref[k],
        model.foot_locs[k],
        model.contacts[k],
        model.J,
        model.sprung_mass,
    )
    x_dot = A * (x - model.x_ref[k]) + B * (u - model.u_ref[k]) + d

    return x_dot
end

function TO.RobotDynamics.discrete_dynamics(
    ::Type{Euler},
    model::TO.RobotDynamics.AbstractModel,
    x::StaticVector,
    u::StaticVector,
    t,
    dt,
)
    k1 = TO.RobotDynamics.dynamics(model, x, u, t) * dt
    x + k1
end

function discrete_jacobian!(
    ::Type{Q},
    ∇f,
    model::AbstractModel,
    z::AbstractKnotPoint{T,N,M},
) where {T,N,M,Q<:Explicit}
    ix, iu, idt = z._x, z._u, N + M + 1
    t = z.t
    fd_aug(s) = discrete_dynamics(Q, model, s[ix], s[iu], t, z.dt)
    ∇f .= ForwardDiff.jacobian(fd_aug, SVector{N + M}(z.z))
    return nothing
end
