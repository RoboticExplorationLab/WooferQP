"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

using ForwardDiff

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
    sprung_mass
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
    sprung_mass
)::SMatrix{n,m,T,n * m} where {T,n,m}
    return ForwardDiff.jacobian(
        (u_var) ->
            NonLinearContinuousDynamics(x, u_var, r, contacts, J, sprung_mass),
        u,
    )
end

function reference_trajectory!(
    x_curr::Vector{T},
    param::ControllerParams,
) where {T<:Number}
    # TODO: integrate the x,y,ψ position from the reference
    α = collect(range(0, 1, length = param.N + 1))
    # interpolate everything but x, y position
    if param.vel_ctrl
        interp_indices = 3:12
    else
        interp_indices = 1:12
    end

    for i = 1:param.N+1
        param.x_ref[interp_indices, i] .=
            (1 - α[i]) * x_curr[interp_indices] +
            α[i] * param.x_des[interp_indices]
    end

    if param.vel_ctrl
        # integrate x, y position
        integ_indices = 1:2
        param.x_ref[integ_indices, :] .=
            repeat(x_curr[integ_indices, 1], 1, param.N + 1) +
            cumsum(param.x_ref[7:8, :] * param.optimizer.dt, dims = 2)
    end
end

function foot_forces!(
    x_curr::AbstractVector{T},
    param::ControllerParams,
) where {T<:Number}
    # x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
    # contacts: 4xN+1 matrix of foot contacts over the planning horizon
    # foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon

    opt = param.optimizer
    N = opt.N
    dt = opt.dt

    for i = 1:N
        opt.x_ref[i][] = SVector{12}(param.x_ref[:, i])

		foot_locs_static = SVector{12}(param.foot_locs[:, i])
		contacts_static = SVector{4}(param.contacts[:, i])

        A_c_i = LinearizedContinuousDynamicsA(opt.x_ref[i][], opt.u_ref[i][], foot_locs_static, contacts_static, opt.J, opt.sprung_mass)
		B_c_i = LinearizedContinuousDynamicsB(opt.x_ref[i][], opt.u_ref[i][], foot_locs_static, contacts_static, opt.J, opt.sprung_mass)
		d_c_i = NonLinearContinuousDynamics(opt.x_ref[i][], opt.u_ref[i][], foot_locs_static, contacts_static, opt.J, opt.sprung_mass)

		# Euler Discretization
		opt.A_d[i][] = oneunit(SMatrix{12,12,T}) + A_c_i*dt
		opt.B_d[i][] = B_c_i*dt
		opt.d_d[i][] = d_c_i*dt + opt.x_ref[i][]

		opt.q[i][] = 2*opt.Q_i*opt.x_ref[i][]
		opt.r[i][] = 2*opt.R_i*opt.u_ref[i][]
    end
    opt.x_ref[N+1][] = SVector{12}(param.x_ref[:, N+1])
	opt.q[N+1][] = 2*opt.Q_f*opt.x_ref[N+1][]

    @time solve!(opt.model)

    param.forces .=
        value.(opt.model, opt.u)[select12(1)]
end

function select12(i)
    return (12*(i-1)+1):(12*(i-1)+12)
end

function select12_3(i, j, k)
    return 12*(i-1) + 3*(j-1) + k
end
