struct OptimizerParams{T, S}
    # discretization length
    dt::T

    # planning horizon length
    N::S

    # Parametron Variables:
    model::Model
    x::Vector{Variable}
    u::Vector{Variable}

    # Parameter Values (to be updated in place)
    x_ref::Vector{SVector{12, T}}
    u_ref::Vector{SVector{12, T}}
    foot_locs::Vector{SVector{12, T}}
    contacts::Vector{SVector{4, T}}
    Q_f::Diagonal
    Q_i::Diagonal
    R_i::Diagonal

    function OptimizerParams(
        dt::T,
        N::S,
        q::Vector{T},
        r::Vector{T};
        n::Integer = 12,
        m::Integer = 12,
    ) where {T<:Number, S<:Integer}
        # initialize model and variables
        model = Model(OSQP.Optimizer(verbose = false))
        x = [Variable(model) for _ = 1:((N+1)*n)]
        u = [Variable(model) for _ = 1:((N)*n)]

        # initialize quadratic cost parameters
        Q_f = Diagonal(SVector{n}(q))
        Q_i = Diagonal(SVector{n}(q))
        R_i = Diagonal(SVector{m}(r))

        x_ref = [@SVector zeros(n) for _ = 1:(N + 1)]
        u_ref = [@SVector zeros(m) for _ = 1:(N)]
        foot_locs = [@SVector zeros(12) for _ = 1:(N)]
        contacts = [@SVector zeros(4) for _ = 1:(N)]

        new{T, S}(dt, N, model, x, u, x_ref, u_ref, foot_locs, contacts, Q_f, Q_i, R_i)
    end
end

function populate_model!(opt::OptimizerParams, μ, min_vert_force, max_vert_force)
    x_ref_param = [Parameter(()->opt.x_ref[i], opt.model) for i=1:(opt.N+1)]
    u_ref_param = [Parameter(()->opt.u_ref[i], opt.model) for i=1:(opt.N)]
    foot_locs_param = [Parameter(()->opt.foot_locs[i], opt.model) for i=1:(opt.N)]
    contacts_param = [Parameter(()->opt.contacts[i], opt.model) for i=1:(opt.N)]

    add_objective!(opt, x_ref_param, u_ref_param)
    add_control_constraints!(opt, μ, min_vert_force, max_vert_force)
    add_dynamics_constraints!(opt, x_ref_param, u_ref_param, foot_locs_param, contacts_param)
end

function add_objective!(opt::OptimizerParams, x_ref, u_ref)
    N = opt.N
    model = opt.model
    x = opt.x
    u = opt.u

    # terminal state cost
    objective = @expression (x[select12(N+1)] - x_ref[N+1])'*opt.Q_f*(x[select12(N+1)] - x_ref[N+1])

    for i=1:N
        # stagewise state penalty
        objective = @expression objective + (x[select12(i)] - x_ref[i])'*opt.Q_i*(x[select12(i)] - x_ref[i])

        #stagewise control penalty
        objective = @expression objective + (u[select12(i)] - u_ref[i])'*opt.R_i*(u[select12(i)] - u_ref[i])
    end

    @objective(model, Minimize, objective)
end

function add_control_constraints!(opt::OptimizerParams, μ::Number, min_vert_force::Number, max_vert_force::Number)
    model = opt.model
    u = opt.u
    N = opt.N

    for i = 1:N
        # Control constraints
        for j = 1:4
            # convert absolute value constraint to linear inequality:
            @constraint(
                model,
                u[select12_3(i, j, 1)] <= μ * u[select12_3(i, j, 3)]
            )
            @constraint(
                model,
                u[select12_3(i, j, 1)] >= -μ * u[select12_3(i, j, 3)]
            )
            @constraint(
                model,
                u[select12_3(i, j, 2)] <= μ * u[select12_3(i, j, 3)]
            )
            @constraint(
                model,
                u[select12_3(i, j, 2)] >= -μ * u[select12_3(i, j, 3)]
            )

            @constraint(model, u[select12_3(i, j, 3)] >= min_vert_force)
            @constraint(model, u[select12_3(i, j, 3)] <= max_vert_force)
        end
    end
end

function add_dynamics_constraints!(opt::OptimizerParams, x_ref, u_ref, foot_locs, contacts)
    model = opt.model
    x = opt.x
    u = opt.u
    J = woofer.inertial.body_inertia
    sprung_mass = woofer.inertial.sprung_mass
    dt = opt.dt
    N = opt.N

    @constraint(
        model,
        x[select12(1)] == x_ref[1]
    )
    for i=1:(N)
        A_d = Parameter(() -> LinearizedDiscreteDynamicsA(x_ref[i](), u_ref[i](), foot_locs[i](), contacts[i](), J, sprung_mass, dt), model)
        B_d = Parameter(() -> LinearizedDiscreteDynamicsB(x_ref[i](), u_ref[i](), foot_locs[i](), contacts[i](), J, sprung_mass, dt), model)
        d_d = Parameter(() -> NonLinearContinuousDynamics(x_ref[i](), u_ref[i](), foot_locs[i](), contacts[i](), J, sprung_mass)*dt + x_ref[i](), model)

        @constraint(
            model,
            x[select12(i + 1)] ==
            A_d * (x[select12(i)] - x_ref[i]) +
            B_d * (u[select12(i)] - u_ref[i]) +
            d_d
        )
    end
end
