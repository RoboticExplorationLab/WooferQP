struct OptimizerParams
    # discretization length
    dt::Float64

    # planning horizon length
    N::Int64

    # Parametron Variables:
    model::Model
    x::Array{Variable,1}
    u::Array{Variable,1}

    # Parameter Values (to be updated in place)
    x_ref_reshaped::Vector{Float64}
    u_ref::Vector{Float64}
    A_d::Array{Array{Float64,2},1}
    B_d::Array{Array{Float64,2},1}
    d_d::Array{Array{Float64,1},1}
    Q::Array{Float64,2}
    R::Array{Float64,2}

    function OptimizerParams(
        dt::T,
        N::Integer,
        q::Vector{T},
        r::Vector{T},
        μ::T,
        min_vert_force::T,
        max_vert_force::T;
        n::Integer = 12,
        m::Integer = 12,
    ) where {T<:Number}
        # initialize model and variables
        model = Model(OSQP.Optimizer(verbose = false))
        x = [Variable(model) for _ = 1:((N+1)*n)]
        u = [Variable(model) for _ = 1:((N)*n)]

        # initialize quadratic cost parameters
        Q = Diagonal(repeat(q, N + 1))
        R = Diagonal(repeat(r, N))

        x_ref_reshaped = zeros((N + 1) * n)
        u_ref = zeros((N) * n)

        A_d = [zeros(n, n) for i = 1:N]
        B_d = [zeros(n, m) for i = 1:N]
        d_d = [zeros(n) for i = 1:N]

        Q_param = Parameter(() -> Q, model)
        R_param = Parameter(() -> R, model)

        x_ref_param = Parameter(() -> x_ref_reshaped, model)
        u_ref_param = Parameter(() -> u_ref, model)

        A_d_param = [Parameter(() -> A_d[i], model) for i = 1:N]
        B_d_param = [Parameter(() -> B_d[i], model) for i = 1:N]
        d_d_param = [Parameter(() -> d_d[i], model) for i = 1:N]

        # objective
        @objective(
            model,
            Minimize,
            transpose(x - x_ref_param) * Q_param * (x - x_ref_param) +
            transpose(u - u_ref_param) * R_param * (u - u_ref_param)
        )

        # define all constraints
        @constraint(model, x[select12(0)] == x_ref_param[select12(0)])
        for i = 0:N-1
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

            # Dynamics constraints
            @constraint(
                model,
                x[select12(i + 1)] ==
                A_d_param[i+1] * (x[select12(i)] - x_ref_param[select12(i)]) +
                B_d_param[i+1] * (u[select12(i)] - u_ref_param[select12(i)]) +
                d_d_param[i+1]
            )
        end

        new(dt, N, model, x, u, x_ref_reshaped, u_ref, A_d, B_d, d_d, Q, R)
    end
end
