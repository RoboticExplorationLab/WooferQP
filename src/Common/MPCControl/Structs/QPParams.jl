struct OptimizerParams{T,S}
    # discretization length
    dt::T

    # planning horizon length
    N::S

    # Parametron Variables:
    model::Model
    x::Vector{Variable}
    u::Vector{Variable}

    Q_f::Diagonal
    Q_i::Diagonal
    R_i::Diagonal

    # Parameter References
    x_ref::Vector{Base.RefValue{SVector{12,T}}}
    u_ref::Vector{Base.RefValue{SVector{12,T}}}

    A_d::Vector{Base.RefValue{SMatrix{12,12,T,144}}}
    B_d::Vector{Base.RefValue{SMatrix{12,12,T,144}}}
    d_d::Vector{Base.RefValue{SVector{12,T}}}

    J::SMatrix{3,3,T,9}
    sprung_mass::T

    function OptimizerParams(
        dt::T,
        N::S,
        q::Vector{T},
        r::Vector{T},
        μ,
        min_vert_force,
        max_vert_force;
        n::Integer = 12,
        m::Integer = 12,
    ) where {T<:Number,S<:Integer}
        # initialize model and variables
        model = Model(OSQP.Optimizer(verbose = false))
        x = [Variable(model) for _ = 1:((N+1)*n)]
        u = [Variable(model) for _ = 1:((N)*n)]

        # initialize quadratic cost parameters
        Q_f = Diagonal(SVector{n}(q))
        Q_i = Diagonal(SVector{n}(q))
        R_i = Diagonal(SVector{m}(r))

        x_ref_param =
            [Parameter(model, val = zero(SVector{12,T})) for _ = 1:(N+1)]
        u_ref_param =
            [Parameter(model, val = zero(SVector{12,T})) for _ = 1:(N)]

        A_d_param =
            [Parameter(model, val = zero(SMatrix{12,12,T})) for _ = 1:(N)]
        B_d_param =
            [Parameter(model, val = zero(SMatrix{12,12,T})) for _ = 1:(N)]
        d_d_param = [Parameter(model, val = zero(SVector{12,T})) for _ = 1:(N)]


        # terminal state cost
        objective = @expression transpose(x[select12(N + 1)] - x_ref_param[N+1]) *
                    Q_f *
                    (x[select12(N + 1)] - x_ref_param[N+1])

        for i = 1:N
            # stagewise state penalty
            objective = @expression objective +
                        transpose(x[select12(i)] - x_ref_param[i]) *
                        Q_i *
                        (x[select12(i)] - x_ref_param[i])

            #stagewise control penalty
            objective = @expression objective +
                        transpose(u[select12(i)] - u_ref_param[i]) *
                        R_i *
                        (u[select12(i)] - u_ref_param[i])
        end

        @objective(model, Minimize, objective)

        # Dynamics constraint:
        @constraint(model, x[select12(1)] == x_ref_param[1])
        for i = 1:(N)
            @constraint(
                model,
                x[select12(i + 1)] ==
                A_d_param[i] * (x[select12(i)] - x_ref_param[i]) +
                B_d_param[i] * (u[select12(i)] - u_ref_param[i]) +
                d_d_param[i]
            )
        end

        # Control Constraints
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

        ref_arr_x = [x_ref_param[i].val for i = 1:(N+1)]
        ref_arr_u = [u_ref_param[i].val for i = 1:(N)]

        ref_arr_A_d = [A_d_param[i].val for i = 1:(N)]
        ref_arr_B_d = [B_d_param[i].val for i = 1:(N)]
        ref_arr_d_d = [d_d_param[i].val for i = 1:(N)]

        new{T,S}(
            dt,
            N,
            model,
            x,
            u,
            Q_f,
            Q_i,
            R_i,
            ref_arr_x,
            ref_arr_u,
            ref_arr_A_d,
            ref_arr_B_d,
            ref_arr_d_d,
            woofer.inertial.body_inertia,
            woofer.inertial.sprung_mass,
        )
    end
end
