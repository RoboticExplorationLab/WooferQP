mutable struct OptimizerParams
	# discretization length
	dt::Float64

	# state and control size
	n::Int64
	m::Int64

	# ALTRO Variables:
	model::Quadruped
	objective::Objective
	constraints::ConstraintList
	problem::Problem
	solver::AugmentedLagrangianSolver

	X0::Vector{Vector{Float64}}
	U0::Vector{Vector{Float64}}

	Q
	R

	function OptimizerParams(
							    dt::T,
							    N::Integer,
							    q::Vector{T},
							    r::Vector{T},
							    x_des::Vector{T},
							    μ::T,
							    min_vert_force::T,
							    max_vert_force::T;
							    n::Integer = 12,
							    m::Integer = 12,
							) where {T<:Number}
		Q = Diagonal(SVector{n}(q))
		R = Diagonal(SVector{m}(r))


		# TODO: better initialization here
		X0 = [zeros(n) for i=1:(N+1)]
		U0 = [zeros(m) for i=1:N]

		model = Quadruped(dt, N)

		constraints = ConstraintList(n,m,N)

		friction = FrictionConstraint(m, μ)
		add_constraint!(constraints, friction, 1:N)

		u_min = @SVector [-Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force]
		u_max = @SVector [Inf, Inf, max_vert_force, Inf, Inf, max_vert_force, Inf, Inf, max_vert_force, Inf, Inf, max_vert_force]
		bound = BoundConstraint(n,m, u_min=u_min, u_max=u_max)
		add_constraint!(constraints, bound, 1:N)

		# objective
		objective = LQRObjective(Q, R, Q, x_des, N)

		tf = dt*N
		problem = Problem(model, objective, x_des, tf, x0=zeros(n), constraints=constraints, integration=DiscreteQuadruped)
		solver = AugmentedLagrangianSolver(problem)

		solve!(solver)

		new(dt, n, m, model, objective, constraints, problem, solver, X0, U0, Q, R)
	end
end
