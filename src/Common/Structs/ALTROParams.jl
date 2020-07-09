mutable struct Quadruped <: AbstractModel
	config::WooferConfig
	x_ref::Vector{Vector{Float64}}
	u_ref::Vector{Vector{Float64}}
	foot_locs::Vector{Vector{Float64}}
	contacts::Vector{Vector{Float64}}
	times::Vector{Float64}
	n::Int64
	m::Int64

	function Quadruped(config::WooferConfig, dt::Float64, N::Int64, n::Int64, m::Int64)
		tf = dt*N
		times = collect(range(0, tf, length=N+1))
		x_ref = [zeros(12) for i=1:(N+1)]
		u_ref = [zeros(12) for i=1:(N+1)]
		foot_locs = [zeros(12) for i=1:(N+1)]
		contacts = [zeros(4) for i=1:(N+1)]

		new(config, x_ref, u_ref, foot_locs, contacts, times, n, m)
	end
end

mutable struct OptimizerParams
	# discretization length
	dt::Float64

	# state and control size
	n::Int64
	m::Int64

	# ALTRO Variables:
	model::Quadruped
	objective::Objective
	constraints::ConstraintSet
	problem::Problem
	solver::ALTROSolver

	X0::Vector{Vector{Float64}}
	U0::Vector{Vector{Float64}}

	Q
	R

	function OptimizerParams(dt::Float64, N::Int64, q::Vector{Float64}, r::Vector{Float64}, x_des::Vector{Float64}; n::Int64=12, m::Int64=12)
		# constants
		μ = 1.0
		min_vert_force = 0.0
		max_vert_force = 133.0

		Q = Diagonal(SVector{n}(q))
		R = Diagonal(SVector{m}(r))


		# TODO: better initialization here
		X0 = [zeros(n) for i=1:(N+1)]
		U0 = [zeros(m) for i=1:N]

		model = Quadruped(woofer, dt, N, n, m)

		constraints = TrajectoryOptimization.ConstraintSet(n,m,N)

		friction = FrictionConstraint(m, μ)
		add_constraint!(constraints, friction, 1:N)

		u_min = [-Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force]
		u_max = [Inf, Inf, max_vert_force, Inf, Inf, max_vert_force, Inf, Inf, max_vert_force, Inf, Inf, max_vert_force]
		bound = BoundConstraint(n,m, u_min=u_min, u_max=u_max)
		add_constraint!(constraints, bound, 1:N)

		# objective
		objective = LQRObjective(Q, R, Q, x_des, N)

		tf = dt*N
		problem = Problem(model, objective, x_des, tf, x0=zeros(n), constraints=constraints)
		solver = ALTROSolver(problem)

		TrajectoryOptimization.solve!(solver)

		new(dt, n, m, model, objective, constraints, problem, solver, X0, U0, Q, R)
	end
end
