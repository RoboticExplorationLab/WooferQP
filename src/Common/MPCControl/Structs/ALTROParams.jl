struct Quadruped{T, S} <: TO.AbstractModel
	x_ref::Vector{SVector}
	u_ref::Vector{SVector}
	foot_locs::Vector{SVector}
	contacts::Vector{SVector}
	times::SVector
	J::SMatrix{3, 3, T}
	sprung_mass::T

	function Quadruped(dt::T, N::S) where {T <: Real, S <: Integer}
		tf = dt*N
		times = SVector{N+1}(collect(range(0, tf, length=N+1)))
		x_ref = [@SVector zeros(T, 12) for i=1:(N+1)]
		u_ref = [@SVector zeros(T, 12) for i=1:(N+1)]
		foot_locs = [@SVector zeros(T, 12) for i=1:(N+1)]
		contacts = [@SVector zeros(T, 4) for i=1:(N+1)]

		new{T, S}(x_ref, u_ref, foot_locs, contacts, times, woofer.inertial.body_inertia, woofer.inertial.sprung_mass)
	end
end

TO.RobotDynamics.state_dim(::Quadruped) = 12
TO.RobotDynamics.control_dim(::Quadruped) = 12

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
		problem = Problem(model, objective, x_des, tf, x0=zeros(n), constraints=constraints)
		solver = ALTROSolver(problem)
		solver.opts.projected_newton = false

		solve!(solver)

		new(dt, n, m, model, objective, constraints, problem, solver, X0, U0, Q, R)
	end
end
