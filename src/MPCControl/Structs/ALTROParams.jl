mutable struct OptimizerParams{T, S, P, A}
	# discretization length
	dt::T

	# state and control size
	n::S
	m::S

	# ALTRO Variables:
	model::Quadruped{T, S}
	objective::Objective
	constraints::ConstraintList
	problem::P
	solver::A

	X0::Vector{SVector{12, T}}
	U0::Vector{SVector{12, T}}

	Q
	R

	function OptimizerParams(
							    dt::T,
							    N::S,
							    q::AbstractVector{T},
							    r::AbstractVector{T},
							    x_des::AbstractVector{T},
							    μ::T,
							    min_vert_force::T,
							    max_vert_force::T;
							    n::Integer = 12,
							    m::Integer = 12,
							) where {T<:Number, S<:Integer}
		Q = Diagonal(SVector{n}(q))
		R = Diagonal(SVector{m}(r))


		# TODO: better initialization here
		u_guess = woofer.inertial.sprung_mass/4*(@SVector [0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1.])

		X0 = [x_des for i=1:(N+1)]
		U0 = [u_guess for i=1:N]

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

		tf = dt*(N-1)
		problem = Problem(model, objective, x_des, tf, x0=zeros(n), constraints=constraints, integration=RD.DiscreteSystemQuadrature)
		solver = ALTROSolver(problem)
		set_options!(solver, projected_newton=false)

		solve!(solver)

		P = typeof(problem)
		A = typeof(solver)

		new{T,S,P,A}(dt, n, m, model, objective, constraints, problem, solver, X0, U0, Q, R)
	end
end