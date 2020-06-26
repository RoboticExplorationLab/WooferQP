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
		u_ref = [zeros(12) for i=1:(N)]
		foot_locs = [zeros(12) for i=1:(N)]
		contacts = [zeros(4) for i=1:(N)]

		new(config, x_ref, u_ref, foot_locs, contacts, times, n, m)
	end
end

struct OptimizerParams
	# discretization length
	dt::Float64

	# planning horizon length
	N::Int64
	# state and control size
	n::Int64
	m::Int64

	# ALTRO Variables:
	model::Quadruped
	objective::Objective
	constraints::ConstraintSet
	X0::Vector{Vector{Float64}}
	U0::Vector{Vector{Float64}}

	Q
	R

	function OptimizerParams(config::WooferConfig, dt::Float64, N::Int64, q::Vector{Float64}, r::Vector{Float64}, x_des::Vector{Float64}; n::Int64=12, m::Int64=12)
		# constants
		μ = 1.0
		min_vert_force = 1
		max_vert_force = 133
		constraints = TrajectoryOptimization.ConstraintSet(n,m,N)

		Q = Diagonal(SVector{n}(q))
		R = Diagonal(SVector{m}(r))

		X0 = [zeros(n) for i=1:(N+1)]
		U0 = [zeros(m) for i=1:N]

		model = Quadruped(config, dt, N, n, m)

		# objective
		objective = LQRObjective(Q, R, Q, x_des, N)

		tf = dt*N

		prob = Problem(model, objective, x_des, tf, constraints=constraints)

		new(dt, N, n, m, model, objective, constraints, X0, U0, Q, R)
	end
end
