struct OptimizerParams
	# discretization length
	dt::Float64

	# planning horizon length
	N::Int64

	# Parametron Variables:
	model::Model
	x::Array{Variable, 1}
	u::Array{Variable, 1}

	# Parameter Values (to be updated in place)
	x_ref_reshaped::Vector{Float64}
	u_ref::Vector{Float64}
	A_d::Array{Array{Float64, 2}, 1}
	B_d::Array{Array{Float64, 2}, 1}
	d_d::Array{Array{Float64, 1}, 1}
	Q::Array{Float64, 2}
	R::Array{Float64, 2}

	function OptimizerParams(dt::Float64, N::Int64, q::Vector{Float64}, r::Vector{Float64}; n::Int64=12, m::Int64=12)
		# initialize model and variables
		model = Model(OSQP.Optimizer(verbose=false))
		x = [Variable(model) for _ = 1:((N+1)*n)]
		u = [Variable(model) for _ = 1:((N)*n)]

		# initialize quadratic cost parameters
		Q = Diagonal(repeat(q, N+1))
		R = Diagonal(repeat(r, N))

		x_ref_reshaped = zeros((N+1)*n)
		u_ref = zeros((N)*n)

		A_d = [zeros(n,n) for i=1:N]
		B_d = [zeros(n,m) for i=1:N]
		d_d = [zeros(n) for i=1:N]

		Q_param = Parameter(()->Q, model)
		R_param = Parameter(()->R, model)

		x_ref_param = Parameter(()->x_ref_reshaped, model)
		u_ref_param = Parameter(()->u_ref, model)

		A_d_param = [Parameter(()->A_d[i], model) for i=1:N]
		B_d_param = [Parameter(()->B_d[i], model) for i=1:N]
		d_d_param = [Parameter(()->d_d[i], model) for i=1:N]

		# constants
		μ = 1.0
		min_vert_force = 1
		max_vert_force = 133

		# objective
		@objective(model, Minimize, transpose(x-x_ref_param)*Q_param*(x-x_ref_param) + transpose(u-u_ref_param)*R_param*(u-u_ref_param))

		# define all constraints
		@constraint(model, x[select12(0)] == x_ref_param[select12(0)])
		for i=0:N-1
			# Control constraints
			for j=1:4
				# convert absolute value constraint to linear inequality:
				@constraint(model, u[select12_3(i,j,1)] <= μ*u[select12_3(i,j,3)])
				@constraint(model, u[select12_3(i,j,1)] >= -μ*u[select12_3(i,j,3)])
				@constraint(model, u[select12_3(i,j,2)] <= μ*u[select12_3(i,j,3)])
				@constraint(model, u[select12_3(i,j,2)] >= -μ*u[select12_3(i,j,3)])

				@constraint(model, u[select12_3(i,j,3)] >= min_vert_force)
				@constraint(model, u[select12_3(i,j,3)] <= max_vert_force)
			end

			# Dynamics constraints
			@constraint(model, x[select12(i+1)] == A_d_param[i+1]*x[select12(i)] + B_d_param[i+1]*u[select12(i)] + d_d_param[i+1])
		end

		new(dt, N, model, x, u, x_ref_reshaped, u_ref, A_d, B_d, d_d, Q, R)
	end
end




struct SwingLegParams
	foot_trajectories::Array{Float64, 2}

	# z distance below robot of foot at middle of swing leg trajectory
	step_height::Float64

	kp_cart::Float64
	kd_cart::Float64

	function SwingLegParams(step_height, wn_cart, zeta_cart)
		foot_trajectories = zeros(12, 4)

		kp_cart = wn_cart^2
		kd_cart = 2*wn_cart*zeta_cart

		new(foot_trajectories, step_height, kp_cart, kd_cart)
	end
end

struct GaitParams
	# add in a cyclic array for the phases?
	num_phases::Int64

	# 4xnum_phase array of contacts for each gait phase
	contact_phases::Array{Int64}
	phase_times::Vector{Float64}

	phase_length::Float64
	alpha::Float64
	beta::Float64
end

GaitParams(num_phases, contact_phases, phase_times) = GaitParams(	num_phases,
																	contact_phases,
																	phase_times,
																	sum(phase_times),
																	0.5,
																	0.5	)

function createTrotGait()
	num_phases = 4
	contact_phases = [	1 1 1 0;
						1 0 1 1;
						1 0 1 1;
						1 1 1 0	]
	phase_times = [0.6, 0.2, 0.6, 0.2]

	return GaitParams(num_phases, contact_phases, phase_times)
end

function createStandingGait()
	return GaitParams(1, [1;1;1;1], [1.0])
end

mutable struct ControllerParams
	# initialize everything once
	mpc_torques
	swing_torques

	prev_phase
	cur_phase

	cur_foot_loc # current foot location calculated by FK
	active_feet # active feet on ground
	active_feet_12 # expanded active feet on ground

	planner_foot_loc # footstep location for FootstepPlanner
	next_foot_loc # actual planned next foot step for SwingLegController

	nom_foot_loc # foot location with all joint angles = 0

	N # mpc number of time steps
	use_lqr # use LQR terminal cost to go in optimization
	vel_ctrl # velocity based control (reference integrates position)

	mpc_update # rate at which mpc_forces are updated
	last_t # last time that foot forces were calculated


	contacts # contact modes over optimization horizon
	foot_locs # body relative foot locations over optimization horizon
	x_ref # state reference trajectory over optimization horizon
	forces # first step of mpc forces
	x_des # desired state for mpc

	optimizer::OptimizerParams
	gait::GaitParams
	swing::SwingLegParams

	function ControllerParams(	N::Int64,
								mpc_update::Float64,
								x_des::Vector{Float64},
								use_lqr::Bool,
								vel_ctrl::Bool,
								optimizer::OptimizerParams,
								gait::GaitParams,
								swing::SwingLegParams	)

		mpc_torques = zeros(12)
		swing_torques = zeros(12)

		prev_phase = 1
		cur_phase = 1

		cur_foot_loc = zeros(12)
		active_feet = zeros(Int64, 4)
		active_feet_12 = zeros(Int64, 12)
		planner_foot_loc = zeros(12)
		next_foot_loc = zeros(12)

		nom_foot_loc = ForwardKinematicsAll(zeros(12))

		# ensures that foot forces are calculated at start
		last_t = -1

		contacts = zeros(Int64, 4, N)
		foot_locs = zeros(12, N)

		x_ref = zeros(12, N+1)

		forces = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*woofer.inertial.sprung_mass*9.81/4

		new(mpc_torques, swing_torques, prev_phase, cur_phase, cur_foot_loc,
			active_feet, active_feet_12, planner_foot_loc, next_foot_loc, nom_foot_loc, N, use_lqr, vel_ctrl, mpc_update,
			last_t, contacts, foot_locs, x_ref, forces, x_des, optimizer, gait, swing)
	end
end
