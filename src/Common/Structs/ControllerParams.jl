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
