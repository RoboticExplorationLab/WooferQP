mutable struct ControllerParams
	# initialize everything once
	mpc_torques
	swing_torques

	prev_phase
	cur_phase
	cur_phase_time # current time relative to beginning of phase
	last_replan_t # last replan time
	replan_update # time between replanning foot placements

	cur_foot_loc # current foot location calculated by FK
	active_feet # active feet on ground
	active_feet_12 # expanded active feet on ground

	trajectory_foot_loc # foot location at beginning of trajectory
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

	function ControllerParams(path::String)
		data = YAML.load(open(path))
		N = data["N"]

		mpc_torques = zeros(12)
		swing_torques = zeros(12)

		prev_phase = 1
		cur_phase = 1
		cur_phase_time = 0.0

		last_replan_t = 0.0
		replan_update = 0.01

		cur_foot_loc = zeros(12)
		active_feet = zeros(Int64, 4)
		active_feet_12 = zeros(Int64, 12)

		trajectory_foot_loc = zeros(12)
		planner_foot_loc = zeros(12)
		next_foot_loc = zeros(12)

		# ensures that foot forces are calculated at start
		last_t = -1

		contacts = zeros(Int64, 4, N+1)
		foot_locs = zeros(12, N+1)

		x_ref = zeros(12, N+1)

		forces = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*woofer.inertial.sprung_mass*9.81/4

		x_des = [0.0; 0.0; data["stance_height"]; zeros(3); data["xy_vel"]; zeros(3); data["omega_z"]]

		nom_foot_loc = ForwardKinematicsAll(zeros(12))
		offset = [1 -1 1 -1]
		Δx = data["foot_dx"]
		Δy = data["foot_dx"]
		for i=1:4
			nom_foot_loc[LegIndexToRange(i)] += [Δx, Δy*offset[i], 0]
		end

		optimizer = OptimizerParams(	data["dynamics_discretization"],
										N,
										data["q"],
										data["r"]
									)

		gait_type = data["gait"]["type"]

		if gait_type == "trot"
			gait = createTrotGait(stance_time=data["gait"]["stance_time"], swing_time=data["gait"]["swing_time"])
		elseif gait_type == "pronk"
			gait = createPronkGait(stance_time=data["gait"]["stance_time"], flight_time=data["gait"]["swing_time"])
		elseif gait_type == "pace"
			gait = createPaceGait(stance_time=data["gait"]["stance_time"], swing_time=data["gait"]["swing_time"])
		else
			gait = createStandingGait()
		end

		swing = SwingLegParams(data["swing"]["step_height"], data["swing"]["omega"], data["swing"]["zeta"])

		new(mpc_torques, swing_torques, prev_phase, cur_phase, cur_phase_time, last_replan_t, replan_update, cur_foot_loc,
			active_feet, active_feet_12, trajectory_foot_loc, planner_foot_loc, next_foot_loc,
			nom_foot_loc, N, data["use_lqr"], data["velocity_control"], data["update_dt"], last_t, contacts, foot_locs, x_ref,
			forces, x_des, optimizer, gait, swing)
	end
end
