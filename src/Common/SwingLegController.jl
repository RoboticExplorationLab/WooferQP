function generateFootTrajectory(v_b::Vector, t0::T, tf::T, i::Int, param::ControllerParams) where {T<:Number}
	#=
	Generate a body relative trajectory for the ith foot via spline interpolation
	This function is called when the phase is switched to a no contact phase

	Updated foot location taken from param.param.next_foot_loc
	=#

	foot_loc_cur = param.cur_foot_loc[LegIndexToRange(i)]

	# eps = 0.05
	eps = 0.0

	# generate cubic spline in x,y to get body relative foot trajectory
	A = [t0^3 	t0^2 	t0 	1;
		 tf^3 	tf^2 	tf 	1;
		 3*t0^2 2*t0^2 	1 	0;
		 3*tf^2 2*tf^2 	1 	0]

	# TODO: add in omega cross term here? probably doesn't matter...
	b_x = [foot_loc_cur[1], param.next_foot_loc[3*(i-1)+1], -v_b[1], -v_b[1]]
	b_y = [foot_loc_cur[2], param.next_foot_loc[3*(i-1)+2], -v_b[2], -v_b[2]]

	# generate cubic spline in z to enforce height constraint and terminal velocity constraint
	A_z =  [t0^3				t0^2				t0				1;
			tf^3				tf^2				tf				1;
			(0.5*(tf+t0))^3		(0.5*(tf+t0))^2		(0.5*(tf+t0))	1;
			3*tf^2 				2*tf^2 				1 				0]

	b_z = [foot_loc_cur[3], param.next_foot_loc[3*(i-1)+3], param.swing.step_height, -eps]

	param.swing.foot_trajectories[1:4, i] 	.= A\b_x
	param.swing.foot_trajectories[5:8, i] 	.= A\b_y
	param.swing.foot_trajectories[9:12, i] .= A_z\b_z
end

function calcSwingTorques(cur_vel::Vector{T}, α::Vector{T}, t::T, i::Integer, param::ControllerParams) where {T<:Number}
	#=
	PD cartesian controller around swing leg trajectory

	puts swing_torques into param
	=#

	t_p = [t^3, t^2, t, 1]
	t_v = [3*t^2, 2*t^2, 1, 0]

	r_des = [dot(param.swing.foot_trajectories[1:4,i], t_p),
			 dot(param.swing.foot_trajectories[5:8,i], t_p),
			 dot(param.swing.foot_trajectories[9:12,i], t_p)]

	v_des = [dot(param.swing.foot_trajectories[1:4,i], t_v),
			 dot(param.swing.foot_trajectories[5:8,i], t_v),
			 dot(param.swing.foot_trajectories[9:12,i], t_v)]

	J = LegJacobian(α, i)

	return J' * (param.swing.kp_cart*(r_des - param.cur_foot_loc[LegIndexToRange(i)]) + param.swing.kd_cart*(v_des - cur_vel))
end
