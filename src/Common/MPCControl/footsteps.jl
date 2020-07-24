function footstep_location(v_b::AbstractVector{T}, ω_z::T, cur_phase::Integer, i::Integer, param::ControllerParams) where {T<:Number}
	# implement body velocity heuristic to get next body relative foot location
	ω = @SVector [0, 0, ω_z]

	v_b_proj = @SVector [v_b[1], v_b[2], T(0.0)]

	next_phase = get_next_phase(cur_phase, param)

	# k = sqrt(param.x_des[3]/9.81)
	k = T(0.0)

	v_des_proj = @SVector [param.x_des[7], param.x_des[8], T(0.0)]

	# TODO: use rotation matrix instead of ω cross term
	next_foot_loc = param.nom_foot_loc[i] +
						param.gait.alpha*param.gait.phase_times[next_phase]*v_b_proj + k*(v_b_proj - v_des_proj)
						param.gait.beta*param.gait.phase_times[next_phase]*Rotations.skew(ω)*param.cur_foot_loc[i]

	return next_foot_loc
end

function foot_history!(t::Number, param::ControllerParams)
	# construct the contact and foot location history for MPC solver
	# inputs:
	#
	# outputs (in param):
	# contacts
	# foot_locs

	t_i = t + param.optimizer.dt

	prev_phase = get_phase(t, param)

	prev_foot_locs = zero(FootstepLocation)
	prev_foot_locs .= param.cur_foot_loc

	# current contact is first column of matrix
	param.contacts[1] = param.gait.contact_phases[prev_phase]

	# cur_foot_loc is first column of matrix
	param.foot_locs[1] .= param.cur_foot_loc

	for i in 2:(param.N+1)
		next_phase = get_phase(t_i, param)

		param.contacts[i] = param.gait.contact_phases[next_phase]

		x_ref_i = param.x_ref[i]
		v_i = @SVector [x_ref_i[7], x_ref_i[8], x_ref_i[9]]
		ω_b_i = @SVector [x_ref_i[10], x_ref_i[11], x_ref_i[12]]

		# rotation matrix from body to inertial
		r = MRP(x_ref_i[4], x_ref_i[5], x_ref_i[6])
		v_b_i = r \ v_i

		for j in 1:4
			if param.gait.contact_phases[prev_phase][j] == 1
				if param.gait.contact_phases[next_phase][j] == 0
					# next foot placement must be planned prior to foot being released
					param.planner_foot_loc[j] = footstep_location(v_b_i, ω_b_i[3], next_phase, j, param)
					prev_foot_locs[j] = param.planner_foot_loc[j]
				else
					# integrate param.x_ref via midpoint to get body relative foot location in the future
					r_dot_i = rdot(prev_foot_locs[j], ω_b_i, v_b_i)
					r_dot_mid = rdot(prev_foot_locs[j] + param.optimizer.dt/2*r_dot_i, ω_b_i, v_b_i)
					prev_foot_locs[j] = prev_foot_locs[j] + param.optimizer.dt*r_dot_mid
				end
			else
				if param.gait.contact_phases[next_phase][j] == 1
					prev_foot_locs[j] = param.planner_foot_loc[j]
				else
					# doesn't matter if foot not in contact
					prev_foot_locs[j] = @SVector zeros(3)
				end
			end
		end

		param.foot_locs[i] .= prev_foot_locs

		t_i += param.optimizer.dt
		prev_phase = next_phase
	end
end

function rdot(r_b::AbstractVector{T}, om_b::AbstractVector{T}, v_b::AbstractVector{T}) where {T<:Number}
	# d/dt of position vector between COM and stationary foot
	return -Rotations.skew(om_b)*r_b - v_b
end
