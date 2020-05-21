function nextFootstepLocation(v_b::Vector{T}, ω_z::T, next_foot_phase::Int, i::Int, param::ControllerParams) where {T<:Number}
	# implement body velocity heuristic to get next body relative foot location
	ω = [0, 0, ω_z]

	xy_select = Diagonal([1, 1, 0])

	# TODO: use rotation matrix instead of ω cross term
	next_foot_loc = param.cur_foot_loc[LegIndexToRange(i)] + param.gait.alpha*param.gait.phase_times[next_foot_phase]*xy_select*v_b +
						param.gait.beta*param.gait.phase_times[next_foot_phase]*SkewSymmetricMatrix(ω)*param.cur_foot_loc[LegIndexToRange(i)]

	return next_foot_loc
end

function constructFootHistory!(t::T, param::ControllerParams) where {T<:Number}
	# construct the contact and foot location history for MPC solver
	# inputs:
	#
	# outputs (in param):
	# contacts
	# foot_locs

	t_i = t + param.optimizer.dt

	prev_phase = getPhase(t, param)

	next_foot_loc = zeros(3)

	prev_foot_locs = zeros(12)

	r_dot = zeros(3)

	prev_foot_locs .= param.cur_foot_loc


	# current contact is first column of matrix
	param.contacts[:,1] .= param.gait.contact_phases[:, prev_phase]

	# cur_foot_loc is first column of matrix
	param.foot_locs[:,1] .= param.cur_foot_loc

	for i in 2:(param.N)
		next_phase = getPhase(t_i, param)

		param.contacts[:, i] .= param.gait.contact_phases[:, next_phase]

		# rotation matrix from body to inertial
		R_b_n = QuatToRotMatrix(ThreeParamToQuat(param.x_ref[4:6, i]))
		v_b_i = R_b_n'*param.x_ref[7:9, i]
		ω_b_i = param.x_ref[10:12, i]

		for j in 1:4
			if param.gait.contact_phases[j, prev_phase] == 1
				if param.gait.contact_phases[j, next_phase] == 0
					# next foot placement must be planned prior to foot being released
					# next_foot_phase = nextPhase(next_phase, param)
					param.planner_foot_loc[LegIndexToRange(j)] .= nextFootstepLocation(v_b_i, ω_b_i[3], nextPhase(next_phase, param), j, param)
					prev_foot_locs[LegIndexToRange(j)] .= param.planner_foot_loc[LegIndexToRange(j)]
				else
					# integrate param.x_ref via midpoint to get body relative foot location in the future
					r_dot_i = rdot(prev_foot_locs[LegIndexToRange(j)], ω_b_i, v_b_i)
					r_dot_mid = rdot(prev_foot_locs[LegIndexToRange(j)] + param.optimizer.dt/2*r_dot_i, ω_b_i, v_b_i)
					prev_foot_locs[LegIndexToRange(j)] .= prev_foot_locs[LegIndexToRange(j)] + param.optimizer.dt*r_dot_mid
				end
			else
				if param.gait.contact_phases[j, next_phase] == 1
					prev_foot_locs[LegIndexToRange(j)] .= param.planner_foot_loc[LegIndexToRange(j)]
				else
					# doesn't matter if foot not in contact
					prev_foot_locs[LegIndexToRange(j)] .= zeros(3)
				end
			end
		end

		param.foot_locs[:, i] .= prev_foot_locs

		t_i += param.optimizer.dt
		prev_phase = next_phase
	end
end

function rdot(r_b::Vector{T}, om_b::Vector{T}, v_b::Vector{T}) where {T<:Number}
	# d/dt of position vector between COM and stationary foot
	return -SkewSymmetricMatrix(om_b)*r_b - v_b
end
