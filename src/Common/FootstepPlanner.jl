@with_kw struct FootstepPlannerParams
	# TODO: make sure this next_foot_loc = swing leg next_foot_loc
	next_foot_locs::Vector{Float64} = zeros(12) # TODO: make this a more intelligent initialization for cases where only two feet start active?
end

function nextFootstepLocation!(foot_loc, cur_foot_loc::Vector{T}, v_b::Vector{T}, ω_z::T, gait::GaitParams, next_foot_phase::Int, i::Int) where {T<:Number}
	# implement body velocity heuristic to get next body relative foot location
  	# r_i = Vector((WOOFER_CONFIG::WooferConfig).HIP_LAYOUT[i, :])
	ω = [0, 0, ω_z]

	# foot_loc[1] = cur_foot_loc[1] + gait.alpha*v_b[1]*gait.phase_times[next_foot_phase] #- gait.beta*r_i[2]*ω_z*gait.phase_times[next_foot_phase]
	# foot_loc[2] = cur_foot_loc[2] + gait.alpha*v_b[2]*gait.phase_times[next_foot_phase] #+ gait.beta*r_i[1]*ω_z*gait.phase_times[next_foot_phase]
	# foot_loc[3] = cur_foot_loc[3]
	xy_select = Diagonal([1, 1, 0])

	foot_loc[LegIndexToRange(i)] .= cur_foot_loc[LegIndexToRange(i)] + gait.alpha*gait.phase_times[next_foot_phase]*xy_select*v_b +
						gait.beta*gait.phase_times[next_foot_phase]*SkewSymmetricMatrix(ω)*cur_foot_loc[LegIndexToRange(i)]
end

function constructFootHistory!(contacts_future::Array{Int, 2}, foot_locs_future::Array{T}, t::T, x_ref::Array{T}, cur_foot_locs::Vector{T}, mpc_config::MPCControllerParams, gait::GaitParams, foot_params::FootstepPlannerParams) where {T<:Number}
	# construct the contact and foot location history for MPC solver
	# cur_foot_locs: 12 dimension vector
	# TODO: make sure this really works

	t_i = t + mpc_config.dt

	prev_phase = getPhase(t, gait)

	next_foot_loc = zeros(3)

	prev_foot_locs = zeros(12)

	r_dot = zeros(3)

	prev_foot_locs .= cur_foot_locs


	# current contact is first column of matrix
	contacts_future[:,1] .= gait.contact_phases[:, prev_phase]

	# cur_foot_loc is first column of matrix
	foot_locs_future[:,1] .= cur_foot_locs

	for i in 2:(mpc_config.N)
		next_phase = getPhase(t_i, gait)

		contacts_future[:, i] .= gait.contact_phases[:, next_phase]

		for j in 1:4
			if gait.contact_phases[j, prev_phase] == 1
				if gait.contact_phases[j, next_phase] == 0
					# next foot placement must be planned prior to foot being released
					# next_foot_phase = nextPhase(next_phase, gait)

					# TODO: add these terms
					# v_b = R*x_est[7:9]
					# ω = R*x_est[10:12]

					nextFootstepLocation!(foot_params.next_foot_locs, prev_foot_locs, x_ref[7:9, i], x_ref[12, i], gait, nextPhase(next_phase, gait), j)
					prev_foot_locs[LegIndexToRange(j)] .= foot_params.next_foot_locs[LegIndexToRange(j)]
				else
					# integrate x_ref via midpoint to get body relative foot location in the future
					rdot!(r_dot, prev_foot_locs[LegIndexToRange(j)], x_ref[10:12, i], x_ref[7:9, i])
					rdot!(r_dot, prev_foot_locs[LegIndexToRange(j)] + mpc_config.dt/2*r_dot, x_ref[10:12, i], x_ref[7:9, i])
					prev_foot_locs[LegIndexToRange(j)] .= prev_foot_locs[LegIndexToRange(j)] + mpc_config.dt*r_dot
				end
			else
				if gait.contact_phases[j, next_phase] == 1
					prev_foot_locs[LegIndexToRange(j)] .= foot_params.next_foot_locs[LegIndexToRange(j)]
				else
					# doesn't matter if foot not in contact
					prev_foot_locs[LegIndexToRange(j)] .= zeros(3)
				end
			end
		end

		foot_locs_future[:, i] .= prev_foot_locs

		t_i += mpc_config.dt
		prev_phase = next_phase
	end
end

function rdot!(r_dot::Vector{T}, r_b::Vector{T}, om_b::Vector{T}, v_b::Vector{T}) where {T<:Number}
	# d/dt of position vector between COM and stationary foot

	A = SkewSymmetricMatrix(om_b)

	r_dot .= -A*r_b - v_b
end
