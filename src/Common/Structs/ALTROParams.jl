struct OptimizerParams
	# discretization length
	dt::Float64

	# planning horizon length
	N::Int64

	# Parametron Variables:
	model::Quadruped
	x::Array{Variable, 1}
	u::Array{Variable, 1}

	# Parameter Values (to be updated in place)
	x_ref_reshaped::Vector{Float64}
	u_ref::Vector{Float64}
	Q::Array{Float64, 2}
	R::Array{Float64, 2}

	function OptimizerParams(dt::Float64, N::Int64, q::Vector{Float64}, r::Vector{Float64}; n::Int64=12, m::Int64=12)

		# constants
		μ = 1.0
		min_vert_force = 1
		max_vert_force = 133
		cons = TO.ConstraintSet(n,m,N)
		add_constraint!(cons, BoundConstraint(n,m, u_min=-10, u_max=10), 1:N-1)


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
			@constraint(model, x[select12(i+1)] == A_d_param[i+1]*(x[select12(i)] - x_ref_param[select12(i)])
			 										+ B_d_param[i+1]*(u[select12(i)] - u_ref_param[select12(i)]) +
													d_d_param[i+1])
		end

		new(dt, N, model, x, u, x_ref_reshaped, u_ref, A_d, B_d, d_d, Q, R)
	end
end
