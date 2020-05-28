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
