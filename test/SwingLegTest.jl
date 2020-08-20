using LinearAlgebra
using Parameters
using Plots

include("../src/Common/SwingLegController.jl")
include("../src/Common/Dynamics.jl")

swing_params = SwingLegParams()
swing_params.next_foot_loc[LegIndexToRange(i)] .= ForwardKinematics([0.1, -0.1, 0.1], i)

i = 1

v_b = zeros(3)
cur_foot_loc = ForwardKinematics(zeros(3), i)

t0 = 0.0
tf = 0.5


generateFootTrajectory(cur_foot_loc, v_b, t0, tf, i, swing_params)

t = collect(t0:0.01:tf)

r_array = zeros(3, size(t,1))
for j=1:size(t,1)
	t_p = [t[j]^3, t[j]^2, t[j], 1]
	r_array[:, j] = [dot(swing_params.foot_trajectories[1:4,i], t_p),
					 dot(swing_params.foot_trajectories[5:8,i], t_p),
					 dot(swing_params.foot_trajectories[9:12,i], t_p)]
end

plot(t, r_array[1, :])
