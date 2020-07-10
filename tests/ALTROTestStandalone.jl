using StaticArrays
using LinearAlgebra
using TrajectoryOptimization
using TimerOutputs

include("../src/Common/Config.jl")
include("../src/Common/Quaternions.jl")
include("../src/Common/Structs/LinearizedFrictionConstraint.jl")
include("../src/Common/Structs/ALTROParams.jl")
include("../src/Common/Structs/SwingLegParams.jl")
include("../src/Common/Structs/GaitParams.jl")
include("../src/Common/Structs/ControllerParams.jl")
include("../src/Common/ALTROSolver.jl")
include("../src/Common/SwingLegController.jl")
include("../src/Common/Gait.jl")
include("../src/Common/FootstepPlanner.jl")
include("../src/Common/Dynamics.jl")

dt = 0.05
mpc_update = 0.001
N = 15

q = [0, 0, 5e4, 1e3, 1e5, 1e3, 1e4, 1e4, 1e2, 1e2, 1e2, 1e4]
r = [1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4]

x_des = [0, 0, 0.32, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 0]
x0 = [0.0, 0, 0.32, 0.0, 0.0, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 0]

use_lqr = false # use lqr in cost to go
vel_ctrl = false # integrate positions, interpolate velocities

swing = SwingLegParams(-0.20, 100, 1)
gait = createTrotGait(stance_time=0.15, swing_time=0.15)

const to = TimerOutput()
optimizer = OptimizerParams(dt, N, q, r, x_des)
param = ControllerParams(N, mpc_update, x_des, use_lqr, vel_ctrl, zeros(12), optimizer, gait, swing)

function update_trajectories(x0, t, param)
	generateReferenceTrajectory!(x0, param)
	constructFootHistory!(t, param)
end
################################################################################

model = Quadruped(woofer, dt, N, n, m)

constraints = TrajectoryOptimization.ConstraintSet(n,m,N)
friction = FrictionConstraint(m, μ)
add_constraint!(constraints, friction, 1:N)

u_min = @SVector [-Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force]
u_max = @SVector [Inf, Inf, max_vert_force, Inf, Inf, max_vert_force, Inf, Inf, max_vert_force, Inf, Inf, max_vert_force]
bound = BoundConstraint(n,m, u_min=u_min, u_max=u_max)
add_constraint!(constraints, bound, 1:N)

objective = LQRObjective(Q, R, Q, x_des, N)

tf = dt*N
problem = Problem(model, objective, x_des, tf, x0=x0, constraints=constraints)
solver = ALTROSolver(problem)

@time TrajectoryOptimization.solve!(solver)

X = states(solver)
U = controls(solver)

x_curr = [0.1, 0.1, 0.25, 0.1, 0.1, 0.1, -0.01, 0.01, 0.1, -0.01, 0.1, 0.2]
trajectories(x_curr, 0.5, param)

for i=1:(N+1)
	model.x_ref[i] .= param.x_ref[:, i]
	# opt.model.u_ref .= [zeros(m) for i=1:(param.N+1)]
	model.contacts[i] .= param.contacts[:, i]
	model.foot_locs[i] .= param.foot_locs[:, i]
end
solver.solver_al.solver_uncon.x0 .= x_curr
solver.opts.projected_newton = false
# # set_initial_state!(solver, x_curr)

TrajectoryOptimization.solve!(solver)
X = states(solver)
U = controls(solver)
