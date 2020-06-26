using StaticArrays
using LinearAlgebra
using TrajectoryOptimization
using Parametron
using OSQP

include("../src/Common/Config.jl")
include("../src/Common/Quaternions.jl")
include("../src/Common/Structs/ALTROParams.jl")
include("../src/Common/Structs/SwingLegParams.jl")
include("../src/Common/Structs/GaitParams.jl")
include("../src/Common/Structs/ControllerParams.jl")
include("../src/Common/ALTROSolver.jl")
include("../src/Common/SwingLegController.jl")
include("../src/Common/Gait.jl")
include("../src/Common/FootstepPlanner.jl")

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
optimizer = OptimizerParams(woofer, dt, N, q, r, x_des)
param = ControllerParams(N, mpc_update, x_des, use_lqr, vel_ctrl, zeros(12), optimizer, gait, swing)

x_ref = zeros(12, N+1)

generateReferenceTrajectory!(x0, param)

contacts = zeros(Int64, 4, N)
foot_locs = zeros(12, N)

α = zeros(12)

cur_foot_loc = ForwardKinematicsAll(α)

t = 0.0

constructFootHistory!(t, param)

solveFootForces!(param)

param.forces
