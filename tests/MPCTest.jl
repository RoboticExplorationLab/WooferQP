using StaticArrays
using Parameters
using LinearAlgebra
using Rotations

include("../src/Common/Config.jl")
include("../src/Common/QPSolverSparse_convex.jl")
include("../src/Common/Dynamics.jl")
include("../src/Common/Gait.jl")
include("../src/Common/FootstepPlanner.jl")
include("../src/Common/Rotations.jl")

dt = 0.05
N = 1

mpc_config = MPCControllerParams(dt, N)
# standingGait = GaitParams(num_phases=1, contact_phases=[1;1;1;1], phase_times=[1.0])
walkingGait = GaitParams()
footstep_config = FootstepPlannerParams()

x_des = [0, 0, 0.32, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 0]
# x0 = [0.0, 0, 0.32, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 0]
x0 = [0.0, 0, 0.32, 0.0, 0.0, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 0]

x_ref = zeros(12, N+1)

generateReferenceTrajectory!(x_ref, x0, x_des, mpc_config)

contacts = zeros(Int64, 4, N)
foot_locs = zeros(12, N)

α = @SVector zeros(3)

r1 = ForwardKinematics(α, 1)
r2 = ForwardKinematics(α, 2)
r3 = ForwardKinematics(α, 3)
r4 = ForwardKinematics(α, 4)


cur_foot_loc = zeros(12)
cur_foot_loc[1:3] = r1
cur_foot_loc[4:6] = r2
cur_foot_loc[7:9] = r3
cur_foot_loc[10:12] = r4

t = 0.0

constructFootHistory!(contacts, foot_locs, t, x_ref, cur_foot_loc, mpc_config, standingGait, footstep_config)

foot_locs
contacts

forces = zeros(12)

# @time solveFootForces!(forces, x0, x_ref, contacts, foot_locs, mpc_config, WOOFER_CONFIG)
(u, x) = solveFootForces!(forces, x_ref, contacts, foot_locs, mpc_config, false)

x[121:132]

u[1:12]
