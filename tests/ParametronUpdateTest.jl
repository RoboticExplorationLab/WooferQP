include("../src/Common/Config.jl")
include("../src/Common/Dynamics.jl")
include("../src/Common/Quaternions.jl")
include("../src/Common/QPSolverSparse_parametron.jl")
include("../src/Common/Gait.jl")
include("../src/Common/FootstepPlanner.jl")

###############################################################################
N = 15
dt = 0.05

mpc_config = MPCControllerParams(dt, N)
standingGait = GaitParams(num_phases=1, contact_phases=[1;1;1;1], phase_times=[1.0])
# walkingGait = GaitParams()
footstep_config = FootstepPlannerParams()

x_des = [0, 0, 0.32, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 0]
# x0 = [0.0, 0, 0.32, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 0]
x0 = [0.0, 0, 0.32, 0.0, 0.0, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 0]

x_ref = zeros(12, N+1)

generateReferenceTrajectory!(x_ref, x0, x_des, mpc_config)

contacts = zeros(Int64, 4, N)
foot_locs = zeros(12, N)

α = zeros(12)

cur_foot_loc = ForwardKinematicsAll(α)

t = 0.0

constructFootHistory!(contacts, foot_locs, t, x_ref, cur_foot_loc, mpc_config, standingGait, footstep_config)

###############################################################################
forces = zeros(12)
solveFootForces!(forces, x_ref, contacts, foot_locs, mpc_config, false)
