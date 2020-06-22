using OSQP
using Parametron
using LinearAlgebra

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
n = 12
m = 12

x_des = [0, 0, 0.32, 0, 0, 0, 0, 0, 0, 0, 0, 0]

model = Model(OSQP.Optimizer())

x = [Variable(model) for _ = 1:((N+1)*n)]
u = [Variable(model) for _ = 1:((N)*n)]

Q = Diagonal(repeat([1e3, 1e3, 5e4, 1e3, 1e3, 1e3, 1e2, 1e2, 1e2, 1, 1, 1e2], N+1))
u_ref_i = zeros(12)

u_ref = repeat(u_ref_i, N)

R = Diagonal(repeat([1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4], N))

x_ref = repeat(x_des, (N+1))

x_ref_param = Parameter(model, val=x_ref)
u_ref_param = Parameter(model, val=u_ref)

@objective(model, Minimize, transpose(x-x_ref_param)*Q*(x-x_ref_param) + transpose(u-u_ref_param)*R*(u-u_ref_param))

# constants
μ = 0.7
min_vert_force = 1
max_vert_force = 133

@constraint(model, x[select12(0)] == x_ref_param[select12(0)])

for i=0:N-1
	# Control constraints
	for j=1:4
		# # use absolute value friction constraint:
		# problem.constraints += abs(u[select12_3(i,j,1)]) <= μ*u[select12_3(i,j,3)]
		# problem.constraints += abs(u[select12_3(i,j,2)]) <= μ*u[select12_3(i,j,3)]

		# convert absolute value constraint to linear inequality:
		@constraint(model, u[select12_3(i,j,1)] <= μ*u[select12_3(i,j,3)])
		@constraint(model, u[select12_3(i,j,1)] >= -μ*u[select12_3(i,j,3)])
		@constraint(model, u[select12_3(i,j,2)] <= μ*u[select12_3(i,j,3)])
		@constraint(model, u[select12_3(i,j,2)] >= -μ*u[select12_3(i,j,3)])

		@constraint(model, u[select12_3(i,j,3)] >= min_vert_force)
		@constraint(model, u[select12_3(i,j,3)] <= max_vert_force)
	end

	A_c_i = LinearizedContinuousDynamicsA(x_ref[select12(i)], woofer.inertial.body_inertia)
	B_c_i = LinearizedContinuousDynamicsB(x_ref[select12(i)], foot_locs[:, i+1], contacts[:, i+1], woofer.inertial.body_inertia)
	d_c_i = NonLinearContinuousDynamics(x_ref[select12(i)], u_ref[select12(i)], foot_locs[:, i+1], woofer.inertial.body_inertia)

	# Euler Discretization
	A_d_i = Array(I + A_c_i*dt)
	B_d_i = Array(B_c_i*dt)
	d_d_i = Array(d_c_i*dt)

	# Dynamics constraints
	@constraint(model, x[select12(i+1)] == A_d_i*x[select12(i)] + B_d_i*u[select12(i)] + d_d_i)
end

solve!(model)

value.(model, u)[1:12]

value.(model, x)[1:12]
