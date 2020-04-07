using Parameters
using LinearAlgebra

include("Dynamics.jl")
include("Config.jl")
include("Utilities.jl")
include("Rotations.jl")

@with_kw mutable struct StateEstimatorParams
    x::SVector{12, Float64}
    P::SMatrix{12, 12, Float64, 144}

    x_::SVector{12, Float64} = @SVector zeros(12)
    P_::SMatrix{12, 12, Float64, 144} = @SMatrix zeros(12, 12)

    y_meas_joint::SVector{12, Float64} = @SVector zeros(12)
    y_meas_mocap::SVector{6, Float64} = @SVector zeros(6)
    
    y_pred_joint::SVector{12, Float64} = @SVector zeros(12)
    y_pred_mocap::SVector{6, Float64} = @SVector zeros(6)

    C_mocap::SMatrix{6, 12, Float64, 72} = SMatrix{6, 12}([I zeros(3, 9); zeros(3, 3) I zeros(3, 6)])

    nu_joint::SVector{12, Float64} = @SVector zeros(12)
    nu_mocap::SVector{6, Float64} = @SVector zeros(6)

    Q::SMatrix{12, 12, Float64, 144}

    R_mocap::SMatrix{6, 6, Float64, 36}
    R_joint::SMatrix{12, 12, Float64, 144}

    g_n::SVector{3, Float64} = SVector{3}(0, 0, 9.81)

    dt::Float64

    net_body_torque::SVector{3, Float64} = @SVector zeros(3) # body frame
    net_body_force::SVector{3, Float64} = @SVector zeros(3) # inertial frame
end

function InitStateEstimator(
    dt,
    x::SVector{12},
    P::SMatrix{12, 12},
    Q::SMatrix{12, 12},
    R_mocap::SMatrix{6, 6},
    R_joint::SMatrix{12, 12}
)
    return StateEstimatorParams(dt = dt, x = x, P = P, Q = Q, R_mocap = R_mocap, R_joint = R_joint)
end

function LinearizedContinuousDynamics(x::SVector{12, T}, J::SMatrix{3, 3}) where T
    x_ϕ = @SVector [x[4], x[5], x[6]]
    x_w = @SVector [x[10], x[11], x[12]]
    dvdv = Diagonal(@SVector ones(3))
    dphidphi = 0.5 * Vmat() * Rmult(UnitQuaternion(x_w)) * [-x_ϕ' / sqrt(abs(1 - x_ϕ' * x_ϕ)); Diagonal(@SVector ones(3))]
    dphidw = 0.5 * Vmat() * Lmult(reconstructq(x_ϕ)) * Vmat()'
    dwdw = -inv(J) * (skew(x_w) * J - skew(J * x_w))

    r1 = [@SMatrix(zeros(3, 6)) dvdv @SMatrix(zeros(3, 3))]
    r2 = [@SMatrix(zeros(3, 3)) dphidphi @SMatrix(zeros(3, 3)) dphidw]
    r3 = @SMatrix(zeros(3, 12))
    r4 = [@SMatrix(zeros(3,9)) dwdw]
    R = [r1; r2; r3; r4]
    return R
end

function ContinuousDynamics(x::SVector{12}, net_body_force::SVector{3}, net_body_torque::SVector{3}, mass, J::SMatrix{3, 3})::SVector{12}
    x_q = reconstructq(@SVector [x[4], x[5], x[6]])
    x_w = @SVector [x[10], x[11], x[12]]
    x_v = @SVector [x[7], x[8], x[9]]
    return [x_v; 0.5 * Vmat() * Lmult(x_q) * Vmat()' * x_w; net_body_force / mass; inv(J) * (net_body_torque - skew(x_w) * J * x_w)]
end

"""
Propogates the dynamics into est_params.x_, est_params.P_
Then does a measurement update using joint kinematics to update est_params.x, est_params.P
"""
function StateEstimatorDynamicsPropogation(
    dt::T,
    joint_pos::SVector{12},
    joint_vel::SVector{12},
    est_params::StateEstimatorParams,
    lqr_forces::SVector{12},
) where T
	x_ϕ = est_params.x[@SVector [4, 5, 6]]
	x_q = reconstructq(x_ϕ)

	# Robot body mass and inertia
	m = (WOOFER_CONFIG::WooferConfig).SPRUNG_MASS
	J = (WOOFER_CONFIG::WooferConfig).INERTIA

    # Convert from continuous system to discrete system
    # TODO handle when estimator errors when the robot falls
	A_discrete = Diagonal(@SVector ones(12)) + LinearizedContinuousDynamics(est_params.x, J) * dt

	# Vector from CoM to foot in body frame
	est_params.net_body_force = @SVector zeros(3)
	est_params.net_body_torque = @SVector zeros(3)
	for i = 1:4
        leg_joint_pos = joint_pos[SLegIndexToRange(i)]
		leg_vec = ForwardKinematics(leg_joint_pos, i)
		lqr_forces_inertial_i = lqr_forces[SLegIndexToRange(i)]

		# Transform lqr force from inertial frame into body frame
		lqr_force_body_i = inv(x_q) * lqr_forces_inertial_i

		# Calculate body-frame torque
		leg_torque_i = cross(leg_vec, lqr_force_body_i)

		est_params.net_body_force = est_params.net_body_force + lqr_forces_inertial_i
		est_params.net_body_torque = est_params.net_body_torque + leg_torque_i
	end

	# Subtract the weight of the robot from the net forces applied by the legs
	est_params.net_body_force = est_params.net_body_force - est_params.g_n * m

	# Propogate dynamics and covariance
	# x_ is the a priori estimate of the true state, x is the a posteriori estimate of the true state
	# P_ is the a priori estimate of the state covariance, P is the a posteriori estimate of the true state covariance
	xdot = ContinuousDynamics(est_params.x, est_params.net_body_force, est_params.net_body_torque, m, J)
	est_params.x_ =  est_params.x + dt * xdot
	est_params.P_ = A_discrete * est_params.P * A_discrete' + est_params.Q

    # Update based on foot velocity measurements
    phi_ = est_params.x_[@SVector [4, 5, 6]]
    q_ = reconstructq(phi_)
    v_ = est_params.x_[@SVector [7, 8, 9]]
    w_ = est_params.x_[@SVector [10, 11, 12]]

    v1 = LegJacobian(joint_pos[SLegIndexToRange(1)], 1) * joint_vel[SLegIndexToRange(1)]
    v2 = LegJacobian(joint_pos[SLegIndexToRange(2)], 2) * joint_vel[SLegIndexToRange(2)]
    v3 = LegJacobian(joint_pos[SLegIndexToRange(3)], 3) * joint_vel[SLegIndexToRange(3)]
    v4 = LegJacobian(joint_pos[SLegIndexToRange(4)], 4) * joint_vel[SLegIndexToRange(4)]
    est_params.y_meas_joint = [v1; v2; v3; v4]

    # calculate the velocity of the CoM in the body-frame
    r_foot1 = ForwardKinematics(joint_pos[SLegIndexToRange(1)], 1)
    r_foot2 = ForwardKinematics(joint_pos[SLegIndexToRange(2)], 2)
    r_foot3 = ForwardKinematics(joint_pos[SLegIndexToRange(3)], 3)
    r_foot4 = ForwardKinematics(joint_pos[SLegIndexToRange(4)], 4)
    x_v_b = inv(q_) * v_
    est_params.y_pred_joint = - [x_v_b; x_v_b; x_v_b; x_v_b] - [cross(w_, r_foot1); cross(w_, r_foot2); cross(w_, r_foot3); cross(w_, r_foot4)]

    # TODO: check that I have the correct rotation
    r = -rotmat(inv(q_))
    rot_mats = [r; r; r; r]
    skews = [skew(r_foot1); skew(r_foot2); skew(r_foot3); skew(r_foot4)]
    C_joint = [@SMatrix(zeros(12, 6)) rot_mats skews]

	# update based on foot measurements
	est_params.nu_joint = est_params.y_meas_joint - est_params.y_pred_joint
	S = C_joint * est_params.P_ * C_joint' + est_params.R_joint
    K = est_params.P_ * C_joint' * inv(S)
    
	# TODO: implement outlier detection that still allows multiple feets at once
	est_params.x = est_params.x_ + K * est_params.nu_joint

	# Joseph Form covariance update
    est_params.P = (I - K * C_joint) * est_params.P_ * (I - K * C_joint)' + K * est_params.R_joint * K'
end

"""
Uses Mocap data to sequentially update the est_params.x and est_params.P
"""
function StateEstimatorMocapUpdate(
	r_hat::SVector{3, T},
    q::UnitQuaternion,
    est_params::StateEstimatorParams,
) where T
    q_hat = vector(q)
    C = est_params.C_mocap
	est_params.y_pred_mocap = C * est_params.x
    est_params.y_meas_mocap = [r_hat; q_hat]
    est_params.nu_mocap = est_params.y_meas_mocap - est_params.y_pred_mocap
	K_mocap = est_params.P * C' * inv(C * est_params.P * C' + est_params.R_mocap)
	est_params.x = est_params.x + K_mocap * est_params.nu_mocap

    # Joseph Form covariance update
	est_params.P = (I - K_mocap * C) * est_params.P * (I - K_mocap * C)' + K_mocap * est_params.R_mocap * K_mocap'
end
