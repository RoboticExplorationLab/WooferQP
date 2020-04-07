using LinearAlgebra
using StaticArrays
using Parameters
using LCMCore

include("../Common/Dynamics.jl")
include("../Common/Config.jl")
include("../Common/StateEstimator.jl")
include("../Common/LQRController.jl")
include("../Common/Rotations.jl")
include("../Common/LCMTypes.jl")

"""
Standalone control and estimation loop that takes in simulation data messages and
mocap messages. Outputs torque command messages
"""
function RunLoop()
    # initialize everything once
    lcm = LCM()

    ψ = 0.0 # yaw angle
    ϕ = 0.0 # pitch angle
    x_des = SVector{12}([0.0, 0.0, 0.30, 0.00, ϕ, ψ, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0])

    lower_dt = 0.001
    mocap_dt = 0.008
    x0 = SVector{12}([0.00, 0.0, 0.30, zeros(9)...])
    P = 1e-3 * Diagonal(@SVector ones(12))
    Q = 1e-1 * Diagonal(@SVector ones(12))
    R_mocap = 1e-2 * Diagonal(@SVector ones(6))
    R_joint = 1e-3 * Diagonal(@SVector ones(12))
    est_params = StateEstimatorParams(
        dt = lower_dt,
        x = x0,
        P = P,
        Q = Q,
        R_mocap = R_mocap,
        R_joint = R_joint,
    )
    last_t = 0.0
    last_mocap_t = 0.0

    lqr_params = initLQRParams(
        dt = lower_dt,
        x0 = x_des,
        rxy_penalty = 1e2,
        rz_penalty = 1e2,
        phi_penalty = 1e2,
        v_penalty = 1,
        w_penalty = 1,
        fxy_penalty = 1e-3,
        fz_penalty = 1e-3,
    )
    # TODO: figure out why lqr_forces is still a fucking Core.Box. Plus there's no documentation on Core.Box. Julia is shit.
    lqr_forces::SVector{12, Float64} = lqr_params.u0

    not_moved = true

    t0 = 0.0 # Time for placeholder LCM messages
    # TODO: Find out why all the LCM message variables are not type stable
    torque_message = TorqueCommand(t0, zeros(SVector{12}))
    torque_message.timestamp = 0.0

    # Adding type assert to SimData(blah) helps with type stability
    sim_message = SimData(zeros(SVector{3}), @SVector([1, 0, 0, 0]), zeros(SVector{3}), zeros(SVector{3}), t0, true, zeros(SVector{12}), zeros(SVector{12}), zeros(SVector{3}), zeros(SVector{3}))::SimData

    odrive_message = ODriveData(t0, zeros(SVector{12}), zeros(SVector{12}))::ODriveData
    last_odrive_message_timestamp = 0.0

    # Set up LCM callbacks and subscribers
    function ODriveCallback(channel::String, received_data::ODriveData)
        odrive_message = received_data
        println("odrive dt: ", odrive_message.timestamp - last_odrive_message_timestamp)
        flush(stdout)
        StateEstimatorDynamicsPropogation(
            odrive_message.timestamp - last_odrive_message_timestamp,
            odrive_message.joint_position,
            odrive_message.joint_velocity,
            est_params,
            lqr_forces,
        )
        last_odrive_message_timestamp = odrive_message.timestamp
    end
    subscribe(lcm, "ODRIVE_DATA", ODriveCallback, ODriveData)

    function SimCallback(channel::String, received_sim_data::SimData)
        sim_message = received_sim_data
    end
    subscribe(lcm, "SIM_DATA", SimCallback, SimData)

    # Note: the actual LCM mocap data is not being stored anywhere, only used for the ekf update
    last_mocap_data_timestamp = 0.0
    function MocapCallback(channel::String, received_mocap_data::MocapData)
        mocap_q = UnitQuaternion(received_mocap_data.orientation...)
        StateEstimatorMocapUpdate(
            received_mocap_data.position,
            mocap_q,
            est_params,
        )
        last_mocap_data_timestamp = received_mocap_data.timestamp
    end
    subscribe(lcm, "MOCAP_DATA", MocapCallback, MocapData)

    println("Starting asynchronous lcm handler")
    @async while true
        handle(lcm)
    end

    println("Press enter to control/estimation loop...")
    readline(stdin)
    println("Starting.")

    loop_count = 0
    while true
        if !sim_message.run_loop
            println("Ending control/estimation loop.")
            break
        end

        # if sim_message.timestamp > 0.5 && not_moved
        #     lqr_params.x0 = SVector{12}([0.0, 0.0, 0.30, 0.2, 0, 0, zeros(6)...])
		#
        #     g = [0, 0, 0, 0, 0, 0, 0, 0, -9.81, 0, 0, 0] * lower_dt
        #     lqr_params.u0 = SVector{12}(pinv(lqr_params.B_d) * ((I - lqr_params.A_d) * lqr_params.x0 - g))
        #     not_moved = false
        # end

        # LQR Balance Controller
        # .= necessary to keep lqr_forces the same type

        lqr_forces = LQRBalance(est_params.x, lqr_params)
        q = reconstructq(est_params.x[@SVector [4, 5, 6]])

        # TODO: figure out why the type assert here is necessary
        lqr_torques = InertialLegForcesToActuatorTorques(
            lqr_forces,
            q,
            odrive_message.joint_position,
        )

        # TOOD: set time
        torque_message.timestamp = 0.0
        torque_message.torque = lqr_torques
        publish(lcm, "TORQUE_COMMAND", encode(torque_message))
		t = time_ns()

        if loop_count % 10 == 0
            # println("force: ", lqr_forces)
            # println("torque: ", lqr_torques)
            flush(stdout)
        end
		last_t = t

        # sleep(0.001)
        loop_count += 1
    end
end
