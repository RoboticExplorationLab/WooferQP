using LinearAlgebra
using StaticArrays
include("Rotations.jl")
include("Utilities.jl")
include("Dynamics.jl")

@with_kw mutable struct LQRParams
    L::SMatrix{12,12, Float64, 144}
    V::SMatrix{12,12, Float64, 144}
    u0::SVector{12, Float64}
    x0::SVector{12, Float64}
    dt::Float64
    A_d::SMatrix{12,12, Float64, 144}
    B_d::SMatrix{12,12, Float64, 144}
    Q::SMatrix{12,12, Float64, 144}
    R::SMatrix{12,12, Float64, 144}
end

function initLQRParams(;dt::AbstractFloat, x0::SVector{12}, rxy_penalty=1e1, rz_penalty=1e2, v_penalty=1, phi_penalty=1e1, w_penalty=1, fxy_penalty=1e-3, fz_penalty=1e-3)::LQRParams
    r1 = ForwardKinematics(@SVector(zeros(3)), 1)
    r2 = ForwardKinematics(@SVector(zeros(3)), 2)
    r3 = ForwardKinematics(@SVector(zeros(3)), 3)
    r4 = ForwardKinematics(@SVector(zeros(3)), 4)
    foot_locs = [r1; r2; r3; r4]

    # linearization point
    x_ϕ = SVector{3}(x0[SVector{3}(4:6)])
    x_q = reconstructq(x_ϕ)

    dvdv = I
    dwdw =  0.5 * Vmat() * Lmult(reconstructq(x_ϕ)) * Vmat()'

    A_c = [zeros(3,6) I zeros(3,3);
           zeros(3,9) dwdw;
           zeros(6, 12)]

    B_c = zeros(12, 12)

    for j = 1:4
        B_c[7:9, LegIndexToRange(j)] .= 1 / WOOFER_CONFIG.SPRUNG_MASS *
                                             Matrix{Float64}(I, 3, 3)

        r_hat = skew(foot_locs[SLegIndexToRange(j)])
        B_c[10:12, SLegIndexToRange(j)] .= inv(WOOFER_CONFIG.INERTIA) * r_hat * rotmat(inv(x_q))
    end

    # First order approximation of zero order hold
    A_d = I + A_c * dt
    B_d = A_d * B_c * dt

    # Discretization via the matrix exponential
    # cont_sys = [A_c B_c;
    #             zeros(12,12) zeros(12,12)]
    # disc_sys = exp(cont_sys * dt)
    # A_d = disc_sys[1:12, 1:12]
    # B_d = disc_sys[1:12, 13:24]

    Q = Diagonal(SVector{12}([  rxy_penalty*ones(2);
                                rz_penalty;
                                phi_penalty*ones(3);
                                v_penalty*ones(3);
                                w_penalty*ones(3)]))
    R = Diagonal(SVector{12}(repeat([fxy_penalty, fxy_penalty, fz_penalty], 4)))
    L = SMatrix{12,12}(dlqr(A_d, B_d, Q, R))
    V = SMatrix{12,12}(dare(A_d, B_d, Q, R))
    # u0 = SVector{12}([0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0] * WOOFER_CONFIG.MASS * 9.81 / 4)

    g = [0, 0, 0, 0, 0, 0, 0, 0, -9.81, 0, 0, 0]*dt
    u0 = SVector{12}(pinv(B_d)*((I - A_d)*x0 - g))

    @show u0

    return LQRParams(
        L = L,
        x0 = x0,
        u0 = u0,
        dt = dt,
        A_d = A_d,
        B_d = B_d,
        Q = Q,
        R = R,
        V = V,
    )
end

function LQRBalance(
    x::SVector{12},
    lqr_params::LQRParams,
)
    forces = lqr_params.u0 - lqr_params.L * (x - lqr_params.x0)

    # make sure feet can't pull on the ground
    forces = setindex(forces, clamp(forces[3], 0.0, Inf), 3)
    forces = setindex(forces, clamp(forces[6], 0.0, Inf), 6)
    forces = setindex(forces, clamp(forces[9], 0.0, Inf), 9)
    forces = setindex(forces, clamp(forces[12], 0.0, Inf), 12)

    return forces
end

function InertialLegForcesToActuatorTorques(forces::SVector{12, T}, q::UnitQuaternion, joint_positions::SVector{12, T})::SVector{12, T} where T
    body_forces = [ inv(q) * forces[SLegIndexToRange(1)];
                    inv(q) * forces[SLegIndexToRange(2)];
                    inv(q) * forces[SLegIndexToRange(3)];
                    inv(q) * forces[SLegIndexToRange(4)]]
    return Force2Torque(-body_forces, joint_positions)
end
