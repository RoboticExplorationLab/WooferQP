using LinearAlgebra
using StaticArrays

include("../src/Common/QuadrupedDynamics.jl")
include("../src/Common/MPCControl/MPCControl.jl")
include("../src/Common/Utilities.jl")
include("../src/Common/Config.jl")

using .QuadrupedDynamics
import .MPCControl

param = MPCControl.ControllerParams(Float64, Int64)

x_est = [0.0, 0.0, 0.28, 0.1, 0.0, 0.2, -0.01, 0.01, 0.05, 0.01, 0.02, 0.03]

MPCControl.reference_trajectory!(x_est, param)
t = 0.0
MPCControl.foot_history!(t, param)
MPCControl.foot_forces!(x_est, param)

param.forces
