using LinearAlgebra
using StaticArrays

include("../src/Common/QuadrupedDynamics.jl")
include("../src/Common/MPCControl/MPCControl.jl")
include("../src/Common/Utilities.jl")
include("../src/Common/Config.jl")

import .MPCControl

function test()
	param = MPCControl.ControllerParams(Float64, Int64)

	t = 0.0

	joint_pos = @SVector zeros(12)
	joint_vel = @SVector zeros(12)
	x_static = @SVector zeros(12)
	actuator_torques = zeros(12)

	@time MPCControl.control!(
		actuator_torques,
		x_static,
		t,
		joint_pos,
		joint_vel,
		param,
	)
end
