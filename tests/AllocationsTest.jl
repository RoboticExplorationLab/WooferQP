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

	q = zeros(27)
	q[4] = 1.0

	q̇ = zeros(26)

	@time MPCControl.control!(
		actuator_torques,
		q,
		q̇,
		t,
		param
	)
end
