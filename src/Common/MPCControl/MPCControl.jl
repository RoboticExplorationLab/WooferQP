module MPCControl

export
	control!

using LinearAlgebra
using Parametron
using OSQP
using StaticArrays
using ..QuadrupedDynamics
import YAML

include("../Config.jl")
include("../Quaternions.jl")
include("../Utilities.jl")

include("Structs/OptimizerParams.jl")
include("Structs/SwingLegParams.jl")
include("Structs/GaitParams.jl")
include("Structs/ControllerParams.jl")
include("control.jl")
include("osqp_solver.jl")
include("swing_leg.jl")
include("gait.jl")
include("footsteps.jl")

end
