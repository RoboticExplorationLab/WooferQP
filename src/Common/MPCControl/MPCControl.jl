module MPCControl

export
	control!

using LinearAlgebra
using StaticArrays
using Rotations
using ..QuadrupedDynamics
using ForwardDiff
import YAML

include("../Config.jl")
include("../Quaternions.jl")
include("../Utilities.jl")

include("find_solver.jl") # uses MPC.yaml to determine which solver is being used and import correct files and packages
include("control.jl")
include("swing_leg.jl")
include("gait.jl")
include("footsteps.jl")

end
