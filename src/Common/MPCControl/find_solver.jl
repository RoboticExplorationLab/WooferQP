const data = YAML.load(open(joinpath(@__DIR__, "MPC.yaml")))

if data["solver"] == "ALTRO"
	const using_altro = true

	using TrajectoryOptimization
	using Altro

	const TO = TrajectoryOptimization

	# include("Structs/LinearModels.jl")
	include("Structs/Integration.jl")
	include("Structs/QuadrupedModel.jl")
	include("Structs/ALTROParams.jl")
	include("Structs/LinearizedFrictionConstraint.jl")
	include("Structs/SwingLegParams.jl")
	include("Structs/GaitParams.jl")
	include("Structs/ControllerParams.jl")

	include("altro_solver.jl")

else
	const using_altro = false

	using Parametron
	using OSQP

	include("Structs/QPParams.jl")
	include("Structs/SwingLegParams.jl")
	include("Structs/GaitParams.jl")
	include("Structs/ControllerParams.jl")

	include("osqp_solver.jl")
end
