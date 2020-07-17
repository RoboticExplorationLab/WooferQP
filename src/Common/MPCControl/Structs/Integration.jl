abstract type Euler <: TO.RobotDynamics.Explicit end

function TO.RobotDynamics.discrete_dynamics(::Type{Euler}, model::TO.RobotDynamics.AbstractModel, x::StaticVector, u::StaticVector, t, dt)
	k1 = TO.RobotDynamics.dynamics(model, x, u, t)*dt
	x + k1
end

abstract type Exponential <: TO.RobotDynamics.Explicit end

# function TO.RobotDynamics.discrete_dynamics(::Type{Exponential}, model::TO.RobotDynamics.LinearTimeInvariantModel, x::StaticVector, u::StaticVector, t, dt)
#
# end
#
# function TO.RobotDynamics.discrete_dynamics(::Type{Exponential}, model::TO.RobotDynamics.LinearTimeVaryingModel, x::StaticVector, u::StaticVector, t, dt)
#
# end
#
# function TO.RobotDynamics.discrete_dynamics(::Type{Exponential}, model::TO.RobotDynamics.AffineTimeInvariantModel, x::StaticVector, u::StaticVector, t, dt)
#
# end
