struct FrictionConstraint{T} <: TrajectoryOptimization.AbstractConstraint{Inequality,Control,4}
	m::Int
	μ::T
	function FrictionConstraint(m::Int, μ::T) where T
		new{T}(m, μ)
	end
end

@inline control_dims(con::FrictionConstraint) = con.m
@inline sense(con::FrictionConstraint) = Inequality

function evaluate(con::FrictionConstraint, u::SVector)
	return @SVector [	norm(u[1:2], 2) - con.μ*u[3],
						norm(u[4:5], 2) - con.μ*u[6],
						norm(u[7:8], 2) - con.μ*u[9],
						norm(u[10:11], 2) - con.μ*u[12]	]
end

# TODO: non-differentiability here, use subgradient when at zero?
# function jacobian(con::FrictionConstraint, u::SVector)
# 	return @SVector [	u[1:2]'/,
# 						]
# end
