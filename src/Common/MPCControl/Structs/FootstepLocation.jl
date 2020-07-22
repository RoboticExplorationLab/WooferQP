mutable struct FootstepLocation{T} <: AbstractArray{T, 1}
	fr::SVector{3, T}
	fl::SVector{3, T}
	br::SVector{3, T}
	bl::SVector{3, T}

	function FootstepLocation(r::AbstractVector{T}) where T
		@assert length(r) == 12

		fr = SVector{3}(r[1], r[2], r[3])
		fl = SVector{3}(r[4], r[5], r[6])
		br = SVector{3}(r[7], r[8], r[9])
		bl = SVector{3}(r[10], r[11], r[12])

		new{T}(fr, fl, br, bl)
	end
end

function footstep_location_from_angles(α::AbstractVector{T}) where T
	r = ForwardKinematicsAll(α)

	return FootstepLocation(r)
end

function Base.getindex(fsl::FootstepLocation, i::Integer)
	if i==1
		return fsl.fr
	elseif i==2
		return fsl.fl
	elseif i==3
		return fsl.br
	elseif i==4
		return fsl.bl
	else
		throw(BoundsError(fsl, i))
	end
end

function Base.setindex!(fsl::FootstepLocation, loc::SVector{3}, i::Integer)
	if i==1
		fsl.fr = loc
	elseif i==2
		fsl.fl = loc
	elseif i==3
		fsl.br = loc
	elseif i==4
		fsl.bl = loc
	else
		throw(BoundsError(fsl, i))
	end
end

Base.convert(::Type{FootstepLocation}, r::Union{Vector, SVector}) = FootstepLocation(r)
Base.eltype(::Type{FootstepLocation{T}}) where T = SVector{3,T}
Base.length(::FootstepLocation) = 4
Base.size(::FootstepLocation) = (4,)
Base.zero(::Type{FootstepLocation}) = FootstepLocation(@SVector zeros(12))
