function LegIndexToRange(i::Int)
    return 3*(i-1)+1:3*(i-1)+3
end

function SLegIndexToRange(i::Int)
    return SVector{3}(3*(i-1)+1:3*(i-1)+3)
end

function Vec12ToSVector(q::Vector)
    @assert length(q) == 12
    return SVector{12}(
        q[1],
        q[2],
        q[3],
        q[4],
        q[5],
        q[6],
        q[7],
        q[8],
        q[9],
        q[10],
        q[11],
        q[12],
    )
end

function Vec3ToSVector(q::Vector)
    @assert length(q) == 3
    return SVector{3}(q[1], q[2], q[3])
end

function Vec4ToSVector(q::Vector)
    @assert length(q) == 4
    return SVector{4}(q[1], q[2], q[3], q[4])
end
