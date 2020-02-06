abstract type dlnSegment end
struct dlnEdge <: dlnSegment end
struct dlnScrew <: dlnSegment end
struct dlnMixed <: dlnSegment end

# Edge
function makeSegment(segType::dlnEdge, slipSys::Integer, data::Matrix{<:Real})
    @views slipPlane = data[1:3, slipSys]
    @views bVec = data[4:6, slipSys]

    edgeSeg = cross(slipPlane, bVec)
    edgeSeg ./= norm(edgeSeg)

    return edgeSeg
end

# Screw
function makeSegment(segType::dlnScrew, slipSys::Integer, data::Matrix{<:Real})
    @views bVec = data[4:6, slipSys]
    screwSeg = bVec ./ norm(bVec)
    return screwSeg
end
