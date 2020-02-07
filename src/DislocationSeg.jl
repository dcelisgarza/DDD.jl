abstract type dlnSegment end
struct dlnEdge <: dlnSegment end
struct dlnScrew <: dlnSegment end
struct dlnMixed <: dlnSegment end

"""
```
makeSegment(segType::dlnEdge, slipSys::Integer, data::Matrix{<:Real})
makeSegment(segType::dlnScrew, slipSys::Integer, data::Matrix{<:Real})
```
Make single dislocation segment of the given type `dlnEdge` or `dlnScrew` from the slip system corresponding to `data[:, slipSys]`.
"""
function makeSegment(segType::dlnEdge, slipSys::Integer, data::Matrix{<:Real})
    @views slipPlane = data[1:3, slipSys]
    @views bVec = data[4:6, slipSys]
    edgeSeg = cross(slipPlane, bVec)
    edgeSeg ./= norm(edgeSeg)
    return edgeSeg
end
function makeSegment(segType::dlnScrew, slipSys::Integer, data::Matrix{<:Real})
    @views bVec = data[4:6, slipSys]
    screwSeg = bVec ./ norm(bVec)
    return screwSeg
end
