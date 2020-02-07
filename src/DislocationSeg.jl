abstract type AbstractDlnSegment end
struct dlnEdge <: AbstractDlnSegment end
struct dlnScrew <: AbstractDlnSegment end
struct dlnMixed <: AbstractDlnSegment end

abstract type AbstractDlnLoop end
struct loopPrism <: AbstractDlnLoop end
struct loopShear <: AbstractDlnLoop end
struct loopDln <: AbstractDlnLoop end

"""
```
makeSegment(segType::dlnEdge, slipSys::Integer, data::Matrix{<:Real})
makeSegment(segType::dlnScrew, slipSys::Integer, data::Matrix{<:Real})
```
Make single dislocation segment of the given type `dlnEdge` or `dlnScrew` from the slip system corresponding to `slipSystems[:, slipSys]`.
"""
function makeSegment(
    segType::dlnEdge,
    slipSys::Integer,
    slipSystems::Matrix{<:Real},
)
    @views slipPlane = slipSystems[1:3, slipSys]
    @views bVec = slipSystems[4:6, slipSys]
    edgeSeg = cross(slipPlane, bVec)
    edgeSeg ./= norm(edgeSeg)
    return edgeSeg
end
function makeSegment(
    segType::dlnScrew,
    slipSys::Integer,
    slipSystems::Matrix{<:Real},
)
    @views bVec = slipSystems[4:6, slipSys]
    screwSeg = bVec ./ norm(bVec)
    return screwSeg
end

function makeLoop!(
    loopType::loopPrism, # Loop type
    network::DislocationNetwork,
    dlnParams::DislocationP,
    slipSystems::Matrix{<:Real},
    nodeSide::Integer = 2, # Nodes per side
)

    paramSlipSys = dlnParams.slipSystems
    numNodes = 4 * nodeSide
    numSlipSystem = length(paramSlipSys)
    numSlipSystem == 0 ? numSlipSystem = 1 : nothing
    lenNumSources = length(dlnParams.numSources)
    lenNumSources == 0 ? numSlipSystem = 1 : nothing
    @assert numSlipSystem == lenNumSources
    if size(network.coord, 1) == 0 || size(network.label, 1) == 0
        network.coord = zeros(numSlipSystem * numNodes, 3)
        network.label = zeros(nodeType, numSlipSystem * numNodes)
    end

    local seg = zeros(3, 2)

    coord = @view network.coord[:, :]
    for i = 1:numSlipSystem
        seg[:, 1] = makeSegment(dlnEdge(), paramSlipSys[i], slipSystems)
        seg[:, 2] = makeSegment(dlnScrew(), paramSlipSys[i], slipSystems)

        idxi = (paramSlipSys[i] - 1) * numNodes
        coord[1+idxi, :] -= seg[:, 1] + seg[:, 2]
        for j = 1:nodeSide
            coord[j + 1 + idxi, :] += seg[:, 1]
            coord[j + 1 + idxi + nodeSide, :] += seg[:, 2]
            coord[j + 1 + idxi + 2 * nodeSide, :] -= seg[:, 1]
            try
                coord[j + 1 + idxi + 3 * nodeSide, :] -= seg[:, 2]
            catch e
            end
        end
    end

    return network
end

# DislocationNetwork(links, bVec, slipPlane, coord, label, numNode = 0,
#     numSeg = 0,
# )

# DislocationP(coreRad, coreRadMag, minSegLen, maxSegLen, minArea, maxArea,
#     maxConnect, remesh, collision, separation, virtualRemesh, edgeDrag,
#     screwDrag, climbDrag, lineDrag, mobility, numSources = 0,
#     slipSystems = 0, distSource = 0.0,
# )
