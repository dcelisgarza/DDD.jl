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
    @views slipPlane = slipSystems[slipSys, 1:3]
    @views bVec = slipSystems[slipSys, 4:6]
    edgeSeg = cross(slipPlane, bVec)
    edgeSeg ./= norm(edgeSeg)
    return edgeSeg
end
function makeSegment(
    segType::dlnScrew,
    slipSys::Integer,
    slipSystems::Matrix{<:Real},
)
    @views bVec = slipSystems[slipSys, 4:6]
    screwSeg = bVec ./ norm(bVec)
    return screwSeg
end

function makeLoop!(
    loopType::loopPrism, # Loop type
    network::DislocationNetwork,
    dlnParams::DislocationP,
    slipSystems::Matrix{<:Real},
    nodeLabels::Vector{nodeType} = nodeType[1; 0; 1; 1; 1; 0; 1; 1],
    nodeSide::Integer = 2, # Nodes per side
)

    # Set up aliases/views.
    numSources = dlnParams.numSources
    paramSlipSys = dlnParams.slipSystems

    # Constants
    numNodes = 4 * nodeSide
    numSlipSystem = length(paramSlipSys)
    numSlipSystem == 0 ? numSlipSystem = 1 : nothing
    lenNumSources = length(numSources)
    lenNumSources == 0 ? numSlipSystem = 1 : nothing
    @assert numSlipSystem == lenNumSources
    @assert length(nodeLabels) == numNodes * numSlipSystem * sum(numSources)

    if network.numNode == 0
        network.links = zeros(numSlipSystem * numNodes, 2)
        network.bVec = zeros(numSlipSystem * numNodes, 3)
        network.slipPlane = zeros(numSlipSystem * numNodes, 3)
        network.coord = zeros(numSlipSystem * numNodes, 3)
        network.label = zeros(nodeType, numSlipSystem * numNodes)
    else
        vcat(network, numNodes * numSlipSystem * sum(numSources))
    end

    links = network.links
    bVec = network.bVec
    slipPlane = network.slipPlane
    coord = network.coord
    label = network.label

    local seg = zeros(3, 2)
    local numNodesTotal::Integer = 0
    local startIdx = findfirst(label .== -1)

    for i = 1:numSlipSystem
        seg[:, 1] = makeSegment(dlnEdge(), paramSlipSys[i], slipSystems)
        seg[:, 2] = makeSegment(dlnScrew(), paramSlipSys[i], slipSystems)

        idx = startIdx + (i - 1) * numNodesTotal
        for j = 1:numSources[i]
            idxj = (j - 1) * nodeSide * 4
            idx += idxj
            coord[idx, :] -= seg[:, 1] + seg[:, 2]
            label[idx] = nodeLabels[mod(idx, nodeSide*4)]
            # Prismatic loop has four sides, do this for each side.
            for k = 1:nodeSide
                coord[idx+k, :] = coord[idx+k-1, :] + seg[:, 1]
                label[idx+k] = nodeLabels[mod(idx+k, nodeSide*4)]
            end
            for k = 1:nodeSide
                coord[idx + k + nodeSide, :] = coord[idx+k-1+nodeSide, :] +
                                               seg[:, 2]
                label[idx + k + nodeSide] = nodeLabels[mod(idx + k, nodeSide*4) + nodeSide]
            end
            for k = 1:nodeSide
                coord[idx + k + 2 * nodeSide, :] = coord[
                    idx+k-1+2*nodeSide,
                    :,
                ] - seg[:, 1]
                label[idx + k + 2 * nodeSide] = nodeLabels[mod(idx + k, nodeSide*4) +
                                                           2 * nodeSide]
            end
            coord[idx + 1 + 3 * nodeSide, :] = coord[idx+3*nodeSide, :] -
                                               seg[:, 2]
            label[idx + 1 + 3 * nodeSide] = nodeLabels[mod(idx + 1, nodeSide*4) + 3 * nodeSide]

            numNodesTotal += nodeSide * 4

            for j = 1:nodeSide*4-1
                links[idx-1+j, :] = [idx - 1 + j; idx + j]
            end
            links[idx-1+nodeSide*4, :] = [idx - 1 + nodeSide * 4; idx]
            bVec[idx:idx-1+nodeSide*4, :] = reshape(
                repeat(slipSystems[paramSlipSys[i], 4:6], inner = nodeSide * 4),
                nodeSide * 4,
                3,
            )
            slipPlane[idx:idx-1+nodeSide*4, :] = reshape(
                repeat(slipSystems[paramSlipSys[i], 1:3], inner = nodeSide * 4),
                nodeSide * 4,
                3,
            )
        end
    end

    network.numNode += numNodesTotal
    network.numSeg += numNodesTotal
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
