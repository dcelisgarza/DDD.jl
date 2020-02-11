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

function limits!(
    lims::Matrix{<:Real},
    segLen::Real=1,
    range::AbstractArray{<:Real} = Float64[-1 -1 -1; 1 1 1].*segLen,
    scale::AbstractArray{<:Real} = Float64[1 1 1; 1 1 1].*segLen,
    buffer::Real = 1.5,
)
    for i = 1:size(lims, 2)
        for j = 1:size(lims, 1)
            lims[j, i] = range[j, i] * scale[j, i] + buffer * segLen
        end
    end
    return lims
end

function scaleCoord(coord::AbstractArray{<:Real}, factor::Real)
    return coord .*= factor
end

function moveCoord(
    coord::AbstractArray{<:Real},
    disp::AbstractVector{<:Real},
    lims::AbstractArray{<:Real},
)
    for i = 1:size(coord, 2)
        for j = 1:size(coord, 1)
            coord[j,i] += lims[1,i] + (lims[2,i] - lims[1,i])*disp[i]
        end
    end
    return coord
end

function makeLoop!(
    loopType::loopPrism, # Loop type
    network::DislocationNetwork,
    dlnParams::DislocationP,
    slipSystems::AbstractArray{<:Real},
    range::AbstractArray{<:Real} = Float64[-1 -1 -1; 1 1 1],
    scale::AbstractArray{<:Real} = Float64[1 1 1; 1 1 1],
    buffer::Real = 1.5,
    nodeLabels::Vector{nodeType} = nodeType[1; 0; 1; 1; 1; 0; 1; 1],
    nodeSide::Integer = 2; # Nodes per side
    dist::Function = f(x...) = Float64[0; 0; 0],
)

    # Aliases/views for dlnParams.
    numSources = dlnParams.numSources
    paramSlipSys = dlnParams.slipSystems
    distSource = dlnParams.distSource



    # Constants
    numNodes = 4 * nodeSide
    numSlipSystem = length(paramSlipSys)
    lenNumSources = length(numSources)
    lenDistSource = length(distSource)
    numLoopsTotal = numSlipSystem * sum(numSources)

    # Error checking.
    @assert numSlipSystem == lenNumSources == lenDistSource
    @assert length(nodeLabels) == numNodes * numSlipSystem * sum(numSources)

    # Allocating memory.
    lenArr = numSlipSystem * numNodes * sum(numSources)
    tmp = findfirst(x -> x == -1, network.label)
    if tmp == nothing
        vcat(network, lenArr)
    elseif length(network.label) - tmp < lenArr
        vcat(network, lenArr - mod(length(network.label) - tmp, lenArr) - 1)
    end

    # Aliases for dlnNetwork.
    links = network.links
    bVec = network.bVec
    slipPlane = network.slipPlane
    coord = network.coord
    label = network.label

    # Local variables for indexing and holding data.
    local seg = zeros(3, 2)
    local lims = zeros(2, 3)
    local numNodesTotal::Integer = 0
    local startIdx = findfirst(x -> x == -1, label)
    local numLoops = 0
    distribution = dist(numLoopsTotal)
    # Loop through slip systems.
    for i = 1:numSlipSystem
        # Make edge and screw segments.
        seg[:, 1] = makeSegment(dlnEdge(), paramSlipSys[i], slipSystems)
        seg[:, 2] = makeSegment(dlnScrew(), paramSlipSys[i], slipSystems)
        # Set limits for dislocations.
        limits!(lims, distSource[i], range, scale, buffer)
        # Set starting index for this slip system.
        idx = startIdx + (i - 1) * numNodesTotal
        idxDist = 1 + (i - 1) * numLoops
        for j = 1:numSources[i]
            # Set index for this source.
            idxj = (j - 1) * numNodes
            idx += idxj

            # First node.
            coord[idx, :] -= seg[:, 1] + seg[:, 2]
            label[idx] = nodeLabels[mod(idx, lenArr)]
            # Prismatic loop necessarily has four sides, do this for each side.
            for k = 1:nodeSide
                coord[idx+k, :] = coord[idx+k-1, :] + seg[:, 1]
                label[idx+k] = nodeLabels[mod(idx + k, lenArr)]
            end
            for k = 1:nodeSide
                coord[idx + k + nodeSide, :] = coord[idx+k-1+nodeSide, :] +
                                               seg[:, 2]
                label[idx + k + nodeSide] = nodeLabels[mod(
                    idx + k,
                    lenArr,
                )+nodeSide]
            end
            for k = 1:nodeSide
                coord[idx + k + 2 * nodeSide, :] = coord[
                    idx+k-1+2*nodeSide,
                    :,
                ] - seg[:, 1]
                label[idx + k + 2 * nodeSide] = nodeLabels[mod(
                    idx + k,
                    lenArr,
                )+2*nodeSide]
            end
            # Last node.
            coord[idx + 1 + 3 * nodeSide, :] = coord[idx+3*nodeSide, :] -
                                               seg[:, 2]
            label[idx + 1 + 3 * nodeSide] = nodeLabels[mod(
                idx + 1,
                lenArr,
            )+3*nodeSide]

            # Scale loop.
            coord[idx:idx+1+3 * nodeSide, :] .= scaleCoord(
                coord[idx:idx+1+3 * nodeSide, :],
                distSource[i],
            )

            # Move loop.
            coord[idx:idx+1+3 * nodeSide, :] .= moveCoord(
                coord[idx:idx+1+3 * nodeSide, :],
                distribution[idxDist, :], lims
            )


            # Increase node total for stride and knowing how many nodes a network has.
            numNodesTotal += numNodes
            numLoops += 1

            # Generate links for this loop.
            for j = 1:nodeSide*4-1
                links[idx-1+j, :] = [idx - 1 + j; idx + j]
            end
            links[idx-1+nodeSide*4, :] = [idx - 1 + numNodes; idx]
            bVec[idx:idx-1+nodeSide*4, :] = reshape(
                repeat(slipSystems[paramSlipSys[i], 4:6], inner = numNodes),
                numNodes,
                3,
            )
            slipPlane[idx:idx-1+nodeSide*4, :] = reshape(
                repeat(slipSystems[paramSlipSys[i], 1:3], inner = numNodes),
                numNodes,
                3,
            )
        end
    end
    # Track number of nodes and links.
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
