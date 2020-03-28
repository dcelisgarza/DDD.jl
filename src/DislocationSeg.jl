function limits!(
    lims::AbstractArray{<:Float64, N1},
    segLen::Float64,
    range::AbstractArray{<:Float64, N2},
    buffer::Float64,
) where {N1, N2}
    @inbounds for i = 1:size(lims, 2)
        for j = 1:size(lims, 1)
            lims[j, i] = range[j, i] + buffer * segLen
        end
    end
    return lims
end

function translatePoints(
    coord::AbstractArray{<:Float64, N1},
    lims::AbstractArray{<:Float64, N1},
    disp::AbstractArray{<:Float64, N2},
) where {N1, N2}
    @inbounds for i = 1:size(coord, 2)
        for j = 1:size(coord, 1)
            coord[j, i] += lims[1, i] + (lims[2, i] - lims[1, i]) * disp[i]
        end
    end
    return coord
end

function loopDistribution(dist::Zeros, n::Integer, args...)
    return zeros(n, 3)
end
function loopDistribution(dist::Rand, n::Integer, args...)
    return rand(n, 3)
end
function loopDistribution(dist::Randn, n::Integer, args...)
    return randn(n, 3)
end
function loopDistribution(dist::Regular, n::Integer, args...)
    error("loopDistribution: regular distribution yet not implemented")
end

function makeConnect(
    links::AbstractArray{Int64, N},
    label::Vector{nodeType},
    maxConnect::Integer = 4,
) where {N}
    iLnk = findall(x -> x != 0, links[:, 1]) # Indices of defined links.
    lenLabel = length(label)
    connectivity = zeros(Int64, lenLabel, 1 + 2 * maxConnect)
    linksConnect = zeros(Int64, size(links, 1), 2)
    @inbounds for i in eachindex(iLnk)
        idx = iLnk[i] # Index of links that contain non zero links.
        # links[idx, :] yields the nodes involved in the link
        n1 = links[idx, 1] # Node 1, it is the row of the coord matrix
        n2 = links[idx, 2] # Node 2

        connectivity[n1, 1] += 1 # Number of other nodes node n1 is connected to.
        connectivity[n2, 1] += 1 # Number of other nodes node n2 is connected to.

        tmp1 = 2 * connectivity[n1, 1]
        tmp2 = 2 * connectivity[n2, 1]

        connectivity[n1, tmp1:(tmp1 + 1)] = [idx, 1] # idx = linkID, 1 = first node in link with linkID
        connectivity[n2, tmp2:(tmp2 + 1)] = [idx, 2] # idx = linkID, 2 = second node in link with linkID

        linksConnect[idx, 1] = connectivity[n1, 1]
        linksConnect[idx, 2] = connectivity[n2, 1]
    end
    return connectivity, linksConnect
end

function makeNetwork!(
    network::DislocationNetwork,
    sources::Union{DislocationLoop, AbstractVector{<:DislocationLoop}},
    maxConnect::Integer = 4,
)
    local nodeTotal::Integer = 0
    local lims = zeros(Float64, 2, 3)
    # Allocate memory.
    @inbounds for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end
    available = findfirst(x -> x == -1, network.label)
    if available == nothing
        malloc(network, nodeTotal)
    else
        check = length(network.label) - (available - 1) - nodeTotal
        check < 0 ? malloc(network, -check) : nothing
    end

    nodeTotal = 0
    initIdx = findfirst(x -> x == -1, network.label)
    initIdx == nothing ? initIdx = 0 : nothing
    @inbounds for i in eachindex(sources)
        # Indices.
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        disp = loopDistribution(sources[i].dist, numLoops)
        limits!(
            lims,
            mean(sources[i].segLen),
            sources[i].range,
            sources[i].buffer,
        )
        for j = 1:numLoops
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            # Prepare to distribute sources.
            network.links[idxi:idxf, :] =
                sources[i].links[1:nodesLoop, :] .+ (nodeTotal + initIdx - 1)
            network.slipPlane[idxi:idxf, :] .=
                sources[i].slipPlane[1:nodesLoop, :]
            network.bVec[idxi:idxf, :] .= sources[i].bVec[1:nodesLoop, :]
            network.coord[idxi:idxf, :] .= sources[i].coord[1:nodesLoop, :]
            # Move loop.
            network.coord[idxi:idxf, :] .= translatePoints(
                sources[i].coord[1:nodesLoop, :],
                lims,
                disp[j, :],
            )
            network.label[idxi:idxf, :] .= sources[i].label[1:nodesLoop, :]
            nodeTotal += nodesLoop
        end
    end

    network.numNode += nodeTotal
    network.numSeg += nodeTotal

    network.connectivity, network.linksConnect =
        makeConnect(network.links, network.label, maxConnect)
    return network
end

function makeNetwork(
    sources::Union{DislocationLoop, AbstractVector{<:DislocationLoop}},
    maxConnect::Integer = 4,
    memBuffer::Integer = 10,
)
    local nodeTotal::Integer = 0
    local lims = zeros(Float64, 2, 3)
    # Allocate memory.
    @inbounds for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end

    local nodeBuffer::Integer = nodeTotal * memBuffer

    network = DislocationNetwork(
        zeros(Int64, nodeBuffer, 2),
        zeros(nodeBuffer, 3),
        zeros(nodeBuffer, 3),
        zeros(nodeBuffer, 3),
        zeros(nodeType, nodeBuffer),
        nodeTotal,
        nodeTotal,
    )

    nodeTotal = 0
    initIdx = findfirst(x -> x == -1, network.label)
    initIdx == nothing ? initIdx = 0 : nothing
    @inbounds for i in eachindex(sources)
        # Indices.
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        disp = loopDistribution(sources[i].dist, numLoops)
        limits!(
            lims,
            mean(sources[i].segLen),
            sources[i].range,
            sources[i].buffer,
        )
        for j = 1:numLoops
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            # Prepare to distribute sources.
            network.links[idxi:idxf, :] =
                sources[i].links[1:nodesLoop, :] .+ (nodeTotal + initIdx - 1)
            network.slipPlane[idxi:idxf, :] .=
                sources[i].slipPlane[1:nodesLoop, :]
            network.bVec[idxi:idxf, :] .= sources[i].bVec[1:nodesLoop, :]
            network.coord[idxi:idxf, :] .= sources[i].coord[1:nodesLoop, :]
            # Move loop.
            network.coord[idxi:idxf, :] .= translatePoints(
                sources[i].coord[1:nodesLoop, :],
                lims,
                disp[j, :],
            )
            network.label[idxi:idxf, :] .= sources[i].label[1:nodesLoop, :]
            nodeTotal += nodesLoop
        end
    end

    network.connectivity, network.linksConnect =
        makeConnect(network.links, network.label, maxConnect)
    return network
end
