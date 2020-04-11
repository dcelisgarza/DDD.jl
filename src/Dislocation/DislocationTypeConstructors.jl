# Zero constructor.
function makeLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen::T3,
    slipSystem::T2,
    _slipPlane::T4,
    _bVec::T4,
    label::T5,
    buffer::T6,
    range::T4,
    dist::T7,
) where {
    T1 <: loopDln,
    T2 <: Int64,
    T3 <: Float64,
    T4 <: AbstractArray{<:Float64, N} where {N},
    T5 <: Vector{nodeType},
    T6 <: Float64,
    T7 <: AbstractDistribution,
}

    nodeTotal = 0
    links = zeros(Int64, nodeTotal, 2)
    coord = zeros(nodeTotal, 3)
    slipPlane = zeros(0, 3)
    bVec = zeros(0, 3)

    return (
        numSides,
        nodeSide,
        numLoops,
        segLen,
        slipSystem,
        links,
        slipPlane,
        bVec,
        coord,
        label,
        buffer,
        range,
        dist,
    )
end

"""
```
makeLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen::T3,
    slipSystem::T2,
    _slipPlane::T4,
    _bVec::T4,
    label::T5,
    buffer::T6,
    range::T7,
    dist::T8,
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T4 <: AbstractArray{<:Float64, N} where {N},
    T5 <: Vector{nodeType},
    T6 <: Float64,
    T7 <: AbstractArray{<:Float64, N} where {N},
    T8 <: AbstractDistribution,
}
```
Constructor function for [`DislocationLoop`](@ref). See [`AbstractDlnStr`](@ref), [`AbstractDlnSeg`](@ref), [`nodeType`](@ref), [`AbstractDistribution`](@ref) for further details.
"""
function makeLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen::T3,
    slipSystem::T2,
    _slipPlane::T4,
    _bVec::T4,
    label::T5,
    buffer::T6,
    range::T7,
    dist::T8,
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T4 <: AbstractArray{<:Float64, N} where {N},
    T5 <: Vector{nodeType},
    T6 <: Float64,
    T7 <: AbstractArray{<:Float64, N} where {N},
    T8 <: AbstractDistribution,
}

    nodeTotal = numSides * nodeSide
    lSegLen = length(segLen)
    @assert length(label) == nodeTotal "makeLoop: All $nodeTotal nodes must be labelled. There are only $(length(label)) labels currently defined."
    @assert lSegLen == nodeTotal "makeLoop: All $nodeTotal segments must have their lengths defined. There are only $lSegLen lengths currently defined."

    _slipPlane = _slipPlane ./ norm(_slipPlane)
    _bVec = _bVec ./ norm(_bVec)

    rotAxis = zeros(eltype(_slipPlane), 3)
    if typeof(loopType) == loopShear
        rotAxis = _slipPlane
    elseif typeof(loopType) == loopPrism
        rotAxis = _bVec
    end

    links = zeros(Int64, nodeTotal, 2)
    coord = zeros(nodeTotal, 3)
    slipPlane = zeros(0, 3)
    bVec = zeros(0, 3)
    seg = zeros(lSegLen, 3)
    @inbounds for i in eachindex(segLen)
        seg[i, :] = makeSegment(segEdge(), _slipPlane, _bVec) .* segLen[i]
    end

    θ = extAngle(numSides)
    rseg = zeros(3)
    @inbounds for i = 1:numSides
        idx = (i - 1) * nodeSide
        rseg = rot3D(
            seg[mod(i - 1, lSegLen) + 1, :],
            rotAxis,
            zeros(3),
            θ * (i - 1),
        )
        for j = 1:nodeSide
            if i == j == 1
                coord[1, :] = zeros(3)
                slipPlane = [slipPlane; _slipPlane']
                bVec = [bVec; _bVec']
                continue
            end
            if idx + j <= nodeTotal
                coord[idx + j, :] += coord[idx + j - 1, :] + rseg
                slipPlane = [slipPlane; _slipPlane']
                bVec = [bVec; _bVec']
            end
        end
    end
    coord .-= mean(coord, dims = 1)

    # Links
    @inbounds for j = 1:(nodeTotal - 1)
        links[j, :] = [j; 1 + j]
    end
    links[nodeTotal, :] = [nodeTotal; 1]

    return (
        numSides,
        nodeSide,
        numLoops,
        segLen,
        slipSystem,
        links,
        slipPlane,
        bVec,
        coord,
        label,
        buffer,
        range,
        dist,
    )
end

"""
```
makeNetwork(
    sources::Union{
            DislocationLoop,
            AbstractVector{<:DislocationLoop}
        }, # Dislocation structures.
    maxConnect::Integer = 4,
    memBuffer::Integer = 10, # Buffer for memory allocation. The code will allocate the total number of nodes times `memBuffer` to reduce dynamic memory allocation during runtime.
    args...;
    checkConsistency::Bool = false, # Check the consistency of the network.
    kw...,
)
```
Constructor for [`DislocationNetwork`](@ref), see [`makeNetwork!`](@ref) for in-place version.
"""
function makeNetwork(
    sources::Union{DislocationLoop, AbstractVector{<:DislocationLoop}},
    maxConnect::Integer = 4,
    args...;
    memBuffer::Integer = 10,
    checkConsistency::Bool = true,
    kw...,
)
    nodeTotal::Integer = 0
    lims = zeros(Float64, 2, 3)
    # Allocate memory.
    @inbounds for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end

    nodeBuffer::Integer = nodeTotal * memBuffer

    network = DislocationNetwork(
        links = zeros(Int64, nodeBuffer, 2),
        slipPlane = zeros(nodeBuffer, 3),
        bVec = zeros(nodeBuffer, 3),
        coord = zeros(nodeBuffer, 3),
        label = zeros(nodeType, nodeBuffer),
        numNode = nodeTotal,
        numSeg = nodeTotal,
    )
    network.maxConnect = maxConnect

    nodeTotal = 0
    initIdx = findfirst(x -> x == -1, network.label)
    initIdx == nothing ? initIdx = 0 : nothing
    @inbounds for i in eachindex(sources)
        # Indices.
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        disp = loopDistribution(sources[i].dist, numLoops, args...; kw...)
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

    getSegmentIdx!(network)
    makeConnect!(network)

    checkConsistency ? checkNetwork(network) : nothing

    return network
end

"""
```
makeNetwork!(
    network::DislocationNetwork,
    sources::Union{
        DislocationLoop,
        AbstractVector{<:DislocationLoop}
    },
    maxConnect::Integer = 4,
    args...;
    checkConsistency::Bool = false,
    kw...,
)
```
In-place constructor for [`DislocationNetwork`](@ref), see [`makeNetwork`](@ref) for constructor.
"""
function makeNetwork!(
    network::DislocationNetwork,
    sources::Union{DislocationLoop, AbstractVector{<:DislocationLoop}},
    maxConnect::Integer = 4,
    args...;
    checkConsistency::Bool = false,
    kw...,
)
    nodeTotal::Integer = 0
    lims = zeros(Float64, 2, 3)
    network.maxConnect = maxConnect
    # Allocate memory.
    @inbounds for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end
    available = findfirst(x -> x == -1, network.label)
    if available == nothing
        push!(network, nodeTotal)
    else
        check = length(network.label) - (available - 1) - nodeTotal
        check < 0 ? push!(network, -check) : nothing
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
        disp = loopDistribution(sources[i].dist, numLoops, args...; kw...)
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

    getSegmentIdx!(network)
    makeConnect!(network)

    checkConsistency ? checkNetwork(network) : nothing

    return network
end

"""
```
loopDistribution(dist<:AbstractDistribution, n::Integer, args...; kw...)
```
Returns `n` points according to the concrete subtype of [`AbstractDistribution`](@ref) given. Overload this function with new concrete subtypes and custom distributions. This and [`limits!`](@ref) are used in [`translatePoints`](@ref) to distribute dislocations in the simulation domain.
"""
function loopDistribution(dist::Zeros, n::Integer, args...; kw...)
    return zeros(n, 3)
end
function loopDistribution(dist::Rand, n::Integer, args...; kw...)
    return rand(n, 3)
end
function loopDistribution(dist::Randn, n::Integer, args...; kw...)
    return randn(n, 3)
end
function loopDistribution(dist::Regular, n::Integer, args...; kw...)
    error("loopDistribution: regular distribution yet not implemented")
end

"""
```
limits!(
    lims::AbstractArray{<:Float64, N1},
    segLen::Float64,
    range::AbstractArray{<:Float64, N2},
    buffer::Float64,
) where {N1, N2}
```
Calculate the spatial limits a dislocation will occupy. This and [`loopDistribution`](@ref) are used in [`translatePoints`](@ref) to distribute dislocations in the simulation domain.
"""
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

"""
```
function translatePoints(
    coord::AbstractArray{<:Float64, N1},
    lims::AbstractArray{<:Float64, N1},
    disp::AbstractArray{<:Float64, N2},
) where {N1, N2}
```
Translate dislocation node coordinates `coord` inside the spatial bounds of `lims` (calculated in [`limits!`](@ref)) according to the displacement `disp` (calculated in [`loopDistribution`](@ref)). Used to distribute sources inside a domain in [`makeNetwork`](@ref) and [`makeNetwork!`](@ref).
"""
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

"""
```
makeConnect!(network::DislocationNetwork)
```
Creates `connectivity` and `linksConnect` matrices in-place. `connectivity` contains the number of other other nodes each node is connected to, up to `maxConnect` other nodes. It also contains the segments in which it's involved. `linksConnect` is the running total of the number of other nodes each node is connected to, will probably be deleted in the future. This is called from [`makeNetwork`](@ref) and [`makeNetwork!`](@ref).
"""
function makeConnect!(network::DislocationNetwork)
    links = network.links
    maxConnect = network.maxConnect

    idx = findall(x -> x != 0, links[:, 1]) # Indices of defined links.
    lenLinks = size(links, 1)
    connectivity = zeros(Int64, lenLinks, 1 + 2 * maxConnect)
    linksConnect = zeros(Int64, lenLinks, 2)
    @inbounds for i in idx
        # links[idx, :] yields the nodes involved in the link
        n1 = links[i, 1] # Node 1, it is the row of the coord matrix
        n2 = links[i, 2] # Node 2

        connectivity[n1, 1] += 1 # Num of other nodes node n1 is connected to.
        connectivity[n2, 1] += 1 # Num of other nodes node n2 is connected to.

        tmp1 = 2 * connectivity[n1, 1]
        tmp2 = 2 * connectivity[n2, 1]

        # i = linkID. 1, 2 are the first and second nodes in link with linkID
        connectivity[n1, tmp1:(tmp1 + 1)] = [i, 1]
        connectivity[n2, tmp2:(tmp2 + 1)] = [i, 2]

        linksConnect[i, 1] = connectivity[n1, 1]
        linksConnect[i, 2] = connectivity[n2, 1]
    end

    network.connectivity = connectivity
    network.linksConnect = linksConnect

    return network
end

"""
```
checkNetwork(network::DislocationNetwork)
```
Checks the validity of the dislocation network. It ensures the following conditions are met by the member variables of `network`:
1. `connectivity` and `links` have the same number of non-zero entries;
1. all entries in `bVec` are non-zero;
1. only the trailing columns of `connectivity` are zeros;
1. consistency between `connectivity` and `links`;
1. `bVec` is conserved among connected nodes;
1. entries in `links` are unique;
1. consistency betwen `connectivity` and `linksConnect`

"""
function checkNetwork(network::DislocationNetwork)
    label = network.label
    idx = findall(x -> x != -1, label)
    links = network.links
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # Connectivity and links must be the same length.
    maximum(connectivity) == length(links[idx, 1]) ? nothing :
    error("Non-empty entries of connectivity should be the same as the non-empty entries of links.")

    bVec = network.bVec
    bSum = zeros(Float64, 3)
    @inbounds for i in idx
        iLinkBuffer = zeros(Int64, 0)
        col = connectivity[i, 1]

        # Check for zero Burgers vectors.
        norm(bVec[i, :]) < eps(eltype(bVec)) ?
        error("Burgers vector must be non-zero.
norm(bVec) = $(norm(bVec[i,:]))") : nothing

        # Trailing columns must be zero.
        sum(connectivity[i, (2 * (col + 1)):end]) != 0 ?
        error("Trailing columns of connectivity must be zero.
        connectivity[$(i), $(2*(col + 1)):end] = $(connectivity[i, (2 * (col + 1)):end]))") :
        nothing

        # Burgers vector conservation.
        bSum .= 0
        # No repeated links.
        neighbours = -ones(Int64, col)
        for j = 1:col
            j2 = 2 * j

            iLink = connectivity[i, j2]     # Link ID.
            pLink = connectivity[i, j2 + 1] # Link position in links.
            bSum += (3 - 2 * pLink) .* bVec[iLink, :] # Running total of burgers vectors.

            neighbours[j] = links[iLink, 3 - connectivity[i, 2 * j + 1]] # Find neighbouring nodes.
            push!(iLinkBuffer, iLink) # Push to link buffer for error printing.

            # Check that node does not appear twice in connectivity.
            for k = (j + 1):col
                connectivity[i, j2] == connectivity[i, 2 * k] ?
                error("Node $(i) cannot appear twice in connectivity, $(connectivity[i, j2])") :
                nothing
            end

            # Check links and connectivity are consistent.
            testNode = links[connectivity[i, j2], connectivity[i, j2 + 1]]
            testNode != i ?
            error("Connectivity and links are not consistent.
connectivity[$(i), :] = $(connectivity[i,:])
links[$(connectivity[i, j2]), :] = [$(links[connectivity[i, j2],1]), $(links[connectivity[i, j2],2])]") :
            nothing
        end

        # Check for Burgers vector conservation.
        abs(dot(bSum, bSum)) > eps(eltype(bSum)) ?
        error("Burgers vector is not conserved, bSum = $(bSum).") : nothing

        # Find unique neighbours.
        uniqueNeighbours = unique(neighbours)
        length(uniqueNeighbours) != length(neighbours) ?
        error("There must be no duplicate links. Each entry in links must correspond to a unique pair of adjacent entries in neighbours.
links = $(iLinkBuffer[1:end-1])
uniqueNeighbours = $(uniqueNeighbours)
neighbours = $(neighbours)") :
        nothing

        # Check connectivity and linksConnect are consistent.
        check = [
            connectivity[links[i, 1], 2 * linksConnect[i, 1]]
            connectivity[links[i, 2], 2 * linksConnect[i, 2]]
        ]
        i != check[1] | i != check[2] ? error("Inconsistent link.
        links[$(i), :] = $(links[i,:])
        linksConnect[$(i),:] = $(linksConnect[i,:])
        connectivity[$(links[i,1]), :] = $(connectivity[links[i,1], :])
        connectivity[$(links[i,2]), :] = $(connectivity[links[i,2], :])") :
        nothing
    end

    return true
end

function getSegmentIdx!(network::DislocationNetwork)
    links = network.links
    label = network.label
    segIdx = zeros(Int64, size(links, 1), 3)
    idx = findall(x -> x != -1, label)
    numSeg::Integer = 0
    for i in idx
        n1 = links[i, 1]
        n2 = links[i, 2]
        (label[n1] == 4 || label[n2] == 4) ? continue : nothing
        numSeg += 1
        segIdx[numSeg, :] = [i, n1, n2]
    end

    network.numSeg = numSeg
    network.segIdx = segIdx

    return network
end
