"""
```
DislocationNetwork(;
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNodeSegConnect::T4 = [0, 0, 0],
    connectivity::T5 = zeros(Int, 0, 0),
    linksConnect::T5 = zeros(Int, 2, 0),
    segIdx::T5 = zeros(Int, 2, 3),
    segForce::T6 = zeros(3, 2, 0),
) where {
    T1 <: AbstractArray{T, N} where {T, N},
    T2 <: AbstractArray{T, N} where {T, N},
    T3 <: AbstractVector{nodeType},
    T4 <: Int,
    T5 <: AbstractArray{Int, N} where {N},
    T6 <: AbstractArray{T, N} where {T, N},
}
```
Keyword constructor for [`DislocationNetwork`](@ref), performs validations but creates dislocation network as provided.
"""
function DislocationNetwork(;
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNode::T4 = zeros(Int, 1),
    numSeg::T4 = zeros(Int, 1),
    maxConnect::T5 = 0,
    connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::T1 = zeros(Int, 2, size(links, 2)),
    segIdx::T1 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links, 2), 0),
) where {T1 <: AbstractArray{T,N} where {T,N},T2 <: AbstractArray{T,N} where {T,N},T3 <: AbstractVector{nodeType},T4 <: AbstractVector{Int},T5 <: Int,T6 <: AbstractArray{T,N} where {T,N},}

    return DislocationNetwork(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        nodeVel,
        nodeForce,
        numNode,
        numSeg,
        maxConnect,
        connectivity,
        linksConnect,
        segIdx,
        segForce,
    )
end

"""
```
DislocationNetwork(
    sources::T1,
    maxConnect::T2 = 4,
    args...;                        # Optional arguments
    memBuffer = nothing,            # Buffer for memory allocation
    checkConsistency::T3 = true,    # Check consistency of generated network
    kw...,                          # Other keyword arguments
) where {
    T1 <: Union{T, AbstractVector{T}} where {T <: DislocationLoop},
    T2 <: Int,
    T3 <: Bool,
}
```
Out of place constructor for [`DislocationNetwork`](@ref). Generates a new dislocation network from previously generated sources.

## Argument Explanation

- `args...` are optional arguments that will be passed on to the `loopDistribution` function which distributes the loops in `sources` according to the type of their `dist` variable.
- `kw...` are optional keyword arguments that will also be passed to `loopDistribution`.
- `memBuffer` is the numerical value for allocating memory in advance, the quantity ``\\textrm{memBuffer} \\times N`` where `N` is the total number of nodes in `sources`, will be the initial number of entries allocated in the matrices that keep the network's data, if it is `nothing` then the number of entries is ``\\textrm{round}(N \\log_{2}(N))``.
"""
function DislocationNetwork(
    sources::T1,
    maxConnect::T2 = 4,
    args...;
    memBuffer = nothing,
    checkConsistency::T3 = true,
    kw...,
) where {T1 <: Union{T,AbstractVector{T},NTuple{N,T} where N} where {T <: DislocationLoop},T2 <: Int,T3 <: Bool,}

    # Initialisation.
    nodeTotal::Int = 0
    lims = zeros(MMatrix{3,2})
    # Calculate node total.
    for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end
    # Memory buffer.
    isnothing(memBuffer) ? nodeBuffer = Int(round(nodeTotal * log2(nodeTotal))) :
    nodeBuffer = nodeTotal * Int(memBuffer)

    # Allocate memory.
    links = zeros(Int, 2, nodeBuffer)
    slipPlane = zeros(3, nodeBuffer)
    bVec = zeros(3, nodeBuffer)
    coord = zeros(3, nodeBuffer)
    label = zeros(nodeType, nodeBuffer)
    nodeVel = zeros(Float64, 3, nodeBuffer)
    nodeForce = zeros(Float64, 3, nodeBuffer)
    numNode = nodeTotal
    numSeg = nodeTotal
    segForce = zeros(Float64, 3, 2, nodeBuffer)

    initIdx = 1
    makeNetwork!(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        sources,
        lims,
        initIdx,
        args...;
        kw...,
    )

    # Calculate number of segments and indexing matrix.
    numSeg, segIdx = getSegmentIdx(links, label)
    # Generate connectivity and linksConnect matrix.
    connectivity, linksConnect = makeConnect(links, maxConnect)

    # Create network.
    network = DislocationNetwork(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        nodeVel,
        nodeForce,
        [numNode],
        [numSeg],
        maxConnect,
        connectivity,
        linksConnect,
        segIdx,
        segForce,
    )

    # Check that the network is generated properly.
    checkConsistency ? checkNetwork(network) : nothing

    return network
end

"""
```
function DislocationNetwork!(
    network::T1,
    sources::T2,
    maxConnect::T3 = 4,
    args...;
    checkConsistency::T3 = true,
    kw...,
) where {
    T1 <: DislocationNetwork,
    T2 <: Union{T, AbstractVector{T}} where {T <: DislocationLoop},
    T3 <: Int,
    T4 <: Bool,
}
```
In-place constructor for [`DislocationNetwork`](@ref). Generates a new dislocation network from already generated sources. If the matrices already in `network` are not large enough to accommodate the additions from `sources`, it will automatically allocate ``\\textrm{round}(N \\log_{2}(N))`` new entries where `N` is the total number of nodes in `sources`.
"""
function DislocationNetwork!(
    network::T1,
    sources::T2,
    maxConnect::T3 = 4,
    args...;
    checkConsistency::T4 = true,
    kw...,
) where {T1 <: DislocationNetwork,T2 <: Union{T,AbstractVector{T},NTuple{N,T} where N} where {T <: DislocationLoop},T3 <: Int,T4 <: Bool,}
    # For comments see DislocationNetwork. It is a 1-to-1 translation except that this one modifies the network in-place.
    
    iszero(network) && return DislocationNetwork(
        sources,
        maxConnect = maxConnect,
        args...;
        checkConsistency = checkConsistency,
        kw...,
    )
    
    @assert network.maxConnect == maxConnect "Maximum connectivity of added network must be equal to that of the existing network."
    
    nodeTotal::Int = 0
    lims = zeros(MMatrix{3,2})
    for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end
    numNode = nodeTotal

    # Allocate memory.
    available = length(findall(x -> x == 0, network.label))
    if nodeTotal > available
        newEntries = Int(round(nodeTotal * log2(nodeTotal)))
        network = push!(network, newEntries)
    end
    
    links = network.links
    slipPlane = network.slipPlane
    bVec = network.bVec
    coord = network.coord
    label = network.label
    
    initIdx::Int = 1
    first = findfirst(x -> x == 0, label)
    isnothing(first) ? initIdx = 1 : initIdx = first
    makeNetwork!(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        sources,
        lims,
        initIdx,
        args...;
        kw...,
    )
    network.numNode[1] += numNode
    
    getSegmentIdx!(network)
    makeConnect!(network)
    
    checkConsistency ? checkNetwork(network) : nothing
    return network
end

function makeNetwork!(
    links,
    slipPlane,
    bVec,
    coord,
    label,
    sources,
    lims,
    initIdx,
    args...;
    kw...,
)
    nodeTotal::Int = 0
    elemT = eltype(coord)
    @inbounds for i in eachindex(sources)
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        # Calculate the normalised displacements for all loops in sources[i] according to their distribution.
        disp = loopDistribution(sources[i].dist, numLoops, args...; kw...)
        # Calculate the real spatial limits of the distributions.
        limits!(lims, mean(sources[i].segLen), sources[i].range, sources[i].buffer)
        # Fill out the data for all loops specified in sources[i].
        for j in 1:numLoops
            # The number of nodes in the loop is nodesLoop, so that's our stride inside sources[i]
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            # Links are numbered sequentially in network so we have to account for previously assigned links.
            links[:, idxi:idxf] .=
                sources[i].links[:, 1:nodesLoop] .+ (nodeTotal + initIdx - 1)
            slipPlane[:, idxi:idxf] = sources[i].slipPlane[:, 1:nodesLoop]
            bVec[:, idxi:idxf] = sources[i].bVec[:, 1:nodesLoop]
            coord[:, idxi:idxf] = sources[i].coord[:, 1:nodesLoop]
            label[idxi:idxf] = sources[i].label[1:nodesLoop]
            # Map the normalised displacements to real space using the real limits and translate the nodes' coordinates accordingly.
            staticDisp = SVector{3,elemT}(disp[1, j], disp[2, j], disp[3, j])
            viewCoord = @view coord[:, idxi:idxf]
            translatePoints!(viewCoord, lims, staticDisp)
            nodeTotal += nodesLoop
        end
    end
    return nothing
end
