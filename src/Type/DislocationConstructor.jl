"""
```
SlipSystem(crystalStruct::T1, slipPlane::T2, bVec::T2) where {T1 <: AbstractCrystalStruct, T2 <: AbstractArray{T, N} where {T, N}}
```
Keyword constructor for [`SlipSystem`](@ref). Throws error if ``\\bm{b} \\not\\perp \\bm{n}`` where ``\\bm{b}`` is the Burgers vector and ``\\bm{n}`` the slip plane.
"""
function SlipSystem(;
    crystalStruct::T1,
    slipPlane::T2,
    bVec::T2
) where {
    T1 <: AbstractCrystalStruct,
    T2 <: AbstractArray{T, N} where {T, N}
}
    return SlipSystem(crystalStruct, slipPlane, bVec)
end

"""
```
function DislocationParameters(coreRad::T1, coreRadMag::T1, minSegLen::T1, maxSegLen::T1, minArea::T1, maxArea::T1, maxConnect::T2, remesh::T3, collision::T3, separation::T3, virtualRemesh::T3, edgeDrag::T1, screwDrag::T1, climbDrag::T1, lineDrag::T1, mobility::T4, parCPU::T3 = false, parGPU::T3 = false) where {T1, T2 <: Int, T3 <: Bool, T4 <: AbstractMobility}
```
Keyword constructor for [`DislocationParameters`](@ref). Validates values and calculates derived quantities.
"""
function DislocationParameters(
    coreRad::T1, 
    coreRadMag::T1,
    minSegLen::T1,
    maxSegLen::T1,
    minArea::T1,
    maxArea::T1,
    edgeDrag::T1,
    screwDrag::T1,
    climbDrag::T1,
    lineDrag::T1,
    maxConnect::T2,
    mobility::T3,
    remesh::T4 = true,
    collision::T4 = true,
    separation::T4 = true,
    virtualRemesh::T4 = true,
    parCPU::T4 = false,
    parGPU::T4 = false
) where {
    T1,
    T2 <: Int,
    T3 <: AbstractMobility,
    T4 <: Bool
}

    coreRad == minSegLen == maxSegLen == 0 ? nothing :
    @assert coreRad < minSegLen < maxSegLen
    minArea == maxArea == 0 ? nothing : @assert minArea < maxArea

    return DislocationParameters(
        coreRad,
        coreRad^2,
        coreRadMag,
        minSegLen,
        maxSegLen,
        minSegLen*2,
        minArea,
        maxArea,
        minArea^2,
        maxArea^2,
        edgeDrag,
        screwDrag,
        climbDrag,
        lineDrag,
        maxConnect,
        mobility,
        remesh,
        collision,
        separation,
        virtualRemesh,
        parCPU,
        parGPU
    )
end
function DislocationParameters(;
    coreRad::T1, 
    coreRadMag::T1,
    minSegLen::T1,
    maxSegLen::T1,
    minArea::T1,
    maxArea::T1,
    edgeDrag::T1,
    screwDrag::T1,
    climbDrag::T1,
    lineDrag::T1,
    maxConnect::T2,
    mobility::T3,
    remesh::T4 = true,
    collision::T4 = true,
    separation::T4 = true,
    virtualRemesh::T4 = true,
    parCPU::T4 = false,
    parGPU::T4 = false
) where {
    T1,
    T2 <: Int,
    T3 <: AbstractMobility,
    T4 <: Bool
}
    return DislocationParameters(
        coreRad, 
        coreRadMag,
        minSegLen,
        maxSegLen,
        minArea,
        maxArea,
        edgeDrag,
        screwDrag,
        climbDrag,
        lineDrag,
        maxConnect,
        mobility,
        remesh,
        collision,
        separation,
        virtualRemesh,
        parCPU,
        parGPU
    )
end

"""
```
DislocationLoop(
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
    T1 <: loopDln,
    T2 <: Int,
    T3 <: Union{T where {T}, AbstractArray{T, N} where {T, N}},
    T4 <: AbstractArray{T, N} where {T, N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution,
}
```
Constructor for a "zero" [`DislocationLoop`](@ref).
"""
function DislocationLoop(
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
    dist::T8
) where {
    T1 <: loopDln,
    T2 <: Int,
    T3 <: Union{T where {T}, AbstractArray{T, N} where {T, N}},
    T4 <: AbstractArray{T, N} where {T, N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution
}

    nodeTotal::Int = 0
    links = zeros(Int, 2, nodeTotal)
    coord = zeros(3, nodeTotal)
    slipPlane = zeros(3, 0)
    bVec = zeros(3, 0)

    return DislocationLoop(
        loopType,
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
        dist
    )
end

"""
```
DislocationLoop(
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
    T2 <: Int,
    T3 <: Union{T where {T}, AbstractArray{T, N} where {T, N}},
    T4 <: AbstractArray{T, N} where {T, N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution,
}
```
Validates inputs and generates a [`DislocationLoop`](@ref) of `loopType` defined by the arguments.
"""
function DislocationLoop(
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
    dist::T8
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int,
    T3 <: Union{T where {T}, AbstractArray{T, N} where {T, N}},
    T4 <: AbstractArray{T, N} where {T, N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution
}

    nodeTotal = numSides * nodeSide # Calculate total number of nodes for memory allocation.
    numSegLen = length(segLen) # Number of segment lengths.

    # Validate input.
    @assert length(label) == nodeTotal "DislocationLoop: All $nodeTotal nodes must be labelled. There are $(length(label)) labels currently defined."
    @assert numSegLen == nodeTotal "DislocationLoop: All $nodeTotal segments must have their lengths defined. There are $numSegLen lengths currently defined."

    # Normalise vectors.
    _slipPlane = _slipPlane ./ norm(_slipPlane)
    _bVec = _bVec ./ norm(_bVec)

    # Pick rotation axis for segments.
    rotAxis = zeros(eltype(_slipPlane), 3)
    # Shear loops rotate around slip plane vector. They have screw, mixed and edge segments.
    if typeof(loopType) == loopShear
        rotAxis = _slipPlane
        # Prismatic loops rotate around Burgers vector. All segments are edge.
    elseif typeof(loopType) == loopPrism
        rotAxis = _bVec
        # Catch all.
    else
        @warn "DislocationLoop: rotation axis for $(typeof(loopType)) not defined, defaulting to prismatic loop."
        rotAxis = _bVec
    end

    # Allocate arrays.
    links = zeros(Int, 2, nodeTotal)
    coord = zeros(3, nodeTotal)
    slipPlane = repeat(_slipPlane, inner = (1, numSegLen))
    bVec = repeat(_bVec, inner = (1, numSegLen))
    seg = zeros(3, numSegLen)

    # Create initial segments.
    for i in eachindex(segLen)
        seg[:, i] = makeSegment(segEdge(), _slipPlane, _bVec) .* segLen[i]
    end

    θ = externalAngle(numSides)  # External angle of a regular polygon with numSides.

    # Loop over polygon's sides.
    for i in 1:numSides
        # Index for side i.
        idx = (i - 1) * nodeSide
        # Rotate segments by external angle of polygon to make polygonal loop.
        rseg = rot3D(seg[:, mod(i - 1, numSegLen) + 1], rotAxis, zeros(3), θ * (i - 1))
        # DO NOT add @simd, this loop works by adding rseg to the previous coordinate to make the loop. Loop over the nodes per side.
        for j in 1:nodeSide
            # Count first node once.
            if i == j == 1
                coord[:, 1] = zeros(3)  # Initial coordinate is on the origin.
                continue
            end
            if idx + j <= nodeTotal
                # Add segment vector to previous coordinate.
                coord[:, idx + j] += coord[:, idx + j - 1] + rseg
            end
        end
    end

    # Find centre of the loop and make it zero.
    meanCoord = mean(coord, dims = 2)
    coord .-= meanCoord

    # Create links matrix.
    for j in 1:(nodeTotal - 1)
        links[:, j] = [j; j + 1]
    end
    links[:, nodeTotal] = [nodeTotal; 1]

    return DislocationLoop(
        loopType,
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
DislocationLoop(;
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
    T2 <: Int,
    T3 <: Union{T where {T}, AbstractArray{T, N} where {T, N}},
    T4 <: AbstractArray{T, N} where {T, N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution,
}
```
Generic keyword constructor for [`DislocationLoop`](@ref). Calls other constructors that dispatch on [`loopType`](@ref).
"""
function DislocationLoop(;
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
    dist::T8
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int,
    T3 <: Union{T where {T}, AbstractArray{T, N} where {T, N}},
    T4 <: AbstractArray{T, N} where {T, N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution
}
    return DislocationLoop(
        loopType,
        numSides,
        nodeSide,
        numLoops,
        segLen,
        slipSystem,
        _slipPlane,
        _bVec,
        label,
        buffer,
        range,
        dist
    )
end

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
    maxConnect::T4 = zeros(Int, 1),
    connectivity::T5 = zeros(Int, 1 + 2 * maxConnect[1], numNode[1]),
    linksConnect::T5 = zeros(Int, 2, numSeg[1]),
    segIdx::T5 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links, 2), 0),
) where {
    T1 <: AbstractArray{T, N} where {T, N},
    T2 <: AbstractArray{T, N} where {T, N},
    T3 <: AbstractVector{nodeType},
    T4 <: AbstractVector{Int},
    T5 <: AbstractArray{Int, N} where {N},
    T6 <: AbstractArray{T, N} where {T, N},
}

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
) where {
    T1 <: Union{T, AbstractVector{T}} where {T <: DislocationLoop},
    T2 <: Int,
    T3 <: Bool,
}

    # Initialisation.
    nodeTotal::Int = 0
    lims = zeros(3, 2)
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

    nodeTotal = 0
    initIdx = 1
    makeNetwork!(links, slipPlane, bVec, coord, label, nodeTotal, sources, lims, initIdx, args...; kw...)

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
        [maxConnect],
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
DislocationNetwork(
    network::T1,
    sources::T2,
    maxConnect::T3 = 4,
    args...;
    checkConsistency::T4 = true,
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
) where {
    T1 <: DislocationNetwork,
    T2 <: Union{T, AbstractVector{T}} where {T <: DislocationLoop},
    T3 <: Int,
    T4 <: Bool,
}
    # For comments see DislocationNetwork. It is a 1-to-1 translation except that this one modifies the network in-place.

    nodeTotal::Int = 0
    lims = zeros(3, 2)
    network.maxConnect[1] = maxConnect
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
    nodeTotal = 0
    first = findfirst(x -> x == 0, label)
    isnothing(first) ? initIdx = 1 : initIdx = first
    makeNetwork!(links, slipPlane, bVec, coord, label, nodeTotal, sources, lims, initIdx, args...; kw...)
    network.numNode[1] += numNode

    getSegmentIdx!(network)
    makeConnect!(network)

    checkConsistency ? checkNetwork(network) : nothing
    return network
end

function makeNetwork!(links, slipPlane, bVec, coord, label, nodeTotal, sources, lims, initIdx, args...; kw...)
    for i in eachindex(sources)
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        disp = loopDistribution(sources[i].dist, numLoops, args...; kw...)
        limits!(lims, mean(sources[i].segLen), sources[i].range, sources[i].buffer)
        for j in 1:numLoops
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            links[:, idxi:idxf] =
                sources[i].links[:, 1:nodesLoop] .+ (nodeTotal + initIdx - 1)
            slipPlane[:, idxi:idxf] .= sources[i].slipPlane[:, 1:nodesLoop]
            bVec[:, idxi:idxf] .= sources[i].bVec[:, 1:nodesLoop]
            coord[:, idxi:idxf] .= sources[i].coord[:, 1:nodesLoop]
            label[idxi:idxf] .= sources[i].label[1:nodesLoop]
            coord[:, idxi:idxf] .=
                translatePoints(sources[i].coord[:, 1:nodesLoop], lims, disp[:, j])
            nodeTotal += nodesLoop
        end
    end
    return nothing
end
