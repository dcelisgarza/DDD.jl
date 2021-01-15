"""
```
SlipSystem(crystalStruct::T1, slipPlane::T2, bVec::T2) 
    where {T1 <: AbstractCrystalStruct,T2}
```
Constructor for [`SlipSystem`](@ref). Checks for orthogonality of the burgers vector and slip plane. It assumes each column corresponds to a slip system.
"""
function SlipSystem(
        crystalStruct::T1,
        slipPlane::T2,
        bVec::T2,
    ) where {T1 <: AbstractCrystalStruct,T2}
    if sum(slipPlane .!= 0) != 0 && sum(bVec .!= 0) != 0
        if ndims(slipPlane) == 1
            @assert isapprox(dot(slipPlane, bVec), 0) "SlipSystem: slip plane, n == $(slipPlane), and Burgers vector, b = $(bVec), must be orthogonal."
        else
            idx = findall(x -> !isapprox(x, 0), vec(sum(slipPlane .* bVec, dims = 1)))
            @assert isempty(idx) "SlipSystem: entries of the slip plane, n[$idx, :] = $(slipPlane[idx,:]), and Burgers vector, b[$idx, :] = $(bVec[idx,:]), are not orthogonal."
        end
    end

    return SlipSystem{T1,T2}(crystalStruct, slipPlane, bVec)
end
"""
```
SlipSystem(; crystalStruct::T1, slipPlane::T2, bVec::T2) 
    where {T1 <: AbstractCrystalStruct,T2}
```
Keyword constructor for [`SlipSystem`](@ref).
"""
function SlipSystem(; crystalStruct::T1, slipPlane::T2, bVec::T2) where {T1 <: AbstractCrystalStruct,T2}
    return SlipSystem(crystalStruct, slipPlane, bVec)
end

"""
```
DislocationParameters(
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
    parGPU::T4 = false,
    slipStepCritLen::T1 = maxSegLen / 2,
    slipStepCritArea::T1 = 0.5 * (slipStepCritLen^2) * sin(2 * π / 360),
) where {T1,T2 <: Int,T3 <: AbstractMobility,T4 <: Bool}
```
The constructor for [`DislocationParameters`](@ref) provides a few default values and calculates derived quantities.
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
    parGPU::T4 = false,
    slipStepCritLen::T1 = maxSegLen / 2,
    slipStepCritArea::T1 = 0.5 * (slipStepCritLen^2) * sin(2 * π / 360),
) where {T1,T2 <: Int,T3 <: AbstractMobility,T4 <: Bool}

    coreRad == minSegLen == maxSegLen == 0 ? nothing :
    @assert coreRad < minSegLen < maxSegLen
    minArea == maxArea == 0 ? nothing : @assert minArea < maxArea

    return DislocationParameters(
        coreRad,
        coreRad^2,
        coreRadMag,
        minSegLen,
        maxSegLen,
        minSegLen * 2,
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
        parGPU,
        slipStepCritLen,
        slipStepCritArea,
    )
end
"""
```
DislocationParameters(;
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
    parGPU::T4 = false,
    slipStepCritLen::T1 = maxSegLen / 2,
    slipStepCritArea::T1 = 0.5 * (slipStepCritLen^2) * sin(2 * π / 360),
) where {T1,T2 <: Int,T3 <: AbstractMobility,T4 <: Bool}
```
Keyword constructor for [`DislocationParameters`](@ref). Calls the positional constructor internally.
"""
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
    parGPU::T4 = false,
    slipStepCritLen::T1 = maxSegLen / 2,
    slipStepCritArea::T1 = 0.5 * (slipStepCritLen^2) * sin(2 * π / 360),
) where {T1,T2 <: Int,T3 <: AbstractMobility,T4 <: Bool}
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
        parGPU,
        slipStepCritLen,
        slipStepCritArea,
)
end

"""
```
DislocationLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen,
    slipSystem::T2,
    _slipPlane,
    _bVec,
    label::T3,
    buffer,
    range,
    dist::T4,
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int,
    T3 <: AbstractVector{nodeType},
    T4 <: AbstractDistribution,
}
```
General constructor for [`DislocationLoop`](@ref).
"""
function DislocationLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen,
    slipSystem::T2,
    _slipPlane,
    _bVec,
    label::T3,
    buffer,
    range,
    dist::T4,
) where {T1 <: AbstractDlnStr,T2 <: Int,T3 <: AbstractVector{nodeType},T4 <: AbstractDistribution,}

    nodeTotal::Int = 0
    links = zeros(MMatrix{2,nodeTotal,Int})
    coord = zeros(MMatrix{3,nodeTotal})
    slipPlane = zeros(MMatrix{3,0})
    bVec = zeros(MMatrix{3,0})

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
    segLen,
    slipSystem::T2,
    _slipPlane::T3,
    _bVec::T3,
    label::T4,
    buffer,
    range,
    dist::T5,
) where {
    T1 <: loopPure,
    T2 <: Int,
    T3 <: AbstractArray{T,N} where {T,N},
    T4 <: AbstractVector{nodeType},
    T5 <: AbstractDistribution,
}
```
Constructor for `loopPure` loops [`DislocationLoop`](@ref).
"""
function DislocationLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen,
    slipSystem::T2,
    _slipPlane::T3,
    _bVec::T3,
    label::T4,
    buffer,
    range,
    dist::T5,
) where {T1 <: loopPure,T2 <: Int,T3 <: AbstractArray{T,N} where {T,N},T4 <: AbstractVector{nodeType},T5 <: AbstractDistribution,}

    nodeTotal = numSides * nodeSide # Calculate total number of nodes for memory allocation.
    numSegLen = length(segLen) # Number of segment lengths.

    # Validate input.
    @assert length(label) == nodeTotal "DislocationLoop: All $nodeTotal nodes must be labelled. There are $(length(label)) labels currently defined."
    @assert numSegLen == nodeTotal "DislocationLoop: All $nodeTotal segments must have their lengths defined. There are $numSegLen lengths currently defined."

    # Normalise vectors.
    elemT = eltype(_slipPlane)
    _slipPlane = _slipPlane / norm(_slipPlane)
    _bVec = _bVec / norm(_bVec)

    # Pick rotation axis for segments.
    # Shear loops rotate around slip plane vector. They have screw, mixed and edge segments.
    if typeof(loopType) == loopShear
        rotAxis = SVector{3,elemT}(_slipPlane[1], _slipPlane[2], _slipPlane[3])
        # Prismatic loops rotate around Burgers vector. All segments are edge.
    else
        rotAxis = SVector{3,elemT}(_bVec[1], _bVec[2], _bVec[3])
        # Catch all.
    end

    # Allocate arrays.
    links = zeros(MMatrix{2,nodeTotal,Int})
    coord = zeros(MMatrix{3,nodeTotal})
    slipPlane = MMatrix{3,nodeTotal}(repeat(_slipPlane, inner = (1, numSegLen)))
    bVec = MMatrix{3,nodeTotal}(repeat(_bVec, inner = (1, numSegLen)))
    seg = zeros(MMatrix{3,numSegLen})

    # Create initial segments.
    staticSlipPlane = SVector{3,elemT}(_slipPlane[1], _slipPlane[2], _slipPlane[3])
    staticBVec = SVector{3,elemT}(_bVec[1], _bVec[2], _bVec[3])
    @inbounds @simd for i in eachindex(segLen)
        seg[:, i] = makeSegment(segEdge(), staticSlipPlane, staticBVec) * segLen[i]
    end

    θ = externalAngle(numSides)  # External angle of a regular polygon with numSides.

    # Loop over polygon's sides.
    origin = SVector{3,elemT}(0, 0, 0)
    @inbounds for i in 1:numSides
        # Index for side i.
        idx = (i - 1) * nodeSide
        # Rotate segments by external angle of polygon to make polygonal loop.
        modIdx = mod(i - 1, numSegLen) + 1
        staticSeg = SVector{3,elemT}(seg[1, modIdx], seg[2, modIdx], seg[3, modIdx])
        rseg = rot3D(staticSeg, rotAxis, origin, θ * (i - 1))
        # DO NOT add @simd, this loop works by adding rseg to the previous coordinate to make the loop. Loop over the nodes per side.
        for j in 1:nodeSide
            # Count first node once.
            if i == j == 1
                coord[:, 1] .= 0 # Initial coordinate is on the origin.
                continue
            end
            if idx + j <= nodeTotal
                # Add segment vector to previous coordinate.
                coord[:, idx + j] += @views coord[:, idx + j - 1] + rseg
            end
        end
    end

    # Find centre of the loop and make it zero.
    meanCoord = mean(coord, dims = 2)
    coord .-= meanCoord

    # Create links matrix.
    @inbounds @simd for j in 1:(nodeTotal - 1)
        links[:, j] .= (j, j + 1)
    end
    links[:, nodeTotal] .= (nodeTotal, 1)

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
DislocationLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen,
    slipSystem::T2,
    _slipPlane::T3,
    _bVec::T3,
    label::T4,
    buffer,
    range,
    dist::T5,
) where {
    T1 <: loopImpure,
    T2 <: Int,
    T3 <: AbstractArray{T,N} where {T,N},
    T4 <: AbstractVector{nodeType},
    T5 <: AbstractDistribution,
}
```
A fallback [`DislocationLoop`](@ref) constructor for other as of yet unimplemented `loopType`.
"""
function DislocationLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen,
    slipSystem::T2,
    _slipPlane::T3,
    _bVec::T3,
    label::T4,
    buffer,
    range,
    dist::T5,
) where {T1 <: loopImpure,T2 <: Int,T3 <: AbstractArray{T,N} where {T,N},T4 <: AbstractVector{nodeType},T5 <: AbstractDistribution,}

    @warn "DislocationLoop: Constructor for $(typeof(loopType)) not defined, defaulting to prismatic loop."
    return DislocationLoop(
        loopPrism(),
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
    T3 <: Union{T where {T},AbstractArray{T,N} where {T,N}},
    T4 <: AbstractArray{T,N} where {T,N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T,N} where {T,N},
    T8 <: AbstractDistribution,
}
```
Keyword [`DislocationLoop`](@ref) constructor calls the positional constructor, which dispatches the appropriate method.
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
    dist::T8,
) where {T1 <: AbstractDlnStr,T2 <: Int,T3 <: Union{T where {T},AbstractArray{T,N} where {T,N}},T4 <: AbstractArray{T,N} where {T,N},T5 <: AbstractVector{nodeType},T6,T7 <: AbstractArray{T,N} where {T,N},T8 <: AbstractDistribution,}
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
DislocationNetwork(
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNode::T4 = zeros(Int, 1),
    numSeg::T4 = zeros(Int, 1),
    maxConnect::T5 = 4,
    connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::T1 = zeros(Int, 2, size(links, 2)),
    segIdx::T1 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links)...),
) where {
    T1 <: AbstractArray{T,N} where {T,N},
    T2 <: AbstractArray{T,N} where {T,N},
    T3 <: AbstractVector{nodeType},
    T4 <: Union{Int,AbstractVector{Int}},
    T5 <: Int,T6 <: AbstractArray{T,N} where {T,N},
}
```
[`DislocationNetwork`](@ref) constructor provides default values, validates inputs and calculates derived quantities.
"""
function DislocationNetwork(
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNode::T4 = zeros(Int, 1),
    numSeg::T4 = zeros(Int, 1),
    maxConnect::T5 = 4,
    connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::T1 = zeros(Int, 2, size(links, 2)),
    segIdx::T1 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links)...),
) where {T1 <: AbstractArray{T,N} where {T,N},T2 <: AbstractArray{T,N} where {T,N},T3 <: AbstractVector{nodeType},T4 <: Union{Int,AbstractVector{Int}},T5 <: Int,T6 <: AbstractArray{T,N} where {T,N},}

    @assert size(links, 1) == size(segForce, 2) == 2
    @assert size(bVec, 1) == size(slipPlane, 1) == size(coord, 1) size(segForce, 1) == 3
    @assert size(links, 2) == size(bVec, 2) == size(slipPlane, 2) == size(segForce, 3)
    @assert size(coord, 2) == length(label)
    @assert length(numNode) == length(numSeg) == 1

    typeof(numNode) <: AbstractVector ? numNodeArr = numNode : numNodeArr = [numNode]
    typeof(numSeg) <: AbstractVector ? numSegArr = numSeg : numSegArr = [numSeg]

    return DislocationNetwork{T1,T2,T3,typeof(numNodeArr),T5,T6}(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        nodeVel,
        nodeForce,
        numNodeArr,
        numSegArr,
        maxConnect,
        connectivity,
        linksConnect,
        segIdx,
        segForce,
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
    numNode::T4 = zeros(Int, 1),
    numSeg::T4 = zeros(Int, 1),
    maxConnect::T5 = 0,
    connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::T1 = zeros(Int, 2, size(links, 2)),
    segIdx::T1 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links, 2), 0),
) where {
    T1 <: AbstractArray{T,N} where {T,N},
    T2 <: AbstractArray{T,N} where {T,N},
    T3 <: AbstractVector{nodeType},T4 <: AbstractVector{Int},
    T5 <: Int,
    T6 <: AbstractArray{T,N} where {T,N},
}
```
Keyword constructor for [`DislocationNetwork`](@ref) calls the positional one.
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
    args...;
    memBuffer = nothing,
    checkConsistency::T3 = true,
    kw...,
) where {
    T1 <: Union{T,
                AbstractVector{T},
                NTuple{N,T} where N
                } where {T <: DislocationLoop},
    T2 <: Int,
    T3 <: Bool,
}
```
The recommended way of creating a network is to use a `source` of type [`DislocationLoop`](@ref). It also accepts arrays and tuples of [`DislocationLoop`](@ref) variables. This automatically generates the network according to the parameters stored in `source` or each of its entries.

# Meaning

* `args...` are optional arguments that will be passed on to the [`loopDistribution`](@ref) function which distributes the loops in `sources` according to the type of their `dist` variable.
* `kw...` are optional keyword arguments that will also be passed to `loopDistribution`.
* `memBuffer` is the numerical value for allocating memory in advance. The quantity ``\\textrm{memBuffer} × N`` where `N` is the total number of nodes in `sources`, will be the initial number of entries allocated in the matrices that keep the network's data. If no `memBuffer` is provided, the number of entries allocated will be ``\\textrm{round}(N \\log_{2}(N))``.
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
DislocationNetwork!(
    network::T1,
    sources::T2,
    maxConnect::T3 = 4,
    args...;
    checkConsistency::T4 = true,
    kw...,
) where {
    T1 <: DislocationNetwork,
    T2 <: Union{T,
                AbstractVector{T},
                NTuple{N,T} where N
                } where {T <: DislocationLoop},
    T3 <: Int,
    T4 <: Bool,
}
```
Adds more `sources` to an existing `network`.
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