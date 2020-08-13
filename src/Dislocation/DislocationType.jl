## Primitives
"""
```
@enum nodeType begin
    none = 0    # Undefined node, value at initialisation.
    intMob = 1  # Internal mobile node.
    intFix = 2  # Internal fixed node.
    srfMob = 3  # Mobile surface node.
    srfFix = 4  # Fixed surface node.
    ext = 5     # External node.
end
```
"""
@enum nodeType begin
    none = 0
    intMob = 1
    intFix = 2
    srfMob = 3
    srfFix = 4
    ext = 5
end

"""
```
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end    # Undefined segment
struct segEdge <: AbstractDlnSeg end    # Edge segment
struct segEdgeN <: AbstractDlnSeg end   # Edge segment
struct segScrew <: AbstractDlnSeg end   # Screw segment
struct segMixed <: AbstractDlnSeg end   # Mixed segment
```
where `segEdge` have ``(\\bm{b} \\perp \\bm{t}) \\perp \\bm{n}``, `segEdgeN` have ``(\\bm{b} \\perp \\bm{t}) \\parallel \\bm{n}``, `segScrew` have ``\\bm{b} \\parallel \\bm{t}``, `segMixed` have ``\\bm{b} \\not\\perp \\bm{t}~ \\&~ \\bm{b} \\not\\perp \\bm{n}`` and ``\\bm{b}`` is the Burgers vector and ``\\bm{n}`` the slip plane.
"""
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end
struct segEdge <: AbstractDlnSeg end
struct segEdgeN <: AbstractDlnSeg end
struct segScrew <: AbstractDlnSeg end
struct segMixed <: AbstractDlnSeg end

"""
```
abstract type AbstractDlnStr end
struct loopDln <: AbstractDlnStr end    # Unclassified loop
struct loopPrism <: AbstractDlnStr end  # Prismatic loop
struct loopShear <: AbstractDlnStr end  # Shear loop
struct loopJog <: AbstractDlnStr end    # Jog
struct loopKink <: AbstractDlnStr end   # Kink
```
"""
abstract type AbstractDlnStr end
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
struct loopMixed <: AbstractDlnStr end
struct loopDln <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end

"""
```
abstract type AbstractDistribution end
struct Zeros <: AbstractDistribution end
struct Rand <: AbstractDistribution end
struct Randn <: AbstractDistribution end
struct Regular <: AbstractDistribution end
```
Distributions for dislocation sources.
"""
abstract type AbstractDistribution end
struct Zeros <: AbstractDistribution end
struct Rand <: AbstractDistribution end
struct Randn <: AbstractDistribution end
struct Regular <: AbstractDistribution end

"""
```
abstract type AbstractMobility end
struct mobBCC <: AbstractMobility end
struct mobFCC <: AbstractMobility end
struct mobHCP <: AbstractMobility end
```
Mobility functions.
"""
abstract type AbstractMobility end
struct mobBCC <: AbstractMobility end
struct mobFCC <: AbstractMobility end
struct mobHCP <: AbstractMobility end

## Structures
"""
```
struct SlipSystem{T1, T2}
    crystalStruct::T1   # Crystal structure
    slipPlane::T2       # Slip plane
    bVec::T2            # Burgers vector
end
```
Structure to store slip systems.
"""
struct SlipSystem{T1, T2}
    crystalStruct::T1
    slipPlane::T2
    bVec::T2
end
"""
```
SlipSystem(;
    crystalStruct::T1,
    slipPlane::T2,
    bVec::T2,
) where {T1 <: AbstractCrystalStruct, T2 <: AbstractArray{T, N} where {T, N}}
```
Keyword constructor for [`SlipSystem`](@ref). Throws error if ``\\bm{b} \\not\\perp \\bm{n}`` where ``\\bm{b}`` is the Burgers vector and ``\\bm{n}`` the slip plane.
"""
@inline function SlipSystem(;
    crystalStruct::T1,
    slipPlane::T2,
    bVec::T2,
) where {T1 <: AbstractCrystalStruct, T2 <: AbstractArray{T, N} where {T, N}}

    if sum(slipPlane .!= 0) != 0 && sum(bVec .!= 0) != 0
        if ndims(slipPlane) == 1
            @assert isapprox(dot(slipPlane, bVec), 0) "SlipSystem: slip plane, n == $(slipPlane), and Burgers vector, b = $(bVec), must be orthogonal."
        else
            idx = findall(x -> !isapprox(x, 0), dimDot(slipPlane, bVec; dim = 1))
            @assert isempty(idx) "SlipSystem: entries of the slip plane, n[$idx, :] = $(slipPlane[idx,:]), and Burgers vector, b[$idx, :] = $(bVec[idx,:]), are not orthogonal."
        end
    end

    return SlipSystem(crystalStruct, slipPlane, bVec)
end

"""
```
struct DislocationP{T1, T2, T3, T4}
    coreRad::T1         # Core radius
    coreRadSq::T1       # Square of core radius
    coreRadMag::T1      # Magnitude of core radius
    minSegLen::T1       # Minimum segment length
    maxSegLen::T1       # Maximum segment length
    twoMinSegLen::T1    # Twice minimum segment length
    minArea::T1         # Minimum area enclosed by 3 segments
    maxArea::T1         # Maximum area enclosed by 3 segments
    minAreaSq::T1       # Squared min area
    maxAreaSq::T1       # Squared max area
    maxConnect::T2      # Maximum connectivity
    remesh::T3          # Remesh flag
    collision::T3       # Collision flag
    separation::T3      # Separation flag
    virtualRemesh::T3   # Virtual remeshing flag
    parCPU::T3          # Parallelise on CPU
    parGPU::T3          # Parallelise on GPU
    edgeDrag::T1        # Drag coefficient edge dislocation
    screwDrag::T1       # Drag coefficient screw dislocation
    climbDrag::T1       # Drag coefficient climb direction
    lineDrag::T1        # Drag coefficient line direction
    mobility::T4        # Mobility law
end
```
Structure to store dislocation parameters.
"""
struct DislocationP{T1, T2, T3, T4}
    coreRad::T1
    coreRadSq::T1
    coreRadMag::T1
    minSegLen::T1
    maxSegLen::T1
    twoMinSegLen::T1
    minArea::T1
    maxArea::T1
    minAreaSq::T1
    maxAreaSq::T1
    maxConnect::T2
    remesh::T3
    collision::T3
    separation::T3
    virtualRemesh::T3
    parCPU::T3
    parGPU::T3
    edgeDrag::T1
    screwDrag::T1
    climbDrag::T1
    lineDrag::T1
    mobility::T4
end
"""
```
DislocationP(;
    coreRad::T1,
    coreRadMag::T1,
    minSegLen::T1,
    maxSegLen::T1,
    minArea::T1,
    maxArea::T1,
    maxConnect::T2,
    remesh::T3,
    collision::T3,
    separation::T3,
    virtualRemesh::T3,
    edgeDrag::T1,
    screwDrag::T1,
    climbDrag::T1,
    lineDrag::T1,
    mobility::T4,
) where {T1, T2 <: Int, T3 <: Bool, T4 <: AbstractMobility}
```
Keyword constructor for [`DislocationP`](@ref). Validates values and calculates derived quantities.
"""
@inline function DislocationP(;
    coreRad::T1,
    coreRadMag::T1,
    minSegLen::T1,
    maxSegLen::T1,
    minArea::T1,
    maxArea::T1,
    maxConnect::T2,
    remesh::T3,
    collision::T3,
    separation::T3,
    virtualRemesh::T3,
    parCPU::T3 = false,
    parGPU::T3 = false,
    edgeDrag::T1,
    screwDrag::T1,
    climbDrag::T1,
    lineDrag::T1,
    mobility::T4,
) where {T1, T2 <: Int, T3 <: Bool, T4 <: AbstractMobility}

    coreRad == minSegLen == maxSegLen == 0 ? nothing :
    @assert coreRad < minSegLen < maxSegLen
    minArea == maxArea == 0 ? nothing : @assert minArea < maxArea
    coreRadSq = coreRad^2

    return DislocationP(
        coreRad,
        coreRadSq,
        coreRadMag,
        minSegLen,
        maxSegLen,
        minSegLen * 2,
        minArea,
        maxArea,
        minArea^2,
        maxArea^2,
        maxConnect,
        remesh,
        collision,
        separation,
        virtualRemesh,
        parCPU,
        parGPU,
        edgeDrag,
        screwDrag,
        climbDrag,
        lineDrag,
        mobility,
    )
end

"""
```
struct DislocationLoop{T1, T2, T3, T4, T5, T6, T7, T8, T9}
    loopType::T1    # Loop type
    numSides::T2    # Number of sides per loop
    nodeSide::T2    # Number of nodes per side
    numLoops::T2    # Number of loops to generate
    segLen::T3      # Segment lengths
    slipSystem::T4  # Slip system
    links::T5       # Links matrix
    slipPlane::T6   # Slip plane for each link
    bVec::T6        # Burgers vector for each link
    coord::T6       # Coordinates of each node
    label::T7       # Label of each node
    buffer::T8      # Mean distance buffer separating each loop centre
    range::T6       # Distribution range of generated loops
    dist::T9        # Distribution of generated loops
end
```
Structure to store dislocation loop.
"""
struct DislocationLoop{T1, T2, T3, T4, T5, T6, T7, T8, T9}
    loopType::T1
    numSides::T2
    nodeSide::T2
    numLoops::T2
    segLen::T3
    slipSystem::T4
    links::T5
    slipPlane::T6
    bVec::T6
    coord::T6
    label::T7
    buffer::T8
    range::T6
    dist::T9
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
@inline function DislocationLoop(;
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

    return DislocationLoop(
        loopType;
        numSides = numSides,
        nodeSide = nodeSide,
        numLoops = numLoops,
        segLen = segLen,
        slipSystem = slipSystem,
        _slipPlane = _slipPlane,
        _bVec = _bVec,
        label = label,
        buffer = buffer,
        range = range,
        dist = dist,
    )

end
"""
```
DislocationLoop(
    loopType::T1;
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
@inline function DislocationLoop(
    loopType::T1;
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
        dist,
    )
end
"""
```
DislocationLoop(
    loopType::T1;
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
@inline function DislocationLoop(
    loopType::T1;
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
    @inbounds @simd for i in eachindex(segLen)
        seg[:, i] = makeSegment(segEdge(), _slipPlane, _bVec) .* segLen[i]
    end

    θ = extAngle(numSides)  # External angle of a regular polygon with numSides.

    # Loop over polygon's sides.
    @inbounds for i in 1:numSides
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
    @inbounds @simd for j in 1:(nodeTotal - 1)
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

# TODO: Make DislocationNetwork immutable to take advantage of Julia 1.5's immutable struct optimisations.
"""
```
mutable struct DislocationNetwork{T1, T2, T3, T4, T5, T6}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeVel::T2
    nodeForce::T2
    numNode::T4         # Number of nodes in network
    numSeg::T4          # Number of segments in network
    maxConnect::T4      # Maximum connections per node
    connectivity::T5    # Connectivity matrix
    linksConnect::T5    # Links involved in connection
    segIdx::T5          # Contains segment index and the nodes of the nodes in said link
    segForce::T6        # Force on each node of each segment
end
```
Structure to store dislocation network.
"""
mutable struct DislocationNetwork{T1, T2, T3, T4, T5, T6}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeVel::T2
    nodeForce::T2
    numNode::T4
    numSeg::T4
    maxConnect::T4
    connectivity::T5
    linksConnect::T5
    segIdx::T5
    segForce::T6
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
    numNode::T4 = 0,
    numSeg::T4 = 0,
    maxConnect::T4 = 0,
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
@inline function DislocationNetwork(;
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNode::T4 = 0,
    numSeg::T4 = 0,
    maxConnect::T4 = 0,
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

    @assert size(links, 1) == size(segForce, 2) == 2
    @assert size(bVec, 1) == size(slipPlane, 1) == size(coord, 1) size(segForce, 1) == 3
    @assert size(links, 2) == size(bVec, 2) == size(slipPlane, 2) == size(segForce, 3)
    @assert size(coord, 2) == length(label)

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
@inline function DislocationNetwork(
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
    @inbounds for i in eachindex(sources)
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
    @inbounds for i in eachindex(sources)
        # Indices.
        idx = 1 + nodeTotal
        nodesLoop = length(sources[i].label)    # Number of nodes in a loop from sources[i].
        numLoops = sources[i].numLoops          # Number loops to generate from sources[i].
        numNodes = numLoops * nodesLoop         # Number of nodes in a loop from sources[i].
        # Calculate spatial distribution for sources[i].
        disp = loopDistribution(sources[i].dist, numLoops, args...; kw...)
        # Calculate spatial limits of the loops generated by sources[i].
        limits!(lims, mean(sources[i].segLen), sources[i].range, sources[i].buffer)
        # Loop through the number of loops generated from sources[i].
        for j in 1:numLoops
            # Initial and final indices of the nodes in sources[i].
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            # Add links from sources[i] to links and adjust by the total number of nodes in the network.
            links[:, idxi:idxf] = sources[i].links[:, 1:nodesLoop] .+ nodeTotal
            # Pass slip plane, burgers vector, coordinates and labels from sources[i] to network.
            slipPlane[:, idxi:idxf] .= sources[i].slipPlane[:, 1:nodesLoop]
            bVec[:, idxi:idxf] .= sources[i].bVec[:, 1:nodesLoop]
            coord[:, idxi:idxf] .= sources[i].coord[:, 1:nodesLoop]
            label[idxi:idxf] .= sources[i].label[1:nodesLoop]
            # Translate loops according to the previously calcualted displacements and limits.
            coord[:, idxi:idxf] .=
                translatePoints(sources[i].coord[:, 1:nodesLoop], lims, disp[:, j])
            # Add nodes in the loop to the total number of nodes.
            nodeTotal += nodesLoop
        end
    end

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
        numNode,
        numSeg,
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
    T2 <: Union{T, AbstractVector{T}} where {T <: DislocationLoop},
    T3 <: Int,
    T4 <: Bool,
}
```
In-place constructor for [`DislocationNetwork`](@ref). Generates a new dislocation network from already generated sources. If the matrices already in `network` are not large enough to accommodate the additions from `sources`, it will automatically allocate ``\\textrm{round}(N \\log_{2}(N))`` new entries where `N` is the total number of nodes in `sources`.
"""
@inline function DislocationNetwork!(
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
    network.maxConnect = maxConnect
    @inbounds for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end

    # Allocate memory.
    available = length(findall(x -> x == 0, network.label))
    if nodeTotal > available
        newEntries = Int(round(nodeTotal * log2(nodeTotal)))
        push!(network, newEntries)
    end

    initIdx::Int = 1
    nodeTotal = 0
    first = findfirst(x -> x == 0, network.label)
    isnothing(first) ? initIdx = 1 : initIdx = first
    @inbounds for i in eachindex(sources)
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        disp = loopDistribution(sources[i].dist, numLoops, args...; kw...)
        limits!(lims, mean(sources[i].segLen), sources[i].range, sources[i].buffer)
        for j in 1:numLoops
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            network.links[:, idxi:idxf] =
                sources[i].links[:, 1:nodesLoop] .+ (nodeTotal + initIdx - 1)
            network.slipPlane[:, idxi:idxf] .= sources[i].slipPlane[:, 1:nodesLoop]
            network.bVec[:, idxi:idxf] .= sources[i].bVec[:, 1:nodesLoop]
            network.coord[:, idxi:idxf] .= sources[i].coord[:, 1:nodesLoop]
            network.label[idxi:idxf] .= sources[i].label[1:nodesLoop]
            network.coord[:, idxi:idxf] .=
                translatePoints(sources[i].coord[:, 1:nodesLoop], lims, disp[:, j])
            nodeTotal += nodesLoop
        end
    end
    network.numNode += nodeTotal

    getSegmentIdx!(network)
    makeConnect!(network)

    checkConsistency ? checkNetwork(network) : nothing

    return network
end
