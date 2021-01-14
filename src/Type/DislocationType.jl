## Dislocation Type Declarations
"""
```
@enum nodeType begin
    none = 0    # Undefined node, value at initialisation.
    intMob = 1  # Internal mobile node.
    intFix = 2  # Internal fixed node.
    srfMob = 3  # Mobile surface node.
    srfFix = 4  # Fixed surface node.
    ext = 5     # External node.
    tmp = 6     # Temporary flag, used during topological operations.
end
```
Different types of nodes behave differently. There are only a finite number of them so an enumerated type provides safety and efficiency. Each value represents a different type of node and therefore its behaviour.

# Meaning
* `none` are uninitialised nodes.
* `intMob` are mobile nodes internal to the convex hull of the domain. They take part in tractions, displacements and dislocation interactions.
* `intFix` are fixed nodes internal to the convex hull of the domain. They participate in the same way as `intMob` nodes except for the fact that their velocities is fixed are zero.
* `srfMob` are mobile nodes that live on the surface of the convex hull of the domain, they are used to track slip steps and therefore participate in the same things as internal nodes but their velocities are restricted to the convex hull surface.
* `srfFix` are fixed surface nodes and have the same characteristics as mobile surface nodes except for having zero velocity.
* `ext` are external nodes that do not participate in dislocation interactions or forces but are used to calculate displacements and track slip steps.
* `tmp` are nodes that are temporarily flagged before they are assigned another type.
"""
@enum nodeType begin
    none = 0
    intMob = 1
    intFix = 2
    srfMob = 3
    srfFix = 4
    virtual = 5
    tmp = 6
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

These types are used to automatically generate segments out of Burgers vectors ``\\bm{b}``, slip planes ``\\bm{n}``, and/or line direction ``\\bm{l}``.

# Meaning
* `segEdge` have ``\\bm{b} ⟂ \\bm{t}`` ,
* `segEdgeN` have ``\\bm{b} ⟂ \\bm{t}`` and ``\\bm{b} ∥ \\bm{n}`` ,
* `segScrew` have ``\\bm{b} ∥ \\bm{t}`` ,
* `segMixed` have none of the above.
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
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
const loopPure = Union{loopPrism,loopShear}
struct loopMixed <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end
const loopImpure = Union{loopMixed,loopJog,loopKink}
const loopDefined = Union{loopPrism,loopShear,loopMixed,loopJog,loopKink}
struct loopDln <: AbstractDlnStr end
```
These types are used to automatically generate dislocation loops for simulation initialisation.

# Meaning
* `loopPrism` are prismatic loops, their Burgers vectors are perpendicular to the their line direction. They are idealised loops that can be automatically generated as n-gons.
* `loopShear` are shear loops, their line direction goes through edge, screw and line segments as the loop goes round. They are idealised loops that can be automatically generated as n-gons.
* `loopPure` are idealised loops.
* `loopMixed` are loops with prismatic and shear character. They have to be hand-made or require a heuristic to automatically generate.
* `loopDln` is a generic loop used for adding methods to Base functions.
* `loopKink` and `loopJog` are structures formed by colliding dislocations. They are not currently used.
* `loopImpure` are non-idealised loops.
* `loopDefined` are defined loop types.
"""
abstract type AbstractDlnStr end
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
const loopPure = Union{loopPrism,loopShear}
struct loopMixed <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end
const loopImpure = Union{loopMixed,loopJog,loopKink}
const loopDefined = Union{loopPrism,loopShear,loopMixed,loopJog,loopKink}
struct loopDln <: AbstractDlnStr end

"""
```
abstract type AbstractDistribution end
struct Zeros <: AbstractDistribution end
struct Rand <: AbstractDistribution end
struct Randn <: AbstractDistribution end
struct Regular <: AbstractDistribution end
```
Spatial distributions for dislocation sources. These are used to automatically generate networks with a given distribution.

# Meaning
* `Zeros` makes the network generation functions place the center of the generated dislocation loops at the origin. This can be used to generate a network and loops can be manually or pseudo-manually distributed in the domain.
* `Rand` makes the network generation functions uniformly distribute the dislocations according to the range and buffer values in the dislocation loop structure.
* `Rand` makes the network generation functions normally distribute the dislocations according to the range and buffer values in the dislocation loop structure.
* `Rand` TBA, will regularly distribute dislocations according to the range, buffer and other args given to the dislocation network generator.
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
Types to dispatch different mobility functions.

# Meaning
* `mobBCC` is used to dispatch the default BCC mobility function.
* `mobFCC` is used to dispatch the default FCC mobility function.
* `mobHCP` is used to dispatch the default HCP mobility function.
"""
abstract type AbstractMobility end
struct mobBCC <: AbstractMobility end
struct mobFCC <: AbstractMobility end
struct mobHCP <: AbstractMobility end

"""
Slip systems.
```
struct SlipSystem{T1, T2}
    crystalStruct::T1   # Crystal structure
    slipPlane::T2       # Slip plane
    bVec::T2            # Burgers vector
end
```
Structure for storing slip systems. 
    
The constructor,
```
SlipSystem(crystalStruct::T1, slipPlane::T2, bVec::T2) 
    where {T1 <: AbstractCrystalStruct,T2}
```
checks for orthogonality of the burgers vector and slip plane. It assumes each column corresponds to a slip system like so,
```
[
    x1  x2  ... xn;
    y1  y2  ... yn;
    z1  z2  ... zn
]
```
where each number corresponds to a slip plane. The keyword constructor,
```
SlipSystem(; crystalStruct::T1, slipPlane::T2, bVec::T2) 
    where {T1 <: AbstractCrystalStruct,T2}
```
simply calls the positional one.
"""
struct SlipSystem{T1,T2}
    crystalStruct::T1
    slipPlane::T2
    bVec::T2

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

        return new{T1,T2}(crystalStruct, slipPlane, bVec)
    end

end
function SlipSystem(; crystalStruct::T1, slipPlane::T2, bVec::T2) where {T1 <: AbstractCrystalStruct,T2}
    return SlipSystem(crystalStruct, slipPlane, bVec)
end

"""
```
struct DislocationParameters{T1,T2,T3,T4}
    coreRad::T1         # Dislocation core radius.
    coreRadSq::T1       # Square of the dislocation core radius.
    coreRadMag::T1      # Magnitude of the core radius (real units for post-processing).
    minSegLen::T1       # Minimum segment length.
    maxSegLen::T1       # Maximum segment length.
    twoMinSegLen::T1    # Twice the minimum segment length.
    minArea::T1         # Minimum area for remeshing.
    maxArea::T1         # Maximum area for remeshing.
    minAreaSq::T1       # Square of the minimum area.
    maxAreaSq::T1       # Square of the maximum area.
    edgeDrag::T1        # Edge drag coefficient.
    screwDrag::T1       # Screw drag coefficient.
    climbDrag::T1       # Climb drag coefficient.
    lineDrag::T1        # Line drag coefficient.
    maxConnect::T2      # Maximum connectivity of nodes.
    mobility::T3        # Dislocation mobility.
    remesh::T4          # Remesh flag.
    collision::T4       # Collision flag.
    separation::T4      # Separation flag.
    virtualRemesh::T4   # Virtual remeshing flag.
    parCPU::T4          # Parallelise on CPU flag.
    parGPU::T4          # Parallelise on GPU flag.
    slipStepCritLen::T1 # Critical length for slip step tracking.
    slipStepCritArea::T1    # Critical area for slip step tracking.
end
```
Contains the dislocation network parameters.

The constructor provides a few default values and calculates any derived quantities.
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

The keyword constructor calls the positional one.
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
"""
struct DislocationParameters{T1,T2,T3,T4}
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
    edgeDrag::T1
    screwDrag::T1
    climbDrag::T1
    lineDrag::T1
    maxConnect::T2
    mobility::T3
    remesh::T4
    collision::T4
    separation::T4
    virtualRemesh::T4
    parCPU::T4
    parGPU::T4
    slipStepCritLen::T1
    slipStepCritArea::T1
end
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
struct DislocationLoop{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    loopType::T1    # Loop type.
    numSides::T2    # Number of sides in the loop.
    nodeSide::T2    # Nodes per side of the loop.
    numLoops::T2    # Number of loops to generate when making the network.
    segLen::T3      # Segment lengths.
    slipSystem::T4  # Slip system.
    links::T5       # Links.
    slipPlane::T6   # Slip planes.
    bVec::T6        # Burgers vectors.
    coord::T6       # Coordinates.
    label::T7       # Node labels.
    buffer::T8      # Buffer for distributions.
    range::T9       # Range for distributions.
    dist::T10       # Distribution.
end
```
Stores dislocation loops and parameters used to generate a dislocation network.

There are a couple of different constructors for DislocationLoop. The first one is used for expanding Base.zero().
```
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
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int,
    T3 <: AbstractVector{nodeType},
    T4 <: AbstractDistribution,
}
```

There is also a constructor for `loopPure` loops.
```
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
) where {
    T1 <: loopPure,
    T2 <: Int,
    T3 <: AbstractArray{T,N} where {T,N},
    T4 <: AbstractVector{nodeType},
    T5 <: AbstractDistribution,
}
```

A fallback for other as of yet unimplemented loop types.
```
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
) where {
    T1 <: loopImpure,
    T2 <: Int,
    T3 <: AbstractArray{T,N} where {T,N},
    T4 <: AbstractVector{nodeType},
    T5 <: AbstractDistribution,
}
```

And a keyword constructor that calls the positional constructor, which dispatches on the appropriate method.
```
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
end
```
"""
struct DislocationLoop{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
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
    range::T9
    dist::T10
end
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
dist,
    )

end
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
        dist,
    )
end







"""
Dislocation network.
```
struct DislocationNetwork{T1, T2, T3, T4, T5, T6}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeVel::T2
    nodeForce::T2
    numNodeSegConnect::T4   # Number of nodes, segments and max connectivity in network
    connectivity::T5        # Connectivity matrix
    linksConnect::T5        # Links involved in connection
    segIdx::T5              # Contains segment index and the nodes of the nodes in said link
    segForce::T6            # Force on each node of each segment
end
```
"""
struct DislocationNetwork{T1,T2,T3,T4,T5,T6}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeVel::T2
    nodeForce::T2
    numNode::T4
    numSeg::T4
    maxConnect::T5
    connectivity::T1
    linksConnect::T1
    segIdx::T1
    segForce::T6

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

        return new{T1,T2,T3,typeof(numNodeArr),T5,T6}(
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
end
