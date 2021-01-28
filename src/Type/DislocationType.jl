"""
```
@enum nodeTypeDln begin
    noneDln = 0    # Undefined node, value at initialisation
    intMobDln = 1  # Internal mobile node
    intFixDln = 2  # Internal fixed node
    srfMobDln = 3  # Mobile surface node
    srfFixDln = 4  # Fixed surface node
    extDln = 5     # External node
    tmpDln = 6     # Temporary flag, used during topological operations
end
```
Different types of nodes behave differently. There are only a finite number of them so an enumerated type provides safety and efficiency. Each value represents a different type of node and therefore its behaviour.

# Meaning
* `noneDln` are uninitialised nodes.
* `intMobDln` are mobile nodes internal to the convex hull of the domain. They take part in tractions, displacements and dislocation interactions.
* `intFixDln` are fixed nodes internal to the convex hull of the domain. They participate in the same way as `intMobDln` nodes except for the fact that their velocities is fixed are zero.
* `srfMobDln` are mobile nodes that live on the surface of the convex hull of the domain, they are used to track slip steps and therefore participate in the same things as internal nodes but their velocities are restricted to the convex hull surface.
* `srfFixDln` are fixed surface nodes and have the same characteristics as mobile surface nodes except for having zero velocity.
* `extDln` are external nodes that do not participate in dislocation interactions or forces but are used to calculate displacements and track slip steps.
* `tmpDln` are nodes that are temporarily flagged before they are assigned another type.
"""
@enum nodeTypeDln begin
    noneDln = 0
    intMobDln = 1
    intFixDln = 2
    srfMobDln = 3
    srfFixDln = 4
    extDln = 5
    tmpDln = 6
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

These types are used to automatically generate segments out of Burgers vectors `b`, slip planes `n`, and/or line direction `l`.

* `segEdge` have `b ⟂ t`,
* `segEdgeN` have `b ⟂ t` and `b ∥ n` ,
* `segScrew` have `b ∥ t` ,
* `segMixed` have noneDln of the above.
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
const loopDefined = Union{loopPure,loopImpure}
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
```
struct SlipSystem{T1,T2,T3}
    crystalStruct::T1
    slipPlane::T2
    bVec::T3
end
```
Stores slip systems. 
"""
struct SlipSystem{T1,T2,T3}
    crystalStruct::T1
    slipPlane::T2
    bVec::T3
end

"""
```
struct DislocationParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21}
    mobility::T1            # Dislocation core radius
    dragCoeffs::T2          # Drag coefficients
    coreRad::T3             # Dislocation core radius
    coreRadSq::T4           # Dislocation core radius squared
    coreRadMag::T5          # Magnitude of core radius
    coreEnergy::T6          # Core energy
    minSegLen::T7           # Minimum segment length
    maxSegLen::T8           # Maximum segment length
    twoMinSegLen::T9        # Minimum segment length times two
    minArea::T10            # Minimum area for remeshing
    maxArea::T11            # Maximum area for remeshing
    minAreaSq::T12          # Minimum area for remeshing squared
    maxAreaSq::T13          # Maximum area for remeshing squared
    slipStepCritLen::T14    # Critical length for slip step tracking
    slipStepCritArea::T15   # Critical area for slip slep tracking
    remesh::T16
    collision::T17
    separation::T18
    virtualRemesh::T19
    parCPU::T20
    parGPU::T21
end
struct DislocationParameters{T1,T2,T3,T4}
    coreRad::T1         
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
Stores the dislocation parameters.
"""
struct DislocationParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21}
    mobility::T1
    dragCoeffs::T2
    coreRad::T3
    coreRadSq::T4
    coreRadMag::T5
    coreEnergy::T6
    minSegLen::T7
    maxSegLen::T8
    twoMinSegLen::T9
    minArea::T10
    maxArea::T11
    minAreaSq::T12
    maxAreaSq::T13
    slipStepCritLen::T14
    slipStepCritArea::T15
    remesh::T16
    collision::T17
    separation::T18
    virtualRemesh::T19
    parCPU::T20
    parGPU::T21
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
Stores a dislocation loop and parameters used to generate a [`DislocationNetwork`](@ref).
"""
struct DislocationLoop{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14}
    loopType::T1
    numSides::T2
    nodeSide::T3
    numLoops::T4
    segLen::T5
    slipSystem::T6
    links::T7
    slipPlane::T8
    bVec::T9
    coord::T10
    label::T11
    buffer::T12
    range::T13
    dist::T14
end
const DislocationLoopCollection = Union{T,AbstractVector{T},NTuple{N,T} where N} where {T <: DislocationLoop}

"""
```
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
end
```
Stores the dislocation network generated from [`DislocationLoop`](@ref).
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
end