"""
Dislocation nodes have labels that change how they are treated by the simulation. There are only given types of nodes so these labels may only take on predefined values and error for anything else.
```
@enum nodeType begin
    undef = 0  # Undefined node, value at initialisation.
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
There are different types of segments. Edge segments are orthogonal to the Burgers vector; screw segments are parallel to the Burgers vector; mixed segments are anything in between. Segment type can be easily inferred during runtime. These are mainly used for multiple dispatch purposes.
```
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end    # Undefined segment.
struct segEdge <: AbstractDlnSeg end    # Edge segment.
struct segEdgeN <: AbstractDlnSeg end   # Edge segment parallel to slip plane.
struct segScrew <: AbstractDlnSeg end   # Screw segment.
struct segMixed <: AbstractDlnSeg end   # Mixed segment.
```
"""
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end
struct segEdge <: AbstractDlnSeg end
struct segEdgeN <: AbstractDlnSeg end
struct segScrew <: AbstractDlnSeg end
struct segMixed <: AbstractDlnSeg end
"""
Dislocation structures have different classifications. Prismatic loops are made up only of edge segments with the same slip system; shear loops are made up of a mixture of segment types with the same slip system; jogs and kinks are steps not contained in the slip plane.
```
abstract type AbstractDlnStr end
struct loopDln <: AbstractDlnStr end    # Unclassified loop.
struct loopPrism <: AbstractDlnStr end  # Prismatic loop.
struct loopShear <: AbstractDlnStr end  # Shear loop.
struct loopJog <: AbstractDlnStr end    # Jog.
struct loopKink <: AbstractDlnStr end   # Kink.
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
struct BCC <: AbstractMobility end
struct FCC <: AbstractMobility end
struct HCP <: AbstractMobility end
```
Mobility functions.
"""
abstract type AbstractMobility end
struct mobBCC <: AbstractMobility end
struct mobFCC <: AbstractMobility end
struct mobHCP <: AbstractMobility end
"""
```
struct SlipSystem{
    T1 <: AbstractString,
    T2 <: AbstractArray{<:Float64, N} where {N},
}
    name::T1
    slipPlane::T2
    bVec::T2
```
Slip systems.
"""
struct SlipSystem{
    T1 <: AbstractCrystalStruct,
    T2 <: AbstractArray{<:Float64, N} where {N},
}
    crystalStruct::T1
    slipPlane::T2
    bVec::T2
    function SlipSystem(; crystalStruct, slipPlane, bVec)
        new{typeof(crystalStruct), typeof(slipPlane)}(
            crystalStruct,
            slipPlane,
            bVec,
        )
    end
end

"""
```
DislocationP{
    T1 <: Float64,
    T2 <: Int64,
    T3 <: Bool,
    T4 <: AbstractMobility,
}
    # Size.
    coreRad::T1     # Core radius.
    coreRadMag::T1  # Magnitude of core Radius.
    # Connectivity.
    minSegLen::T1       # Minimum line length.
    maxSegLen::T1       # Maximum line length.
    minArea::T1         # Minimum area for remeshing.
    maxArea::T1         # Maximum area for remeshing.
    maxConnect::T2      # Maximum number of connections to a node.
    remesh::T3          # Flag for remeshing.
    collision::T3       # Flag for collision handling.
    separation::T3      # Flag for separation handling.
    virtualRemesh::T3   # Flag for virtual remeshing.
    # Mobility.
    edgeDrag::T1    # Drag coefficient edge dislocation.
    screwDrag::T1   # Drag coefficient screw dislocation.
    climbDrag::T1   # Drag coefficient climb.
    lineDrag::T1    # Drag coefficient line.
    mobility::T4    # Mobility law.
```
Dislocation parameters structure. See [`AbstractMobility`](@ref) for more details.
"""
struct DislocationP{
    T1 <: Float64,
    T2 <: Int64,
    T3 <: Bool,
    T4 <: AbstractMobility,
}

    coreRad::T1
    coreRadSq::T1
    coreRadMag::T1
    minSegLen::T1
    maxSegLen::T1
    minArea::T1
    maxArea::T1
    maxConnect::T2
    remesh::T3
    collision::T3
    separation::T3
    virtualRemesh::T3
    edgeDrag::T1
    screwDrag::T1
    climbDrag::T1
    lineDrag::T1
    mobility::T4

    function DislocationP(;
        coreRad,
        coreRadMag,
        minSegLen,
        maxSegLen,
        minArea,
        maxArea,
        maxConnect,
        remesh,
        collision,
        separation,
        virtualRemesh,
        edgeDrag,
        screwDrag,
        climbDrag,
        lineDrag,
        mobility,
    )

        coreRad == minSegLen == maxSegLen == 0 ? nothing :
        @assert coreRad < minSegLen < maxSegLen
        minArea == maxArea == 0 ? nothing : @assert minArea < maxArea
        coreRadSq = coreRad^2
        new{
            typeof(coreRad),
            typeof(maxConnect),
            typeof(remesh),
            typeof(mobility),
        }(
            coreRad,
            coreRadSq,
            coreRadMag,
            minSegLen,
            maxSegLen,
            minArea,
            maxArea,
            maxConnect,
            remesh,
            collision,
            separation,
            virtualRemesh,
            edgeDrag,
            screwDrag,
            climbDrag,
            lineDrag,
            mobility,
        )
    end # constructor
end # DislocationP

"""
```
struct DislocationLoop{
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T4 <: Union{T where {T <: Int64}, AbstractArray{<:Int64, N} where {N}},
    T5 <: AbstractArray{<:Int64, N} where {N},
    T6 <: AbstractArray{<:Float64, N} where {N},
    T7 <: Vector{<:nodeType},
    T8 <: Float64,
    T9 <: AbstractDistribution,
}

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
```
Dislocation loop structure generated via the constructor [`makeLoop`](@ref).
"""
struct DislocationLoop{
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T4 <: Union{T where {T <: Int64}, AbstractArray{<:Int64, N} where {N}},
    T5 <: AbstractArray{<:Int64, N} where {N},
    T6 <: AbstractArray{<:Float64, N} where {N},
    T7 <: Vector{<:nodeType},
    T8 <: Float64,
    T9 <: AbstractDistribution,
}

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

    function DislocationLoop(;
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
        dist, = makeLoop(
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

        new{
            typeof(loopType),
            typeof(numSides),
            typeof(segLen),
            typeof(slipSystem),
            typeof(links),
            typeof(slipPlane),
            typeof(label),
            typeof(buffer),
            typeof(dist),
        }(
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
end

"""
```
DislocationNetwork{
    T1 <: AbstractArray{<:Int64, N} where {N},
    T2 <: AbstractArray{<:Float64, N} where {N},
    T3 <: Vector{nodeType},
    T4 <: Int64,
    T5 <: Integer,
}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeForce::T2
    nodeVel::T2
    numNode::T4 = 0     # Total number of nodes in network.
    numSeg::T4 = 0      # Total number of segs in network.
    maxConnect::T5 = 4  # Maximum connectivity of nodes.
    connectivity::T1
    linksConnect::T1
    segIdx::T1          # segIdx[:,1] is the segment index. Used to find the bVec and slipPlane of a real segment. segIdx[:,2:3] are the indices of the nodes involved in a given link, used to find their coordinates.
```
Dislocation Network structure. See [`DislocationLoop`](@ref), [`makeNetwork`](@ref) and [`makeNetwork!`](@ref) for further details.
"""
mutable struct DislocationNetwork{
    T1 <: AbstractArray{<:Int64, N} where {N},
    T2 <: AbstractArray{<:Float64, N} where {N},
    T3 <: Vector{nodeType},
    T4 <: Int64,
}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeForce::T2
    nodeVel::T2
    numNode::T4
    numSeg::T4
    maxConnect::T4
    connectivity::T1
    linksConnect::T1
    segIdx::T1

    function DislocationNetwork(;
        links,
        slipPlane,
        bVec,
        coord,
        label,
        nodeForce,
        nodeVel,
        numNode = 0,
        numSeg = 0,
        maxConnect = 0,
    )

        @assert size(links, 2) == 2
        @assert size(bVec, 2) == size(slipPlane, 2) == size(coord, 2) == 3
        @assert size(links, 1) == size(bVec, 1) == size(slipPlane, 1)
        @assert size(coord, 1) == size(label, 1)

        new{typeof(links), typeof(bVec), typeof(label), typeof(numNode)}(
            links,
            slipPlane,
            bVec,
            coord,
            label,
            nodeForce,
            nodeVel,
            numNode,
            numSeg,
            maxConnect,
        )
    end # Constructor
end # DislocationNetwork
