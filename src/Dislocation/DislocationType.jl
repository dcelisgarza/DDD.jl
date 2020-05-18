## Primitives
"""
Dislocation nodes have labels that change how they are treated by the simulation. There are only certain types of nodes so these labels may only take on predefined values. Adding new nodes requires adding new entries to the enumerated type.
```
@enum nodeType begin
    none = 0  # Undefined node, value at initialisation.
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

## Structures

"""
```
struct SlipSystem{T1 <: AbstractCrystalStruct, T2 <: AbstractArray{T, N} where {T, N}}
    crystalStruct::T1
    slipPlane::T2
    bVec::T2
end
```
Slip systems.
"""
struct SlipSystem{T1, T2}
    crystalStruct::T1
    slipPlane::T2
    bVec::T2
end
function SlipSystem(;
    crystalStruct::T1,
    slipPlane::T2,
    bVec::T2,
) where {T1 <: AbstractCrystalStruct, T2 <: AbstractArray{T, N} where {T, N}}
    SlipSystem(crystalStruct, slipPlane, bVec)
end

"""
```
DislocationP{T1, T2, T3, T4}
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
struct DislocationP{T1, T2, T3, T4}
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
end
"""
```
```
"""
function DislocationP(;
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
) where {T1 <: Float64, T2 <: Int, T3 <: Bool, T4 <: AbstractMobility}

    coreRad == minSegLen == maxSegLen == 0 ? nothing : @assert coreRad < minSegLen < maxSegLen
    minArea == maxArea == 0 ? nothing : @assert minArea < maxArea
    coreRadSq = coreRad^2

    DislocationP(
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

"""
```
struct DislocationLoop{
    T1 <: AbstractDlnStr,
    T2 <: Int,
    T3 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T4 <: Union{T where {T <: Int}, AbstractArray{<:Int, N} where {N}},
    T5 <: AbstractArray{<:Int, N} where {N},
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
    T6 <: Float64,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution,
}
```
Generic named constructor for dislocation loop. Calls other constructors that dispatch on `loopType`.
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
    T6 <: Float64,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution,
}

    DislocationLoop(
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
    T6 <: Float64,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution,
}

    nodeTotal::Int = 0
    links = zeros(Int, nodeTotal, 2)
    coord = zeros(nodeTotal, 3)
    slipPlane = zeros(0, 3)
    bVec = zeros(0, 3)

    DislocationLoop(
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
    T6 <: Float64,
    T7 <: AbstractArray{T, N} where {T, N},
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

    links = zeros(Int, nodeTotal, 2)
    coord = zeros(nodeTotal, 3)
    slipPlane = zeros(0, 3)
    bVec = zeros(0, 3)
    seg = zeros(lSegLen, 3)
    @inbounds for i in eachindex(segLen)
        seg[i, :] = makeSegment(segEdge(), _slipPlane, _bVec) .* segLen[i]
    end

    θ = extAngle(numSides)
    rseg = zeros(3)
    @inbounds for i in 1:numSides
        idx = (i - 1) * nodeSide
        rseg = rot3D(seg[mod(i - 1, lSegLen) + 1, :], rotAxis, zeros(3), θ * (i - 1))
        for j in 1:nodeSide
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
    @inbounds for j in 1:(nodeTotal - 1)
        links[j, :] = [j; 1 + j]
    end
    links[nodeTotal, :] = [nodeTotal; 1]

    DislocationLoop(
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
DislocationNetwork{
    T1 <: AbstractArray{<:Int, N} where {N},
    T2 <: AbstractArray{<:Float64, N} where {N},
    T3 <: Vector{nodeType},
    T4 <: Int,
    T5 <: Int,
}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    segForce::T2
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
mutable struct DislocationNetwork{T1, T2, T3, T4, T5}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    segForce::T2
    nodeVel::T2
    numNode::T4
    numSeg::T4
    maxConnect::T4
    connectivity::T5
    linksConnect::T5
    segIdx::T5
end # DislocationNetwork
@inline function DislocationNetwork(;
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    segForce::T2,
    nodeVel::T2,
    numNode::T4 = 0,
    numSeg::T4 = 0,
    maxConnect::T4 = 0,
    connectivity::T5 = zeros(Int, 0, 0),
    linksConnect::T5 = zeros(Int, 0, 0),
    segIdx::T5 = zeros(Int, 0, 0),
) where {
    T1 <: AbstractArray{T, N} where {T, N},
    T2 <: AbstractArray{T, N} where {T, N},
    T3 <: AbstractVector{nodeType},
    T4 <: Int,
    T5 <: AbstractArray{Int, N} where {N},
}

    @assert size(links, 2) == 2
    @assert size(bVec, 2) == size(slipPlane, 2) == size(coord, 2) == 3
    @assert size(links, 1) == size(bVec, 1) == size(slipPlane, 1)
    @assert size(coord, 1) == size(label, 1)

    DislocationNetwork(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        segForce,
        nodeVel,
        numNode,
        numSeg,
        maxConnect,
        connectivity,
        linksConnect,
        segIdx,
    )
end # Constructor

"""
```
DislocationNetwork(
    sources::T1,
    maxConnect::T2 = 4,
    args...;
    memBuffer::T2 = 10,
    checkConsistency::T3 = true,
    kw...,
) where {T1 <: Union{DislocationLoop, AbstractVector{<:DislocationLoop}}, T2 <: Int, T3 <: Bool}
```
Constructor for [`DislocationNetwork`](@ref), see [`DislocationNetwork!`](@ref) for in-place version.
"""
@inline function DislocationNetwork(
    sources::T1,
    maxConnect::T2 = 4,
    args...;
    memBuffer::T2 = 10,
    checkConsistency::T3 = true,
    kw...,
) where {T1 <: Union{DislocationLoop, AbstractVector{<:DislocationLoop}}, T2 <: Int, T3 <: Bool}
    nodeTotal::Int = 0
    lims = zeros(2, 3)
    # Allocate memory.
    @inbounds for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end

    nodeBuffer::Int = nodeTotal * memBuffer

    links = zeros(Int, nodeBuffer, 2)
    slipPlane = zeros(nodeBuffer, 3)
    bVec = zeros(nodeBuffer, 3)
    coord = zeros(nodeBuffer, 3)
    label = zeros(nodeType, nodeBuffer)
    segForce = zeros(Float64, nodeBuffer, 3)
    nodeVel = zeros(Float64, nodeBuffer, 3)
    numNode = nodeTotal
    numSeg = nodeTotal

    maxConnect = maxConnect

    nodeTotal = 0
    initIdx = findfirst(x -> x == 0, label)
    initIdx == nothing ? initIdx = 0 : nothing
    @inbounds for i in eachindex(sources)
        # Indices.
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        disp = loopDistribution(sources[i].dist, numLoops, args...; kw...)
        limits!(lims, mean(sources[i].segLen), sources[i].range, sources[i].buffer)
        for j in 1:numLoops
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            # Prepare to distribute sources.
            links[idxi:idxf, :] = sources[i].links[1:nodesLoop, :] .+ (nodeTotal + initIdx - 1)
            slipPlane[idxi:idxf, :] .= sources[i].slipPlane[1:nodesLoop, :]
            bVec[idxi:idxf, :] .= sources[i].bVec[1:nodesLoop, :]
            coord[idxi:idxf, :] .= sources[i].coord[1:nodesLoop, :]
            # Move loop.
            coord[idxi:idxf, :] .=
                translatePoints(sources[i].coord[1:nodesLoop, :], lims, disp[j, :])
            label[idxi:idxf, :] .= sources[i].label[1:nodesLoop, :]
            nodeTotal += nodesLoop
        end
    end

    numSeg, segIdx = getSegmentIdx(links, label)
    connectivity, linksConnect = makeConnect(links, maxConnect)

    network = DislocationNetwork(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        segForce,
        nodeVel,
        numNode,
        numSeg,
        maxConnect,
        connectivity,
        linksConnect,
        segIdx,
    )

    checkConsistency ? checkNetwork(network) : nothing

    return network
end # Constructor

"""
```
DislocationNetwork!(
    network::DislocationNetwork,
    sources::Union{
        DislocationLoop,
        AbstractVector{<:DislocationLoop}
    },
    maxConnect::Int = 4,
    args...;
    checkConsistency::Bool = false,
    kw...,
)
```
In-place constructor for [`DislocationNetwork`](@ref), see [`DislocationNetwork`](@ref) for constructor.
"""
@inline function DislocationNetwork!(
    network::T1,
    sources::T2,
    maxConnect::T3 = 4,
    args...;
    checkConsistency::T4 = false,
    kw...,
) where {
    T1 <: DislocationNetwork,
    T2 <: Union{DislocationLoop, AbstractVector{<:DislocationLoop}},
    T3 <: Int,
    T4 <: Bool,
}
    nodeTotal::Int = 0
    lims = zeros(2, 3)
    network.maxConnect = maxConnect
    # Allocate memory.
    @inbounds for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end
    available = findfirst(x -> x == 0, network.label)
    if available == nothing
        push!(network, nodeTotal)
    else
        check = length(network.label) - (available - 1) - nodeTotal
        check < 0 ? push!(network, -check) : nothing
    end

    nodeTotal = 0
    initIdx = findfirst(x -> x == 0, network.label)
    initIdx == nothing ? initIdx = 0 : nothing
    @inbounds for i in eachindex(sources)
        # Indices.
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        disp = loopDistribution(sources[i].dist, numLoops, args...; kw...)
        limits!(lims, mean(sources[i].segLen), sources[i].range, sources[i].buffer)
        for j in 1:numLoops
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            # Prepare to distribute sources.
            network.links[idxi:idxf, :] =
                sources[i].links[1:nodesLoop, :] .+ (nodeTotal + initIdx - 1)
            network.slipPlane[idxi:idxf, :] .= sources[i].slipPlane[1:nodesLoop, :]
            network.bVec[idxi:idxf, :] .= sources[i].bVec[1:nodesLoop, :]
            network.coord[idxi:idxf, :] .= sources[i].coord[1:nodesLoop, :]
            # Move loop.
            network.coord[idxi:idxf, :] .=
                translatePoints(sources[i].coord[1:nodesLoop, :], lims, disp[j, :])
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
