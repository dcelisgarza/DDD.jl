"""
```
makeSegment(type::AbstractDlnSeg, slipPlane::Vector{T}, bVec::Vector{T}) where {T<:Float64}
```
Make segment depending on the segment type. The type is known at compile type so the compiler will delete unecessary branches.
"""
function makeSegment(
    type::segEdge,
    slipPlane::Vector{T},
    bVec::Vector{T},
) where {T <: Float64}
    edge = cross(slipPlane, bVec)
    return edge ./ norm(edge)
end
function makeSegment(
    type::segEdgeN,
    slipPlane::Vector{T},
    bVec::Vector{T},
) where {T <: Float64}
    return slipPlane ./ norm(slipPlane)
end
function makeSegment(
    type::segScrew,
    slipPlane::Vector{T},
    bVec::Vector{T},
) where {T <: Float64}
    return bVec ./ norm(bVec)
end
"""
```
makeLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segType::T3,
    segLen::T4,
    slipSystem::T2,
    _slipPlane::T5,
    _bVec::T5,
    label::T6,
    buffer::T7,
    range::T5,
    dist::T8,
) where {
    T1 <: loopDln,
    T2 <: Int64,
    T3 <: segNone,
    T4 <: Float64,
    T5 <: AbstractArray{<:Float64, N} where {N},
    T6 <: Vector{nodeType},
    T7 <: Float64,
    T8 <: AbstractDistribution,
}
```

Constructor function for zero loop.
"""
function makeLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segType::T3,
    segLen::T4,
    slipSystem::T2,
    _slipPlane::T5,
    _bVec::T5,
    label::T6,
    buffer::T7,
    range::T5,
    dist::T8,
) where {
    T1 <: loopDln,
    T2 <: Int64,
    T3 <: segNone,
    T4 <: Float64,
    T5 <: AbstractArray{<:Float64, N} where {N},
    T6 <: Vector{nodeType},
    T7 <: Float64,
    T8 <: AbstractDistribution,
}

    nodeTotal = 0
    numSegType = length(segType)
    links = zeros(Int64, nodeTotal, 2)
    coord = zeros(nodeTotal, 3)
    seg = zeros(numSegType, 3)
    slipPlane = zeros(0, 3)
    bVec = zeros(0, 3)

    return numSides,
    nodeSide,
    numLoops,
    segType,
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
end

"""
```
makeLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segType::T3,
    segLen::T4,
    slipSystem::T2,
    _slipPlane::T5,
    _bVec::T5,
    label::T6,
    buffer::T7,
    range::T8,
    dist::T9,
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{
        T where {T <: AbstractDlnSeg},
        AbstractArray{<:AbstractDlnSeg, N} where {N},
    },
    T4 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T5 <: AbstractArray{<:Float64, N} where {N},
    T6 <: Vector{nodeType},
    T7 <: Float64,
    T8 <: AbstractArray{<:Float64, N} where {N},
    T9 <: AbstractDistribution,
}
```
Constructor function for a loop.
"""
function makeLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segType::T3,
    segLen::T4,
    slipSystem::T2,
    _slipPlane::T5,
    _bVec::T5,
    label::T6,
    buffer::T7,
    range::T8,
    dist::T9,
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{
        T where {T <: AbstractDlnSeg},
        AbstractArray{<:AbstractDlnSeg, N} where {N},
    },
    T4 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T5 <: AbstractArray{<:Float64, N} where {N},
    T6 <: Vector{nodeType},
    T7 <: Float64,
    T8 <: AbstractArray{<:Float64, N} where {N},
    T9 <: AbstractDistribution,
}

    nodeTotal = numSides * nodeSide
    lSegLen = length(segLen)
    @assert length(label) == nodeTotal
    @assert length(segType) == 1
    @assert mod(numSides, 2) == 0
    @assert lSegLen == Int(numSides / 2)

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

    return numSides,
    nodeSide,
    numLoops,
    segType,
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
end

"""
```
DislocationLoop{
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{
        T where {T <: AbstractDlnSeg},
        AbstractArray{<:AbstractDlnSeg, N} where {N},
    },
    T4 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T5 <: Union{Int64, AbstractArray{<:Int64, N} where {N}},
    T6 <: AbstractArray{<:Int64, N} where {N},
    T7 <: AbstractArray{<:Float64, N} where {N},
    T8 <: Vector{<:nodeType},
    T9 <: Float64,
    T10 <: AbstractDistribution,
}
    loopType::T1
    numSides::T2
    nodeSide::T2
    numLoops::T2
    segType::T3
    segLen::T4
    slipSystem::T5
    links::T6
    slipPlane::T7
    bVec::T7
    coord::T7
    label::T8
    buffer::T9
    range::T7
    dist::T10
```
Dislocation loop structure.
"""
struct DislocationLoop{
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{
        T where {T <: AbstractDlnSeg},
        AbstractArray{<:AbstractDlnSeg, N} where {N},
    },
    T4 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T5 <: Union{Int64, AbstractArray{<:Int64, N} where {N}},
    T6 <: AbstractArray{<:Int64, N} where {N},
    T7 <: AbstractArray{<:Float64, N} where {N},
    T8 <: Vector{<:nodeType},
    T9 <: Float64,
    T10 <: AbstractDistribution,
}

    loopType::T1
    numSides::T2
    nodeSide::T2
    numLoops::T2
    segType::T3
    segLen::T4
    slipSystem::T5
    links::T6
    slipPlane::T7
    bVec::T7
    coord::T7
    label::T8
    buffer::T9
    range::T7
    dist::T10

    function DislocationLoop(
        loopType,
        numSides,
        nodeSide,
        numLoops,
        segType,
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
        segType,
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
            segType,
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
            typeof(segType),
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
            segType,
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
length(::DislocationLoop) = 1
getindex(x::DislocationLoop, i::Integer) = i == 1 ? x : throw(BoundsError())
eachindex(x::DislocationLoop) = 1
function zero(::Type{DislocationLoop})
    DislocationLoop(
        loopDln(),
        convert(Int64, 0),
        convert(Int64, 0),
        convert(Int64, 0),
        segNone(),
        convert(Float64, 0),
        convert(Int64, 0),
        zeros(Float64, 0, 3),
        zeros(Float64, 0, 3),
        zeros(nodeType, 0),
        convert(Float64, 0),
        zeros(0, 3),
        Zeros(),
    )
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
    links::T1 # Links.
    slipPlane::T2 # Slip planes.
    bVec::T2 # Burgers vectors.
    coord::T2 # Node coordinates.
    label::T3 # Node labels.
    numNode::T4 # Number of dislocations.
    numSeg::T4 # Number of segments.
    maxConnect::T5
    connectivity::T1
    linksConnect::T1
```
Dislocation Network structure.
"""
mutable struct DislocationNetwork{
    T1 <: AbstractArray{<:Int64, N} where {N},
    T2 <: AbstractArray{<:Float64, N} where {N},
    T3 <: Vector{nodeType},
    T4 <: Int64,
    T5 <: Integer,
}
    links::T1 # Links.
    slipPlane::T2 # Slip planes.
    bVec::T2 # Burgers vectors.
    coord::T2 # Node coordinates.
    label::T3 # Node labels.
    numNode::T4 # Number of dislocations.
    numSeg::T4 # Number of segments.
    maxConnect::T5
    connectivity::T1
    linksConnect::T1

    function DislocationNetwork(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        numNode = 0,
        numSeg = 0,
        maxConnect = 4,
    )

        @assert size(links, 2) == 2
        @assert size(bVec, 2) == size(slipPlane, 2) == size(coord, 2) == 3
        @assert size(links, 1) == size(bVec, 1) == size(slipPlane, 1)
        @assert size(coord, 1) == size(label, 1)

        new{
            typeof(links),
            typeof(bVec),
            typeof(label),
            typeof(numNode),
            typeof(maxConnect),
        }(
            links,
            slipPlane,
            bVec,
            coord,
            label,
            numNode,
            numSeg,
            maxConnect,
        )
    end # Constructor
end # DislocationNetwork

zero(::Type{DislocationNetwork}) = DislocationNetwork(
    zeros(Int64, 0, 2),
    zeros(Float64, 0, 3),
    zeros(Float64, 0, 3),
    zeros(Float64, 0, 3),
    zeros(nodeType, 0),
    convert(Int64, 0),
    convert(Int64, 0),
)

function malloc(network::DislocationNetwork, n::Int64)
    network.links = [network.links; zeros(Int64, n, 2)]
    network.slipPlane = [network.slipPlane; zeros(n, 3)]
    network.bVec = [network.bVec; zeros(n, 3)]
    network.coord = [network.coord; zeros(n, 3)]
    network.label = [network.label; zeros(nodeType, n)]
    return network
end

"""
```
struct DislocationP{
    T1 <: Float64,
    T2 <: Int64,
    T3 <: Bool,
    T4 <: AbstractMobility,
}
    # Size.
    coreRad::T1 # Core radius.
    coreRadMag::T1 # Magnitude of core Radius.
    # Connectivity.
    minSegLen::T1 # Minimum line length.
    maxSegLen::T1 # Maximum line length.
    minArea::T1 # Minimum area for remeshing.
    maxArea::T1 # Maximum area for remeshing.
    maxConnect::T2 # Maximum number of connections to a node.
    remesh::T3 # Flag for remeshing.
    collision::T3 # Flag for collision handling.
    separation::T3 # Flag for separation handling.
    virtualRemesh::T3 # Flag for virtual remeshing.
    # Mobility.
    edgeDrag::T1 # Drag coefficient edge dislocation.
    screwDrag::T1 # Drag coefficient screw dislocation.
    climbDrag::T1 # Drag coefficient climb.
    lineDrag::T1 # Drag coefficient line.
    mobility::T4 # Mobility law.
```
Dislocation parameters structure.
"""
struct DislocationP{
    T1 <: Float64,
    T2 <: Int64,
    T3 <: Bool,
    T4 <: AbstractMobility,
}
    # Size.
    coreRad::T1 # Core radius.
    coreRadMag::T1 # Magnitude of core Radius.
    # Connectivity.
    minSegLen::T1 # Minimum line length.
    maxSegLen::T1 # Maximum line length.
    minArea::T1 # Minimum area for remeshing.
    maxArea::T1 # Maximum area for remeshing.
    maxConnect::T2 # Maximum number of connections to a node.
    remesh::T3 # Flag for remeshing.
    collision::T3 # Flag for collision handling.
    separation::T3 # Flag for separation handling.
    virtualRemesh::T3 # Flag for virtual remeshing.
    # Mobility.
    edgeDrag::T1 # Drag coefficient edge dislocation.
    screwDrag::T1 # Drag coefficient screw dislocation.
    climbDrag::T1 # Drag coefficient climb.
    lineDrag::T1 # Drag coefficient line.
    mobility::T4 # Mobility law.
    # Fool-proof constructor.

    function DislocationP(
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

        new{
            typeof(coreRad),
            typeof(maxConnect),
            typeof(remesh),
            typeof(mobility),
        }(
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
    end # constructor
end # DislocationP
