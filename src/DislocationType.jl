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
) where {T<:Float64}
    edge = cross(slipPlane, bVec)
    return edge ./ norm(edge)
end
function makeSegment(
    type::segEdgeN,
    slipPlane::Vector{T},
    bVec::Vector{T},
) where {T<:Float64}
    return slipPlane ./ norm(slipPlane)
end
function makeSegment(
    type::segScrew,
    slipPlane::Vector{T},
    bVec::Vector{T},
) where {T<:Float64}
    return bVec ./ norm(bVec)
end
function makeSegment(
    type::segNone,
    slipPlane::Vector{T},
    bVec::Vector{T},
) where {T<:Float64}
    return [0.0; 0.0; 0.0]
end
"""
```
DislocationLoop{
    T1<:loopSides,
    T2<:Int64,
    T3<:Union{
        T where {T<:AbstractDlnSeg},
        AbstractArray{<:AbstractDlnSeg,N} where {N},
    },
    T4<:Union{T where {T<:Float64},AbstractArray{<:Float64,N} where {N}},
    T5<:Union{Int64,AbstractArray{<:Int64,N} where {N}},
    T6<:AbstractArray{<:Int64,N} where {N},
    T7<:AbstractArray{<:Float64,N} where {N},
    T8<:Vector{<:nodeType},
    T9<:Float64,
    T10<:AbstractArray{<:Float64,N} where {N},
}
    numSides::T1
    nodeSide::T2
    segType::T3
    segLen::T4
    slipSystem::T5
    links::T6
    slipPlane::T7
    bVec::T7
    coord::T7
    label::T8
    numLoops::T2
    buffer::T9
    range::T10
```
Idealised dislocation loop to be used as sources.
"""
function makeLoop(
    loopType::loopDln,
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
    nodeTotal = numSides * nodeSide
    numSegType = length(segType)
    @assert length(label) == nodeTotal
    @assert numSegType == length(segLen) == numSides / 2 == size(_bVec, 1) ==
            size(_slipPlane, 1) == length(slipSystem)
    links = zeros(Int64, nodeTotal, 2)
    coord = zeros(nodeTotal, 3)
    seg = zeros(numSegType, 3)
    slipPlane = zeros(0, 3)
    bVec = zeros(0, 3)

    for i = 1:numSegType
        seg[i, :] = makeSegment(segType[i], _slipPlane[i, :], _bVec[i, :]) *
                    segLen[i]
    end
    # """
    # This if statement is ripe for refactoring but I need properly to think about it because it doesn't seem trivial.
    # """
    if numSides == 4
        # First node.
        coord[1, :] -= seg[1, :] + seg[2, :]
        slipPlane = [slipPlane; _slipPlane[2, :]']
        bVec = [bVec; _bVec[2, :]']
        for k = 1:nodeSide
            coord[1+k, :] = coord[k, :] + seg[1, :]
            slipPlane = [slipPlane; _slipPlane[1, :]']
            bVec = [bVec; _bVec[1, :]']
        end
        for k = 1:nodeSide
            coord[1 + k + nodeSide, :] = coord[k+nodeSide, :] + seg[2, :]
            slipPlane = [slipPlane; _slipPlane[2, :]']
            bVec = [bVec; _bVec[2, :]']
        end
        for k = 1:nodeSide
            coord[1 + k + 2 * nodeSide, :] = coord[k+2*nodeSide, :] - seg[1, :]
            slipPlane = [slipPlane; _slipPlane[1, :]']
            bVec = [bVec; _bVec[1, :]']
        end
        if nodeSide > 1
            # Last node
            coord[2+3*nodeSide, :] = coord[1+3*nodeSide, :] - seg[2, :]
            slipPlane = [slipPlane; _slipPlane[2, :]']
            bVec = [bVec; _bVec[2, :]']
        end
    elseif numSides == 6
        # First node.
        coord[1, :] -= seg[1, :] + seg[2, :] + seg[3, :]
        slipPlane = [slipPlane; _slipPlane[3, :]']
        bVec = [bVec; _bVec[3, :]']
        for k = 1:nodeSide
            coord[1+k, :] = coord[k, :] + seg[1, :]
            slipPlane = [slipPlane; _slipPlane[1, :]']
            bVec = [bVec; _bVec[1, :]']
        end
        for k = 1:nodeSide
            coord[1 + k + nodeSide, :] = coord[k+nodeSide, :] + seg[2, :]
            slipPlane = [slipPlane; _slipPlane[2, :]']
            bVec = [bVec; _bVec[2, :]']
        end
        for k = 1:nodeSide
            coord[1 + k + 2 * nodeSide, :] = coord[k+2*nodeSide, :] + seg[3, :]
            slipPlane = [slipPlane; _slipPlane[3, :]']
            bVec = [bVec; _bVec[3, :]']
        end
        for k = 1:nodeSide
            coord[1 + k + 3 * nodeSide, :] = coord[k+3*nodeSide, :] - seg[1, :]
            slipPlane = [slipPlane; _slipPlane[1, :]']
            bVec = [bVec; _bVec[1, :]']
        end
        for k = 1:nodeSide
            coord[1 + k + 4 * nodeSide, :] = coord[k+4*nodeSide, :] - seg[2, :]
            slipPlane = [slipPlane; _slipPlane[2, :]']
            bVec = [bVec; _bVec[2, :]']
        end
        if nodeSide > 1
            for k = 1:nodeSide-1
                coord[1 + k + 5 * nodeSide, :] = coord[k+5*nodeSide, :] -
                                                 seg[3, :]
                slipPlane = [slipPlane; _slipPlane[3, :]']
                bVec = [bVec; _bVec[3, :]']
            end
        end
        # """
        # # Room for expansion in case we can generate loops with more sides.
        # else
        #     error("more sides for a source loop are undefined")
        # """
    end
    # Links
    for j = 1:nodeTotal-1
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





struct DislocationLoop{
    T0<:AbstractDlnStr,
    T1<:loopSides,
    T2<:Int64,
    T3<:Union{
        T where {T<:AbstractDlnSeg},
        AbstractArray{<:AbstractDlnSeg,N} where {N},
    },
    T4<:Union{T where {T<:Float64},AbstractArray{<:Float64,N} where {N}},
    T5<:Union{Int64,AbstractArray{<:Int64,N} where {N}},
    T6<:AbstractArray{<:Int64,N} where {N},
    T7<:AbstractArray{<:Float64,N} where {N},
    T8<:Vector{<:nodeType},
    T9<:Float64,
    T10<:AbstractArray{<:Float64,N} where {N},
    T11<:AbstractDistribution,
}
    loopType::T0
    numSides::T1
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
    range::T10
    dist::T11
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
            typeof(nodeSide),
            typeof(segType),
            typeof(segLen),
            typeof(slipSystem),
            typeof(links),
            typeof(slipPlane),
            typeof(label),
            typeof(buffer),
            typeof(range),
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
function zero(::Type{DislocationLoop})
    DislocationLoop(
        loopDln(),
        loopSides(4),
        convert(Int64, 1),
        convert(Int64, 0),
        [segNone(), segNone()],
        zeros(Float64, 2),
        zeros(Int64, 2),
        zeros(Float64, 2, 3),
        zeros(Float64, 2, 3),
        zeros(nodeType, 4),
        convert(Float64, 0),
        zeros(2, 3),
        Zeros(),
    )
end

mutable struct DislocationNetwork{
    T1<:AbstractArray{<:Int64,N} where {N},
    T2<:AbstractArray{<:Float64,N} where {N},
    T3<:Vector{nodeType},
    T4<:Int64,
}
    links::T1 # Links.
    slipPlane::T2 # Slip planes.
    bVec::T2 # Burgers vectors.
    coord::T2 # Node coordinates.
    label::T3 # Node labels.
    numNode::T4 # Number of dislocations.
    numSeg::T4 # Number of segments.
    function DislocationNetwork(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        numNode = 0,
        numSeg = 0,
    )
        @assert size(links, 2) == 2
        @assert size(bVec, 2) == size(slipPlane, 2) == size(coord, 2) == 3
        @assert size(links, 1) == size(bVec, 1) == size(slipPlane, 1)
        @assert size(coord, 1) == size(label, 1)
        new{typeof(links),typeof(bVec),typeof(label),typeof(numNode)}(
            links,
            slipPlane,
            bVec,
            coord,
            label,
            numNode,
            numSeg,
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

struct DislocationP{T1<:Float64,T2<:Int64,T3<:Bool,T4<:AbstractMobility}
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
