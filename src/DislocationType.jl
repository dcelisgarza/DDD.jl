"""
```
@enum nodeType begin
    undef = -1
    intMob = 0
    intFix = 1
    srfMob = 2
    srfFix = 3
    ext = 4
end
```
Enumerated type for dislocation nodes.
"""
@enum nodeType begin
    undef = -1
    intMob = 0
    intFix = 1
    srfMob = 2
    srfFix = 3
    ext = 4
end
"""
Overloaded functions for dislocation `nodeType`.
"""
isequal(x::Real, y::nodeType) = isequal(x, Integer(y))
isequal(x::nodeType, y::Real) = isequal(Integer(x), y)
isless(x::Real, y::nodeType) = isless(x, Integer(y))
isless(x::nodeType, y::Real) = isless(Integer(x), y)
==(x::nodeType, y::Real) = isequal(Integer(x), y)
==(x::Real, y::nodeType) = isequal(x, Integer(y))
convert(::Type{nodeType}, x::Real) = nodeType(Integer(x))
zero(::Type{nodeType}) = -1
"""
```
abstract type AbstractDlnSeg end
struct dlnEdge <: AbstractDlnSeg end
struct dlnScrew <: AbstractDlnSeg end
struct dlnMixed <: AbstractDlnSeg end
```
Edge types.
"""
abstract type AbstractDlnSeg end
struct dlnEdge <: AbstractDlnSeg end
struct dlnEdgeN <: AbstractDlnSeg end
struct dlnScrew <: AbstractDlnSeg end
struct dlnMixed <: AbstractDlnSeg end
"""
```
```
"""
function makeSegment(
    type::AbstractDlnSeg,
    slipPlane::Vector{T},
    bVec::Vector{T},
) where {T<:Real}
    if typeof(type) == dlnEdge
        seg = cross(slipPlane, bVec)
    elseif typeof(type) == dlnEdgeN
        seg = slipPlane
    elseif typeof(type) == dlnScrew
        seg = bVec
    else
        error("makeSegment: mixed dislocation segment undefined")
    end
    return seg ./= norm(seg)
end
"""
```
@enum loopSides begin
    four = 4
    six = 6
end
```
Type for number of sides for idealised loops.
"""
@enum loopSides begin
    four = 4
    six = 6
end
"""
Overloaded functions for `loopSides`.
"""
isequal(x::Real, y::loopSides) = isequal(x, Integer(y))
isequal(x::loopSides, y::Real) = isequal(Integer(x), y)
isless(x::Real, y::loopSides) = isless(x, Integer(y))
isless(x::loopSides, y::Real) = isless(Integer(x), y)
==(x::loopSides, y::Real) = isequal(Integer(x), y)
==(x::Real, y::loopSides) = isequal(x, Integer(y))
convert(::Type{loopSides}, x::Real) = loopSides(Integer(x))
*(x::loopSides, y::Real) = *(Int(x), y)
/(x::loopSides, y::Int64) = /(Int(x), y)

"""
```
abstract type AbstractDlnStr end
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
struct loopDln <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end
```
Idealised dislocation structure types.
"""
abstract type AbstractDlnStr end
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
struct loopDln <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end
"""
```
```
Idealised dislocation loop.
"""
mutable struct DislocationLoop{
    T1<:AbstractArray{<:AbstractDlnSeg,N} where {N},
    T2<:loopSides,
    T3<:Integer,
    T4<:Matrix{<:Integer},
    T5<:AbstractArray{<:Real,N} where {N},
    T6<:Vector{<:nodeType},
}
    segType::T1
    numSides::T2
    nodeSide::T3
    links::T4
    slipPlane::T5
    bVec::T5
    coord::T5
    label::T6
    function DislocationLoop(
        segType,
        numSides,
        nodeSide,
        slipplane,
        bvec,
        label,
    )
        nodeTotal = numSides * nodeSide
        numSegType = length(segType)
        @assert length(label) == nodeTotal
        @assert numSegType == numSides / 2 == size(bvec, 1) ==
                size(slipplane, 1)
        links = zeros(Integer, nodeTotal, 2)
        coord = zeros(nodeTotal, 3)
        seg = zeros(numSegType, 3)
        slipPlane = zeros(0, 3)
        bVec = zeros(0, 3)
        for i = 1:numSegType
            seg[i, :] = makeSegment(segType[i], slipplane[i, :], bvec[i, :])
        end
        if numSides == 4
            # First node.
            coord[1, :] -= seg[1, :] + seg[2, :]
            slipPlane = [slipPlane; slipplane[2, :]']
            bVec = [bVec; bvec[2, :]']
            for k = 1:nodeSide
                coord[1+k, :] = coord[k, :] + seg[1, :]
                slipPlane = [slipPlane; slipplane[1, :]']
                bVec = [bVec; bvec[1, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + nodeSide, :] = coord[k+nodeSide, :] + seg[2, :]
                slipPlane = [slipPlane; slipplane[2, :]']
                bVec = [bVec; bvec[2, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + 2 * nodeSide, :] = coord[k+2*nodeSide, :] -
                                                 seg[1, :]
                slipPlane = [slipPlane; slipplane[1, :]']
                bVec = [bVec; bvec[1, :]']
            end
            # Last node
            coord[2+3*nodeSide, :] = coord[1+3*nodeSide, :] - seg[2, :]
            slipPlane = [slipPlane; slipplane[2, :]']
            bVec = [bVec; bvec[2, :]']
        elseif numSides == 6
            # First node.
            coord[1, :] -= seg[1, :] + seg[2, :] + seg[3, :]
            slipPlane = [slipPlane; slipplane[3, :]']
            bVec = [bVec; bvec[3, :]']
            for k = 1:nodeSide
                coord[1+k, :] = coord[k, :] + seg[1, :]
                slipPlane = [slipPlane; slipplane[1, :]']
                bVec = [bVec; bvec[1, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + nodeSide, :] = coord[k+nodeSide, :] + seg[2, :]
                slipPlane = [slipPlane; slipplane[2, :]']
                bVec = [bVec; bvec[2, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + 2 * nodeSide, :] = coord[k+2*nodeSide, :] +
                                                 seg[3, :]
                slipPlane = [slipPlane; slipplane[3, :]']
                bVec = [bVec; bvec[3, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + 3 * nodeSide, :] = coord[k+3*nodeSide, :] -
                                                 seg[1, :]
                slipPlane = [slipPlane; slipplane[1, :]']
                bVec = [bVec; bvec[1, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + 4 * nodeSide, :] = coord[k+4*nodeSide, :] -
                                                 seg[2, :]
                slipPlane = [slipPlane; slipplane[2, :]']
                bVec = [bVec; bvec[2, :]']
            end
            for k = 1:nodeSide-1
                coord[1 + k + 5 * nodeSide, :] = coord[k+5*nodeSide, :] -
                                                 seg[3, :]
                slipPlane = [slipPlane; slipplane[3, :]']
                bVec = [bVec; bvec[3, :]']
            end
        end
        # Links
        for j = 1:nodeTotal-1
            links[j, :] = [j; 1 + j]
        end
        links[nodeTotal, :] = [nodeTotal; 1]

        new{
            typeof(segType),
            typeof(numSides),
            typeof(nodeSide),
            typeof(links),
            typeof(slipPlane),
            typeof(label),
        }(
            segType,
            numSides,
            nodeSide,
            links,
            slipPlane,
            bVec,
            coord,
            label,
        )


    end
end





mutable struct DislocationNetwork{
    T1<:Matrix{<:Integer},
    T2<:Matrix{<:Real},
    T3<:Matrix{<:Real},
    T4<:Vector{nodeType},
    T5<:Integer,
}
    links::T1 # Links.
    slipPlane::T2 # Slip planes.
    bVec::T2 # Burgers vectors.
    coord::T3 # Node coordinates.
    label::T4 # Node labels.
    numNode::T5 # Number of dislocations.
    numSeg::T5 # Number of segments.
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
        new{
            typeof(links),
            typeof(bVec),
            typeof(coord),
            typeof(label),
            typeof(numNode),
        }(
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

function vcat(network::DislocationNetwork, n::Integer)
    network.links = [network.links; zeros(Integer, n, 2)]
    network.slipPlane = [network.slipPlane; zeros(n, 3)]
    network.bVec = [network.bVec; zeros(n, 3)]
    network.coord = [network.coord; zeros(n, 3)]
    network.label = [network.label; zeros(nodeType, n)]
    return network
end






struct DislocationP{
    T1<:AbstractFloat,
    T2<:Integer,
    T3<:Bool,
    T4<:Union{String,Symbol},
    T5<:Union{Integer,AbstractArray{<:Integer}},
    T6<:Union{AbstractFloat,AbstractArray{<:AbstractFloat}},
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
    # Sources
    numSources::T5
    slipSystems::T5
    distSource::T6
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
        numSources = 0,
        slipSystems = 0,
        distSource = 0.0,
    )
        coreRad == minSegLen == maxSegLen == 0 ? nothing :
        @assert coreRad < minSegLen < maxSegLen
        minArea == maxArea == 0 ? nothing : @assert minArea < maxArea
        @assert length(numSources) == length(slipSystems) == length(distSource)
        new{
            typeof(coreRad),
            typeof(maxConnect),
            typeof(remesh),
            typeof(mobility),
            typeof(numSources),
            typeof(distSource),
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
            numSources,
            slipSystems,
            distSource,
        )
    end # constructor
end # DislocationP
