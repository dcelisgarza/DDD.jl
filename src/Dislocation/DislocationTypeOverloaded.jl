isequal(x::Real, y::nodeType) = isequal(x, Int(y))
isequal(x::nodeType, y::Real) = isequal(Int(x), y)
isless(x::Real, y::nodeType) = isless(x, Int(y))
isless(x::nodeType, y::Real) = isless(Int(x), y)
==(x::nodeType, y::Real) = isequal(Int(x), y)
==(x::Real, y::nodeType) = isequal(x, Int(y))
convert(::Type{nodeType}, x::Real) = nodeType(Int(x))
zero(::Type{nodeType}) = 0
getindex(x::nodeType, i::Int) = i == 1 ? Int(x) : throw(BoundsError())
iterate(x::nodeType, i = 1) = (length(x) < i ? nothing : (x[i], i + 1))
length(x::nodeType) = 1

"""
```
makeSegment(type::AbstractDlnSeg, slipPlane::AbstractVector{T}, bVec::AbstractVector{T})
    where {T}
```
Make signle segment depending on the segment type, see [`AbstractDlnSeg`](@ref).
"""
@inline function makeSegment(
    type::segEdge,
    slipPlane::AbstractVector{T},
    bVec::AbstractVector{T},
) where {T}
    edge = cross(slipPlane, bVec)
    return edge ./ norm(edge)
end
@inline function makeSegment(
    type::segEdgeN,
    slipPlane::AbstractVector{T},
    bVec::AbstractVector{T},
) where {T}
    return slipPlane ./ norm(slipPlane)
end
@inline function makeSegment(
    type::segScrew,
    slipPlane::AbstractVector{T},
    bVec::AbstractVector{T},
) where {T}
    return bVec ./ norm(bVec)
end
length(::T) where {T <: AbstractDlnSeg} = 1

# Dislocationloop.
length(::DislocationLoop) = 1
getindex(x::DislocationLoop, i::Int) = i == 1 ? x : throw(BoundsError())
eachindex(x::DislocationLoop) = 1
function zero(::Type{DislocationLoop})
    DislocationLoop(;
        loopType = loopDln(),
        numSides = convert(Int, 0),
        nodeSide = convert(Int, 0),
        numLoops = convert(Int, 0),
        segLen = convert(Float64, 0),
        slipSystem = convert(Int, 0),
        _slipPlane = zeros(0, 3),
        _bVec = zeros(0, 3),
        label = zeros(nodeType, 0),
        buffer = convert(Float64, 0),
        range = zeros(0, 3),
        dist = Zeros(),
    )
end

# DislocationNetwork
zero(::Type{DislocationNetwork}) = DislocationNetwork(
    links = zeros(Int, 0, 2),
    slipPlane = zeros(0, 3),
    bVec = zeros(0, 3),
    coord = zeros(0, 3),
    label = zeros(nodeType, 0),
    segForce = zeros(0, 3),
    nodeVel = zeros(0, 3),
    numNode = convert(Int, 0),
    numSeg = convert(Int, 0),
    maxConnect = convert(Int, 0),
)
function push!(network::DislocationNetwork, n::Int)
    network.links = [network.links; zeros(Int, n, 2)]
    network.slipPlane = [network.slipPlane; zeros(n, 3)]
    network.bVec = [network.bVec; zeros(n, 3)]
    network.coord = [network.coord; zeros(n, 3)]
    network.label = [network.label; zeros(nodeType, n)]
    network.segForce = [network.segForce; zeros(n, 3)]
    network.nodeVel = [network.nodeVel; zeros(n, 3)]
    return network
end
function getindex(network::DislocationNetwork, i::Int...)
    return network.links[i, :],
    network.slipPlane[i, :],
    network.bVec[i, :],
    network.coord[i, :],
    network.label[i, :],
    network.segForce[i, :],
    network.nodeVel[i, :],
    network.connectivity[i, :],
    network.linksConnect[i, :],
    network.segIdx[i, :]
end
