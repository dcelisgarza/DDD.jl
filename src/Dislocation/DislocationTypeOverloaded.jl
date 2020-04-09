isequal(x::Real, y::nodeType) = isequal(x, Integer(y))
isequal(x::nodeType, y::Real) = isequal(Integer(x), y)
isless(x::Real, y::nodeType) = isless(x, Integer(y))
isless(x::nodeType, y::Real) = isless(Integer(x), y)
==(x::nodeType, y::Real) = isequal(Integer(x), y)
==(x::Real, y::nodeType) = isequal(x, Integer(y))
convert(::Type{nodeType}, x::Real) = nodeType(Integer(x))
zero(::Type{nodeType}) = -1
getindex(x::nodeType, i::Integer) = i == 1 ? Int(x) : throw(BoundsError())
iterate(x::nodeType, i=1) = (length(x) < i ? nothing : (x[i], i + 1))
length(x::nodeType) = 1

"""
```
makeSegment(type::AbstractDlnSeg, slipPlane::Vector{T}, bVec::Vector{T})
    where {T<:Float64}
```
Make signle segment depending on the segment type, see [`AbstractDlnSeg`](@ref).
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
length(::T) where {T <: AbstractDlnSeg} = 1

# Dislocationloop.
length(::DislocationLoop) = 1
getindex(x::DislocationLoop, i::Integer) = i == 1 ? x : throw(BoundsError())
eachindex(x::DislocationLoop) = 1
function zero(::Type{DislocationLoop})
    DislocationLoop(
        loopType = loopDln(),
        numSides = convert(Int64, 0),
        nodeSide = convert(Int64, 0),
        numLoops = convert(Int64, 0),
        segType = segNone(),
        segLen = convert(Float64, 0),
        slipSystem = convert(Int64, 0),
        _slipPlane = zeros(Float64, 0, 3),
        _bVec = zeros(Float64, 0, 3),
        label = zeros(nodeType, 0),
        buffer = convert(Float64, 0),
        range = zeros(0, 3),
        dist = Zeros(),
    )
end

# DislocationNetwork
zero(::Type{DislocationNetwork}) = DislocationNetwork(
    links = zeros(Int64, 0, 2),
    slipPlane = zeros(Float64, 0, 3),
    bVec = zeros(Float64, 0, 3),
    coord = zeros(Float64, 0, 3),
    label = zeros(nodeType, 0),
    numNode = convert(Int64, 0),
    numSeg = convert(Int64, 0),
    maxConnect = convert(Int64, 0),
)
function push!(network::DislocationNetwork, n::Int64)
    network.links = [network.links; zeros(Int64, n, 2)]
    network.slipPlane = [network.slipPlane; zeros(n, 3)]
    network.bVec = [network.bVec; zeros(n, 3)]
    network.coord = [network.coord; zeros(n, 3)]
    network.label = [network.label; zeros(nodeType, n)]
    return network
end
