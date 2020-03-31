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

# DislocationNetwork
zero(::Type{DislocationNetwork}) = DislocationNetwork(
    zeros(Int64, 0, 2),
    zeros(Float64, 0, 3),
    zeros(Float64, 0, 3),
    zeros(Float64, 0, 3),
    zeros(nodeType, 0),
    convert(Int64, 0),
    convert(Int64, 0),
)
function push!(network::DislocationNetwork, n::Int64)
    network.links = [network.links; zeros(Int64, n, 2)]
    network.slipPlane = [network.slipPlane; zeros(n, 3)]
    network.bVec = [network.bVec; zeros(n, 3)]
    network.coord = [network.coord; zeros(n, 3)]
    network.label = [network.label; zeros(nodeType, n)]
    return network
end
