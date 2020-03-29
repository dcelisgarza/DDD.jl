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
