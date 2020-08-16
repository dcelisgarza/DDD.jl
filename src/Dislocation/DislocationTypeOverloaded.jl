Base.isequal(x::Real, y::nodeType) = isequal(x, Int(y))
Base.isequal(x::nodeType, y::Real) = isequal(Int(x), y)
Base.:(==)(x::nodeType, y::Real) = ==(Int(x), y)
Base.:(==)(x::Real, y::nodeType) = ==(x, Int(y))
Base.isless(x::Real, y::nodeType) = isless(x, Int(y))
Base.isless(x::nodeType, y::Real) = isless(Int(x), y)
Base.convert(::Type{nodeType}, x::Real) = nodeType(Int(x))
Base.zero(::Type{nodeType}) = 0
Base.getindex(x::nodeType, i::Int) = i == 1 ? Int(x) : throw(BoundsError())
Base.iterate(x::nodeType, i = 1) = (length(x) < i ? nothing : (x[i], i + 1))
Base.length(x::nodeType) = 1
Base.size(x::nodeType) = 1

Base.length(::T) where {T <: AbstractDlnSeg} = 1

# Dislocationloop.
Base.length(::DislocationLoop) = 1
Base.getindex(x::DislocationLoop, i::Int) = i == 1 ? x : throw(BoundsError())
Base.eachindex(x::DislocationLoop) = 1
function Base.zero(::Type{DislocationLoop})
    return DislocationLoop(;
        loopType = loopDln(),
        numSides = convert(Int, 0),
        nodeSide = convert(Int, 0),
        numLoops = convert(Int, 0),
        segLen = convert(Float64, 0),
        slipSystem = convert(Int, 0),
        _slipPlane = zeros(3, 0),
        _bVec = zeros(3, 0),
        label = zeros(nodeType, 0),
        buffer = convert(Float64, 0),
        range = zeros(3, 0),
        dist = Zeros(),
    )
end

# DislocationNetwork
function Base.zero(::Type{DislocationNetwork})
    return DislocationNetwork(
        links = zeros(Int, 2, 0),
        slipPlane = zeros(3, 0),
        bVec = zeros(3, 0),
        coord = zeros(3, 0),
        label = zeros(nodeType, 0),
        nodeVel = zeros(3, 0),
        nodeForce = zeros(3, 0),
        numNodeSegConnect = [convert(Int, 0), convert(Int, 0), convert(Int, 0)],
        segForce = zeros(3, 2, 0),
        linksConnect = zeros(Int, 2, 0),
        segIdx = zeros(Int, 0, 3),
    )
end
function Base.push!(network::DislocationNetwork, n::Int)
    network = DislocationNetwork(;
        links = hcat(network.links, zeros(Int, 2, n)),
        slipPlane = hcat(network.slipPlane, zeros(3, n)),
        bVec = hcat(network.bVec, zeros(3, n)),
        coord = hcat(network.coord, zeros(3, n)),
        label = vcat(network.label, zeros(nodeType, n)),
        nodeVel = hcat(network.nodeVel, zeros(3, n)),
        nodeForce = hcat(network.nodeForce, zeros(3, n)),
        numNodeSegConnect = network.numNodeSegConnect,
        connectivity = zeros(Int, 1 + 2 * network.numNodeSegConnect[3], size(network.coord, 2) + n),
        linksConnect = hcat(network.linksConnect, zeros(Int, 2, n)),
        segIdx = vcat(network.segIdx, zeros(Int, n, 3)),
        segForce = cat(network.segForce, zeros(3, 2, n), dims = 3),
    )

    return network
end
function Base.getindex(network::DislocationNetwork, i::Union{Int, AbstractVector{Int}})
    return network.links[:, i],
    network.slipPlane[:, i],
    network.bVec[:, i],
    network.coord[:, i],
    network.label[i],
    network.nodeVel[:, i],
    network.nodeForce[:, i],
    network.connectivity[:, i],
    network.linksConnect[:, i],
    network.segIdx[i, :],
    network.segForce[:, :, i]
end
