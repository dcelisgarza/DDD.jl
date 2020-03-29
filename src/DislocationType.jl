"""
```
DislocationLoop{
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{
        T where {T <: AbstractDlnSeg},
        AbstractArray{<:AbstractDlnSeg, N} where {N},
    },
    T4 <: Union{
        T where {T <: Float64},
        AbstractArray{<:Float64, N} where {N}
    },
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
    slipSystem::T5  # Slip system/systems of segments.
    links::T6       # Link matrix for dislocation nodes.
    slipPlane::T7   # Slip planes of all segments in loop.
    bVec::T7        # Burgers vector of all segments in loop.
    coord::T7       # Coords of all nodes in loop.
    label::T8
    buffer::T9
    range::T7
    dist::T10
```
Dislocation loop structure generated via the constructor [`makeLoop`](@ref), see its documentation for further details.
"""
struct DislocationLoop{
    T1 <: AbstractDlnStr,
    T2 <: Int64,
    T3 <: Union{
        T where {T <: AbstractDlnSeg},
        AbstractArray{<:AbstractDlnSeg, N} where {N},
    },
    T4 <: Union{T where {T <: Float64}, AbstractArray{<:Float64, N} where {N}},
    T5 <: Union{T where {T <: Int64}, AbstractArray{<:Int64, N} where {N}},
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
# Zero constructor.
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
    loopType::T1,   # Type of loop.
    numSides::T2,   # Number of sides in loop, must be even.
    nodeSide::T2,   # Nodes per side.
    numLoops::T2,   # Number of loops in initial network.
    segType::T3,    # Type of segments in loop (future-proofing).
    segLen::T4,     # Lengths of half the segments in loop, why `numSides` must be even.
    slipSystem::T2, # Index of the loop's slip system.
    _slipPlane::T5, # Slip plane vector.
    _bVec::T5,      # Burgers vector.
    label::T6,      # Labels, length must be equal to `numSides * nodeSide`.
    buffer::T7,     # Buffer distance to prevent overlap.
    range::T8,      # Spatial range of `dist`, defines the space occupied by the loops in the network.
    dist::T9,       # Spatial distribution of loops in `range`.
) where {
    T1 <: AbstractDlnStr, # Dislocation structure.
    T2 <: Int64,
    T3 <: AbstractDlnSeg, # Dislocation segment type.
    T4 <: Union{
        T where {T <: Float64},
        AbstractArray{<:Float64, N} where {N}
    },
    T5 <: AbstractArray{<:Float64, N} where {N},
    T6 <: Vector{nodeType}, # Dislocation node enumerated type.
    T7 <: Float64,
    T8 <: AbstractArray{<:Float64, N} where {N},
    T9 <: AbstractDistribution, # Spatial distribution type, defines the loops' spatial distribution in the network.
}
```
Constructor function for [`DislocationLoop`](@ref). See [`AbstractDlnStr`](@ref), [`AbstractDlnSeg`](@ref), [`nodeType`](@ref), [`AbstractDistribution`](@ref) for further details.
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
    T3 <: AbstractDlnSeg,
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
# Overloaded DislocationLoop functions.
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
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    numNode::T4 = 0     # Total number of nodes in network.
    numSeg::T4 = 0      # Total number of segs in network.
    maxConnect::T5 = 4  # Maximum connectivity of nodes.
    connectivity::T1
    linksConnect::T1
```
Dislocation Network structure. See [`DislocationLoop`](@ref), [`makeNetwork`](@ref), [`makeNetwork!`](@ref) for further details.
"""
mutable struct DislocationNetwork{
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
    numNode::T4
    numSeg::T4
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
        maxConnect = 0,
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

"""
```
makeConnect(
    links::AbstractArray{Int64, N},
    label::Vector{nodeType},
    maxConnect::Integer,
) where {N}
```
Returns tuple of matrices `connectivity` and `linksConnect`. `connectivity` contains the number of other other nodes each node is connected to, up to `maxConnect` other nodes. It also contains the segments in which it's involved. `linksConnect` is the running total of the number of other nodes each node is connected to, will probably be deleted in the future. This is called from [`makeNetwork`](@ref) and [`makeNetwork!`](@ref).
"""
function makeConnect(
    links::AbstractArray{Int64, N},
    label::Vector{nodeType},
    maxConnect::Integer,
) where {N}

    iLnk = findall(x -> x != 0, links[:, 1]) # Indices of defined links.
    lenLabel = length(label)
    connectivity = zeros(Int64, lenLabel, 1 + 2 * maxConnect)
    linksConnect = zeros(Int64, size(links, 1), 2)
    @inbounds for i in eachindex(iLnk)
        idx = iLnk[i] # Index of links that contain non zero links.
        # links[idx, :] yields the nodes involved in the link
        n1 = links[idx, 1] # Node 1, it is the row of the coord matrix
        n2 = links[idx, 2] # Node 2

        connectivity[n1, 1] += 1 # Number of other nodes node n1 is connected to.
        connectivity[n2, 1] += 1 # Number of other nodes node n2 is connected to.

        tmp1 = 2 * connectivity[n1, 1]
        tmp2 = 2 * connectivity[n2, 1]

        connectivity[n1, tmp1:(tmp1 + 1)] = [idx, 1] # idx = linkID, 1 = first node in link with linkID
        connectivity[n2, tmp2:(tmp2 + 1)] = [idx, 2] # idx = linkID, 2 = second node in link with linkID

        linksConnect[idx, 1] = connectivity[n1, 1]
        linksConnect[idx, 2] = connectivity[n2, 1]
    end
    return connectivity, linksConnect
end

"""
```
makeNetwork(
    sources::Union{
                DislocationLoop,
                AbstractVector{<:DislocationLoop}
            }, # Dislocation structures.
    maxConnect::Integer = 4,
    memBuffer::Integer = 10, # Buffer for memory allocation. The code will allocate the total number of nodes times `memBuffer` to reduce dynamic memory allocation during runtime.
)
```
Constructor for [`DislocationNetwork`](@ref), see for in-place version [`makeNetwork!`](@ref).
"""
function makeNetwork(
    sources::Union{DislocationLoop, AbstractVector{<:DislocationLoop}},
    maxConnect::Integer = 4,
    memBuffer::Integer = 10,
)
    nodeTotal::Integer = 0
    lims = zeros(Float64, 2, 3)
    # Allocate memory.
    @inbounds for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end

    nodeBuffer::Integer = nodeTotal * memBuffer

    network = DislocationNetwork(
        zeros(Int64, nodeBuffer, 2),
        zeros(nodeBuffer, 3),
        zeros(nodeBuffer, 3),
        zeros(nodeBuffer, 3),
        zeros(nodeType, nodeBuffer),
        nodeTotal,
        nodeTotal,
    )

    nodeTotal = 0
    initIdx = findfirst(x -> x == -1, network.label)
    initIdx == nothing ? initIdx = 0 : nothing
    @inbounds for i in eachindex(sources)
        # Indices.
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        disp = loopDistribution(sources[i].dist, numLoops)
        limits!(
            lims,
            mean(sources[i].segLen),
            sources[i].range,
            sources[i].buffer,
        )
        for j = 1:numLoops
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            # Prepare to distribute sources.
            network.links[idxi:idxf, :] =
                sources[i].links[1:nodesLoop, :] .+ (nodeTotal + initIdx - 1)
            network.slipPlane[idxi:idxf, :] .=
                sources[i].slipPlane[1:nodesLoop, :]
            network.bVec[idxi:idxf, :] .= sources[i].bVec[1:nodesLoop, :]
            network.coord[idxi:idxf, :] .= sources[i].coord[1:nodesLoop, :]
            # Move loop.
            network.coord[idxi:idxf, :] .= translatePoints(
                sources[i].coord[1:nodesLoop, :],
                lims,
                disp[j, :],
            )
            network.label[idxi:idxf, :] .= sources[i].label[1:nodesLoop, :]
            nodeTotal += nodesLoop
        end
    end

    network.connectivity, network.linksConnect =
        makeConnect(network.links, network.label, maxConnect)
    return network
end
"""
```
makeNetwork!(
    network::DislocationNetwork,
    sources::Union{
                DislocationLoop,
                AbstractVector{<:DislocationLoop}
            },
    maxConnect::Integer = 4,
)
```
In-place constructor for [`DislocationNetwork`](@ref), see [`makeNetwork`](@ref) for constructor.
"""
function makeNetwork!(
    network::DislocationNetwork,
    sources::Union{DislocationLoop, AbstractVector{<:DislocationLoop}},
    maxConnect::Integer = 4,
)
    nodeTotal::Integer = 0
    lims = zeros(Float64, 2, 3)
    # Allocate memory.
    @inbounds for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end
    available = findfirst(x -> x == -1, network.label)
    if available == nothing
        push!(network, nodeTotal)
    else
        check = length(network.label) - (available - 1) - nodeTotal
        check < 0 ? push!(network, -check) : nothing
    end

    nodeTotal = 0
    initIdx = findfirst(x -> x == -1, network.label)
    initIdx == nothing ? initIdx = 0 : nothing
    @inbounds for i in eachindex(sources)
        # Indices.
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        disp = loopDistribution(sources[i].dist, numLoops)
        limits!(
            lims,
            mean(sources[i].segLen),
            sources[i].range,
            sources[i].buffer,
        )
        for j = 1:numLoops
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            # Prepare to distribute sources.
            network.links[idxi:idxf, :] =
                sources[i].links[1:nodesLoop, :] .+ (nodeTotal + initIdx - 1)
            network.slipPlane[idxi:idxf, :] .=
                sources[i].slipPlane[1:nodesLoop, :]
            network.bVec[idxi:idxf, :] .= sources[i].bVec[1:nodesLoop, :]
            network.coord[idxi:idxf, :] .= sources[i].coord[1:nodesLoop, :]
            # Move loop.
            network.coord[idxi:idxf, :] .= translatePoints(
                sources[i].coord[1:nodesLoop, :],
                lims,
                disp[j, :],
            )
            network.label[idxi:idxf, :] .= sources[i].label[1:nodesLoop, :]
            nodeTotal += nodesLoop
        end
    end

    network.numNode += nodeTotal
    network.numSeg += nodeTotal

    network.connectivity, network.linksConnect =
        makeConnect(network.links, network.label, maxConnect)
    return network
end
# Overloaded functions.
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
