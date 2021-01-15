## Dislocation Type Declarations
"""
```
@enum nodeType begin
    none = 0    # Undefined node, value at initialisation.
    intMob = 1  # Internal mobile node.
    intFix = 2  # Internal fixed node.
    srfMob = 3  # Mobile surface node.
    srfFix = 4  # Fixed surface node.
    ext = 5     # External node.
    tmp = 6     # Temporary flag, used during topological operations.
end
```
Different types of nodes behave differently. There are only a finite number of them so an enumerated type provides safety and efficiency. Each value represents a different type of node and therefore its behaviour.

# Meaning
* `none` are uninitialised nodes.
* `intMob` are mobile nodes internal to the convex hull of the domain. They take part in tractions, displacements and dislocation interactions.
* `intFix` are fixed nodes internal to the convex hull of the domain. They participate in the same way as `intMob` nodes except for the fact that their velocities is fixed are zero.
* `srfMob` are mobile nodes that live on the surface of the convex hull of the domain, they are used to track slip steps and therefore participate in the same things as internal nodes but their velocities are restricted to the convex hull surface.
* `srfFix` are fixed surface nodes and have the same characteristics as mobile surface nodes except for having zero velocity.
* `ext` are external nodes that do not participate in dislocation interactions or forces but are used to calculate displacements and track slip steps.
* `tmp` are nodes that are temporarily flagged before they are assigned another type.
"""
@enum nodeType begin
    none = 0
    intMob = 1
    intFix = 2
    srfMob = 3
    srfFix = 4
    virtual = 5
    tmp = 6
end

"""
```
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end    # Undefined segment
struct segEdge <: AbstractDlnSeg end    # Edge segment
struct segEdgeN <: AbstractDlnSeg end   # Edge segment
struct segScrew <: AbstractDlnSeg end   # Screw segment
struct segMixed <: AbstractDlnSeg end   # Mixed segment
```

These types are used to automatically generate segments out of Burgers vectors ``\\bm{b}``, slip planes ``\\bm{n}``, and/or line direction ``\\bm{l}``.

# Meaning
* `segEdge` have ``\\bm{b} ⟂ \\bm{t}`` ,
* `segEdgeN` have ``\\bm{b} ⟂ \\bm{t}`` and ``\\bm{b} ∥ \\bm{n}`` ,
* `segScrew` have ``\\bm{b} ∥ \\bm{t}`` ,
* `segMixed` have none of the above.
"""
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end
struct segEdge <: AbstractDlnSeg end
struct segEdgeN <: AbstractDlnSeg end
struct segScrew <: AbstractDlnSeg end
struct segMixed <: AbstractDlnSeg end

"""
```
abstract type AbstractDlnStr end
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
const loopPure = Union{loopPrism,loopShear}
struct loopMixed <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end
const loopImpure = Union{loopMixed,loopJog,loopKink}
const loopDefined = Union{loopPrism,loopShear,loopMixed,loopJog,loopKink}
struct loopDln <: AbstractDlnStr end
```
These types are used to automatically generate dislocation loops for simulation initialisation.

# Meaning
* `loopPrism` are prismatic loops, their Burgers vectors are perpendicular to the their line direction. They are idealised loops that can be automatically generated as n-gons.
* `loopShear` are shear loops, their line direction goes through edge, screw and line segments as the loop goes round. They are idealised loops that can be automatically generated as n-gons.
* `loopPure` are idealised loops.
* `loopMixed` are loops with prismatic and shear character. They have to be hand-made or require a heuristic to automatically generate.
* `loopDln` is a generic loop used for adding methods to Base functions.
* `loopKink` and `loopJog` are structures formed by colliding dislocations. They are not currently used.
* `loopImpure` are non-idealised loops.
* `loopDefined` are defined loop types.
"""
abstract type AbstractDlnStr end
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
const loopPure = Union{loopPrism,loopShear}
struct loopMixed <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end
const loopImpure = Union{loopMixed,loopJog,loopKink}
const loopDefined = Union{loopPrism,loopShear,loopMixed,loopJog,loopKink}
struct loopDln <: AbstractDlnStr end

"""
```
abstract type AbstractDistribution end
struct Zeros <: AbstractDistribution end
struct Rand <: AbstractDistribution end
struct Randn <: AbstractDistribution end
struct Regular <: AbstractDistribution end
```
Spatial distributions for dislocation sources. These are used to automatically generate networks with a given distribution.

# Meaning
* `Zeros` makes the network generation functions place the center of the generated dislocation loops at the origin. This can be used to generate a network and loops can be manually or pseudo-manually distributed in the domain.
* `Rand` makes the network generation functions uniformly distribute the dislocations according to the range and buffer values in the dislocation loop structure.
* `Rand` makes the network generation functions normally distribute the dislocations according to the range and buffer values in the dislocation loop structure.
* `Rand` TBA, will regularly distribute dislocations according to the range, buffer and other args given to the dislocation network generator.
"""
abstract type AbstractDistribution end
struct Zeros <: AbstractDistribution end
struct Rand <: AbstractDistribution end
struct Randn <: AbstractDistribution end
struct Regular <: AbstractDistribution end

"""
```
abstract type AbstractMobility end
struct mobBCC <: AbstractMobility end
struct mobFCC <: AbstractMobility end
struct mobHCP <: AbstractMobility end
```
Types to dispatch different mobility functions.

# Meaning
* `mobBCC` is used to dispatch the default BCC mobility function.
* `mobFCC` is used to dispatch the default FCC mobility function.
* `mobHCP` is used to dispatch the default HCP mobility function.
"""
abstract type AbstractMobility end
struct mobBCC <: AbstractMobility end
struct mobFCC <: AbstractMobility end
struct mobHCP <: AbstractMobility end

"""
```
struct SlipSystem{T1, T2}
    crystalStruct::T1   # Crystal structure
    slipPlane::T2       # Slip plane
    bVec::T2            # Burgers vector
end
```
Stores slip systems. 
"""
struct SlipSystem{T1,T2}
    crystalStruct::T1
    slipPlane::T2
    bVec::T2
end

"""
```
struct DislocationParameters{T1,T2,T3,T4}
    coreRad::T1         # Dislocation core radius.
    coreRadSq::T1       # Square of the dislocation core radius.
    coreRadMag::T1      # Magnitude of the core radius (real units for post-processing).
    minSegLen::T1       # Minimum segment length.
    maxSegLen::T1       # Maximum segment length.
    twoMinSegLen::T1    # Twice the minimum segment length.
    minArea::T1         # Minimum area for remeshing.
    maxArea::T1         # Maximum area for remeshing.
    minAreaSq::T1       # Square of the minimum area.
    maxAreaSq::T1       # Square of the maximum area.
    edgeDrag::T1        # Edge drag coefficient.
    screwDrag::T1       # Screw drag coefficient.
    climbDrag::T1       # Climb drag coefficient.
    lineDrag::T1        # Line drag coefficient.
    maxConnect::T2      # Maximum connectivity of nodes.
    mobility::T3        # Dislocation mobility.
    remesh::T4          # Remesh flag.
    collision::T4       # Collision flag.
    separation::T4      # Separation flag.
    virtualRemesh::T4   # Virtual remeshing flag.
    parCPU::T4          # Parallelise on CPU flag.
    parGPU::T4          # Parallelise on GPU flag.
    slipStepCritLen::T1 # Critical length for slip step tracking.
    slipStepCritArea::T1    # Critical area for slip step tracking.
end
```
Stores the dislocation parameters.
"""
struct DislocationParameters{T1,T2,T3,T4}
    coreRad::T1
    coreRadSq::T1
    coreRadMag::T1
    minSegLen::T1
    maxSegLen::T1
    twoMinSegLen::T1
    minArea::T1
    maxArea::T1
    minAreaSq::T1
    maxAreaSq::T1
    edgeDrag::T1
    screwDrag::T1
    climbDrag::T1
    lineDrag::T1
    maxConnect::T2
    mobility::T3
    remesh::T4
    collision::T4
    separation::T4
    virtualRemesh::T4
    parCPU::T4
    parGPU::T4
    slipStepCritLen::T1
    slipStepCritArea::T1
end

"""
```
struct DislocationLoop{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    loopType::T1    # Loop type.
    numSides::T2    # Number of sides in the loop.
    nodeSide::T2    # Nodes per side of the loop.
    numLoops::T2    # Number of loops to generate when making the network.
    segLen::T3      # Segment lengths.
    slipSystem::T4  # Slip system.
    links::T5       # Links.
    slipPlane::T6   # Slip planes.
    bVec::T6        # Burgers vectors.
    coord::T6       # Coordinates.
    label::T7       # Node labels.
    buffer::T8      # Buffer for distributions.
    range::T9       # Range for distributions.
    dist::T10       # Distribution.
end
```
Stores dislocation loops and parameters used to generate a [`DislocationNetwork`](@ref).
"""
struct DislocationLoop{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
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
    range::T9
    dist::T10
end

"""
```
struct DislocationNetwork{T1,T2,T3,T4,T5,T6}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeVel::T2
    nodeForce::T2
    numNode::T4
    numSeg::T4
    maxConnect::T5
    connectivity::T1
    linksConnect::T1
    segIdx::T1
    segForce::T6
end
```
Stores the dislocation network generated from [`DislocationLoop`](@ref).

The internal constructor lets users manually define networks,
```
DislocationNetwork(
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNode::T4 = zeros(Int, 1),
    numSeg::T4 = zeros(Int, 1),
    maxConnect::T5 = 4,
    connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::T1 = zeros(Int, 2, size(links, 2)),
    segIdx::T1 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links)...),
) where {
    T1 <: AbstractArray{T,N} where {T,N},
    T2 <: AbstractArray{T,N} where {T,N},
    T3 <: AbstractVector{nodeType},
    T4 <: Union{Int,AbstractVector{Int}},
    T5 <: Int,T6 <: AbstractArray{T,N} where {T,N},
}
```
provides some default values and validates inputs.

Its keyword constructor just calls the structure.
```
DislocationNetwork(;
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNode::T4 = zeros(Int, 1),
    numSeg::T4 = zeros(Int, 1),
    maxConnect::T5 = 0,
    connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::T1 = zeros(Int, 2, size(links, 2)),
    segIdx::T1 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links, 2), 0),
) where {
    T1 <: AbstractArray{T,N} where {T,N},
    T2 <: AbstractArray{T,N} where {T,N},
    T3 <: AbstractVector{nodeType},T4 <: AbstractVector{Int},
    T5 <: Int,
    T6 <: AbstractArray{T,N} where {T,N},
}
```

However, the recommended way of creating a network is to use a variable of type [`DislocationLoop`](@ref) to automatically generate the desired network.
```
DislocationNetwork(
    sources::T1,
    maxConnect::T2 = 4,
    args...;
    memBuffer = nothing,
    checkConsistency::T3 = true,
    kw...,
) where {
    T1 <: Union{T,
                AbstractVector{T},
                NTuple{N,T} where N
                } where {T <: DislocationLoop},
    T2 <: Int,
    T3 <: Bool,
}
```

There is also a mutating version of this constructor that expands the network accordingly.
```
DislocationNetwork!(
    network::T1,
    sources::T2,
    maxConnect::T3 = 4,
    args...;
    checkConsistency::T4 = true,
    kw...,
) where {
    T1 <: DislocationNetwork,
    T2 <: Union{T,
                AbstractVector{T},
                NTuple{N,T} where N
                } where {T <: DislocationLoop},
    T3 <: Int,
    T4 <: Bool,
}
```

## Argument Explanation

- `args...` are optional arguments that will be passed on to the `loopDistribution` function which distributes the loops in `sources` according to the type of their `dist` variable.
- `kw...` are optional keyword arguments that will also be passed to `loopDistribution`.
- `memBuffer` is the numerical value for allocating memory in advance, the quantity ``\\textrm{memBuffer} \\times N`` where `N` is the total number of nodes in `sources`, will be the initial number of entries allocated in the matrices that keep the network's data, if it is `nothing` then the number of entries is ``\\textrm{round}(N \\log_{2}(N))``.
"""
struct DislocationNetwork{T1,T2,T3,T4,T5,T6}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeVel::T2
    nodeForce::T2
    numNode::T4
    numSeg::T4
    maxConnect::T5
    connectivity::T1
    linksConnect::T1
    segIdx::T1
    segForce::T6

    function DislocationNetwork(
        links::T1,
        slipPlane::T2,
        bVec::T2,
        coord::T2,
        label::T3,
        nodeVel::T2,
        nodeForce::T2,
        numNode::T4 = zeros(Int, 1),
        numSeg::T4 = zeros(Int, 1),
        maxConnect::T5 = 4,
        connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
        linksConnect::T1 = zeros(Int, 2, size(links, 2)),
        segIdx::T1 = zeros(Int, size(links, 2), 3),
        segForce::T6 = zeros(3, size(links)...),
    ) where {T1 <: AbstractArray{T,N} where {T,N},T2 <: AbstractArray{T,N} where {T,N},T3 <: AbstractVector{nodeType},T4 <: Union{Int,AbstractVector{Int}},T5 <: Int,T6 <: AbstractArray{T,N} where {T,N},}

        @assert size(links, 1) == size(segForce, 2) == 2
        @assert size(bVec, 1) == size(slipPlane, 1) == size(coord, 1) size(segForce, 1) == 3
        @assert size(links, 2) == size(bVec, 2) == size(slipPlane, 2) == size(segForce, 3)
        @assert size(coord, 2) == length(label)
        @assert length(numNode) == length(numSeg) == 1

        typeof(numNode) <: AbstractVector ? numNodeArr = numNode : numNodeArr = [numNode]
        typeof(numSeg) <: AbstractVector ? numSegArr = numSeg : numSegArr = [numSeg]

        return new{T1,T2,T3,typeof(numNodeArr),T5,T6}(
            links,
            slipPlane,
            bVec,
            coord,
            label,
            nodeVel,
            nodeForce,
            numNodeArr,
            numSegArr,
            maxConnect,
            connectivity,
            linksConnect,
            segIdx,
            segForce,
)
end
end
function DislocationNetwork(;
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNode::T4 = zeros(Int, 1),
    numSeg::T4 = zeros(Int, 1),
    maxConnect::T5 = 0,
    connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::T1 = zeros(Int, 2, size(links, 2)),
    segIdx::T1 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links, 2), 0),
) where {T1 <: AbstractArray{T,N} where {T,N},T2 <: AbstractArray{T,N} where {T,N},T3 <: AbstractVector{nodeType},T4 <: AbstractVector{Int},T5 <: Int,T6 <: AbstractArray{T,N} where {T,N},}

    return DislocationNetwork(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        nodeVel,
        nodeForce,
        numNode,
        numSeg,
        maxConnect,
        connectivity,
        linksConnect,
        segIdx,
segForce,
)
end
function DislocationNetwork(
    sources::T1,
    maxConnect::T2 = 4,
    args...;
    memBuffer = nothing,
    checkConsistency::T3 = true,
    kw...,
) where {T1 <: Union{T,AbstractVector{T},NTuple{N,T} where N} where {T <: DislocationLoop},T2 <: Int,T3 <: Bool,}

    # Initialisation.
    nodeTotal::Int = 0
    lims = zeros(MMatrix{3,2})
    # Calculate node total.
    for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end
    # Memory buffer.
    isnothing(memBuffer) ? nodeBuffer = Int(round(nodeTotal * log2(nodeTotal))) :
    nodeBuffer = nodeTotal * Int(memBuffer)

    # Allocate memory.
    links = zeros(Int, 2, nodeBuffer)
    slipPlane = zeros(3, nodeBuffer)
    bVec = zeros(3, nodeBuffer)
    coord = zeros(3, nodeBuffer)
    label = zeros(nodeType, nodeBuffer)
    nodeVel = zeros(Float64, 3, nodeBuffer)
    nodeForce = zeros(Float64, 3, nodeBuffer)
    numNode = nodeTotal
    numSeg = nodeTotal
    segForce = zeros(Float64, 3, 2, nodeBuffer)

    initIdx = 1
    makeNetwork!(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        sources,
        lims,
        initIdx,
        args...;
        kw...,
    )

    # Calculate number of segments and indexing matrix.
    numSeg, segIdx = getSegmentIdx(links, label)
    # Generate connectivity and linksConnect matrix.
    connectivity, linksConnect = makeConnect(links, maxConnect)

    # Create network.
    network = DislocationNetwork(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        nodeVel,
        nodeForce,
        [numNode],
        [numSeg],
        maxConnect,
        connectivity,
        linksConnect,
        segIdx,
        segForce,
    )

    # Check that the network is generated properly.
    checkConsistency ? checkNetwork(network) : nothing

    return network
end
function DislocationNetwork!(
    network::T1,
    sources::T2,
    maxConnect::T3 = 4,
    args...;
    checkConsistency::T4 = true,
    kw...,
) where {T1 <: DislocationNetwork,T2 <: Union{T,AbstractVector{T},NTuple{N,T} where N} where {T <: DislocationLoop},T3 <: Int,T4 <: Bool,}
    # For comments see DislocationNetwork. It is a 1-to-1 translation except that this one modifies the network in-place.
    
    iszero(network) && return DislocationNetwork(
        sources,
        maxConnect = maxConnect,
        args...;
        checkConsistency = checkConsistency,
        kw...,
    )
    
    @assert network.maxConnect == maxConnect "Maximum connectivity of added network must be equal to that of the existing network."
    
    nodeTotal::Int = 0
    lims = zeros(MMatrix{3,2})
    for i in eachindex(sources)
        nodeTotal += sources[i].numLoops * length(sources[i].label)
    end
    numNode = nodeTotal

    # Allocate memory.
    available = length(findall(x -> x == 0, network.label))
    if nodeTotal > available
        newEntries = Int(round(nodeTotal * log2(nodeTotal)))
        network = push!(network, newEntries)
    end
    
    links = network.links
    slipPlane = network.slipPlane
    bVec = network.bVec
    coord = network.coord
    label = network.label
    
    initIdx::Int = 1
    first = findfirst(x -> x == 0, label)
    isnothing(first) ? initIdx = 1 : initIdx = first
    makeNetwork!(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        sources,
        lims,
        initIdx,
        args...;
        kw...,
    )
    network.numNode[1] += numNode
    
    getSegmentIdx!(network)
    makeConnect!(network)
    
    checkConsistency ? checkNetwork(network) : nothing
    return network
end

function makeNetwork!(
    links,
    slipPlane,
    bVec,
    coord,
    label,
    sources,
    lims,
    initIdx,
    args...;
    kw...,
)
    nodeTotal::Int = 0
    elemT = eltype(coord)
    @inbounds for i in eachindex(sources)
        idx = initIdx + nodeTotal
        nodesLoop = length(sources[i].label)
        numLoops = sources[i].numLoops
        numNodes = numLoops * nodesLoop
        # Calculate the normalised displacements for all loops in sources[i] according to their distribution.
        disp = loopDistribution(sources[i].dist, numLoops, args...; kw...)
        # Calculate the real spatial limits of the distributions.
        limits!(lims, mean(sources[i].segLen), sources[i].range, sources[i].buffer)
        # Fill out the data for all loops specified in sources[i].
        for j in 1:numLoops
            # The number of nodes in the loop is nodesLoop, so that's our stride inside sources[i]
            idxi = idx + (j - 1) * nodesLoop
            idxf = idxi + nodesLoop - 1
            # Links are numbered sequentially in network so we have to account for previously assigned links.
            links[:, idxi:idxf] .=
                sources[i].links[:, 1:nodesLoop] .+ (nodeTotal + initIdx - 1)
            slipPlane[:, idxi:idxf] = sources[i].slipPlane[:, 1:nodesLoop]
            bVec[:, idxi:idxf] = sources[i].bVec[:, 1:nodesLoop]
            coord[:, idxi:idxf] = sources[i].coord[:, 1:nodesLoop]
            label[idxi:idxf] = sources[i].label[1:nodesLoop]
            # Map the normalised displacements to real space using the real limits and translate the nodes' coordinates accordingly.
            staticDisp = SVector{3,elemT}(disp[1, j], disp[2, j], disp[3, j])
            viewCoord = @view coord[:, idxi:idxf]
            translatePoints!(viewCoord, lims, staticDisp)
            nodeTotal += nodesLoop
        end
    end
    return nothing
end