Base.:(==)(x::nodeTypeDln, y::Real) = ==(Int(x), y)
Base.:(==)(x::Real, y::nodeTypeDln) = ==(x, Int(y))
Base.convert(::Type{nodeTypeDln}, x::Real) = nodeTypeDln(Int(x))
Base.zero(::Type{nodeTypeDln}) = nodeTypeDln(0)

Base.getindex(x::DislocationLoop, i) = i == 1 ? x : throw(BoundsError())
Base.eachindex(x::DislocationLoop) = 1

function Base.zero(::Type{DislocationLoop})
    return DislocationLoop(;
        loopType = loopDln(),
        numSides = convert(Int, 0),
        nodeSide = convert(Int, 0),
        numLoops = convert(Int, 0),
        segLen = convert(Float64, 0),
        slipSystemIdx = convert(Int, 0),
        slipSystem = SlipSystem(nothing, zeros(3), zeros(3)),
        label = zeros(nodeTypeDln, 0),
        buffer = convert(Float64, 0),
        range = zeros(3, 0),
        dist = Zeros(),
    )
end

function Base.zero(::Type{DislocationNetwork})
    return DislocationNetwork(;
        links = zeros(Int, 2, 0),
        slipPlane = zeros(3, 0),
        bVec = zeros(3, 0),
        coord = zeros(3, 0),
        label = zeros(nodeTypeDln, 0),
        nodeVel = zeros(3, 0),
        nodeForce = zeros(3, 0),
        numNode = [0],
        numSeg = [0],
        maxConnect = 0,
        segForce = zeros(3, 2, 0),
        linksConnect = zeros(Int, 2, 0),
        segIdx = zeros(Int, 0, 3),
    )
end
function Base.push!(network::DislocationNetwork, n)
    elemT = eltype(network.coord)
    network = DislocationNetwork(;
        links = hcat(network.links, zeros(Int, 2, n)),
        slipPlane = hcat(network.slipPlane, zeros(elemT, 3, n)),
        bVec = hcat(network.bVec, zeros(elemT, 3, n)),
        coord = hcat(network.coord, zeros(elemT, 3, n)),
        label = vcat(network.label, zeros(nodeTypeDln, n)),
        nodeVel = hcat(network.nodeVel, zeros(elemT, 3, n)),
        nodeForce = hcat(network.nodeForce, zeros(elemT, 3, n)),
        numNode = copy(network.numNode),
        numSeg = copy(network.numSeg),
        maxConnect = copy(network.maxConnect),
        connectivity = hcat(
            network.connectivity,
            zeros(Int, size(network.connectivity, 1), n),
        ),
        linksConnect = hcat(network.linksConnect, zeros(Int, 2, n)),
        segIdx = vcat(network.segIdx, zeros(Int, n, 3)),
        segForce = cat(network.segForce, zeros(elemT, 3, 2, n), dims = 3),
    )

    return network
end
function Base.getindex(network::DislocationNetwork, i)
    return @views network.links[:, i],
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
function Base.iszero(network::DislocationNetwork)
    (network.numNode[1] == 0 || network.numSeg[1] == 0) && return true
    return false
end

"""
```
loopDistribution(::Rand, n, args...; kw...)
```
Returns a `3 × n` matrix whose points follow a random uniform distribution. Used by [`DislocationNetwork`](@ref) when the `dist` parameter of [`DislocationLoop`](@ref) is `Rand()`.
"""
loopDistribution(::Rand, n, args...; kw...) = rand(3, n)
"""
```
loopDistribution(::Zeros, n, args...; kw...)
```
Returns a `3 × n` zeros matrix. Used by [`DislocationNetwork`](@ref) when the `dist` parameter of [`DislocationLoop`](@ref) is `Zeros()`.
"""
loopDistribution(::Zeros, n, args...; kw...) = zeros(3, n)
"""
```
loopDistribution(::Rand, n, args...; kw...)
```
Returns a `3 × n` matrix whose points follow a random normal distribution. Used by [`DislocationNetwork`](@ref) when the `dist` parameter of [`DislocationLoop`](@ref) is `Randn()`.
"""
loopDistribution(::Randn, n, args...; kw...) = randn(3, n)
"""
```
loopDistribution(::Rand, n, args...; kw...)
```
Returns a `3 × n` matrix whose points follow are regularly placed. Used by [`DislocationNetwork`](@ref) when the `dist` parameter of [`DislocationLoop`](@ref) is `Regular()`.

!!! note 
    Not yet implemented.
"""
function loopDistribution(::Regular, n, args...; kw...)
    return error("loopDistribution: regular distribution yet not implemented")
end

"""
```
limits!(lims, segLen, range, buffer)
```
Compute the bounding limits within which a [`DislocationLoop`](@ref) will be distributed in a [`DislocationNetwork`](@ref).
"""
function limits!(lims, segLen, range, buffer)
    @inbounds for i in 1:size(lims, 2)
        @simd for j in 1:size(lims, 1)
            lims[j, i] = range[j, i] + buffer * segLen
        end
    end
    return lims
end

"""
```
translatePoints!(coord, lims, disp)
```
Translates coordinates using the limits and displacements calculated by [`limits!`](@ref) and [`loopDistribution`](@ref).
"""
function translatePoints!(coord, lims, disp)
    @inbounds for i in 1:size(coord, 2)
        @simd for j in 1:size(coord, 1)
            coord[j, i] += lims[j, 1] + (lims[j, 2] - lims[j, 1]) * disp[j]
        end
    end
    return coord
end

"""
```
makeSegment(::segEdge, slipPlane, bVec)
```
Return the edge segment `(slipPlane × bVec) / || slipPlane × bVec ||`. Used to create segments for [`DislocationNetwork`](@ref).
"""
function makeSegment(::segEdge, slipPlane, bVec)
    edge = cross(slipPlane, bVec)
    return edge / norm(edge)
end
"""
```
makeSegment(::segEdge, slipPlane, bVec)
```
Return the edge segment `slipPlane / || slipPlane ||`. Used to create segments for [`DislocationNetwork`](@ref).
"""
function makeSegment(::segEdgeN, slipPlane, bVec)
    return slipPlane / norm(slipPlane)
end
"""
```
makeSegment(::segScrew, slipPlane, bVec)
```
Return the screw segment `bVec / || bVec ||`. Used to create segments for [`DislocationNetwork`](@ref).
"""
function makeSegment(::segScrew, slipPlane, bVec)
    return bVec / norm(bVec)
end

"""
```
makeConnect(links, maxConnect)
makeConnect!(network::DislocationNetwork)
```
Creates `connectivity` and `linksConnect` matrices for [`DislocationNetwork`](@ref). 

* `connectivity` contains the number of other other nodes each node is connected to, up to `maxConnect` other nodes. Every `(2i, j)` entry contains a node connected to node `j`. Every `(2i+1, j)` coordinate contains whether that node is the first or second node in the link.
* `linksConnect` relates connections enabled by a link. Analogous to the connectivity of a link.
"""
function genMakeConnect(mutating)
    if !mutating
        name = :makeConnect
        arguments = quote
            (links, maxConnect)
        end
        accumulator = quote
            links = arguments[1]
            maxConnect = arguments[2]
            lenLinks = size(links, 2)
            connectivity = zeros(Int, 1 + 2 * maxConnect, lenLinks)
            linksConnect = zeros(Int, 2, lenLinks)
        end
        retVal = quote
            return connectivity, linksConnect
        end
    else
        name = :makeConnect!
        arguments = quote
            (network::DislocationNetwork)
        end
        accumulator = quote
            network = arguments[1]
            links = network.links
            connectivity = network.connectivity
            linksConnect = network.linksConnect
            connectivity .= 0
            linksConnect .= 0
            lenLinks = size(links, 2)
        end
        retVal = quote
            return nothing
        end
    end

    ex = quote
        function $name(arguments...)
            $accumulator
            # Indices of defined links, undefined links are always at the end so we only need to know the first undefined entry.
            idx = findfirst(x -> x == 0, @view links[1, :])
            isnothing(idx) ? idx = lenLinks : idx -= 1
            # Loop through indices.
            @inbounds @simd for i in 1:idx
                # Node 1, it is the row of the coord matrix.
                n1 = links[1, i]
                n2 = links[2, i]
                # Num of other nodes node n1 is connected to.
                connectivity[1, n1] += 1
                connectivity[1, n2] += 1
                # Next column for n1 in connectivity.
                tmp1 = 2 * connectivity[1, n1]
                tmp2 = 2 * connectivity[1, n2]
                # i = linkID. 1, 2 are the first and second nodes in link with linkID
                connectivity[tmp1:(tmp1 + 1), n1] .= (i, 1)
                connectivity[tmp2:(tmp2 + 1), n2] .= (i, 2)
                # Connectivity of the nodes in link i.
                linksConnect[1, i] = connectivity[1, n1]
                linksConnect[2, i] = connectivity[1, n2]
            end
            return $retVal
        end
    end
    return ex
end
eval(genMakeConnect(true))
eval(genMakeConnect(false))

"""
```
getSegmentIdx(links, label)
```   
Finds the indices of a link and corresponding nodes.

## Returns

`n × 3` matrix is of the form `[i, node1, node2]`. 

* `i` can be used to find the Burgers vector, slip plane and segment forces of segment `i`, eg `bVec[:, i]`. 
* `node1` and `node2` can be used to find the coordinate and velocity of the nodes, eg `l = coord[:, node2] - coord[:, node1]`.
"""
function getSegmentIdx(links, label)
    lenLinks = size(links, 2)
    segIdx = zeros(Int, lenLinks, 3)  # Indexing matrix.

    # Find all defined nodes.
    idx = findfirst(x -> x == 0, @view links[1, :])
    isnothing(idx) ? idx = lenLinks : idx -= 1

    numSeg = 0 # Number of segments.

    # Loop through indices.
    @inbounds for i in 1:idx
        # Nodes.
        n1 = links[1, i]
        n2 = links[2, i]
        # Skip external nodes.
        (label[n1] == extDln || label[n2] == extDln) && continue
        numSeg += 1 # Increment index.
        # Indexing matrix, segment numSeg is made up of link i which is made up from nodes n1 and n2.
        segIdx[numSeg, :] .= (i, n1, n2)
    end

    return numSeg, segIdx
end
"""
```
getSegmentIdx!(network::DislocationNetwork)
```
Mutates the `segIdx` matrix in `network`. Works the same way as [`getSegmentIdx`](@ref).
"""
function getSegmentIdx!(network::DislocationNetwork)
    links = network.links
    label = network.label
    numSeg = network.numSeg
    segIdx = network.segIdx

    lenLinks = size(links, 2)
    lenSegIdx = size(segIdx, 1)
    segIdx .= 0

    idx = findfirst(x -> x == 0, @view links[1, :])
    isnothing(idx) ? idx = lenLinks : idx -= 1
    lNumSeg = 0

    @inbounds for i in 1:idx
        n1 = links[1, i]
        n2 = links[2, i]
        (label[n1] == extDln || label[n2] == extDln) && continue
        lNumSeg += 1
        segIdx[lNumSeg, :] .= (i, n1, n2)
    end

    numSeg[1] = lNumSeg
    return nothing
end

"""
```
checkNetwork(network::DislocationNetwork)
```
Checks the validity of the dislocation network. It ensures the following conditions are met by the member variables of `network`:
1. `connectivity` and `links` have the same number of non-zero entries;
1. all entries in `bVec` are non-zero;
1. only the trailing columns of `connectivity` are zeros;
1. consistency between `connectivity` and `links`;
1. `bVec` is conserved among connected nodes;
1. entries in `links` are unique;
1. consistency betwen `connectivity` and `linksConnect`

"""
function checkNetwork(network::DislocationNetwork)
    links = network.links
    label = network.label
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    idx = findfirst(x -> x == 0, @view links[1, :])
    isnothing(idx) ? idx = size(links, 2) : idx -= 1

    # Max value of connectivity should be the last link.
    maximum(connectivity) == idx ? nothing :
    error(
        "Non-empty entries of connectivity should be the same as the non-empty entries of links.",
    )

    bVec = network.bVec
    elemT = eltype(network.bVec)
    bSum = zeros(3)
    @inbounds for i in 1:idx
        iLinkBuffer = zeros(Int, 0)
        col = connectivity[1, i]

        # Check for zero Burgers vectors.
        staticBVec = SVector{3, elemT}(bVec[1, i], bVec[2, i], bVec[3, i])
        norm(staticBVec) < eps(elemT) ? error("Burgers vector must be non-zero.
                                       norm(bVec) = $(norm(bVec[i,:]))") : nothing

        # Trailing columns must be zero.
        sum(connectivity[(2 * (col + 1)):end, i]) != 0 ?
        error(
            "Trailing columns of connectivity must be zero.
      connectivity[$(i), $(2 * (col + 1)):end] = $(connectivity[i, (2 * (col + 1)):end]))",
        ) : nothing

        # Burgers vector conservation.
        bSum .= 0
        # No repeated links.
        neighbours = -ones(Int, col)
        for j in 1:col
            j2 = 2 * j
            iLink = connectivity[j2, i]     # Link ID.
            colLink = connectivity[j2 + 1, i] # Link position in links.
            colOppLink = 3 - connectivity[2 * j + 1, i] # Opposite column in links
            staticBVec = SVector{3, elemT}(bVec[1, iLink], bVec[2, iLink], bVec[3, iLink])
            bSum .+= (3 - 2 * colLink) * staticBVec # Running total of burgers vectors.

            neighbours[j] = links[colOppLink, iLink] # Find neighbouring nodes.
            push!(iLinkBuffer, iLink) # Push to link buffer for error printing.

            # Check that node does not appear twice in connectivity.
            @simd for k in (j + 1):col
                connectivity[j2, i] == connectivity[2 * k, i] ?
                error(
                    "Node $(i) cannot appear twice in connectivity, $(connectivity[i, j2])",
                ) : nothing
            end

            # Check links and connectivity are consistent.
            testNode = links[connectivity[j2 + 1, i], connectivity[j2, i]]
            testNode != i ?
            error(
                "Connectivity and links are not consistent.
connectivity[:, $(i)] = $(connectivity[:, i])
links[:, $(connectivity[j2, i])] = [$(links[1, connectivity[j2, i]]), $(links[2, connectivity[j2, i]])]",
            ) : nothing
        end

        # Check for Burgers vector conservation.
        abs(dot(bSum, bSum)) > eps(elemT) ?
        error("Burgers vector is not conserved, bSum = $(bSum).") : nothing

        # Find unique neighbours.
        uniqueNeighbours = unique(neighbours)
        length(uniqueNeighbours) != length(neighbours) ?
        error(
            "There must be no duplicate links. Each entry in links must correspond to a unique pair of adjacent entries in neighbours.
links = $(iLinkBuffer[1:end - 1])
uniqueNeighbours = $(uniqueNeighbours)
neighbours = $(neighbours)",
        ) : nothing

        # Check connectivity and linksConnect are consistent.
        check = (
            connectivity[2 * linksConnect[1, i], links[1, i]],
            connectivity[2 * linksConnect[2, i], links[2, i]],
        )
        i != check[1] | i != check[2] ? error("Inconsistent link.
        links[:, $(i)] = $(links[:, i])
        linksConnect[:, $(i)] = $(linksConnect[:, i])
        connectivity[:, $(links[1, i])] = $(connectivity[:, links[1, i]])
        connectivity[:, $(links[2, i])] = $(connectivity[:, links[2, i]])") : nothing
    end

    return true
end
