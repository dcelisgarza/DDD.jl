## Distributions
"""
```
loopDistribution(dist<:AbstractDistribution, n::Int, args...; kw...)
loopDistribution(dist::Zeros, n::Int, args...; kw...) = zeros(3, n)
loopDistribution(dist::Rand, n::Int, args...; kw...) = rand(3, n)
loopDistribution(dist::Randn, n::Int, args...; kw...) = randn(3, n)
loopDistribution(dist::Regular, n::Int, args...; kw...) = error("loopDistribution: regular distribution yet not implemented")
```
Returns a `3 × n` matrix whose points follow the distribution defined by the method dispatching on the concrete subtype of `dist`. When creating custom [`AbstractDistribution`](@ref) subtypes, a corresponding `loopDistribution` method must be created for the custom distribution to be used in the [`DislocationLoop`](@ref) constructors.
"""
loopDistribution(dist::Zeros, n::Int, args...; kw...) = zeros(3, n)
loopDistribution(dist::Rand, n::Int, args...; kw...) = rand(3, n)
loopDistribution(dist::Randn, n::Int, args...; kw...) = randn(3, n)
function loopDistribution(dist::Regular, n::Int, args...; kw...)
    error("loopDistribution: regular distribution yet not implemented")
end

"""
```
limits!(
    lims::T1,
    segLen::T2,
    range::T1,
    buffer::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2}
```
In-place addition of `buffer × segLen` to `range` in order to calculate the limits in which dislocations will exist.
"""
@inline function limits!(
    lims::T1,
    segLen::T2,
    range::T1,
    buffer::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2}

    @inbounds for i in 1:size(lims, 2)
        @simd for j in 1:size(lims, 1)
            lims[j, i] = range[j, i] + buffer * segLen
        end
    end
    return lims
end

"""
```
translatePoints(
    coord::T1,
    lims::T1,
    disp::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: AbstractVector{T} where {T}}
```
Translates coordinates using the limits and displacements calculated by [`limits!`](@ref) and [`loopDistribution`](@ref).
"""
@inline function translatePoints(
    coord::T1,
    lims::T1,
    disp::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: AbstractVector{T} where {T}}

    @inbounds for i in 1:size(coord, 2)
        @simd for j in 1:size(coord, 1)
            coord[j, i] += lims[j, 1] + (lims[j, 2] - lims[j, 1]) * disp[j]
        end
    end
    return coord
end

## Segments
"""
```
makeSegment(
    type::T1,
    slipPlane::T2,
    bVec::T2
) where {T1 <: AbstractDlnSeg, T2 <: AbstractVector{T} where {T}}
```
Make signle segment depending on the segment type, see [`AbstractDlnSeg`](@ref).
"""
@inline function makeSegment(
    type::T1,
    slipPlane::T2,
    bVec::T2,
) where {T1 <: segEdge, T2 <: AbstractVector{T} where {T}}
    edge = cross(slipPlane, bVec)
    return edge ./ norm(edge)
end
@inline function makeSegment(
    type::T1,
    slipPlane::T2,
    bVec::T2,
) where {T1 <: segEdgeN, T2 <: AbstractVector{T} where {T}}
    return slipPlane ./ norm(slipPlane)
end
@inline function makeSegment(
    type::T1,
    slipPlane::T2,
    bVec::T2,
) where {T1 <: segScrew, T2 <: AbstractVector{T} where {T}}
    return bVec ./ norm(bVec)
end

## Auxiliary matrices.
"""
```
makeConnect(
    links::T1,
    maxConnect::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: Int}
```
Creates `connectivity` and `linksConnect` matrices. `connectivity` contains the number of other other nodes each node is connected to, up to `maxConnect` other nodes. Every `(2i, j)` entry contains a node connected to node `j`. Every `(2i+1, j)` coordinate contains whether that node is the first or second node in the link.

The matrix `linksConnect` relates connections enabled by a link. Analogous to the connectivity of a link.
"""
@inline function makeConnect(
    links::T1,
    maxConnect::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: Int}

    lenLinks = size(links, 2)

    # Indices of defined links, undefined links are always at the end so we only need to know the first undefined entry.
    idx = findfirst(x -> x == 0, links[1, :])
    isnothing(idx) ? idx = lenLinks : idx -= 1

    connectivity = zeros(Int, 1 + 2 * maxConnect, lenLinks)
    linksConnect = zeros(Int, 2, lenLinks)

    # Loop through indices.
    for i in 1:idx
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
        connectivity[tmp1:(tmp1 + 1), n1] = [i, 1]
        connectivity[tmp2:(tmp2 + 1), n2] = [i, 2]

        # Connectivity of the nodes in link i.
        linksConnect[1, i] = connectivity[1, n1]
        linksConnect[2, i] = connectivity[1, n2]
    end

    return connectivity, linksConnect
end
"""
```
makeConnect!(network::DislocationNetwork)
```
In-place version of [`makeConnect`](@ref).
"""
@inline function makeConnect!(network::DislocationNetwork)

    # For comments see makeConnect. It is a 1-to-1 translation except that this one modifies the network in-place.
    links = network.links
    connectivity = network.connectivity
    linksConnect = network.linksConnect
    connectivity .= 0
    linksConnect .= 0
    maxConnect = network.numNodeSegConnect[3]
    lenLinks = size(links, 2)
    idx = findfirst(x -> x == 0, links[1, :])
    isnothing(idx) ? idx = lenLinks : idx -= 1
    for i in 1:idx
        n1 = links[1, i]
        n2 = links[2, i]
        connectivity[1, n1] += 1
        connectivity[1, n2] += 1
        tmp1 = 2 * connectivity[1, n1]
        tmp2 = 2 * connectivity[1, n2]
        connectivity[tmp1:(tmp1 + 1), n1] = [i, 1]
        connectivity[tmp2:(tmp2 + 1), n2] = [i, 2]
        linksConnect[1, i] = connectivity[1, n1]
        linksConnect[2, i] = connectivity[1, n2]
    end

    return network
end

"""
```
getSegmentIdx(
    links::T1,
    label::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: AbstractVector{nodeType}}
```
Creates an indexing matrix for quick indexing of dislocation segments for quick access to slip planes, burgers vectors and line vectors. The return `3 × n` matrix is of the form `[i, node1, node2]`. Index `i` can be used to find the Burgers vector, slip plane and segment forces of a segment, eg `bVec[:, i]`. While `node1` and `node2` can be used to find the coordinate and velocity of the nodes, eg `t = coord[:, node2] - coord[:, node1]`.
"""
@inline function getSegmentIdx(
    links::T1,
    label::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: AbstractVector{nodeType}}

    lenLinks = size(links, 2)
    segIdx = zeros(Int, lenLinks, 3)  # Indexing matrix.

    # Find all defined nodes.
    idx = findfirst(x -> x == 0, links[1, :])
    isnothing(idx) ? idx = lenLinks : idx -= 1

    numSeg::Int = 0 # Number of segments.

    # Loop through indices.
    for i in 1:idx
        # Nodes.
        n1 = links[1, i]
        n2 = links[2, i]
        # Skip external nodes.
        (label[n1] == 4 || label[n2] == 4) ? continue : nothing
        numSeg += 1 # Increment index.
        # Indexing matrix, segment numSeg is made up of link i which is made up from nodes n1 and n2.
        segIdx[numSeg, :] = [i, n1, n2]
    end

    return numSeg, segIdx
end
"""
```
getSegmentIdx!(network::DislocationNetwork)
```
In-place version of [`getSegmentIdx`](@ref).
"""
@inline function getSegmentIdx!(network::T1) where {T1 <: DislocationNetwork}

    links = network.links
    label = network.label


    lenLinks = size(links, 2)
    lenSegIdx = size(network.segIdx, 1)
    network.segIdx .= 0

    idx = findfirst(x -> x == 0, links[1, :])
    isnothing(idx) ? idx = lenLinks : idx -= 1
    numSeg::Int = 0
    for i in 1:idx
        n1 = links[1, i]
        n2 = links[2, i]
        (label[n1] == 4 || label[n2] == 4) ? continue : nothing
        numSeg += 1
        network.segIdx[numSeg, :] = [i, n1, n2]
    end
    network.numNodeSegConnect[2] = numSeg
    return network
end

## Check integrity.
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
@inline function checkNetwork(network::T1) where {T1 <: DislocationNetwork}

    links = network.links
    label = network.label
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    idx = findfirst(x -> x == 0, links[1, :])
    isnothing(idx) ? idx = size(links, 2) : idx -= 1

    # Max value of connectivity should be the last link.
    maximum(connectivity) == idx ? nothing :
    error("Non-empty entries of connectivity should be the same as the non-empty entries of links.")

    bVec = network.bVec
    elemT = eltype(network.bVec)
    bSum = zeros(3)
    @inbounds for i in 1:idx
        iLinkBuffer = zeros(Int, 0)
        col = connectivity[1, i]

        # Check for zero Burgers vectors.
        norm(bVec[:, i]) < eps(elemT) ? error("Burgers vector must be non-zero.
                                       norm(bVec) = $(norm(bVec[i,:]))") : nothing

        # Trailing columns must be zero.
        sum(connectivity[(2 * (col + 1)):end, i]) != 0 ?
        error("Trailing columns of connectivity must be zero.
        connectivity[$(i), $(2*(col + 1)):end] = $(connectivity[i, (2 * (col + 1)):end]))") :
        nothing

        # Burgers vector conservation.
        bSum .= 0
        # No repeated links.
        neighbours = -ones(Int, col)
        for j in 1:col
            j2 = 2 * j
            iLink = connectivity[j2, i]     # Link ID.
            colLink = connectivity[j2 + 1, i] # Link position in links.
            colOppLink = 3 - connectivity[2 * j + 1, i] # Opposite column in links
            bSum += (3 - 2 * colLink) * bVec[:, iLink] # Running total of burgers vectors.

            neighbours[j] = links[colOppLink, iLink] # Find neighbouring nodes.
            push!(iLinkBuffer, iLink) # Push to link buffer for error printing.

            # Check that node does not appear twice in connectivity.
            for k in (j + 1):col
                connectivity[j2, i] == connectivity[2 * k, i] ?
                error("Node $(i) cannot appear twice in connectivity, $(connectivity[i, j2])") :
                nothing
            end

            # Check links and connectivity are consistent.
            testNode = links[connectivity[j2 + 1, i], connectivity[j2, i]]
            testNode != i ?
            error("Connectivity and links are not consistent.
connectivity[:, $(i)] = $(connectivity[:, i])
links[:, $(connectivity[j2, i])] = [$(links[1, connectivity[j2, i]]), $(links[2, connectivity[j2, i]])]") :
            nothing
        end

        # Check for Burgers vector conservation.
        abs(dot(bSum, bSum)) > eps(elemT) ?
        error("Burgers vector is not conserved, bSum = $(bSum).") : nothing

        # Find unique neighbours.
        uniqueNeighbours = unique(neighbours)
        length(uniqueNeighbours) != length(neighbours) ?
        error("There must be no duplicate links. Each entry in links must correspond to a unique pair of adjacent entries in neighbours.
links = $(iLinkBuffer[1:end-1])
uniqueNeighbours = $(uniqueNeighbours)
neighbours = $(neighbours)") :
        nothing

        # Check connectivity and linksConnect are consistent.
        check = [
            connectivity[2 * linksConnect[1, i], links[1, i]]
            connectivity[2 * linksConnect[2, i], links[2, i]]
        ]
        i != check[1] | i != check[2] ? error("Inconsistent link.
        links[:, $(i)] = $(links[:, i])
        linksConnect[:, $(i)] = $(linksConnect[:, i])
        connectivity[:, $(links[1, i])] = $(connectivity[:, links[1, i]])
        connectivity[:, $(links[2, i])] = $(connectivity[:, links[2, i]])") : nothing
    end

    return true
end
