## Distributions
"""
```
loopDistribution(dist<:AbstractDistribution, n::Int, args...; kw...)
```
Returns `n` points according to the concrete subtype of [`AbstractDistribution`](@ref) given. Overload this function with new concrete subtypes and custom distributions. This and [`limits!`](@ref) are used in [`translatePoints`](@ref) to distribute dislocations in the simulation domain.
"""
loopDistribution(dist::Zeros, n::Int, args...; kw...) = zeros(3, n)
loopDistribution(dist::Rand, n::Int, args...; kw...) = rand(3, n)
loopDistribution(dist::Randn, n::Int, args...; kw...) = randn(3, n)
function loopDistribution(dist::Regular, n::Int, args...; kw...)
    return error("loopDistribution: regular distribution yet not implemented")
end

"""
```
limits!(
    lims::AbstractArray{<:Float64, N1},
    segLen::Float64,
    range::AbstractArray{<:Float64, N2},
    buffer::Float64,
) where {N1, N2}
```
Calculate the spatial limits a dislocation will occupy. This and [`loopDistribution`](@ref) are used in [`translatePoints`](@ref) to distribute dislocations in the simulation domain.
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
Translate dislocation node coordinates `coord` inside the spatial bounds of `lims` (calculated in [`limits!`](@ref)) according to the displacement `disp` (calculated in [`loopDistribution`](@ref)). Used to distribute sources inside a domain in [`makeNetwork`](@ref) and [`makeNetwork!`](@ref).
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

## Auxiliary matrices.
"""
```
makeConnect(
    links::T1,
    maxConnect::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: AbstractArray{T, N}}
```
Creates `connectivity` and `linksConnect` matrices. `connectivity` contains the number of other other nodes each node is connected to, up to `maxConnect` other nodes. It also contains the segments in which it's involved. `linksConnect` is the connectivity of each link. This is called from [`makeNetwork`](@ref) and [`makeNetwork!`](@ref).
"""
@inline function makeConnect(
    links::T1,
    maxConnect::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: Int}

    # Indices of defined links.
    idx = findall(x -> x != 0, links[1, :])
    lenLinks = size(links, 2)
    connectivity = zeros(Int, 1 + 2 * maxConnect, lenLinks)
    linksConnect = zeros(Int, 2, lenLinks)

    # Loop through indices.
    @inbounds @simd for i in idx
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
    maxConnect = network.maxConnect

    idx = findall(x -> x != 0, links[1,:])
    lenLinks = size(links, 2)
    connectivity = zeros(Int, 1 + 2 * maxConnect, lenLinks)
    linksConnect = zeros(Int, 2, lenLinks)
    @inbounds @simd for i in idx
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

    network.connectivity = connectivity
    network.linksConnect = linksConnect

    return network
end

"""
```
getSegmentIdx(
    links::T1,
    label::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: AbstractVector{nodeType}}
```
"""
@inline function getSegmentIdx(
    links::T1,
    label::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: AbstractVector{nodeType}}

    segIdx = zeros(Int, 3, size(links, 2))  # Indexing matrix.
    idx = findall(x -> x != 0, label)       # Find all defined nodes.
    numSeg::Int = 0 # Number of segments.

    # Loop through indices.
    for i in idx
        # Nodes.
        n1 = links[1, i]
        n2 = links[2, i]
        # Skip external nodes.
        (label[n1] == 4 || label[n2] == 4) ? continue : nothing
        numSeg += 1 # Increment index.
        # Indexing matrix, segment numSeg is made up of link i which is made up from nodes n1 and n2.
        segIdx[:, numSeg] = [i, n1, n2]
    end

    return numSeg, segIdx
end
"""
```
getSegmentIdx!(network::DislocationNetwork)
```
"""
@inline function getSegmentIdx!(network::DislocationNetwork)
    links = network.links
    label = network.label

    segIdx = zeros(Int, 3, size(links, 2))
    idx = findall(x -> x != 0, label)
    numSeg::Int = 0
    for i in idx
        n1 = links[1, i]
        n2 = links[2, i]
        (label[n1] == 4 || label[n2] == 4) ? continue : nothing
        numSeg += 1
        segIdx[:, numSeg] = [i, n1, n2]
    end

    network.numSeg = numSeg
    network.segIdx = segIdx

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
@inline function checkNetwork(network::DislocationNetwork)
    label = network.label
    idx = findall(x -> x != 0, label)
    links = network.links
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # Connectivity and links must be the same length.
    maximum(connectivity) == length(links[1, idx]) ? nothing :
    error("Non-empty entries of connectivity should be the same as the non-empty entries of links.")

    bVec = network.bVec
    bSum = zeros(3)
    @inbounds for i in idx
        iLinkBuffer = zeros(Int, 0)
        col = connectivity[1, i]

        # Check for zero Burgers vectors.
        norm(bVec[:, i]) < eps(eltype(bVec)) ? error("Burgers vector must be non-zero.
                                       norm(bVec) = $(norm(bVec[i,:]))") :
        nothing

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
            pLink = connectivity[j2 + 1, i] # Link position in links.
            bSum += (3 - 2 * pLink) .* bVec[:, iLink] # Running total of burgers vectors.

            neighbours[j] = links[3 - connectivity[2 * j + 1, i], iLink] # Find neighbouring nodes.
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
        abs(dot(bSum, bSum)) > eps(eltype(bSum)) ?
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
