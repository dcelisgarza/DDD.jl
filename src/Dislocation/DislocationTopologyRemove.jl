"""
```
removeNode!(network::DislocationNetwork, nodeGone::Int, lastNode = nothing)
```
In-place remove `nodeGone`. If `nodeGone` is *not* the last node, this replaces the entries corresponding to `nodeGone` with `lastNode`. Else it simply zeros out the entries corresponding to `lastNode`. This also decreases `numNode` by one. If `lastNode` is nothing, this function finds it.

Modifies
```
network.links
network.coord
network.label
network.nodeVel
network.numNode
network.connectivity
```
"""
@inline function removeNode!(
    network::T1,
    nodeGone::T2,
    lastNode = nothing,
) where {T1 <: DislocationNetwork, T2 <: Int}

    links = network.links
    coord = network.coord
    label = network.label
    nodeVel = network.nodeVel
    nodeForce = network.nodeForce
    connectivity = network.connectivity

    isnothing(lastNode) ? lastNode = maximum((network.numNode, 1)) : nothing

    # The if nodeGone is not the last node in the relevant arrays, replace it by lastNode.
    if nodeGone < lastNode
        coord[:, nodeGone] .= coord[:, lastNode]
        label[nodeGone] = label[lastNode]
        nodeVel[:, nodeGone] .= nodeVel[:, lastNode]
        nodeForce[:, nodeGone] .= nodeForce[:, lastNode]
        connectivity[:, nodeGone] .= connectivity[:, lastNode]
        # Change the link
        @simd for j in 1:connectivity[1, nodeGone]
            idx = 2 * j
            links[connectivity[idx + 1, nodeGone], connectivity[idx, nodeGone]] = nodeGone
        end
    end

    # Remove the last node from the arrays since nodeGone was either replaced by lastNode above, or already was the last node. Reduce the number of nodes by one.
    coord[:, lastNode] .= 0
    label[lastNode] = 0
    nodeVel[:, lastNode] .= 0
    nodeForce[:, lastNode] .= 0
    network.numNode -= 1
    connectivity[:, lastNode] .= 0

    return network
end

"""
```
removeConnection!(network::DislocationNetwork, nodeKept::Int, connectGone::Int)
```
In-place remove `connectGone` from the entries corresponding to `nodeKept`. If `connectGone` is not the last connection, it gets replaced by the last connection. Else it simply zeroes out the entries corresponding to `connectGone`. Also reduces the number of connections by one.

Modifies
```
network.connectivity
network.linksConnect
```
"""
@inline function removeConnection!(
    network::T1,
    nodeKept::T2,
    connectGone::T2,
) where {T1 <: DislocationNetwork, T2 <: Int}

    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # Remove entry from linksConnect.
    idx = 2 * connectGone
    link1 = connectivity[idx, nodeKept]
    @assert link1 != 0 "removeConnection!: connection $connectGone for node $nodeKept does not exist."
    link2 = connectivity[idx + 1, nodeKept]
    linksConnect[link2, link1] = 0

    # The last connection to be added to nodeKept is equal to the connectivity of the nodeKept.
    lastConnect = connectivity[1, nodeKept]
    lst = 2 * lastConnect

    # If the connectGone is not the last connection that was made, replace it by the last connection.
    if connectGone < lastConnect
        connectivity[idx:(idx + 1), nodeKept] = connectivity[lst:(lst + 1), nodeKept]
        linksConnect[connectivity[idx + 1, nodeKept], connectivity[idx, nodeKept]] =
            connectGone
    end

    # Change connectivity to reflect that nodeKept has one less connection.
    # Remove the last connection since lastConnect was either replaced by connectGone above, or already was the last connection.
    connectivity[1, nodeKept] -= 1
    connectivity[lst:(lst + 1), nodeKept] .= 0

    return network
end

"""
```
removeLink!(network::DislocationNetwork, linkGone::Int, lastLink = nothing)
```
In-place remove `linkGone` with the last valid link. If `linkGone` is *not* the last link, this replaces the entries corresponding to `linkGone` with `lastLink`. Else it simply zeros out the entries corresponding to `lastLink`. This also decreases `numSeg` by one. If `lastLink` is nothing, this function finds it.

Modifies
```
network.links
network.slipPlane
network.bVec
network.coord
network.label
network.segForce
network.nodeVel
network.numSeg
network.connectivity
network.linksConnect
```
"""
function removeLink!(
    network::T1,
    linkGone::T2,
    lastLink = nothing,
) where {T1 <: DislocationNetwork, T2 <: Int}

    links = network.links
    slipPlane = network.slipPlane
    bVec = network.bVec
    segForce = network.segForce
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    @assert linkGone <= size(links, 2) "removeLink!: link $linkGone not found."

    # Delete linkGone from connectivity for both nodes in the link.
    node1 = links[1, linkGone]
    connectGone1 = linksConnect[1, linkGone]
    removeConnection!(network, node1, connectGone1)

    node2 = links[2, linkGone]
    connectGone2 = linksConnect[2, linkGone]
    removeConnection!(network, node2, connectGone2)

    # Remove link that no longer appears in connectivity and doesn't connect any nodes.
    @assert linksConnect[:, linkGone]' * linksConnect[:, linkGone] == 0 "removeLink!: link $linkGone still has connections and should not be deleted."

    isnothing(lastLink) ? lastLink = maximum((network.numSeg, 1)) : nothing

    # If the linkGone is not the last link, replace it with lastLink.
    if linkGone < lastLink
        links[:, linkGone] = links[:, lastLink]
        slipPlane[:, linkGone] = slipPlane[:, lastLink]
        bVec[:, linkGone] = bVec[:, lastLink]
        segForce[:, :, linkGone] = segForce[:, :, lastLink]
        linksConnect[:, linkGone] = linksConnect[:, lastLink]
        node1 = links[1, linkGone]
        node2 = links[2, linkGone]
        connect1 = 2 * linksConnect[1, linkGone]
        connect2 = 2 * linksConnect[2, linkGone]
        connectivity[connect1, node1] = linkGone
        connectivity[connect2, node2] = linkGone
    end

    # Remove the last link from the arrays since linkGone was either replaced by lastLink above, or already was the last link.
    links[:, lastLink] .= 0
    slipPlane[:, lastLink] .= 0
    bVec[:, lastLink] .= 0
    segForce[:, :, lastLink] .= 0
    network.numSeg -= 1
    linksConnect[:, lastLink] .= 0

    return network
end

"""
```
mergeNode!(network::DislocationNetwork, nodeKept::Int, nodeGone::Int)
```
Merges `nodeGone` into `nodeKept`. After calling this function there are no repeated entries, self-links or double links.
"""
@inline @fastmath function mergeNode!(
    network::T1,
    nodeKept::T2,
    nodeGone::T2,
) where {T1 <: DislocationNetwork, T2 <: Int}

    @assert nodeKept <= network.numNode && nodeGone <= network.numNode "mergeNode!: the node kept after merging, $nodeKept and node removed after merging, $nodeGone, must be in the simulation."

    # Return if both nodes to be merged are the same.
    nodeKept == nodeGone && return 0

    links = network.links
    slipPlane = network.slipPlane
    bVec = network.bVec
    connectivity = network.connectivity
    linksConnect = network.linksConnect
    elemT = eltype(network.bVec)

    nodeKeptConnect = connectivity[1, nodeKept]
    nodeGoneConnect = connectivity[1, nodeGone]
    totalConnect = nodeKeptConnect + nodeGoneConnect

    # Pass connections from nodeGone to nodeKept.
    if size(connectivity, 1) < 2 * totalConnect + 1
        network.connectivity = vcat(
            network.connectivity,
            zeros(
                Int,
                2 * totalConnect + 1 - size(network.connectivity, 1),
                length(network.label),
            ),
        )
        connectivity = network.connectivity
    end
    connectivity[(2 * (nodeKeptConnect + 1)):(2 * totalConnect + 1), nodeKept] =
        connectivity[2:(2 * nodeGoneConnect + 1), nodeGone]
    connectivity[1, nodeKept] = totalConnect

    # Replace nodeGone with nodeKept in links and update linksConnect with the new positions of the links in connectivity.
    @inbounds @simd for i in 1:nodeGoneConnect
        link = connectivity[2 * i, nodeGone]      # Link id for nodeGone
        colLink = connectivity[2 * i + 1, nodeGone] # Position of links where link appears.
        links[colLink, link] = nodeKept         # Replace nodeGone with nodeKept.
        linksConnect[colLink, link] = nodeKeptConnect + i # Increase connectivity of nodeKept.
    end

    # Remove node from network and update index of nodeKept in case it changed.
    lastNode = maximum((network.numNode, 1))
    removeNode!(network, nodeGone, lastNode)
    coord = network.coord
    nodeVel = network.nodeVel
    connectivity = network.connectivity
    nodeKept == lastNode ? nodeKept = nodeGone : nothing

    # Warning: connectivity can be updated by removeLink!().
    # Delete self links of nodeKept.
    i = 1
    @inbounds while i < connectivity[1, nodeKept]
        link = connectivity[2 * i, nodeKept]      # Link id for nodeKept
        colLink = connectivity[2 * i + 1, nodeKept] # Position of links where link appears.
        colNotLink = 3 - colLink # The column of the other node in the position of link in links.
        nodeNotLink = links[colNotLink, link]   # Other node in the link.
        # If nodeKept and nodeNotLink are the same, this is a self-link, therefore delete. Else increase counter.
        nodeKept == nodeNotLink ? removeLink!(network, link) : i += 1
        links = network.links
        connectivity = network.connectivity
    end

    # Warning: connectivity can be updated in the while loop.
    # Delete duplicate links.
    i = 1
    @inbounds while i < connectivity[1, nodeKept]
        link1 = connectivity[2 * i, nodeKept]         # Link id for nodeKept
        colLink1 = connectivity[2 * i + 1, nodeKept]    # Position of links where link appears.
        colNotLink1 = 3 - colLink1 # The column of the other node in the position of link in links.
        nodeNotLink1 = links[colNotLink1, link1]    # Other node in the link.
        # Same but for the next link.
        j = i + 1
        while j <= connectivity[1, nodeKept]
            link2 = connectivity[2 * j, nodeKept]
            colLink2 = connectivity[2 * j + 1, nodeKept]
            colNotLink2 = 3 - colLink2
            nodeNotLink2 = links[colNotLink2, link2]

            # Continue to next iteration if no duplicate links are found.
            if nodeNotLink1 != nodeNotLink2
                j += 1
                continue
            end

            # Fix Burgers vector.
            # If the nodes are on different ends of a link, conservation of Burgers vector in a loop requires we subtract contributions from the two links involved. Else we add them.
            colNotLink1 + colNotLink2 == 3 ? bVec[:, link1] -= bVec[:, link2] :
            bVec[:, link1] += bVec[:, link2]

            # WARNING This calculation is odd. Try using the cross product of the adjacent segments.
            # Fix slip plane.
            # Line direction and velocity of the resultant dislocation.
            t = SVector{3, elemT}(
                coord[1, nodeKept] - coord[1, nodeNotLink1],
                coord[2, nodeKept] - coord[2, nodeNotLink1],
                coord[3, nodeKept] - coord[3, nodeNotLink1],
            )

            v = SVector{3, elemT}(
                nodeVel[1, nodeKept] + nodeVel[1, nodeNotLink1],
                nodeVel[2, nodeKept] + nodeVel[2, nodeNotLink1],
                nodeVel[3, nodeKept] + nodeVel[3, nodeNotLink1],
            )

            # Burgers vector and potential new slip plane.
            b = SVector{3, elemT}(bVec[1, link1], bVec[2, link1], bVec[3, link1])
            n1 = t × b  # For non-screw segments.
            n2 = t × v  # For screw segments.
            if n1 ⋅ n1 > eps(elemT) # non-screw
                slipPlane[:, link1] = n1 / norm(n1)
            elseif n2 ⋅ n2 > eps(elemT) # screw
                slipPlane[:, link1] = n2 / norm(n2)
            end

            # Remove link2 from network and update the index of link1 in case it changed.
            lastLink = maximum((network.numSeg, 1))
            removeLink!(network, link2, lastLink)
            links = network.links
            connectivity = network.connectivity
            link1 == lastLink ? link1 = link2 : nothing

            # If the burgers vector of the new junction is non-zero, continue to the next iteration. Else remove it.
            b = SVector{3, elemT}(bVec[1, link1], bVec[2, link1], bVec[3, link1])
            if isapprox(dot(b, b), 0)
                removeLink!(network, link1)
                links = network.links
                connectivity = network.connectivity
                # If the node that was connected to nodeKept has no connections, remove it from the network and update the index of nodeKept in case it changed.
                if connectivity[1, nodeNotLink1] == 0
                    lastNode = maximum((network.numNode, 1))
                    removeNode!(network, nodeNotLink1, lastNode)
                    coord = network.coord
                    nodeVel = network.nodeVel
                    connectivity = network.connectivity
                    nodeKept == lastNode ? nodeKept = nodeNotLink1 : nothing
                end
            end

            # Since there was a change in the configuration, we check for duplicates again.
            i -= 1
            break
        end
        i += 1
    end

    # If nodeKept has no connections remove it.
    if connectivity[1, nodeKept] == 0
        removeNode!(network, nodeKept)
        nodeKept = 0
    end

    return nodeKept
end

function coarsenNetwork!(
    dlnParams::T1,
    matParams::T2,
    network::T3;
    # mesh::RegularCuboidMesh,
    # dlnFEM::DislocationFEMCorrective;
    parallel::Bool = false,
) where {T1 <: DislocationP, T2 <: MaterialP, T3 <: DislocationNetwork}

    minAreaSq = dlnParams.minAreaSq
    minSegLen = dlnParams.minSegLen
    maxSegLen = dlnParams.maxSegLen

    label = network.label
    links = network.links
    coord = network.coord
    nodeVel = network.nodeVel
    connectivity = network.connectivity
    elemT = eltype(network.coord)

    i = 1
    while i <= network.numNode
        # We only want to coarsen real internal nodes. Else we skip to the next node.
        if !(connectivity[1, i] == 2 && label[i] == 1)
            i += 1
            continue
        end
        link1 = connectivity[2, i]  # Node i in link 1
        link2 = connectivity[4, i]  # Node i in link 2

        colInLink1 = connectivity[3, i] # Column of node i in link 1
        colInLink2 = connectivity[5, i] # Column of node i in link 2

        oppColLink1 = 3 - colInLink1 # Node i is connected via link 1 to the node that is in this column in links.
        oppColLink2 = 3 - colInLink2 # Node i is connected via link 2 to the node that is in this column in links.

        link1_nodeOppI = links[oppColLink1, link1] # Node i is connected to this node as part of link 1.
        link2_nodeOppI = links[oppColLink2, link2] # Node i is connected to this node as part of link 2.

        # We don't want to remesh out segments between two fixed nodes because the nodes by definition do not move and act as a source, thus we skip to the next node.
        if label[link1_nodeOppI] == 2 && label[link2_nodeOppI] == 2
            i += 1
            continue
        end

        # Coordinate of node i
        iCoord = SVector{3, elemT}(coord[1, i], coord[2, i], coord[3, i])
        # Create a triangle formed by the three nodes involved in coarsening.
        coordVec1 =
            SVector{3, elemT}(
                coord[1, link1_nodeOppI],
                coord[2, link1_nodeOppI],
                coord[3, link1_nodeOppI],
            ) - iCoord # Vector between node 1 and the node it's connected to via link 1.
        coordVec2 =
            SVector{3, elemT}(
                coord[1, link2_nodeOppI],
                coord[2, link2_nodeOppI],
                coord[3, link2_nodeOppI],
            ) - iCoord # Vector between node 1 and the node it's connected to via link 2.
        coordVec3 = coordVec2 - coordVec1 # Vector between both nodes connected to iCoord.
        # Lengths of the triangle sides.
        r1 = norm(coordVec1)
        r2 = norm(coordVec2)
        r3 = norm(coordVec3)

        # If coarsening would result in a link whose length is bigger than the maximum allowed we skip ahead to the next node.
        if r3 >= maxSegLen
            i += 1
            continue
        end

        # Half the triangle's perimeter.
        r0 = (r1 + r2 + r3) / 2
        # Heron's formula for the area of a triangle when we know its sides.
        areaSq = r0 * (r0 - r1) * (r0 - r2) * (r0 - r3)

        # Node i velocities.
        iVel = SVector{3, elemT}(nodeVel[1, i], nodeVel[2, i], nodeVel[3, i])
        velVec1 =
            SVector{3, elemT}(
                nodeVel[1, link1_nodeOppI],
                nodeVel[2, link1_nodeOppI],
                nodeVel[3, link1_nodeOppI],
            ) - iVel
        velVec2 =
            SVector{3, elemT}(
                nodeVel[1, link2_nodeOppI],
                nodeVel[2, link2_nodeOppI],
                nodeVel[3, link2_nodeOppI],
            ) - iVel
        velVec3 = velVec2 - velVec1

        # Rate of change of side length with respect to time. We add eps(typeof(r)) to avoid division by zero.
        dr1dt = coordVec1 ⋅ velVec1 / (r1 + eps(elemT))
        dr2dt = coordVec2 ⋅ velVec2 / (r2 + eps(elemT))
        dr3dt = coordVec3 ⋅ velVec3 / (r3 + eps(elemT))
        dr0dt = (dr1dt + dr2dt + dr3dt) / 2

        # Rate of change of triangle area with respect to time.
        dAreaSqdt = dr0dt * (r0 - r1) * (r0 - r2) * (r0 - r3)
        dAreaSqdt += r0 * (dr0dt - dr1dt) * (r0 - r2) * (r0 - r3)
        dAreaSqdt += r0 * (r0 - r1) * (dr0dt - dr2dt) * (r0 - r3)
        dAreaSqdt += r0 * (r0 - r1) * (r0 - r2) * (dr0dt - dr3dt)

        # Coarsen if the area is less than the minimum allowed and shrinking. Or if the length of either of its links is less than the minimum allowed. Else continue to the next node.
        if !(areaSq < minAreaSq && dAreaSqdt < 0 || r1 < minSegLen || r2 < minSegLen)
            i += 1
            continue
        end

        # Merge node i into node link2_nodeOppI.
        nodeMerged = mergeNode!(network, link2_nodeOppI, i)
        getSegmentIdx!(network)
        links = network.links
        coord = network.coord
        label = network.label
        nodeVel = network.nodeVel
        connectivity = network.connectivity

        # If link2_nodeOppI no longer exists there is nothing to calculate and we proceed to the next iteration.
        if nodeMerged == 0
            i += 1
            continue
        end

        for j in 1:connectivity[1, nodeMerged]
            # Find the new link that has been created between nodeMerged and nodeNotMerged.
            linkMerged = connectivity[2 * j, nodeMerged]
            colLinkMerged = connectivity[2 * j + 1, nodeMerged]
            colNotLinkMerged = 3 - colLinkMerged
            nodeNotMerged = links[colNotLinkMerged, linkMerged]

            if nodeNotMerged == i || nodeNotMerged == link1_nodeOppI
                # Calculate segment force for segment linkMerged.
                calcSegForce!(dlnParams, matParams, network, linkMerged)
                # Calculate node velocity.
                nodes = links[:, linkMerged]
                missing, nodeVel[:, nodes] =
                    dlnMobility(dlnParams, matParams, network, nodes)
            end
        end
    end
    return network
end
