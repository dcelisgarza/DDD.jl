## Topology functions

## Deletion
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
    network::DislocationNetwork,
    nodeGone::Int,
    lastNode::Int = nothing,
)
    links = network.links
    coord = network.coord
    label = network.label
    nodeVel = network.nodeVel
    connectivity = network.connectivity

    isnothing(lastNode) ? lastNode = maximum((network.numNode, 1)) : nothing

    # The if nodeGone is not the last node in the relevant arrays, replace it by lastNode.
    if nodeGone < lastNode
        coord[nodeGone, :] .= coord[lastNode, :]
        label[nodeGone] = label[lastNode]
        nodeVel[nodeGone, :] .= nodeVel[lastNode, :]
        connectivity[nodeGone, :] .= connectivity[lastNode, :]
        # Change the link
        @simd for j in 1:connectivity[nodeGone, 1]
            idx = 2 * j
            links[connectivity[nodeGone, idx], connectivity[nodeGone, idx + 1]] = nodeGone
        end
    end

    # Remove the last node from the arrays since nodeGone was either replaced by lastNode above, or already was the last node. Reduce the number of nodes by one.
    coord[lastNode, :] .= 0
    label[lastNode] = 0
    nodeVel[lastNode, :] .= 0
    network.numNode -= 1
    connectivity[lastNode, :] .= 0

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
    network::DislocationNetwork,
    nodeKept::Int,
    connectGone::Int,
)
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # Remove entry from linksConnect.
    idx = 2 * connectGone
    link1 = connectivity[nodeKept, idx]
    @assert link1 != 0 "removeConnection!: connection $connectGone for node $nodeKept does not exist."
    link2 = connectivity[nodeKept, idx + 1]
    linksConnect[link1, link2] = 0

    # The last connection to be added to nodeKept is equal to the connectivity of the nodeKept.
    lastConnect = connectivity[nodeKept, 1]
    lst = 2 * lastConnect

    # If the connectGone is not the last connection that was made, replace it by the last connection.
    if connectGone < lastConnect
        connectivity[nodeKept, idx:(idx + 1)] = connectivity[nodeKept, lst:(lst + 1)]
        linksConnect[connectivity[nodeKept, idx], connectivity[nodeKept, idx + 1]] =
            connectGone
    end

    # Change connectivity to reflect that nodeKept has one less connection.
    # Remove the last connection since lastConnect was either replaced by connectGone above, or already was the last connection.
    connectivity[nodeKept, 1] -= 1
    connectivity[nodeKept, lst:(lst + 1)] .= 0

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
function removeLink!(network::DislocationNetwork, linkGone::Int, lastLink = nothing)
    links = network.links
    slipPlane = network.slipPlane
    bVec = network.bVec
    coord = network.coord
    label = network.label
    segForce = network.segForce
    nodeVel = network.nodeVel
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    @assert linkGone <= size(links, 1) "removeLink!: link $linkGone not found."

    # Delete linkGone from connectivity for both nodes in the link.
    node1 = links[linkGone, 1]
    connectGone1 = linksConnect[linkGone, 1]
    removeConnection!(network, node1, connectGone1)

    node2 = links[linkGone, 2]
    connectGone2 = linksConnect[linkGone, 2]
    removeConnection!(network, node2, connectGone2)

    # Remove link that no longer appears in connectivity and doesn't connect any nodes.
    @assert linksConnect[linkGone, :]' * linksConnect[linkGone, :] == 0 "removeLink!: link $linkGone still has connections and should not be deleted."

    isnothing(lastLink) ? lastLink = maximum((network.numSeg, 1)) : nothing

    # If the linkGone is not the last link, replace it with lastLink.
    if linkGone < lastLink
        links[linkGone, :] = links[lastLink, :]
        slipPlane[linkGone, :] = slipPlane[lastLink, :]
        bVec[linkGone, :] = bVec[lastLink, :]
        segForce[linkGone, :] = segForce[lastLink, :]
        linksConnect[linkGone, :] = linksConnect[lastLink, :]
        node1 = links[linkGone, 1]
        node2 = links[linkGone, 2]
        connect1 = 2 * linksConnect[linkGone, 1]
        connect2 = 2 * linksConnect[linkGone, 2]
        connectivity[node1, connect1] = linkGone
        connectivity[node2, connect2] = linkGone
    end

    # Remove the last link from the arrays since linkGone was either replaced by lastLink above, or already was the last link.
    links[lastLink, :] .= 0
    slipPlane[lastLink, :] .= 0
    bVec[lastLink, :] .= 0
    segForce[lastLink, :] .= 0
    network.numSeg -= 1
    linksConnect[lastLink, :] .= 0

    return network
end

"""
```
mergeNode!(network::DislocationNetwork, nodeKept::Int, nodeGone::Int)
```
Merges `nodeGone` into `nodeKept`. After calling this function there are no repeated entries, self-links or double links.
"""
@inline @fastmath function mergeNode!(
    network::DislocationNetwork,
    nodeKept::Int,
    nodeGone::Int,
)
    @assert nodeKept <= network.numNode && nodeGone <= network.numNode "mergeNode!: the node kept after merging, $nodeKept and node removed after merging, $nodeGone, must be in the simulation."

    # Return if both nodes to be merged are the same.
    nodeKept == nodeGone && return

    links = network.links
    slipPlane = network.slipPlane
    bVec = network.bVec
    coord = network.coord
    label = network.label
    nodeVel = network.nodeVel
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    nodeKeptConnect = connectivity[nodeKept, 1]
    nodeGoneConnect = connectivity[nodeGone, 1]
    totalConnect = nodeKeptConnect + nodeGoneConnect

    # Pass connections from nodeGone to nodeKept.
    if size(connectivity, 2) < 2 * totalConnect + 1
        network.connectivity = hcat(
            network.connectivity,
            zeros(Int, length(label), 2 * totalConnect + 1 - size(network.connectivity, 2)),
        )
        connectivity = network.connectivity
    end
    connectivity[nodeKept, (2 * (nodeKeptConnect + 1)):(2 * totalConnect + 1)] =
        connectivity[nodeGone, 2:(2 * nodeGoneConnect + 1)]
    connectivity[nodeKept, 1] = totalConnect

    # Replace nodeGone with nodeKept in links and update linksConnect with the new positions of the links in connectivity.
    @inbounds @simd for i in 1:nodeGoneConnect
        link = connectivity[nodeGone, 2 * i]      # Link id for nodeGone
        colLink = connectivity[nodeGone, 2 * i + 1] # Position of links where link appears.
        links[link, colLink] = nodeKept         # Replace nodeGone with nodeKept.
        linksConnect[link, colLink] = nodeKeptConnect + i # Increase connectivity of nodeKept.
    end

    # Remove node from network and update index of nodeKept in case it changed.
    lastNode = maximum((network.numNode, 1))
    removeNode!(network, nodeGone, lastNode)
    nodeKept == lastNode ? nodeKept = nodeGone : nothing

    # Warning: connectivity can be updated by removeLink!().
    # Delete self links of nodeKept.
    i = 1
    @inbounds while i < connectivity[nodeKept, 1]
        link = connectivity[nodeKept, 2 * i]      # Link id for nodeKept
        colLink = connectivity[nodeKept, 2 * i + 1] # Position of links where link appears.
        colNotLink = 3 - colLink # The column of the other node in the position of link in links.
        nodeNotLink = links[link, colNotLink]   # Other node in the link.
        # If nodeKept and nodeNotLink are the same, this is a self-link, therefore delete. Else increase counter.
        nodeKept == nodeNotLink ? removeLink!(network, link) : i += 1
    end

    # Warning: connectivity can be updated in the while loop.
    # Delete duplicate links.
    i = 1
    @inbounds while i < connectivity[nodeKept, 1]
        link1 = connectivity[nodeKept, 2 * i]         # Link id for nodeKept
        colLink1 = connectivity[nodeKept, 2 * i + 1]    # Position of links where link appears.
        colNotLink1 = 3 - colLink1 # The column of the other node in the position of link in links.
        nodeNotLink1 = links[link1, colNotLink1]    # Other node in the link.
        # Same but for the next link.
        j = i + 1
        while j <= connectivity[nodeKept, 1]
            link2 = connectivity[nodeKept, 2 * j]
            colLink2 = connectivity[nodeKept, 2 * j + 1]
            colNotLink2 = 3 - colLink2
            nodeNotLink2 = links[link2, colNotLink2]

            # Continue to next iteration if no duplicate links are found.
            if nodeNotLink1 != nodeNotLink2
                j += 1
                continue
            end

            # Fix Burgers vector.
            # If the nodes are on different ends of a link, conservation of Burgers vector in a loop requires we subtract contributions from the two links involved. Else we add them.
            colNotLink1 + colNotLink2 == 3 ? bVec[link1, :] -= bVec[link2, :] :
            bVec[link1, :] += bVec[link2, :]

            # WARNING This calculation is odd. Try using the cross product of the adjacent segments.
            # Fix slip plane.
            # Line direction and velocity of the resultant dislocation.
            t = @SVector [
                coord[nodeKept, 1] - coord[nodeNotLink1, 1],
                coord[nodeKept, 2] - coord[nodeNotLink1, 2],
                coord[nodeKept, 3] - coord[nodeNotLink1, 3],
            ]

            v = @SVector [
                nodeVel[nodeKept, 1] + nodeVel[nodeNotLink1, 1],
                nodeVel[nodeKept, 2] + nodeVel[nodeNotLink1, 2],
                nodeVel[nodeKept, 3] + nodeVel[nodeNotLink1, 3],
            ]

            # Burgers vector and potential new slip plane.
            b = @SVector [bVec[link1, 1], bVec[link1, 2], bVec[link1, 3]]
            n1 = t × b  # For non-screw segments.
            n2 = t × v  # For screw segments.
            if n1 ⋅ n1 > eps(eltype(n1)) # non-screw
                slipPlane[link1, :] = n1 / norm(n1)
            elseif n2 ⋅ n2 > eps(eltype(n2)) # screw
                slipPlane[link1, :] = n2 / norm(n2)
            end

            # Remove link2 from network and update the index of link1 in case it changed.
            lastLink = maximum((network.numSeg, 1))
            removeLink!(network, link2, lastLink)
            link1 == lastLink ? link1 = link2 : nothing

            # If the burgers vector of the new junction is non-zero, continue to the next iteration. Else remove it.
            b = @SVector [bVec[link1, 1], bVec[link1, 2], bVec[link1, 3]]
            if isapprox(dot(b, b), 0)
                removeLink!(network, link1)
                # If the node that was connected to nodeKept has no connections, remove it from the network and update the index of nodeKept in case it changed.
                if connectivity[nodeNotLink1, 1] == 0
                    lastNode = maximum((network.numNode, 1))
                    removeNode!(network, nodeNotLink1, lastNode)
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
    if connectivity[nodeKept, 1] == 0
        removeNode!(network, nodeKept)
        nodeKept = 0
    end

    return nodeKept
end

function splitNode end

function coarsenNetwork!(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    # mesh::RegularCuboidMesh,
    # dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
)
    minAreaSq = dlnParams.minArea^2
    minSegLen = dlnParams.minSegLen
    maxSegLen = dlnParams.maxSegLen

    label = network.label
    links = network.links
    connectivity = network.connectivity
    linksConnect = network.linksConnect
    segForce = network.segForce
    nodeVel = network.nodeVel

    i = 1
    while i <= network.numNode
        # We only want to coarsen real internal nodes. Else we skip to the next node.
        if !(connectivity[i, 1] == 2 && label[i] == 1)
            i += 1
            continue
        end
        link1 = connectivity[i, 2]  # Node i in link 1
        link2 = connectivity[i, 4]  # Node i in link 2

        colInLink1 = connectivity[i, 3] # Column of node i in link 1
        colInLink2 = connectivity[i, 5] # Column of node i in link 2

        oppColLink1 = 3 - colInLink1 # Node i is connected via link 1 to the node that is in this column in links.
        oppColLink2 = 3 - colInLink2 # Node i is connected via link 2 to the node that is in this column in links.

        link1_nodeOppI = links[link1, oppColLink1] # Node i is connected to this node as part of link 1.
        link2_nodeOppI = links[link2, oppColLink2] # Node i is connected to this node as part of link 2.

        # We don't want to remesh out segments between two fixed nodes because the nodes by definition do not move and act as a source, thus we skip to the next node.
        if label[link1_nodeOppI] == 2 && label[link2_nodeOppI] == 2
            i += 1
            continue
        end

        # Coordinate of node i
        iCoord = @SVector [coord[i, 1], coord[i, 2], coord[i, 3]]
        # Create a triangle formed by the three nodes involved in coarsening.
        coordVec1 =
            SVector(
                coord[link1_nodeOppI, 1],
                coord[link1_nodeOppI, 2],
                coord[link1_nodeOppI, 3],
            ) - iCoord # Vector between node 1 and the node it's connected to via link 1.
        coordVec2 =
            SVector(
                coord[link2_nodeOppI, 1],
                coord[link2_nodeOppI, 2],
                coord[link2_nodeOppI, 3],
            ) - iCoord # Vector between node 1 and the node it's connected to via link 2.
        coordVec3 = vec2 - vec1 # Vector between both nodes connected to iCoord.
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
        iVel = @SVector [nodeVel[i, 1], nodeVel[i, 2], nodeVel[i, 3]]
        velVec1 =
            SVector(
                nodeVel[link1_nodeOppI, 1],
                nodeVel[link1_nodeOppI, 2],
                nodeVel[link1_nodeOppI, 3],
            ) - iVel
        velVec2 =
            SVector(
                nodeVel[link2_nodeOppI, 1],
                nodeVel[link2_nodeOppI, 2],
                nodeVel[link2_nodeOppI, 3],
            ) - iVel
        velVec3 = velVec2 - velVec1

        # Rate of change of side length with respect to time. We add eps(typeof(r)) to avoid division by zero.
        dr1dt = coordVec1 ⋅ velVec1 / (r1 + eps(typeof(r1)))
        dr2dt = coordVec2 ⋅ velVec2 / (r2 + eps(typeof(r2)))
        dr3dt = coordVec3 ⋅ velVec3 / (r3 + eps(typeof(r3)))
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

        # If link2_nodeOppI no longer exists there is nothing to calculate and we proceed to the next iteration.
        nodeMerged == 0 ? continue : nothing

        for j in 1:connectivity[nodeMerged, 1]
            # Find the new link that has been created between nodeMerged and nodeNotMerged.
            linkMerged = connectivity[nodeMerged, 2 * j]
            colLinkMerged = connectivity[nodeMerged, 2 * j + 1]
            colNotLinkMerged = 3 - colLinkMerged
            nodeNotMerged = links[linkMerged, colNotLinkMerged]

            if nodeNotMerged == i || nodeNotMerged == link1_nodeOppI
                #remesh 67
                getSegmentIdx!(network)
                # Calculate segment force for segment linkMerged.
                for k in 1:2
                    # Calculate new velocity for the two nodes involved in linkMerged.
                end
            end
        end
    end
    return network
end

## Addition
function splitNode!(
    network::DislocationNetwork,
    splitNode::Int,
    splitConnect::Int,
    midCoord::AbstractVector{T},
    midVel::AbstractVector{T},
) where {T}

    # newNode gets inserted between splitNode and the node it is connected to via the connection splitConnect. We want to take this connection and remove it from splitNode. We then connect splitNode to newNode. Then we assign splitConnect to newNode so that it can connect to the node splitNode used to be connected to, this way we close the loop and the connections move from splitNode -> other, to splitNode -> newNode -> other. We copy the connection to a temporary variable, guaranteeing the data isn't modified by removeConnection!().
    tmpConnect = copy(network.connectivity[splitNode, :])
    # Remove connection from node in preparation of adding the new node.
    removeConnection!(network, splitNode, splitConnect)

    # New node created at the very end.
    network.numNode += 1
    newNode = network.numNode

    # Allocate memory.
    if newNode > length(network.label)
        # A good heuristic for memory allocation is to allocate an extra N log₂(N) entries.
        numNewEntries = Int(round(newNode * log2(newNode)))
        network.coord = vcat(
            network.coord,
            zeros(eltype(network.coord), numNewEntries, size(network.coord, 2)),
        )
        network.label = vcat(network.label, zeros(eltype(network.label), numNewEntries))
        network.nodeVel = vcat(
            network.nodeVel,
            zeros(eltype(network.nodeVel), numNewEntries, size(network.nodeVel, 2)),
        )
        network.connectivity = vcat(
            network.connectivity,
            zeros(
                eltype(network.connectivity),
                numNewEntries,
                size(network.connectivity, 2),
            ),
        )
    end

    links = network.links
    bVec = network.bVec
    network.coord[newNode, :] = midCoord
    network.label[newNode] = nodeType(1)
    network.nodeVel[newNode, :] = midVel

    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # Create connectivity for the new node and update links and linksConnect.
    # This only does the first connection. The rest are made after in case they are needed.
    connectivity[newNode, 1] = 1
    connectivity[newNode, 2:3] = tmpConnect[(2 * splitConnect):(2 * splitConnect + 1)]
    link = connectivity[newNode, 2]
    colLink = connectivity[newNode, 3]
    links[link, colLink] = newNode
    linksConnect[link, colLink] = 1

    # Check if we need a new link between newNode and splitNode for Burgers vector conservation.
    b = @SVector [bVec[link, 1], bVec[link, 2], bVec[link, 3]]
    # If burgers vector is conserved return.
    b ⋅ b == 0 && return network

    # New link created at the end.
    network.numSeg += 1
    newSeg = network.numSeg

    # Allocate memory.
    if newSeg > size(network.links, 1)
        # A good heuristic for memory allocation is to allocate an extra N log₂(N) entries.
        numNewEntries = Int(round(newSeg * log2(newSeg)))
        network.links = vcat(
            network.links,
            zeros(eltype(network.links), numNewEntries, size(network.links, 2)),
        )
        network.slipPlane = vcat(
            network.slipPlane,
            zeros(eltype(network.slipPlane), numNewEntries, size(network.slipPlane, 2)),
        )
        network.bVec = vcat(
            network.bVec,
            zeros(eltype(network.bVec), numNewEntries, size(network.bVec, 2)),
        )
        network.segForce = vcat(
            network.segForce,
            zeros(eltype(network.segForce), numNewEntries, size(network.segForce, 2)),
        )
        network.linksConnect = vcat(
            network.linksConnect,
            zeros(
                eltype(network.linksConnect),
                numNewEntries,
                size(network.linksConnect, 2),
            ),
        )
    end

    bVec = network.bVec
    links = network.links
    slipPlane = network.slipPlane
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # Add new link to network.
    bVec[newSeg, :] = b
    links[newSeg, colLink] = splitNode  # First node in the link is the node that was split.
    colOppLink = 3 - colLink            # Other column in links.
    links[newSeg, colOppLink] = newNode # Second node of the link is the new node.

    # Update connectivity of splitNode.
    newConnect1 = connectivity[splitNode, 1] + 1
    connectivity[splitNode, 1] = newConnect1
    connectivity[splitNode, (2 * newConnect1):(2 * newConnect1 + 1)] = [newSeg, colLink]

    # Update connectivity of newNode.
    newConnect2 = connectivity[newNode, 1] + 1
    connectivity[newNode, 1] = newConnect2
    connectivity[newNode, (2 * newConnect2):(2 * newConnect2 + 1)] = [newSeg, colOppLink]

    # Update linksConnect.
    linksConnect[newSeg, colLink] = newConnect1
    linksConnect[newSeg, colOppLink] = newConnect2

    # Fix slip plane.
    # Line direction and velocity of the resultant dislocation.
    coord = network.coord
    nodeVel = network.nodeVel

    # WARNING This calculation is odd. Try using the cross product of the adjacent segments.
    t = @SVector [
        coord[splitNode, 1] - coord[newNode, 1],
        coord[splitNode, 2] - coord[newNode, 2],
        coord[splitNode, 3] - coord[newNode, 3],
    ]

    v = @SVector [
        nodeVel[splitNode, 1] + nodeVel[newNode, 1],
        nodeVel[splitNode, 2] + nodeVel[newNode, 2],
        nodeVel[splitNode, 3] + nodeVel[newNode, 3],
    ]

    # Potential new slip plane.
    n1 = t × b  # For non-screw segments.
    n2 = t × v  # For screw segments.
    if n1 ⋅ n1 > eps(eltype(n1)) # non-screw
        slipPlane[newSeg, :] = n1 / norm(n1)
    elseif n2 ⋅ n2 > eps(eltype(n2)) # screw
        slipPlane[newSeg, :] = n2 / norm(n2)
    end

    return network
end

function refineNetwork!(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    # mesh::RegularCuboidMesh,
    # dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
)

    maxAreaSq = dlnParams.maxArea^2
    maxSegLen = dlnParams.maxSegLen
    twoMinSegLen = 2 * dlnParams.minSegLen

    label = network.label
    links = network.links
    connectivity = network.connectivity
    linksConnect = network.linksConnect
    segForce = network.segForce
    nodeVel = network.nodeVel
    numNode = network.numNode

    for i in 1:numNode
        if connectivity[i, 1] == 2 && label[i] == 1
            link1 = connectivity[i, 2]  # First connection.
            link2 = connectivity[i, 4]  # Second connection.
            colLink1 = connectivity[i, 3]   # Column where node i is in links of the first connection.
            colLink2 = connectivity[i, 5]   # Column where node i is in links of the second connection.
            oppColLink1 = 3 - colInLink1 # Node i is connected via link 1 to the node that is in this column in links.
            oppColLink2 = 3 - colInLink2 # Node i is connected via link 2 to the node that is in this column in links.
            link1_nodeOppI = links[link1, oppColLink1] # Node i is connected to this node as part of link 1.
            link2_nodeOppI = links[link1, oppColLink2] # Node i is connected to this node as part of link 2.

            # Create triangle formed by the node and its two links.
            iCoord = @SVector [coord[i, 1], coord[i, 2], coord[i, 3]]
            # Side 1
            coordVec1 =
                SVector(
                    coord[link1_nodeOppI, 1],
                    coord[link1_nodeOppI, 2],
                    coord[link1_nodeOppI, 3],
                ) - iCoord
            # Side 2
            coordVec2 =
                SVector(
                    coord[link2_nodeOppI, 1],
                    coord[link2_nodeOppI, 2],
                    coord[link2_nodeOppI, 3],
                ) - iCoord
            # Side 3 (close the triangle)
            coordVec3 = coordVec2 - coordVec1
            # Side lengths.
            r1 = norm(coordVec1)
            r2 = norm(coordVec2)
            r3 = norm(coordVec3)

            # Half the triangle's perimeter.
            r0 = (r1 + r2 + r3) / 2
            # Heron's formula for the area of a triangle when we know its sides.
            areaSq = r0 * (r0 - r1) * (r0 - r2) * (r0 - r3)

            # Check if we have to split the first link.
            if (areaSq > maxAreaSq && r1 >= twoMinSegLen && link1_nodeOppI <= numNode) ||
               r1 > maxSegLen
                midCoord =
                    (
                        SVector(coord[i, 1], coord[i, 2], coord[i, 3]) + SVector(
                            coord[link1_nodeOppI, 1],
                            coord[link1_nodeOppI, 2],
                            coord[link1_nodeOppI, 3],
                        )
                    ) / 2

                midVel =
                    (
                        SVector(nodeVel[i, 1], nodeVel[i, 2], nodeVel[i, 3]) + SVector(
                            nodeVel[link1_nodeOppI, 1],
                            nodeVel[link1_nodeOppI, 2],
                            nodeVel[link1_nodeOppI, 3],
                        )
                    ) / 2

                splitNode!(network, i, 1, midCoord, midVel)

                newNode = network.numNode
                newLink = network.numSeg

                slipPlane[link1, :] == slipPlane[link2, :] ?
                slipPlane[newlink, :] = slipPlane[link1, :] : nothing

                for j in 1:connectivity[newNode, 1]
                    link = connectivity[newNode, 2 * j]
                    colLink = connectivity[newNode, 2 * j + 1]
                    oppColLink = 3 - colLink
                    oldNode = links[link, oppColLink]
                    # remesh 154
                end
                # remesh 159
            end

            # Check if we have to split the second link.
            if (areaSq > maxAreaSq && r2 >= twoMinSegLen && link2_nodeOppI <= numNode) ||
               r2 > maxSegLen
                midCoord =
                    (
                        SVector(coord[i, 1], coord[i, 2], coord[i, 3]) + SVector(
                            coord[link2_nodeOppI, 1],
                            coord[link2_nodeOppI, 2],
                            coord[link2_nodeOppI, 3],
                        )
                    ) / 2

                midVel =
                    (
                        SVector(nodeVel[i, 1], nodeVel[i, 2], nodeVel[i, 3]) + SVector(
                            nodeVel[link2_nodeOppI, 1],
                            nodeVel[link2_nodeOppI, 2],
                            nodeVel[link2_nodeOppI, 3],
                        )
                    ) / 2

                splitNode!(network, i, 2, midCoord, midVel)

                newNode = network.numNode
                newLink = network.numSeg

                slipPlane[link2, :] == slipPlane[link1, :] ?
                slipPlane[newlink, :] = slipPlane[link2, :] : nothing

                for j in 1:connectivity[newNode, 1]
                    link = connectivity[newNode, 2 * j]
                    colLink = connectivity[newNode, 2 * j + 1]
                    oppColLink = 3 - colLink
                    oldNode = links[link, oppColLink]
                    # remesh 133
                end
                # remesh 139
            end

        elseif connectivity[i, 1] > 2 && label[i] == 1
            # remesh 163
        end
    end

    return network
end
#=
end

function remeshNetwork(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    mesh::RegularCuboidMesh,
    dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
) end

function remeshInternal(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    mesh::RegularCuboidMesh,
    dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
)

end

function remeshSurface(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    mesh::RegularCuboidMesh,
    dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
)

end
=#
