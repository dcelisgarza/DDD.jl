"""
```
removeNode!(network::DislocationNetwork, nodeGone, lastNode = nothing)
```
Removes node `nodeGone` from network.
"""
function removeNode!(network::DislocationNetwork, nodeGone, lastNode = nothing)
    links = network.links
    coord = network.coord
    label = network.label
    nodeVel = network.nodeVel
    nodeForce = network.nodeForce
    connectivity = network.connectivity

    isnothing(lastNode) ? lastNode = maximum((network.numNode[1], 1)) : nothing

    # The if nodeGone is not the last node in the relevant arrays, replace it by lastNode.
    if nodeGone < lastNode
        label[nodeGone] = label[lastNode]
        @inbounds @simd for i in 1:3
            coord[i, nodeGone] = coord[i, lastNode]
            nodeVel[i, nodeGone] = nodeVel[i, lastNode]
            nodeForce[i, nodeGone] = nodeForce[i, lastNode]
        end
        @inbounds @simd for i in 1:size(connectivity, 1)
            connectivity[i, nodeGone] = connectivity[i, lastNode]
        end
        # Change the link
        @inbounds @simd for i in 1:connectivity[1, nodeGone]
            idx = 2 * i
            links[connectivity[idx + 1, nodeGone], connectivity[idx, nodeGone]] = nodeGone
        end
    end

    # Remove the last node from the arrays since nodeGone was either replaced by lastNode above, or already was the last node. Reduce the number of nodes by one.
    coord[:, lastNode] .= 0
    label[lastNode] = 0
    nodeVel[:, lastNode] .= 0
    nodeForce[:, lastNode] .= 0
    network.numNode[1] -= 1
    connectivity[:, lastNode] .= 0

    return network
end

"""
```
removeConnection!(network::DislocationNetwork, nodeKept, connectGone)
```
Removes connection `connectGone` from `nodeKept`.
"""
function removeConnection!(network::DislocationNetwork, nodeKept, connectGone)
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
        connectivity[idx, nodeKept] = connectivity[lst, nodeKept]
        connectivity[idx + 1, nodeKept] = connectivity[lst + 1, nodeKept]
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
removeLink!(network::DislocationNetwork, linkGone, lastLink = nothing)
```
Removes link `linkGone` and uses `lastLink` to reoganise the network.
"""
function removeLink!(network::DislocationNetwork, linkGone, lastLink = nothing)
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
    checkLink = SVector{2, Int}(linksConnect[1, linkGone], linksConnect[2, linkGone])
    @assert checkLink ⋅ checkLink == 0 "removeLink!: link $linkGone still has connections and should not be deleted."

    isnothing(lastLink) ? lastLink = maximum((network.numSeg[1], 1)) : nothing

    # If the linkGone is not the last link, replace it with lastLink.
    if linkGone < lastLink
        @inbounds @simd for i in 1:2
            links[i, linkGone] = links[i, lastLink]
            linksConnect[i, linkGone] = linksConnect[i, lastLink]
        end
        @inbounds @simd for i in 1:3
            slipPlane[i, linkGone] = slipPlane[i, lastLink]
            bVec[i, linkGone] = bVec[i, lastLink]
        end
        @inbounds for j in 1:2
            @simd for i in 1:3
                segForce[i, j, linkGone] = segForce[i, j, lastLink]
            end
        end
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
    network.numSeg[1] -= 1
    linksConnect[:, lastLink] .= 0

    return network
end

"""
```
mergeNode!(network::DislocationNetwork, nodeKept, nodeGone)
```
Merges `nodeGone` into `nodeKept`.
"""
function mergeNode!(network::DislocationNetwork, nodeKept, nodeGone)
    @assert nodeKept <= network.numNode[1] && nodeGone <= network.numNode[1] "mergeNode: the node kept after merging, $nodeKept and node removed after merging, $nodeGone, must be in the simulation."

    # Return if both nodes to be merged are the same.
    nodeKept == nodeGone && return 0, network

    elemT = eltype(network.bVec)

    nodeKeptConnect = network.connectivity[1, nodeKept]
    nodeGoneConnect = network.connectivity[1, nodeGone]
    totalConnect = nodeKeptConnect + nodeGoneConnect

    # Pass connections from nodeGone to nodeKept.
    if size(network.connectivity, 1) < 2 * totalConnect + 1
        network = DislocationNetwork(;
            links = network.links,
            slipPlane = network.slipPlane,
            bVec = network.bVec,
            coord = network.coord,
            label = network.label,
            segForce = network.segForce,
            nodeVel = network.nodeVel,
            nodeForce = network.nodeForce,
            numNode = network.numNode,
            numSeg = network.numSeg,
            maxConnect = network.maxConnect,
            linksConnect = network.linksConnect,
            connectivity = vcat(
                network.connectivity,
                zeros(
                    Int,
                    2 * totalConnect + 1 - size(network.connectivity, 1),
                    length(network.label),
                ),
            ),
            segIdx = network.segIdx,
        )
    end

    links = network.links
    slipPlane = network.slipPlane
    bVec = network.bVec
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # This next loop is equivalent to the following expression but doesn't allocate memory.
    # connectivity[(2 * (nodeKeptConnect + 1)):(2 * totalConnect + 1), nodeKept] = connectivity[2:(2 * nodeGoneConnect + 1), nodeGone]
    idx = (2 * (nodeKeptConnect + 1)):(2 * totalConnect + 1)
    lhs = idx[1]
    @inbounds @simd for i in 0:(length(idx) - 1)
        connectivity[lhs + i, nodeKept] = connectivity[2 + i, nodeGone]
    end
    connectivity[1, nodeKept] = totalConnect

    # Replace nodeGone with nodeKept in links and update linksConnect with the new positions of the links in connectivity.
    @inbounds @simd for i in 1:nodeGoneConnect
        link = connectivity[2 * i, nodeGone]      # Link id for nodeGone
        colLink = connectivity[2 * i + 1, nodeGone] # Position of links where link appears.
        links[colLink, link] = nodeKept         # Replace nodeGone with nodeKept.
        linksConnect[colLink, link] = nodeKeptConnect + i # Increase connectivity of nodeKept.
    end

    # Remove node from network and update index of nodeKept in case it changed.
    lastNode = maximum((network.numNode[1], 1))
    removeNode!(network, nodeGone, lastNode)
    coord = network.coord
    nodeVel = network.nodeVel
    connectivity = network.connectivity
    nodeKept == lastNode ? nodeKept = nodeGone : nothing

    # Warning: connectivity can be updated by removeLink!().
    # Delete self links of nodeKept.
    i = 1
    @inbounds while i < connectivity[1, nodeKept]
        link, missing, nodeNotLink = findConnectedNode(network, nodeKept, i) # i connection.

        # If nodeKept and nodeNotLink are the same, this is a self-link, therefore delete. Else increase counter.
        nodeKept == nodeNotLink ? removeLink!(network, link) : i += 1
        links = network.links
        connectivity = network.connectivity
    end

    # Warning: connectivity can be updated in the while loop.
    # Delete duplicate links.
    i = 1
    @inbounds while i < connectivity[1, nodeKept]
        link1, colNotLink1, nodeNotLink1 = findConnectedNode(network, nodeKept, i) # i connection.

        # Same but for the next link.
        j = i + 1
        while j <= connectivity[1, nodeKept]
            link2, colNotLink2, nodeNotLink2 = findConnectedNode(network, nodeKept, j) # i connection.

            # Continue to next iteration if no duplicate links are found.
            if nodeNotLink1 != nodeNotLink2
                j += 1
                continue
            end

            # Fix Burgers vector.
            # If the nodes are on different ends of a link, conservation of Burgers vector in a loop requires we subtract contributions from the two links involved. Else we add them.
            if colNotLink1 + colNotLink2 == 3
                @inbounds @simd for i in 1:3
                    bVec[i, link1] -= bVec[i, link2]
                end
            else
                @inbounds @simd for i in 1:3
                    bVec[i, link1] += bVec[i, link2]
                end
            end

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
            lastLink = maximum((network.numSeg[1], 1))
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
                    lastNode = maximum((network.numNode[1], 1))
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

    return nodeKept, network
end

"""
```
coarsenNetwork!(
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    mesh::AbstractMesh,
    forceDisplacement::ForceDisplacement,
    network::DislocationNetwork,
)
```
Coarsens network such that no links are smaller than the minimum allowable length and so that no two links form triangles with area under the minimum allowed.
"""
function coarsenNetwork!(
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    mesh::AbstractMesh,
    forceDisplacement::ForceDisplacement,
    network::DislocationNetwork,
)
    minAreaSq = dlnParams.minAreaSq
    minSegLen = dlnParams.minSegLen
    maxSegLen = dlnParams.maxSegLen

    label = network.label
    links = network.links
    coord = network.coord
    numNode = network.numNode[1]
    nodeVel = network.nodeVel
    connectivity = network.connectivity
    elemT = eltype(network.coord)

    i = 1
    @inbounds while i <= numNode
        # We only want to coarsen real internal nodes. Else we skip to the next node.
        if !(connectivity[1, i] == 2 && label[i] == intMobDln)
            i += 1
            continue
        end

        missing, missing, link1_nodeOppI = findConnectedNode(network, i, 1) # 1st connection.
        missing, missing, link2_nodeOppI = findConnectedNode(network, i, 2) # 2nd connection.

        # We don't want to remesh out segments between two fixed nodes because the nodes by definition do not move and act as a source, thus we skip to the next node.
        if label[link1_nodeOppI] == intFixDln && label[link2_nodeOppI] == intFixDln
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

        # If coarsening would result in a link whose length is bigger than the maximum allowed, and there r1 and r2 are bigger than the minimum allowed length, we skip to the next node.
        if r3 > maxSegLen && r1 > minSegLen && r2 > minSegLen
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
        nodeMerged, network = mergeNode!(network, link2_nodeOppI, i)
        getSegmentIdx!(network)
        links = network.links
        coord = network.coord
        label = network.label
        nodeVel = network.nodeVel
        numNode = network.numNode[1]
        connectivity = network.connectivity

        # If link2_nodeOppI no longer exists there is nothing to calculate and we proceed to the next iteration.
        if nodeMerged == 0
            i += 1
            continue
        end

        @simd for j in 1:connectivity[1, nodeMerged]
            # Find the new link that has been created between nodeMerged and nodeNotMerged.
            linkMerged = connectivity[2 * j, nodeMerged]
            colLinkMerged = connectivity[2 * j + 1, nodeMerged]
            colNotLinkMerged = 3 - colLinkMerged
            nodeNotMerged = links[colNotLinkMerged, linkMerged]

            if nodeNotMerged == i || nodeNotMerged == link1_nodeOppI
                # Calculate segment force for segment linkMerged.
                calcSegForce!(
                    dlnParams,
                    matParams,
                    mesh,
                    forceDisplacement,
                    network,
                    linkMerged,
                )
                # Calculate node velocity.
                nodes = SVector{2, Int}(links[1, linkMerged], links[2, linkMerged])
                dlnMobility!(dlnParams, matParams, network, nodes)
                nodeVel = network.nodeVel
            end
        end
    end
    return network
end
