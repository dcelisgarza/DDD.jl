"""
```
removeNode!(network::DislocationNetwork, nodeGone::Int, lastNode = missing)
```
In-place remove `nodeGone`. If `nodeGone` is *not* the last node, this replaces the entries corresponding to `nodeGone` with `lastNode`. Else it simply zeros out the entries corresponding to `lastNode`. This also decreases `numNode` by one. If `lastNode` is missing, this function finds it.

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
@inline function removeNode!(network::DislocationNetwork, nodeGone::Int, lastNode = missing)
    links = network.links
    coord = network.coord
    label = network.label
    nodeVel = network.nodeVel
    connectivity = network.connectivity

    if ismissing(lastNode)
        lastNode = maximum((network.numNode, 1))
    end

    # The if nodeGone is not the last node in the relevant arrays, replace it by lastNode.
    if nodeGone < lastNode
        coord[nodeGone, :] .= coord[lastNode, :]
        label[nodeGone] = label[lastNode]
        nodeVel[nodeGone, :] .= nodeVel[lastNode, :]
        connectivity[nodeGone, :] .= connectivity[lastNode, :]
        # Change the link
        @simd for j = 1:connectivity[nodeGone, 1]
            idx = 2 * j
            links[connectivity[nodeGone, idx], connectivity[nodeGone, idx+1]] = nodeGone
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
@inline function removeConnection!(network::DislocationNetwork, nodeKept::Int, connectGone::Int)
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # The last connection to be added to nodeKept is equal to the connectivity of the nodeKept.
    lastConnect = connectivity[nodeKept, 1]
    lst = 2 * lastConnect

    # Remove entry from linksConnect.
    idx = 2 * connectGone
    link1 = connectivity[nodeKept, idx]
    link2 = connectivity[nodeKept, idx+1]
    linksConnect[link1, link2] = 0

    # If the connectGone is not the last connection that was made, replace it by the last connection.
    if connectGone < lastConnect
        connectivity[nodeKept, idx:(idx+1)] = connectivity[nodeKept, lst:(lst+1)]
        linksConnect[connectivity[nodeKept, idx], connectivity[nodeKept, idx+1]] = connectGone
    end

    # Change connectivity to reflect that nodeKept has one less connection.
    # Remove the last connection since lastConnect was either replaced by connectGone above, or already was the last connection.
    connectivity[nodeKept, 1] -= 1
    connectivity[nodeKept, lst:(lst+1)] .= 0

    return network
end

"""
```
removeLink!(network::DislocationNetwork, linkGone::Int, lastLink = missing)
```
In-place remove `linkGone` with the last valid link. If `linkGone` is *not* the last link, this replaces the entries corresponding to `linkGone` with `lastLink`. Else it simply zeros out the entries corresponding to `lastLink`. This also decreases `numSeg` by one. If `lastLink` is missing, this function finds it.

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
function removeLink!(network::DislocationNetwork, linkGone::Int, lastLink = missing)
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

    if ismissing(lastLink)
        lastLink = maximum((network.numSeg, 1))
    end

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
@inline function mergeNode!(network::DislocationNetwork, nodeKept::Int, nodeGone::Int)
    # Return if both nodes to be merged are the same.
    nodeKept == nodeGone && return

    @assert nodeKept <= network.numNode && nodeGone <= network.numNode "mergeNode!: the node kept after merging, $nodeKept and node removed after merging, $nodeGone, must be in the simulation."

    links = network.links
    slipPlane = network.slipPlane
    bVec = network.bVec
    coord = network.coord
    label = network.label
    segForce = network.segForce
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
    connectivity[nodeKept, (2*(nodeKeptConnect+1)):(2*totalConnect+1)] =
        connectivity[nodeGone, 2:(2*nodeGoneConnect+1)]
    connectivity[nodeKept, 1] = totalConnect

    # Replace nodeGone with nodeKept in links and update linksConnect with the new positions of the links in connectivity.
    for i = 1:nodeGoneConnect
        link = connectivity[nodeGone, 2*i]        # Link id for nodeGone
        colLink = connectivity[nodeGone, 2*i+1] # Position of links where link appears.
        links[link, colLink] = nodeKept             # Replace nodeGone with nodeKept.
        linksConnect[link, colLink] = nodeKeptConnect + i # Increase connectivity of nodeKept.
    end

    # Remove node from network and update index of nodeKept in case it changed.
    lastNode = maximum((network.numNode, 1))
    removeNode!(network, nodeGone, lastNode)
    nodeKept == lastNode ? nodeKept = nodeGone : nothing

    # Warning: connectivity can be updated by removeLink!().
    # Delete self links of nodeKept.
    i = 1
    while i < connectivity[nodeKept, 1]
        link = connectivity[nodeKept, 2*i]      # Link id for nodeKept
        colLink = connectivity[nodeKept, 2*i+1] # Position of links where link appears.
        colNotLink = 3 - colLink # The column of the other node in the position of link in links.
        nodeNotLink = links[link, colNotLink]   # Other node in the link.
        # If nodeKept and nodeNotLink are the same, this is a self-link, therefore delete. Else increase counter.
        nodeKept == nodeNotLink ? removeLink!(network, link) : i += 1
    end

    # Warning: connectivity can be updated in the while loop.
    # Delete duplicate links.
    i = 1
    while i < connectivity[nodeKept, 1]
        link1 = connectivity[nodeKept, 2*i]         # Link id for nodeKept
        colLink1 = connectivity[nodeKept, 2*i+1]    # Position of links where link appears.
        colNotLink1 = 3 - colLink1 # The column of the other node in the position of link in links.
        nodeNotLink1 = links[link1, colNotLink1]    # Other node in the link.
        # Same but for the next link.
        j = i + 1
        while j <= connectivity[nodeKept, 1]
            link2 = connectivity[nodeKept, 2*j]
            colLink2 = connectivity[nodeKept, 2*j+1]
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
            if dot(n1, n1) > eps(eltype(n1)) # non-screw
                slipPlane[link1, :] = n1 / norm(n1)
            elseif dot(n2, n2) > eps(eltype(n2)) # screw
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
    connectivity[nodeKept, 1] == 0 ? removeNode!(network, nodeKept) : nothing

    return network
end

function splitNode end

function coarsenMesh(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    # mesh::RegularCuboidMesh,
    # dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
)
    minArea = dlnParams.minArea^2
    minSegLen = dlnParams.minSegLen
    maxSegLen = dlnParams.maxSegLen
    tiny = eps(Float64)

    label = network.label
    idx = findall(x -> x == 1, label)
    links = network.links
    connectivity = network.links
    linksConnect = network.linksConnect
    segForce = network.segForce
    nodeVel = network.nodeVel

    for i in idx
        connectivity[i, 1] != 2 ? continue : nothing
        link1 = connectivity[i, 2]  # Node i in link 1
        link2 = connectivity[i, 4]  # Node i in link 2

        colInLink1 = connectivity[i, 3] # Column of node i in link 1
        colInLink2 = connectivity[i, 5] # Column of node i in link 2

        colNotInLink1 = 3 - colInLink1 # Node i is connected via link 1 to the node that is in this column in links.
        colNotInLink2 = 3 - colInLink2 # Node i is connected via link 2 to the node that is in this column in links.

        link1_nodeNotInLink = links[i, colNotInLink1] # Node i is connected to this node as part of link 1.
        link2_nodeNotInLink = links[i, colNotInLink2] # Node i is connected to this node as part of link 2.

        # We don't want to remesh out segments between two fixed nodes because the nodes by definition do not move and act as a source.
        label[link1_nodeNotInLink] == 2 && label[link2_nodeNotInLink] == 2 ? continue : nothing

        # Coordinate of node i
        iCoord = coord[i, :]
        # Create a triangle formed by the three nodes involved in coarsening.
        coordVec1 = coord[link1_nodeNotInLink, :] - iCoord # Vector between node 1 and the node it's connected to via link 1.
        coordVec2 = coord[link2_nodeNotInLink, :] - iCoord # Vector between node 1 and the node it's connected to via link 2.
        coordVec3 = vec2 - vec1 # Vector between both nodes connected to iCoord.
        # Lengths of the triangle sides.
        r1 = norm(coordVec1)
        r2 = norm(coordVec2)
        r3 = norm(coordVec3)

        # If coarsening would result in a link whose length is bigger than the maximum we don't coarsen.
        r3 >= maxSegLen ? continue : nothing

        # Heron's formula for the area of a triangle when we know its sides.
        # Half the triangle's perimeter.
        r0 = (r1 + r2 + r3) / 2
        areaSq = r0 * (r0 - r1) * (r0 - r2) * (r0 - r3)

        # Node i velocities.
        iVel = nodeVel[i, :]
        velVec1 = nodeVel[link1_nodeNotInLink, :] - iVel
        velVec2 = nodeVel[link2_nodeNotInLink, :] - iVel
        velVec3 = velVec2 - velVec1

        # Rate of change of side length with respect to time. We add eps(Float64) to avoid division by zero.
        dr1dt = dot(coordVec1, velVec1) / (r1 + tiny)
        dr2dt = dot(coordVec2, velVec2) / (r2 + tiny)
        dr3dt = dot(coordVec3, velVec3) / (r3 + tiny)
        dr0dt = (dr1dt + dr2dt + dr3dt) / 2

        # Rate of change of triangle area with respect to time.
        dareaSqdt = dr0dt * (r0 - r1) * (r0 - r2) * (r0 - r3)
        dareaSqdt += r0 * (dr0dt - dr1dt) * (r0 - r2) * (r0 - r3)
        dareaSqdt += r0 * (r0 - r1) * (dr0dt - dr2dt) * (r0 - r3)
        dareaSqdt += r0 * (r0 - r1) * (r0 - r2) * (dr0dt - dr3dt)

        # Coarsen if the area is less than the minimum allowed area and the area is shrinking. Or if the length of link 1 or link 2 less than the minimum area. Else do continue to the next iteration.
        areaSq < minArea && dareaSqdt < 0 || r1 < minSegLen || r2 < minSegLen ? nothing : continue

        # remesh 60
    end

end

#=
function refineMesh(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    mesh::RegularCuboidMesh,
    dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
)

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
