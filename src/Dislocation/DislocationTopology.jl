"""
Replaces node id with the last valid node. Cleans up links and
"""
function removeNode!(network::DislocationNetwork, nodeKept::Int, nodeGone::Int)
    links = network.links
    coord = network.coord
    label = network.label
    nodeVel = network.nodeVel
    connectivity = network.connectivity

    # Find the first undefined node.
    # If there are no undefined nodes, the index of the last defined node is the length of label, else the last defined node is one before the first undefined one.
    firstUndef = findfirst(x -> x == 0, label)
    firstUndef == nothing ? lastNode = length(label) : lastNode = firstUndef - 1

    # The node that is gone is replaced by the last node.
    if nodeGone < lastNode
        coord[nodeGone, :] .= coord[lastNode, :]
        label[nodeGone] .= label[lastNode]
        nodeVel[nodeGone, :] .= nodeVel[lastNode, :]
        connectivity[nodeGone, :] .= connectivity[lastNode, :]
        # Change the link
        for j in 1:connectivity[nodeGone, 1]
            links[connectivity[nodeGone, 2 * j], connectivity[nodeGone, 2 * j + 1]] = nodeGone
        end
    end

    # Zero out the last node.
    coord[lastNode, :] .= 0
    label[lastNode] = 0
    nodeVel[lastNode, :] .= 0
    connectivity[lastNode, :] .= 0
    network.numNode -= 1

    # If the remaining node was the last node, change its index to the removed node.
    nodeKept == lastNode ? nodeKept = nodeGone : nothing

    return network, nodeKept
end

"""
```
removeLink!(network::DislocationNetwork, link::Int)
```
Removes link and its information from network.
"""
function removeLink!(network::DislocationNetwork, link::Int)
    links = network.links
    coord = network.coord
    label = network.label
    nodeVel = network.nodeVel
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    @assert link > size(links, 1) "link $link not found."

    # Delete linkid from connectivity.
    node1 = links[link, 1]
    connectGone1 = linksConnect[link, 1]
    # mergenodes 142


end

"""
Merges and cleans up the information in `network.connectivity` and `network.links` for the nodes that will be merged. This is such that there are no repeated entries, self-links or double links.
"""
function mergeNode!(network::DislocationNetwork, nodeKept::Int, nodeGone::Int)

    # Return if both nodes to be merged are the same.
    nodeKept == nodeGone && return

    links = network.links
    coord = network.coord
    segForce = network.segForce
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    nodeKeptConnect = connectivity[nodeKept, 1]
    nodeGoneConnect = connectivity[nodeGone, 1]
    totalConnect = nodeKeptConnect + nodeGoneConnect

    # Pass connections from nodeGone to nodeKept.
    connectivity[nodeKept, (2 * (nodeKeptConnect + 1)):(2 * (nodeGoneConnect + 1))] =
        connectivity[nodeGone, 2:(2 * nodeGoneConnect + 1)]
    connectivity[nodeKeptConnect, 1] = totalConnect

    # Replace nodeGone with nodeKept in links and update linksConnect with the new positions of the links in connectivity.
    for i in 1:nodeGoneConnect
        link = connectivity[nodeGone, 2 * i]        # Link id for nodeGone
        colLink = connectivity[nodeGone, 2 * i + 1] # Position of links where link appears.
        links[link, colLink] = nodeKept             # Replace nodeGone with nodeKept.
        linksConnect[link, colLink] = nodeKeptConnect + i # Increase connectivity of nodeKept.
    end

    # Remove node and update nodeKept.
    network, nodeKept = removeNode!(network, nodeGone)

    # Delete self links of nodeKept.
    for i in 1:nodeKeptConnect - 1
        link = connectivity[nodeKept, 2 * i]        # Link id for nodeKept
        colLink = connectivity[nodeKept, 2 * i + 1] # Position of links where link appears.
        colNotLink = 3 - colLink # The column of the other node in the position of link in links.
        nodeNotLink = links[link, colNotLink]       # Other node in the link.

        nodeKept == nodeNotLink ? removeLink!(network, link) : continue
        i -= 1
    end
    # mergenodes 46

end



function splitNode end

function coarsenMesh(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    mesh::RegularCuboidMesh,
    dlnFEM::DislocationFEMCorrective;
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
