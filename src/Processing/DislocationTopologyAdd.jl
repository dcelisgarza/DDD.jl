## Addition
function splitNode!(
    network::T1,
    splitNode::T2,
    splitConnect::T2,
    midCoord::T3,
    midVel::T3,
) where {T1 <: DislocationNetwork,T2 <: Int,T3 <: AbstractVector{T} where {T}}

    elemT = eltype(network.bVec)
    # newNode gets inserted between splitNode and the node it is connected to via the connection splitConnect. We want to take this connection and remove it from splitNode. We then connect splitNode to newNode. Then we assign splitConnect to newNode so that it can connect to the node splitNode used to be connected to, this way we close the loop and the connections move from splitNode -> other, to splitNode -> newNode -> other. We copy the connection to a temporary variable, guaranteeing the data isn't modified by removeConnection!().
    tmpConnect = network.connectivity[:, splitNode]
    # Remove connection from node in preparation of adding the new node.
    removeConnection!(network, splitNode, splitConnect)

    # New node created at the very end.
    network.numNode[1] += 1
    newNode = network.numNode[1]

    # Allocate memory.
    if newNode > length(network.label)
        # A good heuristic for memory allocation is to allocate an extra N log₂(N) entries.
        numNewEntries = Int(round(newNode * log2(newNode)))
        network = DislocationNetwork(;
            links = network.links,
            slipPlane = network.slipPlane,
            bVec = network.bVec,
            coord = hcat(
                network.coord,
                zeros(elemT, size(network.coord, 1), numNewEntries),
            ),
            label = vcat(network.label, zeros(nodeType, numNewEntries)),
            nodeVel = hcat(
                network.nodeVel,
                zeros(elemT, size(network.nodeVel, 1), numNewEntries),
            ),
            nodeForce = hcat(
                network.nodeForce,
                zeros(elemT, size(network.nodeForce, 1), numNewEntries),
            ),
            numNode = network.numNode,
            numSeg = network.numSeg,
            maxConnect = network.maxConnect,
            connectivity = hcat(
                network.connectivity,
                zeros(Int, size(network.connectivity, 1), numNewEntries),
            ),
            linksConnect = network.linksConnect,
            segIdx = network.segIdx,
            segForce = network.segForce,
        )
    end

    links = network.links
    bVec = network.bVec
    network.coord[:, newNode] = midCoord
    network.label[newNode] = 1
    network.nodeVel[:, newNode] = midVel
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # Create connectivity for the new node and update links and linksConnect.
    # This only does the first connection. The rest are made after in case they are needed.
    connectivity[1, newNode] = 1
    connectivity[2, newNode] = tmpConnect[(2 * splitConnect)]
    connectivity[3, newNode] = tmpConnect[(2 * splitConnect + 1)]

    link = connectivity[2, newNode]
    colLink = connectivity[3, newNode]
    links[colLink, link] = newNode
    linksConnect[colLink, link] = 1

    # Check if we need a new link between newNode and splitNode for Burgers vector conservation.
    b = SVector{3,elemT}(bVec[1, link], bVec[2, link], bVec[3, link])
    # If burgers vector is conserved return.
    b ⋅ b == 0 && return network

    # New link created at the end.
    network.numSeg[1] += 1
    newSeg = network.numSeg[1]

    # Allocate memory.
    if newSeg > size(network.links, 2)
        # A good heuristic for memory allocation is to allocate an extra N log₂(N) entries.
        numNewEntries = Int(round(newSeg * log2(newSeg)))

        network = DislocationNetwork(;
            links = hcat(network.links, zeros(Int, size(network.links, 1), numNewEntries)),
            slipPlane = hcat(
                network.slipPlane,
                zeros(elemT, size(network.slipPlane, 1), numNewEntries),
            ),
            bVec = hcat(network.bVec, zeros(elemT, size(network.bVec, 1), numNewEntries)),
            coord = network.coord,
            label = network.label,
            nodeVel = network.nodeVel,
            nodeForce = network.nodeForce,
            numNode = network.numNode,
            numSeg = network.numSeg,
            maxConnect = network.maxConnect,
            connectivity = network.connectivity,
            linksConnect = hcat(
                network.linksConnect,
                zeros(Int, size(network.linksConnect, 1), numNewEntries),
            ),
            segIdx = vcat(network.segIdx, zeros(Int, numNewEntries, 3)),
            segForce = cat(
                network.segForce,
                zeros(elemT, size(network.segForce, 1), 2, numNewEntries),
                dims = 3,
            ),
        )
    end

    bVec = network.bVec
    links = network.links
    slipPlane = network.slipPlane
    connectivity = network.connectivity
    linksConnect = network.linksConnect

    # Add new link to network.
    bVec[:, newSeg] = b
    links[colLink, newSeg] = splitNode  # First node in the link is the node that was split.
    colOppLink = 3 - colLink            # Other column in links.
    links[colOppLink, newSeg] = newNode # Second node of the link is the new node.

    # Update connectivity of splitNode.
    newConnect1 = connectivity[1, splitNode] + 1
    connectivity[1, splitNode] = newConnect1
    connectivity[(2 * newConnect1):(2 * newConnect1 + 1), splitNode] .= (newSeg, colLink)

    # Update connectivity of newNode.
    newConnect2 = connectivity[1, newNode] + 1
    connectivity[1, newNode] = newConnect2
    connectivity[(2 * newConnect2):(2 * newConnect2 + 1), newNode] .= (newSeg, colOppLink)

    # Update linksConnect.
    linksConnect[colLink, newSeg] = newConnect1
    linksConnect[colOppLink, newSeg] = newConnect2

    # Fix slip plane.
    # Line direction and velocity of the resultant dislocation.
    coord = network.coord
    nodeVel = network.nodeVel

    # WARNING This calculation's dodgy. Try using the cross product of the adjacent segments.
    t = SVector{3,elemT}(
        coord[1, splitNode] - coord[1, newNode],
        coord[2, splitNode] - coord[2, newNode],
        coord[3, splitNode] - coord[3, newNode],
    )

    v = SVector{3,elemT}(
        nodeVel[1, splitNode] + nodeVel[1, newNode],
        nodeVel[2, splitNode] + nodeVel[2, newNode],
        nodeVel[3, splitNode] + nodeVel[3, newNode],
    )

    # Potential new slip plane.
    n1 = t × b  # For non-screw segments.
    n2 = t × v  # For screw segments.
    if n1 ⋅ n1 > eps(elemT) # non-screw
        slipPlane[:, newSeg] = n1 / norm(n1)
    elseif n2 ⋅ n2 > eps(elemT) # screw
        slipPlane[:, newSeg] = n2 / norm(n2)
    end

    return network
end

function refineNetwork!(
    dlnParams::T1,
    matParams::T2,
    network::T3,
    # mesh::RegularCuboidMesh,
    # dlnFEM::DislocationFEMCorrective;
) where {T1 <: DislocationParameters,T2 <: MaterialParameters,T3 <: DislocationNetwork}

    maxAreaSq = dlnParams.maxAreaSq
    maxSegLen = dlnParams.maxSegLen
    twoMinSegLen = dlnParams.twoMinSegLen

    links = network.links
    coord = network.coord
    label = network.label
    numNode = network.numNode[1]
    nodeVel = network.nodeVel
    connectivity = network.connectivity

    elemT = eltype(network.coord)

    @inbounds for i in 1:numNode
        if connectivity[1, i] == 2 && getNodeType(label[i]) == intMob
            link1 = connectivity[2, i]  # First connection.
            link2 = connectivity[4, i]  # Second connection.
            colLink1 = connectivity[3, i]   # Column where node i is in links of the first connection.
            colLink2 = connectivity[5, i]   # Column where node i is in links of the second connection.
            oppColLink1 = 3 - colLink1 # Node i is connected via link 1 to the node that is in this column in links.
            oppColLink2 = 3 - colLink2 # Node i is connected via link 2 to the node that is in this column in links.
            link1_nodeOppI = links[oppColLink1, link1] # Node i is connected to this node as part of link 1.
            link2_nodeOppI = links[oppColLink2, link2] # Node i is connected to this node as part of link 2.

            # Create triangle formed by the node and its two links.
            iCoord = SVector{3,elemT}(coord[1, i], coord[2, i], coord[3, i])
            # Side 1
            coordVec1 =
                SVector{3,elemT}(
                    coord[1, link1_nodeOppI],
                    coord[2, link1_nodeOppI],
                    coord[3, link1_nodeOppI],
                ) - iCoord
            # Side 2
            coordVec2 =
                SVector{3,elemT}(
                    coord[1, link2_nodeOppI],
                    coord[2, link2_nodeOppI],
                    coord[3, link2_nodeOppI],
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

            # Check if we have to split the second link.
            # Split node the area is greater than the maximum area and if by splitting we would not create two segments smaller than the minimum allowed and if the node is in the simulation. Or if r2 is longer than the maximum length allowed.
            if (areaSq > maxAreaSq && r2 >= twoMinSegLen && link2_nodeOppI <= numNode) ||
               r2 > maxSegLen
                midCoord =
                    SVector{3,elemT}(
                        coord[1, i] + coord[1, link2_nodeOppI],
                        coord[2, i] + coord[2, link2_nodeOppI],
                        coord[3, i] + coord[3, link2_nodeOppI],
                    ) / 2

                midVel =
                    SVector{3,elemT}(
                        nodeVel[1, i] + nodeVel[1, link2_nodeOppI],
                        nodeVel[2, i] + nodeVel[2, link2_nodeOppI],
                        nodeVel[3, i] + nodeVel[3, link2_nodeOppI],
                    ) / 2

                network = splitNode!(network, i, 2, midCoord, midVel)
                getSegmentIdx!(network)
                links = network.links
                slipPlane = network.slipPlane
                coord = network.coord
                label = network.label
                nodeVel = network.nodeVel
                connectivity = network.connectivity

                newNode = network.numNode[1]
                newLink = network.numSeg[1]

                equalSlipPlane = let
                    flag = true
                    @inbounds @simd for j in 1:3
                        flag = flag && isapprox(slipPlane[j, link2], slipPlane[j, link1])
                    end
                    flag
                end
                if equalSlipPlane
                    @inbounds @simd for j in 1:3
                        slipPlane[j, newLink] = slipPlane[j, link2]
                    end
                end

                # Calculate force and mobility for the new node's connectivity.
                j = 1:connectivity[1, newNode]
                link = connectivity[2 * j, newNode]
                colLink = connectivity[2 * j .+ 1, newNode]
                oppColLink = 3 .- colLink
                oldNode = links[oppColLink, link]
                # Calculate segment force for segment link.
                calcSegForce!(dlnParams, matParams, network, link)
                # Calculate old node velocity.
                dlnMobility!(dlnParams, matParams, network, oldNode)
                nodeVel = network.nodeVel

                # Calculate new node velocity.
                dlnMobility!(dlnParams, matParams, network, newNode)
                nodeVel = network.nodeVel
            end

            # Check if we have to split the first link.
            # Split node the area is greater than the maximum area and if by splitting we would not create two segments smaller than the minimum allowed and if the node is in the simulation. Or if r1 is longer than the maximum length allowed.
            if (areaSq > maxAreaSq && r1 >= twoMinSegLen && link1_nodeOppI <= numNode) ||
               r1 > maxSegLen
                midCoord =
                    SVector{3,elemT}(
                        coord[1, i] + coord[1, link1_nodeOppI],
                        coord[2, i] + coord[2, link1_nodeOppI],
                        coord[3, i] + coord[3, link1_nodeOppI],
                    ) / 2

                midVel =
                    SVector{3,elemT}(
                        nodeVel[1, i] + nodeVel[1, link1_nodeOppI],
                        nodeVel[2, i] + nodeVel[2, link1_nodeOppI],
                        nodeVel[3, i] + nodeVel[3, link1_nodeOppI],
                    ) / 2

                network = splitNode!(network, i, 1, midCoord, midVel)
                getSegmentIdx!(network)
                links = network.links
                slipPlane = network.slipPlane
                coord = network.coord
                label = network.label
                nodeVel = network.nodeVel
                connectivity = network.connectivity

                newNode = network.numNode[1]
                newLink = network.numSeg[1]

                equalSlipPlane = let
                    flag = true
                    @inbounds @simd for j in 1:3
                        flag = flag && isapprox(slipPlane[j, link1], slipPlane[j, link2])
                    end
                    flag
                end
                if equalSlipPlane
                    @inbounds @simd for j in 1:3
                        slipPlane[j, newLink] = slipPlane[j, link1]
                    end
                end

                # Calculate force and mobility for the new node's connectivity.
                j = 1:connectivity[1, newNode]
                link = connectivity[2 * j, newNode]
                colLink = connectivity[2 * j .+ 1, newNode]
                oppColLink = 3 .- colLink
                oldNode = links[oppColLink, link]
                # Calculate segment force for segment link.
                calcSegForce!(dlnParams, matParams, network, link)
                # Calculate old node velocity.
                dlnMobility!(dlnParams, matParams, network, oldNode)
                nodeVel = network.nodeVel

                # Calculate new node velocity.
                dlnMobility!(dlnParams, matParams, network, newNode)
                nodeVel = network.nodeVel
            end

        elseif connectivity[1, i] > 2 && getNodeType(label[i]) == intMob
            # Loop through the connections of node i.
            for j in 1:connectivity[1, i]
                # Find the line direction of the link.
                link = connectivity[2 * j, i]
                colLink = connectivity[2 * j + 1, i]
                colOppLink = 3 - colLink
                link_nodeOpp = links[colOppLink, link]
                t = SVector{3,elemT}(
                    coord[1, link_nodeOpp] - coord[1, i],
                    coord[2, link_nodeOpp] - coord[2, i],
                    coord[3, link_nodeOpp] - coord[3, i],
                )
                r1 = norm(t)

                # If the link is smaller than the maximum segment length we skip to the next iteration.
                r1 < maxSegLen ? continue : nothing

                midCoord =
                    SVector{3,elemT}(
                        coord[1, i] + coord[1, link_nodeOpp],
                        coord[2, i] + coord[2, link_nodeOpp],
                        coord[3, i] + coord[3, link_nodeOpp],
                    ) / 2

                midVel =
                    SVector{3,elemT}(
                        nodeVel[1, i] + nodeVel[1, link_nodeOpp],
                        nodeVel[2, i] + nodeVel[2, link_nodeOpp],
                        nodeVel[3, i] + nodeVel[3, link_nodeOpp],
                    ) / 2

                network = splitNode!(network, i, j, midCoord, midVel)
                getSegmentIdx!(network)
                links = network.links
                slipPlane = network.slipPlane
                coord = network.coord
                label = network.label
                nodeVel = network.nodeVel
                connectivity = network.connectivity

                newNode = network.numNode[1]
                newLink = network.numSeg[1]

                @inbounds @simd for k in 1:3
                    slipPlane[k, newLink] = slipPlane[k, link]
                end

                # Calculate force and mobility for the new node's connectivity.
                k = 1:connectivity[1, newNode]
                link = connectivity[2 * k, newNode]
                colLink = connectivity[2 * k .+ 1, newNode]
                colOppLink = 3 .- colLink
                oldNode = links[colOppLink, link]
                # Calculate segment force for segment link.
                calcSegForce!(dlnParams, matParams, network, link)
                # Calculate old node velocity.
                dlnMobility!(dlnParams, matParams, network, oldNode)
                nodeVel = network.nodeVel

                # Calculate new node velocity.
                dlnMobility!(dlnParams, matParams, network, newNode)
                nodeVel = network.nodeVel
            end
        end
    end
    return network
end
