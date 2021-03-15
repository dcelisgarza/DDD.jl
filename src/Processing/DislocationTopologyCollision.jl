function detectCollision(dlnParams::DislocationParameters, network::DislocationNetwork)
    
    collisionDistSq = dlnParams.collisionDistSq
    label = network.label
    numNode = network.numNode[1]
    coord = network.coord
    nodeVel = network.nodeVel
    links = network.links
    connectivity = network.connectivity
    elemT = eltype(coord)
    smallestMinDist::elemT = 0
    collision = false
    collisionType = :null
    n1s1, n2s1, n1s2, n2s2, s1, s2 = 0, 0, 0, 0, 0, 0

    # Hinge condition.
    for i in 1:numNode
        # Only hinge collisions with mobile nodes.
        if !(label[i] == intMobDln || label[i] == srfMobDln)
            continue
        end

        for j in 1:connectivity[1, i]
            # Find connected node.
            link1, rowColLink1, connectedNode1 = findConnectedNode(network, i, j)
            # Make t vector.
            tVec = SVector{3,elemT}(
                coord[1, i] - coord[1, connectedNode1],
                coord[2, i] - coord[2, connectedNode1],
                coord[3, i] - coord[3, connectedNode1]
            )
            # Square of the norm.
            tVecN1 = tVec ⋅ tVec

            # Other connections.
            for k in j + 1:connectivity[1, i]
                link2, rowColLink2, connectedNode2 = findConnectedNode(network, i, k)
                tVec = SVector{3,elemT}(
                        coord[1, i] - coord[1, connectedNode2],
                        coord[2, i] - coord[2, connectedNode2],
                        coord[3, i] - coord[3, connectedNode2]
                    )
                tVecN2 = tVec ⋅ tVec

                # We collapse one of the lines into a point, we want to make sure we collapse the smallest distance so the distance calculation is the actual minimum.
                if tVecN1 < tVecN2
                    link2, rowColLink2, connectedNode2, tVecN2, link1, rowColLink1, connectedNode1, tVecN1 = link1, rowColLink1, connectedNode1, tVecN1, link2, rowColLink2, connectedNode2, tVecN2
                end

                # tVecN1 should now be bigger than tVecN2. i is the hinge node.
                x0 = SVector{3,elemT}(
                        coord[1, i],
                        coord[2, i],
                        coord[3, i]
                    )
                vx0 = SVector{3,elemT}(
                        nodeVel[1, i],
                        nodeVel[2, i],
                        nodeVel[3, i]
                    )

                x1 = SVector{3,elemT}(
                        coord[1, connectedNode1],
                        coord[2, connectedNode1],
                        coord[3, connectedNode1]
                    )
                vx1 = SVector{3,elemT}(
                        nodeVel[1, connectedNode1],
                        nodeVel[2, connectedNode1],
                        nodeVel[3, connectedNode1]
                    )

                # We collapse the other segment into a point.
                y0 = SVector{3,elemT}(
                        coord[1, connectedNode2],
                        coord[2, connectedNode2],
                        coord[3, connectedNode2]
                    )
                vy0 = SVector{3,elemT}(
                        nodeVel[1, connectedNode2],
                        nodeVel[2, connectedNode2],
                        nodeVel[3, connectedNode2]
                    )

                distSq, dDistSqDt, L1, L2 = minimumDistance(
                    x0, x1, y0, y0, vx0, vx1, vy0, vy0
                )

                newSegLen = L1^2 * tVecN1

                if distSq < collisionDistSq && dDistSqDt < -eps(elemT) && newSegLen > collisionDistSq
                    smallestMinDistTmp = 1 / distSq
                    if smallestMinDistTmp > smallestMinDist
                        smallestMinDist = smallestMinDistTmp
                        collision = true
                        collisionType = :hinge
                        n1s1 = i
                        n2s1 = connectedNode1
                        n1s2 = connectedNode2
                        n2s2 = n1s2
                        s1 = link1
                        s2 = link2
                    end
                end
            end
        end
    end

    if collision
        return collision, collisionType, n1s1, n2s1, n1s2, n2s2, s1, s2
    end

    return collision, collisionType, n1s1, n2s1, n1s2, n2s2, s1, s2
end