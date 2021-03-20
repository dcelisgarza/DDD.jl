function detectCollision(dlnParams::DislocationParameters, network::DislocationNetwork, skipSegs)
    
    collisionDistSq = dlnParams.collisionDistSq
    label = network.label
    numNode = network.numNode[1]
    numSeg = network.numSeg[1]
    coord = network.coord
    nodeVel = network.nodeVel
    links = network.links
    bVec = network.bVec
    connectivity = network.connectivity
    elemT = eltype(coord)
    smallestMinDist = zero(elemT)
    collision = false
    collisionType = :null
    n1s1, n2s1, n1s2, n2s2, s1, s2 = 0, 0, 0, 0, 0, 0
    L1, L2 = zero(elemT), zero(elemT)
    mobileNodes = (intMobDln, srfMobDln)

    # Hinge condition.
    for i in 1:numNode
        # Only hinge collisions with mobile nodes.
        if label[i] ∉ mobileNodes
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

                if distSq < collisionDistSq && dDistSqDt < -eps(elemT)
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
        return collision, collisionType, n1s1, n2s1, n1s2, n2s2, s1, s2, L1, L2
    end

    # Collision of unconnected segments, twoline collision.
    for i in 1:numSeg
        n1s1_i = links[1, i]
        n2s1_i = links[2, i]

        label1 = label[n1s1_i]
        label2 = label[n2s1_i]

        if label1 ∉ mobileNodes || label2 ∉ mobileNodes
            continue
        end

        x0 = SVector{3,elemT}(
                coord[1, n1s1_i],
                coord[2, n1s1_i],
                coord[3, n1s1_i]
            )
        vx0 = SVector{3,elemT}(
                nodeVel[1, n1s1_i],
                nodeVel[2, n1s1_i],
                nodeVel[3, n1s1_i]
            )
                
        x1 = SVector{3,elemT}(
                coord[1, n2s1_i],
                coord[2, n2s1_i],
                coord[3, n2s1_i]
            )
        vx1 = SVector{3,elemT}(
                nodeVel[1, n2s1_i],
                nodeVel[2, n2s1_i],
                nodeVel[3, n2s1_i]
            )
        
        bi = SVector{3,elemT}(
                bVec[1, i],
                bVec[2, i],
                bVec[3, i]
            )

        for j in i + 1:numSeg
            n1s2_j = links[1, j]
            n2s2_j = links[2, j]

            # Exclude hinges.
            if n1s1_i != n1s2_j && n1s1_i != n2s2_j && n2s1_i != n1s2_j && n2s1_i != n2s2_j
                label1 = label[n1s2_j]
                label2 = label[n2s2_j]

                if label1 ∉ mobileNodes || label2 ∉ mobileNodes
                    continue
                end
            
                y0 = SVector{3,elemT}(
                    coord[1, n1s2_j],
                    coord[2, n1s2_j],
                    coord[3, n1s2_j]
                )
                vy0 = SVector{3,elemT}(
                    nodeVel[1, n1s2_j],
                    nodeVel[2, n1s2_j],
                    nodeVel[3, n1s2_j]
                )
                    
                y1 = SVector{3,elemT}(
                    coord[1, n2s2_j],
                    coord[2, n2s2_j],
                    coord[3, n2s2_j]
                )
                vy1 = SVector{3,elemT}(
                    nodeVel[1, n2s2_j],
                    nodeVel[2, n2s2_j],
                    nodeVel[3, n2s2_j]
                )
            
                bj = SVector{3,elemT}(
                    bVec[1, j],
                    bVec[2, j],
                    bVec[3, j]
                )

                # Stop superdislocations from forming.
                if ((bi - bj) ⋅ (bi - bj) < eps(elemT) && (x1 - x0) ⋅ (y1 - y0) > eps(elemT)) || ((bi + bj) ⋅ (bi + bj) < eps(elemT) && (x1 - x0) ⋅ (y1 - y0) < eps(elemT))
                    continue
                end

                skip = false
                for k = 1:length(skipSegs)
                    s1, s2 = skipSegs[k]
                    if i == s1 && j == s2 || i == s2 && j == s1
                        skip = true
                        break
                    end
                end

                skip && continue

                distSq, dDistSqDt, L1, L2 = minimumDistance(
                    x0, x1, y0, y1, vx0, vx1, vy0, vy1
                )

                if distSq < collisionDistSq && dDistSqDt < -eps(elemT)
                    smallestMinDistTmp = 1 / distSq
                    if smallestMinDistTmp > smallestMinDist
                        smallestMinDist = smallestMinDistTmp
                        collision = true
                        collisionType = :twoline
                        n1s1 = n1s1_i
                        n2s1 = n2s1_i
                        n1s2 = n1s2_j
                        n2s2 = n1s2_j
                        s1 = i
                        s2 = j
                    end
                end
                
            end
        end
    end

    return collision, collisionType, n1s1, n2s1, n1s2, n2s2, s1, s2, L1, L2
end    