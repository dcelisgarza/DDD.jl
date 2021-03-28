function detectCollision(
    dlnParams::DislocationParameters,
    network::DislocationNetwork,
    skipSegs,
)
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
        label[i] ∉ mobileNodes && continue

        for j in 1:connectivity[1, i]
            # Find connected node.
            link1, rowColLink1, connectedNode1 = findConnectedNode(network, i, j)
            # Make t vector.
            tVec = SVector{3, elemT}(
                coord[1, i] - coord[1, connectedNode1],
                coord[2, i] - coord[2, connectedNode1],
                coord[3, i] - coord[3, connectedNode1],
            )
            # Square of the norm.
            tVecN1 = tVec ⋅ tVec

            # Other connections.
            for k in (j + 1):connectivity[1, i]
                link2, rowColLink2, connectedNode2 = findConnectedNode(network, i, k)
                tVec = SVector{3, elemT}(
                    coord[1, i] - coord[1, connectedNode2],
                    coord[2, i] - coord[2, connectedNode2],
                    coord[3, i] - coord[3, connectedNode2],
                )
                tVecN2 = tVec ⋅ tVec

                # We collapse one of the lines into a point, we want to make sure we collapse the smallest distance so the distance calculation is the actual minimum.
                tVecN1 < tVecN2 ?
                begin
                    link2,
                    rowColLink2,
                    connectedNode2,
                    tVecN2,
                    link1,
                    rowColLink1,
                    connectedNode1,
                    tVecN1 = link1,
                    rowColLink1,
                    connectedNode1,
                    tVecN1,
                    link2,
                    rowColLink2,
                    connectedNode2,
                    tVecN2
                end : nothing

                # tVecN1 should now be bigger than tVecN2. i is the hinge node.
                x0 = SVector{3, elemT}(coord[1, i], coord[2, i], coord[3, i])
                vx0 = SVector{3, elemT}(nodeVel[1, i], nodeVel[2, i], nodeVel[3, i])

                x1 = SVector{3, elemT}(
                    coord[1, connectedNode1],
                    coord[2, connectedNode1],
                    coord[3, connectedNode1],
                )
                vx1 = SVector{3, elemT}(
                    nodeVel[1, connectedNode1],
                    nodeVel[2, connectedNode1],
                    nodeVel[3, connectedNode1],
                )

                # We collapse the other segment into a point.
                y0 = SVector{3, elemT}(
                    coord[1, connectedNode2],
                    coord[2, connectedNode2],
                    coord[3, connectedNode2],
                )
                vy0 = SVector{3, elemT}(
                    nodeVel[1, connectedNode2],
                    nodeVel[2, connectedNode2],
                    nodeVel[3, connectedNode2],
                )

                skip = false
                for k in 1:length(skipSegs)
                    s1Tmp, s2Tmp = skipSegs[k]
                    if link1 == s1Tmp && link2 == s2Tmp || link1 == s2Tmp && link2 == s1Tmp
                        skip = true
                        break
                    end
                end

                skip && continue

                distSq, dDistSqDt, L1, L2 =
                    minimumDistance(x0, x1, y0, y0, vx0, vx1, vy0, vy0)

                if distSq < collisionDistSq && dDistSqDt < -eps(elemT)
                    smallestMinDistTmp = 1 / distSq
                    if smallestMinDistTmp > smallestMinDist
                        smallestMinDist = smallestMinDistTmp
                        collision = true
                        collisionType = :hinge
                        n1s1 = i # Hinge node
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

    collision && return collision, collisionType, n1s1, n2s1, n1s2, n2s2, s1, s2, L1, L2

    # Collision of unconnected segments, twoline collision.
    for i in 1:numSeg
        n1s1_i = links[1, i]
        n2s1_i = links[2, i]

        label1 = label[n1s1_i]
        label2 = label[n2s1_i]

        (label1 ∉ mobileNodes || label2 ∉ mobileNodes) && continue

        x0 = SVector{3, elemT}(coord[1, n1s1_i], coord[2, n1s1_i], coord[3, n1s1_i])
        vx0 = SVector{3, elemT}(nodeVel[1, n1s1_i], nodeVel[2, n1s1_i], nodeVel[3, n1s1_i])

        x1 = SVector{3, elemT}(coord[1, n2s1_i], coord[2, n2s1_i], coord[3, n2s1_i])
        vx1 = SVector{3, elemT}(nodeVel[1, n2s1_i], nodeVel[2, n2s1_i], nodeVel[3, n2s1_i])

        bi = SVector{3, elemT}(bVec[1, i], bVec[2, i], bVec[3, i])

        for j in (i + 1):numSeg
            n1s2_j = links[1, j]
            n2s2_j = links[2, j]

            # Exclude hinges.
            if n1s1_i != n1s2_j && n1s1_i != n2s2_j && n2s1_i != n1s2_j && n2s1_i != n2s2_j
                label1 = label[n1s2_j]
                label2 = label[n2s2_j]

                (label1 ∉ mobileNodes || label2 ∉ mobileNodes) && continue

                y0 = SVector{3, elemT}(coord[1, n1s2_j], coord[2, n1s2_j], coord[3, n1s2_j])
                vy0 = SVector{3, elemT}(
                    nodeVel[1, n1s2_j],
                    nodeVel[2, n1s2_j],
                    nodeVel[3, n1s2_j],
                )

                y1 = SVector{3, elemT}(coord[1, n2s2_j], coord[2, n2s2_j], coord[3, n2s2_j])
                vy1 = SVector{3, elemT}(
                    nodeVel[1, n2s2_j],
                    nodeVel[2, n2s2_j],
                    nodeVel[3, n2s2_j],
                )

                bj = SVector{3, elemT}(bVec[1, j], bVec[2, j], bVec[3, j])

                # Stop superdislocations from forming.
                (
                    (
                        (bi - bj) ⋅ (bi - bj) < eps(elemT) &&
                        (x1 - x0) ⋅ (y1 - y0) > eps(elemT)
                    ) || (
                        (bi + bj) ⋅ (bi + bj) < eps(elemT) &&
                        (x1 - x0) ⋅ (y1 - y0) < eps(elemT)
                    )
                ) && continue

                skip = false
                for k in 1:length(skipSegs)
                    s1Tmp, s2Tmp = skipSegs[k]
                    if i == s1Tmp && j == s2Tmp || i == s2Tmp && j == s1Tmp
                        skip = true
                        break
                    end
                end

                skip && continue

                distSq, dDistSqDt, L1, L2 =
                    minimumDistance(x0, x1, y0, y1, vx0, vx1, vy0, vy1)

                if distSq < collisionDistSq && dDistSqDt < -eps(elemT)
                    smallestMinDistTmp = 1 / distSq
                    if smallestMinDistTmp > smallestMinDist
                        smallestMinDist = smallestMinDistTmp
                        collision = true
                        collisionType = :twoLine
                        n1s1 = n1s1_i
                        n2s1 = n2s1_i
                        n1s2 = n1s2_j
                        n2s2 = n2s2_j
                        s1 = i
                        s2 = j
                    end
                end
            end
        end
    end

    return collision, collisionType, n1s1, n2s1, n1s2, n2s2, s1, s2, L1, L2
end

function resolveCollision(
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    mesh::AbstractMesh,
    forceDisplacement::ForceDisplacement,
    network::DislocationNetwork,
    collisionType,
    n1s1,
    n2s1,
    n1s2,
    n2s2,
    s1,
    s2,
    L1,
    L2,
)
    collisionDist = dlnParams.collisionDist
    collisionDistSq = dlnParams.collisionDistSq
    minSegLenSq = dlnParams.minSegLenSq
    coord = network.coord
    nodeVel = network.nodeVel
    segForce = network.segForce
    connectivity = network.connectivity
    linksConnect = network.linksConnect
    elemT = eltype(coord)
    mergeNode1 = 0

    # Identify closeness to node in segment 1 is the same in both collision types.
    tVec = SVector{3, elemT}(
        coord[1, n2s1] - coord[1, n1s1],
        coord[2, n2s1] - coord[2, n1s1],
        coord[3, n2s1] - coord[3, n1s1],
    )
    tVecSq = tVec ⋅ tVec

    close_n1s1 = ((L1^2 * tVecSq) < minSegLenSq)
    close_n2s1 = (((1 - L1)^2 * tVecSq) < minSegLenSq)

    if collisionType == :twoLine
        if close_n1s1 && L1 <= 0.5
            mergeNode1 = n1s1
        elseif close_n2s1
            mergeNode1 = n2s1
        else
            splitNode = n1s1
            splitConnect = linksConnect[1, s1]

            if close_n1s1
                L1 = collisionDist / sqrt(tVec)
            else
                L1 = 1 - collisionDist / sqrt(tVec)
            end

            midCoord =
                SVector{3, elemT}(coord[1, n1s1], coord[2, n1s1], coord[3, n1s1]) *
                (1 - L1) +
                SVector{3, elemT}(coord[1, n2s1], coord[2, n2s1], coord[3, n2s1]) * L1

            midVel =
                SVector{3, elemT}(nodeVel[1, n1s1], nodeVel[2, n1s1], nodeVel[3, n1s1]) *
                (1 - L1) +
                SVector{3, elemT}(nodeVel[1, n2s1], nodeVel[2, n2s1], nodeVel[3, n2s1]) * L1

            splitNode!(network, splitNode, splitConnect, midCoord, midVel)
            mergeNode1 = network.numNode[1]
        end

        # Second node to merge.
        tVec = SVector{3, elemT}(
            coord[1, n2s2] - coord[1, n1s2],
            coord[2, n2s2] - coord[2, n1s2],
            coord[3, n2s2] - coord[3, n1s2],
        )
        tVecSq = tVec ⋅ tVec

        close_n1s2 = ((L2^2 * tVecSq) < minSegLenSq)
        close_n2s2 = (((1 - L2)^2 * tVecSq) < minSegLenSq)

        if close_n1s2 && L2 <= 0.5
            mergeNode2 = n1s2
        elseif close_n2s2
            mergeNode2 = n2s2
        else
            splitNode = n1s2
            splitConnect = linksConnect[1, s2]

            if close_n1s2
                L2 = collisionDist / sqrt(tVec)
            else
                L2 = 1 - collisionDist / sqrt(tVec)
            end

            midCoord =
                SVector{3, elemT}(coord[1, n1s2], coord[2, n1s2], coord[3, n1s2]) *
                (1 - L2) +
                SVector{3, elemT}(coord[1, n2s2], coord[2, n2s2], coord[3, n2s2]) * L2

            midVel =
                SVector{3, elemT}(nodeVel[1, n1s2], nodeVel[2, n1s2], nodeVel[3, n1s2]) *
                (1 - L2) +
                SVector{3, elemT}(nodeVel[1, n2s2], nodeVel[2, n2s2], nodeVel[3, n2s2]) * L2

            splitNode!(network, splitNode, splitConnect, midCoord, midVel)
            mergeNode2 = network.numNode[1]
        end

        # Power dissipation of unmerged structure.

        # Node 1
        force = SVector(zeros(elemT, 3))
        c = connectivity[1, mergeNode1]
        for i in 1:c
            link = connectivity[2 * i, mergeNode1]
            pos = connectivity[2 * i + 1, mergeNode1]
            force += SVector{3, elemT}(
                segForce[1, 3 - pos, link],
                segForce[2, 3 - pos, link],
                segForce[3, 3 - pos, link],
            )
        end
        nodeVelTmp = SVector{3, elemT}(
            nodeVel[1, mergeNode1],
            nodeVel[2, mergeNode1],
            nodeVel[3, mergeNode1],
        )
        powerPreCollision = 1.05 * nodeVelTmp ⋅ force

        # Node 2
        force = SVector(zeros(elemT, 3))
        c = connectivity[1, mergeNode2]
        for i in 1:c
            link = connectivity[2 * i, mergeNode2]
            pos = connectivity[2 * i + 1, mergeNode2]
            force += SVector{3, elemT}(
                segForce[1, 3 - pos, link],
                segForce[2, 3 - pos, link],
                segForce[3, 3 - pos, link],
            )
        end
        nodeVelTmp = SVector{3, elemT}(
            nodeVel[1, mergeNode2],
            nodeVel[2, mergeNode2],
            nodeVel[3, mergeNode2],
        )
        powerPreCollision += 1.05 * nodeVelTmp ⋅ force

        collisionPoint = findCollisionPoint(network, mergeNode1, mergeNode2)

    elseif collisionType == :hinge
    end

    # calcSegForce(dlnParams, matParams, mesh, forceDisplacement, network, idx)
end

function findCollisionPoint(network::DislocationNetwork, node1, node2)
    bVec = network.bVec
    label = network.label
    coord = network.coord
    nodeVel = network.nodeVel
    connectivity = network.connectivity
    elemT = eltype(coord)

    p1 = SVector{3, elemT}(coord[1, node1], coord[2, node1], coord[3, node1])
    p2 = SVector{3, elemT}(coord[1, node2], coord[2, node2], coord[3, node2])
    label1 = label[node1]
    label2 = label[node2]

    if label1 == intFixDln
        collisionPoint = p1
        return collisionPoint
    elseif label2 == intFixDln
        collisionPoint = p2
        return collisionPoint
    end

    N = 0
    matN = zeros(elemT, MMatrix{3, 3})
    vecN = zeros(elemT, MVector{3})

    N = findPlaneCondition!(matN, vecN, N, network, node1, p1)
    N < 3 ? findPlaneCondition!(matN, vecN, N, network, node2, p2) : nothing

    A = SMatrix{6, 6, elemT}([I(3) matN[:, 1:3]'; matN[:, 1:3] zeros(3, 3)])
    b = [(p1 + p2) / 2; vecN]
    x = A \ b
    collisionPoint = x[1:3]

    return collisionPoint
end

function findPlaneCondition!(matN, vecN, N, network, node, point)
    bVec = network.bVec
    coord = network.coord
    nodeVel = network.nodeVel
    connectivity = network.connectivity
    elemT = eltype(coord)

    newPlaneCondition = 0.875

    for i in 1:connectivity[1, node]
        if N < 3
            link, missing, connectedNode = findConnectedNode(network, node, i)
            tVec =
                point - SVector{3, elemT}(
                    coord[1, connectedNode],
                    coord[2, connectedNode],
                    coord[3, connectedNode],
                )
            tVec = normalize(tVec)

            b = SVector{3, elemT}(bVec[1, link], bVec[2, link], bVec[3, link])
            vel = SVector{3, elemT}(nodeVel[1, node], nodeVel[2, node], nodeVel[3, node])

            n1 = tVec × b
            n2 = tVec × vel

            n1Sq = n1 ⋅ n1
            n2Sq = n2 ⋅ n2

            if n1Sq > eps(elemT)
                plane = normalize(n1)
            elseif n2Sq > eps(elemT)
                plane = normalize(n2)
            end

            if n1Sq > eps(elemT) || n2Sq > eps(elemT)
                if N == 0
                    condition = true
                elseif N == 1
                    condition = (matN[:, 1] ⋅ plane)^2 < newPlaneCondition^2
                else
                    detN = det(
                        SMatrix{3, 3, elemT}(
                            matN[1, 1],
                            matN[2, 1],
                            matN[3, 1],
                            matN[1, 2],
                            matN[2, 2],
                            matN[3, 2],
                            plane[1],
                            plane[2],
                            plane[3],
                        ),
                    )
                    condition = detN^2 < (1 - newPlaneCondition)^2
                end

                if condition
                    N += 1
                    matN[:, N] = plane
                    vecN[N] = plane ⋅ point
                end
            end
        end
    end
    return N
end
