function dlnMobility!(
    mobility::mobBCC,
    dlnParams::DislocationP,
    network::DislocationNetwork,
    σ_PN, # Peierls-Nabarro stress for the bcc material.
    nodeIdx = nothing,
    conList = nothing,
)

    edgeDrag = dlnParams.edgeDrag
    screwDrag = dlnParams.screwDrag
    climbDrag = dlnParams.climbDrag
    lineDrag = dlnParams.lineDrag

    links = network.links
    bVec = network.links
    coord = network.coord
    maxConnect = network.maxConnect
    segForce = network.segForce
    I_st = SMatrix{3, 3}(I)

    # Do it for all nodes if no list is provided.
    if isnothing(nodeIdx)
        numNodes = network.numNodes
        connectivity = network.connectivity
        nodeRange = 1:numNodes
        conList = zeros(maxConnect + 1, numNodes)
        @inbounds @simd for i in nodeRange
            numConnect = connectivity[1, i]
            conList[1, i] = numConnect
            conList[2:(numConnect + 1), i] = 1:numConnect
        end
    else
        numNodes = length(nodeIdx)
        nodeRange = nodeIdx
        @assert size(conList) == (maxConnect + 1, numNodes)
    end

    nodeForce = zeros(3, numNodes)
    nodeVel = zeros(3, numNodes)

    # Loop through nodes.
    for (i, node1) in enumerate(nodeRange)
        totalDrag = SMatrix{3, 3, Float64}(0, 0, 0, 0, 0, 0, 0, 0, 0)
        iNodeForce = SVector{3, Float64}(0, 0, 0)
        numConnect = conList[1, i] # Number of connections.

        # Loop through number of connections.
        for connect in 1:numConnect
            nodeCon = conList[connect + 1, i]
            link = connectivity[2 * nodeCon, node1]
            colLink = connectivity[2 * nodeCon + 1, node1]
            colOppLink = 3 - colLink
            node2 = links[colOppLink, link]

            # Line direction, continue to next iteration if norm is 0 and normalise line direction.
            t = SVector{3, Float64}(
                coord[1, node2] - coord[1, node1],
                coord[2, node2] - coord[2, node1],
                coord[3, node2] - coord[3, node1],
            )
            nT = norm(t)
            # If the segment length is zero, skip to the next iteration.
            nT < eps(Float64) ? continue : nothing
            tN = 1 / nT
            t *= tN # Normalise t.

            # Segment force relevant to node1.
            iNodeConForce = SVector{3, Float64}(
                segForce[1, colOppLink, link],
                segForce[2, colOppLink, link],
                segForce[3, colOppLink, link],
            )
            # If the sress is lower than the Peierls-Nabarro stress, the node behaves as if the force acting on it is zero.
            norm(iNodeConForce) * tN < σ_PN ? iNodeConForce = SVector{3, Float64}(0, 0, 0) :
            nothing
            # Add force from this connection to the total force on the node.
            iNodeForce += iNodeConForce

            # Burgers Vector
            b = SVector{3, Float64}(bVec[1, link], bVec[2, link], bVec[3, link])
            nB = norm(b)
            # Burgers vector families of screw dislocations in BCC, in case there have been collisions and the Burgers vector norm is greater than 1.
            absB = abs.(b)
            absB /= minimum(absB)
            bTypeArr1 = @. abs(1 - absB) < sqrt(eps(Float64))
            bTypeArr2 = @. abs(2 - absB) < sqrt(eps(Float64))
            bTypeArr3 = @. abs(3 - absB) < sqrt(eps(Float64))
            bType1 = sum(flagArr1)
            bType2 = sum(flagArr2)
            bType3 = sum(flagArr3)
            screw = 1 - nB < sqrt(eps(Float64))
            # Normalise Burgers vector.
            b /= nB

            nBnT_2 = nB * nT / 2
            # If it's not a screw-ish dislocation, do edge calculation and continue to next iteration.
            if !(
                bType1 == 3 && screw ||
                bType1 == 2 && bType2 == 1 && screw ||
                bType1 == 1 && bType2 == 1 && bType3 == 1 && screw
            )
                # ξ += |b|| * ||t||/2 * (ξ_climb * I + (ξ_line - ξ_climb) * t ⊗ t)
                totalDrag += nBnT_2 * (climbDrag * I3 + (lineDrag - climbDrag) * t ⊗ t)
                continue
            end

            # cos(θ) = (a ⋅ b)/(||a|| ||b||)
            # Both vectors have been normalised already.
            tDb = t ⋅ b
            cosSqθ = tDb^2
            sinSqθ = 1 - cosSqθ

            # ξ += ||b|| * ||t||/2 (ξ_screw * I + (ξ_line - ξ_screw) * t ⊗ t)
            totalDrag += nBnT_2 * (screwDrag * I3 + (lineDrag - screwDrag) * t ⊗ t)

            # If dislocation segment is pure screw skip to next iteration.
            sinSqθ < eps(Float64) ? continue : nothing

            # P, Q vectors as new local axes.
            p = b × t / tDb
            q = p × t

            # Eqn (112) from Arsenlis et al 2007 MSMSE 15 553
            screwDragSq = screwDrag^2
            # screw_ξ_glide = 1/sqrt(sin^2(θ) / ξ_edge^2 + cos^2(θ) / ξ_screw^2)
            glideDrag = 1 / sqrt(sinSqθ / edgeDrag^2 + cosSqθ / screwDragSq)
            # screw_ξ_climb = sqrt(sin^2(θ) * ξ_climb^2 + cos^2(θ) * ξ_screw^2)
            climbDrag = sqrt(sinSqθ * climbDrag^2 + cosSqθ * screwDragSq)

            # ξ += ||b|| * ||t||/2 *
            #      ((screw_ξ_glide - ξ_screw) * q ⊗ q + (screw_ξ_climb - ξ_screw))
            totalDrag +=
                nBnT_2 * ((glideDrag - screwDrag) * q ⊗ q + (climbDrag - screwDrag) * p ⊗ p)

        end
        # Solve for velocity.
        # ξ v = f
        try
            iNodeVel = totalDrag \ iNodeForce
        catch e
            totalDrag += I3 * maximum(abs.(totalDrag)) * sqrt(eps(Float64)) * 10
            iNodeVel = totalDrag \ iNodeForce
        end

        nodeForce[1, i] = iNodeForce[1]
        nodeForce[2, i] = iNodeForce[2]
        nodeForce[3, i] = iNodeForce[3]
        nodeVel[1, i] = iNodeVel[1]
        nodeVel[2, i] = iNodeVel[2]
        nodeVel[3, i] = iNodeVel[3]
    end
end
