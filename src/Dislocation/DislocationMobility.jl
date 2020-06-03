@inline function dlnMobility(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    nodeIdx = nothing,
)
    return dlnMobility(dlnParams.mobility, dlnParams, matParams, network, nodeIdx)
end
@inline function dlnMobility(
    mobility::mobBCC,
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    nodeIdx = nothing,
)
    # Peierls-Nabarro stress for the bcc material.
    σPN = matParams.σPN
    # Drag coefficients.
    edgeDrag = dlnParams.edgeDrag
    screwDrag = dlnParams.screwDrag
    climbDrag = dlnParams.climbDrag
    lineDrag = dlnParams.lineDrag

    links = network.links
    bVec = network.bVec
    coord = network.coord
    maxConnect = network.maxConnect
    segForce = network.segForce
    connectivity = network.connectivity
    I3 = SMatrix{3, 3}(I)

    # Do it for all nodes if no list is provided.
    if isnothing(nodeIdx)
        numNode = network.numNode
        nodeRange = 1:numNode
    else
        numNode = length(nodeIdx)
        nodeRange = nodeIdx
    end

    nodeForce = zeros(3, numNode)
    nodeVel = zeros(3, numNode)

    # Loop through nodes.
    @fastmath @inbounds for (i, node1) in enumerate(nodeRange)
        totalDrag = SMatrix{3, 3, Float64}(0, 0, 0, 0, 0, 0, 0, 0, 0)
        iNodeForce = SVector{3, Float64}(0, 0, 0)
        iNodeVel = SVector{3, Float64}(0, 0, 0)
        numConnect = connectivity[1, node1]

        # Loop through number of connections.
        for connect in 1:numConnect
            link = connectivity[2 * connect, node1]
            colLink = connectivity[2 * connect + 1, node1]
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
                segForce[1, colLink, link],
                segForce[2, colLink, link],
                segForce[3, colLink, link],
            )
            # If the sress is lower than the Peierls-Nabarro stress, the node behaves as if the force acting on it is zero.
            norm(iNodeConForce) * tN < σPN ? iNodeConForce = SVector{3, Float64}(0, 0, 0) :
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
            bType1 = sum(bTypeArr1)
            bType2 = sum(bTypeArr2)
            bType3 = sum(bTypeArr3)
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
            p = b × t / sqrt(sinSqθ)
            q = p × t

            # Eqn (112) from Arsenlis et al 2007 MSMSE 15 553
            screwDragSq = screwDrag^2
            # screw_ξ_glide = 1/sqrt(sin^2(θ) / ξ_edge^2 + cos^2(θ) / ξ_screw^2)
            ScrewGlideDrag = 1 / sqrt(sinSqθ / edgeDrag^2 + cosSqθ / screwDragSq)
            # screw_ξ_climb = sqrt(sin^2(θ) * ξ_climb^2 + cos^2(θ) * ξ_screw^2)
            ScrewClimbDrag = sqrt(sinSqθ * climbDrag^2 + cosSqθ * screwDragSq)

            # ξ += ||b|| * ||t||/2 *
            #      ((screw_ξ_glide - ξ_screw) * q ⊗ q + (screw_ξ_climb - ξ_screw))
            totalDrag +=
                nBnT_2 * (
                    (ScrewGlideDrag - screwDrag) * q ⊗ q +
                    (ScrewClimbDrag - screwDrag) * p ⊗ p
                )

        end
        # Solve for velocity.
        # ξ v = f
        try
            iNodeVel = totalDrag \ iNodeForce
        catch SingularSystem
            origTotalDrag = totalDrag
            while true
                try
                    totalDrag += I3 * maximum(abs.(origTotalDrag)) * sqrt(eps(Float64))
                    iNodeVel = totalDrag \ iNodeForce
                    break
                catch SingularSystem
                end
            end
        end

        nodeForce[:, i] = iNodeForce
        nodeVel[:, i] = iNodeVel
    end

    return nodeForce, nodeVel
end
