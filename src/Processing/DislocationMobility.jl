"""
Original by Bruce Bromage at the Department of Materials of the University of Oxford, @brucebromage on github.

This is outdated, new capabilities include rotating the frame of reference and better handling of cross-slip.
"""
function dlnMobility(
    dlnParams::DislocationParameters{T1,T2,T3,T4} where {T1 <: mobBCC,T2,T3,T4},
    matParams::MaterialParameters,
    network::DislocationNetwork,
    idx = nothing,
)

    # Peierls-Nabarro stress for the bcc material.
    σPN = matParams.σPN
    # Drag coefficients.
    edgeDrag = dlnParams.dragCoeffs[:edge]
    screwDrag = dlnParams.dragCoeffs[:screw]
    climbDrag = dlnParams.dragCoeffs[:climb]
    lineDrag = dlnParams.dragCoeffs[:line]
    dType = typeof(lineDrag)
    I3 = SMatrix{3,3,dType}(I)

    links = network.links
    bVec = network.bVec
    coord = network.coord
    maxConnect = network.maxConnect
    segForce = network.segForce
    connectivity = network.connectivity
    elemT = eltype(network.bVec)

    # Do it for all nodes if no list is provided.
    if isnothing(idx)
        numNode = network.numNode[1]
        idx = 1:numNode
    else
        numNode = length(idx)
    end

    nodeForce = zeros(elemT, 3, numNode)
    nodeVel = zeros(elemT, 3, numNode)

    # Loop through nodes.
    for (i, node1) in enumerate(idx)
        totalDrag = zeros(SMatrix{3,3,dType})
        iNodeForce = zeros(SVector{3,elemT})
        iNodeVel = zeros(SVector{3,elemT})
        numConnect = connectivity[1, node1]

        # Loop through number of connections.
        for connect in 1:numConnect
            link = connectivity[2 * connect, node1]
            colLink = connectivity[2 * connect + 1, node1]
            colOppLink = 3 - colLink
            node2 = links[colOppLink, link]

            # Line direction, continue to next iteration if norm is 0 and normalise line direction.
            t = SVector{3,elemT}(
                coord[1, node2] - coord[1, node1],
                coord[2, node2] - coord[2, node1],
                coord[3, node2] - coord[3, node1],
            )
            nT = norm(t)
            # If the segment length is zero, skip to the next iteration.
            nT < eps(elemT) ? continue : nothing
            tN = 1 / nT
            t *= tN # Normalise t.

            # Segment force relevant to node1.
            iNodeConForce = SVector{3,elemT}(
                segForce[1, colLink, link],
                segForce[2, colLink, link],
                segForce[3, colLink, link],
            )
            # If the sress is lower than the Peierls-Nabarro stress, the node behaves as if the force acting on it is zero.
            norm(iNodeConForce) * tN < σPN ? iNodeConForce = zeros(SVector{3,elemT}) :
            nothing
            # Add force from this connection to the total force on the node.
            iNodeForce += iNodeConForce

            # Burgers Vector
            b = SVector{3,elemT}(bVec[1, link], bVec[2, link], bVec[3, link])
            nB = norm(b)
            # Burgers vector families of screw dislocations in BCC, in case there have been collisions and the Burgers vector norm is greater than 1.
            absB = abs.(b)
            absB /= minimum(absB)
            sqrtBType = sqrt(eps(elemT))
            bTypeArr1 = @. abs(1 - absB) < sqrtBType
            bTypeArr2 = @. abs(2 - absB) < sqrtBType
            bTypeArr3 = @. abs(3 - absB) < sqrtBType
            bType1 = sum(bTypeArr1)
            bType2 = sum(bTypeArr2)
            bType3 = sum(bTypeArr3)
            screw = 1 - nB < sqrtBType
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
            sinSqθ < eps(typeof(tDb)) ? continue : nothing

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
        iNodeVel = totalDrag \ iNodeForce
        #= 
        try
            iNodeVel = totalDrag \ iNodeForce
        catch SingularSystem
            origTotalDrag = totalDrag
            while true
                try
                    totalDrag += I3 * maximum(abs.(origTotalDrag)) * sqrt(eps(dType))
                    iNodeVel = totalDrag \ iNodeForce
                    break
                catch SingularSystem
                end
            end
        end =#

        nodeForce[:, i] = iNodeForce
        nodeVel[:, i] = iNodeVel
    end

    return nodeForce, nodeVel
end
function dlnMobility!(
    dlnParams::DislocationParameters{T1,T2,T3,T4} where {T1 <: mobBCC,T2,T3,T4},
    matParams::MaterialParameters,
    network::DislocationNetwork,
    idx = nothing,
)

    # Peierls-Nabarro stress for the bcc material.
    σPN = matParams.σPN
    # Drag coefficients.
    edgeDrag = dlnParams.dragCoeffs[:edge]
    screwDrag = dlnParams.dragCoeffs[:screw]
    climbDrag = dlnParams.dragCoeffs[:climb]
    lineDrag = dlnParams.dragCoeffs[:line]
    dType = typeof(lineDrag)
    I3 = SMatrix{3,3,dType}(I)

    links = network.links
    bVec = network.bVec
    coord = network.coord
    nodeVel = network.nodeVel
    nodeForce = network.nodeForce
    maxConnect = network.maxConnect
    segForce = network.segForce
    connectivity = network.connectivity
    elemT = eltype(network.bVec)

    # Do it for all nodes if no list is provided.
    if isnothing(idx)
        numNode = network.numNode[1]
        idx = 1:numNode
    else
        numNode = length(idx)
    end

    # Loop through nodes.
    for node1 in idx
        totalDrag = zeros(SMatrix{3,3,dType})
        iNodeForce = zeros(SVector{3,elemT})
        iNodeVel = zeros(SVector{3,elemT})
        numConnect = connectivity[1, node1]

        # Loop through number of connections.
        for connect in 1:numConnect
            link = connectivity[2 * connect, node1]
            colLink = connectivity[2 * connect + 1, node1]
            colOppLink = 3 - colLink
            node2 = links[colOppLink, link]

            # Line direction, continue to next iteration if norm is 0 and normalise line direction.
            t = SVector{3,elemT}(
                coord[1, node2] - coord[1, node1],
                coord[2, node2] - coord[2, node1],
                coord[3, node2] - coord[3, node1],
            )
            nT = norm(t)
            # If the segment length is zero, skip to the next iteration.
            nT < eps(elemT) ? continue : nothing
            tN = 1 / nT
            t *= tN # Normalise t.

            # Segment force relevant to node1.
            iNodeConForce = SVector{3,elemT}(
                segForce[1, colLink, link],
                segForce[2, colLink, link],
                segForce[3, colLink, link],
            )
            # If the sress is lower than the Peierls-Nabarro stress, the node behaves as if the force acting on it is zero.
            norm(iNodeConForce) * tN < σPN ? iNodeConForce = zeros(SVector{3,elemT}) :
            nothing
            # Add force from this connection to the total force on the node.
            iNodeForce += iNodeConForce

            # Burgers Vector
            b = SVector{3,elemT}(bVec[1, link], bVec[2, link], bVec[3, link])
            nB = norm(b)
            # Burgers vector families of screw dislocations in BCC, in case there have been collisions and the Burgers vector norm is greater than 1.
            absB = abs.(b)
            absB /= minimum(absB)
            sqrtBType = sqrt(eps(elemT))
            bTypeArr1 = @. abs(1 - absB) < sqrtBType
            bTypeArr2 = @. abs(2 - absB) < sqrtBType
            bTypeArr3 = @. abs(3 - absB) < sqrtBType
            bType1 = sum(bTypeArr1)
            bType2 = sum(bTypeArr2)
            bType3 = sum(bTypeArr3)
            screw = 1 - nB < sqrtBType
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
            sinSqθ < eps(typeof(tDb)) ? continue : nothing

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

        iNodeVel = totalDrag \ iNodeForce
        #= 
        try
            iNodeVel = totalDrag \ iNodeForce
        catch SingularSystem
            origTotalDrag = totalDrag
            while true
                try
                    totalDrag += I3 * maximum(abs.(origTotalDrag)) * sqrt(eps(dType))
                    iNodeVel = totalDrag \ iNodeForce
                    break
                catch SingularSystem
                end
            end
        end =#

        nodeForce[:, node1] = iNodeForce
        nodeVel[:, node1] = iNodeVel
    end

    return network
end
