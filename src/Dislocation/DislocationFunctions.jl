"""
```
calcSegForce(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork;
    # mesh::RegularCuboidMesh,
    # dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
)
```
"""
@inline function calcSegForce(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    idx = nothing;
    # mesh::RegularCuboidMesh,
    # dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
)
    isnothing(idx) ? numSeg = network.numSeg : numSeg = length(idx)

    # pkForce = pkForce(mesh, dlnFEM, network)
    selfForce = calcSelfForce(dlnParams, matParams, network, idx)
    segForce = calcSegSegForce(dlnParams, matParams, network, idx; parallel = parallel)

    @inbounds for i in 1:numSeg
        @simd for j in 1:2
            segForce[:, j, i] = selfForce[j][:, i] + segForce[:, j, i]
        end
    end

    return segForce
end
@inline function calcSegForce!(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    idx = nothing;
    # mesh::RegularCuboidMesh,
    # dlnFEM::DislocationFEMCorrective;
    parallel::Bool = true,
)

    # pkForce!(mesh, dlnFEM, network)
    calcSelfForce!(dlnParams, matParams, network, idx)
    calcSegSegForce!(dlnParams, matParams, network, idx; parallel = parallel)

    return network
end

"""
```
calcSelfForce(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
)
```
Calculates the self-interaction force felt by two nodes in a segment. Naturally the forces are equal and opposite to each other.
"""
@inline function calcSelfForce(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    idx = nothing,
)

    μ = matParams.μ
    ν = matParams.ν
    E = matParams.E
    omNuInv = matParams.omνInv
    nuOmNuInv = matParams.νomνInv
    μ4π = matParams.μ4π
    a = dlnParams.coreRad
    aSq = dlnParams.coreRadSq
    bVec = network.bVec
    coord = network.coord
    segIdx = network.segIdx

    # Indices for self force.
    if isnothing(idx)
        # If no index is provided, calculate forces for all segments.
        numSeg = network.numSeg
        idx = 1:numSeg
    else
        # Else, calculate forces only on idx.
        numSeg = length(idx)
    end

    idxBvec = @view segIdx[idx, 1]
    idxNode1 = @view segIdx[idx, 2]
    idxNode2 = @view segIdx[idx, 3]
    # Un normalised segment vectors. Use views for speed.
    bVec = @view bVec[:, idxBvec]
    tVec = @views coord[:, idxNode2] - coord[:, idxNode1]

    selfForceNode2 = zeros(3, numSeg)

    @fastmath @simd for i in eachindex(idx)
        # Finding the norm of each line vector.
        tVecSq = tVec[1, i]^2 + tVec[2, i]^2 + tVec[3, i]^2
        L = sqrt(tVecSq)
        Linv = inv(L)
        # Finding the non-singular norm.
        La = sqrt(tVecSq + aSq)
        tVecI = SVector{3, Float64}(tVec[1, i] * Linv, tVec[2, i] * Linv, tVec[3, i] * Linv)
        bVecI = SVector{3, Float64}(bVec[1, i], bVec[2, i], bVec[3, i])
        # Normalised the dislocation network vector, the sum of all the segment vectors has norm 1.
        # Screw component, scalar projection of bVec onto t.
        bScrew = tVecI ⋅ bVecI
        # Edge component, vector rejection of bVec onto t.
        bEdgeVec = bVecI - bScrew * tVecI
        # Finding the norm squared of each edge component.
        bEdgeSq = bEdgeVec ⋅ bEdgeVec
        #=
        A. Arsenlis et al, Modelling Simul. Mater. Sci. Eng. 15 (2007)
        553?595: gives this expression in appendix A p590
        f^{s}_{43} = -(μ/(4π)) [ t × (t × b)](t ⋅ b) { v/(1-v) ( ln[
        (L_a + L)/a] - 2*(L_a - a)/L ) - (L_a - a)^2/(2La*L) }

        tVec × (tVec × bVec)    = tVec (tVec ⋅ bVec) - bVec (tVec ⋅ tVec)
        = tVec * bScrew - bVec
        = - bEdgeVec
        =#
        # Torsional component of the elastic self interaction force. This is the scalar component of the above equation.
        LaMa = La - a
        # Torsional component of core self interaction.
        tor =
            μ4π *
            bScrew *
            (nuOmNuInv * (log((La + L) / a) - 2 * LaMa * Linv) - LaMa^2 / (2 * La * L))
        torCore = 2 * E * nuOmNuInv * bScrew
        torTot = tor + torCore
        # Longitudinal component of core self interaction.
        lonCore = (bScrew^2 + bEdgeSq * omNuInv) * E
        # Force calculation.
        selfForceNode2[:, i] = torTot * bEdgeVec - lonCore * tVecI
    end
    selfForceNode1 = -selfForceNode2

    return selfForceNode1, selfForceNode2
end
@inline function calcSelfForce!(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    idx = nothing,
)

    μ = matParams.μ
    ν = matParams.ν
    E = matParams.E
    omNuInv = matParams.omνInv
    nuOmNuInv = matParams.νomνInv
    μ4π = matParams.μ4π
    a = dlnParams.coreRad
    aSq = dlnParams.coreRadSq
    bVec = network.bVec
    coord = network.coord
    segIdx = network.segIdx
    segForce = network.segForce

    # Indices for self force.
    if isnothing(idx)
        # If no index is provided, calculate forces for all segments.
        numSeg = network.numSeg
        idx = 1:numSeg
    else
        # Else, calculate forces only on idx.
        numSeg = length(idx)
    end

    # Un normalised segment vectors. Use views for speed.
    idxBvec = @view segIdx[idx, 1]
    idxNode1 = @view segIdx[idx, 2]
    idxNode2 = @view segIdx[idx, 3]
    bVec = @view bVec[:, idxBvec]
    tVec = @views coord[:, idxNode2] - coord[:, idxNode1]

    @fastmath @inbounds @simd for i in eachindex(idx)
        # Finding the norm of each line vector.
        tVecSq = tVec[1, i]^2 + tVec[2, i]^2 + tVec[3, i]^2
        L = sqrt(tVecSq)
        Linv = inv(L)
        # Finding the non-singular norm.
        La = sqrt(tVecSq + aSq)
        tVecI = SVector{3, Float64}(tVec[1, i] * Linv, tVec[2, i] * Linv, tVec[3, i] * Linv)
        bVecI = SVector{3, Float64}(bVec[1, i], bVec[2, i], bVec[3, i])
        # Normalised the dislocation network vector, the sum of all the segment vectors has norm 1.
        # Screw component, scalar projection of bVec onto t.
        bScrew = tVecI ⋅ bVecI
        # Edge component, vector rejection of bVec onto t.
        bEdgeVec = bVecI - bScrew * tVecI
        # Finding the norm squared of each edge component.
        bEdgeSq = bEdgeVec ⋅ bEdgeVec
        #=
        A. Arsenlis et al, Modelling Simul. Mater. Sci. Eng. 15 (2007)
        553?595: gives this expression in appendix A p590
        f^{s}_{43} = -(μ/(4π)) [ t × (t × b)](t ⋅ b) { v/(1-v) ( ln[
        (L_a + L)/a] - 2*(L_a - a)/L ) - (L_a - a)^2/(2La*L) }

        tVec × (tVec × bVec)    = tVec (tVec ⋅ bVec) - bVec (tVec ⋅ tVec)
        = tVec * bScrew - bVec
        = - bEdgeVec
        =#
        # Torsional component of the elastic self interaction force. This is the scalar component of the above equation.
        LaMa = La - a
        # Torsional component of core self interaction.
        tor =
            μ4π *
            bScrew *
            (nuOmNuInv * (log((La + L) / a) - 2 * LaMa * Linv) - LaMa^2 / (2 * La * L))
        torCore = 2 * E * nuOmNuInv * bScrew
        torTot = tor + torCore
        # Longitudinal component of core self interaction.
        lonCore = (bScrew^2 + bEdgeSq * omNuInv) * E

        selfForce = torTot * bEdgeVec - lonCore * tVecI

        segForce[1, 1, idx[i]] += -selfForce[1]
        segForce[2, 1, idx[i]] += -selfForce[2]
        segForce[3, 1, idx[i]] += -selfForce[3]
        segForce[1, 2, idx[i]] += selfForce[1]
        segForce[2, 2, idx[i]] += selfForce[2]
        segForce[3, 2, idx[i]] += selfForce[3]
    end

    return network
end

"""
!!! Note
    This function is based on the SegSegForces function by A. Arsenlis et al. It is optimised for speed and reusability. It has also been locally parallelised.
```
```
 It implements the analytical solution of the force between two dislocation segments. Details are found in Appendix A.1. in ["Enabling Strain Hardening Simulations with Dislocation Dynamics" by A. Arsenlis et al.](https://doi.org/10.1088%2F0965-0393%2F15%2F6%2F001)

At a high level this works by creating a local coordinate frame using the line directions of the dislocation segments and a vector orthogonal to them. The line integrals are then evaluated parametrically utilising this local coordinate. BibTex citation here:

@article{Arsenlis_2007,
	doi = {10.1088/0965-0393/15/6/001},
	url = {https://doi.org/10.1088%2F0965-0393%2F15%2F6%2F001},
	year = 2007,
	month = {jul},
	publisher = {{IOP} Publishing},
	volume = {15},
	number = {6},
	pages = {553--595},
	author = {A Arsenlis and W Cai and M Tang and M Rhee and T Oppelstrup and G Hommes and T G Pierce and V V Bulatov},
	title = {Enabling strain hardening simulations with dislocation dynamics},
	journal = {Modelling and Simulation in Materials Science and Engineering},
	abstract = {Numerical algorithms for discrete dislocation dynamics simulations are investigated for the purpose of enabling strain hardening simulations of single crystals on massively parallel computers. The algorithms investigated include the calculation of forces, the equations of motion, time integration, adaptive mesh refinement, the treatment of dislocation core reactions and the dynamic distribution of data and work on parallel computers. A simulation integrating all these algorithmic elements using the Parallel Dislocation Simulator (ParaDiS) code is performed to understand their behaviour in concert and to evaluate the overall numerical performance of dislocation dynamics simulations and their ability to accumulate percent of plastic strain.}
}
"""
@inline function calcSegSegForce(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    idx = nothing;
    parallel::Bool = false,
)

    # Constants.
    μ = matParams.μ
    μ4π = matParams.μ4π
    μ8π = matParams.μ8π
    μ4πν = matParams.μ4πν
    aSq = dlnParams.coreRadSq
    μ8πaSq = aSq * μ8π
    μ4πνaSq = aSq * μ4πν

    bVec = network.bVec
    coord = network.coord
    segIdx = network.segIdx

    # Un normalised segment vectors. Views for speed.
    numSeg = network.numSeg
    idxBvec = @view segIdx[1:numSeg, 1]
    idxNode1 = @view segIdx[1:numSeg, 2]
    idxNode2 = @view segIdx[1:numSeg, 3]
    bVec = @view bVec[:, idxBvec]
    node1 = @view coord[:, idxNode1]
    node2 = @view coord[:, idxNode2]

    # Calculate segseg forces on every segment.
    if isnothing(idx)
        segSegForce = zeros(3, 2, numSeg)
        if parallel
            # Threadid parallelisation + parallelised reduction.
            TSegSegForce = zeros(Threads.nthreads(), 3, 2, numSeg)
            @fastmath @inbounds @sync for i in 1:numSeg
                Threads.@spawn begin
                    b1 = SVector{3, Float64}(bVec[1, i], bVec[2, i], bVec[3, i])
                    n11 = SVector{3, Float64}(node1[1, i], node1[2, i], node1[3, i])
                    n12 = SVector{3, Float64}(node2[1, i], node2[2, i], node2[3, i])
                    @simd for j in (i + 1):numSeg
                        b2 = SVector{3, Float64}(bVec[1, j], bVec[2, j], bVec[3, j])
                        n21 = SVector{3, Float64}(node1[1, j], node1[2, j], node1[3, j])
                        n22 = SVector{3, Float64}(node2[1, j], node2[2, j], node2[3, j])

                        Fnode1, Fnode2, Fnode3, Fnode4 = calcSegSegForce(
                            aSq,
                            μ4π,
                            μ8π,
                            μ8πaSq,
                            μ4πν,
                            μ4πνaSq,
                            b1,
                            n11,
                            n12,
                            b2,
                            n21,
                            n22,
                        )

                        TSegSegForce[Threads.threadid(), 1, 1, i] += Fnode1[1]
                        TSegSegForce[Threads.threadid(), 2, 1, i] += Fnode1[2]
                        TSegSegForce[Threads.threadid(), 3, 1, i] += Fnode1[3]

                        TSegSegForce[Threads.threadid(), 1, 1, j] += Fnode3[1]
                        TSegSegForce[Threads.threadid(), 2, 1, j] += Fnode3[2]
                        TSegSegForce[Threads.threadid(), 3, 1, j] += Fnode3[3]

                        TSegSegForce[Threads.threadid(), 1, 2, i] += Fnode2[1]
                        TSegSegForce[Threads.threadid(), 2, 2, i] += Fnode2[2]
                        TSegSegForce[Threads.threadid(), 3, 2, i] += Fnode2[3]

                        TSegSegForce[Threads.threadid(), 1, 2, j] += Fnode4[1]
                        TSegSegForce[Threads.threadid(), 2, 2, j] += Fnode4[2]
                        TSegSegForce[Threads.threadid(), 3, 2, j] += Fnode4[3]
                    end
                end
            end

            nthreads = Threads.nthreads()
            TSegSegForce2 = [Threads.Atomic{Float64}(0.0) for i in 1:(numSeg * 3 * 2)]
            TSegSegForce2 = reshape(TSegSegForce2, 3, 2, numSeg)
            @fastmath @inbounds @sync for tid in 1:nthreads
                Threads.@spawn begin
                    start = 1 + ((tid - 1) * numSeg) ÷ nthreads
                    stop = (tid * numSeg) ÷ nthreads
                    domain = start:stop
                    Threads.atomic_add!.(
                        TSegSegForce2[:, :, start:stop],
                        sum(TSegSegForce[:, :, :, start:stop], dims = 1)[1, :, :, :],
                    )
                end
            end
            # This allows type inference and reduces memory allocation.
            segSegForce .= getproperty.(TSegSegForce2, :value)
        elseif parallel == false
            # Serial execution.
            @fastmath @inbounds for i in 1:numSeg
                b1 = SVector{3, Float64}(bVec[1, i], bVec[2, i], bVec[3, i])
                n11 = SVector{3, Float64}(node1[1, i], node1[2, i], node1[3, i])
                n12 = SVector{3, Float64}(node2[1, i], node2[2, i], node2[3, i])
                @simd for j in (i + 1):numSeg
                    b2 = SVector{3, Float64}(bVec[1, j], bVec[2, j], bVec[3, j])
                    n21 = SVector{3, Float64}(node1[1, j], node1[2, j], node1[3, j])
                    n22 = SVector{3, Float64}(node2[1, j], node2[2, j], node2[3, j])

                    Fnode1, Fnode2, Fnode3, Fnode4 = calcSegSegForce(
                        aSq,
                        μ4π,
                        μ8π,
                        μ8πaSq,
                        μ4πν,
                        μ4πνaSq,
                        b1,
                        n11,
                        n12,
                        b2,
                        n21,
                        n22,
                    )
                    segSegForce[1, 1, i] += Fnode1[1]
                    segSegForce[2, 1, i] += Fnode1[2]
                    segSegForce[3, 1, i] += Fnode1[3]

                    segSegForce[1, 1, j] += Fnode3[1]
                    segSegForce[2, 1, j] += Fnode3[2]
                    segSegForce[3, 1, j] += Fnode3[3]

                    segSegForce[1, 2, i] += Fnode2[1]
                    segSegForce[2, 2, i] += Fnode2[2]
                    segSegForce[3, 2, i] += Fnode2[3]

                    segSegForce[1, 2, j] += Fnode4[1]
                    segSegForce[2, 2, j] += Fnode4[2]
                    segSegForce[3, 2, j] += Fnode4[3]
                end
            end
        end
    else # Calculate segseg forces only on segments provided
        lenIdx = length(idx)
        segSegForce = zeros(3, 2, lenIdx)
        @fastmath @inbounds for (k, i) in enumerate(idx)
            b1 = SVector{3, Float64}(bVec[1, i], bVec[2, i], bVec[3, i])
            n11 = SVector{3, Float64}(node1[1, i], node1[2, i], node1[3, i])
            n12 = SVector{3, Float64}(node2[1, i], node2[2, i], node2[3, i])
            for j in 1:numSeg
                i == j ? continue : nothing
                b2 = SVector{3, Float64}(bVec[1, j], bVec[2, j], bVec[3, j])
                n21 = SVector{3, Float64}(node1[1, j], node1[2, j], node1[3, j])
                n22 = SVector{3, Float64}(node2[1, j], node2[2, j], node2[3, j])

                Fnode1, Fnode2, missing, missing = calcSegSegForce(
                    aSq,
                    μ4π,
                    μ8π,
                    μ8πaSq,
                    μ4πν,
                    μ4πνaSq,
                    b1,
                    n11,
                    n12,
                    b2,
                    n21,
                    n22,
                )

                segSegForce[1, 1, k] += Fnode1[1]
                segSegForce[2, 1, k] += Fnode1[2]
                segSegForce[3, 1, k] += Fnode1[3]

                segSegForce[1, 2, k] += Fnode2[1]
                segSegForce[2, 2, k] += Fnode2[2]
                segSegForce[3, 2, k] += Fnode2[3]
            end
        end
    end
    #=
        ## Atomic add parallelisation, slow as heck on a single processor.
        TSegSegForce = [Threads.Atomic{Float64}(0.0) for i = 1:(numSeg * 3 * 2)]
        TSegSegForce = reshape(TSegSegForce, numSeg, 3, 2)
        nthreads = Base.Threads.nthreads()
        @fastmath @inbounds Threads.@threads for tid = 1:nthreads
            start = 1 + ((tid - 1) * numSeg) ÷ nthreads
            stop = (tid * numSeg) ÷ nthreads
            domain = start:stop
            for i = start:stop
                b1 = (bVec[1, i], bVec[2, i], bVec[3, i])
                n11 = (node1[1, i], node1[2, i], node1[3, i])
                n12 = (node2[1, i], node2[2, i], node2[3, i])
                for j = (i + 1):numSeg
                    b2 = (bVec[1 ,j], bVec[2 ,j], bVec[3 ,j])
                    n21 = (node1[1 ,j], node1[2 ,j], node1[3 ,j])
                    n22 = (node2[1 ,j], node2[2 ,j], node2[3 ,j])

                    Fnode1, Fnode2, Fnode3, Fnode4 = calcSegSegForce(
                        aSq,
                        μ4π,
                        μ8π,
                        μ8πaSq,
                        μ4πν,
                        μ4πνaSq,
                        b1,
                        n11,
                        n12,
                        b2,
                        n21,
                        n22,
                    )
                    Threads.atomic_add!.(TSegSegForce[i, :, 1], Fnode1)
                    Threads.atomic_add!.(TSegSegForce[j, :, 1], Fnode3)
                    Threads.atomic_add!.(TSegSegForce[i, :, 2], Fnode2)
                    Threads.atomic_add!.(TSegSegForce[j, :, 2], Fnode4)
                end
            end
        end
        segSegForce = getproperty.(TSegSegForce, :value)
    =#
    return segSegForce
end

@inline function calcSegSegForce!(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
    idx = nothing;
    parallel::Bool = true,
)

    # Constants.
    μ = matParams.μ
    μ4π = matParams.μ4π
    μ8π = matParams.μ8π
    μ4πν = matParams.μ4πν
    aSq = dlnParams.coreRadSq
    μ8πaSq = aSq * μ8π
    μ4πνaSq = aSq * μ4πν

    bVec = network.bVec
    coord = network.coord
    segIdx = network.segIdx

    # Un normalised segment vectors. Views for speed.
    numSeg = network.numSeg
    idxBvec = @view segIdx[1:numSeg, 1]
    idxNode1 = @view segIdx[1:numSeg, 2]
    idxNode2 = @view segIdx[1:numSeg, 3]
    bVec = @view bVec[:, idxBvec]
    node1 = @view coord[:, idxNode1]
    node2 = @view coord[:, idxNode2]
    segForce = network.segForce

    # Calculate segseg forces on every segment.
    if isnothing(idx)
        if parallel
            # Threadid parallelisation + parallelised reduction.
            TSegSegForce = zeros(Threads.nthreads(), 3, 2, numSeg)
            @fastmath @inbounds @sync for i in 1:numSeg
                Threads.@spawn begin
                    b1 = SVector{3, Float64}(bVec[1, i], bVec[2, i], bVec[3, i])
                    n11 = SVector{3, Float64}(node1[1, i], node1[2, i], node1[3, i])
                    n12 = SVector{3, Float64}(node2[1, i], node2[2, i], node2[3, i])
                    @simd for j in (i + 1):numSeg
                        b2 = SVector{3, Float64}(bVec[1, j], bVec[2, j], bVec[3, j])
                        n21 = SVector{3, Float64}(node1[1, j], node1[2, j], node1[3, j])
                        n22 = SVector{3, Float64}(node2[1, j], node2[2, j], node2[3, j])

                        Fnode1, Fnode2, Fnode3, Fnode4 = calcSegSegForce(
                            aSq,
                            μ4π,
                            μ8π,
                            μ8πaSq,
                            μ4πν,
                            μ4πνaSq,
                            b1,
                            n11,
                            n12,
                            b2,
                            n21,
                            n22,
                        )

                        TSegSegForce[Threads.threadid(), 1, 1, i] += Fnode1[1]
                        TSegSegForce[Threads.threadid(), 2, 1, i] += Fnode1[2]
                        TSegSegForce[Threads.threadid(), 3, 1, i] += Fnode1[3]

                        TSegSegForce[Threads.threadid(), 1, 1, j] += Fnode3[1]
                        TSegSegForce[Threads.threadid(), 2, 1, j] += Fnode3[2]
                        TSegSegForce[Threads.threadid(), 3, 1, j] += Fnode3[3]

                        TSegSegForce[Threads.threadid(), 1, 2, i] += Fnode2[1]
                        TSegSegForce[Threads.threadid(), 2, 2, i] += Fnode2[2]
                        TSegSegForce[Threads.threadid(), 3, 2, i] += Fnode2[3]

                        TSegSegForce[Threads.threadid(), 1, 2, j] += Fnode4[1]
                        TSegSegForce[Threads.threadid(), 2, 2, j] += Fnode4[2]
                        TSegSegForce[Threads.threadid(), 3, 2, j] += Fnode4[3]
                    end
                end
            end

            nthreads = Threads.nthreads()
            TSegSegForce2 = [Threads.Atomic{Float64}(0.0) for i in 1:(numSeg * 3 * 2)]
            TSegSegForce2 = reshape(TSegSegForce2, 3, 2, numSeg)
            @fastmath @inbounds @sync for tid in 1:nthreads
                Threads.@spawn begin
                    start = 1 + ((tid - 1) * numSeg) ÷ nthreads
                    stop = (tid * numSeg) ÷ nthreads
                    domain = start:stop
                    Threads.atomic_add!.(
                        TSegSegForce2[:, :, start:stop],
                        sum(TSegSegForce[:, :, :, start:stop], dims = 1)[1, :, :, :],
                    )
                end
            end
            # This allows type inference and reduces memory allocation.
            segForce[:, :, 1:numSeg] += getproperty.(TSegSegForce2, :value)
        else
            # Serial execution.
            @fastmath @inbounds for i in 1:numSeg
                b1 = SVector{3, Float64}(bVec[1, i], bVec[2, i], bVec[3, i])
                n11 = SVector{3, Float64}(node1[1, i], node1[2, i], node1[3, i])
                n12 = SVector{3, Float64}(node2[1, i], node2[2, i], node2[3, i])
                @simd for j in (i + 1):numSeg
                    b2 = SVector{3, Float64}(bVec[1, j], bVec[2, j], bVec[3, j])
                    n21 = SVector{3, Float64}(node1[1, j], node1[2, j], node1[3, j])
                    n22 = SVector{3, Float64}(node2[1, j], node2[2, j], node2[3, j])

                    Fnode1, Fnode2, Fnode3, Fnode4 = calcSegSegForce(
                        aSq,
                        μ4π,
                        μ8π,
                        μ8πaSq,
                        μ4πν,
                        μ4πνaSq,
                        b1,
                        n11,
                        n12,
                        b2,
                        n21,
                        n22,
                    )
                    segForce[1, 1, i] += Fnode1[1]
                    segForce[2, 1, i] += Fnode1[2]
                    segForce[3, 1, i] += Fnode1[3]

                    segForce[1, 1, j] += Fnode3[1]
                    segForce[2, 1, j] += Fnode3[2]
                    segForce[3, 1, j] += Fnode3[3]

                    segForce[1, 2, i] += Fnode2[1]
                    segForce[2, 2, i] += Fnode2[2]
                    segForce[3, 2, i] += Fnode2[3]

                    segForce[1, 2, j] += Fnode4[1]
                    segForce[2, 2, j] += Fnode4[2]
                    segForce[3, 2, j] += Fnode4[3]
                end
            end
        end
    else # Calculate segseg forces only on segments provided
        @fastmath @inbounds for i in idx
            b1 = SVector{3, Float64}(bVec[1, i], bVec[2, i], bVec[3, i])
            n11 = SVector{3, Float64}(node1[1, i], node1[2, i], node1[3, i])
            n12 = SVector{3, Float64}(node2[1, i], node2[2, i], node2[3, i])
            for j in 1:numSeg
                i == j ? continue : nothing
                b2 = SVector{3, Float64}(bVec[1, j], bVec[2, j], bVec[3, j])
                n21 = SVector{3, Float64}(node1[1, j], node1[2, j], node1[3, j])
                n22 = SVector{3, Float64}(node2[1, j], node2[2, j], node2[3, j])

                Fnode1, Fnode2, missing, missing = calcSegSegForce(
                    aSq,
                    μ4π,
                    μ8π,
                    μ8πaSq,
                    μ4πν,
                    μ4πνaSq,
                    b1,
                    n11,
                    n12,
                    b2,
                    n21,
                    n22,
                )

                segForce[1, 1, i] += Fnode1[1]
                segForce[2, 1, i] += Fnode1[2]
                segForce[3, 1, i] += Fnode1[3]

                segForce[1, 2, i] += Fnode2[1]
                segForce[2, 2, i] += Fnode2[2]
                segForce[3, 2, i] += Fnode2[3]
            end
        end
    end
    return network
end
@inline function calcSegSegForce(
    aSq::T1,
    μ4π::T1,
    μ8π::T1,
    μ8πaSq::T1,
    μ4πν::T1,
    μ4πνaSq::T1,
    b1::T2,
    n11::T2,
    n12::T2,
    b2::T2,
    n21::T2,
    n22::T2,
) where {T1, T2 <: AbstractVector{T} where {T}}

    t2 = n22 - n21
    t2N = 1 / norm(t2)
    t2 = t2 * t2N

    t1 = n12 - n11
    t1N = 1 / norm(t1)
    t1 = t1 * t1N

    c = t1 ⋅ t2
    cSq = c * c
    omcSq = 1 - cSq

    if omcSq > sqrt(eps(typeof(omcSq)))

        omcSqI = 1 / omcSq

        # Single cross products.
        t2ct1 = t2 × t1
        t1ct2 = -t2ct1
        b2ct2 = b2 × t2
        b1ct1 = b1 × t1

        # Dot products.
        t2db2 = t2 ⋅ b2
        t2db1 = t2 ⋅ b1
        t1db2 = t1 ⋅ b2
        t1db1 = t1 ⋅ b1

        # Cross dot products.
        t2ct1db2 = t2ct1 ⋅ b2
        t1ct2db1 = t1ct2 ⋅ b1
        b1ct1db2 = b1ct1 ⋅ b2
        b2ct2db1 = b2ct2 ⋅ b1

        # Double cross products.
        t2ct1ct2 = t1 - c * t2
        t1ct2ct1 = t2 - c * t1
        t2cb1ct2 = b1 - t2db1 * t2
        t1cb2ct1 = b2 - t1db2 * t1
        b1ct1ct2 = t2db1 * t1 - c * b1
        b2ct2ct1 = t1db2 * t2 - c * b2

        # Double cross product dot product.
        t2ct1cb1dt1 = t2db1 - t1db1 * c
        t1ct2cb2dt2 = t1db2 - t2db2 * c
        t2ct1cb1db2 = t2db1 * t1db2 - t1db1 * t2db2

        # Integration limits for local coordinates.
        R1 = n21 - n11
        R2 = n22 - n12
        d = (R2 ⋅ t2ct1) * omcSqI

        μ4πd = μ4π * d
        μ4πνd = μ4πν * d
        μ4πνdSq = μ4πνd * d
        μ4πνdCu = μ4πνdSq * d
        μ4πνaSqd = μ4πνaSq * d
        μ8πaSqd = μ8πaSq * d

        lim11 = R1 ⋅ t1
        lim12 = R1 ⋅ t2
        lim21 = R2 ⋅ t1
        lim22 = R2 ⋅ t2

        x1 = (lim12 - c * lim11) * omcSqI
        x2 = (lim22 - c * lim21) * omcSqI
        y1 = (lim11 - c * lim12) * omcSqI
        y2 = (lim21 - c * lim22) * omcSqI

        integ = SegSegInteg(aSq, d, c, cSq, omcSq, omcSqI, x1, y1)
        integ = integ .- SegSegInteg(aSq, d, c, cSq, omcSq, omcSqI, x1, y2)
        integ = integ .- SegSegInteg(aSq, d, c, cSq, omcSq, omcSqI, x2, y1)
        integ = integ .+ SegSegInteg(aSq, d, c, cSq, omcSq, omcSqI, x2, y2)

        # Seg 1, nodes 2-1
        tmp1 = t1db2 * t2db1 + t2ct1cb1db2
        V1 = tmp1 * t2ct1
        V2 = b1ct1 * t1ct2cb2dt2
        V3 = t1ct2 * b1ct1db2 - t1cb2ct1 * t2db1
        V4 = -b1ct1 * t2ct1db2
        V5 = b2ct2ct1 * t2db1 - t1ct2 * b2ct2db1

        tmp1 = μ4πνd * t1ct2db1
        tmp2 = μ4πνd * b2ct2db1
        V7 = μ4πd * V1 - μ4πνd * V2 + tmp1 * b2ct2ct1 + tmp2 * t1ct2ct1

        tmp1 = μ4πν * t2db1
        tmp2 = μ4πν * b2ct2db1
        V8 = μ4π * V5 - tmp1 * b2ct2ct1 + tmp2 * t1ct2

        tmp1 = μ4πν * t1db1
        V9 = -tmp1 * b2ct2ct1 + μ4π * V3 - μ4πν * V4

        tmp1 = μ4πνdCu * t1ct2cb2dt2 * t1ct2db1
        V10 = μ8πaSqd * V1 - μ4πνaSqd * V2 - tmp1 * t1ct2ct1

        tmp1 = μ4πνdSq * t1ct2cb2dt2 * t2db1
        tmp2 = μ4πνdSq * t1ct2cb2dt2 * t1ct2db1
        V11 = μ8πaSq * V5 + tmp1 * t1ct2ct1 - tmp2 * t1ct2

        tmp1 = μ4πνdSq * (t2ct1db2 * t1ct2db1 + t1ct2cb2dt2 * t1db1)
        V12 = μ8πaSq * V3 - μ4πνaSq * V4 + tmp1 * t1ct2ct1

        tmp1 = μ4πνd * (t1ct2cb2dt2 * t1db1 + t2ct1db2 * t1ct2db1)
        tmp2 = μ4πνd * t2ct1db2 * t2db1
        V13 = tmp1 * t1ct2 - tmp2 * t1ct2ct1

        tmp1 = μ4πνd * t1ct2cb2dt2 * t2db1
        V14 = tmp1 * t1ct2

        tmp1 = μ4πνd * t2ct1db2 * t1db1
        V15 = -tmp1 * t1ct2ct1

        tmpVec1 = μ4πν * t2ct1db2 * t1ct2
        V16 = -tmpVec1 * t2db1
        V17 = -tmpVec1 * t1db1

        Fint1 = integ[3] - y2 * integ[1]
        Fint2 = integ[4] - y2 * integ[2]
        Fint3 = integ[6] - y2 * integ[3]
        Fint4 = integ[9] - y2 * integ[7]
        Fint5 = integ[10] - y2 * integ[8]
        Fint6 = integ[12] - y2 * integ[9]
        Fint7 = integ[14] - y2 * integ[10]
        Fint8 = integ[13] - y2 * integ[11]
        Fint9 = integ[17] - y2 * integ[12]
        Fint10 = integ[15] - y2 * integ[13]
        Fint11 = integ[19] - y2 * integ[14]

        Fnode1 =
            (
                V7 * Fint1 +
                V8 * Fint2 +
                V9 * Fint3 +
                V10 * Fint4 +
                V11 * Fint5 +
                V12 * Fint6 +
                V13 * Fint7 +
                V14 * Fint8 +
                V15 * Fint9 +
                V16 * Fint10 +
                V17 * Fint11
            ) * t1N

        Fint1 = y1 * integ[1] - integ[3]
        Fint2 = y1 * integ[2] - integ[4]
        Fint3 = y1 * integ[3] - integ[6]
        Fint4 = y1 * integ[7] - integ[9]
        Fint5 = y1 * integ[8] - integ[10]
        Fint6 = y1 * integ[9] - integ[12]
        Fint7 = y1 * integ[10] - integ[14]
        Fint8 = y1 * integ[11] - integ[13]
        Fint9 = y1 * integ[12] - integ[17]
        Fint10 = y1 * integ[13] - integ[15]
        Fint11 = y1 * integ[14] - integ[19]

        Fnode2 =
            (
                V7 * Fint1 +
                V8 * Fint2 +
                V9 * Fint3 +
                V10 * Fint4 +
                V11 * Fint5 +
                V12 * Fint6 +
                V13 * Fint7 +
                V14 * Fint8 +
                V15 * Fint9 +
                V16 * Fint10 +
                V17 * Fint11
            ) * t1N

        # Seg 2 (nodes 4-3)
        tmp1 = t2db1 * t1db2 + t2ct1cb1db2
        V1 = tmp1 * t1ct2
        V2 = b2ct2 * t2ct1cb1dt1
        V3 = t2ct1 * b1ct1db2 - b1ct1ct2 * t1db2
        V5 = t2cb1ct2 * t1db2 - t2ct1 * b2ct2db1
        V6 = b2ct2 * t1ct2db1

        tmp1 = μ4πνd * t2ct1db2
        tmp2 = μ4πνd * b1ct1db2
        V7 = μ4πd * V1 - μ4πνd * V2 + tmp1 * b1ct1ct2 + tmp2 * t2ct1ct2

        tmp1 = μ4πν * t2db2
        V8 = tmp1 * b1ct1ct2 + μ4π * V5 - μ4πν * V6

        tmp1 = μ4πν * t1db2
        tmp2 = μ4πν * b1ct1db2
        V9 = μ4π * V3 + tmp1 * b1ct1ct2 - tmp2 * t2ct1

        tmp1 = μ4πνdCu * t2ct1cb1dt1 * t2ct1db2
        V10 = μ8πaSqd * V1 - μ4πνaSqd * V2 - tmp1 * t2ct1ct2

        tmp1 = μ4πνdSq * (t1ct2db1 * t2ct1db2 + t2ct1cb1dt1 * t2db2)
        V11 = μ8πaSq * V5 - μ4πνaSq * V6 - tmp1 * t2ct1ct2

        tmp1 = μ4πνdSq * t2ct1cb1dt1 * t1db2
        tmp2 = μ4πνdSq * t2ct1cb1dt1 * t2ct1db2
        V12 = μ8πaSq * V3 - tmp1 * t2ct1ct2 + tmp2 * t2ct1

        tmp1 = μ4πνd * (t2ct1cb1dt1 * t2db2 + t1ct2db1 * t2ct1db2)
        tmp2 = μ4πνd * t1ct2db1 * t1db2
        V13 = tmp1 * t2ct1 - tmp2 * t2ct1ct2

        tmp1 = μ4πνd * t1ct2db1 * t2db2
        V14 = -tmp1 * t2ct1ct2

        tmp1 = μ4πνd * t2ct1cb1dt1 * t1db2
        V15 = tmp1 * t2ct1

        tmpVec1 = μ4πν * t1ct2db1 * t2ct1
        V16 = tmpVec1 * t2db2
        V17 = tmpVec1 * t1db2

        Fint1 = x2 * integ[1] - integ[2]
        Fint2 = x2 * integ[2] - integ[5]
        Fint3 = x2 * integ[3] - integ[4]
        Fint4 = x2 * integ[7] - integ[8]
        Fint5 = x2 * integ[8] - integ[11]
        Fint6 = x2 * integ[9] - integ[10]
        Fint7 = x2 * integ[10] - integ[13]
        Fint8 = x2 * integ[11] - integ[16]
        Fint9 = x2 * integ[12] - integ[14]
        Fint10 = x2 * integ[13] - integ[18]
        Fint11 = x2 * integ[14] - integ[15]

        Fnode3 =
            (
                V7 * Fint1 +
                V8 * Fint2 +
                V9 * Fint3 +
                V10 * Fint4 +
                V11 * Fint5 +
                V12 * Fint6 +
                V13 * Fint7 +
                V14 * Fint8 +
                V15 * Fint9 +
                V16 * Fint10 +
                V17 * Fint11
            ) * t2N

        Fint1 = integ[2] - x1 * integ[1]
        Fint2 = integ[5] - x1 * integ[2]
        Fint3 = integ[4] - x1 * integ[3]
        Fint4 = integ[8] - x1 * integ[7]
        Fint5 = integ[11] - x1 * integ[8]
        Fint6 = integ[10] - x1 * integ[9]
        Fint7 = integ[13] - x1 * integ[10]
        Fint8 = integ[16] - x1 * integ[11]
        Fint9 = integ[14] - x1 * integ[12]
        Fint10 = integ[18] - x1 * integ[13]
        Fint11 = integ[15] - x1 * integ[14]

        Fnode4 =
            (
                V7 * Fint1 +
                V8 * Fint2 +
                V9 * Fint3 +
                V10 * Fint4 +
                V11 * Fint5 +
                V12 * Fint6 +
                V13 * Fint7 +
                V14 * Fint8 +
                V15 * Fint9 +
                V16 * Fint10 +
                V17 * Fint11
            ) * t2N

    else
        Fnode1, Fnode2, Fnode3, Fnode4 = calcParSegSegForce(
            aSq,
            μ4π,
            μ8π,
            μ8πaSq,
            μ4πν,
            μ4πνaSq,
            b1,
            n11,
            n12,
            b2,
            n21,
            n22,
        )
    end

    return Fnode1, Fnode2, Fnode3, Fnode4
end

@inline function calcParSegSegForce(
    aSq::T1,
    μ4π::T1,
    μ8π::T1,
    μ8πaSq::T1,
    μ4πν::T1,
    μ4πνaSq::T1,
    b1::T2,
    n11::T2,
    n12::T2,
    b2::T2,
    n21::T2,
    n22::T2,
) where {T1, T2 <: AbstractVector{T} where {T}}

    flip::Bool = false

    t2 = n22 - n21
    t2N = 1 / norm(t2)
    t2 = t2 * t2N

    t1 = n12 - n11
    t1N = 1 / norm(t1)
    t1 = t1 * t1N

    c = t2 ⋅ t1

    # half of the cotangent of critical θ
    hCotanθc = sqrt((1 - sqrt(eps(typeof(c))) * 1.01) / (sqrt(eps(typeof(c))) * 1.01)) / 2

    # If c is negative we do a swap of n11 and n12 to keep notation consistent and avoid
    if c < 0
        flip = true
        n12, n11 = n11, n12
        t1 = -t1
        b1 = -b1
    end

    # Vector projection and rejection.
    tmp = (n22 - n21) ⋅ t1
    n22m = n21 + tmp * t1
    diff = n22 - n22m
    magDiff = norm(diff)

    tmpVec1 = 0.5 * diff
    tmpVec2 = hCotanθc * magDiff * t1
    n21m = n21 + tmpVec1 + tmpVec2
    n22m = n22m + tmpVec1 - tmpVec2

    # Dot products.
    R = n21m - n11
    Rdt1 = R ⋅ t1

    nd = R - Rdt1 * t1
    ndb1 = nd ⋅ b1
    dSq = nd ⋅ nd
    aSq_dSq = aSq + dSq
    aSq_dSqI = 1 / aSq_dSq

    x1 = n21m ⋅ t1
    x2 = n22m ⋅ t1
    y1 = -n11 ⋅ t1
    y2 = -n12 ⋅ t1

    t1db2 = t1 ⋅ b2
    t1db1 = t1 ⋅ b1
    nddb1 = nd ⋅ b1

    # Cross products.
    b2ct1 = b2 × t1
    b1ct1 = b1 × t1
    ndct1 = nd × t1

    # Cross dot products
    b2ct1db1 = b2ct1 ⋅ b1
    b2ct1dnd = b2ct1 ⋅ nd

    # Double cross products
    b2ct1ct1 = t1db2 * t1 - b2

    integ = ParSegSegInteg(aSq_dSq, aSq_dSqI, x1, y1)
    integ = integ .- ParSegSegInteg(aSq_dSq, aSq_dSqI, x1, y2)
    integ = integ .- ParSegSegInteg(aSq_dSq, aSq_dSqI, x2, y1)
    integ = integ .+ ParSegSegInteg(aSq_dSq, aSq_dSqI, x2, y2)

    tmp = t1db1 * t1db2
    tmpVec1 = tmp * nd
    tmpVec2 = b2ct1dnd * b1ct1
    V1 = μ4πν * (nddb1 * b2ct1ct1 + b2ct1db1 * ndct1 - tmpVec2) - μ4π * tmpVec1

    tmp = (μ4πν - μ4π) * t1db1
    V2 = tmp * b2ct1ct1

    tmp = μ4πν * b2ct1dnd * nddb1
    V3 = -μ8πaSq * tmpVec1 - μ4πνaSq * tmpVec2 - tmp * ndct1

    tmp = μ8πaSq * t1db1
    tmp2 = μ4πν * b2ct1dnd * t1db1
    V4 = -tmp * b2ct1ct1 - tmp2 * ndct1

    # Node 2, n12
    Fint1 = integ[3] - y1 * integ[1]
    Fint2 = integ[6] - y1 * integ[4]
    Fint3 = integ[9] - y1 * integ[7]
    Fint4 = integ[12] - y1 * integ[10]
    Fnode2 = (V1 * Fint1 + V2 * Fint2 + V3 * Fint3 + V4 * Fint4) * t1N

    # Node 1, n11
    Fint1 = y2 * integ[1] - integ[3]
    Fint2 = y2 * integ[4] - integ[6]
    Fint3 = y2 * integ[7] - integ[9]
    Fint4 = y2 * integ[10] - integ[12]
    Fnode1 = (V1 * Fint1 + V2 * Fint2 + V3 * Fint3 + V4 * Fint4) * t1N

    magDiffSq = diff ⋅ diff
    magn21mSq = n21m ⋅ n21m
    magn22mSq = n22m ⋅ n22m

    if magDiffSq > sqrt(eps(typeof(magDiffSq))) * (magn21mSq + magn22mSq)
        nothing, nothing, Fnode1Core, Fnode2Core = calcSegSegForce(
            aSq,
            μ4π,
            μ8π,
            μ8πaSq,
            μ4πν,
            μ4πνaSq,
            b2,
            n21,
            n21m,
            b1,
            n11,
            n12,
        )
        Fnode1 = Fnode1 + Fnode1Core
        Fnode2 = Fnode2 + Fnode2Core

        nothing, nothing, Fnode1Core, Fnode2Core = calcSegSegForce(
            aSq,
            μ4π,
            μ8π,
            μ8πaSq,
            μ4πν,
            μ4πνaSq,
            b2,
            n22m,
            n22,
            b1,
            n11,
            n12,
        )
        Fnode1 = Fnode1 + Fnode1Core
        Fnode2 = Fnode2 + Fnode2Core
    end

    # Segment 2
    # Scalar projection of seg1 (n12-n11) onto t2, not normalised because we need the length.
    tmp = (n12 - n11) ⋅ t2
    # Vector projection of seg 1 to seg 2.
    n12m = n11 + tmp * t2
    # Vector rejection and its magnitude.
    diff = n12 - n12m
    magDiff = norm(diff)

    tmpVec1 = 0.5 * diff
    tmpVec2 = hCotanθc * magDiff * t2
    n11m = n11 + tmpVec1 + tmpVec2
    n12m = n12m + tmpVec1 - tmpVec2

    # Dot products.
    R = n21 - n11m
    Rdt2 = R ⋅ t2

    nd = R - Rdt2 * t2
    dSq = nd ⋅ nd
    aSq_dSq = aSq + dSq
    aSq_dSqI = 1 / aSq_dSq

    x1 = n21 ⋅ t2
    x2 = n22 ⋅ t2
    y1 = -n11m ⋅ t2
    y2 = -n12m ⋅ t2

    t2db2 = t2 ⋅ b2
    t2db1 = t2 ⋅ b1
    nddb2 = nd ⋅ b2

    # Cross products.
    b2ct2 = b2 × t2
    b1ct2 = b1 × t2
    ndct2 = nd × t2

    # Cross dot producs.
    b1ct2db2 = b1ct2 ⋅ b2
    b1ct2dnd = b1ct2 ⋅ nd

    # Double cross products.
    b1ct2ct2 = t2db1 * t2 - b1

    integ = ParSegSegInteg(aSq_dSq, aSq_dSqI, x1, y1)
    integ = integ .- ParSegSegInteg(aSq_dSq, aSq_dSqI, x1, y2)
    integ = integ .- ParSegSegInteg(aSq_dSq, aSq_dSqI, x2, y1)
    integ = integ .+ ParSegSegInteg(aSq_dSq, aSq_dSqI, x2, y2)

    tmp = t2db2 * t2db1
    tmpVec1 = tmp * nd
    tmpVec2 = b1ct2dnd * b2ct2
    V1 = μ4πν * (nddb2 * b1ct2ct2 + b1ct2db2 * ndct2 - tmpVec2) - μ4π * tmpVec1

    tmp = (μ4πν - μ4π) * t2db2
    V2 = tmp * b1ct2ct2

    tmp = μ4πν * b1ct2dnd * nddb2
    V3 = -μ8πaSq * tmpVec1 - μ4πνaSq * tmpVec2 - tmp * ndct2

    tmp = μ8πaSq * t2db2
    tmp2 = μ4πν * b1ct2dnd * t2db2
    V4 = -tmp * b1ct2ct2 - tmp2 * ndct2

    Fint1 = integ[2] - x1 * integ[1]
    Fint2 = integ[5] - x1 * integ[4]
    Fint3 = integ[8] - x1 * integ[7]
    Fint4 = integ[11] - x1 * integ[10]
    Fnode4 = (V1 * Fint1 + V2 * Fint2 + V3 * Fint3 + V4 * Fint4) * t2N

    Fint1 = x2 * integ[1] - integ[2]
    Fint2 = x2 * integ[4] - integ[5]
    Fint3 = x2 * integ[7] - integ[8]
    Fint4 = x2 * integ[10] - integ[11]
    Fnode3 = (V1 * Fint1 + V2 * Fint2 + V3 * Fint3 + V4 * Fint4) * t2N

    magDiffSq = magDiff^2
    magn11mSq = n11m ⋅ n11m
    magn12mSq = n12m ⋅ n12m

    if magDiffSq > sqrt(eps(typeof(magDiffSq))) * (magn11mSq + magn12mSq)
        nothing, nothing, Fnode3Core, Fnode4Core = calcSegSegForce(
            aSq,
            μ4π,
            μ8π,
            μ8πaSq,
            μ4πν,
            μ4πνaSq,
            b1,
            n11,
            n11m,
            b2,
            n21,
            n22,
        )
        Fnode3 = Fnode3 + Fnode3Core
        Fnode4 = Fnode4 + Fnode4Core

        nothing, nothing, Fnode3Core, Fnode4Core = calcSegSegForce(
            aSq,
            μ4π,
            μ8π,
            μ8πaSq,
            μ4πν,
            μ4πνaSq,
            b1,
            n12m,
            n12,
            b2,
            n21,
            n22,
        )
        Fnode3 = Fnode3 + Fnode3Core
        Fnode4 = Fnode4 + Fnode4Core
    end

    # If we flipped the first segment originally, flip the forces round.
    if flip
        Fnode1, Fnode2 = Fnode2, Fnode1
    end

    return Fnode1, Fnode2, Fnode3, Fnode4
end

@inline function ParSegSegInteg(aSq_dSq::T1, aSq_dSqI::T1, x::T1, y::T1) where {T1}

    xpy = x + y
    xmy = x - y
    Ra = sqrt(aSq_dSq + xpy * xpy)
    RaInv = 1 / Ra
    Log_Ra_ypz = log(Ra + xpy)

    tmp = xmy * Ra * aSq_dSqI

    integ1 = Ra * aSq_dSqI
    integ2 = -0.5 * (Log_Ra_ypz - tmp)
    integ3 = -0.5 * (Log_Ra_ypz + tmp)
    integ4 = -Log_Ra_ypz
    integ5 = y * Log_Ra_ypz - Ra
    integ6 = x * Log_Ra_ypz - Ra
    integ7 = aSq_dSqI * (2 * aSq_dSqI * Ra - RaInv)
    integ8 = aSq_dSqI * (tmp - x * RaInv)
    integ9 = -aSq_dSqI * (tmp + y * RaInv)
    integ10 = -aSq_dSqI * xpy * RaInv
    integ11 = RaInv - y * integ10
    integ12 = RaInv - x * integ10

    return integ1,
    integ2,
    integ3,
    integ4,
    integ5,
    integ6,
    integ7,
    integ8,
    integ9,
    integ10,
    integ11,
    integ12
end

@inline function SegSegInteg(
    aSq::T,
    d::T,
    c::T,
    cSq::T,
    omcSq::T,
    omcSqI::T,
    x::T,
    y::T,
) where {T}

    aSq_dSq = aSq + d^2 * omcSq
    xSq = x^2
    ySq = y^2
    Ra = sqrt(aSq_dSq + xSq + ySq + 2 * x * y * c)
    RaInv = 1 / Ra

    Ra_Rd_t1 = Ra + y + x * c
    Ra_Rd_t2 = Ra + x + y * c

    log_Ra_Rd_t1 = log(Ra_Rd_t1)
    xlog_Ra_Rd_t1 = x * log_Ra_Rd_t1

    log_Ra_Rd_t2 = log(Ra_Rd_t2)
    ylog_Ra_Rd_t2 = y * log_Ra_Rd_t2

    RaSq_R_t1_I = RaInv / Ra_Rd_t1
    xRaSq_R_t1_I = x * RaSq_R_t1_I
    xSqRaSq_R_t1_I = x * xRaSq_R_t1_I

    RaSq_R_t2_I = RaInv / Ra_Rd_t2
    yRaSq_R_t2_I = y * RaSq_R_t2_I
    ySqRaSq_R_t2_I = y * yRaSq_R_t2_I

    den = 1 / sqrt(omcSq * aSq_dSq)

    integ1 = -2 * den * atan((1 + c) * (Ra + x + y) * den)

    c_1 = aSq_dSq * integ1
    c_5_6 = (c * Ra - c_1) * omcSqI

    integ2 = (c * log_Ra_Rd_t2 - log_Ra_Rd_t1) * omcSqI
    integ3 = (c * log_Ra_Rd_t1 - log_Ra_Rd_t2) * omcSqI
    integ4 = (c * c_1 - Ra) * omcSqI
    integ5 = ylog_Ra_Rd_t2 + c_5_6
    integ6 = xlog_Ra_Rd_t1 + c_5_6

    c_11_12 = integ1 - c * RaInv
    c_15_18 = c * xRaSq_R_t1_I - RaInv
    x_13_14 = x * c_15_18
    c_19 = c * yRaSq_R_t2_I - RaInv
    y_13_14 = y * c_19
    c_16 = log_Ra_Rd_t2 - (x - c * y) * RaInv - cSq * ySqRaSq_R_t2_I
    z_15_18 = y * c_16
    c_17_19 = log_Ra_Rd_t1 - (y - c * x) * RaInv - cSq * xSqRaSq_R_t1_I

    c15_18_19 = 2 * integ4

    integ7 = (integ1 - xRaSq_R_t1_I - yRaSq_R_t2_I) / (aSq_dSq)
    integ8 = (RaSq_R_t1_I - c * RaSq_R_t2_I) * omcSqI
    integ9 = (RaSq_R_t2_I - c * RaSq_R_t1_I) * omcSqI
    integ10 = (RaInv - c * (xRaSq_R_t1_I + yRaSq_R_t2_I + integ1)) * omcSqI
    integ11 = (xRaSq_R_t1_I + cSq * yRaSq_R_t2_I + c_11_12) * omcSqI
    integ12 = (yRaSq_R_t2_I + cSq * xRaSq_R_t1_I + c_11_12) * omcSqI
    integ13 = (integ3 - x_13_14 + c * (y_13_14 - integ2)) * omcSqI
    integ14 = (integ2 - y_13_14 + c * (x_13_14 - integ3)) * omcSqI
    integ15 = (integ5 - z_15_18 + c * (xSq * c_15_18 - c15_18_19)) * omcSqI
    integ16 = (xSqRaSq_R_t1_I + c * c_16 + 2 * integ2) * omcSqI
    integ17 = (ySqRaSq_R_t2_I + c * c_17_19 + 2 * integ3) * omcSqI
    integ18 = (c15_18_19 - xSq * c_15_18 + c * (z_15_18 - integ5)) * omcSqI
    integ19 = (c15_18_19 - ySq * c_19 + c * (x * c_17_19 - integ6)) * omcSqI

    return integ1,
    integ2,
    integ3,
    integ4,
    integ5,
    integ6,
    integ7,
    integ8,
    integ9,
    integ10,
    integ11,
    integ12,
    integ13,
    integ14,
    integ15,
    integ16,
    integ17,
    integ18,
    integ19
end
