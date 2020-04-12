"""
```
function calcSelfForce(
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
)

    μ = matParams.μ
    ν = matParams.ν
    E = matParams.E
    omNuInv = matParams.omνInv
    nuOmNuInv = matParams.νomνInv
    μ4π = matParams.μ4π
    a = dlnParams.coreRad
    aSq = dlnParams.coreRadSq

    idx = network.segIdx
    coord = network.coord
    numSeg = network.numSeg
    # Un normalised segment vectors.
    bVec = network.bVec[idx[:, 1], :]
    tVec = @. coord[idx[:, 3], :] - coord[idx[:, 2], :]

    # We don't use fused-vectorised operations or dot products because the explicit loop is already 3x faster at 100 dislocations, scales much better in memory and compute time, and can be parallelised more easily. Though the parallelisation overhead isn't worth it unless you are running monstruous simulations.
    L = zeros(numSeg)
    Linv = zeros(numSeg)
    La = zeros(numSeg)
    bScrew = zeros(numSeg)
    bEdgeVec = zeros(numSeg, 3)
    bEdgeSq = zeros(numSeg)
    LaMa = zeros(numSeg)
    tor = zeros(numSeg)
    torCore = zeros(numSeg)
    torTot = zeros(numSeg)
    lonCore = zeros(numSeg)

    @fastmath @inbounds @simd for i = 1:numSeg
        # Finding the norm of each line vector.
        tVecSq = tVec[i, 1]^2 + tVec[i, 2]^2 + tVec[i, 3]^2
        L[i] = sqrt(tVecSq)
        Linv[i] = inv(L[i])
        # Finding the non-singular norm.
        La[i] = sqrt(tVecSq + aSq)
        # Normalised the dislocation network vector, the sum of all the segment vectors has norm 1.
        tVec[i, 1] *= Linv[i]
        tVec[i, 2] *= Linv[i]
        tVec[i, 3] *= Linv[i]
        # Screw component, scalar projection of bVec onto t.
        bScrew[i] =
            bVec[i, 1] * tVec[i, 1] +
            bVec[i, 2] * tVec[i, 2] +
            bVec[i, 3] * tVec[i, 3]
        # Edge component, vector rejection of bVec onto t.
        bEdgeVec[i, 1] = bVec[i, 1] - bScrew[i] * tVec[i, 1]
        bEdgeVec[i, 2] = bVec[i, 2] - bScrew[i] * tVec[i, 2]
        bEdgeVec[i, 3] = bVec[i, 3] - bScrew[i] * tVec[i, 3]
        # Finding the norm squared of each edge component.
        bEdgeSq[i] =
            sqrt(bEdgeVec[i, 1]^2 + bEdgeVec[i, 2]^2 + bEdgeVec[i, 3]^2)
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
        LaMa[i] = La[i] - a
        # Torsional component of core self interaction.
        tor[i] = @. μ4π *
           bScrew[i] *
           (
               nuOmNuInv * (log((La[i] + L[i]) / a) - 2 * LaMa[i] * Linv[i]) -
               LaMa[i]^2 / (2 * La[i] * L[i])
           )
        torCore[i] = 2 * E * nuOmNuInv * bScrew[i]
        torTot[i] = tor[i] + torCore[i]
        # Longitudinal component of core self interaction.
        lonCore[i] = (bScrew[i]^2 + bEdgeSq[i] * omNuInv) * E
    end
    # Force on node 2 = -Force on node 1. This is faster as a fused operation than in the loop. The final assignment uses no memory, julia black magic.
    f2 = @. torTot * bEdgeVec - lonCore * tVec
    f1 = -f2

    return f1, f2
end

"""
!!! Note
    This function is based on the SegSegForces function by A. Arsenlis et al. It is optimised for speed and reusability.
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
    network::DislocationNetwork;
    parallel::Bool = true,
)

    # Constants.
    μ = matParams.μ
    ν = matParams.ν
    μ4π = matParams.μ4π
    μ8π = matParams.μ8π
    μ4πν = matParams.μ4πν
    aSq = dlnParams.coreRadSq
    μ8πaSq = aSq * μ8π
    μ4πνaSq = aSq * μ4πν

    # Un normalised segment vectors.
    numSegs = network.numSeg
    idx = network.segIdx
    coord = network.coord
    bVec = network.bVec[idx[:, 1], :]
    node1 = coord[idx[:, 2], :]
    node2 = coord[idx[:, 3], :]
    SegSegForce = zeros(numSegs, 3, 2)

    if parallel
        # Threadid parallelisation + parallelised reduction.
        # Uses a lot more memory but is faster than the other two.
        TSegSegForce = zeros(Threads.nthreads(), numSegs, 3, 2)
        @fastmath @inbounds Threads.@threads for i = 1:numSegs
            b1 = (bVec[i, 1], bVec[i, 2], bVec[i, 3])
            n11 = (node1[i, 1], node1[i, 2], node1[i, 3])
            n12 = (node2[i, 1], node2[i, 2], node2[i, 3])
            @simd for j = (i + 1):numSegs
                b2 = (bVec[j, 1], bVec[j, 2], bVec[j, 3])
                n21 = (node1[j, 1], node1[j, 2], node1[j, 3])
                n22 = (node2[j, 1], node2[j, 2], node2[j, 3])

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

                TSegSegForce[Threads.threadid(), i, :, 1] .+= Fnode1
                TSegSegForce[Threads.threadid(), j, :, 1] .+= Fnode3
                TSegSegForce[Threads.threadid(), i, :, 2] .+= Fnode2
                TSegSegForce[Threads.threadid(), j, :, 2] .+= Fnode4

            end
        end

        nthreads = Threads.nthreads()
        TSegSegForce2 =
            [Threads.Atomic{Float64}(0.0) for i = 1:(numSegs * 3 * 2)]
        TSegSegForce2 = reshape(TSegSegForce2, numSegs, 3, 2)
        @fastmath @inbounds Threads.@threads for tid = 1:nthreads
            start = 1 + ((tid - 1) * numSegs) ÷ nthreads
            stop = (tid * numSegs) ÷ nthreads
            domain = start:stop
            Threads.atomic_add!.(
                TSegSegForce2[start:stop, :, :],
                sum(TSegSegForce[:, start:stop, :, :], dims = 1)[1, :, :, :],
            )
        end
        # This allows type inference and reduces memory allocation.
        SegSegForce .= getproperty.(TSegSegForce2, :value)
    else
        # Serial execution.
        @fastmath @inbounds for i = 1:numSegs
            b1 = (bVec[i, 1], bVec[i, 2], bVec[i, 3])
            n11 = (node1[i, 1], node1[i, 2], node1[i, 3])
            n12 = (node2[i, 1], node2[i, 2], node2[i, 3])
            @simd for j = (i + 1):numSegs
                b2 = (bVec[j, 1], bVec[j, 2], bVec[j, 3])
                n21 = (node1[j, 1], node1[j, 2], node1[j, 3])
                n22 = (node2[j, 1], node2[j, 2], node2[j, 3])

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
                SegSegForce[i, :, 1] .+= Fnode1
                SegSegForce[j, :, 1] .+= Fnode3
                SegSegForce[i, :, 2] .+= Fnode2
                SegSegForce[j, :, 2] .+= Fnode4
            end
        end
    end

    ### Atomic add parallelisation, slow as heck on a single processor.
    # TSegSegForce = [Threads.Atomic{Float64}(0.0) for i = 1:(numSegs * 3 * 2)]
    # TSegSegForce = reshape(TSegSegForce, numSegs, 3, 2)
    # nthreads = Base.Threads.nthreads()
    # @fastmath @inbounds Threads.@threads for tid = 1:nthreads
    #     start = 1 + ((tid - 1) * numSegs) ÷ nthreads
    #     stop = (tid * numSegs) ÷ nthreads
    #     domain = start:stop
    #     for i = start:stop
    #         b1 = (bVec[i, 1], bVec[i, 2], bVec[i, 3])
    #         n11 = (node1[i, 1], node1[i, 2], node1[i, 3])
    #         n12 = (node2[i, 1], node2[i, 2], node2[i, 3])
    #         for j = (i + 1):numSegs
    #             b2 = (bVec[j, 1], bVec[j, 2], bVec[j, 3])
    #             n21 = (node1[j, 1], node1[j, 2], node1[j, 3])
    #             n22 = (node2[j, 1], node2[j, 2], node2[j, 3])
    #
    #             Fnode1, Fnode2, Fnode3, Fnode4 = calcSegSegForce(
    #                 aSq,
    #                 μ4π,
    #                 μ8π,
    #                 μ8πaSq,
    #                 μ4πν,
    #                 μ4πνaSq,
    #                 b1,
    #                 n11,
    #                 n12,
    #                 b2,
    #                 n21,
    #                 n22,
    #             )
    #             Threads.atomic_add!.(TSegSegForce[i, :, 1], Fnode1)
    #             Threads.atomic_add!.(TSegSegForce[j, :, 1], Fnode3)
    #             Threads.atomic_add!.(TSegSegForce[i, :, 2], Fnode2)
    #             Threads.atomic_add!.(TSegSegForce[j, :, 2], Fnode4)
    #         end
    #     end
    # end
    # SegSegForce = getproperty.(TSegSegForce, :value)

    return SegSegForce
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
) where {T1 <: Float64, T2 <: NTuple{N, <:Float64} where {N}}

    Fnode1::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    Fnode2::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    Fnode3::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    Fnode4::NTuple{3, Float64} = (0.0, 0.0, 0.0)

    t2 = @. n22 - n21
    t2N = 1 / norm(t2)
    t2 = @. t2 * t2N

    t1 = @. n12 - n11
    t1N = 1 / norm(t1)
    t1 = @. t1 * t1N

    c = dot(t1, t2)
    cSq = c * c
    omcSq = 1 - cSq
    if omcSq > eps(Float32)

        omcSqI = 1 / omcSq

        # Single cross products.
        t2ct1 = (
            t2[2] * t1[3] - t2[3] * t1[2],
            t2[3] * t1[1] - t2[1] * t1[3],
            t2[1] * t1[2] - t2[2] * t1[1],
        )

        t1ct2 = .-t2ct1

        b2ct2 = (
            b2[2] * t2[3] - b2[3] * t2[2],
            b2[3] * t2[1] - b2[1] * t2[3],
            b2[1] * t2[2] - b2[2] * t2[1],
        )

        b1ct1 = (
            b1[2] * t1[3] - b1[3] * t1[2],
            b1[3] * t1[1] - b1[1] * t1[3],
            b1[1] * t1[2] - b1[2] * t1[1],
        )

        # Dot products.
        t2db2 = dot(t2, b2)
        t2db1 = dot(t2, b1)
        t1db2 = dot(t1, b2)
        t1db1 = dot(t1, b1)

        # Cross dot products.
        t2ct1db2 = dot(t2ct1, b2)
        t1ct2db1 = dot(t1ct2, b1)
        b1ct1db2 = dot(b1ct1, b2)
        b2ct2db1 = dot(b2ct2, b1)

        # Double cross products.
        t2ct1ct2 = @. t1 - c * t2
        t1ct2ct1 = @. t2 - c * t1
        t2cb1ct2 = @. b1 - t2db1 * t2
        t1cb2ct1 = @. b2 - t1db2 * t1
        b1ct1ct2 = @. t2db1 * t1 - c * b1
        b2ct2ct1 = @. t1db2 * t2 - c * b2

        # Double cross product dot product.
        t2ct1cb1dt1 = t2db1 - t1db1 * c
        t1ct2cb2dt2 = t1db2 - t2db2 * c
        t2ct1cb1db2 = t2db1 * t1db2 - t1db1 * t2db2

        # Integration limits for local coordinates.
        R1 = @. n21 - n11
        R2 = @. n22 - n12
        d = dot(R2, t2ct1) * omcSqI

        μ4πd = μ4π * d
        μ8πd = μ8π * d
        μ4πνd = μ4πν * d
        μ4πνdSq = μ4πνd * d
        μ4πνdCu = μ4πνdSq * d
        μ4πνaSqd = μ4πνaSq * d
        μ8πaSqd = μ8πaSq * d

        lim11 = dot(R1, t1)
        lim12 = dot(R1, t2)
        lim21 = dot(R2, t1)
        lim22 = dot(R2, t2)

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
        V1 = @. tmp1 * t2ct1
        V2 = @. b1ct1 * t1ct2cb2dt2
        V3 = @. t1ct2 * b1ct1db2 - t1cb2ct1 * t2db1
        V4 = @. -b1ct1 * t2ct1db2
        V5 = @. b2ct2ct1 * t2db1 - t1ct2 * b2ct2db1

        tmp1 = μ4πνd * t1ct2db1
        tmp2 = μ4πνd * b2ct2db1
        V7 = @. μ4πd * V1 - μ4πνd * V2 + tmp1 * b2ct2ct1 + tmp2 * t1ct2ct1

        tmp1 = μ4πν * t2db1
        tmp2 = μ4πν * b2ct2db1
        V8 = @. μ4π * V5 - tmp1 * b2ct2ct1 + tmp2 * t1ct2

        tmp1 = μ4πν * t1db1
        V9 = @. -tmp1 * b2ct2ct1 + μ4π * V3 - μ4πν * V4

        tmp1 = μ4πνdCu * t1ct2cb2dt2 * t1ct2db1
        V10 = @. μ8πaSqd * V1 - μ4πνaSqd * V2 - tmp1 * t1ct2ct1

        tmp1 = μ4πνdSq * t1ct2cb2dt2 * t2db1
        tmp2 = μ4πνdSq * t1ct2cb2dt2 * t1ct2db1
        V11 = @. μ8πaSq * V5 + tmp1 * t1ct2ct1 - tmp2 * t1ct2

        tmp1 = μ4πνdSq * (t2ct1db2 * t1ct2db1 + t1ct2cb2dt2 * t1db1)
        V12 = @. μ8πaSq * V3 - μ4πνaSq * V4 + tmp1 * t1ct2ct1

        tmp1 = μ4πνd * (t1ct2cb2dt2 * t1db1 + t2ct1db2 * t1ct2db1)
        tmp2 = μ4πνd * t2ct1db2 * t2db1
        V13 = @. tmp1 * t1ct2 - tmp2 * t1ct2ct1

        tmp1 = μ4πνd * t1ct2cb2dt2 * t2db1
        V14 = @. tmp1 * t1ct2

        tmp1 = μ4πνd * t2ct1db2 * t1db1
        V15 = @. -tmp1 * t1ct2ct1

        tmpVec1 = @. μ4πν * t2ct1db2 * t1ct2
        V16 = @. -tmpVec1 * t2db1
        V17 = @. -tmpVec1 * t1db1

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

        Fnode1 = @. (
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

        Fnode2 = @. (
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
        V1 = @. tmp1 * t1ct2
        V2 = @. b2ct2 * t2ct1cb1dt1
        V3 = @. t2ct1 * b1ct1db2 - b1ct1ct2 * t1db2
        V5 = @. t2cb1ct2 * t1db2 - t2ct1 * b2ct2db1
        V6 = @. b2ct2 * t1ct2db1

        tmp1 = μ4πνd * t2ct1db2
        tmp2 = μ4πνd * b1ct1db2
        V7 = @. μ4πd * V1 - μ4πνd * V2 + tmp1 * b1ct1ct2 + tmp2 * t2ct1ct2

        tmp1 = μ4πν * t2db2
        V8 = @. tmp1 * b1ct1ct2 + μ4π * V5 - μ4πν * V6

        tmp1 = μ4πν * t1db2
        tmp2 = μ4πν * b1ct1db2
        V9 = @. μ4π * V3 + tmp1 * b1ct1ct2 - tmp2 * t2ct1

        tmp1 = μ4πνdCu * t2ct1cb1dt1 * t2ct1db2
        V10 = @. μ8πaSqd * V1 - μ4πνaSqd * V2 - tmp1 * t2ct1ct2

        tmp1 = μ4πνdSq * (t1ct2db1 * t2ct1db2 + t2ct1cb1dt1 * t2db2)
        V11 = @. μ8πaSq * V5 - μ4πνaSq * V6 - tmp1 * t2ct1ct2

        tmp1 = μ4πνdSq * t2ct1cb1dt1 * t1db2
        tmp2 = μ4πνdSq * t2ct1cb1dt1 * t2ct1db2
        V12 = @. μ8πaSq * V3 - tmp1 * t2ct1ct2 + tmp2 * t2ct1

        tmp1 = μ4πνd * (t2ct1cb1dt1 * t2db2 + t1ct2db1 * t2ct1db2)
        tmp2 = μ4πνd * t1ct2db1 * t1db2
        V13 = @. tmp1 * t2ct1 - tmp2 * t2ct1ct2

        tmp1 = μ4πνd * t1ct2db1 * t2db2
        V14 = @. -tmp1 * t2ct1ct2

        tmp1 = μ4πνd * t2ct1cb1dt1 * t1db2
        V15 = @. tmp1 * t2ct1

        tmpVec1 = @. μ4πν * t1ct2db1 * t2ct1
        V16 = @. tmpVec1 * t2db2
        V17 = @. tmpVec1 * t1db2

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

        Fnode3 = @. (
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
        ) .* t2N

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

        Fnode4 = @. (
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
        ) .* t2N

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
) where {T1 <: Float64, T2 <: NTuple{N, <:Float64} where {N}}

    Fnode1::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    Fnode2::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    Fnode3::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    Fnode4::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    flip::Bool = false

    # half of the cotangent of critical θ
    hCotanθc =
        sqrt((1 - sqrt(eps(Float64)) * 1.01) / (sqrt(eps(Float64)) * 1.01)) / 2

    t2 = @. n22 - n21
    t2N = 1 / norm(t2)
    t2 = @. t2 * t2N

    t1 = @. n12 - n11
    t1N = 1 / norm(t1)
    t1 = @. t1 * t1N

    c = dot(t2, t1)
    # If c is negative we do a swap of n11 and n12 to keep notation consistent and avoid
    if c < 0
        flip = true
        n12, n11 = n11, n12
        t1 = .-t1
        b1 = .-b1
    end

    # Vector projection and rejection.
    tmp = dot(n22 .- n21, t1)
    n22m = @. n21 + tmp * t1
    diff = @. n22 - n22m
    magDiff = norm(diff)

    tmpVec1 = @. 0.5 * diff
    tmpVec2 = @. hCotanθc * magDiff * t1
    n21m = @. n21 + tmpVec1 + tmpVec2
    n22m = @. n22m + tmpVec1 - tmpVec2

    # Dot products.
    R = @. n21m - n11
    Rdt1 = dot(R, t1)

    nd = @. R - Rdt1 * t1
    ndb1 = dot(nd, b1)
    dSq = dot(nd, nd)
    aSq_dSq = aSq + dSq
    aSq_dSqI = 1 / aSq_dSq

    x1 = dot(n21m, t1)
    x2 = dot(n22m, t1)
    y1 = -dot(n11, t1)
    y2 = -dot(n12, t1)

    t1db2 = dot(t1, b2)
    t1db1 = dot(t1, b1)
    nddb1 = dot(nd, b1)

    # Cross products.
    b2ct1 = (
        b2[2] * t1[3] - b2[3] * t1[2],
        b2[3] * t1[1] - b2[1] * t1[3],
        b2[1] * t1[2] - b2[2] * t1[1],
    )

    b1ct1 = (
        b1[2] * t1[3] - b1[3] * t1[2],
        b1[3] * t1[1] - b1[1] * t1[3],
        b1[1] * t1[2] - b1[2] * t1[1],
    )

    ndct1 = (
        nd[2] * t1[3] - nd[3] * t1[2],
        nd[3] * t1[1] - nd[1] * t1[3],
        nd[1] * t1[2] - nd[2] * t1[1],
    )

    # Cross dot products
    b2ct1db1 = dot(b2ct1, b1)
    b2ct1dnd = dot(b2ct1, nd)

    # Double cross products
    b2ct1ct1 = @. t1db2 * t1 - b2

    integ = ParSegSegInteg(aSq_dSq, aSq_dSqI, x1, y1)
    integ = integ .- ParSegSegInteg(aSq_dSq, aSq_dSqI, x1, y2)
    integ = integ .- ParSegSegInteg(aSq_dSq, aSq_dSqI, x2, y1)
    integ = integ .+ ParSegSegInteg(aSq_dSq, aSq_dSqI, x2, y2)

    tmp = t1db1 * t1db2
    tmpVec1 = @. tmp * nd
    tmpVec2 = @. b2ct1dnd * b1ct1
    V1 = @. μ4πν * (nddb1 * b2ct1ct1 + b2ct1db1 * ndct1 - tmpVec2) -
            μ4π * tmpVec1

    tmp = (μ4πν - μ4π) * t1db1
    V2 = @. tmp * b2ct1ct1

    tmp = μ4πν * b2ct1dnd * nddb1
    V3 = @. -μ8πaSq * tmpVec1 - μ4πνaSq * tmpVec2 - tmp * ndct1

    tmp = μ8πaSq * t1db1
    tmp2 = μ4πν * b2ct1dnd * t1db1
    V4 = @. -tmp * b2ct1ct1 - tmp2 * ndct1

    # Node 2, n12
    Fint1 = integ[3] - y1 * integ[1]
    Fint2 = integ[6] - y1 * integ[4]
    Fint3 = integ[9] - y1 * integ[7]
    Fint4 = integ[12] - y1 * integ[10]
    Fnode2 = @. (V1 * Fint1 + V2 * Fint2 + V3 * Fint3 + V4 * Fint4) * t1N

    # Node 1, n11
    Fint1 = y2 * integ[1] - integ[3]
    Fint2 = y2 * integ[4] - integ[6]
    Fint3 = y2 * integ[7] - integ[9]
    Fint4 = y2 * integ[10] - integ[12]

    Fnode1 = @. (V1 * Fint1 + V2 * Fint2 + V3 * Fint3 + V4 * Fint4) * t1N

    magDiffSq = dot(diff, diff)
    magn21mSq = dot(n21m, n21m)
    magn22mSq = dot(n22m, n22m)

    if magDiffSq > eps(Float32) * (magn21mSq + magn22mSq)
        missing, missing, Fnode1Core, Fnode2Core = calcSegSegForce(
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
        Fnode1 = @. Fnode1 + Fnode1Core
        Fnode2 = @. Fnode2 + Fnode2Core

        missing, missing, Fnode1Core, Fnode2Core = calcSegSegForce(
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
        Fnode1 = @. Fnode1 + Fnode1Core
        Fnode2 = @. Fnode2 + Fnode2Core
    end

    # Segment 2
    # Scalar projection of seg1 (n12-n11) onto t2, not normalised because we need the length.
    tmp = dot(n12 .- n11, t2)
    # Vector projection of seg 1 to seg 2.
    n12m = @. n11 + tmp * t2
    # Vector rejection and its magnitude.
    diff = @. n12 - n12m
    magDiff = norm(diff)

    tmpVec1 = @. 0.5 * diff
    tmpVec2 = @. hCotanθc * magDiff * t2
    n11m = @. n11 + tmpVec1 + tmpVec2
    n12m = @. n12m + tmpVec1 - tmpVec2

    # Dot products.
    R = @. n21 - n11m
    Rdt2 = dot(R, t2)

    nd = @. R - Rdt2 * t2
    dSq = dot(nd, nd)
    aSq_dSq = aSq + dSq
    aSq_dSqI = 1 / aSq_dSq

    x1 = dot(n21, t2)
    x2 = dot(n22, t2)
    y1 = -dot(n11m, t2)
    y2 = -dot(n12m, t2)

    t2db2 = dot(t2, b2)
    t2db1 = dot(t2, b1)
    nddb2 = dot(nd, b2)

    # Cross products.
    b2ct2 = (
        b2[2] * t2[3] - b2[3] * t2[2],
        b2[3] * t2[1] - b2[1] * t2[3],
        b2[1] * t2[2] - b2[2] * t2[1],
    )

    b1ct2 = (
        b1[2] * t2[3] - b1[3] * t2[2],
        b1[3] * t2[1] - b1[1] * t2[3],
        b1[1] * t2[2] - b1[2] * t2[1],
    )

    ndct2 = (
        nd[2] * t2[3] - nd[3] * t2[2],
        nd[3] * t2[1] - nd[1] * t2[3],
        nd[1] * t2[2] - nd[2] * t2[1],
    )

    # Cross dot producs.
    b1ct2db2 = dot(b1ct2, b2)
    b1ct2dnd = dot(b1ct2, nd)

    # Double cross products.
    b1ct2ct2 = @. t2db1 * t2 - b1

    integ = ParSegSegInteg(aSq_dSq, aSq_dSqI, x1, y1)
    integ = integ .- ParSegSegInteg(aSq_dSq, aSq_dSqI, x1, y2)
    integ = integ .- ParSegSegInteg(aSq_dSq, aSq_dSqI, x2, y1)
    integ = integ .+ ParSegSegInteg(aSq_dSq, aSq_dSqI, x2, y2)

    tmp = t2db2 * t2db1
    tmpVec1 = @. tmp * nd
    tmpVec2 = @. b1ct2dnd * b2ct2
    V1 = @. μ4πν * (nddb2 * b1ct2ct2 + b1ct2db2 * ndct2 - tmpVec2) -
            μ4π * tmpVec2

    tmp = (μ4πν - μ4π) * t2db2
    V2 = @. tmp * b1ct2ct2

    tmp = μ4πν * b1ct2dnd * nddb2
    V3 = @. -μ8πaSq * tmpVec1 - μ4πνaSq * tmpVec2 - tmp * ndct2

    tmp = μ8πaSq * t2db2
    tmp2 = μ4πν * b1ct2dnd * t2db2
    V4 = @. -tmp * b1ct2ct2 - tmp2 * tmp2 * ndct2

    Fint1 = integ[2] - x1 * integ[1]
    Fint2 = integ[5] - x1 * integ[4]
    Fint3 = integ[8] - x1 * integ[7]
    Fint4 = integ[11] - x1 * integ[10]
    Fnode4 = @. (V1 * Fint1 + V2 * Fint2 + V3 * Fint3 + V4 * Fint4) * t2N

    Fint1 = x2 * integ[1] - integ[2]
    Fint2 = x2 * integ[4] - integ[5]
    Fint3 = x2 * integ[7] - integ[8]
    Fint4 = x2 * integ[10] - integ[11]
    Fnode3 = @. (V1 * Fint1 + V2 * Fint2 + V3 * Fint3 + V4 * Fint4) * t2N

    magDiffSq = magDiff^2
    magn11mSq = dot(n11m, n11m)
    magn12mSq = dot(n12m, n12m)

    if magDiffSq > eps(Float32) * (magn11mSq + magn12mSq)
        missing, missing, Fnode3Core, Fnode4Core = calcSegSegForce(
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
        Fnode3 = @. Fnode3 + Fnode3Core
        Fnode4 = @. Fnode4 + Fnode4Core

        missing, missing, Fnode3Core, Fnode4Core = calcSegSegForce(
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
        Fnode3 = @. Fnode3 + Fnode3Core
        Fnode4 = @. Fnode4 + Fnode4Core
    end

    # If we flipped the first segment originally, flip the forces round.
    if flip
        Fnode1, Fnode2 = Fnode2, Fnode1
    end

    return Fnode1, Fnode2, Fnode3, Fnode4
end

@inline function ParSegSegInteg(
    aSq_dSq::T1,
    aSq_dSqI::T1,
    x::T1,
    y::T1,
) where {T1 <: Float64}

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
    aSq::T1,
    d::T1,
    c::T1,
    cSq::T1,
    omcSq::T1,
    omcSqI::T1,
    x::T1,
    y::T1,
) where {T1 <: Float64}

    aSq_dSq = aSq + d^2 * omcSqI
    xSq = y^2
    ySq = x^2
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

    den = 1 / sqrt(omcSqI * aSq_dSq)

    integ1 = -2 * den * atan((1 + c) * (Ra + x + y) * den)

    c_1 = aSq_dSq * integ1
    c_5_6 = (c * Ra - c_1) * omcSqI

    integ2 = (c * log_Ra_Rd_t1 - log_Ra_Rd_t2) * omcSqI
    integ3 = (c * log_Ra_Rd_t2 - log_Ra_Rd_t1) * omcSqI
    integ4 = (c * c_1 - Ra) * omcSqI
    integ5 = ylog_Ra_Rd_t2 + c_5_6
    integ6 = xlog_Ra_Rd_t1 + c_5_6

    c_11_12 = integ1 - c * RaInv
    c_15_18 = c * xRaSq_R_t1_I - RaInv
    x_13_14 = x * c_15_18
    c_19 = c * yRaSq_R_t2_I - RaInv
    y_13_14 = y * c_19
    c_16 = log_Ra_Rd_t1 - (x - c * y) * RaInv - cSq * ySqRaSq_R_t2_I
    z_15_18 = y * c_16
    c_17_19 = log_Ra_Rd_t2 - (y - c * x) * RaInv - cSq * xSqRaSq_R_t1_I

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
