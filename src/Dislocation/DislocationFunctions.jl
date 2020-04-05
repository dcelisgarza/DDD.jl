@inline function calcSelfForce(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
)

    μ = matParams.μ
    ν = matParams.ν
    omNuInv = 1 / (1 - ν)
    E = matParams.E
    a = dlnParams.coreRad

    idx = network.segIdx
    coord = network.coord
    # Un normalised segment vectors.
    bVec = network.bVec[idx[:, 1], :]
    tVec = coord[idx[:, 3], :] - coord[idx[:, 2], :]

    # Finding the norm of each line vector.
    L = dimNorm(tVec; dims = 2)
    Linv = @. 1 / L
    # Finding the non-singular norm.
    La = @. sqrt(L * L + a * a)

    # Normalised the dislocation network vector, the sum of all the segment vectors has norm 1.
    tVec = @. tVec * Linv

    # Screw component, scalar projection of bVec onto t.
    # bScrew[i] = bVec[i,:] ⋅ t[i,:]
    bScrew = dimDot(bVec, tVec; dims = 2)

    # Edge component, vector rejection of bVec onto t.
    # bEdgeVec[i, :] = bVec[i,:] - (bVec[i,:] ⋅ t[i,:]) t
    bEdgeVec = @. bVec - bScrew * tVec
    # Finding the norm squared of each edge component.
    # bEdgeSq[i] = bEdgeVec[i,:] ⋅ bEdgeVec[i,:]
    bEdgeSq = dimDot(bEdgeVec, bEdgeVec; dims = 2)

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
    LaMa = @. La - a
    tor = @. (μ / (4π)) *
       bScrew *
       (
           ν * omNuInv * (log((La + L) / a) - 2 * LaMa * Linv) -
           LaMa^2 / (2 * La * L)
       )

    # Torsional component of core self interaction.
    torCore = @. 2 * E * ν * omNuInv * bScrew
    torTot = @. tor + torCore

    # Longitudinal component of core self interaction.
    lonCore = @. (bScrew^2 + bEdgeSq * omNuInv) * E

    # Force on node 2 = -Force on node 1.
    f2 = @. torTot * bEdgeVec - lonCore * tVec
    f1 = -f2

    return (f1, f2)
end

"""
!!! Note
    This function is based on the SegSegForces function by A. Arsenlis et al. However, has been quite heavily modified to increase its memory and computational efficiency as segment-segment interactions are the bottleneck of 3D dislocation dynamics simulations.
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
)
    μ = matParams.μ
    ν = matParams.ν
    a = dlnParams.coreRad
    # Constants.
    aSq = a^2
    μ4π = μ / (4π)
    μ8π = μ4π / 2
    μ8πaSq = aSq * μ8π
    μ4πν = μ4π / (1 - ν)
    μ4πνaSq = aSq * μ4πν

    # Un normalised segment vectors.
    numSegs = network.numSeg
    idx = network.segIdx
    coord = network.coord
    bVec = transpose(network.bVec[idx[:, 1], :])
    node1 = transpose(coord[idx[:, 2], :])
    node2 = transpose(coord[idx[:, 3], :])

    b1 = zeros(3)
    b2 = zeros(3)
    t1 = zeros(3)
    t2 = zeros(3)
    n11 = zeros(3) # Segment 1, node 1
    n12 = zeros(3) # Segment 1, node 2
    n21 = zeros(3) # Segment 2, node 1
    n22 = zeros(3) # Segment 2, node 2
    t2ct1 = zeros(3)
    t1ct2 = zeros(3)
    b2ct2 = zeros(3)
    b1ct1 = zeros(3)
    t2ct1ct2 = zeros(3)
    t1ct2ct2 = zeros(3)
    t2cb1ct2 = zeros(3)
    t1cb2ct1 = zeros(3)
    b1ct1ct2 = zeros(3)
    b2ct2ct1 = zeros(3)
    integ = zeros(19)
    integt = zeros(19)
    SegSegForcet = zeros(3, 4)
    SegSegForce = zeros(3, 2, numSegs - 1)

    # We're also explcitly using dot products and norms because caching due to variable reuse makes things faster when explicitly typed.
    for i = 1:(numSegs - 1)
        b1 = bVec[:, i]
        n11 = node1[:, i]
        n12 = node2[:, i]
        t1 = @. n12 - n11
        t1N = 1 / sqrt(t1[1] * t1[1] + t1[2] * t1[2] + t1[3] * t1[3])
        t1 = @. t1 * t1N
        for j = (i + 1):numSegs
            b2 = bVec[:, j]
            n21 = node1[:, j]
            n22 = node2[:, j]
            t2 = @. n22 - n21
            t2N = 1 / sqrt(t2[1] * t2[1] + t2[2] * t2[2] + t2[3] * t2[3])
            t2 = @. t2 * t2N

            # t2 × t1
            t2ct1[1] = t2[2] * t1[3] - t2[3] * t1[2]
            t2ct1[2] = t2[3] * t1[1] - t2[1] * t1[3]
            t2ct1[3] = t2[1] * t1[2] - t2[2] * t1[1]
            # b2 × t2
            b2ct2[1] = b2[2] * t2[3] - b2[3] * t2[2]
            b2ct2[2] = b2[3] * t2[1] - b2[1] * t2[3]
            b2ct2[3] = b2[1] * t2[2] - b2[2] * t2[1]
            # b1 × t1
            b1ct1[1] = b1[2] * t1[3] - b1[3] * t1[2]
            b1ct1[2] = b1[3] * t1[1] - b1[1] * t1[3]
            b1ct1[3] = b1[1] * t1[2] - b1[2] * t1[1]

            c = t2[1] * t1[1] + t2[2] * t1[2] + t2[3] * t1[3]
            cSq = c^2
            omcSq = 1 - cSq
            if omcSq > eps(Float64)
                calcSegSegForce(
                    aSq,
                    μ4π,
                    μ8π,
                    μ8πaSq,
                    μ4πν,
                    μ4πνaSq,
                    b1,
                    t1,
                    t1N,
                    n11,
                    n12,
                    b2,
                    t2,
                    t2N,
                    n21,
                    n22,
                    t2ct1,
                    t1ct2,
                    b2ct2,
                    b1ct1,
                    c,
                    cSq,
                    omcSq,
                    integ,
                    integt,
                    SegSegForcet,
                )
            else
                calcSegSegForce(
                    aSq,
                    μ4π,
                    μ8π,
                    μ8πaSq,
                    μ4πν,
                    μ4πνaSq,
                    b1,
                    t1,
                    t1N,
                    n11,
                    n12,
                    b2,
                    t2,
                    t2N,
                    n21,
                    n22,
                    integ,
                    integt,
                    SegSegForcet,
                )
            end
        end
    end
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
    t1::T2,
    t1N::T1,
    n11::T2,
    n12::T2,
    b2::T2,
    t2::T2,
    t2N::T1,
    n21::T2,
    n22::T2,
    t2ct1::T2,
    t1ct2::T2,
    b2ct2::T2,
    b1ct1::T2,
    c::T1,
    cSq::T1,
    omcSq::T1,
    integ::T2,
    integt::T2,
    SegSegForce::T3,
) where {
    T1 <: Float64,
    T2 <: AbstractVector{<:Float64},
    T3 <: AbstractArray{<:Float64, N} where {N},
}

    omcSqI = 1 / omcSq

    R1 = @. n21 - n11
    R2 = @. n22 - n12
    # Caching makes explicit expressions faster than dot products.
    d = (R2[1] * t2ct1[1] + R2[2] * t2ct1[2] + R2[3] * t2ct1[3]) * omcSqI

    μ4πd = μ4π * d
    μ8πd = μ8π * d
    μ4πνd = μ4πν * d
    μ4πνdSq = μ4πνd * d
    μ4πνdCu = μ4πνdSq * d
    μ4πνaSqd = μ4πνaSq * d
    μ8πaSqd = μ8πaSq * d

    lim11 = R1[1] * t1[1] + R1[2] * t1[2] + R1[3] * t1[3]
    lim12 = R1[1] * t2[1] + R1[2] * t2[2] + R1[3] * t2[3]
    lim21 = R2[1] * t1[1] + R2[2] * t1[2] + R2[3] * t1[3]
    lim22 = R2[1] * t2[1] + R2[2] * t2[2] + R2[3] * t2[3]

    t2db2 = t2[1] * b2[1] + t2[2] * b2[2] + t2[3] * b2[3]
    t2db1 = t2[1] * b1[1] + t2[2] * b1[2] + t2[3] * b1[3]
    t1db2 = t1[1] * b2[1] + t1[2] * b2[2] + t1[3] * b2[3]
    t1db1 = t1[1] * b1[1] + t1[2] * b1[2] + t1[3] * b1[3]
    t2ct1db2 = t2ct1[1] * b2[2] + t2ct1[2] * b2[2] + t2ct1[3] * b2[3]
    t1ct2db1 = t1ct2[1] * b1[1] + t1ct2[2] * b1[2] + t1ct2[3] * b1[3]

    x1 = (lim12 - c * lim11) * omcSqI
    x2 = (lim22 - c * lim21) * omcSqI
    y1 = (lim11 - c * lim12) * omcSqI
    y2 = (lim21 - c * lim22) * omcSqI

    integ = SegSegInteg(aSq, d, c, cSq, omcSq, omcSqI, x1, y1, integt)
    integ -= SegSegInteg(aSq, d, c, cSq, omcSq, omcSqI, x1, y2, integt)
    integ -= SegSegInteg(aSq, d, c, cSq, omcSq, omcSqI, x2, y1, integt)
    integ += SegSegInteg(aSq, d, c, cSq, omcSq, omcSqI, x2, y2, integt)

    #=
        μ4πd
        μ8πd
        μ4πνd
        μ4πνdSq
        μ4πνdCu
        μ4πνaSqd
        μ8πaSqd
        aSq
        μ4π
        μ8π
        μ8πaSq
        μ4πν
        μ4πνaSq
    =#

    # t1 × t2
    t1ct2 = -t2ct1
    # b2 × t2
    b2ct2[1] = b2[2] * t2[3] - b2[3] * t2[2]
    b2ct2[2] = b2[3] * t2[1] - b2[1] * t2[3]
    b2ct2[3] = b2[1] * t2[2] - b2[2] * t2[1]
    # b1 × t1
    b1ct1[1] = b1[2] * t1[3] - b1[3] * t1[2]
    b1ct1[2] = b1[3] * t1[1] - b1[1] * t1[3]
    b1ct1[3] = b1[1] * t1[2] - b1[2] * t1[1]
    # t2 ⋅ b2
    t2db2 = t2[1] * b2[1] + t2[2] * b2[2] + t2[3] * b2[3]
    # t2 ⋅ b1
    t2db1 = t2[1] * b1[1] + t2[2] * b1[2] + t2[3] * b1[3]
    # t1 ⋅ b2
    t1db2 = t1[1] * b2[1] + t1[2] * b2[2] + t1[3] * b2[3]
    # t1 ⋅ b1
    t1db1 = t1[1] * b1[1] + t1[2] * b1[2] + t1[3] * b1[3]
    # (t2 × t1) ⋅ b2
    t2ct1db2 = t2ct1[1] * b2[1] + t2ct1[2] * b2[2] + t2ct1[3] * b2[3]
    # (t1 × t2) ⋅ b1
    t1ct2db1 = t1ct2[1] * b1[1] + t1ct2[2] * b1[2] + t1ct2[3] * b1[3]
    # (b1 × t1) ⋅ b2
    b1ct1db2 = b1ct1[1] * b2[1] + b1ct1[2] * b2[2] + b1ct1[3] * b2[3]
    # (b2 × t2) ⋅ b1
    b2ct2db1 = b2ct2[1] * b1[1] + b2ct2[2] * b1[2] + b2ct2[3] * b1[3]

    t2ct1ct2 = @. t1 - c * t2
    t1ct2ct1 = @. t2 - c * t1
    t2cb1ct2 = @. b1 - t2db1 * t2
    t1cb2ct1 = @. b2 - t1db2 * t1
    b1ct1ct2 = @. t2db1 * t1 - c * b1
    b2ct2ct1 = @. t1db2 * t2 - c * b2
    t2ct1cb1dt1 = t2db1 - t1db1 * c
    t1ct2cb2dt2 = t1db2 - t2db2 * c
    t2ct1cb1db2 = t2db1 * t1db2 - t1db1 * t2db2

    # Seg 2 (nodes 4-3)
    tmp1 = t2db1 * t1db2 + t2ct1cb1db2
    I00a = @. tmp1 * t1ct2
    I00b = @. b2ct2 * t2ct1cb1dt1
    I01a = @. t2ct1 * b1ct1db2 - b1ct1ct2 * t1db2
    I10a = @. t2cb1ct2 * t1db2 - t2ct1 * b2ct2db1
    I10b = @. b2ct2 * t1ct2db1

    tmp1 = μ4πνd * t2ct1db2
    tmp2 = μ4πνd * b1ct1db2
    V[:, 1] = @. μ4πd * I00a - μ4πνd * I00b + tmp1 * b1ct1ct2 + tmp2 * t2ct1ct2


    tmp1 = μ4πν * t2db2
    V[:, 2] = @. tmp1 * b1ct1ct2 + μ4π * I10a - μ4πν * I10b

    tmp1 = μ4πν * t1db2
    tmp2 = μ4πν * b1ct1db2
    V[:, 3] = @. μ4π * I01a + tmp1 * b1ct1ct2 - tmp2 * t2ct1

    tmp1 = μ4πνdCu * t2ct1cb1dt1 * t2ct1db2
    V[:, 4] = @. μ8πaSqd * I00a - μ4πνaSqd * I00b - tmp1 * t2ct1ct2

    tmp1 = μ4πνdSq * (t1ct2db1 * t2ct1db2 + t2ct1cb1dt1 * t2db2)
    V[:, 5] = @. μ8πaSq * I10a - a2m4pn * I10b - tmp1 * t2ct1ct2

    tmp1 = μ4πνdSq * t2ct1cb1dt1 * t1db2
    tmp2 = μ4πνdSq * t2ct1cb1dt1 * t2ct1db2
    V[:, 6] = @. μ8πaSq * I01a - tmp1 * t2ct1ct2 + tmp2 * t2ct1

    tmp1 = μ4πνd * (t2ct1cb1dt1 * t2db2 + t1ct2db1 * t2ct1db2)
    tmp2 = μ4πνd * t1ct2db1 * t1db2
    V[:, 7] = @. tmp1 * t2ct1 - tmp2 * t2ct1ct2

    tmp1 = μ4πνd * t1ct2db1 * t2db2
    V[:, 8] = @. -tmp1 * t2ct1ct2

    tmp1 = μ4πνd * t2ct1cb1dt1 * t1db2
    V[:, 9] = @. tmp1 * t2ct1

    tmp1 = μ4πν * t1ct2db1 * t2db2
    V[:, 10] = @. tmp1 * t2ct1

    tmp1 = (μ4πν * t1ct2db1 * t1db2)
    V[:, 11] = @. tmp1 * t2ct1

    Fint[1] = integ[2] - x1 * integ[1]
    Fint[2] = integ[5] - x1 * integ[2]
    Fint[3] = integ[4] - x1 * integ[3]
    Fint[4] = integ[8] - x1 * integ[7]
    Fint[5] = integ[11] - x1 * integ[8]
    Fint[6] = integ[10] - x1 * integ[9]
    Fint[7] = integ[13] - x1 * integ[10]
    Fint[8] = integ[16] - x1 * integ[11]
    Fint[9] = integ[14] - x1 * integ[12]
    Fint[10] = integ[18] - x1 * integ[13]
    Fint[11] = integ[15] - x1 * integ[14]

    fp4 = sum(V.*Fint',dims=2) .* t2N


    Fint[1] = x2 * integ[1] - integ[2]
    Fint[2] = x2 * integ[2] - integ[5]
    Fint[3] = x2 * integ[3] - integ[4]
    Fint[4] = x2 * integ[7] - integ[8]
    Fint[5] = x2 * integ[8] - integ[11]
    Fint[6] = x2 * integ[9] - integ[10]
    Fint[7] = x2 * integ[10] - integ[13]
    Fint[8] = x2 * integ[11] - integ[16]
    Fint[9] = x2 * integ[12] - integ[14]
    Fint[10] = x2 * integ[13] - integ[18]
    Fint[11] = x2 * integ[14] - integ[15]

    fp3x =
        (
            V[:, 1] * Fint[1] +
            V[:, 2] * Fint[2] +
            V[:, 3] * Fint[3] +
            V[:, 4] * Fint[4] +
            V[:, 5] * Fint[5] +
            V[:, 6] * Fint[6] +
            V[:, 7] * Fint[7] +
            V[:, 8] * Fint[8] +
            V[:, 9] * Fint[9] +
            V[:, 10] * Fint[10] +
            V[:, 11] * Fint[11]
        ) * t2N

    fp3y =
        (
            I_003y * Fint[1] +
            I_103y * Fint[2] +
            I_013y * Fint[3] +
            I_005y * Fint[4] +
            I_105y * Fint[5] +
            I_015y * Fint[6] +
            I_115y * Fint[7] +
            I_205y * Fint[8] +
            I_025y * Fint[9] +
            I_215y * Fint[10] +
            I_125y * Fint[11]
        ) * t2N

    fp3z =
        (
            I_003z * Fint[1] +
            I_103z * Fint[2] +
            I_013z * Fint[3] +
            I_005z * Fint[4] +
            I_105z * Fint[5] +
            I_015z * Fint[6] +
            I_115z * Fint[7] +
            I_205z * Fint[8] +
            I_025z * Fint[9] +
            I_215z * Fint[10] +
            I_125z * Fint[11]
        ) * t2N

    # seg 1, nodes 2-1
    tmp1 = t1db2 * t2db1 + t2ct1cb1db2

    I00ax = tmp1 * t2ct1
    I00ay = tmp1 * tctpy
    I00az = tmp1 * tctpz

    I00bx = b1ct1 * t1ct2cb2dt2
    I00by = bpctpy * t1ct2cb2dt2
    I00bz = bpctpz * t1ct2cb2dt2

    tmp1 = μ4πνd * t1ct2db1
    tmp2 = μ4πνd * b2ct2db1

    V[:, 1] = μ4πd * I00ax - μ4πνd * I00bx + tmp1 * b2ct2ct1 + tmp2 * t1ct2ct1
    I_003y = μ4πd * I00ay - μ4πνd * I00by + tmp1 * bctctpy + tmp2 * tpctctpy
    I_003z = μ4πd * I00az - μ4πνd * I00bz + tmp1 * bctctpz + tmp2 * tpctctpz

    tmp1 = μ4πνdCu * t1ct2cb2dt2 * t1ct2db1

    V[:, 4] = μ8πaSqd * I00ax - μ4πνaSqd * I00bx - tmp1 * t1ct2ct1
    I_005y = μ8πaSqd * I00ay - μ4πνaSqd * I00by - tmp1 * tpctctpy
    I_005z = μ8πaSqd * I00az - μ4πνaSqd * I00bz - tmp1 * tpctctpz

    I01ax = t1ct2 * b1ct1db2 - t1cb2ct1 * t2db1
    I01ay = tpcty * b1ct1db2 - tpcbctpy * t2db1
    I01az = tpctz * b1ct1db2 - tpcbctpz * t2db1

    I01bx = -b1ct1 * t2ct1db2
    I01by = -bpctpy * t2ct1db2
    I01bz = -bpctpz * t2ct1db2

    tmp1 = μ4πν * t1db1

    V[:, 3] = -tmp1 * b2ct2ct1 + μ4π * I01ax - μ4πν * I01bx
    I_013y = -tmp1 * bctctpy + μ4π * I01ay - μ4πν * I01by
    I_013z = -tmp1 * bctctpz + μ4π * I01az - μ4πν * I01bz

    tmp1 = μ4πνdSq * (t2ct1db2 * t1ct2db1 + t1ct2cb2dt2 * t1db1)

    V[:, 6] = μ8πaSq * I01ax - a2m4pn * I01bx + tmp1 * t1ct2ct1
    I_015y = μ8πaSq * I01ay - a2m4pn * I01by + tmp1 * tpctctpy
    I_015z = μ8πaSq * I01az - a2m4pn * I01bz + tmp1 * tpctctpz

    I10ax = b2ct2ct1 * t2db1 - t1ct2 * b2ct2db1
    I10ay = bctctpy * t2db1 - tpcty * b2ct2db1
    I10az = bctctpz * t2db1 - tpctz * b2ct2db1

    tmp1 = μ4πν * t2db1
    tmp2 = μ4πν * b2ct2db1

    V[:, 2] = μ4π * I10ax - tmp1 * b2ct2ct1 + tmp2 * t1ct2
    I_103y = μ4π * I10ay - tmp1 * bctctpy + tmp2 * tpcty
    I_103z = μ4π * I10az - tmp1 * bctctpz + tmp2 * tpctz

    tmp1 = μ4πνdSq * t1ct2cb2dt2 * t2db1
    tmp2 = μ4πνdSq * t1ct2cb2dt2 * t1ct2db1

    V[:, 5] = μ8πaSq * I10ax + tmp1 * t1ct2ct1 - tmp2 * t1ct2
    I_105y = μ8πaSq * I10ay + tmp1 * tpctctpy - tmp2 * tpcty
    I_105z = μ8πaSq * I10az + tmp1 * tpctctpz - tmp2 * tpctz

    tmp1 = (μ4πνd * t2ct1db2 * t1db1)

    V[:, 9] = -tmp1 * t1ct2ct1
    I_025y = -tmp1 * tpctctpy
    I_025z = -tmp1 * tpctctpz

    tmp1 = (μ4πνd * t1ct2cb2dt2 * t2db1)

    V[:, 8] = tmp1 * t1ct2
    I_205y = tmp1 * tpcty
    I_205z = tmp1 * tpctz

    tmp1 = μ4πνd * (t1ct2cb2dt2 * t1db1 + t2ct1db2 * t1ct2db1)
    tmp2 = μ4πνd * t2ct1db2 * t2db1

    V[:, 7] = tmp1 * t1ct2 - tmp2 * t1ct2ct1
    I_115y = tmp1 * tpcty - tmp2 * tpctctpy
    I_115z = tmp1 * tpctz - tmp2 * tpctctpz

    tmp1 = (μ4πν * t2ct1db2 * t1db1)

    V[:, 11] = -tmp1 * t1ct2
    I_125y = -tmp1 * tpcty
    I_125z = -tmp1 * tpctz

    tmp1 = (μ4πν * t2ct1db2 * t2db1)

    V[:, 10] = -tmp1 * t1ct2
    I_215y = -tmp1 * tpcty
    I_215z = -tmp1 * tpctz

    Fint[1] = integ[3] - y2 * integ[1]
    Fint[2] = integ[4] - y2 * integ[2]
    Fint[3] = integ[6] - y2 * integ[3]
    Fint[4] = integ[9] - y2 * integ[7]
    Fint[5] = integ[10] - y2 * integ[8]
    Fint[6] = integ[12] - y2 * integ[9]
    Fint[7] = integ[14] - y2 * integ[10]
    Fint[8] = integ[13] - y2 * integ[11]
    Fint[9] = integ[17] - y2 * integ[12]
    Fint[10] = integ[15] - y2 * integ[13]
    Fint[11] = integ[19] - y2 * integ[14]

    fp1x =
        (
            V[:, 1] * Fint[1] +
            V[:, 2] * Fint[2] +
            V[:, 3] * Fint[3] +
            V[:, 4] * Fint[4] +
            V[:, 5] * Fint[5] +
            V[:, 6] * Fint[6] +
            V[:, 7] * Fint[7] +
            V[:, 8] * Fint[8] +
            V[:, 9] * Fint[9] +
            V[:, 10] * Fint[10] +
            V[:, 11] * Fint[11]
        ) * t1N

    fp1y =
        (
            I_003y * Fint[1] +
            I_103y * Fint[2] +
            I_013y * Fint[3] +
            I_005y * Fint[4] +
            I_105y * Fint[5] +
            I_015y * Fint[6] +
            I_115y * Fint[7] +
            I_205y * Fint[8] +
            I_025y * Fint[9] +
            I_215y * Fint[10] +
            I_125y * Fint[11]
        ) * t1N

    fp1z =
        (
            I_003z * Fint[1] +
            I_103z * Fint[2] +
            I_013z * Fint[3] +
            I_005z * Fint[4] +
            I_105z * Fint[5] +
            I_015z * Fint[6] +
            I_115z * Fint[7] +
            I_205z * Fint[8] +
            I_025z * Fint[9] +
            I_215z * Fint[10] +
            I_125z * Fint[11]
        ) * t1N

    Fint[1] = y1 * integ[1] - integ[3]
    Fint[2] = y1 * integ[2] - integ[4]
    Fint[3] = y1 * integ[3] - integ[6]
    Fint[4] = y1 * integ[7] - integ[9]
    Fint[5] = y1 * integ[8] - integ[10]
    Fint[6] = y1 * integ[9] - integ[12]
    Fint[7] = y1 * integ[10] - integ[14]
    Fint[8] = y1 * integ[11] - integ[13]
    Fint[9] = y1 * integ[12] - integ[17]
    Fint[10] = y1 * integ[13] - integ[15]
    Fint[11] = y1 * integ[14] - integ[19]

    fp2x =
        (
            V[:, 1] * Fint[1] +
            V[:, 2] * Fint[2] +
            V[:, 3] * Fint[3] +
            V[:, 4] * Fint[4] +
            V[:, 5] * Fint[5] +
            V[:, 6] * Fint[6] +
            V[:, 7] * Fint[7] +
            V[:, 8] * Fint[8] +
            V[:, 9] * Fint[9] +
            V[:, 10] * Fint[10] +
            V[:, 11] * Fint[11]
        ) * t1N

    fp2y =
        (
            I_003y * Fint[1] +
            I_103y * Fint[2] +
            I_013y * Fint[3] +
            I_005y * Fint[4] +
            I_105y * Fint[5] +
            I_015y * Fint[6] +
            I_115y * Fint[7] +
            I_205y * Fint[8] +
            I_025y * Fint[9] +
            I_215y * Fint[10] +
            I_125y * Fint[11]
        ) * t1N

    fp2z =
        (
            I_003z * Fint[1] +
            I_103z * Fint[2] +
            I_013z * Fint[3] +
            I_005z * Fint[4] +
            I_105z * Fint[5] +
            I_015z * Fint[6] +
            I_115z * Fint[7] +
            I_205z * Fint[8] +
            I_025z * Fint[9] +
            I_215z * Fint[10] +
            I_125z * Fint[11]
        ) * t1N

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
    t1::T2,
    t1N::T1,
    n11::T2,
    n12::T2,
    b2::T2,
    t2::T2,
    t2N::T1,
    n21::T2,
    n22::T2,
    integ::T2,
    integt::T2,
    SegSegForce::T3,
) where {
    T1 <: Float64,
    T2 <: AbstractVector{<:Float64},
    T3 <: AbstractArray{<:Float64, N} where {N},
}
    return SegSegForce
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
    integ::T2,
) where {T1 <: Float64, T2 <: AbstractVector{<:Float64}}

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

    integ[1] = -2 * den * atan((1 + c) * (Ra + x + y) * den)

    c_1 = aSq_dSq * integ[1]
    c_5_6 = (c * Ra - c_1) * omcSqI

    integ[2] = (c * log_Ra_Rd_t1 - log_Ra_Rd_t2) * omcSqI
    integ[3] = (c * log_Ra_Rd_t2 - log_Ra_Rd_t1) * omcSqI
    integ[4] = (c * c_1 - Ra) * omcSqI
    integ[5] = ylog_Ra_Rd_t2 + c_5_6
    integ[6] = xlog_Ra_Rd_t1 + c_5_6

    c_11_12 = integ[1] - c * RaInv
    c_15_18 = c * xRaSq_R_t1_I - RaInv
    x_13_14 = x * c_15_18
    c_19 = c * yRaSq_R_t2_I - RaInv
    y_13_14 = y * c_19
    c_16 = log_Ra_Rd_t1 - (x - c * y) * RaInv - cSq * ySqRaSq_R_t2_I
    z_15_18 = y * c_16
    c_17_19 = log_Ra_Rd_t2 - (y - c * x) * RaInv - cSq * xSqRaSq_R_t1_I

    c15_18_19 = 2 * integ[4]

    integ[7] = (integ[1] - xRaSq_R_t1_I - yRaSq_R_t2_I) / (aSq_dSq)
    integ[8] = (RaSq_R_t1_I - c * RaSq_R_t2_I) * omcSqI
    integ[9] = (RaSq_R_t2_I - c * RaSq_R_t1_I) * omcSqI
    integ[10] = (RaInv - c * (xRaSq_R_t1_I + yRaSq_R_t2_I + integ[1])) * omcSqI
    integ[11] = (xRaSq_R_t1_I + cSq * yRaSq_R_t2_I + c_11_12) * omcSqI
    integ[12] = (yRaSq_R_t2_I + cSq * xRaSq_R_t1_I + c_11_12) * omcSqI
    integ[13] = (integ[3] - x_13_14 + c * (y_13_14 - integ[2])) * omcSqI
    integ[14] = (integ[2] - y_13_14 + c * (x_13_14 - integ[3])) * omcSqI
    integ[15] = (integ[5] - z_15_18 + c * (xSq * c_15_18 - c15_18_19)) * omcSqI
    integ[16] = (xSqRaSq_R_t1_I + c * c_16 + 2 * integ[2]) * omcSqI
    integ[17] = (ySqRaSq_R_t2_I + c * c_17_19 + 2 * integ[3]) * omcSqI
    integ[18] = (c15_18_19 - xSq * c_15_18 + c * (z_15_18 - integ[5])) * omcSqI
    integ[19] = (c15_18_19 - ySq * c_19 + c * (x * c_17_19 - integ[6])) * omcSqI

    return integ
end
