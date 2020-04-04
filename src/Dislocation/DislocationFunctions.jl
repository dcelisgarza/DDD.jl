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
```
```
Analytical solution of the force between two dislocation segments. Details are found in Appendix A.1. in ["Enabling Strain Hardening Simulations with Dislocation Dynamics" by A. Arsenlis et al.](https://doi.org/10.1088%2F0965-0393%2F15%2F6%2F001)

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
    SegSegForce = zeros(3, numSegs)

    for i in 1:numSegs
        b1 = bVec[:, i]
        n11 = node1[:, i]
        n12 = node2[:, i]
        t1 = @. n12 - n11
        t1N = 1 / norm(t1)
        t1 = @. t1 * t1N
        for j = (i + 1):numSegs
            b2 = bVec[:, j]
            n21 = node1[:, j]
            n22 = node2[:, j]
            t2 = @. n22 - n21
            t2N = 1 / norm(t2)
            t2 = @. t2 * t2N

            t2ct1 = cross(t2, t1)
            c = dot(t2, t1)
            cSq = c^2
            omcSq = 1 - cSq
            if omcSq > sqrt(eps(Float64))
                SegSegForce[:, j] = calcSegSegForce(
                    a,
                    μ,
                    ν,
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
                    c,
                    cSq,
                    omcSq,
                )
            else
                SegSegForce[:, j] = calcSegSegForce(
                    a,
                    μ,
                    ν,
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
                )
            end
        end
    end
    return transpose(SegSegForce)
end
@inline function calcSegSegForce(
    a::T1,
    μ::T1,
    ν::T1,
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
    c::T1,
    cSq::T1,
    omcSq::T1,
) where {T1 <: Float64, T2 <: AbstractVector{<:Float64}}
    SegSegForce = zeros(3)

    omcSqI = 1 / omcSq

    R1 = @. n21 - n11
    R2 = @. n22 - n12
    d = dot(R2, t2ct1) * omcSqI

    lim11 = dot(R1, t1)
    lim12 = dot(R1, t2)
    lim21 = dot(R2, t1)
    lim22 = dot(R2, t2)

    x1 = (lim12 - c * lim11) * omcSqI
    x2 = (lim22 - c * lim21) * omcSqI

    y1 = (lim11 - c * lim12) * omcSqI
    y2 = (lim21 - c * lim22) * omcSqI

    calcSegSegIntegrals(a, d, c, cSq, omcSq, omcSqI, x1, y1)
    calcSegSegIntegrals(a, d, c, cSq, omcSq, omcSqI, x1, y2)
    calcSegSegIntegrals(a, d, c, cSq, omcSq, omcSqI, x2, y1)
    calcSegSegIntegrals(a, d, c, cSq, omcSq, omcSqI, x2, y2)

    return SegSegForce
end

@inline function calcSegSegForce(
    a::T1,
    μ::T1,
    ν::T1,
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
) where {T1 <: Float64, T2 <: AbstractVector{<:Float64}}
    SegSegForce = zeros(3)
    return SegSegForce
end

@inline function calcSegSegIntegrals(a, d, c, cSq, omcSq, omcSqI, x1, y1) end
