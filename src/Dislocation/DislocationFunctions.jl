function calcSelfForce(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
)
    mu = matParams.μ
    nu = matParams.ν
    E = matParams.E
    a = dlnParams.coreRad
    omNuInv = 1 / (1 - nu)

    # Un normalised segment vectors.
    bVec = getSegBvec(network)
    tVec = calcSegVector(network)

    # Finding the norm of the network's total segment length.
    len = sqrt.(sum(tVec .* tVec, dims = 2))
    lenInv = 1 ./ len

    # With the core radius.
    lenA = sqrt.(sum(len .* len) + a * a)
    lenAInv = 1 ./ lenA

    # Normalised the dislocation network vector, the sum of all the segment vectors has norm 1.
    t = segVec .* lenInv

    # Screw component, project bVec onto t.
    sB = sum(bVec .* t, dims = 2)
    sBSq = sB .* sB

    # Edge component, rejection of bVec onto t.
    eB = bVec .- sB .* t
    eBSq = sum(eB, eB, dims = 2)

    # A. Arsenlis et al, Modelling Simul. Mater. Sci. Eng. 15 (2007)
    # 553?595: gives this expression in appendix A p590
    # f^{s}_{43} = -(mu/(4π)) [ t × (t × b)](t ⋅ b) { v/(1-v) ( ln[
    # (L_a + L)/a] - 2*(L_a - a)/L ) - (L_a - a)^2/(2La*L) }

end
