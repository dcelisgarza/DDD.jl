function calcSelfForce(
    dlnParams::DislocationP,
    matParams::MaterialP,
    network::DislocationNetwork,
)
    #=
        This is equivalent to doing the dot product over each entry of the vectors.
        c[i,:] =  sum(a .* b, dims = 2) = a[i,:] ⋅ b[i,:]
    =#
    μ = matParams.μ
    ν = matParams.ν
    E = matParams.E
    a = dlnParams.coreRad

    # Un normalised segment vectors.
    bVec = getSegBvec(network)
    tVec = calcSegVector(network)

    # Finding the norm of each line vector.
    L = sqrt.(sum(tVec .* tVec, dims = 2))
    # With the core radius.
    La = sqrt.(sum(len .* len) + a * a)

    # Normalised the dislocation network vector, the sum of all the segment vectors has norm 1.
    tVec = @. segVec / len

    # Screw component, scalar projection of bVec onto t.
    # bScrew[i] = bVec[i,:] ⋅ t[i,:]
    bScrew = sum(bVec .* tVec, dims = 2)
    # Finding the norm of the projection squared.
    bScrewSq = @. bScrew * bScrew

    # Edge component, vector rejection of bVec onto t.
    # bEdgeVec[i, :] = bVec[i,:] - (bVec[i,:] ⋅ t[i,:]) t
    bEdgeVec = @. bVec - bScrew * tVec
    # Finding the norm squared of each edge component.
    # bEdgeSq[i] = bEdgeVec[i,:] ⋅ bEdgeVec[i,:]
    bEdgeSq = sum(bEdgeVec .* bEdgeVec, dims = 2)

    #=
        A. Arsenlis et al, Modelling Simul. Mater. Sci. Eng. 15 (2007)
        553?595: gives this expression in appendix A p590
        f^{s}_{43} = -(μ/(4π)) [ t × (t × b)](t ⋅ b) { v/(1-v) ( ln[
        (L_a + L)/a] - 2*(L_a - a)/L ) - (L_a - a)^2/(2La*L) }

    #   tVec × (tVec × bVec)    = tVec (tVec ⋅ bVec) - bVec (tVec ⋅ tVec)
                                = tVec * bScrew - bVec
                                = - bEdgeVec
    =#
    # Torsional component of the elastic self interaction force. This is the scalar component of the above equation.
    tor = @. (μ / (4π)) *
       bScrew *
       (
           ν / (1 - ν) * (log((La + L) / a) - 2 * (La - a) / a) -
           (La - a) * (La - a) / (2 * La * L)
       )

    # Torsional component of core self interaction.
    torCore = @. 2 * E * ν / (1 - ν) * bScrew
    torTot = @. tor + tCore

    # Longitudinal component of core self interaction.
    lonCore = @. (bScrewSq + bEdgeSq / (1 - ν)) * E

    # Force on node 2 = -Force on node 1.
    f2 = @. torTot * bEdgeVec - lonCore * tVec
    f1 = -f2

    return f1, f2
end
