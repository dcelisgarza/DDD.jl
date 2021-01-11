
function buildMesh(matParams::T1, femParams::T2) where {T1 <: MaterialParameters,T2 <: FEMParameters}
    return buildMesh(femParams.type, matParams, femParams)
end
function buildMesh(::T1, matParams::T2, femParams::T3) where {T1 <: DispatchRegularCuboidMesh,T2 <: MaterialParameters,T3 <: FEMParameters}
    return RegularCuboidMesh(femParams.order, matParams, femParams)
end

"""
E Tarleton edmund.tarleton@materials.ox.ac.uk
3D FEM code using linear 8 node element with 8 integration pts (2x2x2) per element.

   4.-------.3
   |\\       \\
   | \\      |\\
   1.-\\---- .2\\
    \\  \\   \\ \\
     \\ 8.--------.7
      \\ |     \\ |
       \\|      \\|
        5.--------.6

rectangular domain.
(y going in to the screen) note this is rotated about x axis w.r.t. local (s1,s2,s3) system.
 -------------------------------------------------

^z                   (mx,my)
|
|
|------>x-----------------------------------------
"""
function RegularCuboidMesh(order::T1, matParams::T2, femParams::T3) where {T1 <: LinearElement, T2 <: MaterialParameters, T3 <: FEMParameters}

    μ = matParams.μ
    ν = matParams.ν
    νomνInv = matParams.νopνInv # ν/(1+ν)
    dx = femParams.dx
    dy = femParams.dy
    dz = femParams.dz
    mx = femParams.mx
    my = femParams.my
    mz = femParams.mz

    dxType = typeof(dx)
    mxType = typeof(mx)

    # Element width, height, depth
    w = dx / mx
    h = dy / my
    d = dz / mz

    numElem = mx * my * mz
    mx1 = mx + 1
    my1 = my + 1
    mz1 = mz + 1
    numNode = mx1 * my1 * mz1

    coord = zeros(dxType, 3, numNode)           # Node coordinates.
    connectivity = zeros(mxType, 8, numElem)    # Element connectivity.
    localK = zeros(dxType, 24, 24)          # K for an element.
    nx = Array{SMatrix{3,8,dxType}}(undef, 8) # Gauss nodes real space.
    B = zeros(dxType, 6, 24, 8)

    nodeEl = 1:8 # Local node numbers.
    dofLocal = Tuple(Iterators.flatten((3 * (nodeEl .- 1) .+ 1, 3 * (nodeEl .- 1) .+ 2, 3 * (nodeEl .- 1) .+ 3)))
    
    V1 = zeros(dxType, numElem * 24^2)
    V2 = zeros(dxType, numElem * 24^2)
    V3 = zeros(dxType, numElem * 24^2)

    # For a regular cuboid mesh this is predefined.
    vertices = SMatrix{3,8,dxType}(
                                0, 0, 0,
                                dx, 0, 0,
                                0, dy, 0,
                                dx, dy, 0,
                                0, 0, dz,
                                dx, 0, dz,
                                0, dy, dz,
                                dx, dy, dz
                            )
    E = 2 * μ * (1 + ν)                 # Young's modulus
    Lame = E * νomνInv / (1 - 2 * ν)    # Lamé's relation
    CDiag = Lame + 2 * μ
    
    # Stiffness tensor.
    C = SMatrix{6,6,dxType}(
        CDiag,  Lame,   Lame,   0,  0,  0,
        Lame,   CDiag,  Lame,   0,  0,  0,
        Lame,   Lame,   CDiag,  0,  0,  0,
        0,      0,      0,      μ,  0,  0,
        0,      0,      0,      0,  μ,  0,
        0,      0,      0,      0,  0,  μ
    )
    
    # Local nodes Gauss Quadrature. Using 8 nodes per element.
    p = 1 / sqrt(3)
    gaussNodes = SMatrix{3,8,dxType}(
        -p, -p, -p,
         p, -p, -p,
         p,  p, -p,
        -p,  p, -p,
        -p, -p,  p,
         p, -p,  p,
         p,  p,  p,
        -p,  p,  p,
    )

    # Shape functions and their derivatives.
    N = shapeFunction(LinearQuadrangle3D(), gaussNodes[1, :], gaussNodes[2, :], gaussNodes[3, :])
    dNdS = shapeFunctionDeriv(LinearQuadrangle3D(), gaussNodes[1, :], gaussNodes[2, :], gaussNodes[3, :])

    # Fill node coordinates.
    @inbounds for k in 1:mz1
        for j in 1:my1
            @simd for i in 1:mx1
                globalNode = i + (k - 1) * mx1 + (j - 1) * mx1 * mz1
                coord[1, globalNode] = (i - 1) * w
                coord[2, globalNode] = (j - 1) * h
                coord[3, globalNode] = (k - 1) * d
            end
        end
    end

    # Fill element connectivity.
    @inbounds for k in 1:mz
        for j in 1:my
            @simd for i in 1:mx
                globalElem = i + (k - 1) * mx + (j - 1) * mx * mz
                connectivity[1, globalElem] = i + (k - 1) * mx1 + mx1 * mz1 * j
                connectivity[2, globalElem] = connectivity[1, globalElem] + 1
                connectivity[4, globalElem] = i + k * mx1 + mx1 * mz1 * j
                connectivity[3, globalElem] = connectivity[4, globalElem] + 1
                connectivity[5, globalElem] = i + (k - 1) * mx1 + mx1 * mz1 * (j - 1)
                connectivity[6, globalElem] = connectivity[5, globalElem] + 1
                connectivity[8, globalElem] = i + k * mx1 + mx1 * mz1 * (j - 1)
                connectivity[7, globalElem] = connectivity[8, globalElem] + 1
            end
        end
    end

    localCoord = SMatrix{3,8,dxType}(
        coord[1, connectivity[1, 1]], coord[2, connectivity[1, 1]], coord[3, connectivity[1, 1]], 
        coord[1, connectivity[2, 1]], coord[2, connectivity[2, 1]], coord[3, connectivity[2, 1]], 
        coord[1, connectivity[3, 1]], coord[2, connectivity[3, 1]], coord[3, connectivity[3, 1]], 
        coord[1, connectivity[4, 1]], coord[2, connectivity[4, 1]], coord[3, connectivity[4, 1]], 
        coord[1, connectivity[5, 1]], coord[2, connectivity[5, 1]], coord[3, connectivity[5, 1]], 
        coord[1, connectivity[6, 1]], coord[2, connectivity[6, 1]], coord[3, connectivity[6, 1]], 
        coord[1, connectivity[7, 1]], coord[2, connectivity[7, 1]], coord[3, connectivity[7, 1]], 
        coord[1, connectivity[8, 1]], coord[2, connectivity[8, 1]], coord[3, connectivity[8, 1]]            
    )

    localdNdS = SMatrix{8,3,dxType}(
        dNdS[1][1, 1], dNdS[1][1, 2], dNdS[1][1, 3], dNdS[1][1, 4], dNdS[1][1, 5], dNdS[1][1, 6], dNdS[1][1, 7], dNdS[1][1, 8],
        dNdS[1][2, 1], dNdS[1][2, 2], dNdS[1][2, 3], dNdS[1][2, 4], dNdS[1][2, 5], dNdS[1][2, 6], dNdS[1][2, 7], dNdS[1][2, 8],
        dNdS[1][3, 1], dNdS[1][3, 2], dNdS[1][3, 3], dNdS[1][3, 4], dNdS[1][3, 5], dNdS[1][3, 6], dNdS[1][3, 7], dNdS[1][3, 8],
    )

    J = localCoord * localdNdS
    detJ = det(J)
    invJ = inv(J)

    @inbounds @simd for q in nodeEl
        nx[q] = invJ' * dNdS[q]
    end

    # Gauss quadrature nodes.
    @inbounds for q in nodeEl
        # Nodes in element.
        @simd for a in nodeEl
            idx = (a - 1) * 3
            B[1, idx + 1, q] = nx[q][1, a]
            B[1, idx + 2, q] = 0
            B[1, idx + 3, q] = 0

            B[2, idx + 1, q] = 0
            B[2, idx + 2, q] = nx[q][2, a]
            B[2, idx + 3, q] = 0

            B[3, idx + 1, q] = 0
            B[3, idx + 2, q] = 0
            B[3, idx + 3, q] = nx[q][3, a]

            B[4, idx + 1, q] = nx[q][2, a]
            B[4, idx + 2, q] = nx[q][1, a]
            B[4, idx + 3, q] = 0

            B[5, idx + 1, q] = nx[q][3, a]
            B[5, idx + 2, q] = 0
            B[5, idx + 3, q] = nx[q][1, a]

            B[6, idx + 1, q] = 0
            B[6, idx + 2, q] = nx[q][3, a]
            B[6, idx + 3, q] = nx[q][2, a]
        end
    end

    @inbounds @simd for q in nodeEl
        localK += B[:, :, q]' * C * B[:, :, q]
    end

    localK = localK * detJ
    cntr = 0
    @inbounds for p in 1:numElem
        globalNode = connectivity[nodeEl, p]; # Global node numbers
        dofGlobal = Tuple(Iterators.flatten((3 * (globalNode .- 1) .+ 1, 3 * (globalNode .- 1) .+ 2, 3 * (globalNode .- 1) .+ 3))) # Global degree of freedom
        for i = 1:24
            @simd for j = 1:24
                cntr += 1
                V1[cntr] = dofGlobal[i]
                V2[cntr] = dofGlobal[j]
                V3[cntr] = localK[dofLocal[i], dofLocal[j]]
            end
        end
    end

    numNode3 = 3 * numNode
    globalK = sparse(V1, V2, V3, numNode3, numNode3)

    return RegularCuboidMesh(
        order,
        vertices,
        C,
        dx,
        dy,
        dz,
        mx,
        my,
        mz,
        numElem,
        numNode,
        w,
        h,
        d,
        B,
        coord,
        connectivity,
        globalK,
    )
    
end
function RegularCuboidMesh(; order::T1, matParams::T2, femParams::T3) where {T1<:AbstractElementOrder, T2 <: MaterialParameters, T3 <: FEMParameters}
    return RegularCuboidMesh(order, matParams, femParams)
end