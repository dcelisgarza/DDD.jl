"""
```
FEMParameters(; 
    type::AbstractMesh,
    order::AbstractElementOrder,
    model::AbstractModel,
    dx, dy, dz, mx, my, mz
)
```
Creates [`FEMParameters`](@ref).
"""
function FEMParameters(; type::AbstractMesh, order::AbstractElementOrder, model::AbstractModel, dx, dy, dz, mx, my, mz)
    return FEMParameters(type, order, model, dx, dy, dz, mx, my, mz)
end

"""
```
ForceDisplacement(; uTilde, uHat, u, fTilde, fHat, f)
```
Creates [`ForceDisplacement`](@ref).
"""
function ForceDisplacement(; uTilde, uHat, u, fTilde, fHat, f)
    return ForceDisplacement(uTilde, uHat, u, fTilde, fHat, f)
end

"""
```
buildMesh(
    matParams::MaterialParameters, 
    femParams::FEMParameters{F1,F2,F3,F4,F5} where {F1<:DispatchRegularCuboidMesh,F2,F3,F4,F5}
)
```
Creates a [`RegularCuboidMesh`](@ref).
"""
function buildMesh(matParams::MaterialParameters, femParams::FEMParameters{F1,F2,F3,F4,F5} where {F1<:DispatchRegularCuboidMesh,F2,F3,F4,F5})
    return RegularCuboidMesh(matParams, femParams)
end

"""
```
RegularCuboidMesh(
    matParams::MaterialParameters,
    femParams::FEMParameters{F1,F2,F3,F4,F5} where {F1<:DispatchRegularCuboidMesh,F2<:LinearElement,F3,F4,F5}
)
```
Created by: E. Tarleton `edmund.tarleton@materials.ox.ac.uk`

3D FEM code using linear 8 node element with 8 integration pts (2x2x2) per element.

```
   4.-------.3
   | \\       |\\
   |  \\      | \\    my
   1.--\\---- .2 \\
    \\   \\     \\  \\
     \\  8.--------.7
      \\  |      \\ |  mz
       \\ |       \\|
         5.--------.6
             mx
y   ^z
 ↖  |
  \\ |
   \\|---->x
```
!!! note
    This is rotated about the `x` axis w.r.t. to the local `(s1, s2, s3)` system. [`calc_σHat`](@ref) uses custom linear shape functions that eliminate the need for a Jacobian, speeding up the calculation. We keep the Jacobian in this function because it's only run at simulation initialisation so it's not performance critical, which lets us use standard shape functions.
"""
function RegularCuboidMesh(
    matParams::MaterialParameters,
    femParams::FEMParameters{F1,F2,F3,F4,F5} where {F1<:DispatchRegularCuboidMesh,F2<:LinearElement,F3,F4,F5}
)
    μ = matParams.μ
    ν = matParams.ν
    νomνInv = matParams.νopνInv # ν/(1+ν)
    order = femParams.order
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
    mxm1 = mx - 1
    mym1 = my - 1
    mzm1 = mz - 1
    mx1mz1 = mx1 * mz1
    mymx1mz1 = my * mx1mz1
    mzmx1 = mz * mx1
    numNode = mx1 * my1 * mz1

    # Length scale.
    scale = SVector{3,dxType}(dx, dy, dz) * cbrt(dx * dy * dz)

    # For a regular cuboid mesh this is predefined. Making a volumetric polytope allows us to check if points lie inside the volume simply by doing [x, y, z] ∈ vertices, or checking if they do not belong by doing [x, y, z] ∉ vertices. The symbols are typed as \in and \notin + Tab.
    vtx = SMatrix{3,8,dxType}(
                0, 0, 0,
                dx, 0, 0,
                0, dy, 0,
                dx, dy, 0,
                0, 0, dz,
                dx, 0, dz,
                0, dy, dz,
                dx, dy, dz
            )

    vertices = VPolytope(vtx)
    
    # Faces as defined by the vertices.
    faces = SMatrix{4,6,mxType}(
        2, 1, 4, 3, # xy plane @ min z
        5, 6, 7, 8, # xy plane @ max z
        1, 2, 5, 6, # xz plane @ min y
        4, 3, 8, 7, # xz plane @ max y
        3, 1, 7, 5, # yz plane @ min x
        2, 4, 6, 8, # yz plane @ max x
    )

    faceMidPt = mean(vtx[:, faces], dims = 2)[:,1,:]
    
    # Face normal of the corresponding face.
    faceNorm = SMatrix{3,6,dxType}(
         0,  0, -1,
         0,  0,  1,
         0, -1,  0,
         0,  1,  0,
        -1,  0,  0,
         1,  0,  0,
    )

    coord = zeros(dxType, 3, numNode)           # Node coordinates.
    connectivity = zeros(mxType, 8, numElem)    # Element connectivity.
    
    # Nodes corresponding to the vertices.
    cornerNode = (
                x0y0z0 = 1,
                x1y0z0 = mx1,
                x0y1z0 = 1 + mymx1mz1,
                x1y1z0 = mx1 + mymx1mz1,
                x0y0z1 = 1 + mzmx1,
                x1y0z1 = mx1mz1,
                x0y1z1 = 1 + mymx1mz1 + mzmx1,
                x1y1z1 = mx1 + mymx1mz1 + mzmx1
            )

    # Nodes corresponding to the edges.
    edgeNode = (x_y0z0 = zeros(mxType, mxm1), x_y1z0 = zeros(mxType, mxm1), x_y0z1 = zeros(mxType, mxm1), x_y1z1 = zeros(mxType, mxm1), 
                y_x0z0 = zeros(mxType, mym1), y_x1z0 = zeros(mxType, mym1), y_x0z1 = zeros(mxType, mym1), y_x1z1 = zeros(mxType, mym1),
                z_x0y0 = zeros(mxType, mzm1), z_x1y0 = zeros(mxType, mzm1), z_x0y1 = zeros(mxType, mzm1), z_x1y1 = zeros(mxType, mzm1))
    # Nodes corresponding to the faces.
    faceNode = (xy_z0 = zeros(mxType, mxm1 * mym1), xy_z1 = zeros(mxType, mxm1 * mym1),
                xz_y0 = zeros(mxType, mxm1 * mzm1), xz_y1 = zeros(mxType, mxm1 * mzm1),
                yz_x0 = zeros(mxType, mym1 * mzm1), yz_x1 = zeros(mxType, mym1 * mzm1))

    localK = zeros(dxType, 24, 24)  # K for an element.
    B = zeros(dxType, 6, 24, 8)     # Jacobian matrix.

    nodeEl = 1:8 # Local node numbers.
    dofLocal = Tuple(Iterators.flatten((3 * (nodeEl .- 1) .+ 1, 3 * (nodeEl .- 1) .+ 2, 3 * (nodeEl .- 1) .+ 3)))
    
    # TODO
    V1 = zeros(dxType, numElem * 24^2)
    V2 = zeros(dxType, numElem * 24^2)
    V3 = zeros(dxType, numElem * 24^2)

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
    @inbounds @simd for k in 1:mz1
        km1 = k - 1
        for j in 1:my1
            jm1 = j - 1
            for i in 1:mx1
                globalNode = i + km1 * mx1 + jm1 * mx1mz1
                coord[1, globalNode] = (i - 1) * w
                coord[2, globalNode] = jm1 * h
                coord[3, globalNode] = km1 * d
            end
        end
    end

    # Fill element connectivity.
    @inbounds @simd for k in 1:mz
        km1 = k - 1
        kmx1 = k * mx1
        km1mx1 = km1 * mx1
        for j in 1:my
            jm1 = j - 1
            mx1mz1j = mx1mz1 * j
            mx1mz1jm1 = mx1mz1 * jm1
            for i in 1:mx
                globalElem = i + km1 * mx + jm1 * mx * mz
                connectivity[1, globalElem] = i + km1mx1 + mx1mz1j
                connectivity[2, globalElem] = connectivity[1, globalElem] + 1
                connectivity[4, globalElem] = i + kmx1 + mx1mz1j
                connectivity[3, globalElem] = connectivity[4, globalElem] + 1
                connectivity[5, globalElem] = i + km1mx1 + mx1mz1jm1
                connectivity[6, globalElem] = connectivity[5, globalElem] + 1
                connectivity[8, globalElem] = i + kmx1 + mx1mz1jm1
                connectivity[7, globalElem] = connectivity[8, globalElem] + 1
            end
        end
    end

    # Edge node along x
    @inbounds @simd for i in 1:mxm1
        ip1 = i + 1
        edgeNode[:x_y0z0][i] = ip1            # x_y0z0
        edgeNode[:x_y1z0][i] = ip1 + mymx1mz1 # x_y1z0
        edgeNode[:x_y0z1][i] = ip1 + mzmx1    # x_y0z1
        edgeNode[:x_y1z1][i] = ip1 + mzmx1 + mymx1mz1 # x_y1z1
    end
    # Edge node along y
    @inbounds @simd for i in 1:mym1
        imx1mz1 = i * mx1mz1
        edgeNode[:y_x0z0][i] = 1 + imx1mz1    # y_x0z0
        edgeNode[:y_x1z0][i] = mx1 + imx1mz1  # y_x1z0
        edgeNode[:y_x0z1][i] = 1 + imx1mz1 + mzmx1    # y_x0z1
        edgeNode[:y_x1z1][i] = mx1 + imx1mz1 + mzmx1  # y_x1z1
    end
    # Edge node along z
    @inbounds @simd for i in 1:mzm1
        ip1 = i + 1
        ip1mx1 = ip1 * mx1
        imx1 = i * mx1
        edgeNode[:z_x0y0][i] = 1 + imx1 # z_x0y0
        edgeNode[:z_x1y0][i] = ip1mx1   # z_x1y0
        edgeNode[:z_x0y1][i] = 1 + mymx1mz1 + imx1   # z_x0y1
        edgeNode[:z_x1y1][i] = mymx1mz1 + ip1mx1 # z_x1y1
    end
    # Face node xy
    @inbounds @simd for j in 1:mym1
        jm1 = j - 1
        jmx1mz1 = j * mx1mz1
        mxm1jm1 = mxm1 * jm1
        jmx1mz1p1 = jmx1mz1 + 1
        for i in 1:mxm1
            ipmxm1jm1 = i + mxm1jm1
            faceNode[:xy_z0][ipmxm1jm1] = jmx1mz1p1 + i  # xy_z0
            faceNode[:xy_z1][ipmxm1jm1] = jmx1mz1p1 + i + mzmx1  # xy_z1
        end
    end
    # Face node xz
    @inbounds @simd for j in 1:mzm1
        jm1 = j - 1
        jmx1 = j * mx1
        mxm1jm1 = mxm1 * jm1
        jmx1p1 = 1 + jmx1
        for i in 1:mxm1
            ipmxm1jm1 = i + mxm1jm1
            faceNode[:xz_y0][ipmxm1jm1] = jmx1p1 + i # xz_y0
            faceNode[:xz_y1][ipmxm1jm1] = jmx1p1 + i + mymx1mz1 # xz_y1
        end
    end
    # Face node yz
    @inbounds @simd for j in 1:mym1
        jm1 = j - 1
        jmx1mz1 = j * mx1mz1
        for i in 1:mzm1
            imx1 = i * mx1
            jmx1mz1pimx1 = jmx1mz1 + imx1
            faceNode[:yz_x0][i + jm1*mzm1] = 1 + jmx1mz1pimx1    # yz_x0
            faceNode[:yz_x1][i + jm1*mzm1] = mx1 + jmx1mz1pimx1  # yz_x1
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
    invJ = inv(J)'

    # Gauss quadrature nodes.
    @inbounds @simd for q in nodeEl
        nx = invJ * dNdS[q]
        # Nodes in element.
        for a in nodeEl
            idx = (a - 1) * 3
            idx1 = idx + 1
            idx2 = idx + 2
            idx3 = idx + 3
            B[1, idx1, q] = nx[1, a]
            B[1, idx2, q] = 0
            B[1, idx3, q] = 0

            B[2, idx1, q] = 0
            B[2, idx2, q] = nx[2, a]
            B[2, idx3, q] = 0

            B[3, idx1, q] = 0
            B[3, idx2, q] = 0
            B[3, idx3, q] = nx[3, a]

            B[4, idx1, q] = nx[2, a]
            B[4, idx2, q] = nx[1, a]
            B[4, idx3, q] = 0

            B[5, idx1, q] = nx[3, a]
            B[5, idx2, q] = 0
            B[5, idx3, q] = nx[1, a]

            B[6, idx1, q] = 0
            B[6, idx2, q] = nx[3, a]
            B[6, idx3, q] = nx[2, a]
    end
    end

    @inbounds @simd for q in nodeEl
        localK += B[:, :, q]' * C * B[:, :, q]
    end

    localK = localK * detJ
    cntr = 0
    @inbounds @simd for p in 1:numElem
        globalNode = connectivity[nodeEl, p]; # Global node numbers
        dofGlobal = Tuple(Iterators.flatten((3 * (globalNode .- 1) .+ 1, 3 * (globalNode .- 1) .+ 2, 3 * (globalNode .- 1) .+ 3))) # Global degree of freedom
        for i = 1:24
            for j = 1:24
                cntr += 1
                V1[cntr] = dofGlobal[i]
                V2[cntr] = dofGlobal[j]
                V3[cntr] = localK[dofLocal[i], dofLocal[j]]
            end
    end
    end

    numNode3 = 3 * numNode
    globalK = sparse(V1, V2, V3, numNode3, numNode3)
    droptol!(globalK, eps(dxType))


    return RegularCuboidMesh(
            order,
            dx,
            dy,
            dz,
            mx,
            my,
            mz,
            w,
            h,
            d,
            scale,
            numElem,
            numNode,
            C,
            vertices,
            faces,
            faceNorm,
            faceMidPt,
            cornerNode,
            edgeNode,
            faceNode,
            coord,
            connectivity,
            globalK,
        )
end

"""
```
BoundaryNode(; type, index, node)
```
Create [`BoundaryNode`](@ref).
"""
function BoundaryNode(; type, index, node)
    return BoundaryNode(type, index, node)
end

"""
```
Boundaries(
    ::FEMParameters{T1,T2,T3,T4,T5} where {T1,T2,T3<:CantileverLoad,T4,T5},
    femMesh::RegularCuboidMesh; 
    kw...
)
```
Creates [`Boundaries`](@ref) for loading a hexahedral cantilever.
"""
function Boundaries(
    ::FEMParameters{T1,T2,T3,T4,T5} where {T1,T2,T3<:CantileverLoad,T4,T5},
    femMesh::RegularCuboidMesh; 
    kw...
)
    numNode = femMesh.numNode
    numNode3 = numNode * 3
    faceNorm = femMesh.faceNorm
    cornerNode = femMesh.cornerNode
    edgeNode = femMesh.edgeNode
    faceNode = femMesh.faceNode
    K = femMesh.K
    
    if !haskey(kw, :uGamma)
        uGamma = BoundaryNode(;
                type = [cornerFE; cornerFE; cornerFE; cornerFE; edgeFE; edgeFE; edgeFE; edgeFE; faceFE],
                index = [:x0y0z0; :x0y1z0; :x0y0z1; :x0y1z1; :y_x0z0; :y_x0z1; :z_x0y0; :z_x0y1; :yz_x0],
                node = [cornerNode[:x0y0z0]; cornerNode[:x0y1z0]; cornerNode[:x0y0z1]; cornerNode[:x0y1z1]; edgeNode[:y_x0z0]; edgeNode[:y_x0z1]; edgeNode[:z_x0y0]; edgeNode[:z_x0y1]; faceNode[:yz_x0]]
            )
    else
        uGamma = kw[:uGamma]
    end

    if !haskey(kw, :mGamma)
        mGamma = BoundaryNode(;
                type = [cornerFE; cornerFE; edgeFE],
                index = [:x1y0z1; :x1y1z1; :y_x1z1],
                node = [cornerNode[:x1y0z1]; cornerNode[:x1y1z1]; edgeNode[:y_x1z1]]
            )
    else
        mGamma = kw[:mGamma]
    end

    if !haskey(kw, :tGamma)
        # Indices
        uGammaLbl = uGamma.type
        mGammaLbl = mGamma.type
        uCorner = findall(x -> x == cornerFE, uGammaLbl)
        uEdge = findall(x -> x == edgeFE, uGammaLbl)
        uFace = findall(x -> x == faceFE, uGammaLbl)
        mCorner = findall(x -> x == cornerFE, mGammaLbl)
        mEdge = findall(x -> x == edgeFE, mGammaLbl)
        mFace = findall(x -> x == faceFE, mGammaLbl)

        uGammaIdx = uGamma.index
        mGammaIdx = mGamma.index
        muCorner = union(uGamma.index[uCorner], mGamma.index[mCorner])
        muEdge = union(uGamma.index[uEdge], mGamma.index[mEdge])
        muFace = union(uGamma.index[uFace], mGamma.index[mFace])

        tCorner = setdiff(keys(cornerNode), muCorner)
        tEdge = setdiff(keys(edgeNode), muEdge)
        tFace = setdiff(keys(faceNode), muFace)
        tGamma = BoundaryNode(;
            type = [fill(cornerFE, length(tCorner)); fill(edgeFE, length(tEdge)); fill(faceFE, length(tFace))],
            index = [tCorner; tEdge; tFace],
            node = setdiff([cornerNode...; collect(Iterators.flatten(edgeNode)); collect(Iterators.flatten(faceNode))], [uGamma.node; mGamma.node])
        )
    else 
        tGamma = kw[:tGamma]
    end
    
    uGammaNode = uGamma.node
    mGammaNode = mGamma.node
    
    !haskey(kw, "uDofs") ? uDofs = copy(reshape(reshape([3 * uGammaNode .- 2; 3 * uGammaNode .- 1; 3 * uGammaNode], :, 3)', :)) : uDofs = kw["uDofs"]
    
    !haskey(kw, "tDofs") ? tDofs = setdiff(1:numNode, uDofs) : tDofs = kw["tDofs"]

    !haskey(kw, "mDofs") ? mDofs = [] : mDofs = setdiff(1:numNode, (uDofs, tDofs))

    C = cholesky(Hermitian(K[tDofs, tDofs]))

    noExit = 5

    boundaries = Boundaries(;
        noExit = noExit,
        uGamma = uGamma,
        tGamma = tGamma,
        mGamma = mGamma,
        uDofs = uDofs,
        tDofs = tDofs,
        mDofs = mDofs,
        tK = C
    )
    
    forceDisplacement = ForceDisplacement(;
        uTilde = spzeros(numNode3), 
        uHat = spzeros(numNode3), 
        u = spzeros(numNode3), 
        fTilde = spzeros(numNode3), 
        fHat = spzeros(numNode3), 
        f = spzeros(numNode3)
    )
    return boundaries, forceDisplacement
end

"""
```
Boundaries(; noExit, uGamma, tGamma, mGamma, uDofs, tDofs, mDofs, tK)
```
Creates [`Boundaries`](@ref).
"""
function Boundaries(; noExit, uGamma, tGamma, mGamma, uDofs, tDofs, mDofs, tK)
    uGammaDln = collect(Iterators.flatten(vcat(uGamma.node, mGamma.node)))
    tGammaDln = collect(Iterators.flatten(vcat(tGamma.node, mGamma.node)))
    uDofsDln = copy(reshape(reshape([3 * uGammaDln .- 2; 3 * uGammaDln .- 1; 3 * uGammaDln], :, 3)', :))
    tDofsDln = copy(reshape(reshape([3 * tGammaDln .- 2; 3 * tGammaDln .- 1; 3 * tGammaDln], :, 3)', :))
    return Boundaries(noExit, uGammaDln, tGammaDln, uDofsDln, tDofsDln, uGamma, tGamma, mGamma, uDofs, tDofs, mDofs, tK)
end