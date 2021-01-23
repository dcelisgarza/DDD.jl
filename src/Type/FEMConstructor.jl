"""
```
FEMParameters(
    type::T1,
    order::T2,
    model::T3,
    dx::T4,
    dy::T4,
    dz::T4,
    mx::T5,
    my::T5,
    mz::T5
) where {
    T1 <: AbstractMesh,
    T2 <: AbstractModel,
    T3 <: AbstractElementOrder,
    T4 <: AbstractFloat,
    T5 <: Integer
}
```
Constructor for [`FEMParameters`](@ref).
"""
function FEMParameters(type::T1, order::T2, model::T3, dx::T4, dy::T4, dz::T4, mx::T5, my::T5, mz::T5) where {T1 <: AbstractMesh,T2 <: AbstractElementOrder,T3 <: AbstractModel, T4 <: AbstractFloat,T5 <: Integer}
    FEMParameters{T1,T2,T3,T4,T5}(type, order, model, dx, dy, dz, mx, my, mz)
end
"""
```
FEMParameters(;
    type::T1,
    order::T2,
    model::T3,
    dx::T4,
    dy::T4,
    dz::T4,
    mx::T5,
    my::T5,
    mz::T5
) where {
    T1 <: AbstractMesh,
    T2 <: AbstractModel,
    T3 <: AbstractElementOrder,
    T4 <: AbstractFloat,
    T5 <: Integer
}
```
Keyword constructor for [`FEMParameters`](@ref). Calls the positional constructor.
"""
function FEMParameters(; type::T1, order::T2, model::T3, dx::T4, dy::T4, dz::T4, mx::T5, my::T5, mz::T5) where {T1 <: AbstractMesh,T2 <: AbstractElementOrder,T3<:AbstractModel, T4 <: AbstractFloat,T5 <: Integer}
    return FEMParameters(type, order, model, dx, dy, dz, mx, my, mz)
end
"""
```
Boundaries(; uGamma, tGamma, mGamma, uDofs, tDofs, mDofs)
```
Boundary condition keyword constructor.
"""
function Boundaries(; uGamma, tGamma, mGamma, uDofs, tDofs, mDofs)
    return Boundaries(uGamma, tGamma, mGamma, uDofs, tDofs, mDofs)
end
"""
```
buildMesh(matParams::T1, femParams::T2) 
    where {T1 <: MaterialParameters,T2 <: FEMParameters}
```
Internally calls another `buildMesh` that dispatches based on [`FEMParameters`](@ref).`type`.
"""
function buildMesh(matParams::T1, femParams::T2) where {T1 <: MaterialParameters,T2 <: FEMParameters}
    return buildMesh(femParams.type, matParams, femParams)
end
"""
```
buildMesh(::T1, matParams::T2, femParams::T3) 
    where {T1 <: DispatchRegularCuboidMesh,T2 <: MaterialParameters,T3 <: FEMParameters}
```
Builds a [`RegularCuboidMesh`](@ref) by dispatching on `femParams.type == DispatchRegularCuboidMesh()`.
"""
function buildMesh(::T1, matParams::T2, femParams::T3) where {T1 <: DispatchRegularCuboidMesh,T2 <: MaterialParameters,T3 <: FEMParameters}
    return RegularCuboidMesh(femParams.order, matParams, femParams)
end

"""
```
RegularCuboidMesh(
    order::T1, 
    matParams::T2, 
    femParams::T3
) where {
    T1 <: LinearElement,
    T2 <: MaterialParameters,
    T3 <: FEMParameters
}
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
function RegularCuboidMesh(order::T1, matParams::T2, femParams::T3) where {T1 <: LinearElement,T2 <: MaterialParameters,T3 <: FEMParameters}

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
    mxm1 = mx - 1
    mym1 = my - 1
    mzm1 = mz - 1
    mx1mz1 = mx1 * mz1
    mymx1mz1 = my * mx1mz1
    mzmx1 = mz * mx1
    numNode = mx1 * my1 * mz1

    coord = zeros(dxType, 3, numNode)           # Node coordinates.
    connectivity = zeros(mxType, 8, numElem)    # Element connectivity.
    
    # Nodes corresponding to the vertices.
    cornerNode = SVector{8,mxType}(1, mx1, 1 + mymx1mz1, mx1 + mymx1mz1, 1 + mzmx1, mx1mz1, 
                                    1 + mymx1mz1 + mzmx1, mx1 + mymx1mz1 + mzmx1)
    # Nodes corresponding to the edges.
    edgeNode = (zeros(mxType, mxm1), zeros(mxType, mym1), zeros(mxType, mxm1), zeros(mxType, mym1), 
                zeros(mxType, mxm1), zeros(mxType, mym1), zeros(mxType, mxm1), zeros(mxType, mym1),
                zeros(mxType, mzm1), zeros(mxType, mzm1), zeros(mxType, mzm1), zeros(mxType, mzm1))
    # Nodes corresponding to the faces.
    faceNode = (zeros(mxType, mxm1 * mym1), zeros(mxType, mxm1 * mym1),
                zeros(mxType, mxm1 * mzm1), zeros(mxType, mxm1 * mzm1),
                zeros(mxType, mym1 * mzm1), zeros(mxType, mym1 * mzm1))

    localK = zeros(dxType, 24, 24)  # K for an element.
    B = zeros(dxType, 6, 24, 8)     # Jacobian matrix.

    nodeEl = 1:8 # Local node numbers.
    dofLocal = Tuple(Iterators.flatten((3 * (nodeEl .- 1) .+ 1, 3 * (nodeEl .- 1) .+ 2, 3 * (nodeEl .- 1) .+ 3)))
    
    V1 = zeros(dxType, numElem * 24^2)
    V2 = zeros(dxType, numElem * 24^2)
    V3 = zeros(dxType, numElem * 24^2)

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
        edgeNode[1][i] = ip1
        edgeNode[3][i] = ip1 + mymx1mz1
        edgeNode[5][i] = ip1 + mzmx1
        edgeNode[7][i] = ip1 + mzmx1 + mymx1mz1
    end
    # Edge node along y
    @inbounds @simd for i in 1:mym1
        imx1mz1 = i * mx1mz1
        edgeNode[2][i] = 1 + imx1mz1
        edgeNode[4][i] = mx1 + imx1mz1
        edgeNode[6][i] = 1 + imx1mz1 + mzmx1
        edgeNode[8][i] = mx1 + imx1mz1 + mzmx1
    end
    # Edge node along z
    @inbounds @simd for i in 1:mzm1
        ip1 = i + 1
        ip1mx1 = ip1 * mx1
        imx1 = i * mx1
        edgeNode[9][i] = 1 + imx1
        edgeNode[10][i] = ip1mx1
        edgeNode[11][i] = 1 + mymx1mz1 + imx1
        edgeNode[12][i] = mymx1mz1 + ip1mx1
    end
    # Face node xy
    @inbounds @simd for j in 1:mym1
        jm1 = j - 1
        jmx1mz1 = j * mx1mz1
        mxm1jm1 = mxm1 * jm1
        jmx1mz1p1 = jmx1mz1 + 1
        for i in 1:mxm1
            ipmxm1jm1 = i + mxm1jm1
            faceNode[1][ipmxm1jm1] = jmx1mz1p1 + i
            faceNode[2][ipmxm1jm1] = jmx1mz1p1 + i + mzmx1
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
            faceNode[3][ipmxm1jm1] = jmx1p1 + i
            faceNode[4][ipmxm1jm1] = jmx1p1 + i + mymx1mz1
        end
    end

    @inbounds @simd for j in 1:mym1
        jm1 = j - 1
        jmx1mz1 = j * mx1mz1
        for i in 1:mzm1
            imx1 = i * mx1
            jmx1mz1pimx1 = jmx1mz1 + imx1
            faceNode[5][i + jm1*mzm1] = 1 + jmx1mz1pimx1
            faceNode[6][i + jm1*mzm1] = mx1 + jmx1mz1pimx1
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
        vertices,
        faces,
        faceMidPt,
        faceNorm,
        C,
        dx,
        dy,
        dz,
        scale,
        mx,
        my,
        mz,
        numElem,
        numNode,
        w,
        h,
        d,
        coord,
        connectivity,
        cornerNode,
        edgeNode,
        faceNode,
        globalK,
)
    
end
"""
```
RegularCuboidMesh(;
    order::T1, 
    matParams::T2, 
    femParams::T3
) where {
    T1 <: AbstractElementOrder,
    T2 <: MaterialParameters,
    T3 <: FEMParameters
}
```
Keyword constructor for [`RegularCuboidMesh`](@ref).
"""
function RegularCuboidMesh(; order::T1, matParams::T2, femParams::T3) where {T1 <: AbstractElementOrder,T2 <: MaterialParameters,T3 <: FEMParameters}
    return RegularCuboidMesh(order, matParams, femParams)
end

"""
```
Boundaries(femParams::T1, femMesh::T2, args...; kw...) where {T1 <: FEMParameters,T2 <: AbstractMesh}
```
Predefined boundaries.
"""
function Boundaries(femParams::T1, femMesh::T2, args...; kw...) where {T1 <: FEMParameters,T2 <: AbstractMesh}
    return Boundaries(femParams.model, femMesh, args...; kw...)
end
"""
```
Boundaries(::T1, femMesh::T2, args...; kw...) where {T1 <: CantileverLoad,T2 <: RegularCuboidMesh}
```
Cantilever loading boundaries.
"""
function Boundaries(::T1, femMesh::T2, args...; kw...) where {T1 <: CantileverLoad,T2 <: RegularCuboidMesh}
    numNode = femMesh.numNode
    faceNorm = femMesh.faceNorm
    cornerNode = femMesh.cornerNode
    edgeNode = femMesh.edgeNode
    faceNode = femMesh.faceNode

    !haskey(kw, "uGamma") ? uGamma = [cornerNode[1]; cornerNode[3]; cornerNode[5]; cornerNode[7]; edgeNode[2]; edgeNode[6]; edgeNode[9]; edgeNode[11]; faceNode[5]] : uGamma = kw["uGamma"]

    !haskey(kw, "mGamma") ? mGamma = [cornerNode[6]; cornerNode[8]; edgeNode[8]] : mGamma = kw["mGamma"]

    !haskey(kw, "tGamma") ? tGamma = setdiff([cornerNode; edgeNode; faceNode], [uGamma; mGamma]) : tGamma = kw["tGamma"]
    
    !haskey(kw, "uDofs") ? uDofs = [3 * uGamma .- 2; 3 * uGamma .- 1; 3 * uGamma; 3 * mGamma] : uDofs = kw["uDofs"]
    
    !haskey(kw, "tDofs") ? tDofs = setdiff(1:numNode, uDofs) : tDofs = kw["tDofs"]

    !haskey(kw, "mDofs") ? mDofs = nothing : mDofs = setdiff(1:numNode, (uDofs, tDofs))

    return Boundaries(; uGamma = uGamma, tGamma = tGamma, mGamma = mGamma, uDofs = uDofs, tDofs = tDofs, mDofs = mDofs)
end