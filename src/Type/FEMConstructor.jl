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
function FEMParameters(;
    type::AbstractMesh,
    order::AbstractElementOrder,
    model::AbstractModel,
    dx,
    dy,
    dz,
    mx,
    my,
    mz,
)
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
function buildMesh(
    matParams::MaterialParameters,
    femParams::FEMParameters{
        F1,
        F2,
        F3,
        F4,
        F5,
    } where {F1 <: DispatchRegularCuboidMesh, F2, F3, F4, F5},
)
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

    https://www.mcs.anl.gov/uploads/cels/papers/P1573A.pdf
"""
function RegularCuboidMesh(
    matParams::MaterialParameters,
    femParams::FEMParameters{
        F1,
        F2,
        F3,
        F4,
        F5,
    } where {F1 <: DispatchRegularCuboidMesh, F2 <: LinearElement, F3, F4, F5},
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
    scale = SVector{3, dxType}(dx, dy, dz) * cbrt(dx * dy * dz)

    # For a regular cuboid mesh this is predefined. Making a volumetric polytope allows us to check if points lie inside the volume simply by doing [x, y, z] ∈ vertices, or checking if they do not belong by doing [x, y, z] ∉ vertices. The symbols are typed as \in and \notin + Tab.
    vtx = SMatrix{3, 8, dxType}(
        0, 0, 0,
        dx, 0, 0,
        dx, dy, 0,
        0, dy, 0,
        0, 0, dz,
        dx, 0, dz,
        dx, dy, dz,
        0, dy, dz,
    )

    vertices = VPolytope(vtx)

    # Faces as defined by the vertices.
    faces = SMatrix{4, 6, mxType}(
        1, 2, 6, 5, # xz plane @ min y
        2, 3, 7, 6, # yz plane @ max x
        3, 4, 8, 7, # xz plane @ max y
        4, 1, 5, 8, # yz plane @ min x
        4, 3, 2, 1, # xy plane @ min z
        5, 6, 7, 8, # xy plane @ max z
    )

    faceMidPt = mean(vtx[:, faces], dims = 2)[:, 1, :]

    # Face normal of the corresponding face.
    faceNorm =
        SMatrix{3, 6, dxType}(
            0, -1, 0,   # xz plane @ min y
            1, 0, 0,    # yz plane @ max x
            0, 1, 0,    # xz plane @ max y
            -1, 0, 0,   # yz plane @ min x
            0, 0, -1,   # xy plane @ min z
            0, 0, 1,    # xy plane @ max z
        )

    coord = zeros(dxType, 3, numNode)           # Node coordinates.
    connectivity = zeros(mxType, 8, numElem)    # Element connectivity.

    cornerNode = SVector{8, Symbol}(
        :x0y0z0,
        :x1y0z0,
        :x1y1z0,
        :x0y1z0,
        :x0y0z1,
        :x1y0z1,
        :x1y1z1,
        :x0y1z1,
    )
    edgeNode = SVector{12, Symbol}(
        :x_y0z0,
        :y_x1z0,
        :x_y1z0,
        :y_x0z0,
        :x_y0z1,
        :y_x1z1,
        :x_y1z1,
        :y_x0z1,
        :z_x0y0,
        :z_x1y0,
        :z_x1y1,
        :z_x0y1,
    )
    faceNode = SVector{6, Symbol}(
        :xz_y0,
        :yz_x1,
        :xz_y1,
        :yz_x0,
        :xy_z0,
        :xy_z1,
    )
    surfNode = (
        # Corners
        x0y0z0 = 1,
        x1y0z0 = mx1,
        x1y1z0 = mx1 * my1,
        x0y1z0 = mx1 * my + 1,
        x0y0z1 = 1 + mx1 * my1 * mz,
        x1y0z1 = mx1 + mx1 * my1 * mz,
        x1y1z1 = mx1 * my1 * mz1,
        x0y1z1 = mx1 * my + 1 + mx1 * my1 * mz,
        # Edges
        x_y0z0 = zeros(mxType, mxm1),
        y_x1z0 = zeros(mxType, mym1),
        x_y1z0 = zeros(mxType, mxm1),
        y_x0z0 = zeros(mxType, mym1),
        x_y0z1 = zeros(mxType, mxm1),
        y_x1z1 = zeros(mxType, mym1),
        x_y1z1 = zeros(mxType, mxm1),
        y_x0z1 = zeros(mxType, mym1),
        z_x0y0 = zeros(mxType, mzm1),
        z_x1y0 = zeros(mxType, mzm1),
        z_x1y1 = zeros(mxType, mzm1),
        z_x0y1 = zeros(mxType, mzm1),
        # Faces
        xz_y0 = zeros(mxType, mxm1 * mzm1),
        yz_x1 = zeros(mxType, mym1 * mzm1),
        xz_y1 = zeros(mxType, mxm1 * mzm1),
        yz_x0 = zeros(mxType, mym1 * mzm1),
        xy_z0 = zeros(mxType, mxm1 * mym1),
        xy_z1 = zeros(mxType, mxm1 * mym1),
    )
    surfNodeArea = (
        # Corners
        x0y0z0 = (h * d + w * d + w * h) / 4,
        x1y0z0 = (h * d + w * d + w * h) / 4,
        x1y1z0 = (h * d + w * d + w * h) / 4,
        x0y1z0 = (h * d + w * d + w * h) / 4,
        x0y0z1 = (h * d + w * d + w * h) / 4,
        x1y0z1 = (h * d + w * d + w * h) / 4,
        x1y1z1 = (h * d + w * d + w * h) / 4,
        x0y1z1 = (h * d + w * d + w * h) / 4,
        # Edges
        x_y0z0 = (h + d) * w / 2,
        y_x1z0 = (w + d) * h / 2,
        x_y1z0 = (h + d) * w / 2,
        y_x0z0 = (w + d) * h / 2,
        x_y0z1 = (h + d) * w / 2,
        y_x1z1 = (w + d) * h / 2,
        x_y1z1 = (h + d) * w / 2,
        y_x0z1 = (w + d) * h / 2,
        z_x0y0 = (h + w) * d / 2,
        z_x1y0 = (h + w) * d / 2,
        z_x1y1 = (h + w) * d / 2,
        z_x0y1 = (h + w) * d / 2,
        # Surfaces
        xz_y0 = w * d,
        yz_x1 = h * d,
        xz_y1 = w * d,
        yz_x0 = h * d,
        xy_z0 = w * h,
        xy_z1 = w * h,
    )
    
    surfNodeNorm = (
        # Corners
        x0y0z0 = normalize(faceNorm[:, 4] + faceNorm[:, 1] + faceNorm[:, 5]),
        x1y0z0 = normalize(faceNorm[:, 2] + faceNorm[:, 1] + faceNorm[:, 5]),
        x1y1z0 = normalize(faceNorm[:, 2] + faceNorm[:, 3] + faceNorm[:, 5]),
        x0y1z0 = normalize(faceNorm[:, 4] + faceNorm[:, 3] + faceNorm[:, 5]),
        x0y0z1 = normalize(faceNorm[:, 4] + faceNorm[:, 1] + faceNorm[:, 6]),
        x1y0z1 = normalize(faceNorm[:, 2] + faceNorm[:, 1] + faceNorm[:, 6]),
        x1y1z1 = normalize(faceNorm[:, 2] + faceNorm[:, 3] + faceNorm[:, 6]),
        x0y1z1 = normalize(faceNorm[:, 4] + faceNorm[:, 3] + faceNorm[:, 6]),
        # Edges
        x_y0z0 = normalize(faceNorm[:, 1] + faceNorm[:, 5]),
        y_x1z0 = normalize(faceNorm[:, 2] + faceNorm[:, 5]),
        x_y1z0 = normalize(faceNorm[:, 3] + faceNorm[:, 5]),
        y_x0z0 = normalize(faceNorm[:, 4] + faceNorm[:, 5]),
        x_y0z1 = normalize(faceNorm[:, 1] + faceNorm[:, 6]),
        y_x1z1 = normalize(faceNorm[:, 2] + faceNorm[:, 6]),
        x_y1z1 = normalize(faceNorm[:, 3] + faceNorm[:, 6]),
        y_x0z1 = normalize(faceNorm[:, 4] + faceNorm[:, 6]),
        z_x0y0 = normalize(faceNorm[:, 4] + faceNorm[:, 1]),
        z_x1y0 = normalize(faceNorm[:, 2] + faceNorm[:, 1]),
        z_x1y1 = normalize(faceNorm[:, 2] + faceNorm[:, 3]),
        z_x0y1 = normalize(faceNorm[:, 4] + faceNorm[:, 3]),
        # Faces
        xz_y0 = faceNorm[:, 1],
        yz_x1 = faceNorm[:, 2],
        xz_y1 = faceNorm[:, 3],
        yz_x0 = faceNorm[:, 4],
        xy_z0 = faceNorm[:, 5],
        xy_z1 = faceNorm[:, 6],
    )

    nodeEl = 1:8 # Local node numbers.
    dofLocal = Tuple(
        Iterators.flatten((
            3 * (nodeEl .- 1) .+ 1,
            3 * (nodeEl .- 1) .+ 2,
            3 * (nodeEl .- 1) .+ 3,
        )),
    )

    E = 2 * μ * (1 + ν)                 # Young's modulus
    Lame = E * νomνInv / (1 - 2 * ν)    # Lamé's relation
    CDiag = Lame + 2 * μ

    # Stiffness tensor.
    C = SMatrix{6, 6, dxType}(
        CDiag,
        Lame,
        Lame,
        0,
        0,
        0,
        Lame,
        CDiag,
        Lame,
        0,
        0,
        0,
        Lame,
        Lame,
        CDiag,
        0,
        0,
        0,
        0,
        0,
        0,
        μ,
        0,
        0,
        0,
        0,
        0,
        0,
        μ,
        0,
        0,
        0,
        0,
        0,
        0,
        μ,
    )

    # Local nodes Gauss Quadrature. Using 8 nodes per element.
    p = 1 / sqrt(3)

    gaussNodes = SMatrix{3, 8, dxType}(
        -p, -p, -p, # x0y0z0
        p, -p, -p,  # x1y0z0
        p, p, -p,   # x1y1z0
        -p, p, -p,  # x0y1z0
        -p, -p, p,  # x0y0z1
        p, -p, p,   # x1y0z1
        p, p, p,    # x1y1z1
        -p, p, p,   # x0y1z1
    )

    # Shape functions and their derivatives.
    N = shapeFunction(
        LinearQuadrangle3D(),
        gaussNodes[1, :],
        gaussNodes[2, :],
        gaussNodes[3, :],
    )
    dNdS = shapeFunctionDeriv(
        LinearQuadrangle3D(),
        gaussNodes[1, :],
        gaussNodes[2, :],
        gaussNodes[3, :],
    )

    # Fill node coordinates.
    @inbounds @simd for k in 1:mz1
        km1 = k - 1
        for j in 1:my1
            jm1 = j - 1
            for i in 1:mx1
                globalNode = i + jm1 * mx1 + km1 * mx1 * my1
                coord[1, globalNode] = (i - 1) * w
                coord[2, globalNode] = jm1 * h
                coord[3, globalNode] = km1 * d
            end
        end
    end

    # Fill element connectivity.
    @inbounds @simd for k in 1:mz
        km1 = k - 1
        for j in 1:my
            jm1 = j - 1
            for i in 1:mx
                globalElem = i + jm1 * mx + km1 * mx * my
                connectivity[1, globalElem] = i + jm1*mx1 + km1*mx1*my1
                connectivity[2, globalElem] = connectivity[1, globalElem] + 1
                connectivity[4, globalElem] = i + j*mx1 + km1*mx1*my1
                connectivity[3, globalElem] = connectivity[4, globalElem] + 1
                connectivity[5, globalElem] = i + jm1*mx1 + k*mx1*my1
                connectivity[6, globalElem] = connectivity[5, globalElem] + 1
                connectivity[8, globalElem] = i + j*mx1 + k*mx1*my1
                connectivity[7, globalElem] = connectivity[8, globalElem] + 1
            end
        end
    end

    # Edge node along x
    @inbounds @simd for i in 1:mxm1
        surfNode[:x_y0z0][i] = i + 1                              # x_y0z0
        surfNode[:x_y1z0][i] = i + 1 + mx1 * my                   # x_y1z0
        surfNode[:x_y0z1][i] = i + 1 + mx1 * my1 * mz             # x_y0z1
        surfNode[:x_y1z1][i] = i + 1 + mx1 * my1 * mz + mx1 * my  # x_y1z1
    end
    # Edge node along y
    @inbounds @simd for i in 1:mym1
        surfNode[:y_x0z0][i] = i * mx1 + 1                      # y_x0z0
        surfNode[:y_x1z0][i] = i * mx1 + mx1                    # y_x1z0
        surfNode[:y_x0z1][i] = i * mx1 + 1 + mx1 * my1 * mz     # y_x0z1
        surfNode[:y_x1z1][i] = i * mx1 + mx1 + mx1 * my1 * mz   # y_x1z1
    end
    # Edge node along z
    @inbounds @simd for i in 1:mzm1
        surfNode[:z_x0y0][i] = i * mx1 * my1 + 1                # z_x0y0
        surfNode[:z_x1y0][i] = i * mx1 * my1 + mx1              # z_x1y0
        surfNode[:z_x0y1][i] = i * mx1 * my1 + 1 + mx1 * my     # z_x0y1
        surfNode[:z_x1y1][i] = i * mx1 * my1 + mx1 + mx1 * my   # z_x1y1
    end

    # Face node xy
    @inbounds @simd for j in 1:mym1
        jm1 = j - 1
        for i in 1:mxm1
            surfNode[:xy_z0][i + mxm1 * jm1] = j * mx1 + 1 + i                  # xy_z0
            surfNode[:xy_z1][i + mxm1 * jm1] = j * mx1 + 1 + i + mx1 * my1 * mz # xy_z1
        end
    end

    # Face node xz
    @inbounds @simd for j in 1:mzm1
        jm1 = j - 1
        for i in 1:mxm1
            surfNode[:xz_y0][i + mxm1 * jm1] = j * mx1 * my1 + 1 + i            # xz_y0
            surfNode[:xz_y1][i + mxm1 * jm1] = j * mx1 * my1 + 1 + i + mx1 * my # xz_y1
        end
    end
    
    # Face node yz
    @inbounds @simd for j in 1:mzm1
        jm1 = j - 1
        for i in 1:mym1
            surfNode[:yz_x0][i + mym1 * jm1] = j * mx1 * my1 + 1 + i * mx1  # yz_x0
            surfNode[:yz_x1][i + mym1 * jm1] = j * mx1 * my1 + mx1 + i * mx1  # yz_x0
        end
    end

    localCoord = SMatrix{3, 8, dxType}(
        coord[1, connectivity[1, 1]],
        coord[2, connectivity[1, 1]],
        coord[3, connectivity[1, 1]],
        coord[1, connectivity[2, 1]],
        coord[2, connectivity[2, 1]],
        coord[3, connectivity[2, 1]],
        coord[1, connectivity[3, 1]],
        coord[2, connectivity[3, 1]],
        coord[3, connectivity[3, 1]],
        coord[1, connectivity[4, 1]],
        coord[2, connectivity[4, 1]],
        coord[3, connectivity[4, 1]],
        coord[1, connectivity[5, 1]],
        coord[2, connectivity[5, 1]],
        coord[3, connectivity[5, 1]],
        coord[1, connectivity[6, 1]],
        coord[2, connectivity[6, 1]],
        coord[3, connectivity[6, 1]],
        coord[1, connectivity[7, 1]],
        coord[2, connectivity[7, 1]],
        coord[3, connectivity[7, 1]],
        coord[1, connectivity[8, 1]],
        coord[2, connectivity[8, 1]],
        coord[3, connectivity[8, 1]],
    )

    localdNdS = SMatrix{8, 3, dxType}(
        dNdS[1][1, 1],
        dNdS[1][1, 2],
        dNdS[1][1, 3],
        dNdS[1][1, 4],
        dNdS[1][1, 5],
        dNdS[1][1, 6],
        dNdS[1][1, 7],
        dNdS[1][1, 8],
        dNdS[1][2, 1],
        dNdS[1][2, 2],
        dNdS[1][2, 3],
        dNdS[1][2, 4],
        dNdS[1][2, 5],
        dNdS[1][2, 6],
        dNdS[1][2, 7],
        dNdS[1][2, 8],
        dNdS[1][3, 1],
        dNdS[1][3, 2],
        dNdS[1][3, 3],
        dNdS[1][3, 4],
        dNdS[1][3, 5],
        dNdS[1][3, 6],
        dNdS[1][3, 7],
        dNdS[1][3, 8],
    )

    J = localCoord * localdNdS
    detJ = det(J)
    invJ = inv(J)'

    B = zeros(dxType, 6, 24, 8)     # Jacobian matrix.
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
            # B[1, idx2, q] = 0
            # B[1, idx3, q] = 0

            # B[2, idx1, q] = 0
            B[2, idx2, q] = nx[2, a]
            # B[2, idx3, q] = 0

            # B[3, idx1, q] = 0
            # B[3, idx2, q] = 0
            B[3, idx3, q] = nx[3, a]

            B[4, idx1, q] = nx[2, a]
            B[4, idx2, q] = nx[1, a]
            # B[4, idx3, q] = 0

            B[5, idx1, q] = nx[3, a]
            # B[5, idx2, q] = 0
            B[5, idx3, q] = nx[1, a]

            # B[6, idx1, q] = 0
            B[6, idx2, q] = nx[3, a]
            B[6, idx3, q] = nx[2, a]
        end
    end

    localK = zeros(dxType, 24, 24)  # K for an element.
    @inbounds @simd for q in nodeEl
        localK .+= B[:, :, q]' * C * B[:, :, q]
    end

    localK = localK * detJ
    cntr = 0
    V1 = zeros(dxType, numElem * 24^2)
    V2 = zeros(dxType, numElem * 24^2)
    V3 = zeros(dxType, numElem * 24^2)
    @inbounds @simd for p in 1:numElem
        globalNode = connectivity[nodeEl, p] # Global node numbers
        dofGlobal = Tuple(
            Iterators.flatten((
                3 * (globalNode .- 1) .+ 1,
                3 * (globalNode .- 1) .+ 2,
                3 * (globalNode .- 1) .+ 3,
            )),
        ) # Global degree of freedom
        for i in 1:24
            for j in 1:24
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

    # Surface elements
    mxmy = mx * my
    mxmz = mx * mz
    mymz = my * mz
    surfElemNode = zeros(mxType, 2 * (mxmy + mxmz + mymz), 4)
    cntr = 0
    
    @inbounds @simd for i in 1:size(faces, 2)
        label = vec(connectivity[faces[:, i], :]')
        
        if i == 1 || i == 3
            xyz = 2
            dimCoord = @view coord[xyz, label]
            i == 1 ? lim = minimum(dimCoord) : lim = maximum(dimCoord)
        elseif i == 2 || i == 4
            xyz = 1
            dimCoord = @view coord[xyz, label]
            i == 4 ? lim = minimum(dimCoord) : lim = maximum(dimCoord)
        else
            xyz = 3
            dimCoord = @view coord[xyz, label]
            i == 5 ? lim = minimum(dimCoord) : lim = maximum(dimCoord)
        end

        idx = findall(x -> x ≈ lim, dimCoord)
        n = div(length(idx), 4)

        label = label[idx]

        surfElemNode[(1:n) .+ cntr, :] = reshape(label, :, 4)
        cntr += n
    end
    

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
        surfNode,
        surfNodeArea,
        surfNodeNorm,
        surfElemNode,
        coord,
        connectivity,
        globalK,
    )
end

"""
```
BoundaryNode(; index, node)
```
Create [`BoundaryNode`](@ref).
"""
function BoundaryNode(; index, node)
    return BoundaryNode(index, node)
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
    ::FEMParameters{T1, T2, T3, T4, T5} where {T1, T2, T3 <: CantileverLoad, T4, T5},
    femMesh::RegularCuboidMesh;
    kw...,
)
    numNode = femMesh.numNode
    numNode3 = numNode * 3
    faceNorm = femMesh.faceNorm
    surfNode = femMesh.surfNode
    K = femMesh.K

    if !haskey(kw, :uGamma)
        uIndex = SVector{9, Symbol}(
            [
                :x0y0z0
                :x0y1z0
                :x0y0z1
                :x0y1z1
                :y_x0z0
                :y_x0z1
                :z_x0y0
                :z_x0y1
                :yz_x0
            ],
        )
        uNodes = collect(Iterators.flatten([surfNode[uIndex[i]] for i in 1:length(uIndex)]))
        uGamma = BoundaryNode(; index = uIndex, node = uNodes)
    else
        uGamma = kw[:uGamma]
    end

    if !haskey(kw, :mGamma)
        mIndex = SVector{3, Symbol}([:x1y0z1; :x1y1z1; :y_x1z1])
        mNodes = collect(Iterators.flatten([surfNode[mIndex[i]] for i in 1:length(mIndex)]))
        mGamma = BoundaryNode(; index = mIndex, node = mNodes)
    else
        mGamma = kw[:mGamma]
    end

    if !haskey(kw, :tGamma)
        tIndex = SVector{14, Symbol}(
            Symbol.(setdiff(string.(keys(surfNode)), string.([uIndex; mIndex]))),
        )
        tNodes = collect(Iterators.flatten([surfNode[tIndex[i]] for i in 1:length(tIndex)]))
        tGamma = BoundaryNode(; index = tIndex, node = tNodes)
    else
        tGamma = kw[:tGamma]
    end

    uGammaNode = uGamma.node
    mGammaNode = mGamma.node
    tGammaNode = tGamma.node

    !haskey(kw, "uDofs") ?
    uDofs =
        sort!([3 * uGammaNode .- 2; 3 * uGammaNode .- 1; 3 * uGammaNode; 3 * mGammaNode]) :
    uDofs = kw["uDofs"]

    !haskey(kw, "tDofs") ? tDofs = sort!(setdiff(1:numNode3, uDofs)) : tDofs = kw["tDofs"]

    C = cholesky(Hermitian(K[tDofs, tDofs]))

    noExit = 5

    boundaries = Boundaries(;
        noExit = noExit,
        uGamma = uGamma,
        tGamma = tGamma,
        mGamma = mGamma,
        uDofs = uDofs,
        tDofs = tDofs,
        mDofs = nothing,
        tK = C,
    )

    forceDisplacement = ForceDisplacement(;
        uTilde = spzeros(numNode3),
        uHat = spzeros(numNode3),
        u = spzeros(numNode3),
        fTilde = spzeros(numNode3),
        fHat = spzeros(numNode3),
        f = spzeros(numNode3),
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
    uGammaDln = sort!([uGamma.node; mGamma.node])
    tGammaDln = sort!([tGamma.node; mGamma.node])
    uDofsDln = sort!([3 * uGammaDln .- 2; 3 * uGammaDln .- 1; 3 * uGammaDln])
    tDofsDln = sort!([3 * tGammaDln .- 2; 3 * tGammaDln .- 1; 3 * tGammaDln])
    return Boundaries(
        noExit,
        uGammaDln,
        tGammaDln,
        uDofsDln,
        tDofsDln,
        uGamma,
        tGamma,
        mGamma,
        uDofs,
        tDofs,
        mDofs,
        tK,
    )
end
