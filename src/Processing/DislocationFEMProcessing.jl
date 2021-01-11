"""
```
calc_σHat(
    mesh::RegularCuboidMesh,
    dlnFEM::DislocationFEMCorrective,
    x0::AbstractArray{T, N} where {T, N},
)
```
Calculate the reaction from a dislocation.
"""
function calc_σHat(
    mesh::RegularCuboidMesh,
    forceDisplacement::ForceDisplacement,
    x0::AbstractVector{T} where {T},
)

    # Unroll structure.
    mx = mesh.mx        # num elem in x
    my = mesh.my        # num elem in y
    mz = mesh.mz        # num elem in z
    wInv = 1 / mesh.w    # 1 / width
    hInv = 1 / mesh.h    # 1 / height
    dInv = 1 / mesh.d    # 1 / depth
    C = mesh.C
    connectivity = mesh.connectivity
    coord = mesh.coord
    elemT = eltype(coord)

    uHat = forceDisplacement.uHat

    x, y, z = x0

    # Find element index closest to the coordinate.
    i::Int = max(ceil(x * wInv), 1)
    j::Int = max(ceil(y * hInv), 1)
    k::Int = max(ceil(z * dInv), 1)

    # Calculate index of the elements.
    idx = i + (k - 1) * mx + (j - 1) * mx * mz

    # Diametrically opposed points in a cubic element.
    n1 = connectivity[1, idx]
    n7 = connectivity[7, idx]

    # Find element midpoints.
    xc = 0.5 * (coord[1, n1] + coord[1, n7])
    yc = 0.5 * (coord[2, n1] + coord[2, n7])
    zc = 0.5 * (coord[3, n1] + coord[3, n7])

    # Setting up Jacobian.
    ds1dx = 2 * wInv
    ds2dy = 2 * hInv
    ds3dz = 2 * dInv

    s1 = (x - xc) * ds1dx
    s2 = (y - yc) * ds2dy
    s3 = (z - zc) * ds3dz

    pm1 = (-1, 1, 1, -1, - 1, 1, 1, -1)
    pm2 = (1, 1, 1, 1, -1, -1, -1, -1)
    pm3 = (-1, -1, 1, 1, -1, -1, 1, 1)

    ds1dx *= 0.125
    ds2dy *= 0.125
    ds3dz *= 0.125

    B = zeros(elemT, 6, 24)
    U = MVector{24,elemT}(zeros(24))
    @inbounds @simd for i in 1:8
        dNdS1 = pm1[i] * (1 + pm2[i] * s2) * (1 + pm3[i] * s3) * ds1dx
        dNdS2 = (1 + pm1[i] * s1) * pm2[i] * (1 + pm3[i] * s3) * ds2dy
        dNdS3 = (1 + pm1[i] * s1) * (1 + pm2[i] * s2) * pm3[i] * ds3dz
        # Indices calculated once for performance.
        idx1 = 3 * i
        idx2 = 3 * (i - 1)
        # Linear index of the z-coordinate of node i of the idx'th FE element.
        idx3 = 3 * connectivity[i, idx]

        # Constructing the Jacobian for node i.
        B[1, idx2 + 1] = dNdS1
        B[2, idx2 + 2] = dNdS2
        B[3, idx2 + 3] = dNdS3

        B[4, idx2 + 1] = B[2, idx2 + 2]
        B[4, idx2 + 2] = B[1, idx2 + 1]

        B[5, idx2 + 1] = B[3, idx1 + 0]
        B[5, idx1 + 0] = B[1, idx2 + 1]

        B[6, idx2 + 2] = B[3, idx1 + 0]
        B[6, idx1 + 0] = B[2, idx2 + 2]

        # Building uhat for the nodes of the finite element closest to the point of interest. From idx3, the finite element is the i2'th element and the the node we're looking at is the j'th node. The node index is idx3 = label[i2,j].
        U[idx1 - 2] = uHat[idx3 - 2]
        U[idx1 - 1] = uHat[idx3 - 1]
        U[idx1 - 0] = uHat[idx3 - 0]
    end

    # Isotropic stress tensor in vector form.
    # σ_vec[1] = σ_xx
    # σ_vec[2] = σ_yy
    # σ_vec[3] = σ_zz
    # σ_vec[4] = σ_xy = σ_yx
    # σ_vec[5] = σ_xz = σ_xz
    # σ_vec[6] = σ_yz = σ_zy
    # B*U transforms U from the nodes of the closest finite element with index i2=idx[i] to the point of interest [s1, s2, s3].
    σ_vec = C * B * U
    σ = SMatrix{3,3,elemT}(
        σ_vec[1], σ_vec[4], σ_vec[5],
        σ_vec[4], σ_vec[2], σ_vec[6],
        σ_vec[5], σ_vec[6], σ_vec[3],
    )

    return σ
end
