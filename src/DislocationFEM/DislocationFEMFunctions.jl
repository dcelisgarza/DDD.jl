"""
```
```
Calculate the Peach-Koehler force on segments.
"""
function calcPKForce(
    mesh::RegularCuboidMesh,
    dlnFEM::DislocationFEMCorrective,
    network::DislocationNetwork,
)
    # Un normalised segment vectors.
    numSeg = network.numSeg
    idx = network.segIdx
    coord = network.coord
    bVec = network.bVec[idx[:, 1], :]
    node1 = coord[idx[:, 2], :]
    node2 = coord[idx[:, 3], :]
    tVec = node2 - node1
    midNode = 0.5 * (node1 + node2)
    PKForce = zeros(numSeg, 3)

    @simd for i = 1:numSeg
        x0 = midNode[i, :]
        b = bVec[i, :]
        t = tVec[i, :]
        σ_hat = calc_σ_hat(mesh, dlnFEM, x0)
        PKForce[i, :] = (σ_hat * b) × t
    end

    return PKForce
end

"""
```
```
Calculate the reaction from a dislocation.
"""
function calc_σ_hat(
    mesh::RegularCuboidMesh,
    dlnFEM::DislocationFEMCorrective,
    x0::AbstractArray{<:Float64, N} where {N},
)

    # These are just aliases to reduce verbosity, they don't impact performance.
    mx = mesh.numElem[1]    # num elem in x
    my = mesh.numElem[2]    # num elem in y
    mz = mesh.numElem[3]    # num elem in z
    w = 1 / mesh.sizeElem[1]  # 1 / width
    h = 1 / mesh.sizeElem[2]  # 1 / height
    d = 1 / mesh.sizeElem[3]  # 1 / depth
    sizeMesh = mesh.sizeMesh
    C = mesh.stiffTensor
    label = mesh.label
    coord = mesh.coord

    Uhat = dlnFEM.uHat

    x = x0[1]
    y = x0[2]
    z = x0[3]

    # Stress tensor.
    # σ[1,:] = σ_xx
    # σ[2,:] = σ_yy
    # σ[3,:] = σ_zz
    # σ[4,:] = σ_xy = σ_yx
    # σ[5,:] = σ_xz = σ_xz
    # σ[6,:] = σ_yz = σ_zy
    σ = zeros(6)
    B = zeros(6, 24)
    U = zeros(24)

    # Find element index closest to the coordinate.
    i = maximum(ceil(x * w), 1)
    j = maximum(ceil(y * h), 1)
    k = maximum(ceil(z * d), 1)

    # Calculate index of the elements.
    idx = i + (k - 1) * mx + (j - 1) * mx * mz

    # Diametrically opposed points in a cubic element.
    n1 = label[idx, 1]
    n7 = label[idx, 7]
    # Find element midpoints.
    xc = 0.5 * (coord[n1, 1] + coord[n7, 1])
    yc = 0.5 * (coord[n1, 2] + coord[n7, 2])
    zc = 0.5 * (coord[n1, 3] + coord[n7, 3])
    # Setting up Jacobian.
    #=
        # The code is this but simplified for performance.
        a = 0.5*w
        b = 0.5*h
        c = 0.5*d
        s1 = (x0[:,1]-xc[:,1])/a
        s2 = (x0[:,2]-xc[:,2])/h
        s3 = (x0[:,3]-xc[:,3])/d
        ds1dx = 1/a
        ds2dy = 1/h
        ds3dz = 1/d
    =#
    ds1dx = 2 * w
    dx2dy = 2 * h
    dx3dz = 2 * d
    s1 = ds1dx * (x - xc)
    s2 = dx2dy * (y - yc)
    s3 = dx3dz * (z - zc)

    dNdS = shapeFunctionDeriv(LinearQuadrangle3D(), s1, s2, s3)
    dNdS[:, 1] .*= ds1dx
    dNdS[:, 2] .*= ds2dy
    dNdS[:, 3] .*= ds3dz

    @inbounds for i = 1:size(dNdS, 1)
        # Indices calculated once for performance.
        idx1 = 3 * i
        idx2 = 3 * (i - 1)
        # label[a, b] is the index of the node b, in FE element a.
        idx3 = 3 * label[i]
        # Constructing the Jacobian for node i.
        B[1, idx2 + 1] = dNdS[i, 1]
        B[2, idx2 + 2] = dNdS[i, 2]
        B[3, idx2 + 3] = dNdS[i, 3]
        B[4, idx2 + 1] = B[2, idx2 + 2]
        B[4, idx2 + 2] = B[1, idx2 + 1]
        B[5, idx2 + 1] = B[3, idx1 + 0]
        B[5, idx1 + 0] = B[1, idx2 + 1]
        B[6, idx2 + 2] = B[3, idx1 + 0]
        B[6, idx1 + 0] = B[2, idx2 + 2]
        # Building uhat for the nodes of the finite element closest to the point of interest. From idx3, the finite element is the i2'th element and the the node we're looking at is the j'th node. The node index is idx3 = label[i2,j].
        U[idx - 2] = uHat[idx3 - 2]
        U[idx - 1] = uHat[idx3 - 1]
        U[idx - 0] = uHat[idx3 - 0]
    end
    # B*U transforms U from the nodes of the closest finite element with index i2=idx[i] to the point of interest [s1, s2, s3].
    σ = C * B * U
end
