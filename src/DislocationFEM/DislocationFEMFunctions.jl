"""
```
```
Calculate the Peach-Koehler force on segments.
"""
function pkForce(
    dlnFE::DislocationFEMCorrective,
    mesh::RegularCuboidMesh,
    segments::AbstractArray{<:Float64, N} where {N},
) end

"""
```
```
"""
function hatStress(
    mesh::RegularCuboidMesh,
    dlnFEM::DislocationFEMCorrective,
    x0::AbstractArray{<:Float64, N} where {N},
)

    # These are just aliases to reduce verbosity, they don't impact performance.
    mx = mesh.numElem[1]    # num elem in x
    my = mesh.numElem[2]    # num elem in y
    mz = mesh.numElem[3]    # num elem in z
    w = mesh.sizeElem[1]    # width
    h = mesh.sizeElem[2]    # height
    d = mesh.sizeElem[3]    # depth
    sizeMesh = mesh.SizeMesh
    C = mesh.stiffTensor
    label = mesh.label
    coord = mesh.coord
    Uhat = dlnFEM.uHat

    # Find element index closest to the coordinate.
    i = ceil(x0[:, 1] / w)
    j = ceil(x0[:, 2] / h)
    k = ceil(x0[:, 3] / d)
    # Correct index if indices end up being outside of the simulation domain
    i[findall(x -> x == 0, i)] .= 1
    i[findall(x -> x > mx, i)] .= mx
    j[findall(x -> x == 0, j)] .= 1
    j[findall(x -> x > my, j)] .= my
    k[findall(x -> x == 0, k)] .= 1
    k[findall(x -> x > mz, k)] .= mz

    # Calculate index of the elements.
    idx = i + (k - 1) * mx + (j - 1) * mx * mz

    # Diametrically opposed points in a cubic element.
    opp = [label[idx, 1], label[idx, 7]]
    # Find element midpoints.
    xc =
        [
            (coord[opp[:, 1], 1] + coord[opp[:, 2], 1]),
            (coord[opp[:, 1], 2] + coord[opp[:, 2], 2]),
            (coord[opp[:, 1], 3] + coord[opp[:, 2], 3]),
        ] * 0.5

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
    ds1dx = 2 / w
    dx2dy = 2 / h
    dx3dz = 2 / d
    s1 = ds1dx * (x0[:, 1] - xc[:, 1])
    s2 = dx2dy * (x0[:, 2] - xc[:, 2])
    s3 = dx3dz * (x0[:, 3] - xc[:, 3])
    # Sign vectors.
    sv1 = [-1; 1; 1; -1; -1; 1; 1; -1]
    sv2 = [1; 1; 1; 1; -1; -1; -1; -1]
    sv3 = [-1; -1; 1; 1; -1; -1; 1; 1]

    # Stress tensor.
    # sigma[1,:] = σ_xx
    # sigma[2,:] = σ_yy
    # sigma[3,:] = σ_zz
    # sigma[4,:] = σ_xy = σ_yx
    # sigma[5,:] = σ_xz = σ_xz
    # sigma[6,:] = σ_yz = σ_zy
    if ndims(s1) > 1
        sigma = zeros(6, length(s1))
    else
        sigma = zeros(6)
    end

    dNds1 = zeros(8)
    dNds2 = zeros(8)
    dNds3 = zeros(8)
    B = zeros(6, 24)
    U = zeros(24)

    ds1dx /= 8
    dx2dy /= 8
    dx3dz /= 8
    # Traverse s1, s2 and s3 along their indices in the fastest way, works for any dimensionality. Realistically we'll only use 3D (volume), 2D (surface), 1D (line) or 0D (point) 'meshes'.
    @inbounds for i in eachindex(s1)
        # Index for the surface element nearest the point of interest.
        i2 = idx[i]
        # Traverse the shape functions.
        for j in eachindex(sv1)
            # Indices calculated once for performance.
            idx1 = 3 * j
            idx2 = 3 * (j - 1)
            # label[a, b] is the index of the node b, in FE element a.
            idx3 = 3 * label[i2, j]
            # Derivatives of the shape functions for node j and point i.
            dNds1[j] = sv1[j] * (1 + sv2[j] * s2[i]) * (1 + sv3[j] * s3[i])
            dNds2[j] = sv2[j] * (1 + sv1[j] * s1[i]) * (1 + sv3[j] * s3[i])
            dNds3[j] = sv3[j] * (1 + sv1[j] * s1[i]) * (1 + sv2[j] * s2[i])
            # Constructing the Jacobian for node j.
            B[1, idx2 + 1] = dNds1[j] * ds1dx
            B[2, idx2 + 2] = dNds2[j] * ds2dy
            B[3, idx2 + 3] = dNds3[j] * ds3dz
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
        # B*U transforms U from the nodes of the closest finite element with index i2=idx[i] to the point of interest [s1[i], s2[i], s3[i]].
        sigma[:, i] = C * B * U
    end
end
