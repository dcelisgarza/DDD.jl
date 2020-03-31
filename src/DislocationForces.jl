"""
```
```
Calculate the Peach-Koehler force on segments.
"""
function pkForce(
    dlnFE::DislocationFEM,
    mesh::CuboidMesh,
    segments::AbstractArray{<:Float64, N} where {N},
) end

"""
```
```
"""
function hatStress(
    mesh::RegularCuboidMesh,
    x0::AbstractArray{<:Float64, N} where {N},
)

    # These are just aliases to reduce verbosity, they don't impact performance.
    mx = mesh.numElem[1]    # num elem in x
    my = mesh.numElem[2]    # num elem in y
    mz = mesh.numElem[3]    # num elem in z
    w = mesh.sizeElem[1]    # width
    h = mesh.sizeElem[2]    # height
    d = mesh.sizeElem[3]    # depth
    C = mesh.stiffTensor
    label = mesh.label
    coord = mesh.coord

    # Find element index closest to the coordinate.
    i = ceil(x0[:, 1] / w)
    j = ceil(x0[:, 2] / h)
    k = ceil(x0[:, 3] / d)
    # Correct index if the result is under 1. Julia is 1-indexed, if using floor we'd have to correct the last index.
    i[findall(x -> x == 0, i)] .= 1
    j[findall(x -> x == 0, i)] .= 1
    k[findall(x -> x == 0, i)] .= 1

    # Calculate index of the elements.
    idx = i + (k - 1) * mx + (j - 1) * mx * mz

    # Find element midpoints.
    lbl = [label[idx, 1]; label[idx, 7]]
    midCoord = [
        0.5 * (coord[lbl[1], 1] + coord[lbl[2], 1]),
        0.5 * (coord[lbl[1], 2] + coord[lbl[2], 2]),
        0.5 * (coord[lbl[1], 3] + coord[lbl[2], 3]),
    ]

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
    # Derivatives of the shape functions.
    if ndims(s1) > 1
        dNds1 = zeros(8, nPoints)
        dNds2 = zeros(8, nPoints)
        dNds3 = zeros(8, nPoints)
        B = zeros(6, 24, nPoints)
    else
        dNds1 = zeros(8)
        dNds2 = zeros(8)
        dNds3 = zeros(8)
        B = zeros(6, 24)
    end

    ds1dx /= 8
    dx2dy /= 8
    dx3dz /= 8
    @inbounds for i in eachindex(s1)
        for j in eachindex(sv1)
            dNds1[j, i] = sv1[j] * (1 + sv2[j] * s2[i]) * (1 + sv3[j] * s3[i])
            dNds2[j, i] = sv2[j] * (1 + sv1[j] * s1[i]) * (1 + sv3[j] * s3[i])
            dNds3[j, i] = sv3[j] * (1 + sv1[j] * s1[i]) * (1 + sv2[j] * s2[i])

            idx1 = 3 * j
            idx2 = 3 * (j - 1)

            B[1, idx2 + 1, i] = dNds1[j, i] * ds1dx
            B[2, idx2 + 2, i] = dNds2[j, i] * ds2dy
            B[3, idx2 + 3, i] = dNds3[j, i] * ds3dz

            B[4, idx2 + 1, i] = B[2, idx2 + 2, i]
            B[4, idx2 + 2, i] = B[1, idx2 + 1, i]

            B[5, idx2 + 1, i] = B[3, idx1, i]
            B[5, idx1, i] = B[1, idx2 + 1, i]

            B[6, idx2 + 2, i] = B[3, 3 * a, i]
            B[6, idx1, i] = B[2, idx2 + 2, i]
        end
    end
end
