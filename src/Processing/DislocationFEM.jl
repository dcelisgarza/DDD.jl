"""
```
calc_σTilde(
    x0,
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    network::DislocationNetwork,
    idx = nothing,
)
```
Computes the stress tensor, `̃σ`, induced by dislocations on points `x0`.

## Returns

```
σxx = σ[1, :]
σyy = σ[2, :]
σzz = σ[3, :]
σxy = σ[4, :]
σxz = σ[5, :]
σyz = σ[6, :]
```
"""
function calc_σTilde(
    x0,
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    network::DislocationNetwork,
    idx = nothing,
)
    # Constants
    μ4π = matParams.μ4π
    μ8π = matParams.μ8π
    μ4πν = matParams.μ4πν
    νμ4πν = matParams.νμ4πν
    a2 = dlnParams.coreRadSq
    a2μ8π = a2 * μ8π

    bVec = network.bVec
    coord = network.coord
    segIdx = network.segIdx
    elemT = eltype(network.bVec)

    if isnothing(idx)
        numSeg = network.numSeg[1]
        idx = 1:numSeg
    end

    idxBvec = @view segIdx[idx, 1]
    idxNode1 = @view segIdx[idx, 2]
    idxNode2 = @view segIdx[idx, 3]
    bVec = @view bVec[:, idxBvec]
    coord1 = @view coord[:, idxNode1]
    coord2 = @view coord[:, idxNode2]

    numPoints = div(length(x0), 3)
    numPoints == 1 ? σ = zeros(6) : σ = zeros(6, numPoints)
    
    # Loop over segments.
    @inbounds @simd for i in eachindex(idx)
        b = SVector{3,elemT}(bVec[1, i], bVec[2, i], bVec[3, i])
        n1 = SVector{3,elemT}(coord1[1, i], coord1[2, i], coord1[3, i])
        n2 = SVector{3,elemT}(coord2[1, i], coord2[2, i], coord2[3, i])
        t = normalize(n2 - n1)
        # Loop over grid points.
        for j in 1:numPoints
            x = SVector{3,elemT}(x0[1, j], x0[2, j], x0[3, j])

            R = x - n1
            Rdt = R ⋅ t
            nd = R - Rdt * t
            d2 = nd ⋅ nd
            
            s1 = -Rdt
            s2 = -(x - n2) ⋅ t
            a2_d2 = a2 + d2
            a2d2inv = 1 / a2_d2

            Ra = sqrt(a2_d2 + s1 * s1)
            Rainv = 1 / Ra
            Ra3inv = Rainv * Rainv * Rainv
            sRa3inv = s1 * Ra3inv

            s_03a = s1 * Rainv * a2d2inv
            s_13a = -Rainv
            s_05a = (2 * s_03a + sRa3inv) * a2d2inv
            s_15a = -Ra3inv
            s_25a = s_03a - sRa3inv

            Ra = sqrt(a2_d2 + s2 * s2)
            Rainv = 1 / Ra
            Ra3inv = Rainv * Rainv * Rainv
            sRa3inv = s2 * Ra3inv

            s_03b = s2 * Rainv * a2d2inv
            s_13b = -Rainv
            s_05b = (2 * s_03b + sRa3inv) * a2d2inv
            s_15b = -Ra3inv
            s_25b = s_03b - sRa3inv

            s_03 = s_03b - s_03a
            s_13 = s_13b - s_13a
            s_05 = s_05b - s_05a
            s_15 = s_15b - s_15a
            s_25 = s_25b - s_25a

            txb = t × b
            dxb = nd × b

            dxbdt = dxb ⋅ t

            dmdxx = nd[1] * nd[1]
            dmdyy = nd[2] * nd[2]
            dmdzz = nd[3] * nd[3]
            dmdxy = nd[1] * nd[2]
            dmdyz = nd[2] * nd[3]
            dmdxz = nd[1] * nd[3]

            tmtxx = t[1] * t[1]
            tmtyy = t[2] * t[2]
            tmtzz = t[3] * t[3]
            tmtxy = t[1] * t[2]
            tmtyz = t[2] * t[3]
            tmtxz = t[1] * t[3]

            tmdxx = 2 * t[1] * nd[1]
            tmdyy = 2 * t[2] * nd[2]
            tmdzz = 2 * t[3] * nd[3]
            tmdxy = t[1] * nd[2] + t[2] * nd[1]
            tmdyz = t[2] * nd[3] + t[3] * nd[2]
            tmdxz = t[1] * nd[3] + t[3] * nd[1]

            tmtxbxx = 2 * t[1] * txb[1]
            tmtxbyy = 2 * t[2] * txb[2]
            tmtxbzz = 2 * t[3] * txb[3]
            tmtxbxy = t[1] * txb[2] + t[2] * txb[1]
            tmtxbyz = t[2] * txb[3] + t[3] * txb[2]
            tmtxbxz = t[1] * txb[3] + t[3] * txb[1]

            dmtxbxx = 2 * nd[1] * txb[1]
            dmtxbyy = 2 * nd[2] * txb[2]
            dmtxbzz = 2 * nd[3] * txb[3]
            dmtxbxy = nd[1] * txb[2] + nd[2] * txb[1]
            dmtxbyz = nd[2] * txb[3] + nd[3] * txb[2]
            dmtxbxz = nd[1] * txb[3] + nd[3] * txb[1]

            tmdxbxx = 2 * t[1] * dxb[1]
            tmdxbyy = 2 * t[2] * dxb[2]
            tmdxbzz = 2 * t[3] * dxb[3]
            tmdxbxy = t[1] * dxb[2] + t[2] * dxb[1]
            tmdxbyz = t[2] * dxb[3] + t[3] * dxb[2]
            tmdxbxz = t[1] * dxb[3] + t[3] * dxb[1]

            common = μ4πν * dxbdt

            I_03xx = common + μ4πν * dmtxbxx - μ4π * tmdxbxx
            I_03yy = common + μ4πν * dmtxbyy - μ4π * tmdxbyy
            I_03zz = common + μ4πν * dmtxbzz - μ4π * tmdxbzz
            I_03xy = μ4πν * dmtxbxy - μ4π * tmdxbxy
            I_03yz = μ4πν * dmtxbyz - μ4π * tmdxbyz
            I_03xz = μ4πν * dmtxbxz - μ4π * tmdxbxz

            I_13xx = -νμ4πν * tmtxbxx
            I_13yy = -νμ4πν * tmtxbyy
            I_13zz = -νμ4πν * tmtxbzz
            I_13xy = -νμ4πν * tmtxbxy
            I_13yz = -νμ4πν * tmtxbyz
            I_13xz = -νμ4πν * tmtxbxz

            I_05xx = common * (a2 + dmdxx) - a2μ8π * tmdxbxx
            I_05yy = common * (a2 + dmdyy) - a2μ8π * tmdxbyy
            I_05zz = common * (a2 + dmdzz) - a2μ8π * tmdxbzz
            I_05xy = common * dmdxy - a2μ8π * tmdxbxy
            I_05yz = common * dmdyz - a2μ8π * tmdxbyz
            I_05xz = common * dmdxz - a2μ8π * tmdxbxz

            I_15xx = a2μ8π * tmtxbxx - common * tmdxx
            I_15yy = a2μ8π * tmtxbyy - common * tmdyy
            I_15zz = a2μ8π * tmtxbzz - common * tmdzz
            I_15xy = a2μ8π * tmtxbxy - common * tmdxy
            I_15yz = a2μ8π * tmtxbyz - common * tmdyz
            I_15xz = a2μ8π * tmtxbxz - common * tmdxz

            I_25xx = common * tmtxx
            I_25yy = common * tmtyy
            I_25zz = common * tmtzz
            I_25xy = common * tmtxy
            I_25yz = common * tmtyz
            I_25xz = common * tmtxz

            σ[1, j] += I_03xx * s_03 + I_13xx * s_13 + I_05xx * s_05 + I_15xx * s_15 + I_25xx * s_25
            σ[2, j] += I_03yy * s_03 + I_13yy * s_13 + I_05yy * s_05 +  I_15yy * s_15 + I_25yy * s_25
            σ[3, j] += I_03zz * s_03 + I_13zz * s_13 + I_05zz * s_05 + I_15zz * s_15 + I_25zz * s_25
            σ[4, j] += I_03xy * s_03 + I_13xy * s_13 + I_05xy * s_05 + I_15xy * s_15 + I_25xy * s_25
            σ[5, j] += I_03xz * s_03 + I_13xz * s_13 + I_05xz * s_05 + I_15xz * s_15 + I_25xz * s_25      
            σ[6, j] += I_03yz * s_03 + I_13yz * s_13 + I_05yz * s_05 + I_15yz * s_15 + I_25yz * s_25
        end
    end
    return σ
end

"""
```
calc_σTilde!(
    σ,
    x0,
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    network::DislocationNetwork,
idx = nothing,
)
```
In-place computation of the stress tensor,`̃σ`, induced by dislocations on points `x0`.

## Returns

```
σxx = σ[1, :]
σyy = σ[2, :]
σzz = σ[3, :]
σxy = σ[4, :]
σxz = σ[5, :]
σyz = σ[6, :]
```
"""
function calc_σTilde!(
    σ,
    x0,
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    network::DislocationNetwork,
    idx = nothing,
)
    @assert ndims(σ) == ndims(x0) && size(σ, 2) == size(x0, 2) && size(σ, 1) == 6
    # Constants
    μ4π = matParams.μ4π
    μ8π = matParams.μ8π
    μ4πν = matParams.μ4πν
    νμ4πν = matParams.νμ4πν
    a2 = dlnParams.coreRadSq
    a2μ8π = a2 * μ8π
    
    bVec = network.bVec
    coord = network.coord
    segIdx = network.segIdx
    elemT = eltype(network.bVec)
    σ .= zero(elemT)

    if isnothing(idx)
        numSeg = network.numSeg[1]
    idx = 1:numSeg
    end

    idxBvec = @view segIdx[idx, 1]
    idxNode1 = @view segIdx[idx, 2]
    idxNode2 = @view segIdx[idx, 3]
    bVec = @view bVec[:, idxBvec]
    coord1 = @view coord[:, idxNode1]
    coord2 = @view coord[:, idxNode2]

    numPoints = div(length(x0), 3)

    # Loop over segments.
    @inbounds @simd for i in eachindex(idx)
        b = SVector{3,elemT}(bVec[1, i], bVec[2, i], bVec[3, i])
        n1 = SVector{3,elemT}(coord1[1, i], coord1[2, i], coord1[3, i])
        n2 = SVector{3,elemT}(coord2[1, i], coord2[2, i], coord2[3, i])
        t = normalize(n2 - n1)
        # Loop over grid points.
        for j in 1:numPoints
            x = SVector{3,elemT}(x0[1, j], x0[2, j], x0[3, j])

            R = x - n1
            Rdt = R ⋅ t
            nd = R - Rdt * t
            d2 = nd ⋅ nd
            
            s1 = -Rdt
            s2 = -(x - n2) ⋅ t
            a2_d2 = a2 + d2
            a2d2inv = 1 / a2_d2

            Ra = sqrt(a2_d2 + s1 * s1)
            Rainv = 1 / Ra
            Ra3inv = Rainv * Rainv * Rainv
            sRa3inv = s1 * Ra3inv

            s_03a = s1 * Rainv * a2d2inv
            s_13a = -Rainv
            s_05a = (2 * s_03a + sRa3inv) * a2d2inv
            s_15a = -Ra3inv
            s_25a = s_03a - sRa3inv

            Ra = sqrt(a2_d2 + s2 * s2)
            Rainv = 1 / Ra
            Ra3inv = Rainv * Rainv * Rainv
            sRa3inv = s2 * Ra3inv

            s_03b = s2 * Rainv * a2d2inv
            s_13b = -Rainv
            s_05b = (2 * s_03b + sRa3inv) * a2d2inv
            s_15b = -Ra3inv
            s_25b = s_03b - sRa3inv

            s_03 = s_03b - s_03a
            s_13 = s_13b - s_13a
            s_05 = s_05b - s_05a
            s_15 = s_15b - s_15a
            s_25 = s_25b - s_25a

            txb = t × b
            dxb = nd × b

            dxbdt = dxb ⋅ t

            dmdxx = nd[1] * nd[1]
            dmdyy = nd[2] * nd[2]
            dmdzz = nd[3] * nd[3]
            dmdxy = nd[1] * nd[2]
            dmdyz = nd[2] * nd[3]
            dmdxz = nd[1] * nd[3]

            tmtxx = t[1] * t[1]
            tmtyy = t[2] * t[2]
            tmtzz = t[3] * t[3]
            tmtxy = t[1] * t[2]
            tmtyz = t[2] * t[3]
            tmtxz = t[1] * t[3]

            tmdxx = 2 * t[1] * nd[1]
            tmdyy = 2 * t[2] * nd[2]
            tmdzz = 2 * t[3] * nd[3]
            tmdxy = t[1] * nd[2] + t[2] * nd[1]
            tmdyz = t[2] * nd[3] + t[3] * nd[2]
            tmdxz = t[1] * nd[3] + t[3] * nd[1]

            tmtxbxx = 2 * t[1] * txb[1]
            tmtxbyy = 2 * t[2] * txb[2]
            tmtxbzz = 2 * t[3] * txb[3]
            tmtxbxy = t[1] * txb[2] + t[2] * txb[1]
            tmtxbyz = t[2] * txb[3] + t[3] * txb[2]
            tmtxbxz = t[1] * txb[3] + t[3] * txb[1]

            dmtxbxx = 2 * nd[1] * txb[1]
            dmtxbyy = 2 * nd[2] * txb[2]
            dmtxbzz = 2 * nd[3] * txb[3]
            dmtxbxy = nd[1] * txb[2] + nd[2] * txb[1]
            dmtxbyz = nd[2] * txb[3] + nd[3] * txb[2]
            dmtxbxz = nd[1] * txb[3] + nd[3] * txb[1]

            tmdxbxx = 2 * t[1] * dxb[1]
            tmdxbyy = 2 * t[2] * dxb[2]
            tmdxbzz = 2 * t[3] * dxb[3]
            tmdxbxy = t[1] * dxb[2] + t[2] * dxb[1]
            tmdxbyz = t[2] * dxb[3] + t[3] * dxb[2]
            tmdxbxz = t[1] * dxb[3] + t[3] * dxb[1]

            common = μ4πν * dxbdt

            I_03xx = common + μ4πν * dmtxbxx - μ4π * tmdxbxx
            I_03yy = common + μ4πν * dmtxbyy - μ4π * tmdxbyy
            I_03zz = common + μ4πν * dmtxbzz - μ4π * tmdxbzz
            I_03xy = μ4πν * dmtxbxy - μ4π * tmdxbxy
            I_03yz = μ4πν * dmtxbyz - μ4π * tmdxbyz
            I_03xz = μ4πν * dmtxbxz - μ4π * tmdxbxz

            I_13xx = -νμ4πν * tmtxbxx
            I_13yy = -νμ4πν * tmtxbyy
            I_13zz = -νμ4πν * tmtxbzz
            I_13xy = -νμ4πν * tmtxbxy
            I_13yz = -νμ4πν * tmtxbyz
            I_13xz = -νμ4πν * tmtxbxz

            I_05xx = common * (a2 + dmdxx) - a2μ8π * tmdxbxx
            I_05yy = common * (a2 + dmdyy) - a2μ8π * tmdxbyy
            I_05zz = common * (a2 + dmdzz) - a2μ8π * tmdxbzz
            I_05xy = common * dmdxy - a2μ8π * tmdxbxy
            I_05yz = common * dmdyz - a2μ8π * tmdxbyz
            I_05xz = common * dmdxz - a2μ8π * tmdxbxz

            I_15xx = a2μ8π * tmtxbxx - common * tmdxx
            I_15yy = a2μ8π * tmtxbyy - common * tmdyy
            I_15zz = a2μ8π * tmtxbzz - common * tmdzz
            I_15xy = a2μ8π * tmtxbxy - common * tmdxy
            I_15yz = a2μ8π * tmtxbyz - common * tmdyz
            I_15xz = a2μ8π * tmtxbxz - common * tmdxz

            I_25xx = common * tmtxx
            I_25yy = common * tmtyy
            I_25zz = common * tmtzz
            I_25xy = common * tmtxy
            I_25yz = common * tmtyz
            I_25xz = common * tmtxz

            σ[1, j] += I_03xx * s_03 + I_13xx * s_13 + I_05xx * s_05 + I_15xx * s_15 + I_25xx * s_25
            σ[2, j] += I_03yy * s_03 + I_13yy * s_13 + I_05yy * s_05 +  I_15yy * s_15 + I_25yy * s_25
            σ[3, j] += I_03zz * s_03 + I_13zz * s_13 + I_05zz * s_05 + I_15zz * s_15 + I_25zz * s_25
            σ[4, j] += I_03xy * s_03 + I_13xy * s_13 + I_05xy * s_05 + I_15xy * s_15 + I_25xy * s_25
            σ[5, j] += I_03xz * s_03 + I_13xz * s_13 + I_05xz * s_05 + I_15xz * s_15 + I_25xz * s_25      
            σ[6, j] += I_03yz * s_03 + I_13yz * s_13 + I_05yz * s_05 + I_15yz * s_15 + I_25yz * s_25
        end
end
    return σ
end

"""
Calculates displacements from dislocations.
Bruce Bromage, bruce.bromage@materials.ox.ac.uk
Github: @brucebromage

[Calculating dislocation displacements on the surface of a volume](https://iopscience.iop.org/article/10.1088/1361-651X/aae404/pdf)
B Bromage and E Tarleton
Published 29 October 2018 • © 2018 IOP Publishing Ltd
Modelling and Simulation in Materials Science and Engineering, Volume 26,
Number 8

```latex
@article{bromage2018calculating,
  title={Calculating dislocation displacements on the surface of a volume},
  author={Bromage, B and Tarleton, E},
  journal={Modelling and Simulation in Materials Science and Engineering},
  volume={26},
  number={8},
pages={085007},
  year={2018},
  publisher={IOP Publishing}
}
```
"""
function calc_uTilde!(
    forceDisplacement::ForceDisplacement,
    mesh::AbstractMesh,
    boundaries::Boundaries,
    matParams::MaterialParameters,
    network::DislocationNetwork
)

    uTilde = forceDisplacement.uTilde
    uTilde .= zero(eltype(uTilde))

    C = mean(mesh.vertices.vertices)
    faces = mesh.faces
    faceNorm = mesh.faceNorm
    faceMidPt = mesh.faceMidPt
    coordFE = mesh.coord

    # FE nodes with applied displacements.
    uNodes = boundaries.uGammaDln
    uDofs = boundaries.uDofsDln

    # Coordinates of the FE nodes with applied displacements.
    uCoord = @view coordFE[:, uNodes]

    numNode = network.numNode[1]
    numSeg = network.numSeg[1]
    label = network.label
    links = network.links
    bVec = network.bVec
    coordDln = network.coord
    elemT = eltype(coordDln)
    
    tmpArr = zeros(elemT, 3)

    if typeof(mesh) <: AbstractRegularCuboidMesh
        # We only need to check half the normals for a cuboid.
    normals = faceNorm[:, 1:2:5]
    else
        @warn "calc_uTilde! could use fewer normals for mesh of type $(typeof(mesh)), see $(@__LINE__)"
        normals = faceNorm
    end
    
    surfNodeCons = spzeros(Int, numNode)
    @inbounds for i in 1:numSeg
        link1 = links[1, i]
        link2 = links[2, i]
        # Find external nodes connected to a surface node, stores it and skips the virtual segment.
        if (label[link1] == srfMobDln || label[link1] == srfFixDln) && label[link2] == extDln
            surfNodeCons[link2] = link1
            continue
        elseif (label[link2] == srfMobDln || label[link2] == srfFixDln) && label[link1] == extDln
            surfNodeCons[link1] = link2
            continue
        end

        intSeg::Bool = true

        # Coords of first node of the segment.
        A = SVector{3,elemT}(coordDln[1, link1], coordDln[2, link1], coordDln[3, link1])
        # Coords of second node of the segment.
        B = SVector{3,elemT}(coordDln[1, link2], coordDln[2, link2], coordDln[3, link2])
        # Segment's burgers vector.
        b = SVector{3,elemT}(coordDln[1, i], coordDln[2, i], coordDln[3, i])

        # Find external segments.
        if label[link1] == extDln || label[link2] == extDln
            intSeg = false
            # Find the nearest intersect with the mesh to find the point on the mesh's surface the node was projected out from.
            distMinA::elemT = Inf
            distMinB::elemT = Inf
            distMinC::elemT = Inf
            distMinTmpA::elemT = 0
            distMinTmpB::elemT = 0
            intersectA = SVector{3,elemT}(Inf, Inf, Inf)
            intersectB = SVector{3,elemT}(Inf, Inf, Inf)
            for j in 1:size(normals, 2)
                distMinTmpA, intersectTmpA, missing = findIntersectVolume(mesh, normals[:, j], A, tmpArr)
                distMinTmpB, intersectTmpB, missing = findIntersectVolume(mesh, normals[:, j], B, tmpArr)
                if distMinTmpA < distMinA
                    distMinA = distMinTmpA
                    intersectA = intersectTmpA
                end
                if distMinTmpB < distMinB
                    distMinB = distMinTmpB
                    intersectB = intersectTmpB
                end
            end

            # In case the projection was out of a corner or edge, there will be no intersect so we have to find it iteratively using the planes of the surfaces.
            # We find the intersects with the planes that get us closer to the centre of the domain.
            if isinf(distMinA)
                Atmp = A
                # Loop through the faces of the solid.
                for i in 1:size(faceNorm, 2)
                    # Calculate the distance from the node to the plane of face i.
                    distAtmp = ((Atmp + faceMidPt[:, i]) ⋅ faceNorm[:, i])^2
                    # Find the intersect of the node with the plane of face i (nodes are external so we have to flip the sign of the face normal)
                    intersectATmp = Atmp - faceNorm[:, i] * sqrt(distAtmp)
                    # Calculate the distance from the domain's centroid.
                    distTmpC = norm(intersectATmp - C)
                    # If the distance from the intersect and the centroid is smaller than the minimum, we move the coordinate to the intersect. This ensures we get closer to the centroid after every iteration, eventually we will intersect the domain.
                    if distTmpC <= distMinC
                        distMinC = distTmpC
                        Atmp = intersectATmp
                    end
                end
                intersectA = Atmp

            end

            if isinf(distMinB)
                Btmp = B
                for i in 1:size(faceNorm, 2)
                    distBtmp = ((Btmp + faceMidPt[:, i]) ⋅ faceNorm[:, i])^2
                    intersectBTmp = Btmp - faceNorm[:, i] * sqrt(distBtmp)
                    distTmpC = norm(intersectBTmp - C)
                    if distTmpC <= distMinC
                        distMinC = distTmpC
                        Btmp = intersectBTmp
                    end
                end
                intersectB = Btmp

            end

            Aprime = intersectA # Point where A was projected from.
            Bprime = intersectB # Point where B was projected from.

            # Check if A or B are connected to a surface point. If any of them are, change the intersect to that point.
            surfConA = surfNodeCons[link1]
            surfConB = surfNodeCons[link2]

            surfConA != 0 ? Aprime = SVector{3,elemT}(coordDln[1, surfConA], coordDln[2, surfConA], coordDln[3, surfConA]) : nothing
            surfConB != 0 ? Bprime = SVector{3,elemT}(coordDln[1, surfConB], coordDln[2, surfConB], coordDln[3, surfConB]) : nothing
            
            # Closure point for external loop.
            Cprime = (A + Bprime) * 0.5
        end

        # Use Barnett triangles to calculate displacements using a closure point.
        if intSeg
            calcDisplacementDislocationTriangle!(uTilde, uDofs, matParams, A, B, C, b, uCoord)
        else
            # For external loops close internal loop with a normal closure point and calculate the displacement for the external loop with an auxiliary closure point.
            calcDisplacementDislocationTriangle!(uTilde, uDofs, Aprime, Bprime, C, b, uCoord)
            calcDisplacementDislocationTriangle!(uTilde, uDofs, Aprime, A, Cprime, b, uCoord)
            calcDisplacementDislocationTriangle!(uTilde, uDofs, A, B, Cprime, b, uCoord)
            calcDisplacementDislocationTriangle!(uTilde, uDofs, B, Bprime, Cprime, b, uCoord)
            calcDisplacementDislocationTriangle!(uTilde, uDofs, Bprime, Aprime, Cprime, b, uCoord) 
        end
    end

    return nothing
end
function calc_uTilde(
    forceDisplacement::ForceDisplacement,
    mesh::AbstractMesh,
    boundaries::Boundaries,
    matParams::MaterialParameters,
    network::DislocationNetwork
)

    C = mean(mesh.vertices.vertices)
    scale = mesh.scale
    faces = mesh.faces
    faceNorm = mesh.faceNorm
    faceMidPt = mesh.faceMidPt
    coordFE = mesh.coord

    # FE nodes with applied displacements.
    uNodes = boundaries.uGammaDln
    uDofs = 1:length(boundaries.uDofsDln)
    uTilde = zeros(length(uDofs))


    # Coordinates of the FE nodes with applied displacements.
    uCoord = @view coordFE[:, uNodes]

    numNode = network.numNode[1]
    numSeg = network.numSeg[1]
    label = network.label
    links = network.links
    bVec = network.bVec
    coordDln = network.coord
    elemT = eltype(coordDln)
    
    tmpArr = zeros(elemT, 3)

    if typeof(mesh) <: AbstractRegularCuboidMesh
        # We only need to check half the normals for a cuboid.
    normals = faceNorm[:, 1:2:5]
    else
        @warn "calc_uTilde! could use fewer normals for mesh of type $(typeof(mesh)), see $(@__LINE__)"
        normals = faceNorm
    end
    
    surfNodeCons = spzeros(Int, numNode)
    @inbounds for i in 1:numSeg
        link1 = links[1, i]
        link2 = links[2, i]
        # Find external nodes connected to a surface node, stores it and skips the virtual segment.
        if (label[link1] == srfMobDln || label[link1] == srfFixDln) && label[link2] == extDln
            surfNodeCons[link2] = link1
            continue
        elseif (label[link2] == srfMobDln || label[link2] == srfFixDln) && label[link1] == extDln
            surfNodeCons[link1] = link2
            continue
        end

        intSeg::Bool = true

        # Coords of first node of the segment.
        A = SVector{3,elemT}(coordDln[1, link1], coordDln[2, link1], coordDln[3, link1])
        # Coords of second node of the segment.
        B = SVector{3,elemT}(coordDln[1, link2], coordDln[2, link2], coordDln[3, link2])
        # Segment's burgers vector.
        b = SVector{3,elemT}(coordDln[1, i], coordDln[2, i], coordDln[3, i])

        # Find external segments.
        faceA = 0
        faceB = 0
        if label[link1] == extDln || label[link2] == extDln
            intSeg = false

            # Find the nearest intersect with the mesh to find the point on the mesh's surface the node was projected out from.
            distMinA::elemT = Inf
            distMinB::elemT = Inf
            distMinC::elemT = Inf

            distMinTmpA::elemT = 0
            distMinTmpB::elemT = 0

            intersectA = SVector{3,elemT}(Inf, Inf, Inf)
            intersectB = SVector{3,elemT}(Inf, Inf, Inf)
            for j in 1:size(normals, 2)
                distMinTmpA, intersectTmpA, faceA = findIntersectVolume(mesh, normals[:, j], A, tmpArr)
                distMinTmpB, intersectTmpB, faceB = findIntersectVolume(mesh, normals[:, j], B, tmpArr)
                if distMinTmpA < distMinA
                    distMinA = distMinTmpA
                    intersectA = intersectTmpA
                end
                if distMinTmpB < distMinB
                    distMinB = distMinTmpB
                    intersectB = intersectTmpB
                end
            end

            # In case the projection was out of a corner or edge, there will be no intersect so we have to find it iteratively using the planes of the surfaces.
            # We find the intersects with the planes that get us closer to the centre of the domain.
            if isinf(distMinA)
                Atmp = A
                # Loop through the faces of the solid.
                for i in 1:size(faceNorm, 2)
                    # Calculate the distance from the node to the plane of face i.
                    distAtmp = ((Atmp + faceMidPt[:, i]) ⋅ faceNorm[:, i])^2
                    # Find the intersect of the node with the plane of face i (nodes are external so we have to flip the sign of the face normal)
                    intersectATmp = Atmp - faceNorm[:, i] * sqrt(distAtmp)
                    # Calculate the distance from the domain's centroid.
                    distTmpC = norm(intersectATmp - C)
                    # If the distance from the intersect and the centroid is smaller than the minimum, we move the coordinate to the intersect. This ensures we get closer to the centroid after every iteration, eventually we will intersect the domain.
                    if distTmpC <= distMinC
                        distMinC = distTmpC
                        Atmp = intersectATmp
                    end
                end
                intersectA = Atmp

            end

            if isinf(distMinB)
                Btmp = B
                for i in 1:size(faceNorm, 2)
                    distBtmp = ((Btmp + faceMidPt[:, i]) ⋅ faceNorm[:, i])^2
                    intersectBTmp = Btmp - faceNorm[:, i] * sqrt(distBtmp)
                    distTmpC = norm(intersectBTmp - C)
                    if distTmpC <= distMinC
                        distMinC = distTmpC
                        Btmp = intersectBTmp
                    end
                end
                intersectB = Btmp

            end

            Aprime = intersectA # Point where A was projected from.
            Bprime = intersectB # Point where B was projected from.

            # Check if A or B are connected to a surface point. If any of them are, change the intersect to that point.
            surfConA = surfNodeCons[link1]
            surfConB = surfNodeCons[link2]

            surfConA != 0 ? Aprime = SVector{3,elemT}(coordDln[1, surfConA], coordDln[2, surfConA], coordDln[3, surfConA]) : nothing
            surfConB != 0 ? Bprime = SVector{3,elemT}(coordDln[1, surfConB], coordDln[2, surfConB], coordDln[3, surfConB]) : nothing
            
            # Closure point for external loop.
            Cprime = (A + Bprime) * 0.5
        end

        # Use Barnett triangles to calculate displacements using a closure point.
        if intSeg
            calcDisplacementDislocationTriangle!(uTilde, uDofs, matParams, A, B, C, b, uCoord)
        else
            # For external loops close internal loop with a normal closure point and calculate the displacement for the external loop with an auxiliary closure point.
            calcDisplacementDislocationTriangle!(uTilde, uDofs, Aprime, Bprime, C, b, uCoord)
            calcDisplacementDislocationTriangle!(uTilde, uDofs, Aprime, A, Cprime, b, uCoord)
            calcDisplacementDislocationTriangle!(uTilde, uDofs, A, B, Cprime, b, uCoord)
            calcDisplacementDislocationTriangle!(uTilde, uDofs, B, Bprime, Cprime, b, uCoord)
            calcDisplacementDislocationTriangle!(uTilde, uDofs, Bprime, Aprime, Cprime, b, uCoord)
        end
    end
    return uTilde
end

"""
```
calcDisplacementDislocationTriangle!(uTilde, uDofs, matParams, A, B, C, b, P)
```
Written by F.Hofmann 9/7/09
Routine to compute the displacement field for a triangular dislocation loop ABC at point P.

Modified by F.Hofmann 5/11/18

# Inputs

- `A, B, C`: three column vectors defining the nodes of the dislocation loop.
- `P`: column vector with coordinates of the point at which the displacement is evaluated in dimension 1. Then in dimension 2 this is a list of points at which to evaluate the field.
- `b`: column vector with 3 burgers vector components.

Translated by Daniel Celis Garza.
"""
function calcDisplacementDislocationTriangle!(uTilde, uDofs, matParams, A, B, C, b, P)
    elemT = eltype(A)

    omνInv8π = matParams.omνInv8π
    om2νomνInv8π = matParams.om2νomνInv8π

    # Segment tangent vectors t.
    AB = B - A
    BC = C - B
    CA = A - C
    missing, AB = safeNorm(AB)
    missing, BC = safeNorm(BC)
    missing, CA = safeNorm(CA)

    @inbounds @simd for i in 1:size(P, 2)
        # Local R vectors.
        Pi = SVector{3,elemT}(P[1, i], P[2, i], P[3, i])
        
        RA = A - Pi
        RB = B - Pi
        RC = C - Pi

        RAn, RA = safeNorm(RA)
        RBn, RB = safeNorm(RB)
        RCn, RC = safeNorm(RC)  

        # Compute fs
        fAB = (b × AB) * log((RBn / RAn) * (1 + RB ⋅ AB) / (1 + RA ⋅ AB))
        fBC = (b × BC) * log((RCn / RBn) * (1 + RC ⋅ BC) / (1 + RB ⋅ BC))
        fCA = (b × CA) * log((RAn / RCn) * (1 + RA ⋅ CA) / (1 + RC ⋅ CA))

        # Compute gs
        RAdRB = clamp(RA ⋅ RB, -1, 1)
        RBdRC = clamp(RB ⋅ RC, -1, 1)
        RCdRA = clamp(RC ⋅ RA, -1, 1)
        gAB = (b ⋅ (RA × RB)) * (RA + RB) / (1 + RAdRB)
        gBC = (b ⋅ (RB × RC)) * (RB + RC) / (1 + RBdRC)
        gCA = (b ⋅ (RC × RA)) * (RC + RA) / (1 + RCdRA)

        # Compute solid angle
        θa = acos(RAdRB)
        θb = acos(RBdRC)
        θc = acos(RCdRA)
        θ = (θa + θb + θc) / 2

        factor = max(tan(θ / 2) * tan((θ - θa) / 2) * tan((θ - θb) / 2)  * tan((θ - θc) / 2), 0)
        ω = 4 * atan(sqrt(factor))
        ω = -sign(RA ⋅ (CA × AB)) * ω

        u = -b * ω / (4 * π) - om2νomνInv8π * (fAB + fBC + fCA) + omνInv8π * (gAB + gBC + gCA)
        any(x -> isnan(x), u) ? println(i, "b = ", b, " ") : nothing
        idx = 3 * i
        uTilde[uDofs[idx - 2]] += u[1]
        uTilde[uDofs[idx - 1]] += u[2]
        uTilde[uDofs[idx]] += u[3]
    end

    return nothing
end
function calcDisplacementDislocationTriangle!(uTilde, uDofs, A, B, C, b, P)
    elemT = eltype(A)

    # Segment tangent vectors t.
    AB = B - A
    CA = A - C
    missing, AB = safeNorm(AB)
    missing, CA = safeNorm(CA)

    for i in 1:size(P, 2)
        # Local R vectors.
        Pi = SVector{3,elemT}(P[1, i], P[2, i], P[3, i])
        
        RA = A - Pi
        RB = B - Pi
        RC = C - Pi

        missing, RA = safeNorm(RA)
        missing, RB = safeNorm(RB)
        missing, RC = safeNorm(RC)  

        # Compute gs
        RAdRB = clamp(RA ⋅ RB, -1, 1)
        RBdRC = clamp(RB ⋅ RC, -1, 1)
        RCdRA = clamp(RC ⋅ RA, -1, 1)

        # Compute solid angle
        θa = acos(RAdRB)
        θb = acos(RBdRC)
        θc = acos(RCdRA)
        θ = (θa + θb + θc) / 2

        factor = max(tan(θ / 2) * tan((θ - θa) / 2) * tan((θ - θb) / 2)  * tan((θ - θc) / 2), 0)
        ω = 4 * atan(sqrt(factor))
        ω = -sign(RA ⋅ (CA × AB)) * ω

        u = -b * ω / (4 * π)

        idx = 3 * i
        uTilde[uDofs[idx - 2]] += u[1]
        uTilde[uDofs[idx - 1]] += u[2]
        uTilde[uDofs[idx]] += u[3]
end

    return nothing
end

