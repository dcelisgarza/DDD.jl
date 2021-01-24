"""
```
calc_σTilde(
    x0,
    dlnParams::T1,
    matParams::T2,
    network::T3,
    idx = nothing,
) where {T1 <: DislocationParameters,T2 <: MaterialParameters,T3 <: DislocationNetwork}
```
Calculate stress on points `x0`.
"""
function calc_σTilde(
    x0,
    dlnParams::T1,
    matParams::T2,
    network::T3,
    idx = nothing,
) where {T1 <: DislocationParameters,T2 <: MaterialParameters,T3 <: DislocationNetwork}

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
    σ = zeros(6, numPoints)

    # Loop over grid points.
    @inbounds @simd for j in 1:numPoints
        x = SVector{3,elemT}(x0[1, j], x0[2, j], x0[3, j])
        # Loop over segments.
        for i in eachindex(idx)
            b = SVector{3,elemT}(bVec[1, i], bVec[2, i], bVec[3, i])
            n1 = SVector{3,elemT}(coord1[1, i], coord1[2, i], coord1[3, i])
            n2 = SVector{3,elemT}(coord2[1, i], coord2[2, i], coord2[3, i])
            t = normalize(n2 - n1)

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
            σ[5, j] += I_03yz * s_03 + I_13yz * s_13 + I_05yz * s_05 + I_15yz * s_15 + I_25yz * s_25
            σ[6, j] += I_03xz * s_03 + I_13xz * s_13 + I_05xz * s_05 + I_15xz * s_15 + I_25xz * s_25      
        end
    end
    return σ
end