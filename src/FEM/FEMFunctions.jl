function shapeFunction(
    shape::LinearQuadrangle3D,
    x::Float64,
    y::Float64,
    z::Float64,
)
    # N[n](x,y,z) := shape function n.
    N = zeros(8)

    N[1] = (1 - x) * (1 - y) * (1 - z)
    N[2] = (1 + x) * (1 - y) * (1 - z)
    N[3] = (1 + x) * (1 + y) * (1 - z)
    N[4] = (1 - x) * (1 + y) * (1 - z)
    N[5] = (1 - x) * (1 - y) * (1 + z)
    N[6] = (1 + x) * (1 - y) * (1 + z)
    N[7] = (1 + x) * (1 + y) * (1 + z)
    N[8] = (1 - x) * (1 + y) * (1 + z)

    N .*= 0.125

    return N
end

function shapeFunction(
    shape::LinearQuadrangle3D,
    x::AbstractVector{<:Float64},
    y::AbstractVector{<:Float64},
    z::AbstractVector{<:Float64},
)
    @assert length(x) == length(y) == length(z)
    # N[n, p](x_vec,y_vec,z_vec) := shape function n for point p.
    N = zeros(8, length(x))
    for i in eachindex(x)
        N[1, i] = (1 - x[i]) * (1 - y[i]) * (1 - z[i])
        N[2, i] = (1 + x[i]) * (1 - y[i]) * (1 - z[i])
        N[3, i] = (1 + x[i]) * (1 + y[i]) * (1 - z[i])
        N[4, i] = (1 - x[i]) * (1 + y[i]) * (1 - z[i])
        N[5, i] = (1 - x[i]) * (1 - y[i]) * (1 + z[i])
        N[6, i] = (1 + x[i]) * (1 - y[i]) * (1 + z[i])
        N[7, i] = (1 + x[i]) * (1 + y[i]) * (1 + z[i])
        N[8, i] = (1 - x[i]) * (1 + y[i]) * (1 + z[i])
    end

    N .*= 0.125

    return N
end

function shapeFunctionDeriv(
    shape::LinearQuadrangle3D,
    x::Float64,
    y::Float64,
    z::Float64,
)
    # dNdS[n, x](x,y,z) := x'th derivative of shape function n.
    # dNdS[n, x](x,y,z) = dN[a, b] / dx
    dNdS = zeros(8, 3)

    dNdS[1, 1] = -(1 - y) * (1 - z)
    dNdS[2, 1] = (1 - y) * (1 - z)
    dNdS[3, 1] = (1 + y) * (1 - z)
    dNdS[4, 1] = -(1 + y) * (1 - z)
    dNdS[5, 1] = -(1 - y) * (1 + z)
    dNdS[6, 1] = (1 - y) * (1 + z)
    dNdS[7, 1] = (1 + y) * (1 + z)
    dNdS[8, 1] = -(1 + y) * (1 + z)

    dNdS[1, 2] = -(1 - x) * (1 - z)
    dNdS[2, 2] = -(1 + x) * (1 - z)
    dNdS[3, 2] = (1 + x) * (1 - z)
    dNdS[4, 2] = (1 - x) * (1 - z)
    dNdS[5, 2] = -(1 - x) * (1 + z)
    dNdS[6, 2] = -(1 + x) * (1 + z)
    dNdS[7, 2] = (1 + x) * (1 + z)
    dNdS[8, 2] = (1 - x) * (1 + z)

    dNdS[1, 3] = -(1 - x) * (1 - y)
    dNdS[2, 3] = -(1 + x) * (1 - y)
    dNdS[3, 3] = -(1 + x) * (1 + y)
    dNdS[4, 3] = -(1 - x) * (1 + y)
    dNdS[5, 3] = (1 - x) * (1 - y)
    dNdS[6, 3] = (1 + x) * (1 - y)
    dNdS[7, 3] = (1 + x) * (1 + y)
    dNdS[8, 3] = (1 - x) * (1 + y)

    dNdS .*= 0.125

    return dNdS
end

function shapeFunctionDeriv(
    shape::LinearQuadrangle3D,
    x::AbstractVector{<:Float64},
    y::AbstractVector{<:Float64},
    z::AbstractVector{<:Float64},
)
    @assert length(x) == length(y) == length(z)
    # dNdS[n, x, p](x,y,z) := x'th derivative of shape function n for point p.
    # dNdS[n, x, p](x,y,z) = dN[a, b, p] / dx
    dNdS = zeros(8, 3, length(x))

    for i in eachindex(x)
        dNdS[1, 1, i] = -(1 - y[i]) * (1 - z[i])
        dNdS[2, 1, i] = (1 - y[i]) * (1 - z[i])
        dNdS[3, 1, i] = (1 + y[i]) * (1 - z[i])
        dNdS[4, 1, i] = -(1 + y[i]) * (1 - z[i])
        dNdS[5, 1, i] = -(1 - y[i]) * (1 + z[i])
        dNdS[6, 1, i] = (1 - y[i]) * (1 + z[i])
        dNdS[7, 1, i] = (1 + y[i]) * (1 + z[i])
        dNdS[8, 1, i] = -(1 + y[i]) * (1 + z[i])

        dNdS[1, 2, i] = -(1 - x[i]) * (1 - z[i])
        dNdS[2, 2, i] = -(1 + x[i]) * (1 - z[i])
        dNdS[3, 2, i] = (1 + x[i]) * (1 - z[i])
        dNdS[4, 2, i] = (1 - x[i]) * (1 - z[i])
        dNdS[5, 2, i] = -(1 - x[i]) * (1 + z[i])
        dNdS[6, 2, i] = -(1 + x[i]) * (1 + z[i])
        dNdS[7, 2, i] = (1 + x[i]) * (1 + z[i])
        dNdS[8, 2, i] = (1 - x[i]) * (1 + z[i])

        dNdS[1, 3, i] = -(1 - x[i]) * (1 - y[i])
        dNdS[2, 3, i] = -(1 + x[i]) * (1 - y[i])
        dNdS[3, 3, i] = -(1 + x[i]) * (1 + y[i])
        dNdS[4, 3, i] = -(1 - x[i]) * (1 + y[i])
        dNdS[5, 3, i] = (1 - x[i]) * (1 - y[i])
        dNdS[6, 3, i] = (1 + x[i]) * (1 - y[i])
        dNdS[7, 3, i] = (1 + x[i]) * (1 + y[i])
        dNdS[8, 3, i] = (1 - x[i]) * (1 + y[i])
    end

    dNdS .*= 0.125

    return dNdS
end
