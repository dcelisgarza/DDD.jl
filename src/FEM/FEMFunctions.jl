"""
```
shapeFunction(shape<:AbstractShapeFunction, x, y, z)
```
Returns the shape functions of type `typeof(shape) <: AbstractShapeFunction`. If `x,y,z` are floats returns a vector of length `N`, different shape functons have different numbers of nodes. If given vectors, returns an array of size `(N, length(x))`.
!!! note
    All coordinate vectors must be of equal length.
[`shapeFunctionDeriv`](@ref) are the 1st order derivatives of the shape functions.
"""
@inline function shapeFunction(shape::T, x, y, z) where {T <: LinearQuadrangle3D}

    # N[n](x,y,z) := shape function n.
    omx = 1 - x
    omy = 1 - y
    omz = 1 - z
    opx = 1 + x
    opy = 1 + y
    opz = 1 + z

    N =
        SVector(
            omx * omy * omz,
            opx * omy * omz,
            opx * opy * omz,
            omx * opy * omz,
            omx * omy * opz,
            opx * omy * opz,
            opx * opy * opz,
            omx * opy * opz,
        ) / 8

    return N
end
@inline function shapeFunction(
    shape::T1,
    x::T2,
    y::T2,
    z::T2,
) where {T1 <: LinearQuadrangle3D, T2 <: AbstractVector{T} where {T}}

    @assert length(x) == length(y) == length(z)
    # N[n, p](x_vec,y_vec,z_vec) := shape function n for point p.
    N = zeros(8, length(x))
    omx = 0.0
    omy = 0.0
    omz = 0.0
    opx = 0.0
    opy = 0.0
    opz = 0.0
    @inbounds @simd for i in eachindex(x)
        omx = 1 - x[i]
        omy = 1 - y[i]
        omz = 1 - z[i]
        opx = 1 + x[i]
        opy = 1 + y[i]
        opz = 1 + z[i]

        N[1, i] = omx * omy * omz
        N[2, i] = opx * omy * omz
        N[3, i] = opx * opy * omz
        N[4, i] = omx * opy * omz
        N[5, i] = omx * omy * opz
        N[6, i] = opx * omy * opz
        N[7, i] = opx * opy * opz
        N[8, i] = omx * opy * opz
    end

    N /= 8

    return N
end

"""
```
shapeFunctionDeriv(shape<:AbstractShapeFunction, x, y, z)
```
Returns the first order derivative of the shape functions, [`shapeFunction`](@ref), of type `typeof(shape) <: AbstractShapeFunction`. If `x,y,z` are floats returns a 2D array of size `(N, 3)`. If given vectors, returns a 3D array of size `(N, 3, length(x))`.
!!! note
    All coordinate vectors must be of equal length.
"""
@inline function shapeFunctionDeriv(shape::T, x, y, z) where {T <: LinearQuadrangle3D}

    # dNdS[n, x](x,y,z) := x'th derivative of shape function n.
    # dNdS[n, x](x,y,z) = dN[a, b] / dx
    omx = 1 - x
    omy = 1 - y
    omz = 1 - z
    opx = 1 + x
    opy = 1 + y
    opz = 1 + z

    dNdS =
        SMatrix{8, 3}(
            -omy * omz,
            omy * omz,
            opy * omz,
            -opy * omz,
            -omy * opz,
            omy * opz,
            opy * opz,
            -opy * opz,
            -omx * omz,
            -opx * omz,
            opx * omz,
            omx * omz,
            -omx * opz,
            -opx * opz,
            opx * opz,
            omx * opz,
            -omx * omy,
            -opx * omy,
            -opx * opy,
            -omx * opy,
            omx * omy,
            opx * omy,
            opx * opy,
            omx * opy,
        ) / 8

    return dNdS
end
@inline function shapeFunctionDeriv(
    shape::T1,
    x::T2,
    y::T2,
    z::T2,
) where {T1 <: LinearQuadrangle3D, T2 <: AbstractVector{T} where {T}}

    @assert length(x) == length(y) == length(z)
    # dNdS[n, x, p](x,y,z) := x'th derivative of shape function n for point p.
    # dNdS[n, x, p](x,y,z) = dN[a, b, p] / dx
    dNdS = zeros(8, 3, length(x))
    omx = 0.0
    omy = 0.0
    omz = 0.0
    opx = 0.0
    opy = 0.0
    opz = 0.0
    @inbounds @simd for i in eachindex(x)
        omx = 1 - x[i]
        omy = 1 - y[i]
        omz = 1 - z[i]
        opx = 1 + x[i]
        opy = 1 + y[i]
        opz = 1 + z[i]

        dNdS[1, 1, i] = -omy * omz
        dNdS[2, 1, i] = omy * omz
        dNdS[3, 1, i] = opy * omz
        dNdS[4, 1, i] = -opy * omz
        dNdS[5, 1, i] = -omy * opz
        dNdS[6, 1, i] = omy * opz
        dNdS[7, 1, i] = opy * opz
        dNdS[8, 1, i] = -opy * opz

        dNdS[1, 2, i] = -omx * omz
        dNdS[2, 2, i] = -opx * omz
        dNdS[3, 2, i] = opx * omz
        dNdS[4, 2, i] = omx * omz
        dNdS[5, 2, i] = -omx * opz
        dNdS[6, 2, i] = -opx * opz
        dNdS[7, 2, i] = opx * opz
        dNdS[8, 2, i] = omx * opz

        dNdS[1, 3, i] = -omx * omy
        dNdS[2, 3, i] = -opx * omy
        dNdS[3, 3, i] = -opx * opy
        dNdS[4, 3, i] = -omx * opy
        dNdS[5, 3, i] = omx * omy
        dNdS[6, 3, i] = opx * omy
        dNdS[7, 3, i] = opx * opy
        dNdS[8, 3, i] = omx * opy
    end

    dNdS .*= 0.125

    return dNdS
end
