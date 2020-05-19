"""
```
shapeFunction(shape<:AbstractShapeFunction, x, y, z)
```
Returns the shape functions of type `typeof(shape) <: AbstractShapeFunction`. If `x,y,z` are floats returns a vector of length `N`, different shape functons have different numbers of nodes. If given vectors, returns an array of size `(N, length(x))`.
!!! note
    All coordinate vectors must be of equal length.
[`shapeFunctionDeriv`](@ref) are the 1st order derivatives of the shape functions.
"""
@inline function shapeFunction(shape::LinearQuadrangle3D, x, y, z)
    # N[n](x,y,z) := shape function n.
    N = zeros(8)
    omx = 1 - x
    omy = 1 - y
    omz = 1 - z
    opx = 1 + x
    opy = 1 + y
    opz = 1 + z

    N[1] = omx * omy * omz
    N[2] = opx * omy * omz
    N[3] = opx * opy * omz
    N[4] = omx * opy * omz
    N[5] = omx * omy * opz
    N[6] = opx * omy * opz
    N[7] = opx * opy * opz
    N[8] = omx * opy * opz

    N .*= 0.125

    return N
end
@inline function shapeFunction(
    shape::LinearQuadrangle3D,
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
) where {T}
    @assert length(x) == length(y) == length(z)
    # N[n, p](x_vec,y_vec,z_vec) := shape function n for point p.
    N = zeros(8, length(x))
    omx::T = 0
    omy::T = 0
    omz::T = 0
    opx::T = 0
    opy::T = 0
    opz::T = 0
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

    N .*= 0.125

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
@inline function shapeFunctionDeriv(shape::LinearQuadrangle3D, x, y, z)
    # dNdS[n, x](x,y,z) := x'th derivative of shape function n.
    # dNdS[n, x](x,y,z) = dN[a, b] / dx
    dNdS = zeros(8, 3)
    omx = 1 - x
    omy = 1 - y
    omz = 1 - z
    opx = 1 + x
    opy = 1 + y
    opz = 1 + z

    dNdS[1, 1] = -omy * omz
    dNdS[2, 1] = omy * omz
    dNdS[3, 1] = opy * omz
    dNdS[4, 1] = -opy * omz
    dNdS[5, 1] = -omy * opz
    dNdS[6, 1] = omy * opz
    dNdS[7, 1] = opy * opz
    dNdS[8, 1] = -opy * opz

    dNdS[1, 2] = -omx * omz
    dNdS[2, 2] = -opx * omz
    dNdS[3, 2] = opx * omz
    dNdS[4, 2] = omx * omz
    dNdS[5, 2] = -omx * opz
    dNdS[6, 2] = -opx * opz
    dNdS[7, 2] = opx * opz
    dNdS[8, 2] = omx * opz

    dNdS[1, 3] = -omx * omy
    dNdS[2, 3] = -opx * omy
    dNdS[3, 3] = -opx * opy
    dNdS[4, 3] = -omx * opy
    dNdS[5, 3] = omx * omy
    dNdS[6, 3] = opx * omy
    dNdS[7, 3] = opx * opy
    dNdS[8, 3] = omx * opy

    dNdS .*= 0.125

    return dNdS
end
@inline function shapeFunctionDeriv(
    shape::LinearQuadrangle3D,
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
) where {T}
    @assert length(x) == length(y) == length(z)
    # dNdS[n, x, p](x,y,z) := x'th derivative of shape function n for point p.
    # dNdS[n, x, p](x,y,z) = dN[a, b, p] / dx
    dNdS = zeros(8, 3, length(x))
    omx::T = 0
    omy::T = 0
    omz::T = 0
    opx::T = 0
    opy::T = 0
    opz::T = 0
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
