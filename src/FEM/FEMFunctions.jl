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
        SVector{8, typeof(x)}(
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

    numPoints = length(x)
    xType = eltype(x)
    @assert numPoints == length(y) == length(z)
    # N[n, p](x_vec,y_vec,z_vec) := shape function n for point p.
    N = Array{SVector{8, xType}}(undef, numPoints)

    @inbounds @simd for i in eachindex(x)
        omx = 1 - x[i]
        omy = 1 - y[i]
        omz = 1 - z[i]
        opx = 1 + x[i]
        opy = 1 + y[i]
        opz = 1 + z[i]

        N[i] = SVector{8, xType}(
            omx * omy * omz,
            opx * omy * omz,
            opx * opy * omz,
            omx * opy * omz,
            omx * omy * opz,
            opx * omy * opz,
            opx * opy * opz,
            omx * opy * opz,
        ) / 8
    end

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
        SMatrix{3, 8, typeof(x)}(
            -omy * omz,
            -omx * omz,
            -omx * omy,
            omy * omz,
            -opx * omz,
            -opx * omy,
            opy * omz,
            opx * omz,
            -opx * opy,
            -opy * omz,
            omx * omz,
            -omx * opy,
            -omy * opz,
            -omx * opz,
            omx * omy,
            omy * opz,
            -opx * opz,
            opx * omy,
            opy * opz,
            opx * opz,
            opx * opy,
            -opy * opz,
            omx * opz,
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

    numPoints = length(x)
    xType = eltype(x)
    @assert numPoints == length(y) == length(z)
    # dNdS[n, x, p](x,y,z) := x'th derivative of shape function n for point p.
    # dNdS[n, x, p](x,y,z) = dN[a, b, p] / dx
    dNdS = Array{SMatrix{3, 8, xType}}(undef, numPoints)#zeros(8, 3, length(x))

    @inbounds @simd for i in eachindex(x)
        omx = 1 - x[i]
        omy = 1 - y[i]
        omz = 1 - z[i]
        opx = 1 + x[i]
        opy = 1 + y[i]
        opz = 1 + z[i]

        dNdS[i] = SMatrix{3, 8, xType}(
            -omy * omz,
            -omx * omz,
            -omx * omy,
            omy * omz,
            -opx * omz,
            -opx * omy,
            opy * omz,
            opx * omz,
            -opx * opy,
            -opy * omz,
            omx * omz,
            -omx * opy,
            -omy * opz,
            -omx * opz,
            omx * omy,
            omy * opz,
            -opx * opz,
            opx * omy,
            opy * opz,
            opx * opz,
            opx * opy,
            -opy * opz,
            omx * opz,
            omx * opy,
        ) / 8
    end

    return dNdS
end
