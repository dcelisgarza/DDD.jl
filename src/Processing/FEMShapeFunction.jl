"""
```
shapeFunction(::LinearQuadrangle3D, x, y, z)
```
Computes the linear shape functions `N[1:8]` for an `(x, y, z)` point on a 3D linear quadrangle.
"""
function shapeFunction(::LinearQuadrangle3D, x, y, z)
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
"""
```
shapeFunction(::LinearQuadrangle3D, x::AbstractVector, y::AbstractVector, z::AbstractVector)
```
Computes the linear shape functions `N[1:8][p]` for `(x, y, z)` point `p` on a 3D linear quadrangle.
"""
function shapeFunction(
    ::LinearQuadrangle3D,
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractVector,
)
    numPoints = length(x)
    xType = eltype(x)
    @assert numPoints == length(y) == length(z)
    # N[n, p](x_vec,y_vec,z_vec) := shape function n for point p.
    N = Array{SVector{8, xType}}(undef, numPoints)

    for i in eachindex(x)
        omx = 1 - x[i]
        omy = 1 - y[i]
        omz = 1 - z[i]
        opx = 1 + x[i]
        opy = 1 + y[i]
        opz = 1 + z[i]

        N[i] =
            SVector{8, xType}(
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
Computes the derivatives of the linear shape functions `N[1:3, 1:8]` for an `(x, y, z)` point on a 3D linear quadrangle.

## Returns
```
dNdS[x, n](x,y,z) := x'th derivative of shape function n.
```
"""
function shapeFunctionDeriv(::LinearQuadrangle3D, x, y, z)
    # dNdS[x, n](x,y,z) := x'th derivative of shape function n.
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
"""
```
shapeFunctionDeriv(shape<:AbstractShapeFunction, x::AbstractVector, y::AbstractVector, z::AbstractVector)
```
Computes the derivatives of the linear shape functions `N[1:3, 1:8]` for an `(x, y, z)` point on a 3D linear quadrangle.

## Returns
```
dNdS[x, n, p](x,y,z) := x'th derivative of shape function n for point p.
```
"""
function shapeFunctionDeriv(
    ::LinearQuadrangle3D,
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractVector,
)
    numPoints = length(x)
    xType = eltype(x)
    @assert numPoints == length(y) == length(z)

    # dNdS[x, n, p](x,y,z) := x'th derivative of shape function n for point p.
    dNdS = Array{SMatrix{3, 8, xType}}(undef, numPoints)

    for i in eachindex(x)
        omx = 1 - x[i]
        omy = 1 - y[i]
        omz = 1 - z[i]
        opx = 1 + x[i]
        opy = 1 + y[i]
        opz = 1 + z[i]

        dNdS[i] =
            SMatrix{3, 8, xType}(
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
