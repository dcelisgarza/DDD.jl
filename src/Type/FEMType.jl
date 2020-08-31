"""
Mesh types.
```
abstract type AbstractMesh end
```
"""
abstract type AbstractMesh end

"""
Shape function types.
```
abstract type AbstractShapeFunction end
abstract type AbstractShapeFunction3D <: AbstractShapeFunction end
abstract type AbstractShapeFunction2D <: AbstractShapeFunction end
struct LinearQuadrangle3D <:AbstractShapeFunction3D end
struct LinearQuadrangle2D <:AbstractShapeFunction2D end
```
"""
abstract type AbstractShapeFunction end
abstract type AbstractShapeFunction3D <: AbstractShapeFunction end
abstract type AbstractShapeFunction2D <: AbstractShapeFunction end
struct LinearQuadrangle3D <: AbstractShapeFunction3D end
struct LinearQuadrangle2D <: AbstractShapeFunction2D end

"""
Cuboid mesh.
```
struct RegularCuboidMesh{
    T1 <: AbstractArray{T3, N} where {T3, N},
    T2 <: AbstractArray{T4, N} where {T4, N},
}
    label::T1
    sizeElem::T2
    sizeMesh::T2
    stiffTensor::T2
    coord::T2
    vertices::T2
end
```
"""
struct RegularCuboidMesh{
    T1 <: AbstractArray{T3, N} where {T3, N},
    T2 <: AbstractArray{T4, N} where {T4, N},
}
    label::T1
    sizeElem::T2
    sizeMesh::T2
    stiffTensor::T2
    coord::T2
    vertices::T2
end