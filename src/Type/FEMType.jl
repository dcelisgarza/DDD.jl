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
abstract type AbstractElementOrder end
struct LinearElement <: AbstractElementOrder end
abstract type AbstractShapeFunction end
abstract type AbstractShapeFunction3D <: AbstractShapeFunction end
abstract type AbstractShapeFunction2D <: AbstractShapeFunction end
struct LinearQuadrangle3D <: AbstractShapeFunction3D end
struct LinearQuadrangle2D <: AbstractShapeFunction2D end

struct FEMParams{T1, T2, T3}
    order::T1
    dx::T2
    dy::T2
    dz::T2
    mx::T3
    my::T3
    mz::T3
end
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
struct RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9} <: AbstractMesh
    order::T1
    vertices::T2
    C::T3
    dx::T4
    dy::T4
    dz::T4
    mx::T5
    my::T5
    mz::T5
    w::T4
    h::T4
    d::T4
    B::T6
    coord::T7
    connectivity::T8
    K::T9
end
