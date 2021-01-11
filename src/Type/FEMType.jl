"""
Mesh types.
```
abstract type AbstractMesh end
```
"""
abstract type AbstractMesh end
abstract type AbstractRegularCuboidMesh <: AbstractMesh end
struct DispatchRegularCuboidMesh <: AbstractRegularCuboidMesh end

"""
Element orders.
"""
abstract type AbstractElementOrder end
struct LinearElement <: AbstractElementOrder end

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

struct FEMParameters{T1,T2,T3,T4}
    type::T1
    order::T2
    dx::T3
    dy::T3
    dz::T3
    mx::T4
    my::T4
    mz::T4
    function FEMParameters(type::T1, order::T2, dx::T3, dy::T3, dz::T3, mx::T4, my::T4, mz::T4) where {T1 <: AbstractMesh,T2 <: AbstractElementOrder,T3 <: AbstractFloat,T4 <: Integer}
        new{T1,T2,T3,T4}(type, order, dx, dy, dz, mx, my, mz)
    end
end
function FEMParameters(; type::T1, order::T2, dx::T3, dy::T3, dz::T3, mx::T4, my::T4, mz::T4) where {T1 <: AbstractMesh,T2 <: AbstractElementOrder,T3 <: AbstractFloat,T4 <: Integer}
    return FEMParameters(type, order, dx, dy, dz, mx, my, mz)
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
struct RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9} <: AbstractRegularCuboidMesh
    order::T1
    vertices::T2
    C::T3
    dx::T4
    dy::T4
    dz::T4
    mx::T5
    my::T5
    mz::T5
    numElem::T5
    numNode::T5
    w::T4
    h::T4
    d::T4
    B::T6
    coord::T7
    connectivity::T8
    K::T9
end

struct ForceDisplacement{T1,T2,T3,T4}
    u::T1       # Displacement.
    f::T2       # Force.
    uHat::T3    # Corrective displacement.
    fHat::T4    # Force displacement.
end