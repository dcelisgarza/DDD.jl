"""
```
abstract type AbstractMesh end
abstract type AbstractRegularCuboidMesh <: AbstractMesh end
struct DispatchRegularCuboidMesh <: AbstractRegularCuboidMesh end
```
FE mesh types for dispatch.
"""
abstract type AbstractMesh end
abstract type AbstractRegularCuboidMesh <: AbstractMesh end
struct DispatchRegularCuboidMesh <: AbstractRegularCuboidMesh end

"""
```
abstract type AbstractElementOrder end
struct LinearElement <: AbstractElementOrder end
```
Element orders for dispatch.
"""
abstract type AbstractElementOrder end
struct LinearElement <: AbstractElementOrder end

"""
```
abstract type AbstractShapeFunction end
abstract type AbstractShapeFunction3D <: AbstractShapeFunction end
abstract type AbstractShapeFunction2D <: AbstractShapeFunction end
struct LinearQuadrangle3D <:AbstractShapeFunction3D end
struct LinearQuadrangle2D <:AbstractShapeFunction2D end
```
Shape function types for dispatch.
"""
abstract type AbstractShapeFunction end
abstract type AbstractShapeFunction3D <: AbstractShapeFunction end
abstract type AbstractShapeFunction2D <: AbstractShapeFunction end
struct LinearQuadrangle3D <: AbstractShapeFunction3D end
struct LinearQuadrangle2D <: AbstractShapeFunction2D end

"""
```
struct FEMParameters{T1,T2,T3,T4}
    type::T1    # Mesh type.
    order::T2   # Element order.
    dx::T3      # Size in x.
    dy::T3      # Size in y.
    dz::T3      # Size in z.
    mx::T4      # Number of elements in x.
    my::T4      # Number of elements in y.
    mz::T4      # Number of elements in z.
end
```
Stores the FE parameters.
"""
struct FEMParameters{T1,T2,T3,T4}
    type::T1
    order::T2
    dx::T3
    dy::T3
    dz::T3
    mx::T4
    my::T4
    mz::T4
end
"""
```
struct RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9} <: AbstractRegularCuboidMesh
    order::T1           # Element order.
    vertices::T2        # Vertices.
    C::T3               # Stiffness tensor.
    dx::T4              # Size in x.  
    dy::T4              # Size in y.
    dz::T4              # Size in z.
    mx::T5              # Number of elements in x.
    my::T5              # Number of elements in y.
    mz::T5              # Number of elements in z.
    numElem::T5         # Total number of elements.
    numNode::T5         # Total number of nodes.
    w::T4               # Element width (size in x).
    h::T4               # Element height (size in y).
    d::T4               # Element depth (size in z).
    B::T6               # Jacobian matrix.
    coord::T7           # Node coordinates.
    connectivity::T8    # Node connectivity.
    K::T9               # Stiffness matrix.
end
```
Stores data for a regular cuboid mesh.
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

"""
```
struct ForceDisplacement{T1,T2,T3,T4}
    u::T1       # Displacement.
    f::T2       # Force.
    uHat::T3    # Corrective displacement.
    fHat::T4    # Corrective force.
end
```
Stores displacements and forces on the FE nodes.
"""
struct ForceDisplacement{T1,T2,T3,T4}
    u::T1
    f::T2
    uHat::T3
    fHat::T4
end