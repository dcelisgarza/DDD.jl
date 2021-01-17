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
struct RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11} <: AbstractRegularCuboidMesh
    order::T1           # Element order.
    vertices::T2        # Vertices.
    faces::T3           # Faces.
    faceNorm::T4        # Face normals.
    C::T5               # Stiffness tensor.
    dx::T6              # Size in x.  
    dy::T6              # Size in y.
    dz::T6              # Size in z.
    mx::T7              # Number of elements in x.
    my::T7              # Number of elements in y.
    mz::T7              # Number of elements in z.
    numElem::T7         # Total number of elements.
    numNode::T7         # Total number of nodes.
    w::T6               # Element width (size in x).
    h::T6               # Element height (size in y).
    d::T6               # Element depth (size in z).
    B::T8               # Jacobian matrix.
    coord::T9           # Node coordinates.
    connectivity::T10   # Node connectivity.
    K::T11              # Stiffness matrix.
end
```
Stores data for a regular cuboid mesh.
"""
struct RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11} <: AbstractRegularCuboidMesh
    order::T1
    vertices::T2
    faces::T3
    faceNorm::T4
    C::T5
    dx::T6
    dy::T6
    dz::T6
    mx::T7
    my::T7
    mz::T7
    numElem::T7
    numNode::T7
    w::T6
    h::T6
    d::T6
    B::T8
    coord::T9
    connectivity::T10
    K::T11
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