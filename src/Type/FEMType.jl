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
abstract type AbstractModel end
abstract type AbstractCantileverBend <: AbstractModel end
struct CantileverLoad <: AbstractCantileverBend end
```
"""
abstract type AbstractModel end
abstract type AbstractCantileverBend <: AbstractModel end
struct CantileverLoad <: AbstractCantileverBend end
"""
```
struct FEMParameters{T1,T2,T3,T4}
    type::T1    # Mesh type.
    order::T2   # Element order.
    model::T3   # M
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
struct FEMParameters{T1,T2,T3,T4,T5}
    type::T1
    order::T2
    model::T3
    dx::T4
    dy::T4
    dz::T4
    mx::T5
    my::T5
    mz::T5
end
"""
```
struct RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15} <: AbstractRegularCuboidMesh
    order::T1           # Element order.
    vertices::T2        # Vertices.
    faces::T3           # Faces.
    faceMidPt::T4       # Face midpoints.
    faceNorm::T4        # Face normals.
    C::T5               # Stiffness tensor.
    dx::T6              # Size in x.  
    dy::T6              # Size in y.
    dz::T6              # Size in z.
    scale::T7           # Length scale.
    mx::T8              # Number of elements in x.
    my::T8              # Number of elements in y.
    mz::T8              # Number of elements in z.
    numElem::T8         # Total number of elements.
    numNode::T8         # Total number of nodes.
    w::T6               # Element width (size in x).
    h::T6               # Element height (size in y).
    d::T6               # Element depth (size in z).
    coord::T9           # Node coordinates.
    connectivity::T10   # Node connectivity.
    corners::T11        # Corner labels.
    edges::T12          # Edge labels.
    faces::T13          # Face labels.
    K::T14              # Stiffness matrix.
end
```
Stores data for a regular cuboid mesh.
"""
struct RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14} <: AbstractRegularCuboidMesh
    order::T1
    vertices::T2
    faces::T3
    faceMidPt::T4
    faceNorm::T4
    C::T5
    dx::T6
    dy::T6
    dz::T6
    scale::T7
    mx::T8
    my::T8
    mz::T8
    numElem::T8
    numNode::T8
    w::T6
    h::T6
    d::T6
    coord::T9
    connectivity::T10
    cornerNode::T11
    edgeNode::T12
    faceNode::T13
    K::T14
end

"""
```
mutable struct ForceDisplacement{T1,T2,T3,T4}
    uTilde::T1  # Dislocation displacements.
    uHat::T2    # Corrective displacements.
    u::T3       # Displacement.
    fTilde::T4  # Dislocation tractions.
    fHat::T5    # Corrective tractions.
    f::T6       # Force.
end
```
Stores displacements and forces on the FE nodes.
"""
struct ForceDisplacement{T1,T2,T3,T4,T5,T6}
    uTilde::T1
    uHat::T2
    u::T3
    fTilde::T4
    fHat::T5
    f::T6
end

"""
```
struct Boundaries{T1,T2,T3,T4,T5,T6}
    uGamma::T1  # Nodes with displacement boundaries.
    tGamma::T2  # Nodes with traction boundaries.
    mGamma::T3  # Nodes with displacement and traction boundaries.
    uDofs::T4   # Degrees of freedom with specified displacements.
    tDofs::T5   # Degrees of feedom with specified tractions.
    mDofs::T6   # Degrees of feedom with specified displacements and tractions.
end
```
Stores the nodes and degrees of freedom upon which the different boundary conditions are applied.
"""
struct Boundaries{T1,T2,T3,T4,T5,T6}
    uGamma::T1
    tGamma::T2
    mGamma::T3
    uDofs::T4
    tDofs::T5
    mDofs::T6
end