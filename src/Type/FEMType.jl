"""
```
@enum nodeTypeFE begin
    noneFE = 0  # Uninitialised node
    corner = 1  # Corner node
    edge = 2    # Edge node
    face = 3    # Face node
    intFE = 4   # Internal node
end
```
Finite element node type.
"""
@enum nodeTypeFE begin
    noneFE = 0
    cornerFE = 1
    edgeFE = 2
    faceFE = 3
    intFE = 4
end
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
Finite element orders for dispatch.
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
Abstract types for dispatching different models.
"""
abstract type AbstractModel end
abstract type AbstractCantileverBend <: AbstractModel end
struct CantileverLoad <: AbstractCantileverBend end

"""
```
struct FEMParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9}
    type::T1
    order::T2
    model::T3
    dx::T4
    dy::T5
    dz::T6
    mx::T7
    my::T8
    mz::T9
end
```
Stores the finite element parameters.
"""
struct FEMParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9}
    type::T1
    order::T2
    model::T3
    dx::T4
    dy::T5
    dz::T6
    mx::T7
    my::T8
    mz::T9
end
"""
```
struct RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21,T22,T23,T24} <: AbstractRegularCuboidMesh
    order::T1           # Element order
    dx::T2              # Size in x
    dy::T3              # Size in y
    dz::T4              # Size in z
    mx::T5              # Elements in x
    my::T6              # Elements in y
    mz::T7              # Elements in z
    w::T8               # Width
    h::T9               # Height
    d::T10              # Depth
    scale::T11          # Mesh scale
    numElem::T12        # Number of elements
    numNode::T13        # Number of nodes
    C::T14              # Stiffness tensor
    vertices::T15       # Vertices
    faces::T16          # Faces
    faceNorm::T17       # Face normals
    faceMidPt::T18      # Face mid-points
    cornerNode::T19     # Corner nodes set
    edgeNode::T20       # Edge nodes set
    faceNode::T21       # Face node set
    coord::T22          # Node coordinates
    connectivity::T23   # Node connectivity
    K::T24              # Stiffness matrix
end
```
Stores regular a cuboid mesh.
"""
struct RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21,T22,T23,T24} <: AbstractRegularCuboidMesh
    order::T1
    dx::T2
    dy::T3
    dz::T4
    mx::T5
    my::T6
    mz::T7
    w::T8
    h::T9
    d::T10
    scale::T11
    numElem::T12
    numNode::T13
    C::T14
    vertices::T15
    faces::T16
    faceNorm::T17
    faceMidPt::T18
    cornerNode::T19
    edgeNode::T20
    faceNode::T21
    coord::T22
    connectivity::T23
    K::T24
end

"""
```
struct ForceDisplacement{T1,T2,T3,T4}
    uTilde::T1  # Dislocation displacements.
    uHat::T2    # Corrective displacements.
    u::T3       # Displacement.
    fTilde::T4  # Dislocation tractions.
    fHat::T5    # Corrective tractions.
    f::T6       # Force.
end
```
Stores displacements and forces applied on the FE nodes.
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
struct Boundaries{T1,T2,T3,T4,T5,T6,T7}
    uGamma::T1  # Nodes with displacement boundaries.
    tGamma::T2  # Nodes with traction boundaries.
    mGamma::T3  # Nodes with displacement and traction boundaries.
    uDofs::T4   # Degrees of freedom with specified displacements.
    tDofs::T5   # Degrees of feedom with specified tractions.
    mDofs::T6   # Degrees of feedom with specified displacements and tractions.
    tK::T7      # Stiffness matrix of traction degrees of freedom.
end
```
Stores the nodes and degrees of freedom upon which the different boundary conditions are applied.
"""
struct Boundaries{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}
    uGammaDln::T1
    tGammaDln::T2
    uDofsDln::T3
    tDofsDln::T4
    uGamma::T5
    tGamma::T6
    mGamma::T7
    uDofs::T8
    tDofs::T9
    mDofs::T10
    tK::T11
end