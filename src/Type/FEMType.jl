"""
```
@enum nodeTypeFE begin
    noneFE = 0
    corner = 1  # Corner node
    edge = 2    # Edge node
    face = 3    # Face node
    intFE = 4   # Internal node
end
```
Finite element node type.

## Instances

- `noneFE = 0`: Uninitialised
- `corner = 1`: Corner
- `edge = 2`: Edge
- `face = 3`: Face
- `intFE = 4`: Internal
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

## Fields

- `type`: Mesh type
- `order`: Element order
- `model`: Experimental model
- `dx`: Dimension in x
- `dy`: Dimension in y
- `dz`: Dimension in z
- `mx`: Elements in x
- `my`: Elements in y
- `mz`: Elements in z
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
```
Stores regular a cuboid mesh.

## Fields

- `order`: Element order
- `dx`: Size in x
- `dy`: Size in y
- `dz`: Size in z
- `mx`: Elements in x
- `my`: Elements in y
- `mz`: Elements in z
- `w`: Width
- `h`: Height
- `d`: Depth
- `scale`: Mesh scale
- `numElem`: Number of elements
- `numNode`: Number of nodes
- `C`: Stiffness tensor
- `vertices`: Vertices
- `faces`: Faces
- `faceNorm`: Face normals
- `faceMidPt`: Face mid-points
- `cornerNode`: Corner nodes set
- `edgeNode`: Edge nodes set
- `faceNode`: Face node set
- `coord`: Node coordinates
- `connectivity`: Node connectivity
- `K`: Stiffness matrix
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
struct ForceDisplacement{T1,T2,T3,T4,T5,T6}
    uTilde::T1
    uHat::T2
    u::T3
    fTilde::T4
    fHat::T5
    f::T6
end
```
Stores displacements and forces applied on the FE nodes.

## Fields

- `uTilde`: Dislocation displacements
- `uHat`: Corrective displacements
- `u`: Displacements
- `fTilde` Dislocation forces
- `fHat`: Corrective forces
- `f`: Forces
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
```
Stores the nodes and degrees of freedom upon which the different boundary conditions are applied.

## Fields

- `uGammaDln`: Nodes on which dislocation displacements are calculated
- `tGammaDln`: Nodes on which dislocation tractions are calculated
- `uDofsDln`: Degrees of freedom on which dislocation displacements are calculated
- `tDofsDln`: Degrees of freedom on which dislocation tractions are calculated
- `uGamma`: Nodes with displacement boundaries
- `tGamma`: Nodes with traction boundaries
- `mGamma`: Nodes with displacement and traction boundaries
- `uDofs`: Degrees of freedom with specified displacements
- `tDofs`: Degrees of feedom with specified tractions
- `mDofs0`: Degrees of feedom with specified displacements and tractions
- `tK1`: Stiffness matrix of traction degrees of freedom
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

struct BoundaryNode{T1,T2,T3}
    type::T1
    index::T2
    node::T3
end