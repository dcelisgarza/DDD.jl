## IO

### Input

```@docs
loadJSON
loadDislocationLoop
loadMaterialParameters
loadFEMParameters
loadBoundaries
loadForceDisplacement
loadIntegrationParameters
loadSlipSystem
loadDislocationParameters
loadNetwork
loadIntegrationTime
loadParameters
```

### Output

```@docs
saveJSON
```

## Types

```@docs
nodeTypeDln
AbstractDlnSeg
AbstractDlnStr
AbstractDistribution
AbstractMobility
SlipSystem{T1,T2,T3}
DislocationParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21}
DislocationLoop{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14}
DislocationLoopCollection
DislocationNetwork{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14}
nodeTypeFE
AbstractMesh
AbstractElementOrder
AbstractShapeFunction
AbstractModel
FEMParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9}
RegularCuboidMesh{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21,T22,T23,T24}
ForceDisplacement{T1,T2,T3,T4}
Boundaries{T1,T2,T3,T4,T5,T6,T7}
AbstractIntegrator
IntegrationParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
IntegrationTime{T1,T2,T3}
AbstractCrystalStruct
MaterialParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13}
```

## Constructors

```@docs
SlipSystem(;
    crystalStruct::AbstractCrystalStruct,
    slipPlane::AbstractArray,
    bVec::AbstractArray
)
DislocationParameters(;
    mobility::AbstractMobility,
    dragCoeffs = (edge = 1.0, screw = 2.0, climb = 1e9),
    coreRad = 1.0,
    coreRadMag = 1.0,
    coreEnergy = 1 / (4 * π) * log(coreRad / 0.1),
    minSegLen = 2 * coreRad,
    maxSegLen = 20 * coreRad,
    minArea = coreRad^2 / sqrt(2),
    maxArea = 100 * minArea,
    slipStepCritLen = maxSegLen / 2,
    slipStepCritArea = 0.5 * (slipStepCritLen^2) * sind(1),
    remesh = true,
    collision = true,
    separation = true,
    virtualRemesh = true,
    parCPU = false,
    parGPU = false,
)
DislocationLoop(
    loopType::AbstractDlnStr,
    numSides,
    nodeSide,
    numLoops,
    segLen,
    slipSystem,
    _slipPlane,
    _bVec,
    label,
    buffer,
    range,
    dist,
)
DislocationLoop(
    loopType::loopPure,
    numSides,
    nodeSide,
    numLoops,
    segLen,
    slipSystem,
    _slipPlane::AbstractArray,
    _bVec::AbstractArray,
    label::AbstractVector{nodeTypeDln},
    buffer,
    range,
    dist::AbstractDistribution,
)
DislocationLoop(
    loopType::loopImpure,
    numSides,
    nodeSide,
    numLoops,
    segLen,
    slipSystem,
    _slipPlane::AbstractArray,
    _bVec::AbstractArray,
    label::AbstractVector{nodeTypeDln},
    buffer,
    range,
    dist::AbstractDistribution,
)
DislocationLoop(;
    loopType::AbstractDlnStr,
    numSides,
    nodeSide,
    numLoops,
    segLen,
    slipSystem,
    _slipPlane,
    _bVec,
    label,
    buffer,
    range,
    dist,
)
DislocationNetwork(;
    links::AbstractArray,
    slipPlane::AbstractArray,
    bVec::AbstractArray,
    coord::AbstractArray,
    label::AbstractVector{nodeTypeDln},
    nodeVel::AbstractArray,
    nodeForce::AbstractArray,
    numNode = length(label),
    numSeg = size(links, 2),
    maxConnect = 4,
    connectivity::AbstractArray = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::AbstractArray = zeros(Int, 2, size(links, 2)),
    segIdx::AbstractArray = zeros(Int, size(links, 2), 3),
    segForce::AbstractArray = zeros(3, 2, size(links, 2)),
)
DislocationNetwork(
    sources::DislocationLoopCollection,
    maxConnect = 4,
    args...;
    memBuffer = nothing,
    checkConsistency = true,
    kw...,
)
DislocationNetwork!
makeNetwork!
FEMParameters(; 
    type::AbstractMesh,
    order::AbstractElementOrder,
    model::AbstractModel,
    dx, dy, dz, mx, my, mz
)
ForceDisplacement(; uTilde, uHat, u, fTilde, fHat, f)
Boundaries(; uGamma, tGamma, mGamma, uDofs, tDofs, mDofs, tK)
Boundaries(
    ::FEMParameters{T1,T2,T3,T4,T5} where {T1,T2,T3<:CantileverLoad,T4,T5},
    femMesh::RegularCuboidMesh; 
    kw...
)
buildMesh
RegularCuboidMesh(
    matParams::MaterialParameters,
    femParams::FEMParameters{F1,F2,F3,F4,F5} where {F1<:DispatchRegularCuboidMesh,F2<:LinearElement,F3,F4,F5}
)
IntegrationParameters(;
    method::AbstractIntegrator,
    tmin = 0.0,
    tmax = 1e13,
    dtmin = 1e-3,
    dtmax = Inf,
    abstol = 1e-6,
    reltol = 1e-6,
    maxchange = 1.2,
    exponent = 20.0,
    maxiter = 10,
)
IntegrationTime(; dt = 0.0, time = 0.0, step = 0)
MaterialParameters(;
    crystalStruct::AbstractCrystalStruct,
    μ = 1.0,
    μMag = 1.0,
    ν = 0.5,
    σPN = 0.0,
)
```

## Processing

```@docs
calc_σTilde
calc_σTilde!
calcSegForce
calcSegForce!
calc_σHat
calcPKForce
calcPKForce!
calcSelfForce
calcSelfForce!
calcSegSegForce
calcSegSegForce!
dlnMobility
dlnMobility!
splitNode!
refineNetwork!
removeNode!
removeConnection!
removeLink!
mergeNode!
coarsenNetwork!
findIntersectVolume
makeSurfaceNode!
remeshSurfaceNetwork!
coarsenVirtualNetwork!
shapeFunction
shapeFunctionDeriv
deriv!
integrate!
```

## Utility

```@docs
loopDistribution
limits!
translatePoints!
makeSegment
makeConnect
makeConnect!
getSegmentIdx
getSegmentIdx!
checkNetwork
```

## Post Processing

```@docs
plotNodes
plotNodes!
plotFEDomain
```

## Misc support functions

```@docs
compStruct
internalAngle
externalAngle
rot3D
⊗
linePlaneIntersect
gausslegendre(n::Integer, a, b)
```
