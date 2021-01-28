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
```

## Constructors

```@docs
SlipSystem
DislocationParameters
DislocationLoop
DislocationNetwork
DislocationNetwork!
makeNetwork!
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
