# Constructors

```@docs
SlipSystem(crystalStruct::T1, slipPlane::T2, bVec::T2) where {T1 <: AbstractCrystalStruct,T2}
```

```@docs
SlipSystem(;crystalStruct::T1, slipPlane::T2, bVec::T2) where {T1 <: AbstractCrystalStruct,T2}
```

```@docs
DislocationParameters(
    coreRad::T1,
    coreRadMag::T1,
    minSegLen::T1,
    maxSegLen::T1,
    minArea::T1,
    maxArea::T1,
    edgeDrag::T1,
    screwDrag::T1,
    climbDrag::T1,
    lineDrag::T1,
    maxConnect::T2,
    mobility::T3,
    remesh::T4 = true,
    collision::T4 = true,
    separation::T4 = true,
    virtualRemesh::T4 = true,
    parCPU::T4 = false,
    parGPU::T4 = false,
    slipStepCritLen::T1 = maxSegLen / 2,
    slipStepCritArea::T1 = 0.5 * (slipStepCritLen^2) * sin(2 * π / 360),
) where {T1,T2 <: Int,T3 <: AbstractMobility,T4 <: Bool}
```

```@docs
DislocationParameters(;
    coreRad::T1,
    coreRadMag::T1,
    minSegLen::T1,
    maxSegLen::T1,
    minArea::T1,
    maxArea::T1,
    edgeDrag::T1,
    screwDrag::T1,
    climbDrag::T1,
    lineDrag::T1,
    maxConnect::T2,
    mobility::T3,
    remesh::T4 = true,
    collision::T4 = true,
    separation::T4 = true,
    virtualRemesh::T4 = true,
    parCPU::T4 = false,
    parGPU::T4 = false,
    slipStepCritLen::T1 = maxSegLen / 2,
    slipStepCritArea::T1 = 0.5 * (slipStepCritLen^2) * sin(2 * π / 360),
) where {T1,T2 <: Int,T3 <: AbstractMobility,T4 <: Bool}
```

```@docs
DislocationLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen,
    slipSystem::T2,
    _slipPlane,
    _bVec,
    label::T3,
    buffer,
    range,
    dist::T4,
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int,
    T3 <: AbstractVector{nodeTypeDln},
    T4 <: AbstractDistribution,
}
```

```@docs
DislocationLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen,
    slipSystem::T2,
    _slipPlane::T3,
    _bVec::T3,
    label::T4,
    buffer,
    range,
    dist::T5,
) where {
    T1 <: loopPure,
    T2 <: Int,
    T3 <: AbstractArray{T,N} where {T,N},
    T4 <: AbstractVector{nodeTypeDln},
    T5 <: AbstractDistribution,
}
```

```@docs
DislocationLoop(
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen,
    slipSystem::T2,
    _slipPlane::T3,
    _bVec::T3,
    label::T4,
    buffer,
    range,
    dist::T5,
) where {
    T1 <: loopImpure,
    T2 <: Int,
    T3 <: AbstractArray{T,N} where {T,N},
    T4 <: AbstractVector{nodeTypeDln},
    T5 <: AbstractDistribution,
}
```

```@docs
DislocationLoop(;
    loopType::T1,
    numSides::T2,
    nodeSide::T2,
    numLoops::T2,
    segLen::T3,
    slipSystem::T2,
    _slipPlane::T4,
    _bVec::T4,
    label::T5,
    buffer::T6,
    range::T7,
    dist::T8,
) where {
    T1 <: AbstractDlnStr,
    T2 <: Int,
    T3 <: Union{T where {T},AbstractArray{T,N} where {T,N}},
    T4 <: AbstractArray{T,N} where {T,N},
    T5 <: AbstractVector{nodeTypeDln},
    T6,
    T7 <: AbstractArray{T,N} where {T,N},
    T8 <: AbstractDistribution,
}
```

```@docs
DislocationNetwork(
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNode::T4 = zeros(Int, 1),
    numSeg::T4 = zeros(Int, 1),
    maxConnect::T5 = 4,
    connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::T1 = zeros(Int, 2, size(links, 2)),
    segIdx::T1 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links)...),
) where {
    T1 <: AbstractArray{T,N} where {T,N},
    T2 <: AbstractArray{T,N} where {T,N},
    T3 <: AbstractVector{nodeTypeDln},
    T4 <: Union{Int,AbstractVector{Int}},
    T5 <: Int,T6 <: AbstractArray{T,N} where {T,N},
}
```

```@docs
DislocationNetwork(;
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    nodeForce::T2,
    numNode::T4 = zeros(Int, 1),
    numSeg::T4 = zeros(Int, 1),
    maxConnect::T5 = 0,
    connectivity::T1 = zeros(Int, 1 + 2 * maxConnect, length(label)),
    linksConnect::T1 = zeros(Int, 2, size(links, 2)),
    segIdx::T1 = zeros(Int, size(links, 2), 3),
    segForce::T6 = zeros(3, size(links, 2), 0),
) where {
    T1 <: AbstractArray{T,N} where {T,N},
    T2 <: AbstractArray{T,N} where {T,N},
    T3 <: AbstractVector{nodeTypeDln},T4 <: AbstractVector{Int},
    T5 <: Int,
    T6 <: AbstractArray{T,N} where {T,N},
}
```

```@docs
DislocationNetwork(
    sources::T1,
    maxConnect::T2 = 4,
    args...;
    memBuffer = nothing,
    checkConsistency::T3 = true,
    kw...,
) where {
    T1 <: Union{T,
                AbstractVector{T},
                NTuple{N,T} where N
                } where {T <: DislocationLoop},
    T2 <: Int,
    T3 <: Bool,
}
```

```@docs
DislocationNetwork!(
    network::T1,
    sources::T2,
    maxConnect::T3 = 4,
    args...;
    checkConsistency::T4 = true,
    kw...,
) where {
    T1 <: DislocationNetwork,
    T2 <: Union{T,
                AbstractVector{T},
                NTuple{N,T} where N
                } where {T <: DislocationLoop},
    T3 <: Int,
    T4 <: Bool,
}
```

```@docs
makeNetwork!
```

### Utility

```@docs
loopDistribution
```

```@docs
limits!
```

```@docs
translatePoints!
```

```@docs
makeSegment
```

```@docs
makeConnect
```

```@docs
makeConnect!
```

```@docs
getSegmentIdx
```

```@docs
getSegmentIdx!
```

```@docs
checkNetwork
```
