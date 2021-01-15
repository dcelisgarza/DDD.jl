"""
```
SlipSystem(crystalStruct::T1, slipPlane::T2, bVec::T2) 
    where {T1 <: AbstractCrystalStruct,T2}
```
Constructor for [`SlipSystem`](@ref). Checks for orthogonality of the burgers vector and slip plane. It assumes each column corresponds to a slip system.
"""
function SlipSystem(
        crystalStruct::T1,
        slipPlane::T2,
        bVec::T2,
    ) where {T1 <: AbstractCrystalStruct,T2}
    if sum(slipPlane .!= 0) != 0 && sum(bVec .!= 0) != 0
        if ndims(slipPlane) == 1
            @assert isapprox(dot(slipPlane, bVec), 0) "SlipSystem: slip plane, n == $(slipPlane), and Burgers vector, b = $(bVec), must be orthogonal."
        else
            idx = findall(x -> !isapprox(x, 0), vec(sum(slipPlane .* bVec, dims = 1)))
            @assert isempty(idx) "SlipSystem: entries of the slip plane, n[$idx, :] = $(slipPlane[idx,:]), and Burgers vector, b[$idx, :] = $(bVec[idx,:]), are not orthogonal."
        end
    end

    return SlipSystem{T1,T2}(crystalStruct, slipPlane, bVec)
end
"""
```
SlipSystem(; crystalStruct::T1, slipPlane::T2, bVec::T2) 
    where {T1 <: AbstractCrystalStruct,T2}
```
Keyword constructor for [`SlipSystem`](@ref).
"""
function SlipSystem(; crystalStruct::T1, slipPlane::T2, bVec::T2) where {T1 <: AbstractCrystalStruct,T2}
    return SlipSystem(crystalStruct, slipPlane, bVec)
end

"""
```
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
The constructor for [`DislocationParameters`](@ref) provides a few default values and calculates derived quantities.
"""
function DislocationParameters(
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

    coreRad == minSegLen == maxSegLen == 0 ? nothing :
    @assert coreRad < minSegLen < maxSegLen
    minArea == maxArea == 0 ? nothing : @assert minArea < maxArea

    return DislocationParameters(
        coreRad,
        coreRad^2,
        coreRadMag,
        minSegLen,
        maxSegLen,
        minSegLen * 2,
        minArea,
        maxArea,
        minArea^2,
        maxArea^2,
        edgeDrag,
        screwDrag,
        climbDrag,
        lineDrag,
        maxConnect,
        mobility,
        remesh,
        collision,
        separation,
        virtualRemesh,
        parCPU,
        parGPU,
        slipStepCritLen,
        slipStepCritArea,
    )
end
"""
```
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
Keyword constructor for [`DislocationParameters`](@ref). Calls the positional constructor internally.
"""
function DislocationParameters(;
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
    return DislocationParameters(
        coreRad,
        coreRadMag,
        minSegLen,
        maxSegLen,
        minArea,
        maxArea,
        edgeDrag,
        screwDrag,
        climbDrag,
        lineDrag,
        maxConnect,
        mobility,
        remesh,
        collision,
        separation,
        virtualRemesh,
        parCPU,
        parGPU,
        slipStepCritLen,
        slipStepCritArea,
)
end

"""
```
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
    T3 <: AbstractVector{nodeType},
    T4 <: AbstractDistribution,
}
```
General constructor for [`DislocationLoop`](@ref).
"""
function DislocationLoop(
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
) where {T1 <: AbstractDlnStr,T2 <: Int,T3 <: AbstractVector{nodeType},T4 <: AbstractDistribution,}

    nodeTotal::Int = 0
    links = zeros(MMatrix{2,nodeTotal,Int})
    coord = zeros(MMatrix{3,nodeTotal})
    slipPlane = zeros(MMatrix{3,0})
    bVec = zeros(MMatrix{3,0})

    return DislocationLoop(
        loopType,
        numSides,
        nodeSide,
        numLoops,
        segLen,
        slipSystem,
        links,
        slipPlane,
        bVec,
        coord,
        label,
        buffer,
        range,
        dist
    )
end
"""
```
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
    T4 <: AbstractVector{nodeType},
    T5 <: AbstractDistribution,
}
```
Constructor for `loopPure` loops [`DislocationLoop`](@ref).
"""
function DislocationLoop(
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
) where {T1 <: loopPure,T2 <: Int,T3 <: AbstractArray{T,N} where {T,N},T4 <: AbstractVector{nodeType},T5 <: AbstractDistribution,}

    nodeTotal = numSides * nodeSide # Calculate total number of nodes for memory allocation.
    numSegLen = length(segLen) # Number of segment lengths.

    # Validate input.
    @assert length(label) == nodeTotal "DislocationLoop: All $nodeTotal nodes must be labelled. There are $(length(label)) labels currently defined."
    @assert numSegLen == nodeTotal "DislocationLoop: All $nodeTotal segments must have their lengths defined. There are $numSegLen lengths currently defined."

    # Normalise vectors.
    elemT = eltype(_slipPlane)
    _slipPlane = _slipPlane / norm(_slipPlane)
    _bVec = _bVec / norm(_bVec)

    # Pick rotation axis for segments.
    # Shear loops rotate around slip plane vector. They have screw, mixed and edge segments.
    if typeof(loopType) == loopShear
        rotAxis = SVector{3,elemT}(_slipPlane[1], _slipPlane[2], _slipPlane[3])
        # Prismatic loops rotate around Burgers vector. All segments are edge.
    else
        rotAxis = SVector{3,elemT}(_bVec[1], _bVec[2], _bVec[3])
        # Catch all.
    end

    # Allocate arrays.
    links = zeros(MMatrix{2,nodeTotal,Int})
    coord = zeros(MMatrix{3,nodeTotal})
    slipPlane = MMatrix{3,nodeTotal}(repeat(_slipPlane, inner = (1, numSegLen)))
    bVec = MMatrix{3,nodeTotal}(repeat(_bVec, inner = (1, numSegLen)))
    seg = zeros(MMatrix{3,numSegLen})

    # Create initial segments.
    staticSlipPlane = SVector{3,elemT}(_slipPlane[1], _slipPlane[2], _slipPlane[3])
    staticBVec = SVector{3,elemT}(_bVec[1], _bVec[2], _bVec[3])
    @inbounds @simd for i in eachindex(segLen)
        seg[:, i] = makeSegment(segEdge(), staticSlipPlane, staticBVec) * segLen[i]
    end

    θ = externalAngle(numSides)  # External angle of a regular polygon with numSides.

    # Loop over polygon's sides.
    origin = SVector{3,elemT}(0, 0, 0)
    @inbounds for i in 1:numSides
        # Index for side i.
        idx = (i - 1) * nodeSide
        # Rotate segments by external angle of polygon to make polygonal loop.
        modIdx = mod(i - 1, numSegLen) + 1
        staticSeg = SVector{3,elemT}(seg[1, modIdx], seg[2, modIdx], seg[3, modIdx])
        rseg = rot3D(staticSeg, rotAxis, origin, θ * (i - 1))
        # DO NOT add @simd, this loop works by adding rseg to the previous coordinate to make the loop. Loop over the nodes per side.
        for j in 1:nodeSide
            # Count first node once.
            if i == j == 1
                coord[:, 1] .= 0 # Initial coordinate is on the origin.
                continue
            end
            if idx + j <= nodeTotal
                # Add segment vector to previous coordinate.
                coord[:, idx + j] += @views coord[:, idx + j - 1] + rseg
            end
        end
    end

    # Find centre of the loop and make it zero.
    meanCoord = mean(coord, dims = 2)
    coord .-= meanCoord

    # Create links matrix.
    @inbounds @simd for j in 1:(nodeTotal - 1)
        links[:, j] .= (j, j + 1)
    end
    links[:, nodeTotal] .= (nodeTotal, 1)

    return DislocationLoop(
        loopType,
        numSides,
        nodeSide,
        numLoops,
        segLen,
        slipSystem,
        links,
        slipPlane,
        bVec,
        coord,
        label,
        buffer,
        range,
        dist,
    )
end
"""
```
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
    T4 <: AbstractVector{nodeType},
    T5 <: AbstractDistribution,
}
```
A fallback [`DislocationLoop`](@ref) constructor for other as of yet unimplemented `loopType`.
"""
function DislocationLoop(
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
) where {T1 <: loopImpure,T2 <: Int,T3 <: AbstractArray{T,N} where {T,N},T4 <: AbstractVector{nodeType},T5 <: AbstractDistribution,}

    @warn "DislocationLoop: Constructor for $(typeof(loopType)) not defined, defaulting to prismatic loop."
    return DislocationLoop(
        loopPrism(),
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
end
"""
```
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
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T,N} where {T,N},
    T8 <: AbstractDistribution,
}
end
```
Keyword [`DislocationLoop`](@ref) constructor calls the positional constructor, which dispatches the appropriate method.
"""
function DislocationLoop(;
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
) where {T1 <: AbstractDlnStr,T2 <: Int,T3 <: Union{T where {T},AbstractArray{T,N} where {T,N}},T4 <: AbstractArray{T,N} where {T,N},T5 <: AbstractVector{nodeType},T6,T7 <: AbstractArray{T,N} where {T,N},T8 <: AbstractDistribution,}
    return DislocationLoop(
        loopType,
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
        dist
    )

end