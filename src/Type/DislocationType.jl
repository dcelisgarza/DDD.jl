## Dislocation Type Declarations
"""
Since there is only a finite number of node types making the node type an enumerated type provides safety.
```
@enum nodeType begin
    none = 0    # Undefined node, value at initialisation.
    intMob = 1  # Internal mobile node.
    intFix = 2  # Internal fixed node.
    srfMob = 3  # Mobile surface node.
    srfFix = 4  # Fixed surface node.
    ext = 5     # External node.
end
```
"""
@enum nodeType begin
    none = 0
    intMob = 1
    intFix = 2
    srfMob = 3
    srfFix = 4
    ext = 5
end

"""
Dislocation segment types.
```
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end    # Undefined segment
struct segEdge <: AbstractDlnSeg end    # Edge segment
struct segEdgeN <: AbstractDlnSeg end   # Edge segment
struct segScrew <: AbstractDlnSeg end   # Screw segment
struct segMixed <: AbstractDlnSeg end   # Mixed segment
```
where `segEdge` have ``(\\bm{b} \\perp \\bm{t}) \\perp \\bm{n}``, `segEdgeN` have ``(\\bm{b} \\perp \\bm{t}) \\parallel \\bm{n}``, `segScrew` have ``\\bm{b} \\parallel \\bm{t}``, `segMixed` have ``\\bm{b} \\not\\perp \\bm{t}~ \\&~ \\bm{b} \\not\\perp \\bm{n}`` and ``\\bm{b}`` is the Burgers vector and ``\\bm{n}`` the slip plane.
"""
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end
struct segEdge <: AbstractDlnSeg end
struct segEdgeN <: AbstractDlnSeg end
struct segScrew <: AbstractDlnSeg end
struct segMixed <: AbstractDlnSeg end

"""
Dislocation loop types.
```
abstract type AbstractDlnStr end
struct loopDln <: AbstractDlnStr end    # Unclassified loop
struct loopPrism <: AbstractDlnStr end  # Prismatic loop
struct loopShear <: AbstractDlnStr end  # Shear loop
struct loopJog <: AbstractDlnStr end    # Jog
struct loopKink <: AbstractDlnStr end   # Kink
```
"""
abstract type AbstractDlnStr end
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
struct loopMixed <: AbstractDlnStr end
struct loopDln <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end

"""
Initial statistical distributions of dislocations in the domain.
```
abstract type AbstractDistribution end
struct Zeros <: AbstractDistribution end
struct Rand <: AbstractDistribution end
struct Randn <: AbstractDistribution end
struct Regular <: AbstractDistribution end
```
Distributions for dislocation sources.
"""
abstract type AbstractDistribution end
struct Zeros <: AbstractDistribution end
struct Rand <: AbstractDistribution end
struct Randn <: AbstractDistribution end
struct Regular <: AbstractDistribution end

"""
Dislocation mobility types.
```
abstract type AbstractMobility end
struct mobBCC <: AbstractMobility end
struct mobFCC <: AbstractMobility end
struct mobHCP <: AbstractMobility end
```
Mobility functions.
"""
abstract type AbstractMobility end
struct mobBCC <: AbstractMobility end
struct mobFCC <: AbstractMobility end
struct mobHCP <: AbstractMobility end

"""
Slip systems.
```
struct SlipSystem{T1, T2}
    crystalStruct::T1   # Crystal structure
    slipPlane::T2       # Slip plane
    bVec::T2            # Burgers vector
end
```
Structure to store slip systems.
"""
struct SlipSystem{T1, T2}
    crystalStruct::T1
    slipPlane::T2
    bVec::T2

    function SlipSystem(crystalStruct::T1, slipPlane::T2, bVec::T2) where {
        T1 <: AbstractCrystalStruct,
        T2 <: AbstractArray{T, N} where {T, N}
    }
        if sum(slipPlane .!= 0) != 0 && sum(bVec .!= 0) != 0
            if ndims(slipPlane) == 1
                @assert isapprox(dot(slipPlane, bVec), 0) "SlipSystem: slip plane, n == $(slipPlane), and Burgers vector, b = $(bVec), must be orthogonal."
            else
                idx = findall(x -> !isapprox(x, 0), vec(sum(slipPlane .* bVec, dims = 1)))
                @assert isempty(idx) "SlipSystem: entries of the slip plane, n[$idx, :] = $(slipPlane[idx,:]), and Burgers vector, b[$idx, :] = $(bVec[idx,:]), are not orthogonal."
            end
        end

        new{T1, T2}(crystalStruct, slipPlane, bVec)
    end
    
end

"""
Dislocation loop parameters.
```
struct DislocationParameters{T1, T2, T3, T4}
    coreRad::T1         # Core radius
    coreRadSq::T1       # Square of core radius
    coreRadMag::T1      # Magnitude of core radius
    minSegLen::T1       # Minimum segment length
    maxSegLen::T1       # Maximum segment length
    twoMinSegLen::T1    # Twice minimum segment length
    minArea::T1         # Minimum area enclosed by 3 segments
    maxArea::T1         # Maximum area enclosed by 3 segments
    minAreaSq::T1       # Squared min area
    maxAreaSq::T1       # Squared max area
    edgeDrag::T1        # Drag coefficient edge dislocation
    screwDrag::T1       # Drag coefficient screw dislocation
    climbDrag::T1       # Drag coefficient climb direction
    lineDrag::T1        # Drag coefficient line direction
    maxConnect::T2      # Maximum connectivity
    remesh::T3          # Remesh flag
    collision::T3       # Collision flag
    separation::T3      # Separation flag
    virtualRemesh::T3   # Virtual remeshing flag
    parCPU::T3          # Parallelise on CPU
    parGPU::T3          # Parallelise on GPU
    mobility::T4        # Mobility law
end
```
"""
struct DislocationParameters{T1, T2, T3, T4}
    coreRad::T1
    coreRadSq::T1
    coreRadMag::T1
    minSegLen::T1
    maxSegLen::T1
    twoMinSegLen::T1
    minArea::T1
    maxArea::T1
    minAreaSq::T1
    maxAreaSq::T1
    edgeDrag::T1
    screwDrag::T1
    climbDrag::T1
    lineDrag::T1
    maxConnect::T2
    mobility::T3
    remesh::T4
    collision::T4
    separation::T4
    virtualRemesh::T4
    parCPU::T4
    parGPU::T4
end

"""
Dislocation loop.
```
struct DislocationLoop{T1, T2, T3, T4, T5, T6, T7, T8, T9}
    loopType::T1    # Loop type
    numSides::T2    # Number of sides per loop
    nodeSide::T2    # Number of nodes per side
    numLoops::T2    # Number of loops to generate
    segLen::T3      # Segment lengths
    slipSystem::T4  # Slip system
    links::T5       # Links matrix
    slipPlane::T6   # Slip plane for each link
    bVec::T6        # Burgers vector for each link
    coord::T6       # Coordinates of each node
    label::T7       # Label of each node
    buffer::T8      # Mean distance buffer separating each loop centre
    range::T6       # Distribution range of generated loops
    dist::T9        # Distribution of generated loops
end
```
"""
struct DislocationLoop{T1, T2, T3, T4, T5, T6, T7, T8, T9}
    loopType::T1
    numSides::T2
    nodeSide::T2
    numLoops::T2
    segLen::T3
    slipSystem::T4
    links::T5
    slipPlane::T6
    bVec::T6
    coord::T6
    label::T7
    buffer::T8
    range::T6
    dist::T9
end

"""
Dislocation network.
```
struct DislocationNetwork{T1, T2, T3, T4, T5, T6}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeVel::T2
    nodeForce::T2
    numNodeSegConnect::T4   # Number of nodes, segments and max connectivity in network
    connectivity::T5        # Connectivity matrix
    linksConnect::T5        # Links involved in connection
    segIdx::T5              # Contains segment index and the nodes of the nodes in said link
    segForce::T6            # Force on each node of each segment
end
```
"""
struct DislocationNetwork{T1, T2, T3, T4, T5, T6}
    links::T1
    slipPlane::T2
    bVec::T2
    coord::T2
    label::T3
    nodeVel::T2
    nodeForce::T2
    numNode::T4
    numSeg::T4
    maxConnect::T5
    connectivity::T1
    linksConnect::T1
    segIdx::T1
    segForce::T6

    function DislocationNetwork(
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
        T1 <: AbstractArray{T, N} where {T, N},
        T2 <: AbstractArray{T, N} where {T, N},
        T3 <: AbstractVector{nodeType},
        T4 <: AbstractVector{Int},
        T5 <: Int,
        T6 <: AbstractArray{T, N} where {T, N},
    }

        @assert size(links, 1) == size(segForce, 2) == 2
        @assert size(bVec, 1) == size(slipPlane, 1) == size(coord, 1) size(segForce, 1) == 3
        @assert size(links, 2) == size(bVec, 2) == size(slipPlane, 2) == size(segForce, 3)
        @assert size(coord, 2) == length(label)

        new{T1, T2, T3, T4, T5, T6}(links, slipPlane, bVec, coord, label, nodeVel, nodeForce, numNode, numSeg, maxConnect, connectivity, linksConnect, segIdx, segForce)
    end
end