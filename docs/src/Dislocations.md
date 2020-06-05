# Dislocation Generation

## Types, Structs and Constructors

Dislocations are described, generated and validated by custom types, structures and functions. By subtyping the provided types with new concrete types, users can define functions which dispatch specifically on their new types while minimising the need for code rewrites, as multiple dispatch takes care of everything during JIT compilation. Structures have not had their default constructors overwritten, we provide custom constructors whose use is recommended instead.

In discrete dislocation dynamics, dislocations are described by nodes connected by segments. The nodes are labelled according to their type, which is used by the software to decide how they are treated. However, labels are discrete variables, so they cannot take on any value. Additionally, accidentally using non-existent node types may produce silent and difficult to track errors. It is also impractical to validate node types at runtime. We solve these issues by defining a custom enumerated type, which not only limits possible values but informs users and developers of what the values represent.
```@docs
nodeType
```
Of course, mislabelling a node with an erroneous but defined value may still occur. Preventing such bugs is the task of users and developers, however the problem may is eased by the self-descriptive nature of enumerated types.

Dislocations also have different idealised segment types which are characterised by the relationship between the line direction and Burgers vector. These idealised types are used in the code for in loop generation. We've defined a few common types, most of which are not currently used but may prove useful in the future, for example in the statistical analysis of the dislocation network.
```@docs
AbstractDlnSeg
```

Dislocation loops are idealised as having different classifications. Prismatic loops are made up only of edge segments generally with the same slip system; shear loops are made up of a mixture of segment types with the same slip system; jogs and kinks are steps not contained in the slip plane. These idealisations can be used to automate loop generation with minimal rewriting via multiple dispatch. We provide the following types for such a purpose.
```@docs
AbstractDlnStr
```

Generating dislocation networks often involves distributing the initial loops within the simulation domain in a particular way. We again define custom structures that enable us to make use of multiple dispatch.
```@docs
AbstractDistribution
```

Mobility functions describe how dislocations move within a material. There are many variations of such functions and users may want/need to use different functions for different purposes. Creating concrete mobility types lets users define a function for their new concrete mobility type and carry on with their lives.
```@docs
AbstractMobility
```

As previously mentioned, idealised dislocation segments live on slip systems, which are pairings of slip plane and Burgers vector. These can be stored in the following structure.
```@docs
SlipSystem
```
As slip systems are defined with respect to pure edge dislocations, we recommended users use the keyword constructor as it validates the orthogonality of paired Burgers vector and slip plane.
```@docs
SlipSystem(;
    crystalStruct::T1,
    slipPlane::T2,
    bVec::T2,
) where {T1 <: AbstractCrystalStruct, T2 <: AbstractArray{T, N} where {T, N}}
```

The simulation requires certain parameters pertaining to the dislocation network being modelled. These values control certain aspects of the simulation and are stored in the following structure.
```@docs
DislocationP
```
We recommend the use of the keyword constructor as it performs sanity checks on various parameters, ensuring a hierarchy of values is maintained.
```@docs
DislocationP(;
    coreRad::T1,
    coreRadMag::T1,
    minSegLen::T1,
    maxSegLen::T1,
    minArea::T1,
    maxArea::T1,
    maxConnect::T2,
    remesh::T3,
    collision::T3,
    separation::T3,
    virtualRemesh::T3,
    edgeDrag::T1,
    screwDrag::T1,
    climbDrag::T1,
    lineDrag::T1,
    mobility::T4,
) where {T1, T2 <: Int, T3 <: Bool, T4 <: AbstractMobility}
```

Dislocation loops form the basis of a network, we provide a structure to store these loops, whether idealised or otherwise.
```@docs
DislocationLoop
```
Again we recommend the use of the key-word only constructor as it provides an interface to call specific constructors which dispatch on `loopType` and generate the loop automatically. The default constructor should only be used for truly custom loops, however creating a concrete subtype of [`AbstractDlnStr`](@ref) and a constructor which dispatches on this new type is highly recommended since it facilitates testing, reproducibility and seamlessly integrates with existing infrastructure.
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
    T3 <: Union{T where {T}, AbstractArray{T, N} where {T, N}},
    T4 <: AbstractArray{T, N} where {T, N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution,
}
```
The first concrete [`DislocationLoop`](@ref) constructor is the "zero" constructor. The way Julia's multiple dispatch works is by dispatching on the most specific method for the inputs provided. Therefore this constructor will be called whenever the loop type is `loopDln()`.
```@docs
DislocationLoop(
    loopType::T1;
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
    T1 <: loopDln,
    T2 <: Int,
    T3 <: Union{T where {T}, AbstractArray{T, N} where {T, N}},
    T4 <: AbstractArray{T, N} where {T, N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution,
}
```
We also provide a catch-all constructor that generates shear or prismatic loops depending on whether `loopType` is `loopShear()` or `loopPrism()`. It acts as a fallback for other loop types but generates prismatic loops. Such behaviour can be overridden by defining new methods which dispatch on more a specific `loopType`.
```@docs
DislocationLoop(
    loopType::T1;
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
    T3 <: Union{T where {T}, AbstractArray{T, N} where {T, N}},
    T4 <: AbstractArray{T, N} where {T, N},
    T5 <: AbstractVector{nodeType},
    T6,
    T7 <: AbstractArray{T, N} where {T, N},
    T8 <: AbstractDistribution,
}
```
The `dist` parameter refers to a concrete subtype of [`AbstractDistribution`](@ref). When defining a new distribution it is important to define a new version of the [`loopDistribution`](@ref) function.
```@docs
loopDistribution
```

The dislocation loops contain all the data relevant to a single loop. This data is then used to populate a dislocation network, which is a mutable structure because it evolves over time.
```@docs
DislocationNetwork
```
Again we provide a keyword constructor which performs some sanity checks and loads the data into the structure. This is the constructor to use when loading data from a previously generated network.
```@docs
DislocationNetwork(;
    links::T1,
    slipPlane::T2,
    bVec::T2,
    coord::T2,
    label::T3,
    nodeVel::T2,
    numNode::T4 = 0,
    numSeg::T4 = 0,
    maxConnect::T4 = 0,
    connectivity::T5 = zeros(Int, 0, 0),
    linksConnect::T5 = zeros(Int, 2, 0),
    segIdx::T5 = zeros(Int, 2, 3),
    segForce::T6 = zeros(3, 2, 0),
) where {
    T1 <: AbstractArray{T, N} where {T, N},
    T2 <: AbstractArray{T, N} where {T, N},
    T3 <: AbstractVector{nodeType},
    T4 <: Int,
    T5 <: AbstractArray{Int, N} where {N},
    T6 <: AbstractArray{T, N} where {T, N},
}
```
However, if the aim is to generate a new network then use the following constructor.
```@docs
DislocationNetwork(
    sources::T1,
    maxConnect::T2 = 4,
    args...;
    memBuffer = nothing,
    checkConsistency::T3 = true,
    kw...,
) where {
    T1 <: Union{T, AbstractVector{T}} where {T <: DislocationLoop},
    T2 <: Int,
    T3 <: Bool,
}
```
If adding to an existing network, use the mutating (also called in-place) constructor.
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
    T2 <: Union{T, AbstractVector{T}} where {T <: DislocationLoop},
    T3 <: Int,
    T4 <: Bool,
}
```

Dislocation network constructors use a few internal functions to distribute loops about the domain as well as create auxiliary matrices and verify the integrity of the generated network. As previously mentioned, [`loopDistribution`](@ref) is used to generate points from a particular distribution. These points must be scaled and adjusted by limits generated the `limits` function.
```@docs
limits!(
    lims::T1,
    segLen::T2,
    range::T1,
    buffer::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2}
```
The limits, together with the aforementioned distributions are used to translate coordinates with the `translatePoints` function.
```@docs
translatePoints(
    coord::T1,
    lims::T1,
    disp::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: AbstractVector{T} where {T}}
```

In order to traverse the network, it is useful to define a few auxiliary matrices containing relational information about nodes and links. These are created by the `makeConnect` functions.
```@docs
makeConnect(
    links::T1,
    maxConnect::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: Int}

makeConnect!(network::DislocationNetwork)
```

It's also useful to define another matrix for indexing segments quickly, this matrix is defined by the `getSegmentIdx` functions.
```@docs
getSegmentIdx(
    links::T1,
    label::T2,
) where {T1 <: AbstractArray{T, N} where {T, N}, T2 <: AbstractVector{nodeType}}

getSegmentIdx!(network::DislocationNetwork)
```

The validity of the network can be checked by `checkNetwork`.
```@docs
checkNetwork(network::DislocationNetwork)
```
