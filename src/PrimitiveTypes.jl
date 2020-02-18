"""
```
@enum nodeType begin
    undef = -1
    intMob = 0
    intFix = 1
    srfMob = 2
    srfFix = 3
    ext = 4
end
```
Enumerated type for dislocation nodes.
"""
@enum nodeType begin
    undef = -1
    intMob = 0
    intFix = 1
    srfMob = 2
    srfFix = 3
    ext = 4
end
"""
Overloaded functions for dislocation `nodeType`.
"""
isequal(x::Real, y::nodeType) = isequal(x, Integer(y))
isequal(x::nodeType, y::Real) = isequal(Integer(x), y)
isless(x::Real, y::nodeType) = isless(x, Integer(y))
isless(x::nodeType, y::Real) = isless(Integer(x), y)
==(x::nodeType, y::Real) = isequal(Integer(x), y)
==(x::Real, y::nodeType) = isequal(x, Integer(y))
convert(::Type{nodeType}, x::Real) = nodeType(Integer(x))
zero(::Type{nodeType}) = -1
"""
```
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end
struct segEdge <: AbstractDlnSeg end
struct segEdgeN <: AbstractDlnSeg end
struct segScrew <: AbstractDlnSeg end
struct segMixed <: AbstractDlnSeg end
```
Dislocation segment types.
"""
abstract type AbstractDlnSeg end
struct segNone <: AbstractDlnSeg end
struct segEdge <: AbstractDlnSeg end
struct segEdgeN <: AbstractDlnSeg end
struct segScrew <: AbstractDlnSeg end
struct segMixed <: AbstractDlnSeg end
length(::T) where {T<:AbstractDlnSeg} = 1
"""
```
abstract type AbstractDlnStr end
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
struct loopDln <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end
```
Idealised dislocation structure types. Unsure of how to use them.
"""
abstract type AbstractDlnStr end
struct loopPrism <: AbstractDlnStr end
struct loopShear <: AbstractDlnStr end
struct loopMixed <: AbstractDlnStr end
struct loopDln <: AbstractDlnStr end
struct loopJog <: AbstractDlnStr end
struct loopKink <: AbstractDlnStr end
"""
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
```
abstract type AbstractMaterial end
struct BCC <:AbstractMaterial end
struct FCC <:AbstractMaterial end
struct HCP <:AbstractMaterial end
```
Crystal structures.
"""
abstract type AbstractMaterial end
struct BCC <: AbstractMaterial end
struct FCC <: AbstractMaterial end
struct HCP <: AbstractMaterial end
"""
```
abstract type AbstractMobility end
struct BCC <:AbstractMobility end
struct FCC <:AbstractMobility end
struct HCP <:AbstractMobility end
```
Mobility functions.
"""
abstract type AbstractMobility end
struct mobBCC <: AbstractMobility end
struct mobFCC <: AbstractMobility end
struct mobHCP <: AbstractMobility end
"""
```
abstract type AbstractIntegrator end
struct CustomTrapezoid <:AbstractIntegrator end
```
Integrator types.
"""
abstract type AbstractIntegrator end
struct CustomTrapezoid <: AbstractIntegrator end
