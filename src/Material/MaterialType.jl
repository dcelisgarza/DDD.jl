"""
```
abstract type AbstractCrystalStruct end
struct BCC <: AbstractCrystalStruct end
struct FCC <: AbstractCrystalStruct end
struct HCP <: AbstractCrystalStruct end
```
Crystal structure types.
"""
abstract type AbstractCrystalStruct end
struct BCC <: AbstractCrystalStruct end
struct FCC <: AbstractCrystalStruct end
struct HCP <: AbstractCrystalStruct end

"""
```
struct MaterialP{T1, T2}
    μ::T1
    μMag::T1
    ν::T1
    E::T1
    crystalStruct::T2
end
```
"""
struct MaterialP{T1, T2}
    μ::T1
    μMag::T1
    ν::T1
    E::T1
    omνInv::T1
    νomνInv::T1
    μ4π::T1
    μ8π::T1
    μ4πν::T1
    crystalStruct::T2
    σPN::T1
end
@inline function MaterialP(; μ, μMag, ν, E, crystalStruct::AbstractCrystalStruct, σPN = 0.0)
    omνInv = 1 / (1 - ν)
    νomνInv = ν * omνInv
    μ4π = μ / (4π)
    μ8π = μ4π / 2
    μ4πν = μ4π * omνInv
    return MaterialP(μ, μMag, ν, E, omνInv, νomνInv, μ4π, μ8π, μ4πν, crystalStruct, σPN)
end

# function zero(::Type{MaterialP})
#     return MaterialP(0.0, 0.0, 0.0, :empty)
# end
