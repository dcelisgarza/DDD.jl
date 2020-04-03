"""
```
abstract type AbstractCrystalStruct end
struct BCC <: AbstractCrystalStruct end
struct FCC <: AbstractCrystalStruct end
struct HCP <: AbstractCrystalStruct end
```
Crystal structures.
"""
abstract type AbstractCrystalStruct end
struct BCC <: AbstractCrystalStruct end
struct FCC <: AbstractCrystalStruct end
struct HCP <: AbstractCrystalStruct end
"""
```
MaterialP{T1 <: Float64, T2 <: AbstractCrystalStruct}
    μ::T1
    μMag::T1
    ν::T1
    E::T1
    crystalStruct::T2
```
"""
struct MaterialP{T1 <: Float64, T2 <: AbstractCrystalStruct}
    μ::T1
    μMag::T1
    ν::T1
    E::T1
    crystalStruct::T2

    function MaterialP(μ, μMag, ν, E, crystalStruct)
        new{typeof(μ), typeof(crystalStruct)}(μ, μMag, ν, E, crystalStruct)
    end # constructor
end # MaterialP

# function zero(::Type{MaterialP})
#     return MaterialP(0.0, 0.0, 0.0, :empty)
# end
