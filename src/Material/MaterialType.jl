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
    omνInv::T1
    νomνInv::T1
    μ4π::T1
    μ8π::T1
    μ4πν::T1
    crystalStruct::T2

    function MaterialP(; μ, μMag, ν, E, crystalStruct)
        omνInv = 1 / (1 - ν)
        νomνInv = ν * omνInv
        μ4π = μ / (4π)
        μ8π = μ4π / 2
        μ4πν = μ4π * omνInv
        new{typeof(μ), typeof(crystalStruct)}(
            μ,
            μMag,
            ν,
            E,
            omνInv,
            νomνInv,
            μ4π,
            μ8π,
            μ4πν,
            crystalStruct,
        )
    end # constructor
end # MaterialP

# function zero(::Type{MaterialP})
#     return MaterialP(0.0, 0.0, 0.0, :empty)
# end
