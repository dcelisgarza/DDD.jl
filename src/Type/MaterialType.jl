"""
Dislocation structure types.
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
Material parameters.
```
struct MaterialParameters{T1, T2}
    μ::T1
    μMag::T1
    ν::T1
    E::T1
    crystalStruct::T2
end
```
"""
struct MaterialParameters{T1, T2}
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