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
struct MaterialParameters{T1,T2}
    crystalStruct::T1   # Crystal structure.
    μ::T2               # Shear modulus.
    μMag::T2            # Magnitude of shear modulus.
    ν::T2               # Poisson ratio.
    E::T2               # Young's modulus.
    omνInv::T2          # 1 / (1 - ν)
    opνInv::T2          # 1 / (1 + ν)
    νomνInv::T2         # ν / (1 - ν)
    νopνInv::T2         # v / (1 + ν)
    μ4π::T2             # μ / (4π)
    μ8π::T2             # μ / (8π)
    μ4πν::T2            # μ / (4π (1 - ν))
    σPN::T2             # Peierls-Nabarro stress.
end
```
Store material parameters.
"""
struct MaterialParameters{T1,T2}
    crystalStruct::T1
    μ::T2
    μMag::T2
    ν::T2
    E::T2
    omνInv::T2
    opνInv::T2
    νomνInv::T2
    νopνInv::T2
    μ4π::T2
    μ8π::T2
    μ4πν::T2
νμ4πν::T2
    σPN::T2
end
