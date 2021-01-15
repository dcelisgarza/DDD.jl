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
    μ::T1               # Shear modulus.
    μMag::T1            # Magnitude of shear modulus.
    ν::T1               # Poisson ratio.
    E::T1               # Young's modulus.
    omνInv::T1          # 1 / (1 - ν)
    opνInv::T1          # 1 / (1 + ν)
    νomνInv::T1         # ν / (1 - ν)
    νopνInv::T1         # v / (1 + ν)
    μ4π::T1             # μ / (4π)
    μ8π::T1             # μ / (8π)
    μ4πν::T1            # μ / (4π (1 - ν))
    crystalStruct::T2   # Crystal structure.
    σPN::T1             # Peierls-Nabarro stress.
end
```
Store material parameters.
"""
struct MaterialParameters{T1,T2}
    μ::T1
    μMag::T1
    ν::T1
    E::T1
    omνInv::T1
    opνInv::T1
    νomνInv::T1
    νopνInv::T1
    μ4π::T1
    μ8π::T1
    μ4πν::T1
    crystalStruct::T2
    σPN::T1
end
