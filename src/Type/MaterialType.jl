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
struct MaterialParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15}
    crystalStruct::T1   # Crystal structure.
    μ::T2               # Shear modulus.
    μMag::3             # Magnitude of shear modulus.
    ν::T4               # Poisson ratio.
    E::T5               # Young's modulus.
    omνInv::T6          # 1 / (1 - ν)
    opνInv::T7          # 1 / (1 + ν)
    νomνInv::T8         # ν / (1 - ν)
    νopνInv::T9         # v / (1 + ν)
    μ4π::T10            # μ / (4π)
    μ8π::T11            # μ / (8π)
    μ4πν::T12           # μ / (4π (1 - ν))
    omνInv8π::T13       # 1 / (8π (1 - ν))
    om2νomνInv8π::T14   # (1 - 2 * ν) / (8π (1 - ν))
    σPN::T15            # Peierls-Nabarro stress.
end
```
Store material parameters.
"""
struct MaterialParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15}
    crystalStruct::T1
    μ::T2
    μMag::T3
    ν::T4
    omνInv::T5
    opνInv::T6
    νomνInv::T7
    νopνInv::T8
    μ4π::T9
    μ8π::T10
    μ4πν::T11
    νμ4πν::T12
    omνInv8π::T13
    om2νomνInv8π::T14
    σPN::T15
end
