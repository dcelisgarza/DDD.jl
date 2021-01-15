"""
```
MaterialParameters(
    μ::T1,
    μMag::T1,
    ν::T1,
    E::T1,
    crystalStruct::T2,
    σPN::T1 = 0,
) where {T1,T2 <: AbstractCrystalStruct}
```
Constructor for [`MaterialParameters`](@ref) automatically calculates derived quantities and sets the Peierls-Nabarro stress to a default of `σPN = 0`.
"""
function MaterialParameters(
    μ::T1,
    μMag::T1,
    ν::T1,
    E::T1,
    crystalStruct::T2,
    σPN::T1 = 0,
) where {T1,T2 <: AbstractCrystalStruct}
    omνInv = 1 / (1 - ν)
    opνInv = 1 / (1 + ν)
    νomνInv = ν * omνInv
    νopνInv = ν * opνInv
    μ4π = μ / (4π)
    μ8π = μ4π / 2
    μ4πν = μ4π * omνInv
    return MaterialParameters(
        μ,
        μMag,
        ν,
        E,
        omνInv,
        opνInv,
        νomνInv,
        νopνInv,
        μ4π,
        μ8π,
        μ4πν,
        crystalStruct,
        σPN,
    )
end
"""
```
MaterialParameters(;
    μ::T1,
    μMag::T1,
    ν::T1,
    E::T1,
    crystalStruct::T2,
    σPN::T1 = 0,
) where {T1,T2 <: AbstractCrystalStruct}
```
Keyword constructor for [`MaterialParameters`](@ref).
"""
function MaterialParameters(;
    μ::T1,
    μMag::T1,
    ν::T1,
    E::T1,
    crystalStruct::T2,
    σPN::T1 = 0.0,
) where {T1,T2 <: AbstractCrystalStruct}
    return MaterialParameters(μ, μMag, ν, E, crystalStruct, σPN)
end
