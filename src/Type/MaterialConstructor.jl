"""
```
MaterialParameters(;
    crystalStruct::AbstractCrystalStruct,
    μ = 1.0,
    μMag = 1.0,
    ν = 0.5,
    σPN = 0.0,
)
```
Creates [`MaterialParameters`](@ref).
"""
function MaterialParameters(;
    crystalStruct::AbstractCrystalStruct,
    μ = 1.0,
    μMag = 1.0,
    ν = 0.5,
    σPN = 0.0,
)
    omνInv = 1 / (1 - ν)
    opνInv = 1 / (1 + ν)
    νomνInv = ν * omνInv
    νopνInv = ν * opνInv
    μ4π = μ / (4π)
    μ8π = μ4π / 2
    μ4πν = μ4π * omνInv
    νμ4πν = ν * μ4πν

    return MaterialParameters(
        crystalStruct,
        μ,
        μMag,
        ν,
        omνInv,
        opνInv,
        νomνInv,
        νopνInv,
        μ4π,
        μ8π,
        μ4πν,
        νμ4πν,
        σPN,
    )
end
