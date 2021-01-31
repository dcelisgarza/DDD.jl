"""
```
MaterialParameters(;
    crystalStruct::AbstractCrystalStruct,
    μ = 1.0,
    μMag = 1.0,
    ν = 0.2,
    σPN = 0.0,
)
```
Creates [`MaterialParameters`](@ref).
"""
function MaterialParameters(;
    crystalStruct::AbstractCrystalStruct,
    μ = 1.0,
    μMag = 1.0,
    ν = 0.2,
    σPN = 0.0,
)
    omνInv = 1 / (1 - ν)
    opνInv = 1 / (1 + ν)
    νomνInv = ν * omνInv        # ν / (1 - ν)
    νopνInv = ν * opνInv        # ν / (1 + ν)
    μ4π = μ / (4π)              # μ / (4π)
    μ8π = μ4π / 2               # μ / (8π)
    μ4πν = μ4π * omνInv         # μ / (1 - ν)
    νμ4πν = ν * μ4πν            # νμ / (1 - ν)
    omνInv8π = omνInv / (8π)    # 1 / (8π (1 - ν))
    om2νomνInv8π = (1 - 2ν) * omνInv8π # (1 - 2ν) / (8π (1 - ν))

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
        omνInv8π,
        om2νomνInv8π,
        σPN,
    )
end
