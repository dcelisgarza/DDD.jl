function MaterialParameters(μ, μMag, ν, E, crystalStruct::AbstractCrystalStruct, σPN = 0.0)
    omνInv = 1 / (1 - ν)
    νomνInv = ν * omνInv
    μ4π = μ / (4π)
    μ8π = μ4π / 2
    μ4πν = μ4π * omνInv
    return MaterialParameters(μ, μMag, ν, E, omνInv, νomνInv, μ4π, μ8π, μ4πν, crystalStruct, σPN)
end
function MaterialParameters(; μ, μMag, ν, E, crystalStruct::AbstractCrystalStruct, σPN = 0.0)
    return MaterialParameters(μ, μMag, ν, E, crystalStruct, σPN)
end