function MaterialParameters(
    μ::T1,
    μMag::T1,
    ν::T1,
    E::T1,
    crystalStruct::T2,
    σPN::T1 = 0.0,
) where {T1, T2 <: AbstractCrystalStruct}
    omνInv = 1 / (1 - ν)
    νomνInv = ν * omνInv
    μ4π = μ / (4π)
    μ8π = μ4π / 2
    μ4πν = μ4π * omνInv
    return MaterialParameters(
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
        σPN,
    )
end
function MaterialParameters(;
    μ::T1,
    μMag::T1,
    ν::T1,
    E::T1,
    crystalStruct::T2,
    σPN::T1 = 0.0,
) where {T1, T2 <: AbstractCrystalStruct}
    return MaterialParameters(μ, μMag, ν, E, crystalStruct, σPN)
end
