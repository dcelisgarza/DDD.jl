struct MaterialP{T1 <: Float64, T2 <: AbstractCrystalStruct}
    μ::T1
    μMag::T1
    ν::T1
    crystalStruct::T2

    function MaterialP(μ, μMag, ν, crystalStruct)
        new{typeof(μ), typeof(crystalStruct)}(μ, μMag, ν, crystalStruct)
    end # constructor
end # MaterialP

# function zero(::Type{MaterialP})
#     return MaterialP(0.0, 0.0, 0.0, :empty)
# end
