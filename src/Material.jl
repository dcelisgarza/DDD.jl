"""
    MaterialP{
        T1<:Union{AbstractFloat,Vector{<:AbstractFloat}},
        T2<:Union{String,Vector{String}}
    }(
        μ::T1,
        μMag::T1,
        ν::T1,
        crystalStruct::T2
    )

This structure defines material parameters. Allows for the definition of multiple materials in a single invocation by giving the arguments as Vectors (1D arrays) of the relevant type.

# Arguments

```
μ               # Young's modulus.
μMag           # Magnitude of Young's modulus.
ν               # Poisson's ratio.
crystalStruct  # Crystal structure.
```

# Examples

## Correct Declaration

### Single material
```jldoctest
julia> sample_material = MaterialP(1.0, 1e5, 0.28, "bcc")
MaterialP{Float64,String}(1.0, 100000.0, 0.28, "bcc")
```

### Multiple (two) materials
```jldoctest
julia> sample_material = MaterialP([1.0, 0.8], [1e5, 1e6], [0.28, 0.31], ["bcc", "fcc"])
MaterialP{Array{Float64,1},Array{String,1}}([1.0, 0.8], [100000.0, 1.0e6], [0.28, 0.31], ["bcc", "fcc"])
```

## Incorrect declaration
```jldoctest
julia> material = MaterialP(1.0, [1e5, 2e6], 0.28, "bcc")
ERROR: AssertionError: size(μ) == size(μMag) == size(ν) == n
```

## Immutability
```jldoctest
julia> sample_material = MaterialP(1.0, 1e5, 0.28, "bcc")
MaterialP{Float64,String}(1.0, 100000.0, 0.28, "bcc")

julia> sample_material.μ = 0.7
ERROR: setfield! immutable struct of type MaterialP cannot be changed

julia> sample_material = MaterialP(0.7, 1e5, 0.28, "bcc")
MaterialP{Float64,String}(0.7, 100000.0, 0.28, "bcc")
```
"""
struct MaterialP{T1<:Float64,T2<:AbstractMaterial}
    μ::T1
    μMag::T1
    ν::T1
    crystalStruct::T2

    function MaterialP(μ, μMag, ν, crystalStruct)
        if typeof(crystalStruct) == Symbol
            @eval crystalStruct = $crystalStruct
        end
        new{typeof(μ),typeof(crystalStruct)}(μ, μMag, ν, crystalStruct)
    end # constructor
end # MaterialP

# function zero(::Type{MaterialP})
#     return MaterialP(0.0, 0.0, 0.0, :empty)
# end
