# Constructors

```@docs
MaterialParameters(
    μ::T1,
    μMag::T1,
    ν::T1,
    E::T1,
    crystalStruct::T2,
    σPN::T1 = 0,
) where {T1,T2 <: AbstractCrystalStruct}
```

```@docs
MaterialParameters(;
    μ::T1,
    μMag::T1,
    ν::T1,
    E::T1,
    crystalStruct::T2,
    σPN::T1 = 0,
) where {T1,T2 <: AbstractCrystalStruct}
```