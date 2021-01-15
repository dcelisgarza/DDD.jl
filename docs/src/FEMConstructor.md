# Constructors

```@docs
FEMParameters(type::T1, order::T2, dx::T3, dy::T3, dz::T3, mx::T4, my::T4, mz::T4) where {T1 <: AbstractMesh,T2 <: AbstractElementOrder,T3 <: AbstractFloat,T4 <: Integer}
```

```@docs
FEMParameters(; type::T1, order::T2, dx::T3, dy::T3, dz::T3, mx::T4, my::T4, mz::T4) where {T1 <: AbstractMesh,T2 <: AbstractElementOrder,T3 <: AbstractFloat,T4 <: Integer}
```

```@docs
buildMesh
```

```@docs
RegularCuboidMesh(
    order::T1, 
    matParams::T2, 
    femParams::T3
) where {
    T1 <: LinearElement,
    T2 <: MaterialParameters,
    T3 <: FEMParameters
}
```

```@docs
RegularCuboidMesh(;
    order::T1, 
    matParams::T2, 
    femParams::T3
) where {
    T1 <: AbstractElementOrder,
    T2 <: MaterialParameters,
    T3 <: FEMParameters
}
```
