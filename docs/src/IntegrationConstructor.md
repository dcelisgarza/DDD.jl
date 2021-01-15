# Constructors

```@docs
IntegrationParameters(
    method::T1,
    tmin::T2,
    tmax::T2,
    dtmin::T2 = 1e-3,
    dtmax::T2 = Inf,
    abstol::T2 = 1e-6,
    reltol::T2 = 1e-6,
    maxchange::T2 = 1.2,
    exponent::T2 = 20.0,
    maxiter::T3 = 10,
) where {
    T1 <: AbstractIntegrator,
    T2 <: AbstractFloat,
    T3 <: Integer
}
```

```@docs
IntegrationParameters(;
    method::T1,
    tmin::T2,
    tmax::T2,
    dtmin::T2 = 1e-3,
    dtmax::T2 = Inf,
    abstol::T2 = 1e-6,
    reltol::T2 = 1e-6,
    maxchange::T2 = 1.2,
    exponent::T2 = 20.0,
    maxiter::T3 = 10,
) where {
    T1 <: AbstractIntegrator,
    T2 <: AbstractFloat,
    T3 <: Integer
}
```

```@docs
IntegrationTime(dt::T1, time::T1, step::T2) 
    where {T1 <: AbstractFloat,T2 <: Integer}
```

```@docs
IntegrationTime(; dt::T1, time::T1, step::T2) 
    where {T1 <: AbstractFloat,T2 <: Integer}
```
