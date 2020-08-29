"""
Integrator type.
```
abstract type AbstractIntegrator end
struct CustomTrapezoid <:AbstractIntegrator end
```
Integrator types.
"""
abstract type AbstractIntegrator end
struct CustomTrapezoid <: AbstractIntegrator end

"""
Integration parameters.
```
struct IntegrationParameters{T1, T2, T3}
    dt::T1
    tmin::T1
    tmax::T1
    method::T2
    abstol::T1
    reltol::T1
    time::T1
    step::T3
end
```
"""
struct IntegrationParameters{T1, T2, T3}
    method::T1
    tmin::T2
    tmax::T2
    dtmin::T2
    dtmax::T2
    abstol::T2
    reltol::T2
    maxchange::T2
    exponent::T2
    maxiter::T3
end

"""
Structure to keep track of integration time.
```
mutable struct IntegrationTime{T1, T2}
    dt::T1
    time::T1
    step::T2
end
```
"""
mutable struct IntegrationTime{T1, T2}
    dt::T1
    time::T1
    step::T2
end