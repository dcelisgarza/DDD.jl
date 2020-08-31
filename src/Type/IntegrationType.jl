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

    function IntegrationParameters(method::T1, tmin::T2, tmax::T2, dtmin::T2 = 1e-3, dtmax::T2 = Inf, abstol::T2 = 1e-6, reltol::T2 = 1e-6, maxchange::T2 = 1.2, exponent::T2 = 20.0, maxiter::T3 = 10) where {T1 <: AbstractIntegrator, T2 <: AbstractFloat, T3 <: Integer}
        new{T1, T2, T3}(method, tmin, tmax, dtmin, dtmax, abstol, reltol, maxchange, exponent, maxiter)
    end
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
struct IntegrationTime{T1, T2}
    dt::T1
    time::T1
    step::T2
    function IntegrationTime(dt::T1, time::T1, step::T2) where {T1 <: AbstractFloat, T2 <: Integer}
        new{T1, T2}(dt, time, step)
    end
end