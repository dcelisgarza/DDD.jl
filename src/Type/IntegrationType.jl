"""
```
abstract type AbstractIntegrator end
struct AdaptiveEulerTrapezoid <: AbstractIntegrator end
```
Integrator types for dispatch.
"""
abstract type AbstractIntegrator end
struct AdaptiveEulerTrapezoid <: AbstractIntegrator end

"""
```
struct IntegrationParameters{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    method::T1
    tmin::T2
    tmax::T3
    dtmin::T4
    dtmax::T5
    abstol::T6
    reltol::T7
    maxchange::T8
    exponent::T9
    maxiter::T10
end
```
Stores integration parameters.
"""
struct IntegrationParameters{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10}
    method::T1
    tmin::T2
    tmax::T3
    dtmin::T4
    dtmax::T5
    abstol::T6
    reltol::T7
    maxchange::T8
    exponent::T9
    maxiter::T10
end

"""
```
struct IntegrationTime{T1,T2,T3}
    dt::T1      # Current time step.
    time::T2    # Current simulation time.
    step::T3    # Current simulation step.
end
```
Store integration time and steps.
"""
struct IntegrationTime{T1, T2, T3}
    dt::T1
    time::T2
    step::T3
end
