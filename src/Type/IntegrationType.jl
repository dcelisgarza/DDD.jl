"""
```
abstract type AbstractIntegrator end
struct AdaptiveEulerTrapezoid <: AbstractIntegrator end
```
Integrator types.
"""
abstract type AbstractIntegrator end
struct AdaptiveEulerTrapezoid <: AbstractIntegrator end

"""
```
struct IntegrationParameters{T1, T2, T3}
    dt::T1      # Initial timestep.
    tmin::T1    # Minimum timestep.
    tmax::T1    # Maximum timestep.
    method::T2  # Integration method.
    abstol::T1  # Absolute tolerance.
    reltol::T1  # Relative tolerance.
    time::T1    # Starting simulation time.
    step::T3    # Starting simulation step.
end
```
Storage of integration parameters.
"""
struct IntegrationParameters{T1,T2,T3}
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
```
struct IntegrationTime{T1, T2}
    dt::T1      # Current time step.
    time::T1    # Current simulation time.
    step::T2    # Current simulation step.
end
```
Store integration time and steps.
"""
struct IntegrationTime{T1,T2}
    dt::T1
    time::T1
    step::T2
end
