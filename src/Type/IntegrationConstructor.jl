"""
```
IntegrationParameters(;
    method::AbstractIntegrator,
    tmin = 0.0,
    tmax = 1e13,
    dtmin = 1e-3,
    dtmax = Inf,
    abstol = 1e-6,
    reltol = 1e-6,
    maxchange = 1.2,
    exponent = 20.0,
    maxiter = 10,
)
```
Creates [`IntegrationParameters`](@ref).
"""
function IntegrationParameters(;
    method::AbstractIntegrator,
    tmin = 0.0,
    tmax = 1e13,
    dtmin = 1e-3,
    dtmax = Inf,
    abstol = 1e-6,
    reltol = 1e-6,
    maxchange = 1.2,
    exponent = 20.0,
    maxiter = 10,
)
    return IntegrationParameters(
            method,
            tmin,
            tmax,
            dtmin,
            dtmax,
            abstol,
            reltol,
            maxchange,
            exponent,
            maxiter,
        )
end

"""
```
IntegrationTime(; dt = 0.0, time = 0.0, step = 0)
```
Creates [`IntegrationTime`](@ref).
"""
function IntegrationTime(; dt = 0.0, time = 0.0, step = 0)
    return IntegrationTime(dt, time, step)
end
