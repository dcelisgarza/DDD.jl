"""
```
abstract type AbstractIntegrator end
struct CustomTrapezoid <:AbstractIntegrator end
```
Integrator types.
"""
abstract type AbstractIntegrator end
struct CustomTrapezoid <: AbstractIntegrator end

"""
```
mutable struct IntegrationP{T1, T2, T3}
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
This structure contains the integration parameters for the simulation.
"""
struct IntegrationP{T1, T2}
    tmin::T1
    tmax::T1
    method::T2
    abstol::T1
    reltol::T1
    dt0::T1
end
function IntegrationP(;
    tmin,
    tmax,
    method,
    abstol = 1e-6,
    reltol = 1e-6,
    dt0 = 1.0,
)
    return IntegrationP(tmin, tmax, method, abstol, reltol, dt0)
end
mutable struct IntegrationVar{T1, T2}
    dt::T1
    time::T1
    step::T2
end
function IntegrationVar(;dt = 0.0, time = 0.0, step = 0)
    return IntegrationVar(dt, time, step)
end

function deriv(du, u, p, t)
end

# function zero(::Type{IntegrationP})
#     return IntegrationP(0.0, zeros(2), :empty)
# end

# function CustomTrapezoid(
#     Network::DislocationNetwork,
#     Material::MaterialP,
#     Mesh::CuboidMesh,
#     Integration::IntegrationP,
# )
#
# end
