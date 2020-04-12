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
IntegrationP{
    T1 <: Float64,
    T2 <: AbstractIntegrator,
    T3 <: Int64,
}
    dt::T1
    tmin::T1
    tmax::T1
    method::T2
    abstol::T1
    reltol::T1
    time::T1
    step::T3
```
This structure contains the integration parameters for the simulation.
"""
mutable struct IntegrationP{
    T1 <: Float64,
    T2 <: AbstractIntegrator,
    T3 <: Int64,
}
    dt::T1
    tmin::T1
    tmax::T1
    method::T2
    abstol::T1
    reltol::T1
    time::T1
    step::T3

    function IntegrationP(;
        dt,
        tmin,
        tmax,
        method,
        abstol = 1e-6,
        reltol = 1e-6,
        time = 0.0,
        step = 0,
    )
        new{typeof(dt), typeof(method), typeof(step)}(
            dt,
            tmin,
            tmax,
            method,
            abstol,
            reltol,
            time,
            step,
        )
    end
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
