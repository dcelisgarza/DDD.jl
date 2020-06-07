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
function IntegrationP(; tmin, tmax, method, abstol = 1e-6, reltol = 1e-6, dt0 = 1.0)
    return IntegrationP(tmin, tmax, method, abstol, reltol, dt0)
end
mutable struct IntegrationVar{T1, T2}
    dt::T1
    time::T1
    step::T2
end
function IntegrationVar(; dt = 0.0, time = 0.0, step = 0)
    return IntegrationVar(dt, time, step)
end

function deriv!(dlnParams::T1, matParams::T2, network::T3; parallel = false)

    calcSegForce!(dlnParams, matParams, network; parallel = parallel)
    dlnMobility!(dlnParams, matParams, network)

    links = network.links
    coord = network.coord
    nodeVel = network.nodeVel
    label = network.label
    numNode = network.numNode

    # Fixed nodes do not move, internal fixed, surface fixed, external.
    idxFixed = findall(x -> x == 2 || x == 4 || x == 5, label)
    !isnothing(idxFixed) ? nodeVel[:, idxFixed] = [0, 0, 0] : nothing

    # Make surface velocity zero along the surface normal.
    idxSurf = findall(x -> x == 3, label)
    @inbounds for node in idxSurf
        # Find the links where node appears.
        nodeLink1 = links[1, :] .== node
        nodeLink2 = links[2, :] .== node

        # Find the node opposite to node in all its links.
        nodeOppLink1 = links[2, nodeLink1]
        nodeOppLink2 = links[1, nodeLink2]

        # Labels of the nodes connected to node.
        conLabel = label[(nodeOppLink1, nodeOppLink2)]
        # Find only those that are external.
        extIdx = findall(x -> x == 5, conLabel)

        # Correct only if node is connected to external nodes.
        if !isnothing(extIdx)
            numExt = length(extIdx)
            conNodes = coord[:, extIdx]
            # If there is more than one external connection, calculate the mean normal.
            if numExt > 1
                conNodes = mean(conNodes, dims = 2)
                normalise!(conNodes)
            end
            # Vector rejection to remove velocity in the direction of the surface normal.
            nodeVel[:, node] -= (nodeVel[:, node] ⋅ conNodes) * conNodes
        end
    end

    return network
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
