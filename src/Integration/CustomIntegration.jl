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
struct IntegrationP{T1, T2, T3}
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
function IntegrationP(;
    method,
    tmin,
    tmax,
    dtmin = 1e-3,
    dtmax = Inf,
    abstol = 1e-6,
    reltol = 1e-6,
    maxchange = 1.2,
    exponent = 20,
    maxiter = 10,
)
    return IntegrationP(
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
mutable struct IntegrationVar{T1, T2}
    dt::T1
    time::T1
    step::T2
end
function IntegrationVar(; dt = 0.0, time = 0.0, step = 0)
    return IntegrationVar(dt, time, step)
end

function deriv!(
    dlnParams::T1,
    matParams::T2,
    network::T3;
    parallel::Bool = false,
) where {T1 <: DislocationP, T2 <: MaterialP, T3 <: DislocationNetwork}

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

function integrate!(
    intParams::T1,
    intVars::T2,
    dlnParams::T3,
    matParams::T4,
    network::T5;
    parallel::Bool = false,
) where {
    T1 <: IntegrationP,
    T2 <: IntegrationVar,
    T3 <: DislocationP,
    T4 <: MaterialP,
    T5 <: DislocationNetwork,
}
    return integrate!(
        intParams.method,
        intParams,
        intVars,
        dlnParams,
        matParams,
        network;
        parallel = false,
    )
    # Calculate current velocity.
    deriv!(dlnParams, matParams, network; parallel = parallel)
    # Store current position and velocity.
    initCoord = copy(network.coord)
    initVel = copy(network.nodeVel)
end

function integrate!(
    intMethod::CustomTrapezoid,
    intParams::T1,
    intVars::T2,
    dlnParams::T3,
    matParams::T4,
    network::T5;
    parallel::Bool = false,
) where {
    T1 <: IntegrationP,
    T2 <: IntegrationVar,
    T3 <: DislocationP,
    T4 <: MaterialP,
    T5 <: DislocationNetwork,
}

    numNode = network.numNode
    numNode == 0 && return network

    # Calculate current velocity.
    deriv!(dlnParams, matParams, network; parallel = parallel)

    # Store current position and velocity.
    idx = 1:numNode
    initCoord = copy(network.coord[:, idx])
    initVel = copy(network.nodeVel[:, idx])

    # Preallocate arrays.
    sizeCoord = size(initCoord)
    distance = zeros(sizeCoord)
    err = zeros(sizeCoord)

    dtmin = intParams.dtmin
    dtmax = intParams.dtmax
    abstol = intParams.abstol
    reltol = intParams.reltol
    maxchange = intParams.maxchange
    exponent = intParams.exponent
    maxiter = intParams.maxiter

    dt = intVars.dt
    dtOld = zero(typeof(dt))
    dtOldGood::Bool = false

    count = 0
    for i in 1:maxiter
        # Advance coordinates with forward euler.
        network.coord[:, idx] .= initCoord + initVel * dt

        # Calculate new velocity from new coords.
        deriv!(dlnParams, matParams, network; parallel = parallel)
        coord = network.coord
        nodeVel = network.nodeVel

        # Calculate the distance moved.
        distance .= coord[:, idx] - initCoord
        maxDist = maximum(abs.(distance))

        # Calculate the error with Euler trapezoid method.
        err .= distance - (nodeVel[:, idx] + initVel) / 2 * dt
        maxErr = maximum(abs.(err))

        # Tentative new timestep.
        factor =
            maxchange *
            (1 / (1 + (maxchange^exponent - 1) * (maxErr / reltol)))^(1 / exponent)
        dtNew = dt * factor

        # Check if errors are under the tolerances.
        if maxDist < abstol && maxErr < reltol
            # If the errors are under the tolerances, increase the time step.
            dtOldGood = true
            dtOld = dt
            dt = min(dtNew, dtmax)
        else
            # If the errors are over the tolerances, check if the previous time step is good.
            if dtOldGood
                # If it was, use the previous time step.
                network.coord[:, idx] .= initCoord + initVel * dtOld
                deriv!(dlnParams, matParams, network; parallel = parallel)
                break
            else
                # If it wasn't, make the timestep smaller.
                dt = min(dtNew, dt / 2)
            end
        end
        # Break if dt is less than the minimum allowed.
        dt <= dtmin && break
    end

    intVars.dt = dt
    intVars.time += dt
    intVars.step += 1

    return intVars, network
end
