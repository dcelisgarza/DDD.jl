#=
function deriv!(
    dlnParams::T1,
    matParams::T2,
    network::T3,
) where {T1 <: DislocationParameters, T2 <: MaterialParameters, T3 <: DislocationNetwork}

    calcSegForce!(dlnParams, matParams, network)
    dlnMobility!(dlnParams, matParams, network)

    links = network.links
    coord = network.coord
    nodeVel = network.nodeVel
    label = network.label
    numNode = network.numNode[1]

    # Fixed nodes do not move, internal fixed, surface fixed, external.
    idxFixed = findall(x -> x == 2 || x == 4 || x == 5, label)
    !isnothing(idxFixed) ? nodeVel[:, idxFixed] .= [0, 0, 0] : nothing

    # Make surface velocity zero along the surface normal.
    idxSurf = findall(x -> x == 3, label)
    for node in idxSurf
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
                normalize!(conNodes)
            end
            # Vector rejection to remove velocity in the direction of the surface normal.
            nodeVel[:, node] -= (nodeVel[:, node] ⋅ conNodes) * conNodes
        end
    end

    return nothing
end

function integrate!(
    intMethod::CustomTrapezoid,
    intParams::T1,
    intVars::T2,
    dlnParams::T3,
    matParams::T4,
    network::T5,
) where {
    T1 <: IntegrationParameters,
    T2 <: IntegrationTime,
    T3 <: DislocationParameters,
    T4 <: MaterialParameters,
    T5 <: DislocationNetwork,
}

    numNode = network.numNode[1]
    numNode == 0 && return network

    # Calculate current velocity.
    deriv!(dlnParams, matParams, network)

    # Store current position and velocity.
    idx = 1:numNode
    initCoord = network.coord[:, idx]
    initVel = network.nodeVel[:, idx]

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
        network.coord[:, idx] = initCoord + initVel * dt

        # Calculate new velocity from new coords.
        deriv!(dlnParams, matParams, network)
        coord = network.coord
        nodeVel = network.nodeVel

        # Calculate the distance moved.
        distance = @views coord[:, idx] - initCoord
        maxDist = maximum(abs.(distance))

        # Calculate the error with Euler trapezoid method.
        err = @views distance - (nodeVel[:, idx] + initVel) / 2 * dt
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
                network.coord[:, idx] = initCoord + initVel * dtOld
                deriv!(dlnParams, matParams, network)
                break
            else
                # If it wasn't, make the timestep smaller.
                dt = min(dtNew, dt / 2)
            end
        end
        # Break if dt is less than the minimum allowed.
        dt <= dtmin && break
    end

    intVars = IntegrationTime(dt, time + dt, step + 1)

    return intVars
end
=#
