"""
```
deriv!(
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    mesh::AbstractMesh,
    forceDisplacement::ForceDisplacement,
    network::DislocationNetwork,
)
```
Computes the nodal velocities of a network.
"""
function deriv!(
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    mesh::AbstractMesh,
    forceDisplacement::ForceDisplacement,
    network::DislocationNetwork,
)
    calcSegForce!(dlnParams, matParams, mesh, forceDisplacement, network)
    dlnMobility!(dlnParams, matParams, network)

    links = network.links
    coord = network.coord
    nodeVel = network.nodeVel
    label = network.label
    numNode = network.numNode[1]

    # Fixed nodes do not move, internal fixed, surface fixed, external.
    fixedNode = (intFixDln, srfFixDln, extDln)
    idxFixed = findall(x -> x ∈ fixedNode, label)
    !isnothing(idxFixed) ? nodeVel[:, idxFixed] .= 0 : nothing

    # Make surface velocity zero along the surface normal.
    idxSurf = findall(x -> x == srfMobDln, label)
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
        extIdx = findall(x -> x == extDln, conLabel)

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
            nodeVel[:, node] -= @views (nodeVel[:, node] ⋅ conNodes) * conNodes
        end
    end

    return nothing
end

"""
```
integrate!(
    intParams::IntegrationParameters{T1,T2,T3} where {T1 <: AdaptiveEulerTrapezoid,T2,T3},
    intVars::IntegrationTime,
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    mesh::AbstractMesh,
    forceDisplacement::ForceDisplacement,
    network::DislocationNetwork,
)
```
Integrates nodal velocities using a time-adaptive Euler-Trapezoid method.
"""
function integrate!(
    intParams::IntegrationParameters{
        T1,
        T2,
        T3,
    } where {T1 <: AdaptiveEulerTrapezoid, T2, T3},
    intVars::IntegrationTime,
    dlnParams::DislocationParameters,
    matParams::MaterialParameters,
    mesh::AbstractMesh,
    forceDisplacement::ForceDisplacement,
    network::DislocationNetwork,
)
    numNode = network.numNode[1]
    numNode == 0 && return network

    # Calculate current velocity.
    deriv!(dlnParams, matParams, mesh, forceDisplacement, network)

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

    dt = max(intVars.dt, dtmin)
    time = intVars.time
    step = intVars.step

    dtOld = zero(typeof(dt))
    dtOldGood::Bool = false
    convergent::Bool = false

    counter = 1
    while !convergent
        # Advance coordinates with forward euler.
        network.coord[:, idx] = initCoord + initVel * dt

        # Calculate new velocity from new coords.
        deriv!(dlnParams, matParams, mesh, forceDisplacement, network)
        coord = network.coord
        nodeVel = network.nodeVel

        # Calculate the distance moved.
        distance = @views coord[:, idx] - initCoord
        maxDist = maximum(abs.(distance))

        # Calculate the error with Euler trapezoid method.
        err = @views distance - (nodeVel[:, idx] + initVel) / 2 * dt
        maxErr = maximum(abs.(err))

        if dt <= dtmin
            counter = maxiter + 1
            # Check if errors are under the tolerances.
        elseif maxDist < abstol && maxErr < reltol
            dtOld = dt
            factor =
                maxchange *
                (1 / (1 + (maxchange^exponent - 1) * (maxErr / reltol)))^(1 / exponent)
            # If the errors are under the tolerances, increase the time step.
            dtOldGood = true
            dt = min(dt * factor, dtmax)
            counter += 1
            # If the errors are over the tolerances, check if the previous time step is good.
        elseif dtOldGood
            dt = dtOld
            counter = maxiter
            # If it wasn't, make the timestep smaller.
        else
            dt = dt / 2
        end

        if counter > maxiter || dt == dtmax
            convergent = true
        end
    end

    intVars = IntegrationTime(dt, time + dt, step + 1)

    return intVars
end
