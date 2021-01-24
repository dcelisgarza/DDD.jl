function findIntersectVolume(mesh::AbstractMesh, l, l0, tmpArr)
    vertices = mesh.vertices
    faceMidPt = mesh.faceMidPt
    faceNorm = mesh.faceNorm
    numFaces = size(faceMidPt, 2)
    elemT = eltype(l)

    distMin::elemT = Inf
    face = 0
    # Declare new node.
    newNode = SVector{3,elemT}(0, 0, 0)
    for i in 1:numFaces
        oldMin = distMin
        intersect = linePlaneIntersect(faceNorm[:, i], faceMidPt[:, i], l, l0)
        if isnothing(intersect)
            continue
        else
            tmpArr[:] = intersect[:]
            if tmpArr ∉ vertices
                continue
            end
        end
        distMin = min(sum((intersect - l0).^2), oldMin)
        if distMin < oldMin
            # The new node will be placed at the point where the line direction intersects the nearest face, which is face i.
            newNode = intersect
            face = i
        end
    end
    return distMin, newNode, face
end

function remeshSurfaceNetwork!(
    mesh::AbstractMesh, 
    network::DislocationNetwork
)
    scale = mesh.scale
    vertices = mesh.vertices
    faceMidPt = mesh.faceMidPt
    faceNorm = mesh.faceNorm
    numFaces = size(faceMidPt, 2)
    
    label = network.label
    links = network.links
    coord = network.coord
    elemT = eltype(coord)
    numNode = network.numNode[1]
    connectivity = network.connectivity

    # Findall internal nodes.
    idx = findall(x -> x == intMob, label)
    numInt = length(idx)
    # Find the location of the nodes that are newly outside the domain, P.
    nodeCoord = zeros(elemT, 3)
    @inbounds @simd for i in 1:numInt
        idxi = idx[i]
        for j in 1:3
            nodeCoord[j] = coord[j, idxi]
        end
        if nodeCoord ∉ vertices
            label[i] = tmp
        end
    end

    for node1 in 1:numNode
        # We only want to check connections of mobile internal or external nodes.
        label[node1] == intMob || label[node1] == ext ? nothing : continue
        # Find the number of connections to the node.
        numCon = connectivity[1, node1]
        # Loop through all connections of node1.
        for j in 1:numCon
            idx = 2 * j
            # Link where node1 appears for its connection j.
            linkId = connectivity[idx, node1]
            # Column where node1 appears.
            colLink = connectivity[idx + 1, node1]
            # Neigbour node from j'th link node1 appears.
            node2 = links[3 - colLink, linkId]

            if label[node2] == tmp
                network = makeSurfaceNode!(mesh, network, node1, node2, j)
                getSegmentIdx!(network)
            end
        end
    end

    label = network.label
    coord = network.coord
    nodeVel = network.nodeVel
    numNode = network.numNode[1]
    intersectArr = zeros(elemT, 3)

    for i in 1:numNode
        # We only want to change the temporary nodes.
        label[i] != tmp ? continue : nothing

        # Assume nodes moved with linear velocity and crossed the surface. We use this velocity as the vector to move them back.
        vel = SVector{3,elemT}(nodeVel[1, i], nodeVel[2, i], nodeVel[3, i])
        velN = norm(vel)
        iszero(velN) ? throw(ErrorException("norm of the nodal velocity must be greater than zero")) : nothing
        vel = vel / velN
        l0 = SVector{3,elemT}(coord[1, i], coord[2, i], coord[3, i])
        
        distMin, newNode, face = findIntersectVolume(mesh, vel, l0, intersectArr)

        # If there is no intersect find the nearest plane to project it out of.
        if isinf(distMin)
            # D = |(x0 + p0) ⋅ n/||n||, where n := plane normal, p0 a point on the plane, p0 ⋅ n = d, from the plane equation ax + by + cz = d, x0 is a point in space.
            distances = ((l0[1] .+ faceMidPt[1, :]) .* faceNorm[1, :] .+ (l0[2] .+ faceMidPt[2, :]) .* faceNorm[2, :] .+ (l0[3] .+ faceMidPt[3, :]) .* faceNorm[3, :]).^2
            distMin, face = findmin(distances)
            newNode = l0 + faceNorm[:, face] * distMin
        end
        # If there is an intersect, move the node to the intersect and we will project it out of the face it intersected.
        coord[:, i] = newNode

        for j in 1:3
            coord[j, i] = coord[j, i] + faceNorm[j, face] * scale[j]
        end
        
        label[i] = ext
    end

    # Find connected surface nodes and merge them.
    links = network.links
    connectivity = network.connectivity
    node1 = 1
    while node1 < numNode
        if label[node1] != srfMob || label[node1] != srfFix
            node1 += 1
            continue
        end

        # Find number of connections to node.
        check::Bool = false
        numCon = connectivity[1, node1]
        # Loop through the number of connections.
        for j in 1:numCon
            idx = 2 * j
            linkId = connectivity[idx, node1]
            colLink = connectivity[idx + 1, node1]
            node2 = links[3 - colLink, linkId]
            # Check if the node connected is a surface node aso we can merge node1 into node2.
            if label[node2] == srfMob || label[node2] == srfFix
                missing, network = mergeNode!(network, node2, node1)
                getSegmentIdx!(network)
                links = network.links
                label = network.label
                numNode = network.numNode[1]
                connectivity = network.connectivity
                node1 = 1
                check = true
                break
            end
        end
        check ? continue : nothing
        node1 += 1
    end

    # Find surface nodes that are only connected to virtual nodes and project them to be external too.
    numNode = network.numNode[1]
    label = network.label
    coord = network.coord
    connectivity = network.connectivity
    for node1 in 1:numNode
        label[node1] == srfMob || label[node1] == srfFix ? continue : nothing

        # Number of external connections.
        numExtCon = 0
        numCon = connectivity[1, node1]
        for j in 1:numCon
            linkId = connectivity[idx, node1]
            colLink = connectivity[idx + 1, node1]
            node2 = links[3 - colLink, linkId]

            if label[node2] == ext
                numExtCon += 1
            end
        end

        # If the surface node is only connected to external nodes move away from the surface.
        if numCon == numExtCon
            # Check if node1 intersects a surface by projecting along the surface normal. Do this for all surfaces and find the minimum distance to a surface. Save the surface node1 intercepts and the coordinate where it does so.
            l0 = SVector{3,elemT}(coord[1, node1], coord[2, node1], coord[3, node1])
            face = 0
            distMin::elemT = Inf
            newNode = SVector{3,elemT}(0, 0, 0)
            for j in 1:numFaces
                oldDist = distMin
                distMin, newNodeMin, faceMin = findIntersectVolume(mesh, faceNorm[:, j], l0, intersectArr)
                distMin = min(distMin, oldDist)
                if distMin < oldDist
                    face = faceMin
                    newNode = newNodeMin
                end
            end

            # If there is no intersect with a surface find the nearest plane to project out of.
            if isinf(distMin)
                # D = |(x0 + p0) ⋅ n/||n||, where n := plane normal, p0 a point on the plane, p0 ⋅ n = d, from the plane equation ax + by + cz = d, x0 is a point in space.
                distances = ((l0[1] .+ faceMidPt[1, :]) .* faceNorm[1, :] .+ (l0[2] .+ faceMidPt[2, :]) .* faceNorm[2, :] .+ (l0[3] .+ faceMidPt[3, :]) .* faceNorm[3, :]).^2
                distMin, face = findmin(distances)
                # Point where the node intersects with the plane.
                newNode = l0 + faceNorm[:, face] * distMin
            end

            # Move the node to the intersect and we will project it out of the face it intersected.
            coord[:, node1] = newNode

            for j in 1:3
                coord[j, node1] = coord[j, node1] + faceNorm[j, face] * scale[j]
            end

            label[node1] = ext
        end
    end

    # Ensure surface nodes are only connected to one virtual node.
    node1 = 1
    while node1 < numNode
        # Only check surface nodes.
        if label[node1] != srfMob || label[node1] != srfFix
            node1 += 1
            continue
        end

        # Only check surface nodes connected to more than 2 node.
        if connectivity[1, node1] <= 2
            node1 += 1
            continue
        end

        # Go through the connections checking if there is more than one connection to an external node.
        for j in 1:connectivity[1, node1] - 1
            idx1 = 2 * j
            linkId1 = connectivity[idx1, node1]
            colLink1 = connectivity[idx1 + 1, node1]
            node2 = links[3 - colLink1, linkId1]

            # If the first connection is not external go to the next iteration.
            label[node2] != ext ? continue : nothing

            for k in j + 1:connectivity[1, node1]
                idx2 = 2 * k
                linkId2 = connectivity[idx2, node1]
                colLink2 = connectivity[idx2 + 1, node1]
                node3 = links[3 - colLink2, linkId2]
                # If the next connection is not external go to the next iteration.
                label[node2] != ext ? continue : nothing

                missing, network = mergeNode!(network, node3, node1)
                getSegmentIdx!(network)
                links = network.links
                label = network.label
                numNode = network.numNode[1]
                connectivity = network.connectivity
                # We merged a node, so we want to repeat the outer loop.
                node1 = 1
                break
            end
            # If we get here it means the inner loop found one external connection to an external node but the inner loop didn't find another, so we can cut the outer loop as it would just retread the steps of the inner one.
            break
        end
        node1 += 1
    end

    # Make surface nodes if necessary.
    numNode = network.numNode[1]
    links = network.links
    connectivity = network.connectivity
    for node1 in 1:numNode
        label[node1] != ext ? continue : nothing

        numCon = connectivity[1, node1]
        for j in 1:numCon
            idx = 2 * j
            # Link where node1 appears for its connection j.
            linkId = connectivity[idx, node1]
            # Column where node1 appears.
            colLink = connectivity[idx + 1, node1]
            # Neigbour node from j'th link node1 appears.
            node2 = links[3 - colLink, linkId]

            if label[node2] == intMob || label[node2] == intFix
                network = makeSurfaceNode!(mesh, network, node1, node2, j)
                getSegmentIdx!(network)
            end
        end    
    end

    return network
end

function makeSurfaceNode!(
    mesh::AbstractMesh, 
    network::DislocationNetwork,
    node1::Integer,
    node2::Integer,
    idx::Integer
)
    vertices = mesh.vertices
    faceMidPt = mesh.faceMidPt
    faceNorm = mesh.faceNorm
    numFaces = size(faceMidPt, 2)
    coord = network.coord
    elemT = eltype(coord)

    coordNode1 = SVector{3,elemT}(coord[1, node1], coord[2, node1], coord[3, node1])
    coordNode2 = SVector{3,elemT}(coord[1, node2], coord[2, node2], coord[3, node2])
    t = coordNode2 - coordNode1
    intersectArr = zeros(elemT, 3)
    distMin, newNode, missing = findIntersectVolume(mesh, t, coordNode1, intersectArr)

    isinf(distMin) && return network
    # We set the surface node velocity to zero. It will be calculated later.
    network = splitNode!(network, node1, idx, newNode, SVector{3,elemT}(0, 0, 0))
    # Set the label of the new node to be a surface node.
    network.label[network.numNode[1]] = srfMob

    return network
end

"""
```
coarsenVirtualNetwork!(
    dlnParams::T1,
    network::T2,
) where {T1 <: DislocationParameters, T2 <: DislocationNetwork}
```
Check whether virtual nodes can be eliminated based on:
1) If they are not connected to any surface nodes
2) If they are not due to an angle change in the simulated volume surface

Bruce Bromage, Github @brucebromage
Michromechanical Testing Group
Department of Materials, University of Oxford
bruce.bromage@materials.ox.ac.uk
May 2017

Adapted Jan 2021 Daniel Celis Garza, Github @dcelisgarza
"""
    function coarsenVirtualNetwork!(
    dlnParams::T1,
    network::T2,
) where {T1 <: DislocationParameters,T2 <: DislocationNetwork}

        critLen = dlnParams.slipStepCritLen
        critArea = dlnParams.slipStepCritArea
        label = network.label
        links = network.links
        coord = network.coord
        numNode = network.numNode[1]
        connectivity = network.connectivity
        elemT = eltype(coord)

        i = 1
        while i <= numNode
        # Only find virtual nodes with two connections.
            if label[i] == ext && connectivity[1, i] == 2
            # This is where node i appears in connectivity.
                node1 = connectivity[2, i] # Link where node i appears first.
                linkCol1 = 3 - connectivity[3, i] # Column of links where it appears.
                node2 = connectivity[4, i] # Link where node i appears second.
                linkCol2 = 3 - connectivity[5, i] # Column of links where it appears.

                linkNode1 = links[node1, linkCol1] # First node connected to target node.
                linkNode2 = links[node2, linkCol2] # Second node connected to target node.

            # Only if both nodes are virtual.
                if label[linkNode1] == label[linkNode2] == ext
                # Coordinate of node i.
                    iCoord = SVector{3,elemT}(coord[1, i], coord[2, i], coord[3, i])

                # Vector of link 1.
                    coordVec1 =
                    SVector{3,elemT}(
                        coord[1, linkNode1],
                        coord[2, linkNode1],
                        coord[3, linkNode1],
                    ) - iCoord
                    normVec1 = norm(coordVec1)

                # Vector of link 2.
                    coordVec2 =
                    SVector{3,elemT}(
                        coord[1, linkNode2],
                        coord[2, linkNode2],
                        coord[3, linkNode2],
                    ) - iCoord
                    normVec2 = norm(coordVec2)

                # Angle between vectors.
                    θ = acos(norm(coordVec1 ⋅ coordVec2) / (normVec1 * normVec2))

                # Area formed by the triangle formed by link 1 and link 2.
                    area = 0.5 * normVec1 * normVec2 * sin(θ)

                # If the length of link 1 and the angle change are below the critical size, merge node i to the first connected node.
                    if normVec1 < critLen && area < critArea
                        nothing, network = mergeNode!(network, linkNode1, i)
                        getSegmentIdx!(network)
                        links = network.links
                        coord = network.coord
                        label = network.label
                        numNode = network.numNode[1]
                        connectivity = network.connectivity
                    # If the length of link 1 and the angle change are below the critical size, merge node i to the second connected node.
                    elseif normVec2 < critLen && area < critArea
                        nothing, network = mergeNode!(network, linkNode2, i)
                        getSegmentIdx!(network)
                        links = network.links
                        coord = network.coord
                        label = network.label
                        numNode = network.numNode[1]
                        connectivity = network.connectivity
                    else
                        i += 1
                    end
                else
                    i += 1
                end
            else
                i += 1
            end

        end

    end
