"""
```
findIntersectVolume(mesh::AbstractMesh, l, l0, tmpArr)
```
Find the shortest intersecting distance between a vector `l` passing through the point `l0` and an [`AbstractMesh`](@ref).
"""
function findIntersectVolume(mesh::AbstractMesh, l, l0)
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
        if isnothing(intersect) || any(isinf.(intersect))
            continue
        else
            if Array(intersect) ∉ vertices
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

"""
```
makeSurfaceNode!(mesh::AbstractMesh,  network::DislocationNetwork, node1, node2, idx)
```
Creates a surface node between `node1` and `node2`, using connection `idx` of `node1` of a [`DislocationNetwork`](@ref) on the surface of an [`AbstractMesh`](@ref)
"""
function makeSurfaceNode!(mesh::AbstractMesh,  network::DislocationNetwork, node1, node2, idx)
    vertices = mesh.vertices
    faceMidPt = mesh.faceMidPt
    faceNorm = mesh.faceNorm
    numFaces = size(faceMidPt, 2)
    coord = network.coord
    elemT = eltype(coord)

    coordNode1 = SVector{3,elemT}(coord[1, node1], coord[2, node1], coord[3, node1])
    coordNode2 = SVector{3,elemT}(coord[1, node2], coord[2, node2], coord[3, node2])
    t = coordNode2 - coordNode1
    distMin, newNode, missing = findIntersectVolume(mesh, t, coordNode1)

    isinf(distMin) && return network
    # We set the surface node velocity to zero. It will be calculated later.
    network = splitNode!(network, node1, idx, newNode, SVector{3,elemT}(0, 0, 0))
    # Set the label of the new node to be a surface node.
    network.label[network.numNode[1]] = srfMobDln

    return network
end

"""
```
remeshSurfaceNetwork!(mesh::AbstractMesh, boundaries::Boundaries, network::DislocationNetwork)
```
Remeshes a [`DislocationNetwork`](@ref)'s nodes on the surface of an [`AbstractMesh`](@ref).
"""
function remeshSurfaceNetwork!(mesh::AbstractMesh, boundaries::Boundaries, network::DislocationNetwork)
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
    intNodes = (intMobDln, intFixDln)
    srfNodes = (srfMobDln, srfFixDln)

    # Find internal nodes that are newly exited and mark them as such.
    for i in 1:numNode
        label[i] ∉ intNodes && continue
        coord[:, i] ∈ vertices && continue
        label[i] = tmpDln
    end

    # Generate surface nodes by finding which internal nodes are connected to newly exited nodes.
    for node1 in 1:numNode
        # Only look for nodes that are internal.
        label[node1] ∉ intNodes && continue
        # Find the number of connections to the node.
        numCon = connectivity[1, node1]
        # Loop through all connections of node1 to see if it is connected to a newly exited node.
        for j in 1:numCon
            missing, missing, node2 = findConnectedNode(network, node1, j) # j connection.

            # If node1 is connected to a newly exited node 2, create surface node.
            label[node2] == tmpDln ? network = makeSurfaceNode!(mesh, network, node1, node2, j) : nothing
        end
    end

    label = network.label
    coord = network.coord
    nodeVel = network.nodeVel
    numNode = network.numNode[1]

    # Project newly exited nodes out to pseudo-infinity.
    for i in 1:numNode
        # We only want to change nodes flagged as having just exited.
        label[i] != tmpDln && continue

        # Assume nodes moved with linear velocity and crossed the surface. We use this velocity as the vector to move them back onto the surface and project from there.
        vel = SVector{3,elemT}(nodeVel[1, i], nodeVel[2, i], nodeVel[3, i])
        velN = norm(vel)
        iszero(velN) ? throw(ErrorException("norm of the nodal velocity must be greater than zero")) : nothing
        vel = vel / velN
        l0 = SVector{3,elemT}(coord[1, i], coord[2, i], coord[3, i])
        
        distMin, newNode, face = findIntersectVolume(mesh, vel, l0)

        # If there is no intersect, find the nearest plane to project it out of.
        if isinf(distMin)
            # D = (x0 + p0) ⋅ n/||n||, where n := plane normal, p0 a point on the plane, p0 ⋅ n = d, from the plane equation ax + by + cz = d, x0 is a point in space.
            distances = ((l0[1] .+ faceMidPt[1, :]) .* faceNorm[1, :] .+ (l0[2] .+ faceMidPt[2, :]) .* faceNorm[2, :] .+ (l0[3] .+ faceMidPt[3, :]) .* faceNorm[3, :]).^2
            distMin, face = findmin(distances)
            newNode = l0 + faceNorm[:, face] * distMin
        end
        # If there is an intersect, move the node to the intersect and we will project it out of the face it intersected.
        coord[:, i] = newNode

        if face ∉ boundaries.noExit
            for j in 1:3
                coord[j, i] = coord[j, i] + faceNorm[j, face] * scale[j]
            end
            
            label[i] = extDln
        else
            label[i] = srfFixDln
        end
    end

    # Find connected surface nodes and merge them.
    links = network.links
    connectivity = network.connectivity
    node1 = 1
    while node1 < numNode
            if label[node1] ∉ srfNodes
        node1 += 1
            continue
        end
        
        # Find number of connections to node.
        check = false
        numCon = connectivity[1, node1]
        # Loop through the number of connections.
        for j in 1:numCon
            missing, missing, node2 = findConnectedNode(network, node1, j) # j connection.
            
            # Check if the node connected is a surface node aso we can merge node1 into node2.
            if label[node2] ∈ srfNodes
                missing, network = mergeNode!(network, node2, node1)
                links = network.links
                label = network.label
                numNode = network.numNode[1]
                connectivity = network.connectivity
                node1 = 1
                check = true
                break
            end
        end
        check ? continue : node1 += 1
    end

    # Find surface nodes that are only connected to virtual nodes and project them to be external too.
    numNode = network.numNode[1]
    label = network.label
    coord = network.coord
    connectivity = network.connectivity
    for node1 in 1:numNode
        label[node1] ∉ srfNodes && continue
            
        # Number of external connections.
        numExtCon = 0
        numCon = connectivity[1, node1]
        for j in 1:numCon
            missing, missing, node2 = findConnectedNode(network, node1, j) # j connection.
            label[node2] == extDln ? numExtCon += 1 : nothing
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
            distMin, newNodeMin, faceMin = findIntersectVolume(mesh, faceNorm[:, j], l0)
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

            label[node1] = extDln
        end
    end

    # Ensure surface nodes are only connected to one virtual node.
    node1 = 1
    while node1 < numNode
        # Only check surface nodes with more than two connections.
        if label[node1] ∉ srfNodes || connectivity[1, node1] <= 2
            node1 += 1
            continue
        end

        # Go through the connections checking if there is more than one connection to an external node.
        for j in 1:connectivity[1, node1] - 1
            missing, missing, node2 = findConnectedNode(network, node1, j) # j connection.

            # If the first connection is not external go to the next iteration.
            label[node2] != extDln && continue

            for k in j + 1:connectivity[1, node1]
                missing, missing, node3 = findConnectedNode(network, node1, k) # k connection.
                
                # If the next connection is not external go to the next iteration.
                label[node3] != extDln && continue

                missing, network = mergeNode!(network, node3, node1)
                links = network.links
                label = network.label
                numNode = network.numNode[1]
                connectivity = network.connectivity
                # We merged a node, so we want to repeat the outer loop.
                node1 = 0
                break
            end
            # If we get here it means the inner loop found one external connection to an external node but the inner loop didn't find another, so we can cut the outer loop as it would just retread the steps of the inner one.
            break
        end
        node1 += 1
    end

    # Make surface nodes if necessary. Only want to check external nodes.
    numNode = network.numNode[1]
    links = network.links
    connectivity = network.connectivity
    for node1 in 1:numNode
        label[node1] != extDln && continue

        numCon = connectivity[1, node1]
        for j in 1:numCon
            missing, missing, node2 = findConnectedNode(network, node1, j) # j connection.
            # If the node2 internal, create a surface node where the segment intersects the mesh.
            label[node2] ∈ intNodes ? network = makeSurfaceNode!(mesh, network, node1, node2, j) : nothing
        end    
    end

    return network
end

"""
```
coarsenVirtualNetwork!(dlnParams::DislocationParameters, network::DislocationNetwork)
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
function coarsenVirtualNetwork!(dlnParams::DislocationParameters, network::DislocationNetwork)
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
        if label[i] == extDln && connectivity[1, i] == 2
            missing, missing, linkNode1 = findConnectedNode(network, i, 1) # 1st connection.
            missing, missing, linkNode2 = findConnectedNode(network, i, 2) # 2nd connection.

            # Only if both nodes are virtual.
            if label[linkNode1] == label[linkNode2] == extDln
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
