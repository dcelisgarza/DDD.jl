function remeshSurfaceNetwork!()
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
        if getNodeType(label[i]) == ext && connectivity[1, i] == 2
            # This is where node i appears in connectivity.
            node1 = connectivity[2, i] # Link where node i appears first.
            linkCol1 = 3 - connectivity[3, i] # Column of links where it appears.
            node2 = connectivity[4, i] # Link where node i appears second.
            linkCol2 = 3 - connectivity[5, i] # Column of links where it appears.

            linkNode1 = links[node1, linkCol1] # First node connected to target node.
            linkNode2 = links[node2, linkCol2] # Second node connected to target node.

            # Only if both nodes are virtual.
            if getNodeType(label[linkNode1]) == getNodeType(label[linkNode2]) == ext
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
