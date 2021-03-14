## Topology functions
"""
```
findConnectedNode(network::DislocationNetwork, node, connection)
```
Use a `node`'s `connection` to find the other node in the link. Returns a tuple with the link corresponding to the `connection` of `node`, the row of links the connected node appears, and the connected node. 

## Examples
We want to find the node connected to `node` 5 via `connection` 3.

```julia
julia> link, rowColLink, connectedNode = findConnectedNode(network, 5, 10)
julia> link
10
julia> oppRowLink
1
julia> connectedNode
25
```
Means that the link corresponding to `connection` 3 of `node` 5 is `network.links[:, 10]`; the node connected to `node` 5 is found on `network.links[oppRowLink, 10]` and is node 25, so the node we're looking for is `network.links[oppRowLink, 10] == 25` and the node whose connection we are looking for is `network.links[3-oppRowLink, 10] == 5`. In this case `oppRowLink == 1`, so through its 3rd connection, node 5 is the second node on link 10, where it is connected to node 25.
"""
function findConnectedNode(network::DislocationNetwork, node, connection)
    links = network.links
    connectivity = network.connectivity

    link = connectivity[2 * connection, node] # Find link that corresponds to node's connection.
    colLink = connectivity[2 * connection .+ 1, node] # Find whether node corresponds to the first or second entry of the link corresponding to connection.
    rowColLink = 3 .- colLink # The other node corresponds to the other entry in the link.
    connectedNode = links[rowColLink, link] # Connected node.

    return link, rowColLink, connectedNode
end
include("DislocationTopologyRemove.jl")
include("DislocationTopologyAdd.jl")
include("DislocationTopologySpecial.jl")
include("DislocationTopologyCollision.jl")