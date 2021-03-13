## Topology functions
function findConnectedNode(network::DislocationNetwork, node, connection)
    links = network.links
    connectivity = network.connectivity

    link = connectivity[2 * connection, node] # Find link that corresponds to node's connection.
    colLink = connectivity[2 * connection .+ 1, node] # Find whether node corresponds to the first or second entry of the link corresponding to connection.
    oppColLink = 3 .- colLink # The other node corresponds to the other entry in the link.
    connectedNode = links[oppColLink, link] # Connected node.

    return link, oppColLink, connectedNode
end
include("DislocationTopologyRemove.jl")
include("DislocationTopologyAdd.jl")
include("DislocationTopologySpecial.jl")