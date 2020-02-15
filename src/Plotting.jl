function plotNodes(network::DislocationNetwork, args...; kw...)
    idx = idxLabel(network, -1; condition = !=)
    coord = network.coord
    fig = plot()
    for i in idx
        n1 = network.links[i, 1]
        n2 = network.links[i, 2]
        plot!(
            coord[[n1, n2], 1],
            coord[[n1, n2], 2],
            coord[[n1, n2], 3],
            args...;
            kw...,
        )
        # quiver needs to be implemented in Plots.jl but we can use python.
        """
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        """
    end
    return fig
end

function plotNodes!(fig, network::DislocationNetwork, args...; kw...)
    idx = idxLabel(network, -1; condition = !=)
    coord = network.coord
    for i in idx
        n1 = network.links[i, 1]
        n2 = network.links[i, 2]
        plot!(
            coord[[n1, n2], 1],
            coord[[n1, n2], 2],
            coord[[n1, n2], 3],
            args...;
            kw...,
        )
        # quiver needs to be implemented in Plots.jl but we can use python.
        """
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        """
    end
    return display()
end

function plotNodes(loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != -1, loop.label)
    coord = loop.coord
    fig = plot()
    for i in idx
        n1 = loop.links[i, 1]
        n2 = loop.links[i, 2]
        plot!(
            coord[[n1, n2], 1],
            coord[[n1, n2], 2],
            coord[[n1, n2], 3],
            args...;
            kw...,
        )
        # quiver needs to be implemented in Plots.jl but we can use python.
        """
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        """
    end
    return fig
end

function plotNodes!(fig, loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != -1, loop.label)
    coord = loop.coord
    for i in idx
        n1 = loop.links[i, 1]
        n2 = loop.links[i, 2]
        plot!(
            coord[[n1, n2], 1],
            coord[[n1, n2], 2],
            coord[[n1, n2], 3],
            args...;
            kw...,
        )
        # quiver needs to be implemented in Plots.jl but we can use python.
        """
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        """
    end
    return fig
end
