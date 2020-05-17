"""
```
plotNodes(network::DislocationNetwork, args...; kw...)
```
Plots dislocation network as nodes connected by segments. Returns a new figure. See [`plotNodes!`](@ref) for mutating version.
"""
function plotNodes(network::DislocationNetwork, args...; kw...)
    idx = findall(x -> x != 0, network.label)
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
        #=
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        =#
    end
    return fig
end

"""
```
plotNodes!(fig, network::DislocationNetwork, args...; kw...)
```
Updates figure to plot dislocation network as nodes connected by segments. See [`plotNodes`](@ref) for non-mutating version.
"""
function plotNodes!(fig, network::DislocationNetwork, args...; kw...)
    idx = findall(x -> x != 0, network.label)
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
        #=
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        =#
    end
    return fig
end

"""
```
plotNodes(loop::DislocationLoop, args...; kw...)
```
Plots dislocation network as nodes connected by segments. Returns a new figure. See [`plotNodes!`](@ref) for mutating version.
"""
function plotNodes(loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != 0, loop.label)
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
        #=
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        =#
    end
    return fig
end

"""
```
plotNodes!(fig, loop::DislocationLoop, args...; kw...)
```
Updates figure to plot dislocation network as nodes connected by segments. See [`plotNodes`](@ref) for non-mutating version.
"""
function plotNodes!(fig, loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != 0, loop.label)
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
        #=
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        =#
    end
    return fig
end

#=
function plotNodesMakie(network::DislocationNetwork, args...; kw...)
    idx = findall(x -> x != 0, loop.label)
    coord = network.coord
    fig = Scene()
    meshscatter!(coord, args...; kw...)
    for i in idx
        n1 = network.links[i, 1]
        n2 = network.links[i, 2]
        lines!(
            coord[[n1, n2], 1],
            coord[[n1, n2], 2],
            coord[[n1, n2], 3],
            args...;
            kw...,
        )
    end
    return fig
end

function plotNodesMakie!(fig, network::DislocationNetwork, args...; kw...)
    idx = findall(x -> x != 0, loop.label)
    coord = network.coord
    meshscatter!(coord, args...; kw...)
    for i in idx
        n1 = network.links[i, 1]
        n2 = network.links[i, 2]
        lines!(
            coord[[n1, n2], 1],
            coord[[n1, n2], 2],
            coord[[n1, n2], 3],
            args...;
            kw...,
        )
    end
    return fig
end
=#
