"""
```
plotNodes(obj::Union{DislocationLoop, DislocationNetwork}, args...; kw...)
```
Plots dislocation network as nodes connected by segments. Returns a new figure. See [`plotNodes!`](@ref) for mutating version.
"""
function plotNodes(network::DislocationNetwork, args...; kw...)
    idx = findall(x -> x != 0, network.links[1, :])
    coord = network.coord
    fig = plot()
    for i in idx
        n1 = network.links[1, i]
        n2 = network.links[2, i]
        plot!(coord[1, [n1, n2]], coord[2, [n1, n2]], coord[3, [n1, n2]], args...; kw...)
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
plotNodes!(fig, obj::Union{DislocationLoop, DislocationNetwork}, args...; kw...)
```
Updates figure to plot dislocation network as nodes connected by segments. See [`plotNodes`](@ref) for non-mutating version.
"""
function plotNodes!(fig, network::DislocationNetwork, args...; kw...)
    idx = findall(x -> x != 0, network.links[1, :])
    coord = network.coord
    for i in idx
        n1 = network.links[1, i]
        n2 = network.links[2, i]
        plot!(coord[1, [n1, n2]], coord[2, [n1, n2]], coord[3, [n1, n2]], args...; kw...)
        #=
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        =#
    end
    return fig
end
function plotNodes(loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != 0, loop.links[1, :])
    coord = loop.coord
    fig = plot()
    for i in idx
        n1 = loop.links[1, i]
        n2 = loop.links[2, i]
        plot!(coord[1, [n1, n2]], coord[2, [n1, n2]], coord[3, [n1, n2]], args...; kw...)
        #=
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        =#
    end
    return fig
end
function plotNodes!(fig, loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != 0, loop.links[1, :])
    coord = loop.coord
    for i in idx
        n1 = loop.links[1, i]
        n2 = loop.links[2, i]
        plot!(coord[1, [n1, n2]], coord[2, [n1, n2]], coord[3, [n1, n2]], args...; kw...)
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
        n1 = network.links[1, i]
        n2 = network.links[2, i]
        lines!(
            coord[1, [n1, n2]],
            coord[2, [n1, n2]],
            coord[3, [n1, n2]],
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
        n1 = network.links[1, i]
        n2 = network.links[2, i]
        lines!(
            coord[1, [n1, n2]],
            coord[2, [n1, n2]],
            coord[3, [n1, n2]],
            args...;
            kw...,
        )
    end
    return fig
end
=#
