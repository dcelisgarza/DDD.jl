"""
```
plotNodes(obj::Union{DislocationLoop, DislocationNetwork}, args...; kw...)
```
Plots dislocation network as nodes connected by segments. Returns a new figure. See [`plotNodes!`](@ref) for mutating version.
"""
function plotNodes(network::DislocationNetwork, args...; kw...)
    idx = findall(x -> x != 0, @view network.links[1, :])
    coord = network.coord
    links = network.links
    elemT = eltype(links)
    fig = plot()
    for i in idx
        n1n2 = SVector{2, elemT}(links[1, i], links[2, i])
        plot!(coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
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
    idx = findall(x -> x != 0, @view network.links[1, :])
    coord = network.coord
    links = network.links
    elemT = eltype(links)
    for i in idx
        n1n2 = SVector{2, elemT}(links[1, i], links[2, i])
        plot!(coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        #=
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        =#
    end
    return fig
end
function plotNodes(loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != 0, @view loop.links[1, :])
    coord = loop.coord
    links = loop.links
    elemT = eltype(links)
    fig = plot()
    for i in idx
        n1n2 = SVector{2, elemT}(links[1, i], links[2, i])
        plot!(coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        #=
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...)
        =#
    end
    return fig
end
function plotNodes!(fig, loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != 0, @view loop.links[1, :])
    coord = loop.coord
    links = loop.links
    elemT = eltype(links)
    for i in idx
        n1n2 = SVector{2, elemT}(links[1, i], links[2, i])
        plot!(coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
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
    idx = findall(x -> x != 0, @view loop.label)
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
    idx = findall(x -> x != 0, @view loop.label)
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
