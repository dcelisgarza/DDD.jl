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
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        # quiver needs to be implemented in Plots.jl but we can use python.
        #= 
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end

    return fig
end

function plotNodes(mesh::T1, network::DislocationNetwork, args...; kw...) where {T1 <: AbstractMesh}
    idx = findall(x -> x != 0, @view network.links[1, :])
    coord = network.coord
    links = network.links
    elemT = eltype(links)
    fig = plot()
    for i in idx
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        # quiver needs to be implemented in Plots.jl but we can use python.
        #= 
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end

    # Plot FEM domain.
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, 8)
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, 8)
    face1 = SVector{5,Int}(1, 2, 4, 3, 1)
    face2 = SVector{5,Int}(5, 6, 8, 7, 5)
    surf1 = vertices[:, face1]
    surf2 = vertices[:, face2]
    plot!(fig, surf1[1, :], surf1[2, :], surf1[3, :], linecolor = :black, linewidth = 2)
    plot!(fig, surf2[1, :], surf2[2, :], surf2[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [1, 5]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [2, 6]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [3, 7]];
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [4, 8]];
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    
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
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        plot!(coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        #= 
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
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
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        #= 
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end
    return fig
end
function plotNodes!(fig, loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != 0, @view loop.links[1, :])
    coord = loop.coord
    links = loop.links
    elemT = eltype(links)
    for i in idx
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        #= 
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end
    return fig
end
function plotNodes!(fig, mesh::T1, network::DislocationNetwork, args...; kw...) where {T1 <: AbstractMesh}
    idx = findall(x -> x != 0, @view network.links[1, :])
    coord = network.coord
    links = network.links
    elemT = eltype(links)
    for i in idx
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        # quiver needs to be implemented in Plots.jl but we can use python.
        #= 
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end

    # Plot FEM domain.
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, 8)
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, 8)
    face1 = SVector{5,Int}(1, 2, 4, 3, 1)
    face2 = SVector{5,Int}(5, 6, 8, 7, 5)
    surf1 = vertices[:, face1]
    surf2 = vertices[:, face2]
    plot!(fig, surf1[1, :], surf1[2, :], surf1[3, :], linecolor = :black, linewidth = 2)
    plot!(fig, surf2[1, :], surf2[2, :], surf2[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [1, 5]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [2, 6]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [3, 7]];
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [4, 8]];
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    
    return fig
end