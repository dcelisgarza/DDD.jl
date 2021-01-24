"""
```
plotNodes(obj::Union{DislocationLoop, DislocationNetwork}, args...; kw...)
```
Plots dislocation network as nodes connected by segments. Returns a new figure. See [`plotNodes!`](@ref) for mutating version.
"""
function plotNodes(network::DislocationNetwork, args...; kw...)
    label = network.label
    coord = network.coord
    links = network.links
    numNode = network.numNode[1]
    elemT = eltype(links)
    fig = plot()
    for i in 1:numNode
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        label[n1n2[1]] == ext || label[n1n2[2]] == ext ? continue : nothing
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        # quiver needs to be implemented in Plots.jl but we can use python.
        #= 
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
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
    label = network.label
    coord = network.coord
    links = network.links
    numNode = network.numNode[1]
    elemT = eltype(links)
    for i in 1:numNode
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        label[n1n2[1]] == ext || label[n1n2[2]] == ext ? continue : nothing
        plot!(coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        #= 
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end
    return fig
end

function plotNodes(mesh::T1, network::DislocationNetwork, args...; kw...) where {T1 <: AbstractMesh}
    label = network.label
    coord = network.coord
    links = network.links
    numNode = network.numNode[1]
    elemT = eltype(links)
    fig = plot()
    for i in 1:numNode
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        label[n1n2[1]] == ext || label[n1n2[2]] == ext ? continue : nothing
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        # quiver needs to be implemented in Plots.jl but we can use python.
        #= 
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end

    # Plot FEM domain.
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, :)
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, :)
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
    
    xlims = (minimum(vertices[1, :]) - 0.2 * maximum(vertices[1, :]), maximum(vertices[1, :]) + 0.2 * maximum(vertices[1, :]))
    ylims = (minimum(vertices[2, :]) - 0.2 * maximum(vertices[2, :]), maximum(vertices[2, :]) + 0.2 * maximum(vertices[2, :]))
    zlims = (minimum(vertices[3, :]) - 0.2 * maximum(vertices[3, :]), maximum(vertices[3, :]) + 0.2 * maximum(vertices[1, :]))
    plot!(fig, xlims = xlims, ylims = ylims, zlims = zlims)
    return fig
end

function plotNodes!(fig, mesh::T1, network::DislocationNetwork, args...; kw...) where {T1 <: AbstractMesh}
    label = network.label
    coord = network.coord
    links = network.links
    numNode = network.numNode[1]
    elemT = eltype(links)
    for i in 1:numNode
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        label[n1n2[1]] == ext || label[n1n2[2]] == ext ? continue : nothing
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        # quiver needs to be implemented in Plots.jl but we can use python.
        #= 
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end

    # Plot FEM domain.
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, :)
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, :)
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
    
    xlims = (minimum(vertices[1, :]), maximum(vertices[1, :]))
    ylims = (minimum(vertices[2, :]), maximum(vertices[2, :]))
    zlims = (minimum(vertices[3, :]), maximum(vertices[3, :]))
    plot!(fig, xlims = xlims, ylims = ylims, zlims = zlims)

    return fig
end

function plotNodes(loop::DislocationLoop, args...; kw...)
    idx = findall(x -> x != 0, loop.label)
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
    idx = findall(x -> x != 0, loop.label)
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

function plotFEDomain(mesh::T) where {T <: AbstractMesh}
    cornerNode = mesh.cornerNode
    edgeNode = mesh.edgeNode
    faceNode = mesh.faceNode
    coord = mesh.coord

    fig = scatter(coord[1, cornerNode], coord[2, cornerNode], coord[3, cornerNode], markershape = :diamond, markersize = 3, label = "Corners")
    @inbounds for i in 1:length(edgeNode)
        scatter!(fig, coord[1, edgeNode[i]], coord[2, edgeNode[i]], coord[3, edgeNode[i]], markershape = :circle, markersize = 3, label = "Edge $i")
    end
    @inbounds for i in 1:length(faceNode)
        scatter!(fig, coord[1, faceNode[i]], coord[2, faceNode[i]], coord[3, faceNode[i]], markershape = :square, markersize = 3, label = "Face $i")
    end
    return fig
end