"""
```
plotNodes(network::DislocationNetwork, args...; kw...)
```
Plots [`DislocationNetwork`](@ref).
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
        label[n1n2[1]] == extDln || label[n1n2[2]] == extDln && continue
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        # quiver needs to be implemented in Plots.jl but we can use python.
        #= 
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end
    plot!(fig; xlabel = "x", ylabel = "y", zlabel = "z")

    return fig
end
"""
```
plotNodes!(fig, network::DislocationNetwork, args...; kw...)
```
In-place plots [`DislocationNetwork`](@ref).
"""
function plotNodes!(fig, network::DislocationNetwork, args...; kw...)
    label = network.label
    coord = network.coord
    links = network.links
    numNode = network.numNode[1]
    elemT = eltype(links)
    for i in 1:numNode
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        label[n1n2[1]] == extDln || label[n1n2[2]] == extDln && continue
        plot!(coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        #= 
        # quiver needs to be implemented in Plots.jl but we can use python.
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end
    plot!(fig; xlabel = "x", ylabel = "y", zlabel = "z")
    return fig
end

"""
```
plotNodes(mesh::AbstractMesh, network::DislocationNetwork, args...; kw...)
```
Plots [`DislocationNetwork`](@ref) inside [`AbstractMesh`](@ref).
"""
function plotNodes(mesh::AbstractMesh, network::DislocationNetwork, args...; kw...)
    label = network.label
    coord = network.coord
    links = network.links
    numNode = network.numNode[1]
    elemT = eltype(links)
    fig = plot()
    for i in 1:numNode
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        label[n1n2[1]] == extDln || label[n1n2[2]] == extDln ? continue : nothing
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        # quiver needs to be implemented in Plots.jl but we can use python.
        #= 
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end

    # Plot FEM domain.
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, :)
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, :)
    face1 = [mesh.faces[:, 5]; mesh.faces[1, 5]]
    face2 = [mesh.faces[:, 6]; mesh.faces[1, 6]]
    surf1 = vertices[:, face1]
    surf2 = vertices[:, face2]
    plot!(fig, surf1[1, :], surf1[2, :], surf1[3, :], linecolor = :black, linewidth = 2)
    plot!(fig, surf2[1, :], surf2[2, :], surf2[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [1, 5]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [2, 6]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [3, 7]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [4, 8]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)

    xlims = (
        minimum(vertices[1, :]) - 0.2 * maximum(vertices[1, :]),
        maximum(vertices[1, :]) + 0.2 * maximum(vertices[1, :]),
    )
    ylims = (
        minimum(vertices[2, :]) - 0.2 * maximum(vertices[2, :]),
        maximum(vertices[2, :]) + 0.2 * maximum(vertices[2, :]),
    )
    zlims = (
        minimum(vertices[3, :]) - 0.2 * maximum(vertices[3, :]),
        maximum(vertices[3, :]) + 0.2 * maximum(vertices[1, :]),
    )
    plot!(fig, xlims = xlims, ylims = ylims, zlims = zlims)
    plot!(fig; xlabel = "x", ylabel = "y", zlabel = "z")
    return fig
end

"""
```
plotNodes!(fig, mesh::AbstractMesh, network::DislocationNetwork, args...; kw...)
```
In-place plots [`DislocationNetwork`](@ref) inside [`AbstractMesh`](@ref).
"""
function plotNodes!(
    fig,
    mesh::T1,
    network::DislocationNetwork,
    args...;
    kw...,
) where {T1 <: AbstractMesh}
    label = network.label
    coord = network.coord
    links = network.links
    numNode = network.numNode[1]
    elemT = eltype(links)
    for i in 1:numNode
        n1n2 = SVector{2,elemT}(links[1, i], links[2, i])
        label[n1n2[1]] == extDln || label[n1n2[2]] == extDln ? continue : nothing
        plot!(fig, coord[1, n1n2], coord[2, n1n2], coord[3, n1n2], args...; kw...)
        # quiver needs to be implemented in Plots.jl but we can use python.
        #= 
        lVec = coord[n2, :] - coord[n1, :]
        quiver!([coord[n1,1]], [coord[n1,2]], [coord[n1,3]], args...; quiver=([lVec[1]], [lVec[2]], [lVec[3]]), kw...) =#
    end

    # Plot FEM domain.
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, :)
    vertices = reshape(collect(Iterators.flatten(mesh.vertices.vertices)), 3, :)
    face1 = [mesh.faces[:, 5]; mesh.faces[1, 5]]
    face2 = [mesh.faces[:, 6]; mesh.faces[1, 6]]
    surf1 = vertices[:, face1]
    surf2 = vertices[:, face2]
    plot!(fig, surf1[1, :], surf1[2, :], surf1[3, :], linecolor = :black, linewidth = 2)
    plot!(fig, surf2[1, :], surf2[2, :], surf2[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [1, 5]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [2, 6]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [3, 7]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)
    side = vertices[:, [4, 8]]
    plot!(fig, side[1, :], side[2, :], side[3, :], linecolor = :black, linewidth = 2)

    xlims = (minimum(vertices[1, :]), maximum(vertices[1, :]))
    ylims = (minimum(vertices[2, :]), maximum(vertices[2, :]))
    zlims = (minimum(vertices[3, :]), maximum(vertices[3, :]))
    plot!(fig, xlims = xlims, ylims = ylims, zlims = zlims)
    plot!(fig; xlabel = "x", ylabel = "y", zlabel = "z")

    return fig
end

"""
```
plotNodes(loop::DislocationLoop, args...; kw...)
```
Plots [`DislocationLoop`](@ref).
"""
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
"""
```
plotNodes!(fig, loop::DislocationLoop, args...; kw...)
```
In-place plots [`DislocationLoop`](@ref).
"""
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
    plot!(fig; xlabel = "x", ylabel = "y", zlabel = "z")
    return fig
end

"""
```
plotFEDomain(mesh::AbstractMesh)
```
Plots corners, edges and surfaces of [`AbstractMesh`](@ref).
"""
function plotFEDomain(mesh::AbstractMesh, args...; kw...)
    surfNode = mesh.surfNode
    cornerNode = mesh.cornerNode
    edgeNode = mesh.edgeNode
    faceNode = mesh.faceNode
    coord = mesh.coord

    cornerIdx =
        collect(Iterators.flatten([surfNode[cornerNode[i]] for i in 1:length(cornerNode)]))
    fig = scatter(
        coord[1, cornerIdx],
        coord[2, cornerIdx],
        coord[3, cornerIdx],
        markershape = :diamond,
        markersize = 3,
        label = "Corners",
        args...;
    kw...,
    )

    @inbounds for i in 1:length(edgeNode)
        scatter!(
            fig,
            coord[1, surfNode[edgeNode[i]]],
            coord[2, surfNode[edgeNode[i]]],
            coord[3, surfNode[edgeNode[i]]],
            markershape = :circle,
            markersize = 3,
            label = "Edge $i",
            args...;
            kw...,
        )
    end
    @inbounds for i in 1:length(faceNode)
        scatter!(
            fig,
            coord[1, surfNode[faceNode[i]]],
            coord[2, surfNode[faceNode[i]]],
            coord[3, surfNode[faceNode[i]]],
            markershape = :square,
            markersize = 3,
            label = "Face $i",
            args...;
            kw...,
        )
    end
    plot!(fig; xlabel = "x", ylabel = "y", zlabel = "z")
    return fig
end

function plotBoundaries(boundaries::Boundaries, mesh::RegularCuboidMesh, args...; kw...)
    uGamma = boundaries.uGamma.node
    mGamma = boundaries.mGamma.node
    tGamma = boundaries.tGamma.node
    feCoord = mesh.coord

    fig = scatter(
        feCoord[1, uGamma],
        feCoord[2, uGamma],
        feCoord[3, uGamma];
        markershape = :square,
        markercolor = :red,
        markersize = 3,
        label = "Fixed",
    )
    scatter!(
        fig,
        feCoord[1, mGamma],
        feCoord[2, mGamma],
        feCoord[3, mGamma];
        markershape = :diamond,
        markercolor = :cyan,
        markersize = 3,
        label = "Load",
    )
    scatter!(
        fig,
        feCoord[1, tGamma],
        feCoord[2, tGamma],
        feCoord[3, tGamma];
        markershape = :circle,
        markercolor = :black,
        markersize = 3,
        label = "Traction",
    )
    plot!(fig, args...; xlabel = "x", ylabel = "y", zlabel = "z", kw...)
    return fig
end
