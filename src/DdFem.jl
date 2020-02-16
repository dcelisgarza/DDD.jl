"""
Cuboid mesh structure. Incomplete.
"""
mutable struct CuboidMesh{
    T1<:AbstractArray{<:Float64,N} where {N},
    T2<:AbstractArray{<:Float64,N} where {N},
}
    elem::T1
    vertices::T2
    coord::T2
    label::T1

    function CuboidMesh(
        elem,
        vertices = nothing,
        coord = nothing,
        label = nothing,
    )
        if vertices == nothing
            vertices = zeros(typeof(elem[1]), 8, 3)
        else
            @assert length(vertices) == 3
        end
        vertices[2:2:end, 1] .= elem[1]
        vertices[3:4, 2] .= elem[2]
        vertices[7:end, 2] .= elem[2]
        vertices[5:end, 3] .= elem[3]

        if !(coord == nothing && label == nothing)
            @assert length(elem) == size(coord, 2)
            @assert length(label) == size(coord, 1)
            new{typeof(elem),typeof(vertices)}(elem, vertices, coord, label)
        else
            new{typeof(elem),typeof(vertices)}(elem, vertices)
        end
    end #constructor
end #CuboidMesh

# function zero(::Type{CuboidMesh})
#     return CuboidMesh(zeros(3))
# end
