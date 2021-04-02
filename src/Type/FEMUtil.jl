Base.:(==)(x::nodeTypeFE, y::Real) = ==(Int(x), y)
Base.:(==)(x::Real, y::nodeTypeFE) = ==(x, Int(y))
Base.convert(::Type{nodeTypeFE}, x::Real) = nodeTypeFE(Int(x))
Base.zero(::Type{nodeTypeFE}) = nodeTypeFE(0)
function findSurfaceElemCoord(mesh::AbstractMesh, idx = nothing)
    isnothing(idx) ? idx = (1:size(mesh.surfElemNode, 1)) : nothing
    return reshape(mesh.coord[:, vec(mesh.surfElemNode[idx,:])], :, 4)
end
