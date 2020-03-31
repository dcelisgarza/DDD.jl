struct RegularCuboidMesh{
    T1 <: AbstractArray{<:Int64, N} where {N},
    T2 <: AbstractArray{<:Float64, N} where {N},
}
    numElem::T1
    sizeElem::T2
    sizeMesh::T2
    stiffTensor::T2
    label::T1
    coord::T2
end
