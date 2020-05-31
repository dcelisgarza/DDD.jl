"""
Related functions: [`coordLbl`](@ref)
```
idxLabel(network::DislocationNetwork, label::Int; condition::Function = ==)
```
Find indices for dislocations whose `label` meets `condition(x, label)`.
"""
@inline function idxLabel(network::DislocationNetwork, label::Int; condition::Function = ==)
    return findall(x -> condition(x, label), network.label)
end
"""
Related functions: [`idxLabel`](@ref)
```
coordLbl(network::DislocationNetwork, label::Int)
```
Get coordinates for the nodes with a given label (node type).
"""
@inline function coordLbl(network::DislocationNetwork, label::Int)
    return network.coord[:, idxLabel(network, label)]
end
"""
```
coordIdx(network::DislocationNetwork,
    index::Union{Int,AbstractArray{<:Int, N}) where {N}
```
Get coordinates for the node(s) that with the `index` or vector of indices
provided.
"""
@inline function coordIdx(
    network::DislocationNetwork,
    index::Union{Int, AbstractArray{<:Int, N}},
) where {N}
    return network.coord[:, index]
end
"""
Related functions: [`dataCond`](@ref)
```
idxCond(network::DislocationNetwork, fieldname::Symbol,
    args...; condition::Function)
```
Find index/indices whose `fieldname` meets the `condition(fieldname, args...)` where condition can be any function that uses `fieldname` and `args` to make a comparison.
```
idxCond(data::Union{
        AbstractArray{<:Real,N1},
        AbstractArray{<:nodeType, {N2}}
    },
    val::Real; condition::Function = ==) where {N1,N2}
```
Find index/indices of data that meet the condition(data, val).
```
idxCond(network::DislocationNetwork, fieldname::Symbol, val::Real;
    condition::Function = ==)
```
Find index/indices of node whose `fieldname`, `(:links, :bVec, :slipPlane, :coord, :label, :numNode, :numSeg)`, meets `condition(fieldname, val)`. If `fieldName` is multidimensional it will return `CartesianIndex`.
```
idxCond(network::DislocationNetwork, fieldname::Symbol, idxComp::Int,
    val::Real; condition::Function = ==)
```
Find index/indices of node whose `fieldname` meets `condition(fieldname[:, idxComp], val)`. It errors if the fieldname provided does not have a column `idxComp`.
"""
@inline function idxCond(
    network::DislocationNetwork,
    fieldname::Symbol,
    condition::Function,
    args...,
)
    return findall(x -> condition(x, args...), getproperty(network, fieldname))
end
@inline function idxCond(
    data::Union{AbstractArray{<:Real, N1}, AbstractArray{<:nodeType, N2}},
    val::Real;
    condition::Function = ==,
) where {N1, N2}
    return findall(x -> condition(x, val), data)
end
@inline function idxCond(
    network::DislocationNetwork,
    fieldname::Symbol,
    val::Real;
    condition::Function = ==,
)
    return findall(x -> condition(x, val), getproperty(network, fieldname))
end
@inline function idxCond(
    network::DislocationNetwork,
    fieldname::Symbol,
    idxComp::Int,
    val::Real;
    condition::Function = ==,
)
    return findall(x -> condition(x, val), getproperty(network, fieldname)[idxComp, :])
end
"""
Related functions: `idxCond`
```
dataCond(network::DislocationNetwork, dataField::Symbol, val::Real;
    condition::Function = ==)
```
Get the data whose `dataField`, `(:links, :bVec, :slipPlane, :coord, :label, :numNode, :numSeg)`, meets the `condition(dataField, val)`. If `dataField` is multidimensional it will use `CartesianIndex` to find it.
```
dataCond(network::DislocationNetwork, dataField::Symbol, idxComp::Int,
    val::Real; condition::Function = ==)
```
Get the data whose `dataField`, meets the `condition(dataField[:,idxComp], val)`. It errors if the fieldname provided does not have a column `idxComp`.
```
dataCond(network::DislocationNetwork, dataField::Symbol, condField::Symbol,
    val::Real; condition::Function = ==)
```
Get the data from `dataField` that corresponds to the `condField` (from the same pool as `dataField`) that meets the `condition(condField, val)`. If `condField` is multidimensional it will use `CartesianIndex` in its search. `dataField` and `condField` must have the same number of rows.
```
dataCond(network::DislocationNetwork, dataField::Symbol, condField::Symbol,
    idxComp::Int, val::Real; condition::Function = ==)
```
Get the data from `dataField` that corresponds to the `condField` that meets the `condition(condField[:, idxComp], val)`. `dataField` and `condField` must have the same number of rows.
"""
@inline function dataCond(
    network::DislocationNetwork,
    dataField::Symbol,
    val::Real;
    condition::Function = ==,
)
    data = getproperty(network, dataField)
    idx = idxCond(data, val; condition = condition)
    return data[idx]
end
@inline function dataCond(
    network::DislocationNetwork,
    dataField::Symbol,
    idxComp::Int,
    val::Real;
    condition::Function = ==,
)
    data = getproperty(network, dataField)
    idx = idxCond(data[idxComp, :], val; condition = condition)
    return data[:, idx]
end
@inline function dataCond(
    network::DislocationNetwork,
    dataField::Symbol,
    condField::Symbol,
    val::Real;
    condition::Function = ==,
)
    data = getproperty(network, dataField)
    cond = getproperty(network, condField)
    @assert size(data, 2) == size(cond, 1) ||
            size(data, 1) == size(cond, 2) ||
            size(data) == size(cond) "Number of rows of both fields must be equal."
    idx = idxCond(cond, val; condition = condition)
    return if ndims(cond) > 1
        return data[idx]
    else
        return data[:, idx]
    end
end
@inline function dataCond(
    network::DislocationNetwork,
    dataField::Symbol,
    condField::Symbol,
    idxComp::Int,
    val::Real;
    condition::Function = ==,
)
    data = getproperty(network, dataField)
    cond = getproperty(network, condField)
    @assert size(data, 1) == size(cond, 1) "Number of rows of both fields must be equal."
    idx = idxCond(cond[idxComp, :], val; condition = condition)
    return data[:, idx]
end
