"""
```
JSON.lower(
    t::T,
) where {
    T <: Union{
        AbstractCrystalStruct,
        AbstractMobility,
        AbstractIntegrator,
        AbstractDlnSeg,
        AbstractDlnStr,
        AbstractDistribution,
    },
}

JSON.lower(t::nodeType)
```
Extensions to `JSON.lower` for custom types. Allows these variables to be serialised properly.
"""
function JSON.lower(
    t::T,
) where {
    T <: Union{
        AbstractCrystalStruct,
        AbstractMobility,
        AbstractIntegrator,
        AbstractDlnSeg,
        AbstractDlnStr,
        AbstractDistribution,
    },
}
    return string(t)
end
JSON.lower(t::nodeType) = Int(t)

"""
```
save(filename::AbstractString, args...; mode::AbstractString = "w")
```
Wrapper for `JSON.print` to a file, `args` are the variables or structures you want to save.
"""
function save(filename::AbstractString, args...; mode::AbstractString = "w")
    return open(filename, mode) do io
        return JSON.print(io, args)
    end
end
