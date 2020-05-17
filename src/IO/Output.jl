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
} = string(t)
JSON.lower(t::nodeType) = Int(t)

"""
```
save(filename::AbstractString, args...; mode::AbstractString = "w")
```
Wrapper for `JSON.print` to a file, `args` are the variables or structures you want to save.
"""
function save(filename::AbstractString, args...; mode::AbstractString = "w")
    open(filename, mode) do io
        JSON.print(io, args)
    end
end
