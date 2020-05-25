#=
    This is a place for code that was written and is now unused and probably won't be reused. Placed here just in case.
=#

"""
```
makeType(type::Any, supertype::Any)
```
Macro to make a `type` that is a subtype of `supertype`.
"""
macro makeType(type::Any, supertype::Any)
    type = Symbol(type)
    supertype = Symbol(supertype)
    ex = quote
        struct $type <: $supertype end
    end
    return esc(ex)
end
"""
```
string_as_varname_macro(s::AbstractString, v::Any)
```
Turn a string `s` into a variable whose name is the string `s` and value is `v`.
"""
macro string_as_varname_macro(s::AbstractString, v::Any)
    s = Symbol(s)
    return esc(:($s = $v))
end

"""
```
translateEnum(
    valType::DataType,
    dict::Dict{T1, T2},
    key::T3,
) where {T1, T2, T3}
```
Translates the string name of enumerated types to the actual Julia type.
"""
@inline function translateEnum(
    valType::DataType,
    dict::Dict{T1,T2},
    key::T3,
) where {T1,T2,T3}
    instanceDict = makeInstanceDict(valType)
    @inbounds @simd for i in eachindex(dict[key])
        dict[key][i] = instanceDict[dict[key][i]]
    end
    return dict
end

"""
```
translateEnum(
    valType::DataType,
    dict::Dict{T1, T2},
    key::T3,
) where {T1, T2, T3}
```
Translates the string name of enumerated types to the actual Julia type.
"""
@inline function translateEnum(
    valType::DataType,
    dict::Dict{T1,T2},
    key::T3,
) where {T1,T2,T3}
    instanceDict = makeInstanceDict(valType)
    @inbounds @simd for i in eachindex(dict[key])
        dict[key][i] = instanceDict[dict[key][i]]
    end
    return dict
end
