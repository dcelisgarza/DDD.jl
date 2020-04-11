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
    esc(ex)
end
"""
```
string_as_varname_macro(s::AbstractString, v::Any)
```
Turn a string `s` into a variable whose name is the string `s` and value is `v`.
"""
macro string_as_varname_macro(s::AbstractString, v::Any)
    s = Symbol(s)
    esc(:($s = $v))
end
