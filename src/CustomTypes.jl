"""
Overloaded functions for strings.
"""
function size(::String)
    return ()
end

"""
Overloaded functions for symbols.
"""
function size(::Symbol)
    return ()
end

macro makeType(type::Any, supertype::Any)
    type = Symbol(type)
    supertype = Symbol(supertype)
    ex = quote
        struct $type <: $supertype end
    end
    esc(ex)
end

function compStruct(arg1, arg2)
    @assert typeof(arg1) == typeof(arg2)

    names = fieldnames(typeof(arg1))
    for i in names
        result = getproperty(arg1, i) == getproperty(arg1, i)
        if result == false return false end
    end
    return true
end

macro string_as_varname_macro(s::AbstractString, v::Any)
    s = Symbol(s)
    esc(:($s = $v))
end
