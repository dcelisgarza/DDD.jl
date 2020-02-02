"""
Author: Daniel Celis Garza
Date: 2020/01/23

Module for custom types and function overloads.

# Defines
    Point3D, overloaded functions (a::Point3D, b::Number)
# Exports
    Point3D
# Imports
    Base: +, -, *, /, ^
"""
module CustomTypes
import Base: +, -, *, /, ^, zero, size

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

# macro string_as_varname_macro(s::AbstractString, v::Any)
#     s = Symbol(s)
#     esc(:($s = $v))
# end
end # module
