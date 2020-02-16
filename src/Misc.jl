"""
```
inclusiveComparison(data, args...)
```
Compare data to a tuple, return `true` if it is equal to any arg, `false` if it is not equal to any.
"""
function inclusiveComparison(data, args...) :: Bool
    for i = 1:length(args)
        if data == args[i]
            return true
        end
    end
    return false
end
"""
```
compStruct(arg1, arg2)
```
Compares structures to see if they are equal.
"""
function compStruct(arg1, arg2)
    @assert typeof(arg1) == typeof(arg2)

    names = fieldnames(typeof(arg1))
    for i in names
        result = getproperty(arg1, i) == getproperty(arg2, i)
        if result == false
            return false
        end
    end
    return true
end
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
"""
```
intAngle(n::Int64)
```
Calculates the interior angle of a regular polygon with `n` sides.
"""
function intAngle(n::Int64)
    return (n - 2) * π / n
end
"""
```
extAngle(n::Int64)
```
Calculates the exterior angle of a regular polygon with `n` sides.
"""
function extAngle(n::Int64)
    return π - intAngle(n)
end
"""
```
rot3D(xyz::Vector{<:Float64}, uvw::Vector{Float64}, abc::Vector{<:Float64}, θ::Float64)
```
Rotate point `xyz` about the line with direction vector `uvw` that crosses the point `abc` by the angle `θ`.
"""
function rot3D(
    xyz::AbstractVector{<:Float64},
    uvw::AbstractVector{<:Float64},
    abc::AbstractVector{<:Float64},
    θ::Float64,
)
    #https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions
    isapprox(norm(uvw), 1.0) ? nothing : uvw ./= norm(uvw)

    local costheta = cos(θ)
    local onemcostheta = 1 - costheta
    local sintheta = sin(θ)
    local xyzDOTuvw = dot(xyz, uvw)

    return [
        (abc[1] * (uvw[2]^2 + uvw[3]^2) -
         uvw[1] * (abc[2] * uvw[2] + abc[3] * uvw[3] - xyzDOTuvw)) *
        onemcostheta +
        xyz[1] * costheta +
        sintheta * (-abc[3] * uvw[2] + abc[2] * uvw[3] - uvw[3] * xyz[2] +
         uvw[2] * xyz[3])

        (abc[2] * (uvw[1]^2 + uvw[3]^2) -
         uvw[2] * (abc[1] * uvw[1] + abc[3] * uvw[3] - xyzDOTuvw)) *
        onemcostheta +
        xyz[2] * costheta +
        sintheta *
        (abc[3] * uvw[1] - abc[1] * uvw[3] + uvw[3] * xyz[1] - uvw[1] * xyz[3])

        (abc[3] * (uvw[1]^2 + uvw[2]^2) -
         uvw[3] * (abc[1] * uvw[1] + abc[2] * uvw[2] - xyzDOTuvw)) *
        onemcostheta +
        xyz[3] * costheta +
        sintheta * (-abc[2] * uvw[1] + abc[1] * uvw[2] - uvw[2] * xyz[1] +
         uvw[1] * xyz[2])
    ]
end
