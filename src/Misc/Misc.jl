#=
    Miscelaneous code that carries out general tasks.
=#

"""
```
function makeTypeDict(valType::DataType)
```
Inputs contain strings that correspond to DDD data types. This function atuomatically creates a dictionary for all concrete subtypes of a given `valType`.

# Examples

```@jldoctest
julia> abstract type MyAbstractType end
julia> struct MyStruct1 <: MyAbstractType end
julia> struct MyStruct2 <: MyAbstractType end
julia> makeTypeDict(MyAbstractType)
Dict{String,Any} with 4 entries:
  "DDD.MyStruct1()" => MyStruct1()
  "DDD.MyStruct2()" => MyStruct2()
  "MyStruct1()"     => MyStruct1()
  "MyStruct2()"     => MyStruct2()
```
"""
function makeTypeDict(valType::DataType)
    subTypes = subtypes(valType)
    dict = Dict{String, Any}()
    for subType in subTypes
        strSubType = string(subType) * "()"
        push!(dict, strSubType => eval(subType()))
        if strSubType[1:4] == "DDD."
            push!(dict, strSubType[5:end] => eval(subType()))
        else
            push!(dict, "DDD." * strSubType => eval(subType()))
        end
    end
    return dict
end

"""
```
inclusiveComparison(data, args...)
```
Compare data to a tuple, return `true` if it is equal to any arg, `false` if it is not equal to any.

# Examples

```@jldoctest
julia> inclusiveComparison("f", 1,4,5,"f")
true
julia> inclusiveComparison(23.246, 1.5, 4, 5, "f")
false
```
"""
function inclusiveComparison(data, args...)::Bool
    for i in eachindex(args)
        if data == args[i]
            return true
        end
    end
    return false
end

"""
```
compStruct(arg1, arg2; verbose::Bool=false)
```
Function that compares values of the fields of two variables `arg1` and `arg2` with the same structure. If `verbose = true`, it will print which fields are different from each other.

# Examples

```@jldoctest
julia> struct MyStruct1; x; end
julia> test1 = MyStruct1(1)
MyStruct1(1)
julia> test2 = MyStruct1(5)
MyStruct1(5)
julia> compStruct(test1, test2; verbose = true)
Structures differ in field: x.
false
julia> compStruct(1, 1; verbose = true)
true
julia> compStruct(1, [1]; verbose = true)
false
```
"""
function compStruct(arg1, arg2; verbose::Bool = false)
    if typeof(arg1) != typeof(arg2)
        @warn "compStruct: Variables have different types:
            typeof(arg1) = $(typeof(arg1))
            typeof(arg2) = $(typeof(arg2))"
        return false
    end
    names = fieldnames(typeof(arg1))
    flag::Bool = true
    for i in names
        result = getproperty(arg1, i) == getproperty(arg2, i)
        if result == false
            flag = false
            verbose ? println("Structures differ in field: $(i).") : nothing
        end
    end
    return flag
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
rot3D(xyz::Vector{<:Real}, uvw::Vector{<:Real}, abc::Vector{<:Real},
    θ::Real)
```
Rotate point `xyz` about the line with direction vector `uvw` that crosses the point `abc` by the angle `θ`. Further details found [here](https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions).

# Examples

```@jldoctest
julia> rot3D([1;1;1],[1;0;0],[0;0;0],π/2)
3-element Array{Float64,1}:
  1.0
 -0.9999999999999999
  1.0
julia> rot3D([1;1;1],[1;0;0],[0;0;0],-π/2)
3-element Array{Float64,1}:
1.0
1.0
-0.9999999999999999
julia> rot3D([1;1;1],[1;0;0],[0;0;0],π)
3-element Array{Float64,1}:
  1.0
 -1.0000000000000002
 -0.9999999999999999
```
"""
function rot3D(
    xyz::AbstractVector{<:Real},
    uvw::AbstractVector{<:Real},
    abc::AbstractVector{<:Real},
    θ::Real,
)

    isapprox(norm(uvw), 1.0) ? nothing : uvw ./= norm(uvw)

    cosθ = cos(θ)
    onemcosθ = 1 - cosθ
    sintheta = sin(θ)
    xyzDOTuvw = dot(xyz, uvw)

    return [
        (
            abc[1] * (uvw[2] * uvw[2] + uvw[3] * uvw[3]) -
            uvw[1] * (abc[2] * uvw[2] + abc[3] * uvw[3] - xyzDOTuvw)
        ) * onemcosθ +
        xyz[1] * cosθ +
        sintheta * (
            -abc[3] * uvw[2] + abc[2] * uvw[3] - uvw[3] * xyz[2] +
            uvw[2] * xyz[3]
        )
        (
            abc[2] * (uvw[1] * uvw[1] + uvw[3] * uvw[3]) -
            uvw[2] * (abc[1] * uvw[1] + abc[3] * uvw[3] - xyzDOTuvw)
        ) * onemcosθ +
        xyz[2] * cosθ +
        sintheta *
        (abc[3] * uvw[1] - abc[1] * uvw[3] + uvw[3] * xyz[1] - uvw[1] * xyz[3])
        (
            abc[3] * (uvw[1] * uvw[1] + uvw[2] * uvw[2]) -
            uvw[3] * (abc[1] * uvw[1] + abc[2] * uvw[2] - xyzDOTuvw)
        ) * onemcosθ +
        xyz[3] * cosθ +
        sintheta * (
            -abc[2] * uvw[1] + abc[1] * uvw[2] - uvw[2] * xyz[1] +
            uvw[1] * xyz[2]
        )
    ]
end