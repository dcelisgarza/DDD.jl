dictkeys(d::Dict) = (collect(Symbol.(keys(d)))...,)
dictvalues(d::Dict) = (collect(values(d))...,)
namedtuple(d::Dict{T1,T2}) where {T1,T2} = NamedTuple{dictkeys(d)}(dictvalues(d))

"""
```
makeInstanceDict(valType::DataType)
```
Make a dictionary of enumerated variable instances. Helps in translating JSON files.
"""
function makeInstanceDict(valType::DataType)
    insts = instances(valType)
    dict = Dict{String,valType}()
    for inst in insts
        push!(dict, string(inst) => inst)
    end
    return dict
end

"""
```
subTypeTree(t; dict = Dict(), level = 1, cutoff = 0)
```
Create subtype dictionary. Adapted from https://github.com/JuliaLang/julia/issues/24741
"""
function subTypeTree(t; dict = Dict(), level = 1, cutoff = 0)
    level > cutoff ? push!(dict, t => supertype(t)) : nothing
    for s in subtypes(t)
        subTypeTree(s; level = level + 1, dict = dict, cutoff = cutoff)
    end
    return dict
end

"""
```
makeTypeDict(valType::DataType)
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
function makeTypeDict(valType::DataType; cutoff = 1)
    primitive = subTypeTree(valType; cutoff = cutoff)
    dict = Dict{String,Any}()
    for (key, val) in primitive
        strSubType = string(key) * "()"
        try
            push!(dict, strSubType => key())
            if strSubType[1:4] == "DDD."
                push!(dict, strSubType[5:end] => key())
            else
                push!(dict, "DDD." * strSubType => key())
            end
        catch
        end
    end
    return dict
end

"""
```
inclusiveComparison(data, args...)::Bool
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
compStruct(arg1, arg2; verbose::Bool = false)
```
Compares values of the fields of two variables `arg1` and `arg2` with the same structure. If `verbose = true`, it will print which fields are different from each other.

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
        !verbose ? nothing :
        @warn "compStruct: Variables have different types:\n\ttypeof(arg1) = $(typeof(arg1))\n\ttypeof(arg2) = $(typeof(arg2))"
        return false
    end
    names = fieldnames(typeof(arg1))
    flag::Bool = true
    for name in names
        result = getproperty(arg1, name) == getproperty(arg2, name)
        if result == false
            flag = false
            !verbose ? nothing : @warn "Structures differ in field: $(name)."
        end
    end
    return flag
end

"""
```
internalAngle(n::Int)
```
Compute the interior angle of a regular polygon with `n` sides.
"""
internalAngle(n::Int) = (n - 2) * π / n

"""
```
externalAngle(n::Int)
```
Compute the exterior angle of a regular polygon with `n` sides.
"""
externalAngle(n::Int) = π - internalAngle(n)

"""
```
rot3D(xyz, uvw, abc, θ)
```
Rotate point `xyz` about the vector `uvw` that crosses point `abc` by the angle `θ`. Further details found [here](https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions).

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
function rot3D(xyz, uvw, abc, θ)
    uvw = normalize(uvw)

    cosθ = cos(θ)
    onemcosθ = 1 - cosθ
    sinθ = sin(θ)
    xyzDOTuvw = dot(xyz, uvw)

    return SVector{3,typeof(xyzDOTuvw)}(
        (
            abc[1] * (uvw[2] * uvw[2] + uvw[3] * uvw[3]) -
            uvw[1] * (abc[2] * uvw[2] + abc[3] * uvw[3] - xyzDOTuvw)
        ) * onemcosθ +
        xyz[1] * cosθ +
        sinθ * (-abc[3] * uvw[2] + abc[2] * uvw[3] - uvw[3] * xyz[2] + uvw[2] * xyz[3]),
        (
            abc[2] * (uvw[1] * uvw[1] + uvw[3] * uvw[3]) -
            uvw[2] * (abc[1] * uvw[1] + abc[3] * uvw[3] - xyzDOTuvw)
        ) * onemcosθ +
        xyz[2] * cosθ +
        sinθ * (abc[3] * uvw[1] - abc[1] * uvw[3] + uvw[3] * xyz[1] - uvw[1] * xyz[3]),
        (
            abc[3] * (uvw[1] * uvw[1] + uvw[2] * uvw[2]) -
            uvw[3] * (abc[1] * uvw[1] + abc[2] * uvw[2] - xyzDOTuvw)
        ) * onemcosθ +
        xyz[3] * cosθ +
        sinθ * (-abc[2] * uvw[1] + abc[1] * uvw[2] - uvw[2] * xyz[1] + uvw[1] * xyz[2]),
    )
end

"""
```
⊗(x::AbstractVector, y::AbstractVector)
```
Tensor product.
"""
⊗(x::AbstractVector{T1}, y::AbstractVector{T2}) where {T1,T2} = x * y'

"""
```
linePlaneIntersect(n::T, p0::T, l::T, l0::T) where {T <: AbstractVector}
```
Finds the intersect between a line and a plane. `n` is the plane normal, `p0` is a point on the plane, `l` is the line vector, `l0` is a point on a line.
"""
function linePlaneIntersect(n, p0, l, l0)
    den = l ⋅ n
    num = (p0 - l0) ⋅ n

    if den ≈ 0
        num ≈ 0 && return Inf
        return nothing
    end

    d = num / den
    intersect = l0 + d * l
    return intersect
end

"""
```
gausslegendre(n::Integer, a, b)
```
Compute Gauss-Legendre quadrature points and weights for the interval `[a, b]`.
"""
function FastGaussQuadrature.gausslegendre(n::Integer, a, b)
    bma = (b - a) * 0.5
    bpa = (b + a) * 0.5
    x, w = gausslegendre(n)
    @inbounds @simd for i in 1:n
        x[i] = bma * x[i] + bpa
        w[i] = bma * w[i]
    end
    return x, w
end