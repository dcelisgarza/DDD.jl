function compStruct(arg1, arg2)
    @assert typeof(arg1) == typeof(arg2)

    names = fieldnames(typeof(arg1))
    for i in names
        result = getproperty(arg1, i) == getproperty(arg1, i)
        if result == false
            return false
        end
    end
    return true
end

macro makeType(type::Any, supertype::Any)
    type = Symbol(type)
    supertype = Symbol(supertype)
    ex = quote
        struct $type <: $supertype end
    end
    esc(ex)
end

macro string_as_varname_macro(s::AbstractString, v::Any)
    s = Symbol(s)
    esc(:($s = $v))
end

function intAngle(numSides::Integer)
    return (numSides - 2) * π / numSides
end

function extAngle(numSides::Integer)
    return π - intAngle(numSides)
end

function arbRot3D(
    abc::Vector{<:Real},
    uvw::Vector{Real},
    vec::Vector{<:Real},
    θ::Real,
)
    #https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions
    local costheta = cos(θ)
    local onemcostheta = 1 - costheta
    local sintheta = sin(θ)
    local dotxyzuvw = dot(xyz, uvw)

    return [
        (abc[1] * (uvw[2]^2 + uvw[3]^2) -
         uvw[1] * (abc[2] * uvw[2] + abc[3] * uvw[3] - dotxyzuvw)) *
        onemcostheta +
        xyz[1] * costheta +
        sintheta * (-abc[3] * uvw[2] + abc[2] * uvw[3] - uvw[3] * xyz[2] +
         uvw[2] * xyz[3])

        (abc[2] * (uvw[1]^2 + uvw[3]^2) -
         uvw[2] * (abc[1] * uvw[1] + abc[3] * uvw[3] - dotxyzuvw)) *
        onemcostheta +
        xyz[2] * costheta +
        sintheta *
        (abc[3] * uvw[1] - abc[1] * uvw[3] + uvw[3] * xyz[1] - uvw[1] * xyz[3])

        (abc[3] * (uvw[1]^2 + uvw[2]^2) -
         uvw[3] * (abc[1] * uvw[1] + abc[2] * uvw[2] - dotxyzuvw)) *
        onemcostheta +
        xyz[3] * costheta +
        sintheta * (-abc[2] * uvw[1] + abc[1] * uvw[2] - uvw[2] * xyz[1] +
         uvw[1] * xyz[2])
    ]
end
