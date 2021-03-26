Base.:(==)(x::nodeTypeFE, y::Real) = ==(Int(x), y)
Base.:(==)(x::Real, y::nodeTypeFE) = ==(x, Int(y))
Base.convert(::Type{nodeTypeFE}, x::Real) = nodeTypeFE(Int(x))
Base.zero(::Type{nodeTypeFE}) = nodeTypeFE(0)
