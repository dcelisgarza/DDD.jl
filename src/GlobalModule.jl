"""
Author: Daniel Celis Garza
Date: 2020/01/23

Global module includes `CustomTypes.jl`, `DislocationParams.jl`, `MaterialParams.jl`
"""
module GlobalModule
include("CustomTypes.jl")
include("DislocationBase.jl")
include("MaterialBase.jl")
include("DdFemBase.jl")
include("CustomIntegration.jl")
include("Input.jl")
include("Output.jl")
end # module
