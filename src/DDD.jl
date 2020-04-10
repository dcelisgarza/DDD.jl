module DDD

using LinearAlgebra,
    Plots,
    Statistics,
    InteractiveUtils,
    JSON

import Base:
    zero, isequal, isless, convert, ==, *, /, length, getindex, eachindex, push!, iterate

include("./Misc/Misc.jl")
export makeTypeDict
include("./Integration/CustomIntegration.jl")
include("./Material/MaterialBase.jl")
include("./Dislocation/DislocationBase.jl")
export nodeType
export SlipSystem, AbstractDlnStr, loopPrism, loopShear, DislocationLoop
export AbstractDistribution, Zeros, Rand
include("./FEM/FEMBase.jl")
include("./DislocationFEM/DislocationFEMBase.jl")
include("./IO/IOBase.jl")
export load, loadSlipSystem, loadDislocationP, loadMaterialP, loadIntegrationP
export save
include("./PostProcessing/Plotting.jl")
end # module
