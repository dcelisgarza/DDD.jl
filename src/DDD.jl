module DDD

using LinearAlgebra, Plots, Statistics, InteractiveUtils, JSON

import Base: zero, isequal, isless, convert, ==, *, /, length, getindex
import Base: eachindex, push!, iterate

include("./Misc/Misc.jl")
# Miscelaneous.
export makeTypeDict
include("./Integration/CustomIntegration.jl")
include("./Material/MaterialBase.jl")
include("./Dislocation/DislocationBase.jl")
# Distributions.
export AbstractDistribution, Zeros, Rand
# Dislocation types.
export nodeType, SlipSystem, AbstractDlnStr, loopPrism, loopShear,
export DislocationLoop, DislocationNetwork
# Dislocation functions.
export makeNetwork, makeNetwork!
include("./FEM/FEMBase.jl")
include("./DislocationFEM/DislocationFEMBase.jl")
include("./IO/IOBase.jl")
# Imports.
export load, loadDislocationP, loadMaterialP, loadIntegrationP
export loadSlipSystem, loadDislocationLoop, loadParams
export loadDislocationLoop
# Export.
export save
include("./PostProcessing/Plotting.jl")
export plotNodes, plotNodes!
end # module
