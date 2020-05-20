module DDD

using LinearAlgebra, Plots, Statistics, InteractiveUtils, JSON

import Base: zero, isequal, isless, convert, ==, *, /, length, getindex
import Base: eachindex, push!, iterate

include("./Misc/Misc.jl")
# Miscelaneous.
export makeTypeDict, compStruct, intAngle, extAngle, rot3D
include("./Integration/CustomIntegration.jl")
include("./Material/MaterialBase.jl")
include("./Dislocation/DislocationBase.jl")
# Distributions.
export AbstractDistribution, Zeros, Rand, Randn, Regular
# Dislocation types.
export nodeType, SlipSystem, AbstractDlnStr, loopPrism, loopShear
export DislocationLoop, DislocationNetwork, DislocationNetwork!
# Dislocation functions.
export checkNetwork, loopDistribution, calcSelfForce
export calcSegSegForce, calcParSegSegForce, calcSegForce
# Topology
export removeNode!, mergeNode!
include("./FEM/FEMBase.jl")
export shapeFunction, shapeFunctionDeriv, LinearQuadrangle3D
include("./DislocationFEM/DislocationFEMBase.jl")
include("./IO/IOBase.jl")
# Imports.
export load, loadDislocationP, loadMaterialP, loadIntegrationP
export loadSlipSystem, loadDislocationLoop, loadParams
export loadDislocationLoop, loadNetwork
# Export.
export save
include("./PostProcessing/Plotting.jl")
export plotNodes, plotNodes!
end # module
