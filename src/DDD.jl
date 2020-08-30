module DDD

using LinearAlgebra, Plots, Statistics, InteractiveUtils, JSON, StaticArrays

# Miscelaneous.
include("./Misc/Misc.jl")
export makeTypeDict, compStruct, intAngle, externalAngle, rot3D, makeInstanceDict
export inclusiveComparison, ⊗

include("./Type/TypeBase.jl")
export BCC, FCC, HCP
export MaterialParameters
export nodeType, SlipSystem, DislocationParameters, DislocationLoop, DislocationNetwork, DislocationNetwork!, checkNetwork
export IntegrationParameters, IntegrationTime
export Rand, Randn, Zeros, Regular, loopDistribution
# include("./Material/MaterialBase.jl")
# export MaterialParameters, AbstractCrystalStruct, BCC, FCC, HCP
# include("./Dislocation/DislocationBase.jl")
# export DislocationParameters
# # Distributions.
# export AbstractDistribution, Zeros, Rand, Randn, Regular, mobBCC
# # Dislocation types.
# export nodeType, SlipSystem, AbstractDlnStr, loopPrism, loopShear
# export DislocationLoop, DislocationNetwork, DislocationNetwork!
# export getSegmentIdx!, getSegmentIdx, makeConnect!, makeConnect
# # Dislocation functions.
# export checkNetwork, loopDistribution, calcSelfForce, calcSelfForce!
# export calcSegSegForce, calcSegSegForce!, calcParSegSegForce, calcSegForce, calcSegForce!
# export dlnMobility, dlnMobility!
# # Topology
# export removeNode!, mergeNode, splitNode, coarsenNetwork, refineNetwork
# include("./FEM/FEMBase.jl")
# export shapeFunction, shapeFunctionDeriv, LinearQuadrangle3D
# include("./DislocationFEM/DislocationFEMBase.jl")
# include("./Integration/CustomIntegration.jl")
# export IntegrationParameters, IntegrationTime, CustomTrapezoid, integrate!
# include("./PostProcessing/Plotting.jl")
# export plotNodes, plotNodes!
include("./IO/IOBase.jl")
# Imports.
export load, loadDislocationParameters, loadMaterialParameters, loadIntegrationParameters
export loadSlipSystem, loadDislocationLoop, loadParams
export loadDislocationLoop, loadNetwork, loadIntegrationTime
# Export.
export save
end # module
