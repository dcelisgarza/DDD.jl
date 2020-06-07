module DDD

using LinearAlgebra, Plots, Statistics, InteractiveUtils, JSON, StaticArrays

include("./Misc/Misc.jl")
# Miscelaneous.
export makeTypeDict, compStruct, intAngle, extAngle, rot3D, makeInstanceDict
export translateEnum, inclusiveComparison, dimDot, dimNorm, ⊗
include("./Material/MaterialBase.jl")
export MaterialP, AbstractCrystalStruct, BCC, FCC, HCP
include("./Dislocation/DislocationBase.jl")
export DislocationP
# Distributions.
export AbstractDistribution, Zeros, Rand, Randn, Regular, mobBCC
# Dislocation types.
export nodeType, SlipSystem, AbstractDlnStr, loopPrism, loopShear
export DislocationLoop, DislocationNetwork, DislocationNetwork!
export getSegmentIdx!, getSegmentIdx, makeConnect!, makeConnect
# Dislocation functions.
export checkNetwork, loopDistribution, calcSelfForce, calcSelfForce!
export calcSegSegForce, calcSegSegForce!, calcParSegSegForce, calcSegForce, calcSegForce!
export dlnMobility, dlnMobility!
# Topology
export removeNode!, mergeNode!, splitNode!, coarsenNetwork!, refineNetwork!
include("./FEM/FEMBase.jl")
export shapeFunction, shapeFunctionDeriv, LinearQuadrangle3D
include("./DislocationFEM/DislocationFEMBase.jl")
include("./Integration/CustomIntegration.jl")
export IntegrationP, IntegrationVar, CustomTrapezoid, integrate!
include("./PostProcessing/Plotting.jl")
export plotNodes, plotNodes!
include("./IO/IOBase.jl")
# Imports.
export load, loadDislocationP, loadMaterialP, loadIntegrationP
export loadSlipSystem, loadDislocationLoop, loadParams
export loadDislocationLoop, loadNetwork, loadIntegrationVar
# Export.
export save
end # module
