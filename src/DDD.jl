module DDD

using LinearAlgebra, Plots, Statistics, InteractiveUtils, JSON, StaticArrays

import Base: zero, isequal, isless, convert, ==, *, /, length, getindex
import Base: eachindex, push!, iterate, size, show

include("./Misc/Misc.jl")
# Miscelaneous.
export makeTypeDict, compStruct, intAngle, extAngle, rot3D, makeInstanceDict
export translateEnum, inclusiveComparison, dimDot, dimNorm
include("./Integration/CustomIntegration.jl")
export IntegrationP, CustomTrapezoid
include("./Material/MaterialBase.jl")
export MaterialP, AbstractCrystalStruct, BCC, FCC, HCP
include("./Dislocation/DislocationBase.jl")
export DislocationP
# Distributions.
export AbstractDistribution, Zeros, Rand, Randn, Regular, mobBCC
# Dislocation types.
export nodeType, SlipSystem, AbstractDlnStr, loopPrism, loopShear
export DislocationLoop, DislocationNetwork, DislocationNetwork!
export getSegmentIdx!, makeConnect!
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
