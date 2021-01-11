module DDD

using LinearAlgebra,
    SparseArrays, Plots, Statistics, InteractiveUtils, JSON, StaticArrays, FileIO

# Miscelaneous.
include("./Misc/Misc.jl")
export makeTypeDict, compStruct, intAngle, externalAngle, rot3D, makeInstanceDict
export inclusiveComparison, ⊗

include("./Type/TypeBase.jl")
export AbstractCrystalStruct,
    BCC,
    FCC,
    HCP,
    loopDln,
    loopPrism,
    loopShear,
    AbstractShapeFunction,
    AbstractShapeFunction2D,
    AbstractShapeFunction3D,
    LinearQuadrangle2D,
    LinearQuadrangle3D,
    LinearElement,
    AbstractMesh,
    AbstractRegularCuboidMesh,
    DispatchRegularCuboidMesh,
    RegularCuboidMesh,
    buildMesh,
    FEMParameters,
    ForceDisplacement
export MaterialParameters
export nodeType,
    SlipSystem,
    DislocationParameters,
    DislocationLoop,
    DislocationNetwork,
    DislocationNetwork!,
    checkNetwork,
    getSegmentIdx,
    getSegmentIdx!,
    makeConnect,
    makeConnect!
export IntegrationParameters, IntegrationTime
export Rand, Randn, Zeros, Regular, loopDistribution

include("./Processing/ProcessingBase.jl")
export calcSegForce,
    calcSegForce!,
    calcSelfForce,
    calcSelfForce!,
    calcSegSegForce,
    calcSegSegForce!,
    calc_σHat
export dlnMobility, dlnMobility!
export mergeNode!, splitNode!, coarsenNetwork!, refineNetwork!
export shapeFunction, shapeFunctionDeriv

include("./PostProcessing/Plotting.jl")
export plotNodes, plotNodes!

include("./IO/IOBase.jl")
export loadJSON,
    loadDislocationParametersJSON, loadMaterialParametersJSON, loadIntegrationParametersJSON
export loadSlipSystemJSON, loadDislocationLoopJSON, loadParametersJSON
export loadDislocationLoopJSON, loadNetworkJSON, loadIntegrationTimeJSON
export loadFEMParametersJSON
export saveJSON

end # module
