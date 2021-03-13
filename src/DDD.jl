module DDD

using LinearAlgebra,
    SparseArrays, SuiteSparse, Plots, Statistics, InteractiveUtils, JSON, StaticArrays, FileIO, LazySets, FastGaussQuadrature

# Miscelaneous.
include("./Misc/Misc.jl")
export compStruct, internalAngle, externalAngle, rot3D
export ⊗, linePlaneIntersect, gausslegendre, safeNorm

include("./Type/TypeBase.jl")
export AbstractDlnSeg,
    AbstractCrystalStruct,
    BCC,
    FCC,
    HCP,
    loopDln,
    loopPrism,
    loopShear,
    loopPure,
    loopImpure,
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
    ForceDisplacement,
    ForceDisplacementDot,
    AbstractMobility,
    mobBCC,
    mobFCC,
    mobHCP,
    AbstractDlnStr,
    AbstractDistribution,
    limits!,
    translatePoints!,
    makeNetwork!,
    makeSegment,
    AbstractElementOrder,
    AbstractIntegrator,
    BoundaryNode,
    Boundaries,
    AbstractModel,
    AbstractCantileverBend,
    CantileverLoad
export MaterialParameters
export nodeTypeDln,
    nodeTypeFE,
    SlipSystem,
    DislocationParameters,
    DislocationLoop,
    DislocationLoopCollection,
    DislocationNetwork,
    DislocationNetwork!,
    checkNetwork,
    getSegmentIdx,
    getSegmentIdx!,
    makeConnect,
    makeConnect!,
    integrate!
export IntegrationParameters, IntegrationTime, AbstractIntegrator, AdaptiveEulerTrapezoid
export Rand, Randn, Zeros, Regular, loopDistribution

include("./Processing/ProcessingBase.jl")
export calcSegForce,
    calcSegForce!,
    calcSelfForce,
    calcSelfForce!,
    calcSegSegForce,
    calcSegSegForce!,
    calc_σHat,
    calcPKForce,
    calcPKForce!,
    remeshSurfaceNetwork!,
    calc_σTilde,
    calc_σTilde!,
    findIntersectVolume,
    calc_uTilde,
    calc_uTilde!,
    calcDisplacementDislocationTriangle!
export dlnMobility, dlnMobility!
export findConnectedNode, mergeNode!, splitNode!, coarsenNetwork!, refineNetwork!, makeSurfaceNode!, coarsenVirtualNetwork!
export shapeFunction, shapeFunctionDeriv, deriv!

include("./PostProcessing/Plotting.jl")
export plotNodes, plotNodes!, plotFEDomain

include("./IO/IOBase.jl")
export loadJSON,
    loadDislocationParameters, loadMaterialParameters, loadIntegrationParameters, loadBoundaries, loadForceDisplacement
export loadSlipSystem, loadDislocationLoop, loadParameters
export loadDislocationLoop, loadNetwork, loadIntegrationTime
export loadFEMParameters
export saveJSON

end # module
