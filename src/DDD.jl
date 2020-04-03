module DDD

using CSV,
    DataFrames,
    LinearAlgebra,
    DelimitedFiles,
    Plots,
    Statistics,
    InteractiveUtils
import Base:
    zero, isequal, isless, convert, ==, *, /, length, getindex, eachindex, push!

include("./Misc/Misc.jl")
include("./Integration/CustomIntegration.jl")
include("./Material/MaterialBase.jl")
include("./Dislocation/DislocationBase.jl")
include("./FEM/FEMBase.jl")
include("./DislocationFEM/DislocationFEMBase.jl")
include("./IO/IOBase.jl")
include("./PostProcessing/Plotting.jl")
export shapeFunction, shapeFunctionDeriv
export inclusiveComparison
export compStruct, intAngle, extAngle, rot3D, makeTypeDict
export getSegVector
export nodeType, loopSides, AbstractDlnSeg, segNone, segEdge, segEdgeN, segScrew
export segMixed, AbstractDlnStr, loopPrism, loopShear, loopMixed, loopDln
export AbstractDistribution, Zeros, Rand, Randn, Regular
export AbstractCrystalStruct, BCC, FCC, HCP
export AbstractMobility, mobBCC, mobFCC, mobHCP
export AbstractIntegrator, CustomTrapezoid
export AbstractMesh, RegularCuboidMesh, AbstractShapeFunction
export AbstractShapeFunction3D, AbstractShapeFunction2D
export LinearQuadrangle3D, LinearQuadrangle2D



export hatStress
export DislocationP, DislocationNetwork
export coordLbl, coordIdx, idxLabel, idxCond, dataCond
export makeSegment,
    makeLoop,
    DislocationLoop,
    makeNetwork,
    makeNetwork!,
    loopDistribution,
    makeConnect,
    checkNetwork,
    getSegmentIdx!
export loadCSV, loadParams, loadDln, loadSlipSys

export saveParams

export plotNodes, plotNodes!

end # module
