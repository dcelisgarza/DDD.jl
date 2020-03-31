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
include("./Dislocation/DislocationBase.jl")
include("./Material/MaterialBase.jl")
include("./IO/IOBase.jl")
include("./PostProcessing/Plotting.jl")


export inclusiveComparison, compStruct, intAngle, extAngle, rot3D, makeTypeDict
export nodeType, loopSides, AbstractDlnSeg, segNone, segEdge, segEdgeN, segScrew
export segMixed, AbstractDlnStr, loopPrism, loopShear, loopMixed, loopDln
export AbstractDistribution, Zeros, Rand, Randn, Regular
export AbstractCrystalStruct, BCC, FCC, HCP
export AbstractMobility, mobBCC, mobFCC, mobHCP
export AbstractIntegrator, CustomTrapezoid, AbstractMesh
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
