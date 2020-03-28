module DDD

using CSV,
    DataFrames,
    LinearAlgebra,
    DelimitedFiles,
    Plots,
    Statistics,
    InteractiveUtils
import Base:
    zero, isequal, isless, convert, ==, *, /, length, getindex, eachindex

include("Misc.jl")
export inclusiveComparison, compStruct, intAngle, extAngle, rot3D, makeTypeDict

include("PrimitiveTypes.jl")
export nodeType, loopSides, AbstractDlnSeg, segNone, segEdge, segEdgeN, segScrew
export segMixed, AbstractDlnStr, loopPrism, loopShear, loopMixed, loopDln
export AbstractDistribution, Zeros, Rand, Randn, Regular
export AbstractCrystalStruct, BCC, FCC, HCP
export AbstractMobility, mobBCC, mobFCC, mobHCP
export AbstractIntegrator, CustomTrapezoid

include("DislocationBase.jl")
export DislocationP, DislocationNetwork, malloc
export coordLbl, coordIdx, idxLabel, idxCond, dataCond
export makeSegment,
    makeLoop,
    DislocationLoop,
    makeNetwork,
    makeNetwork!,
    loopDistribution,
    makeConnect

include("Material.jl")

include("DdFem.jl")

include("CustomIntegration.jl")

include("Input.jl")
export loadCSV, loadParams, loadDln, loadSlipSys

include("Output.jl")
export saveParams

include("Plotting.jl")
export plotNodes, plotNodes!

end # module
