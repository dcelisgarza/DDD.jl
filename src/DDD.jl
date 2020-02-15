module DDD

using CSV, DataFrames, LinearAlgebra, DelimitedFiles, Plots, Statistics
import Base: zero, isequal, isless, convert, ==, *, /, length

include("Misc.jl")
export compStruct, intAngle, extAngle

include("PrimitiveTypes.jl")
export nodeType,
       loopSides,
       AbstractDlnSegment,
       segNone,
       segEdge,
       segEdgeN,
       segScrew,
       segMixed,
       AbstractDlnStr,
       loopPrism,
       loopShear,
       loopMixed,
       AbstractDistribution,
       Zeros,
       Rand,
       Randn,
       Regular,
       AbstractMaterial,
       BCC,
       FCC,
       HCP,
       AbstractMobility,
       mobBCC,
       mobFCC,
       mobHCP,
       AbstractIntegrator,
       CustomTrapezoid

include("DislocationBase.jl")
export DislocationP, DislocationNetwork, malloc
export coordLbl, coordIdx, idxLabel, idxCond, dataCond
export makeSegment, makeLoop!, DislocationLoop, makeNetwork!, loopDistribution

include("Material.jl")

include("DdFem.jl")

include("CustomIntegration.jl")

include("Input.jl")
export loadCSV, loadParams, loadDln

include("Output.jl")
export saveParams

include("Plotting.jl")
export plotNodes, plotNodes!

end # module
