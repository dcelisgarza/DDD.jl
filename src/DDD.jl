module DDD

using CSV, DataFrames, LinearAlgebra, DelimitedFiles, Plots
import Base: zero, isequal, isless, convert, ==, *, /, length

include("CustomTypes.jl")
export compStruct, intAngle, extAngle

include("DislocationBase.jl")
export DislocationP, DislocationNetwork, nodeType, loopSides
export coordLbl, coordIdx, idxLabel, idxCond, dataCond
export AbstractDlnSegment,
       dlnNone,
       dlnEdge,
       dlnEdgeN,
       dlnScrew,
       dlnMixed,
       AbstractDlnStr,
       loopPrism,
       loopShear,
       loopMixed,
       makeSegment,
       makeLoop!,
       DislocationLoop

include("Material.jl")

include("DdFem.jl")

include("CustomIntegration.jl")

include("Input.jl")
export loadCSV, loadParams, cleanFieldDf

include("Output.jl")
export saveParams

include("Plotting.jl")
export plotNodes, plotNodes!

end # module
