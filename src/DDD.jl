module DDD

using CSV, DataFrames, LinearAlgebra, DelimitedFiles, Plots
import Base: zero, isequal, isless, convert, ==, *, /

include("CustomTypes.jl")
export compStruct, intAngle, extAngle

include("DislocationBase.jl")
export DislocationP, DislocationNetwork, nodeType, loopSides
export coordLbl, coordIdx, idxLabel, idxCond, dataCond
export AbstractDlnSegment,
       dlnEdge,
       dlnEdgeN,
       dlnScrew,
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
export plotNodes

end # module
