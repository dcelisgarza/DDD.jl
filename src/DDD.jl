module DDD

using CSV, DataFrames, LinearAlgebra, DelimitedFiles
import Base: zero, isequal, isless, convert, ==, vcat

include("CustomTypes.jl")
export compStruct, intAngle, extAngle

include("DislocationBase.jl")
export DislocationP, DislocationNetwork, nodeType
export coordLbl, coordIdx, idxLabel, idxCond, dataCond
export AbstractDlnSegment,
       dlnEdge,
       dlnScrew,
       AbstractDlnLoop,
       loopPrism,
       loopShear,
       loopMixed,
       makeSegment,
       makeLoop!

include("Material.jl")

include("DdFem.jl")

include("CustomIntegration.jl")

include("Input.jl")
export loadCSV, loadParams, cleanFieldDf

include("Output.jl")
export saveParams

end # module
