module DDD

using CSV, DataFrames, LinearAlgebra, DelimitedFiles
import Base: zero, isequal, isless, convert, ==

include("CustomTypes.jl")
export compStruct

include("DislocationBase.jl")
export DislocationP, DislocationNetwork, nodeType
export coordLbl, coordIdx, idxLabel, idxCond, dataCond
export dlnSegment, dlnEdge, dlnScrew, makeSegment

include("Material.jl")

include("DdFem.jl")

include("CustomIntegration.jl")

include("Input.jl")
export loadCSV, loadParams

include("Output.jl")
export saveParams

end # module
