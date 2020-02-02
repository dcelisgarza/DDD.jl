include("GlobalModule.jl")
using CSV, DataFrames
using .GlobalModule.CustomTypes
using .GlobalModule.DislocationBase
using .GlobalModule.MaterialBase
using .GlobalModule.CustomIntegration
using .GlobalModule.DdFemBase
using .GlobalModule.Input
using .GlobalModule.Output

inFilename = "C:/Users/Daniel Celis Garza/.julia/dev/DDD/inputs/simParams/sample"
outFilename = "C:/Users/Daniel Celis Garza/.julia/dev/DDD/outputs/simParams/sample"


dlnParams, matParams, intParams = loadParams(inFilename)
# df = loadParams(inFilename)
# dlnNetwork = genSourcesBCC(dlnParams)
saveParams(dlnParams, matParams, intParams, outFilename; delim = ',')
# df1 = CSV.read(inFilename*"Params"*".csv";header=0,copycols=true)
# df2 = CSV.read(outFilename*"Params"*".csv";header=0,copycols=false)


dlnParams2, matParams2, intParams2 = loadParams(outFilename)
