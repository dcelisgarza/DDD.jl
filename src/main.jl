using DDD

# cd(@__DIR__)
inFilename = "../inputs/simParams/sampleParams.csv"
outFilename = "../outputs/simParams/sampleParams.csv"

dlnParams, matParams, intParams = loadParams(inFilename)
saveParams(dlnParams, matParams, intParams, outFilename; delim = ',')
