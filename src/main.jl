using DDD

cd(@__DIR__)
inFileParams = "../inputs/simParams/sampleParams.csv"
inFileSlipSys = "../data/slipSystems/bcc.csv"
inFileDln = "../inputs/dln/sampleDln.csv"

outFileParams = "../outputs/simParams/sampleParams.csv"

dlnParams, matParams, intParams, slipSystems, sources =
    loadParams(inFileParams, inFileSlipSys, inFileDln)
saveParams(dlnParams, matParams, intParams, outFileParams; delim = ',')
