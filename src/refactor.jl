using DDD
# cd(@__DIR__)
fileDislocationParameters = "DDD/inputs/simParams/sampleDislocationParameters.JSON"
fileMaterialParameters = "DDD/inputs/simParams/sampleMaterialParameters.JSON"
fileIntegrationParameters = "DDD/inputs/simParams/sampleIntegrationParameters.JSON"
fileSlipSystem = "DDD/data/slipSystems/BCC.JSON"
fileDislocationLoop = "DDD/inputs/dln/samplePrismShear.JSON"
fileIntVar = "DDD/inputs/simParams/sampleIntegrationTime.JSON"
dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
    fileDislocationParameters,
    fileMaterialParameters,
    fileIntegrationParameters,
    fileSlipSystem,
    fileDislocationLoop,
)

network = DislocationNetwork(dislocationLoop)
network = DislocationNetwork!(network, dislocationLoop)
network.numSeg
network.numNode
network.coord