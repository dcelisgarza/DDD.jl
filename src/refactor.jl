using DDD
# cd(@__DIR__)
fileDislocationParameters = "DDD/inputs/simParams/sampleDislocationParameters.JSON"
fileMaterialParameters = "DDD/inputs/simParams/sampleMaterialParameters.JSON"
fileIntegrationParameters = "DDD/inputs/simParams/sampleIntegrationParameters.JSON"
fileSlipSystem = "DDD/data/slipSystems/BCC.JSON"
fileDislocationLoop = "DDD/inputs/dln/sampleDislocation.JSON"
fileIntVar = "DDD/inputs/simParams/sampleIntegrationTime.JSON"
# dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
#     fileDislocationParameters,
#     fileMaterialParameters,
#     fileIntegrationParameters,
#     fileSlipSystem,
#     fileDislocationLoop,
# )





slipSystem = SlipSystem(DDD.BCC(), [1; 0; 1], [-1; 0; 1])
# slipSystemName = SlipSystem(; crystalStruct = DDD.BCC(), slipPlane = [1; 0; 1], bVec = [-1; 0; 1])

# params = DislocationParameters(0.25, 0.25^2, 156., 11., 22., 11. * 2, 122., 444., 122.0^2, 444.0^2, 1.0, 1e-6, 1e6, 1e-9, 4, DDD.mobBCC(), true, true, true, true, true, true)
# paramsAuto = DislocationParameters(0.25, 156., 11., 22., 122., 444., 1.0, 1e-6, 1e6, 1e-9, 4, DDD.mobBCC())
# paramsName = DislocationParameters(; coreRad = 0.25, coreRadMag = 156., minSegLen = 11., maxSegLen = 22., minArea = 122., maxArea = 444., edgeDrag = 1.0, screwDrag = 1e-6, climbDrag = 1e6, lineDrag = 1e-9, maxConnect = 4, mobility = DDD.mobBCC())

# loop = DislocationLoop(DDD.loopJog(), 5, 2, 2, 100., 5, [1 2], [1 0 1], [-1 0 1], [0 0 0], [DDD.nodeType(1)], 0., [0 10; 0 10; 0 10], DDD.Zeros())


dictDislocationLoop = load(fileDislocationLoop)
dislocationLoop = loadDislocationLoop(dictDislocationLoop, slipSystem)

network = DislocationNetwork(dislocationLoop)
network = DislocationNetwork!(network, dislocationLoop)
network.numSeg
network.numNode
network.label