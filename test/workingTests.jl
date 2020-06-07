using Revise, BenchmarkTools, Plots, LinearAlgebra
# using plotlyjs
# https://github.com/sglyon/ORCA.jl/issues/8#issuecomment-629049679
# https://stackoverflow.com/a/48509426

using DDD
cd(@__DIR__)
plotlyjs()

fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
fileSlipSystem = "../data/slipSystems/BCC.JSON"
fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
    fileDislocationP,
    fileMaterialP,
    fileIntegrationP,
    fileSlipSystem,
    fileDislocationLoop,
)
network = DislocationNetwork(dislocationLoop, memBuffer = 1)
DislocationNetwork!(network, dislocationLoop)

shearDecagon = DislocationLoop(
    loopShear();
    numSides = 10,
    nodeSide = 1,
    numLoops = 1,
    segLen = [300; 700; 1100; 1500; 1900; 1900; 1500; 1100; 700; 300],
    slipSystem = 4,
    _slipPlane = slipSystems.slipPlane[:, 4],
    _bVec = slipSystems.bVec[:, 4],
    label = nodeType[1; 1; 1; 1; 1; 1; 1; 1; 1; 1],
    buffer = 0.0,
    range = Float64[0 0; 0 0; 0 0],
    dist = Zeros(),
)

network = DislocationNetwork(shearDecagon)
network.coord[:, 11] = vec(mean(network2.coord, dims = 2))
network.label[11] = 1
network.numNode = 11
network.links[:, 11] = [11; 2]
network.links[:, 12] = [11; 4]
network.links[:, 13] = [11; 5]
network.links[:, 14] = [11; 10]
network.bVec[:, 11:14] .= network.bVec[:, 1]
network.slipPlane[:, 11:14] .= network.slipPlane[:, 1]
makeConnect!(network)
getSegmentIdx!(network)
plotlyjs()
fig1 = plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)
coarsenNetwork!(dlnParams, matParams, network)
fig1 = plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)
refineNetwork!(dlnParams, matParams, network)
fig1 = plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)
refineNetwork!(dlnParams, matParams, network)
fig1 = plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)
