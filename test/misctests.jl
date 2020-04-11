using Revise
using DDD
cd(@__DIR__)

fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
fileSlipSystem = "../data/slipSystems/SlipSystems.JSON"
fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
    fileDislocationP,
    fileMaterialP,
    fileIntegrationP,
    fileSlipSystem,
    fileDislocationLoop,
)
network = makeNetwork(dislocationLoop; memBuffer = 1)

pentagon = DislocationLoop(
    loopType = loopPrism(),
    numSides = 5,
    nodeSide = 1,
    numLoops = 10,
    segLen = 10 * ones(5),
    slipSystem = 3,
    _slipPlane = slipSystems.slipPlane[3, :],
    _bVec = slipSystems.bVec[3, :],
    label = nodeType[0; 0; 0; 0; 0],
    buffer = 0.0,
    range = Float64[-100 -100 -100; 100 100 100],
    dist = Rand(),
)

makeNetwork!(network, pentagon)
using Makie

function plotNodesMakie(network::DislocationNetwork, args...; kw...)
    idx = findall(x -> x != -1, network.label)
    coord = network.coord
    fig = Scene()
    meshscatter!(coord, args...; kw...)
    for i in idx
        n1 = network.links[i, 1]
        n2 = network.links[i, 2]
        lines!(
            coord[[n1, n2], 1],
            coord[[n1, n2], 2],
            coord[[n1, n2], 3],
            args...;
            kw...,
        )
    end
    return fig
end

fig = plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :black,
    markercolor = :black,
    legend = false,
)

scene1 = plotNodesMakie(
    network,
    linewidth = 2,
    markersize = 0.5,
    strokecolor = :black,
    color = :black,
)

# self = calcSelfForce(dlnParams, matParams, network)
# @btime calcSelfForce(dlnParams, matParams, network)
# par = calcSegSegForce(dlnParams, matParams, network; parallel = true)
# @btime calcSegSegForce(dlnParams, matParams, network; parallel = true)

fileDump = "../outputs/dump/sampleDump.JSON"
save(fileDump, dlnParams, matParams, intParams, slipSystems, dislocationLoop)
simulation = load(fileDump)
dlnParams2 = loadDislocationP(simulation[1])
matParams2 = loadMaterialP(simulation[2])
intParams2 = loadIntegrationP(simulation[3])
slipSystems2 = loadSlipSystem(simulation[4])
dislocationLoop2 = zeros(DislocationLoop, length(simulation[5]))
for i in eachindex(dislocationLoop2)
    dislocationLoop2[i] = loadDislocationLoop(simulation[5][i], slipSystems2)
end
compStruct(dlnParams, dlnParams2)
compStruct(matParams, matParams2)
compStruct(intParams, intParams2)
compStruct(slipSystems, slipSystems2)
compStruct(dislocationLoop, dislocationLoop2)

fileDislocationNetwork = "../outputs/dln/sampleNetwork.JSON"
save(fileDislocationNetwork, network)
networkDict = load(fileDislocationNetwork)

network2 = loadNetwork(fileDislocationNetwork)
compStruct(network, network2; verbose=true)
network
















network2
