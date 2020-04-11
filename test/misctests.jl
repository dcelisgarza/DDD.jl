using Revise, BenchmarkTools
using DDD
cd(@__DIR__)

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

self = calcSelfForce(dlnParams, matParams, network)
@time calcSelfForce(dlnParams, matParams, network)
par = calcSegSegForce(dlnParams, matParams, network; parallel = true)
@time calcSegSegForce(dlnParams, matParams, network; parallel = true)

μ = matParams.μ
ν = matParams.ν
μ4π = matParams.μ4π
μ8π = matParams.μ8π
μ4πν = matParams.μ4πν
aSq = dlnParams.coreRadSq
μ8πaSq = aSq * μ8π
μ4πνaSq = aSq * μ4πν

idx = network.segIdx
coord = network.coord
bVec = network.bVec[idx[:, 1], :]
node1 = coord[idx[:, 2], :]
node2 = coord[idx[:, 3], :]

b1 = (bVec[1, 1], bVec[1, 2], bVec[1, 3])
n11 = (node1[1, 1], node1[1, 2], node1[1, 3])
n12 = (node2[1, 1], node2[1, 2], node2[1, 3])

b2 = (bVec[2, 1], bVec[2, 2], bVec[2, 3])
n21 = (node1[2, 1], node1[2, 2], node1[2, 3])
n22 = (node2[2, 1], node2[2, 2], node2[2, 3])


Fnode1, Fnode2, Fnode3, Fnode4 = calcSegSegForce(
    aSq,
    μ4π,
    μ8π,
    μ8πaSq,
    μ4πν,
    μ4πνaSq,
    b1,
    n11,
    n12,
    b2,
    n21,
    n22,
)

b1 = (bVec[1, 1], bVec[1, 2], bVec[1, 3])
n11 = (0.0, 0.0, 0.0)
n12 = (1.0, 0.0, 0.0)

b2 = (bVec[2, 1], bVec[2, 2], bVec[2, 3])
n21 = (0.0, 1.00000001, 0.000000001)
n22 = (2.0, 1.0, 0.0)


Fnode1, Fnode2, Fnode3, Fnode4 = calcParSegSegForce(
    aSq,
    μ4π,
    μ8π,
    μ8πaSq,
    μ4πν,
    μ4πνaSq,
    b1,
    n11,
    n12,
    b2,
    n21,
    n22,
)
