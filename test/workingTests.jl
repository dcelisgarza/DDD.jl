using Revise, BenchmarkTools, Plots
# using plotlyjs
# https://github.com/sglyon/ORCA.jl/issues/8#issuecomment-629049679
# https://stackoverflow.com/a/48509426

using DDD
cd(@__DIR__)

# using Makie
# function plotNodesMakie(network::DislocationNetwork, args...; kw...)
#     idx = findall(x -> x != 0, network.label)
#     coord = network.coord
#     fig = Scene()
#     meshscatter!(coord, args...; kw...)
#     for i in idx
#         n1 = network.links[i, 1]
#         n2 = network.links[i, 2]
#         lines!(
#             coord[[n1, n2], 1],
#             coord[[n1, n2], 2],
#             coord[[n1, n2], 3],
#             args...;
#             kw...,
#         )
#     end
#     return fig
# end

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
dislocationLoop[1]
network = DislocationNetwork(dislocationLoop, memBuffer = 1)
DislocationNetwork!(network, dislocationLoop)

pentagon = DislocationLoop(
    loopPrism();
    numSides = 5,
    nodeSide = 1,
    numLoops = 2,
    segLen = 10 * ones(5),
    slipSystem = 4,
    _slipPlane = slipSystems.slipPlane[:, 4],
    _bVec = slipSystems.bVec[:, 4],
    label = nodeType[1; 2; 1; 2; 1],
    buffer = 0.0,
    range = Float64[-100 100; -100 100; -100 100],
    dist = Rand(),
)

using LinearAlgebra
mean(pentagon.coord, dims = 2)
network = DislocationNetwork(pentagon; memBuffer = 1)
# DislocationNetwork!(network, pentagon)
network.coord = [
    -38.9526 -15.6253 -70.0006
    -44.416 -7.63875 -72.5237
    -51.8751 -6.78528 -65.9181
    -51.0216 -14.2443 -59.3125
    -43.0351 -19.7078 -61.8356
    -20.9647 26.5869 85.4687
    -26.4281 34.5734 82.9456
    -33.8872 35.4269 89.5512
    -33.0337 27.9679 96.1568
    -25.0472 22.5044 93.6336
];
network2 = deepcopy(network)
mergeNode!(network2, 1, 3)

@time mergeNode!(network2, 1, 10)

@allocated mergeNode!(network2, 1, 3)

@allocated mergeNode!(network2, 1, 3)

@allocated mergeNode!(network2, 1, 3)

using Plots
plotlyjs()
fig = plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :red,
    markercolor = :red,
    legend = false,
    size = (750, 750),
)
plotNodes!(
    fig,
    network2,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)

plotNodes!(
    fig,
    network1,
    m = 1,
    l = 3,
    linecolor = :red,
    markercolor = :red,
    legend = false,
)

scene1 = plotNodesMakie(
    network,
    linewidth = 2,
    markersize = 0.5,
    strokecolor = :black,
    color = :black,
)
using BenchmarkTools
self = calcSelfForce(dlnParams, matParams, network)
@btime calcSelfForce(dlnParams, matParams, network)

idx = [2, 5, 7]
selfIdx = calcSelfForce(dlnParams, matParams, network, idx)
isapprox(self[1][:, idx], selfIdx[1])

self[2][:, idx] == selfIdx[2]
@btime calcSelfForce(dlnParams, matParams, network, idx)

@allocated calcSelfForce(dlnParams, matParams, network)
@btime calcSelfForce(dlnParams, matParams, network)

ser = calcSegSegForce(dlnParams, matParams, network; parallel = false)
@btime calcSegSegForce(dlnParams, matParams, network; parallel = false)
@allocated calcSegSegForce(dlnParams, matParams, network; parallel = false)

ser = calcSegSegForce(dlnParams, matParams, network; parallel = false)
@btime calcSegSegForce(dlnParams, matParams, network; parallel = false)

ser[:, 2, idx]
serIdx[:, 2, :]

@test ser[2][:, idx] == serIdx[2]

idx = rand(1:(network.numNode))
serIdx = calcSegSegForce(dlnParams, matParams, network, idx; parallel = false)
@btime calcSegSegForce(dlnParams, matParams, network, idx; parallel = false)
@test ser[1][:, idx] == serIdx[1]
@test ser[2][:, idx] == serIdx[2]

par = calcSegSegForce(dlnParams, matParams, network; parallel = true)
@allocated calcSegSegForce(dlnParams, matParams, network; parallel = true)
@btime calcSegSegForce(dlnParams, matParams, network; parallel = true)

tser = calcSegForce(dlnParams, matParams, network; parallel = false)
@allocated calcSegForce(dlnParams, matParams, network; parallel = false)
@btime calcSegForce(dlnParams, matParams, network; parallel = false)

tpar = calcSegForce(dlnParams, matParams, network; parallel = true)
@allocated calcSegForce(dlnParams, matParams, network; parallel = true)
@btime calcSegForce(dlnParams, matParams, network; parallel = true)

isapprox.(par, ser)
isapprox.(tser, tpar)
isapprox.(tser .- ser, self)
isapprox.(tpar .- ser, self)

@btime calcSelfForce(dlnParams, matParams, network)
@btime calcSegSegForce(dlnParams, matParams, network; parallel = false)
@btime calcSegForce(dlnParams, matParams, network; parallel = false)

remo .- tot
@btime calcSegSegForce(dlnParams, matParams, network; parallel = false)

par = calcSegSegForce(dlnParams, matParams, network; parallel = true)
ser = calcSegSegForce(dlnParams, matParams, network; parallel = false)

maximum.(par .- ser)
minimum.(par .- ser)

idx = network.segIdx
coord = network.coord
bVec = network.bVec[idx[:, 1], :]
node1 = coord[idx[:, 2], :]
node2 = coord[idx[:, 3], :]

@btime tuple = (coord[idx[:, 2], :], coord[idx[:, 3], :])
@btime mean((node1, node2))
@btime 0.5 * (node1 + node2)

test = [1 1 1; 2 2 2; 3 3 3]
test2 = [1; 2; 3]
test[1, :]
test * test2

test .* test2
cross(vec(sum(test .* test2, dims = 1)), test2)

a = tot[1][:, :]
b = tot[2][:, :]
@btime 0.5 * (a + b)
@btime 0.5 * (tot[1][:, :] + tot[2][:, :])
@btime test = (tot[1][:, :], tot[2][:, :])
@btime mean(test)
@btime mean(tot[:])

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

Fnode1, Fnode2, Fnode3, Fnode4 =
    calcSegSegForce(aSq, μ4π, μ8π, μ8πaSq, μ4πν, μ4πνaSq, b1, n11, n12, b2, n21, n22)

b1 = (bVec[1, 1], bVec[1, 2], bVec[1, 3])
n11 = (0.0, 0.0, 0.0)
n12 = (1.0, 0.0, 0.0)

b2 = (bVec[2, 1], bVec[2, 2], bVec[2, 3])
n21 = (0.0, 1.00000001, 0.000000001)
n22 = (2.0, 1.0, 0.0)

Fnode1, Fnode2, Fnode3, Fnode4 =
    calcParSegSegForce(aSq, μ4π, μ8π, μ8πaSq, μ4πν, μ4πνaSq, b1, n11, n12, b2, n21, n22)

x = rand(10000000, 3)'
xt = copy(transpose(x))
y = transpose(rand(10000000, 3))
yt = copy(transpose(y))
function foo(x)
    for i in 1:size(x, 1)
        (x[i, :] .* i .- i) .^ i
    end
end
@time foo(x)
@time foo(xt)
@time foo(y)
@time foo(yt)
