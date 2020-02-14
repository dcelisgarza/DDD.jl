using DDD
using Test, Plots

import DelimitedFiles: readdlm
plotlyjs()
cd(@__DIR__)
params = "../inputs/simParams/sampleParams.csv"
dlnParams, matParams, intParams = loadParams(params)
network = DislocationNetwork(
    zeros(Integer, 50, 2),
    zeros(50, 3),
    zeros(50, 3),
    zeros(50, 3),
    zeros(nodeType, 50),
    0,
    0,
)
slipsys = "../data/slipSystems/bcc.csv"
slipSystems = readdlm(slipsys, ',')
function dist(int::Integer)
    return rand(int, 3)
end


rang = Float64[-1 -1 -1; 1 1 1]
scale = Float64[1 1 1; 1 1 1] .* dlnParams.distSource * 10
makeLoop!(
    loopPrism(),
    network,
    dlnParams,
    slipSystems,
    rang,
    scale;
    dist = dist,
)
makeLoop!(
    loopShear(),
    network,
    dlnParams,
    slipSystems,
    rang,
    scale;
    dist = dist,
)
plot(
    network.coord[:, 1],
    network.coord[:, 2],
    network.coord[:, 3],
    m = 3,
    l = 3,
)
idxLabel(network, -1; condition = !==)
plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :red,
    markercolor = :red,
    legend = false,
)


# plot!(zeros(n), zeros(n), 1:n, w=10)

using Statistics, LinearAlgebra
abs(mean(network.coord)) < eps(Float64) * maximum(abs.(network.coord[:, :]))


network = DislocationNetwork(
    zeros(Integer, 8, 2),
    zeros(8, 3),
    zeros(8, 3),
    zeros(8, 3),
    zeros(nodeType, 8),
    0,
    0,
)
rang = Float64[0 0 0; 0 0 0]
scale = Float64[1 1 1; 1 1 1]
makeLoop!(loopShear(), network, dlnParams, slipSystems, rang, scale, 0.0)

fig = plot()
testSlip = Float64[ 1 1 1 -1  1 0;
                    1 1 1 -1  0 1;
                    1 1 1  0 -1 1]
shearLoop = DislocationLoop(
    [dlnEdge(), dlnScrew()],
    loopSides(4),
    2,
    [testSlip[1, 1:3]'; testSlip[1, 1:3]'],
    [testSlip[1, 4:6]'; testSlip[1, 4:6]'],
    nodeType[1; 0; 1; 1; 1; 0; 1; 1],
)
plotNodes!(
    fig,
    shearLoop,
    m = 1,
    l = 3,
    linecolor = :red,
    markercolor = :red,
    legend = false,
)
prismaticLoop = DislocationLoop(
    [dlnEdge(), dlnEdgeN()],
    loopSides(4),
    2,
    [testSlip[1, 1:3]'; testSlip[1, 1:3]'],
    [testSlip[1, 4:6]'; testSlip[1, 4:6]'],
    nodeType[1; 0; 1; 1; 1; 0; 1; 1],
)
plotNodes!(
    fig,
    prismaticLoop,
    m = 1,
    l = 3,
    linecolor = :orange,
    markercolor = :orange,
    legend = false,
)
prismaticLoop61 = DislocationLoop(
    [dlnEdge(), dlnEdge(), dlnEdge()],
    loopSides(6),
    3,
    testSlip[:,1:3],
    testSlip[:,4:6],
    nodeType[1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0],
)
plotNodes!(
    fig,
    prismaticLoop61,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)

prismaticLoop62 = DislocationLoop(
    [dlnEdgeN(), dlnEdge(), dlnEdge()],
    loopSides(6),
    3,
    [slipSystems[1, 1:3]'; slipSystems[2, 1:3]'; slipSystems[3, 1:3]'],
    [slipSystems[1, 4:6]'; slipSystems[2, 4:6]'; slipSystems[3, 4:6]'],
    nodeType[1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0],
)
plotNodes!(
    fig,
    prismaticLoop62,
    m = 1,
    l = 3,
    linecolor = :indigo,
    markercolor = :indigo,
    legend = false,
)

prismaticLoop62 = DislocationLoop(
    [dlnEdgeN(), dlnEdgeN(), dlnEdge()],
    loopSides(6),
    3,
    [slipSystems[1, 1:3]'; slipSystems[2, 1:3]'; slipSystems[3, 1:3]'],
    [slipSystems[1, 4:6]'; slipSystems[2, 4:6]'; slipSystems[3, 4:6]'],
    nodeType[1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0],
)
plotNodes!(
    fig,
    prismaticLoop62,
    m = 1,
    l = 3,
    linecolor = :purple,
    markercolor = :purple,
    legend = false,
)

prismaticLoop63 = DislocationLoop(
    [dlnEdgeN(), dlnEdgeN(), dlnEdgeN()],
    loopSides(6),
    3,
    [slipSystems[1, 1:3]'; slipSystems[2, 1:3]'; slipSystems[3, 1:3]'],
    [slipSystems[1, 4:6]'; slipSystems[2, 4:6]'; slipSystems[3, 4:6]'],
    nodeType[1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0],
)
plotNodes!(
    fig,
    prismaticLoop63,
    m = 1,
    l = 3,
    linecolor = :cyan,
    markercolor = :cyan,
    legend = false,
)
