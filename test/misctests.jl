using DDD
using Test, Plots
using BenchmarkTools
import DelimitedFiles: readdlm
plotlyjs()
cd(@__DIR__)
params = "../inputs/simParams/sampleParams.csv"
slipsys = "../data/slipSystems/bcc.csv"
source = "../inputs/dln/samplePrismatic.csv"
dlnParams, matParams, intParams, slipSystems, loops = loadParams(
    params,
    slipsys,
    source,
)
network = DislocationNetwork(
    zeros(Int64, 0, 2),
    zeros(0, 3),
    zeros(0, 3),
    zeros(0, 3),
    zeros(nodeType, 0),
    0,
    0,
)
makeNetwork!(network, loops)
fig = plot()
plotNodes!(
    fig,
    network,
    m = 1,
    l = 3,
    linecolor = :black,
    markercolor = :black,
    legend = false,
)

fig = plot()
plot(loops[1].coord[:,1], loops[1].coord[:,2], loops[1].coord[:,3], m=3, l=3)
using Statistics
mean(loops[1].coord, dims=1)
plotNodes!(
    fig,
    loops[1],
    m = 1,
    l = 3,
    linecolor = :black,
    markercolor = :black,
    legend = false,
)
plotNodes!(
    fig,
    loops[2],
    m = 1,
    l = 3,
    linecolor = :red,
    markercolor = :red,
    legend = false,
)
plotNodes!(
    fig,
    loops[3],
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)
plotNodes!(
    fig,
    loops[4],
    m = 1,
    l = 3,
    linecolor = :orange,
    markercolor = :orange,
    legend = false,
)

loop = DislocationLoop(
    loopDln(),
    4,
    convert(Int64, 1),
    convert(Int64, 1),
    [segNone(), segNone()],
    zeros(Float64, 2),
    zeros(Int64, 2),
    zeros(Float64, 2, 3),
    zeros(Float64, 2, 3),
    zeros(nodeType, 4),
    convert(Float64, 0),
    zeros(2, 3),
    Zeros(),
)
network = DislocationNetwork(
    zeros(Int64, 0, 2),
    zeros(0, 3),
    zeros(0, 3),
    zeros(0, 3),
    zeros(nodeType, 0),
    0,
    0,
)
makeNetwork!(network, loop)
connectivity, linksConnect = makeConnect(network, dlnParams)
