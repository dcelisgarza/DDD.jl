using DDD
using Test, Plots
using BenchmarkTools
import DelimitedFiles: readdlm
plotlyjs()
cd(@__DIR__)
params = "../inputs/simParams/sampleParams.csv"
slipsys = "../data/slipSystems/bcc.csv"
source = "../inputs/dln/sampleDlnRand.csv"
dlnParams, matParams, intParams, slipSystems, loops = loadParams(
    params,
    slipsys,
    source,
)
network = DislocationNetwork(
    zeros(Int64, 712, 2),
    zeros(712, 3),
    zeros(712, 3),
    zeros(712, 3),
    zeros(nodeType, 712),
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
