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
# function dist(int::Integer)
#     return rand(int, 3)
# end
#
#
# rang = Float64[-1 -1 -1; 1 1 1]
# scale = Float64[1 1 1; 1 1 1] .* dlnParams.distSource * 10
# makeLoop!(
#     loopPrism(),
#     network,
#     dlnParams,
#     slipSystems,
#     rang,
#     scale;
#     dist = dist,
# )
# makeLoop!(
#     loopShear(),
#     network,
#     dlnParams,
#     slipSystems,
#     rang,
#     scale;
#     dist = dist,
# )
# plot(
#     network.coord[:, 1],
#     network.coord[:, 2],
#     network.coord[:, 3],
#     m = 3,
#     l = 3,
# )
# idxLabel(network, -1; condition = !==)
# plotNodes(
#     network,
#     m = 1,
#     l = 3,
#     linecolor = :red,
#     markercolor = :red,
#     legend = false,
# )
#
#
# # plot!(zeros(n), zeros(n), 1:n, w=10)
#
# using Statistics, LinearAlgebra
# abs(mean(network.coord)) < eps(Float64) * maximum(abs.(network.coord[:, :]))
#
#
# network = DislocationNetwork(
#     zeros(Integer, 8, 2),
#     zeros(8, 3),
#     zeros(8, 3),
#     zeros(8, 3),
#     zeros(nodeType, 8),
#     0,
#     0,
# )
# rang = Float64[0 0 0; 0 0 0]
# scale = Float64[1 1 1; 1 1 1]
# makeLoop!(loopShear(), network, dlnParams, slipSystems, rang, scale, 0.0)
zeroloop = zeros(DislocationLoop, 2)




fig = plot()
testSlip = Float64[
    1 1 1 -1 1 0
    1 1 1 -1 0 1
    1 1 1 0 -1 1
]
shearLoop = DislocationLoop(
    loopSides(4),
    2,
    [segEdge(), segScrew()],
    [1.0, 1],
    [testSlip[1, 1:3]'; testSlip[1, 1:3]'],
    [testSlip[1, 4:6]'; testSlip[1, 4:6]'],
    nodeType[1; 0; 1; 1; 1; 0; 1; 1],
    1,
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
    loopSides(4),
    2,
    [segEdge(), segEdgeN()],
    [1.0, 1.0],
    [testSlip[1, 1:3]'; testSlip[1, 1:3]'],
    [testSlip[1, 4:6]'; testSlip[1, 4:6]'],
    nodeType[1; 0; 1; 1; 1; 0; 1; 1],
    1,
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
    loopSides(6),
    3,
    [segEdge(), segEdge(), segEdge()],
    [1.0, 1, 5],
    testSlip[:, 1:3],
    testSlip[:, 4:6],
    nodeType[1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0],
    1,
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
    loopSides(6),
    3,
    [segEdgeN(), segScrew(), segEdge()],
    [1.0, 1, 1],
    [slipSystems[1, 1:3]'; slipSystems[2, 1:3]'; slipSystems[3, 1:3]'],
    [slipSystems[1, 4:6]'; slipSystems[2, 4:6]'; slipSystems[3, 4:6]'],
    nodeType[1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0],
    1,
)
plotNodes!(
    fig,
    prismaticLoop62,
    m = 1,
    l = 3,
    linecolor = :brown,
    markercolor = :brown,
    legend = false,
)

prismaticLoop62 = DislocationLoop(
    loopSides(6),
    3,
    [segEdgeN(), segEdgeN(), segEdge()],
    [1.0, 1, 1],
    [slipSystems[1, 1:3]'; slipSystems[2, 1:3]'; slipSystems[3, 1:3]'],
    [slipSystems[1, 4:6]'; slipSystems[2, 4:6]'; slipSystems[3, 4:6]'],
    nodeType[1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0],
    1,
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

using DataFrames

slipsys = "../data/slipSystems/bcc.csv"
slipSystems = readdlm(slipsys, ',')
df = loadCSV("../inputs/dln/sampleDln.csv"; header = 1, transpose = true)
difLoops = nrow(df)
loops = zeros(DislocationLoop, difLoops)
dict = Dict(
    "segEdge" => segEdge(),
    "segEdgeN" => segEdgeN(),
    "segScrew" => segScrew(),
    "segMixed" => segMixed(),
)
for i = 1:difLoops
    st = split.(df[i, :segType], ";")
    segType = [dict[st[i]] for i = 1:length(st)]
    sl = split.(df[i, :segLen], ";")
    segLen = parse.(Float64, sl)
    ss = split.(df[i, :slipSystem], ";")
    slipSystem = parse.(Int, ss)
    lbl = split.(df[i,:label],";")
    label = convert.(nodeType,parse.(Int, lbl))
    loops[i] = DislocationLoop(
        loopSides(df[i, :numSides]),
        df[i, :nodeSide],
        segType,
        segLen,
        slipSystems[slipSystem, 1:3],
        slipSystems[slipSystem, 4:6],
        label,
        df[i,:numLoops]
    )
end
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

label = convert.(nodeType,label)

ss = split.(df[1, :slipSystem], ";")
ss = parse.(Int, ss)
slipSystems[ss,1:3]

stype = split.(df[1, :segType], ";")

stype = [dict[stype[i]] for i = 1:length(stype)]
