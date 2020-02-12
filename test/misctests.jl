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
