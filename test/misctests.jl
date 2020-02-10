using DDD
using Test, DelimitedFiles

cd(@__DIR__)
params = "../inputs/simParams/sampleParams.csv"
dlnParams, matParams, intParams = loadParams(params)
network = DislocationNetwork(
    zeros(Integer, 0, 2),
    zeros(0, 3),
    zeros(0, 3),
    zeros(0, 3),
    zeros(nodeType, 0),
    0,
    0,
)

slipsys = "../data/slipSystems/bcc.csv"
slipSystems = readdlm(slipsys, ',')
network = makeLoop!(loopPrism(), network, dlnParams, slipSystems)

using Plots
plotlyjs()
n = 100
ts = range(0, stop=8π, length=n)
x = ts .* map(cos, ts)
y = (0.1ts) .* map(sin, ts)
z = 1:n
plot(x, y, z, zcolor=reverse(z), m=(10, 0.8, :blues, Plots.stroke(0)), leg=false, cbar=true, w=5)
plot!(zeros(n), zeros(n), 1:n, w=10)

plot(network.coord[:,1],network.coord[:,2],network.coord[:,3],m=3,l=3)
