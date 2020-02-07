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

network2 = DislocationNetwork(
    zeros(Integer, 0, 2),
    zeros(0, 3),
    zeros(0, 3),
    zeros(0, 3),
    zeros(nodeType, 0),
    0,
    0,
)

slipsys = "../data/slipSystems/bcc.csv"
slipSystems = readdlm(slipsys, ',', Float64)
network = makeLoop!(loopPrism(), network, dlnParams, slipSystems)

network2 = makeLoop!(loopPrism(), network2, dlnParams, slipSystems)

intAngle(6) / π * 180
