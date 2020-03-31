using DDD
cd(@__DIR__)
params = "../inputs/simParams/sampleParams.csv"
slipsys = "../data/slipSystems/bcc.csv"
source = "../inputs/dln/samplePrismaticShear.csv"
dlnParams, matParams, intParams, slipSystems, loops =
    loadParams(params, slipsys, source)

network = DislocationNetwork(
    zeros(Int64, 3, 2),
    zeros(3, 3),
    zeros(3, 3),
    zeros(3, 3),
    zeros(nodeType, 3),
    0,
    0,
)
makeNetwork!(network, loops)
network2 = makeNetwork(loops,4,1)
compStruct(network, network2; verbose = false)


using Test, Plots
gr()
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
plot(
    loops[1].coord[:, 1],
    loops[1].coord[:, 2],
    loops[1].coord[:, 3],
    m = 3,
    l = 3,
)
using Statistics
mean(loops[1].coord, dims = 1)
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
# connectivity, linksConnect = makeConnect(network, dlnParams)

include("../src/PlottingMakie.jl")

plotNodesMakie(network; markersize = 0.6, linewidth = 3)

fig = Scene()
plotNodesMakie!(fig, network; markersize = 0.6, linewidth = 3)

# dlnSegTypes = subtypes(AbstractDlnSeg)
function makeTypeDict(valType::DataType)
    subTypes = subtypes(valType)
    dict = Dict{String, Any}()

    for i in eachindex(subTypes)
        push!(dict, string(subTypes[i]()) => eval(subTypes[i]()))
    end

    return dict
end

dict = makeTypeDict(AbstractDlnSeg)
test = subtypes(AbstractDlnSeg)
push!(dict, string(test[1]()) => eval(test[1]()))

dlnTypes = Dict(
    "loopPrism()" => loopPrism(),
    "loopShear()" => loopShear(),
    "loopMixed()" => loopMixed(),
    "loopDln()" => loopDln(),
    "DDD.loopPrism()" => loopPrism(),
    "DDD.loopShear()" => loopShear(),
    "DDD.loopMixed()" => loopMixed(),
    "DDD.loopDln()" => loopDln(),
)
