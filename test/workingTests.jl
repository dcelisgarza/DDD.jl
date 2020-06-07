using Revise, BenchmarkTools, Plots, LinearAlgebra
# using plotlyjs
# https://github.com/sglyon/ORCA.jl/issues/8#issuecomment-629049679
# https://stackoverflow.com/a/48509426

using DDD
cd(@__DIR__)
plotlyjs()

fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
fileSlipSystem = "../data/slipSystems/BCC.JSON"
fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
fileIntVar = "../inputs/simParams/sampleIntegrationTime.JSON"
dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
    fileDislocationP,
    fileMaterialP,
    fileIntegrationP,
    fileSlipSystem,
    fileDislocationLoop,
)
intVars = loadIntegrationVar(fileIntVar)
# network = DislocationNetwork(dislocationLoop, memBuffer = 1)
# DislocationNetwork!(network, dislocationLoop)

shearDecagon = DislocationLoop(
    loopShear();
    numSides = 10,
    nodeSide = 1,
    numLoops = 1,
    segLen = [300; 700; 1100; 1500; 1900; 1900; 1500; 1100; 700; 300],
    slipSystem = 4,
    _slipPlane = slipSystems.slipPlane[:, 4],
    _bVec = slipSystems.bVec[:, 4],
    label = nodeType[1; 1; 1; 1; 1; 1; 1; 1; 1; 1],
    buffer = 0.0,
    range = Float64[0 0; 0 0; 0 0],
    dist = Zeros(),
)

network = DislocationNetwork(shearDecagon)
network.coord[:, 11] = vec(mean(network.coord, dims = 2))
network.label[11] = 1
network.numNode = 11
network.links[:, 11] = [11; 2]
network.links[:, 12] = [11; 4]
network.links[:, 13] = [11; 5]
network.links[:, 14] = [11; 10]
network.bVec[:, 11:14] .= network.bVec[:, 1]
network.slipPlane[:, 11:14] .= network.slipPlane[:, 1]
makeConnect!(network)
getSegmentIdx!(network)
plotlyjs()
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)

network2 = deepcopy(network)
intVars2 = deepcopy(intVars)
numSeg = network.numSeg
intVars2
intParams
network2.coord[:, 1:numSeg] - network.coord[:, 1:numSeg]
fig2 = plotNodes(
    network2,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)


network2 = deepcopy(network)
refineNetwork!(dlnParams, matParams, network2)

foo(intParams, intVars, dlnParams, matParams, network2)
function foo(intParams, intVars, dlnParams, matParams, network)
    # coarsenNetwork!(dlnParams, matParams, network)
    # refineNetwork!(dlnParams, matParams, network)
    integrate!(intParams, intVars, dlnParams, matParams, network)
end
function baar(intParams, intVars, dlnParams, matParams, network)

    network2 = deepcopy(network)
    intVars2 = deepcopy(intVars)

    anim = @animate for i in 1:500


        fig = plotNodes(
            network2,
            m = 3,
            l = 2,
            linecolor = :blue,
            markercolor = :blue,
            legend = false,
        )
        foo(intParams, intVars2, dlnParams, matParams, network2)
        # if mod(i, 10) == 0
            # plotNodes!(
            #     fig,
            #     network2,
            #     m = 1,
            #     l = 3,
            #     linecolor = :blue,
            #     markercolor = :blue,
            #     legend = false,
            #     show=true
            # )
            # plot!(fig)
        # end
    end every 5
    gif(anim, "../examples/integ.gif")
end

baar(intParams, intVars, dlnParams, matParams, network)

# for i in 1:1000
#     integrate!(intParams, intVars2, dlnParams, matParams, network2)
#     plotNodes!(
#         fig2,
#         network2,
#         m = 1,
#         l = 3,
#         linecolor = :blue,
#         markercolor = :blue,
#         legend = false,
#     )
# end

coarsenNetwork!(dlnParams, matParams, network)
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)
refineNetwork!(dlnParams, matParams, network)
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)
refineNetwork!(dlnParams, matParams, network)
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)
