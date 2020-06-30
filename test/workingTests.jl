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
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)

using Serialization, BSON

@time serialize("test2.jls", (network, dlnParams, matParams, intParams, slipSystems, dislocationLoop))
@time save("test3.json", (network, dlnParams, matParams, intParams, slipSystems, dislocationLoop))
@time bson("test4.bson", a=(network, dlnParams, matParams, intParams, slipSystems, dislocationLoop))

networkOUT1, dlnParamsOUT1, matParamsOUT1 = open("test.txt", "r") do io
    deserialize(io)
end

networkOUT1
dlnParamsOUT1

networkOUT2 = deserialize("test2.txt")

compStruct(networkOUT1, network)
compStruct(networkOUT2, network)

## Remeshing and integration
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

function foo(intParams, intVars, dlnParams, matParams, network)
    coarsenNetwork!(dlnParams, matParams, network)
    refineNetwork!(dlnParams, matParams, network)
    integrate!(intParams, intVars, dlnParams, matParams, network)
end
network2 = deepcopy(network)
@allocated foo(intParams, intVars, dlnParams, matParams, network2)
@allocated foo(intParams, intVars, dlnParams, matParams, network2)

@allocated foo(intParams, intVars, dlnParams, matParams, network2)
@allocated foo(intParams, intVars, dlnParams, matParams, network2)
gr()
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
            # camera=(60,30),
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
        # println(network2.numNode)
    end every 10

    gif(anim, "uau4.gif")
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



prisPentagon = DislocationLoop(
          loopPrism();    # Prismatic loop, all segments are edge segments.
          numSides = 5,   # 5-sided loop.
          nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
          numLoops = 20,  # Number of loops of this type to generate when making a network.
          segLen = 500 * ones(5),  # Length of each segment between nodes, equal to the number of nodes.
          slipSystem = 2, # Slip System (assuming slip systems are stored in a file, this is the index).
          _slipPlane = slipSystems.slipPlane[:, 2],  # Slip plane of the segments.
          _bVec = slipSystems.bVec[:, 2],            # Burgers vector of the segments.
          label = nodeType[1; 1; 1; 1; 1],    # Node labels, has to be equal to the number of nodes.
          buffer = 0.0,   # Buffer to increase the dislocation spread.
          range = Float64[          # Distribution range
                        -5000 5000; # xmin, xmax
                        -5000 5000; # ymin, ymax
                        -5000 5000  # zmin, zmax
                      ],
          dist = Rand(),  # Loop distribution.
      )

prismHeptagon = DislocationLoop(
        loopPrism();    # Shear loop
        numSides = 7,
        nodeSide = 1,   # 3 nodes per side, it devides the side into equal segments.
        numLoops = 20,
        segLen = 700 * ones(7),  # The hexagon's side length is 10, each segment is 10/3.
        slipSystem = 1,
        _slipPlane = slipSystems.slipPlane[:, 1],
        _bVec = slipSystems.bVec[:, 1],
        label = nodeType[1; 1; 1; 1; 1; 2; 1],
        buffer = 0.0,
        range = Float64[
                      -5000 5000;
                      -5000 5000;
                      -5000 5000
                    ],
        dist = Rand(),
    )

network = DislocationNetwork(
          [prismHeptagon, prisPentagon]; # Dispatch type, bespoke functions dispatch on this.
          memBuffer = 1 # Buffer for memory allocation.
       )
gr()


fig = plotNodes(
        network,
        m = 3,
        l = 3,
        linecolor = :blue,
        markercolor = :blue,
        legend = false,
        # camera=(60,30),
        # size=(400,400)
      )

baar(intParams, intVars, dlnParams, matParams, network)
