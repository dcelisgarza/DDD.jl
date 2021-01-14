using DDD
matParams = 
##
using Revise
using Plots, BenchmarkTools, LinearAlgebra
using DDD

# Define material parameters.
MaterialParameters =
    MaterialParameters(; μ = 1.0, μMag = 1.45e5, ν = 0.28, E = 1.0, crystalStruct = BCC())
# Define dislocation parameters.
DislocationParameters = DislocationParameters(;
    coreRad = 90.78,
    coreRadMag = 1e-3,
    minSegLen = 300.0,
    maxSegLen = 1500.0,
    minArea = 30 * 300.0,
    maxArea = 20 * 30 * 300.0,
    maxConnect = 4,
    remesh = true,
    collision = true,
    separation = true,
    virtualRemesh = true,
    edgeDrag = 1.0,
    screwDrag = 2.0,
    climbDrag = 1e10,
    lineDrag = 0.0,
    mobility = mobBCC(),
)
# Define integration parameters
IntegrationParameters = IntegrationParameters(;
    dt = 1e3,
    tmin = 0.0,
    tmax = 1e10,
    method = CustomTrapezoid(),
    abstol = 1e-6,
    reltol = 1e-6,
    time = 0.0,
    step = 0,
)
# Define the slip system.
slipSystems =
    SlipSystem(crystalStruct = BCC(), slipPlane = [1.0; 1.0; 1.0], bVec = [1.0; -1.0; 0.0])
# dot(a,b) == a ⋅ b == a' * b == 0
isapprox(slipSystems.slipPlane ⋅ slipSystems.bVec, 0) ==
isapprox(dot(slipSystems.slipPlane, slipSystems.bVec), 0) ==
isapprox(slipSystems.slipPlane' * slipSystems.bVec, 0)
# Define loops.
# Create initial structure.
prisPentagon = DislocationLoop(
    loopPrism();    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = 10 * ones(5),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane,  # Slip plane of the segments.
    _bVec = slipSystems.bVec,            # Burgers vector of the segments.
    label = nodeType[1; 2; 1; 2; 1],    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = Float64[          # Distribution range
        -100 100 # xmin, xmax
        -100 100 # ymin, ymax
        -100 100  # zmin, zmax
    ],
    dist = Rand(),  # Loop distribution.
)
shearHexagon = DislocationLoop(
    loopShear();    # Shear loop
    numSides = 6,
    nodeSide = 3,   # 3 nodes per side, it devides the side into equal segments.
    numLoops = 20,
    segLen = 10 * ones(3 * 6) / 3,  # The hexagon's side length is 10, each segment is 10/3.
    slipSystem = 1,
    _slipPlane = slipSystems.slipPlane,
    _bVec = slipSystems.bVec,
    label = nodeType[1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2],
    buffer = 0.0,
    range = Float64[
        -100 100
        -100 100
        -100 100
    ],
    dist = Rand(),
)
network = DislocationNetwork([shearHexagon, prisPentagon], memBuffer = 1)

plotlyjs()
fig1 = plotNodes(
    shearHexagon,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)
plotNodes!(
    fig1,
    prisPentagon,
    m = 1,
    l = 3,
    linecolor = :red,
    markercolor = :red,
    legend = false,
)
plot!(fig1, camera = (100, 35), size = (400, 400))

fig2 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)
plot!(fig2, camera = (110, 40), size = (400, 400))
# cd(@__DIR__)
# savefig(fig1, "loops.png")
# savefig(fig2, "network.png")
