using Revise
using Plots, BenchmarkTools, LinearAlgebra, Interact
using DDD

# Define material parameters.
materialP = MaterialP(; μ = 1.0, μMag = 1.45e5, ν = 0.28, E = 1.0, crystalStruct = BCC())
# Define dislocation parameters.
dislocationP = DislocationP(;
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
integrationP = IntegrationP(;
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
slipSystem = SlipSystem(crystalStruct = BCC(), slipPlane = [1.0; 1.0; 1.0], bVec = [1.0; -1.0; 0.0])
# dot(a,b) == a ⋅ b == a' * b == 0
isapprox(slipSystem.slipPlane ⋅ slipSystem.bVec, 0) ==
isapprox(dot(slipSystem.slipPlane, slipSystem.bVec), 0) ==
isapprox(slipSystem.slipPlane' * slipSystem.bVec, 0)
# Define loops.
# Create initial structure.
prisPentagon = DislocationLoop(
    loopPrism();    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = 10 * ones(5),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystem.slipPlane,  # Slip plane of the segments.
    _bVec = slipSystem.bVec,            # Burgers vector of the segments.
    label = nodeType[1; 2; 1; 2; 1],    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to move away from the minimum limits of the range.
    range = Float64[-100 -100 -100; 100 100 100],   # Distribution range.
    dist = Rand(),  # Loop distribution.
)
shearHexagon = DislocationLoop(
    loopShear();    # Shear loop
    numSides = 6,
    nodeSide = 3,   # 3 nodes per side, it devides the side into equal segments.
    numLoops = 20,
    segLen = 10 * ones(3 * 6) / 3,  # The side length is 10, each segment is a third of 10/3.
    slipSystem = 1,
    _slipPlane = slipSystem.slipPlane,
    _bVec = slipSystem.bVec,
    label = nodeType[1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2],
    buffer = 0.0,
    range = Float64[-100 -100 -100; 100 100 100],
    dist = Rand(),
)
network = DislocationNetwork([shearHexagon, prisPentagon]; memBuffer = 1)

plotlyjs()
fig1 = plotNodes(
    shearHexagon,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
    size = (750, 750),
)
plotNodes!(fig1, prisPentagon, m = 1, l = 3, linecolor = :red, markercolor = :red, legend = false)
plot!(fig1, camera=(100,35))

fig2 = plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
    size = (750, 750),
)
plot!(fig2, camera=(110,40))
cd(@__DIR__)
savefig(fig1, "loops.svg")
savefig(fig2, "network.svg")
