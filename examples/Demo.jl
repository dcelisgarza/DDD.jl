using DDD, Plots, LinearAlgebra, Statistics, SparseArrays, StaticArrays, BenchmarkTools
plotlyjs()

## Parameters
dlnParams = DislocationParameters(;
    mobility = mobBCC(),
    dragCoeffs = (edge = 1, screw = 1e-1, climb = 1e9, line = 1e-5),
    coreRad = 0.015 * sqrt(2),
    minSegLen = 0.15 * sqrt(2),
    maxSegLen = 1.5 * sqrt(2),
    coreEnergy = 1 / (4 * π) * log(0.015 * sqrt(2) / 3.5e-5),
    coreRadMag = 3.5e-4,
);

matParams = MaterialParameters(; crystalStruct = BCC(), μ = 1, μMag = 80e3, ν = 0.25);

femParamsC = FEMParameters(;
    type = DispatchRegularCuboidMesh(),
    order = LinearElement(),
    model = CantileverLoad(),
    dx = 23.0,
    dy = 17.0,
    dz = 13.0,
    mx = 7,
    my = 5,
    mz = 3,
);

slipSystems = SlipSystem(;
    crystalStruct = BCC(),
    slipPlane = Float64[-1 1; 1 -1; 0 0],
    bVec = Float64[1 1; 1 1; 1 -1],
);

intParams = IntegrationParameters(;
    method = AdaptiveEulerTrapezoid(),
    abstol = dlnParams.collisionDist / 2,
    reltol = dlnParams.collisionDist / 2,
);

intTime = IntegrationTime(; dt = 0.0, time = 0.0);

# See fieldnames
fieldnames(typeof(matParams))

# # Browse docs/get variable info
# ?DislocationParameters

# ?femParamsC

## Build mesh & boundary conditions
meshC = buildMesh(matParams, femParamsC)
# Plot mesh
figFE = plotFEDomain(meshC; camera = (10, -2))
# savefig(figFE, "mesh.svg")

# Generate cantilever boundary conditions.
boundaryC, forceDispC = Boundaries(femParamsC, meshC)
# Plot boundaries
figCBC = plotBoundaries(boundaryC, meshC)
# savefig(figCBC, "cantilever.svg")

# Pillar loading
femParamsP = FEMParameters(;
    type = DispatchRegularCuboidMesh(),
    order = LinearElement(),
    model = PillarLoad(),
    dx = 17.0,
    dy = 13.0,
    dz = 23.0,
    mx = 5,
    my = 3,
    mz = 7,
);
meshP = buildMesh(matParams, femParamsP);
boundaryP, forceDispP = Boundaries(femParamsP, meshP);
figPBC = plotBoundaries(boundaryP, meshP)
# savefig(figPBC, "pillar.svg")

## Dislocations generation
dx, dy, dz = femParamsC.dx, femParamsC.dy, femParamsC.dz
segLen = (dlnParams.minSegLen + dlnParams.maxSegLen) / 2

prismOct = DislocationLoop(;
    loopType = loopPrism(),
    numSides = 8,
    nodeSide = 1,
    numLoops = 5,
    segLen = segLen * ones(8),
    slipSystemIdx = 1,
    slipSystem = slipSystems,
    label = nodeTypeDln.(ones(Int, 8)),
    buffer = 0,
    range = [0 dx; 0 dy; 0 dz],
    dist = Rand(),
)
shearPent = DislocationLoop(;
    loopType = loopShear(),
    numSides = 5,
    nodeSide = 2,
    numLoops = 5,
    segLen = segLen * ones(10),
    slipSystemIdx = 2,
    slipSystem = slipSystems,
    label = nodeTypeDln.(ones(Int, 10)),
    buffer = 0,
    range = [0 dx; 0 dy; 0 dz],
    dist = Rand(),
)
network = DislocationNetwork((prismOct, shearPent))

dlnFig = plotNodes(
    prismOct,
    m = 2,
    l = 3,
    linecolor = :blue,
    markershape = :circle,
    markercolor = :blue,
    legend = false,
)
plotNodes!(
    dlnFig,
    shearPent,
    m = 2,
    l = 3,
    linecolor = :red,
    markershape = :square,
    markercolor = :red,
    legend = false,
    camera = (10, 20),
)
# savefig(dlnFig, "loops.svg")

networkFig = plotNodes(
    network,
    m = 2,
    l = 3,
    linecolor = :blue,
    markershape = :circle,
    markercolor = :blue,
    legend = false,
)

networkFEFig = plotNodes(
    meshC,
    network,
    m = 2,
    l = 3,
    linecolor = :blue,
    markershape = :circle,
    markercolor = :blue,
    legend = false,
    camera = camera = (-15, 25),
)
# savefig(networkFEFig, "networkPreRem.svg")

# Add more to network
network = DislocationNetwork!(network, prismOct);

## Calculations
# Peach-Köhler force. We don't have a function to apply loading yet so we just generate a random applied stress.
forceDispC.uHat .= rand(-50:50, length(forceDispC.uHat)) .* rand(length(forceDispC.uHat))

fPK = calcPKForce(meshC, forceDispC, network)
fPK69 = calcPKForce(meshC, forceDispC, network, 69)
fPK420 = calcPKForce(meshC, forceDispC, network, [4, 20])
fPK[:, 69] ≈ fPK69
fPK[:, [4, 20]] ≈ fPK420

# Self forces on nodes 1 and 2 of the segments.
sf = calcSelfForce(dlnParams, matParams, network)
sf42 = calcSelfForce(dlnParams, matParams, network, 42)
sf619 = calcSelfForce(dlnParams, matParams, network, [6, 1, 9])
sf[1][:, 42] ≈ sf42[1] && sf[2][:, 42] ≈ sf42[2]
sf[1][:, [6, 1, 9]] ≈ sf619[1] && sf[2][:, [6, 1, 9]] ≈ sf619[2]

# Remote forces.
rf = calcSegSegForce(dlnParams, matParams, network)
rf13 = calcSegSegForce(dlnParams, matParams, network, 13)
rf666 = calcSegSegForce(dlnParams, matParams, network, [66, 6])
rf[:, :, 13] ≈ rf13
rf[:, :, [66, 6]] ≈ rf666

# Total forces.
f = calcSegForce(dlnParams, matParams, meshC, forceDispC, network)
f57 = calcSegForce(dlnParams, matParams, meshC, forceDispC, network, 57)
f316 = calcSegForce(dlnParams, matParams, meshC, forceDispC, network, [31, 6])
f[:, :, 57] ≈ f57
f[:, :, [31, 6]] ≈ f316

calcSegForce!(dlnParams, matParams, meshC, forceDispC, network)
f ≈ network.segForce[:, :, 1:network.numSeg[1]]

# Mobility
# ! This mobility law is outdated.
nodeForce, nodeVel = dlnMobility(dlnParams, matParams, network)
dlnMobility!(dlnParams, matParams, network)
nodeVel ≈ network.nodeVel[:, 1:network.numNode[1]]
nodeForce ≈ network.nodeForce[:, 1:network.numNode[1]]

# Remeshing
# Uncomment following 3 lines to add more loops to network.
# network = DislocationNetwork!(network, (prismOct, shearPent))
# calcSegForce!(dlnParams, matParams, meshC, forceDispC, network)
# dlnMobility!(dlnParams, matParams, network)

# Remesh the surface.
network = remeshSurfaceNetwork!(meshC, boundaryC, network)
# Coarsen the internal network.
network = coarsenNetwork!(dlnParams, matParams, meshC, forceDispC, network)
# Refine the internal network.
network = refineNetwork!(dlnParams, matParams, meshC, forceDispC, network)
remFEFig = plotNodes(
    meshC,
    network,
    m = 2,
    l = 3,
    linecolor = :blue,
    markershape = :circle,
    markercolor = :blue,
    legend = false,
    camera = (-15, 25),
)
idx = findall(x -> x ∈ (3, 4), network.label)
scatter!(
    network.coord[1, idx],
    network.coord[2, idx],
    network.coord[3, idx],
    m = 2,
    markercolor = :red,
)
# savefig(remFEFig, "networkPostRem.svg")

#= 
    collision == true if a collision was found false otherwise.
    collisionType == :null if no collision was found, 
        :hinge if two segments sharing a node collide,
        :twoLine if two unconnected segments collide.
    n1s1 == node 1 segment 1, n2s1 == node 2 segment 1
    n1s2 == node 1 segment 2, n2s2 == node 2 segment 2
    s1 == segment 1, s2 == segment 2
    L1, L2 == normalised position along the segments where collision will occur =#
collision, collisionType, n1s1, n2s1, n1s2, n2s2, s1, s2, L1, L2 =
    detectCollision(dlnParams, network)

# Dislocation displacements on the displacement nodes.
uTilde = calc_uTilde(meshC, boundaryC, matParams, network)
calc_uTilde!(forceDispC, meshC, boundaryC, matParams, network)
uTilde ≈ forceDispC.uTilde[boundaryC.uDofsDln]

reshape(uTilde, :, 3)' == forceDispC.uTilde[reshape(boundaryC.uDofsDln, :, 3)']
meshC.coord[:, boundaryC.uGammaDln] == meshC.coord[reshape(boundaryC.uDofsDln, :, 3)']

coordDisp = copy(meshC.coord[:, boundaryC.uGammaDln])
coordDisp = coordDisp + forceDispC.uTilde[reshape(boundaryC.uDofsDln, :, 3)']

deformedMesh = deepcopy(meshC)
deformedMesh.coord[reshape(boundaryC.uDofsDln, :, 3)'] +=
    20 * forceDispC.uTilde[reshape(boundaryC.uDofsDln, :, 3)']
plotBoundaries(boundaryC, deformedMesh; camera = (-20, 20))
# savefig("dispBound.svg")

scatter!(coordDisp[1, :], coordDisp[2, :], coordDisp[3, :])

meshC.coord[reshape(boundaryC.uDofsDln, :, 3)'] - coordDisp
meshC.connectivity

boundaryC.uGammaDln

# Numeric tractions q = 1
N, A = getTGammaDlnNormsArea(boundaryC, meshC)
points = meshC.coord[:, boundaryC.tGammaDln]
σTilde = calc_σTilde(points, dlnParams, matParams, network)
T = calcNumericTractions(σTilde, N, A)

# Gauss quadrature

# Field point stresses
numNode = 102;
numSeg = numNode - 1;
len = numNode + 1;
xrange = range(-100, 100, length = len);
yrange = range(-100, 100, length = len);

X = ones(length(yrange)) .* xrange';
Y = ones(length(xrange))' .* yrange;
Z = zeros(length(xrange))' .* zeros(len);
points = [X[:]'; Y[:]'; Z[:]'];
X, Y, Z = nothing, nothing, nothing;

a = 5.0
# We need some dummy parameters.
matParams = MaterialParameters(; crystalStruct = BCC(), μ = 1.0, μMag = 1.0, ν = 0.28);
dlnParams = DislocationParameters(;
    coreRad = a,
    coreRadMag = 1.0,
    minSegLen = a + 2,
    maxSegLen = a + 3,
    minArea = a + 1,
    maxArea = a + 2,
    mobility = mobBCC(),
);

# Screw dislocation.

l = Float64[0; 0; 1]
b = Float64[0; 0; 1]
n = b × l

links = zeros(Int, 2, numSeg)
slipPlane = zeros(3, numSeg)
bVec = zeros(3, numSeg)
coord = [zeros(len)'; zeros(len)'; xrange']
label = zeros(nodeTypeDln, len)
nodeVel = similar(coord)
nodeForce = similar(coord)
for i in 1:numSeg
    links[:, i] .= (i, i + 1)
    slipPlane[:, i] = n
    bVec[:, i] = b
end

network = DislocationNetwork(;
    links = links,
    slipPlane = slipPlane,
    bVec = bVec,
    coord = coord,
    label = label,
    nodeVel = nodeVel,
    nodeForce = nodeForce,
)
makeConnect!(network)
getSegmentIdx!(network)

# Out of place calculation.
stress = reshape(calc_σTilde(points, dlnParams, matParams, network), 6, len, :)

σ = zeros(6, size(points, 2))
# In-place calculation.
calc_σTilde!(σ, points, dlnParams, matParams, network)
σ = reshape(σ, 6, len, :)

println("σ ≈ stress : ", σ ≈ stress)

figXY = contourf(xrange, yrange, σ[5, :, :], levels = 30)
figXZ = contourf(xrange, yrange, σ[6, :, :], levels = 30)

# Edge dislocation.

l = Float64[0; 0; 1]
b = Float64[1; 0; 0]
n = b × l

links = zeros(Int, 2, numSeg)
slipPlane = zeros(3, numSeg)
bVec = zeros(3, numSeg)
coord = [zeros(len)'; zeros(len)'; xrange']
label = zeros(nodeTypeDln, len)
nodeVel = similar(coord)
nodeForce = similar(coord)
for i in 1:numSeg
    links[:, i] .= (i, i + 1)
    slipPlane[:, i] = n
    bVec[:, i] = b
end

network = DislocationNetwork(;
    links = links,
    slipPlane = slipPlane,
    bVec = bVec,
    coord = coord,
    label = label,
    nodeVel = nodeVel,
    nodeForce = nodeForce,
)
makeConnect!(network)
getSegmentIdx!(network)

# Out of place calculation.
stress = reshape(calc_σTilde(points, dlnParams, matParams, network), 6, len, :)

σ = zeros(6, size(points, 2))
# In-place calculation.
calc_σTilde!(σ, points, dlnParams, matParams, network)
σ = reshape(σ, 6, len, :)

println("σ ≈ stress : ", σ ≈ stress)

figXX = contourf(xrange, yrange, σ[1, :, :], levels = 30)
figYY = contourf(xrange, yrange, σ[2, :, :], levels = 30)
figZZ = contourf(xrange, yrange, σ[3, :, :], levels = 30)
figXY = contourf(xrange, yrange, σ[4, :, :], levels = 30)
