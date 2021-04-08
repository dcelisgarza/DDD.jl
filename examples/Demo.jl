using DDD, Plots, LinearAlgebra, Statistics, SparseArrays, StaticArrays, BenchmarkTools
plotlyjs()

dlnParams = DislocationParameters(;
    mobility = mobBCC(),
    dragCoeffs = (edge = 1, screw = 1e-1, climb = 1e9, line = 1e-5),
    coreRad = 0.015 * sqrt(2),
    minSegLen = 0.15 * sqrt(2),
    maxSegLen = 1.5 * sqrt(2),
    coreEnergy = 1 / (4 * π) * log(0.015 * sqrt(2) / 3.5e-5),
    coreRadMag = 3.5e-4,
)
# (:mobility, :dragCoeffs, :coreRad, :coreRadSq, :coreRadMag, :coreEnergy, :minSegLen, :maxSegLen, :twoMinSegLen, :minSegLenSq, :minArea, :maxArea, :minAreaSq, :maxAreaSq, :collisionDist, :collisionDistSq, :slipStepCritLen, :slipStepCritArea, :remesh, :collision, :separation, :virtualRemesh, :parCPU, :parGPU)

matParams = MaterialParameters(; crystalStruct = BCC(), μ = 1, μMag = 80e3, ν = 0.25)
# (:crystalStruct, :μ, :μMag, :ν, :omνInv, :opνInv, :νomνInv, :νopνInv, :μ4π, :μ8π, :μ4πν, :νμ4πν, :omνInv8π, :om2νomνInv8π, :σPN)

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
)
# (:type, :order, :model, :dx, :dy, :dz, :mx, :my, :mz)

slipSystems = SlipSystem(;
    crystalStruct = BCC(),
    slipPlane = Float64[-1 1; 1 -1; 0 0],
    bVec = Float64[1 1; 1 1; 1 -1],
)
# (:crystalStruct, :slipPlane, :bVec)

intParams = IntegrationParameters(;
    method = AdaptiveEulerTrapezoid(),
    abstol = dlnParams.collisionDist / 2,
    reltol = dlnParams.collisionDist / 2,
)
# (:method, :tmin, :tmax, :dtmin, :dtmax, :abstol, :reltol, :maxchange, :exponent, :maxiter)

intTime = IntegrationTime(; dt = 0.0, time = 0.0)
# (:crystalStruct, :slipPlane, :bVec)

meshC = buildMesh(matParams, femParamsC)
# (:order, :dx, :dy, :dz, :mx, :my, :mz, :w, :h, :d, :scale, :numElem, :numNode, :C, :vertices, :faces, :faceNorm, :faceMidPt, :cornerNode, :edgeNode, :faceNode, :surfNode, :surfNodeArea, :surfNodeNorm, :surfElemNode, :coord, :connectivity, :K)
figFE = plotFEDomain(meshC)

boundaryC, forceDispC = Boundaries(femParamsC, meshC)
# (:noExit, :uGammaDln, :tGammaDln, :uDofsDln, :tDofsDln, :uGamma, :tGamma, :mGamma, :uDofs, :tDofs, :mDofs, :tK)
# (:uTilde, :uHat, :u, :fTilde, :fHat, :f)

figCBC = plotBoundaries(boundaryC, meshC)

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
# (:numNode, :numSeg, :maxConnect, :label, :links, :connectivity, :linksConnect, :slipPlane, :segIdx, :bVec, :coord, :nodeVel, :nodeForce, :segForce)

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
    markerstyle = :square,
    markercolor = :red,
    legend = false,
)

networkFE = plotNodes(
    meshC,
    network,
    m = 2,
    l = 3,
    linecolor = :blue,
    markershape = :circle,
    markercolor = :blue,
    legend = false,
)

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
f316 = calcSegForce(
    dlnParams,
    matParams,
    meshC,
    forceDispC,
    network,
    [31, 6],
)
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

# Field point stresses
points = meshC.coord[:, boundaryC.tGammaDln]
σTilde = calc_σTilde(points, dlnParams, matParams, network)
σTilde2 = zeros(6, length(boundaryC.tGammaDln))
calc_σTilde!(σTilde2 , points, dlnParams, matParams, network)
σTilde ≈ σTilde2

# Dislocation displacements on the displacement nodes.
uTilde = calc_uTilde(meshC, boundaryC, matParams, network)
calc_uTilde!(forceDispC, meshC, boundaryC, matParams, network)
uTilde ≈ forceDispC.uTilde[boundaryC.uDofsDln]

# Numeric tractions q = 1

# Gauss quadrature




femParamsPillar = FEMParameters(;
    type = DispatchRegularCuboidMesh(),
    order = LinearElement(),
    model = PillarLoad(),
    dx = 17.0,
    dy = 13.0,
    dz = 23.0,
    mx = 5,
    my = 3,
    mz = 7,
)
meshC = buildMesh(matParams, femParamsPillar)
pillarLoadBC, forceDispC = Boundaries(femParamsPillar, meshC)
figPBC = plotBoundaries(pillarLoadBC, meshC)