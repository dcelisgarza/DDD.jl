using DDD, Plots, LinearAlgebra, Statistics, SparseArrays, StaticArrays, BenchmarkTools
plotlyjs()

dlnParams = DislocationParameters(;
    mobility = mobBCC(),
    dragCoeffs = (edge = 1, screw = 1e-1, climb = 1e9),
    coreRad = 0.015 * sqrt(2),
    minSegLen = 0.15 * sqrt(2),
    maxSegLen = 1.5 * sqrt(2),
    coreEnergy = 1 / (4 * π) * log(0.015 * sqrt(2) / 3.5e-5),
    coreRadMag = 3.5e-4,
)
# (:mobility, :dragCoeffs, :coreRad, :coreRadSq, :coreRadMag, :coreEnergy, :minSegLen, :maxSegLen, :twoMinSegLen, :minSegLenSq, :minArea, :maxArea, :minAreaSq, :maxAreaSq, :collisionDist, :collisionDistSq, :slipStepCritLen, :slipStepCritArea, :remesh, :collision, :separation, :virtualRemesh, :parCPU, :parGPU)

matParams = MaterialParameters(;
    crystalStruct = BCC(),
    μ = 1,
    μMag = 80e3,
    ν = 0.25,
)
# (:crystalStruct, :μ, :μMag, :ν, :omνInv, :opνInv, :νomνInv, :νopνInv, :μ4π, :μ8π, :μ4πν, :νμ4πν, :omνInv8π, :om2νomνInv8π, :σPN)

femParams = FEMParameters(;
    type = DispatchRegularCuboidMesh(),
    order = LinearElement(),
    model = CantileverLoad(),
    dx = 43.0,
    dy = 37.0,
    dz = 31.0,
    mx = 27,
    my = 23,
    mz = 19
)
# (:type, :order, :model, :dx, :dy, :dz, :mx, :my, :mz)

slipSystems = SlipSystem(;
    crystalStruct = BCC(),
    slipPlane = Float64[-1 1 ;1 -1;0 0],
    bVec = Float64[1 1;1 1;1 -1]
)
# (:crystalStruct, :slipPlane, :bVec)

intParams = IntegrationParameters(;
    method = AdaptiveEulerTrapezoid(),
    abstol = dlnParams.collisionDist / 2,
    reltol = dlnParams.collisionDist / 2,
)
# (:method, :tmin, :tmax, :dtmin, :dtmax, :abstol, :reltol, :maxchange, :exponent, :maxiter)

intTime = IntegrationTime(;
    dt = 0.0,
    time = 0.0
)
# (:crystalStruct, :slipPlane, :bVec)

regularCuboidMesh = buildMesh(matParams, femParams)
# (:order, :dx, :dy, :dz, :mx, :my, :mz, :w, :h, :d, :scale, :numElem, :numNode, :C, :vertices, :faces, :faceNorm, :faceMidPt, :cornerNode, :edgeNode, :faceNode, :surfNode, :surfNodeArea, :surfNodeNorm, :surfElemNode, :coord, :connectivity, :K)
figFE = plotFEDomain(regularCuboidMesh)

cantileverBC, forceDisplacement = Boundaries(femParams, regularCuboidMesh)
# (:noExit, :uGammaDln, :tGammaDln, :uDofsDln, :tDofsDln, :uGamma, :tGamma, :mGamma, :uDofs, :tDofs, :mDofs, :tK)
# (:uTilde, :uHat, :u, :fTilde, :fHat, :f)

feCoord = regularCuboidMesh.coord
uGamma = cantileverBC.uGamma.node
mGamma = cantileverBC.mGamma.node
tGamma = cantileverBC.tGamma.node

function plotBoundaries(::FEMParameters{T1,T2,T3} where {T1,T2,T3 <: CantileverLoad}, boundaries::Boundaries, mesh::RegularCuboidMesh)
    uGamma = boundaries.uGamma.node
    mGamma = boundaries.mGamma.node
    tGamma = boundaries.tGamma.node
    feCoord = mesh.coord

    fig = scatter(feCoord[1, uGamma], feCoord[2, uGamma], feCoord[3, uGamma]; markershape = :square, markercolor = :red, markersize = 3, label = "Fixed End")
    scatter!(fig, feCoord[1, mGamma], feCoord[2, mGamma], feCoord[3, mGamma]; markershape = :diamond, markercolor = :cyan, markersize = 3, label = "Loaded Edge")
    scatter!(fig, feCoord[1, tGamma], feCoord[2, tGamma], feCoord[3, tGamma]; markershape = :circle, markercolor = :black, markersize = 3, label = "Traction Surface")
    return fig
end
plotBoundaries(femParams, cantileverBC, regularCuboidMesh)

bcFig = scatter(feCoord[1, uGamma], feCoord[2, uGamma], feCoord[3, uGamma]; markershape = :square, markercolor = :red, markersize = 3, label = "Fixed End")
scatter!(bcFig, feCoord[1, mGamma], feCoord[2, mGamma], feCoord[3, mGamma]; markershape = :diamond, markercolor = :cyan, markersize = 3, label = "Loaded Edge")
scatter!(bcFig, feCoord[1, tGamma], feCoord[2, tGamma], feCoord[3, tGamma]; markershape = :circle, markercolor = :black, markersize = 3, label = "Traction Surface")

dx, dy, dz = femParams.dx, femParams.dy, femParams.dz
segLen = (dlnParams.minSegLen + dlnParams.maxSegLen) / 2

prismOct = DislocationLoop(;
    loopType = loopPrism(),
    numSides = 8,
    nodeSide = 1,
    numLoops = 20,
    segLen = segLen * ones(8),
    slipSystemIdx = 1,
    slipSystem = slipSystems,
    label = nodeTypeDln.(ones(Int, 8)),
    buffer = 0,
    range = [0 dx; 0 dy; 0 dz],
    dist = Rand()
)

shearPent = DislocationLoop(;
    loopType = loopShear(),
    numSides = 5,
    nodeSide = 2,
    numLoops = 30,
    segLen = segLen * ones(10),
    slipSystemIdx = 2,
    slipSystem = slipSystems,
    label = nodeTypeDln.(ones(Int, 10)),
    buffer = 0,
    range = [0 dx; 0 dy; 0 dz],
    dist = Rand()
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
    regularCuboidMesh,
    network,
    m = 2,
    l = 3,
    linecolor = :blue,
    markershape = :circle,
    markercolor = :blue,
    legend = false,
)
