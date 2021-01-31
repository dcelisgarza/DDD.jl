using DDD, Test, SparseArrays, LinearAlgebra, StaticArrays

@testset "Integral tests" begin
    dlnParams = DislocationParameters(; mobility = mobBCC())
    matParams = MaterialParameters(; crystalStruct = BCC())
    femParams = FEMParameters(; 
                            type = DispatchRegularCuboidMesh(), 
                            order = LinearElement(), 
                            model = CantileverLoad(), 
                            dx = 1013.0, dy = 1987.0, dz = 2999.0,
                            mx = 19, my = 23, mz = 29
                        )
    intParams = IntegrationParameters(; method = AdaptiveEulerTrapezoid())
    slipSystem = SlipSystem(; crystalStruct = BCC(), slipPlane = Float64[-1;1;0], bVec = Float64[1;1;1])
    dx, dy, dz = femParams.dx, femParams.dy, femParams.dz
    segLen = (dx + dy + dz) / 30
    prismLoop = DislocationLoop(;
        loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
        numSides = 8,   # 8-sided loop.
        nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 1,   # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystem = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        _slipPlane = slipSystem.slipPlane[:, 1],  # Slip plane of the segments.
        _bVec = slipSystem.bVec[:, 1],            # Burgers vector of the segments.
        label = SVector{8,nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3,2,Float64}(dx / 2, dy / 2, dz / 2, dx / 2, dy / 2, dz / 2),  # Distribution range
        dist = Zeros(),  # Loop distribution.
    )
    shearLoop = DislocationLoop(;
        loopType = loopShear(),    # Prismatic loop, all segments are edge segments.
        numSides = 8,   # 8-sided loop.
        nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 1,   # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystem = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        _slipPlane = slipSystem.slipPlane[:, 1],  # Slip plane of the segments.
        _bVec = slipSystem.bVec[:, 1],            # Burgers vector of the segments.
        label = SVector{8,nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3,2,Float64}(dx / 2, dy / 2, dz / 2, dx / 2, dy / 2, dz / 2),  # Distribution range
        dist = Zeros(),  # Loop distribution.
    )
    network = DislocationNetwork((prismLoop, shearLoop))
    regularCuboidMesh = buildMesh(matParams, femParams)
    cantileverBC, forceDisplacement = Boundaries(femParams, regularCuboidMesh)

    coordFE = regularCuboidMesh.coord
    uGamma = cantileverBC.uGamma[:node]
    uDofsDln = cantileverBC.uDofsDln
    numSeg = network.numSeg[1]

    segForce = calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    @test isapprox(segForce, network.segForce[:,:,1:numSeg])

    idx = 1:2:16
    segForceIdx = calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network, idx)
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(segForceIdx, network.segForce[:,:,idx])

    σTilde = calc_σTilde(coordFE[:, uGamma], dlnParams, matParams, network)
    σTilde2 = zeros(size(σTilde))
    calc_σTilde!(σTilde2, coordFE[:, uGamma], dlnParams, matParams, network)
    @test isapprox(σTilde, σTilde2)

    uTilde = calc_uTilde(forceDisplacement, regularCuboidMesh, cantileverBC, matParams, network)
    calc_uTilde!(forceDisplacement, regularCuboidMesh, cantileverBC, matParams, network)
    @test isapprox(uTilde, forceDisplacement.uTilde[uDofsDln])

    idx = [5771;6852;7419;7436;6903;5822;5255;5238;7442;6852;6298;5213;5232;5822;6376;7461]
    uHatIdx = unique(regularCuboidMesh.connectivity[:, idx]) * 3
    uHatDofs = [uHatIdx .- 2; uHatIdx .- 1; uHatIdx]
    randUHat = sprand(length(uHatDofs), 0.5) * 10
    forceDisplacement.uHat[uHatDofs] = randUHat

    segForce = calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    @test isapprox(segForce, network.segForce[:,:,1:numSeg])

    idx = 1:2:16
    segForceIdx = calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network, idx)
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(segForceIdx, network.segForce[:,:,idx])

    network2 = deepcopy(network)
    remeshSurfaceNetwork!(regularCuboidMesh, network2)
    @test compStruct(network2, network)

    prismLoopIntersect = DislocationLoop(;
        loopType = loopShear(),    # Prismatic loop, all segments are edge segments.
        numSides = 8,   # 8-sided loop.
        nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 1,   # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystem = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        _slipPlane = slipSystem.slipPlane[:, 1],  # Slip plane of the segments.
        _bVec = slipSystem.bVec[:, 1],            # Burgers vector of the segments.
        label = SVector{8,nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3,2,Float64}(segLen / 2, segLen / 2, segLen / 2, segLen / 2, segLen / 2, segLen / 2),  # Distribution range
        dist = Zeros(),  # Loop distribution.
    )
    DislocationNetwork!(network2, prismLoopIntersect)
    numNode = network2.numNode[1]

    network2.nodeVel[:, 1:numNode] .= rand(3, numNode)

    remeshSurfaceNetwork!(regularCuboidMesh,  network2)

    network2.numNode[1] == numNode + 2
    label = network2.label
    ext = findall(x -> x == 5, label)
    surf = findall(x -> x == 3, label)
    @test length(ext) == 5
    @test length(surf) == 2
    @test compStruct(network2, network) == false
end