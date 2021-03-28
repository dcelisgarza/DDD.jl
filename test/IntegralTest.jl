using DDD, Test, SparseArrays, LinearAlgebra, StaticArrays
cd(@__DIR__)

@testset "Integral tests" begin
    dlnParams = DislocationParameters(; mobility = mobBCC())
    matParams = MaterialParameters(; crystalStruct = BCC())
    femParams = FEMParameters(;
        type = DispatchRegularCuboidMesh(),
        order = LinearElement(),
        model = CantileverLoad(),
        dx = 1013.0,
        dy = 1987.0,
        dz = 2999.0,
        mx = 19,
        my = 23,
        mz = 29,
    )
    intParams = IntegrationParameters(; method = AdaptiveEulerTrapezoid())
    slipSystem = SlipSystem(;
        crystalStruct = BCC(),
        slipPlane = Float64[-1; 1; 0],
        bVec = Float64[1; 1; 1],
    )
    dx, dy, dz = femParams.dx, femParams.dy, femParams.dz
    segLen = (dx + dy + dz) / 30
    prismLoop = DislocationLoop(;
        loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
        numSides = 8,   # 8-sided loop.
        nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 1,   # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystemIdx = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        slipSystem = slipSystem,  # Slip plane of the segments.
        label = SVector{8, nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3, 2, Float64}(dx / 2, dy / 2, dz / 2, dx / 2, dy / 2, dz / 2),  # Distribution range
        dist = Zeros(),  # Loop distribution.
    )
    shearLoop = DislocationLoop(;
        loopType = loopShear(),    # Prismatic loop, all segments are edge segments.
        numSides = 8,   # 8-sided loop.
        nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 1,   # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystemIdx = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        slipSystem = slipSystem,  # Slip plane of the segments.
        label = SVector{8, nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3, 2, Float64}(dx / 2, dy / 2, dz / 2, dx / 2, dy / 2, dz / 2),  # Distribution range
        dist = Zeros(),  # Loop distribution.
    )
    network = DislocationNetwork([prismLoop, shearLoop])
    regularCuboidMesh = buildMesh(matParams, femParams)
    cantileverBC, forceDisplacement = Boundaries(femParams, regularCuboidMesh)

    coordFE = regularCuboidMesh.coord
    uGamma = cantileverBC.uGamma.node
    uDofsDln = cantileverBC.uDofsDln
    numSeg = network.numSeg[1]

    segForce =
        calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    @test isapprox(segForce, network.segForce[:, :, 1:numSeg])

    idx = 1:2:16
    segForceIdx = calcSegForce(
        dlnParams,
        matParams,
        regularCuboidMesh,
        forceDisplacement,
        network,
        idx,
    )
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(segForceIdx, network.segForce[:, :, idx])

    σTilde = calc_σTilde(coordFE[:, uGamma], dlnParams, matParams, network)
    σTilde2 = zeros(size(σTilde))
    calc_σTilde!(σTilde2, coordFE[:, uGamma], dlnParams, matParams, network)
    @test isapprox(σTilde, σTilde2)

    uTilde = calc_uTilde(regularCuboidMesh, cantileverBC, matParams, network)
    calc_uTilde!(forceDisplacement, regularCuboidMesh, cantileverBC, matParams, network)
    @test isapprox(uTilde, forceDisplacement.uTilde[uDofsDln])

    idx = [
        5771
        6852
        7419
        7436
        6903
        5822
        5255
        5238
        7442
        6852
        6298
        5213
        5232
        5822
        6376
        7461
    ]
    uHatIdx = unique(regularCuboidMesh.connectivity[:, idx]) * 3
    uHatDofs = [uHatIdx .- 2; uHatIdx .- 1; uHatIdx]
    randUHat = sprand(length(uHatDofs), 0.5)
    randIdx = findall(!iszero, randUHat)
    randUHat[randIdx] .= 1:(length(randIdx) * 0.01)
    forceDisplacement.uHat[uHatDofs] = randUHat

    segForce =
        calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    @test isapprox(segForce, network.segForce[:, :, 1:numSeg])

    idx = 1:2:16
    segForceIdx = calcSegForce(
        dlnParams,
        matParams,
        regularCuboidMesh,
        forceDisplacement,
        network,
        idx,
    )
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(segForceIdx, network.segForce[:, :, idx])

    network2 = deepcopy(network)
    remeshSurfaceNetwork!(regularCuboidMesh, cantileverBC, network2)
    getSegmentIdx!(network2)
    @test compStruct(network2, network)

    prismLoopIntersect = DislocationLoop(;
        loopType = loopShear(),    # Prismatic loop, all segments are edge segments.
        numSides = 8,   # 8-sided loop.
        nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 1,   # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystemIdx = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        slipSystem = slipSystem,  # Slip plane of the segments.
        label = SVector{8, nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3, 2, Float64}(
            segLen / 2,
            segLen / 2,
            segLen / 2,
            segLen / 2,
            segLen / 2,
            segLen / 2,
        ),  # Distribution range
        dist = Zeros(),  # Loop distribution.
    )
    DislocationNetwork!(network2, prismLoopIntersect)
    numNode = network2.numNode[1]

    network2.nodeVel[1, 1:numNode] .= [
        0.848409575718539,
        0.6652242160082342,
        0.25113206533518206,
        0.17655640207336298,
        0.33256805603630246,
        0.9969648349174243,
        0.12851456965799923,
        0.05591976789414588,
        0.29757306575603804,
        0.2723507312786162,
        0.8874071844807814,
        0.4429058105879311,
        0.8868742664817155,
        0.4598309743323694,
        0.03931969223629794,
        0.12373984299000429,
        0.03861761480269377,
        0.7603101804979915,
        0.4543587651493177,
        0.27460243590637345,
        0.3,
        0.5434785594380296,
        0.5395078480115256,
        0.8654239757909008,
    ]
    network2.nodeVel[2, 1:numNode] .= [
        0.8636080299591222,
        0.37959594965548926,
        0.9920015465135401,
        0.702809412405236,
        0.14450700182749543,
        0.4816990411158426,
        0.2905765122000159,
        0.319012344824527,
        0.05952526906892097,
        0.4392216093355392,
        0.8752415268798401,
        0.3592278150943291,
        0.8529140139462188,
        0.0993328662955082,
        0.14584336379153573,
        0.43739168117035354,
        0.3338215613113067,
        0.08518695747920568,
        0.7462534959446294,
        0.5272909447878065,
        0.5,
        0.7878515928217127,
        0.5106406121544798,
        0.22360673480968973,
    ]
    network2.nodeVel[3, 1:numNode] .= [
        0.38025510658061124,
        0.028912486176383645,
        0.24395899994533154,
        0.8387909328083565,
        0.041107268905112626,
        0.6138145136425359,
        0.5325226796794968,
        0.21895716122452558,
        0.27194053775758875,
        0.6585257914167739,
        0.6004564313623306,
        0.10944691655246075,
        0.7382692437252569,
        0.9260924108803184,
        0.46106334186849685,
        0.5657795344686234,
        0.09638787333766241,
        0.2854907820375585,
        0.04887167777640866,
        0.850652937369534,
        0.2,
        0.7285145955076195,
        0.39460521446188235,
        0.0044684000305827976,
    ]

    remeshSurfaceNetwork!(regularCuboidMesh, cantileverBC, network2)
    getSegmentIdx!(network2)

    network2.numNode[1] == numNode + 2
    label = network2.label
    ext = findall(x -> x == 5, label)
    surf = findall(x -> x == 3, label)
    @test length(ext) == 3
    @test length(surf) == 1
    @test compStruct(network2, network) == false

    numSeg = network2.numSeg[1]
    segForce =
        calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network2)
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network2)
    @test isapprox(segForce, network2.segForce[:, :, 1:numSeg])

    idx = 1:2:16
    segForceIdx = calcSegForce(
        dlnParams,
        matParams,
        regularCuboidMesh,
        forceDisplacement,
        network2,
        idx,
    )
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network2, idx)
    @test isapprox(segForceIdx, network2.segForce[:, :, idx])

    σTilde = calc_σTilde(coordFE[:, uGamma], dlnParams, matParams, network2)
    σTilde2 = zeros(size(σTilde))
    calc_σTilde!(σTilde2, coordFE[:, uGamma], dlnParams, matParams, network2)
    @test isapprox(σTilde, σTilde2)

    uTilde = calc_uTilde(regularCuboidMesh, cantileverBC, matParams, network2)
    calc_uTilde!(forceDisplacement, regularCuboidMesh, cantileverBC, matParams, network2)
    @test isapprox(uTilde, forceDisplacement.uTilde[uDofsDln])

    idx = [
        5771
        6852
        7419
        7436
        6903
        5822
        5255
        5238
        7442
        6852
        6298
        5213
        5232
        5822
        6376
        7461
    ]
    uHatIdx = unique(regularCuboidMesh.connectivity[:, idx]) * 3
    uHatDofs = [uHatIdx .- 2; uHatIdx .- 1; uHatIdx]
    randUHat = sprand(length(uHatDofs), 0.5)
    randIdx = findall(!iszero, randUHat)
    randUHat[randIdx] .= 1:(length(randIdx) * 0.01)
    forceDisplacement.uHat[uHatDofs] = randUHat

    segForce =
        calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network2)
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network2)
    @test isapprox(segForce, network2.segForce[:, :, 1:numSeg])

    idx = 1:2:16
    segForceIdx = calcSegForce(
        dlnParams,
        matParams,
        regularCuboidMesh,
        forceDisplacement,
        network2,
        idx,
    )
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network2, idx)
    @test isapprox(segForceIdx, network2.segForce[:, :, idx])

    network3 = deepcopy(network2)
    remeshSurfaceNetwork!(regularCuboidMesh, cantileverBC, network3)
    getSegmentIdx!(network3)

    compStruct(network3, network2)
    coarsenVirtualNetwork!(dlnParams, network3)
    @test compStruct(network3, network2)
end
