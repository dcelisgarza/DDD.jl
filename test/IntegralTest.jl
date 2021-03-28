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
        0.7148901889745141,
        0.10526964060096944,
        0.7038040749450563,
        0.9942853358602315,
        0.5209345527205054,
        0.9870544902707334,
        0.30684645103929964,
        0.7051655060748223,
        0.07334165035468887,
        0.27273647674203194,
        0.1160202967848678,
        0.6301327760714026,
        0.3819425782046564,
        0.45831256554373434,
        0.12112688382206094,
        0.5381334236550053,
        0.17928196607698998,
        0.022709432453792866,
        0.44323022551317925,
        0.6438541879202242,
        0.47413802450796005,
        0.003779619119881117,
        0.959979232067655,
        0.48097298829509527,
    ]
    network2.nodeVel[2, 1:numNode] .= [
        0.5651131933578388,
        0.10791774525913067,
        0.63487059263634,
        0.5163123641541374,
        0.22151119293431365,
        0.27669947756931856,
        0.19942512197253914,
        0.4059776094900336,
        0.08264799460993633,
        0.6767456320347092,
        0.9736199703806871,
        0.8091298360198254,
        0.5635227594205467,
        0.2431143863899523,
        0.785351623353673,
        0.4611200680355141,
        0.7200163333834613,
        0.48110624588973216,
        0.20950933276652828,
        0.05190375487521082,
        0.03686221869674755,
        0.3730038293946101,
        0.30874254232031184,
        0.4146316224254616,
    ]
    network2.nodeVel[3, 1:numNode] .= [
        0.8558270563635098,
        0.6437634767707383,
        0.6190435004230979,
        0.08784253136806419,
        0.8632503495066168,
        0.11429209621406189,
        0.12221799840096659,
        0.33849307986842714,
        0.6403947185869319,
        0.3637313311773851,
        0.9546414393908917,
        0.5662235687559833,
        0.8775219339303044,
        0.21926659450593444,
        0.5741266736432413,
        0.4348897841453008,
        0.44219356603365756,
        0.23163772406660454,
        0.25640994295825226,
        0.03538648273719436,
        0.7287992956844731,
        0.11886293329970665,
        0.9888538124016337,
        0.2812829069918916,
    ]

    remeshSurfaceNetwork!(regularCuboidMesh, cantileverBC, network2)
    getSegmentIdx!(network2)

    network2.numNode[1] == numNode + 2
    label = network2.label
    ext = findall(x -> x == 5, label)
    surf = findall(x -> x == 3, label)
    @test length(ext) == 5
    @test length(surf) == 2
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
