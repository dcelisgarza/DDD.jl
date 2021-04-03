using DDD
using Test, SparseArrays, StaticArrays
cd(@__DIR__)

@testset "Forces" begin
    fileDislocationParameters = "./testData/sampleDislocationParameters.json"
    fileMaterialParameters = "./testData/sampleMaterialParameters.json"
    fileFEMParameters = "./testData/sampleFEMParameters.json"
    fileIntegrationParameters = "./testData/sampleIntegrationParameters.json"
    fileSlipSystem = "./testData/BCC.json"
    fileDislocationLoop = "./testData/samplePrismShear.json"
    dlnParams, matParams, femParams, intParams, slipSystems, dislocationLoop =
        loadParameters(
            fileDislocationParameters,
            fileMaterialParameters,
            fileFEMParameters,
            fileIntegrationParameters,
            fileSlipSystem,
            fileDislocationLoop,
        )
    regularCuboidMesh = buildMesh(matParams, femParams)
    numFEMNode = regularCuboidMesh.numNode
    dx = regularCuboidMesh.dx
    dy = regularCuboidMesh.dy
    dz = regularCuboidMesh.dz
    segLen = dx / 8

    f = spzeros(3 * numFEMNode)
    f[[112, 118, 133, 141, 213, 244, 262, 272, 317]] = [
        0.43048187784858616,
        0.22724536603830137,
        0.4340867899691503,
        0.6660863546953892,
        0.30358515797696106,
        0.2945958951093859,
        0.7278367502911502,
        0.7095924334694701,
        0.1642050526375538,
    ]
    fHat = spzeros(3 * numFEMNode)
    fHat[[32, 48, 55, 88, 138, 148, 191, 230, 253, 335]] = [
        0.09706224225842108,
        0.07773687633248638,
        0.13682398802299178,
        0.4752286167553166,
        0.7423196193496164,
        0.8286077556473421,
        0.7023632196408749,
        0.9813639162461198,
        0.5296701796678411,
        0.5523797553266823,
    ]
    u = spzeros(3 * numFEMNode)
    u[[30, 127, 195, 221, 316, 325, 338, 348, 370]] = [
        0.8792592573507609,
        0.8430664083925272,
        0.4711050560756602,
        0.4860071865093816,
        0.7905698600135145,
        0.39047211692578077,
        0.6545538020629462,
        0.5446700211111557,
        0.8865721648558644,
    ]
    uHat = spzeros(3 * numFEMNode)
    uHat[[91, 126, 130, 195, 217, 226, 229, 256, 281, 293, 309, 342]] = [
        0.5231621885339968,
        0.5771429489788034,
        0.7151190318538345,
        0.7283662326812077,
        0.6314274719472075,
        0.9814688915693632,
        0.5672795171250207,
        0.002712918060655989,
        0.1788941754890383,
        0.188299784057536,
        0.8489027048214433,
        0.029995302953659708,
    ]
    forceDisplacement =
        ForceDisplacement(nothing, uHat * 1000, u * 1000, nothing, fHat * 1000, f * 1000)

    prismSquare = DislocationLoop(;
        loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
        numSides = 4,   # 5-sided loop.
        nodeSide = 2,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 2,  # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystemIdx = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        slipSystem = slipSystems,  # Slip plane of the segments.
        label = SVector{8, nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3, 2, Float64}(
            0 + segLen,
            0 + segLen,
            0 + segLen,
            dx - segLen,
            dy - segLen,
            dz - segLen,
        ),  # Distribution range
        dist = Rand(),  # Loop distribution.
    )

    shearSquare = DislocationLoop(;
        loopType = loopShear(),    # Prismatic loop, all segments are edge segments.
        numSides = 4,   # 5-sided loop.
        nodeSide = 2,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 2,  # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystemIdx = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        slipSystem = slipSystems,  # Slip plane of the segments.
        label = SVector{8, nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3, 2, Float64}(
            0 + segLen,
            0 + segLen,
            0 + segLen,
            dx - segLen,
            dy - segLen,
            dz - segLen,
        ),  # Distribution range
        dist = Rand(),  # Loop distribution.
    )
    network = DislocationNetwork([prismSquare, shearSquare])
    numSeg = network.numSeg[1]
    network.coord[:, 1:numSeg] = [
        1201.5837679702274 1283.233426063 1424.6547823003095 1566.076138537619 1484.4264804448464 1402.7768223520739 1261.3554661147643 1119.9341098774548 304.6214867750916 386.2711448678642 527.6925011051737 669.1138573424832 587.4641992497106 505.81454115693805 364.39318491962854 222.97182868231897 1270.8257124306385 1352.4753705234111 1467.9454243613363 1583.4154781992615 1501.765820106489 1420.1161620137163 1304.6461081757911 1189.176054337866 462.2031932442291 543.8528513370018 659.3229051749269 774.792959012852 693.1433009200794 611.4936428273068 496.02358898938166 380.5535351514565
        696.4782710631935 614.8286129704209 756.2499692077304 897.6713254450399 979.3209835378125 1060.970641630585 919.5492853932756 778.1279291559661 302.5604363129223 220.91077822014972 362.33213445745923 503.75349069476874 585.4031487875413 667.0528068803139 525.6314506430044 384.210094405695 1055.3880794417687 973.7384213489961 858.2683675110709 742.7983136731457 824.4479717659183 906.0976298586909 1021.5676836966161 1137.0377375345413 1195.4446644351433 1113.7950063423707 998.3249525044455 882.8548986665203 964.5045567592929 1046.1542148520655 1161.6242686899907 1277.094322527916
        1219.6759003965913 1382.9752165821365 1382.9752165821365 1382.9752165821365 1219.6759003965913 1056.376584211046 1056.376584211046 1056.376584211046 1604.4473849682922 1767.7467011538374 1767.7467011538374 1767.7467011538374 1604.4473849682922 1441.148068782747 1441.148068782747 1441.148068782747 951.679209138293 1114.9785253238383 999.5084714859131 884.038417647988 720.7391014624428 557.4397852768976 672.9098391148227 788.3798929527478 1431.0547935274708 1594.354109713016 1478.8840558750908 1363.4140020371656 1200.1146858516204 1036.8153696660752 1152.2854235040004 1267.7554773419254
    ]

    fPKTest = [
        0.0 0.0 -0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        -0.0 0.0 0.0
        0.0298244946801427 -3.935141923601883 -1.982483209141013
        -0.4565798368140385 0.4565798368140385 40.759564106652746
        11.942933454784617 -11.942933454784617 50.64827964686894
        0.0 0.0 -0.0
        0.0 0.0 -0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        -0.0 0.0 0.0
        -0.0 0.0 0.0
        0.0 -0.0 0.0
        0.0 -0.0 0.0
        0.0 0.0 -0.0
        34.47756746101095 25.45149478892432 -4.513036336043308
        -13.183687913414278 -7.188653151905362 -5.995034761508917
        -7.178181870597174 -0.03305076155172504 -7.1451311090454555
        -0.0 0.0 0.0
        -41.951995497071024 -22.03609089968592 9.957952298692552
        25.62507197894962 -0.6852147795546628 26.310286758504304
        60.72042574006913 -1.1681574909891737 61.88858323105836
        108.2078469749699 72.74142517615954 -17.73321089940518
        0.0 0.0 -0.0
        0.0 0.0 -0.0
        0.0 0.0 -0.0
        -0.0 0.0 0.0
        0.015530976306194433 0.027503198516044464 0.005986111104925015
        0.006498597861991796 0.016097702989461442 -0.00959910512746965
        0.0 -0.0 0.0
        0.0 0.0 -0.0
    ]

    fPK = calcPKForce(regularCuboidMesh, forceDisplacement, network)
    @test isapprox(fPK', fPKTest)

    network.segForce .= 0
    calcPKForce!(regularCuboidMesh, forceDisplacement, network)
    @test isapprox(fPK, network.segForce[:, 1, 1:numSeg] * 2)
    @test isapprox(network.segForce[:, 1, 1:numSeg], network.segForce[:, 2, 1:numSeg])

    idx = 5
    fPK = calcPKForce(regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(fPK, fPKTest[idx, :])

    network.segForce .= 0
    calcPKForce!(regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(fPK, network.segForce[:, 1, idx] * 2)
    @test isapprox(network.segForce[:, 1, idx], network.segForce[:, 2, idx])

    idx = [5; 10; 9; 23]
    fPK = calcPKForce(regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(fPK', fPKTest[idx, :])

    network.segForce .= 0
    calcPKForce!(regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(fPK, network.segForce[:, 1, idx] * 2)
    @test isapprox(network.segForce[:, 1, idx], network.segForce[:, 2, idx])

    idx = 15:32
    fPK = calcPKForce(regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(fPK', fPKTest[idx, :])

    network.segForce .= 0
    calcPKForce!(regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(fPK, network.segForce[:, 1, idx] * 2)
    @test isapprox(network.segForce[:, 1, idx], network.segForce[:, 2, idx])

    dlnParamsPar = DislocationParameters(;
        coreRad = dlnParams.coreRad,
        coreRadMag = dlnParams.coreRadMag,
        coreEnergy = dlnParams.coreEnergy,
        minSegLen = dlnParams.minSegLen,
        maxSegLen = dlnParams.maxSegLen,
        minArea = dlnParams.minArea,
        maxArea = dlnParams.maxArea,
        remesh = dlnParams.remesh,
        collision = dlnParams.collision,
        separation = dlnParams.separation,
        virtualRemesh = dlnParams.virtualRemesh,
        parCPU = true,
        parGPU = dlnParams.parGPU,
        dragCoeffs = dlnParams.dragCoeffs,
        mobility = dlnParams.mobility,
    )

    pentagon = DislocationLoop(;
        loopType = loopPrism(),
        numSides = 5,
        nodeSide = 1,
        numLoops = 2,
        segLen = 10 * ones(5),
        slipSystemIdx = 4,
        slipSystem = slipSystems,
        label = nodeTypeDln[1; 2; 1; 2; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork(pentagon; memBuffer = 1)
    network.coord[:, 6:end] .+= [10; 10; 10]
    for i in eachindex(network.segForce)
        network.segForce[i] = i
    end
    for i in eachindex(network.nodeVel)
        network.nodeVel[i] = -i
    end

    selfForce = calcSelfForce(dlnParams, matParams, network)

    f1 = [
        -1.109241905058071 0.758809517064033 -0.350432387994039
        -0.118537684538429 1.035981587065609 0.917443902527180
        1.035981587065609 -0.118537684538429 0.917443902527180
        0.758809517064033 -1.109241905058071 -0.350432387994039
        -0.567011514533144 -0.567011514533144 -1.134023029066286
        -1.109241905058071 0.758809517064033 -0.350432387994039
        -0.118537684538429 1.035981587065609 0.917443902527180
        1.035981587065609 -0.118537684538429 0.917443902527180
        0.758809517064033 -1.109241905058071 -0.350432387994039
        -0.567011514533144 -0.567011514533144 -1.134023029066286
    ]
    f2 = [
        1.109241905058071 -0.758809517064033 0.350432387994039
        0.118537684538429 -1.035981587065609 -0.917443902527180
        -1.035981587065609 0.118537684538429 -0.917443902527180
        -0.758809517064033 1.109241905058071 0.350432387994039
        0.567011514533144 0.567011514533144 1.134023029066286
        1.109241905058071 -0.758809517064033 0.350432387994039
        0.118537684538429 -1.035981587065609 -0.917443902527180
        -1.035981587065609 0.118537684538429 -0.917443902527180
        -0.758809517064033 1.109241905058071 0.350432387994039
        0.567011514533144 0.567011514533144 1.134023029066286
    ]
    @test isapprox(selfForce[1], f1')
    @test isapprox(selfForce[2], f2')

    idx = rand(1:(network.numNode[1]), Int(network.numNode[1] / 2))
    selfIdx = calcSelfForce(dlnParams, matParams, network, idx)
    @test isapprox(selfForce[1][:, idx], selfIdx[1])
    @test isapprox(selfForce[2][:, idx], selfIdx[2])

    idx = rand(1:(network.numNode[1]))
    selfIdx = calcSelfForce(dlnParams, matParams, network, idx)
    @test isapprox(selfForce[1][:, idx], selfIdx[1])
    @test isapprox(selfForce[2][:, idx], selfIdx[2])

    pentagon = DislocationLoop(;
        loopType = loopPrism(),
        numSides = 5,
        nodeSide = 1,
        numLoops = 2,
        segLen = 10 * ones(5),
        slipSystemIdx = 4,
        slipSystem = slipSystems,
        label = nodeTypeDln[1; 2; 1; 2; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork(pentagon, memBuffer = 1)
    remoteForceSer = calcSegSegForce(dlnParams, matParams, network)
    remoteForcePar = calcSegSegForce(dlnParamsPar, matParams, network)
    f1 = [
        0.000330019407456 0.001179598001408 0.001509617408864
        0.001578610251880 -0.000645615383196 0.000932994868684
        0.000645615383195 -0.001578610251879 -0.000932994868684
        -0.001179598001402 -0.000330019407454 -0.001509617408857
        -0.001374647041125 0.001374647041125 -0.000000000000000
        0.000330019407456 0.001179598001408 0.001509617408865
        0.001578610251881 -0.000645615383196 0.000932994868685
        0.000645615383195 -0.001578610251879 -0.000932994868684
        -0.001179598001403 -0.000330019407455 -0.001509617408857
        -0.001374647041125 0.001374647041125 -0.000000000000000
    ]
    f2 = [
        0.000330019407454 0.001179598001403 0.001509617408857
        0.001578610251880 -0.000645615383196 0.000932994868684
        0.000645615383196 -0.001578610251880 -0.000932994868684
        -0.001179598001408 -0.000330019407456 -0.001509617408864
        -0.001374647041126 0.001374647041126 0.000000000000000
        0.000330019407454 0.001179598001402 0.001509617408856
        0.001578610251878 -0.000645615383195 0.000932994868683
        0.000645615383196 -0.001578610251881 -0.000932994868685
        -0.001179598001408 -0.000330019407456 -0.001509617408864
        -0.001374647041126 0.001374647041126 0.000000000000000
    ]
    @test isapprox(remoteForceSer[:, 1, :], f1')
    @test isapprox(remoteForceSer[:, 2, :], f2')
    @test isapprox(remoteForcePar[:, 1, :], f1')
    @test isapprox(remoteForcePar[:, 2, :], f2')

    pentagon = DislocationLoop(;
        loopType = loopPrism(),
        numSides = 5,
        nodeSide = 1,
        numLoops = 2,
        segLen = 10 * ones(5),
        slipSystemIdx = 4,
        slipSystem = slipSystems,
        label = nodeTypeDln[1; 2; 1; 2; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork(pentagon, memBuffer = 1)
    network.coord[:, 6:end] .+= [20; 20; 20]
    remoteForceSer = calcSegSegForce(dlnParams, matParams, network)
    remoteForcePar = calcSegSegForce(dlnParamsPar, matParams, network)
    f1 = [
        0.000217825304625 0.000839481257576 0.001128278735983
        0.001146987202523 -0.000492397907541 0.000704213490464
        0.000532132124482 -0.001268053797081 -0.000764723861474
        -0.000989719230921 -0.000261682028865 -0.001314774875872
        -0.001103824816939 0.001093373550371 0.000005225633284
        0.000258286008162 0.000976109417854 0.001296054439544
        0.001312903402718 -0.000551465905248 0.000792350410865
        0.000507633513349 -0.001180265566845 -0.000725717309154
        -0.000831056838483 -0.000215895485883 -0.001116145458036
        -0.001051183335307 0.001060779799870 -0.000004798232281
    ]
    f2 = [
        0.000215895485882 0.000831056838483 0.001116145458036
        0.001180265566845 -0.000507633513349 0.000725717309154
        0.000551465905248 -0.001312903402718 -0.000792350410864
        -0.000976109417853 -0.000258286008162 -0.001296054439544
        -0.001060779799870 0.001051183335308 0.000004798232281
        0.000261682028865 0.000989719230921 0.001314774875872
        0.001268053797081 -0.000532132124482 0.000764723861474
        0.000492397907540 -0.001146987202523 -0.000704213490464
        -0.000839481257576 -0.000217825304625 -0.001128278735983
        -0.001093373550372 0.001103824816939 -0.000005225633284
    ]
    @test isapprox(remoteForceSer[:, 1, :], f1')
    @test isapprox(remoteForceSer[:, 2, :], f2')
    @test isapprox(remoteForcePar[:, 1, :], f1')
    @test isapprox(remoteForcePar[:, 2, :], f2')

    idx = rand(1:(network.numSeg[1]), Int(network.numSeg[1] / 2))
    serIdx = calcSegSegForce(dlnParams, matParams, network, idx)
    @test isapprox(remoteForceSer[:, :, idx], serIdx)
    @test isapprox(remoteForcePar[:, :, idx], serIdx)

    idx = rand(1:(network.numSeg[1]))
    serIdx = calcSegSegForce(dlnParams, matParams, network, idx)
    @test isapprox(remoteForceSer[:, :, idx], serIdx)
    @test isapprox(remoteForcePar[:, :, idx], serIdx)

    hexagonPris = DislocationLoop(;
        loopType = loopPrism(),
        numSides = 6,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(6),
        slipSystemIdx = 1,
        slipSystem = slipSystems,
        label = nodeTypeDln[1; 2; 1; 2; 1; 2],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    hexagonShear = DislocationLoop(;
        loopType = loopShear(),
        numSides = 6,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(6),
        slipSystemIdx = 1,
        slipSystem = slipSystems,
        label = nodeTypeDln[1; 2; 1; 2; 1; 2],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork([hexagonPris, hexagonShear], memBuffer = 1)
    remoteForceSer = calcSegSegForce(dlnParams, matParams, network)
    remoteForcePar = calcSegSegForce(dlnParamsPar, matParams, network)

    f1 = [
        0.000148330200377 0.000867173203941 -0.001163833604695
        -0.001162438253910 0.000146934849592 -0.000868568554726
        -0.001310620554548 -0.000720386254088 0.000295117150230
        -0.000148330200377 -0.000867173203940 0.001163833604693
        0.001162438253908 -0.000146934849592 0.000868568554724
        0.001310620554548 0.000720386254089 -0.000295117150230
        0.000487619084780 0.000552951524211 -0.000501376468106
        -0.000840679364846 -0.000186252085257 -0.000213189755230
        -0.001337886532722 -0.000747434694617 0.000295225919052
        -0.000487619084782 -0.000552951524212 0.000501376468101
        0.000840679364847 0.000186252085255 0.000213189755231
        0.001337886532723 0.000747434694617 -0.000295225919053
    ]
    f2 = [
        0.000146934849592 0.000868568554725 -0.001162438253909
        -0.001163833604693 0.000148330200376 -0.000867173203941
        -0.001310620554545 -0.000720386254089 0.000295117150228
        -0.000146934849592 -0.000868568554727 0.001162438253910
        0.001163833604695 -0.000148330200376 0.000867173203942
        0.001310620554545 0.000720386254088 -0.000295117150228
        0.000480410434831 0.000546521015274 -0.000507348104804
        -0.000843647543328 -0.000196923065664 -0.000210680448989
        -0.001337886532724 -0.000747434694618 0.000295225919053
        -0.000480410434829 -0.000546521015273 0.000507348104808
        0.000843647543328 0.000196923065666 0.000210680448988
        0.001337886532722 0.000747434694618 -0.000295225919052
    ]
    @test isapprox(remoteForceSer[:, 1, :], f1')
    @test isapprox(remoteForceSer[:, 2, :], f2')
    @test isapprox(remoteForcePar[:, 1, :], f1')
    @test isapprox(remoteForcePar[:, 2, :], f2')

    PKForce = calcPKForce(regularCuboidMesh, forceDisplacement, network)
    selfForce = calcSelfForce(dlnParams, matParams, network)
    remoteForceSer = calcSegSegForce(dlnParams, matParams, network)
    remoteForcePar = calcSegSegForce(dlnParamsPar, matParams, network)
    sumForceSer = (
        selfForce[1] .+ remoteForceSer[:, 1, :] + 0.5 * PKForce,
        selfForce[2] .+ remoteForceSer[:, 2, :] + 0.5 * PKForce,
    )
    sumForcePar = (
        selfForce[1] .+ remoteForcePar[:, 1, :] + 0.5 * PKForce,
        selfForce[2] .+ remoteForcePar[:, 2, :] + 0.5 * PKForce,
    )

    totalForceSer =
        calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    totalForcePar =
        calcSegForce(dlnParamsPar, matParams, regularCuboidMesh, forceDisplacement, network)

    @test isapprox(totalForceSer[:, 1, :], sumForceSer[1])
    @test isapprox(totalForcePar[:, 2, :], sumForcePar[2])
    @test isapprox(totalForceSer[:, 1, :], totalForcePar[:, 1, :])
    @test isapprox(totalForceSer[:, 2, :], totalForcePar[:, 2, :])

    idx = rand(1:(network.numSeg[1]), Int(network.numSeg[1] / 2))
    totalForceIdx = calcSegForce(
        dlnParams,
        matParams,
        regularCuboidMesh,
        forceDisplacement,
        network,
        idx,
    )
    totalForceSer =
        calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    totalForcePar =
        calcSegForce(dlnParamsPar, matParams, regularCuboidMesh, forceDisplacement, network)
    @test isapprox(totalForceSer[:, :, idx], totalForceIdx)
    @test isapprox(totalForcePar[:, :, idx], totalForceIdx)

    idx = rand(1:(network.numSeg[1]))
    totalForceIdx = calcSegForce(
        dlnParams,
        matParams,
        regularCuboidMesh,
        forceDisplacement,
        network,
        idx,
    )
    totalForceSer =
        calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    totalForcePar =
        calcSegForce(dlnParamsPar, matParams, regularCuboidMesh, forceDisplacement, network)
    @test isapprox(totalForceSer[:, :, idx], totalForceIdx)
    @test isapprox(totalForcePar[:, :, idx], totalForceIdx)

    # In-place functions.
    numSeg = network.numSeg[1]
    idx = rand(1:numSeg)
    self = calcSelfForce(dlnParams, matParams, network)
    selfIdx = calcSelfForce(dlnParams, matParams, network, idx)
    @test isapprox(self[1][:, idx], selfIdx[1])
    @test isapprox(self[2][:, idx], selfIdx[2])
    network.segForce .= 0
    calcSelfForce!(dlnParams, matParams, network)
    @test isapprox(self[1], network.segForce[:, 1, 1:numSeg])
    isapprox(self[2], network.segForce[:, 2, 1:numSeg])
    network.segForce .= 0
    calcSelfForce!(dlnParams, matParams, network, idx)
    @test isapprox(self[1][:, idx], network.segForce[:, 1, idx])
    @test isapprox(self[2][:, idx], network.segForce[:, 2, idx])

    ser = calcSegSegForce(dlnParams, matParams, network)
    network.segForce .= 0
    calcSegSegForce!(dlnParams, matParams, network)
    @test isapprox(network.segForce[:, :, 1:numSeg], ser)
    network.segForce .= 0
    calcSegSegForce!(dlnParams, matParams, network, idx)
    @test isapprox(network.segForce[:, :, idx], ser[:, :, idx])

    par = calcSegSegForce(dlnParamsPar, matParams, network)
    network.segForce .= 0
    calcSegSegForce!(dlnParamsPar, matParams, network)
    @test isapprox(network.segForce[:, :, 1:numSeg], par)
    network.segForce .= 0
    calcSegSegForce!(dlnParamsPar, matParams, network, idx)
    @test isapprox(network.segForce[:, :, idx], par[:, :, idx])
    @test isapprox(par, ser)

    tot = calcSegForce(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    @test isapprox(tot[:, 1, :], self[1] + ser[:, 1, :] + 0.5 * PKForce)
    @test isapprox(tot[:, 2, :], self[2] + ser[:, 2, :] + 0.5 * PKForce)
    network.segForce .= 0
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network)
    @test isapprox(tot, network.segForce[:, :, 1:numSeg])

    tot = calcSegForce(
        dlnParams,
        matParams,
        regularCuboidMesh,
        forceDisplacement,
        network,
        idx,
    )
    @test isapprox(tot[:, 1, :], self[1][:, idx] + ser[:, 1, idx] + 0.5 * PKForce[:, idx])
    @test isapprox(tot[:, 2, :], self[2][:, idx] + ser[:, 2, idx] + 0.5 * PKForce[:, idx])
    network.segForce .= 0
    calcSegForce!(dlnParams, matParams, regularCuboidMesh, forceDisplacement, network, idx)
    @test isapprox(tot, network.segForce[:, :, idx])

    network.coord .+= dx * dy * dz
    PKForce = calcPKForce(regularCuboidMesh, forceDisplacement, network)
    @test isapprox(PKForce, zeros(3, numSeg))
end

@testset "Calculating sigma_hat" begin
    fileDislocationParameters = "../inputs/simParams/sampleDislocationParameters.json"
    fileMaterialParameters = "../inputs/simParams/sampleMaterialParameters.json"
    fileFEMParameters = "../inputs/simParams/sampleFEMParameters.json"
    fileIntegrationParameters = "../inputs/simParams/sampleIntegrationParameters.json"
    fileSlipSystem = "../data/slipSystems/BCC.json"
    fileDislocationLoop = "../inputs/dln/samplePrismShear.json"
    fileIntVar = "../inputs/simParams/sampleIntegrationTime.json"
    dlnParams, matParams, femParams, intParams, slipSystems, dislocationLoop =
        loadParameters(
            fileDislocationParameters,
            fileMaterialParameters,
            fileFEMParameters,
            fileIntegrationParameters,
            fileSlipSystem,
            fileDislocationLoop,
        )

    regularCuboidMesh = buildMesh(matParams, femParams)
    numFEMNode = regularCuboidMesh.numNode

    f = spzeros(3 * numFEMNode)
    f[[112, 118, 133, 141, 213, 244, 262, 272, 317]] = [
        0.43048187784858616,
        0.22724536603830137,
        0.4340867899691503,
        0.6660863546953892,
        0.30358515797696106,
        0.2945958951093859,
        0.7278367502911502,
        0.7095924334694701,
        0.1642050526375538,
    ]
    fHat = spzeros(3 * numFEMNode)
    fHat[[32, 48, 55, 88, 138, 148, 191, 230, 253, 335]] = [
        0.09706224225842108,
        0.07773687633248638,
        0.13682398802299178,
        0.4752286167553166,
        0.7423196193496164,
        0.8286077556473421,
        0.7023632196408749,
        0.9813639162461198,
        0.5296701796678411,
        0.5523797553266823,
    ]
    u = spzeros(3 * numFEMNode)
    u[[30, 127, 195, 221, 316, 325, 338, 348, 370]] = [
        0.8792592573507609,
        0.8430664083925272,
        0.4711050560756602,
        0.4860071865093816,
        0.7905698600135145,
        0.39047211692578077,
        0.6545538020629462,
        0.5446700211111557,
        0.8865721648558644,
    ]
    uHat = spzeros(3 * numFEMNode)
    uHat[[67, 68, 69, 195, 217, 226, 229, 256, 281, 293, 309, 342]] = [
        0.5231621885339968,
        0.5771429489788034,
        0.7151190318538345,
        0.7283662326812077,
        0.6314274719472075,
        0.9814688915693632,
        0.5672795171250207,
        0.002712918060655989,
        0.1788941754890383,
        0.188299784057536,
        0.8489027048214433,
        0.029995302953659708,
    ]

    forceDisplacement =
        ForceDisplacement(nothing, uHat * 1000, u * 1000, nothing, fHat * 1000, f * 1000)
    regularCuboidMesh.coord[:, 91]
    σHat = calc_σHat(regularCuboidMesh, forceDisplacement, [1800.0, 1400.0, 200.0])
    σHatTest = [
        -1.5186579125763586 -0.6876907109455002 -1.2291546581681494
        -0.6876907109455002 -1.5861338631323667 -0.35243484259464375
        -1.2291546581681494 -0.35243484259464375 -0.8481461758746458
    ]
    @test isapprox(σHat, σHatTest)

    σHat = calc_σHat(regularCuboidMesh, forceDisplacement, [1893.0, 1356.0, 345.0])
    σHatTest = [
        -0.5560725333509386 -0.169126060231533 -0.9758794760722914
        -0.169126060231533 -0.4428115859951268 0.11892368539106829
        -0.9758794760722914 0.11892368539106829 -0.5401850945259058
    ]
    @test isapprox(σHat, σHatTest)
end
