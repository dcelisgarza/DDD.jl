using DDD
using Test
cd(@__DIR__)

@testset "Forces" begin
    fileDislocationParameters = "../inputs/simParams/sampleDislocationParameters.JSON"
    fileMaterialParameters = "../inputs/simParams/sampleMaterialParameters.JSON"
    fileIntegrationParameters = "../inputs/simParams/sampleIntegrationParameters.JSON"
    fileSlipSystem = "../data/slipSystems/BCC.JSON"
    fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
    dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
        fileDislocationParameters,
        fileMaterialParameters,
        fileIntegrationParameters,
        fileSlipSystem,
        fileDislocationLoop,
    )
    dlnParamsPar = DislocationParameters(;
        coreRad = dlnParams.coreRad,
        coreRadMag = dlnParams.coreRadMag,
        minSegLen = dlnParams.minSegLen,
        maxSegLen = dlnParams.maxSegLen,
        minArea = dlnParams.minArea,
        maxArea = dlnParams.maxArea,
        maxConnect = dlnParams.maxConnect,
        remesh = dlnParams.remesh,
        collision = dlnParams.collision,
        separation = dlnParams.separation,
        virtualRemesh = dlnParams.virtualRemesh,
        parCPU = true,
        parGPU = dlnParams.parGPU,
        edgeDrag = dlnParams.edgeDrag,
        screwDrag = dlnParams.screwDrag,
        climbDrag = dlnParams.climbDrag,
        lineDrag = dlnParams.lineDrag,
        mobility = dlnParams.mobility,
    )

    pentagon = DislocationLoop(
        loopPrism();
        numSides = 5,
        nodeSide = 1,
        numLoops = 2,
        segLen = 10 * ones(5),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1],
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

    idx = rand(1:(network.numNodeSegConnect[1]), Int(network.numNodeSegConnect[1] / 2))
    selfIdx = calcSelfForce(dlnParams, matParams, network, idx)
    @test isapprox(selfForce[1][:, idx], selfIdx[1])
    @test isapprox(selfForce[2][:, idx], selfIdx[2])

    idx = rand(1:(network.numNodeSegConnect[1]))
    selfIdx = calcSelfForce(dlnParams, matParams, network, idx)
    @test isapprox(selfForce[1][:, idx], selfIdx[1])
    @test isapprox(selfForce[2][:, idx], selfIdx[2])

    pentagon = DislocationLoop(
        loopPrism();
        numSides = 5,
        nodeSide = 1,
        numLoops = 2,
        segLen = 10 * ones(5),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1],
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

    pentagon = DislocationLoop(
        loopPrism();
        numSides = 5,
        nodeSide = 1,
        numLoops = 2,
        segLen = 10 * ones(5),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1],
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

    idx = rand(1:(network.numNodeSegConnect[2]), Int(network.numNodeSegConnect[2] / 2))
    serIdx = calcSegSegForce(dlnParams, matParams, network, idx)
    @test isapprox(remoteForceSer[:, :, idx], serIdx)
    @test isapprox(remoteForcePar[:, :, idx], serIdx)

    idx = rand(1:(network.numNodeSegConnect[2]))
    serIdx = calcSegSegForce(dlnParams, matParams, network, idx)
    @test isapprox(remoteForceSer[:, :, idx], serIdx)
    @test isapprox(remoteForcePar[:, :, idx], serIdx)

    hexagonPris = DislocationLoop(
        loopPrism();
        numSides = 6,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(6),
        slipSystem = 1,
        _slipPlane = slipSystems.slipPlane[:, 1],
        _bVec = slipSystems.bVec[:, 1],
        label = nodeType[1; 2; 1; 2; 1; 2],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    hexagonShear = DislocationLoop(
        loopShear();
        numSides = 6,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(6),
        slipSystem = 1,
        _slipPlane = slipSystems.slipPlane[:, 1],
        _bVec = slipSystems.bVec[:, 1],
        label = nodeType[1; 2; 1; 2; 1; 2],
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

    selfForce = calcSelfForce(dlnParams, matParams, network)
    remoteForceSer = calcSegSegForce(dlnParams, matParams, network)
    remoteForcePar = calcSegSegForce(dlnParamsPar, matParams, network)
    sumForceSer =
        (selfForce[1] .+ remoteForceSer[:, 1, :], selfForce[2] .+ remoteForceSer[:, 2, :])
    sumForcePar =
        (selfForce[1] .+ remoteForcePar[:, 1, :], selfForce[2] .+ remoteForcePar[:, 2, :])
    totalForceSer = calcSegForce(dlnParams, matParams, network)
    totalForcePar = calcSegForce(dlnParamsPar, matParams, network)

    @test isapprox(totalForceSer[:, 1, :], sumForceSer[1])
    @test isapprox(totalForcePar[:, 2, :], sumForcePar[2])
    @test isapprox(totalForceSer[:, 1, :], totalForcePar[:, 1, :])
    @test isapprox(totalForceSer[:, 2, :], totalForcePar[:, 2, :])

    idx = rand(1:(network.numNodeSegConnect[2]), Int(network.numNodeSegConnect[2] / 2))
    totalForceIdx = calcSegForce(dlnParams, matParams, network, idx)
    totalForceSer = calcSegForce(dlnParams, matParams, network)
    totalForcePar = calcSegForce(dlnParamsPar, matParams, network)
    @test isapprox(totalForceSer[:, :, idx], totalForceIdx)
    @test isapprox(totalForcePar[:, :, idx], totalForceIdx)

    idx = rand(1:(network.numNodeSegConnect[2]))
    totalForceIdx = calcSegForce(dlnParams, matParams, network, idx)
    totalForceSer = calcSegForce(dlnParams, matParams, network)
    totalForcePar = calcSegForce(dlnParamsPar, matParams, network)
    @test isapprox(totalForceSer[:, :, idx], totalForceIdx)
    @test isapprox(totalForcePar[:, :, idx], totalForceIdx)

    # In-place functions.
    numSeg = network.numNodeSegConnect[2]
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

    tot = calcSegForce(dlnParams, matParams, network)
    @test isapprox(tot[:, 1, :], self[1] + ser[:, 1, :])
    @test isapprox(tot[:, 2, :], self[2] + ser[:, 2, :])
    network.segForce .= 0
    calcSegForce!(dlnParams, matParams, network)
    @test isapprox(tot, network.segForce[:, :, 1:numSeg])

    tot = calcSegForce(dlnParams, matParams, network, idx)
    @test isapprox(tot[:, 1, :], self[1][:, idx] + ser[:, 1, idx])
    @test isapprox(tot[:, 2, :], self[2][:, idx] + ser[:, 2, idx])
    network.segForce .= 0
    calcSegForce!(dlnParams, matParams, network, idx)
    @test isapprox(tot, network.segForce[:, :, idx])
end
