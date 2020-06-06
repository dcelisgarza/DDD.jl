using DDD
using Test

cd(@__DIR__)
@testset "BCC mobility" begin
    fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
    fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
    fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
    fileSlipSystem = "../data/slipSystems/BCC.JSON"
    fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
    dlnParams, matParams, intParams, slipSystems, missing = loadParams(
        fileDislocationP,
        fileMaterialP,
        fileIntegrationP,
        fileSlipSystem,
        fileDislocationLoop,
    )

    prismPentagon = DislocationLoop(
        loopPrism();
        numSides = 5,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(5),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )

    shearHexagon = DislocationLoop(
        loopShear();
        numSides = 6,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(6),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork([shearHexagon, prismPentagon], memBuffer = 1)
    calcSegForce!(dlnParams, matParams, network)

    force, vel = dlnMobility(dlnParams, matParams, network)
    testVel = [
        -0.036175006034684 -0.036175006034656 0.088938030619749
        0.060538472465831 0.060538472465831 0.121076825150312
        0.071350296347816 0.071350296347788 -0.018587478143454
        0.036175184446079 0.036175184446051 -0.088937984731268
        -0.060538377421600 -0.060538377421600 -0.121076879461304
        -0.071350353768921 -0.071350353768893 0.018587492151485
        -0.000039076066951 -0.000039076038383 0.000039076070642
        0.000014926772741 0.000014926761829 -0.000014926738200
        0.000048298670216 0.000048298634903 -0.000048298652559
        0.000014926727728 0.000014926716815 -0.000014926751357
        -0.000039076030376 -0.000039076001807 0.000039075998116
    ]

    testForce = [
        -0.228724254571166 -0.226639000305794 0.900128776494974
        0.686441227083849 0.686441227083850 1.372882454167666
        0.677022353942198 0.674937099676816 -0.003532577752966
        0.228724254571162 0.226639000305798 -0.900128776495018
        -0.686441227083830 -0.686441227083827 -1.372882454167693
        -0.677022353942161 -0.674937099676790 0.003532577752984
        -0.543143464917538 1.326707393542668 0.784736210352455
        0.991807802754387 0.277588328715755 1.268948329254914
        1.156114370762920 -1.155148398080229 -0.000482986341346
        -0.277289793905604 -0.991509267944234 -1.269246864065067
        -1.327488914694218 0.542361943765987 -0.783954689200902
    ]

    @test isapprox(vel', testVel, rtol = 1e-5)
    @test isapprox(force', testForce)

    idx = rand(1:(network.numNode))
    forceIdx, velIdx = dlnMobility(dlnParams, matParams, network, idx)
    @test isapprox(force[:, idx], forceIdx)
    @test isapprox(vel[:, idx], velIdx)

    idx = rand(1:(network.numNode), 5)
    forceIdx, velIdx = dlnMobility(dlnParams, matParams, network, idx)
    @test isapprox(force[:, idx], forceIdx)
    @test isapprox(vel[:, idx], velIdx)
end
