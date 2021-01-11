using DDD
using Test, StaticArrays, SparseArrays
cd(@__DIR__)

@testset "Calculating sigma_hat" begin
    fileDislocationParameters = "../inputs/simParams/sampleDislocationParameters.json"
    fileMaterialParameters = "../inputs/simParams/sampleMaterialParameters.json"
    fileFEMParameters = "../inputs/simParams/sampleFEMParameters.json"
    fileIntegrationParameters = "../inputs/simParams/sampleIntegrationParameters.json"
    fileSlipSystem = "../data/slipSystems/BCC.json"
    fileDislocationLoop = "../inputs/dln/samplePrismShear.json"
    fileIntVar = "../inputs/simParams/sampleIntegrationTime.json"
    dlnParams, matParams, femParams, intParams, slipSystems, dislocationLoop =
        loadParametersJSON(
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
    f[[112, 118, 133, 141, 213, 244, 262, 272, 317]] = 
    [
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
    fHat[[32, 48, 55, 88, 138, 148, 191, 230, 253, 335]] = 
    [
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
    u[[30, 127, 195, 221, 316, 325, 338, 348, 370]] = 
    [
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
    uHat[[91, 126, 130, 195, 217, 226, 229, 256, 281, 293, 309, 342]] = 
    [
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

    forceDisplacement = ForceDisplacement(u * 1000, f * 1000, uHat * 1000, fHat * 1000)

    @time σHat = calc_σHat(regularCuboidMesh, forceDisplacement, [1575.0, 985.0, 1341.0])
    σHatTest = [
        -0.012540380057364 -0.084737139530472 0
        -0.084737139530472 -0.032246691576078 0.015024315519587
        0 0.015024315519587 -0.012540380057364
    ]
    @test isapprox(σHat, σHatTest)
end
