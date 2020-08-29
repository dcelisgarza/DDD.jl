using DDD
using Test

cd(@__DIR__)
@testset "Compare loading input and output parameters" begin
    # Load and create.
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
    fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
    fileIntegTime = "../inputs/simParams/sampleIntegrationTime.JSON"
    integTime = loadIntegrationTime(fileIntegTime)
    network = DislocationNetwork(dislocationLoop, memBuffer = 1)
    # Dump simulation.
    paramDump = "../outputs/simParams/sampleDump.JSON"
    save(
        paramDump,
        dlnParams,
        matParams,
        intParams,
        slipSystems,
        dislocationLoop,
        integTime,
    )
    networkDump = "../outputs/dln/sampleNetwork.JSON"
    save(networkDump, network)
    # Reload simulation.
    simulation = load(paramDump)
    dlnParams2 = loadDislocationParameters(simulation[1])
    matParams2 = loadMaterialParameters(simulation[2])
    intParams2 = loadIntegrationParameters(simulation[3])
    slipSystems2 = loadSlipSystem(simulation[4])
    dislocationLoop2 = zeros(DislocationLoop, length(simulation[5]))
    for i in eachindex(dislocationLoop2)
        dislocationLoop2[i] = loadDislocationLoop(simulation[5][i], slipSystems2)
    end
    network2 = loadNetwork(networkDump)
    integTime2 = loadIntegrationTime(simulation[6])
    # Ensure the data is all the same.
    @test compStruct(dlnParams, dlnParams2; verbose = true)
    @test compStruct(matParams, matParams2; verbose = true)
    @test compStruct(intParams, intParams2; verbose = true)
    @test compStruct(slipSystems, slipSystems2; verbose = true)
    @test compStruct(dislocationLoop, dislocationLoop2; verbose = true)
    @test compStruct(network, network2; verbose = true)
    @test compStruct(integTime, integTime2; verbose = true)
end
