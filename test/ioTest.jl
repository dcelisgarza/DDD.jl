using DDD
using Test, FileIO

cd(@__DIR__)
@testset "Compare loading input and output parameters" begin
    # Load and create.
    fileDislocationParameters = "./testData/sampleDislocationParameters.json"
    fileMaterialParameters = "./testData/sampleMaterialParameters.json"
    fileIntegrationParameters = "./testData/sampleIntegrationParameters.json"
    fileSlipSystem = "./testData/BCC.json"
    fileDislocationLoop = "./testData/samplePrismShear.json"
    dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParametersJSON(
        fileDislocationParameters,
        fileMaterialParameters,
        fileIntegrationParameters,
        fileSlipSystem,
        fileDislocationLoop,
    )
    fileDislocationLoop = "./testData/samplePrismShear.json"
    fileIntegTime = "./testData/sampleIntegrationTime.json"
    integTime = loadIntegrationTimeJSON(fileIntegTime)
    network = DislocationNetwork(dislocationLoop, memBuffer = 1)
    # Dump simulationJSON.
    paramDumpJSON = "./testData/sampleDump.json"
    paramDumpJLD2 = "./testData/sampleDump.jld2"
    saveJSON(
        paramDumpJSON,
        dlnParams,
        matParams,
        intParams,
        slipSystems,
        dislocationLoop,
        integTime,
    )
    save(
        paramDumpJLD2,
        "dlnParams",
        dlnParams,
        "matParams",
        matParams,
        "intParams",
        intParams,
        "slipSystems",
        slipSystems,
        "dislocationLoop",
        dislocationLoop,
        "integTime",
        integTime,
    )
    networkDumpJSON = "./testData/sampleNetwork.json"
    networkDumpJLD2 = "./testData/sampleNetwork.jld2"
    saveJSON(networkDumpJSON, network)
    save(networkDumpJLD2, "network", network)
    # Reload simulationJSON.
    simulationJSON = loadJSON(paramDumpJSON)
    dlnParams2 = loadDislocationParametersJSON(simulationJSON[1])
    matParams2 = loadMaterialParametersJSON(simulationJSON[2])
    intParams2 = loadIntegrationParametersJSON(simulationJSON[3])
    slipSystems2 = loadSlipSystemJSON(simulationJSON[4])
    dislocationLoop2 = zeros(DislocationLoop, length(simulationJSON[5]))
    for i in eachindex(dislocationLoop2)
        dislocationLoop2[i] = loadDislocationLoopJSON(simulationJSON[5][i], slipSystems2)
    end
    network2 = loadNetworkJSON(networkDumpJSON)
    integTime2 = loadIntegrationTimeJSON(simulationJSON[6])
    # Reload simulationJLD2.
    dlnParams3, matParams3, intParams3, slipSystems3, dislocationLoop3, integTime3 = load(
        paramDumpJLD2,
        "dlnParams",
        "matParams",
        "intParams",
        "slipSystems",
        "dislocationLoop",
        "integTime",
    )
    network3 = load(networkDumpJLD2, "network")
    # Ensure the data is all the same.
    @test compStruct(dlnParams, dlnParams2; verbose = true)
    @test compStruct(matParams, matParams2; verbose = true)
    @test compStruct(intParams, intParams2; verbose = true)
    @test compStruct(slipSystems, slipSystems2; verbose = true)
    @test compStruct(dislocationLoop, dislocationLoop2; verbose = true)
    @test compStruct(network, network2; verbose = true)
    @test compStruct(integTime, integTime2; verbose = true)

    @test compStruct(dlnParams, dlnParams3; verbose = true)
    @test compStruct(matParams, matParams3; verbose = true)
    @test compStruct(intParams, intParams3; verbose = true)
    @test compStruct(slipSystems, slipSystems3; verbose = true)
    @test compStruct(dislocationLoop, dislocationLoop3; verbose = true)
    @test compStruct(network, network3; verbose = true)
    @test compStruct(integTime, integTime3; verbose = true)

    rm(paramDumpJSON; force = true)
    rm(paramDumpJLD2; force = true)
    rm(networkDumpJSON; force = true)
    rm(networkDumpJLD2; force = true)
end
