using DDD
using Test, FileIO

cd(@__DIR__)
@testset "Compare loading input and output parameters" begin
    # Load and create.
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

    missing, missing, missing, missing, missing, dislocationLoop2 = loadParameters(
        fileDislocationParameters,
        fileMaterialParameters,
        fileFEMParameters,
        fileIntegrationParameters,
        fileSlipSystem,
        "./testData/sampleDislocation.json",
    )

    @test typeof(dislocationLoop2) <: DislocationLoop

    fileDislocationLoop = "./testData/samplePrismShear.json"
    fileIntegTime = "./testData/sampleIntegrationTime.json"
    integTime = loadIntegrationTime(fileIntegTime)
    network = DislocationNetwork(dislocationLoop, memBuffer = 1)
    # Dump simulationJSON.
    paramDumpJSON = "./testData/sampleDump.json"
    paramDumpJLD2 = "./testData/sampleDump.jld2"
    saveJSON(
        paramDumpJSON,
        dlnParams,
        matParams,
        femParams,
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
        "femParams",
        femParams,
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
    dlnParams2 = loadDislocationParameters(simulationJSON[1])
    matParams2 = loadMaterialParameters(simulationJSON[2])
    femParams2 = loadFEMParameters(simulationJSON[3])
    intParams2 = loadIntegrationParameters(simulationJSON[4])
    slipSystems2 = loadSlipSystem(simulationJSON[5])
    dislocationLoop2 = zeros(DislocationLoop, length(simulationJSON[6]))
    for i in eachindex(dislocationLoop2)
        dislocationLoop2[i] = loadDislocationLoop(simulationJSON[6][i], slipSystems2)
    end
    network2 = loadNetwork(networkDumpJSON)
    integTime2 = loadIntegrationTime(simulationJSON[7])
    # Reload simulationJLD2.
    dlnParams3,
    matParams3,
    femParams3,
    intParams3,
    slipSystems3,
    dislocationLoop3,
    integTime3 = load(
        paramDumpJLD2,
        "dlnParams",
        "matParams",
        "femParams",
        "intParams",
        "slipSystems",
        "dislocationLoop",
        "integTime",
    )
    network3 = load(networkDumpJLD2, "network")
    # Ensure the data is all the same.
    @test compStruct(dlnParams, dlnParams2; verbose = true)
    @test compStruct(matParams, matParams2; verbose = true)
    @test compStruct(femParams, femParams2; verbose = true)
    @test compStruct(intParams, intParams2; verbose = true)
    @test compStruct(slipSystems, slipSystems2; verbose = true)
    @test compStruct(dislocationLoop, dislocationLoop2; verbose = true)
    @test compStruct(network, network2; verbose = true)
    @test compStruct(integTime, integTime2; verbose = true)

    @test compStruct(dlnParams, dlnParams3; verbose = true)
    @test compStruct(matParams, matParams3; verbose = true)
    @test compStruct(femParams, femParams3; verbose = true)
    @test compStruct(intParams, intParams3; verbose = true)
    @test compStruct(slipSystems, slipSystems3; verbose = true)
    @test compStruct(dislocationLoop, dislocationLoop3; verbose = true)
    @test compStruct(network, network3; verbose = true)
    @test compStruct(integTime, integTime3; verbose = true)

    sampleRegCubMeshDump = "./testData/sampleRegCubMesh.jld2"
    regularCuboidMesh = buildMesh(matParams, femParams)
    save(sampleRegCubMeshDump, "mesh", regularCuboidMesh)
    regularCuboidMesh2 = load(sampleRegCubMeshDump, "mesh")
    @test compStruct(regularCuboidMesh, regularCuboidMesh2; verbose = true)

    rm(paramDumpJSON; force = true)
    rm(paramDumpJLD2; force = true)
    rm(networkDumpJSON; force = true)
    rm(networkDumpJLD2; force = true)
    rm(sampleRegCubMeshDump; force = true)
end
