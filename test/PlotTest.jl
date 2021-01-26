using DDD
using Test
using Plots
gr()
cd(@__DIR__)

@testset "Plot nodes" begin
    # Load and create.
    fileSlipSystem = "./testData/BCC.json"
    fileDislocationLoop = "./testData/sampleDislocation.json"
    dictSlipSystem = loadJSON(fileSlipSystem)
    slipSystems = loadSlipSystemJSON(dictSlipSystem)
    # There can be multiple dislocations per simulation parameters.
    dictDislocationLoop = loadJSON(fileDislocationLoop)
    loops = zeros(DislocationLoop, length(dictDislocationLoop))
    if typeof(dictDislocationLoop) <: AbstractArray
        loops = zeros(DislocationLoop, length(dictDislocationLoop))
        for i in eachindex(loops)
            loops[i] = loadDislocationLoopJSON(dictDislocationLoop[i], slipSystems)
        end
    else
        loops = loadDislocationLoopJSON(dictDislocationLoop, slipSystems)
    end
    network = DislocationNetwork(loops)

    function sumNodes(loops)
        totalNodes = 0
        for i in eachindex(loops)
            totalNodes += loops[i].numSides * loops[i].nodeSide * loops[i].numLoops
        end
        return totalNodes
    end
    totalNodes = sumNodes(loops)
    # Modifying current plot.
    fig = plot()
    plotNodes!(fig, network, m = 1, l = 3, legend = false)
    @test fig.n == totalNodes
    # Creating new plot.
    fig2 = plotNodes(network, m = 1, l = 3, legend = false)
    @test fig2.n == totalNodes

    function plotLoops(fig, loops)
        totalNodes = 0
        return for i in eachindex(loops)
            plotNodes!(fig, loops[i], m = 1, l = 3, legend = false)
            totalNodes += loops[i].numSides * loops[i].nodeSide * loops[i].numLoops
            @test fig.n == totalNodes
        end
    end
    function plotLoops(loops)
        return for i in eachindex(loops)
            fig = plotNodes(loops[i], m = 1, l = 3, legend = false)
            @test fig.n == loops[i].numSides * loops[i].nodeSide * loops[i].numLoops
        end
    end
    fig = plot()
    plotLoops(fig, loops)
    plotLoops(loops)

    DictMaterialParameters = loadJSON("./testData/sampleMaterialParameters.json")
    DictFEMParameters = loadJSON("./testData/sampleFEMParameters.json")
    matParams = loadMaterialParametersJSON(DictMaterialParameters)
    femParams = loadFEMParametersJSON(DictFEMParameters)
    regularCuboidMesh = buildMesh(matParams, femParams)

    fig = plot()
    plotNodes!(fig, regularCuboidMesh, network, m = 1, l = 3, legend = false)
    @test fig.n == totalNodes + 6
    # Creating new plot.
    fig2 = plotNodes(regularCuboidMesh, network, m = 1, l = 3, legend = false)
    @test fig2.n == totalNodes + 6

    fig3 = plotFEDomain(regularCuboidMesh)
    @test fig3.n == 19
end
