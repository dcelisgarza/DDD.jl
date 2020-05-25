using DDD
using Test
using Plots
gr()
cd(@__DIR__)

@testset "Plot nodes" begin
    # Load and create.
    fileSlipSystem = "../data/slipSystems/SlipSystems.JSON"
    fileDislocationLoop = "../inputs/dln/sampleDislocation.JSON"
    dictSlipSystem = load(fileSlipSystem)
    slipSystems = loadSlipSystem(dictSlipSystem[1])
    # There can be multiple dislocations per simulation parameters.
    dictDislocationLoop = load(fileDislocationLoop)
    loops = zeros(DislocationLoop, length(dictDislocationLoop))
    for i in eachindex(loops)
        loops[i] = loadDislocationLoop(dictDislocationLoop[i], slipSystems)
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
end
