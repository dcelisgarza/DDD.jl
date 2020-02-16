using DDD
using Test
using Plots
gr()
cd(@__DIR__)

@testset "Plot nodes" begin
    params = "../inputs/simParams/sampleParams.csv"
    slipsys = "../data/slipSystems/bcc.csv"
    source = "../inputs/dln/sampleDln.csv"
    dlnParams, matParams, intParams, slipSystems, loops = loadParams(
        params,
        slipsys,
        source,
    )
    network = zero(DislocationNetwork)
    makeNetwork!(network, loops)
    function sumNodes(loops)
        totalNodes = 0
        for i in eachindex(loops)
            totalNodes += loops[i].numSides * loops[i].nodeSide *
                          loops[i].numLoops
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
        for i in eachindex(loops)
            plotNodes!(fig, loops[i], m = 1, l = 3, legend = false)
            totalNodes += loops[i].numSides * loops[i].nodeSide *
                          loops[i].numLoops
            @test fig.n == totalNodes
        end
    end
    function plotLoops(loops)
        for i in eachindex(loops)
            fig = plotNodes(loops[i], m = 1, l = 3, legend = false)
            @test fig.n == loops[i].numSides * loops[i].nodeSide *
                           loops[i].numLoops
        end
    end
    fig = plot()
    plotLoops(fig, loops)
    plotLoops(loops)
end
