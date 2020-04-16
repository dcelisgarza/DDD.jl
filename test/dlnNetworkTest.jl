using DDD
using Test

import LinearAlgebra: dot, cross, norm
import Statistics: mean
import Random: seed!

import DDD:
    segEdge,
    segEdgeN,
    segScrew,
    makeSegment,
    idxLabel,
    idxCond,
    dataCond,
    coordIdx,
    coordLbl,
    inclusiveComparison
cd(@__DIR__)

@testset "Generate single segments" begin
    fileSlipSystem = "../data/slipSystems/SlipSystems.JSON"
    dictSlipSystem = load(fileSlipSystem)
    slipSystems = loadSlipSystem(dictSlipSystem[1])
    slipSysInt = 1
    slipPlane = slipSystems.slipPlane[1, :]
    bVec = slipSystems.bVec[1, :]
    edge = makeSegment(segEdge(), slipPlane, bVec)
    edgeN = makeSegment(segEdgeN(), slipPlane, bVec)
    screw = makeSegment(segScrew(), slipPlane, bVec)
    @test abs(dot(edge, screw)) < eps(Float64)
    @test abs(dot(edgeN, screw)) < eps(Float64)
    @test abs(dot(edge, bVec)) < eps(Float64)
    @test abs(dot(edgeN, bVec)) < eps(Float64)
    @test isapprox(edge, cross(slipPlane, bVec) ./ norm(cross(slipPlane, bVec)))
    @test isapprox(norm(edge), norm(screw))
    @test isapprox(norm(edge), 1.0)
end

@testset "Dislocation indexing functions" begin
    cnd = [==, >=, <=, <, >, !=]
    numNode = 10
    numSeg = 20
    links = zeros(Int64, numSeg, 2)
    bVec = zeros(numSeg, 3)
    slipPlane = zeros(numSeg, 3)
    coord = zeros(numNode, 3)
    label = zeros(nodeType, numNode)
    lenLinks = size(links, 1)
    [links[i, :] = [i, i + lenLinks] for i = 1:lenLinks]
    [bVec[i, :] = [i, i + lenLinks, i + 2 * lenLinks] for i = 1:lenLinks]
    [slipPlane[i, :] = -[i, i + lenLinks, i + 2 * lenLinks] for i = 1:lenLinks]
    lenLabel = length(label)
    [label[i + 1] = mod(i, 6) - 1 for i = 0:(lenLabel - 1)]
    [
        coord[i, :] = convert.(Float64, [i, i + lenLabel, i + 2 * lenLabel])
        for i = 1:length(label)
    ]
    network = DislocationNetwork(
        links = links,
        slipPlane = slipPlane,
        bVec = bVec,
        coord = coord,
        label = label,
        segForce = zeros(size(coord)),
        numNode = convert(Int64, numNode),
        numSeg = convert(Int64, numSeg),
        maxConnect = convert(Int64, 0),
    )
    @test isequal(network.label[1], -1)
    @test isequal(-1, network.label[1])
    @test -1 == network.label[1]
    @test -2.0 < network.label[1]
    rnd = rand(-1:numNode)
    @test idxLabel(network, rnd) == findall(x -> x == rnd, label)
    @test coordLbl(network, rnd) == coord[findall(x -> x == rnd, label), :]
    @test idxCond(network, :label, inclusiveComparison, 0, 1) == [2; 3; 8; 9]
    @test idxCond(network, :label, rnd; condition = <=) ==
          findall(x -> x <= rnd, label)
    rnd = rand(1:numSeg)
    @test idxCond(network, :bVec, rnd; condition = cnd[1]) ==
          findall(x -> cnd[1](x, rnd), bVec)
    @test dataCond(network, :slipPlane, rnd; condition = cnd[2]) ==
          slipPlane[findall(x -> cnd[2](x, rnd), slipPlane)]
    rnd = rand(-1:1)
    @test dataCond(network, :slipPlane, :bVec, rnd; condition = cnd[3]) ==
          slipPlane[findall(x -> cnd[3](x, rnd), bVec)]
    rnd = rand(1:numNode)
    @test dataCond(network, :coord, :label, rnd; condition = cnd[4]) ==
          coord[findall(x -> cnd[4](x, rnd), label), :]
    col = rand(1:3)
    @test idxCond(network, :bVec, col, rnd; condition = cnd[5]) ==
          findall(x -> cnd[5](x, rnd), bVec[:, col])
    @test dataCond(network, :bVec, col, rnd; condition = cnd[6]) ==
          bVec[findall(x -> cnd[6](x, rnd), bVec[:, col]), :]
    rnd = rand(-1:1)
    @test dataCond(network, :slipPlane, :bVec, col, rnd; condition = cnd[1]) ==
          slipPlane[findall(x -> cnd[1](x, rnd), bVec[:, col]), :]
    rnd = rand(1:numNode)
    @test coordIdx(network, rnd) == coord[rnd, :]
    rnd = rand(1:numNode, numNode)
    @test coordIdx(network, rnd) == coord[rnd, :]
end

@testset "Loop generation" begin

    slipfile = "../data/slipSystems/SlipSystems.JSON"
    loopfile = "../inputs/dln/sampleDislocation.JSON"

    dictSlipSystem = load(slipfile)
    slipSystems = loadSlipSystem(dictSlipSystem[1])

    dictDislocationLoop = load(loopfile)
    loops = zeros(DislocationLoop, length(dictDislocationLoop))
    for i in eachindex(loops)
        loops[i] = loadDislocationLoop(dictDislocationLoop[i], slipSystems)
    end
    # Check that the midpoint of the loops is at (0,0,0)
    for i in eachindex(loops)
        @test mean(loops[i].coord) <
              maximum(abs.(loops[i].coord)) * eps(Float64)
    end
    # Populate a dislocation network with the loops.
    # Test one branch of memory allocation.
    network = zero(DislocationNetwork)
    makeNetwork!(network, loops[1])
    @test network.numNode ==
          loops[1].numSides * loops[1].nodeSide * loops[1].numLoops
    # Test other branch of memory allocation.
    network = zero(DislocationNetwork)
    makeNetwork!(network, loops)
    function sumNodes(loops)
        totalNodes = 0
        for i in eachindex(loops)
            totalNodes +=
                loops[i].numSides * loops[i].nodeSide * loops[i].numLoops
        end
        return totalNodes
    end
    # Check that the memory was allocated correctly. Only need to check the first and last, they are transfered sequentially so if both pass, the rest have to have been transfered correctly.
    totalNodes = sumNodes(loops)
    @test totalNodes ==
          network.numNode ==
          network.numSeg ==
          size(network.links, 1) ==
          size(network.slipPlane, 1) ==
          size(network.bVec, 1) ==
          size(network.coord, 1) ==
          size(network.label, 1)
    # Check that the first loop was transfered correctly.
    nodeLoop = loops[1].numSides * loops[1].nodeSide * loops[1].numLoops
    @test network.links[1:nodeLoop, :] == loops[1].links
    @test network.slipPlane[1:nodeLoop, :] == loops[1].slipPlane
    @test network.bVec[1:nodeLoop, :] == loops[1].bVec
    @test network.coord[1:nodeLoop, :] == loops[1].coord
    @test network.label[1:nodeLoop] == loops[1].label
    # Check that the last loop was transfered correctly.
    nodeLoop = loops[end].numSides * loops[end].nodeSide * loops[end].numLoops
    @test network.links[(1 + end - nodeLoop):end, :] ==
          loops[end].links .+ (totalNodes - nodeLoop)
    @test network.slipPlane[(1 + end - nodeLoop):end, :] == loops[end].slipPlane
    @test network.bVec[(1 + end - nodeLoop):end, :] == loops[end].bVec
    @test network.coord[(1 + end - nodeLoop):end, :] == loops[end].coord
    @test network.label[(1 + end - nodeLoop):end] == loops[end].label
    network2 = DislocationNetwork(
        links = zeros(Int64, 1, 2),
        slipPlane = zeros(1, 3),
        bVec = zeros(1, 3),
        coord = zeros(1, 3),
        label = zeros(nodeType, 1),
        segForce = zeros(1, 3),
        numNode = convert(Int64, 0),
        numSeg = convert(Int64, 0),
        maxConnect = convert(Int64, 0),
    )
    makeNetwork!(network2, loops)
    @test compStruct(network, network2)
    network3 = makeNetwork(loops)
    @test network.links == network3.links[1:totalNodes, :]
    @test network.slipPlane == network3.slipPlane[1:totalNodes, :]
    @test network.bVec == network3.bVec[1:totalNodes, :]
    @test network.coord == network3.coord[1:totalNodes, :]
    @test network.label == network3.label[1:totalNodes]
    @test network.numNode == network3.numNode
    @test network.numSeg == network3.numSeg
    # Test distributions.
    n = 5
    seed!(1234)
    randArr = loopDistribution(Rand(), n)
    seed!(1234)
    test = rand(n, 3)
    @test randArr == test
    seed!(1234)
    randArr = loopDistribution(Randn(), n)
    seed!(1234)
    test = randn(n, 3)
    @test randArr == test
    @test_throws ErrorException loopDistribution(Regular(), n)
    @test checkNetwork(network)
    network.label[1] = 4
    @test !compStruct(network, network3)
end

@testset "Overloaded type functions" begin
    var = segEdge()
    @test length(var) == 1
    @test zero(nodeType) == -1
    loop = zero(DislocationLoop)
    @test length(loop) == 1
    @test loop[1] == loop
    @test_throws BoundsError loop[2]
    @test eachindex(loop[1]) == 1
end
