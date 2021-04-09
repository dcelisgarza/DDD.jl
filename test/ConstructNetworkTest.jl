using DDD
using Test

import LinearAlgebra: dot, cross, norm
import Statistics: mean
import Random: seed!

import DDD: segEdge, segEdgeN, segScrew, makeSegment, inclusiveComparison
cd(@__DIR__)

@testset "Generate single segments" begin
    testSlip = [-1.0, 0.0, 1.0]
    testBVec = [1.0, 1.0, -1.0]
    @test_throws AssertionError SlipSystem(;
        crystalStruct = BCC(),
        slipPlane = testSlip,
        bVec = testBVec,
    )
    testSlip = hcat(testSlip, [-1.0, 0.0, 1.0])
    testBVec = hcat(testBVec, [1.0, 1.0, 1.0])
    @test_throws AssertionError SlipSystem(;
        crystalStruct = BCC(),
        slipPlane = testSlip,
        bVec = testBVec,
    )
    fileSlipSystem = "./testData/BCC.json"
    dictSlipSystem = loadJSON(fileSlipSystem)
    slipSystems = loadSlipSystem(dictSlipSystem)
    slipSysInt = 1
    slipPlane = slipSystems.slipPlane[:, 1]
    bVec = slipSystems.bVec[:, 1]
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
    network = DislocationNetwork(;
        links = rand(1:10, 2, 10),
        slipPlane = rand(3, 10),
        bVec = rand(3, 10),
        coord = rand(3, 10),
        label = nodeTypeDln.(rand(0:5, 10)),
        segForce = rand(3, 2, 10),
        nodeVel = rand(3, 10),
        nodeForce = rand(3, 10),
        numNode = [10],
        numSeg = [10],
        maxConnect = 10,
    )
    getSegmentIdx!(network)
    makeConnect!(network)

    idx = rand(1:10)
    @test network[idx] == (
        network.links[:, idx],
        network.slipPlane[:, idx],
        network.bVec[:, idx],
        network.coord[:, idx],
        network.label[idx],
        network.nodeVel[:, idx],
        network.nodeForce[:, idx],
        network.connectivity[:, idx],
        network.linksConnect[:, idx],
        network.extSeg[idx],
        network.segIdx[idx, :],
        network.segForce[:, :, idx],
    )

    idx = rand(1:10, 5)
    @test network[idx] == (
        network.links[:, idx],
        network.slipPlane[:, idx],
        network.bVec[:, idx],
        network.coord[:, idx],
        network.label[idx],
        network.nodeVel[:, idx],
        network.nodeForce[:, idx],
        network.connectivity[:, idx],
        network.linksConnect[:, idx],
        network.extSeg[idx],
        network.segIdx[idx, :],
        network.segForce[:, :, idx],
    )
end

@testset "Loop generation" begin
    slipfile = "./testData/BCC.json"
    loopfile = "./testData/samplePrismShear.json"

    dictSlipSystem = loadJSON(slipfile)
    slipSystems = loadSlipSystem(dictSlipSystem)

    dictDislocationLoop = loadJSON(loopfile)

    if typeof(dictDislocationLoop) <: AbstractArray
        loops = zeros(DislocationLoop, length(dictDislocationLoop))
        for i in eachindex(loops)
            loops[i] = loadDislocationLoop(dictDislocationLoop[i], slipSystems)
        end
    else
        loops = loadDislocationLoop(dictDislocationLoop, slipSystems)
    end

    # Check that the midpoint of the loops is at (0,0,0)
    for i in eachindex(loops)
        @test mean(loops[i].coord) < maximum(abs.(loops[i].coord)) * eps(Float64)
    end
    # Populate a dislocation network with the loops.
    # Test one branch of memory allocation.
    network = zero(DislocationNetwork)
    network = DislocationNetwork!(network, loops[1])
    @test network.numNode[1] == loops[1].numSides * loops[1].nodeSide * loops[1].numLoops
    # Test other branch of memory allocation.
    network = zero(DislocationNetwork)
    network = DislocationNetwork!(network, loops)
    function sumNodes(loops)
        totalNodes = 0
        for i in eachindex(loops)
            totalNodes += loops[i].numSides * loops[i].nodeSide * loops[i].numLoops
        end
        return totalNodes
    end
    # Check that the memory was allocated correctly. Only need to check the first and last, they are transfered sequentially so if both pass, the rest have to have been transfered correctly.
    totalNodes = sumNodes(loops)
    @test Int(round(totalNodes * log2(totalNodes))) ==
          Int(round(network.numNode[1] * log2(network.numNode[1]))) ==
          Int(round(network.numSeg[1] * log2(network.numSeg[1]))) ==
          size(network.links, 2) ==
          size(network.slipPlane, 2) ==
          size(network.bVec, 2) ==
          size(network.coord, 2) ==
          size(network.label, 1)
    # Check that the first loop was transfered correctly.
    nodeLoop1 = loops[1].numSides * loops[1].nodeSide * loops[1].numLoops
    @test network.links[:, 1:size(loops[1].links, 2)] == loops[1].links
    @test network.slipPlane[:, 1:size(loops[1].links, 2)] == loops[1].slipPlane
    @test network.bVec[:, 1:size(loops[1].links, 2)] == loops[1].bVec
    @test network.coord[:, 1:size(loops[1].links, 2)] == loops[1].coord
    @test network.label[1:size(loops[1].links, 2)] == loops[1].label
    # Check that the last loop was transfered correctly.
    nodeLoop2 = loops[end].numSides * loops[end].nodeSide * loops[end].numLoops
    @test network.links[:, (totalNodes - size(loops[end].links, 2) + 1):totalNodes] ==
          loops[end].links .+ nodeLoop1 .+ nodeLoop2 .- size(loops[end].links, 2)
    @test network.slipPlane[:, (totalNodes - size(loops[end].links, 2) + 1):totalNodes] ==
          loops[end].slipPlane
    @test network.bVec[:, (totalNodes - size(loops[end].links, 2) + 1):totalNodes] ==
          loops[end].bVec
    @test network.coord[:, (totalNodes - size(loops[end].links, 2) + 1):totalNodes] ==
          loops[end].coord
    @test network.label[(totalNodes - size(loops[end].links, 2) + 1):totalNodes] ==
          loops[end].label
    network2 = DislocationNetwork(;
        links = zeros(Int, 2, 1),
        slipPlane = zeros(3, 1),
        bVec = zeros(3, 1),
        coord = zeros(3, 1),
        label = zeros(nodeTypeDln, 1),
        segForce = zeros(3, 2, 1),
        nodeVel = zeros(3, 1),
        nodeForce = zeros(3, 1),
        numNode = [0],
        numSeg = [0],
        maxConnect = 0,
    )
    network2 = DislocationNetwork!(network2, loops)
    @test compStruct(network, network2, verbose = true)
    network3 = DislocationNetwork(loops)
    @test network.links[:, 1:totalNodes] == network3.links[:, 1:totalNodes]
    @test network.slipPlane[:, 1:totalNodes] == network3.slipPlane[:, 1:totalNodes]
    @test network.bVec[:, 1:totalNodes] == network3.bVec[:, 1:totalNodes]
    @test network.coord[:, 1:totalNodes] == network3.coord[:, 1:totalNodes]
    @test network.label[1:totalNodes] == network3.label[1:totalNodes]
    @test network.numNode == network3.numNode
    @test network.numSeg == network3.numSeg
    # Test distributions.
    n = 5
    seed!(1234)
    randArr = loopDistribution(Rand(), n)
    seed!(1234)
    test = rand(3, n)
    @test randArr == test
    seed!(1234)
    randArr = loopDistribution(Randn(), n)
    seed!(1234)
    test = randn(3, n)
    @test randArr == test
    @test_throws ErrorException loopDistribution(Regular(), n)
    @test checkNetwork(network)
    network.label[1] = 4
    @test !compStruct(network, network3)

    import DDD: loopKink
    loopType = loopKink
    @test_logs (
        :warn,
        "DislocationLoop: Constructor for $loopType not defined, defaulting to prismatic loop.",
    ) DislocationLoop(;
        loopType = loopKink(),
        numSides = 4,
        nodeSide = 1,
        numLoops = 1,
        segLen = ones(4),
        slipSystemIdx = 4,
        slipSystem = slipSystems,
        label = nodeTypeDln[1; 1; 1; 1],
        buffer = 0.0,
        range = Float64[0 0; 0 0; 0 0],
        dist = Zeros(),
    )

    @test !iszero(network)
    networkZero = zero(DislocationNetwork)
    @test iszero(networkZero)
    networkZero = DislocationNetwork!(networkZero, loops[1]; maxConnect = 0)
    @test !iszero(networkZero)

    network = DislocationNetwork(loops[1], memBuffer = 1)
    numNode = network.numNode[1]
    numSeg = network.numSeg[1]
    network = DislocationNetwork!(network, loops[1])
    network = DislocationNetwork!(network, loops[1])
    @test network.numNode[1] == 3 * numNode
    @test network.numSeg[1] == 3 * numSeg
end

@testset "Overloaded type functions" begin
    @test nodeTypeDln(1) == 1
    @test 1 == nodeTypeDln(1)
    @test convert(nodeTypeDln, 2) == nodeTypeDln(2)
    @test zero(nodeTypeDln) == nodeTypeDln(0)
end
