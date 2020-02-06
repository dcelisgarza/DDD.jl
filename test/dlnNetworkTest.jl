using DDD
using Test

import DelimitedFiles: readdlm
import LinearAlgebra: dot, cross, norm
cd(@__DIR__)

@testset "Generate single segments" begin
    inFilename = "../data/slipSystems/bcc.csv"
    data = readdlm(inFilename, ',', Float64)
    slipSysInt = 5
    slipSystem = data[:, slipSysInt]
    edge = makeSegment(dlnEdge(), slipSysInt, data)
    screw = makeSegment(dlnScrew(), slipSysInt, data)
    @test abs(dot(edge, screw)) < eps(Float64)
    @test abs(dot(edge, slipSystem[4:6])) < eps(Float64)
    @test isapprox(
        edge,
        cross(slipSystem[1:3], slipSystem[4:6]) ./
        norm(cross(slipSystem[1:3], slipSystem[4:6])),
    )
    @test isapprox(norm(edge), norm(screw))
    @test isapprox(norm(edge), 1.0)
end

@testset "Dislocation indexing functions" begin
    numNode = 10
    numSeg = 20
    links = zeros(Integer, numSeg, 2)
    bVec = zeros(Integer, numSeg, 3)
    slipPlane = zeros(Integer, numSeg, 3)
    coord = zeros(numNode, 3)
    label = zeros(nodeType, numNode)
    lenLinks = size(links, 1)
    [links[i, :] = [i, i + lenLinks] for i = 1:lenLinks]
    [bVec[i, :] = [i, i + lenLinks, i + 2 * lenLinks] for i = 1:lenLinks]
    [slipPlane[i, :] = -[i, i + lenLinks, i + 2 * lenLinks] for i = 1:lenLinks]
    lenLabel = length(label)
    [label[i+1] = mod(i, 6) - 1 for i = 0:lenLabel-1]
    [coord[i, :] = convert.(Float64, [i, i + lenLabel, i + 2 * lenLabel]) for i = 1:length(label)]
    network = DislocationNetwork(
        links,
        bVec,
        slipPlane,
        coord,
        label,
        numNode,
        numSeg,
    )
    rnd = rand(-1:numNode)
    @test idxLabel(network, rnd) == findall(x -> x == rnd, label)
    @test coordLbl(network, rnd) == coord[findall(x -> x == rnd, label), :]
    @test idxCond(network, :label, rnd; condition = <=) == findall(
        x -> x <= rnd,
        label,
    )
    rnd = rand(1:numSeg)
    cnd = rand([==, >=, <=, <, >, !==])
    @test idxCond(network, :bVec, rnd; condition = cnd) == findall(
        x -> cnd(x, rnd),
        bVec,
    )
    cnd = rand([==, >=, <=, <, >, !==])
    @test dataCond(
        network,
        :slipPlane,
        rnd;
        condition = cnd,
    ) == slipPlane[findall(x -> cnd(x, rnd), slipPlane)]
    col = rand(1:3)
    cnd = rand([==, >=, <=, <, >, !==])
    @test idxCond(network, :bVec, col, rnd; condition = cnd) == findall(
        x -> cnd(x, rnd),
        bVec[:, col],
    )
    cnd = rand([==, >=, <=, <, >, !==])
    @test dataCond(network, :bVec, col, rnd; condition = cnd) == bVec[
        findall(x -> cnd(x, rnd), bVec[:, col]),
        :,
    ]
    rnd = rand(1:numNode)
    @test coordIdx(network, rnd) == coord[rnd, :]
    rnd = rand(1:numNode, numNode)
    @test coordIdx(network, rnd) == coord[rnd, :]
end
