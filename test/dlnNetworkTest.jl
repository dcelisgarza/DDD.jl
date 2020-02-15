using DDD
using Test

import DelimitedFiles: readdlm
import LinearAlgebra: dot, cross, norm
import Statistics: mean
cd(@__DIR__)

@testset "Generate single segments" begin
    inFilename = "../data/slipSystems/bcc.csv"
    data = readdlm(inFilename, ',', Float64)
    slipSysInt = 1
    slipSystem = data[slipSysInt, :]
    edge = makeSegment(segEdge(), slipSysInt, data)
    screw = makeSegment(segScrew(), slipSysInt, data)
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
    cnd = [==, >=, <=, <, >, !=]
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
        slipPlane,
        bVec,
        coord,
        label,
        numNode,
        numSeg,
    )
    @test isequal(network.label[1], -1)
    @test isequal(-1, network.label[1])
    @test -1 == network.label[1]
    @test -2.0 < network.label[1]
    rnd = rand(-1:numNode)
    @test idxLabel(network, rnd) == findall(x -> x == rnd, label)
    @test coordLbl(network, rnd) == coord[findall(x -> x == rnd, label), :]
    @test idxCond(network, :label, rnd; condition = <=) == findall(
        x -> x <= rnd,
        label,
    )
    rnd = rand(1:numSeg)
    @test idxCond(network, :bVec, rnd; condition = cnd[1]) == findall(
        x -> cnd[1](x, rnd),
        bVec,
    )
    @test dataCond(
        network,
        :slipPlane,
        rnd;
        condition = cnd[2],
    ) == slipPlane[findall(x -> cnd[2](x, rnd), slipPlane)]
    rnd = rand(-1:1)
    @test dataCond(
        network,
        :slipPlane,
        :bVec,
        rnd;
        condition = cnd[3],
    ) == slipPlane[findall(x -> cnd[3](x, rnd), bVec)]
    rnd = rand(1:numNode)
    @test dataCond(network, :coord, :label, rnd; condition = cnd[4]) == coord[
        findall(x -> cnd[4](x, rnd), label),
        :,
    ]
    col = rand(1:3)
    @test idxCond(network, :bVec, col, rnd; condition = cnd[5]) == findall(
        x -> cnd[5](x, rnd),
        bVec[:, col],
    )
    @test dataCond(network, :bVec, col, rnd; condition = cnd[6]) == bVec[
        findall(x -> cnd[6](x, rnd), bVec[:, col]),
        :,
    ]
    rnd = rand(-1:1)
    @test dataCond(
        network,
        :slipPlane,
        :bVec,
        col,
        rnd;
        condition = cnd[1],
    ) == slipPlane[findall(x -> cnd[1](x, rnd), bVec[:, col]), :]
    rnd = rand(1:numNode)
    @test coordIdx(network, rnd) == coord[rnd, :]
    rnd = rand(1:numNode, numNode)
    @test coordIdx(network, rnd) == coord[rnd, :]
end

@testset "Loop generation" begin
    slipfile = "../data/slipSystems/bcc.csv"
    loopfile = "../inputs/dln/sampleDln.csv"
    slipSystems = readdlm(slipfile, ',')
    df = loadCSV(loopfile; header = 1, transpose = true)
    loops = loadDln(df, slipSystems)
    for i in eachindex(loops)
        @test mean(loops[i].coord) < maximum(abs.(loops[i].coord)) *
                                     eps(Float64)
    end
end
