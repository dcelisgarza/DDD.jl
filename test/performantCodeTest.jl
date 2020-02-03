using DDD
using Test, Traceur, DelimitedFiles

import LinearAlgebra: norm
cd(@__DIR__)
# Filenames
inFilename = "../data/slipSystems/bcc.csv"
slipSystems = readdlm(inFilename, ',', Float64)
slipSysInt = 5

@testset "Performant vector operations" begin
    @test isempty(Traceur.warnings(
        () -> makeSegment(dlnEdge(), slipSysInt, slipSystems),
        modules = [DDD],
    ))
    @test isempty(Traceur.warnings(
        () -> makeSegment(dlnScrew(), slipSysInt, slipSystems),
        modules = [DDD],
    ))
end


numNode = 5
numSeg = 7

links = zeros(Integer, numSeg, 2)
bVec = zeros(numSeg, 3)
slipPlane = zeros(numSeg, 3)

coord = zeros(numNode, 3)
label = zeros(Integer, numNode)

lenLinks = size(links, 1)
[links[i, :] = convert.(Float64, [i, i + lenLinks]) for i = 1:lenLinks]
[bVec[i, :] = convert.(Float64, [i, i + lenLinks, i + 2 * lenLinks]) for i = 1:lenLinks]
[slipPlane[i, :] = -convert.(Float64, [i, i + lenLinks, i + 2 * lenLinks]) for i = 1:lenLinks]

lenLabel = length(label)
[label[i] = i for i = 1:lenLabel]
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

@testset "Performant indexing operations" begin
    @test isempty(Traceur.warnings(() -> getIndex(network, 5), modules = [DDD]))
    @test isempty(Traceur.warnings(() -> getCoord(network, 2), modules = [DDD]))
    @test isempty(Traceur.warnings(
        () -> getCoord(network, [1; 3; 5]),
        modules = [DDD],
    ))
end

@testset "Optimised dynamic indexing operations" begin
    @test length(Traceur.warnings(
        () -> getIndex(network, :label, >=, 3),
        modules = [DDD],
    )) == 2
    @test length(Traceur.warnings(
        () -> getIndex(network, :coord, 2, >=, 8.0),
        modules = [DDD],
    )) == 3
    @test length(Traceur.warnings(
        () -> getData(network, :coord, :coord, 2, >=, 8.0),
        modules = [DDD],
    )) == 5
end
