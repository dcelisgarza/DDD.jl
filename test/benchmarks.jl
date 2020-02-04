using DDD
using Test, Traceur, DelimitedFiles, BenchmarkTools
import LinearAlgebra: norm

cd(@__DIR__)

# Segment.
inFilename = "../data/slipSystems/bcc.csv"
slipSystems = readdlm(inFilename, ',', Float64)
slipSysInt = 5
makeSegment(dlnEdge(), slipSysInt, slipSystems)
makeSegment(dlnScrew(), slipSysInt, slipSystems)

@benchmark makeSegment(dlnEdge(), slipSysInt, slipSystems)
@benchmark makeSegment(dlnScrew(), slipSysInt, slipSystems)
@trace(makeSegment(dlnEdge(), slipSysInt, slipSystems), modules = [DDD])
@trace(makeSegment(dlnEdge(), slipSysInt, slipSystems), modules = [DDD])

# Dislocation
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
[label[i] = -i for i = 1:lenLabel]
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
@benchmark DislocationNetwork(
    links,
    bVec,
    slipPlane,
    coord,
    label,
    numNode,
    numSeg,
)
@trace(
    DislocationNetwork(links, bVec, slipPlane, coord, label, numNode, numSeg),
    modules = [DDD],
)

# Indexing
idxLabel(network, 5)
idxCond(network, :label, -3; condition = >=)
idx = idxCond(network, :bVec, 10; condition = >=)
lidx = LinearIndices(bVec)[idx]
LinearIndices(idx)
lidx[1:7:end]
LinearIndices(bVec)[lidx]
dataCond(network, :bVec, 10; condition = >=)
idxCond(network, :bVec, 2, 14.; condition = ==)
dataCond(network, :bVec, 2, 14.; condition = ==)


coordLbl(network, -2)
coordIdx(network, 3)
coordIdx(network, [1; 3; 5])







dataCond(network, :links, :bVec, 2, 1; condition = ==)

@benchmark idxLabel(network, 5)
@benchmark idxCond(network, :label, -3; condition = >=)
@benchmark idxCond(network, :bVec, 3; condition = >=)
@benchmark idxCond(network, :bVec, 1, 5.0; condition = <)
@benchmark coordLbl(network, -2)
@benchmark coordIdx(network, 3)
@benchmark coordIdx(network, [1; 3; 5])
@benchmark dataCond(network, :coord, 2, 8.0; condition = >=)
@benchmark dataCond(network, :coord, :coord, 2, 8.0; condition = >=)
@benchmark dataCond(network, :coord, :label, 1, 8.0; condition = >=)
@benchmark dataCond(network, :links, :bVec, 2, 1; condition = ==)

@benchmark getIndex(network, 5)
@benchmark getIndex(network, :label, >=, 3)
@benchmark getIndex(network, :coord, 2, >=, 8.0)
@benchmark getCoord(network, 2)
@benchmark getCoord(network, [1; 3; 5])
@benchmark getData(network, :coord, :coord, 2, >=, 8.0)
@trace(getIndex(network, 5), modules = [DDD])
@trace(coordIdx(network, 2), modules = [DDD])
@trace(coordIdx(network, [1; 3; 5]), modules = [DDD])

# @testset "Performant indexing operations" begin
#     @test isempty(Traceur.warnings(() -> getIndex(network, 5), modules = [DDD]))
#     @test isempty(Traceur.warnings(() -> getCoord(network, 2), modules = [DDD]))
#     @test isempty(Traceur.warnings(
#         () -> getCoord(network, [1; 3; 5]),
#         modules = [DDD],
#     ))
# end

# @testset "Optimised dynamic indexing operations" begin
#     @test length(Traceur.warnings(
#         () -> getIndex(network, :label, >=, 3),
#         modules = [DDD],
#     )) == 2
#     @test length(Traceur.warnings(
#         () -> getIndex(network, :coord, 2, >=, 8.0),
#         modules = [DDD],
#     )) == 3
#     @test length(Traceur.warnings(
#         () -> getData(network, :coord, :coord, 2, >=, 8.0),
#         modules = [DDD],
#     )) == 5
# end
