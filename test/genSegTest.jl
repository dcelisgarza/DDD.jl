using DDD
using Test

import LinearAlgebra: dot, cross, norm
cd(@__DIR__)
# Filenames
inFilename = "../data/slipSystems/bcc.csv"
df = loadCSV(inFilename; header = 1)

slipSysInt = 5
slipSystem = convert.(Float64, Vector(df[slipSysInt, :]))
edge = makeSegment!(dlnEdge(), slipSysInt, df)
screw = makeSegment!(dlnScrew(), slipSysInt, df)

edgeDscrew = dot(edge, screw)
edgeDbVec = dot(edge, slipSystem[4:6])
edgeTest = cross(slipSystem[1:3], slipSystem[4:6])
edgeTest ./= norm(edgeTest)
normEdge = norm(edge)
normScrew = norm(screw)

@testset "Generate single segments" begin
    @test isapprox(edgeDscrew, 0.0)
    @test isapprox(edgeDbVec, 0.0)
    @test isapprox(edge, edgeTest)
    @test isapprox(normEdge, normScrew)
    @test isapprox(normEdge, 1.0)
end
