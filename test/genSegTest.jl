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

@testset "Generate single segments" begin
    @test isapprox(dot(edge, screw), 0.0)
    @test isapprox(dot(edge, slipSystem[4:6]), 0.0)
    @test isapprox(
        edge,
        cross(slipSystem[1:3], slipSystem[4:6]) ./
        norm(cross(slipSystem[1:3], slipSystem[4:6])),
    )
    @test isapprox(norm(edge), norm(screw))
    @test isapprox(norm(edge), 1.0)
end
