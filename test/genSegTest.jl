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
