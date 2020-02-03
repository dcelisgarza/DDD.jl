using DDD
using Test

import DelimitedFiles: readdlm
import LinearAlgebra: dot, cross, norm
cd(@__DIR__)
# Filenames
inFilename = "../data/slipSystems/bcc.csv"
data = readdlm(inFilename,',',Float64)

slipSysInt = 5
slipSystem = data[:,slipSysInt]
edge = makeSegment(dlnEdge(), slipSysInt, data)
screw = makeSegment(dlnScrew(), slipSysInt, data)

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
