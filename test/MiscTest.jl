using DDD
using Test

cd(@__DIR__)

@testset "Geometry" begin
    @test isapprox(intAngle(3), π / 3)
    @test isapprox(intAngle(4), π / 2)
    @test isapprox(intAngle(6), 2π / 3)
    @test isapprox(extAngle(3), 2π / 3)
    @test isapprox(extAngle(4), π / 2)
    @test isapprox(extAngle(6), π / 3)
end
