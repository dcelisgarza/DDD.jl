using DDD
using Test

cd(@__DIR__)

@testset "Geometry" begin
    arr = Int[3; 4; 6]
    @test isapprox(intAngle(arr[1]), π / 3)
    @test isapprox(intAngle(arr[2]), π / 2)
    @test isapprox(intAngle(arr[3]), 2π / 3)
    @test isapprox(extAngle(arr[1]), 2π / 3)
    @test isapprox(extAngle(arr[2]), π / 2)
    @test isapprox(extAngle(arr[3]), π / 3)
    xyz = [1.0, 0.0, 0.0]
    θ = pi / 2
    uvw = [0.0, 5.0, 0.0]
    abc = [0.0, 0.0, 0.0]
    p = rot3D(xyz, uvw, abc, θ)
    @test isapprox(p, [0.0, 0.0, -1.0])
    uvw = [0.0, 0.0, 20.0]
    abc = [0.0, 0.0, 0.0]
    p = rot3D(xyz, uvw, abc, θ)
    @test isapprox(p, [0.0, 1.0, 0.0])
    uvw = [1.0, 0.0, 0.0]
    abc = [0.0, 0.0, 0.0]
    p = rot3D(xyz, uvw, abc, θ)
    @test isapprox(p, xyz)
    xyz = [-23.0, 29.0, -31.0]
    uvw = [11.0, -13.0, 17.0]
    abc = [-2.0, 5.0, 7.0]
    θ = 37 / 180 * pi
    p = rot3D(xyz, uvw, abc, θ)
    @test isapprox(p, [-21.1690, 31.0685, -30.6029]; atol = 1e-4)
    @test compStruct(1, 1.2) == false
end
