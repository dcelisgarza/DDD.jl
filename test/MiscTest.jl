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

@testset "Auxiliary" begin
    dict = Dict(
        "intFix" => nodeType(2),
        "none" => nodeType(0),
        "intMob" => nodeType(1),
        "srfFix" => nodeType(4),
        "ext" => nodeType(5),
        "srfMob" => nodeType(3),
    )

    @test makeInstanceDict(nodeType) == dict
    data = rand(5)
    @test inclusiveComparison(data[rand(1:5)], data...)
    @test !inclusiveComparison(data, data[rand(1:5)] * 6)

    arr = [1 2 3; 5 7 11]
    @test dimDot(arr, arr, dim = 1) == [26, 53, 130]
    @test dimDot(arr, arr, dim = 2) == [14, 195]
    @test dimNorm(arr, dim = 1) == [sqrt(1^2 + 5^2), sqrt(2^2 + 7^2), sqrt(3^2 + 11^2)]
    @test dimNorm(arr, dim = 2) == [sqrt(1^2+2^2+3^2), sqrt(5^2+7^2+11^2)]
end
