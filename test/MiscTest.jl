using DDD
using Test

cd(@__DIR__)

@testset "Geometry" begin
    arr = Int[3; 4; 6]
    @test isapprox(intAngle(arr[1]), π / 3)
    @test isapprox(intAngle(arr[2]), π / 2)
    @test isapprox(intAngle(arr[3]), 2π / 3)
    @test isapprox(externalAngle(arr[1]), 2π / 3)
    @test isapprox(externalAngle(arr[2]), π / 2)
    @test isapprox(externalAngle(arr[3]), π / 3)
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

    planenorm = Float64[0, 0, 1]
    planepnt  = Float64[0, 0, 5]
    raydir = Float64[0, -1, -2]
    raypnt = Float64[0,  0, 10]

    ψ = linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
    @test isapprox(ψ,  [0, -2.5, 5.0])

    planenorm = Float64[0, 2, 1]
    planepnt  = Float64[0, 0, 5]
    raydir = Float64[0, -1, -2]
    raypnt = Float64[0,  0, 10]
    ψ = linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
    @test isapprox(ψ, [0.0, -1.25, 7.5])

    planenorm = Float64[0, 0, 1]
    planepnt  = Float64[0, 0, 5]
    raydir = Float64[0, 1, 2]
    raypnt = Float64[0,  0, 10]

    ψ = linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
    @test isapprox(ψ,  [0, -2.5, 5.0])

    planenorm = Float64[0, 0, 1]
    planepnt  = Float64[0, 0, 5]
    raydir = Float64[0, 1, -2]
    raypnt = Float64[0,  0, 10]

    ψ = linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
    @test isapprox(ψ,  [0, 2.5, 5.0])

    planenorm = Float64[0, 0, 1]
    planepnt  = Float64[0, 0, 5]
    raydir = Float64[0, 1, 0]
    raypnt = Float64[0,  0, 5]
    ψ = linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
    @test isinf(ψ)

    planenorm = Float64[0, 0, 1]
    planepnt  = Float64[0, 0, 5]
    raydir = Float64[0, 1, 0]
    raypnt = Float64[0,  0, 6]
    ψ = linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
    @test isnothing(ψ)
end

@testset "Auxiliary" begin
    dict = Dict(
        "intFix" => nodeType(2),
        "none" => nodeType(0),
        "intMob" => nodeType(1),
        "srfFix" => nodeType(4),
        "virtual" => nodeType(5),
        "srfMob" => nodeType(3),
        "tmp" => nodeType(6)
    )

    @test makeInstanceDict(nodeType) == dict
    data = rand(5)
    @test inclusiveComparison(data[rand(1:5)], data...)
    @test !inclusiveComparison(data, data[rand(1:5)] * 6)
end
