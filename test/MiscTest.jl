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
    xyz = [1., 0., 0.]
    θ = pi/2
    uvw = [0.,5.,0.]
    abc = [0.,0.,0.]
    p = rot3D(xyz,uvw,abc,θ)
    @test isapprox(p,[0.,0.,-1.])
    uvw = [0.,0.,20.]
    abc = [0.,0.,0.]
    p = rot3D(xyz,uvw,abc,θ)
    @test isapprox(p,[0.,1.,0.])
    uvw = [1.,0.,0.]
    abc = [0.,0.,0.]
    p = rot3D(xyz,uvw,abc,θ)
    @test isapprox(p,xyz)
    xyz = [-23., 29., -31.]
    uvw = [11.,-13.,17.]
    abc = [-2.,5.,7.]
    θ = 37/180*pi
    p = rot3D(xyz,uvw,abc,θ)
    @test isapprox(p, [-21.1690, 31.0685, -30.6029]; atol = 1e-4)
end
