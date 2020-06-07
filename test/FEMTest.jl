using DDD
using Test, StaticArrays
cd(@__DIR__)
@testset "Shape functions" begin
    points = Float64[
        0 0 0
        1 0 0
        1 1 0
        0 1 0
        0 0 1
        1 0 1
        1 1 1
        0 1 1
        0.5 0 0
        0.5 0.5 0
        0 0.5 0
        0 0 0.5
        0.5 0 0.5
        0.5 0.5 0.5
        0 0.5 0.5
    ]
    Nall = shapeFunction(LinearQuadrangle3D(), points[:, 1], points[:, 2], points[:, 3])
    @test all(isapprox.(sum.(Nall), 1))

    dNdSall =
        shapeFunctionDeriv(LinearQuadrangle3D(), points[:, 1], points[:, 2], points[:, 3])
    checkSum = sum.(dNdSall, dims = 2)
    for i in 1:size(points, 1)
        @test vec(checkSum[i]) == SVector(0.0, 0.0, 0.0)

        N = shapeFunction(LinearQuadrangle3D(), points[i, 1], points[i, 2], points[i, 3])
        @test isapprox(Nall[i], N)

        dNdS = shapeFunctionDeriv(
            LinearQuadrangle3D(),
            points[i, 1],
            points[i, 2],
            points[i, 3],
        )

        @test isapprox(dNdSall[i], dNdS)
    end
end
