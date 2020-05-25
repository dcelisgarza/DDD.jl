using DDD
using Test
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
    @test any(isapprox.(sum(Nall, dims = 1), 1))

    dNdSall =
        shapeFunctionDeriv(LinearQuadrangle3D(), points[:, 1], points[:, 2], points[:, 3])

    for i in 1:size(points, 1)
        @test isapprox.(sum(dNdSall[:, :, i], dims = 1), 0) == Bool[1 1 1]

        N = shapeFunction(LinearQuadrangle3D(), points[i, 1], points[i, 2], points[i, 3])
        @test isapprox(Nall[:, i], N)

        dNdS = shapeFunctionDeriv(
            LinearQuadrangle3D(),
            points[i, 1],
            points[i, 2],
            points[i, 3],
        )

        @test isapprox(dNdSall[:, :, i], dNdS)
    end
end
