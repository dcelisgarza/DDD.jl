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
    Nall = shapeFunction(
        LinearQuadrangle3D(),
        points[:, 1],
        points[:, 2],
        points[:, 3],
    )
    @test any(isapprox.(sum(Nall, dims = 1), 1))

    dNdSall = shapeFunctionDeriv(
        LinearQuadrangle3D(),
        points[:, 1],
        points[:, 2],
        points[:, 3],
    )

    for i = 1:size(points, 1)
        @test isapprox.(sum(dNdSall[:, :, i], dims = 1), 0) == Bool[1 1 1]

        N = shapeFunction(
            LinearQuadrangle3D(),
            points[i, 1],
            points[i, 2],
            points[i, 3],
        )
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

@testset "Compare loading input and output parameters" begin
    # Filenames
    inFileParams = "../inputs/simParams/sampleParams.csv"
    inFileSlipSys = "../data/slipSystems/bcc.csv"
    inFileDln = "../inputs/dln/samplePrismaticShear.csv"
    outFileParams = "../outputs/simParams/sampleParams.csv"
    # Load parameters.
    dlnParams, matParams, intParams, slipSystems, sources =
        loadParams(inFileParams, inFileSlipSys, inFileDln)
    saveParams(dlnParams, matParams, intParams, outFileParams; delim = ',')
    dlnParams2, matParams2, intParams2, slipSystems, sources =
        loadParams(outFileParams, inFileSlipSys, inFileDln)
    @test compStruct(dlnParams, dlnParams2)
    @test compStruct(matParams, matParams2)
    @test compStruct(intParams, intParams2)
end