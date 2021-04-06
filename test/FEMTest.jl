using DDD
using Test, StaticArrays, SparseArrays, LinearAlgebra
cd(@__DIR__)

@testset "nodeFE type" begin
    @test nodeTypeFE(0) == 0
    @test 1 == nodeTypeFE(1)
    @test nodeTypeFE(3) == 3.0
    @test zero(nodeTypeFE) == 0
    @test convert(nodeTypeFE, 3.0) == nodeTypeFE(3)
end

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

@testset "Checking arbitrary points of the mesh and connectivity" begin
    DictMaterialParameters = loadJSON("./testData/sampleMaterialParameters.json")
    matParams = loadMaterialParameters(DictMaterialParameters)
    femParams = FEMParameters(;
        type = DispatchRegularCuboidMesh(),
        order = LinearElement(),
        model = CantileverLoad(),
        dx = Float64(1009),
        dy = Float64(1013),
        dz = Float64(1019),
        mx = 11,
        my = 13,
        mz = 17,
    )
    regularCuboidMesh = buildMesh(matParams, femParams)

    @test typeof(regularCuboidMesh) <: AbstractRegularCuboidMesh &&
          typeof(femParams.type) <: AbstractRegularCuboidMesh
    @test regularCuboidMesh.order == femParams.order

    testVertices = [
        0 0 0
        1009 0 0
        1009 1013 0
        0 1013 0
        0 0 1019
        1009 0 1019
        1009 1013 1019
        0 1013 1019
    ]
    testC = [
        3.272727272727273 1.272727272727273 1.272727272727273 0 0 0
        1.272727272727273 3.272727272727273 1.272727272727273 0 0 0
        1.272727272727273 1.272727272727273 3.272727272727273 0 0 0
        0 0 0 1.000000000000000 0 0
        0 0 0 0 1.000000000000000 0
        0 0 0 0 0 1.000000000000000
    ]
    vertices =
        reshape(collect(Iterators.flatten(regularCuboidMesh.vertices.vertices)), 3, 8)
    faces = regularCuboidMesh.faces
    normals = regularCuboidMesh.faceNorm

    @test isapprox(vertices', testVertices)

    faceCoord = vertices[:, faces]
    p = faceCoord[:, 2, :] - faceCoord[:, 1, :]
    q = faceCoord[:, 4, :] - faceCoord[:, 1, :]
    n = reshape(
        collect(Iterators.flatten([normalize(p[:, i] × q[:, i]) for i in 1:size(p, 2)])),
        3,
        6,
    )
    @test isapprox(n, normals)

    mx = regularCuboidMesh.mx
    my = regularCuboidMesh.my
    mz = regularCuboidMesh.mz
    faceNorm = regularCuboidMesh.faceNorm
    lbl = regularCuboidMesh.surfElemNode
    coord = regularCuboidMesh.coord
    for (i, val) in enumerate([
        1,
        1 + mx * mz,
        1 + mx * mz + my * mz,
        1 + 2 * mx * mz + my * mz,
        1 + 2 * mx * mz + 2 * my * mz,
        1 + 2 * mx * mz + 2 * my * mz + mx * my,
    ])
        s = coord[:, lbl[:, val]]
        p = SVector{3, eltype(s)}(s[:, 2] - s[:, 1])
        q = SVector{3, eltype(s)}(s[:, 4] - s[:, 1])
        @test normalize(p × q) ≈ faceNorm[:, i]
    end

    @test isapprox(regularCuboidMesh.C, testC)

    @test regularCuboidMesh.dx == femParams.dx
    @test regularCuboidMesh.dy == femParams.dy
    @test regularCuboidMesh.dz == femParams.dz
    @test regularCuboidMesh.mx == femParams.mx
    @test regularCuboidMesh.my == femParams.my
    @test regularCuboidMesh.mz == femParams.mz
    @test regularCuboidMesh.numElem ==
          regularCuboidMesh.mx * regularCuboidMesh.my * regularCuboidMesh.mz
    @test regularCuboidMesh.numNode ==
          (regularCuboidMesh.mx + 1) *
          (regularCuboidMesh.my + 1) *
          (regularCuboidMesh.mz + 1)

    @test regularCuboidMesh.dx / regularCuboidMesh.mx == regularCuboidMesh.w
    @test regularCuboidMesh.dy / regularCuboidMesh.my == regularCuboidMesh.h
    @test regularCuboidMesh.dz / regularCuboidMesh.mz == regularCuboidMesh.d

    idx = [
        9 3
        9 18
        6 5
        5 20
        10 9
    ]

    testCoord1 = 1.0e+02 * [
        7.338181818181819 0 0
        1.834545454545455 0 0
    ]
    testCoord2 = [733.8181818181819 0.0 0.0; 458.6363636363637 77.92307692307692 0.0]
    testCoord3 = 1.0e+02 * [
        4.586363636363637 0 0
        3.669090909090909 0 0
    ]
    testCoord4 = [366.90909090909093 0.0 0.0; 642.0909090909091 77.92307692307692 0.0]
    testCoord5 = 1.0e+02 * [
        8.255454545454546 0 0
        7.338181818181819 0 0
    ]
    coord = regularCuboidMesh.coord
    @test isapprox(coord[:, idx[1, :]]', testCoord1)
    @test isapprox(coord[:, idx[2, :]]', testCoord2)
    @test isapprox(coord[:, idx[3, :]]', testCoord3)
    @test isapprox(coord[:, idx[4, :]]', testCoord4)
    @test isapprox(coord[:, idx[5, :]]', testCoord5)

    testCon = [
        2508 2519 2352 2339 2507
        1518 1529 1362 1349 1517
        2750 2761 2594 2581 2749
        2750 2761 2594 2581 2749
        2703 2714 2547 2534 2702
    ]
    connectivity = regularCuboidMesh.connectivity

    idxCon = [6, 8, 3, 1, 5]
    idxNode = [2002, 1149, 2201, 2201, 2158]
    @test connectivity[idxCon, idxNode]' == testCon

    KTest = [
        -27.845980957491278
        -45.099576129125296
        -64.81800121385578
        -7.379079254079251
        3.5419580419580443
        29.51631701631701
        5.676247771836005
        -9.901489001393466
        -22.704991087344027
        7.379079254079253
    ]
    K = regularCuboidMesh.K
    droptol!(K, 1e-14)

    idxK = [
        CartesianIndex(2194, 2695)
        CartesianIndex(8235, 7734)
        CartesianIndex(5636, 5672)
        CartesianIndex(7785, 8326)
        CartesianIndex(7018, 7524)
        CartesianIndex(3367, 2868)
        CartesianIndex(6289, 5753)
        CartesianIndex(8368, 7903)
        CartesianIndex(724, 764)
        CartesianIndex(5121, 4654)
    ]
    @test isapprox(K[idxK], KTest)
end

@testset "Boundary conditions" begin
    DictMaterialParameters = loadJSON("./testData/sampleMaterialParameters.json")
    matParams = loadMaterialParameters(DictMaterialParameters)
    femParams = FEMParameters(;
        type = DispatchRegularCuboidMesh(),
        order = LinearElement(),
        model = CantileverLoad(),
        dx = Float64(1009),
        dy = Float64(1013),
        dz = Float64(1019),
        mx = 11,
        my = 13,
        mz = 17,
    )

    regularCuboidMesh = buildMesh(matParams, femParams)
    dx = regularCuboidMesh.dx
    dz = regularCuboidMesh.dz
    faces = regularCuboidMesh.faces
    coord = regularCuboidMesh.coord
    connectivity = regularCuboidMesh.connectivity
    surfNode = regularCuboidMesh.surfNode
    coord = regularCuboidMesh.coord

    cantileverBC, forceDisplacement = Boundaries(femParams, regularCuboidMesh)
    uGamma = cantileverBC.uGamma.node
    mGamma = cantileverBC.mGamma.node
    left = findall(x -> x == 0, coord[1, :])

    loadEdge1 = findall(x -> x ≈ dx, coord[1, :])
    loadEdge2 = findall(x -> x ≈ dz, coord[3, :])
    loadEdge = intersect(loadEdge1, loadEdge2)
    @test norm(coord[:, left]) ≈ norm(coord[:, uGamma])
    @test norm(coord[:, loadEdge]) ≈ norm(coord[:, mGamma])
end
