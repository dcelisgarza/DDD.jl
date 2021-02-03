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
    for i = 1:size(points, 1)
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
        0 1013 0
        1009 1013 0
        0 0 1019
        1009 0 1019
        0 1013 1019
        1009 1013 1019
    ]
    testC = [
        3.272727272727273 1.272727272727273 1.272727272727273 0 0 0
        1.272727272727273 3.272727272727273 1.272727272727273 0 0 0
        1.272727272727273 1.272727272727273 3.272727272727273 0 0 0
        0 0 0 1.000000000000000 0 0
        0 0 0 0 1.000000000000000 0
        0 0 0 0 0 1.000000000000000
    ]
    vertices = reshape(collect(Iterators.flatten(regularCuboidMesh.vertices.vertices)), 3, 8)
    faces = regularCuboidMesh.faces
    normals = regularCuboidMesh.faceNorm

    @test isapprox(vertices', testVertices)
    
    faceCoord = vertices[:, faces]
    p = faceCoord[:, 2, :] - faceCoord[:, 1, :]
    q = faceCoord[:, 3, :] - faceCoord[:, 1, :]
    n = reshape(collect(Iterators.flatten([normalize(p[:,i] × q[:,i]) for i in 1:size(p, 2)])), 3, 6)
    @test isapprox(n, normals)

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
    testCoord2 = 1.0e+02 * [
        7.338181818181819 0 0
        4.586363636363637 0 0.599411764705882
    ]
    testCoord3 = 1.0e+02 * [
        4.586363636363637 0 0
        3.669090909090909 0 0
    ]
    testCoord4 = 1.0e+02 * [
        3.669090909090909 0 0
        6.420909090909091 0 0.599411764705882
    ]
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
        1966 1977 2194 2181 1965
        2795 2806 3023 3010 2794
        682 693 910 897 681
        1998 2009 2226 2213 1997
        1125 1136 1353 1340 1124
    ]
    connectivity = regularCuboidMesh.connectivity

    idxCon = [6, 8, 3, 1, 5]
    idxNode = [1703, 2430, 592, 1732, 976]
    connectivity[idxCon, idxNode]' == testCon

    KTest = [
        29.516317016317018,
        -27.845980957491285,
        -73.225643849016691,
        -19.731828987678465,
        -0.284596060433909,
        32.172011818817367,
        1.042355371900828,
        -7.379079254079253,
        -5.676247771836008,
        -0.885489510489512,
    ]
    K = regularCuboidMesh.K
    droptol!(K, 1e-14)

    idxK = [
        CartesianIndex(3435, 3400)
        CartesianIndex(1108, 1069)
        CartesianIndex(8973, 9009)
        CartesianIndex(3019, 3664)
        CartesianIndex(6670, 6706)
        CartesianIndex(5488, 4840)
        CartesianIndex(675, 29)
        CartesianIndex(1042, 357)
        CartesianIndex(5710, 5096)
        CartesianIndex(1918, 1275)
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
    cornerNode = regularCuboidMesh.cornerNode
    edgeNode = regularCuboidMesh.edgeNode
    faceNode = regularCuboidMesh.faceNode
    coord = regularCuboidMesh.coord

    numNode = regularCuboidMesh.numNode
    numNode3 = numNode * 3
    uGamma = (type = nodeTypeFE(1),# Type
                idx = :x0y0z0, # Index
                node = cornerNode[:x0y0z0])
    tGamma = (type = nodeTypeFE(2),# Type
                idx = :x_y0z1, # Index
                node = edgeNode[:x_y0z1]) 
    mGamma = (type = nodeTypeFE(3),# Type
                idx = :xy_z0, # Index
                node = faceNode[:xy_z0])
    testGamma, missing = Boundaries(femParams, regularCuboidMesh; uGamma = uGamma, tGamma = tGamma, mGamma = mGamma)
    @test isequal(uGamma, testGamma.uGamma)
    @test isequal(tGamma, testGamma.tGamma)
    @test isequal(mGamma, testGamma.mGamma)

    cantileverBC, forceDisplacement = Boundaries(femParams, regularCuboidMesh)
    uGamma = cantileverBC.uGamma[:node]
    mGamma = cantileverBC.mGamma[:node]
    left = findall(x -> x == 0, coord[1, :])

    loadEdge1 = findall(x -> x ≈ dx, coord[1, :])
    loadEdge2 = findall(x -> x ≈ dz, coord[3, :])
    loadEdge = intersect(loadEdge1, loadEdge2)
    @test norm(coord[:, left]) ≈ norm(coord[:, uGamma])
    @test norm(coord[:, loadEdge]) ≈ norm(coord[:, mGamma])
end