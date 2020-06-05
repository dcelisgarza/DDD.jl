using DDD
using Test, Statistics
cd(@__DIR__)

@testset "Merge nodes" begin
    fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
    fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
    fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
    fileSlipSystem = "../data/slipSystems/SlipSystems.JSON"
    fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
    dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
        fileDislocationP,
        fileMaterialP,
        fileIntegrationP,
        fileSlipSystem,
        fileDislocationLoop,
    )

    pentagon = DislocationLoop(
        loopPrism();
        numSides = 5,
        nodeSide = 1,
        numLoops = 2,
        segLen = 10 * ones(5),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork(pentagon; memBuffer = 1)
    network.coord[:, 6:end] .+= [10; 10; 10]
    factor = rand()

    for j in 1:size(network.segForce, 1)
        for i in 1:size(network.segForce, 3)
            network.segForce[j, 1, i] = i + (j - 1) * size(network.segForce, 3)
            network.segForce[j, 2, i] = network.segForce[j, 1, i] * factor
            network.nodeVel[j, i] = -i - (j - 1) * size(network.segForce, 3)
        end
    end

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 1, 1)
    @test compStruct(networkTest, network)

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 1, 2)
    @test !compStruct(networkTest, network)
    links = [
        2 6
        1 3
        3 4
        4 5
        5 1
        6 7
        7 8
        8 9
        9 2
        0 0
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0 0 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0 0 0
    ]
    coord =
        1.0e+02 * [
            -0.971749965049807 -1.069074863996580 -1.040824829046386
            -0.830925136003420 -0.928250034950193 -0.859175170953614
            -1.060150095500755 -0.939849904499245 -1.000000000000000
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.871749965049807 -0.969074863996580 -0.940824829046386
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            -0.960150095500755 -0.839849904499245 -0.900000000000000
            -0.885559421232031 -0.848384617786012 -0.833944039018043
            0 0 0
        ]
    label = [
        1
        1
        1
        2
        1
        1
        2
        1
        2
        0
    ]
    segForce = [
        10 20 30
        2 12 22
        3 13 23
        4 14 24
        5 15 25
        6 16 26
        7 17 27
        8 18 28
        9 19 29
        0 0 0
    ]
    nodeVel = [
        -1 -11 -21
        -10 -20 -30
        -3 -13 -23
        -4 -14 -24
        -5 -15 -25
        -6 -16 -26
        -7 -17 -27
        -8 -18 -28
        -9 -19 -29
        0 0 0
    ]
    connectivity = [
        2 2 1 5 2 0 0 0 0
        2 9 2 1 1 0 0 0 0
        2 2 2 3 1 0 0 0 0
        2 3 2 4 1 0 0 0 0
        2 4 2 5 1 0 0 0 0
        2 6 1 1 2 0 0 0 0
        2 6 2 7 1 0 0 0 0
        2 7 2 8 1 0 0 0 0
        2 8 2 9 1 0 0 0 0
        0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        2 2
        1 1
        2 1
        2 1
        2 2
        1 1
        2 1
        2 1
        2 1
        0 0
    ]

    @test isapprox(networkTest.links, links')
    @test isapprox(networkTest.slipPlane, slipPlane')
    @test isapprox(networkTest.bVec, bVec')
    @test isapprox(networkTest.coord, coord')
    @test isapprox(Int.(networkTest.label), label)
    @test isapprox(networkTest.segForce[:, 1, :], segForce')
    @test isapprox(networkTest.segForce[:, 2, :], factor * segForce')
    @test isapprox(networkTest.nodeVel, nodeVel')
    @test networkTest.numNode == 9
    @test networkTest.numSeg == 9
    @test isapprox(networkTest.connectivity, connectivity')
    @test isapprox(networkTest.linksConnect, linksConnect')

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 1, 3)
    @test !compStruct(networkTest, network)
    links = [
        2 3
        3 6
        1 4
        4 5
        5 1
        6 7
        7 8
        8 2
        0 0
        0 0
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0 0 0
        0 0 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0 0 0
        0 0 0
    ]
    coord =
        1.0e+02 * [
            -0.971749965049807 -1.069074863996580 -1.040824829046386
            -0.885559421232031 -0.848384617786012 -0.833944039018043
            -0.830925136003420 -0.928250034950193 -0.859175170953614
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.871749965049807 -0.969074863996580 -0.940824829046386
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            -0.960150095500755 -0.839849904499245 -0.900000000000000
            0 0 0
            0 0 0
        ]
    label = [
        1
        2
        1
        2
        1
        1
        2
        1
        0
        0
    ]
    segForce = [
        9 19 29
        10 20 30
        3 13 23
        4 14 24
        5 15 25
        6 16 26
        7 17 27
        8 18 28
        0 0 0
        0 0 0
    ]
    nodeVel = [
        -1 -11 -21
        -9 -19 -29
        -10 -20 -30
        -4 -14 -24
        -5 -15 -25
        -6 -16 -26
        -7 -17 -27
        -8 -18 -28
        0 0 0
        0 0 0
    ]
    connectivity = [
        2 3 1 5 2 0 0 0 0
        2 8 2 1 1 0 0 0 0
        2 1 2 2 1 0 0 0 0
        2 3 2 4 1 0 0 0 0
        2 4 2 5 1 0 0 0 0
        2 6 1 2 2 0 0 0 0
        2 6 2 7 1 0 0 0 0
        2 7 2 8 1 0 0 0 0
        0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        2 1
        2 2
        1 1
        2 1
        2 2
        1 1
        2 1
        2 1
        0 0
        0 0
    ]

    @test isapprox(networkTest.links, links')
    @test isapprox(networkTest.slipPlane, slipPlane')
    @test isapprox(networkTest.bVec, bVec')
    @test isapprox(networkTest.coord, coord')
    @test isapprox(Int.(networkTest.label), label)
    @test isapprox(networkTest.segForce[:, 1, :], segForce')
    @test isapprox(networkTest.segForce[:, 2, :], factor * segForce')
    @test isapprox(networkTest.nodeVel, nodeVel')
    @test networkTest.numNode == 8
    @test networkTest.numSeg == 8
    @test isapprox(networkTest.connectivity, connectivity')
    @test isapprox(networkTest.linksConnect, linksConnect')

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 10, 7)
    @test !compStruct(networkTest, network)
    links = [
        1 2
        2 3
        3 4
        4 5
        5 1
        6 7
        7 8
        8 6
        0 0
        0 0
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0 0 0
        0 0 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0 0 0
        0 0 0
    ]
    coord =
        1.0e+02 * [
            -0.971749965049807 -1.069074863996580 -1.040824829046386
            -1.051615382213988 -1.014440578767969 -1.066055960981957
            -1.060150095500755 -0.939849904499245 -1.000000000000000
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.885559421232031 -0.848384617786012 -0.833944039018043
            -0.830925136003420 -0.928250034950193 -0.859175170953614
            -0.960150095500755 -0.839849904499245 -0.900000000000000
            0 0 0
            0 0 0
        ]
    label = [
        1
        2
        1
        2
        1
        2
        1
        1
        0
        0
    ]
    segForce = [
        1 11 21
        2 12 22
        3 13 23
        4 14 24
        5 15 25
        9 19 29
        7 17 27
        8 18 28
        0 0 0
        0 0 0
    ]
    nodeVel = [
        -1 -11 -21
        -2 -12 -22
        -3 -13 -23
        -4 -14 -24
        -5 -15 -25
        -9 -19 -29
        -10 -20 -30
        -8 -18 -28
        0 0 0
        0 0 0
    ]
    connectivity = [
        2 1 1 5 2 0 0 0 0
        2 1 2 2 1 0 0 0 0
        2 2 2 3 1 0 0 0 0
        2 3 2 4 1 0 0 0 0
        2 4 2 5 1 0 0 0 0
        2 8 2 6 1 0 0 0 0
        2 6 2 7 1 0 0 0 0
        2 7 2 8 1 0 0 0 0
        0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        1 1
        2 1
        2 1
        2 1
        2 2
        2 1
        2 1
        2 1
        0 0
        0 0
    ]
    @test isapprox(networkTest.links, links')
    @test isapprox(networkTest.slipPlane, slipPlane')
    @test isapprox(networkTest.bVec, bVec')
    @test isapprox(networkTest.coord, coord')
    @test isapprox(Int.(networkTest.label), label)
    @test isapprox(networkTest.segForce[:, 1, :], segForce')
    @test isapprox(networkTest.segForce[:, 2, :], factor * segForce')
    @test isapprox(networkTest.nodeVel, nodeVel')
    @test networkTest.numNode == 8
    @test networkTest.numSeg == 8
    @test isapprox(networkTest.connectivity, connectivity')
    @test isapprox(networkTest.linksConnect, linksConnect')

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 7, 10)
    @test !compStruct(networkTest, network)
    links = [
        1 2
        2 3
        3 4
        4 5
        5 1
        6 7
        7 8
        8 6
        0 0
        0 0
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0 0 0
        0 0 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0 0 0
        0 0 0
    ]
    coord =
        1.0e+02 * [
            -0.971749965049807 -1.069074863996580 -1.040824829046386
            -1.051615382213988 -1.014440578767969 -1.066055960981957
            -1.060150095500755 -0.939849904499245 -1.000000000000000
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.885559421232031 -0.848384617786012 -0.833944039018043
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            -0.960150095500755 -0.839849904499245 -0.900000000000000
            0 0 0
            0 0 0
        ]
    label = [
        1
        2
        1
        2
        1
        2
        2
        1
        0
        0
    ]
    segForce = [
        1 11 21
        2 12 22
        3 13 23
        4 14 24
        5 15 25
        9 19 29
        7 17 27
        8 18 28
        0 0 0
        0 0 0
    ]
    nodeVel = [
        -1 -11 -21
        -2 -12 -22
        -3 -13 -23
        -4 -14 -24
        -5 -15 -25
        -9 -19 -29
        -7 -17 -27
        -8 -18 -28
        0 0 0
        0 0 0
    ]
    connectivity = [
        2 1 1 5 2 0 0 0 0
        2 1 2 2 1 0 0 0 0
        2 2 2 3 1 0 0 0 0
        2 3 2 4 1 0 0 0 0
        2 4 2 5 1 0 0 0 0
        2 8 2 6 1 0 0 0 0
        2 6 2 7 1 0 0 0 0
        2 7 2 8 1 0 0 0 0
        0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        1 1
        2 1
        2 1
        2 1
        2 2
        2 1
        2 1
        2 1
        0 0
        0 0
    ]
    @test isapprox(networkTest.links, links')
    @test isapprox(networkTest.slipPlane, slipPlane')
    @test isapprox(networkTest.bVec, bVec')
    @test isapprox(networkTest.coord, coord')
    @test isapprox(Int.(networkTest.label), label)
    @test isapprox(networkTest.segForce[:, 1, :], segForce')
    @test isapprox(networkTest.segForce[:, 2, :], factor * segForce')
    @test isapprox(networkTest.nodeVel, nodeVel')
    @test networkTest.numNode == 8
    @test networkTest.numSeg == 8
    @test isapprox(networkTest.connectivity, connectivity')
    @test isapprox(networkTest.linksConnect, linksConnect')

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 1, 10)
    @test !compStruct(networkTest, network)
    links = [
        1 2
        2 3
        3 4
        4 5
        5 1
        6 7
        7 8
        8 9
        9 1
        1 6
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
    ]
    coord =
        1.0e+02 * [
            -0.971749965049807 -1.069074863996580 -1.040824829046386
            -1.051615382213988 -1.014440578767969 -1.066055960981957
            -1.060150095500755 -0.939849904499245 -1.000000000000000
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.871749965049807 -0.969074863996580 -0.940824829046386
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            -0.960150095500755 -0.839849904499245 -0.900000000000000
            -0.885559421232031 -0.848384617786012 -0.833944039018043
            0 0 0
        ]
    label = [
        1
        2
        1
        2
        1
        1
        2
        1
        2
        0
    ]
    segForce = [
        1 11 21
        2 12 22
        3 13 23
        4 14 24
        5 15 25
        6 16 26
        7 17 27
        8 18 28
        9 19 29
        10 20 30
    ]
    nodeVel = [
        -1 -11 -21
        -2 -12 -22
        -3 -13 -23
        -4 -14 -24
        -5 -15 -25
        -6 -16 -26
        -7 -17 -27
        -8 -18 -28
        -9 -19 -29
        0 0 0
    ]
    connectivity = [
        4 1 1 5 2 9 2 10 1
        2 1 2 2 1 0 0 0 0
        2 2 2 3 1 0 0 0 0
        2 3 2 4 1 0 0 0 0
        2 4 2 5 1 0 0 0 0
        2 6 1 10 2 0 0 0 0
        2 6 2 7 1 0 0 0 0
        2 7 2 8 1 0 0 0 0
        2 8 2 9 1 0 0 0 0
        0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        1 1
        2 1
        2 1
        2 1
        2 2
        1 1
        2 1
        2 1
        2 3
        4 2
    ]
    @test isapprox(networkTest.links, links')
    @test isapprox(networkTest.slipPlane, slipPlane')
    @test isapprox(networkTest.bVec, bVec')
    @test isapprox(networkTest.coord, coord')
    @test isapprox(Int.(networkTest.label), label)
    @test isapprox(networkTest.segForce[:, 1, :], segForce')
    @test isapprox(networkTest.segForce[:, 2, :], factor * segForce')
    @test isapprox(networkTest.nodeVel, nodeVel')
    @test networkTest.numNode == 9
    @test networkTest.numSeg == 10
    @test isapprox(networkTest.connectivity, connectivity')
    @test isapprox(networkTest.linksConnect, linksConnect')

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 10, 1)
    @test !compStruct(networkTest, network)
    links = [
        1 2
        2 3
        3 4
        4 5
        5 1
        6 7
        7 8
        8 9
        9 1
        1 6
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
    ]
    coord =
        1.0e+02 * [
            -0.830925136003420 -0.928250034950193 -0.859175170953614
            -1.051615382213988 -1.014440578767969 -1.066055960981957
            -1.060150095500755 -0.939849904499245 -1.000000000000000
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.871749965049807 -0.969074863996580 -0.940824829046386
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            -0.960150095500755 -0.839849904499245 -0.900000000000000
            -0.885559421232031 -0.848384617786012 -0.833944039018043
            0 0 0
        ]
    label = [
        1
        2
        1
        2
        1
        1
        2
        1
        2
        0
    ]
    segForce = [
        1 11 21
        2 12 22
        3 13 23
        4 14 24
        5 15 25
        6 16 26
        7 17 27
        8 18 28
        9 19 29
        10 20 30
    ]
    nodeVel = [
        -10 -20 -30
        -2 -12 -22
        -3 -13 -23
        -4 -14 -24
        -5 -15 -25
        -6 -16 -26
        -7 -17 -27
        -8 -18 -28
        -9 -19 -29
        0 0 0
    ]
    connectivity = [
        4 9 2 10 1 1 1 5 2
        2 1 2 2 1 0 0 0 0
        2 2 2 3 1 0 0 0 0
        2 3 2 4 1 0 0 0 0
        2 4 2 5 1 0 0 0 0
        2 6 1 10 2 0 0 0 0
        2 6 2 7 1 0 0 0 0
        2 7 2 8 1 0 0 0 0
        2 8 2 9 1 0 0 0 0
        0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        3 1
        2 1
        2 1
        2 1
        2 4
        1 1
        2 1
        2 1
        2 1
        2 2
    ]
    @test isapprox(networkTest.links, links')
    @test isapprox(networkTest.slipPlane, slipPlane')
    @test isapprox(networkTest.bVec, bVec')
    @test isapprox(networkTest.coord, coord')
    @test isapprox(Int.(networkTest.label), label)
    @test isapprox(networkTest.segForce[:, 1, :], segForce')
    @test isapprox(networkTest.segForce[:, 2, :], factor * segForce')
    @test isapprox(networkTest.nodeVel, nodeVel')
    @test networkTest.numNode == 9
    @test networkTest.numSeg == 10
    @test isapprox(networkTest.connectivity, connectivity')
    @test isapprox(networkTest.linksConnect, linksConnect')

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 10, 1)
    mergeNode!(networkTest, 1, 3)
    @test !compStruct(networkTest, network)
    links = [
        3 1
        1 6
        1 4
        4 5
        5 1
        6 7
        7 2
        2 3
        0 0
        0 0
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0 0 0
        0 0 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0 0 0
        0 0 0
    ]
    coord =
        1.0e+02 * [
            -0.830925136003420 -0.928250034950193 -0.859175170953614
            -0.960150095500755 -0.839849904499245 -0.900000000000000
            -0.885559421232031 -0.848384617786012 -0.833944039018043
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.871749965049807 -0.969074863996580 -0.940824829046386
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            0 0 0
            0 0 0
            0 0 0
        ]
    label = [
        1
        1
        2
        2
        1
        1
        2
        0
        0
        0
    ]
    segForce = [
        9 19 29
        10 20 30
        3 13 23
        4 14 24
        5 15 25
        6 16 26
        7 17 27
        8 18 28
        0 0 0
        0 0 0
    ]
    nodeVel = [
        -10 -20 -30
        -8 -18 -28
        -9 -19 -29
        -4 -14 -24
        -5 -15 -25
        -6 -16 -26
        -7 -17 -27
        0 0 0
        0 0 0
        0 0 0
    ]
    connectivity = [
        4 1 2 2 1 3 1 5 2 0 0 0 0
        2 7 2 8 1 0 0 0 0 0 0 0 0
        2 8 2 1 1 0 0 0 0 0 0 0 0
        2 3 2 4 1 0 0 0 0 0 0 0 0
        2 4 2 5 1 0 0 0 0 0 0 0 0
        2 6 1 2 2 0 0 0 0 0 0 0 0
        2 6 2 7 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        2 1
        2 2
        3 1
        2 1
        2 4
        1 1
        2 1
        2 1
        0 0
        0 0
    ]
    @test isapprox(networkTest.links, links')
    @test isapprox(networkTest.slipPlane, slipPlane')
    @test isapprox(networkTest.bVec, bVec')
    @test isapprox(networkTest.coord, coord')
    @test isapprox(Int.(networkTest.label), label)
    @test isapprox(networkTest.segForce[:, 1, :], segForce')
    @test isapprox(networkTest.segForce[:, 2, :], factor * segForce')
    @test isapprox(networkTest.nodeVel, nodeVel')
    @test networkTest.numNode == 7
    @test networkTest.numSeg == 8
    @test isapprox(networkTest.connectivity, connectivity')
    @test isapprox(networkTest.linksConnect, linksConnect')

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 10, 1)
    mergeNode!(networkTest, 1, 3)
    mergeNode!(networkTest, 7, 3)
    @test !compStruct(networkTest, network)
    links = [
        3 1
        1 2
        1 4
        4 5
        5 1
        2 3
        0 0
        0 0
        0 0
        0 0
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0 0 0
        0 0 0
        0 0 0
        0 0 0
    ]
    coord =
        1.0e+02 * [
            -0.830925136003420 -0.928250034950193 -0.859175170953614
            -0.871749965049807 -0.969074863996580 -0.940824829046386
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            0 0 0
            0 0 0
            0 0 0
            0 0 0
            0 0 0
        ]
    label = [
        1
        1
        2
        2
        1
        0
        0
        0
        0
        0
    ]
    segForce = [
        9 19 29
        10 20 30
        3 13 23
        4 14 24
        5 15 25
        6 16 26
        0 0 0
        0 0 0
        0 0 0
        0 0 0
    ]
    nodeVel = [
        -10 -20 -30
        -6 -16 -26
        -7 -17 -27
        -4 -14 -24
        -5 -15 -25
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
    ]
    connectivity = [
        4 1 2 2 1 3 1 5 2 0 0 0 0
        2 6 1 2 2 0 0 0 0 0 0 0 0
        2 6 2 1 1 0 0 0 0 0 0 0 0
        2 3 2 4 1 0 0 0 0 0 0 0 0
        2 4 2 5 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        2 1
        2 2
        3 1
        2 1
        2 4
        1 1
        0 0
        0 0
        0 0
        0 0
    ]
    @test isapprox(networkTest.links, links')
    @test isapprox(networkTest.slipPlane, slipPlane')
    @test isapprox(networkTest.bVec, bVec')
    @test isapprox(networkTest.coord, coord')
    @test isapprox(Int.(networkTest.label), label)
    @test isapprox(networkTest.segForce[:, 1, :], segForce')
    @test isapprox(networkTest.segForce[:, 2, :], factor * segForce')
    @test isapprox(networkTest.nodeVel, nodeVel')
    @test networkTest.numNode == 5
    @test networkTest.numSeg == 6
    @test isapprox(networkTest.connectivity, connectivity')
    @test isapprox(networkTest.linksConnect, linksConnect')

    networkTest = deepcopy(network)
    mergeNode!(networkTest, 10, 1)
    mergeNode!(networkTest, 1, 3)
    mergeNode!(networkTest, 7, 3)
    mergeNode!(networkTest, 5, 2)
    @test !compStruct(networkTest, network)
    links = [
        3 1
        2 3
        1 4
        4 2
        0 0
        0 0
        0 0
        0 0
        0 0
        0 0
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
    ]
    coord =
        1.0e+02 * [
            -0.830925136003420 -0.928250034950193 -0.859175170953614
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            0 0 0
            0 0 0
            0 0 0
            0 0 0
            0 0 0
            0 0 0
        ]
    label = [
        1
        1
        2
        2
        0
        0
        0
        0
        0
        0
    ]
    segForce = [
        9 19 29
        6 16 26
        3 13 23
        4 14 24
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
    ]
    nodeVel = [
        -10 -20 -30
        -5 -15 -25
        -7 -17 -27
        -4 -14 -24
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
    ]
    connectivity = [
        2 1 2 3 1 0 0 0 0 0 0 0 0
        2 4 2 2 1 0 0 0 0 0 0 0 0
        2 2 2 1 1 0 0 0 0 0 0 0 0
        2 3 2 4 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        2 1
        2 1
        2 1
        2 1
        0 0
        0 0
        0 0
        0 0
        0 0
        0 0
    ]
    @test isapprox(networkTest.links, links')
    @test isapprox(networkTest.slipPlane, slipPlane')
    @test isapprox(networkTest.bVec, bVec')
    @test isapprox(networkTest.coord, coord')
    @test isapprox(Int.(networkTest.label), label)
    @test isapprox(networkTest.segForce[:, 1, :], segForce')
    @test isapprox(networkTest.segForce[:, 2, :], factor * segForce')
    @test isapprox(networkTest.nodeVel, nodeVel')
    @test networkTest.numNode == 4
    @test networkTest.numSeg == 4
    @test isapprox(networkTest.connectivity, connectivity')
    @test isapprox(networkTest.linksConnect, linksConnect')
end

@testset "Split node" begin
    fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
    fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
    fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
    fileSlipSystem = "../data/slipSystems/SlipSystems.JSON"
    fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
    dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
        fileDislocationP,
        fileMaterialP,
        fileIntegrationP,
        fileSlipSystem,
        fileDislocationLoop,
    )

    pentagon = DislocationLoop(
        loopPrism();
        numSides = 5,
        nodeSide = 1,
        numLoops = 2,
        segLen = 10 * ones(5),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork(pentagon; memBuffer = 1)
    network.coord[:, 6:end] .+= [10; 10; 10]
    factor = rand()
    for j in 1:size(network.segForce, 1)
        for i in 1:size(network.segForce, 3)
            network.segForce[j, 1, i] = i + (j - 1) * size(network.segForce, 3)
            network.segForce[j, 2, i] = network.segForce[j, 1, i] * factor
            network.nodeVel[j, i] = -i - (j - 1) * size(network.segForce, 3)
        end
    end

    networkTest = deepcopy(network)
    midCoord = vec(mean(networkTest.coord, dims = 2))
    midVel = vec(mean(networkTest.nodeVel, dims = 2))
    numNode = networkTest.numNode
    numSeg = networkTest.numSeg
    newEntries = Int(round(11 * log2(11)))
    splitNode!(networkTest, 1, 1, midCoord, midVel)
    links = [
        11 2
        2 3
        3 4
        4 5
        5 1
        6 7
        7 8
        8 9
        9 10
        10 6
        1 11
    ]

    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0.815779215265140 -0.437524113804341 0.378255101460798
    ]

    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
    ]

    label = [
        1
        2
        1
        2
        1
        1
        2
        1
        2
        1
        1
    ]
    coord =
        1e2 * [
            -0.971749965049807 -1.069074863996580 -1.040824829046386
            -1.051615382213988 -1.014440578767969 -1.066055960981957
            -1.060150095500755 -0.939849904499245 -1.000000000000000
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.871749965049807 -0.969074863996580 -0.940824829046386
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            -0.960150095500755 -0.839849904499245 -0.900000000000000
            -0.885559421232031 -0.848384617786012 -0.833944039018043
            -0.830925136003420 -0.928250034950193 -0.859175170953614
            -0.950000000000000 -0.950000000000000 -0.950000000000000
        ]

    nodeVel = [
        -1.000000000000000 -11.000000000000000 -21.000000000000000
        -2.000000000000000 -12.000000000000000 -22.000000000000000
        -3.000000000000000 -13.000000000000000 -23.000000000000000
        -4.000000000000000 -14.000000000000000 -24.000000000000000
        -5.000000000000000 -15.000000000000000 -25.000000000000000
        -6.000000000000000 -16.000000000000000 -26.000000000000000
        -7.000000000000000 -17.000000000000000 -27.000000000000000
        -8.000000000000000 -18.000000000000000 -28.000000000000000
        -9.000000000000000 -19.000000000000000 -29.000000000000000
        -10.000000000000000 -20.000000000000000 -30.000000000000000
        -5.500000000000000 -15.500000000000000 -25.500000000000000
    ]

    connectivity = [
        2 5 2 11 1 0 0 0 0
        2 1 2 2 1 0 0 0 0
        2 2 2 3 1 0 0 0 0
        2 3 2 4 1 0 0 0 0
        2 4 2 5 1 0 0 0 0
        2 6 1 10 2 0 0 0 0
        2 6 2 7 1 0 0 0 0
        2 7 2 8 1 0 0 0 0
        2 8 2 9 1 0 0 0 0
        2 9 2 10 1 0 0 0 0
        2 1 1 11 2 0 0 0 0
    ]

    linksConnect = [
        1 1
        2 1
        2 1
        2 1
        2 1
        1 1
        2 1
        2 1
        2 1
        2 2
        2 2
    ]

    @test networkTest.links[:, 1:11] == links'
    @test isapprox(networkTest.slipPlane[:, 1:11], slipPlane')
    @test isapprox(networkTest.bVec[:, 1:11], bVec')
    @test networkTest.label[1:11] == label
    @test isapprox(networkTest.coord[:, 1:11], coord')
    @test isapprox(networkTest.nodeVel[:, 1:11], nodeVel')
    @test networkTest.connectivity[:, 1:11] == connectivity'
    @test networkTest.linksConnect[:, 1:11] == linksConnect'
    @test size(networkTest.links, 2) ==
          size(networkTest.slipPlane, 2) ==
          size(networkTest.bVec, 2) ==
          length(networkTest.label) ==
          size(networkTest.coord, 2) ==
          size(networkTest.segForce, 3) ==
          size(networkTest.nodeVel, 2) ==
          size(networkTest.connectivity, 2) ==
          size(networkTest.linksConnect, 2) ==
          numNode + newEntries
    @test networkTest.numNode == numNode + 1
    @test networkTest.numSeg == numSeg + 1

    midCoord += [-5, 11, -7]
    midVel += [-6, 1, 9]
    numNode = networkTest.numNode
    numSeg = networkTest.numSeg
    splitNode!(networkTest, 8, 1, midCoord, midVel)

    links = [
        11 2
        2 3
        3 4
        4 5
        5 1
        6 7
        7 12
        8 9
        9 10
        10 6
        1 11
        12 8
    ]

    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0.815779215265140 -0.437524113804341 0.378255101460798
        -0.589338571143760 0.784066911661008 0.194728340517247
    ]

    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
    ]

    label = [
        1
        2
        1
        2
        1
        1
        2
        1
        2
        1
        1
        1
    ]

    coord =
        1.0e+02 * [
            -0.971749965049807 -1.069074863996580 -1.040824829046386
            -1.051615382213988 -1.014440578767969 -1.066055960981957
            -1.060150095500755 -0.939849904499245 -1.000000000000000
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -0.930925136003420 -1.028250034950193 -0.959175170953614
            -0.871749965049807 -0.969074863996580 -0.940824829046386
            -0.951615382213988 -0.914440578767969 -0.966055960981957
            -0.960150095500755 -0.839849904499245 -0.900000000000000
            -0.885559421232031 -0.848384617786012 -0.833944039018043
            -0.830925136003420 -0.928250034950193 -0.859175170953614
            -0.950000000000000 -0.950000000000000 -0.950000000000000
            -1.000000000000000 -0.840000000000000 -1.020000000000000
        ]

    nodeVel = [
        -1.000000000000000 -11.000000000000000 -21.000000000000000
        -2.000000000000000 -12.000000000000000 -22.000000000000000
        -3.000000000000000 -13.000000000000000 -23.000000000000000
        -4.000000000000000 -14.000000000000000 -24.000000000000000
        -5.000000000000000 -15.000000000000000 -25.000000000000000
        -6.000000000000000 -16.000000000000000 -26.000000000000000
        -7.000000000000000 -17.000000000000000 -27.000000000000000
        -8.000000000000000 -18.000000000000000 -28.000000000000000
        -9.000000000000000 -19.000000000000000 -29.000000000000000
        -10.000000000000000 -20.000000000000000 -30.000000000000000
        -5.500000000000000 -15.500000000000000 -25.500000000000000
        -11.500000000000000 -14.500000000000000 -16.500000000000000
    ]

    connectivity = [
        2 5 2 11 1 0 0 0 0
        2 1 2 2 1 0 0 0 0
        2 2 2 3 1 0 0 0 0
        2 3 2 4 1 0 0 0 0
        2 4 2 5 1 0 0 0 0
        2 6 1 10 2 0 0 0 0
        2 6 2 7 1 0 0 0 0
        2 8 1 12 2 0 0 0 0
        2 8 2 9 1 0 0 0 0
        2 9 2 10 1 0 0 0 0
        2 1 1 11 2 0 0 0 0
        2 7 2 12 1 0 0 0 0
    ]

    linksConnect = [
        1 1
        2 1
        2 1
        2 1
        2 1
        1 1
        2 1
        1 1
        2 1
        2 2
        2 2
        2 2
    ]
    # Check that no more memory was allocated.
    @test size(networkTest.links, 2) ==
          size(networkTest.slipPlane, 2) ==
          size(networkTest.bVec, 2) ==
          length(networkTest.label) ==
          size(networkTest.coord, 2) ==
          size(networkTest.segForce, 3) ==
          size(networkTest.nodeVel, 2) ==
          size(networkTest.connectivity, 2) ==
          size(networkTest.linksConnect, 2) ==
          48
    @test networkTest.numNode == numNode + 1
    @test networkTest.numSeg == numSeg + 1
    @test networkTest.links[:, 1:12] == links'
    @test isapprox(networkTest.slipPlane[:, 1:12], slipPlane')
    @test isapprox(networkTest.bVec[:, 1:12], bVec')
    @test networkTest.label[1:12] == label
    @test isapprox(networkTest.coord[:, 1:12], coord')
    @test isapprox(networkTest.nodeVel[:, 1:12], nodeVel')
    @test networkTest.connectivity[:, 1:12] == connectivity'
    @test networkTest.linksConnect[:, 1:12] == linksConnect'
    @test_throws AssertionError splitNode!(networkTest, 8, 3, midCoord, midVel)
end

@testset "Coarsen network" begin
    fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
    fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
    fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
    fileSlipSystem = "../data/slipSystems/SlipSystems.JSON"
    fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
    dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
        fileDislocationP,
        fileMaterialP,
        fileIntegrationP,
        fileSlipSystem,
        fileDislocationLoop,
    )

    prismPentagon = DislocationLoop(
        loopPrism();
        numSides = 5,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(5),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )

    shearHexagon = DislocationLoop(
        loopShear();
        numSides = 6,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(6),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork([shearHexagon, prismPentagon], memBuffer = 1)

    network2 = deepcopy(network)
    calcSegForce!(dlnParams, matParams, network2)
    coarsenNetwork!(dlnParams, matParams, network2)
    links = [
        6 2
        2 3
        3 4
        4 6
        5 1
        1 7
        7 8
        8 5
        0 0
        0 0
        0 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0 0 0
        0 0 0
        0 0 0
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        0 0 0
        0 0 0
        0 0 0
    ]
    coord =
        1.0e+02 * [
            -0.985559421232031 -0.948384617786012 -0.933944039018043
            -1.040824829046386 -1.040824829046386 -1.081649658092773
            -1.070412414523193 -1.070412414523193 -0.990824829046386
            -1.029587585476807 -1.029587585476807 -0.909175170953614
            -1.060150095500755 -0.939849904499245 -1.000000000000000
            -0.929587585476807 -0.929587585476807 -1.009175170953614
            -0.971749965049807 -1.069074863996580 -1.040824829046386
            -1.051615382213988 -1.014440578767969 -1.066055960981957
            0 0 0
            0 0 0
            0 0 0
        ]
    label = [
        2
        2
        1
        2
        1
        1
        1
        2
        0
        0
        0
    ]
    nodeVel = [
        -0.000006689570426 -0.000006689583194 0.000006689512325
        0.094210513111597 0.094210513111611 0.225800116243039
        0 0 0
        0.082778555488260 0.082778555488232 -0.181030113207538
        0 0 0
        -0.197155144338323 -0.197155144338295 0.050824276836052
        -0.000032128097977 -0.000032128031107 0.000032128094986
        0 0 0
        0 0 0
        0 0 0
        0 0 0
    ]
    segForce1 = [
        -0.731165101601140 -0.729070538387561 -0.845685760271632
        -0.107679944143240 -0.108364876814072 1.137684251104435
        0.568002311233386 0.566611169671008 1.133727803147232
        0.577108924635657 0.577103498327851 -0.577838384605370
        1.036545887927879 -0.119085496493695 0.916735911922166
        0.117275438015306 -1.035927969236551 -0.917667535670576
        -1.109169644845321 0.759306566573759 -0.349584829798629
        -0.117507614989553 1.035899543796946 0.917669635538588
        0 0 0
        0 0 0
        0 0 0
    ]
    segForce2 = [
        0.728639551563188 0.730718822192315 0.847032229650163
        0.109020042708812 0.108325930005808 -1.137260380900199
        -0.566020717832909 -0.567411859395269 -1.134318254985340
        -0.577594899151688 -0.577589472843882 0.576866435573308
        -1.035418014951076 0.117989143835423 -0.918151164384458
        -0.119798333631734 1.036035282908816 0.917220387684654
        1.109315417743939 -0.758311215081191 0.351278693716326
        0.119568482835041 -1.036062901586534 -0.917218898263512
        0 0 0
        0 0 0
        0 0 0
    ]
    connectivity = [
        2 5 2 6 1 0 0 0 0
        2 1 2 2 1 0 0 0 0
        2 2 2 3 1 0 0 0 0
        2 3 2 4 1 0 0 0 0
        2 8 2 5 1 0 0 0 0
        2 4 2 1 1 0 0 0 0
        2 7 1 6 2 0 0 0 0
        2 7 2 8 1 0 0 0 0
        0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0
    ]
    linksConnect = [
        2 1
        2 1
        2 1
        2 1
        2 1
        2 2
        1 1
        2 1
        0 0
        0 0
        0 0
    ]

    @test isapprox(network2.links', links)
    @test isapprox(network2.slipPlane', slipPlane)
    @test isapprox(network2.bVec', bVec)
    @test isapprox(network2.coord', coord)
    @test network2.label == label
    @test isapprox(network2.nodeVel', nodeVel, rtol = 1e-6)
    @test network2.numNode == 8
    @test network2.numSeg == 8
    @test isapprox(network2.segForce[:, 1, :]', segForce1)
    @test isapprox(network2.segForce[:, 2, :]', segForce2)
    @test isapprox(network2.connectivity', connectivity)
    @test isapprox(network2.linksConnect', linksConnect)
end

@testset "Refine network" begin
    fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
    fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
    fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
    fileSlipSystem = "../data/slipSystems/SlipSystems.JSON"
    fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
    dlnParams, matParams, intParams, slipSystems, missing = loadParams(
        fileDislocationP,
        fileMaterialP,
        fileIntegrationP,
        fileSlipSystem,
        fileDislocationLoop,
    )

    prismPentagon = DislocationLoop(
        loopPrism();
        numSides = 5,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(5),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )

    shearHexagon = DislocationLoop(
        loopShear();
        numSides = 6,
        nodeSide = 1,
        numLoops = 1,
        segLen = 10 * ones(6),
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 2; 1; 2; 1; 1],
        buffer = 0.0,
        range = Float64[-100 100; -100 100; -100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork([shearHexagon, prismPentagon], memBuffer = 1)
    network2 = deepcopy(network)
    refineNetwork!(dlnParams, matParams, network2)
    @test compStruct(network, network2)
end

@testset "Coarsen and Refine network" begin
    fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
    fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
    fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
    fileSlipSystem = "../data/slipSystems/SlipSystems.JSON"
    fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
    dlnParams, matParams, intParams, slipSystems, missing = loadParams(
        fileDislocationP,
        fileMaterialP,
        fileIntegrationP,
        fileSlipSystem,
        fileDislocationLoop,
    )
    shearDecagon = DislocationLoop(
        loopShear();
        numSides = 10,
        nodeSide = 1,
        numLoops = 1,
        segLen = [300; 700; 1100; 1500; 1900; 1900; 1500; 1100; 700; 300],#,300; 700; 1100; 1500; 1900
        slipSystem = 4,
        _slipPlane = slipSystems.slipPlane[:, 4],
        _bVec = slipSystems.bVec[:, 4],
        label = nodeType[1; 1; 1; 1; 1; 1; 1; 1; 1; 1],
        buffer = 0.0,
        range = Float64[0 0; 0 0; 0 0],
        dist = Zeros(),
    )

    network = DislocationNetwork(shearDecagon, memBuffer = 1)
    calcSegForce!(dlnParams, matParams, network)
    dlnMobility(dlnParams, matParams, network)
    coarsenNetwork!(dlnParams, matParams, network)
    refineNetwork!(dlnParams, matParams, network)
    links = [
        11 2
        2 3
        3 4
        12 5
        13 6
        6 7
        7 8
        8 9
        14 10
        10 1
        1 11
        4 12
        5 13
        9 14
    ]
    slipPlane = [
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
        -0.707106781186547 0.707106781186547 0
    ]
    bVec = [
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
        0.577350269189626 0.577350269189626 -0.577350269189626
    ]
    label = [
        1
        1
        1
        1
        1
        1
        1
        1
        1
        1
        1
        1
        1
        1
    ]
    coord =
        1.0e+03 * [
            0.398589381337062 0.398589381337062 -2.768451689518856
            -0.070157063645514 -0.070157063645514 -2.993292834846829
            -0.812930298558462 -0.812930298558462 -2.666833276870797
            -1.447335912588014 -1.447335912588014 -1.464727194290905
            -1.464584433169251 -1.464584433169251 0.435116214276106
            -0.688912681287911 -0.688912681287911 1.986459718038786
            0.315543986531893 0.315543986531893 2.468262172313013
            1.058317221444842 1.058317221444842 2.141802614336981
            1.357096625276723 1.357096625276723 1.280844500762801
            0.877843003306892 0.877843003306892 -0.743803594378028
            0.164216158845774 0.164216158845774 -2.880872262182843
            -1.455960172878632 -1.455960172878632 -0.514805490007400
            -1.076748557228581 -1.076748557228581 1.210787966157446
            1.117469814291808 1.117469814291808 0.268520453192387
        ]
    nodeVel = [
        0.000222065664539 0.000222065664539 0.002382975226818
        -0.001873207466708 -0.001873207466708 0.001493537050782
        0 0 0
        0.000651728292191 0.000651728292191 -0.000426938702841
        0.000449561250821 0.000449561250821 -0.000204569234931
        0.000411949365045 0.000411949365045 -0.000288275324690
        0 0 0
        -0.000704586698816 -0.000704586698816 -0.000660901785972
        -0.000686592757211 -0.000686592757211 0.000004734099786
        0.000066419652455 0.000066419652427 0.000280595407503
        0.000052846123310 0.000052846123310 0.000025348406879
        -0.000082703071257 -0.000082703071229 0.009109351081861
        -0.000068503691857 -0.000068503691857 -0.000137007195931
        0.000000000000114 0.000000000000142 0.000369912943938
    ]
    connectivity = [
        2 10 2 11 1 0 0 0 0
        2 1 2 2 1 0 0 0 0
        2 2 2 3 1 0 0 0 0
        2 3 2 12 1 0 0 0 0
        2 4 2 13 1 0 0 0 0
        2 5 2 6 1 0 0 0 0
        2 6 2 7 1 0 0 0 0
        2 7 2 8 1 0 0 0 0
        2 8 2 14 1 0 0 0 0
        2 9 2 10 1 0 0 0 0
        2 1 1 11 2 0 0 0 0
        2 4 1 12 2 0 0 0 0
        2 5 1 13 2 0 0 0 0
        2 9 1 14 2 0 0 0 0
    ]
    linksConnect = [
        1 1
        2 1
        2 1
        1 1
        1 1
        2 1
        2 1
        2 1
        1 1
        2 1
        2 2
        2 2
        2 2
        2 2
    ]
    segForce1 = [
        -0.766768378154913 -0.766768378154913 -0.708641892768090
        -0.718396487134110 -0.718396487134110 0.225564777238983
        -0.212518734977064 -0.212518734977064 1.069836712554321
        0.293298957155174 0.293298957155174 1.259960158642018
        0.599880026908032 0.599880026908032 1.101154516691398
        0.786643951604992 0.786643951604992 0.625769121566676
        0.718235597298696 0.718235597298696 -0.226296902616202
        0.035877823434913 0.035877823434913 -1.188271431333817
        -0.484980547736341 -0.484980547736341 -1.200959726167859
        -0.523671777310287 -0.523671777310287 -1.182642555933405
        -0.769213753534785 -0.769213753534785 -0.698445707240829
        0.319365221259123 0.319365221259123 1.260433465537571
        0.613190147666069 0.613190147666069 1.087844395933361
        -0.527086373979404 -0.527086373979404 -1.181026020373306
    ]
    segForce2 = [
        0.744054365041219 0.744054365041219 0.803349762031084
        0.778598719650432 0.778598719650432 0.048384052884796
        0.400271787527061 0.400271787527061 -0.871665199548957
        -0.191415295624645 -0.191415295624645 -1.258110172032725
        -0.519489638440782 -0.519489638440782 -1.181544905158647
        -0.722099487828575 -0.722099487828575 -0.894892354277212
        -0.780230986160347 -0.780230986160347 -0.055811642895184
        -0.221361739665044 -0.221361739665044 1.059533986299295
        0.413834851518675 0.413834851518675 1.234641463435401
        0.330611424565010 0.330611424565010 1.274041030101837
        0.745279830964200 0.745279830964200 0.798240085205638
        -0.202020359230095 -0.202020359230095 -1.258302737016806
        -0.534801549457458 -0.534801549457458 -1.166232994141972
        0.398175155005969 0.398175155005969 1.242055064140309
    ]
    numNode = 14
    numSeg = 14
    @test network.numNode == numNode
    @test network.numSeg == numSeg
    @test network.links[:, 1:numSeg]' == links
    @test isapprox(network.slipPlane[:, 1:numSeg]', slipPlane)
    @test isapprox(network.bVec[:, 1:numSeg]', bVec)
    @test network.label[1:numSeg] == label
    @test isapprox(network.coord[:, 1:numSeg]', coord)
    @test isapprox(network.nodeVel[:, 1:9]', nodeVel[1:9, :], atol = 1e-8)
    @test isapprox(network.nodeVel[:, 10:numNode]', nodeVel[10:end, :], atol = 1e-2)
    @test network.connectivity[:, 1:numSeg]' == connectivity
    @test network.linksConnect[:, 1:numSeg]' == linksConnect
    @test isapprox(network.segForce[:, 1, 1:numSeg]', segForce1)
    @test isapprox(network.segForce[:, 2, 1:numSeg]', segForce2)
end
