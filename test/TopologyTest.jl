using Revise
using DDD
using Test, Statistics
cd(@__DIR__)

# @testset "Merge nodes" begin
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

for j in 1:size(network.segForce, 1)
    for i in 1:size(network.segForce, 2)
        network.segForce[j, i] = i + (j-1)*size(network.segForce, 2)
        network.nodeVel[i] = -i - (j-1)*size(network.segForce, 2)
    end
end

network.segForce'
for i in eachindex(network.segForce)
    network.segForce[i] = i
end
for i in eachindex(network.nodeVel)
    network.nodeVel[i] = -i
end
network.segForce'
network.nodeVel

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

segForce










network.segForce




networkTest.segForce


@test isapprox(networkTest.links, links')
@test isapprox(networkTest.slipPlane, slipPlane')
@test isapprox(networkTest.bVec, bVec')
@test isapprox(networkTest.coord, coord')
@test isapprox(Int.(networkTest.label), label)
@test isapprox(networkTest.segForce, segForce')
@test isapprox(networkTest.nodeVel, nodeVel')
@test networkTest.numNode == 9
@test networkTest.numSeg == 9
@test isapprox(networkTest.connectivity, connectivity)
@test isapprox(networkTest.linksConnect, linksConnect)

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
@test isapprox(networkTest.links, links)
@test isapprox(networkTest.slipPlane, slipPlane)
@test isapprox(networkTest.bVec, bVec)
@test isapprox(networkTest.coord, coord)
@test isapprox(Int.(networkTest.label), label)
@test isapprox(networkTest.segForce, segForce)
@test isapprox(networkTest.nodeVel, nodeVel)
@test networkTest.numNode == 8
@test networkTest.numSeg == 8
@test isapprox(networkTest.connectivity, connectivity)
@test isapprox(networkTest.linksConnect, linksConnect)

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
@test isapprox(networkTest.links, links)
networkTest.slipPlane
@test isapprox(networkTest.slipPlane, slipPlane)
@test isapprox(networkTest.bVec, bVec)
@test isapprox(networkTest.coord, coord)
@test isapprox(Int.(networkTest.label), label)
@test isapprox(networkTest.segForce, segForce)
@test isapprox(networkTest.nodeVel, nodeVel)
@test networkTest.numNode == 8
@test networkTest.numSeg == 8
@test isapprox(networkTest.connectivity, connectivity)
@test isapprox(networkTest.linksConnect, linksConnect)

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
@test isapprox(networkTest.links, links)
@test isapprox(networkTest.slipPlane, slipPlane)
@test isapprox(networkTest.bVec, bVec)
@test isapprox(networkTest.coord, coord)
@test isapprox(Int.(networkTest.label), label)
@test isapprox(networkTest.segForce, segForce)
@test isapprox(networkTest.nodeVel, nodeVel)
@test networkTest.numNode == 8
@test networkTest.numSeg == 8
@test isapprox(networkTest.connectivity, connectivity)
@test isapprox(networkTest.linksConnect, linksConnect)

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
@test isapprox(networkTest.links, links)
@test isapprox(networkTest.slipPlane, slipPlane)
@test isapprox(networkTest.bVec, bVec)
@test isapprox(networkTest.coord, coord)
@test isapprox(Int.(networkTest.label), label)
@test isapprox(networkTest.segForce, segForce)
@test isapprox(networkTest.nodeVel, nodeVel)
@test networkTest.numNode == 9
@test networkTest.numSeg == 10
@test isapprox(networkTest.connectivity, connectivity)
@test isapprox(networkTest.linksConnect, linksConnect)

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
@test isapprox(networkTest.links, links)
@test isapprox(networkTest.slipPlane, slipPlane)
@test isapprox(networkTest.bVec, bVec)
@test isapprox(networkTest.coord, coord)
@test isapprox(Int.(networkTest.label), label)
@test isapprox(networkTest.segForce, segForce)
@test isapprox(networkTest.nodeVel, nodeVel)
@test networkTest.numNode == 9
@test networkTest.numSeg == 10
@test isapprox(networkTest.connectivity, connectivity)
@test isapprox(networkTest.linksConnect, linksConnect)

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
@test isapprox(networkTest.links, links)
@test isapprox(networkTest.slipPlane, slipPlane)
@test isapprox(networkTest.bVec, bVec)
@test isapprox(networkTest.coord, coord)
@test isapprox(Int.(networkTest.label), label)
@test isapprox(networkTest.segForce, segForce)
@test isapprox(networkTest.nodeVel, nodeVel)
@test networkTest.numNode == 7
@test networkTest.numSeg == 8
@test isapprox(networkTest.connectivity, connectivity)
@test isapprox(networkTest.linksConnect, linksConnect)

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
@test isapprox(networkTest.links, links)
@test isapprox(networkTest.slipPlane, slipPlane)
@test isapprox(networkTest.bVec, bVec)
@test isapprox(networkTest.coord, coord)
@test isapprox(Int.(networkTest.label), label)
@test isapprox(networkTest.segForce, segForce)
@test isapprox(networkTest.nodeVel, nodeVel)
@test networkTest.numNode == 5
@test networkTest.numSeg == 6
@test isapprox(networkTest.connectivity, connectivity)
@test isapprox(networkTest.linksConnect, linksConnect)

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
@test isapprox(networkTest.links, links)
@test isapprox(networkTest.slipPlane, slipPlane)
@test isapprox(networkTest.bVec, bVec)
@test isapprox(networkTest.coord, coord)
@test isapprox(Int.(networkTest.label), label)
@test isapprox(networkTest.segForce, segForce)
@test isapprox(networkTest.nodeVel, nodeVel)
@test networkTest.numNode == 4
@test networkTest.numSeg == 4
@test isapprox(networkTest.connectivity, connectivity)
@test isapprox(networkTest.linksConnect, linksConnect)
# end

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
        _slipPlane = slipSystems.slipPlane[4, :],
        _bVec = slipSystems.bVec[4, :],
        label = nodeType[1; 2; 1; 2; 1],
        buffer = 0.0,
        range = Float64[-100 -100 -100; 100 100 100],
        dist = Zeros(),
    )
    network = DislocationNetwork(pentagon; memBuffer = 1)
    network.coord[6:end, :] .+= [10 10 10]
    for i in eachindex(network.segForce)
        network.segForce[i] = i
    end
    for i in eachindex(network.nodeVel)
        network.nodeVel[i] = -i
    end

    networkTest = deepcopy(network)
    midCoord = vec(mean(networkTest.coord, dims = 1))
    midVel = vec(mean(networkTest.nodeVel, dims = 1))
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

    @test networkTest.links[1:11, :] == links
    @test isapprox(networkTest.slipPlane[1:11, :], slipPlane)
    @test isapprox(networkTest.bVec[1:11, :], bVec)
    @test networkTest.label[1:11] == label
    @test isapprox(networkTest.coord[1:11, :], coord)
    @test isapprox(networkTest.nodeVel[1:11, :], nodeVel)
    @test networkTest.connectivity[1:11, :] == connectivity
    @test networkTest.linksConnect[1:11, :] == linksConnect
    @test size(networkTest.links, 1) ==
          size(networkTest.slipPlane, 1) ==
          size(networkTest.bVec, 1) ==
          length(networkTest.label) ==
          size(networkTest.coord, 1) ==
          size(networkTest.segForce, 1) ==
          size(networkTest.nodeVel, 1) ==
          size(networkTest.connectivity, 1) ==
          size(networkTest.linksConnect, 1) ==
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
    @test size(networkTest.links, 1) ==
          size(networkTest.slipPlane, 1) ==
          size(networkTest.bVec, 1) ==
          length(networkTest.label) ==
          size(networkTest.coord, 1) ==
          size(networkTest.segForce, 1) ==
          size(networkTest.nodeVel, 1) ==
          size(networkTest.connectivity, 1) ==
          size(networkTest.linksConnect, 1) ==
          48
    @test networkTest.numNode == numNode + 1
    @test networkTest.numSeg == numSeg + 1
    @test networkTest.links[1:12, :] == links
    @test isapprox(networkTest.slipPlane[1:12, :], slipPlane)
    @test isapprox(networkTest.bVec[1:12, :], bVec)
    @test networkTest.label[1:12] == label
    @test isapprox(networkTest.coord[1:12, :], coord)
    @test isapprox(networkTest.nodeVel[1:12, :], nodeVel)
    @test networkTest.connectivity[1:12, :] == connectivity
    @test networkTest.linksConnect[1:12, :] == linksConnect
    @test_throws AssertionError splitNode!(networkTest, 8, 3, midCoord, midVel)
end
