##
using Plots
plotlyjs()
##
using BenchmarkTools, LinearAlgebra, StaticArrays, SparseArrays, DDD
cd(@__DIR__)
fileDislocationParameters = "../inputs/simParams/sampleDislocationParameters.json"
fileMaterialParameters = "../inputs/simParams/sampleMaterialParameters.json"
fileFEMParameters = "../inputs/simParams/sampleFEMParameters.json"
fileIntegrationParameters = "../inputs/simParams/sampleIntegrationParameters.json"
fileSlipSystem = "../data/slipSystems/BCC.json"
fileDislocationLoop = "../inputs/dln/samplePrismShear.json"
fileIntVar = "../inputs/simParams/sampleIntegrationTime.json"
dlnParams, matParams, femParams, intParams, slipSystems, dislocationLoop =
    loadParametersJSON(
        fileDislocationParameters,
        fileMaterialParameters,
        fileFEMParameters,
        fileIntegrationParameters,
        fileSlipSystem,
        fileDislocationLoop,
    )
intVars = loadIntegrationTimeJSON(fileIntVar)
# network = DislocationNetwork(dislocationLoop)
dx, dy, dz = femParams.dx, femParams.dy, femParams.dz

segLen = rand() * (dx * dy * dz) / 1e8

@btime regularCuboidMesh = buildMesh(matParams, femParams)
numFEMNode = regularCuboidMesh.numNode

f = sparsevec(
    [112, 118, 133, 141, 213, 244, 262, 272, 317, 3 * numFEMNode],
    [
        0.43048187784858616,
        0.22724536603830137,
        0.4340867899691503,
        0.6660863546953892,
        0.30358515797696106,
        0.2945958951093859,
        0.7278367502911502,
        0.7095924334694701,
        0.1642050526375538,
        0,
    ],
)
fHat = sparsevec(
    [32, 48, 55, 88, 138, 148, 191, 230, 253, 335, 3 * numFEMNode],
    [
        0.09706224225842108,
        0.07773687633248638,
        0.13682398802299178,
        0.4752286167553166,
        0.7423196193496164,
        0.8286077556473421,
        0.7023632196408749,
        0.9813639162461198,
        0.5296701796678411,
        0.5523797553266823,
        0,
    ],
)
u = sparsevec(
    [30, 127, 195, 221, 316, 325, 338, 348, 370, 3 * numFEMNode],
    [
        0.8792592573507609,
        0.8430664083925272,
        0.4711050560756602,
        0.4860071865093816,
        0.7905698600135145,
        0.39047211692578077,
        0.6545538020629462,
        0.5446700211111557,
        0.8865721648558644,
        0,
    ],
)
uHat = sparsevec(
    [91, 126, 130, 195, 217, 226, 229, 256, 281, 293, 309, 342, 3 * numFEMNode],
    [
        0.5231621885339968,
        0.5771429489788034,
        0.7151190318538345,
        0.7283662326812077,
        0.6314274719472075,
        0.9814688915693632,
        0.5672795171250207,
        0.002712918060655989,
        0.1788941754890383,
        0.188299784057536,
        0.8489027048214433,
        0.029995302953659708,
        0,
    ],
)
forceDisplacement = ForceDisplacement(u * 1000, f * 1000, uHat * 1000, fHat * 1000)


σHat = calc_σHat(regularCuboidMesh, forceDisplacement, [1575.0, 985.0, 1341.0])
σHatTest = [
    -0.023035166204661 -0.155651908782923 0
    -0.155651908782923 -0.059233284526271 -0.015024315519587
    0 -0.015024315519587 -0.023035166204661
]
isapprox(σHat, σHatTest)

σHat = calc_σHat(regularCuboidMesh, forceDisplacement, [1893.0, 408.0, 1782.0])
σHatTest = [
    -0.607540206946205 0 -0.972551012187583
    0 -0.607540206946205 -0.265466730367529
    -0.972551012187583 -0.265466730367529 -1.562246246433098
]
isapprox(σHat, σHatTest)

@btime σHat = calc_σHat(regularCuboidMesh, forceDisplacement, [1575.0, 985.0, 1341.0])



Bold = copy(B)
isapprox(Bold, B)
##
prismPentagon = DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = segLen * SVector{5}(ones(5)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 2, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 2],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 2],            # Burgers vector of the segments.
    label = SVector{5,nodeType}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(
        0 + segLen,
        0 + segLen,
        0 + segLen,
        dx - segLen,
        dy - segLen,
        dz - segLen,
    ),  # Distribution range
    dist = Rand(),  # Loop distribution.
)
prismHeptagon = DislocationLoop(;
    loopType = loopPrism(),    # Shear loop
    numSides = 7,
    nodeSide = 1,   # 3 nodes per side, it devides the side into equal segments.
    numLoops = 20,
    segLen = segLen * SVector{7}(ones(7)),  # The hexagon's side length is 10, each segment is 10/3.
    slipSystem = 1,
    _slipPlane = slipSystems.slipPlane[:, 1],
    _bVec = slipSystems.bVec[:, 1],
    label = SVector{7,nodeType}(1, 1, 1, 1, 1, 2, 1),
    buffer = 0,
    range = SMatrix{3,2,Float64}(
        0 + segLen,
        0 + segLen,
        0 + segLen,
        dx - segLen,
        dy - segLen,
        dz - segLen,
    ),  # Distribution range
    dist = Rand(),
)

network = DislocationNetwork((prismHeptagon, prismPentagon))
using Plots
gr()
plotNodes(
    regularCuboidMesh,
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,
)
plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,
)
##




##
nodeEl = 1:8 # Local node numbers.
dofLocal = Tuple(Iterators.flatten((3 * (nodeEl .- 1) .+ 1,
    3 * (nodeEl .- 1) .+ 2,
    3 * (nodeEl .- 1) .+ 3,)))

test = [
    -0.006220084679281 0.006220084679281 -0.006220084679281
    0.006220084679281 0.001666666666667 -0.001666666666667
    0.001666666666667 0.000446581987385 0.001666666666667
    -0.001666666666667 0.001666666666667 0.006220084679281
    -0.001666666666667 -0.006220084679281 -0.001666666666667
    0.001666666666667 -0.001666666666667 -0.000446581987385
    0.000446581987385 -0.000446581987385 0.000446581987385
    -0.000446581987385 -0.001666666666667 0.001666666666667
]
isapprox(nx', test)
##

@time constructMesh(matParams, dx, dy, dz, mx, my, mz)



x = coord[1, connect[1, :]]
y = coord[2, connect[1, :]]
z = coord[3, connect[1, :]]

plotlyjs()
figure = scatter(x, y, z)

figure = scatter(mesh[1, :], mesh[2, :], mesh[3, :])

realCoord = Array{SMatrix{3,8}}(undef, 8)
size(N[1])
size(dNdS[1])
dNdS[1]


p = 1 / sqrt(3)
gaussNodes = SMatrix{3,8}(
    -p,
    -p,
    -p,
    p,
    -p,
    -p,
    p,
    p,
    -p,
    -p,
    p,
    -p,
    -p,
    -p,
    p,
    p,
    -p,
    p,
    p,
    p,
    p,
    -p,
    p,
    p,
)
typeof(gaussNodes[1, :]) <: AbstractVector
##
test = normalize(SVector{3}(rand(3)))
normalize!(test)
println("ASDF")
##


using Polyhedra
import GLPK
lib = DefaultLibrary{Float64}(GLPK.Optimizer)
DefaultLibrary

P1 = polyhedron(vrep([
    -1.9 -1.7
    -1.8 0.5
    1.7 0.7
    1.9 -0.3
    0.9 -1.1
]))

vrep([
    -1.9 -1.7
    -1.8 0.5
    1.7 0.7
    1.9 -0.3
    0.9 -1.1
])


using PyCall
using Conda
using QHull
dx = 2000
dy = 2000
dz = 2000

vertices = [
    0 0 0
    dx 0 0
    0 dy 0
    dx dy 0
    0 0 dz
    dx 0 dz
    0 dy dz
    dx dy dz
]
vertices = [0, 0, 0; dx,0, 0; 0,dy, 0; dx,dy, 0; 0,0, dz; dx, 0, dz; 0, dy, dz; dx, dy, dz];


test = points(vrep(vertices))
P2 = polyhedron(vrep(vertices))
removevredundancy!(P2)

@show hrep(P2)
@show vrep(P2)

in(P2, polyhedron(vrep([-10.0 -5 -6])))


ininterior(vrep(vertices), vrep([1 2 3]))
using Plots
plot(P2, color = "red", alpha = 0.2)
convexhull(P2)

using QHull
P2

hrep(P2)

p = rand(10, 2)
ch = chull(p)
ch.points         # original points
ch.vertices       # indices to line segments forming the convex hull
ch.simplices      # the simplexes forming the convex hull
show(ch)
##
remoteForce = calcSegSegForce(dlnParams, matParams, network)
selfForce = calcSegForce(dlnParams, matParams, network)
calcSegForce!(dlnParams, matParams, network)
isapprox(network.segForce[:, :, 1:network.numSeg[1]], selfForce)

@btime calcSegForce!(dlnParams, matParams, network)

using BenchmarkTools

a = rand(6, 6)
b = rand(6)

function foo(a, b)
    a[2:3, 3] = @view b[3:4]
    return a
end
foo(a, b)
function bar(a, b)
    a[2, 3] = b[3]
    a[4, 3] = b[4]
    return a
end
bar(a, b)

a = rand(10, 10)
b = rand(10, 10)
i = rand(1:10)
i2 = rand(1:10)
j = rand(1:10)
j2 = rand(1:10)

function comparing(a, b, i, j, i2, j2)
    a[:, i] == b[:, j] ? a[:, i2] = b[:, j2] : nothing
end
comparing(a, b, i, j, i2, j2)
function comparing2(a, b, i, j, i2, j2)
    isapprox(a[:, i], b[:, j]) ? a[:, i2] = b[:, j2] : nothing
end
comparing2(a, b, i, j, i2, j2)
function comparing2(a, b, i, j, i2, j2)
    isapprox(a[:, i], b[:, j]) ? a[:, i2] = b[:, j2] : nothing
end
comparing2(a, b, i, j, i2, j2)
function comparing3(a, b, i, j, i2, j2)

    equalSlipPlane = let
        flag = true
        for k = 1:3
            flag = flag && isapprox(a[k, i], b[k, j])
        end
        flag
    end
    equalSlipPlane ? for i = 1:3
        a[i, i2] = b[i, j2]
    end : nothing
end
comparing3(a, b, i, j, i2, j2)

[i for i = 1:3]

@btime foo(a, b)
@btime bar(a, b)
@btime comparing(a, b, i, j, i2, j2)
@btime comparing2(a, b, i, j, i2, j2)
@btime comparing3(a, b, i, j, i2, j2)

@btime comparing(a, a, i, i, i2, j2)
@btime comparing2(a, a, i, i, i2, j2)
@btime comparing3(a, a, i, i, i2, j2)
a[1:3, i2]
b[1:3, j2]

c = rand(10000, 10000)
d = rand(10000, 10000)
e = rand(10000, 10000)

function sendit(c, d, e)
    c[500:5000, 10000] .= 0
    d[500:5000, 10000] .= 0
    e[500:5000, 10000] .= 0
    return c, d, e
end
sendit(c, d, e)
function sendit2(c, d, e)
    @inbounds @simd for i = 500:5000
        c[i, 10000] = 0
        d[i, 10000] = 0
        e[i, 10000] = 0
    end
    return c, d, e
end
sendit2(c, d, e)

@btime sendit(c, d, e)
@btime sendit2(c, d, e)

for j in 1:2, i in 1:3
    println(i, j)
end

test = rand(10, 5)
using StaticArrays, LinearAlgebra
wat = SVector{3,Float64}(3, 3, 3)

function watanabe(a, b)
    a[3:5, 4] = b / norm(b)
    return a
end

@btime watanabe(test, wat)

##
prismPentagon = DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = 500 * SVector{5}(ones(5)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 2, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 2],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 2],            # Burgers vector of the segments.
    label = SVector{5,nodeType}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
    dist = Rand(),  # Loop distribution.
)

prismHeptagon = DislocationLoop(;
    loopType = loopPrism(),    # Shear loop
    numSides = 7,
    nodeSide = 1,   # 3 nodes per side, it devides the side into equal segments.
    numLoops = 20,
    segLen = 700 * SVector{7}(ones(7)),  # The hexagon's side length is 10, each segment is 10/3.
    slipSystem = 1,
    _slipPlane = slipSystems.slipPlane[:, 1],
    _bVec = slipSystems.bVec[:, 1],
    label = SVector{7,nodeType}(1, 1, 1, 1, 1, 2, 1),
    buffer = 0.0,
    range = SMatrix{3,2,Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
    dist = Rand(),
)

network = DislocationNetwork((prismHeptagon, prismPentagon))

##
plotlyjs()
fig1 = plotNodes(
    network;
    m = 1,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,
)

@allocated plotNodes(
    network;
    m = 1,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,
)

@btime plotNodes(
    network;
    m = 1,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,
)

fig2 = plotNodes(
    prismHeptagon;
    m = 1,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,
)

##
dx = 2000
dy = 2000
dz = 2000
numNode = 100
numSeg = numNode - 1
links = zeros(Int, 2, numSeg)
slipPlane = zeros(3, numSeg)
bVec = zeros(3, numNode - 1)
coord = zeros(3, numNode)
label = zeros(nodeType, numNode)

coord[1, :] .= dx / 8
coord[2, :] .= dy / 2
coord[3, :] = range(0, dz, length = numNode)
b = Float64[1; 1; 1]
n = Float64[-1; 1; 0]

for i = 1:(numSeg - 1)
    links[:, i] .= (i, i + 1)
    bVec[:, i] = b
    slipPlane[:, i] = n
end
links = hcat(links, [numNode + 1, 1], [numNode + 2, numNode])
bVec = hcat(bVec, b, b)
slipPlane = hcat(slipPlane, n, n)
coord = hcat(coord, [0; 0; -1e3 * dz], [0; 0; 1e3 * dz])
append!(label, nodeType[3, 3])

test = DislocationNetwork(
    links,
    slipPlane,
    bVec,
    coord,
    label,
    zeros(size(coord)),
    zeros(size(coord)),
    numNode + 2,
    numSeg + 2,
)
makeConnect!(test)
getSegmentIdx!(test)

test2 = DislocationNetwork(
    links,
    slipPlane,
    bVec,
    coord,
    label,
    zeros(size(coord)),
    zeros(size(coord)),
    [numNode + 2],
    [numSeg + 2],
)
makeConnect!(test2)
getSegmentIdx!(test2)

compStruct(test, test2, verbose = true)

test.nodeVel
plotlyjs()
fig1 = plotNodes(
    test;
    m = 3,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,
    xlims = (0, dx),
    ylims = (0, dy),
    zlims = (0, dz),
)

##
shearDecagon = DislocationLoop(;
    loopType = loopShear(),
    numSides = 10,
    nodeSide = 1,
    numLoops = 1,
    segLen = [300; 700; 1100; 1500; 1900; 1900; 1500; 1100; 700; 300],
    slipSystem = 4,
    _slipPlane = slipSystems.slipPlane[:, 4],
    _bVec = slipSystems.bVec[:, 4],
    label = nodeType[1; 1; 1; 1; 1; 1; 1; 1; 1; 1],
    buffer = 0.0,
    range = Float64[0 0; 0 0; 0 0],
    dist = Zeros(),
)

network = DislocationNetwork(shearDecagon, memBuffer = 1)

network = DislocationNetwork(shearDecagon)
network.coord[:, 11] = vec(mean(network.coord, dims = 2))
network.label[11] = 1
network.numNode[1] = 11
network.links[:, 11] = [11; 2]
network.links[:, 12] = [11; 4]
network.links[:, 13] = [11; 5]
network.links[:, 14] = [11; 10]
network.bVec[:, 11:14] .= network.bVec[:, 1]
network.slipPlane[:, 11:14] .= network.slipPlane[:, 1]
makeConnect!(network)
getSegmentIdx!(network)

network.label
network.numNode[1]
network.numSeg[1]

@btime refineNetwork!(dlnParams, matParams, network)
@btime coarsenNetwork!(dlnParams, matParams, network)

refineNetwork!(dlnParams, matParams, network)

@allocated refineNetwork!(dlnParams, matParams, network)

@code_warntype refineNetwork!(dlnParams, matParams, network)

network.label
network.numNode[1]
network.numSeg[1]

##
plotlyjs()
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)

using JSON3, StructTypes, FileIO
StructTypes.StructType(::Type{<:DislocationNetwork}) = StructTypes.Struct()
StructTypes.StructType(::Type{nodeType}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{nodeType}) = Int
open("test.json3", "w") do io
    return JSON3.write(io, network)
end
newvar = open("test.json3", "r") do io
    return JSON3.read(io)
end

links = reshape(newvar["links"], 2, :)
slipPlane = reshape(newvar["slipPlane"], 3, :)
bVec = reshape(newvar["bVec"], 3, :)
coord = reshape(newvar["coord"], 3, :)
label = nodeType.(newvar["label"])
nodeVel = reshape(newvar["nodeVel"], 3, :)
nodeForce = reshape(newvar["nodeForce"], 3, :)
numNode = newvar["numNode"]
numSeg = newvar["numSeg"]
maxConnect = newvar["maxConnect"]
connectivity = reshape(newvar["connectivity"], 1 + 2 * newvar["maxConnect"], :)
linksConnect = reshape(newvar["linksConnect"], 2, :)
segIdx = reshape(newvar["segIdx"], :, 3)
segForce = reshape(newvar["segForce"], 3, 2, :)

network2 = DislocationNetwork(;
    links = copy(reshape(newvar["links"], 2, :)),
    slipPlane = copy(reshape(newvar["slipPlane"], 3, :)),
    bVec = copy(reshape(newvar["bVec"], 3, :)),
    coord = copy(reshape(newvar["coord"], 3, :)),
    label = copy(nodeType.(newvar["label"])),
    nodeVel = copy(Float64.(reshape(newvar["nodeVel"], 3, :))),
    nodeForce = copy(Float64.(reshape(newvar["nodeForce"], 3, :))),
    numNode = copy(newvar["numNode"]),
    numSeg = copy(newvar["numSeg"]),
    maxConnect = copy(newvar["maxConnect"]),
    connectivity = copy(reshape(newvar["connectivity"], 1 + 2 * newvar["maxConnect"], :)),
    linksConnect = copy(reshape(newvar["linksConnect"], 2, :)),
    segIdx = copy(reshape(newvar["segIdx"], :, 3)),
    segForce = copy(Float64.(reshape(newvar["segForce"], 3, 2, :))),
)

compStruct(network, network2)
FileIO.save("test.jld2", "network", network, "shearDecagon", [shearDecagon, shearDecagon])
a, b = FileIO.load("test.jld2", "network", "shearDecagon")
b
network3, shearDecagon2 = FileIO.loadJSON("test.jld2", "network", "shearDecagon")
compStruct(network, network3)
compStruct(shearDecagon, shearDecagon2)

network.segForce

abstract type Vehicle end

struct Car <: Vehicle
    type::String
    make::String
    model::String
    seatingCapacity::Int
    topSpeed::Float64
end

struct Truck <: Vehicle
    type::String
    make::String
    model::String
    payloadCapacity::Float64
end

StructTypes.StructType(::Type{Vehicle}) = StructTypes.AbstractType()
StructTypes.StructType(::Type{Car}) = StructTypes.Struct()
StructTypes.StructType(::Type{Truck}) = StructTypes.Struct()
StructTypes.subtypekey(::Type{Vehicle}) = :type
StructTypes.subtypes(::Type{Vehicle}) = (car = Car, truck = Truck)

car = JSON3.read(
    """
{
    "type": "car",
    "make": "Mercedes-Benz",
    "model": "S500",
    "seatingCapacity": 5,
    "topSpeed": 250.1
}""",
    Vehicle,
)

using Serialization, BSON

@time serialize(
    "test2.jls",
    (network, dlnParams, matParams, intParams, slipSystems, dislocationLoop),
)
@time saveJSON(
    "test3.json",
    (network, dlnParams, matParams, intParams, slipSystems, dislocationLoop),
)
@time bson(
    "test4.bson",
    a = (network, dlnParams, matParams, intParams, slipSystems, dislocationLoop),
)

networkOUT1, dlnParamsOUT1, matParamsOUT1 = open("test.txt", "r") do io
    return deserialize(io)
end

networkOUT1
dlnParamsOUT1

networkOUT2 = deserialize("test2.txt")

compStruct(networkOUT1, network)
compStruct(networkOUT2, network)

## Remeshing and integration
network2 = deepcopy(network)
intVars2 = deepcopy(intVars)
numSeg = network.numNodeSegConnect[2]
intVars2
intParams
network2.coord[:, 1:numSeg] - network.coord[:, 1:numSeg]
fig2 = plotNodes(
    network2,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)

function foo(dlnParams, matParams, network)
    return network = refineNetwork!(dlnParams, matParams, network)
end
function bar(dlnParams, matParams, network)
    return network = coarsenNetwork!(dlnParams, matParams, network)
end
foo(dlnParams, matParams, network)
bar(dlnParams, matParams, network)

@btime foo(dlnParams, matParams, network)
@btime bar(dlnParams, matParams, network)

numSeg
network.numNodeSegConnect[2]

function foo(intParams, intVars, dlnParams, matParams, network)
    network = coarsenNetwork!(dlnParams, matParams, network)
    network = refineNetwork!(dlnParams, matParams, network)
    return integrate!(intParams, intVars, dlnParams, matParams, network)
end
network2 = deepcopy(network)
@allocated foo(intParams, intVars, dlnParams, matParams, network2)
@allocated foo(intParams, intVars, dlnParams, matParams, network2)

@allocated foo(intParams, intVars, dlnParams, matParams, network2)
@allocated foo(intParams, intVars, dlnParams, matParams, network2)
gr()
function baar(intParams, intVars, dlnParams, matParams, network)

    network2 = deepcopy(network)
    intVars2 = deepcopy(intVars)

    anim = @animate for i = 1:500
        fig = plotNodes(
            network2,
            m = 3,
            l = 2,
            linecolor = :blue,
            markercolor = :blue,
            legend = false,
            # camera=(60,30),
        )
        foo(intParams, intVars2, dlnParams, matParams, network2)
        # if mod(i, 10) == 0
        # plotNodes!(
        #     fig,
        #     network2,
        #     m = 1,
        #     l = 3,
        #     linecolor = :blue,
        #     markercolor = :blue,
        #     legend = false,
        #     show=true
        # )
        # plot!(fig)
        # end
        # println(network2.numNode)
    end every 10

    return gif(anim, "uau4.gif")
end

baar(intParams, intVars, dlnParams, matParams, network)

# for i in 1:1000
#     integrate!(intParams, intVars2, dlnParams, matParams, network2)
#     plotNodes!(
#         fig2,
#         network2,
#         m = 1,
#         l = 3,
#         linecolor = :blue,
#         markercolor = :blue,
#         legend = false,
#     )
# end
network2 = deepcopy(network)
@btime network3 = coarsenNetwork!(dlnParams, matParams, network2)
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)
network = refineNetwork!(dlnParams, matParams, network)
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)
network = refineNetwork!(dlnParams, matParams, network)
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)

calcSegSegForce(dlnParams, matParams, network)

network.numSeg[1]
network.label

network = DislocationNetwork!(network, [prismHeptagon, prismPentagon])
prismHeptagon

network.label
network.numSeg[1]
size(network.connectivity)
size(network.links)
findall(x -> x != 0, prismHeptagon.links[1, :])
gr()
fig = plotNodes(
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
    # camera=(60,30),
    # size=(400,400)
)

dlnParamsPar = DislocationParameters(;
    coreRad = dlnParams.coreRad,
    coreRadMag = dlnParams.coreRadMag,
    minSegLen = dlnParams.minSegLen,
    maxSegLen = dlnParams.maxSegLen,
    minArea = dlnParams.minArea,
    maxArea = dlnParams.maxArea,
    maxConnect = dlnParams.maxConnect,
    remesh = dlnParams.remesh,
    collision = dlnParams.collision,
    separation = dlnParams.separation,
    virtualRemesh = dlnParams.virtualRemesh,
    parCPU = true,
    parGPU = dlnParams.parGPU,
    edgeDrag = dlnParams.edgeDrag,
    screwDrag = dlnParams.screwDrag,
    climbDrag = dlnParams.climbDrag,
    lineDrag = dlnParams.lineDrag,
    mobility = dlnParams.mobility,
)

remoteForceSer = calcSegSegForce(dlnParams, matParams, network)
remoteForcePar = calcSegSegForce(dlnParamsPar, matParams, network)
isapprox(remoteForceSer, remoteForcePar)
calcSegSegForce!(dlnParams, matParams, network)
isapprox(remoteForceSer, network.segForce[:, :, 1:(network.numNodeSegConnect[2])])
network.segForce .= 0
calcSegSegForce!(dlnParamsPar, matParams, network)
isapprox(remoteForcePar, network.segForce[:, :, 1:(network.numNodeSegConnect[2])])
network.segForce .= 0

@btime calcSegSegForce(dlnParams, matParams, network)
@btime calcSegSegForce!(dlnParams, matParams, network)
network.segForce .= 0
@btime calcSegSegForce(dlnParamsPar, matParams, network)
@btime calcSegSegForce!(dlnParamsPar, matParams, network)
network.segForce .= 0

@code_warntype calcSegForce(dlnParamsPar, matParams, network)
using StaticArrays
@code_warntype calcSegSegForce(
    0.5,
    0.2,
    0.55,
    0.70,
    0.13,
    0.4,
    SVector(0, 2, 3),
    SVector(1, -5, 6),
    SVector(1, 1, 1),
    SVector(-1, 5, 2),
    SVector(1, 2, 3),
    SVector(1, 2, -3),
)

network.segForce[:, 1, 1] = [1, 2, 3]

network.segForce
remoteForcePar
remoteForceSer
baar(intParams, intVars, dlnParams, matParams, network)

##
prismPentagonSlow = DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = 500 * ones(5),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 2, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 2],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 2],            # Burgers vector of the segments.
    label = nodeType[1; 1; 1; 1; 1],    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = Float64[          # Distribution range
        -5000 5000 # xmin, xmax
        -5000 5000 # ymin, ymax
        -5000 5000  # zmin, zmax
    ],
    dist = Rand(),  # Loop distribution.
)
prismHeptagonSlow = DislocationLoop(;
    loopType = loopPrism(),    # Shear loop
    numSides = 7,
    nodeSide = 1,   # 3 nodes per side, it devides the side into equal segments.
    numLoops = 20,
    segLen = 700 * ones(7),  # The hexagon's side length is 10, each segment is 10/3.
    slipSystem = 1,
    _slipPlane = slipSystems.slipPlane[:, 1],
    _bVec = slipSystems.bVec[:, 1],
    label = nodeType[1; 1; 1; 1; 1; 2; 1],
    buffer = 0.0,
    range = Float64[
        -5000 5000
        -5000 5000
        -5000 5000
    ],
    dist = Rand(),
)

prismPentagonFast = DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = 500 * SVector{5}(ones(5)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 2, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 2],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 2],            # Burgers vector of the segments.
    label = SVector{5,nodeType}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
    dist = Rand(),  # Loop distribution.
)
prismHeptagonFast = DislocationLoop(;
    loopType = loopPrism(),    # Shear loop
    numSides = 7,
    nodeSide = 1,   # 3 nodes per side, it devides the side into equal segments.
    numLoops = 20,
    segLen = 700 * SVector{7}(ones(7)),  # The hexagon's side length is 10, each segment is 10/3.
    slipSystem = 1,
    _slipPlane = slipSystems.slipPlane[:, 1],
    _bVec = slipSystems.bVec[:, 1],
    label = SVector{7,nodeType}(1, 1, 1, 1, 1, 2, 1),
    buffer = 0.0,
    range = SMatrix{3,2,Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
    dist = Rand(),
)

@btime DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = 500 * SVector{5}(ones(5)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 2, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 2],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 2],            # Burgers vector of the segments.
    label = SVector{5,nodeType}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
    dist = Rand(),  # Loop distribution.
)

@btime DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = 500 * ones(5),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 2, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 2],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 2],            # Burgers vector of the segments.
    label = nodeType[1; 1; 1; 1; 1],    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = Float64[          # Distribution range
        -5000 5000 # xmin, xmax
        -5000 5000 # ymin, ymax
        -5000 5000  # zmin, zmax
    ],
    dist = Rand(),  # Loop distribution.
)

sourcesFast = (prismHeptagonFast, prismPentagonFast)
sourcesSlow = (prismHeptagonSlow, prismPentagonSlow)

@btime DislocationNetwork(sourcesFast)
@btime DislocationNetwork(sourcesSlow)

@btime DislocationNetwork(prismHeptagonSlow)
@btime DislocationNetwork(prismHeptagonFast)
@btime DislocationNetwork(prismPentagonSlow)
@btime DislocationNetwork(prismPentagonFast)

@btime DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = 500 * ones(5),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 2, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 2],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 2],            # Burgers vector of the segments.
    label = nodeType[1; 1; 1; 1; 1],    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = Float64[          # Distribution range
        -5000 5000 # xmin, xmax
        -5000 5000 # ymin, ymax
        -5000 5000  # zmin, zmax
    ],
    dist = Rand(),  # Loop distribution.
)
@btime DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 5,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 20,  # Number of loops of this type to generate when making a network.
    segLen = 500 * SVector{5}(ones(5)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 2, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 2],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 2],            # Burgers vector of the segments.
    label = SVector{5,nodeType}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
    dist = Rand(),  # Loop distribution.
)

@btime DislocationLoop(;
    loopType = loopPrism(),    # Shear loop
    numSides = 7,
    nodeSide = 1,   # 3 nodes per side, it devides the side into equal segments.
    numLoops = 20,
    segLen = 700 * ones(7),  # The hexagon's side length is 10, each segment is 10/3.
    slipSystem = 1,
    _slipPlane = slipSystems.slipPlane[:, 1],
    _bVec = slipSystems.bVec[:, 1],
    label = nodeType[1; 1; 1; 1; 1; 2; 1],
    buffer = 0.0,
    range = Float64[
        -5000 5000
        -5000 5000
        -5000 5000
    ],
    dist = Rand(),
)
@btime DislocationLoop(;
    loopType = loopPrism(),    # Shear loop
    numSides = 7,
    nodeSide = 1,   # 3 nodes per side, it devides the side into equal segments.
    numLoops = 20,
    segLen = 700 * SVector{7}(ones(7)),  # The hexagon's side length is 10, each segment is 10/3.
    slipSystem = 1,
    _slipPlane = slipSystems.slipPlane[:, 1],
    _bVec = slipSystems.bVec[:, 1],
    label = SVector{7,nodeType}(1, 1, 1, 1, 1, 2, 1),
    buffer = 0.0,
    range = SMatrix{3,2,Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
    dist = Rand(),
)

##
@btime DislocationNetwork((prismHeptagon, prismPentagon))
@btime DislocationNetwork!(network, (prismHeptagon, prismPentagon))
@btime calcSegForce(dlnParams, matParams, network)
@btime calcSegForce!(dlnParams, matParams, network)
