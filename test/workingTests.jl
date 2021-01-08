##
using BenchmarkTools, Plots, LinearAlgebra, StaticArrays
using DDD

cd(@__DIR__)
fileDislocationParameters = "../inputs/simParams/sampleDislocationParameters.json"
fileMaterialParameters = "../inputs/simParams/sampleMaterialParameters.json"
fileIntegrationParameters = "../inputs/simParams/sampleIntegrationParameters.json"
fileSlipSystem = "../data/slipSystems/BCC.json"
fileDislocationLoop = "../inputs/dln/samplePrismShear.json"
fileIntVar = "../inputs/simParams/sampleIntegrationTime.json"
dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParametersJSON(
    fileDislocationParameters,
    fileMaterialParameters,
    fileIntegrationParameters,
    fileSlipSystem,
    fileDislocationLoop,
)
intVars = loadIntegrationTimeJSON(fileIntVar)

network = DislocationNetwork(dislocationLoop)
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
        for k in 1:3
            flag = flag && isapprox(a[k, i], b[k, j])
        end
        flag
    end
    equalSlipPlane ? for i in 1:3
        a[i, i2] = b[i, j2]
    end : nothing
end
comparing3(a, b, i, j, i2, j2)

[i for i in 1:3]

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
    @inbounds @simd for i in 500:5000
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
wat = SVector{3, Float64}(3, 3, 3)

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
    label = SVector{5, nodeType}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3, 2, Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
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
    label = SVector{7, nodeType}(1, 1, 1, 1, 1, 2, 1),
    buffer = 0.0,
    range = SMatrix{3, 2, Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
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

for i in 1:(numSeg - 1)
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

    anim = @animate for i in 1:500
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
    label = SVector{5, nodeType}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3, 2, Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
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
    label = SVector{7, nodeType}(1, 1, 1, 1, 1, 2, 1),
    buffer = 0.0,
    range = SMatrix{3, 2, Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
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
    label = SVector{5, nodeType}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3, 2, Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
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
    label = SVector{5, nodeType}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0.0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3, 2, Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
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
    label = SVector{7, nodeType}(1, 1, 1, 1, 1, 2, 1),
    buffer = 0.0,
    range = SMatrix{3, 2, Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
    dist = Rand(),
)

##
@btime DislocationNetwork((prismHeptagon, prismPentagon))
@btime DislocationNetwork!(network, (prismHeptagon, prismPentagon))
@btime calcSegForce(dlnParams, matParams, network)
@btime calcSegForce!(dlnParams, matParams, network)
