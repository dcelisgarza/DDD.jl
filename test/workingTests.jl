using Revise, BenchmarkTools, Plots, LinearAlgebra
# using plotlyjs
# https://github.com/sglyon/ORCA.jl/issues/8#issuecomment-629049679
# https://stackoverflow.com/a/48509426

using DDD
cd(@__DIR__)
plotlyjs()

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
# network = DislocationNetwork(dislocationLoop, memBuffer = 1)
# DislocationNetwork!(network, dislocationLoop)

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
    return network = refineNetwork(dlnParams, matParams, network)
end
function bar(dlnParams, matParams, network)
    return network = coarsenNetwork(dlnParams, matParams, network)
end
foo(dlnParams, matParams, network)
bar(dlnParams, matParams, network)

@btime foo(dlnParams, matParams, network)
@btime bar(dlnParams, matParams, network)

numSeg
network.numNodeSegConnect[2]

function foo(intParams, intVars, dlnParams, matParams, network)
    coarsenNetwork(dlnParams, matParams, network)
    refineNetwork(dlnParams, matParams, network)
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
@btime network3 = coarsenNetwork(dlnParams, matParams, network2)
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)
refineNetwork(dlnParams, matParams, network)
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)
refineNetwork(dlnParams, matParams, network)
fig1 =
    plotNodes(network, m = 1, l = 3, linecolor = :blue, markercolor = :blue, legend = false)

prisPentagon = DislocationLoop(
    loopPrism();    # Prismatic loop, all segments are edge segments.
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

prismHeptagon = DislocationLoop(
    loopPrism();    # Shear loop
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

using Random
Random.seed!(1337)
network = DislocationNetwork(
    [prismHeptagon, prisPentagon]; # Dispatch type, bespoke functions dispatch on this.
    memBuffer = 1, # Buffer for memory allocation.
)

@time for _ in 1:50
    global network = DislocationNetwork!(
        network,
        [prismHeptagon, prisPentagon]; # Dispatch type, bespoke functions dispatch on this.
        memBuffer = 1, # Buffer for memory allocation.
    )
end

network.numNodeSegConnect
size(network.connectivity)
size(network.links)

network.coord[:, (network.numNodeSegConnect[2] - 1):(network.numNodeSegConnect[2] + 1)]
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

struct test{T1, T2, T3}
    a::T1
    b::T2
    c::T3
end

var = test(1, 2.5, Float64[1 2 3; 4 5 6])

var
var.c[:] = Float64[5, 6, 7, 9]

var.a .= 2
var.c[1] = -5
var.c[2:3] = [-10, 25]

var.c
var.c .= vcat(var.c, [7 8 9])

var.c .= resize!(reshape(var.c, :), 9)
var.c

var.c
var.c .= [-1, -1, -1, -1]

var.c

resize!(var.c, 2)
var.c

mutable struct mutate_me
    a::Array{Int, 2}
    b::Vector{Int}
end

struct immutate_me
    a::Array{Int, 2}
    b::Vector{Int}
end

function foo!(variable, condition)
    if condition
        # variable.a = vcat(variable.a, zeros(Int, size(variable.a)))
        push!(vec(variable.a), vec(zeros(Int, size(variable.a))))
    end
    variable.a .= LinearIndices(variable.a)
    return nothing
end

function foo(variable, condition)
    if condition
        # variable = immutate_me(vcat(variable.a, zeros(Int, size(variable.a))))
        push!(vec(variable.a), vec(zeros(Int, size(variable.a))))
    end
    variable.a .= LinearIndices(variable.a)
    return variable
end

mutating_var = mutate_me([0 0 0; 0 0 0], [0, 0, 0])
immutating_var = immutate_me([0 0 0; 0 0 0], [0, 0, 0])

# Always works
foo!(mutating_var, rand(Bool))
mutating_var

# Works when no resizing is needed, obviously
foo!(immutating_var, rand(Bool))
immutating_var

# immutating_var changes when no resizing is needed, does't work otherwise
foo(immutating_var, rand(Bool))

# immutating_var always changes
immutating_var = foo(immutating_var, rand(Bool))

function watanabe(var)
    append!(var.b, zeros(Int, 2 * length(var.b)))
    var.b .= LinearIndices(var.b)
    return nothing
end

mutating_var = mutate_me([0 0 0; 0 0 0], [0, 0, 0])
immutating_var = immutate_me([0 0 0; 0 0 0], [0, 0, 0])
watanabe(mutating_var)
watanabe(immutating_var)
mutating_var
immutating_var

prod(size(mutating_var.a))

# Methods to implement for linear arrays
# https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array-1
struct LinearArray{T, N} <: AbstractArray{T, N} end
struct VectorisedArray{T, N} <: AbstractArray{T, N}
    size::NTuple{N, Int}
    data::Vector{T}
end

Base.size(V::VectorisedArray) = V.size
Base.IndexStyle(::Type{<:VectorisedArray}) = IndexLinear()
Base.getindex(V::VectorisedArray, i::Int) = V.data[i]
function Base.getindex(V::VectorisedArray, I::Vararg{Int, N}) where {N}
    idx = I[1]

    for i in 2:N
        idx += prod(V.size[1:(i - 1)]) * (I[i] - 1)
    end

    return V.data[idx]
end
function Base.setindex!(V::VectorisedArray, v, i::Int)
    return V.data[i] = v
end
function Base.setindex!(V::VectorisedArray, v, I::Vararg{Int, N}) where {N}
    return V.data[I] = v
end
function Base.hcat(V::VectorisedArray, i::Int)
    size = zeros(Int, length(V.size))
    size .= V.size
    size[end] += i
    append!(V.data, zeros(eltype(V.data), prod(V.size) * i))
    V = VectorisedArray(Tuple(size), V.data)
    return V
end

a = [5, 6]

b = (7, 8)
c = zeros(Int, length(b))

c .= b
vec(b)
Tuple(a)
dim = (3, 2)
I = (2, 2)
data = [1, 2, 3, 4, 5, 6]

data[I[1] + dim[2] * (I[2] - 1)]

I[1] + dim[2] * (I[2] - 1) + dim[2] * dim[3] * (I[3] - 1)

println(wak...)

sum(wak[end:-1:1] .* (test[1:end] .- 1))
sum(I[end:-1:1] * (V.dim[1:end] - 1))

VectorisedArray((4, 2), [1, 2, 3, 4, 5, 6, 7, 8])
size([1 2])
test = VectorisedArray((4, 2), [1, 2, 3, 4, 5, 6, 7, 8])

test = push!(test, 5)

test[2, 1]

test = (3, 2)
test[1]

VectorisedArray{T, N}(init, I...) where {T, N} = 1

test = (5, 6)
println(test...)

var = LinearArray{Int, 2}(undef, 2, 5)
size(var)

struct foo{T1, T2}
    a::T1
    b::T2
end

function bar(foo)

    a = foo.a
    b = foo.b

    a .= -1
    a[1, 1:2] = [50, -50]
    b .-= 30
    return nothing
end

test = foo([4 -5; 5 2; 6 57], [5])
test

bar(test)

test

test.a

test.a .= 0

test
