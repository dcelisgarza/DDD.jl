using Revise
using DDD
cd(@__DIR__)




slipSystem = load("../inputs/simParams/SlipSystems.JSON")
dislocation = load("../inputs/simParams/sampleDislocation.JSON")
slipSystem = loadSlipSystem(slipSystem[1])
loadDislocationLoop(
    dislocation[1],
    slipSystem,
)
























params = "../inputs/simParams/sampleParams.csv"
slipsys = "../data/slipSystems/bcc.csv"
source = "../inputs/dln/samplePrismaticShear.csv"
using LinearAlgebra
dlnParams, matParams, intParams, slipSystems, loops =
    loadParams(params, slipsys, source)
network = makeNetwork(loops; memBuffer = 2)
makeNetwork!(network, loops)

# self = calcSelfForce(dlnParams, matParams, network)
# @btime calcSelfForce(dlnParams, matParams, network)
# par = calcSegSegForce(dlnParams, matParams, network; parallel = true)
# @btime calcSegSegForce(dlnParams, matParams, network; parallel = true)

output = "../outputs/dln/sampleDln.JSON"
# , intParams, slipSystems, loops
save(output, dlnParams, matParams, intParams, slipSystems, loops)
dict = load(output)

zip(dict[4][:])

compStruct(intParams, intParams2)
eval(dict[1]["mobility"][1:(end - 2)])

dict[end][1] = translateEnum(nodeType, dict[end][1], "label")
dict[end][2] = translateEnum(nodeType, dict[end][2], "label")

replaceslip = readdlm(slipsys, ',')

ss = SlipSystem(
    name = "BCC",
    slipPlane = replaceslip[:, 1:3],
    bVec = replaceslip[:, 4:6],
)

"../data/slipSystems/SlipSystems.JSON"
save("../data/slipSystems/SlipSystems.JSON", ss)

# import JSON: lower
# # JSON.lower(t::T) where {T<:AbstractCrystalStruct} = string(t)
# # JSON.lower(t::T) where {T<:AbstractMobility} = string(t)
# # JSON.lower(t::T) where {T<:AbstractIntegrator} = string(t)
# # JSON.lower(t::T) where {T<:AbstractDlnSeg} = string(t)
# # JSON.lower(t::T) where {T<:AbstractDlnStr} = string(t)
# # JSON.lower(t::T) where {T<:AbstractDistribution} = string(t)
# JSON.lower(t::T) where {T<:Union{AbstractCrystalStruct, AbstractMobility, AbstractIntegrator, AbstractDlnSeg, AbstractDlnStr, AbstractDistribution}} = string(t)

testing = "../outputs/dln/test3.JSON"
open(testing, "w") do io
    JSON.print(io, loops)
end
JSON.parsefile(testing)

for i = 1:length(in["label"])
    in["label"][i] = labels[in["label"][i]]
end

save(output, network, dlnParams, matParams)
load(output)

insts = instances(nodeType)
labels = Dict{String, nodeType}()
for inst in insts
    push!(labels, string(inst) => inst)
end
labels

labels = dict["label"]
first = labels
string(insts[1])

parse(first[1])

f(x::nodeType) = "I'm a Fruit with value: $(Int(x))"
f(intMob)

function subtypetree(t, level = 1, dict = Dict())
    push!(dict, t => supertype(t))
    for s in subtypes(t)
        subtypetree(s, level + 1, dict)
    end
    return dict
end
abstract type test <: AbstractShapeFunction2D end
wakanda = subTypeTree(AbstractShapeFunction)

wakanda = makeTypeDict(AbstractShapeFunction; cutoff = 3)
for (key, val) in wakanda
    println(key)
end
supertype(test)
dict

inst = instances(nodeType)
inst[1]

lbl = dict["label"]
label1 = lbl[1]
eval(:label1)

var = nodeType(instances(lbl[1]))
var.value

names = fieldnames(DislocationNetwork)
df = DataFrame(var = Any[], val = Any[])
push!(df, (:links, getfield(network, :links)))
for name in names
    push!(df, (name, getfield(network, name)))
end
df[4, 2]

readdlm(output)
fieldnames(DislocationNetwork)

saveNetwork(network, output)
size(network.segIdx)
data = CSV.read(output; delim = ',')
string = chop(data[4, 2], head = 1, tail = 1)
pieces = split(string, ";")
map(pieces) do piece
    pieces2 = split(string, " ")
    map(pieces2) do piece2
        parse(Float64, piece2)
    end
end

coord == network.coord
import Base: getindex, iterate, length

@btime calcSelfForce(dlnParams, matParams, network)
Juno.@profiler for i = 1:10000
    calcSelfForce(dlnParams, matParams, network)
end

# @benchmark calcSelfForce(dlnParams, matParams, network)
@btime calcSegSegForce(dlnParams, matParams, network; parallel = true)
dsp1 = serial - par1
maximum(dsp1)
minimum(dsp1)

serial = calcSegSegForce(dlnParams, matParams, network; parallel = false)
@btime calcSegSegForce(dlnParams, matParams, network; parallel = false)

network2 = makeNetwork(loops; memBuffer = 1)
par2 = calcSegSegForce(dlnParams, matParams, network1)
@time calcSegSegForce(dlnParams, matParams, network)
dsp2 = serial - par2
maximum(dsp2)
minimum(dsp2)

diff = par1 - par2
diff = segseg - segsegFT
maximum(diff)
minimum(diff)

segsegT = calcSegSegForce(dlnParams, matParams, network)

diffFT_T = segsegFT - segsegFT

mean(segseg)
maximum(segseg)
minimum(segseg)

Juno.@profiler (
    for i = 1:1000
        calcSegSegForce(dlnParams, matParams, network)
    end
)

isapprox(sum(segseg) - 2 * eps(Float64), 0)
Profile.clear()

for i = 1:10000
    @profile calcSelfForce(dlnParams, matParams, network)
end

@benchmark calcSelfForce(dlnParams, matParams, network)

@benchmark calcSelfForce(dlnParams, matParams, network)

using Makie
function plotNodesMakie(network::DislocationNetwork, args...; kw...)
    idx = idxLabel(network, -1; condition = !=)
    coord = network.coord
    fig = Scene()
    meshscatter!(coord, args...; kw...)
    for i in idx
        n1 = network.links[i, 1]
        n2 = network.links[i, 2]
        lines!(
            coord[[n1, n2], 1],
            coord[[n1, n2], 2],
            coord[[n1, n2], 3],
            args...;
            kw...,
        )
    end
    return fig
end

function plotNodesMakie!(fig, network::DislocationNetwork, args...; kw...)
    idx = idxLabel(network, -1; condition = !=)
    coord = network.coord
    meshscatter!(coord, args...; kw...)
    for i in idx
        n1 = network.links[i, 1]
        n2 = network.links[i, 2]
        lines!(
            coord[[n1, n2], 1],
            coord[[n1, n2], 2],
            coord[[n1, n2], 3],
            args...;
            kw...,
        )
    end
    return fig
end

segVec = getSegVector(network)

sqrt.(sum(segVec .* segVec, dims = 2))

segVec[:, :] .* segVec[:, :]

sum(segVec[1, :] .* segVec[1, :], dims = 1)
norm(segVec, 2)

scene1 = plotNodesMakie(
    network,
    linewidth = 2,
    markersize = 0.1,
    strokecolor = :green,
    color = :orange,
)
plotNodesMakie!(
    scene1,
    network2,
    linewidth = 2,
    markersize = 0.1,
    strokecolor = :blue,
    color = :blue,
)

trythis = rand(100000, 3)
trythis2 = rand(100000, 3)

dimDot(trythis, trythis)
@benchmark dimNorm(trythis; dims = 2)
@benchmark sqrt.(sum(trythis .* trythis, dims = 2))
#=
using Plots
params = "../inputs/simParams/sampleParams.csv"
slipsys = "../data/slipSystems/bcc.csv"
source = "../inputs/dln/samplePrismaticShear.csv"
dlnParams, matParams, intParams, slipSystems, loops =
    loadParams(params, slipsys, source)

network = DislocationNetwork(
    zeros(Int64, 3, 2),
    zeros(3, 3),
    zeros(3, 3),
    zeros(3, 3),
    zeros(nodeType, 3),
    0,
    0,
)
makeNetwork!(network, loops)

points = Float64[
    0 0 0
    1 0 0
    1 1 0
    0 1 0
    0 0 1
    1 0 1
    1 1 1
    0 1 1
    0.5 0.5 0.25
]
Nall = shapeFunction(
    LinearQuadrangle3D(),
    points[:, 1],
    points[:, 2],
    points[:, 3],
)

N1 = shapeFunction(
    LinearQuadrangle3D(),
    points[2, 1],
    points[2, 2],
    points[2, 3],
)
Nall[:,2] == N1
isapprox(Nall[:,2], N1)

dNdS = shapeFunctionDeriv(
    LinearQuadrangle3D(),
    points[:, 1],
    points[:, 2],
    points[:, 3],
)

isapprox.(sum(dNdS[:,1,1]), 0)
for i = 1:size(points,1)
    # for j = 1:3
        @test isapprox.(sum(dNdS[:,:,i], dim=2), 0)
    # end
end
sum(dNdS[:,:,1], dims=1)


isapprox.(sum(dNdS[:,:,1], dims=1), 0)


gr()
fig = plot()
plotNodes!(
    fig,
    network,
    m = 1,
    l = 3,
    linecolor = :black,
    markercolor = :black,
    legend = false,
)

fig = plot()
plot(
    loops[1].coord[:, 1],
    loops[1].coord[:, 2],
    loops[1].coord[:, 3],
    m = 3,
    l = 3,
)
using Statistics
mean(loops[1].coord, dims = 1)
plotNodes!(
    fig,
    loops[1],
    m = 1,
    l = 3,
    linecolor = :black,
    markercolor = :black,
    legend = false,
)
plotNodes!(
    fig,
    loops[2],
    m = 1,
    l = 3,
    linecolor = :red,
    markercolor = :red,
    legend = false,
)
plotNodes!(
    fig,
    loops[3],
    m = 1,
    l = 3,
    linecolor = :blue,
    markercolor = :blue,
    legend = false,
)
plotNodes!(
    fig,
    loops[4],
    m = 1,
    l = 3,
    linecolor = :orange,
    markercolor = :orange,
    legend = false,
)

network = DislocationNetwork(
    zeros(Int64, 0, 2),
    zeros(0, 3),
    zeros(0, 3),
    zeros(0, 3),
    zeros(nodeType, 0),
    0,
    0,
)
makeNetwork!(network, loop)
# connectivity, linksConnect = makeConnect(network, dlnParams)

include("../src/PlottingMakie.jl")

plotNodesMakie(network; markersize = 0.6, linewidth = 3)

fig = Scene()
plotNodesMakie!(fig, network; markersize = 0.6, linewidth = 3)

# dlnSegTypes = subtypes(AbstractDlnSeg)
function makeTypeDict(valType::DataType)
    subTypes = subtypes(valType)
    dict = Dict{String, Any}()

    for i in eachindex(subTypes)
        push!(dict, string(subTypes[i]()) => eval(subTypes[i]()))
    end

    return dict
end

dict = makeTypeDict(AbstractDlnSeg)
test = subtypes(AbstractDlnSeg)
push!(dict, string(test[1]()) => eval(test[1]()))

dlnTypes = Dict(
    "loopPrism()" => loopPrism(),
    "loopShear()" => loopShear(),
    "loopMixed()" => loopMixed(),
    "loopDln()" => loopDln(),
    "DDD.loopPrism()" => loopPrism(),
    "DDD.loopShear()" => loopShear(),
    "DDD.loopMixed()" => loopMixed(),
    "DDD.loopDln()" => loopDln(),
)

test = rand(3, 3)
var = zeros(6, length(test))
for i in eachindex(test)
    var[:, i] .= i
end
sxx = reshape(var[1, :], size(test))

network2 = makeNetwork(loops, 4, 1)
compStruct(network, network2; verbose = false)
compStruct(loops, loops; verbose = false)
test1 = MyStruct1(1)

xyz = ones(3);
uvw = [-0.5; 2.0; 1.0];
abc = zeros(3);

using BenchmarkTools
test = rand(30000, 3)
x = test[:, 1]
y = test[:, 2]
z = test[:, 3]
shapeFunction(LinearQuadrangle3D(), x, y, z)
shapeFunctionDeriv(LinearQuadrangle3D(), x, y, z)
@benchmark shapeFunction(LinearQuadrangle3D(), x, y, z)
@benchmark shapeFunctionDeriv(LinearQuadrangle3D(), x, y, z)



@benchmark shapeFunction(LinearQuadrangle3D(), x, y, z)
@benchmark shapeFunctionDeriv(LinearQuadrangle3D(), x, y, z)


@benchmark shapeFunction(
    LinearQuadrangle3D(),
    test[:, 1],
    test[:, 2],
    test[:, 3],
)

@benchmark shapeFunctionDeriv(LinearQuadrangle3D(), x, y, z)
@benchmark shapeFunctionDeriv(
    LinearQuadrangle3D(),
    test[:, 1],
    test[:, 2],
    test[:, 3],
)
=#
