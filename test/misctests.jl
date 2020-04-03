using DDD
using Test
cd(@__DIR__)

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

params = "../inputs/simParams/sampleParams.csv"
slipsys = "../data/slipSystems/bcc.csv"
source = "../inputs/dln/sampleDln.csv"
using LinearAlgebra
dlnParams, matParams, intParams, slipSystems, loops =
    loadParams(params, slipsys, source)
network = makeNetwork(loops; memBuffer = 1)
segVec = getSegVector(network)

sqrt.(sum(segVec.*segVec, dims=2))



segVec[:,:].*segVec[:,:]

sum(segVec[1,:].*segVec[1,:],dims=1)
norm(segVec,2)


scene1 = plotNodesMakie(network, linewidth = 2, markersize = 0.1, strokecolor =:green, color=:orange)
plotNodesMakie!(scene1, network2, linewidth = 2, markersize = 0.1, strokecolor =:blue, color=:blue)

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
