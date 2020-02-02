include("../src/GlobalModule.jl")
using BenchmarkTools
using .GlobalModule.CustomTypes
using .GlobalModule.DislocationBase
using .GlobalModule.MaterialBase
using .GlobalModule.CustomIntegration
using .GlobalModule.DdFemBase

"""
Author: Daniel Celis Garza
Date: 2020/01/23

Testing MaterialPParams module.
"""

material = MaterialP(1.0, 1e5, 0.28, "bcc")
dislocation = DislocationP(
    0.5,
    1e5,
    5,
    6.0,
    50.3,
    75.3,
    4,
    true,
    true,
    true,
    true,
    1.0,
    2.0,
    3.0,
    4.0,
    :bcc,
)
elem = ones(3)
mesh = CuboidMesh(elem)
tspan = zeros(2)
tspan = [0.0, 1000.0]

# integrator = IntegrationP(0.3, tspan, CustomTrapezoid)

println(dislocation)
println(fieldnames(typeof(dislocation)))

n_dln = 0
n_seg = 0
links = zeros(Integer, 500, 2)
b_vec = zeros(500, 3)
n_vec = zeros(500, 3)
coord = zeros(500, 3)
label = zeros(Integer, 500)
network = DislocationNetwork(links, b_vec, n_vec, coord, label, n_dln, n_seg)
n_dln = nothing
n_seg = nothing
links = nothing
b_vec = nothing
n_vec = nothing

label[1:2:end] .= 5

idx1 = getIndex(network, :label, ==, 5)
data = getData(network, :coord, :coord, 3, ==, 0)
# Implement the array conditions
test = [1 2 3; -4 2 6; 1 5 6; 5 2 9]
index = findall(x -> x .== 5, test)
lindex = LinearIndices(test)[index]




# idx2 = getIndex(network, 0)
# data1 = coord[findall(x -> x == 0, label), :]
# data2 = getCoord(network, 0)
# data3 = coord[getIndex(network, :label, ==, 0), :]
# data4 = coord[getIndex(network, 0), :]
# data5 = coord[idx1, :]
# println(data1 == data2 == data3 == data4)
# println(idx1 == idx2)
# println(data2 == coord[idx1, :])
# println(data2 == coord[idx2, :])

# println("Index")
# @benchmark findall(x -> x == 5, network.label)
# @benchmark getIndex(network, :label, ==, 0)
# @benchmark getIndex(network, 0)
#
# println("Data")
# @benchmark coord[findall(x -> x == 0, label), :]
# @benchmark getCoord(network, 0)
# @benchmark coord[getIndex(network, :label, ==, 0), :]
# @benchmark coord[getIndex(network, 0), :]
# @benchmark coord[idx1, :]
