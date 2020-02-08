using DDD
using Test, DelimitedFiles

cd(@__DIR__)
params = "../inputs/simParams/sampleParams.csv"
dlnParams, matParams, intParams = loadParams(params)
network = DislocationNetwork(
    zeros(Integer, 0, 2),
    zeros(0, 3),
    zeros(0, 3),
    zeros(0, 3),
    zeros(nodeType, 0),
    0,
    0,
)

slipsys = "../data/slipSystems/bcc.csv"
slipSystems = readdlm(slipsys, ',')
network = makeLoop!(loopPrism(), network, dlnParams, slipSystems)

# intAngle(6) / π * 180
#
# using BenchmarkTools, Traceur
#
# function travA(A::Matrix{<:Float64})
#     # A[:,1] .= 1
#     # A[:,2] .= 2
#     # A[:,3] .= 3
#     for j = 1:size(A,2)
#         for i = 1:size(A,1)
#             result = A[i,j]
#         end
#         # A[i,1] = 1
#         # A[i,2] = 2
#         # A[i,3] = 3
#     end
# end
#
#
# function travB(B::Matrix{<:Float64})
#     # B[1,:] .= 1
#     # B[2,:] .= 2
#     # B[3,:] .= 3
#     for i = 1:size(B,2)
#         for j = 1:size(B,1)
#             result = B[j,i]
#         end
#         # B[:,i] = 1:3
#         # B[2,i] = 2
#         # B[3,i] = 3
#     end
# end
#
# A = rand(1000000,3)
# B = rand(3,1000000)
# @benchmark travA(A)
# @benchmark travB(B)
#
#
# idx = rand(1:size(A,1), 10000)
# @benchmark A[idx,:]
# @benchmark B[:,idx]
# @benchmark B[:,idx]
