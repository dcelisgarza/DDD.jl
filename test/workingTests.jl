##
using DDD, StaticArrays, SparseArrays, LinearAlgebra
dlnParams = DislocationParameters(; mobility = mobBCC())
matParams = MaterialParameters(; crystalStruct = BCC())
femParams = FEMParameters(; 
                        type = DispatchRegularCuboidMesh(), 
                        order = LinearElement(), 
                        model = CantileverLoad(), 
                        dx = 1013.0, dy = 1987.0, dz = 2999.0,
                        mx = 3, my = 5, mz = 7
                    )
slipSystem = SlipSystem(; crystalStruct = BCC(), slipPlane = Float64[-1;1;0], bVec = Float64[1;1;1])
intParams = IntegrationParameters(; method = AdaptiveEulerTrapezoid())
intTime = IntegrationTime()


length(fieldnames(typeof(dlnParams)))
length(fieldnames(typeof(matParams)))
length(fieldnames(typeof(femParams)))
length(fieldnames(typeof(slipSystem)))
length(fieldnames(typeof(intParams)))

nodeTypeDln.(ones(Int, 8))
dx, dy, dz = femParams.dx, femParams.dy, femParams.dz
segLen = (dx + dy + dz) / 30
prismLoop = DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 8,   # 8-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 1,   # Number of loops of this type to generate when making a network.
    segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystemIdx = 1, # Slip System index (assuming slip systems are stored in a file, this is the index).
    slipSystem = slipSystem,  # Slip system.
    label = SVector{8,nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(dx / 2, dy / 2, dz / 2, dx / 2, dy / 2, dz / 2),  # Distribution range
    dist = Zeros(),  # Loop distribution.
)
shearLoop = DislocationLoop(;
    loopType = loopShear(),    # Prismatic loop, all segments are edge segments.
    numSides = 8,   # 8-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 10,   # Number of loops of this type to generate when making a network.
    segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystemIdx = 1, # Slip System index (assuming slip systems are stored in a file, this is the index).
    slipSystem = slipSystem,  # Slip system.
    label = SVector{8,nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(0, 0, 0, dx, dy, dz),  # Distribution range
    dist = Rand(),  # Loop distribution.
)
network = DislocationNetwork([prismLoop, shearLoop])
regularCuboidMesh = buildMesh(matParams, femParams)
cantileverBC, forceDisplacement = Boundaries(femParams, regularCuboidMesh)

cornerNode = regularCuboidMesh.cornerNode
edgeNode = regularCuboidMesh.edgeNode
faceNode = regularCuboidMesh.faceNode
uGamma = (type = nodeTypeFE(1),# Type
                idx = :x0y0z0, # Index
                node = cornerNode[:x0y0z0])
tGamma = (type = [nodeTypeFE(2)],# Type
                idx = :x_y0z1, # Index
                node = edgeNode[:x_y0z1]) 
mGamma = (type = [nodeTypeFE(3)],# Type
                idx = :xy_z0, # Index
                node = faceNode[:xy_z0])
testGamma, testForceDisp = Boundaries(femParams, regularCuboidMesh; uGamma = uGamma, tGamma = tGamma, mGamma = mGamma)
testGamma.tGamma
saveJSON("test.json", testGamma)
testGammaDict = loadJSON("test.json")
testGamma2 = loadBoundaries(testGammaDict)
for i in fieldnames(typeof(testGamma2))
    isnothing(getproperty(testGamma2, i)) ? continue : nothing
    println(isequal(getproperty(testGamma2, i), getproperty(testGamma, i)))
end



testGamma.uGamma

fig = plotNodes(
    regularCuboidMesh,
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,
    dpi = 480,
)
savefig("remeshed.pdf")
plotNodes!(
    fig,
    shearLoop,
    m = 1,
    l = 3,
    linecolor = :red,
    marker = :star,
    markercolor = :red,
    legend = false,
    dpi = 480,
)
savefig("shearOct.pdf")

plotFEDomain(regularCuboidMesh, dpi = 480)
savefig("fenodeset.pdf")


uGammaTest = [1
     5
     9
    13
    17
    21
    25
    29
    33
    37
    41
    45
    49
    53
    57
    61
    65
    69
    73
    77
    81
    85
    89
    93
    97
   101
   105
   109
   113
   117
   121
   125
   129
   133
   137
   141
   145
   149
   153
   157
   161
   165
   169
   173
   177
   181
   185
   189
    32
    64
    96
   128
   160
   192]

isempty(setdiff(uGammaTest, cantileverBC.uGammaDln))

using Plots
plotlyjs()



uTilde = calc_uTilde(forceDisplacement, regularCuboidMesh, cantileverBC, matParams, network)

ux = uTilde[1:3:end]
uy = uTilde[2:3:end]
uz = uTilde[3:3:end]
uxTest = [0.002902379237464
  -0.000825320267641
  -0.000216505389301
   0.000600933838415
   0.002962973902505
   0.002054587382829
   0.000532260383111
  -0.000540416450265
  -0.000513437090367
  -0.000490402528635
  -0.000006896717830
   0.000474259103626
   0.005397183215013
   0.009359752463691
   0.011706897280988
   0.007291869665819
   0.001849610098638
   0.000042381888479
  -0.001702117485825
  -0.002843414555045
  -0.002410738110906
   0.000539089616592
   0.001820044610570
   0.001234703455380
   0.006602728364944
   0.015312299260891
   0.025589237191566
   0.012678168452899
  -0.000080901942204
  -0.000928424965431
   0.004779683506410
   0.012116230740317
   0.016235807189654
  -0.004825885604473
  -0.005475614327638
  -0.001426544494001
   0.000383569748001
  -0.002632875658310
  -0.015472962888402
  -0.023681134911467
  -0.002317495197571
  -0.000044720241262
  -0.001875362848391
  -0.005706341419847
  -0.010385105064809
  -0.002953997587868
   0.002119801972975
   0.001235749845661
   0.000825318245526
  -0.002902377246890
   0.000540411473077
  -0.000532267836520
  -0.002054593129104
  -0.002962975005666]
uyTest = [0.006736418683634
   0.001427864907960
  -0.000568837442163
  -0.004853668130073
   0.005429278910802
   0.003267939449470
   0.001853527265136
   0.001574961743085
  -0.000909892991336
  -0.001196498288706
  -0.002159487676931
  -0.003753398175197
   0.012104479012635
   0.020456492891772
   0.025313726190326
   0.015967912028855
   0.004258394441267
   0.000136951612747
   0.002102127166782
   0.001900891347573
  -0.002757058791196
  -0.010749182282816
  -0.011816328549622
  -0.008061685925799
   0.011133663353635
   0.024101263652828
   0.039102875634676
   0.021008262032718
   0.001506249617503
  -0.001057542236559
   0.006594654891793
   0.015672178731245
   0.031102310536085
   0.009726722788471
  -0.002387545889418
  -0.001786153091243
   0.003369486742511
   0.007604874919956
   0.014348185530808
   0.001131034211979
  -0.006296249645573
  -0.003745640177948
   0.002919184183457
   0.005486126324909
   0.004463346231076
  -0.009841525638686
  -0.012246278262990
  -0.006950703681434
  -0.001427865129897
  -0.006736420211361
  -0.001574961561103
  -0.001853530003290
  -0.003267946489236
  -0.005429285071188]
uzTest = [0.007979227585252
  -0.001595216413407
   0.001236157444620
  -0.005673318105667
   0.008293770884076
   0.005993892765416
   0.002011169419053
  -0.000819980038978
   0.001104091301591
  -0.000522872096892
  -0.003424476808622
  -0.005491397737416
   0.011382400765936
   0.013338935578350
   0.008721285792825
   0.001445577151523
   0.000718316355596
   0.001430179748307
  -0.002556418925075
  -0.003080950976256
  -0.002751213481589
  -0.004929546030899
  -0.007995530649524
  -0.007572519527238
   0.014145265788676
   0.021928019175189
   0.018591784539152
   0.002491290422773
   0.002776834409128
   0.002296013176721
   0.010859267736399
   0.019333926283518
   0.020285675852899
   0.007236945304440
   0.004653694540224
   0.000570477123928
   0.002222698048389
   0.000587231913281
  -0.001747377517760
   0.000180577548373
  -0.005817246538319
  -0.005196819047295
  -0.002546929151370
  -0.005515336569039
  -0.006145943990649
  -0.006972817101953
  -0.011041482330800
  -0.008537489721315
   0.001595214917576
  -0.007979222870157
   0.000819972886257
  -0.002011182465812
  -0.005993903435130
  -0.008293772390987]

isapprox(ux, uxTest, rtol = 5e-5)
isapprox(uy, uyTest, rtol = 5e-5)
isapprox(uz, uzTest, rtol = 5e-5)

calc_uTilde!(forceDisplacement, regularCuboidMesh, cantileverBC, matParams, network)

isapprox(uTilde[1:3:end], forceDisplacement.uTilde[3 * uGammaDln .- 2])
isapprox(uTilde[2:3:end], forceDisplacement.uTilde[3 * uGammaDln .- 1])
isapprox(uTilde[3:3:end], forceDisplacement.uTilde[3 * uGammaDln])

uTildeTest = [0.002902379237464
   0.006736418683634
   0.007979227585252
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.005397183215013
   0.012104479012635
   0.011382400765936
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.009359752463691
   0.020456492891772
   0.013338935578350
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.011706897280988
   0.025313726190326
   0.008721285792825
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.007291869665819
   0.015967912028855
   0.001445577151523
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.001849610098638
   0.004258394441267
   0.000718316355596
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.000042381888479
   0.000136951612747
   0.001430179748307
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000216505389301
  -0.000568837442163
   0.001236157444620
                   0
                   0
                   0
                   0
                   0
                   0
   0.000825318245526
  -0.001427865129897
   0.001595214917576
   0.002962973902505
   0.005429278910802
   0.008293770884076
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.006602728364944
   0.011133663353635
   0.014145265788676
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.015312299260891
   0.024101263652828
   0.021928019175189
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.025589237191566
   0.039102875634676
   0.018591784539152
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.012678168452899
   0.021008262032718
   0.002491290422773
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000080901942204
   0.001506249617503
   0.002776834409128
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000928424965431
  -0.001057542236559
   0.002296013176721
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000513437090367
  -0.000909892991336
   0.001104091301591
                   0
                   0
                   0
                   0
                   0
                   0
   0.000540411473077
  -0.001574961561103
   0.000819972886257
   0.002054587382829
   0.003267939449470
   0.005993892765416
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.004779683506410
   0.006594654891793
   0.010859267736399
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.012116230740317
   0.015672178731245
   0.019333926283518
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.016235807189654
   0.031102310536085
   0.020285675852899
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.004825885604473
   0.009726722788471
   0.007236945304440
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.005475614327638
  -0.002387545889418
   0.004653694540224
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.001426544494001
  -0.001786153091243
   0.000570477123928
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000490402528635
  -0.001196498288706
  -0.000522872096892
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000532267836520
  -0.001853530003290
  -0.002011182465812
   0.000532260383111
   0.001853527265136
   0.002011169419053
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.000383569748001
   0.003369486742511
   0.002222698048389
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.002632875658310
   0.007604874919956
   0.000587231913281
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.015472962888402
   0.014348185530808
  -0.001747377517760
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.023681134911467
   0.001131034211979
   0.000180577548373
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.002317495197571
  -0.006296249645573
  -0.005817246538319
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000044720241262
  -0.003745640177948
  -0.005196819047295
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000006896717830
  -0.002159487676931
  -0.003424476808622
                   0
                   0
                   0
                   0
                   0
                   0
  -0.002054593129104
  -0.003267946489236
  -0.005993903435130
  -0.000540416450265
   0.001574961743085
  -0.000819980038978
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.001875362848391
   0.002919184183457
  -0.002546929151370
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.005706341419847
   0.005486126324909
  -0.005515336569039
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.010385105064809
   0.004463346231076
  -0.006145943990649
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.002953997587868
  -0.009841525638686
  -0.006972817101953
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.002119801972975
  -0.012246278262990
  -0.011041482330800
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.001235749845661
  -0.006950703681434
  -0.008537489721315
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.000474259103626
  -0.003753398175197
  -0.005491397737416
                   0
                   0
                   0
                   0
                   0
                   0
  -0.002962975005666
  -0.005429285071188
  -0.008293772390987
  -0.000825320267641
   0.001427864907960
  -0.001595216413407
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.001702117485825
   0.002102127166782
  -0.002556418925075
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.002843414555045
   0.001900891347573
  -0.003080950976256
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.002410738110906
  -0.002757058791196
  -0.002751213481589
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.000539089616592
  -0.010749182282816
  -0.004929546030899
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.001820044610570
  -0.011816328549622
  -0.007995530649524
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.001234703455380
  -0.008061685925799
  -0.007572519527238
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.000600933838415
  -0.004853668130073
  -0.005673318105667
                   0
                   0
                   0
                   0
                   0
                   0
  -0.002902377246890
  -0.006736420211361
  -0.007979222870157]

isapprox(uTildeTest, forceDisplacement.uTilde, rtol = 5e-5)






norm(forceDisplacement.uTilde)
norm(uTilde)
forceDisplacement.uTilde


##
using Plots
# plotlyjs()
gr()
##
using BenchmarkTools, LinearAlgebra, StaticArrays, SparseArrays, DDD, LazySets, FastGaussQuadrature
cd(@__DIR__)
fileDislocationParameters = "../inputs/simParams/sampleDislocationParameters.json"
fileMaterialParameters = "../inputs/simParams/sampleMaterialParameters.json"
fileFEMParameters = "../inputs/simParams/sampleFEMParameters.json"
fileIntegrationParameters = "../inputs/simParams/sampleIntegrationParameters.json"
fileSlipSystem = "../data/slipSystems/BCC.json"
fileDislocationLoop = "../inputs/dln/samplePrismShear.json"
fileIntVar = "../inputs/simParams/sampleIntegrationTime.json"
dlnParams, matParams, femParams, intParams, slipSystems, dislocationLoop =
    loadParameters(
        fileDislocationParameters,
        fileMaterialParameters,
        fileFEMParameters,
        fileIntegrationParameters,
        fileSlipSystem,
        fileDislocationLoop,
    )
intVars = loadIntegrationTime(fileIntVar)
# network = DislocationNetwork(dislocationLoop)
femParams = FEMParameters(
                    femParams.type, 
                    femParams.order,
                    femParams.model,
                    57.,
                    43.,
                    37.,
                    23,
                    27,
                    13
                )
##
regularCuboidMesh = buildMesh(matParams, femParams)
cantileverBC, forceDisplacement = Boundaries(femParams, regularCuboidMesh)

calc_uTilde!(forceDisplacement, regularCuboidMesh, cantileverBC, matParams, network)

old = copy(forceDisplacement.uTilde)

array = copy(forceDisplacement.uTilde)
inPlace = copy(forceDisplacement.uTilde)
isequal(old,forceDisplacement.uTilde)

array
forceDisplacement.uTilde

numSeg = network.numSeg[1]
network.nodeVel[:, 1:numSeg] = rand(3, numSeg)
network = remeshSurfaceNetwork!(regularCuboidMesh,  network)



array - inPlace

length(a)
size(a)
size(b)
size(c)
# saveJSON("cantileverBC.json", cantileverBC)
# cantileverDict = loadJSON("cantileverBC.json")
# loadBoundaries(cantileverDict)
# cantileverDict["tGamma"]

cantileverBC.tK[:P]
tDofs = cantileverBC.tDofs


cantileverBC.tK \ A ≈ K[tDofs,tDofs] \ A
@time cantileverBC.tK \ A 
@time K[tDofs,tDofs] \ A


forceDisplacement.u[3 * cantileverBC.mGamma[:node]] .= -1
saveJSON("forceDisplacement.json", forceDisplacement)
forceDisplacement
save("forceDisp.jld2", "forceDisplacement", forceDisplacement, "cantileverBC", cantileverBC)

Matrix(cantileverBC.tK)
serialize("out.out", cantileverBC.tK)

A = deserialize("out.out")
loadedBC.tK

save("out.jld2", "test", cantileverBC.tK)

using JLD2
@save "out.jld2" cantileverBC.tK
regularCuboidMesh.K

loadedForceDisp = load("forceDisp.jld2", "forceDisplacement")

isequal(loadedForceDisp.fHat, forceDisplacement.fHat)

ldict = loadForceDisplacement(loadJSON("forceDisplacement.json"))
sparse(ldict["fTilde"])


dict = loadJSON("forceDisplacement.json")
sparse(Float64.(dict["u"]))


isequal

loadedBC = loadBoundaries(cantileverDict)

loadedBC.tGamma[:node]




##

dx, dy, dz = regularCuboidMesh.dx, regularCuboidMesh.dy, regularCuboidMesh.dz
segLen = dx / 5
prismSquare = DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 8,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 10,  # Number of loops of this type to generate when making a network.
    segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 1],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 1],            # Burgers vector of the segments.
    label = SVector{8,nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(0, 0, 0, dx, dy, dz),  # Distribution range
    dist = Rand(),  # Loop distribution.
)

shearSquare = DislocationLoop(;
    loopType = loopShear(),    # Prismatic loop, all segments are edge segments.
    numSides = 8,   # 5-sided loop.
    nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 10,  # Number of loops of this type to generate when making a network.
    segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 5, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 1],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 1],            # Burgers vector of the segments.
    label = SVector{8,nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(0, 0, 0, dx, dy, dz),  # Distribution range
    dist = Rand(),  # Loop distribution.
)
network = DislocationNetwork((shearSquare, prismSquare))
##
uGamma = cantileverBC.uGamma[:node]
coord = regularCuboidMesh.coord
nodeY = regularCuboidMesh.my + 1
nodeZ = regularCuboidMesh.mz + 1

y = range(0, dy, length = nodeY)
z = range(0, dz, length = nodeZ)

length(y)
@which Boundaries(femParams, regularCuboidMesh)
σ =  reshape(calc_σTilde(coord[:, cantileverBC.uGamma], dlnParams, matParams, network), 6, nodeY, :)
contourf(z, y, σ[1, :, :])
contourf(z, y, σ[2, :, :])
contourf(z, y, σ[3, :, :])
contourf(z, y, σ[4, :, :])
contourf(z, y, σ[5, :, :])
contourf(z, y, σ[6, :, :])

##
numNode = 301
numSeg = numNode - 1
len = numNode + 1
xrange = range(-300, 300, length = len)
yrange = range(-300, 300, length = len)

X = ones(length(yrange)) .* xrange'
Y = ones(length(xrange))' .* yrange
Z = zeros(length(xrange))' .* zeros(len)
points = [X[:]'; Y[:]'; Z[:]']

l = Float64[0; 0; 1]
b = Float64[0; 0; 1]
n = b × l
a = 5 * norm(b)

links = zeros(Int, 2, numSeg)
slipPlane = zeros(3, numSeg)
bVec = zeros(3, numSeg)
coord = [zeros(len)'; zeros(len)'; xrange']
label = zeros(nodeTypeDln, len)
nodeVel = similar(coord)
nodeForce = similar(coord)
for i in 1:numSeg
    links[:, i] .= (i, i + 1)
    slipPlane[:, i] = n
    bVec[:, i] = b
end

matParams = MaterialParameters(;
        crystalStruct = BCC(),
        μ = 1.0,
        μMag = 1.0,
        ν = 0.28,
        E = 0.1,
    )
dlnParams = DislocationParameters(;
        coreRad = a,
        coreRadMag = 1.,
        minSegLen = a + 2,
        maxSegLen = a + 3,
        minArea = a + 1,
        maxArea = a + 2,
        edgeDrag = 1.,
        screwDrag = 1.,
        climbDrag = 1.,
        lineDrag = 1.,
        maxConnect = 4,
        mobility = mobBCC(),
    )
network = DislocationNetwork(;
        links = links,
        slipPlane = slipPlane,
        bVec = bVec,
        coord = coord,
        label = label,
        nodeVel = nodeVel,
        nodeForce = nodeForce,
    )
makeConnect!(network)
getSegmentIdx!(network)

stress = reshape(calc_σTilde(points, dlnParams, matParams, network), 6, len, :)

σ = zeros(6, size(points, 2))
calc_σTilde!(σ, points, dlnParams, matParams, network)
σ =  reshape(σ, 6, len, :)

isequal(σ, stress)

plt = contourf(xrange, yrange, stress[1,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[2,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[3,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[4,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[5,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[6,:,:], levels = 30)
display(plt)


##
numNode = 301
numSeg = numNode - 1
len = numNode + 1
xrange = range(-300, 300, length = len)
yrange = range(-300, 300, length = len)

X = ones(length(yrange)) .* xrange'
Y = ones(length(xrange))' .* yrange
Z = zeros(length(xrange))' .* zeros(len)
points = [X[:]'; Y[:]'; Z[:]']

l = Float64[0; 0; 1]
b = Float64[1; 0; 0]
n = b × l
a = 5 * norm(b)

links = zeros(Int, 2, numSeg)
slipPlane = zeros(3, numSeg)
bVec = zeros(3, numSeg)
coord = [zeros(len)'; zeros(len)'; xrange']
label = zeros(nodeTypeDln, len)
nodeVel = similar(coord)
nodeForce = similar(coord)
for i in 1:numSeg
    links[:, i] .= (i, i + 1)
    slipPlane[:, i] = n
    bVec[:, i] = b
end

matParams = MaterialParameters(;
        crystalStruct = BCC(),
        μ = 1.0,
        μMag = 1.0,
        ν = 0.28,
        E = 0.1,
    )
dlnParams = DislocationParameters(;
        coreRad = a,
        coreRadMag = 1.,
        minSegLen = a + 2,
        maxSegLen = a + 3,
        minArea = a + 1,
        maxArea = a + 2,
        edgeDrag = 1.,
        screwDrag = 1.,
        climbDrag = 1.,
        lineDrag = 1.,
        maxConnect = 4,
        mobility = mobBCC(),
    )
network = DislocationNetwork(;
        links = links,
        slipPlane = slipPlane,
        bVec = bVec,
        coord = coord,
        label = label,
        nodeVel = nodeVel,
        nodeForce = nodeForce,
    )
makeConnect!(network)
getSegmentIdx!(network)

stress = reshape(calc_σTilde(points, dlnParams, matParams, network), 6, len, :)

σ = zeros(6, size(points, 2))
calc_σTilde!(σ, points, dlnParams, matParams, network)
σ =  reshape(σ, 6, len, :)

isequal(σ, stress)

plt = contourf(xrange, yrange, stress[1,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[2,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[3,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[4,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[5,:,:], levels = 30)
display(plt)
plt = contourf(xrange, yrange, stress[6,:,:], levels = 30)
display(plt)

##


matParams = MaterialParameters(; crystalStruct = BCC(), μ = 1., μMag = 1., ν = 0.305, E = 1.)
##
regularCuboidMesh.vertices
dx, dy, dz
fig1 = plotNodes(
    regularCuboidMesh, 
    network,
    m = 1,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,
)
# png(fig1, "original")
xlims!(0, dx)
ylims!(0, dy)
zlims!(0, dz)
plot!(aspect_ratio = (1, 0, 1))

network.nodeVel[:, 1:network.numNode[1]] .= rand(3, network.numNode[1])
network2 = deepcopy(network)
network2 = remeshSurfaceNetwork!(regularCuboidMesh,  network2)

label = network2.label
coord = network2.coord
surface = findall(x -> x == 3, label)
external = findall(x -> x == 5, label)
temporary = findall(x -> x == 6, label)
regularCuboidMesh.faceMidPt[:, 6]

plotNodes!(
    fig1,
    regularCuboidMesh, 
    network2,
    m = 2,
    l = 3,
    linecolor = :orange,
    marker = :circle,
    markercolor = :orange,
    legend = false,
)
regularCuboidMesh.faceMidPt[1, :] .* regularCuboidMesh.faceNorm[1, :] .+ regularCuboidMesh.faceMidPt[2, :] .* regularCuboidMesh.faceNorm[2, :] .+ regularCuboidMesh.faceMidPt[3, :] .* regularCuboidMesh.faceNorm[3, :]
fig2 = plotNodes(
    regularCuboidMesh, 
    network2,
    m = 1,
    l = 3,
    linecolor = :black,
    marker = :circle,
    markercolor = :black,
    legend = false,
)

scatter!(coord[1, surface], coord[2, surface],coord[3, surface], markersize = 2, markercolor = :cyan)

# network2.coord[:, external]
# scatter!(fig1, coord[1, surface], coord[2, surface], coord[3, surface], markersize = 2, markercolor = :red)
# scatter!(fig1, coord[1, external], coord[2, external], coord[3, external], markersize = 2, markercolor = :black)





# scatter!(coord[1, temporary], coord[2, temporary],coord[3, temporary], markersize = 3, markercolor = :green)



maximum(network2.connectivity[1,:])





##



using Plots
vertices = reshape(collect(Iterators.flatten(regularCuboidMesh.vertices.vertices)), 3, 8)
faces = regularCuboidMesh.faces
normals = regularCuboidMesh.faceNorm

# https://stackoverflow.com/questions/47021821/julia-flattening-array-of-array-tuples
using RecursiveArrayTools

regularCuboidMesh.vertices.vertices
norm(regularCuboidMesh.vertices, [1,2,3])



faceCoord = vertices[:, faces]

p = faceCoord[:, 2, :] - faceCoord[:, 1, :]
q = faceCoord[:, 3, :] - faceCoord[:, 1, :]

n = reshape(collect(Iterators.flatten([normalize(p[:,i] × q[:,i]) for i in 1:size(p, 2)])), 3, 6)
isapprox(n, normals)

faceCoord[:, 3, 5]
faceCoord[:, 2, 5]
faceCoord[:, 1, 5]
p[:, 1] × q[:, 1]


regularCuboidMesh.vertices.vertices[1:2]

plotNodes!(regularCuboidMesh, network)


# Infinite domain.
infDom = VPolytope([
    [Inf, Inf, Inf],
    [Inf, Inf, Inf],
    [Inf, Inf, Inf],
    [Inf, Inf, Inf],
    [Inf, Inf, Inf],
    [Inf, Inf, Inf],
    [Inf, Inf, Inf],
    [Inf, Inf, Inf],
    ])

[Inf,Inf,Inf] ∉ infDom

##
# Surface remeshing
using Polyhedra
regularCuboidMesh.vertices.vertices

vertices = [0.0 0.0 0.0;
2000.0 0.0 0.0;
0.0 2000.0 0.0;
2000.0 2000.0 0.0;
0.0 0.0 2000.0;
2000.0 0.0 2000.0;
0.0 2000.0 2000.0;
2000.0 2000.0 2000.0]
vRep = vrep(vertices)
poly = polyhedron(vRep)
npoints(poly)
hrep(poly)
allhalfspaces(poly.hrep)

# @time begin
vertices = regularCuboidMesh.vertices
plotNodes(regularCuboidMesh, network,
m = 1,
l = 3,
linecolor = :blue,
marker = :circle,
markercolor = :blue,
legend = false,)

# Find which nodes are outside the domain.
coord = network.coord
label = network.label
# Findall internal nodes.
idx = findall(x -> x == 1, label)
# Find the location of the nodes that are outside the domain, P.
indices = map((x, y, z) -> eltype(vertices)[x, y, z] ∉ vertices, coord[1, idx], coord[2, idx], coord[3, idx])
outside = findall(indices)
# Change the label of the nodes outside to a temporary flag for newly exited nodes.
# label[outside] .= 6

# Plot nodes outside the domain.
scatter!(coord[1, outside], coord[2, outside], coord[3, outside], markersize = 3)
# end



##
# Define plane
planenorm = Float64[0, 0, 1]
planepnt  = Float64[0, 0, 5]
 
# Define ray
raydir = Float64[0, -1, -2]
raypnt = Float64[0,  0, 10]
 
planenorm = Float64[0, 0, 1]
planepnt  = Float64[0, 0, 5]
raydir = Float64[0, 1, 0]
raypnt = Float64[0,  0, 5]
ψ = linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
isinf(ψ)

planenorm = Float64[0, 0, 1]
planepnt  = Float64[0, 0, 5]
raydir = Float64[0, 1, 0]
raypnt = Float64[0, 0, 6]
ψ = linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
isnothing(ψ)



@btime linePlaneIntersect(planenorm, planepnt, raydir, raypnt)



planenorm = Float64[0, 0, 1]
planepnt  = Float64[0, 0, 5]
raydir = Float64[0, 0, 1]
raypnt = Float64[0,  0, 5]

ψ = linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
@test isapprox(ψ,  [0, -2.5, 5.0])

# Define plane
planenorm = SVector(0, 0, 1)
planepnt  = SVector(0, 0, 5)
 
# Define ray
raydir = SVector(0, -1, 0)
raypnt = SVector(0,  0, 10)

@btime lineplanecollision(planenorm, planepnt, raydir, raypnt)
@btime linePlaneIntersect(planenorm, planepnt, raydir, raypnt)


ψ = lineplanecollision(planenorm, planepnt, raydir, raypnt)
linePlaneIntersect(planenorm, planepnt, raydir, raypnt)
println("Intersection at $ψ")





##
dx, dy, dz = femParams.dx, femParams.dy, femParams.dz
numFEMNode = regularCuboidMesh.numNode
segLen = rand() * (dx * dy * dz) / 1e8
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
forceDisplacement = ForceDisplacement(nothing, uHat, u, nothing, fHat, f)


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

segLen = dx / 10
prismSquare = DislocationLoop(;
    loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
    numSides = 4,   # 5-sided loop.
    nodeSide = 2,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 2,  # Number of loops of this type to generate when making a network.
    segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 1],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 1],            # Burgers vector of the segments.
    label = SVector{8,nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(0 + segLen, 0 + segLen, 0 + segLen, dx - segLen, dy - segLen, dz - segLen),  # Distribution range
    dist = Rand(),  # Loop distribution.
)

shearSquare = DislocationLoop(;
    loopType = loopShear(),    # Prismatic loop, all segments are edge segments.
    numSides = 4,   # 5-sided loop.
    nodeSide = 2,   # One node per side, if 1 nodes will be in the corners.
    numLoops = 2,  # Number of loops of this type to generate when making a network.
    segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
    slipSystem = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
    _slipPlane = slipSystems.slipPlane[:, 1],  # Slip plane of the segments.
    _bVec = slipSystems.bVec[:, 1],            # Burgers vector of the segments.
    label = SVector{8,nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
    buffer = 0,   # Buffer to increase the dislocation spread.
    range = SMatrix{3,2,Float64}(0 + segLen, 0 + segLen, 0 + segLen, dx - segLen, dy - segLen, dz - segLen),  # Distribution range
    dist = Rand(),  # Loop distribution.
)
network = DislocationNetwork((prismSquare, shearSquare))
plotNodes(regularCuboidMesh, network,
    m = 1,
    l = 3,
    linecolor = :blue,
    marker = :circle,
    markercolor = :blue,
    legend = false,)
fPK = calcPKForce(regularCuboidMesh, forceDisplacement, network)

numSeg = network.numSeg[1]
println([network.links[:, 1:numSeg]' network.bVec[:,1:numSeg]' network.slipPlane[:, 1:numSeg]'])
println([network.coord[:, 1:numSeg]' zeros(numSeg)])

println(network.links[:, 1:numSeg])
println(network.bVec[:, 1:numSeg])
println(network.slipPlane[:, 1:numSeg])

println(network.coord[:, 1:numSeg])



fPKTest = 1.0e+02 *
[
    0.983568207059471   0.550732784894686  -0.216417711082393
  -0.074266349477605   0.074266349477605  -0.636911320120312
   0.155045597251722  -0.155045597251722  -0.044881839670221
   0.101431295661336  -0.067753859578523  -0.084592577619930
   0.273213406343103   0.068697355372783  -0.102258025485160
  -0.071275154999800   0.071275154999800  -0.338908217279899
   0.022167783169158  -0.022167783169158  -1.136846650257031
   1.436820246111006   0.473392983022663  -0.481713631544171
  -0.168015348058644  -0.409635686533833  -0.120810169237595
                   0                   0                   0
  -0.127152152429410   0.127152152429410  -0.532135533886545
  -0.607064449348708  -0.376897289465542   0.115083579941583
  -1.327884941063532  -0.715166073499301   0.306359433782116
   0.441147725403124  -0.441147725403124   1.360442463840005
   0.255278751159306  -0.255278751159306   0.748442624130778
                   0                   0                   0
  -0.367714422335496  -0.154013286390333   0.106850567972581
   0.045819274971325  -0.035698104354998   0.081517379326323
  -0.066266665225592  -0.050212580764005  -0.016054084461586
  -0.090887767612372  -0.032124956563013   0.029381405524679
                   0                   0                   0
                   0                   0                   0
                   0                   0                   0
  -0.095974787977392  -0.060606642986924   0.017684072495234
  -0.911537345550952  -0.433723676128685   0.238906834711134
   0.166970234231067  -0.039020251746418   0.205990485977485
  -0.040569130620300  -0.070322668905174   0.029753538284874
  -0.080046044506237  -0.059386350692248   0.010329846906995
   0.000155309763062   0.000275031985160   0.000059861111049
   0.000064985978620   0.000160977029895  -0.000095991051275
                   0                   0                   0
                   0                   0                   0
]

isapprox(fPK', fPKTest)




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
    label = SVector{5,nodeTypeDln}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
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
    label = SVector{7,nodeTypeDln}(1, 1, 1, 1, 1, 2, 1),
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
    label = SVector{5,nodeTypeDln}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
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
    label = SVector{7,nodeTypeDln}(1, 1, 1, 1, 1, 2, 1),
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
label = zeros(nodeTypeDln, numNode)

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
append!(label, nodeTypeDln[3, 3])

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
    label = nodeTypeDln[1; 1; 1; 1; 1; 1; 1; 1; 1; 1],
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
StructTypes.StructType(::Type{nodeTypeDln}) = StructTypes.NumberType()
StructTypes.numbertype(::Type{nodeTypeDln}) = Int
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
label = nodeTypeDln.(newvar["label"])
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
    label = copy(nodeTypeDln.(newvar["label"])),
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
    label = nodeTypeDln[1; 1; 1; 1; 1],    # Node labels, has to be equal to the number of nodes.
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
    label = nodeTypeDln[1; 1; 1; 1; 1; 2; 1],
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
    label = SVector{5,nodeTypeDln}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
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
    label = SVector{7,nodeTypeDln}(1, 1, 1, 1, 1, 2, 1),
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
    label = SVector{5,nodeTypeDln}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
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
    label = nodeTypeDln[1; 1; 1; 1; 1],    # Node labels, has to be equal to the number of nodes.
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
    label = nodeTypeDln[1; 1; 1; 1; 1],    # Node labels, has to be equal to the number of nodes.
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
    label = SVector{5,nodeTypeDln}(1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
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
    label = nodeTypeDln[1; 1; 1; 1; 1; 2; 1],
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
    label = SVector{7,nodeTypeDln}(1, 1, 1, 1, 1, 2, 1),
    buffer = 0.0,
    range = SMatrix{3,2,Float64}(-5000, -5000, -5000, 5000, 5000, 5000),  # Distribution range
    dist = Rand(),
)

##
@btime DislocationNetwork((prismHeptagon, prismPentagon))
@btime DislocationNetwork!(network, (prismHeptagon, prismPentagon))
@btime calcSegForce(dlnParams, matParams, network)
@btime calcSegForce!(dlnParams, matParams, network)
