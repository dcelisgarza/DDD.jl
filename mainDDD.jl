using DDD, Plots
cd(@__DIR__)
plotlyjs()

fileDislocationParameters = "./inputs/simParams/sampleDislocationParameters.json"
fileMaterialParameters = "./inputs/simParams/sampleMaterialParameters.json"
fileFEMParameters = "./inputs/simParams/sampleFEMParameters.json"
fileIntegrationParameters = "./inputs/simParams/sampleIntegrationParameters.json"
fileSlipSystem = "./data/slipSystems/BCC.json"
fileDislocationLoop = "./inputs/dln/FRSources.json"
fileIntVar = "./inputs/simParams/sampleIntegrationTime.json"

# function main(fileDislocationParameters, fileMaterialParameters, fileFEMParameters, fileIntegrationParameters, fileSlipSystem, fileDislocationLoop, fileIntVar)
dlnParams, matParams, femParams, intParams, slipSystems, dislocationLoop = loadParametersJSON(fileDislocationParameters, fileMaterialParameters, fileFEMParameters, fileIntegrationParameters, fileSlipSystem, fileDislocationLoop)
intVars = loadIntegrationTimeJSON(fileIntVar)

mesh = buildMesh(matParams, femParams)
boundaries, forceDisplacement = Boundaries(femParams, mesh)
network = DislocationNetwork(dislocationLoop)
network.bVec[:,:] *= sqrt(3) / 2
network.label[1:8]
network.coord[:,1:8]

calcSegForce!(dlnParams, matParams, mesh, forceDisplacement, network)
dlnMobility!(dlnParams, matParams, network)
network = remeshSurfaceNetwork!(mesh, network)

plotFEDomain(mesh)
plotNodes(mesh, network, m = 1, l = 3, linecolor = :blue, marker = :circle, markercolor = :blue, legend = false)

counter = 0
@time begin
    while counter < 40# intVars.time < intParams.tmax
        intVars = integrate!(intParams, intVars, dlnParams, matParams, mesh, forceDisplacement, network)
        # plotNodes(mesh, network, m = 1, l = 3, linecolor = :blue, marker = :circle, markercolor = :blue, legend = false)
        counter += 1
    end
end
intVars.time
intVars.dt
intParams.reltol
plotNodes(mesh, network, m = 1, l = 3, linecolor = :blue, marker = :circle, markercolor = :blue, legend = false)
network.numNode[1]
# end

# main(fileDislocationParameters, fileMaterialParameters, fileFEMParameters, fileIntegrationParameters, fileSlipSystem, fileDislocationLoop, fileIntVar)