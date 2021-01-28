using DDD
using Test, LinearAlgebra
cd(@__DIR__)

@testset "Calculating sigma_tilde" begin
    numNode = 50
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
    )
    dlnParams = DislocationParameters(;
        coreRad = a,
        coreRadMag = 1.,
        minSegLen = a + 2,
        maxSegLen = a + 3,
        minArea = a + 1,
        maxArea = a + 2,
        dragCoeffs = (edge = 1., screw = 1., climb = 1., line = 1.),
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

    @test isequal(σ, stress)

    ##
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
    )
    dlnParams = DislocationParameters(;
        coreRad = a,
        coreRadMag = 1.,
        minSegLen = a + 2,
        maxSegLen = a + 3,
        minArea = a + 1,
        maxArea = a + 2,
        dragCoeffs = (edge = 1., screw = 1., climb = 1., line = 1.),
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

    @test isequal(σ, stress)
end