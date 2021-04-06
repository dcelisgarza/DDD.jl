using DDD
using Test, LinearAlgebra, StaticArrays
cd(@__DIR__)

@testset "Calculating sigma_tilde" begin
    numNode = 2
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

    matParams = MaterialParameters(; crystalStruct = BCC(), μ = 1.0, μMag = 1.0, ν = 0.28)
    dlnParams = DislocationParameters(;
        coreRad = a,
        coreRadMag = 1.0,
        minSegLen = a + 2,
        maxSegLen = a + 3,
        minArea = a + 1,
        maxArea = a + 2,
        dragCoeffs = (edge = 1.0, screw = 1.0, climb = 1.0, line = 1.0),
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

    σTilde = calc_σTilde(points, dlnParams, matParams, network)'

    test_σTilde = [
        0.000248130106066 0.000035465440391 -0.000011334371408 0.000035426063059 -0.000045958542661 -0.000059335461657
        0 0 0 -0.000260417766102 0.000030180624794 0
        -0.000248130106066 -0.000035465440391 0.000011334371408 0.000035426063059 -0.000045958542661 0.000059335461657
        0.000260598566573 0.000390644731208 0.000015657608042 0 -0.000207853177578 0
        0 0 0 0 -0.001871493087898 0
        -0.000260598566573 -0.000390644731208 -0.000015657608042 0 -0.000207853177578 0
        0.000248130106066 0.000035465440391 -0.000011334371408 -0.000035426063059 -0.000045958542661 0.000059335461657
        0 0 0 0.000260417766102 0.000030180624794 0
        -0.000248130106066 -0.000035465440391 0.000011334371408 -0.000035426063059 -0.000045958542661 -0.000059335461657
    ]

    @test isapprox(σTilde[:, 1:4], test_σTilde[:, 1:4])
    @test isapprox(σTilde[:, 6], test_σTilde[:, 5])
    @test isapprox(σTilde[:, 5], test_σTilde[:, 6])

    σ = zeros(6, size(points, 2))
    calc_σTilde!(σ, points, dlnParams, matParams, network)
    @test isequal(σ', σTilde)

    ##

    numNode = 2
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

    matParams = MaterialParameters(; crystalStruct = BCC(), μ = 1.0, μMag = 1.0, ν = 0.28)
    dlnParams = DislocationParameters(;
        coreRad = a,
        coreRadMag = 1.0,
        minSegLen = a + 2,
        maxSegLen = a + 3,
        minArea = a + 1,
        maxArea = a + 2,
        dragCoeffs = (edge = 1.0, screw = 1.0, climb = 1.0, line = 1.0),
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

    σTilde = calc_σTilde(points, dlnParams, matParams, network)'

    test_σTilde =
        1.0e-03 * [
            0 0 0 0 -0.076573455482457 0.076573455482457
            0 0 0 0 -0.187565879762850 0
            0 0 0 0 -0.076573455482457 -0.076573455482457
            0 0 0 0 0 0.187565879762850
            0 0 0 0 0 0
            0 0 0 0 0 -0.187565879762850
            0 0 0 0 0.076573455482457 0.076573455482457
            0 0 0 0 0.187565879762850 0
            0 0 0 0 0.076573455482457 -0.076573455482457
        ]

    @test isapprox(σTilde[:, 1:4], test_σTilde[:, 1:4])
    @test isapprox(σTilde[:, 6], test_σTilde[:, 5])
    @test isapprox(σTilde[:, 5], test_σTilde[:, 6])

    σ = zeros(6, size(points, 2))
    calc_σTilde!(σ, points, dlnParams, matParams, network)
    @test isequal(σ', σTilde)

    ##

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

    matParams = MaterialParameters(; crystalStruct = BCC(), μ = 1.0, μMag = 1.0, ν = 0.28)
    dlnParams = DislocationParameters(;
        coreRad = a,
        coreRadMag = 1.0,
        minSegLen = a + 2,
        maxSegLen = a + 3,
        minArea = a + 1,
        maxArea = a + 2,
        dragCoeffs = (edge = 1.0, screw = 1.0, climb = 1.0, line = 1.0),
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
    σ = reshape(σ, 6, len, :)

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

    matParams = MaterialParameters(; crystalStruct = BCC(), μ = 1.0, μMag = 1.0, ν = 0.28)
    dlnParams = DislocationParameters(;
        coreRad = a,
        coreRadMag = 1.0,
        minSegLen = a + 2,
        maxSegLen = a + 3,
        minArea = a + 1,
        maxArea = a + 2,
        dragCoeffs = (edge = 1.0, screw = 1.0, climb = 1.0, line = 1.0),
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
    σ = reshape(σ, 6, len, :)

    @test isequal(σ, stress)
end

@testset "Calculate uTilde" begin
    dlnParams = DislocationParameters(; mobility = mobBCC())
    matParams = MaterialParameters(; crystalStruct = BCC(), ν = 0.28)
    femParams = FEMParameters(;
        type = DispatchRegularCuboidMesh(),
        order = LinearElement(),
        model = CantileverLoad(),
        dx = 1013.0,
        dy = 1987.0,
        dz = 2999.0,
        mx = 3,
        my = 5,
        mz = 7,
    )
    intParams = IntegrationParameters(; method = AdaptiveEulerTrapezoid())
    slipSystem = SlipSystem(;
        crystalStruct = BCC(),
        slipPlane = Float64[-1; 1; 0],
        bVec = Float64[1; 1; 1],
    )
    dx, dy, dz = femParams.dx, femParams.dy, femParams.dz
    segLen = (dx + dy + dz) / 30
    prismLoop = DislocationLoop(;
        loopType = loopPrism(),    # Prismatic loop, all segments are edge segments.
        numSides = 8,   # 8-sided loop.
        nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 1,   # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystemIdx = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        slipSystem = slipSystem,  # Slip system of the segments.
        label = SVector{8, nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3, 2, Float64}(dx / 2, dy / 2, dz / 2, dx / 2, dy / 2, dz / 2),  # Distribution range
        dist = Zeros(),  # Loop distribution.
    )
    shearLoop = DislocationLoop(;
        loopType = loopShear(),    # Prismatic loop, all segments are edge segments.
        numSides = 8,   # 8-sided loop.
        nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
        numLoops = 1,   # Number of loops of this type to generate when making a network.
        segLen = segLen * SVector{8}(ones(8)),  # Length of each segment between nodes, equal to the number of nodes.
        slipSystemIdx = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
        slipSystem = slipSystem,  # Slip plane of the segments.
        label = SVector{8, nodeTypeDln}(1, 1, 1, 1, 1, 1, 1, 1),    # Node labels, has to be equal to the number of nodes.
        buffer = 0,   # Buffer to increase the dislocation spread.
        range = SMatrix{3, 2, Float64}(dx / 2, dy / 2, dz / 2, dx / 2, dy / 2, dz / 2),  # Distribution range
        dist = Zeros(),  # Loop distribution.
    )
    network = DislocationNetwork([prismLoop, shearLoop])
    regularCuboidMesh = buildMesh(matParams, femParams)
    cantileverBC, forceDisplacement = Boundaries(femParams, regularCuboidMesh)

    uGammaTest = [
        1
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
        192
    ]

    @test [32, 64, 96, 128, 160] == setdiff(uGammaTest, cantileverBC.uGammaDln)

    uTilde = calc_uTilde(regularCuboidMesh, cantileverBC, matParams, network)

    ux = uTilde[1:3:end]
    uy = uTilde[2:3:end]
    uz = uTilde[3:3:end]
    uxTest = [
        0.002902351345488887,
        -0.0008253133519966701,
        -0.00021650827733668473,
        0.0006009438041774093,
        0.0029629460677653143,
        0.0020545681539812758,
        0.0005322552345154338,
        -0.0005404119399346637,
        -0.0005134326512485413,
        -0.0004903899090375116,
        -6.880892240796615e-6,
        0.0004742726937890663,
        0.005397126562750102,
        0.009359641926982321,
        0.01170673090883372,
        0.007291725493998129,
        0.0018495485834692132,
        4.236557903560154e-5,
        -0.0017021045567632155,
        -0.002843394458582297,
        -0.0024107193135325805,
        0.0005391062331057765,
        0.0018200619637690372,
        0.0012347172382416167,
        0.006602657554255273,
        0.00477962525180953,
        0.00038355589475357044,
        -0.00187535023767425,
        0.015312108727395595,
        0.012116022126695285,
        -0.002632912972606305,
        -0.0057063039169111505,
        0.025588846219305697,
        0.016235303324262124,
        -0.015472906634952174,
        -0.010385046987745625,
        0.012677799867875676,
        -0.004826705167431574,
        -0.023681161774827538,
        -0.0029539571317446283,
        -8.098032109580808e-5,
        -0.005475618765255888,
        -0.002317412841145468,
        0.0021198441533540995,
        -0.0009284282909694205,
        -0.0014265217871674203,
        -4.468666199580377e-5,
        0.0012357735270278093,
        0.0008253133519967209,
        -0.0029023513454888464,
        0.0005404119399346326,
        -0.000532255234515522,
        -0.002054568153981224,
        -0.0029629460677652792,
    ]
    uyTest = [
        0.0067363815836566756,
        0.0014278511601506676,
        -0.0005688368686884211,
        -0.004853663242082534,
        0.0054292511205919925,
        0.0032679235999618647,
        0.001853515275120392,
        0.0015749470920581972,
        -0.0009098880504148402,
        -0.0011964949197497293,
        -0.002159488471687042,
        -0.003753398090388903,
        0.012104398173005865,
        0.020456327309446968,
        0.02531347315436703,
        0.015967697461884905,
        0.004258310757917565,
        0.000136935016830809,
        0.0021021084169944824,
        0.0019008757293132049,
        -0.002757046866739879,
        -0.010749158149593692,
        -0.011816314972861354,
        -0.008061677900843734,
        0.011133590727917276,
        0.006594615229743162,
        0.0033694655216900955,
        0.0029191581465615618,
        0.024101063034655975,
        0.01567205085188615,
        0.0076048333077830595,
        0.005486079417588696,
        0.03910246222571273,
        0.03110200505064155,
        0.01434813388039752,
        0.004463311985544097,
        0.0210078769595474,
        0.009726316969711102,
        0.001130961200306315,
        -0.009841537292647575,
        0.0015061696382455902,
        -0.0023875791696053553,
        -0.006296297689713076,
        -0.01224629630206807,
        -0.001057544898685026,
        -0.0017861516064492313,
        -0.0037456489215073716,
        -0.0069507076726078865,
        -0.0014278511601506416,
        -0.006736381583656756,
        -0.0015749470920582325,
        -0.0018535152751203109,
        -0.0032679235999618314,
        -0.005429251120592,
    ]
    uzTest = [
        0.007979193931819731,
        -0.0015951911425738472,
        0.0012361385228138791,
        -0.005673325438189695,
        0.008293744049109469,
        0.005993888724568718,
        0.0020111921045511444,
        -0.0008199479417571705,
        0.0011040570890579328,
        -0.0005229167518047081,
        -0.00342451549514798,
        -0.005491419587598859,
        0.011382337565374375,
        0.013338843818961211,
        0.008721215722078795,
        0.0014455787776552964,
        0.0007183293235734062,
        0.0014301694225654616,
        -0.002556389674620718,
        -0.0030809269872969602,
        -0.002751202465258223,
        -0.0049295376454319825,
        -0.007995531449452647,
        -0.007572525652961138,
        0.014145193432864162,
        0.010859230173295977,
        0.0022227278814157146,
        -0.0025468808611998463,
        0.021927866951300702,
        0.019333802591945715,
        0.0005872825065838683,
        -0.005515269881330632,
        0.018591615331351674,
        0.020285444264391116,
        -0.0017472723394227748,
        -0.006145908832257655,
        0.002491294207535403,
        0.007236926045437346,
        0.0001805379885031523,
        -0.006972815701468629,
        0.0027768377925153826,
        0.00465361358956461,
        -0.005817365214834877,
        -0.011041518738134993,
        0.0022959756942656615,
        0.0005704038917998867,
        -0.00519688720088853,
        -0.008537520866297306,
        0.0015951911425738617,
        -0.007979193931819718,
        0.0008199479417571756,
        -0.0020111921045511297,
        -0.005993888724568653,
        -0.008293744049109446,
    ]
    @test isapprox(ux, uxTest, rtol = 5e-5)
    @test isapprox(uy, uyTest, rtol = 5e-5)
    @test isapprox(uz, uzTest, rtol = 5e-5)

    calc_uTilde!(forceDisplacement, regularCuboidMesh, cantileverBC, matParams, network)

    uGammaDofs = cantileverBC.uDofsDln
    @test isapprox(uTilde, forceDisplacement.uTilde[uGammaDofs])

    uTildeTest = [
        0.002902351345488887,
        -0.0017021045567632155,
        0.012677799867875676,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0014278511601506676,
        0.0019008757293132049,
        0.009726316969711102,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.0015951911425738472,
        -0.0030809269872969602,
        0.007236926045437346,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.00021650827733668473,
        -0.0024107193135325805,
        -0.023681161774827538,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.0005688368686884211,
        -0.002757046866739879,
        0.001130961200306315,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0067363815836566756,
        0.0021021084169944824,
        0.0210078769595474,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0029629460677653143,
        0.0018200619637690372,
        -8.098032109580808e-5,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.0005134326512485413,
        0.00038355589475357044,
        -0.0009284282909694205,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.0009098880504148402,
        0.0033694655216900955,
        -0.001057544898685026,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0011040570890579328,
        0.0022227278814157146,
        0.0022959756942656615,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.0004903899090375116,
        -0.00187535023767425,
        -0.0014265217871674203,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0005322552345154338,
        0.006602657554255273,
        -0.002317412841145468,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0054292511205919925,
        -0.011816314972861354,
        0.0015061696382455902,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.0011964949197497293,
        0.0029191581465615618,
        -0.0017861516064492313,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.0005229167518047081,
        -0.0025468808611998463,
        0.0005704038917998867,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -6.880892240796615e-6,
        0.015312108727395595,
        -4.468666199580377e-5,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.002159488471687042,
        0.024101063034655975,
        -0.0037456489215073716,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.001853515275120392,
        0.011133590727917276,
        -0.006296297689713076,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.008293744049109469,
        -0.007995531449452647,
        0.0027768377925153826,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.00342451549514798,
        0.021927866951300702,
        -0.00519688720088853,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0004742726937890663,
        0.012116022126695285,
        0.0012357735270278093,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.003753398090388903,
        0.01567205085188615,
        -0.0069507076726078865,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.005491419587598859,
        0.019333802591945715,
        -0.008537520866297306,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0020111921045511444,
        0.014145193432864162,
        -0.005817365214834877,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0020545681539812758,
        0.0012347172382416167,
        -0.005475618765255888,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.005397126562750102,
        -0.002632912972606305,
        0.0008253133519967209,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.012104398173005865,
        0.0076048333077830595,
        -0.0014278511601506416,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.011382337565374375,
        0.0005872825065838683,
        0.0015951911425738617,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.009359641926982321,
        -0.0057063039169111505,
        -0.0029023513454888464,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.0005404119399346637,
        0.00477962525180953,
        0.0021198441533540995,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0032679235999618647,
        -0.008061677900843734,
        -0.0023875791696053553,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.020456327309446968,
        0.005486079417588696,
        -0.006736381583656756,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.013338843818961211,
        -0.005515269881330632,
        -0.007979193931819718,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.01170673090883372,
        0.025588846219305697,
        0.0005404119399346326,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.02531347315436703,
        0.03910246222571273,
        -0.0015749470920582325,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0015749470920581972,
        0.006594615229743162,
        -0.01224629630206807,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.005993888724568718,
        -0.007572525652961138,
        0.00465361358956461,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.008721215722078795,
        0.018591615331351674,
        0.0008199479417571756,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.007291725493998129,
        0.016235303324262124,
        -0.000532255234515522,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.015967697461884905,
        0.03110200505064155,
        -0.0018535152751203109,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0014455787776552964,
        0.020285444264391116,
        -0.0020111921045511297,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.0008199479417571705,
        0.010859230173295977,
        -0.011041518738134993,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.007979193931819731,
        -0.002556389674620718,
        0.002491294207535403,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0018495485834692132,
        -0.015472906634952174,
        -0.002054568153981224,
        0.0012361385228138791,
        -0.002751202465258223,
        0.0001805379885031523,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0007183293235734062,
        -0.0017472723394227748,
        -0.005993888724568653,
        0.0006009438041774093,
        0.0005391062331057765,
        -0.0029539571317446283,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        4.236557903560154e-5,
        -0.010385046987745625,
        -0.0029629460677652792,
        -0.004853663242082534,
        -0.010749158149593692,
        -0.009841537292647575,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.000136935016830809,
        0.004463311985544097,
        -0.005429251120592,
        -0.005673325438189695,
        -0.0049295376454319825,
        -0.006972815701468629,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0014301694225654616,
        -0.006145908832257655,
        -0.008293744049109446,
        -0.0008253133519966701,
        -0.002843394458582297,
        -0.004826705167431574,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.004258310757917565,
        0.01434813388039752,
        -0.0032679235999618314,
    ]
    @test isapprox(uTildeTest, forceDisplacement.uTilde, rtol = 5e-5)
end
