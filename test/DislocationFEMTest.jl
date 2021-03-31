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

    @test isempty(setdiff(uGammaTest, cantileverBC.uGammaDln))

    uTilde = calc_uTilde(regularCuboidMesh, cantileverBC, matParams, network)

    ux = uTilde[1:3:end]
    uy = uTilde[2:3:end]
    uz = uTilde[3:3:end]
    uxTest = [
        0.002902351345488887,
        0.005397126562750102,
        0.009359641926982321,
        0.01170673090883372,
        0.007291725493998129,
        0.0018495485834692132,
        4.236557903560154e-5,
        -0.00021650827733668473,
        0.0008253133519967209,
        0.0029629460677653143,
        0.006602657554255273,
        0.015312108727395595,
        0.025588846219305697,
        0.012677799867875676,
        -8.098032109580808e-5,
        -0.0009284282909694205,
        -0.0005134326512485413,
        0.0005404119399346326,
        0.0020545681539812758,
        0.00477962525180953,
        0.012116022126695285,
        0.016235303324262124,
        -0.004826705167431574,
        -0.005475618765255888,
        -0.0014265217871674203,
        -0.0004903899090375116,
        -0.000532255234515522,
        0.0005322552345154338,
        0.00038355589475357044,
        -0.002632912972606305,
        -0.015472906634952174,
        -0.023681161774827538,
        -0.002317412841145468,
        -4.468666199580377e-5,
        -6.880892240796615e-6,
        -0.002054568153981224,
        -0.0005404119399346637,
        -0.00187535023767425,
        -0.0057063039169111505,
        -0.010385046987745625,
        -0.0029539571317446283,
        0.0021198441533540995,
        0.0012357735270278093,
        0.0004742726937890663,
        -0.0029629460677652792,
        -0.0008253133519966701,
        -0.0017021045567632155,
        -0.002843394458582297,
        -0.0024107193135325805,
        0.0005391062331057765,
        0.0018200619637690372,
        0.0012347172382416167,
        0.0006009438041774093,
        -0.0029023513454888464,
    ]
    uyTest = [
        0.0067363815836566756,
        0.012104398173005865,
        0.020456327309446968,
        0.02531347315436703,
        0.015967697461884905,
        0.004258310757917565,
        0.000136935016830809,
        -0.0005688368686884211,
        -0.0014278511601506416,
        0.0054292511205919925,
        0.011133590727917276,
        0.024101063034655975,
        0.03910246222571273,
        0.0210078769595474,
        0.0015061696382455902,
        -0.001057544898685026,
        -0.0009098880504148402,
        -0.0015749470920582325,
        0.0032679235999618647,
        0.006594615229743162,
        0.01567205085188615,
        0.03110200505064155,
        0.009726316969711102,
        -0.0023875791696053553,
        -0.0017861516064492313,
        -0.0011964949197497293,
        -0.0018535152751203109,
        0.001853515275120392,
        0.0033694655216900955,
        0.0076048333077830595,
        0.01434813388039752,
        0.001130961200306315,
        -0.006296297689713076,
        -0.0037456489215073716,
        -0.002159488471687042,
        -0.0032679235999618314,
        0.0015749470920581972,
        0.0029191581465615618,
        0.005486079417588696,
        0.004463311985544097,
        -0.009841537292647575,
        -0.01224629630206807,
        -0.0069507076726078865,
        -0.003753398090388903,
        -0.005429251120592,
        0.0014278511601506676,
        0.0021021084169944824,
        0.0019008757293132049,
        -0.002757046866739879,
        -0.010749158149593692,
        -0.011816314972861354,
        -0.008061677900843734,
        -0.004853663242082534,
        -0.006736381583656756,
    ]
    uzTest = [
        0.007979193931819731,
        0.011382337565374375,
        0.013338843818961211,
        0.008721215722078795,
        0.0014455787776552964,
        0.0007183293235734062,
        0.0014301694225654616,
        0.0012361385228138791,
        0.0015951911425738617,
        0.008293744049109469,
        0.014145193432864162,
        0.021927866951300702,
        0.018591615331351674,
        0.002491294207535403,
        0.0027768377925153826,
        0.0022959756942656615,
        0.0011040570890579328,
        0.0008199479417571756,
        0.005993888724568718,
        0.010859230173295977,
        0.019333802591945715,
        0.020285444264391116,
        0.007236926045437346,
        0.00465361358956461,
        0.0005704038917998867,
        -0.0005229167518047081,
        -0.0020111921045511297,
        0.0020111921045511444,
        0.0022227278814157146,
        0.0005872825065838683,
        -0.0017472723394227748,
        0.0001805379885031523,
        -0.005817365214834877,
        -0.00519688720088853,
        -0.00342451549514798,
        -0.005993888724568653,
        -0.0008199479417571705,
        -0.0025468808611998463,
        -0.005515269881330632,
        -0.006145908832257655,
        -0.006972815701468629,
        -0.011041518738134993,
        -0.008537520866297306,
        -0.005491419587598859,
        -0.008293744049109446,
        -0.0015951911425738472,
        -0.002556389674620718,
        -0.0030809269872969602,
        -0.002751202465258223,
        -0.0049295376454319825,
        -0.007995531449452647,
        -0.007572525652961138,
        -0.005673325438189695,
        -0.007979193931819718,
    ]

    @test isapprox(ux, uxTest, rtol = 5e-5)
    @test isapprox(uy, uyTest, rtol = 5e-5)
    @test isapprox(uz, uzTest, rtol = 5e-5)

    calc_uTilde!(forceDisplacement, regularCuboidMesh, cantileverBC, matParams, network)

    uGammaDln = cantileverBC.uGammaDln
    @test isapprox(uTilde[1:3:end], forceDisplacement.uTilde[3 * uGammaDln .- 2])
    @test isapprox(uTilde[2:3:end], forceDisplacement.uTilde[3 * uGammaDln .- 1])
    @test isapprox(uTilde[3:3:end], forceDisplacement.uTilde[3 * uGammaDln])

    uTildeTest = [
        0.002902379237464
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
        -0.007979222870157
    ]

    @test isapprox(uTildeTest, forceDisplacement.uTilde, rtol = 5e-5)
end
