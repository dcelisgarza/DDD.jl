# Using DDD.jl

## Initialisation

Before running a simulation we need to initialise the simulation. For this example, we will use the keyword initialisers because they automatically calculate derived quantities, perform input validations, provide default values, and are make for self-documenting code.

Dislocations live in a material, as such we need a few constants that describe it. These are encapsulated in the immutable <sup>[1](#1)</sup> structure `MaterialP`. Note that we use unicode to denote variables as per convention written `\mu -> μ` and `\nu -> ν`. Here we create a basic material.
```julia
julia> materialP = MaterialP(;
          μ = 1.0,                  # Shear modulus.
          μMag = 145e3,             # Shear modulus magnitude.
          ν = 0.28,                 # Poisson ratio.
          E = 1.0,                  # Young's modulus, MPa.
          crystalStruct = BCC(),    # Crystal structure.
          σPN = 0.0                 # Peierls-Nabarro stress for the material.
        )
MaterialP{Float64,BCC}(1.0, 145000.0, 0.28, 1.0, 1.3888888888888888, 0.3888888888888889, 0.07957747154594767, 0.039788735772973836, 0.11052426603603843, BCC(), 0.0)
```
Note that a few extra constants have been automatically calculated by the constructor. We find what these correspond to using the `fieldnames()` on the type of `materialP`, which is `MaterialP`.
```julia
julia> fieldnames(typeof(materialP))
(:μ, :μMag, :ν, :E, :σPN, :omνInv, :νomνInv, :μ4π, :μ8π, :μ4πν, :crystalStruct)
```
Where `omνInv = 1/(1-ν)`, `νomνInv = v/(1-ν)`, `μ4π = μ/(4π)`, `μ8π = μ/(8π)`, `μ4πν = μ/[4π(1-ν)]`. These precomputed variables are used in various places and are there to avoid recalculating them later.

Our dislocations also have certain constant characteristics that are encapsulated in their own immutable structure, `DislocationP`. These parameters are somewhat arbitrary as long as they approximately hold certain proportions.
```julia
julia> dislocationP = DislocationP(;
          coreRad = 90.0,       # Dislocation core radius, referred to as a.
          coreRadMag = 3.2e-4,  # Magnitude of the core radius.
          minSegLen = 320.0,    # Minimum segment length.
          maxSegLen = 1600.0,   # Maximum segment length.
          minArea = 45000.0,    # Minimum allowable area enclosed by two segments.
          maxArea = 20*45000.0, # Maximum allowable area enclosed by two segments.
          maxConnect = 4,       # Maximum number of connections a node can have.
          remesh = true,        # Flag for remeshing.
          collision = true,     # Flag for collision checking.
          separation = true,    # Flag for node separation.
          virtualRemesh = true, # Flag for remeshing virtual nodes.
          edgeDrag = 1.0,       # Drag coefficient for edge segments.
          screwDrag = 2.0,      # Drag coefficient for screw segments.
          climbDrag = 1e10,     # Drag coefficient along the climb direction.
          lineDrag = 0.0,       # Drag coefficient along the line direction.
          mobility = mobBCC(),  # Mobility type for mobility function specialisation.
        )
DislocationP{Float64,Int64,Bool,mobBCC}(90.0, 8100.0, 0.00032, 320.0, 1600.0, 45000.0, 900000.0, 4, true, true, true, true, 1.0, 2.0, 1.0e10, 0.0, mobBCC())
```

The integration parameters are placed into the following mutable structure.
```julia
julia> integrationP = IntegrationP(;
          dt = 1e3,
          tmin = 0.0,
          tmax = 1e10,
          method = CustomTrapezoid(),
          abstol = 1e-6,
          reltol = 1e-6,
          time = 0.0,
          step = 0,
        )
IntegrationP{Float64,CustomTrapezoid,Int64}(1000.0, 0.0, 1.0e10, CustomTrapezoid(), 1.0e-6, 1.0e-6, 0.0, 0)
```

>[!WARNING]
>
>`IntegrationP` will undergo revisions. Probably be split into two, or perhaps eliminated completely in order to use/extend the state of the art `DifferentialEquations.jl` framework.

Within a given material, we have multiple slip systems, which can be loaded into their own immutable structure. Here we only define a single slip system, but we have the capability of adding more by making the `slipPlane` and `bVec` arguments `n × 3` matrices rather than vectors.
```julia
julia> slipSystems = SlipSystem(;
          crystalStruct = BCC(),
          slipPlane = [1.0; 1.0; 1.0],  # Slip plane.
          bVec = [1.0; -1.0; 0.0]       # Burgers vector.
       )
SlipSystem{BCC,Array{Float64,1}}(BCC(), [1.0, 1.0, 1.0], [1.0, -1.0, 0.0])
```
>[!WARNING]
>
>This may change to perform validity checks regarding the relationship between burgers vector and slip plane.

We also need dislocation sources. We make use of Julia's type system to create standard functions for loop generation. We provide a way of easily and quickly generating loops whose segments inhabit the same slip system. However, new `DislocationLoop()` methods can be made by subtyping `AbstractDlnStr`, and dispatching on the new type. One may of also course also use the default constructor and build the initial structures manually.

Here we make a regular pentagonal prismatic dislocation loop, and a regular hexagonal prismatic dislocation loop. Note that the segments may be of arbitrary length, but having asymmetric sides may result in a very ugly and irregular dislocations that may be unphysical or may end up remeshing once the simulation gets under way. As such we recommend making the segment lengths symmetric.
```julia
julia> prisPentagon = DislocationLoop(
          loopPrism();    # Prismatic loop, all segments are edge segments.
          numSides = 5,   # 5-sided loop.
          nodeSide = 1,   # One node per side, if 1 nodes will be in the corners.
          numLoops = 20,  # Number of loops of this type to generate when making a network.
          segLen = 10 * ones(5),  # Length of each segment between nodes, equal to the number of nodes.
          slipSystem = 1, # Slip System (assuming slip systems are stored in a file, this is the index).
          _slipPlane = slipSystems.slipPlane,  # Slip plane of the segments.
          _bVec = slipSystems.bVec,            # Burgers vector of the segments.
          label = nodeType[1; 2; 1; 2; 1],    # Node labels, has to be equal to the number of nodes.
          buffer = 0.0,   # Buffer to increase the dislocation spread.
          range = Float64[          # Distribution range
                        -100 100; # xmin, xmax          
                        -100 100; # ymin, ymax
                        -100 100  # zmin, zmax
                      ],
          dist = Rand(),  # Loop distribution.
      )
DislocationLoop{loopPrism,Int64,Array{Float64,1},Int64,Array{Int64,2},Array{Float64,2},Array{nodeType,1},Float64,Rand}(loopPrism(), 5, 1, 20, [10.0, 10.0, 10.0, 10.0, 10.0], 1, [1 2 … 4 5; 2 3 … 5 1], [0.5773502691896258 0.5773502691896258 … 0.5773502691896258 0.5773502691896258; 0.5773502691896258 0.5773502691896258 … 0.5773502691896258 0.5773502691896258; 0.5773502691896258 0.5773502691896258 … 0.5773502691896258 0.5773502691896258], [0.7071067811865475 0.7071067811865475 … 0.7071067811865475 0.7071067811865475; -0.7071067811865475 -0.7071067811865475 … -0.7071067811865475 -0.7071067811865475; 0.0 0.0 … 0.0 0.0], [-1.932030909139515 4.820453044614565 … -1.785143053581134 -6.014513813778146; -1.932030909139515 4.820453044614565 … -1.785143053581134 -6.014513813778146; -8.055755266097462 -5.087941102678986 … 8.123251093712414 0.10921054317980072], nodeType[DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob], 0.0, [-100.0 100.0; -100.0 100.0; -100.0 100.0], Rand())

julia> shearHexagon = DislocationLoop(
          loopShear();    # Shear loop
          numSides = 6,
          nodeSide = 3,   # 3 nodes per side, it devides the side into equal segments.
          numLoops = 20,
          segLen = 10 * ones(3 * 6) / 3,  # The hexagon's side length is 10, each segment is 10/3.
          slipSystem = 1,
          _slipPlane = slipSystems.slipPlane,
          _bVec = slipSystems.bVec,
          label = nodeType[1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2],
          buffer = 0.0,
          range = Float64[
                        -100 100;
                        -100 100;
                        -100 100
                      ],
          dist = Rand(),
      )
DislocationLoop{loopShear,Int64,Array{Float64,1},Int64,Array{Int64,2},Array{Float64,2},Array{nodeType,1},Float64,Rand}(loopShear(), 6, 3, 20, [3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335], 1, [1 2 … 17 18; 2 3 … 18 1], [0.5773502691896258 0.5773502691896258 … 0.5773502691896258 0.5773502691896258; 0.5773502691896258 0.5773502691896258 … 0.5773502691896258 0.5773502691896258; 0.5773502691896258 0.5773502691896258 … 0.5773502691896258 0.5773502691896258], [0.7071067811865475 0.7071067811865475 … 0.7071067811865475 0.7071067811865475; -0.7071067811865475 -0.7071067811865475 … -0.7071067811865475 -0.7071067811865475; 0.0 0.0 … 0.0 0.0], [5.443310539518175 6.8041381743977185 … 1.3608276348795458 4.082482904638633; -6.804138174397717 -5.443310539518174 … -6.804138174397715 -8.164965809277255; 1.3608276348795436 -1.3608276348795432
… 5.443310539518167 4.082482904638622], nodeType[DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix], 0.0, [-100.0 100.0; -100.0 100.0; -100.0 100.0], Rand())
```
The dislocation loops will be centred about the origin, but the `range`, `buffer` and `dist` parameters will distribute the dislocations about the simulation domain once the dislocation network is generated. The type of `dist` must be a concrete subtype of `AbstractDistribution` and `loopDistribution()` method should dispatch on this concrete subtype. If a non-suported distribution is required, you only need to create a concrete subtype of `AbstractDistribution` and a new method of `loopDistribution()` to dispatch on the new type. This is all the reworking needed, since multiple dispatch will take care of any new distributions when generating the dislocation network.

Note also the array of `nodeType`, this is an enumerated type which ensures node types are limited to only those supported by the model, lowers memory footprint and increases performance.

We can then plot our loops to see our handy work. We use `plotlyjs()` because it provides a nice interactive experience, but of course, since this is Julia any plotting backend will work. Note that since they have the same slip system but one is a shear and the other a prismatic loop, they are orthogonal to each other.
```julia
julia> using Plots
julia> plotlyjs()
julia> fig1 = plotNodes(
          shearHexagon,
          m = 1,
          l = 3,
          linecolor = :blue,
          markercolor = :blue,
          legend = false,
        )
julia> plotNodes!(fig1, prisPentagon, m = 1, l = 3,
                  linecolor = :red, markercolor = :red, legend = false)
julia> plot!(fig1, camera=(100,35), size=(400,400))
```
![loops](/examples/loops.png)

After generating our primitive loops, we can create a network using either a vector of dislocation loops or a single dislocation loop. The network may also be created manually, and new constructor methods may be defined for bespoke cases. For our purposes, we use the constructor that dispatches on `Union{DislocationLoop, AbstractVector{<:DislocationLoop}}`, meaning a single variable whose type is `DislocationLoop` or a vector of them. Here we use a vector with both our loop structures.

Since the dislocation network is a constantly evolving entity, this necessarily means this is a mutable structure.
```julia
julia> network = DislocationNetwork(
          [shearHexagon, prisPentagon]; # Dispatch type, bespoke functions dispatch on this.
          memBuffer = 1 # Buffer for memory allocation.
       )
 DislocationNetwork{Array{Int64,2},Array{Float64,2},Array{nodeType,1},Int64,Array{Int64,2},Array{Float64,3}}([1 2 … 459 460; 2 3 … 460 456], [0.5773502691896258 0.5773502691896258 … 0.5773502691896258 0.5773502691896258; 0.5773502691896258 0.5773502691896258 … 0.5773502691896258 0.5773502691896258; 0.5773502691896258 0.5773502691896258 … 0.5773502691896258 0.5773502691896258], [0.7071067811865475 0.7071067811865475 … 0.7071067811865475 0.7071067811865475; -0.7071067811865475 -0.7071067811865475 … -0.7071067811865475 -0.7071067811865475; 0.0 0.0 … 0.0 0.0], [46.39761283211718 47.75844046699673 … -57.496002894414175 -61.72537365461118; -49.937613036904054 -48.57678540202451 … -59.28380906452193 -63.51317982471895; 41.67317243370178 38.9515171639427 … -2.1285556467706765 -10.14259619730329], nodeType[DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix  …  DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intMob, DDD.intFix, DDD.intMob, DDD.intFix, DDD.intMob], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], 460, 460, 4, [2 2 … 2 2; 1 1 … 458 459; … ; 0 0 … 0 0; 0 0 … 0 0], [1 2 … 2 2; 1 1 … 1 2], [1 1 2; 2 2 3; … ; 459 459 460; 460 460 456],
 [0.0 0.0; 0.0 0.0; 0.0 0.0]
 [0.0 0.0; 0.0 0.0; 0.0 0.0]
 [0.0 0.0; 0.0 0.0; 0.0 0.0]
 ...
 [0.0 0.0; 0.0 0.0; 0.0 0.0]
 [0.0 0.0; 0.0 0.0; 0.0 0.0]
 [0.0 0.0; 0.0 0.0; 0.0 0.0])
```
This method automatically takes the previously defined loops and scatters them according to the parameters provided in the `DislocationLoop` structure. Furthermore, the `memBuffer` defaults to 10. The number of entries allocated for the matrices is the total number of nodes in the network times `memBuffer`. Here we allocate enough memory for all the nodes but no more. Since julia is dynamic we can allocate memory when needed. However for performance reasons it is advisable to minimise memory management as much as possible.

This function will also automatically calculate other quantities to keep track of the network's links, nodes and segments.
```julia
julia> fieldnames(typeof(network))
(:links, :slipPlane, :bVec, :coord, :label, :nodeVel, :numNode, :numSeg, :maxConnect, :connectivity, :linksConnect, :segIdx, :segForce)
```

We can also view our network with `plotNodes`.
```julia
julia> fig2 = plotNodes(
          network,
          m = 1,
          l = 3,
          linecolor = :blue,
          markercolor = :blue,
          legend = false,
        )
julia> plot!(fig2, camera=(110,40), size=(400,400))
```
![network](/examples/network.png)

<a name="1">1</a>: Immutability is translated into code performance.

## IO

The package provides a way to load and save its parameters using `JSON` files. While this is *not* the most performant format for IO, it is a popular and portable, web-friendly file format that is very human readable (and therefore easy to manually create). It also produces smaller file sizes both compressed and uncompressed. Which is why it is so popular for online data sharing.

`JSON` files are representations of dictionaries with `(key, value)` pairs, which are analogous to the `(key, value)` pair of structures. This makes it so any changes to any structure will automatically be taken care of by the `JSON` library. Arrays are recursively linearised into vectors of vectors using the calling language's preferred storage order. This means arrays preserve their shape and dimensionality regardless of whether the inputting or outputting language stores arrays in column- or row-major order.

They also have the added advantage of being designed for sending over the web, so they have small compressed and uncompressed file sizes, smaller than `BSON`, `HDF5` and its variants such as `MAT` and `JLD2` formats. They also aren't plagued by the portability issues these other formats have, as well as being generally easier to read, create and work with, while being within a factor of ≈3 of the aforementioned filetypes. They do however, allocate more memory than the others while being created, but the end result is half the size or less of these other file formats.

### Sample JSON File

This is a sample `JSON` file for a dislocation loop. They can be compactified by editors to decrease storage space by removing unnecessary line breaks and spaces. Here we show a somewhat longified view which is very human readable and trivially easy to create manually. Note that arrays are recursively linearised as vectors of vectors, where the linearisation follows the calling language's memory order. This means arrays will keep their shape and dimensionality regardless of the language that opens the JSON file.
```JSON
[
  {
    "loopType": "DDD.loopPrism()",
    "numSides": 4,
    "nodeSide": 2,
    "numLoops": 1,
    "segLen": [1, 1, 1, 1, 1, 1, 1, 1],
    "slipSystem": 1,
    "label": [2, 1, 2, 1, 2, 1, 2, 1],
    "buffer": 0,
    "range": [[0, 0, 0], [0, 0, 0]],
    "dist": "DDD.Zeros()"
  }
]
```
This file describes an array, denoted by the `[]` at the top and bottom of the file, of a structure denoted by the `{}` on the second and penultimate lines. We could remove the `[]` but having all files be represent arrays (even if they are of length 1) simplifies users' and developers' lives by letting the same IO functions work for every case.

The keys are on the left side of the colon and the values on the right. This would get loaded to a dictionary with the same `(key, value)` pair shown here. Since the keys are the structure's field names and the values their value, everything can be easily matched to the constructor function.

Also, with since `JSON` files represent dictionaries, users can add irrelevant data to the file either as extra entries to the top level array, or into the structure definition without breaking their ability to use the file.

### Initialisation, Data Dump, and Reloading

One can load all their parameters at once like so.
```julia
fileDislocationP = "../inputs/simParams/sampleDislocationP.JSON"
fileMaterialP = "../inputs/simParams/sampleMaterialP.JSON"
fileIntegrationP = "../inputs/simParams/sampleIntegrationP.JSON"
fileSlipSystem = "../data/slipSystems/SlipSystems.JSON"
fileDislocationLoop = "../inputs/dln/samplePrismShear.JSON"
dlnParams, matParams, intParams, slipSystems, dislocationLoop = loadParams(
    fileDislocationP,
    fileMaterialP,
    fileIntegrationP,
    fileSlipSystem,
    fileDislocationLoop,
)
```
which not only loads the data but returns the aforementioned structures. If there is a single file holding all the parameters, then all the filenames would be the same, but nothing else would change as the file would be loaded into a large dictionary and only the relevant `(key, value)` pairs are used in each case.

Users may also load individual structures as follows.
```julia
dictDislocationP = load(fileDislocationP)
dislocationP = loadDislocationP(dictDislocationP[1])

dictMaterialP = load(fileMaterialP)
materialP = loadMaterialP(dictMaterialP[1])

dictIntegrationP = load(fileIntegrationP)
integrationP = loadIntegrationP(dictIntegrationP[1])

dictSlipSystem = load(fileSlipSystem)
slipSystems = loadSlipSystem(dictSlipSystem[1])

# There can be multiple dislocation types per simulation.
dictDislocationLoop = load(fileDislocationLoop)
dislocationLoop = zeros(DislocationLoop, length(dictDislocationLoop))
for i in eachindex(dislocationLoop)
    dislocationLoop[i] = loadDislocationLoop(dictDislocationLoop[i], slipSystems)
end
```
The reason why all the dictionary arguments of the `load<struct_name>` all have an index is that the files specify an array, as that keeps all files consistent with each other, particularly when saving more than one variable or constant in a single file. Individually loading files like this is useful when recovering previous save states where the data was dumped into a single file, as shown here.
```julia
# Dump simulation parameters into a single file. Creates an array where each entry is one of the structs.
paramDump = "../outputs/simParams/sampleDump.JSON"
save(paramDump, dlnParams, matParams, intParams, slipSystems, dislocationLoop)

# Dump network data into a separate file.
networkDump = "../outputs/dln/sampleNetwork.JSON"
save(networkDump, network)

# Reload parameters.
simulation = load(paramDump)
dlnParams2 = loadDislocationP(simulation[1])
matParams2 = loadMaterialP(simulation[2])
intParams2 = loadIntegrationP(simulation[3])
slipSystems2 = loadSlipSystem(simulation[4])
dislocationLoop2 = zeros(DislocationLoop, length(simulation[5]))
for i in eachindex(dislocationLoop2)
    dislocationLoop2[i] = loadDislocationLoop(simulation[5][i], slipSystems2)
end

# Reload network.
network2 = loadNetwork(networkDump)
```
The reason why `loadNetwork()` is different from the others is that the other values are constants, so for record keeping they would only need to be saved once per simulation. The network might have to be saved at multiple times so it gets a simpler function that calls `load()` internally.

### Against the Unbridled Pursuit of Performance

For the sake of open, reproducible and portable science it is recommended users utilise `JSON` or a standard delimited file format for their IO. If IO is a performance bottleneck these are some incremental steps one should take to improve it before creating a custom IO format. Beware that your mileage may vary when using other IO formats, some may not be fully mature yet others may be abandoned in favour of better implementations.

1. Use buffered IO.
1. Use Julia's in-built task and asyncronous functionality via `tasks` and `async` for either multiple IO streams or an asyncronous IO process while the other threads/cores carry on with the simulation.
1. Use [`BSON`](https://github.com/JuliaIO/BSON.jl), `JSON`'s binary counterpart, though this may break compatibility with other systems, particularly those with different word size and architecture. Furthermore, the binary nature of `BSON` may be used to inject code into a programme so should only be used for self-generated files.
1. Use [`JLD2`](https://github.com/JuliaIO/JLD2.jl). Though the package currently says it has remained largely untested in the wild.
1. Use `DelimitedFiles`.
1. Use binary streams.
1. Use [`Parquet`](https://github.com/JuliaIO/Parquet.jl)
1. Create your own format and IO stream.

TO BE WRITTEN: HOW TO EXTEND METHODS TO EXPAND FUNCTIONALITY

# TODO

## Working Objectives
- [x] IO
  - [x] Input validation
    - [ ] Sensible input generators
  - [ ] Performance
    - [ ] Compression
    - [ ] Asyncronicity
- [ ] Topology functions
  - [ ] Internal Remeshing
    - [ ] Coarsen mesh
      - [x] Merge nodes
        - [ ] Test all edge cases
    - [ ] Refine mesh
      - [x] Split nodes
    - [ ] Surface remeshing
    - [ ] Virtual node remeshing
- [x] Self-segment force
- [x] Seg-seg force
  - [ ] Test tiny segment edge case
  - [ ] Distributed and gpu parallelisation
- [ ] PK force
  - [x] Implementation
  - [ ] Tests
- [ ] Post processing
  - [x] Plot nodes
    - [ ] Asyncronicity
  - [ ] Plot recipe
  - [ ] Statistical analysis
- [ ] Mobility function
  - [ ] BCC
  - [ ] FCC
- [ ] Integration
  - [ ] Refactor integrator structures
  - [ ] Look into using [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) for structure and perhaps use/extension of methods
  - [ ] Make integrator
- [ ] Couple to FEM, perhaps use a package from [JuliaFEM](http://www.juliafem.org/).
  - [ ] Mesh and FE matrices generation
  - [ ] Boundary conditions
    - [ ] Neuman
    - [ ] Dirichlet
  - [ ] Displacements
    - [ ] Parallelisation
  - [ ] Tractions
    - [ ] Parallelisation

### Tentative Objectives

- [ ] Custom 3-vec type, place x,y,z coordinates in contiguous memory instead of columns, ie [x1 y1 z1; x2 y2 z2] -> [x1;y1;z1;x2;y2;z2], have to define custom array type, `getindex(arr, (a,b)) = arr[3*(a-1)+b]`, out of bounds and all the rest. Watch [this](https://www.youtube.com/watch?v=jS9eouMJf_Y).
- [ ] Keep an eye on [JuliaIO](https://github.com/JuliaIO), [JuliaFEM](https://github.com/JuliaFEM/), [SciML](https://github.com/sciml) because their methods might be useful.