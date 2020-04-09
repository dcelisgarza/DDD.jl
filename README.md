# DDD

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dcelisgarza.github.io/DDD.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dcelisgarza.github.io/DDD.jl/dev)
[![Build Status](https://travis-ci.com/dcelisgarza/DDD.jl.svg?branch=master)](https://travis-ci.com/dcelisgarza/DDD.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/dcelisgarza/DDD.jl?svg=true)](https://ci.appveyor.com/project/dcelisgarza/DDD-jl)
[![Codecov](https://codecov.io/gh/dcelisgarza/DDD.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dcelisgarza/DDD.jl)
[![Coveralls](https://coveralls.io/repos/github/dcelisgarza/DDD.jl/badge.svg?branch=master)](https://coveralls.io/github/dcelisgarza/DDD.jl?branch=master)

New generation of 3D Discrete Dislocation Dynamics codes.

Dislocation dynamics is a complex field with an enormous barrier to entry. The aim of this project is to create a codebase that is:

- Easy to use.
- Easy to maintain.
- Easy to develop for.
- Modular.
- Idiot proof.
- Well documented and tested.
- Performant.
- Easily parallelisable.

## Current TODO:
- [ ] Custom 3-vec type, place x,y,z coordinates in contiguous memory instead of columns, ie [x1 y1 z1; x2 y2 z2] -> [x1;y1;z1;x2;y2;z2], have to define custom array type, `getindex(arr, (a,b)) = arr[3*(a-1)+b]`, out of bounds and all the rest. Watch [this](https://www.youtube.com/watch?v=jS9eouMJf_Y).
- [x] Generate docs
  - [x] Documented Misc
  - [ ] Upload docs
- [ ] Optimise
- [ ] Specialised integrator
  - [ ] Perhaps later make use of [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) for their stepping and event handling routines.
- [ ] Calculate segment-segment interactions.
- [ ] Mobility laws
  - [ ] BCC
  - [ ] FCC
- [ ] Topology operations
  - [ ] Split
  - [ ] Merge
- [ ] Couple to FEM, perhaps use a package from [JuliaFEM](http://www.juliafem.org/).
  - [ ] Boundary conditions
    - [ ] Neuman
    - [ ] Dirichlet
  - [ ] Displacements
  - [ ] Tractions
