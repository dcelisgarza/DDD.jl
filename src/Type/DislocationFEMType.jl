"""
Corrective stress and force.
```
struct DislocationFEMCorrective{T1 <: AbstractArray{T2, N} where {T2, N}}
    uHat::T1
    fHat::T1
end
```
"""
struct DislocationFEMCorrective{T1 <: AbstractArray{T, N} where {T, N}}
    uHat::T1
    fHat::T1
end