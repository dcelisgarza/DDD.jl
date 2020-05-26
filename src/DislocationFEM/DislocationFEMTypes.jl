"""
```
mutable struct DislocationFEMCorrective{T1 <: AbstractArray{T2, N} where {T2, N}}
    uHat::T1 # U_hat
    fHat::T1 # F_hat
end
```
Corrective stress and force.
"""
mutable struct DislocationFEMCorrective{T1 <: AbstractArray{T, N} where {T, N}}
    uHat::T1 # U_hat
    fHat::T1 # F_hat
end
