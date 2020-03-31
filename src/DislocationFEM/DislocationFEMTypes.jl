struct DislocationFEMCorrective{T1 <: AbstractArray{<:Float64, N} where {N}}
    uHat::T1 # U_hat
    fHat::T1 # F_hat
end
