function IntegrationParameters(; method::T1, tmin::T2, tmax::T2, dtmin::T2 = 1e-3, dtmax::T2 = Inf, abstol::T2 = 1e-6, reltol::T2 = 1e-6, maxchange::T2 = 1.2, exponent::T2 = 20.0, maxiter::T3 = 10) where {T1 <: AbstractIntegrator, T2 <: AbstractFloat, T3 <: Integer}
    return IntegrationParameters(method, tmin, tmax, dtmin, dtmax, abstol, reltol, maxchange, exponent, maxiter)
end
function IntegrationTime(; dt::T1 = 0.0, time::T1 = 0.0, step::T2 = 0) where {T1 <: AbstractFloat, T2 <: Integer}
    return IntegrationTime(dt, time, step)
end