function limits!(
    lims::AbstractArray{<:Float64, N1},
    segLen::Float64,
    range::AbstractArray{<:Float64, N2},
    buffer::Float64,
) where {N1, N2}
    @inbounds for i = 1:size(lims, 2)
        for j = 1:size(lims, 1)
            lims[j, i] = range[j, i] + buffer * segLen
        end
    end
    return lims
end

function translatePoints(
    coord::AbstractArray{<:Float64, N1},
    lims::AbstractArray{<:Float64, N1},
    disp::AbstractArray{<:Float64, N2},
) where {N1, N2}
    @inbounds for i = 1:size(coord, 2)
        for j = 1:size(coord, 1)
            coord[j, i] += lims[1, i] + (lims[2, i] - lims[1, i]) * disp[i]
        end
    end
    return coord
end

function loopDistribution(dist::Zeros, n::Integer, args...)
    return zeros(n, 3)
end
function loopDistribution(dist::Rand, n::Integer, args...)
    return rand(n, 3)
end
function loopDistribution(dist::Randn, n::Integer, args...)
    return randn(n, 3)
end
function loopDistribution(dist::Regular, n::Integer, args...)
    error("loopDistribution: regular distribution yet not implemented")
end
