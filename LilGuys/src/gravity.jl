import Interpolations: linear_interpolation, Line


"""
given a centered snapshot, returns a interpolated potential a a function 
of r
"""
function calc_radial_Φ(snap::Snapshot)
    rs = calc_r(snap.pos)
    # work inside out
    idx = sortperm(rs)
    rs = rs[idx]
    m = snap.m[idx]

    N = length(snap)
    M_in = 0
    Φs = zeros(N)

    for i in 1:N
        Φ_out = calc_Φ(m[i:end], rs[i:end])
        Φ_in = calc_Φ(M_in, rs[i])
        Φs[i] = Φ_out + Φ_in

        M_in += m[i]
    end


    lerp =  linear_interpolation(rs, Φs, extrapolation_bc=Line())

    function Φ(r)
        if r < rs[1]
            return lerp(r)
        else
            return calc_Φ(M_in, r)
        end
    end

    return Φ
end

function calc_Φ(snap::Snapshot, pos)
    return calc_Φ(snap.pos, snap.m, pos)
end

function calc_Φ(pos::Matrix{T}, m::Vector{T}, x) where T <: Real
    rs = calc_r(pos .- x)
    return calc_Φ(m, rs)
end

function calc_Φ(m::Vector{T}, rs::Vector{T}) where T <: Real
    return -G * sum(m ./ rs)
end

function calc_Φ(m::Real, R::Real)
    return -G * m / R
end

function calc_F_grav(snap::Snapshot, pos)
    return calc_F_grav(snap.pos, snap.m, pos)
end

function calc_F_grav(pos::Matrix{F}, m::Vector{F}, x)
    dr = x .- pos 
    rs = calc_r(dr)
    r_hat = dr ./ rs
    f =  -G * sum(m .* r_hat ./ rs.^2, dims=2)
    f[:, rs .== 0] .= 0
    return f
end
