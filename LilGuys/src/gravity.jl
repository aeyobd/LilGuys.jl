import Interpolations: linear_interpolation, Line


"""
given a centered snapshot, returns a interpolated potential a a function 
of r
"""
function calc_radial_Φ(positions::Matrix{T}, masses) where T <: Real
    radii = calc_r(positions)
    # work inside out
    idx = sortperm(rs)
    rs = rs[idx]
    ms = masses[idx]

    N = length(ms)
    M_in = cumsum(ms)

    Φs_out = zeros(F, N)
    for i in 2:N
        Φs_out[i-1] = calc_Φ(m[i:end], rs[i:end])
    end
    Φ_cen = calc_Φ(m, rs)

    return r -> _interpolated_Φ(r, rs, Ms_in, Φs_out, Φ_cen)
end

function _interpolated_Φ(r, rs, Ms_in, Φs_out, Φ_cen)
    if r < rs[1]
        return Φ_cen
    end
    idx = searchsortedlast(rs, r)
    Φ_in = calc_Φ(M_in[idx], r)
    Φ = Φ_out[idx] + Φ_in
    return Φ
end



"""
Gravitational potential due to particles in snapshot at position x_vec
"""
function calc_Φ(snap::Snapshot, x_vec)
    return calc_Φ(snap.x_vec, snap.m, x_vec)
end


"""
Gravitational potential due to ensemble of masses at positions evaluated at x_vec
"""
function calc_Φ(positions::Matrix{T}, masses::Vector{T}, x_vec) where T <: Real
    radii = calc_r(positions .- x_vec)
    return calc_Φ(masses, radii)
end



"""
Potential due to a collection of masses at given radii, or equivalently
potential inside centred shells of masses and radii
"""
function calc_Φ(masses::Vector{T}, radii::Vector{T}) where T <: Real
    return sum(Φ_point.(masses, radii))
end


"""
One point potential law (-Gm/r)
"""
function calc_Φ(mass::Real, radius::Real)
    if R == 0
        return -Inf
    end
    return -G * mass / radius
end


"""
Force of gravity from masses at given positions evaluated at x_vec
"""
function calc_F_grav(positions::Matrix{F}, masses::Vector{F}, x_vec)
    dr = x_vec .- positions
    rs = calc_r(dr)
    r_hat = dr ./ rs
    force =  sum(calc_F_point.(masses, rs) .* r_hat, dims=2)
    force[:, rs .== 0] .= 0  # remove divergences
    return force
end


function calc_F_grav(snap::Snapshot, x_vec)
    return calc_F_grav(snap.x_vec, snap.m, x_vec)
end


function calc_F_point(mass::Real, radius::Real)
    if radius == 0
        return 0
    end
    return -G * mass / radius^2
end
