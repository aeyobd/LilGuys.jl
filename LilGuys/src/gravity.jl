# This file containes methods to evaluate the gravitational 
# potential and force



"""
given a centered snapshot, returns a interpolated potential a a function 
of r
"""
function calc_radial_Φ(radii::AbstractVector{T}, masses) where T <: Real
    # work inside out
    idx = sortperm(radii)
    rs_sorted = radii[idx]
    ms_sorted = masses[idx]

    N = length(ms_sorted)
    Ms_in = cumsum(ms_sorted)

    Φs_out = zeros(F, N)
    for i in 2:N
        Φs_out[i-1] = calc_Φ(ms_sorted[i:end], rs_sorted[i:end])
    end
    Φ_cen = calc_Φ(ms_sorted, rs_sorted)

    return r -> _interpolated_Φ(r, rs_sorted, Ms_in, Φs_out, Φ_cen)
end



"""
    calc_discrete_spherical_Φ(radii, masses)

Given a collection of masses at given radii,
returns the potential at each radius.

The potential is calculated as
```math
Φ(r) = -G M(r) / r - \\int_r^\\infty G dm/dr(r') / r' dr'
```
"""
function calc_radial_discrete_Φ(masses::AbstractVector{T}, radii::AbstractVector{T}) where T <: Real
    # work inside out
    idx = sortperm(radii)
    rs_sorted = radii[idx]
    ms_sorted = masses[idx]

    M_in = cumsum(ms_sorted)

    Φ_shells = calc_Φ.(ms_sorted, rs_sorted)[2:end]

    # include each shell outside the current one...
    Φ_out = Vector{T}(undef, length(masses))
    Φ_out[1:end-1] .= reverse(cumsum(reverse(Φ_shells)))
    Φ_out[end] = 0

    Φ_in = calc_Φ.(M_in, rs_sorted)

    return Φ_out .+ Φ_in
end



function calc_radial_Φ(masses::Vector{T}, positions::Matrix{T}) where T <: Real
    radii = calc_r(positions)
    return calc_radial_Φ(radii, masses)
end

function calc_radial_Φ(snap::Snapshot)
    return calc_radial_Φ(snap.positions, snap.masses)
end

function _interpolated_Φ(r, rs, Ms_in, Φs_out, Φ_cen)
    if r < rs[1]
        return Φ_cen
    end
    idx = searchsortedlast(rs, r)
    Φ_in = calc_Φ(Ms_in[idx], r)
    Φ = Φs_out[idx] + Φ_in
    return Φ
end



"""
Gravitational potential due to particles in snapshot at position x_vec
"""
function calc_Φ(snap::Snapshot, x_vec)
    return calc_Φ(snap.masses, snap.positions, x_vec)
end


"""
Gravitational potential due to ensemble of masses at positions evaluated at x_vec
"""
function calc_Φ(masses::Vector{T}, positions::Matrix{T}, x_vec) where T <: Real
    radii = calc_r(positions .- x_vec)
    return calc_Φ(masses, radii)
end



"""
Potential due to a collection of masses at given radii, or equivalently
potential inside centred shells of masses and radii
"""
function calc_Φ(masses::Vector{T}, radii::Vector{T}) where T <: Real
    return sum(calc_Φ.(masses, radii))
end


"""
One point potential law (-Gm/r)
"""
function calc_Φ(mass::Real, radius::Real)
    if radius == 0
        return -Inf
    end
    return -G * mass / radius
end


"""
Force of gravity from masses at given positions evaluated at x_vec
"""
function calc_F_grav(masses::Vector{F}, positions::Matrix{F}, x_vec)
    dr = x_vec .- positions
    rs = calc_r(dr)
    r_hat = dr ./ rs
    force =  sum(calc_F_point.(masses, rs) .* r_hat, dims=2)
    force[:, rs .== 0] .= 0  # remove divergences
    return force
end


function calc_F_grav(snap::Snapshot, x_vec)
    return calc_F_grav(snap.masses, snap.positions, x_vec)
end


function calc_F_point(mass::Real, radius::Real)
    if radius == 0
        return 0
    end
    return -G * mass / radius^2
end
