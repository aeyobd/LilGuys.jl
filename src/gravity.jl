# This file containes methods to evaluate the gravitational 
# potential and force



struct DistributionFunction{F} 
    ρ::Vector{F}
    ψ::Vector{F}
    d2ρ_dψ2::LinInterp
end


"""
    DistributionFunction(ρ, ψ, r)

Given three vectors for the density, potential and radius,
creates a distribution function object by taking gradients wrt r.
Calling the object with a binding energy ϵ returns the distribution function at
that energy.
"""
function DistributionFunction(ρ::AbstractVector{<:Real}, ψ::AbstractVector{<:Real}, r::AbstractVector{<:Real}; force_positive=false)

    if !issorted(r) | !issorted(-ψ)
        throw(ArgumentError("arrays must be sorted"))
    end

    ρ1 = gradient(ρ, r)
    ρ2 = gradient(ρ1, r)
    ψ1 = gradient(ψ, r)
    ψ2 = gradient(ψ1, r)

    d2ρ_dψ2 = @. ψ1^-2 * ρ2 - ψ1^-3 * ρ1 * ψ2

    if force_positive
        d2ρ_dψ2 = abs.(d2ρ_dψ2)
    end
    d2_interp = lerp(reverse(ψ), reverse(d2ρ_dψ2))

    F = Float64 # TODO: make this more general
    return DistributionFunction{F}(ρ, ψ, d2_interp)
end


"""
Returns the distribution function at a given energy ϵ
"""
function (df::DistributionFunction)(ϵ)
    if df.ψ[1] < ϵ || ϵ < df.ψ[end]
        throw(DomainError("ϵ must be between the minimum and maximum value of ψ: $(df.ψ[1]) < ϵ < $(df.ψ[end]); got $ϵ"))
    end

    f_integrand(ψ) = 1/(√8*π^2) * df.d2ρ_dψ2(ψ) /√(ϵ - ψ)
    f = integrate(f_integrand, 0, ϵ)

    return f
end



"""
given a centered snapshot, returns a interpolated potential a a function 
of r
"""
function calc_radial_Φ(radii::AbstractVector{T}, masses::AbstractVector) where {T <: Real}
    # work inside out
    idx = sortperm(radii)
    rs_sorted = radii[idx]
    ms_sorted = masses[idx]

    N = length(ms_sorted)
    Ms_in = cumsum(ms_sorted)

    Φs_out = zeros(T, N)

    for i in 1:N-1
        Φs_out[i] = calc_Φ(rs_sorted[i+1:end], ms_sorted[i+1:end])
    end
    Φ_cen = calc_Φ(rs_sorted, ms_sorted)

    return r -> _interpolated_Φ(r, rs_sorted, Ms_in, Φs_out, Φ_cen)
end


function calc_radial_Φ(positions::AbstractMatrix{<:Real}, masses::AbstractVector{<:Real}) 
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
    Φ_in = calc_Φ(r, Ms_in[idx])
    Φ = Φs_out[idx] + Φ_in
    return Φ
end




"""
    calc_discrete_spherical_Φ(positions, masses)

Given a collection of masses at given radii,
returns the potential at each radius.

The potential is calculated as
```math
Φ(r) = -G M(r) / r - \\int_r^\\infty G dm/dr(r') / r' dr'
```
"""
function calc_radial_discrete_Φ(radii::AbstractVector{T}, masses::AbstractVector) where T <: Real
    # work inside out
    idx = sortperm(radii)
    rs_sorted = radii[idx]
    ms_sorted = masses[idx]

    M_in = cumsum(ms_sorted)

    Φ_shells = calc_Φ.(rs_sorted, ms_sorted)[2:end]

    # include each shell outside the current one...
    Φ_out = Vector{T}(undef, length(masses))
    Φ_out[1:end-1] .= reverse(cumsum(reverse(Φ_shells)))
    Φ_out[end] = 0

    Φ_in = calc_Φ.(rs_sorted, M_in)

    Φ = Φ_out .+ Φ_in
    return Φ[invperm(idx)]
end


function calc_radial_discrete_Φ(positions::Matrix{T}, masses::AbstractVector) where T <: Real
    radii = calc_r(positions)
    return calc_radial_discrete_Φ(radii, masses)
end


function calc_radial_discrete_Φ(snap::Snapshot)
    return calc_radial_discrete_Φ(calc_r(snap), snap.masses)
end



"""
    calc_Φ(snap, x_vec)

Gravitational potential due to particles in snapshot at position x_vec
"""
function calc_Φ(snap::Snapshot, x_vec)
    return calc_Φ(snap.positions, snap.masses, x_vec)
end

function calc_Φ(snap::Snapshot)
    Φ = Vector{F}(undef, length(snap))

    N = length(snap)
    for i in 1:N
        j = 1:N .!= i
        positions = snap.positions[:, j]
        masses = snap.masses[j]

        Φ[i] = calc_Φ(positions, masses, snap.positions[:, i])
    end

    return Φ
end


"""
    calc_Φ(masses, positions, x_vec)

Gravitational potential due to ensemble of masses at positions evaluated at x_vec
"""
function calc_Φ(positions::Matrix{T}, masses::Vector{T}, x_vec) where T <: Real
    radii = calc_r(positions .- x_vec)
    return calc_Φ(radii, masses)
end



"""
    calc_Φ(radii, masses)

Potential due to a collection of masses at given radii, or equivalently
potential inside centred shells of masses and radii
"""
function calc_Φ(radii::Vector{T}, masses::Vector{T}) where T <: Real
    return sum(calc_Φ.(radii, masses))
end


"""
    calc_Φ(radius, mass)

One point potential law (-Gm/r)
"""
function calc_Φ(radius::Real, mass::Real)
    if radius == 0
        return -Inf
    end
    return -G * mass / radius
end



"""
    calc_F_grav(positions, masses, x_vec)

Force of gravity from masses at given positions evaluated at x_vec
"""
function calc_F_grav(positions::AbstractMatrix, masses::AbstractVector, x_vec)
    dr = x_vec .- positions
    rs = calc_r(dr)
    r_hat = dr ./ rs
    force =  sum(calc_F_point.(rs, masses) .* r_hat, dims=2)
    force[:, rs .== 0] .= 0  # remove divergences
    return force
end



function calc_F_grav(snap::Snapshot, x_vec)
    return calc_F_grav(snap.masses, snap.positions, x_vec)
end


"""
    calc_F_point(radius, mass)

Force of gravity from a point mass at given radius
    
F = -G m / r^2
"""
function calc_F_point(radius::Real, mass::Real)
    if radius == 0
        return 0
    end
    return -G * mass / radius^2
end
