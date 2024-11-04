import Base: @kwdef

import LinearAlgebra: diag, dot, norm, normalize, ⋅, × 

import TOML



F = Float64


"""
An observed 2D density profile
"""
@kwdef mutable struct StellarProfile
    r_units::String
    distance::F = NaN
    mass_units::String = ""
    sigma_v::F = NaN

    log_r::Vector{F}
    log_r_bins::Vector{F}
    counts::Vector{F}

    mass_in_annulus::Vector{F}
    mass_in_annulus_err::Vector{F}

    M_in::Vector{F}
    M_in_err::Vector{F}

    Sigma::Vector{F}
    Sigma_err::Vector{F}
    Sigma_m::Vector{F}
    Sigma_m_err::Vector{F}
    log_Sigma::Vector{F}
    log_Sigma_err::Vector{F}

    Gamma::Vector{F}
    Gamma_err::Vector{F}
    Gamma_max::Vector{F}
    Gamma_max_err::Vector{F}

    time::F = NaN
end



"""
    StellarProfile(radii; arguments...)

Calculate the properties of a density profile given the radii `rs` and the units of the radii `r_units`.


# Arguments
- `weights`: The weights of the radii.
- `bins`: Passed to Arya.histogram. Bins in log r for histogram.
- `normalization`: The normalization of the profile.
    - :mass: normalizes the profile by the total mass.
    - :central: normalizes the profile by the mass within a central radius, r_centre
    - :none: no normalization.
- `r_centre`: The central radius for normalization if `normalization=:central`.
- `distance`: distance (optional)
- `sigma_v`: velocity dispersion (just as an annotation)
- `r_units`: The units of the radii, entirely for self-documentation currently.
"""
function StellarProfile(rs; 
        weights=nothing, 
        bins=nothing, 
        normalization=:mass,
        r_centre=30,
        r_units="", 
        distance=NaN,
        sigma_v=NaN,
    )

    if weights === nothing
        weights = ones(length(rs))
    end

    if any(rs .< 0)
        throw(DomainError("Radii must be positive."))
    end

    if length(rs) < 2
        throw(ArgumentError("Radii must have at least 2 elements."))
    end
    log_r_bin, values, err = histogram(log10.(rs), bins, weights=weights, normalization=:none) # counting histograms
    err[isnan.(err)] .= 0
    log_r = midpoints(log_r_bin)
    δ_log_r = diff(log_r_bin) ./ 2


    mass_per_annulus = values .± err
    _, counts, _ = histogram(log10.(rs), log_r_bin, normalization=:none)

    # approximate error for no particles in bin
    err[isnan.(err)] .= 1

    err[counts .== 0] .= 1

    if normalization == :mass
        mass_per_annulus = mass_per_annulus ./ sum(value.(mass_per_annulus))
    elseif normalization == :central
        filt_cen = log_r .< log10(r_centre)
        Σ_m_cen = sum(value.(mass_per_annulus[filt_cen])) ./ (π * r_centre^2)
        mass_per_annulus = mass_per_annulus ./ Σ_m_cen
    elseif normalization == :none
        mass_per_annulus = mass_per_annulus
    else
        error("normalization not implemented: $normalization")
    end

    Σ = calc_Σ_from_hist(log_r_bin, mass_per_annulus)
    M_in = cumsum(mass_per_annulus)
    Σ_m = calc_Σ_mean_from_hist(log_r_bin, M_in)
    Γ = calc_Γ(log_r, Σ)
    Γ_max = calc_Γ_max(Σ, Σ_m)

    log_Σ = log10.(Σ)
    log_Σ[Σ .== 0] .= NaN


    prof = StellarProfile(
        r_units = r_units,
        distance = distance,
        sigma_v = sigma_v,
        log_r = value.(log_r),
        log_r_bins = log_r_bin,
        counts = counts,
        M_in = value.(M_in),
        M_in_err = nan_uncertainty(M_in),
        mass_in_annulus = value.(mass_per_annulus),
        mass_in_annulus_err = nan_uncertainty(mass_per_annulus),
        Sigma = value.(Σ),
        Sigma_err = nan_uncertainty(Σ),
        Sigma_m = value.(Σ_m),
        Sigma_m_err = nan_uncertainty(Σ_m),
        log_Sigma = value.(log_Σ),
        log_Sigma_err = nan_uncertainty(log_Σ),
        Gamma = value.(Γ),
        Gamma_err = nan_uncertainty(Γ),
        Gamma_max = value.(Γ_max),
        Gamma_max_err = nan_uncertainty(Γ_max)
    )

    return prof
end


function StellarProfile(snap::Snapshot;
        r_units = "kpc",
        eye_vector=[0, 0, 1],
        kwargs...
    )

    if r_units == "kpc"
        x_vec = normalize([1, 0, 0] × eye_vector)
        y_vec = normalize(x_vec × eye_vector)

        N = size(snap.positions, 2)
        ra = [dot(snap.positions[:, i], x_vec) for i in 1:N]
        dec = [dot(snap.positions[:, i], y_vec) for i in 1:N]

        weights = snap.weights
        ra0 = snap.x_cen ⋅ x_vec
        dec0 = snap.x_cen ⋅ y_vec
    else 
        projected = to_gaia(snap, add_centre=true)
        ra = projected.ra
        dec = projected.dec
        ra0 = projected.ra[1]
        dec0 = projected.dec[1]
        weights = projected.weights
    end

    xi, eta = to_tangent(ra, dec, ra0, dec0)
    r = sqrt.(xi .^ 2 + eta .^ 2)

    return StellarProfile(r; r_units=r_units, weights=weights, kwargs...)
end


function Base.print(io::IO, prof::StellarProfile)
    TOML.print(io, struct_to_dict(prof))
end

function StellarProfile(filename::String; kwargs...)
    t = TOML.parsefile(filename)
    t = merge(t, kwargs)
    t = dict_to_tuple(t)

    return StellarProfile(;t...)
end



"""
    nan_uncertainty(v)

Returns the uncertainty of a value `v` with zero uncertainties as NaN.
"""
function nan_uncertainty(v::AbstractVector{<:Measurement})
    x = uncertainty.(v)
    x[x .== 0] .= NaN
    return x
end



"""
    ellipticity_to_aspect(ellipticity)

Converts the ellipticity to the aspect ratio (b/a) of an ellipse.
"""
function ellipticity_to_aspect(ellipticity::Real)
    return 1 - ellipticity
end



# Density methods
"""
    calc_Σ_from_hist(log_r_bin, mass_per_annulus)

Calculate the surface density given the radii `log_r_bin` and the mass per annuli `mass_per_annulus`.
"""
function calc_Σ_from_hist(log_r_bin::AbstractVector{<:Real}, mass_per_annulus::AbstractVector{<:Real})
    r = 10. .^ log_r_bin
    As = π * diff(r .^ 2)

    Σ = mass_per_annulus ./ As 
	return Σ 
end


"""
    calc_Σ_from_1D_density(log_r_bin, density1d)

Calculate the surface density given the radii `log_r_bin` and the 1D density `density1d`.
"""
function calc_Σ_from_1D_density(log_r_bin::AbstractVector{<:Real}, density1d::AbstractVector{<:Real})
    r = 10 .^ log_r_bin
    As = π * diff(r .^ 2)

    Σ = density ./ (2π * log(10) * r .^ 2)
    return Σ
end


"""
    calc_Γ(log_rs, Σs)

Calculate the logarithmic slope of the density profile given the radii `log_rs` and the surface densities `Σs`.
"""
function calc_Γ(log_rs::AbstractVector{<:Real}, Σs::AbstractVector{<:Real})
    d_log_r = gradient(log_rs)
    d_log_Σ = gradient(log10.(Σs))

    return d_log_Σ ./ d_log_r
end


"""
    calc_Σ_mean_from_hist(log_r_bin, M_in)

Calculates the mean surface density from the interior mass to each bin
"""
function calc_Σ_mean_from_hist(log_r_bin::AbstractVector{<:Real}, M_in::AbstractVector{<:Real})
    r = 10 .^ log_r_bin[2:end]
    Areas = @. π * r^2
    Σ_bar = M_in ./ Areas
    return Σ_bar
end


"""
Given the surface density Σ and the mean surface density Σ_m, calculates the maximum slope of the density profile.
""" 
function calc_Γ_max(Σ, Σ_m)
    Γ_max = @. 2*(1 - Σ / Σ_m)
    Γ_max[value.(Σ) .== 0] .= NaN
    return Γ_max
end



# Utility functions


"""
    calc_r_ell_sky(ra, dec, a, b, PA; weights=nothing, centre="mean", units="arcmin")

Given a set of sky coordinates (ra, dec), computes the elliptical radius of each point with respect to the centre of the ellipse defined by the parameters (a, b, PA).
"""
function calc_r_ell_sky(ra, dec, a, b, PA; weights=nothing,
        centre="mean",
        units="arcmin"
    )
    ra0, dec0 = calc_centre2D(ra, dec, centre, weights)

    x, y = to_tangent(ra, dec, ra0, dec0)

    r_ell = calc_r_ell(x, y, a, b, PA)


    if units == "arcmin"
        r_ell = 60r_ell
    elseif units == "arcsec"
        r_ell = 3600r_ell
    elseif units == "deg"
        r_ell = 1r_ell
    else
        error("units not implemented: $units")
    end

    return r_ell
end


function calc_r_ell_sky(ra, dec, ell, PA; kwargs...)
    aspect = ellipticity_to_aspect(ell)
    b = sqrt(aspect)
    a = 1/b
    return calc_r_ell_sky(ra, dec, a, b, PA; kwargs...)
end


function calc_r_ell_sky(ra, dec; kwargs...)
    return calc_r_ell_sky(ra, dec, 0, 0; kwargs...)
end


"""
    calc_r_ell(x, y, a, [b, ]PA)

computes the elliptical radius of a point (x, y) with respect to the center (0, 0) and the ellipse parameters (a, b, PA).
If using sky coordinates, x and y should be tangent coordinates.

Note that the position angle is the astronomy definition, i.e. measured from the North to the East (clockwise) in xi / eta.
"""
function calc_r_ell(x, y, args...)
    x_p, y_p = shear_points_to_ellipse(x, y, args...)

	r_sq = @. (x_p)^2 + (y_p)^2
	return sqrt.(r_sq)
end



"""
    shear_points_to_ellipse(x, y, a, b, PA)

Transforms x and y into the sheared rotated frame of the ellipse.
Position angle is measured from y axis in the direction of positive x.
(North to East on sky).
"""
function shear_points_to_ellipse(x, y, a, b, PA)
    @assert_same_size x y

    if a <= 0 || b <= 0
        throw(DomainError("a and b must be positive. Got a = $a, b = $b."))
    end

    # rotate
    θ = @. deg2rad(PA - 90)
    x_p = @. x * cos(θ) + -y * sin(θ)
    y_p = @. x * sin(θ) + y * cos(θ)

    # scale
    x_p = x_p ./ a
    y_p = y_p ./ b

    return x_p, y_p
end


function shear_points_to_ellipse(x, y, ell, PA)
    aspect = ellipticity_to_aspect(ell)
    b = sqrt(aspect)
    a = 1/b
    return shear_points_to_ellipse(x, y, a, b, PA)
end




"""
    spherical_mean(ra, dec, weights=nothing)

Calculates the spherical mean of a set of sky coordinates (ra, dec) with optional weights. The mean is calculated by taking the average of the unit vectors of the points, so avoids problems with periodicity and spherical geometry.
"""
function spherical_mean(ra::AbstractVector{<:Real}, dec::AbstractVector{<:Real}, weights::Union{Nothing, AbstractVector{<:Real}}=nothing)
    if length(ra) != length(dec)
        throw(DimensionMismatch("ra and dec must have the same length."))
    end

    pos = unit_vector(ra, dec)
    if weights === nothing
        mean_pos = centroid(pos)
    else
        mean_pos = centroid(pos, weights)
    end

    ra0, dec0, r = cartesian_to_sky(mean_pos...)
    if r == 0
        throw(DomainError("Mean position is at the origin."))
    end
    return ra0, dec0
end



"""
    calc_centre2D(ra, dec, centre_method, weights=nothing)

Calculates the centre of vectors of ra and dec (with weights) using one of the
following methods:
- "mean": calculates the mean of the vectors
- Tuple{ra0, dec0}: uses the given ra0 and dec0 as the centre (bypasses calculation)
This is a utility function.

"""
function calc_centre2D(ra::AbstractVector{<:Real}, dec::AbstractVector{<:Real}, centre_method, weights=nothing)
    if centre_method == "mean"
        ra0, dec0 = spherical_mean(ra, dec, weights)
    elseif centre_method isa Tuple
        ra0, dec0 = centre_method
    else
        throw(ArgumentError("centre method not implemented: $centre_method"))
    end

    return ra0, dec0
end



"""
    to_orbit_coords(ra, dec, ra0, dec0, PA) 

Given the position angle of an orbit vector at (ra0, dec0) calculates a rotated
sky frame centred on RA, DEC and with the x-axis aligned with the orbit.
The orbital position angle can be found in the orbit properties file from
analyze_orbit.jl notebook.
"""
function to_orbit_coords(ra, dec, ra0::Real, dec0::Real, PA::Real)
	# want to rotate to dec, ra of centre, then rotate 
	α = deg2rad(ra0)
	δ = deg2rad(dec0)
	ϖ = deg2rad(90 - PA)
	Rmat = Rx_mat(-ϖ) * Ry_mat(δ) * Rz_mat(-α) 

	coords = unit_vector(ra, dec)
	coords =  Rmat * coords
	ra, dec, _ = cartesian_to_sky(coords[1, :], coords[2, :], coords[3, :])

	ra .-= 360 * (ra .> 180)

	ra, dec
end



function to_orbit_coords(ra::Real, dec::Real, ra0::Real, dec0::Real, PA::Real)
    ra, dec =  to_orbit_coords([ra], [dec], ra0, dec0, PA)
    return ra[1], dec[1]
end


function scale(prof::StellarProfile, r_scale, m_scale)
    return StellarProfile(
        r_units = prof.r_units,
        distance = prof.distance,
        log_r = prof.log_r .+ log10(r_scale),
        log_r_bins = prof.log_r_bins .+ log10(r_scale),
        counts = prof.counts,
        M_in = prof.M_in * m_scale,
        M_in_err = prof.M_in_err * m_scale,
        mass_in_annulus = prof.mass_in_annulus * m_scale,
        mass_in_annulus_err = prof.mass_in_annulus_err * m_scale,
        Sigma = prof.Sigma * m_scale / r_scale^2,
        Sigma_err = prof.Sigma_err * m_scale / r_scale^2,
        Sigma_m = prof.Sigma_m * m_scale / r_scale^2,
        Sigma_m_err = prof.Sigma_m_err * m_scale / r_scale^2,
        log_Sigma = prof.log_Sigma .+ log10(m_scale) .- 2log10(r_scale),
        log_Sigma_err = prof.log_Sigma_err .+ log10(m_scale) .- 2log10(r_scale),
        Gamma = prof.Gamma,
        Gamma_err = prof.Gamma_err,
        Gamma_max = prof.Gamma_max,
        Gamma_max_err = prof.Gamma_max_err
    )
end
