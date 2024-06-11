import Base: @kwdef

import LinearAlgebra: diag

import LilGuys as lguys

import Arya: value, err, histogram
import StatsBase: weights, mean
import TOML


# import Polyhedra

using Measurements

F = Float64

@kwdef mutable struct ObsProfile
    log_r::Vector{F}
    log_r_bins::Vector{F}
    r_units::String
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
end


"""
general method to convert a struct to a dictionary
"""
function struct_to_dict(S)
    return Dict(key=>getfield(S, key) for key in fieldnames(typeof(S)))
end


function dict_to_tuple(D)
    return NamedTuple((Symbol(key), value) for (key, value) in D)
end


function Base.print(io::IO, prof::ObsProfile)
    TOML.print(io, struct_to_dict(prof))
end


function ObsProfile(filename)
    t = dict_to_tuple(TOML.parsefile(filename))
    return ObsProfile(;t...)
end





"""
    r_ell_max(xi, eta, ecc, PA)

returns the maximum radius of an ellipse within the convex boundary defined by the points.

Not implemented
"""
function r_ell_max(xi, eta, ecc, PA)
    error("not implemented")
end



function convex_hull(xi, eta)
    ps = zip(xi, eta)
    p = Polyhedra.convexhull(ps...)
    b = Polyhedra.planar_hull(p)
    return b.points.points
end



function mean_centre(ra, dec, mass)
	return mean(ra, weights(mass)), mean(dec, weights(mass))
end



function calc_centre2D(ra, dec, mass, centre_method)
	if centre_method == "mean"
		ra0, dec0 = mean_centre(ra, dec, mass)
    else
        error("centre method not implemented: $centre_method")
	end

	return ra0, dec0
end



# Density methods


# Σ = h.values ./ (2π * log(10) * r .^ 2) # is equivalent because of derivative of log r
"""
    calc_Σ(log_r_bin, mass_per_annulus)

Calculate the surface density given the radii `log_r_bin` and the mass per annuli `mass_per_annulus`.
"""
function calc_Σ(log_r_bin, mass_per_annulus)
    r = 10 .^ log_r_bin
    As = π * diff(r .^ 2)

    Σ = mass_per_annulus ./ As 
	return Σ 
end


"""
    calc_Γ(log_rs, Σs)

Calculate the logarithmic slope of the density profile given the radii `log_rs` and the surface densities `Σs`.
"""
function calc_Γ(log_rs, Σs)
	d_log_r = lguys.gradient(log_rs)
	d_log_Σ = lguys.gradient(log10.(Σs))

	return d_log_Σ ./ d_log_r
end


function calc_M_in(log_r_bin, mass_per_annulus)
    M_in = cumsum(mass_per_annulus)
    return M_in 
end

function calc_Σ_mean(log_r_bin, M_in)
    r = 10 .^ log_r_bin[2:end]
	Areas = @. π * r^2
	Σ_bar = M_in ./ Areas
	return Σ_bar
end



function calc_Γ_max(Σ, Σ_m)
    Γ_max = @. 2*(1 - Σ / Σ_m)
    Γ_max[value.(Σ) .== 0] .= NaN
    return Γ_max
end


"""
    calc_properties(rs, r_units; weights=nothing, bins=20, normalization="mass")

Calculate the properties of a density profile given the radii `rs` and the units of the radii `r_units`.

"""
function calc_properties(rs; r_units="", weights=nothing, bins=20, normalization=true)
    if weights === nothing
        weights = ones(length(rs))
    end


    h = histogram(log10.(rs), bins, weights=weights, normalization=:count)
    log_r_bin = h.bins
    h.err[isnan.(h.err)] .= 0

    mass_per_annulus = h.values .± h.err
    counts = h.values

    if normalization
        mass_per_annulus = mass_per_annulus ./ sum(mass_per_annulus)
    end

    log_r = lguys.midpoint(log_r_bin)
    δ_log_r = diff(log_r_bin) ./ 2
    log_r = log_r

    Σ = calc_Σ(log_r_bin, mass_per_annulus)

    M_in = calc_M_in(log_r_bin, mass_per_annulus)
    Σ_m = calc_Σ_mean(log_r_bin, M_in)

    Γ = calc_Γ(log_r, Σ)
    Γ_max = calc_Γ_max(Σ, Σ_m)


	log_Σ = log10.(Σ)
    log_Σ[Σ .== 0] .= NaN


    prof = ObsProfile(
        r_units = r_units,
        log_r = value.(log_r),
        log_r_bins = log_r_bin,
        counts = counts,
        M_in = value.(M_in),
        M_in_err = err.(M_in),
        mass_in_annulus = value.(mass_per_annulus),
        mass_in_annulus_err = err.(mass_per_annulus),
        Sigma = value.(Σ),
        Sigma_err = err.(Σ),
        Sigma_m = value.(Σ_m),
        Sigma_m_err = err.(Σ_m),
        log_Sigma = value.(log_Σ),
        log_Sigma_err = err.(log_Σ),
        Gamma = value.(Γ),
        Gamma_err = err.(Γ),
        Gamma_max = value.(Γ_max),
        Gamma_max_err = err.(Γ_max)
    )

    return prof
end
