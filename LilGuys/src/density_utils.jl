import Base: @kwdef

import LinearAlgebra: diag

import LilGuys as lguys

import Arya: value, err, histogram
import StatsBase: weights, mean
import TOML
using Measurements

F = Float64

@kwdef mutable struct ObsProfile
    log_r::Vector{F}
    log_r_bins::Vector{F}
    r_units::String
    counts::Vector{F}
    mass::Vector{F}
    mass_err::Vector{F}

    Sigma::Vector{F}
    Sigma_err::Vector{F}
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
calculates the elliptical radii given the ra and dec of stars
"""
function calc_radii(ra, dec; weights=nothing, centre_method="mean", ecc=0, PA=0)
    if weights === nothing
        weights = ones(length(ra))
    end

	ra0, dec0 = calc_centre(ra, dec, weights, centre_method)
	xi, eta  = lguys.to_tangent(ra, dec, ra0, dec0)

	b = sqrt(1 - ecc)
	a = 1/b
	r_ell = lguys.calc_r_ell(xi, eta, a, b, PA-90)

	r_max = sqrt(maximum(xi .^ 2 + eta .^ 2))
    r_cut = r_max * sqrt(1 - ecc)
    filt = r_ell .< r_cut

	Ncut = sum(map(!, filt))
	println("cut $Ncut stars on edge ")

    r_ell, filt
end



function add_r_ell!(stars, ecc, PA)
    r_ell, filt = lguys.calc_r_ell(stars.xi, stars.eta, a, b, PA)
    stars[:, "r_ell"] = r_ell;
    stars = stars[filt, :]
end



function mean_centre(ra, dec, mass)
	return mean(ra, weights(mass)), mean(dec, weights(mass))
end



function calc_centre(ra, dec, mass, centre_method)
	if centre_method == "mean"
		ra0, dec0 = mean_centre(ra, dec, mass)
    else
        error("centre method not implemented: $centre_method")
	end

	return ra0, dec0
end



# Density methods


# Σ = h.values ./ (2π * log(10) * r .^ 2) # is equivalent because of derivative of log r
function calc_Σ(log_r_bin, hist, hist_err)
    As = π * diff((10 .^ log_r_bin) .^ 2)

    Σ = hist ./ As 

    δ_Σ = hist_err ./ As
	return Σ 
end


function calc_Σ_mean(log_r_bin, hist)
    r = 10 .^ lguys.midpoint(log_r_bin)
    counts = cumsum(hist .* diff(log_r_bin))
	Areas = @. π * r^2
	Σ_bar = counts ./ Areas
	return Σ_bar
end


function calc_Γ(log_rs, Σs)
	d_log_r = lguys.gradient(log_rs)
	d_log_Σ = lguys.gradient(log10.(Σs))

	return d_log_Σ ./ d_log_r
end


function calc_M_in(hist)
    M_in = cumsum(hist)
    return M_in
end


function calc_Γ_max(Σ, Σ_m)
    return @. 2*(1 - Σ / Σ_m)
end


"""
    calc_properties(rs, r_units; weights=nothing, bins=20, normalization="mass")

Calculate the properties of a density profile given the radii `rs` and the units of the radii `r_units`.

"""
function calc_properties(rs, r_units; weights=nothing, bins=20, normalization="mass")
    if weights === nothing
        weights = ones(length(rs))
    end
    h = histogram(log10.(rs), bins, weights=weights)
    log_r_bin = h.bins
    δ_hist = h.err
    hist = h.values .± h.err

    log_r = lguys.midpoint(log_r_bin)
    δ_log_r = diff(log_r_bin) ./ 2
    log_r = log_r .± δ_log_r

    Σ = calc_Σ(log_r_bin, hist, δ_hist)
    Σ_m = calc_Σ_mean(log_r_bin, hist)
    Γ = calc_Γ(log_r, Σ)
    Γ_max = calc_Γ_max(Σ, Σ_m)

    M_in = cumsum(hist)

	log_Σ = log10.(Σ)


    prof = ObsProfile(
        r_units = r_units,
        log_r = value.(log_r),
        log_r_bins = log_r_bin,
        counts = value.(hist),
        mass = value.(M_in),
        mass_err= err.(M_in),
        Sigma = value.(Σ),
        Sigma_err = err.(Σ),
        log_Sigma = value.(log_Σ),
        log_Sigma_err = err.(log_Σ),
        Gamma = value.(Γ),
        Gamma_err = err.(Γ),
        Gamma_max = value.(Γ_max),
        Gamma_max_err = err.(Γ_max)
    )

    return prof
end
