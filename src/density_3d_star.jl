"""
Represents a stellar density profile in 3D,
"""

@kwdef struct StellarProfile3D
    "1D averaged velocity dispersion"
    sigma_vx::F 

    "the stellar density"
    rho::Vector{F}

    "error in the density"
    rho_err::Vector{F}

    "the log of the radius of the bin"
    log_r::Vector{F}

    "the log radii used for bins"
    log_r_bins::Vector{F}

    "the mass inside the shell defined by the log_r_bins"
    mass_in_shell::Vector{F}
    mass_in_shell_err::Vector{F}

    quantiles::Vector{F}
    r_quantile::Vector{F}
    time::F = NaN
end


function StellarProfile3D(snap; bins=nothing, quantiles=[0.01, 0.03, 0.1, 0.5, 0.9, 0.99])
    log_r_snap = log10.(calc_r(snap))

    bins, hist, err = histogram(log_r_snap, bins, weights=snap.weights)
    log_r_bins = bins
    mass_in_shell = hist 
    mass_in_shell_err  = err

    Nb = length(mass_in_shell)
    counts = zeros(Nb)

    for i in 1:length(Nb)
        filt = log_r_snap .>= log_r_bins[i]
        filt .&= log_r_snap .< log_r_bins[i+1]
        w = snap.weights[filt]
        counts[i] = effective_sample_size(w)
    end
        

    log_r = midpoints(log_r_bins)
    r_bins = 10 .^ log_r_bins
    r = 10 .^ log_r
    rel_err = mass_in_shell_err ./ mass_in_shell

    rho = calc_ρ_from_hist(r_bins, mass_in_shell)
    rho_err = rho .* rel_err # TODO: does not include binning error


    M_in = cumsum(mass_in_shell)
    M_in_err = 1 ./ sqrt.(cumsum(counts)) .* M_in

    σ_vx = calc_σv_x(snap)
    σ_v_1d = calc_σv_1d(snap)
    @info "σ_vx = $σ_vx, σ_v_1d = $σ_v_1d"

    r_quantile = find_r_quantile_star(snap, quantiles)
    return StellarProfile3D(
        sigma_vx=σ_v_1d, 
        rho=rho,
        rho_err=rho_err,
        log_r=log_r,
        log_r_bins=log_r_bins,
        mass_in_shell=mass_in_shell,
        mass_in_shell_err=mass_in_shell_err,
        quantiles=quantiles,
        r_quantile=r_quantile,
        time = snap.time,
    )
end



function scale(prof::StellarProfile3D, r_scale, v_scale, m_scale)

    return StellarProfile3D(
        sigma_vx=prof.sigma_vx * v_scale,
        rho=prof.rho * m_scale / r_scale^3,
        rho_err=prof.rho_err * m_scale / r_scale^3,
        log_r=prof.log_r .+ log10(r_scale),
        log_r_bins=prof.log_r_bins .+ log10(r_scale),
        mass_in_shell=prof.mass_in_shell * m_scale,
        mass_in_shell_err=prof.mass_in_shell_err * m_scale,
        quantiles=prof.quantiles,
        r_quantile=prof.r_quantile * r_scale,
        time=prof.time,
    )
end


"""
    find_r_quantile_star(snap, p)

Find the radius(ii) containing a given fraction of the total stellar mass of the system
`p` may be a real or a vector of real representing the fractions.
"""
function find_r_quantile_star(snap::Snapshot, p)
    r = calc_r(snap)
    m = snap.weights ./ sum(snap.weights)

    return quantile(r, m, p)
end



"""
    calc_σv_x(snap)

Calculate the velocity dispersion in the x-direction
"""
function calc_σv_x(snap)
	vs = get_v_x(snap)
	w = snap.weights
    return std(vs, w)
end


"""
    calc_σv_1d(snap)

Calculate the average 1D orthoganal velocity dispersion
"""
function calc_σv_1d(snap)
    v = calc_v(snap)
    w = snap.weights
    σ = sqrt( sum(v .^ 2 .* w) / sum(w) )

    return σ / √3
end
