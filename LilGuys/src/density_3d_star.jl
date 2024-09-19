
@kwdef struct StellarProfile3D
    sigma_vx::F 

    rho::Vector{F}
    rho_err::Vector{F}
    log_r::Vector{F}
    log_r_bins::Vector{F}

    mass_in_shell::Vector{F}
    mass_in_shell_err::Vector{F}
end


function StellarProfile3D(snap; bins=nothing)
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
    return StellarProfile3D(
        sigma_vx=σ_vx, 
        rho=rho,
        rho_err=rho_err,
        log_r=log_r,
        log_r_bins=log_r_bins,
        mass_in_shell=mass_in_shell,
        mass_in_shell_err=mass_in_shell_err,
    )
end




"""
    find_r_containing_Mstar(snap, p)

Find the radius(ii) containing a given fraction of the total stellar mass of the system
"""
function find_r_containing_Mstar(snap, p)
    r = lguys.calc_r(snap)
    m = snap.weights ./ sum(snap.weights)

    return quantile(r, m, p)
end



function calc_σv_x(snap)
	vs = get_v_x(snap)
	w = snap.weights
    return std(vs, w)
end
