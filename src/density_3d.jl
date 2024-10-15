import TOML
import LinearAlgebra: ×


"""
A struct representing a 3-dimensional 
mass profile (likely of a snapshot).
All properties are in code units.
"""
@kwdef struct MassProfile3D
    """ Total Energy """
    E::F

    """ Total Kinetic energy """
    K::F

    """ Total Potential energy """
    W::F

    """ Total Angular momentum vector wrt centre"""
    L::Vector{F}

    """ Maximum circular velocity """
    v_circ_max::F
    """ Radius of maximum circular velocity"""
    r_circ_max::F

    "Number of bound particles"
    N_bound::Int

    "Log radius of bin middles"
    log_r::Vector{F}
    "Log radius of bin edges"
    log_r_bins::Vector{F}
    "Number of particles in each bin"
    counts::Vector{Int}

    "Mass in each shell"
    mass_in_shell::Vector{F}
    mass_in_shell_err::Vector{F}

    "Mass enclosed within each shell (up to exterior bin edge)"
    M_in::Vector{F}
    M_in_err::Vector{F}

    "Mean density in bin"
    rho::Vector{F}
    "Density error in each shell"
    rho_err::Vector{F}

    "The radius for each value of v_circ"
    r_circ::Vector{F}
    "Circular velocity"
    v_circ::Vector{F}
    "Circular velocity error"
    v_circ_err::Vector{F}
    "The number of particles within each r_circ"
    n_circ::Vector{F}

    "Circular crossing time / dynamical time at bin"
    t_circ::Vector{F}

    "the radius at which q of the total mass is enclosed where q is given in `quantiles`"
    r_quantile::Vector{F}

    "Quantiles used for M_quantile"
    quantiles::Vector{F}

    time::F = NaN
end

function Base.print(io::IO, prof::MassProfile3D)
    TOML.print(io, struct_to_dict(prof))
end

function MassProfile3D(filename::String)
    t = dict_to_tuple(TOML.parsefile(filename))
    return MassProfile3D(;t...)
end


"""
    MassProfile3D(snap; bins=100, filt_bound=true)

Given a snapshot, computes the density, circular velocity, and energy; 
returning an MassProfile3D object.

`bins` may be an integer, array, or function and is passed to `histogram`
"""
function MassProfile3D(snap::Snapshot;
        bins=nothing,
        filt_bound=true,
        quantiles=[0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999]
    )

    if filt_bound
        filt = get_bound(snap)
        snap = snap[filt]
    end

    N_bound = length(snap)

    K = calc_K_tot(snap)

    if !isnothing(snap.Φs)
        W = calc_W_tot(snap)
        E = W + K
    else
        E = NaN
        W = NaN
    end

    L = calc_L_tot(snap)

    log_r_snap = log10.(calc_r(snap))

    bins, hist, err = histogram(log_r_snap, bins, weights=snap.masses)
    log_r_bins = bins
    mass_in_shell = hist 
    mass_in_shell_err  = err

    _, counts, _ = histogram(log_r_snap, bins)

    log_r = midpoints(log_r_bins)
    r_bins = 10 .^ log_r_bins
    r = 10 .^ log_r
    rel_err = mass_in_shell_err ./ mass_in_shell

    rho = calc_ρ_from_hist(r_bins, mass_in_shell)
    rho_err = rho .* rel_err # TODO: does not include binning error

    M_in = cumsum(mass_in_shell)
    M_in_err = 1 ./ sqrt.(cumsum(counts)) .* M_in

    skip = min(round(Int, default_n_per_bin(log_r_snap) / 1.5), 200)
    @info "circ vel bins $skip"
    r_circ, v_circ, n_circ = calc_v_circ(snap, skip=skip)
    v_circ_err = v_circ ./ sqrt.(n_circ)
    t_circ = 2π * r_circ ./ v_circ

    fit = fit_v_r_circ_max(r_circ, v_circ)
    v_circ_max = fit.v_circ_max
    r_circ_max = fit.r_circ_max

    r_quantile = 10 .^ quantile(log_r_snap, quantiles)

    return MassProfile3D(
        E=E,
        K=K,
        W=W,
        L=L,
        log_r=log_r,
        log_r_bins=log_r_bins,
        counts=counts,
        mass_in_shell=mass_in_shell,
        mass_in_shell_err=mass_in_shell_err,
        M_in=M_in,
        M_in_err=M_in_err,
        rho=rho,
        rho_err=rho_err,
        r_circ=r_circ,
        v_circ=v_circ,
        v_circ_err=v_circ_err,
        t_circ=t_circ,
        n_circ=n_circ,
        v_circ_max=v_circ_max,
        r_circ_max=r_circ_max,
        N_bound=N_bound,
        r_quantile=r_quantile,
        quantiles=quantiles,
        time=snap.time,
    )
end



"""
    calc_ρ_from_hist(bins, counts)
"""
function calc_ρ_from_hist(bins::AbstractVector{<:Real}, counts::AbstractVector{<:Real}) 
    if length(bins) != length(counts) + 1
        throw(DimensionMismatch("Bins must be one longer than counts, got sizes $(length(bins)), $(length(counts))"))
    end

    volumes = 4π/3 * diff(bins .^ 3)
    return counts ./ volumes
end



""" 
	calc_v_rad(snap)

returns the radial velocities relative to the snapshot centre in code units
"""
function calc_v_rad(snap)
    return calc_v_rad(snap.positions, snap.velocities, x_cen=snap.x_cen, v_cen=snap.v_cen)
end

"""
    calc_v_rad(positions, velocities; x_cen=zeros(3), v_cen=zeros(3))

Calculates the radial velocities relative to x_cen, v_cen.
"""
function calc_v_rad(positions::AbstractMatrix{<:Real}, velocities::AbstractMatrix{<:Real}; x_cen=zeros(3), v_cen=zeros(3))

    @assert_same_size positions velocities
    @assert_3vector positions

    x_vec = positions .- x_cen
    v_vec = velocities .- v_cen

    # normalize
    x_hat = x_vec ./ calc_r(x_vec)'

    # dot product
    v_rad = sum(x_hat .* v_vec, dims=1)

    # matrix -> vector
    v_rad = dropdims(v_rad, dims=1)
    
    return v_rad 
end



"""
    calc_v_circ(snap; filter_bound=true, skip=10)

Returns a list of the sorted radii and circular velocity from a snapshot for the given centre.
Removes unbound particles first if `filter_bound`, and skips every `skip` particles (i.e. each radius contains n times skip particles inclusive)
"""
function calc_v_circ(snap::Snapshot; filter_bound::Bool=true, skip::Integer=10)
    r = calc_r(snap)
    m = snap.masses

    if filter_bound
        filt = get_bound(snap)

        r = r[filt]
        m = m[filt]
    end

    idx = sortperm(r)
    m = snap.masses[idx]
    r = r[idx]
    M = cumsum(m)

    idx_skip = skip:skip:length(r)
    r = r[idx_skip]
    M = M[idx_skip]
    v_circ = calc_v_circ.(r, M)

    return r, v_circ, idx_skip
end



"""
    calc_v_circ_max(r, v_circ; q=80, p0=[6., 30.])

Fits the maximum circular velocity of a rotation curve assuming a NFW
profile. Returns the parameters of the fit and the range of radii used.
The fit is done on the top `q` quantile of the circular velocity, as th
tails often extend strangely
"""
function fit_v_r_circ_max(r::AbstractArray{<:Real}, v_circ::AbstractArray{<:Real}; q=0.9, p0=nothing)
    if p0 == nothing
        idx = argmax(v_circ)
        p0 = [r[idx], v_circ[idx]]
    end

    if length(r) != length(v_circ)
        throw(DimensionMismatch("r and v_circ must have the same length. Got $(length(r)) and $(length(v_circ))"))
    elseif length(r) < 2
        throw(ArgumentError("r and v_circ must have at least 2 elements"))
    end


    filt = v_circ .> quantile(v_circ, q)

    local fit, converged

    if sum(filt) < 2
        converged = false
        fit = nothing
    else
        try
            fit = curve_fit(_v_circ_max_model, r[filt], v_circ[filt], p0)
            converged = fit.converged
        catch ArgumentError
            converged = false
            fit = nothing
        end
    end


    if !converged
        @warn "Fit did not converge, using simple maximum."
        idx = argmax(v_circ)
        v_circ_max = v_circ[idx]
        r_circ_max = r[idx]
        r_min = NaN
        r_max = NaN
    else
        r_circ_max = fit.param[1]
        v_circ_max = fit.param[2]
        r_min = minimum(r[filt]) 
        r_max = maximum(r[filt])
    end

    return (; 
        r_circ_max=r_circ_max,
        v_circ_max=v_circ_max,
        fit=fit,
        converged=converged,
        r_min=r_min,
        r_max=r_max
       )
end


"""
    fit_v_r_max(snap; kwargs...)

Fits circular velocity of snapshot
"""
function fit_v_r_circ_max(snap; kwargs...)
    r, v_circ = calc_v_circ(snap)
    return fit_v_r_circ_max(r, v_circ; kwargs...)
end



@doc raw"""
NFW circular velocity model.

```math
v_circ(r) = v_circ_max/β * sqrt( (log(1 + x) - x/(1+x)) / x )
```

where x = r / r_circ_max * α_nfw
and α_nfw, β are the solution so v_circ_max_model(1) = 1.
"""
function _v_circ_max_model(r, param)
	r_circ_max, v_circ_max = param

	x = r ./ r_circ_max .* α_nfw
    β = 0.46499096281742197
    return @. v_circ_max / β  * sqrt(A_NFW(x) / x)
end



"""
    calc_M_in(out, radius; idxs=(1:10:length(output)))

Calculates the number of bound particles within a given radius for a set of snapshots.
"""
function calc_M_in(output::Output, radius::Real; idxs=(1:10:length(output)))
	N = length(idxs)
	M = Vector{Float64}(undef, N)
	
	for i in eachindex(idxs)
		snap = output[idxs[i]]
        M[i] = calc_M_in(snap, radius)
	end

	return M
end


"""
    calc_M_in(snap::Snapshot, radius::Real)

Calculates the number of bound particles within a given radius for a snapshot.
"""
function calc_M_in(snap::Snapshot, radius::Real)
    filt = get_bound(snap)
    rs = calc_r(snap)[filt]
    return sum(rs .<= radius)
end
