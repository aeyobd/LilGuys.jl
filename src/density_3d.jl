import TOML
import LinearAlgebra: ×, ⋅


# TODO: split MassProfile3D in half


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
        filt = bound_particles(snap)
        snap = snap[filt]
    end

    N_bound = length(snap)
    if N_bound == 0
        throw(ArgumentError("No bound particles in snapshot"))
    end

    K = kinetic_energy(snap)

    if !isnothing(snap.potential)
        W = potential_energy(snap)
        E = W + K
    else
        E = NaN
        W = NaN
    end

    L = angular_momentum(snap)

    log_r_snap = log10.(radii(snap))

    bins, hist, err = histogram(log_r_snap, bins, weights=snap.masses, errors=:weighted)
    log_r_bins = bins
    mass_in_shell = hist 
    mass_in_shell_err  = err

    _, counts, _ = histogram(log_r_snap, bins)

    log_r = midpoints(log_r_bins)
    r_bins = 10 .^ log_r_bins
    r = 10 .^ log_r
    rel_err = mass_in_shell_err ./ mass_in_shell

    rho = density_from_hist(r_bins, mass_in_shell)
    rho_err = rho .* rel_err # TODO: does not include binning error

    M_in = cumsum(mass_in_shell)
    M_in_err = 1 ./ sqrt.(cumsum(counts)) .* M_in

    skip = min(round(Int, default_n_per_bin(log_r_snap) / 1.5), 200)
    @info "circ vel bins $skip"
    r_circ, v_c, n_circ = v_circ(snap, skip=skip)
    v_circ_err = v_c ./ sqrt.(n_circ)
    t_circ = 2π * r_circ ./ v_c

    fit = fit_v_r_circ_max(r_circ, v_c)
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
        v_circ=v_c,
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
    density_from_hist(bins, counts)
"""
function density_from_hist(bins::AbstractVector{<:Real}, counts::AbstractVector{<:Real}) 
    if length(bins) != length(counts) + 1
        throw(DimensionMismatch("Bins must be one longer than counts, got sizes $(length(bins)), $(length(counts))"))
    end

    volumes = 4π/3 * diff(bins .^ 3)
    return counts ./ volumes
end



""" 
	radial_velocities(snap)

returns the radial velocities relative to the snapshot centre in code units
"""
function radial_velocities(snap)
    return radial_velocities(snap.positions, snap.velocities, x_cen=snap.x_cen, v_cen=snap.v_cen)
end

"""
    radial_velocities(positions, velocities; x_cen=zeros(3), v_cen=zeros(3))

Calculates the radial velocities relative to x_cen, v_cen.
"""
function radial_velocities(positions::AbstractMatrix{<:Real}, velocities::AbstractMatrix{<:Real}; x_cen=zeros(3), v_cen=zeros(3))

    @assert_same_size positions velocities
    @assert_3vector positions

    x_vec = positions .- x_cen
    v_vec = velocities .- v_cen

    # normalize
    x_hat = x_vec ./ radii(x_vec)'

    # dot product
    v_r = sum(x_hat .* v_vec, dims=1)

    # matrix -> vector
    v_r = dropdims(v_r, dims=1)
    
    return v_r 
end



"""
    v_circ(snap; filter_bound=:recursive_1d, skip=10)

Returns a list of the sorted radii and circular velocity from a snapshot for the given centre.
Skips every `skip` particles (i.e. each radius contains n times skip particles inclusive)
filter_bound may be
- :recursive_1d: recursively removes unbound particles, recomputing the potential each time assuming spherical symmetry
- :recursive: like recursive_1d but computes the full potential
- :simple: only removes unbound particles once
- :false: does not filter particles
"""
function v_circ(snap::Snapshot; filter_bound=:recursive_1d, skip::Integer=10)
    r = radii(snap)
    m = snap.masses

    if filter_bound != :false
        if filter_bound == :simple
            filt = bound_particles(snap)
        elseif filter_bound == :recursive_1d
            filt = bound_particles_recursive_1D(snap)
        else
            throw(ArgumentError("Unknown filter_bound: $filter_bound"))
        end

        r = r[filt]
        m = m[filt]
    end

    if sum(r) == 0
        @warn "No bound particles in snapshot"
        return zeros(0), zeros(0), zeros(0)
    end

    idx = sortperm(r)
    m = snap.masses[idx]
    r = r[idx]
    M = cumsum(m)

    idx_skip = skip:skip:length(r)
    r = r[idx_skip]
    M = M[idx_skip]
    v_c = v_circ.(r, M)

    return r, v_c, idx_skip
end



"""
    fit_v_r_circ_max(r, v_circ; q=80, p0=[6., 30.])

Fits the maximum circular velocity of a rotation curve assuming a NFW
profile. Returns the parameters of the fit and the range of radii used.
The fit is done on the top `q` quantile of the circular velocity, as th
tails often extend strangely
"""
function fit_v_r_circ_max(r::AbstractArray{<:Real}, v_circ::AbstractArray{<:Real}; q=0.9, p0=nothing)
    if length(r) != length(v_circ)
        throw(DimensionMismatch("r and v_circ must have the same length. Got $(length(r)) and $(length(v_circ))"))
    elseif length(r) < 2
        @warn "r and v_circ must have at least 2 elements"
        return (; r_circ_max=NaN, v_circ_max=NaN, fit=nothing, converged=false, r_min=NaN, r_max=NaN)
    end

    if p0 == nothing
        idx = argmax(v_circ)
        p0 = [r[idx], v_circ[idx]]
    end


    filt = v_circ .> quantile(v_circ, q)

    local fit, converged

    if sum(filt) < 2
        converged = false
        fit = nothing
    else
        try
            fit = Interface.LsqFit.curve_fit(_v_circ_max_model, r[filt], v_circ[filt], p0)
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
    r, v_c = v_circ(snap)
    return fit_v_r_circ_max(r, v_c; kwargs...)
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
    M_in(snap::Snapshot, radius::Real)

Calculates the number of bound particles within a given radius for a snapshot.
"""
function M_in(snap::Snapshot, radius::Real)
    filt = bound_particles(snap)
    rs = radii(snap)[filt]
    return sum(rs .<= radius)
end


"""
    β_prof(snap; r_bins)

Computes the velocity dispersion and anisotropy profiles
for the snapshot assuming spherical symmetry.
Returns a tuple or radii bins, velocity dispersion, and anisotropy
"""
function β_prof(snap; r_bins)
    vel = snap.velocities .- snap.v_cen
    pos = snap.positions .- snap.x_cen
    N = size(vel, 2)

    rs = radii(snap)
    r_hat = pos ./ rs'
    
    ϕ_hat = hcat([r_hat[:, i] × [0,0,1] for i in 1:N]...)
    ϕ_hat ./= radii(ϕ_hat)'
    θ_hat = hcat([r_hat[:, i] × ϕ_hat[:, i] for i in 1:N]...)
    θ_hat ./= radii(θ_hat)'

    v_r = [vel[:, i] ⋅ r_hat[:, i] for i in 1:N]
    v_ϕ = [vel[:, i] ⋅ ϕ_hat[:, i] for i in 1:N]
    v_θ = [vel[:, i] ⋅ θ_hat[:, i] for i in 1:N]


    Nb = length(r_bins) - 1
    βs = zeros(Nb)
    σ3d = zeros(Nb)

    for i in 1:Nb
        filt = rs .> r_bins[i]
        filt .&= rs .<= r_bins[i+1]

        σ2_r = std(v_r[filt])^2
        σ2_t = std(v_θ[filt])^2 + std(v_ϕ[filt])^2

        βs[i] = 1 - σ2_t / 2σ2_r
        σ3d[i] = σ2_t^2 + σ2_r^2
    end

    return σ3d, βs
end


"""
    cartesian_to_cylindrical(x, y, z)

Converts cartesian coordinates to cylindrical coordinates
"""
function cartesian_to_cylindrical(x::T, y::T, z::T) where T <: Union{Real, AbstractVector}
    R = sqrt.(x.^2 + y.^2)
    z = z
    ϕ = atan.(y, x)

    return r, z, ϕ
end
