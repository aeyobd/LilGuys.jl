import TOML
import LinearAlgebra: ×, ⋅


"""
    MassProfile(<keyword arguments>)

A struct representing a 3-dimensional 
mass profile and cummulative / summative properties
(of a snapshot).

All properties are in code units.
"""
@kwdef struct MassProfile
    "Radius for v_circ"
    radii::Vector{F}

    "The number of particles within each radius"
    counts::Vector{F}

    "Enclosed mass"
    M_in::Vector{Measurement{F}}

    "Snapshot time"
    time::F = NaN

    "Additional annotations"
    annotations::Dict{String, Any} = Dict{String, Any}()
end

radii(prof::MassProfile) = prof.radii
masses(prof::MassProfile) = prof.mass
circular_velocity(prof::MassProfile) = @. √(G*prof.M_in / prof.radii)


"""
    MassScalars

A struct representing properties of an entire gravitation system
(like energies, angular momentum, and maximum circular velocity).
"""
@kwdef struct MassScalars
    # scalars
    """ Total Energy """
    E::F

    "Total Kinetic energy"
    K::F

    "Total Potential energy"
    W::F

    "Total Angular momentum vector wrt centre"
    L::Vector{F}

    "Maximum circular velocity"
    v_circ_max::F

    "Radius of maximum circular velocity"
    r_circ_max::F

    "Number of bound particles"
    N_bound::Int

    "bound mass"
    bound_mass::F

    "Snapshot time"
    time::F = NaN
end





function Base.print(io::IO, prof::MassProfile)
    TOML.print(io, struct_to_dict(prof))
end


function MassProfile(filename::String)
    t = dict_to_tuple(TOML.parsefile(filename))
    return MassProfile(;t...)
end



"""
    MassProfile(snap; bins=100, filt_bound=true)

Given a snapshot, computes the density, circular velocity, and energy; 
returning an MassProfile object.

`bins` may be an integer, array, or function and is passed to `histogram`
"""
function MassProfile(snap::Snapshot;
        bins=nothing,
        filt_bound=:recursive_1D,
        quantiles=[0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999],
        skip = nothing,
    )

    r = radii(snap)
    m = snap.masses

    if filt_bound != :false
        filt = bound_particles(snap, method=filt_bound)
        r = r[filt]
        m = m[filt]
    end

    if length(r) == 0
        throw(DomainError("No bound particles in snapshot"))
    end

    idx = sortperm(r)
    m = snap.masses[idx]
    r = r[idx]
    M = cumsum(m)


    if isnothing(skip)
        skip = min(round(Int, default_n_per_bin(r) / 1.5), 200)
        @info "MassProfile binsize: $skip"
    end

    idx_skip = skip:skip:length(r)
    counts = sqrt.(idx_skip)
    r = r[idx_skip]
    M = M[idx_skip]
    M_err = @. 1/sqrt(counts) * M


    return MassProfile(
        radii = r,
        M_in = Measurement.(M, M_err),
        counts = counts,
        time = snap.time,
    )
end


function MassScalars(snap, prof::MassProfile)
    K = kinetic_energy(snap)

    if !isnothing(snap.potential)
        W = potential_energy(snap)
        E = W + K
    else
        E = NaN
        W = NaN
    end

    L = angular_momentum(snap)

    fit = fit_v_r_circ_max(radii(prof), middle.(circular_velocity(prof)))
    v_circ_max = fit.v_circ_max
    r_circ_max = fit.r_circ_max


    filt = bound_particles(snap)
    N_bound = sum(filt)
    bound_mass = sum(snap.masses[filt])

    return MassScalars(
        # scalars
        v_circ_max=v_circ_max,
        r_circ_max=r_circ_max,
        E=E,
        K=K,
        W=W,
        L=L,
        N_bound=N_bound,
        bound_mass = bound_mass
   )
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

