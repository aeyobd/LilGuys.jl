import DensityEstimators: histogram
import StatsBase: percentile
import TOML


"""
An simulated 3D density profile. 
All properties are in code units.
"""
@kwdef struct ObsProfile3D
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
    r_circ_max::F
    N_bound::Int

    log_r::Vector{F}
    log_r_bins::Vector{F}
    counts::Vector{Int}

    mass_in_shell::Vector{F}
    mass_in_shell_err::Vector{F}

    M_in::Vector{F}
    M_in_err::Vector{F}

    rho::Vector{F}
    rho_err::Vector{F}

    v_circ::Vector{F}
    v_circ_err::Vector{F}

    t_circ::Vector{F}
end


"""
A collection of 3D density profiles.
"""
struct Profiles3D <: AbstractVector{ObsProfile3D}
    snapshot_index::Vector{Int}
    profiles::Vector{ObsProfile3D}
end

function Base.size(profs::Profiles3D)
    N =  length(profs.profiles)

    return (N,)
end



function Base.getindex(profs::Profiles3D, i::Int)
    return profs.profiles[i]
end


function Profiles3D(filename::String)
    h5open(filename, "r") do f
        snapshot_index = read(f, "snapshot_index")

        profiles = ObsProfile3D[]

        for i in eachindex(snapshot_index)
            df = Dict()

            for key in fieldnames(ObsProfile3D)
                ndims = length(size(read(f, string(key))))
                if ndims == 1
                    data = f[string(key)][i]
                else
                    data = f[string(key)][:, i]
                end
                df[key] = data
            end

            prof = ObsProfile3D(;df...)

            push!(profiles, prof)
        end
        return Profiles3D(snapshot_index, profiles)
    end
end


function save(filename::String, profs::Profiles3D)
    h5open(filename, "w") do f
        write(f, "snapshot_index", profs.snapshot_index)

        for key in fieldnames(ObsProfile3D)
            data = getfield.(profs.profiles, key)
            if eltype(data) <: Real
                write(f, string(key), data)
            elseif eltype(data) <: AbstractVector
                write(f, string(key), hcat(data...))
            end
        end
    end
end



function Base.print(io::IO, prof::ObsProfile3D)
    TOML.print(io, struct_to_dict(prof))
end


function ObsProfile3D(filename::String)
    t = dict_to_tuple(TOML.parsefile(filename))
    return ObsProfile3D(;t...)
end

"""
    calc_profile(snap)
"""
function calc_profile(snap::Snapshot;
        bins=100,
        filt_bound=true,
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

    h = histogram(log10.(calc_r(snap)), bins, weights=snap.masses, normalization=:none)
    log_r_bins = h.bins
    log_r = midpoints(log_r_bins)
    r = 10 .^ log_r


    counts = histogram(log10.(calc_r(snap)), bins, normalization=:none).values

    mass_in_shell = h.values
    mass_in_shell_err = h.err

    V = 4π/3 * diff((10 .^ log_r_bins) .^ 3)
    rho = mass_in_shell ./ V
    rho_err = mass_in_shell_err ./ V


    M_in = cumsum(mass_in_shell)
    M_in_err = cumsum(counts) .^ -0.5 .* M_in

    # circular velocity is mass inclusive
    r_right = 10 .^ log_r_bins[2:end]
    v_circ = calc_v_circ.(r_right, M_in)
    v_circ_err = zeros(length(v_circ))
    t_circ = r_right ./ v_circ
    t_circ_err = zeros(length(t_circ))

    fit = fit_v_r_circ_max(r, v_circ)
    v_circ_max = fit.v_circ_max
    r_circ_max = fit.r_circ_max

    return ObsProfile3D(
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
        v_circ=v_circ,
        v_circ_err=v_circ_err,
        t_circ=t_circ,
        v_circ_max=v_circ_max,
        r_circ_max=r_circ_max,
        N_bound=N_bound
    )
end



"""

sorts a snapshot by radius from 0
"""
function sort_by_r(snap::Snapshot)
    return snap[sortperm(calc_r(snap))]
end



"""
    calc_m_hist(r, r_bins[, masses])

Calculates the density profile given a set of particles located at `r` with masses `masses` by binning into `r_bins`.
"""
function calc_ρ_hist(r::AbstractVector{T}, bins::AbstractVector; weights=nothing) where T <: Real
    if weights == nothing
        weights = ones(length(r))
    end

    counts = histogram(r, bins, weights=weights, normalization=:none).values

    Vs = 4π/3 * diff(bins .^ 3)
    return bins, counts ./ Vs
end



function calc_ρ_hist(r::AbstractVector{T}, bins::Int; weights=nothing, equal_width=false) where T <: Real
    if equal_width
        x1 = minimum(r)
        x2 = maximum(r)
        x = LinRange(x1, x2, bins)
        bins = 10 .^ x
    else
        bins = percentile(r, LinRange(0, 100, bins+1))
    end
    return calc_ρ_hist(r, bins; weights=weights)
end


function calc_ρ_hist(r::AbstractVector{T}; weights=nothing) where T <: Real
    r_bins = round(Int64, 0.1 * sqrt(length(r)))
    return calc_ρ_hist(r, r_bins, weights=weights)
end


function calc_ρ_hist(snap::Snapshot, bins; weights=snap.masses, x_cen=snap.x_cen)
    r = calc_r(snap, x_cen)
    return calc_ρ_hist(r, bins; weights=weights)
end



""" 
	calc_v_rad(snap)

returns the radial velocities relative to the snapshot centre in code units
"""
function calc_v_rad(snap)
	x_vec = snap.positions .- snap.x_cen
	v_vec = snap.velocities .- snap.v_cen

	# normalize
	x_hat = x_vec ./ calc_r(x_vec)'

	# dot product
	v_rad = sum(x_hat .* v_vec, dims=1)

	# matrix -> vector
	v_rad = dropdims(v_rad, dims=1)
	
	return v_rad 
end



"""
    calc_v_circ(snap; x_cen, filter_bound)

Returns a list of the sorted radii and circular velocity from a snapshot for the given centre.
"""
function calc_v_circ(snap::Snapshot; x_cen=snap.x_cen, filter_bound=true)
    r = calc_r(snap.positions .- x_cen)
    m = snap.masses[sortperm(r)]
    r = sort(r)
    M = cumsum(m)

    if filter_bound
        filt = get_bound(snap)

        r = r[filt]
        M = M[filt]
    end
    return r, calc_v_circ.(r, M)
end



"""
    calc_v_circ_max(r, v_circ; percen=80, p0=[6., 30.])

Fits the maximum circular velocity of a rotation curve assuming a NFW
profile. Returns the parameters of the fit and the range of radii used.
"""
function fit_v_r_circ_max(r, v_circ; percen=80, p0=[6., 30.])
    filt = v_circ .> percentile(v_circ, percen)

    local fit, converged
    try
        fit = curve_fit(v_circ_max_model, r[filt], v_circ[filt], p0)
        converged = fit.converged
    catch ArgumentError
        converged = false
        fit = nothing
    end


    if converged == false
        @warn "Fit did not converge, using simple maximum."
        idx = argmax(v_circ)
        v_circ_max = v_circ[idx]
        r_circ_max = r[idx]
    else
        r_circ_max = fit.param[1]
        v_circ_max = fit.param[2]
    end

    return (; 
        r_circ_max=r_circ_max,
        v_circ_max=v_circ_max,
        r_min = minimum(r[filt]), 
        r_max = maximum(r[filt]),
        fit=fit,
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
function v_circ_max_model(r, param)
	Rmx,Vmx = param
	x = r ./ Rmx .* α_nfw
	inner = @. (nm.log(1+x) - x/(1+x)) / x
	return @. Vmx / 0.46499096281742197 * nm.sqrt(inner)
end



"""
    calc_M_h(output, radius; idxs=(1:10:length(output)))

Calculates the number of particles within a given radius for a set of snapshots.
"""
function get_M_h(output::Output, radius; idxs=(1:10:length(output)))
	N = length(idxs)
	M = Vector{Float64}(undef, N)
	
	for i in eachindex(idxs)
		snap = output[idxs[i]]
		ϵ = calc_ϵ(snap)
		filt = ϵ .> 0
		
		rs = calc_r(snap[filt])
		filt2 = rs .< radius

		M[i] = sum(filt2)
	end

	return M
end


"""
    to_gaia(snap; params...)

Converts a snapshot to a Gaia-like DataFrame.


# Arguments
- `p_min::Float=1e-20`: Minimum probability (relative to maximum) for stars to be included
- `filt_bound::Bool=false`: If true, only bound particles are included
- `add_centre::Bool=true`: If true, adds the centre of the snapshot as the first observation (with zero weight)
## Arguments for `to_sky`
- `SkyFrame::CoordinateFrame=ICRS`: The frame to convert the observations to
- `invert_velocity::Bool=false`: If true, the velocity is inverted. Useful for when time is backwards (e.g. orbital analysis)

"""
function to_gaia(snap::Snapshot; 
        p_min=1e-20, filt_bound=false, add_centre=true, 
        set_to_distance=nothing,
        kwargs...)

    add_weights = !(snap.weights isa ConstVector)

    if !add_weights
        filt = trues(length(snap))
    else
        filt = snap.weights .>= p_min * maximum(snap.weights)
    end

    @info "excluded stellar mass: $(sum(snap.weights[.!filt]) / sum(snap.weights))"

    if filt_bound
        filt_b = get_bound(snap)
        @info "unbound stellar mass: $(sum(snap.weights[.!filt_b]) / sum(snap.weights))"
        filt .&= get_bound(snap)
    end

    @info "Converting $(sum(filt)) particles to Gaia-like DataFrame"

    if set_to_distance != nothing
        @assert !(snap.x_cen[1] == snap.x_cen[2] == snap.x_cen[3] == 0) "Snapshot is not centred"

        f = GalactocentricFrame()
        sun_vec = [f.d, f.z_sun, 0]
        r_vec = snap.x_cen .- sun_vec
        vec_final = r_vec ./ norm(r_vec) * set_to_distance .+ sun_vec
        vec_shift = vec_final .- snap.x_cen
        snap.x_cen .+= vec_shift
        snap.positions .+= vec_shift

        offset = calc_r(snap.x_cen, sun_vec) - set_to_distance
        println(vec_final, r_vec, sun_vec)
        @assert abs(offset) < 1e-3 "Offset is too large: $offset"
    end

    observations = to_sky(snap[filt]; kwargs...)
    obs_cen = phase_to_sky(snap.x_cen, snap.v_cen; kwargs...)

    df = to_frame(observations)
    df[!, :index] = snap.index[filt]


    if add_weights
        df[!, :weights] = snap.weights[filt]
    end

    if add_centre
        df_cen = to_frame([obs_cen])
        df_cen[!, :index] = [0]
        if add_weights
            df_cen[!, :weights] = [0]
        end
        df = vcat(df_cen, df)
    end

    df[!, :xi], df[!, :eta] = to_tangent(df.ra, df.dec, obs_cen.ra, obs_cen.dec)

    df[!, :r_ell] = 60*calc_r_ell(df.xi, df.eta, 0, 0)

    return df
end


"""
    to_sky(snap::Snapshot, invert_velocity=false, verbose=false, SkyFrame=ICRS)

Returns a list of observations based on snapshot particles. 

# Arguments
- `verbose::Bool=false`: If true, prints the progress
- `SkyFrame::CoordinateFrame=ICRS`: The frame to convert the observations to
- `invert_velocity::Bool=false`: If true, the velocity is inverted. Useful for when time is backwards (e.g. orbital analysis)
- `kwargs...`: Additional arguments for `phase_to_sky`
"""
function to_sky(snap::Snapshot; 
        verbose::Bool=false,
        SkyFrame=ICRS,
        kwargs...
    )
    observations = Vector{SkyFrame}(undef, length(snap))

    for i in 1:length(snap)
        if verbose
            print("converting $(i)/($(length(snap))\r")
        end

        pos = snap.positions[:, i]
        vel = snap.velocities[:, i]
        obs = phase_to_sky(pos, vel; SkyFrame, kwargs...)
        observations[i] = obs
    end


    return observations
end



"""
    phase_to_sky(pos, vel; invert_velocity=false, SkyFrame=ICRS)
Converts a phase space position and velocity (i.e. simulation point) to a sky observation.

# Arguments
- `invert_velocity::Bool=false`: If true, the velocity is inverted. Useful for when time is backward(e.g. orbital analysis)
- `SkyFrame::CoordinateFrame=ICRS`: The frame to convert the observations to
"""
function phase_to_sky(pos::Vector{F}, vel::Vector{F}; invert_velocity=false, SkyFrame=ICRS)
    if invert_velocity
        vel *=-1
    end

    gc = Galactocentric(pos*R2KPC, vel*V2KMS)
    obs = transform(SkyFrame, gc)
    return obs
end


"""
    to_frame(obs::AbstractVector{CoordinateFrame})

Converts a list of observations to a DataFrame with the same column names
"""
function to_frame(obs::AbstractVector{T}) where T<:CoordinateFrame
    cols = propertynames(obs[1])
    df = DataFrame()
    for col in cols
        val = getproperty.(obs, col)
        
        if val[1] isa Union{Real, String}
            df[!, Symbol(col)] = val
        else
            @info "Skipping column $col"
        end 
    end

    return df
end



