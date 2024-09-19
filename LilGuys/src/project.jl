

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
    println(names(df))
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


