

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
- `kwargs...`: Additional arguments for `phase_to_sky`
- `set_to_distance::Union{Nothing, Real}=nothing`: If set, the snapshot is shifted to this distance from the sun
"""
function to_gaia(snap::Snapshot; 
        p_min=1e-20, filt_bound=false, add_centre=true, 
        set_to_distance=nothing,
        filt_wrong_hemisphere=false,
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
        filt .&= filt_b
    end

    if set_to_distance != nothing
        snap = shift_snapshot_to_distance(snap, set_to_distance)
    end

    @info "transforming coordinate frames"
    observations = to_sky(snap[filt]; kwargs...)
    obs_cen = phase_to_sky(snap.x_cen, snap.v_cen; kwargs...)

    @info "converting $(sum(filt)) particles to Gaia-like DataFrame"
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
    df[!, :r_ell] = 60*calc_R_ell(df.xi, df.eta, 0, 0)


    if filt_wrong_hemisphere
        filt_nan = isnan.(df.xi) .| isnan.(df.eta)
        @info "excluded $(sum(filt_nan)) observations on wrong hemisphere"
        df = df[.!filt_nan, :]
    end

    return df
end

function shift_snapshot_to_distance(snap::Snapshot, set_to_distance::Real)
    snap = deepcopy(snap)

    if snap.x_cen[1] == snap.x_cen[2] == snap.x_cen[3] == 0
        @warn "Snapshot centre is at the origin"
    end

    f = GalactocentricFrame()
    ρ = asin(f.z_sun/f.d)
    sun_vec = [-f.d*cos(ρ), 0, f.d*sin(ρ)]

    r_vec = snap.x_cen .- sun_vec
    vec_final = r_vec ./ norm(r_vec) * set_to_distance .+ sun_vec
    vec_shift = vec_final .- snap.x_cen

    @info "shifting snapshot $(vec_shift) to $(vec_final) "
    snap.x_cen .+= vec_shift
    snap.positions .+= vec_shift

    offset = radii(snap.x_cen, sun_vec) - set_to_distance

    @assert abs(offset) < 1e-3 "Offset is too large: $offset"

    return snap
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
            @info "converting $(i)/($(length(snap))\r"
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



# Utility functions


"""
    calc_R_ell_sky(ra, dec, a, b, PA; weights=nothing, centre="mean", units="arcmin")

Given a set of sky coordinates (ra, dec), computes the elliptical radius of each point with respect to the centre of the ellipse defined by the parameters (a, b, PA).
"""
function calc_R_ell_sky(ra, dec, a, b, PA; weights=nothing,
        centre="mean",
        units="arcmin"
    )
    ra0, dec0 = calc_centre2D(ra, dec, centre, weights)

    x, y = to_tangent(ra, dec, ra0, dec0)

    R_ell = calc_R_ell(x, y, a, b, PA)


    if units == "arcmin"
        R_ell = 60R_ell
    elseif units == "arcsec"
        R_ell = 3600R_ell
    elseif units == "deg"
        R_ell = 1R_ell
    else
        error("units not implemented: $units")
    end

    return R_ell
end


function calc_R_ell_sky(ra, dec, ell, PA; kwargs...)
    aspect = ellipticity_to_aspect(ell)
    b = sqrt(aspect)
    a = 1/b
    return calc_R_ell_sky(ra, dec, a, b, PA; kwargs...)
end


function calc_R_ell_sky(ra, dec; kwargs...)
    return calc_R_ell_sky(ra, dec, 0, 0; kwargs...)
end


"""
    calc_R_ell(x, y, a, [b, ]PA)

computes the elliptical radius of a point (x, y) with respect to the center (0, 0) and the ellipse parameters (a, b, PA).
If using sky coordinates, x and y should be tangent coordinates.

Note that the position angle is the astronomy definition, i.e. measured from the North to the East (clockwise) in xi / eta.
"""
function calc_R_ell(x, y, args...)
    x_p, y_p = shear_points_to_ellipse(x, y, args...)

    R_sq = @. (x_p)^2 + (y_p)^2
    return sqrt.(R_sq)
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
    to_orbit_coords(ra, dec, ra0, dec0, PA; unit=:degree) 

Given the position angle of an orbit vector at (ra0, dec0) calculate a rotated
sky frame centred on RA, DEC and with the x-axis aligned with the orbit.
The orbital position angle can be found in the orbit properties file from
analyze_orbit.jl notebook. Returns coordinates in the units specified by `unit`
"""
function to_orbit_coords(ra, dec, ra0::Real, dec0::Real, PA::Real; unit=:degree)
	# want to rotate to dec, ra of centre, then rotate 
	α = deg2rad(ra0)
	δ = deg2rad(dec0)
	ϖ = deg2rad(90 - PA)
	Rmat = Rx_mat(-ϖ) * Ry_mat(δ) * Rz_mat(-α) 

	coords = unit_vector(ra, dec)
	coords =  Rmat * coords
	ra, dec, _ = cartesian_to_sky(coords[1, :], coords[2, :], coords[3, :])

	ra .-= 360 * (ra .> 180)

    if unit == :arcmin
        return 60ra, 60dec
    elseif unit == :degree
        return ra, dec
    else
        @error "unknown unit $unit"
    end
end



function to_orbit_coords(ra::Real, dec::Real, ra0::Real, dec0::Real, PA::Real)
    ra, dec =  to_orbit_coords([ra], [dec], ra0, dec0, PA)
    return ra[1], dec[1]
end

