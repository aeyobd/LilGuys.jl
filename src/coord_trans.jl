import Base: @kwdef


"""Convert a coordinate frame to ICRS"""
function to_icrs(coord::CoordinateFrame)::ICRS
    throw(MethodError(to_icrs, (coord, )))
end

"""Convert a coordinate frame from ICRS"""
function from_icrs(frame::Type{<:CoordinateFrame}, coord::ICRS)::CoordinateFrame
    throw(MethodError(from_icrs, (frame, coord)))
end


"""
    transform(T, obs)

Transform a coordinate of type `obs` to an object of type `T`.

Note, all kwargs passed to the `from_icrs` implementation for the given types.
"""
function transform(::Type{T}, obs::CoordinateFrame; kwargs...) where {T<:CoordinateFrame}
    return from_icrs(T, to_icrs(obs); kwargs...)
end


# identity
function transform(::Type{T}, obs::T) where {T}
    return obs
end


function to_icrs(icrs::ICRS)
    return icrs
end

function from_icrs(T::Type{<:ICRS}, icrs::ICRS)
    return icrs
end


# automatic cartesian conversion
function transform(::Type{<:Cartesian{T, <:Real}}, obs::T) where {T<:CoordinateFrame}
    return Cartesian(obs)
end

# automatic cartesian conversion
function transform(::Type{T}, obs::Cartesian{T}) where {T<:CoordinateFrame}
    return to_sky(obs)
end


# to make ICRS behave with multiple dispatch
function transform(::Type{<:ICRS}, obs::Cartesian{ICRS})
    return to_sky(obs)
end

function transform(::Type{<:Cartesian{ICRS}}, obs::ICRS)
    return Cartesian(obs)
end


function transform(::Type{<:ICRS}, obs::ICRS)
    return obs
end

function transform(::Type{<:Cartesian{ICRS}}, obs::Cartesian{ICRS})
    return obs
end

function to_icrs(cart::Cartesian{ICRS, <:Real})
    return to_sky(cart)
end



# ========= Cartesian =========
#

function to_icrs(coord::Cartesian{<:T, <:Real}; kwargs...) where {T<:CoordinateFrame}
    return to_icrs(to_sky(coord); kwargs...)
end


function from_icrs(::Type{<:Cartesian{T, <:Real}}, icrs::ICRS; kwargs...) where {T<:CoordinateFrame}
    return Cartesian(from_icrs(T, icrs; kwargs...))
end



function Cartesian{T, F}(sc::T) where {T<:AbstractSkyCoord, F<:Real}
    position = sky_to_cartesian(sc.ra, sc.dec, sc.distance)
    velocity = _observation_to_cartesian_velocity(sc)
    return Cartesian{T, F}(position, velocity, sc)
end


function Cartesian(sc::T) where {T<:AbstractSkyCoord}
    F = typeof(sc.ra)
    return Cartesian{T, F}(sc)
end

function Cartesian{<:T}(sc::T) where {T<:AbstractSkyCoord}
    F = typeof(sc.ra)
    return Cartesian{T, F}(sc)
end

function to_sky(coord::Cartesian{T, <:Real}; kwargs...) where {T<:CoordinateFrame}
    ra, dec, distance = cartesian_to_sky(coord.x, coord.y, coord.z)
    pmra, pmdec, radial_velocity = _cartesian_to_observation_velocity(coord)

    return T(ra=ra, dec=dec, distance=distance, pmra=pmra, pmdec=pmdec, 
             radial_velocity=radial_velocity; kwargs...)
end


function _cartesian_to_observation_position(cart::Cartesian)
    x, y, z = cart.x, cart.y, cart.z

    return cartesian_to_sky(x, y, z)
end


function _cartesian_to_observation_velocity(cart::Cartesian)
    x, y, z = cart.x, cart.y, cart.z
    v_x, v_y, v_z = cart.v_x, cart.v_y, cart.v_z

    R = sqrt(x^2 + y^2)
    r = sqrt(x^2 + y^2 + z^2)

    v_r = (x*v_x + y*v_y + z*v_z) / r
    v_α = (x*v_y - y*v_x) / (r*R)
    v_δ = v_z*R/r^2 - (x*v_x + y*v_y) * z/(r^2 * R)

    return [v_α/kms_per_kpc_mas_per_yr, v_δ/kms_per_kpc_mas_per_yr, v_r]
end




function _observation_to_cartesian_velocity(obs::AbstractSkyCoord)
    rv = obs.radial_velocity
    α = obs.ra
    δ = obs.dec
    v_α_cosδ = obs.pmra * kms_per_kpc_mas_per_yr * obs.distance
    v_δ = obs.pmdec  * kms_per_kpc_mas_per_yr * obs.distance
    d = obs.distance

    vx = (rv * cosd(α) * cosd(δ) 
          - sind(α) * v_α_cosδ 
          - cosd(α) * sind(δ) * v_δ
         )

    vy = (rv * sind(α) * cosd(δ) 
          + cosd(α) * v_α_cosδ 
          - sind(α) * sind(δ) * v_δ
         )
    vz = (rv * sind(δ) 
          + cosd(δ) * v_δ)

    return [vx, vy, vz]
end


# ========= Galactocentric transformations =========


function from_icrs(::Type{<:Galactocentric}, icrs::ICRS{F}; frame=default_gc_frame) where {F<:Real}
    cart = Cartesian(icrs)
    return _icrs_cart_to_galcen(cart; frame=frame)
end


function transform(::Type{<:Galactocentric}, cart::Cartesian{<:ICRS, <:Real}; frame=default_gc_frame)
    return _icrs_cart_to_galcen(cart; frame=frame)
end

function _icrs_cart_to_galcen(cart::Cartesian{<:ICRS, F}; frame=default_gc_frame) where {F}
    x_gc = _heliocen_to_galcen_position(cart, frame=frame)
    v_gc = _heliocen_to_galcen_velocity(cart, frame=frame)

    return Galactocentric{F}(x_gc..., v_gc...; frame=frame)
end

function to_icrs(galcen::Galactocentric)::ICRS
    cart = _galcen_to_icrs_cart(galcen)
    return to_sky(cart)
end


function transform(::Type{<:Cartesian{ICRS, <:Real}}, galcen::Galactocentric)
    return _galcen_to_icrs_cart(galcen)
end

function _galcen_to_icrs_cart(galcen::Galactocentric)
    x_icrs = _galcen_to_heliocen_position(galcen)
    v_icrs = _galcen_to_heliocen_velocity(galcen)

    return Cartesian{ICRS}(x_icrs, v_icrs)
end



# ========= GSR transformations =========


function from_icrs(::Type{<:GSR}, icrs::ICRS{F}; kwargs...) where {F<:Real}
    distance1 = isnan(icrs.distance) ? 1 : icrs.distance
    pmra1 = isnan(icrs.pmra) ? 0 : icrs.pmra
    pmdec1 = isnan(icrs.pmdec) ? 0 : icrs.pmdec

    icrs2 = ICRS(ra=icrs.ra, dec=icrs.dec, distance=distance1, pmra=pmra1,
                 pmdec=pmdec1, radial_velocity=icrs.radial_velocity)

    gc = from_icrs(Galactocentric, icrs2; kwargs...)
    gsr = transform(GSR, _gc_to_gsr_cart(gc))

    if isnan(icrs.distance) 
        distance = NaN
    else
        distance = gsr.distance
    end
    if isnan(icrs.distance) || isnan(icrs.pmra)
        pmra = NaN
    else
        pmra = gsr.pmra
    end

    if isnan(icrs.distance) || isnan(icrs.pmdec)
        pmdec = NaN
    else
        pmdec = gsr.pmdec
    end

    return GSR(ra=gsr.ra, dec=gsr.dec, distance=distance, pmra=pmra, pmdec=pmdec, radial_velocity=gsr.radial_velocity)
end


function from_icrs(::Type{<:Cartesian{GSR}}, icrs::ICRS; kwargs...) 
    gc = from_icrs(Galactocentric, icrs)
    return _gc_to_gsr_cart(gc)
end

function _gc_to_gsr_cart(galcen::Galactocentric)
    x_gsr = _galcen_to_heliocen_position(galcen)
    v_gsr = _galcen_to_gsr_velocity(galcen)

    return Cartesian{GSR}(x_gsr, v_gsr)
end



function transform(::Type{<:Galactocentric}, cart::Cartesian{<:GSR, <:Real}; frame=default_gc_frame)
    x_gc = _heliocen_to_galcen_position(cart; frame)
    v_gc = _gsr_to_galcen_velocity(cart; frame)

    return Galactocentric(x_gc, v_gc, frame=frame)
end


function to_icrs(gsr::GSR)::ICRS
    distance1 = isnan(gsr.distance) ? 1 : gsr.distance
    pmra1 = isnan(gsr.pmra) ? 0 : gsr.pmra
    pmdec1 = isnan(gsr.pmdec) ? 0 : gsr.pmdec
    gsr2 = GSR(ra=gsr.ra, dec=gsr.dec, distance=distance1, pmra=pmra1, pmdec=pmdec1, radial_velocity=gsr.radial_velocity)

    cart = Cartesian(gsr2)
    gc = transform(Galactocentric, cart)

    icrs = to_icrs(gc)

    if isnan(gsr.distance) 
        distance = NaN
    else
        distance = icrs.distance
    end
    if isnan(gsr.distance) || isnan(gsr.pmra)
        pmra = NaN
    else
        pmra = icrs.pmra
    end
    if isnan(gsr.distance) || isnan(gsr.pmdec)
        pmdec = NaN
    else
        pmdec = icrs.pmdec
    end

    return ICRS(ra=icrs.ra, dec=icrs.dec, distance=distance, pmra=pmra, pmdec=pmdec, radial_velocity=icrs.radial_velocity)
end


function to_icrs(gsr::Cartesian{<:GSR, <:Real})::ICRS
    gc = transform(Galactocentric, gsr)
    return to_icrs(gc)
end


function transform(::Type{<:Cartesian{GSR, <:Real}}, galcen::Galactocentric)
    return _gc_to_gsr_cart(galcen)
end




function _heliocen_to_galcen_position(x_vec::Vector{F}; frame=default_gc_frame) where {F<:Real}
    sun_gc = [-1, 0, 0] .* frame.d

    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    r_vec_prime = R_mat*x_vec .+ sun_gc

    x_gc = H_mat * r_vec_prime
    return x_gc
end


function _heliocen_to_galcen_position(cart::Cartesian; frame=default_gc_frame)
    return _heliocen_to_galcen_position(position(cart); frame=frame)
end


function _gsr_to_galcen_velocity(v_vec::Vector{F}; frame=default_gc_frame) where {F<:Real}
    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    v_gc = H_mat * (R_mat * v_vec)
    return v_gc 
end


function _heliocen_to_galcen_velocity(v_vec::Vector{F}; frame=default_gc_frame) where {F<:Real}
    v_gc = _gsr_to_galcen_velocity(v_vec; frame=frame)
    return v_gc .+ frame.v_sun
end


function _gsr_to_galcen_velocity(cart::Cartesian{<:GSR}; frame=default_gc_frame)
    return _gsr_to_galcen_velocity(velocity(cart); frame=frame)
end



function _heliocen_to_galcen_velocity(cart::Cartesian{<:ICRS}; frame=default_gc_frame)
    return _heliocen_to_galcen_velocity(velocity(cart); frame=frame)
end



function _galcen_to_heliocen_position(x_vec::Vector{F}; frame=default_gc_frame) where {F<:Real}
    sun_gc = [-1, 0, 0] .* frame.d

    R_mat = _coordinate_R_inv(frame)
    H_mat = _coordinate_H_inv(frame)

    r_vec_prime = H_mat * x_vec
    x_icrs = R_mat * (r_vec_prime .- sun_gc)
    return x_icrs
end


function _galcen_to_heliocen_position(galcen::Galactocentric)
    return _galcen_to_heliocen_position(position(galcen), frame=galcen.frame)
end


function _galcen_to_gsr_velocity(v_vec::Vector{F}; frame=default_gc_frame) where {F<:Real}
    R_mat = _coordinate_R_inv(frame)
    H_mat = _coordinate_H_inv(frame)

    v_icrs = R_mat * (H_mat * (v_vec )) # .- frame.v_sun))
    return v_icrs
end

function _galcen_to_heliocen_velocity(v_vec::Vector{F}; frame=default_gc_frame) where {F<:Real}
    return _galcen_to_gsr_velocity(v_vec .- frame.v_sun, frame=frame) 
end

function _galcen_to_gsr_velocity(galcen::Galactocentric)
    return _galcen_to_gsr_velocity(velocity(galcen), frame=galcen.frame)
end

function _galcen_to_heliocen_velocity(galcen::Galactocentric)
    return _galcen_to_heliocen_velocity(velocity(galcen), frame=galcen.frame)
end




