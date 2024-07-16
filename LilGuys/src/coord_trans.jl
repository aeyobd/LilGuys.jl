import Base: @kwdef



"""
    transform(T, obs)

Transforms a coordinate of type `obs` to an object of type `T`.


Notes: if not overridden, transform first attempts to convert into ICRS for unknown conversions.
All coordinate types should define an initializer which takes a ICRS object and a method `to_icrs` which
"""
function transform(::Type{T}, obs::CoordinateFrame; kwargs...) where T<:CoordinateFrame
    return T(to_icrs(obs); kwargs...)
end


# required implementation
function transform(::Type{T}, obs::ICRS) where T<:CoordinateFrame
    return T(obs)
end

# required implementation
function transform(::Type{ICRS}, obs::CoordinateFrame)
    return to_icrs(obs)
end


# identity
function transform(::Type{T}, obs::T) where T
    return obs
end


# automatic cartesian conversion
function transform(::Type{Cartesian{T}}, obs::T) where T <: CoordinateFrame
    return Cartesian(obs)
end

# automatic cartesian conversion
function transform(::Type{T}, obs::Cartesian{T}) where T <: CoordinateFrame
    return to_sky(obs)
end


# to make ICRS behave with multiple dispatch
function transform(::Type{ICRS}, obs::Cartesian{ICRS})
    return to_sky(obs)
end

function transform(::Type{Cartesian{ICRS}}, obs::ICRS)
    return Cartesian(obs)
end


function transform(::Type{ICRS}, obs::ICRS)
    return obs
end

function transform(::Type{Cartesian{ICRS}}, obs::Cartesian{ICRS})
    return obs
end

function to_icrs(cart::Cartesian{ICRS})
    return to_sky(cart)
end



# ========= Cartesian =========

function to_sky(coord::Cartesian{T}; kwargs...) where T <: CoordinateFrame
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



function Cartesian{T}(sc::T) where T<:SkyCoord
    position = sky_to_cartesian(sc.ra, sc.dec, sc.distance)
    velocity = _observation_to_cartesian_velocity(sc)
    return Cartesian{T}(position, velocity, sc)
end


function Cartesian(sc::T) where T<:SkyCoord
    return Cartesian{T}(sc)
end


function _observation_to_cartesian_velocity(obs::SkyCoord)
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


function Galactocentric(icrs::ICRS; kwargs...)
    cart = Cartesian(icrs)
    return Galactocentric(cart; kwargs...)
end


function Galactocentric(cart::Cartesian{ICRS}; frame=default_gc_frame)
    x_gc = _heliocen_to_galcen_position(cart, frame)
    v_gc = _heliocen_to_galcen_velocity(cart, frame)

    return Galactocentric(x_gc, v_gc; frame=frame)
end


function to_icrs(galcen::Galactocentric)::ICRS
    cart = transform(Cartesian{ICRS}, galcen)
    return to_sky(cart)
end


function transform(::Type{Cartesian{ICRS}}, galcen::Galactocentric; frame=default_gc_frame)
    x_icrs = _galcen_to_heliocen_position(galcen)
    v_icrs = _galcen_to_heliocen_velocity(galcen)

    return Cartesian{ICRS}(x_icrs, v_icrs)
end



# ========= GSR transformations =========


function GSR(icrs::ICRS; kwargs...)
    cart = Cartesian{GSR}(icrs; kwargs...)
    return to_sky(cart)
end


function Cartesian{GSR}(icrs::ICRS; kwargs...)
    gc = Galactocentric(icrs)
    return transform(Cartesian{GSR}, gc; kwargs...)
end


function transform(::Type{Galactocentric}, cart::Cartesian{GSR}, frame=default_gc_frame)
    x_gc = _heliocen_to_galcen_position(cart, frame)
    v_gc = _gsr_to_galcen_velocity(cart, frame)

    return Galactocentric(x_gc, v_gc, frame=frame)
end


function to_icrs(gsr::GSR)::ICRS
    cart = Cartesian(gsr)
    return to_icrs(cart)
end


function to_icrs(gsr::Cartesian{GSR})::ICRS
    gc = transform(Galactocentric, gsr)
    return to_icrs(gc)
end


function transform(::Type{Cartesian{GSR}}, galcen::Galactocentric)
    x_gsr = _galcen_to_heliocen_position(galcen)
    v_gsr = _galcen_to_gsr_velocity(galcen)

    return Cartesian{GSR}(x_gsr, v_gsr)
end










function _heliocen_to_galcen_position(x_vec::Vector{F}, frame=GalactocentricFrame())
    sun_gc = [-1, 0, 0] .* frame.d

    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    r_vec_prime = R_mat*x_vec .+ sun_gc

    x_gc = H_mat * r_vec_prime
    return x_gc
end


function _heliocen_to_galcen_position(cart::Cartesian, frame=default_gc_frame)
    return _heliocen_to_galcen_position(position_of(cart), frame)
end


function _gsr_to_galcen_velocity(v_vec::Vector{F}, frame=default_gc_frame)
    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    v_gc = H_mat * (R_mat * v_vec)
    return v_gc 
end


function _heliocen_to_galcen_velocity(v_vec::Vector{F}, frame=default_gc_frame)
    v_gc = _gsr_to_galcen_velocity(v_vec, frame)
    return v_gc .+ frame.v_sun
end


function _gsr_to_galcen_velocity(cart::Cartesian{GSR}, frame=default_gc_frame)
    return _gsr_to_galcen_velocity(velocity_of(cart), frame)
end

function _heliocen_to_galcen_velocity(cart::Cartesian{ICRS}, frame=default_gc_frame)
    return _heliocen_to_galcen_velocity(velocity_of(cart), frame)
end



function _galcen_to_heliocen_position(x_vec::Vector{F}, frame=default_gc_frame)
    sun_gc = [-1, 0, 0] .* frame.d

    R_mat = _coordinate_R_inv(frame)
    H_mat = _coordinate_H_inv(frame)

    r_vec_prime = H_mat * x_vec
    x_icrs = R_mat * (r_vec_prime .- sun_gc)
    return x_icrs
end


function _galcen_to_heliocen_position(galcen::Galactocentric)
    return _galcen_to_heliocen_position(position_of(galcen), galcen.frame)
end


function _galcen_to_gsr_velocity(v_vec::Vector{F}, frame=default_gc_frame)
    R_mat = _coordinate_R_inv(frame)
    H_mat = _coordinate_H_inv(frame)

    v_icrs = R_mat * (H_mat * (v_vec )) # .- frame.v_sun))
    return v_icrs
end

function _galcen_to_heliocen_velocity(v_vec::Vector{F}, frame=default_gc_frame)
    return _galcen_to_gsr_velocity(v_vec .- frame.v_sun, frame) 
end

function _galcen_to_gsr_velocity(galcen::Galactocentric)
    return _galcen_to_gsr_velocity(velocity_of(galcen), galcen.frame)
end

function _galcen_to_heliocen_velocity(galcen::Galactocentric)
    return _galcen_to_heliocen_velocity(velocity_of(galcen), galcen.frame)
end


"""Galactocentric rotation matrix"""
function _coordinate_R(frame::GalactocentricFrame)
    η = deg2rad(frame.η)
    α = deg2rad(frame.ra)
    δ = deg2rad(frame.dec)

    return Rx_mat(η) * Ry_mat(δ) * Rz_mat(-α)
end


"""Inverse Galactocentric rotation matrix"""
function _coordinate_R_inv(frame::GalactocentricFrame)
    η = deg2rad(frame.η)
    α = deg2rad(frame.ra)
    δ = deg2rad(frame.dec)

    return Rz_mat(α) * Ry_mat(-δ) * Rx_mat(-η)
end


"""Galactocentric height rotation matrix"""
function _coordinate_H(frame::GalactocentricFrame)
    θ = asin(frame.z_sun / frame.d)
    return Ry_mat(θ)
end


"""Inverse Galactocentric height rotation matrix"""
function _coordinate_H_inv(frame::GalactocentricFrame)
    θ = asin(frame.z_sun / frame.d)
    return Ry_mat(-θ)
end

