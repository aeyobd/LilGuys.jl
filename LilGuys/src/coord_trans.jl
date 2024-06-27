import Base: @kwdef

const kpc_mas_yr = 4.740470463533348 # km/s


@kwdef struct GalactocentricFrame
    d::F = 8.122 # kpc
    ra::F = 266.4051 # Sag A* ra in deg
    dec::F = -28.936175 # Sag A* dec in deg
    η::F = 58.5986320306 # degrees
    z_sun::F = 0.0208 # kpc, solar height above galactic midplane
    v_sun::Vector{F} =  [12.9, 245.6, 7.78]
end

function gc_rest_frame()
    return GalactocentricFrame(; v_sun=[0., 0., 0.])
end



function transform(::Type{T}, obs::T) where T
    return obs
end


function transform(::Type{PhasePoint{T}}, obs::SkyCoord{T}) where T
    x, y, z = _observation_to_cartesian_position(obs)
    if any(isnan, [obs.pm_ra, obs.pm_dec, obs.radial_velocity])
        v_x, v_y, v_z = NaN, NaN, NaN
    else
        v_x, v_y, v_z = _observation_to_cartesian_velocity(obs)
    end

    return PhasePoint{T}(x, y, z, v_x, v_y, v_z)
end


function transform(::Type{SkyCoord{T}}, cart::PhasePoint{T}) where T
    ra, dec, r = _cartesian_to_observation_position(cart)
    if any(isnan, [cart.v_x, cart.v_y, cart.v_z])
        pm_ra, pm_dec, radial_velocity = NaN, NaN, NaN
    else
        pm_ra, pm_dec, radial_velocity = _cartesian_to_observation_velocity(cart)
    end

    return SkyCoord{T}(ra, dec, r, pm_ra, pm_dec, radial_velocity)
end


function transform(::Type{Galactocentric}, cart::ICRS_Cartesian, frame=GalactocentricFrame())
    x_gc = _cartesian_to_galcen_position(cart, frame)
    v_gc = _cartesian_to_galcen_velocity(cart, frame)

    return Galactocentric(x_gc, v_gc)
end


function transform(::Type{ICRS_Cartesian}, galcen::Galactocentric, frame=GalactocentricFrame())
    x_icrs = _galcen_to_cartesian_position(galcen, frame)
    v_icrs = _galcen_to_cartesian_velocity(galcen, frame)

    return ICRS_Cartesian(x_icrs, v_icrs)
end



function transform(::Type{Galactocentric}, cart::HelioRest_Cartesian, frame=gc_rest_frame())
    x_gc = _cartesian_to_galcen_position(cart, frame)
    v_gc = _cartesian_to_galcen_velocity(cart, frame)

    return Galactocentric(x_gc, v_gc)
end


function transform(::Type{HelioRest_Cartesian}, galcen::Galactocentric, frame=gc_rest_frame())
    x_icrs = _galcen_to_cartesian_position(galcen, frame)
    v_icrs = _galcen_to_cartesian_velocity(galcen, frame)

    return HelioRest_Cartesian(x_icrs, v_icrs)
end


function transform(::Type{Galactocentric}, obs::ICRS, frame=GalactocentricFrame())
    cart = transform(ICRS_Cartesian, obs)
    return transform(Galactocentric, cart, frame)
end



function transform(::Type{ICRS}, galcen::Galactocentric, frame=GalactocentricFrame())
    cart = transform(ICRS_Cartesian, galcen, frame)
    return transform(ICRS, cart)
end


function transform(::Type{Galactocentric}, obs::HelioRest, frame=gc_rest_frame())
    cart = transform(HelioRest_Cartesian, obs)
    return transform(Galactocentric, cart, frame)
end



function transform(::Type{HelioRest}, galcen::Galactocentric, frame=gc_rest_frame())
    cart = transform(HelioRest_Cartesian, galcen, frame)
    return transform(HelioRest, cart)
end


function transform(::Type{HelioRest_Cartesian}, obs::ICRS_Cartesian)
    gc = transform(Galactocentric, obs)
    return transform(HelioRest_Cartesian, gc)
end


function transform(::Type{ICRS_Cartesian}, obs::HelioRest_Cartesian)
    gc = transform(Galactocentric, obs)
    return transform(ICRS_Cartesian, gc)
end


function transform(::Type{ICRS}, obs::HelioRest)
    gc = transform(Galactocentric, obs)
    return transform(ICRS, gc)
end


function transform(::Type{HelioRest}, obs::ICRS)
    gc = transform(Galactocentric, obs)
    return transform(HelioRest, gc)
end




function _observation_to_cartesian_position(obs::SkyCoord)
    return obs.distance * unit_vector(obs.ra, obs.dec)
end

function _observation_to_cartesian_velocity(obs::SkyCoord)
    rv = obs.radial_velocity
    α = obs.ra
    δ = obs.dec
    v_α_cosδ = obs.pm_ra * kpc_mas_yr * obs.distance
    v_δ = obs.pm_dec  * kpc_mas_yr * obs.distance
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




function _cartesian_to_observation_position(cart::PhasePoint)
    x, y, z = cart.x, cart.y, cart.z

    R = sqrt(x^2 + y^2)
    r = sqrt(x^2 + y^2 + z^2)

    ra = mod(atand(y, x), 360) # atan2 -> 0 to 360
    dec = atand(z / R) # atan -> -90 to 90

    return [ra, dec, r]
end


function _cartesian_to_observation_velocity(cart::PhasePoint)
    x, y, z = cart.x, cart.y, cart.z
    v_x, v_y, v_z = cart.v_x, cart.v_y, cart.v_z

    R = sqrt(x^2 + y^2)
    r = sqrt(x^2 + y^2 + z^2)

    v_r = (x*v_x + y*v_y + z*v_z) / r
    v_α = (x*v_y - y*v_x) / (r*R)
    v_δ = v_z*R/r^2 - (x*v_x + y*v_y) * z/(r^2 * R)

    return [v_α/kpc_mas_yr, v_δ/kpc_mas_yr, v_r]
end

function _cartesian_to_galcen_position(x_vec::Vector{F}, frame=GalactocentricFrame())
    sun_gc = [-1, 0, 0] .* frame.d

    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    r_vec_prime = R_mat*x_vec .+ sun_gc

    x_gc = H_mat * r_vec_prime
    return x_gc
end

function _cartesian_to_galcen_position(cart::PhasePoint, frame=GalactocentricFrame())
    return _cartesian_to_galcen_position(cart.position, frame)
end

function _cartesian_to_galcen_velocity(v_vec::Vector{F}, frame=GalactocentricFrame())
    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    v_gc = H_mat * (R_mat * v_vec)
    v_gc .+= frame.v_sun
    return v_gc 
end

function _cartesian_to_galcen_velocity(cart::PhasePoint, frame=GalactocentricFrame())
    return _cartesian_to_galcen_velocity(cart.velocity, frame)
end



function _galcen_to_cartesian_position(x_vec::Vector{F}, frame=GalactocentricFrame())
    sun_gc = [-1, 0, 0] .* frame.d

    R_mat = _coordinate_R_inv(frame)
    H_mat = _coordinate_H_inv(frame)

    r_vec_prime = H_mat * x_vec
    x_icrs = R_mat * (r_vec_prime .- sun_gc)
    return x_icrs
end

function _galcen_to_cartesian_position(galcen::Galactocentric, frame=GalactocentricFrame())
    return _galcen_to_cartesian_position(galcen.position, frame)
end


function _galcen_to_cartesian_velocity(v_vec::Vector{F}, frame=GalactocentricFrame())
    R_mat = _coordinate_R_inv(frame)
    H_mat = _coordinate_H_inv(frame)

    v_icrs = R_mat * (H_mat * (v_vec .- frame.v_sun))
    return v_icrs
end


function _galcen_to_cartesian_velocity(galcen::Galactocentric, frame=GalactocentricFrame())
    return _galcen_to_cartesian_velocity(galcen.velocity, frame)
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


"""
    rand_coord(obs::ICRS, err::ICRS)

Generate a random coordinate based on the observed coordinate and its error
assumed to be normally distributed.
"""
function rand_coord(obs::ICRS, err::ICRS)
    return ICRS(
        ra = obs.ra,
        dec = obs.dec,
        pm_ra = obs.pm_ra + randn() * err.pm_ra,
        pm_dec = obs.pm_dec + randn() * err.pm_dec,
        radial_velocity = obs.radial_velocity + randn() * err.radial_velocity,
        distance = obs.distance + randn() * err.distance,
       )
end


"""
    rand_coords(obs::ICRS, err::ICRS, N::Int)

Generate N random coordinates based on the observed coordinate and its error
assumed to be normally distributed.
"""
function rand_coords(obs::ICRS, err::ICRS, N::Int)
    return [rand_coord(obs, err) for _ in 1:N]
end
