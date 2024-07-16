import Base: @kwdef


function transform(::Type{T}, obs::T) where T
    return obs
end


function transform(::Type{PhasePoint{T}}, obs::SkyCoord{T}) where T
    x, y, z = _observation_to_cartesian_position(obs)
    if any(isnan, [obs.pmra, obs.pmdec, obs.radial_velocity])
        v_x, v_y, v_z = NaN, NaN, NaN
    else
        v_x, v_y, v_z = _observation_to_cartesian_velocity(obs)
    end

    return PhasePoint{T}(x, y, z, v_x, v_y, v_z)
end


function transform(::Type{SkyCoord{T}}, cart::PhasePoint{T}) where T
    ra, dec, r = _cartesian_to_observation_position(cart)
    if any(isnan, [cart.v_x, cart.v_y, cart.v_z])
        pmra, pmdec, radial_velocity = NaN, NaN, NaN
    else
        pmra, pmdec, radial_velocity = _cartesian_to_observation_velocity(cart)
    end

    return SkyCoord{T}(ra, dec, r, pmra, pmdec, radial_velocity)
end


function transform(::Type{Galactocentric}, cart::ICRS_Cartesian, frame=default_gc_frame)
    x_gc = _heliocen_to_galcen_position(cart, frame)
    v_gc = _heliocen_to_galcen_velocity(cart, frame)

    return Galactocentric(x_gc, v_gc)
end


function transform(::Type{ICRS_Cartesian}, galcen::Galactocentric, frame=default_gc_frame)
    x_icrs = _galcen_to_heliocen_position(galcen)
    v_icrs = _galcen_to_heliocen_velocity(galcen)

    return ICRS_Cartesian(x_icrs, v_icrs)
end



function transform(::Type{Galactocentric}, cart::GSR_Cartesian, frame=default_gc_frame)
    x_gc = _heliocen_to_galcen_position(cart, frame)
    v_gc = _gsr_to_galcen_velocity(cart, frame)

    return Galactocentric(x_gc, v_gc)
end


function transform(::Type{GSR_Cartesian}, galcen::Galactocentric)
    x_gsr = _galcen_to_heliocen_position(galcen)
    v_gsr = _galcen_to_gsr_velocity(galcen)

    return GSR_Cartesian(x_gsr, v_gsr)
end


function transform(::Type{Galactocentric}, obs::ICRS, frame=default_gc_frame)
    cart = transform(ICRS_Cartesian, obs)
    return transform(Galactocentric, cart, frame)
end



function transform(::Type{ICRS}, galcen::Galactocentric, frame=default_gc_frame)
    cart = transform(ICRS_Cartesian, galcen, frame)
    return transform(ICRS, cart)
end


function transform(::Type{Galactocentric}, obs::GSR)
    cart = transform(GSR_Cartesian, obs)
    return transform(Galactocentric, cart, default_gc_frame)
end



function transform(::Type{GSR}, galcen::Galactocentric)
    cart = transform(GSR_Cartesian, galcen)
    return transform(GSR, cart)
end


function transform(::Type{GSR_Cartesian}, obs::ICRS_Cartesian)
    gc = transform(Galactocentric, obs)
    return transform(GSR_Cartesian, gc)
end


function transform(::Type{ICRS_Cartesian}, obs::GSR_Cartesian)
    gc = transform(Galactocentric, obs)
    return transform(ICRS_Cartesian, gc)
end


function transform(::Type{ICRS}, obs::GSR)
    gc = transform(Galactocentric, obs)
    return transform(ICRS, gc)
end


function transform(::Type{GSR}, obs::ICRS)
    gc = transform(Galactocentric, obs)
    return transform(GSR, gc)
end




function _observation_to_cartesian_position(obs::SkyCoord)
    return obs.distance * unit_vector(obs.ra, obs.dec)
end

function _observation_to_cartesian_velocity(obs::SkyCoord)
    rv = obs.radial_velocity
    α = obs.ra
    δ = obs.dec
    v_α_cosδ = obs.pmra * kpc_mas_yr * obs.distance
    v_δ = obs.pmdec  * kpc_mas_yr * obs.distance
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

    return cartesian_to_sky(x, y, z)
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


function _heliocen_to_galcen_position(x_vec::Vector{F}, frame=GalactocentricFrame())
    sun_gc = [-1, 0, 0] .* frame.d

    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    r_vec_prime = R_mat*x_vec .+ sun_gc

    x_gc = H_mat * r_vec_prime
    return x_gc
end


function _heliocen_to_galcen_position(cart::PhasePoint, frame=default_gc_frame)
    return _heliocen_to_galcen_position(get_position(cart), frame)
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


function _gsr_to_galcen_velocity(cart::GSR_Cartesian, frame=default_gc_frame)
    return _gsr_to_galcen_velocity(get_velocity(cart), frame)
end

function _heliocen_to_galcen_velocity(cart::ICRS_Cartesian, frame=default_gc_frame)
    return _heliocen_to_galcen_velocity(get_velocity(cart), frame)
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
    return _galcen_to_heliocen_position(get_position(galcen), galcen.frame)
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
    return _galcen_to_gsr_velocity(get_velocity(galcen), galcen.frame)
end

function _galcen_to_heliocen_velocity(galcen::Galactocentric)
    return _galcen_to_heliocen_velocity(get_velocity(galcen), galcen.frame)
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

