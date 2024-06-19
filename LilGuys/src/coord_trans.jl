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


"""
    to_tangent(α, δ, α_0, δ_0)

Computes the tangent plane coordinates of a point (α, δ) with respect to a reference point (α_0, δ_0).
"""
function to_tangent(α, δ, α_0, δ_0)
	denom = @. (sind(δ) * sind(δ_0) 
		+ cosd(δ) * cosd(δ_0) * cosd(α - α_0)
	)
	
	eta_num = @. (sind(δ_0) * cosd(δ) * cosd(α-α_0)
		-cosd(δ_0) * sind(δ) 
	)
	
	xi_num = @. cosd(δ) * sind(α - α_0)
	
	xi = @. rad2deg(xi_num/denom)
	eta = @. -rad2deg(eta_num / denom)

	return xi, eta
end




function _unit_vector(ra::F, dec::F)
    x = cosd(dec) * cosd(ra)
    y = cosd(dec) * sind(ra)
    z = sind(dec)
    return [x, y, z]
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
    return obs.distance * _unit_vector(obs.ra, obs.dec)
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


function _coordinate_H(θ::F)
    H = [
        cosd(θ) 0 sind(θ)
        0 1 0
        -sind(θ) 0 cosd(θ)
    ]
    return H
end

function _coordinate_H(frame::GalactocentricFrame)
    θ = asind(frame.z_sun / frame.d)
    return _coordinate_H(θ)
end

function _coordinate_H_inv(frame::GalactocentricFrame)
    θ = asind(frame.z_sun / frame.d)
    return _coordinate_H(-θ)
end


function Rx_mat(u::Real)
    c = cos(u)
    s = sin(u)
    return [1  0  0
            0  c  s
            0 -s  c]
end


function Ry_mat(v::Real)
    c = cos(u)
    s = sin(u)
    return [c  0  s
            0  1  0
           -s  0  c]
end

function Rz_mat(w::Real)
    c = cos(u)
    s = sin(u)

    return [c -s  0
            s  c  0
            0  0  1]
end


function R_mat(u::Real, v::Real, w::Real)
    return Rz_mat(w) * Ry_mat(v) * Rx_mat(u)
end


function _coordinate_R1(δ::F)
    R1 = [cosd(δ) 0 sind(δ)
        0 1 0
        -sind(δ) 0 cosd(δ)]
    return R1
end

function _coordinate_R2(α::F)
    R2 =[cosd(α) sind(α) 0
          -sind(α) cosd(α) 0
          0 0 1]
    return R2
end

function _coordinate_R3(η::F)
    R3 = [1 0 0
        0 cosd(η) sind(η)
        0 -sind(η) cosd(η)]
    return R3
end

function _coordinate_R(frame::GalactocentricFrame)
    return _coordinate_R3(frame.η) * _coordinate_R1(frame.dec) * _coordinate_R2(frame.ra)
end

function _coordinate_R_inv(frame::GalactocentricFrame)
    return _coordinate_R2(-frame.ra) * _coordinate_R1(-frame.dec) * _coordinate_R3(-frame.η)
end






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

function rand_coords(obs::ICRS, err::ICRS, N::Int)
    return [rand_coord(obs, err) for _ in 1:N]
end
