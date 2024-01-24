using Printf
import Base: @kwdef

# Todo: rewrite astropy coordinates with naitive Julia

@kwdef struct PhasePoint
    x::F
    y::F
    z::F
    v_x::F
    v_y::F
    v_z::F
end

@kwdef struct Observation
    ra::F
    dec::F
    distance::F
    pm_ra::F
    pm_dec::F
    radial_velocity::F
end


function PhasePoint(pos::Vector{F}, vel::Vector{F})
    if length(pos) != 3
        error("position must be a 3-vector")
    end
    if length(vel) != 3
        error("velocity must be a 3-vector")
    end
    return PhasePoint(pos..., vel...)
end


const kpc_mas_yr = 4.740470463533348 # km/s


@kwdef struct GalactocentricFrame
    d::F = 8.122 # kpc
    ra::F = 266.4051 # Sag A* ra in deg
    dec::F = -28.936175 # Sag A* dec in deg
    η::F = 58.5986320306 # degrees
    z_sun::F = 0.0208 # kpc, solar height above galactic midplane
    v_sun::Vector{F} =  [12.9, 245.6, 7.78]
end


function Base.getproperty(pp::PhasePoint, sym::Symbol)
    if sym == :position
        return [pp.x, pp.y, pp.z]
    elseif sym == :velocity
        return [pp.v_x, pp.v_y, pp.v_z]
    elseif sym ∈ fieldnames(PhasePoint)
        return getfield(pp, sym)
    else
        error("$(sym) is not a valid field")
    end
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

function _xyz_icrs(ra, dec, d)
    x = d * cosd(dec) * cosd(ra)
    y = d * cosd(dec) * sind(ra)
    z = d * sind(dec)
    return [x, y, z]
end


function _v_xyz_icrs(obs::Observation)
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


function _xyz_icrs(obs::Observation)
    pos = _xyz_icrs(obs.ra, obs.dec, obs.distance)
    vel = _v_xyz_icrs(obs)
    return pos, vel
end


function _from_xyz_icrs(x, y, z)
    R = sqrt(x^2 + y^2)

    ra = atand(y / x)
    dec = atand(z / R)
    if y < 0
        ra += 180
    end
    d = sqrt(x^2 + y^2 + z^2)
    return [ra, dec, d]
end


function _from_v_xyz_icrs(pp::PhasePoint)
end


function to_galcen(obs::Observation, frame=GalactocentricFrame())
    x_vec, v_vec = _xyz_icrs(obs)

    x_gc = _to_galcen(x_vec, frame)
    v_gc = _to_galcen_v(v_vec, frame)

    return PhasePoint(x_gc / R0, v_gc / V0)
end


function _to_galcen(x_vec::Vector{F}, frame=GalactocentricFrame())
    sun_gc = [-1, 0, 0] .* frame.d

    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    r_vec_prime = R_mat*x_vec .+ sun_gc

    x_gc = H_mat * r_vec_prime
    return x_gc
end


function _to_galcen_v(v_vec::Vector{F}, frame=GalactocentricFrame())
    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    v_gc = H_mat * (R_mat * v_vec)
    v_gc .+= frame.v_sun
    return v_gc 
end

function _from_galcen_v(v_vec::Vector{F}, frame=GalactocentricFrame())
    R_mat = _coordinate_R(frame)
    H_mat = _coordinate_H(frame)

    v_icrs = R_mat * (H_mat * v_vec)
    v_icrs .-= frame.v_sun
    return v_icrs
end

"""
Transforms a PhasePoint in Galactocentric coordinates to ICRS xyz coordinates
Returns ([x, y, z], [vx, vy, vz]) in kpc and km/s
"""
function _galactocentric_to_icrs_cartesian(phase::PhasePoint, frame=GalactocentricFrame())


end


"""
Returns a list of observations based on snapshot particles
"""
function to_sky(snap::Snapshot; invert_velocity::Bool=false, verbose::Bool=false)
    observations = Observation[]

    for i in 1:length(snap)
        if verbose
            print("converting $(i)/($(length(snap))\r")
        end
        pos = snap.positions[:, i]
        vel = snap.velocities[:, i]
        if invert_velocity
            vel *=-1
        end
        phase = PhasePoint(pos, vel)
        obs = to_sky(phase)
        push!(observations, obs)
    end
    return observations
end



function rand_coord(obs::Observation, err::Observation)
    return Observation(
        ra = obs.ra,
        dec = obs.dec,
        pm_ra = obs.pm_ra + randn() * err.pm_ra,
        pm_dec = obs.pm_dec + randn() * err.pm_dec,
        radial_velocity = obs.radial_velocity + randn() * err.radial_velocity,
        distance = obs.distance + randn() * err.distance,
       )
end

function rand_coords(obs::Observation, err::Observation, N::Int)
    return [rand_coord(obs, err) for _ in 1:N]
end


function Base.show(io::IO, pp::PhasePoint)
    x = pp.position * R0
    v = pp.velocity * V0
    @printf io "(%4.2f, %4.2f, %4.2f) kpc, " x...
    @printf io "(%4.2f, %4.2f, %4.2f) km/s" v...
    return io
end


function Base.show(io::IO, obs::Observation)
    @printf io "observation at (%4.2f, %4.2f)"  obs.ra obs.dec
    return io
end




