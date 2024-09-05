using Printf
import Base: +, *
import TOML


"""An abstract type for a (sky) coordinate frame"""
abstract type CoordinateFrame end


""" A cartesian coordinate representation.

Requires properties x, y, z, v_x, v_y, v_z
"""
abstract type AbstractCartesian <: CoordinateFrame end


""" A spherical (sky) coordinate representation.

Requires properties ra, dec, distance, pmra, pmdec, radial_velocity"""
abstract type AbstractSkyCoord <: CoordinateFrame end


"""
A type representing a 3D point
"""
struct Point3D{F}
    x::F
    y::F
    z::F

    function Point3D(x, y, z)
        promoted = promote(x, y, z)
        F = typeof(promoted[1])
        return new{F}(promoted...)
    end
end


"""
A type representing a 6D phase space point
"""
struct Point6D{F}
    x::F
    y::F
    z::F
    v_x::F
    v_y::F
    v_z::F

    function Point6D(x, y, z, v_x, v_y, v_z)
        promoted = promote(x, y, z, v_x, v_y, v_z)
        F = typeof(promoted[1])
        return new{F}(promoted...)
    end
end


function Base.getproperty(pp::AbstractCartesian, name::Symbol)
    if :coord ∈ fieldnames(typeof(pp)) && name ∈ fieldnames(typeof(getfield(pp, :coord)))
        return getfield(pp.coord, name)
    else
        return getfield(pp, name)
    end
end


"""
A type representing a sky coordinate (i.e. point on sphere with optional velocity information)
"""
@kwdef struct SkyCoord{F} <: AbstractSkyCoord
    """ Right ascension in degrees """
    ra::F

    """ Declination in degrees """
    dec::F

    """ distance in kpc """
    distance::F = NaN

    """ proper motion in right ascension in mas/yr times cos(dec)"""
    pmra::F = NaN

    """ Proper motion in declination in mas/yr """
    pmdec::F = NaN

    """ Radial velocity in km/s """
    radial_velocity::F = NaN


    function SkyCoord(ra, dec, distance=1, pmra=NaN, pmdec=NaN, radial_velocity=NaN)
        promoted = promote(ra, dec, distance, pmra, pmdec, radial_velocity)
        F = typeof(promoted[1])
        return new{F}(promoted...)
    end
end


function Base.getproperty(coord::AbstractSkyCoord, name::Symbol)
    if :coord ∈ fieldnames(typeof(coord)) && name ∈ fieldnames(typeof(getfield(coord, :coord)))
        return getfield(coord.coord, name)
    else
        return getfield(coord, name)
    end
end


"""
A type representing a Galactic coordinate frame specification
"""
@kwdef struct GalactocentricFrame <: CoordinateFrame
    """ distance from sun to centre of galaxy in kpc"""
    d::F = 8.122  # ± 0.033
    
    """ICRS ra of  galactic centre in degrees"""
    ra::F = 266.4051 

    """ICRS dec of galactic centre in degrees"""
    dec::F = -28.936175 

    """ roll for the galacic frame in degrees """
    η::F = 58.5986320306 # degrees, exact?

    """ height of the sun above the galactic midplane in kpc"""
    z_sun::F = 0.0208 # ± 0.0003

    """ solar motion wrt galactic standard of rest in km/s"""
    v_sun::Vector{F} =  [12.9, 245.6, 7.78] # ± [3.0, 1.4, 0.08]
end



@doc raw"""
The default GC frame as defined by Astropy V4.0

- J2000 'ra' = 17h45m37.224s, 'dec' = -28°56'10.23" (appendix), 'pmracosdec=-3.151±0.018 mas/yr', 'pmdec=-5.547±0.026' (Table 2). from 'https://ui.adsabs.harvard.edu/abs/2004ApJ...616..872R',
- 'distance' = 8.122 ± 0.033 kpc, solar radial velocity = 11 + 1.9 ± 3 km/s.  from  'https://ui.adsabs.harvard.edu/abs/2018A%26A...615L..15G',
- 'z_sun' = 20.8±0.3pc from 'https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.1417B',
- 'roll': Not implemented


'v_sun' is slightly more complicated, relying on parameters from above and using  the procedure described in 'https://ui.adsabs.harvard.edu/abs/2018RNAAS...2..210D'. This results in 
- 'v_sun' = [-12.9 ± 3.0, 245.6 ± 1.4, 7.78 ± 0.08] km/s

basically, the calculation is
```math
v_x = d * (pmra * cos(η) + pmdec * sin(η)) 
v_y = d * (pmra * sin(η) - pmdec * cos(η))
```
"""
const default_gc_frame = GalactocentricFrame()



""" 
ICRS frame
"""
struct ICRS{F} <: AbstractSkyCoord
    coord::SkyCoord{F}

    function ICRS{F}(; ra, dec, distance=NaN, pmra=NaN, pmdec=NaN, radial_velocity=NaN)
        return new{F}(SkyCoord(ra, dec, distance, pmra, pmdec, radial_velocity))
    end
end

function ICRS(; ra, dec, distance=NaN, pmra=NaN, pmdec=NaN, radial_velocity=NaN)
    F = typeof(promote(ra, dec, distance, pmra, pmdec, radial_velocity)[1])
    return ICRS{F}(ra=ra, dec=dec, distance=distance, pmra=pmra, pmdec=pmdec, radial_velocity=radial_velocity)
end



""" 
Galactic standard of rest frame. I.E. ICRS minus the solar velocity.
"""
struct GSR{F} <: AbstractSkyCoord
    coord::SkyCoord{F}
    frame::GalactocentricFrame

    function GSR(; ra, dec, distance=NaN, pmra=NaN, pmdec=NaN, radial_velocity=NaN, frame=default_gc_frame)
        return new{F}(SkyCoord(ra, dec, distance, pmra, pmdec, radial_velocity), frame)
    end
end




""" A Galactocentric point """
struct Galactocentric{F} <: AbstractCartesian
    coord::Point6D{F}
    frame::GalactocentricFrame

    function Galactocentric(x, y, z, v_x=NaN, v_y=NaN, v_z=NaN, frame=default_gc_frame)
        return new{F}(Point6D(x, y, z, v_x, v_y, v_z), frame)
    end
end


function Galactocentric(; x, y, z, v_x, v_y, v_z, frame::GalactocentricFrame = default_gc_frame)
    return Galactocentric(x, y, z, v_x, v_y, v_z, frame)
end


function Galactocentric(pos::Vector{<:Real}, vel::Vector{<:Real} = [NaN, NaN, NaN]; frame::GalactocentricFrame = default_gc_frame)
    if length(pos) != 3
        error("position must be a 3-vector")
    end
    if length(vel) != 3
        error("velocity must be a 3-vector")
    end
    return Galactocentric(pos..., vel..., frame)
end



function Base.show(io::IO, coord::AbstractSkyCoord) 
    print(io, "$(typeof(coord)) at ")
    @printf io "(%4.2f, %4.2f) deg" coord.ra coord.dec
    if coord.distance !== NaN
        @printf io ", "
        @printf io "d = %4.2f kpc" coord.distance
    end
    if coord.pmra !== NaN && coord.pmdec !== NaN
        @printf io ", "
        @printf io "μ = (%4.2f, %4.2f) mas/yr" coord.pmra coord.pmdec
    end

    if coord.radial_velocity !== NaN
        @printf io ", "
        @printf io "v_los = %4.2f km/s" coord.radial_velocity
    end

    return io
end




"""
    Cartesian

A Cartesian representation of a given skycoord type
"""
struct Cartesian{T<:CoordinateFrame, F<:Real} <: AbstractCartesian 
    coord::Point6D{F}
    skycoord::Union{T, Nothing}

    function Cartesian{T, F}(x, y, z, v_x, v_y, v_z, coord=nothing) where {T<:CoordinateFrame, F<:Real}
        return new{T, F}(Point6D(x, y, z, v_x, v_y, v_z), coord)
    end
end

function Cartesian{T, F}(; x, y, z, v_x=NaN, v_y=NaN, v_z=NaN, coord=nothing) where {T<:CoordinateFrame, F<:Real}
    return Cartesian{T, F}(x, y, z, v_x, v_y, v_z, coord)
end

function Cartesian{T}(; x, y, z, v_x=NaN, v_y=NaN, v_z=NaN, coord=nothing) where {T<:CoordinateFrame}
    F = typeof(promote(x, y, z, v_x, v_y, v_z)[1])
    return Cartesian{T, F}(x, y, z, v_x, v_y, v_z, coord)
end


function Cartesian{T, F}(pos::AbstractVector{F}, vel::AbstractVector{F} = [NaN, NaN, NaN], coord=nothing) where {T<:CoordinateFrame, F<:Real}
    if length(pos) != 3
        error("position must be a 3-vector")
    end
    if length(vel) != 3
        error("velocity must be a 3-vector")
    end
    return Cartesian{T, F}(pos..., vel..., coord)
end


function Cartesian{T}(pos::AbstractVector{F}, vel::AbstractVector{F} = [NaN, NaN, NaN], coord=nothing) where {T<:CoordinateFrame, F<:Real}
    return Cartesian{T, F}(pos, vel, coord)
end



"""
    position_of(pp::AbstractCartesian)

Returns [x, y, z], the position of phase point.
"""
function position_of(pp::AbstractCartesian)
    return [pp.x, pp.y, pp.z]
end


"""
    velocity_of(pp::AbstractPhasePoint)

Returns [v_x, v_y, v_z], the velocity of phase point.
"""
function velocity_of(pp::AbstractCartesian)
    return [pp.v_x, pp.v_y, pp.v_z]
end



function Base.show(io::IO, pp::Cartesian{T}) where {T}
    print(io, "$T point at ")
    @printf io "(%4.2f, %4.2f, %4.2f) kpc, " pp.x pp.y pp.z

    if pp.v_x !== NaN && pp.v_y !== NaN && pp.v_z !== NaN
        @printf io "(%4.2f, %4.2f, %4.2f) km/s" pp.v_x pp.v_y pp.v_z
    end
    return io
end






"""
    rand_coord(obs::ICRS, err::ICRS)

Generate a random coordinate based on the observed coordinate and its error
assumed to be normally distributed.
"""
function rand_coord(obs::ICRS, err::ICRS)
    return ICRS(
        ra = obs.ra + randn() * err.ra,
        dec = obs.dec + randn() * err.dec,
        pmra = obs.pmra + randn() * err.pmra,
        pmdec = obs.pmdec + randn() * err.pmdec,
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



function coord_from_file(filename::String)
    args = TOML.parsefile(filename)

    labels = ["ra", "dec", "pmra", "pmdec", "radial_velocity", "distance"]

    kwargs = Dict{Symbol, Float64}()
    for label in labels
        if haskey(args, label)
            kwargs[Symbol(label)] = args[label]
        end
    end

    return ICRS{Float64}(;kwargs...)
end

function coord_err_from_file(filename::String)
    args = TOML.parsefile(filename)

    labels = ["ra_err", "dec_err", "pmra_err", "pmdec_err", "radial_velocity_err", "distance_err"]

    kwargs = Dict{Symbol, Float64}()

    for label in labels
        if haskey(args, label)
            k = replace(label, "_err" => "")
            kwargs[Symbol(k)] = args[label]
        end
    end

    return ICRS{Float64}(;kwargs...)
end
