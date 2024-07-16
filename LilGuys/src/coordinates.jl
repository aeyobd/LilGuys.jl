using Printf
import Base: +, *

"""An abstract type for a (sky) coordinate frame"""
abstract type CoordinateFrame end
abstract type AbstractPhasePoint <: CoordinateFrame end



@kwdef struct GalactocentricFrame
    """ distance from sun to centre of galaxy in kpc"""
    d::F = 8.122 
    
    """ICRS ra of  galactic centre in degrees"""
    ra::F = 266.4051 

    """ICRS dec of galactic centre in degrees"""
    dec::F = -28.936175 

    """ roll for the galacic frame in degrees """
    η::F = 58.5986320306 # degrees

    """ height of the sun above the galactic midplane in kpc"""
    z_sun::F = 0.0208 

    """ solar motion wrt galactic standard of rest in km/s"""
    v_sun::Vector{F} =  [12.9, 245.6, 7.78]
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

@kwdef struct PhasePoint{T} <: AbstractPhasePoint
    x::F
    y::F
    z::F
    v_x::F = NaN
    v_y::F = NaN
    v_z::F = NaN
end



function PhasePoint{T}(pos::Vector{F}, vel::Vector{F} = [NaN, NaN, NaN]) where T
    if length(pos) != 3
        error("position must be a 3-vector")
    end
    if length(vel) != 3
        error("velocity must be a 3-vector")
    end
    return PhasePoint{T}(pos..., vel...)
end


"""
    get_position(pp::AbstractPhasePoint)

Returns [x, y, z], the position of phase point.
"""
function get_position(pp::AbstractPhasePoint)
    return [pp.x, pp.y, pp.z]
end


"""
    get_velocity(pp::AbstractPhasePoint)

Returns [v_x, v_y, v_z], the velocity of phase point.
"""
function get_velocity(pp::AbstractPhasePoint)
    return [pp.v_x, pp.v_y, pp.v_z]
end


""" 
A sky coordinate in an arbitrary frame
"""
@kwdef struct SkyCoord{T} <: CoordinateFrame
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
end



""" A Galactocentric point """
@kwdef struct Galactocentric <: AbstractPhasePoint
    x::F
    y::F
    z::F
    v_x::F
    v_y::F
    v_z::F
    frame::GalactocentricFrame = default_gc_frame
end

function Galactocentric(pos::Vector{F}, vel::Vector{F} = [NaN, NaN, NaN], frame::GalactocentricFrame = default_gc_frame) 
    if length(pos) != 3
        error("position must be a 3-vector")
    end
    if length(vel) != 3
        error("velocity must be a 3-vector")
    end
    return Galactocentric(pos..., vel..., frame)
end



""" A point in the simulation frame """
const SimPoint = PhasePoint{:SimPoint}

""" A point in the ICRS frame """
const ICRS_Cartesian = PhasePoint{:ICRS} 

""" A point in the GSR frame """
const GSR_Cartesian = PhasePoint{:GSR}

""" Galactic standard of rest """
const GSR = SkyCoord{:GSR}

""" International Celestial Reference System """
const ICRS = SkyCoord{:ICRS}
    


function Base.show(io::IO, pp::PhasePoint{T}) where T
    print(io, "$T point at ")
    @printf io "(%4.2f, %4.2f, %4.2f) kpc, " pp.x pp.y pp.z

    if pp.v_x !== NaN && pp.v_y !== NaN && pp.v_z !== NaN
        @printf io "(%4.2f, %4.2f, %4.2f) km/s" pp.v_x pp.v_y pp.v_z
    end
    return io
end


function Base.show(io::IO, coord::SkyCoord{T}) where T
    print(io, "$T at ")
    @printf io "(%4.2f, %4.2f) deg" coord.ra coord.dec
    if coord.distance !== NaN
        @printf ", "
        @printf io "d = %4.2f kpc" coord.distance
    end
    if coord.pmra !== NaN && coord.pmdec !== NaN
        @printf ", "
        @printf io "μ = (%4.2f, %4.2f) mas/yr" coord.pmra coord.pmdec
    end

    if coord.radial_velocity !== NaN
        @printf ", "
        @printf io "v_r = %4.2f km/s" coord.radial_velocity
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
