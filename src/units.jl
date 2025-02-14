"""Internal float unit """
F = Float64

"""Optionally specified vector (union vector and nothing)"""
OptVector = Union{Vector{F}, Nothing}

"""Optionally specified matrix (union matrix{F} and nothing)"""
OptMatrix = Union{Matrix{F}, Nothing}


"""Gravitational constent in internal units"""
const G = 1.

"""Internal mass to solar masses"""
const M2MSUN = 1e10 #Msun

"""Internal length scale to kpc"""
const R2KPC = 1. # kpc

"""Internal time unit to Gyr"""
const T2GYR = 4.715e-3 # Gyr pm 0.02 Myr

"""Internal velocity unit to km/s"""
const V2KMS = 207.4 # km/s (pm 1 fron G uncertainty)


"""The tangental velocity in km/s of an object at 1 kpc moving with a proper motion of 1 mas / year. """
const kms_per_kpc_mas_per_yr = 4.740470463533348 


const kpc_per_Gyr_per_kms = 1.0227 # TODO double check this one
# source IAU
const SECONDS_PER_YEAR = 31_557_600 # seconds; exact IAU
const M_PER_PC = 3.085677581491367e+16 # meters; exact IAU
const M_PER_AU = 149_597_870_700 # meters; exact IAU


"""Number of arcmin in radians"""
const ARCMIN_PER_RAD = (60 * 180) / π  # exact; mathematical


"""
Calculates the physical diameter given the angular diameter and distance.

TODO: could also use Unitful to be more general
"""
function arcmin_to_kpc(arcmin::Real, distance::Real)
    return arcmin * distance  / ARCMIN_PER_RAD
end


"""
Converts a physical length to a sky angular diameter in arcminutes
"""
function kpc_to_arcmin(length::Real, distance::Real)
    return length / distance * ARCMIN_PER_RAD
end


"""
Converts a proper motion (mas/yr) at a given distance (kpc) to the tangental velocity (km/s)
"""
function pm_to_kms(pm::Real, distance::Real)
    return pm * distance * kms_per_kpc_mas_per_yr
end


"""
    dm_to_dist(dm::Real)

Converts a distance modulus to a distance in kpc
"""
function dm_to_dist(dm::Real)
    return 10 .^ (dm / 5 + 1 - 3)
end
