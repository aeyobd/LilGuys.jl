# TODO: could use Unitful here to be more general


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

const SECONDS_PER_YEAR = 31_557_600 # seconds; exact IAU
const M_PER_PC = 3.085677581491367e+16 # meters; exact IAU
const M_PER_AU = 149_597_870_700 # meters; exact IAU


"""Number of arcmin in radians"""
const ARCMIN_PER_RAD = (60 * 180) / π  # exact; mathematical


"""
    arcmin2kpc(arcmin::Real, distance::Real)

Convert an angular scale in arcmin to a physical length in kpc assuming a distance in kpc.
"""
function arcmin2kpc(arcmin::Real, distance::Real)
    return arcmin * distance  / ARCMIN_PER_RAD
end


"""
    kpc2arcmin(length::Real, distance::Real)

Convert a physical length in kpc to a sky angular diameter in arcminutes.
"""
function kpc2arcmin(length::Real, distance::Real)
    return length / distance * ARCMIN_PER_RAD
end


"""
    pm2kms(proper_motion::Real, distance::Real)

Convert a proper motion in mas/yr at a distance in kpc to the tangental velocity in km/s.
"""
function pm2kms(proper_motion::Real, distance::Real)
    return proper_motion * distance * kms_per_kpc_mas_per_yr
end


"""
    kms2pm(velocity::Real, distance::Real)

Convert a tangental velocity in km/s at a distance in kpc to a proper motion in mas/yr.
"""
function kms2pm(velocity::Real, distance::Real)
    return velocity / distance / kms_per_kpc_mas_per_yr
end


"""
    kpc2dm(dm::Real)

Convert a distance in kpc to a distance modulus.
"""
function kpc2dm(dist::Real)
    return 5 * (log10(dist)  - 1 + 3)
end


"""
    dm2kpc(dm::Real)

Convert a distance modulus to a distance in kpc
"""
function dm2kpc(dm::Real)
    return 10 .^ (dm / 5 + 1 - 3)
end
