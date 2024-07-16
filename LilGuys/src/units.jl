const G = 1.
const M2MSUN = 1e10 #Msun
const R2KPC = 1. # kpc
const T2GYR = 4.715e-3 # Gyr pm 0.02 Myr
const V2KMS = 207.4 # km/s (pm 1 fron G uncertainty)

const kms_per_kpc_mas_per_yr = 4.740470463533348 

F = Float64
OptVector = Union{Vector{F}, Nothing}
OptMatrix = Union{Matrix{F}, Nothing}


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

function pm_to_kms(pm::Real, distance::Real)
    return pm * distance * kms_per_kpc_mas_per_yr
end
