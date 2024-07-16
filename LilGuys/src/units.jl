const G = 1.
const M2MSUN = 1e10 #Msun
const R2KPC = 1. # kpc
const T2GYR = 4.715e-3 # Gyr pm 0.02 Myr
const V2KMS = 207.4 # km/s (pm 1 fron G uncertainty)


F = Float64
OptVector = Union{Vector{F}, Nothing}
OptMatrix = Union{Matrix{F}, Nothing}


# source IAU
const SECONDS_PER_YEAR = 31_557_600 # seconds; exact IAU
const M_PER_PC = 3.085677581491367e+16 # meters; exact IAU
const M_PER_AU = 149_597_870_700 # meters; exact IAU


"""Number of arcmin in radians"""
const ARCMIN_PER_RAD = π / (60 * 180) # exact; mathematical


"""
Calculates the physical diameter given the angular diameter and distance.

TODO: could also use Unitful to be more general
"""
function arcmin_to_kpc(arcmin::Real, distance::Real)
    return arcmin * ARCMIN_IN_RAD * distance 
end


"""
Converts a physical length to a sky angular diameter in arcminutes
"""
function kpc_to_arcmin(length::Real, distance::Real)
    return length / distance / arcmin_IN_rad
end
