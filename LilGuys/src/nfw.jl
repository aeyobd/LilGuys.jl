import Roots: find_zero


@doc raw"""
    NFW(M_s, r_s)

A Navarro-Frenk-White profile. The density profile is given by

```math
ρ(r) = \frac{M_s}{4π r_s^3} \frac{1}{x (1 + x)^2}
```

where M_s is the total mass, and r_s is the scale radius.
The profile may be specified in terms of 
- M_s, r_s. The scale mass and scale radius
- M200, c. The 
- M200, r_s
- V_circ_max, r_circ_max
"""
struct NFW <: AbstractProfile
    M_s::Float64
    r_s::Float64
    c::Union{Nothing, Float64}
end

function NFW(; M_s=nothing, r_s=nothing, c=nothing, M200=nothing, R200=nothing,
    V_circ_max=nothing, r_circ_max=nothing)

    if (M200 !== nothing) 
        R200 = calc_R200(M200)
        if r_s !== nothing
            c = r_s / R200
        elseif c !== nothing
            r_s = R200 / c
        end
        M_s = M200 / A_NFW(c)
    end


    if (V_circ_max !== nothing) && (r_circ_max !== nothing)
        M_s = V_circ_max^2 * r_circ_max / G / A_NFW(α_nfw)
        r_s = r_circ_max / α_nfw
    end

    if (M_s === nothing) || (r_s === nothing)
        error("Either M_s and r_s must be given, or M200 and R200, or V_circ_max and r_circ_max")
    end


    if c == nothing
        c = calc_c(NFW(M_s, r_s, nothing))
    end
    return NFW(M_s, r_s, c)
end





@doc raw"""the ratio between scale radius and radius of maximum circular velocity for a NFW halo



Solution to 
```math
0 = -A(x)/x^2 + 1/(x+1)^2
```
"""
const α_nfw = 2.1625815870646098348565536696032645735

"""
renormalized hubble constant from Plank collaboration (2018)
"""
const h_hubble = 0.674
const ρ_crit = 277.5366*h_hubble^2 / 1e10 # code units, 10^10 M_sun / kpc^3


function get_ρ_s(profile::NFW)
    M_s, r_s = profile.M_s, profile.r_s
    return M_s / (4π * r_s^3)
end

function calc_ρ(profile::NFW, r::Real)
    x = r / profile.r_s
    ρ_s = get_ρ_s(profile)
    return (ρ_s / 3) / (x * (1 + x)^2)
end


function calc_M(profile::NFW, r::Real)
    x = r / profile.r_s
    return profile.M_s * A_NFW(x)
end


@doc raw"""
NFW dimensionless scale function. Useful for many other formulae.

```math
A_NFW(c) = \log(1 + c) - \frac{c}{1 + c}
```
"""
function A_NFW(c::Real)
    return log(1 + c) - c / (1 + c)
end


function calc_Φ(profile::NFW, r::Real)
    x = r / profile.r_s
    return -G * profile.M_s / profile.r_s * log(1 + x) / x
end


function calc_V_circ_max(profile::NFW)
    r_max = calc_r_circ_max(profile)
    M = calc_M(profile, r_max)

    return sqrt(G * M / r_max)
end


function calc_r_circ_max(profile::NFW)
    return α_nfw * profile.r_s
end



"""
NFW concentration parameter = r_s / r_200
"""
function calc_c(profile::NFW; tol=1e-3)
    f(c) = 200ρ_crit - lguys.calc_ρ_mean(profile, c * profile.r_s)
    c = find_zero(f, 10)

    if abs(f(c)) > tol
        error("failed to solve for c")
    end
    
    return c
end


function calc_c(M200::Real, R200::Real)
    return R200 / calc_R200(M200)
end


"""
Mean density instide radius
"""
function calc_ρ_mean(profile::NFW, r::Real)
    return calc_M(profile, r) / (4π/3 * r^3)
end


"""
The virial radius, i.e. the radius where the mean inner density is 200 times 
"""
function calc_R200(profile::NFW)
    return profile.r_s * profile.c
end


function calc_M200(profile::NFW)
    return A_NFW(profile.c) * profile.M_s
end


function calc_R200(M200::Real)
    return (3 * M200 / (4π * 200 * ρ_crit))^(1/3)
end


