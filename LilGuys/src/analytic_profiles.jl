"""A series of methods describing common mass profiles"""


import QuadGK: quadgk 
import Base: @kwdef
import SpecialFunctions: besselk
import ForwardDiff: derivative
import TOML
import Roots: find_zero


abstract type AbstractProfile end
abstract type SphericalProfile <: AbstractProfile end

Base.Broadcast.broadcastable(p::AbstractProfile) = Ref(p)


const RECOGNIZED_PROFILES = (
    :Plummer,
    :Exp2D,
    :Exp3D,
    :LogCusp2D,
    :KingProfile,
    :NFW,
    :TruncNFW,
)


"""
    load_profile(params)

Loads in the instance of a Profile from the parameter file or dictionary. 

The file/dictionary may be structured as
```toml
[<name>]
# r_s = 1
# other kwargs
```
or 

```toml
[profile]
[profile.<name>]
# r_s = 1
# other kwargs
```

Throws an error if more than one top-level key is present in the profiles, or
if the profile is not in `RECOGNIZED_PROFILES`.
"""
function load_profile(params::Dict{String, <:Any})
    if "profile" ∈ keys(params)
        return load_profile(params["profile"])
    end

    if length(keys(params)) != 1
        throw(ArgumentError("Only one profile can be loaded at a time"))
    end

    profile_class = first(keys(params))
    kwargs = params[profile_class]

    profile_class = Symbol(profile_class)

    if profile_class ∉ RECOGNIZED_PROFILES
        throw(KeyError("Profile $profile_class not recognized. Recognized profiles are $RECOGNIZED_PROFILES"))
    end

    profile_class = getproperty(LilGuys, profile_class)

    profile = profile_class(;dict_to_tuple(kwargs)...)

	return profile
end


function load_profile(filename::String)
    params = TOML.parsefile(filename)
    return load_profile(params)
end



"""
    calc_ρ(profile, r)

Calculates the 3D density of the profile at radius r
"""
function calc_ρ end


"""
    calc_Σ(profile, R)

Calculates the projected 2D surface density of the profile at radius R
"""
function calc_Σ end


"""
    calc_M(profile, r)

Calculates the mass enclosed within 3D radius r
"""
function calc_M end



@kwdef struct Plummer <: SphericalProfile
    M::Float64 = 1
    r_s::Float64 = 1
end



"""
    Exp2D(M, R_s)

An exponential disk profile in 2D. The density profile is given by

Σ(r) = Σ_0 * exp(-r/R_s)

where M is the total mass and R_s is the scale radius, and Σ_0 = M / (2π * R_s^2)
"""
@kwdef struct Exp2D <: SphericalProfile
    M::Float64 = 1
    R_s::Float64 = 1
end


"""
    Exp3D(M, r_s)

An exponential disk profile in 3D. The density profile is given by

ρ(r) = ρ_0 * exp(-r/r_s)

where M is the total mass and r_s is the scale radius, and ρ_0 = M / (8π * r_s^3)
"""
@kwdef struct Exp3D <: SphericalProfile
    M::Float64 = 1
    r_s::Float64 = 1
end


@doc raw"""
    LogCusp2D(M, R_s)

A logarithmic cusp profile in 2D. The density profile is given by
```math
\rho(r) = 
```
"""
@kwdef struct LogCusp2D <: SphericalProfile
    M::Float64 = 1
    R_s::Float64 = 1
end



@doc raw"""
    KingProfile(; R_t[, R_s, M, k])

A King profile. The density profile is given by

```math
\Sigma(R) = k [ (1+(R/R_s)^2)^(-1/2) - (1+(R_t/R_s)^2)^(-1/2) ]^2
```

where `k` is the density scaling, `R_s` is the core radius, and `R_t` is
the tidal radius. 

# Parameters
- `R_t::Real`: the tidal radius
- `R_s::Real=1`: the core radius
- `M::Real`: the total mass. One of `M` or `k` must be specified
- `k::Real`: the density scaling. One of `M` or `k` must be specified

"""
struct KingProfile <: SphericalProfile
    k::Float64
    R_s::Float64
    R_t::Float64 
end

function KingProfile(;kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)

    if :R_s ∉ keys(kwargs)
        kwargs[:R_s] = 1
    end

    if :R_t ∉ keys(kwargs)
        throw(ArgumentError("R_t must be specified"))
    end

    if :M ∈ keys(kwargs)
        k = _get_king_k(kwargs[:M], kwargs[:R_s], kwargs[:R_t])
        pop!(kwargs, :M)
        kwargs[:k] = k
    elseif :k ∈ keys(kwargs)
        k = kwargs[:k]
    else
        k = 1
    end

    return KingProfile(kwargs[:k], kwargs[:R_s], kwargs[:R_t])
end


function get_Σ_s(profile::Exp2D)
    return profile.M / (2π * profile.R_s^2)
end

function calc_R_h(profile::Exp2D)
    return 1.67835 * profile.R_s
end

function calc_r_h(profile::Exp2D)
    r_h_over_r_s = 2.2235172865036716
    return r_h_over_r_s * profile.R_s
end


function calc_ρ(profile::Exp2D, r::Real)
    Σ_s = get_Σ_s(profile)
    ρ_s = Σ_s / (π * profile.R_s)

    x = r / profile.R_s
    return ρ_s * besselk(0, x)
end


function calc_Σ(profile::Exp2D, R::Real)
    Σ_s = get_Σ_s(profile)
    x = R / profile.R_s
    return Σ_s * exp(-x)
end

# function calc_M(profile::Exp2D, r::Float64)
#     M, R_s = profile.M, profile.R_s
#     return 2/3 π x (3 π L_1(x) K_2(x) + (3 π L_2(x) - 4 x) K_1(x)) + constant
#     L_n is the modified struve function
#     K_n is the modified second kind of bessel function
# end



function get_ρ_s(profile::Exp3D)
    return profile.M / (8π * profile.r_s^3)
end


function calc_ρ(profile::Exp3D, r::Real)
    M, r_s = profile.M, profile.r_s
    ρ_s = get_ρ_s(profile)
    return ρ_s * exp(-r/r_s)
end


function calc_Σ(profile::Exp3D, r::Real)
    M, r_s = profile.M, profile.r_s
    ρ_s = get_ρ_s(profile)
    Σ0 = 2ρ_s * r_s
    x = r/r_s
    return Σ0 * x * besselk(1, x)
end

function calc_M(profile::Exp3D, r::Real)
    M, r_s = profile.M, profile.r_s
    x = r / r_s
    return M * 1/2 * (2 - (x^2 + 2*x + 2)*exp(-x))
end


function calc_ρ(profile::LogCusp2D, r::Real)
    M, r_s = profile.M, profile.R_s
    x = r / r_s
    ρ_s = get_ρ_s(profile)
    return ρ_s * exp(-x) / x
end


function get_ρ_s(profile::LogCusp2D)
    M, r_s = profile.M, profile.R_s
    return M / (4π * r_s^3)
end

function calc_Σ(profile::LogCusp2D, r::Real)
    M, r_s = profile.M, profile.R_s
    ρ_s = get_ρ_s(profile)
    Σ_s = ρ_s * 2 * r_s
    x = r / r_s
    return Σ_s * besselk(0, x)
end

function calc_M(profile::LogCusp2D, r::Real)
    M, r_s = profile.M, profile.R_s
    x = r / r_s
    return M * (1 - exp(-x) - x*exp(-x))
end


function calc_Σ(profile::KingProfile, r::Real)
    if r > profile.R_t
        return 0
    end
    r_s, r_t = profile.R_s, profile.R_t
    k = profile.k
    x = r / r_s
    x_t = r_t / r_s
    return k * (
        (1 + x^2)^(-1/2) 
        - (1 + x_t^2)^(-1/2)
       )^2
end



function calc_ρ(profile::KingProfile, r::Real)
    if r > profile.R_t
        return 0
    end
    if r < 0
        return NaN
    end

    r_s, r_t = profile.R_s, profile.R_t
    k = profile.k
    # equation 27 in King 1962
    x = sqrt( (1+(r/r_s)^2)/(1+(r_t/r_s)^2) )
    K = k / (π*r_s * (1 + (r_t/r_s)^2)^(3/2))
    return K / x^2 * ( acos(x)/x - √(1-x^2) )
end


"""
Given the total mass, scale radius, and tidal radius of a king profile, 
returns the scaling constant k.
"""
function _get_king_k(M::Real, r_s::Real, r_t::Real)
    x_t = r_t / r_s

    M1 = π*r_s^2 * (
            log(1+x_t^2) - 4*( √(1+x_t^2) - 1)/√(1+x_t^2) + x_t^2/(1+x_t^2)
       )

    return M / M1
end


function calc_M_2D(profile::KingProfile, R::Real)
    r_s, r_t = profile.R_s, profile.R_t

    k = profile.k
    x = R / r_s
    x_t = r_t / r_s

    if x > x_t
        x = x_t
    end

    return π*r_s^2*k * (
        log(1+x^2) - 4*( √(1+x^2) - 1)/√(1+x_t^2) + x^2/(1+x_t^2)
       )
end


function calc_M(profile::SphericalProfile, r::Real)
    if r < 0
        return NaN
    end

    return calc_M_from_ρ(profile, r)
end


function calc_M_2D(profile::SphericalProfile, R::Real)
    return calc_M_2D_from_Σ(profile, R)
end

function calc_M_2D_from_Σ(profile::SphericalProfile, R::Real)
    return quadgk(R -> 2π * R * calc_Σ(profile, R), 0, R)[1]
end

function calc_M_from_ρ(profile::SphericalProfile, r::Real)
    return quadgk(r -> 4π * r^2 * calc_ρ(profile, r), 0, r)[1]
end


function calc_Σ_from_ρ(profile::SphericalProfile, R::Real)
    integrand(r) = calc_ρ(profile, r) * r / sqrt(r^2 - R^2)
    return 2*quadgk(integrand, R, Inf)[1]
end


"""
    calc_v_circ(profile, r)

Calculates the circular velocity at radius r
"""
function calc_v_circ(profile::SphericalProfile, r)
    return sqrt(G * calc_M(profile, r) / r)
end


"""
calculate the 3D density profile of a profile at radius r
"""
function calc_ρ_from_Σ(profile::SphericalProfile, r::Real)
    integrand(R) = derivative(R->calc_Σ(profile, R), R) / sqrt(R^2 - r^2)
    return -1/π * quadgk(integrand, r, Inf)[1]
end

"""
Approximate scale radius
"""
function get_r_s(profile::KingProfile)
    return profile.R_s / sqrt(2)
end

function calc_R_h(profile::KingProfile)
    r_s = get_r_s(profile)
    M0 = get_M_tot(profile)
    return find_zero(R -> calc_M_2D(profile, R) / M0 - 1/2, r_s)
end


function calc_r_h(profile::KingProfile)
    M0 = get_M_tot(profile)
    r_s = get_r_s(profile)
    return find_zero(r -> calc_M(profile, r) / M0 - 1/2, r_s)
end


function calc_r_h(profile::SphericalProfile)
    M0 = get_M_tot(profile)
    r_s = 1.0
    return find_zero(r -> calc_M(profile, r) / M0 - 1/2, r_s)
end

function calc_R_h(profile::SphericalProfile)
    r_s = 1.0 # initial guess
    M0 = get_M_tot(profile)
    return find_zero(R -> calc_M_2D(profile, R) / M0 - 1/2, r_s)
end


function get_M_tot(profile::SphericalProfile)
    if :M in fieldnames(typeof(profile))
        return profile.M
    else
        return calc_M(profile, Inf)
    end
end


function calc_Φ_from_ρ(profile::SphericalProfile, r::Real; integrate_M=false)
    integrand(r) = calc_ρ(profile, r) * r

    if integrate_M
        M_in = calc_M_from_ρ(profile, r)
    else
        M_in = calc_M(profile, r)
    end

    Φ_in = -G * M_in / r

    Φ_out = -4π * G * quadgk(integrand, r, Inf)[1]

    return Φ_in + Φ_out
end

