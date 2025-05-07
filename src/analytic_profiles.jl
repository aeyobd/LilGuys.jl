"""A series of methods describing common mass profiles"""


import Base: @kwdef
import SpecialFunctions: besselk, gamma, gamma_inc
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
    :ExpCusp,
    :KingProfile,
    :NFW,
    :TruncNFW,
    :CoredNFW,
    :Sersic,
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
    if !isfile(filename)
        throw(ArgumentError("File $filename does not exist"))
    end

    params = TOML.parsefile(filename)
    return load_profile(params)
end



"""
    density(profile, r)

Calculates the 3D density of the profile at radius r
"""
function density end


"""
    surface_density(profile, R)

Calculates the projected 2D surface density of the profile at radius R
"""
function surface_density end


"""
    mass(profile, r)

Calculates the mass enclosed within 3D radius r
"""
function mass end



"""
    scale(profile, radius_scale[, mass_scale])

Scales the profile by a factor in radius (and optionally mass).
"""
function scale(profile::AbstractProfile, radius_scale::Real, mass_scale::Real=1)
    raise(NotImplementedError()) 
end



@doc raw"""
    Plummer(M, r_s)

A Plummer profile. The density profile is given by

``
\Sigma(r) = \frac{M}{4\pi r_s^2} \left(1 + \frac{r^2}{r_s^2}\right)^{-2}
``

where M is the total mass and r_s is the scale radius. 
"""
@kwdef struct Plummer{F<:Real} <: SphericalProfile
    M::F = 1.
    r_s::F = 1.
end




@doc raw"""
    Exp2D(M, R_s)

An exponential disk profile in 2D. The density profile is given by

``
\Sigma = \Sigma_0 \, \exp(-R/R_s)
``

where M is the total mass and R_s is the scale radius, and Σ_0 = M / (2π * R_s^2)
"""
@kwdef struct Exp2D{F<:Real} <: SphericalProfile
    M::F = 1
    R_s::F = 1
end

function Exp2D(M::Real, R_s::Real)
    return Exp2D(promote(M, R_s)...)
end


@doc raw"""
    Exp3D(M, r_s)

An exponential disk profile in 3D. The density profile is given by

``  
\rho = \rho_0 \, \exp(-r/r_s)
``

where M is the total mass and r_s is the scale radius, and density_0 = M / (8π * r_s^3)
"""
@kwdef struct Exp3D{F<:Real} <: SphericalProfile
    M::F = 1
    r_s::F = 1
end

function Exp3D(M::Real, r_s::Real)
    return Exp3D(promote(M, r_s)...)
end


@doc raw"""
    LogCusp2D(M, R_s)

A logarithmic cusp profile in 2D. The density profile is given by

``math
\Sigma(r) = \Sigma_0\, 
``
"""
@kwdef struct LogCusp2D{F<:Real} <: SphericalProfile
    M::F = 1
    R_s::F = 1
end

function LogCusp2D(M::Real, R_s::Real)
    return LogCusp2D(promote(M, R_s)...)
end


@doc raw"""
    KingProfile(; R_t[, R_s, M, k, c])

A King profile. The density profile is given by

```math
\Sigma(R) = k [ (1+(R/R_s)^2)^{-1/2} - (1+(R_t/R_s)^2)^{-1/2} ]^2
```

where `k` is the density scaling, `R_s` is the core radius, and `R_t` is
the tidal radius. 

# Parameters
- `R_t::Real`: the tidal radius
- `R_s::Real=1`: the core radius
- `c::Real`: the concentration (i.e. R_t / R_s). One of `c` or `R_t` must be specified
- `M::Real`: the total mass. One of `M` or `k` must be specified
- `k::Real`: the density scaling. One of `M` or `k` must be specified

"""
struct KingProfile{F<:Real} <: SphericalProfile
    k::F
    R_s::F
    R_t::F 
end

function KingProfile(k::Real, R_s::Real, R_t::Real)
    return KingProfile(promote(k, R_s, R_t)...)
end

@doc raw"""
    ExpCusp(M, r_s)

An exponential cusp profile in 3D. The density profile is given by

```math
\rho = \rho_0 \exp(-r/r_s) / (r/r_s)
```
where $M_{\rm tot} = 4\pi \rho_0 r_s^3$ is the total mass.
This is the asymptotic result of an NFW profile under heavy
tidal stripping.
"""
@kwdef struct ExpCusp{F<:Real} <: SphericalProfile
    M::F = 1
    r_s::F = 1
end

function ExpCusp(M::Real, r_s::Real)
    return ExpCusp(promote(M, r_s)...)
end

function KingProfile(;kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)

    if :R_s ∉ keys(kwargs)
        kwargs[:R_s] = 1
    end

    if :c ∈ keys(kwargs)
        kwargs[:R_t] = kwargs[:c] * kwargs[:R_s]
        pop!(kwargs, :c)
    end

    if :R_t ∉ keys(kwargs)
        throw(ArgumentError("R_t or c must be specified"))
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


# Calculations
#
#

function density(profile::Plummer, r::Real)
    M, r_s = profile.M, profile.r_s
    x = r / r_s
    ρ0 = get_ρ0(profile)
    return ρ0 * (1 + x^2)^(-5/2)
end

function get_ρ0(profile::Plummer)
    return 3*profile.M / (4π * profile.r_s^3)
end

function get_Σ0(profile::Plummer)
    return profile.M / (π * profile.r_s^2)
end

function surface_density(profile::Plummer, R::Real)
    M, r_s = profile.M, profile.r_s
    x = R / r_s
    Σ0 = get_Σ0(profile)
    return Σ0 * (1 + x^2)^(-2)
end
    

function mass(profile::Plummer, r::Real)
    M, r_s = profile.M, profile.r_s
    x = r / r_s
    return M * x^3 / (1 + x^2)^(3/2)
end


function mass_2D(profile::Plummer, R::Real)
    M, r_s = profile.M, profile.r_s
    x = R / r_s
    return M * x^2 / (1 + x^2)
end

function R_h(profile::Plummer)
    return profile.r_s
end

function r_h(profile::Plummer)
    return (2^(2/3) - 1)^(-1/2) * profile.r_s
end

function scale(profile::Plummer, radius_scale::Real, mass_scale::Real=1)
    return Plummer(profile.M * mass_scale, profile.r_s * radius_scale)
end


############ Exp2D 
function get_Σ_s(profile::Exp2D)
    return profile.M / (2π * profile.R_s^2)
end

function R_h(profile::Exp2D)
    α = 1.6783469900166605 # solution to 1/2 = 1 - (1 + x) exp(-x)
    return α * profile.R_s
end

function r_h(profile::Exp2D)
    r_h_over_r_s = 2.2235172865036716
    return r_h_over_r_s * profile.R_s
end


function density(profile::Exp2D, r::Real)
    Σ_s = get_Σ_s(profile)
    ρ_s = Σ_s / (π * profile.R_s)

    x = r / profile.R_s
    return ρ_s * besselk(0, x)
end


function surface_density(profile::Exp2D, R::Real)
    Σ_s = get_Σ_s(profile)
    x = R / profile.R_s
    return Σ_s * exp(-x)
end

function scale(profile::Exp2D, radius_scale::Real, mass_scale::Real)
    return Exp2D(profile.M * mass_scale, profile.R_s * radius_scale)
end

# function mass(profile::Exp2D, r::Float64)
#     M, R_s = profile.M, profile.R_s
#     return 2/3 π x (3 π L_1(x) K_2(x) + (3 π L_2(x) - 4 x) K_1(x)) + constant
#     L_n is the modified struve function
#     K_n is the modified second kind of bessel function
# end


############ Exp3D

function get_ρ_s(profile::Exp3D)
    return profile.M / (8π * profile.r_s^3)
end


function density(profile::Exp3D, r::Real)
    M, r_s = profile.M, profile.r_s
    ρ_s = get_ρ_s(profile)
    return ρ_s * exp(-r/r_s)
end


function surface_density(profile::Exp3D, r::Real)
    M, r_s = profile.M, profile.r_s
    ρ_s = get_ρ_s(profile)
    Σ0 = 2ρ_s * r_s
    x = r/r_s
    return Σ0 * x * besselk(1, x)
end

function mass(profile::Exp3D, r::Real)
    M, r_s = profile.M, profile.r_s
    x = r / r_s
    return M * 1/2 * (2 - (x^2 + 2*x + 2)*exp(-x))
end

function scale(profile::Exp3D, radius_scale::Real, mass_scale::Real)
    return Exp3D(profile.M * mass_scale, profile.r_s * radius_scale)
end


############ LogCusp2D

function density(profile::LogCusp2D, r::Real)
    M, r_s = profile.M, profile.R_s
    x = r / r_s
    ρ_s = get_ρ_s(profile)
    return ρ_s * exp(-x) / x
end


function get_ρ_s(profile::LogCusp2D)
    M, r_s = profile.M, profile.R_s
    return M / (4π * r_s^3)
end

function surface_density(profile::LogCusp2D, r::Real)
    M, r_s = profile.M, profile.R_s
    ρ_s = get_ρ_s(profile)
    Σ_s = ρ_s * 2 * r_s
    x = r / r_s
    return Σ_s * besselk(0, x)
end

function mass(profile::LogCusp2D, r::Real)
    M, r_s = profile.M, profile.R_s
    x = r / r_s
    return M * (1 - exp(-x) - x*exp(-x))
end

function scale(profile::LogCusp2D, radius_scale::Real, mass_scale::Real)
    return LogCusp2D(
        M = profile.M * mass_scale, 
        R_s = profile.R_s * radius_scale
    )
end


############ ExpCusp

function density(profile::ExpCusp, r::Real)
    r_s = profile.r_s
    x = r / r_s
    ρ_s = get_ρ_s(profile)
    return ρ_s * exp(-x) / x
end


function get_ρ_s(profile::ExpCusp)
    M, r_s = profile.M, profile.r_s
    return M / (4π * r_s^3)
end


function mass(profile::ExpCusp, r::Real)
    M, r_s = profile.M, profile.r_s
    x = r / r_s
    return M * (1 - (1 + x) * exp(-x))
end

function r_circ_max(profile::ExpCusp)
    α = 1.7932821329007609 # root of x^2 + x + 1 - exp(x)
    return α * profile.r_s
end

function v_circ_max(profile::ExpCusp)
    rm = r_circ_max(profile)
    return v_circ(profile, rm)
end

function r_h(profile::ExpCusp)
    α = 1.6783469900166605 # solution to 1/2 = 1 - (1 + x) exp(-x)
    return α * profile.r_s
end


function scale(profile::ExpCusp, radius_scale::Real, mass_scale::Real)
    return ExpCusp(
        M = profile.M * mass_scale, 
        r_s = profile.r_s * radius_scale
    )
end

############ KingProfile

function surface_density(profile::KingProfile, r::Real)
    if r > profile.R_t
        return 0.
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



function density(profile::KingProfile, r::Real)
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


function mass_2D(profile::KingProfile, R::Real)
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


function scale(profile::KingProfile, radius_scale::Real, mass_scale::Real)
    return KingProfile(
        k = profile.k * mass_scale / radius_scale^2, # k is a 2D density
        R_s = profile.R_s * radius_scale,
        R_t = profile.R_t * radius_scale
    )
end



function mass(profile::SphericalProfile, r::Real)
    if r < 0
        return NaN
    end

    return mass_from_density(profile, r)
end


function mass_2D(profile::SphericalProfile, R::Real)
    return mass_from_density_2D(profile, R)
end

function mass_from_density_2D(profile::SphericalProfile, R::Real)
    return integrate(R -> 2π * R * surface_density(profile, R), 0, R)
end

function mass_from_density(profile::SphericalProfile, r::Real)
    return integrate(r -> 4π * r^2 * density(profile, r), 0, r)
end


function surface_density_from_density(profile::SphericalProfile, R::Real)
    integrand(r) = density(profile, r) * r / sqrt(r^2 - R^2)
    return 2*integrate(integrand, R, Inf)
end

function surface_density(profile::SphericalProfile, R::Real)
    return surface_density_from_density(profile, R)
end



"""
    v_circ(profile, r)

Calculates the circular velocity at radius r
"""
function v_circ(profile::SphericalProfile, r)
    return sqrt(G * mass(profile, r) / r)
end


@doc raw"""
calculate the 3D density profile of a profile at radius r provided the 
profile has a 2D surface density profile implemented.
Abel inversion formula:

``
\rho(r) = -1/π \int_r^\infty dR dΣ/dR / \sqrt{R^2 - r^2}
``
"""
function density_from_surface_density(profile::SphericalProfile, r::Real)
    integrand(R) = derivative(R->surface_density(profile, R), R) / sqrt(R^2 - r^2)
    return -1/π * integrate(integrand, r, Inf)
end

"""
Approximate scale radius
"""
function get_r_s(profile::KingProfile)
    return profile.R_s / sqrt(2)
end


"""
Mean density instide radius
"""
function mean_density(profile::SphericalProfile, r::Real)
    return mass(profile, r) / (4π/3 * r^3)
end

function R_h(profile::KingProfile)
    r_s = get_r_s(profile)
    M0 = mass(profile)
    return find_zero(R -> mass_2D(profile, R) / M0 - 1/2, r_s)
end


function r_h(profile::KingProfile)
    M0 = mass(profile)
    r_s = get_r_s(profile)
    return find_zero(r -> mass(profile, r) / M0 - 1/2, r_s)
end




@doc raw"""
    Sersic(M, r_s, n)

A Sérsic profile. The density profile is given by

``math
\Sigma(r) = \Sigma_h \exp(-b_n ((r/r_h)^{1/n} - 1))
``

where $\Sigma_h$ is the surface density at the half-light radius $r_h$.
$b_n$ is a function of $n$, which is calculated as the solution to $\gamma(2n, b_n) = 1/2\Gamma(2n)$ where $\gamma$ is the lower incomplete gamma function and $\Gamma$ is the gamma function.
"""
struct Sersic{F<:Real} <: SphericalProfile
    Σ_h::F 
    R_h::F
    n::F
    _b_n::F
end

_sersic_allowed_kwargs = (:Σ_h, :R_h, :n, :_b_n)

function Sersic(Σ_h::Real, R_h::Real, n::Real, _b_n::Real)
    return Sersic(promote(Σ_h, R_h, n, _b_n)...)
end

function Sersic(;kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)

    if any(k -> k ∉ _sersic_allowed_kwargs, keys(kwargs))
        throw(ArgumentError("Only the following kwargs are allowed: $_sersic_allowed_kwargs"))
    end

    R_h = get(kwargs, :R_h, 1)
    n = get(kwargs, :n, 1)
    Σ_h = get(kwargs, :Σ_h, 1)
    if :_b_n in keys(kwargs)
        _b_n = kwargs[:_b_n]
    else
        _b_n = b_n(n)
    end

    return Sersic(Σ_h, R_h, n, _b_n)
end


function b_n(n::Real)
    return find_zero(b_n -> gamma_inc(2n, b_n)[1] - 1/2, guess_b_n(n))
end

@doc raw"""
    guess_b_n(n)

Returns a guess for the value of $b_n$ given $n$.
"""
function guess_b_n(n::Real)
    if n > 0.36
        return 2n - 1/3 + 4/(405n) + 46/(25515n^2) + 131/(1148175n^3) - 2194697/(30690717750n^4)
    else 
        return 0.5
    end
end


function surface_density(profile::Sersic, R::Real)
    Σ_h, R_h, n = profile.Σ_h, profile.R_h, profile.n
    return Σ_h * exp(-profile._b_n * ((R/R_h)^(1/n) - 1))
end




# General methods

function r_h(profile::SphericalProfile)
    M0 = mass(profile)
    r_s = 1.0
    return find_zero(r -> mass(profile, r) / M0 - 1/2, r_s)
end

function R_h(profile::SphericalProfile)
    r_s = 1.0 # initial guess
    M0 = mass(profile)
    return find_zero(R -> mass_2D(profile, R) / M0 - 1/2, r_s)
end


function mass(profile::SphericalProfile)
    if :M in fieldnames(typeof(profile))
        return profile.M
    else
        return mass(profile, Inf)
    end
end

function potential(profile::SphericalProfile, r::Real)
    return potential_from_density(profile, r)
end



function potential_from_density(profile::SphericalProfile, r::Real; integrate_M=false)
    integrand(r) = density(profile, r) * r

    if integrate_M
        M_in = mass_from_density(profile, r)
    else
        M_in = mass(profile, r)
    end

    Φ_in = -G * M_in / r

    Φ_out = -4π * G * integrate(integrand, r, Inf)

    return Φ_in + Φ_out
end



"""
    σv_star_mean(dm_profile, stellar_profile)

Given a dark matter (mass) profile and a stellar profile, assuming isotropic and spherical system,
returns the light-weighted line of sight (1D) velocity dispersion for the system (in code units).
"""
function σv_star_mean(dm_profile::SphericalProfile, stellar_profile::SphericalProfile; R_min=0, R_max=Inf)
    ρ(r) = density(stellar_profile, r)
    m(r) = mass(dm_profile, r)
    integrand(r) = ρ(r) * G * m(r) * 4π/3 * r

    weighted_σ2 = (
            integrate(integrand, R_min, R_max)
    )

    Mtot = mass(stellar_profile, R_max)

    sqrt(weighted_σ2 / Mtot)
end


