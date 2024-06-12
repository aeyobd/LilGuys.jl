"""A series of methods describing common mass profiles"""


import QuadGK: quadgk 
import Base: @kwdef
import SpecialFunctions: besselk
import ForwardDiff: derivative

abstract type AbstractProfile end

Base.Broadcast.broadcastable(p::AbstractProfile) = Ref(p)

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


"""
    calc_V_circ(profile, r)

Calculates the circular velocity at radius r
"""
function calc_V_circ end


"""
    ABC_Profile(α, β, γ, M, r_s)

The general ABC density profile (not implemented yet)
"""
@kwdef struct ABC_Profile <: AbstractProfile
    α::Float64
    β::Float64
    γ::Float64
    M::Float64 = 1
    r_s::Float64 = 1
end


@kwdef struct Plummer <: AbstractProfile
    a::Float64
    M::Float64 = 1
    r_s::Float64 = 1
end


@kwdef struct Sersic <: AbstractProfile
    n::Float64
    M::Float64 = 1
    r_s::Float64 = 1
end



"""
    Exp2D(M, R_s)

An exponential disk profile in 2D. The density profile is given by

Σ(r) = Σ_0 * exp(-r/R_s)

where M is the total mass and R_s is the scale radius, and Σ_0 = M / (2π * R_s^2)
"""
@kwdef struct Exp2D <: AbstractProfile
    M::Float64 = 1
    R_s::Float64 = 1
end


"""
    Exp3D(M, r_s)

An exponential disk profile in 3D. The density profile is given by

ρ(r) = ρ_0 * exp(-r/r_s)

where M is the total mass and r_s is the scale radius, and ρ_0 = M / (8π * r_s^3)
"""
@kwdef struct Exp3D <: AbstractProfile
    M::Float64 = 1
    r_s::Float64 = 1
end


@doc raw"""
    LogCusp2D(M, R_s)

A logarithmic cusp profile in 2D. The density profile is given by
"""
@kwdef struct LogCusp2D <: AbstractProfile
    M::Float64 = 1
    R_s::Float64 = 1
end


@doc raw"""
    NFW(M_s, r_s)

A Navarro-Frenk-White profile. The density profile is given by

```math
ρ(r) = \frac{M_s}{4π r_s^3} \frac{1}{x (1 + x)^2}
```

where M_s is the total mass, and r_s is the scale radius.
"""
@kwdef struct NFW <: AbstractProfile
    M_s::Float64 = 1
    r_s::Float64 = 1
end


@doc raw"""
    KingProfile(M, R_s, R_t)

A King profile. The density profile is given by

```math
ρ(R) = ρ_0 * [ (1+(R/R_s)^2)^(-1/2) - (1+(R_t/R_s)^2)^(-1/2) ]^2
```

where M is the (approximate) total mass, `R_s` is the core radius, and `R_t` is
the tidal radius. 

"""
@kwdef struct KingProfile <: AbstractProfile
    M::Float64 
    R_s::Float64 
    R_t::Float64 
end


function get_Σ_s(profile::Exp2D)
    return profile.M / (2π * profile.R_s^2)
end

function get_R_h(profile::Exp2D)
    return 1.67835 * profile.R_s
end

function get_r_h(profile::Exp2D)
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

function M(profile::Exp3D, r::Real)
    M, r_s = profile.M, profile.r_s
    x = r / r_s
    return -4π * M *(x^2 + 2*x + 2)*exp(-x)
end


function calc_ρ(profile::LogCusp2D, r::Real)
    M, r_s = profile.M, profile.r_s
    x = r / r_s
    ρ_s = get_ρ_s(profile)
    return ρ_s * exp(-x) / x
end


function get_ρ_s(profile::LogCusp2D)
    M, r_s = profile.M, profile.r_s
    return M / (4π * r_s^2)
end

function calc_Σ(profile::LogCusp2D, r::Real)
    M, r_s = profile.M, profile.r_s
    ρ_s = get_ρ_s(profile)
    Σ_s = ρ_s * 2 * r_s
    x = r / r_s
    return Σ_s * besselk(0, x)
end


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


function A_NFW(c::Real)
    return log(1 + c) - c / (1 + c)
end

function calc_Φ(profile::NFW, r::Real)
    x = r / profile.r_s
    return -G * profile.M_s / profile.r_s * log(1 + x) / x
end


function calc_ρ(profile::ABC_Profile, r::Real)
    α, β, γ = profile.α, profile.β, profile.γ
    return r^(-γ) * (1 + r^a) ^ ((γ - β)/α)
end


function calc_Σ(profile::KingProfile, r::Real)
    if r > profile.R_t
        return 0
    end
    M, r_s, r_t = profile.M, profile.R_s, profile.R_t
    Σ_s = M / (2π * r_s^2)
    x = r / r_s
    x_t = r_t / r_s
    return Σ_s * (
        (1 + x^2)^(-1/2) 
        - (1 + x_t^2)^(-1/2)
       )^2
end


function calc_ρ(profile::KingProfile, r::Real)
    if r > profile.R_t
        return 0
    end
    M, Rs, Rt = profile.M, profile.R_s, profile.R_t
    Σ_s = M / (2π * Rs^2)

    ρ = -1/2 * (
            π*sqrt(Rs^2 + r^2)*Rs^4*sqrt((Rs^2 + Rt^2)/Rs^2) 
            - 4*Rs^5
            - 4*Rs^3*r^2
       )/(Rs^2 * (
            Rs^4*sqrt((Rs^2 + Rt^2)/Rs^2) 
            + 2*Rs^2*r^2*sqrt((Rs^2 + Rt^2)/Rs^2)
            + r^4*sqrt((Rs^2 + Rt^2)/Rs^2)
           ))

    return (-1/π * Σ_s)* ρ
end


function calc_M(profile::AbstractProfile, r::Real)
    return quadgk(r -> 4π * r^2 * calc_ρ(profile, r), 0, r)[1]
end


function calc_M_2D(profile::AbstractProfile, R::Real)
    return quadgk(R -> 2π * R * calc_Σ(profile, R), 0, R)[1]
end


function calc_Σ_from_ρ(profile::AbstractProfile, R::Real)
    integrand(r) = calc_ρ(profile, r) * r / sqrt(r^2 - R^2)
    return 2*quadgk(integrand, R, Inf)[1]
end


function calc_V_circ(profile::AbstractProfile, r)
    return sqrt(G * calc_M(profile, r) / r)
end


"""
calculate the 3D density profile of a profile at radius r
"""
function calc_ρ_from_Σ(profile::AbstractProfile, r::Real)
    integrand(R) = derivative(R->calc_Σ(profile, R), R) / sqrt(R^2 - r^2)
    return -1/π * quadgk(integrand, r, Inf)[1]
end

