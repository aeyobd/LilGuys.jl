"""A series of methods describing common mass profiles"""


import QuadGK: quadgk 
import Base: @kwdef
import SpecialFunctions: besselk

abstract type AbstractProfile end

Base.Broadcast.broadcastable(p::AbstractProfile) = Ref(p)


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


@kwdef struct LogCusp2D <: AbstractProfile
    M::Float64 = 1
    r_s::Float64 = 1
end


@kwdef struct NFW <: AbstractProfile
    M_s::Float64 = 1
    r_s::Float64 = 1
end


"""
    KingProfile(M, r_s, r_t)

A King profile. The density profile is given by

ρ(r) = ρ_0 * [ (1+(r/r_s)^2)^(-1/2) - (1+(r_t/r_s)^2)^(-1/2) ]^2

where M is the (approximate) total mass, r_s is the core radius, and r_t is the tidal radius. 

"""
@kwdef struct KingProfile <: AbstractProfile
    M::Float64 
    r_s::Float64 
    r_t::Float64 
end


function get_Σ_s(profile::Exp2D)
    return profile.M / (2π * profile.R_s^2)
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
    if r > profile.r_t
        return 0
    end
    M, r_s, r_t = profile.M, profile.r_s, profile.r_t
    Σ_s = M / (2π * r_s^2)
    x = r / r_s
    x_t = r_t / r_s
    return Σ_s * (
        (1 + x^2)^(-1/2) 
        - (1 + x_t^2)^(-1/2)
       )^2
end



function calc_M(profile::AbstractProfile, r::Real)
    return quadgk(r -> 4π * r^2 * calc_ρ(profile, r), 0, r)[1]
end

function calc_Σ(profile::AbstractProfile, r::Real)
    return quadgk(r -> 2π * r * calc_ρ(profile, r), 0, r)[1]
end

function calc_V_circ(profile::AbstractProfile, r)
    return sqrt(G * calc_M(profile, r) / r)
end
