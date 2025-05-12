abstract type GeneralNFW <: SphericalProfile end
import SpecialFunctions: expint



"""
The virial radius, i.e. the radius where the mean inner density is 200 times 
"""
function solve_R200(profile::GeneralNFW; δ=200, tol=1e-3)
    f(r) = mean_density(profile, r) - δ*ρ_crit

    R200 = find_zero(f, [0.1profile.r_s, 1000profile.r_s])

    if abs(f(R200)) > tol
        error("failed to solve for c")
    end
    
    return R200
end

function R200(profile::GeneralNFW)
    return solve_R200(profile)
end

function M200(profile::GeneralNFW)
    return mass(profile, R200(profile))
end


"""
NFW concentration parameter = r_s / r_200
"""
function concentration(profile::GeneralNFW)
    r200 = solve_R200(profile)
    c = r200 / profile.r_s

    return c
end


function concentration(M200::Real, r_s::Real)
    return R200(M200) / r_s
end

function solve_r_circ_max(profile::GeneralNFW)
    r_max = find_zero(r -> 4π * r * density(profile, r) - mass(profile, r)/r^2, [0.0001profile.r_s, 1000profile.r_s])
    return r_max
end

function r_circ_max(profile::GeneralNFW)
    return solve_r_circ_max(profile)
end

function v_circ_max(profile::GeneralNFW)
    r_max = r_circ_max(profile)
    return v_circ(profile, r_max)
end

@doc raw"""
    NFW(M_s, r_s)

A Navarro-Frenk-White profile. The density profile is given by

```math
ρ(r) = \frac{M_s}{4π r_s^3} \frac{1}{x (1 + x)^2}
```

where M_s is the total mass, and r_s is the scale radius.
The profile may be specified in terms of 
- M_s, r_s
- M200, c
- M200, r_s
- v_circ_max, r_circ_max
where v_circ_max is the maximum circular velocity, and r_circ_max is the radius at which it occurs; M200 is the mass within the virial radius, and c is the concentration parameter. 
"""
struct NFW <: GeneralNFW
    M_s::Float64
    r_s::Float64
    c::Union{Nothing, Float64}
end

function NFW(; c=nothing, kwargs...)
    arg_names = Set(keys(kwargs))
    valid_kwargs = [:M_s, :r_s, :M200, :v_circ_max, :r_circ_max]

    for arg in arg_names
        if !(arg ∈ valid_kwargs)
            throw(ArgumentError("Invalid keyword argument, $arg"))
        end
    end

    if arg_names == Set([:M_s, :r_s])
        M_s, r_s = kwargs[:M_s], kwargs[:r_s]
    elseif (arg_names == Set([:M200]) ) && (c !== nothing)
        M200 = kwargs[:M200]
        M_s, r_s = _NFW_from_M200_c(M200, c)
    elseif arg_names == Set([:M200, :r_s])
        M200, r_s = kwargs[:M200], kwargs[:r_s]
        M_s, r_s = _NFW_from_M200_r_s(M200, r_s; c=c)
        c = R200(M200) / r_s
    elseif arg_names == Set([:v_circ_max, :r_circ_max])
        v_circ_max, r_circ_max = kwargs[:v_circ_max], kwargs[:r_circ_max]
        M_s, r_s = _NFW_from_v_circ_max_r_circ_max(v_circ_max, r_circ_max)
    else
        throw(ArgumentError("Invalid keyword argument combination: $arg_names"))
    end

    if c === nothing
        c = concentration(NFW(M_s, r_s, nothing))
    end

    return NFW(M_s, r_s, c)
end


function _NFW_from_M200_c(M200, c)
    Rvir = R200(M200)
    r_s = Rvir / c
    M_s = M200 / A_NFW(c)
    return M_s, r_s
end

function _NFW_from_v_circ_max_r_circ_max(v_circ_max, r_circ_max)
    M_s = v_circ_max^2 * r_circ_max / G / A_NFW(α_nfw)
    r_s = r_circ_max / α_nfw
    return M_s, r_s
end


function _NFW_from_M200_r_s(M200, r_s; c=nothing)
    if c === nothing
        c = R200(M200) / r_s
    end
    M_s = M200 / A_NFW(c)
    return M_s, r_s
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
    V_s = 4π/3 * r_s^3 
    return M_s / V_s
end

function density(profile::NFW, r::Real)
    x = r / profile.r_s
    ρ_s = get_ρ_s(profile)
    return (ρ_s / 3) / (x * (1 + x)^2)
end


function mass(profile::NFW, r::Real)
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
    if c < 0
        return NaN
    elseif c === Inf
        return Inf
    end
    return log1p(c) - c / (1 + c)
end


function potential(profile::NFW, r::Real)
    Φ_0 = -G * profile.M_s / profile.r_s
    x = r / profile.r_s

    if (r < 0 )
        throw(DomainError(r, "r must be positive"))
    elseif x == 0 
        return Φ_0 # limiting case
    elseif x === Inf
        return 0.
    end

    return Φ_0 * log1p(x) / x
end


function v_circ_max(profile::NFW)
    r_max = r_circ_max(profile)
    M_max = mass(profile, r_max)

    return sqrt(G * M_max / r_max)
end


function r_circ_max(profile::NFW)
    return α_nfw * profile.r_s
end




"""
The virial radius, i.e. the radius where the mean inner density is 200 times 
"""
function R200(profile::NFW)
    return profile.r_s * profile.c
end


function M200(profile::NFW)
    return A_NFW(profile.c) * profile.M_s
end


function R200(M200::Real)
    return (3 * M200 / (4π * 200 * ρ_crit))^(1/3)
end


@doc raw"""
A truncated NFW profile

The density profile is given by

```math
\rho(r) = \rho_{NFW}(r) \exp(-r/r_t)
```
"""
struct TruncNFW <: GeneralNFW
    M_s::Float64
    r_s::Float64
    r_t::Float64
    c::Union{Nothing, Float64}
end 


function TruncNFW(; r_t=nothing, trunc=nothing, kwargs...)
    nfw = NFW(; kwargs...)

    if trunc !== nothing && r_t !== nothing
        throw(ArgumentError("Cannot specify both trunc and r_t"))
    elseif trunc !== nothing
        r_t = trunc * nfw.r_s
    elseif r_t !== nothing
        # pass
    else
        throw(ArgumentError("Must specify either trunc or r_t"))
    end

    return TruncNFW(nfw.M_s, nfw.r_s, r_t, nfw.c)
end


function density(profile::TruncNFW, r::Real)
    nfw = NFW(profile.M_s, profile.r_s, profile.c)
    return density(nfw, r) * exp(-(r/profile.r_t))
end

function mass(profile::TruncNFW, r::Real)
    x = r / profile.r_s
    t = profile.r_t / profile.r_s
    B = (1+1/t)*exp(1/t)
    return profile.M_s * ( B * expinti(-(x+1)/t) + exp(-x/t)/(1+x) - B*expinti(-1/t) - 1 )
end


function mass(profile::TruncNFW)
    t = profile.r_t / profile.r_s
    A = -(1+1/t)*exp(1/t)*expinti(-1/t) - 1

    return profile.M_s * A
end


function potential(profile::TruncNFW, r::Real)
    Φ_in = - G * mass(profile, r) / r

    Φ0 = -G * profile.M_s / profile.r_s
    x = r / profile.r_s
    t = profile.r_t / profile.r_s
    Φ_out = Φ0 * (1/(1+x)*exp(-x/t) + expinti(-(1+x)/t)*exp(1/t)/t)

    return Φ_in + Φ_out
end



@doc raw"""
    CoredNFW(; r_c, kwargs...)

@peñarrubia+2012 implementaation of a cored NFW profile. The density profile is given by

```math
\rho(r) = \rho_s (r_c / r_s + r / r_s)^{-1} (1 + r / r_s)^{-2}
```

where r_c is the core radius and r_s is the scale radius. 
"""
struct CoredNFW <: GeneralNFW
    M_s::Float64
    r_s::Float64
    r_c::Float64
    r_t::Float64
    c::Union{Nothing, Float64}
end

function CoredNFW(; r_c, r_s, M_s, c=nothing, r_t=100r_s)
    if c === nothing
        c = concentration(CoredNFW(M_s, r_s, r_c, r_t, nothing))
    end

    return CoredNFW(M_s, r_s, r_c, r_t, c)
end

function get_ρ_s(profile::CoredNFW)
    M_s, r_s = profile.M_s, profile.r_s
    V_s = 4π/3 * r_s^3 
    return M_s / V_s
end


function density(profile::CoredNFW, r::Real)
    r_c, r_s = profile.r_c, profile.r_s
    ρ_s = get_ρ_s(profile)
    return ρ_s/3 * (r_c / r_s + r / r_s)^(-1) * (1 + r / r_s)^(-2) * exp(-r/profile.r_t)
end

function mass(profile::CoredNFW, r::Real)
    # result from sagemath, maybe I will do this integral one day
    c = profile.r_c
    s = profile.r_s
    t = profile.r_t
    ρ_s = get_ρ_s(profile)
    Ei(x) = -expint(-x)
    M(r) =  4π * (ρ_s/3) * (
             (c^2*r*s^3 + c^2*s^4)*t*Ei(-(c + r)/t)*exp(c/t)
             - (c*s^5 - s^6)*t * exp(-r/t)
             - (
                (2*c*r*s^4 - s^6 + (2*c - r)*s^5)*t*Ei(-(r + s)/t) 
                + (c*r*s^5 - s^7 + (c - r)*s^6)*Ei(-(r + s)/t)
                )*exp(s/t)
         )/(
                     (c^2*r - (2*c - r)*s^2 + s^3 + (c^2 - 2*c*r)*s)*t
                )

    return M(r) - M(0)
end


function mass(profile::CoredNFW)
    return mass(profile, 300*profile.r_t)
end

