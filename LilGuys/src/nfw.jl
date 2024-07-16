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
- v_circ_max, r_circ_max
"""
struct NFW <: AbstractProfile
    M_s::Float64
    r_s::Float64
    c::Union{Nothing, Float64}
end

function NFW(; M_s=nothing, r_s=nothing, c=nothing, M200=nothing, R200=nothing,
    v_circ_max=nothing, r_circ_max=nothing)

    if (M200 !== nothing) 
        R200 = calc_R200(M200)
        if r_s !== nothing
            c = r_s / R200
        elseif c !== nothing
            r_s = R200 / c
        end
        M_s = M200 / A_NFW(c)
    end


    if (v_circ_max !== nothing) && (r_circ_max !== nothing)
        M_s = v_circ_max^2 * r_circ_max / G / A_NFW(α_nfw)
        r_s = r_circ_max / α_nfw
    end

    if (M_s === nothing) || (r_s === nothing)
        error("Either M_s and r_s must be given, or M200 and R200, or v_circ_max and r_circ_max")
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
    V_s = 4π/3 * r_s^3 
    return M_s / V_s
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


function calc_v_circ_max(profile::NFW)
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
    f(c) = 200ρ_crit - calc_ρ_mean(profile, c * profile.r_s)
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



module Ludlow 
    export solve_rmax, c_ludlow

    using ..LilGuys
    import Roots: find_zero

	h = 0.674 #pm 0.05
	Ω_m0 = 0.315
	Ω_Λ0 = 1 - Ω_m0
	σ8 = 0.811

	n_s = 0.965
	
	ρ_crit = 277.5366*h^2 # M⊙/kpc = (3H^2 / 8πG)

	const G = 4.30091e-6 # km^2/s^2 kpc / M⊙

	Ω_Λ(z) = Ω_Λ0 / (Ω_Λ0 + Ω_m0 * (1+z)^3)
	Ω_m(z) = 1 - Ω_Λ(z)

	Ψ(z) = Ω_m(z)^(4/7) - Ω_Λ(z) + (1 + Ω_m(z)/2) * (1 + Ω_Λ(z)/70)
	D(z) = Ω_m(z) / Ω_m0 * Ψ(0) / Ψ(z) * (1+z)^-1

	ξ(M) = 1/(h * M/1e10) # with M in M⊙
	σ(M, z) = D(z) * 22.26*ξ(M)^0.292 / (1 + 1.53*ξ(M)^0.275 + 3.36*ξ(M)^0.198)

	c_0(z) = 3.395 * (1+z)^-0.215
	β(z) = 0.307 * (1+z)^0.540
	γ_1(z) = 0.628  * (1+z)^-0.047
	γ_2(z) = 0.317 * (1+z)^-0.893
	a(z) = 1/(1+z)
	ν_0(z) = 1/D(z) * (4.135 - 0.564/a(z) - 0.210/a(z)^2 + 0.0557/a(z)^3 - 0.00348/a(z)^4)
	δ_sc = 1.686
	ν(M, z) = δ_sc / σ(M, z)


    """
    Calculates the approximate concentration of a halo given the initial mass and redshift using the n-body fit from @ludlow2016
    """
    function c_ludlow(M, z)
        x = ν(M * M2MSUN, z)/ν_0(z)
        result = c_0(z)
        result *= x^(-γ_1(z))
        result *= (1 + x^(1/β(z)))^(-β(z) * (γ_2(z) - γ_1(z)))
        return result
    end

    """
        solve_M200_c(Vcmax, δlogc=0; interval=[0.001, 1000])

    Solves for the mass and concentration of a halo given the maximum circular velocity. Calls Roots.find_zero to solve for the mass on the given interval and can apply a multiplicative factor `δlogc` to the concentration parameter.

    See also [`c_ludlow`](@ref), [`solve_rmax`](@ref)
    """
    function solve_M200_c(Vcmax, δlogc=0; interval=[0.001, 1000])
        dc = 10 ^ (0 + δlogc)

        f(M200) = LilGuys.calc_v_circ_max(NFW(M200=M200, c=dc * c_ludlow(M200, 0.))) - Vcmax

        M200 = find_zero(f, interval)
        return M200, c_ludlow(M200, 0.) * dc
    end

    """
    solve_rm(Vcmax, δlogc=0; kwargs...)

    Solves for the radius of maximum circular velocity given the maximum circular velocity
    """
    function solve_rmax(Vcmax, δlogc=0; kwargs...)
        M200, c = solve_M200_c(Vcmax, δlogc; kwargs...)
        return calc_r_circ_max(NFW(M200=M200, c=c))
    end

end
