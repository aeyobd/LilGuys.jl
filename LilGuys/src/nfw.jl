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
struct NFW <: SphericalProfile
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
        c = calc_R200(M200) / r_s
    elseif arg_names == Set([:v_circ_max, :r_circ_max])
        v_circ_max, r_circ_max = kwargs[:v_circ_max], kwargs[:r_circ_max]
        M_s, r_s = _NFW_from_v_circ_max_r_circ_max(v_circ_max, r_circ_max)
    else
        throw(ArgumentError("Invalid keyword argument combination: $arg_names"))
    end

    if c === nothing
        c = calc_c(NFW(M_s, r_s, nothing))
    end

    return NFW(M_s, r_s, c)
end


function _NFW_from_M200_c(M200, c)
    R200 = calc_R200(M200)
    r_s = R200 / c
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
        c = calc_R200(M200) / r_s
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
    if c < 0
        return NaN
    elseif c === Inf
        return Inf
    end
    return log1p(c) - c / (1 + c)
end


function calc_Φ(profile::NFW, r::Real)
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
    c = find_zero(f, [0.1, 1000])

    if abs(f(c)) > tol
        error("failed to solve for c")
    end
    
    return c
end


function calc_c(M200::Real, r_s::Real)
    return calc_R200(M200) / r_s
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


@doc raw"""
A truncated NFW profile

The density profile is given by

```math
\rho(r) = \rho_{NFW}(r) \exp(-r/r_t)
```
"""
struct TruncNFW <: SphericalProfile
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


function calc_ρ(profile::TruncNFW, r::Real)
    nfw = NFW(profile.M_s, profile.r_s, profile.c)
    return calc_ρ(nfw, r) * exp(-(r/profile.r_t))
end

function calc_M(profile::TruncNFW, r::Real)
    x = r / profile.r_s
    t = profile.r_t / profile.r_s
    B = (1+1/t)*exp(1/t)
    return profile.M_s * ( B * expinti(-(x+1)/t) + exp(-x/t)/(1+x) - B*expinti(-1/t) - 1 )
end


function calc_M_tot(profile::TruncNFW)
    t = profile.r_t / profile.r_s
    A = -(1+1/t)*exp(1/t)*expinti(-1/t) - 1

    return profile.M_s * A
end


function calc_Φ(profile::TruncNFW, r::Real)
    Φ_in = - G * calc_M(profile, r) / r

    Φ0 = -G * profile.M_s / profile.r_s
    x = r / profile.r_s
    t = profile.r_t / profile.r_s
    Φ_out = Φ0 * (1/(1+x)*exp(-x/t) + expinti(-(1+x)/t)*exp(1/t)/t)

    return Φ_in + Φ_out
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
        c_ludlow(M, z)

        Calculates the approximate concentration of a halo given the initial
        mass (M200) and redshift using the n-body fit from @ludlow2016
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

    Solves for the mass and concentration of a halo given the maximum circular
    velocity. Calls Roots.find_zero to solve for the mass on the given interval
    and can apply a multiplicative factor `δlogc` to the concentration
    parameter.

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

    Solves for the radius of maximum circular velocity given the maximum
    circular velocity in code units and an optional offset to the concentration mass
    relation, `δlogc`. Calls `solve_M200_c` to solve for the mass and
    concentration of the halo, passing along kwargs.
    """
    function solve_rmax(v_circ_max, δlogc=0; kwargs...)
        M200, c = solve_M200_c(v_circ_max, δlogc; kwargs...)
        return calc_r_circ_max(NFW(M200=M200, c=c))
    end

end
