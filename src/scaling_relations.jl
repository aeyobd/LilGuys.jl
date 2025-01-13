"""
    M_s_from_vel_fattahi(v_max)

Given the maximum circular velocity (in code units), returns the stellar mass (in code units).
Uses the fit from @fattahi2018.
"""
function M_s_from_vel_fattahi(v_max)
	ν = v_max * V2KMS / 50 # km/s
	α = 3.36
	γ = -2.4
	m_0 = 3e-2 # 1e10 Msun
	return m_0 * ν^α * exp(-ν^γ)
end


"""
    vel_from_M_s_fattahi(Mstar)

Inverse of `M_s_from_vel_fattahi`. Given the stellar mass of a galaxy, returns the maximum circular velocity of the halo (in code units).
"""
function vel_from_M_s_fattahi(Mstar) # code units
    return find_zero(ν -> M_s_from_vel_fattahi(ν) - Mstar, 0.1)
end




@doc raw"""
    v_circ_tidal_EN21(x)

Returns the relative v_circ_max of a halo undergoing tidal evolution as a function of r_circ_max / r_circ_max_initial.

The relationship is from Errani & Navarro (2021), equation 5:
``
V / V_0 = 2^α\ x^β \ (1 + x^2)^(-α)
``

where $x = r_circ_max / r_circ_max_initial$, and $α = 0.4$, $β = 0.65$.
"""
function v_circ_EN21(x)
    α = 0.4
    β = 0.65
    y = @. 2^α * x^β * (1 + x^2)^(-α)
end


"""
    EN21_tidal_track(r_circ_max_0, v_circ_max_0; x_min = 0.1, n_points = 1000)

Computes the tidal track as described in Errani & Navarro (2021) for a halo with initial r_circ_max and v_circ_max. 

See also [`v_circ_EN21`](@ref). 

Returns the r_circ_max and v_circ_max of the halo as it evolves under tidal forces over a range from r_circ_max_initial to r_circ_max * x_min.
"""
function EN21_tidal_track(r_circ_max_0, v_circ_max_0; x_min = 0.1, n_points = 1000)
    x = 10 .^ LinRange(0, log10(x_min), n_points)
    v_circ_max = v_circ_max_0 * v_circ_EN21(x)
    r_circ_max = r_circ_max_0 * x
    return r_circ_max, v_circ_max
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
