#!/usr/bin/env julia

import LilGuys as lguys

using DataFrames
using HDF5

import Statistics: quantile
import TOML
import QuadGK: quadgk

import DensityEstimators

using ArgParse

include("script_utils.jl")


function get_args()
    s = ArgParseSettings(
        description="""Creates a N-body model by using
Eddington inversion on the density profile.
""",
    )

    @add_arg_table s begin
        "input"
            help="Input halo"
            required=true
        "-o", "--output_file"
            help="output file"
            required=true
        "-n", "--number"
            help="number of stars to create"
        "--num-energy-bins"
            help="number of energy bins"
            default=100
            
    end

    args = parse_args(s)

    params = load_params(args["input"])

    return params
end


function main(args)
    halo = lguys.load_profile(params["input"])

    r = 10 .^ LinRange(-5, 5, 10_000)
    Ψ = -lguys.calc_Φ.(halo, r)
    ρ = lguys.calc_ρ.(halo, r)
    M = lguys.calc_M.(halo, r)
    M_inv = lguys.lerp(M ./ M[end], r)

    # distribution functions
    E = make_energy_bins(ψ, params)
    f_dm = lguys.DistributionFunction(ρ, ψ, r)
    dP_dr = 

    f = lguys.lerp(E, f_dm)
    if any(f_dm .< 0)
        @warn "negative distribution function"
    end

    # sample
    
    N = params["number"]
    radii = M_inv.(rand(N))
    Φs = lguys.calc_Φ.(halo, radii)
    # velocity via rejection sampling
    v_max = sqrt.(2Φs)

    v = v_max .* rand(N)
    while true
        v = v_max .* rand(N)
        if all(rand(N) .< lguys.calc_f(v, f, E))
            break
        end
    end
end


function sample_r(M_inv)
    M_inv(rand())
end

function sample_rv(M_inv, p_v)
    r = sample_r(M_inv)
    v = p_v(r) * rand()
    return r, v
end


"""
loads in the parameterfile and will 
automatically populate the centres file and output files
"""
function load_params(paramname)
    params = TOML.parsefile(paramname); 

    if "output_file" ∉ keys(params)
            params["output_file"] = splitext(paramname)[1]
    end
    
    if "mock_file" ∉ keys(params)
            params["mock_file"] = splitext(paramname)[1] * "_mock_stars.fits"
    end
    
    params
end


function make_energy_bins(ψ, params)
	E_max = ψ[1]
	E_min = ψ[end]
	E = LinRange(E_min, E_max, params["num-energy-bins"] + 1)
end



"""
Given the stellar mass function M_s, returns total missing stars
"""
function print_missing(radii, r_bins, profile)
    r_h = lguys.calc_r_h(profile)
    @info "r_h = $r_h"
    @info " $(sum(radii .< r_h)) stars within (3D) half-light radius"

    M_s_tot = lguys.get_M_tot(profile)
    M_s(r) = lguys.calc_M(profile, r) / M_s_tot

	N_s_out = 1 - M_s(r_bins[end])
	N_s_in = M_s(r_bins[1])
	@info "missing $N_s_out stars outside last bin"
	@info "missing $N_s_in stars inside first bin"
end



if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    run_script_with_output(main, args)
end
