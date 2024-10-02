### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ eaa321be-5d93-11ef-3e83-5b622d34bc3d
begin
	import Pkg; Pkg.activate()

	using Arya
	using CairoMakie
	
	using LilGuys
end

# ╔═╡ cf73ae3a-65a3-447f-a407-1c62cde4601c
md"""
# Synthetic recovery
By creating random samples of stars, what is our accuracy with simple density recovery?
"""

# ╔═╡ 14a0f47b-872b-4049-9d65-876c0beb73b8
import DensityEstimators: bins_min_width_equal_number

# ╔═╡ 869a05a9-efac-4810-99cd-6364dc7201c7
f_Σ(x) = exp(-x) / 2π

# ╔═╡ 060851b5-50c4-4c02-aab0-48206fde07f6
Γ_exp(x) = -x 

# ╔═╡ b4d30aef-95a6-4579-9965-54922272b3f4
N = 10000

# ╔═╡ 0108c0de-6826-44e6-8896-28fe6d48a58d
M_exp(x) = N * (1 - exp(-x) - x*exp(-x))

# ╔═╡ d38e8e46-a9b4-4251-a106-3582e9faeb88
radii = LilGuys.sample_Σ(f_Σ, N)

# ╔═╡ 8f5ebb2d-1c17-4e03-93f5-947fda8b35a8
hist(radii)

# ╔═╡ 8c461d77-2100-457c-8430-59a9b6e6ddef
prof = LilGuys.StellarProfile(radii, normalization=:none, bins= bins_min_width_equal_number(log10.(radii), dx_min=0.05, N_per_bin_min=50))

# ╔═╡ cee96228-eb9d-4819-8d6b-34d79266bfb7
let
	fig, ax = FigAxis(xlabel="log r", ylabel="log Σ")

	errscatter!(prof.log_r, prof.log_Sigma, yerr=prof.log_Sigma_err)

	x = LinRange(-2, 1, 100)
	y = f_Σ.(10 .^ x)

	lines!(x, log10.(y * N))
	fig
end

# ╔═╡ 603bbf97-dc7f-4a2b-a554-e4bfad4b5ef2
let
	fig, ax = FigAxis(xlabel="log r", ylabel="M(r)")

	errscatter!(prof.log_r_bins[2:end], prof.M_in, yerr=prof.M_in_err)

	x = LinRange(-2, 1, 100)
	y = M_exp.(10 .^ x)

	lines!(x, y)
	fig
end

# ╔═╡ dedc1ab1-250c-4df4-b55f-9f1dac0e9d4b
let
	fig, ax = FigAxis(xlabel="log r", ylabel="Γ")

	errscatter!(prof.log_r, prof.Gamma, yerr=prof.Gamma_err)

	x = LinRange(-2, 1, 100)
	y = Γ_exp.(10 .^ x)

	lines!(x, y)
	fig
end

# ╔═╡ bad9ee8e-7338-40d8-a40a-4cfb6d3ad023
function calc_χ2(x, x_exp, xerr)
	return sum(
		@. (x-x_exp)^2 / xerr^2
	) / length(x)
end

# ╔═╡ 40a68cd5-6d7d-4421-8d04-718804c14adf
calc_χ2(prof.Sigma, N * f_Σ.(10 .^ prof.log_r), prof.Sigma_err)

# ╔═╡ Cell order:
# ╟─cf73ae3a-65a3-447f-a407-1c62cde4601c
# ╠═eaa321be-5d93-11ef-3e83-5b622d34bc3d
# ╠═14a0f47b-872b-4049-9d65-876c0beb73b8
# ╠═869a05a9-efac-4810-99cd-6364dc7201c7
# ╠═0108c0de-6826-44e6-8896-28fe6d48a58d
# ╠═060851b5-50c4-4c02-aab0-48206fde07f6
# ╠═b4d30aef-95a6-4579-9965-54922272b3f4
# ╠═d38e8e46-a9b4-4251-a106-3582e9faeb88
# ╠═8f5ebb2d-1c17-4e03-93f5-947fda8b35a8
# ╠═8c461d77-2100-457c-8430-59a9b6e6ddef
# ╠═cee96228-eb9d-4819-8d6b-34d79266bfb7
# ╠═603bbf97-dc7f-4a2b-a554-e4bfad4b5ef2
# ╠═dedc1ab1-250c-4df4-b55f-9f1dac0e9d4b
# ╠═bad9ee8e-7338-40d8-a40a-4cfb6d3ad023
# ╠═40a68cd5-6d7d-4421-8d04-718804c14adf
