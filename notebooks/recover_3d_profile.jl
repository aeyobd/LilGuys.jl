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

By creating random samples of dark matter particles (only with positions) in an NFW halo, how closely does our methods resemble the real profile?

"""

# ╔═╡ 14a0f47b-872b-4049-9d65-876c0beb73b8
import DensityEstimators: bins_min_width_equal_number

# ╔═╡ 0994c436-fb69-4056-89e8-e8aa1f3bb603
halo = LilGuys.TruncNFW(M_s = 1.5, r_s=5, trunc=100)

# ╔═╡ 869a05a9-efac-4810-99cd-6364dc7201c7
f_ρ(x) = calc_ρ(halo, x)

# ╔═╡ 0108c0de-6826-44e6-8896-28fe6d48a58d
M_exp(x) = calc_M(halo, x)

# ╔═╡ 060851b5-50c4-4c02-aab0-48206fde07f6
Γ_exp(x) = -x 

# ╔═╡ b4d30aef-95a6-4579-9965-54922272b3f4
N = 10_000

# ╔═╡ d38e8e46-a9b4-4251-a106-3582e9faeb88
radii = LilGuys.sample_ρ(f_ρ, N, log_r=LinRange(-5, 5, 10000))

# ╔═╡ 8f5ebb2d-1c17-4e03-93f5-947fda8b35a8
hist(radii)

# ╔═╡ 5eacb7c5-edcb-4afa-b605-0b9f7cd71d8e
snap = Snapshot(radii' .* LilGuys.rand_unit(N), zeros(3, N), ones(N) * M_exp(maximum(radii)) / N)

# ╔═╡ 1b420832-8f58-44e6-add5-b20e324e951b
bins = bins_min_width_equal_number(log10.(radii), dx_min=0.05, N_per_bin_min=100)#[2:end-1]

# ╔═╡ 8c461d77-2100-457c-8430-59a9b6e6ddef
prof = LilGuys.MassProfile3D(snap, bins=bins)

# ╔═╡ cee96228-eb9d-4819-8d6b-34d79266bfb7
let
	fig, ax = FigAxis(xlabel="log r", ylabel=L"\log\,\rho")

	log_x = prof.log_r
	x = 10 .^ log_x
	
	errscatter!(log_x, log10.(prof.rho))

	xpred = LinRange(-1, 3.5, 100)
	ypred = f_ρ.(10 .^ xpred)

	lines!(xpred, log10.(ypred), color=COLORS[2])
	hidexdecorations!(ax, grid=false)

	res = ( prof.rho .- f_ρ.(x)) ./ f_ρ.(x)

	ax2 = Axis(fig[2, 1], ylabel="residual", xlabel="log r")

	errscatter!(log_x, res, yerr=prof.rho_err ./ prof.rho)
	hlines!(0, color=:black)

	rowsize!(fig.layout, 2, Relative(0.3))
	linkxaxes!(ax, ax2)
	fig
end

# ╔═╡ 603bbf97-dc7f-4a2b-a554-e4bfad4b5ef2
let
	fig, ax = FigAxis(xlabel="log r", ylabel="M(r)")

	log_x = prof.log_r_bins[2:end]
	x = 10 .^ log_x
	
	errscatter!(log_x, prof.M_in, yerr=prof.M_in_err)

	xpred = LinRange(-1, 3.5, 100)
	ypred = M_exp.(10 .^ xpred)

	lines!(xpred, ypred, color=COLORS[2])
	hidexdecorations!(ax, grid=false)

	res = (prof.M_in .- M_exp.(x)) ./ prof.M_in

	ax2 = Axis(fig[2, 1], 
		ylabel="residual", xlabel="log r",
		limits=(nothing, nothing, -0.03, 0.03)
	)

	errscatter!(log_x, res, yerr=prof.M_in_err ./ prof.M_in)
	hlines!(0, color=:black)

	rowsize!(fig.layout, 2, Relative(0.3))
	linkxaxes!(ax, ax2)
	fig
end

# ╔═╡ 4b02472a-701b-4d3f-8fdf-0ae3a8b49832
let
	fig, ax = FigAxis(xlabel="log r", ylabel="v circ")

	log_x = prof.log_r_bins[2:end]
	x = 10 .^ log_x
	
	errscatter!(log_x, prof.v_circ, yerr=prof.v_circ_err)

	xpred = LinRange(-1, 3.5, 100)
	ypred = calc_v_circ.(halo, 10 .^ xpred)

	lines!(xpred, ypred, color=COLORS[2])
	hidexdecorations!(ax, grid=false)

	res = (prof.v_circ .- calc_v_circ.(halo, x)) ./ prof.v_circ

	ax2 = Axis(fig[2, 1], 
		ylabel="residual", xlabel="log r",
		limits=(nothing, nothing, -0.03, 0.03)
	)

	errscatter!(log_x, res, yerr=prof.v_circ_err ./ prof.v_circ)
	hlines!(0, color=:black)

	rowsize!(fig.layout, 2, Relative(0.3))
	linkxaxes!(ax, ax2)
	fig
end

# ╔═╡ 94eaa3a2-fa51-4032-b0a5-6cd8a06a149f
cumsum(prof.mass_in_shell)

# ╔═╡ bad9ee8e-7338-40d8-a40a-4cfb6d3ad023
function calc_χ2(x, x_exp, xerr)
	return sum(
		@. (x-x_exp)^2 / xerr^2
	) / length(x)
end

# ╔═╡ 40a68cd5-6d7d-4421-8d04-718804c14adf
calc_χ2(prof.rho[2:end-1], f_ρ.(10 .^ prof.log_r[2:end-1]), prof.rho_err[2:end-1])

# ╔═╡ cfacda96-2c99-4c1b-b86b-65c9d4383ca1
calc_χ2(prof.M_in, M_exp.(10 .^ prof.log_r_bins[2:end]), prof.M_in_err)

# ╔═╡ 0a14c5d1-aeff-4c0a-81a9-4b258b35b26c
calc_χ2(prof.v_circ, calc_v_circ.(halo, 10 .^ prof.log_r_bins[2:end]), prof.v_circ_err)

# ╔═╡ Cell order:
# ╟─cf73ae3a-65a3-447f-a407-1c62cde4601c
# ╠═eaa321be-5d93-11ef-3e83-5b622d34bc3d
# ╠═14a0f47b-872b-4049-9d65-876c0beb73b8
# ╠═0994c436-fb69-4056-89e8-e8aa1f3bb603
# ╠═869a05a9-efac-4810-99cd-6364dc7201c7
# ╠═0108c0de-6826-44e6-8896-28fe6d48a58d
# ╠═060851b5-50c4-4c02-aab0-48206fde07f6
# ╠═b4d30aef-95a6-4579-9965-54922272b3f4
# ╠═d38e8e46-a9b4-4251-a106-3582e9faeb88
# ╠═8f5ebb2d-1c17-4e03-93f5-947fda8b35a8
# ╠═5eacb7c5-edcb-4afa-b605-0b9f7cd71d8e
# ╠═1b420832-8f58-44e6-add5-b20e324e951b
# ╠═8c461d77-2100-457c-8430-59a9b6e6ddef
# ╠═cee96228-eb9d-4819-8d6b-34d79266bfb7
# ╠═603bbf97-dc7f-4a2b-a554-e4bfad4b5ef2
# ╠═4b02472a-701b-4d3f-8fdf-0ae3a8b49832
# ╠═94eaa3a2-fa51-4032-b0a5-6cd8a06a149f
# ╠═bad9ee8e-7338-40d8-a40a-4cfb6d3ad023
# ╠═40a68cd5-6d7d-4421-8d04-718804c14adf
# ╠═cfacda96-2c99-4c1b-b86b-65c9d4383ca1
# ╠═0a14c5d1-aeff-4c0a-81a9-4b258b35b26c
