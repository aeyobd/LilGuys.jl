### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 43307c80-298c-11ef-3cfe-afb7123432ca
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	import LilGuys as lguys

	using Arya
end

# ╔═╡ 9ecf3fea-823f-4c2f-82cc-f33d5f3ff787
md"""
# Random samples & Recovery

To ensure visually that my density methods work well, this notebook creates
a randomly generated sample of stars using a spherically symmetric king profile in 3d, and compares the resulting density profiles we calculate with the expected analytic profile. Especially as the de-projection formulae can be a little more complex, this provides nice peace of mind that these formula are valid for painting on stars.
"""

# ╔═╡ ffb3c6d0-a9da-41ae-91f0-5a6fa2c1bdf9
import LilGuys.Plots as LP

# ╔═╡ ce17db49-ab20-4f6b-bb10-c08e970ff156
md"""
## The random sample

This code just creates the random sample of stars. The library call used here in LilGuys is  inverse-cdf sampling, i.e. the probability distribution is calculated at a number of points and then is linearly interpolated using an inverse function.
"""

# ╔═╡ 1c55d5eb-1065-4ef0-bb97-ada81c195432
N = 10000

# ╔═╡ d586226b-37dc-4f56-b632-d8d411848401
R_s=π / 100

# ╔═╡ 97ae40b8-5b48-4828-a139-549d8cc9a124
profile = lguys.KingProfile(R_s=R_s, c=5, M=1)

# ╔═╡ 1cc8aaaa-b20b-401c-be78-6162da324d7e
radii = lguys.sample_ρ(x -> lguys.calc_ρ.(profile, x), N)

# ╔═╡ 194f9836-58cf-43ec-9689-f30bcc1feba1
pos = radii' .* lguys.rand_unit(N)

# ╔═╡ 18e0fe22-54d0-4f94-bb87-2f4c5a30a145
vel = zeros(size(pos));

# ╔═╡ b59b12b8-ffc4-4c4a-840c-6cce3160818c
md"""
The histogram below simply verifies that the sampling worked reasonably well
"""

# ╔═╡ 444ba26c-1f99-467e-a304-3879190dc03a
hist(log10.(radii),
	axis=(; xlabel="log radius", ylabel="counts")
)

# ╔═╡ b0dc752d-e206-426f-9ed6-84107230fce1
md"""
# Testing the 3D profile
"""

# ╔═╡ 07e7f1a4-8cea-4832-8d61-5a0e6ef5f48d
begin 
	snap = lguys.Snapshot(pos, vel, ones(N) / N)
	snap.weights = snap.masses
end

# ╔═╡ 3f172c56-8906-4e5e-a984-6a1a985ebeaf
prof_3d = lguys.MassProfile3D(snap)

# ╔═╡ b2ea83f1-375c-41c5-8eab-4a436e4c46bb
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		yscale=log10,
		xlabel="log radius / kpc",
		ylabel=L"$\log\;\rho$ [code units]",
		limits=(-3, -0.5, 1, 1e5)
	)

	scatter!(prof_3d.log_r, prof_3d.rho)

	x = LinRange(-3, 0, 1000)
	y = lguys.calc_ρ.(profile, 10 .^ x)
	
	vlines!(log10.(R_s))
	vlines!(log10(lguys.calc_r_h(profile)))
	r_h1 = lguys.quantile(radii, 0.5)
	
	vlines!(log10(r_h1))

	
	lines!(x, y)
	fig
end

# ╔═╡ 0ed5ce67-55ed-44db-89ca-0c27b63a2a26
md"""
## 2D orthoganal projection
"""

# ╔═╡ 1befbfe1-65d4-4473-95aa-6cc456714a11
R = @. sqrt(pos[1, :]^2 + pos[2, :]^2)

# ╔═╡ 6cce5fb9-c12f-4af6-b188-a79c53c37fa5
prof_2d_ortho = lguys.StellarProfile(R, bins=100)

# ╔═╡ 32e87484-508c-4280-8d38-c9275ba70c74
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		limits=(-3, -0.5, -3, 4),
		xlabel = "R / kpc",
		ylabel = L"$\Sigma$ / [code units]"
	)

	scatter!(prof_2d_ortho.log_r, prof_2d_ortho.log_Sigma)

	x = 10 .^ LinRange(-3, 2.2, 1000)
	y = lguys.calc_Σ.(profile, x)
	
	lines!(log10.(x), log10.(y), color=COLORS[2])
	fig
end

# ╔═╡ baa4da2c-5f8c-4a58-ac7b-f3a77fc8f164
md"""
## Sky Projection
"""

# ╔═╡ a524219c-9d19-4858-a475-eed124f9c8a1
md"""
for simplicity, we do not shift the location of the snapshot, so the stars should appear on the sky at the location of Sagitarius A*
"""

# ╔═╡ 27aa2cd1-b365-4701-b698-fb2fb7d0f3eb
distance = 8.122 # distance to sag A* in kpc

# ╔═╡ 3637933e-7454-44da-a396-5d9e9e47dffd
gaia = lguys.to_gaia(snap, add_centre=false)

# ╔═╡ b8a32ed2-c7ea-4d0d-a9b4-59e97a31db61
md"""
Should be centred at Sag A*
"""

# ╔═╡ f8f890a1-2dc6-49f5-8756-bac581af478c
ra0, dec0 = lguys.calc_centre2D(gaia.ra, gaia.dec, "mean")

# ╔═╡ f5ff4ee3-c437-4152-82ce-b8747dd4a6a3
xi, eta = lguys.to_tangent(gaia.ra, gaia.dec, ra0, dec0)

# ╔═╡ 7dd90b94-1ea2-475d-9742-6eddab8baaac
r_ell = lguys.calc_r_ell(xi, eta, 1, 1, 0) * 60 # arcmins

# ╔═╡ a75a7819-3d25-4455-bd05-f8482b79d2e5
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel="xi / degrees", ylabel="eta / degrees",
		aspect=DataAspect(),
		xgridvisible=false,
		ygridvisible=false,
	)

	
	scatter!(xi, eta, alpha=0.2, color=:black, markersize=3)

	fig
end

# ╔═╡ 24fd05e6-1d01-4909-b348-21192b26861c
props_sky = lguys.StellarProfile(r_ell, weights=snap.masses)

# ╔═╡ 3ab11ec8-1899-486c-97c5-0133607309c8
R_s_arcmin = lguys.kpc_to_arcmin(R_s, distance)

# ╔═╡ 94d54949-1118-40fe-9e68-62b716bbc1ae
prof_sky = lguys.KingProfile(R_s=R_s_arcmin, c=5, M=1)

# ╔═╡ 132b0c72-27fd-48e1-86b6-6a5291a2c145
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel="log r / degrees",
		ylabel=L"\log \Sigma",
		limits=(-1, 2, -6, -2)
	)

	scatter!(props_sky.log_r, props_sky.log_Sigma)

	x = 10 .^ LinRange(-4, 2.7, 1000)
	y = lguys.calc_Σ.(prof_sky, x)
	
	lines!(log10.(x), log10.(y), color=COLORS[2])
	fig
end

# ╔═╡ ede57bf6-c009-4c77-a439-e41121c3beca
filename = "test_sky_recon.fits"

# ╔═╡ 9965c814-bf78-4b0b-9818-778029786afa
lguys.write_fits(filename, gaia)

# ╔═╡ 909c3ec2-f947-454f-bea8-095e196a2fab
md"""
## Automated stellar profile
"""

# ╔═╡ 4661b78d-2ea0-45c3-8e2e-c1d5bcf3d7f7
md"""
run 
```
stellar_profile.jl test_sky_recon.fits --mass-column weights
```
on above and then fit_profile
"""

# ╔═╡ 772bcaa1-9e99-471e-bbd0-7316a4ac65b9
prof_recon = lguys.StellarProfile("test_sky_recon_profile.toml")

# ╔═╡ edfb1e71-90b7-407f-bfdc-8bb3d0dba95e
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel="log r / arcmin",
		ylabel=L"\log \Sigma",
		limits=(-1, 2, -6, -2)
	)

	errscatter!(prof_recon.log_r, prof_recon.log_Sigma, yerr=prof_recon.log_Sigma_err)

	x = 10 .^ LinRange(-1, 3, 1000)
	y = lguys.calc_Σ.(prof_sky, x)
	
	lines!(log10.(x), log10.(y), color=COLORS[2])
	
	fig
end

# ╔═╡ Cell order:
# ╟─9ecf3fea-823f-4c2f-82cc-f33d5f3ff787
# ╠═43307c80-298c-11ef-3cfe-afb7123432ca
# ╠═ffb3c6d0-a9da-41ae-91f0-5a6fa2c1bdf9
# ╠═ce17db49-ab20-4f6b-bb10-c08e970ff156
# ╠═1c55d5eb-1065-4ef0-bb97-ada81c195432
# ╠═d586226b-37dc-4f56-b632-d8d411848401
# ╠═97ae40b8-5b48-4828-a139-549d8cc9a124
# ╠═1cc8aaaa-b20b-401c-be78-6162da324d7e
# ╠═194f9836-58cf-43ec-9689-f30bcc1feba1
# ╠═18e0fe22-54d0-4f94-bb87-2f4c5a30a145
# ╟─b59b12b8-ffc4-4c4a-840c-6cce3160818c
# ╠═444ba26c-1f99-467e-a304-3879190dc03a
# ╠═b0dc752d-e206-426f-9ed6-84107230fce1
# ╠═07e7f1a4-8cea-4832-8d61-5a0e6ef5f48d
# ╠═3f172c56-8906-4e5e-a984-6a1a985ebeaf
# ╠═b2ea83f1-375c-41c5-8eab-4a436e4c46bb
# ╟─0ed5ce67-55ed-44db-89ca-0c27b63a2a26
# ╠═1befbfe1-65d4-4473-95aa-6cc456714a11
# ╠═6cce5fb9-c12f-4af6-b188-a79c53c37fa5
# ╠═32e87484-508c-4280-8d38-c9275ba70c74
# ╟─baa4da2c-5f8c-4a58-ac7b-f3a77fc8f164
# ╠═a524219c-9d19-4858-a475-eed124f9c8a1
# ╠═27aa2cd1-b365-4701-b698-fb2fb7d0f3eb
# ╠═3637933e-7454-44da-a396-5d9e9e47dffd
# ╟─b8a32ed2-c7ea-4d0d-a9b4-59e97a31db61
# ╠═f8f890a1-2dc6-49f5-8756-bac581af478c
# ╠═f5ff4ee3-c437-4152-82ce-b8747dd4a6a3
# ╠═7dd90b94-1ea2-475d-9742-6eddab8baaac
# ╠═a75a7819-3d25-4455-bd05-f8482b79d2e5
# ╠═24fd05e6-1d01-4909-b348-21192b26861c
# ╠═3ab11ec8-1899-486c-97c5-0133607309c8
# ╠═94d54949-1118-40fe-9e68-62b716bbc1ae
# ╠═132b0c72-27fd-48e1-86b6-6a5291a2c145
# ╠═ede57bf6-c009-4c77-a439-e41121c3beca
# ╠═9965c814-bf78-4b0b-9818-778029786afa
# ╟─909c3ec2-f947-454f-bea8-095e196a2fab
# ╟─4661b78d-2ea0-45c3-8e2e-c1d5bcf3d7f7
# ╠═772bcaa1-9e99-471e-bbd0-7316a4ac65b9
# ╠═edfb1e71-90b7-407f-bfdc-8bb3d0dba95e
