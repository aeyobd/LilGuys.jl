using Test
using Makie # need for makie ext
using PythonCall # needed for fits IO
import Agama
import LilGuys as lguys

using LilGuys

using HDF5

import Random
import Distributions: Chisq, quantile

import Logging

Random.seed!(1234)


macro test_χ2(measurements, errors, expected, p=0.9995)
    return quote
        using Test

        y1 = $(esc(measurements))
        y2 = $(esc(expected))
        err = $(esc(errors))

        @assert length(y1) == length(y2) == length(err)

        χ2 = ((y1 .- y2) ./ err).^2
        filt = @. isfinite(y1) & isfinite(y2) & isfinite(err)
        χ2 = sum(χ2[filt])
        ndof = sum(filt) - 1
        
        chimax = quantile(Chisq(ndof), $p)
        # A simple test for consistency: chi2 should not exceed a chosen threshold, say 3 times ndof
        if !(χ2 < chimax)
            @warn "χ2 = $χ2 > $chimax"
        end

        @test χ2 < chimax 
        χ2
    end
end

macro test_χ2(measurements, expected)
    p = 0.9995
    return quote
        using Test
        import LilGuys: middle, error_interval

        y1 = $(esc(measurements))
        y2 = $(esc(expected))

        @assert length(y1) == length(y2)

        err = maximum.(error_interval.(y1))
        y1 = middle.(y1)

        χ2 = ((y1 .- y2) ./ err).^2
        filt = @. isfinite(y1) & isfinite(y2) & isfinite(err)
        χ2 = sum(χ2[filt])
        ndof = sum(filt) - 1
        
        chimax = quantile(Chisq(ndof), $p)
        # A simple test for consistency: chi2 should not exceed a chosen threshold, say 3 times ndof
        if !(χ2 < chimax)
            @warn "χ2 = $χ2 > $chimax"
        end

        @test χ2 < chimax 
        χ2
    end
end


function snap_from_density(halo, N=10_000; log_r_interp=LinRange(-5, 5, 1000))
    M_0 = lguys.mass(halo)
    ρ(r) = lguys.density(halo, r)
    r = lguys.sample_density(ρ, N, log_r=LinRange(-5, 5, 10_000))

    mass = M_0/N  * (1 .+ 0.0randn(N))
    M = sum(mass)

    snap = lguys.Snapshot(positions=r' .* lguys.rand_unit(N), velocities=zeros(3, N), masses=mass, index=1:N, header=lguys.make_default_header(1, N), potential=-ones(N))

end
