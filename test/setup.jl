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
