
import Interpolations: linear_interpolation, Line

function ρ_E21(x, x0, r_s, r_t, c)
end

"""
given a centered snapshot, returns a interpolated potential a a function 
of r"""
function calc_radial_Φ(snap::Snapshot)
    # work outside in 
    rs = sort(r(snap.pos))
    m = snap.m[1:end]

    N = length(snap)
    M_inside = cumsum(m)

    Φ_inside = -G * M_inside ./ rs
    Φ_outside = -G * [ sum(snap.m[i] ./ rs[i:end]) for i in 1:N]
    push!(Φ_outside, 0)

    pushfirst!(rs, 0)
    pushfirst!(Φ_inside, 0)

    Φs = Φ_inside .+ Φ_outside

    return linear_interpolation(rs, Φs, extrapolation_bc=Line())
end
