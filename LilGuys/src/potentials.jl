
struct HernquistParams
    M::F
    a::F
end


function Φ_hernquist(r::F, params)
    return -G * params.M / (r + params.a)
end

function Φ_hernquist(x_vec::Vector{F}, params)
    r = calc_r(x_vec)
    return Φ_hernquist(r, params)
end


function a_hernquist(r::F, params)
    return -G * params.M / (r + params.a)^2
end


function a_hernquist(x_vec::Vector{F}, hernquist_params)
    r = calc_r(x_vec)
    r_hat = x_vec ./ r
    return -r_hat * a_hernquist(r, hernquist_params)
end



struct NFWParams
    M_s::F
    r_s::F
end

function NFWParams(; M200::F, c::F)
    M_s = M200 / A_nfw(c)
    return NFWParams(M_s, r_s)
end


"""
The nfw scaling function A(x) = log(1 + x) - x / (1 + x)
"""
function A_nfw(x::F)
    return log(1 + x) - x / (1 + x)
end


function a_nfw(r, nfw_params)
    GM = G * nfw_params.M_s
    r_s = nfw_params.r_s
    x = r / r_s
    return -GM / r^2 * A_nfw(x) / x^2
end

function Φ_nfw(r, nfw_params)
    GM = G * nfw_params.M_s
    r_s = nfw_params.r_s
    x = r / r_s
    return -GM / r * log(1 + x) / x
end


function a_miyamoto(r, miyamoto_params)

end

function Φ_miyamoto(r, miyamoto_params)

end

