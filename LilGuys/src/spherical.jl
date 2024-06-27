"""
    to_tangent(α, δ, α_0, δ_0)

Computes the tangent plane coordinates of a point (α, δ) with respect to a reference point (α_0, δ_0).
"""
function to_tangent(α, δ, α_0, δ_0)
	denom = @. (sind(δ) * sind(δ_0) 
		+ cosd(δ) * cosd(δ_0) * cosd(α - α_0)
	)
	
	eta_num = @. (sind(δ_0) * cosd(δ) * cosd(α-α_0)
		-cosd(δ_0) * sind(δ) 
	)
	
	xi_num = @. cosd(δ) * sind(α - α_0)
	
	xi = @. rad2deg(xi_num/denom)
	eta = @. -rad2deg(eta_num / denom)

	return xi, eta
end



"""
    unit_vector(ra::F, dec::F)

Returns the unit vector pointing at the posision (ra, dec) in degrees on the sky.
"""
function unit_vector(ra::F, dec::F)
    x = cosd(dec) * cosd(ra)
    y = cosd(dec) * sind(ra)
    z = sind(dec)
    return [x, y, z]
end


function unit_vector(ra::Vector{F}, dec::Vector{F}) where F
    x = cosd.(dec) .* cosd.(ra)
    y = cosd.(dec) .* sind.(ra)
    z = sind.(dec)
    return [x y z]
end


function cartesian_to_sky(x, y, z)

    R = sqrt(x^2 + y^2)
    r = sqrt(x^2 + y^2 + z^2)

    ra = mod(atand(y, x), 360) # atan2 -> 0 to 360
    dec = atand(z / R) # atan -> -90 to 90

    return [ra, dec, r]
end


function cartesian_to_sky(x::Vector{F}, y::Vector{F}, z::Vector{F}) where F
    R = sqrt.(x.^2 + y.^2)
    r = sqrt.(x.^2 + y.^2 + z.^2)

    ra = mod.(atand.(y, x), 360) # atan2 -> 0 to 360
    dec = atand.(z ./ R) # atan -> -90 to 90

    return [ra dec r]
end


"""
    Rx_mat(u)

rotation matrix around x axis (radians)
"""
function Rx_mat(u::Real)
    c = cos(u)
    s = sin(u)
    return [1  0  0
            0  c  s
            0 -s  c]
end


"""
    Ry_mat(v)

rotation matrix around y axis (radians)
"""
function Ry_mat(v::Real)
    c = cos(v)
    s = sin(v)
    return [c  0  s
            0  1  0
           -s  0  c]
end


"""
    Rz_mat(w)

rotation matrix around z axis (radians)
"""
function Rz_mat(w::Real)
    c = cos(w)
    s = sin(w)

    return [c -s  0
            s  c  0
            0  0  1]
end


function R_mat(u::Real, v::Real, w::Real)
    return Rz_mat(w) * Ry_mat(v) * Rx_mat(u)
end


