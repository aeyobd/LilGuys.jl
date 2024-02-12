
function ϕ_eff(positions, velocities, masses; h=0.08)
	N = length(masses)
	ϕ = zeros(N)

	for i in 1:N
		if i % 100 == 0
			print("\r $i / $N")
		end

		ϕ_g = @. -masses / sqrt(calc_r(positions[:, i], positions)^2 + h^2)
		ϕ_v =  1/N * 1/2 * calc_r(velocities[:, i], velocities).^2
			
		ϕ[i] = sum(min.(ϕ_g .+ ϕ_v, 0))
		
	end
	return ϕ

end

function calc_ρ_eff(x_vec, v_vec, positions, velocities, masses; h)
    N = length(masses)

    ϕ_g = @. -masses / sqrt(calc_r(x_vec, positions)^2 + h^2)
    ϕ_v =  1/N * 1/2 * calc_r(v_vec, velocities).^2
    ϕs = @. min(ϕ_g + ϕ_v, 0)
    return sum(ϕs)
end
