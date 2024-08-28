function cal_Δz(z₊ₕ)
  n = length(z₊ₕ)
  Δz = zeros(n)
  Δz[1] = 0 - z₊ₕ[1]
  for i in 2:n
    Δz[i] = z₊ₕ[i-1] - z₊ₕ[i]
  end
  Δz
end

# Δz: 1:N
# Δz₊ₕ: 0:N
function cal_Δz₊ₕ(z)
  z₊ₕ = C2F(z)
  n = length(z)
  Δz₊ₕ = zeros(n)

  for i in 2:n
    Δz₊ₕ[i-1] = z[i-1] - z[i]
  end
  Δz₊ₕ[end] = z[end] - z₊ₕ[end]
  Δz₊ₕ
end

# Δz = [2, 4, 6, 10]
# z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
# cal_Δz(z₊ₕ) == Δz
# cal_Δz₊ₕ(z, z₊ₕ) == Δz₊ₕ
