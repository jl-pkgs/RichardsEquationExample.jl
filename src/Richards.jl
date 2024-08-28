include("Soil.jl")


function RichardsEquation(du, u, p::Soil, t)
  (; n, Q, z, Δz, K, ψ, ψ0, param) = p

  @. K = van_genuchten_K(u; param)
  @. ψ = van_genuchten_ψ(u; param)

  Q0 = -K[1] * ((ψ0 - ψ[1]) / (0 - z[1]) + 1)
  @inbounds for i in 1:n-1
    K₊ₕ = (K[i] + K[i+1]) / 2
    Δz₊ₕ = z[i] - z[i+1]
    Q[i] = -K₊ₕ * ((ψ[i] - ψ[i+1]) / Δz₊ₕ + 1)
  end

  Q[n] = -K[n]
  du[1] = -(Q0 - Q[1]) / (0 - z[1] * 2)

  @inbounds for i in 2:n
    du[i] = -(Q[i-1] - Q[i]) / Δz[i]
  end
end


# Q[i-1]
function cal_Q(i::Int, K::Vector{FT}, ψ::Vector{FT}, z::Vector{FT}) where {FT}
  K₊ₕ = (K[i] + K[i-1]) / 2
  Δz₊ₕ = z[i-1] - z[i]
  Q = -K₊ₕ * ((ψ[i-1] - ψ[i]) / Δz₊ₕ + 1)
  Q
end

function RichardsEquation_low(du, u, p::Soil, t)
  (; n, z, z₊ₕ, K, ψ, ψ0, param) = p

  @. K = van_genuchten_K(u; param)
  @. ψ = van_genuchten_ψ(u; param)

  # @inbounds 
  @inbounds for i in 2:n-1
    Q_up = cal_Q(i, K, ψ, z) # 多计算了1次
    Q_down = cal_Q(i + 1, K, ψ, z)
    Δz = z₊ₕ[i-1] - z₊ₕ[i] # 
    du[i] = -(Q_up - Q_down) / Δz
  end

  ## boundary
  i = 1
  K₊ₕ_up = K[1]
  Q_up = -K₊ₕ_up * ((ψ0 - ψ[i]) / (0 - z[i]) + 1)
  Q_down = cal_Q(i + 1, K, ψ, z)
  Δz = 0 - z₊ₕ[i]
  du[i] = -(Q_up - Q_down) / Δz

  i = n
  Q_up = cal_Q(i, K, ψ, z)
  Q_down = -K[i]
  Δz = (z[i-1] - z[i])
  du[i] = -(Q_up - Q_down) / Δz
end
