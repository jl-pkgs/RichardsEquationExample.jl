using DifferentialEquations
using HydroTools
using Plots
using Parameters


# Function to calculate hydraulic conductivity from water content
function van_genuchten_K(θ; param)
  (; θs, θr, Ks, α, n, m) = param
  Se = (θ - θr) / (θs - θr)
  # Se = clamp(Se, 0, 1)
  # effective_saturation = Se^0.5
  # term = (1 - (1 - Se^(1 / m))^m)^2
  # return Ks * effective_saturation * term

  if Se <= 1
    # Special case for:
    # - `soil_texture = 1`: Haverkamp et al. (1977) sand
    # - `soil_texture = 2`: Yolo light clay
    ψ = van_genuchten_ψ(θ; param)
    Ks * 1.175e6 / (1.175e6 + abs(ψ)^4.74) # Haverkamp et al. (1977) sand
  # Ks * 124.6 / (124.6 + abs(ψ)^1.77)   # Yolo light clay
  else
    Ks
  end
end

# Function to calculate pressure head psi from water content
function van_genuchten_ψ(θ; param)
  (; θs, θr, α, n, m) = param
  if θ <= θr
    return -Inf  # Return a very high positive number indicating very dry conditions
  elseif θ >= θs
    return 0  # Saturated condition, psi is zero
  else
    return -(1 / α) * (((θs - θr) / (θ - θr))^(1 / m) - 1)^(1 / n)
  end
end

# Define the function representing the system of ODEs for soil moisture transport
function cal_Q(i::Int, K::Vector{FT}, ψ::Vector{FT}, z::Vector{FT}) where {FT}
  K₊ₕ = (K[i] + K[i-1]) / 2
  Δz₊ₕ = z[i-1] - z[i]
  Q = -K₊ₕ * ((ψ[i-1] - ψ[i]) / Δz₊ₕ + 1)
  Q
end



function RichardsEquation_ode(du, u, p::Soil, t)
  (; n, z, z₊ₕ, K, ψ, ψ0, param) = p

  # @inbounds for i in 1:n
  #   K[i] = van_genuchten_K(u[i]; param)
  #   ψ[i] = van_genuchten_ψ(u[i]; param)
  # end
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
  Δz = (z[i-1] - z[i]) / 2
  du[i] = -(Q_up - Q_down) / Δz
end



@with_kw mutable struct Soil{FT}
  n::Int = 10
  z::Vector{FT} = zeros(FT, n)
  z₊ₕ::Vector{FT} = zeros(FT, n)
  Δ::Vector{FT} = zeros(FT, n)
  Δz₊ₕ::Vector{FT} = zeros(FT, n)

  u::Vector{FT} = zeros(FT, n) # θ
  Q::Vector{FT} = zeros(FT, n)
  K::Vector{FT} = zeros(FT, n)
  ψ::Vector{FT} = zeros(FT, n)
  ψ0::FT = FT(0.0)

  param::NamedTuple = (; θs=0.287, θr=0.075, Ks=34 / 3600, α=0.027, n=3.96, m=1)
end

# van Genuchten parameters
# θs = 0.287  # Saturated water content
# θr = 0.075  # Residual water content
# α = 0.027  # 1/cm
param = (; θs=0.287, θr=0.075, Ks=34 / 3600, α=0.027, n=3.96, m=1)
θ0 = 0.267
ψ0 = van_genuchten_ψ(θ0; param)

n = 150
dz = ones(n) # Δz₊ₕ


p = Soil{Float64}(; n=150, ψ0, z, z₊ₕ, Δz₊ₕ)
# (; z, z₊ₕ, K, ψ, ψ0, param) = p

dt = 5
ntim = 0.8 * 3600 / dt
