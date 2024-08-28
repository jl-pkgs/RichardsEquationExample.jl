using DifferentialEquations
using HydroTools
using Plots
using Parameters
include("Soil_depth.jl")


@with_kw mutable struct Soil{FT}
  n::Int = 10
  z::Vector{FT} = zeros(FT, n)
  z₊ₕ::Vector{FT} = zeros(FT, n)
  Δz::Vector{FT} = zeros(FT, n)
  Δz₊ₕ::Vector{FT} = zeros(FT, n)

  u::Vector{FT} = ones(FT, n) .* 0.1# θ
  Q::Vector{FT} = zeros(FT, n)
  K::Vector{FT} = zeros(FT, n)
  ψ::Vector{FT} = zeros(FT, n)
  ψ0::FT = FT(0.0)

  param::NamedTuple = (; θs=0.287, θr=0.075, Ks=34 / 3600, α=0.027, n=3.96, m=1)
end

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
