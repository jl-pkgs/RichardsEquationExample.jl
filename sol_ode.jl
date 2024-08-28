using DifferentialEquations
using HydroTools
using Plots

# Function to calculate hydraulic conductivity from water content
function van_genuchten_K(θ; param)
  (; θs, θr, Ks, α, n, m) = param
  Se = (θ - θr) / (θs - θr)
  Se = clamp(Se, 0, 1)
  effective_saturation = Se^0.5
  term = (1 - (1 - Se^(1 / m))^m)^2
  return Ks * effective_saturation * term
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

function soil_moisture_transport(du, u, p, t)
  n = length(u)
  (; dz, dt, θ0, z, z₊ₕ, K, ψ, ψ0, param) = p

  K .= van_genuchten_K.(u; param)
  ψ .= van_genuchten_ψ.(u; param)

  # @inbounds 
  @inbounds for i in 2:n-1
    Q_up = cal_Q(i, K, ψ, z)
    Q_down = cal_Q(i + 1, K, ψ, z)
    Δz = z₊ₕ[i-1] - z₊ₕ[i] # 
    du[i] = -(Q_up - Q_down) / Δz
  end

  ## boundary
  i = 1
  K₊ₕ_up = K[1]
  Q_up = -K₊ₕ_up * ((ψ0 - ψ[i]) / (0 - z[i]/2) + 1)
  Q_down = cal_Q(i + 1, K, ψ, z)
  Δz = 0 - z₊ₕ[i]
  du[i] = -(Q_up - Q_down) / Δz

  i = n
  Q_up = cal_Q(i, K, ψ, z)
  Q_down = -K[i]
  Δz = (z[i-1] - z[i]) / 2
  du[i] = -(Q_up - Q_down) / Δz
end


# van Genuchten parameters
θs = 0.287  # Saturated water content
θr = 0.075  # Residual water content
α = 0.027  # 1/cm
n = 3.96
m = 1
Ks = 34 / 3600
param = (; θs, θr, Ks, α, n, m)

θ0 = 0.267
ψ0 = van_genuchten_ψ(θ0; param)

n = 150
dz = ones(n)
z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)

dt = 5
ntim = 0.8 * 3600 / dt

# n = 100
u0 = fill(0.1, n) |> collect # Example initial soil moisture profile
tspan = (0.0, 0.8 * 3600)  # Time span for the simulation
p = (; dz, dt, θ0, ψ0, z, z₊ₕ, K=zeros(n), ψ=zeros(n), param=param)


prob = ODEProblem(soil_moisture_transport, u0, tspan, p);
@time sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200);
# plot(sol)

begin
  gr(; framestyle=:box)
  _u = sol.u[end]
  ψ = van_genuchten_ψ.(_u; param)
  p1 = plot(sol.u[end], z; xlabel="θ", ylabel="z", label="θ")
  # p2 = plot(ψ, z; xlabel="ψ", ylabel="z", label="ψ")
  # plot(p1, p2)
end

begin
  fig = plot()
  for i in 5:length(sol.u)-1
    _u = sol.u[i+1]
    _t = sol.t[i+1]
    plot!(fig, _u, z; label="t = $_t")
  end
  # _u = cat(sol.u..., dims=2)
  fig
  # plot(sol)
end
