include("src/Richards.jl")
# Define the function representing the system of ODEs for soil moisture transport

# van Genuchten parameters
# θs = 0.287  # Saturated water content
# θr = 0.075  # Residual water content
# α = 0.027  # 1/cm
param = (; θs=0.287, θr=0.075, Ks=34 / 3600, α=0.027, n=3.96, m=1)
θ0 = 0.267
ψ0 = van_genuchten_ψ(θ0; param)

n = 150
Δz = ones(n) # Δz₊ₕ
z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)



begin
  using Ipaper
  set_seed(1)
  θs = 0.287
  θr = 0.075
  u1 = rand(n) .* (θs - θr) .+ θr
  u2 = deepcopy(u1)

  p1 = Soil{Float64}(; n=150, ψ0, z, z₊ₕ, Δz, Δz₊ₕ, u=u1)
  du1 = zeros(n)
  RichardsEquation_ode_V1(du1, p1.u, p1, 0.0)

  p2 = Soil{Float64}(; n=150, ψ0, z, z₊ₕ, Δz, Δz₊ₕ, u=u2)
  du2 = zeros(n)
  RichardsEquation_ode_V2(du2, p2.u, p2, 0.0)

  du1 == du2
end
# (; z, z₊ₕ, K, ψ, ψ0, param) = p

