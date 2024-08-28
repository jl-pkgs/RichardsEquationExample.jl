includet("src/Richards.jl")
using JLD2, Test
bonan = load("data/output_bonan.jld2")


param = (; θs=0.287, θr=0.075, Ks=34 / 3600, α=0.027, n=3.96, m=1)
θ0 = 0.267
ψ0 = van_genuchten_ψ(θ0; param)

n = 150
Δz = ones(n) # Δz₊ₕ
z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)


function solve_ode()
  u0 = fill(0.1, n) |> collect # Example initial soil moisture profile
  tspan = (0.0, 0.8 * 3600)  # Time span for the simulation
  p = Soil{Float64}(; n=150, ψ0, z, z₊ₕ, Δz, Δz₊ₕ)
  prob = ODEProblem(RichardsEquation, u0, tspan, p)
  solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
end


begin
  @time sol = solve_ode()

  @testset "Richards" begin
    max_error = maximum(abs, sol.u[end] .- bonan["θ"])
    @test max_error < 1e-3
  end
end

@profview for i = 1:10
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200);
end

# begin
#   gr(; framestyle=:box)
#   _u = sol.u[end]
#   ψ = van_genuchten_ψ.(_u; param)
#   p1 = plot(sol.u[end], z; xlabel="θ", ylabel="z", label="θ")
#   # p2 = plot(ψ, z; xlabel="ψ", ylabel="z", label="ψ")
#   # plot(p1, p2)
#   Plots.savefig("soil_moisture_profile.png")
# end

# begin
#   fig = plot()
#   for i in 5:length(sol.u)-1
#     _u = sol.u[i+1]
#     _t = sol.t[i+1]
#     plot!(fig, _u, z; label="t = $_t")
#   end
#   # _u = cat(sol.u..., dims=2)
#   fig
#   # plot(sol)
# end
