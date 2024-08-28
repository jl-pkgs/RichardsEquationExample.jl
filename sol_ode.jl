include("main_ode.jl")

# n = 100
u0 = fill(0.1, n) |> collect # Example initial soil moisture profile
tspan = (0.0, 0.8 * 3600)  # Time span for the simulation


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
