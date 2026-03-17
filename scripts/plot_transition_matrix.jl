using Observables
using CurrentSheetTestParticle
using PitchAngleAnalysis
using CairoMakie

include("fields.jl")
includet("fields.jl")

params = Observable((v = 1.0, B0 = 10.0, a = 0.1, By = 0.0))

fig = Figure(; size = (800, 800))
Label(fig[0, 1], "Bₓ = B₀ (z/L) / (1 + a (z/L)²)²"; tellwidth = false)
Label(fig[0, 2], "Bₓ = B₀ tanh(z/L) / (1 + a (z/L)²)²"; tellwidth = false)
@lift begin
    save_everystep = false
    tspan = (0.0, 1024.0)
    init_kwargs = (; Nw = 64, Nϕ = 32)
    cfg1 = BConfig(; B0 = $params.B0, a = $params.a, By = $params.By)
    cfg2 = BConfig2(; B0 = $params.B0, a = $params.a, By = $params.By)
    p1 = ProblemParamsBase(; v = $params.v, B = cfg1, init_kwargs)
    p2 = ProblemParamsBase(; v = $params.v, B = cfg2, init_kwargs)
    sols1 = solve_params(p1; save_everystep, tspan, load = true)[1]
    sols2 = solve_params(p2; save_everystep, tspan, load = true)[1]
    plot_transition_matrix!(fig[1, 1], sols1)
    plot_transition_matrix!(fig[1, 2], sols2)
    plot_mm_conservation!(fig[2, 1], sols1)
    plot_mm_conservation!(fig[2, 2], sols2)
end
fig

# %%
using MakieBake

bake_html(
    params => (
        (@o _.v) => [0.1, 1.0, 10.0],
        (@o _.B0) => [1.0, 10.0],
        (@o _.a) => [0.0, 0.01, 0.05, 0.1, 0.2, 0.5],
        (@o _.By) => [-1.0, 0.0, 1.0],
    );
    blocks = [fig],
    outdir = "./web/transition_matrices"
)
