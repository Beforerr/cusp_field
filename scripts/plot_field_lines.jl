using LinearAlgebra: norm
using CairoMakie, MakieBake

fieldline_slope(c, z::Number) = c(z)[1:2] ./ c(z)[3]

function cumulative_trapezoid(xs, ys)
    out = zeros(Float64, length(xs))
    for i in 2:length(xs)
        out[i] = out[i - 1] + (xs[i] - xs[i - 1]) * (ys[i] + ys[i - 1]) / 2
    end
    return out
end

function fieldline_profile(c; zmax = 5.0, n = 401)
    nside = (n + 1) ÷ 2
    zpos = collect(range(0.0, zmax, length = nside))
    zneg = collect(range(0.0, -zmax, length = nside))
    z = vcat(reverse(zneg[2:end]), zpos)

    slopes_pos = fieldline_slope.(Ref(c), zpos)
    slopes_neg = fieldline_slope.(Ref(c), zneg)

    xpos = cumulative_trapezoid(zpos, first.(slopes_pos))
    ypos = cumulative_trapezoid(zpos, last.(slopes_pos))
    xneg = cumulative_trapezoid(zneg, first.(slopes_neg))
    yneg = cumulative_trapezoid(zneg, last.(slopes_neg))

    κ = norm.(curvature.(c, z))
    x = vcat(reverse(xneg[2:end]), xpos)
    y = vcat(reverse(yneg[2:end]), ypos)
    return (; z, x, y, κ)
end

fieldline1 = @lift fieldline_profile(BConfig(; B0 = $params.B0, a = $params.a, By = $params.By); zmax = 5.0)
fieldline2 = @lift fieldline_profile(BConfig2(; B0 = $params.B0, a = $params.a, By = $params.By); zmax = 5.0)

z1 = lift(data -> data.z, fieldline1)
x1 = lift(data -> data.x, fieldline1)
y1 = lift(data -> data.y, fieldline1)
κ1 = lift(data -> data.κ, fieldline1)

z2 = lift(data -> data.z, fieldline2)
x2 = lift(data -> data.x, fieldline2)
y2 = lift(data -> data.y, fieldline2)
κ2 = lift(data -> data.κ, fieldline2)

fieldline_fig = Figure(; size = (900, 600))

ax11 = Axis(fieldline_fig[1, 1], xlabel = "z/L", ylabel = "position", title = "Field line through origin")
lines!(ax11, z1, x1, color = Makie.wong_colors()[1], label = "x(z)")
lines!(ax11, z1, y1, color = Makie.wong_colors()[2], linestyle = :dash, label = "y(z)")
scatter!(ax11, [0.0], [0.0], color = :black, markersize = 8)
axislegend(ax11, position = :rb)

ax12 = Axis(fieldline_fig[1, 2], xlabel = "z/L", ylabel = "position", title = "Field line through origin")
lines!(ax12, z2, x2, color = Makie.wong_colors()[1], label = "x(z)")
lines!(ax12, z2, y2, color = Makie.wong_colors()[2], linestyle = :dash, label = "y(z)")
scatter!(ax12, [0.0], [0.0], color = :black, markersize = 8)
axislegend(ax12, position = :rb)

ax21 = Axis(fieldline_fig[2, 1], xlabel = "z/L", ylabel = "κ(z)", title = "Curvature profile")
lines!(ax21, z1, κ1, color = Makie.wong_colors()[3])

ax22 = Axis(fieldline_fig[2, 2], xlabel = "z/L", ylabel = "κ(z)", title = "Curvature profile")
lines!(ax22, z2, κ2, color = Makie.wong_colors()[3])

Label(fieldline_fig[0, 1], "Bₓ = B₀ (z/L) / (1 + a (z/L)²)²"; tellwidth = false)
Label(fieldline_fig[0, 2], "Bₓ = B₀ tanh(z/L) / (1 + a (z/L)²)²"; tellwidth = false)

fieldline_fig

bake_html(
    params => (
        (@o _.B0) => [1.0, 10.0],
        (@o _.a) => [0.0, 0.05, 0.1, 0.2, 0.5],
        (@o _.By) => [0.0, 0.5, 1.0],
    );
    blocks = [fieldline_fig],
    outdir = "./web/fieldlines"
)
