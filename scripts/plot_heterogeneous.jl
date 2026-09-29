#!/usr/bin/env julia
# Plot only the standalone heterogeneous record; no baseline figures are changed.
include(joinpath(@__DIR__, "figure_style.jl"))

root = normpath(joinpath(@__DIR__, ".."))
datadir = isempty(ARGS) ? joinpath(root, "data", "heterogeneous") : abspath(ARGS[1])
outdir = length(ARGS) < 2 ? joinpath(root, "outputs", "figures") : abspath(ARGS[2])
mkpath(outdir)
d = load_results(joinpath(datadir, "heterogeneous.jld2"))
long = load_results(joinpath(datadir, "long_time.jld2"))
long["limit"]["T"] == 1.0 || error("The long-time record lacks the limiting interface at t=1")
initial_radius(th) = 1 + 0.12cos(3th) + 0.06sin(2th) + 0.04cos(5th + 0.3)
material(x, y) = 1 + 0.25x + 0.10y^2

# Drawn at the paper's scale, like the other figures, so that the colour bar
# and fonts match them at \textwidth.
fig = Figure(size=(WIDTH, 245), figure_padding=(4, 8, 4, 4))
left = fig[1, 1] = GridLayout()
right = fig[1, 2] = GridLayout()
axg = Axis(left[1, 1]; aspect=DataAspect(), title="(a) Heterogeneous recession",
    xlabel=L"x_1", ylabel=L"x_2", xticks=[-1, 0, 1], yticks=[-1, 0, 1])
xx = range(-1.30, 1.30, length=601)
yy = range(-1.30, 1.30, length=601)
field = [hypot(x, y) <= initial_radius(atan(y, x)) ? material(x, y) : NaN for x in xx, y in yy]
# i0 has no critical point, so its extremes over the initial solid are attained on
# the interface; they are also the extremes of the space-time figure.
lo, hi = extrema(material.(d["R_initial"] .* cos.(d["theta"]), d["R_initial"] .* sin.(d["theta"])))
hm = heatmap!(axg, xx, yy, field; colormap=banded_diverging(lo, hi; center=1.0),
    colorrange=(lo, hi), rasterize=3)
th = vcat(d["theta"], 2π)
closed(v) = vcat(v, v[1])
rinit = closed(d["R_initial"])
initplot = lines!(axg, rinit .* cos.(th), rinit .* sin.(th);
    color=TEX_BLUE, linewidth=1.2, linestyle=:dash)
rad = closed(long["limit"]["R"])
finalplot = lines!(axg, rad .* cos.(th), rad .* sin.(th); color=TEX_RED, linewidth=1.7)
limits!(axg, -1.30, 1.30, -1.30, 1.30)
Colorbar(left[1, 2], hm; label=L"i_0(x_1,x_2)", ticks=0.8:0.1:1.2)
Legend(left[2, 1], [initplot, finalplot], [L"t=0", L"t=1"];
    orientation=:horizontal, tellwidth=false, tellheight=true)
colsize!(left, 1, Aspect(1, 1.0))
colgap!(left, 6)

axe = Axis(right[1, 1]; title="(b) Conductivity convergence", xlabel=L"\kappa",
    ylabel="Error at T = 0.5", xscale=log10, yscale=log10,
    xticks=([2.0^-6, 2.0^-4, 2.0^-2], [L"2^{-6}", L"2^{-4}", L"2^{-2}"]),
    yticks=([1e-5,1e-4,1e-3,1e-2], [L"10^{-5}",L"10^{-4}",L"10^{-3}",L"10^{-2}"]))
k = d["kappas"]
plots = []
for (fieldname, color, marker) in (("deviation", TEX_BLUE, :circle), ("remainder", TEX_RED, :rect))
    for m in 1:2
        push!(plots, scatterlines!(axe, k, d["errors/$fieldname"][m,:];
            color, marker, linestyle=m==1 ? :solid : :dash, markersize=5, linewidth=1.2))
    end
end
xguide = [k[end], k[end-2]]
lines!(axe, xguide, 0.003 .* (xguide ./ xguide[1]); color=GREY, linewidth=0.8, linestyle=:dot)
text!(axe, 0.036, 0.010; text="slope 1", fontsize=10.5, color=GREY, rotation=0.20)
lines!(axe, xguide, 0.6e-5 .* (xguide ./ xguide[1]).^2; color=GREY, linewidth=0.8, linestyle=:dot)
text!(axe, 0.035, 1.3e-5; text="slope 2", fontsize=10.5, color=GREY, rotation=0.38)
limits!(axe, 0.012, 0.32, 4e-6, 0.035)
Legend(right[2, 1], plots,
    [L"E_0", L"E_1", L"E_0^{\mathrm{corr}}", L"E_1^{\mathrm{corr}}"]; orientation=:horizontal,
    nbanks=2, tellwidth=false, tellheight=true, colgap=10)
colgap!(fig.layout, 22)
rowgap!(left, 5)
rowgap!(right, 5)

pngdir = normpath(outdir) == joinpath(root, "figures", "paper") ?
    joinpath(root, "figures", "previews") : outdir
mkpath(pngdir)
save(joinpath(outdir, "fig_heterogeneous.pdf"), fig; pt_per_unit=0.75)
save(joinpath(pngdir, "fig_heterogeneous.png"), fig; px_per_unit=2)
println("Wrote heterogeneous PDF to ", outdir, " and PNG to ", pngdir)
