#!/usr/bin/env julia
# Regenerate the numerical figures of the paper from stored experiment data.
#   julia --project=. scripts/make_figures.jl [F6 F7] [--data DIR] [--out DIR]
# No solver run is needed; each figure takes seconds.
include(joinpath(@__DIR__, "figure_style.jl"))

include(joinpath(@__DIR__, "cli.jl"))
const OPTIONS = parse_cli(ARGS; output="figures", data=true,
    usage="julia --project=. scripts/make_figures.jl [F6 F7] [--data DIR] [--out DIR] [--quick]")
const RESDIR = OPTIONS.datadir
const FIGDIR = OPTIONS.outdir
loadE(name) = load_results(joinpath(RESDIR, name * ".jld2"))
klog(κ) = round(Int, -log2(κ))

# ===========================================================================
# F6 — the limit and its corrector (flagship)
# ===========================================================================
function figure_corrector()
    E5 = loadE("E5")
    get(E5, "i_star", nothing) == 1.0 || error("Current figures require constant-current E5 data")
    kappas = E5["disk/kappas"]
    fig = Figure(size=(WIDTH, 245))
    dev = E5["disk/dev_norms"]; rem = E5["disk/rem_norms"]
    markers = [:circle, :rect, :utriangle, :diamond]
    ax = Axis(fig[1, 1]; xlabel=L"\kappa", ylabel=L"\|R_\kappa-\tilde{R}_0\|_{H^m}", xscale=log10, yscale=log10)
    pl = [scatterlines!(ax, kappas, dev[m, :]; color=ORDER_COLORS[m], marker=markers[m], linestyle=LINESTYLES[m], linewidth=1.2) for m in 1:4]
    g1 = slope_guide!(ax, kappas, dev[1, end] * 0.35, kappas[end], 1; linestyle=:dot)
    ax = Axis(fig[1, 2]; xlabel=L"\kappa", ylabel=L"\|R_\kappa-\tilde{R}_0-\kappa \tilde{R}_1\|_{H^m}", xscale=log10, yscale=log10)
    for m in 1:4
        scatterlines!(ax, kappas, rem[m, :]; color=ORDER_COLORS[m], marker=markers[m], linestyle=LINESTYLES[m], linewidth=1.2)
    end
    g2 = slope_guide!(ax, kappas, rem[1, end] * 0.35, kappas[end], 2; linestyle=:dash)
    Legend(fig[2, 1:2], [pl..., g1, g2], [L"m=0", L"m=1", L"m=2", L"m=3", "slope 1", "slope 2"]; orientation=:horizontal, tellheight=true, tellwidth=false, nbanks=1)
    rowgap!(fig.layout, 7); colgap!(fig.layout, 25)
    savefig(fig, "fig_corrector"; dir=FIGDIR)
end
# ===========================================================================
# F7 — bulk field and the absence of a boundary layer
# ===========================================================================
function figure_bulk()
    E5 = loadE("E5")
    get(E5, "i_star", nothing) == 1.0 || error("Current figures require constant-current E5 data")
    nr = E5["nr"]
    kappas = E5["disk/kappas"]
    u0 = E5["disk/limit/u0"]
    # One row: scaled differences on four moving domains with a common scale.
    sel = [κ for κ in (1.0, 0.125, 0.015625, 0.00390625) if κ in kappas]
    diffs = [(E5["disk/$κ/phi"] .- u0) ./ κ for κ in sel]
    dmin = minimum(minimum, diffs); dmax = maximum(maximum, diffs)
    # Bright diverging colours keep the zero potential difference near white;
    # truncating the symmetric map to [dmin, dmax] keeps zero white and equal slopes.
    M = max(-dmin, dmax)
    DIFFMAP = cgrad([get(ELECTRO, t) for t in range((1 + dmin / M) / 2, (1 + dmax / M) / 2, length=256)])
    # The colour bar spans the extremes reached; the tiny pad keeps the extreme nodes inside the end bands.
    dmin -= 1e-6 * M; dmax += 1e-6 * M
    fig = Figure(size=(WIDTH, 165))
    tc = nothing
    for (i, κ) in enumerate(sel)
        k = klog(κ)
        ax = Axis(fig[1, i]; aspect=DataAspect(), title=k == 0 ? L"\kappa=1" : L"\kappa=2^{-%$k}")
        tc = field_panel!(ax, E5["disk/$κ/R"], diffs[i], nr; colorrange=(dmin, dmax), colormap=DIFFMAP, contours=false, extend=nothing)
    end
    Colorbar(fig[1, length(sel)+1], tc; label=L"(\phi_\kappa-\tilde{u}_0)/\kappa")
    colgap!(fig.layout, 6)
    savefig(fig, "fig_bulk"; dir=FIGDIR)
end
const FIGURES = Dict("F6" => figure_corrector, "F7" => figure_bulk)

function main()
    todo = isempty(OPTIONS.names) ? ["F6", "F7"] : OPTIONS.names
    for name in todo
        haskey(FIGURES, name) || error("Unknown figure $name; choose F6 or F7")
    end
    for name in todo
        println("Figure $name")
        FIGURES[name]()
    end
end
main()
