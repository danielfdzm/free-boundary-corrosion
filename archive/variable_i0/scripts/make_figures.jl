#!/usr/bin/env julia
# Regenerate the numerical figures of the paper from stored experiment data.
#   julia --project=. scripts/make_figures.jl [F1 F4 F6 F7] [--data DIR] [--out DIR]
# No solver run is needed; each figure takes seconds.
include(joinpath(@__DIR__, "figure_style.jl"))

include(joinpath(@__DIR__, "cli.jl"))
const OPTIONS = parse_cli(ARGS; output="figures", data=true,
    usage="julia --project=. scripts/make_figures.jl [F1 F4 F6 F7] [--data DIR] [--out DIR] [--quick]")
const RESDIR = OPTIONS.datadir
const FIGDIR = OPTIONS.outdir
loadE(name) = load_results(joinpath(RESDIR, name * ".jld2"))
klog(κ) = round(Int, -log2(κ))

# ===========================================================================
# F1 — hero: the model in one picture (flower under M2, kappa = 1/2)
# ===========================================================================
function figure_hero()
    E2 = loadE("E2")
    ns = E2["lifetime/flower/nsnapshots"]
    nr = E2["lifetime/flower/nr"]
    snaps = [(E2["lifetime/flower/snapshot/$i/t"], E2["lifetime/flower/snapshot/$i/R"], E2["lifetime/flower/snapshot/$i/phi"]) for i in 1:ns]
    vmax = ceil(maximum(maximum(abs, s[3]) for s in snaps) * 20) / 20
    fig = Figure(size=(WIDTH, 210))
    tc = nothing
    theta = theta_grid(length(snaps[1][2]))
    thc = vcat(theta, theta[1])
    previous_styles = (:dash, :dot)
    for (i, (t, R, phi)) in enumerate(snaps)
        ax = Axis(fig[1, i]; aspect=DataAspect(), title=L"t = %$(round(t, digits=2))")
        tc = field_panel!(ax, R, phi, nr; colorrange=(-vmax, vmax), streamlines=true, ncontours=8)
        # interfaces at the previous snapshot times, as discontinuous curves
        for j in 1:i-1
            Rp = vcat(snaps[j][2], snaps[j][2][1])
            lines!(ax, Rp .* cos.(thc), Rp .* sin.(thc); color=SOLIDDARK, linestyle=previous_styles[j], linewidth=1.1)
        end
        limits!(ax, -2.05, 2.05, -2.05, 2.05)
    end
    Colorbar(fig[1, ns+1], tc; label=L"\phi_\kappa", ticks=WilkinsonTicks(5))
    colgap!(fig.layout, 4)
    savefig(fig, "fig_hero"; dir=FIGDIR)
end

# ===========================================================================
# F4 — transfer factor and sharpness (flagship)
# ===========================================================================
function figure_transfer()
    E3 = loadE("E3")
    kappas = E3["kappas"]; ks = E3["ks"]; ns = E3["ns"]; a = E3["a"]
    id = 1                                   # delta = 0.05: the regime of the lemma
    tau = E3["tau0"][:, :, id]; mu = E3["mu"][:, :, id]
    lam = E3["lambda"]
    lin = [κ * l / (a + κ * l) for κ in kappas, l in lam]
    ln2n = log2.(ns); ln2k = float.(ks)
    lamf(n) = n * (1 - (1 / 2)^(2n)) / (1 + (1 / 2)^(2n))          # lambda_n on the unit circle in B_2
    fig = Figure(size=(WIDTH, 420))
    ax3 = Axis3(fig[1, 1]; azimuth=-1.05, elevation=0.33, aspect=(1, 1, 0.7), xlabel=L"\log_2 n", ylabel=L"\log_2(1/\kappa)", zlabel=L"\tau",
        xticks=(ln2n, string.(Int.(ln2n))), yticks=(0:2:ks[end], string.(0:2:ks[end])), zticks=[0, 0.5, 1], protrusions=(50, 20, 10, 25),
        xlabeloffset=30, ylabeloffset=35, zlabeloffset=35)
    surface!(ax3, ln2n, ln2k, permutedims(tau); colormap=:viridis, colorrange=(0, 1), rasterize=4)
    wireframe!(ax3, ln2n, ln2k, permutedims(lin); color=(:black, 0.4), linewidth=0.4)
    contour!(ax3, ln2n, ln2k, permutedims(tau); levels=[0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8], colormap=:viridis, colorrange=(0, 1), linewidth=0.7,
        transformation=(:xy, -0.03))
    panel_label!(fig, (1, 1), "(a)")
    # (b) ratio tau/mu on the corridor scale, centred at 1
    ratio = tau ./ mu
    cmap = cgrad([DEEP, WATER, NEUTRAL, SOLIDBASE, SOLIDDARK], [0.0, 0.3, 0.6, 0.8, 1.0])
    ax = Axis(fig[1, 2]; xlabel=L"\log_2 n", ylabel=L"\log_2(1/\kappa)", xticks=(ln2n, string.(Int.(ln2n))), yticks=(0:2:ks[end], string.(0:2:ks[end])))
    hm = heatmap!(ax, ln2n, ln2k, permutedims(ratio); colormap=cmap, colorrange=(0.25, 1.5))
    xlo, xhi, ylo, yhi = -0.5, log2(ns[end]) + 0.5, -0.5, ks[end] + 0.5
    limits!(ax, xlo, xhi, ylo, yhi)
    nn = range(xlo, xhi, length=400)
    kk = [-log2(a / lamf(2.0^x)) for x in nn]
    keep = (kk .>= ylo) .& (kk .<= yhi)
    lines!(ax, nn[keep], kk[keep]; color=:white, linewidth=1.6)
    xl = 4.3
    text!(ax, xl, -log2(a / lamf(2.0^xl)) - 0.9; text=L"\kappa\lambda_n = a", color=:white, align=(:left, :top), fontsize=FS_SMALL)
    Colorbar(fig[1, 3], hm; label=L"\tau/\mu_{\kappa,n}", ticks=([0.25, 0.5, 0.75, 1.0, 1.25, 1.5], ["1/4", "0.5", "0.75", "1", "1.25", "3/2"]))
    panel_label!(fig, (1, 2), "(b)")
    # (c) tau vs kappa at fixed n
    ax = Axis(fig[2, 1:3]; xlabel=L"\kappa", ylabel=L"\tau(\kappa,n)", xscale=log10, yscale=log10)
    sel = [n for n in (1, 8, 64) if n in ns]
    for (i, n) in enumerate(sel)
        jn = findfirst(==(n), ns)
        lines!(ax, kappas, lin[:, jn]; color=(GREY, 0.9), linestyle=:dash, linewidth=0.9, label=i == 1 ? L"\kappa\lambda_n/(a+\kappa\lambda_n)" : nothing)
        scatterlines!(ax, kappas, tau[:, jn]; color=kappa_color(3i; kmax=9), marker=[:circle, :rect, :utriangle][i], label=L"n=%$n", linewidth=0.8)
    end
    slope_guide!(ax, kappas, tau[end, 1] * 0.45, kappas[end], 1; label="slope 1", linestyle=:dot)
    Legend(fig[3, 1:3], ax; orientation=:horizontal, tellheight=true, tellwidth=false, nbanks=1)
    panel_label!(fig, (2, 1), "(c)")
    rowsize!(fig.layout, 1, Relative(0.6))
    rowgap!(fig.layout, 6)
    savefig(fig, "fig_transfer"; dir=FIGDIR)
end
# ===========================================================================
# F6 — the limit and its corrector (flagship)
# ===========================================================================
function figure_corrector()
    E5 = loadE("E5")
    theta = E5["theta"]; nt = length(theta)
    kappas = E5["disk/kappas"]
    κmin = kappas[end]; kmax = klog(κmin)
    tl = E5["disk/limit/times"]
    R0h = E5["disk/limit/Rhist"]; Wh = E5["disk/limit/Whist"]
    Rk = E5["disk/$κmin/Rhist"]; tk = E5["disk/$κmin/times"]
    idx = [findfirst(t -> isapprox(t, s; atol=1e-9), tk) for s in tl]
    Z = (Rk[:, idx] .- R0h) ./ κmin
    step = max(1, nt ÷ 128)
    th_c = vcat(theta[1:step:end], 2π); Zc = vcat(Z[1:step:end, :], Z[1:1, :]); Wc = vcat(Wh[1:step:end, :], Wh[1:1, :])
    vmax = maximum(abs, Zc)
    fig = Figure(size=(WIDTH, 450))
    ax3 = Axis3(fig[1, 1:2]; azimuth=-1.0, elevation=0.36, aspect=(1.3, 1, 0.5), xlabel=L"\theta", ylabel=L"t", zlabel=L"(R_\kappa-\tilde{R}_0)/\kappa",
        xticks=PI_TICKS, protrusions=(40, 10, 10, 20), zlabeloffset=55, ylabeloffset=35, xlabeloffset=22)
    surface!(ax3, th_c, tl, Zc; colormap=ELECTRO, colorrange=(-vmax, vmax), rasterize=4)
    wireframe!(ax3, th_c[1:2:end], tl[1:2:end], Wc[1:2:end, 1:2:end]; color=(:black, 0.3), linewidth=0.35)
    contour!(ax3, th_c, tl, Zc; levels=9, colormap=ELECTRO, colorrange=(-vmax, vmax), linewidth=0.6, transformation=(:xy, -1.3 * vmax))
    Colorbar(fig[1, 3]; limits=(-vmax, vmax), colormap=ELECTRO, label=L"(R_\kappa-\tilde{R}_0)/\kappa,\ \kappa=2^{-%$kmax}", height=Relative(0.65))
    panel_label!(fig, (1, 1), "(a)")
    dev = E5["disk/dev_norms"]; rem = E5["disk/rem_norms"]
    markers = [:circle, :rect, :utriangle, :diamond]
    ax = Axis(fig[2, 1]; xlabel=L"\kappa", ylabel=L"\|R_\kappa-\tilde{R}_0\|_{H^m}\ \text{at}\ T", xscale=log10, yscale=log10)
    pl = [scatterlines!(ax, kappas, dev[m, :]; color=kappa_color(2(m - 1); kmax=8), marker=markers[m], linestyle=LINESTYLES[m], linewidth=0.8) for m in 1:4]
    g1 = slope_guide!(ax, kappas, dev[1, end] * 0.35, kappas[end], 1; linestyle=:dot)
    panel_label!(fig, (2, 1), "(b)")
    ax = Axis(fig[2, 2:3]; xlabel=L"\kappa", ylabel=L"\|R_\kappa-\tilde{R}_0-\kappa \tilde{R}_1\|_{H^m}\ \text{at}\ T", xscale=log10, yscale=log10)
    for m in 1:4
        scatterlines!(ax, kappas, rem[m, :]; color=kappa_color(2(m - 1); kmax=8), marker=markers[m], linestyle=LINESTYLES[m], linewidth=0.8)
    end
    g2 = slope_guide!(ax, kappas, rem[1, end] * 0.35, kappas[end], 2; linestyle=:dash)
    panel_label!(fig, (2, 2), "(c)")
    Legend(fig[3, 1:3], [pl..., g1, g2], [L"m=0", L"m=1", L"m=2", L"m=3", "slope 1", "slope 2"]; orientation=:horizontal, tellheight=true, tellwidth=false, nbanks=1)
    rowsize!(fig.layout, 1, Relative(0.5))
    rowgap!(fig.layout, 8); colgap!(fig.layout, 10)
    savefig(fig, "fig_corrector"; dir=FIGDIR)
end
# ===========================================================================
# F7 — bulk field and the absence of a boundary layer
# ===========================================================================
function figure_bulk()
    E5 = loadE("E5")
    nt = E5["nt"]; nr = E5["nr"]
    kappas = E5["disk/kappas"]
    κ1 = kappas[1]
    phi1 = E5["disk/$κ1/phi"]; R1 = E5["disk/$κ1/R"]
    u0 = E5["disk/limit/u0"]; R0 = E5["disk/limit/R"]
    vmax = ceil(max(maximum(abs, phi1), maximum(abs, u0)) * 20) / 20
    fig = Figure(size=(WIDTH, 590))
    ga = fig[1, 1] = GridLayout()
    ax = Axis(ga[1, 1]; aspect=DataAspect(), title=L"\phi_\kappa\ \text{on}\ \Omega_\kappa(T),\ \kappa=1")
    field_panel!(ax, R1, phi1, nr; colorrange=(-vmax, vmax))
    Label(ga[1, 1, TopLeft()], "(a)"; font=:bold, fontsize=FS, padding=(0, 6, 4, 0), halign=:right)
    ax = Axis(ga[1, 2]; aspect=DataAspect(), title=L"\tilde{u}_0 = \mathcal{P}_{\Gamma_{*,0}(T)}\phi_{\mathrm{eq}}\ \text{on}\ \Omega_0(T)")
    tc = field_panel!(ax, R0, u0, nr; colorrange=(-vmax, vmax))
    Label(ga[1, 2, TopLeft()], "(b)"; font=:bold, fontsize=FS, padding=(0, 6, 4, 0), halign=:right)
    Colorbar(ga[1, 3], tc; label=L"\phi")
    colgap!(ga, 6)
    sel = [κ for κ in (1.0, 0.125, 0.015625) if κ in kappas]
    diffs = [(E5["disk/$κ/phi"] .- u0) ./ κ for κ in sel]
    dmax = ceil(maximum(maximum(abs, d) for d in diffs[2:end]) * 20) / 20
    # Symmetric slate-blue / ivory / muted-copper map, matching the geometry figure.
    DIFFMAP = cgrad([colorant"#526D82", colorant"#A9BDC8", colorant"#F7F5F0",
        colorant"#D5B5A6", colorant"#9A6048"], [0.0, 0.25, 0.5, 0.75, 1.0])
    gc = fig[2, 1] = GridLayout()
    tc = nothing
    for (i, κ) in enumerate(sel)
        ax = Axis(gc[1, i]; aspect=DataAspect(), title=L"\kappa=2^{-%$(klog(κ))}")
        tc = field_panel!(ax, E5["disk/$κ/R"], diffs[i], nr; colorrange=(-dmax, dmax), colormap=DIFFMAP, contours=false)
        i == 1 && Label(gc[1, 1, TopLeft()], "(c)"; font=:bold, fontsize=FS, padding=(0, 6, 4, 0), halign=:right)
    end
    Colorbar(gc[1, length(sel)+1], tc; label=L"(\phi_\kappa-\tilde{u}_0)/\kappa")
    colgap!(gc, 6)
    ax = Axis(fig[3, 1]; xlabel=L"\text{distance to } \Gamma_{*,\kappa}(T)", ylabel=L"(\phi_\kappa-\tilde{u}_0)/\kappa")
    for κ in kappas
        k = klog(κ)
        (k in (0, 3, 6, 8)) || continue
        prof = E5["disk/$κ/profiles"]; dist = E5["disk/$κ/profile_distances"]
        for (ir, ls) in zip((1, 2), (:solid, :dash))
            lines!(ax, dist[:, ir], prof[:, ir]; color=kappa_color(k; kmax=8), linestyle=ls, linewidth=1.0, label=ir == 1 ? L"\kappa=2^{-%$k}" : nothing)
        end
    end
    text!(ax, 0.98, 0.95; text=L"\text{solid: ray } \theta=0,\ \text{dashed: ray } \theta=\pi/2", space=:relative, align=(:right, :top), fontsize=FS_SMALL)
    axislegend(ax; position=:rc, nbanks=4)
    panel_label!(fig, (3, 1), "(d)")
    rowsize!(fig.layout, 1, Relative(0.43)); rowsize!(fig.layout, 2, Relative(0.30)); rowsize!(fig.layout, 3, Relative(0.27))
    rowgap!(fig.layout, 4)
    savefig(fig, "fig_bulk"; dir=FIGDIR)
end
const FIGURES = Dict("F1" => figure_hero, "F4" => figure_transfer,
    "F6" => figure_corrector, "F7" => figure_bulk)

function main()
    todo = isempty(OPTIONS.names) ? ["F1", "F4", "F6", "F7"] : OPTIONS.names
    for name in todo
        haskey(FIGURES, name) || error("Unknown figure $name; choose F1, F4, F6 or F7")
    end
    for name in todo
        println("Figure $name")
        FIGURES[name]()
    end
end
main()
