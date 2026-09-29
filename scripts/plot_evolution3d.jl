#!/usr/bin/env julia
# Space-time visualization of an actual stored planar-interface trajectory.
# No simulation is performed and the stored numerical record is never changed.
include(joinpath(@__DIR__, "figure_style.jl"))

const EVOLUTION_ROOT = normpath(joinpath(@__DIR__, ".."))
# Drawn at the paper's scale for \includegraphics[width=0.7\textwidth], so that
# fonts and the colorbar match the other figures.
const EVOLUTION_SIZE = (round(Int, 0.7WIDTH), 345)

function main()
    length(ARGS) <= 2 || error("Usage: plot_evolution3d.jl [DATA_DIR] [OUTPUT_DIR]")
    datadir = isempty(ARGS) ? joinpath(EVOLUTION_ROOT, "data", "heterogeneous") : abspath(ARGS[1])
    outdir = length(ARGS) < 2 ? joinpath(EVOLUTION_ROOT, "outputs", "figures") : abspath(ARGS[2])
    record_path = joinpath(datadir, "long_time.jld2")
    d = load_results(record_path)
    run = d["coupled"]
    run["kappa"] == 0.25 || error("The supplied record lacks the kappa=1/4 trajectory")
    theta, times, T = d["theta"], run["times"], run["T"]
    R = run["Rhist"]
    @assert size(R) == (length(theta), length(times))
    @assert times[1] == 0 && times[end] == T == 2.5
    @assert d["i0_formula"] == "1+0.25*x+0.10*y^2"
    # Close the angular seam. Every time slice is an actual saved Heun state.
    th = vcat(theta, 2π)
    radius = vcat(R, R[1:1, :])
    X = radius .* cos.(th)
    Y = radius .* sin.(th)
    Z = ones(length(th)) * times'
    current = @. 1 + 0.25X + 0.10Y^2
    # For the display of i0 only, the surface is sampled on a four times finer
    # angular grid: each saved interface is refined by trigonometric interpolation,
    # which is its Fourier representation in the scheme. The lines use the saved nodes.
    nfine = 4length(theta)
    thfine = vcat(theta_grid(nfine), 2π)
    Rfine = reduce(hcat, [resample(R[:, j], nfine) for j in eachindex(times)])
    Rfine = vcat(Rfine, Rfine[1:1, :])
    Xfine = Rfine .* cos.(thfine)
    Yfine = Rfine .* sin.(thfine)
    Zfine = ones(length(thfine)) * times'
    # The colour scale of Figure 4(a): 24 bands from LaTeX's blue at the minimum
    # reached, through cyan, white at i0 = 1, yellow and orange, to LaTeX's red at the maximum.
    lo, hi = extrema(current)
    bands = banded_diverging(lo, hi; center=1.0, low=COOL_ELECTRO, high=WARM_ELECTRO)
    view_azimuth = -1.10
    # CairoMakie does not depth-test separate line plots against the surface.
    # Restrict the intermediate contour overlays to outward-facing portions;
    # the complete initial and final slices remain represented by the surface.
    front = falses(size(X))
    for j in eachindex(times)
        dR = vcat(fourier_derivative(R[:,j]), fourier_derivative(R[:,j])[1])
        nx = @. radius[:,j]*cos(th) + dR*sin(th)
        ny = @. radius[:,j]*sin(th) - dR*cos(th)
        front[:,j] .= nx .* cos(view_azimuth) .+ ny .* sin(view_azimuth) .> 0
    end
    mask_line(v, visible) = [visible[i] ? v[i] : NaN for i in eachindex(v)]

    fig = Figure(size=EVOLUTION_SIZE, figure_padding=(2, 4, 2, 2), backgroundcolor=:white)

    zt = 0:0.5:T
    ax = Axis3(fig[1, 1]; azimuth=view_azimuth, elevation=0.38,
        aspect=(1, 1, 1.08), perspectiveness=0.20, viewmode=:fitzoom,
        xlabel=L"x_1", ylabel=L"x_2", zlabel="Time t",
        xticks=[-1,0,1], yticks=[-1,0,1], zticks=(collect(zt), [@sprintf("%g", t) for t in zt]),
        xlabeloffset=22, ylabeloffset=24, zlabeloffset=31,
        xgridcolor=(:black, .065), ygridcolor=(:black, .065), zgridcolor=(:black, .055),
        xypanelcolor=(:white, 0), xzpanelcolor=(:white, 0), yzpanelcolor=(:white, 0),
        xspinecolor_1=(:black,.34), yspinecolor_1=(:black,.34), zspinecolor_1=(:black,.34),
        xspinecolor_2=(:black,.14), yspinecolor_2=(:black,.14), zspinecolor_2=(:black,.14),
        xspinecolor_3=(:black,.14), yspinecolor_3=(:black,.14), zspinecolor_3=(:black,.14),
        protrusions=(30, 46, 0, 0))
    limits!(ax, -1.28, 1.28, -1.28, 1.28, 0, 1.02T)
    shell = surface!(ax, Xfine, Yfine, Zfine; color=@.(1 + 0.25Xfine + 0.10Yfine^2),
        colormap=bands, colorrange=(lo, hi), shading=NoShading, rasterize=6)
    # The rings are measured interface states, not a decorative grid.
    for t in 0.25:0.25:T-0.25
        j = argmin(abs.(times .- t))
        lines!(ax, mask_line(X[:,j], front[:,j]), Y[:,j], Z[:,j];
            color=(:black,.55), linewidth=.75)
    end
    # Selected fixed angular coordinates reveal the taper of the radial graph.
    for j in 1:64:length(theta)
        lines!(ax, mask_line(X[j,:], front[j,:]), Y[j,:], Z[j,:];
            color=(:black,.15), linewidth=.55)
    end
    initial = lines!(ax, mask_line(X[:,1], front[:,1]), Y[:,1], Z[:,1];
        color=:black, linewidth=1.9, linestyle=:dash)
    final = lines!(ax, X[:,end], Y[:,end], Z[:,end]; color=:black, linewidth=2.3)
    Colorbar(fig[1, 2], shell; label=L"i_0(x_1,x_2)", ticks=0.8:0.1:1.2,
        height=Relative(.62), tellheight=false)
    tfinal = @sprintf("%g", T)
    Legend(fig[2, 1], [initial, final], [L"\Gamma_{\kappa}(0)", L"\Gamma_{\kappa}(%$tfinal)"];
        orientation=:horizontal, tellwidth=false, colgap=16)
    rowgap!(fig.layout, 9)
    colgap!(fig.layout, 4)

    # Publishing into figures/paper keeps its PNG in figures/previews.
    # Fresh/default reproduction outputs keep PDF and PNG together.
    pngdir = normpath(outdir) == joinpath(EVOLUTION_ROOT, "figures", "paper") ?
        joinpath(EVOLUTION_ROOT, "figures", "previews") : outdir
    mkpath(outdir); mkpath(pngdir)
    pdf = joinpath(outdir, "fig_evolution3d.pdf")
    png = joinpath(pngdir, "fig_evolution3d.png")
    save(pdf, fig; pt_per_unit=.75)
    save(png, fig; px_per_unit=3)
    println("Wrote ", pdf, " and ", png)
end

main()
