#!/usr/bin/env julia
# Space-time visualization of an actual stored planar-interface trajectory.
# No simulation is performed and the stored numerical record is never changed.
using CairoMakie, FreeBoundaryNumerics, LinearAlgebra
CairoMakie.activate!(type="pdf")

const EVOLUTION_ROOT = normpath(joinpath(@__DIR__, ".."))
const EVOLUTION_CYAN = colorant"#00B8F0"
const EVOLUTION_BLUE = colorant"#246BFD"
const EVOLUTION_VIOLET = colorant"#7A3FFC"
const EVOLUTION_MAGENTA = colorant"#E83EAC"
const EVOLUTION_ORANGE = colorant"#FF8A1F"
const EVOLUTION_YELLOW = colorant"#FFD43B"
const EVOLUTION_INK = colorant"#17324D"
const EVOLUTION_MAP = cgrad([EVOLUTION_CYAN, EVOLUTION_BLUE, EVOLUTION_VIOLET,
    EVOLUTION_MAGENTA, EVOLUTION_ORANGE, EVOLUTION_YELLOW])

function main()
    length(ARGS) <= 2 || error("Usage: plot_evolution3d.jl [DATA_DIR] [OUTPUT_DIR]")
    datadir = isempty(ARGS) ? joinpath(EVOLUTION_ROOT, "data", "heterogeneous") : abspath(ARGS[1])
    outdir = length(ARGS) < 2 ? joinpath(EVOLUTION_ROOT, "outputs", "figures") : abspath(ARGS[2])
    record_path = joinpath(datadir, "heterogeneous.jld2")
    d = load_results(record_path)
    idx = findfirst(==(0.25), d["kappas"])
    isnothing(idx) && error("The supplied record lacks the kappa=1/4 trajectory")
    run = d["coupled"][idx]
    theta, times = d["theta"], run["times"]
    R = run["Rhist"]
    @assert size(R) == (length(theta), length(times))
    @assert times[1] == 0 && times[end] == 0.5
    @assert d["i0_formula"] == "1+0.25*x+0.10*y^2"
    # Close the angular seam. Every time slice is an actual saved Heun state.
    th = vcat(theta, 2π)
    radius = vcat(R, R[1:1, :])
    X = radius .* cos.(th)
    Y = radius .* sin.(th)
    Z = ones(length(th)) * times'
    current = @. 1 + 0.25X + 0.10Y^2
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

    set_theme!(Theme(fontsize=15,
        Axis3=(xticklabelsize=13, yticklabelsize=13, zticklabelsize=13,
            xlabelsize=18, ylabelsize=18, zlabelsize=16),
        Colorbar=(ticklabelsize=13, labelsize=15),
        Legend=(framevisible=false, labelsize=13, patchsize=(23,10), colgap=20)))
    fig = Figure(size=(680, 495), figure_padding=(12, 22, 8, 13), backgroundcolor=:white)
    Label(fig[0, 1], "Evolution of a heterogeneous interface";
        fontsize=19, font=:bold, color=EVOLUTION_INK, tellwidth=false)
    Label(fig[1, 1], L"(x_1,x_2,t),\quad \kappa=1/4";
        fontsize=16, color=EVOLUTION_INK, tellwidth=false)

    ax = Axis3(fig[2, 1]; azimuth=view_azimuth, elevation=0.38,
        aspect=(1, 1, 1.08), perspectiveness=0.20, viewmode=:fit,
        xlabel=L"x_1", ylabel=L"x_2", zlabel="Time t",
        xticks=[-1,0,1], yticks=[-1,0,1], zticks=([0,.1,.2,.3,.4,.5], ["0","0.1","0.2","0.3","0.4","0.5"]),
        xlabeloffset=19, ylabeloffset=19, zlabeloffset=37,
        xgridcolor=(:black, .065), ygridcolor=(:black, .065), zgridcolor=(:black, .055),
        xypanelcolor=(:white, 0), xzpanelcolor=(:white, 0), yzpanelcolor=(:white, 0),
        xspinecolor_1=(:black,.34), yspinecolor_1=(:black,.34), zspinecolor_1=(:black,.34),
        xspinecolor_2=(:black,.14), yspinecolor_2=(:black,.14), zspinecolor_2=(:black,.14),
        xspinecolor_3=(:black,.14), yspinecolor_3=(:black,.14), zspinecolor_3=(:black,.14),
        protrusions=(36, 55, 15, 10))
    limits!(ax, -1.28, 1.28, -1.28, 1.28, 0, .51)
    shell = surface!(ax, X, Y, Z; color=current, colormap=EVOLUTION_MAP,
        colorrange=(.70,1.35), shading=NoShading, rasterize=4)
    # The white rings are measured interface states, not a decorative grid.
    for t in .05:.05:.45
        j = argmin(abs.(times .- t))
        lines!(ax, mask_line(X[:,j], front[:,j]), Y[:,j], Z[:,j];
            color=(:white,.76), linewidth=1.1)
    end
    # Selected fixed angular coordinates reveal the taper of the radial graph.
    for j in 1:64:length(theta)
        lines!(ax, mask_line(X[j,:], front[j,:]), Y[j,:], Z[j,:];
            color=(EVOLUTION_INK,.15), linewidth=.65)
    end
    initial = lines!(ax, mask_line(X[:,1], front[:,1]), Y[:,1], Z[:,1];
        color=EVOLUTION_BLUE, linewidth=2.3)
    final = lines!(ax, X[:,end], Y[:,end], Z[:,end]; color=EVOLUTION_MAGENTA, linewidth=2.8)
    # A fine dark rim preserves the endpoint geometry over the bright surface.
    lines!(ax, X[:,end], Y[:,end], Z[:,end]; color=(EVOLUTION_INK,.55), linewidth=.7)
    Colorbar(fig[2, 2], shell; label=L"i_0(x_1,x_2)", width=13,
        ticks=[.75,1.0,1.25], height=Relative(.66), tellheight=false)
    Legend(fig[3, 1], [initial, final], [L"\Gamma_{\kappa}(0)", L"\Gamma_{\kappa}(0.5)"];
        orientation=:horizontal, tellwidth=false)
    Label(fig[4, 1], "Height represents time; every horizontal slice is a planar interface.";
        fontsize=12.5, color=EVOLUTION_INK, tellwidth=false)
    rowgap!(fig.layout, 4)
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
