# House style shared by all figures: palette locked to the paper's TikZ colours,
# Computer Modern math via LaTeXStrings, 10pt text, no top/right spines.
using CairoMakie, GeometryBasics, Contour, JLD2, Printf, Statistics, LinearAlgebra
using FreeBoundaryNumerics
CairoMakie.activate!(type="pdf")

const WATER = colorant"#125EB9"
const DEEP = colorant"#0B3C7A"
const SOLIDDARK = colorant"#703416"
const SOLIDBASE = colorant"#CD7534"
const NEUTRAL = colorant"#F5F5F2"
const WATERFILL = colorant"#78B9F0"
const ELECTRO = cgrad([DEEP, WATER, NEUTRAL, SOLIDBASE, SOLIDDARK])
const GREY = colorant"#6E6E6E"
const WIDTH = 566.0          # px units; with pt_per_unit = 0.75 this is 5.9 in = \textwidth
const FS = 13.3              # 10 pt
const FS_SMALL = 12.0        # 9 pt
const LINESTYLES = [:solid, :dash, :dot, :dashdot, (:dash, :dense), (:dot, :dense), (:dashdot, :dense), :solid, :dash]

set_theme!(Theme(
    fontsize=FS,
    Axis=(xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false,
        xticklabelsize=FS_SMALL, yticklabelsize=FS_SMALL, titlesize=FS, titlefont=:regular),
    Axis3=(xticklabelsize=FS_SMALL, yticklabelsize=FS_SMALL, zticklabelsize=FS_SMALL,
        xspinesvisible=true, yspinesvisible=true, zspinesvisible=true),
    Legend=(framevisible=false, labelsize=FS_SMALL, rowgap=0, patchsize=(18, 8)),
    Colorbar=(ticklabelsize=FS_SMALL, labelsize=FS, size=9),
    Lines=(linewidth=1.3,),
    Scatter=(markersize=6,),
))

"Colour of the family member kappa = 2^-k on the viridis scale (k in [0, kmax])."
kappa_color(k; kmax=8) = get(cgrad(:viridis), 0.92 * k / kmax)
kappa_style(k) = LINESTYLES[mod1(round(Int, k) + 1, length(LINESTYLES))]

"Colorbar for a kappa family, ticked 1, 1/4, 1/16, 1/64, 1/256."
function kappa_colorbar!(pos; kmax=8, label=L"\kappa", vertical=true)
    ticks = 0:2:kmax
    labels = [k == 0 ? "1" : "1/$(2^k)" for k in ticks]
    Colorbar(pos; colormap=cgrad(:viridis)[range(0, 0.92, length=64)], limits=(0, kmax), ticks=(collect(ticks), labels),
        label=label, vertical=vertical, flipaxis=vertical)
end

panel_label!(fig, pos, txt) = Label(fig[pos..., TopLeft()], txt; font=:bold, fontsize=FS, padding=(0, 6, 4, 0), halign=:right)

const PI_TICKS = ([0, π / 2, π, 3π / 2, 2π], ["0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"])

savefig(fig, name; dir) = (mkpath(dir); save(joinpath(dir, name * ".pdf"), fig; pt_per_unit=0.75); println("  wrote ", joinpath(dir, name * ".pdf")))

"Slope guide through (x0, y0) with the given slope on log-log axes."
function slope_guide!(ax, x, y0, x0, slope; color=GREY, label=nothing, linestyle=:dash)
    xs = [minimum(x), maximum(x)]
    ys = y0 .* (xs ./ x0) .^ slope
    lines!(ax, xs, ys; color=color, linestyle=linestyle, linewidth=0.9, label=label)
end

"Contour lines of a P1 field on the fitted annular mesh, computed on the logical (theta, eta) grid and mapped to the plane."
function logical_contours(R::Vector{Float64}, phi::Vector{Float64}, nt::Int, nr::Int, levels; B=2.0)
    Z = zeros(nt + 1, nr + 1)
    for k in 0:nr, j in 1:nt
        Z[j, k+1] = phi[k*nt+j]
    end
    Z[nt+1, :] .= Z[1, :]
    th = collect(0:nt) .* (2π / nt)
    eta = collect(0:nr) ./ nr
    Rext = vcat(R, R[1])
    segs = Vector{Vector{Point2f}}()
    for cl in Contour.levels(Contour.contours(th, eta, Z, levels))
        for line in Contour.lines(cl)
            xs, ys = Contour.coordinates(line)
            pts = Point2f[]
            for (t, e) in zip(xs, ys)
                j = clamp(floor(Int, t / (2π / nt)) + 1, 1, nt)
                frac = t / (2π / nt) - (j - 1)
                Rt = (1 - frac) * Rext[j] + frac * Rext[j+1]
                r = Rt + e * (B - Rt)
                push!(pts, Point2f(r * cos(t), r * sin(t)))
            end
            push!(segs, pts)
        end
    end
    return segs
end

"Piecewise-constant gradient of the P1 field as a function of the position (zero outside the electrolyte)."
function p1_gradient_field(R::Vector{Float64}, phi::Vector{Float64}, x::Vector{Float64}, y::Vector{Float64}, mesh)
    nt, nr, B = mesh.nt, mesh.nr, mesh.B
    dth = 2π / nt
    return function (p)
        px, py = p[1], p[2]
        th = mod(atan(py, px), 2π)
        j = clamp(floor(Int, th / dth) + 1, 1, nt)
        jp = j == nt ? 1 : j + 1
        frac = th / dth - (j - 1)
        Rt = (1 - frac) * R[j] + frac * R[jp]
        r = hypot(px, py)
        eta = (r - Rt) / (B - Rt)
        (eta <= 0.002 || eta >= 0.995) && return Point2f(0, 0)
        k = clamp(floor(Int, eta * nr), 0, nr - 1)
        u = eta * nr - k
        a = k * nt + j; b = k * nt + jp; c = (k + 1) * nt + jp; d = (k + 1) * nt + j
        tri = u < frac ? (a, c, b) : (a, d, c)
        x1, y1, x2, y2, x3, y3 = x[tri[1]], y[tri[1]], x[tri[2]], y[tri[2]], x[tri[3]], y[tri[3]]
        det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        b1, b2, b3 = y2 - y3, y3 - y1, y1 - y2
        c1, c2, c3 = x3 - x2, x1 - x3, x2 - x1
        gx = (phi[tri[1]] * b1 + phi[tri[2]] * b2 + phi[tri[3]] * b3) / det
        gy = (phi[tri[1]] * c1 + phi[tri[2]] * c2 + phi[tri[3]] * c3) / det
        return Point2f(gx, gy)
    end
end

"Mesh coordinates and triangulation for the field plots."
function field_mesh(R::Vector{Float64}, nr::Int; B=2.0)
    nt = length(R)
    mesh = AnnulusMesh(nt, nr; B=B)
    x, y = node_coordinates(mesh, R)
    return mesh, x, y
end

"Filled field on the fitted mesh with the solid interior flat, the interface and the wall."
function field_panel!(ax, R::Vector{Float64}, phi::Vector{Float64}, nr::Int; colorrange, colormap=ELECTRO, levels=24,
        contours=true, ncontours=12, interface=true, wall=true, streamlines=false, B=2.0)
    mesh, x, y = field_mesh(R, nr; B=B)
    nt = mesh.nt
    tc = tricontourf!(ax, x, y, phi; triangulation=mesh.tri, levels=range(colorrange[1], colorrange[2], length=levels + 1),
        colormap=colormap, extendlow=:auto, extendhigh=:auto)
    if contours
        lv = range(colorrange[1], colorrange[2], length=ncontours + 2)[2:end-1]
        for seg in logical_contours(R, phi, nt, nr, collect(lv); B=B)
            lines!(ax, seg; color=(:white, 0.75), linewidth=0.5)
        end
    end
    if streamlines
        f = p1_gradient_field(R, phi, x, y, mesh)
        streamplot!(ax, f, -B .. B, -B .. B; color=p -> RGBAf(1, 1, 1, 0.45), linewidth=0.55, arrow_size=0,
            gridsize=(28, 28), density=0.9, stepsize=0.01, maxsteps=400)
    end
    th = mesh.theta
    poly!(ax, Point2f.(R .* cos.(th), R .* sin.(th)); color=colorant"#EFECE6", strokecolor=SOLIDDARK, strokewidth=interface ? 1.6 : 0)
    if wall
        tw = range(0, 2π, length=361)
        lines!(ax, B .* cos.(tw), B .* sin.(tw); color=WATER, linestyle=:dash, linewidth=0.9)
    end
    hidedecorations!(ax); hidespines!(ax)
    return tc
end

fmt_sci(x) = x == 0 ? "0" : (@sprintf("%.2e", x) |> s -> begin
    m, e = split(s, "e"); ee = parse(Int, e)
    "\$$(m)\\times10^{$(ee)}\$"
end)
fmt_ord(x) = isnan(x) ? "---" : @sprintf("%.2f", x)
