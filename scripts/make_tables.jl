#!/usr/bin/env julia
# Current conductivity-sweep and coupled-refinement tables from stored data.
#   julia --project=. scripts/make_tables.jl [--quick]
using FreeBoundaryNumerics, JLD2, Printf
include(joinpath(@__DIR__, "cli.jl"))
const OPTIONS = parse_cli(ARGS; output="tables", data=true,
    usage="julia --project=. scripts/make_tables.jl [--data DIR] [--out DIR] [--quick]")
isempty(OPTIONS.names) || error("make_tables.jl takes no positional arguments")
const RESDIR = OPTIONS.datadir
const TABDIR = OPTIONS.outdir
mkpath(TABDIR)

function sci(x)
    x == 0 && return "\$0\$"
    s = @sprintf("%.2e", x)
    m, e = split(s, "e")
    return "\$$(m)\\times10^{$(parse(Int, e))}\$"
end
ord(x) = isnan(x) ? "---" : @sprintf("\$%.2f\$", x)
kap(κ) = κ == 1 ? "\$1\$" : "\$2^{-$(round(Int, -log2(κ)))}\$"

function table_sweep()
    E5 = load_results(joinpath(RESDIR, "E5.jld2"))
    get(E5, "i_star", nothing) == 1.0 || error("Current tables require constant-current E5 data")
    kappas = E5["disk/kappas"]
    # Exact continuum references only: the matched-reference displacement agrees
    # to three digits and the matched norms are plotted in the corrector figure.
    devf = E5["disk/dev_norms_fine"]; remf = E5["disk/rem_norms_fine"]
    exc = E5["disk/excess"]
    o = (v) -> observed_orders(kappas, v)
    io = IOBuffer()
    println(io, "\\begin{tabular}{@{}lrrrrrr@{}}")
    println(io, "\\toprule")
    println(io, "\$\\kappa\$ & \$\\|R_\\kappa-R_0\\|_{H^0}\$ & ord. & \$\\|\\rho_\\kappa\\|_{H^0}\$ & ord. & \$\\int_{\\Gamma_\\kappa(T)} i_d-i_*|\\Gamma_\\kappa(T)|\$ & ord.\\\\")
    println(io, "\\midrule")
    od0, or0, oex = o(devf[1, :]), o(remf[1, :]), o(exc)
    for i in eachindex(kappas)
        println(io, kap(kappas[i]), " & ", sci(devf[1, i]), " & ", ord(od0[i]), " & ", sci(remf[1, i]), " & ", ord(or0[i]),
            " & ", sci(exc[i]), " & ", ord(oex[i]), "\\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    write(joinpath(TABDIR, "conductivity_sweep.tex"), String(take!(io)))
    println("  wrote conductivity_sweep.tex")
end

# Historical frozen-domain tables are retained only in archive/variable_i0.
table_sweep()
coupled_refinement_path = joinpath(RESDIR, "E5_coupled_refinement.jld2")
if isfile(coupled_refinement_path)
    include(joinpath(@__DIR__, "coupled_refinement_table.jl"))
    write_coupled_refinement_table(load_results(coupled_refinement_path), TABDIR)
    println("  wrote coupled_refinement.tex and coupled_refinement.csv")
end
