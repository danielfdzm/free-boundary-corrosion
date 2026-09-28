#!/usr/bin/env julia
# Frozen expansion, conductivity sweep, and refinement-budget tables from stored data.
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

function table_frozen()
    E3 = load_results(joinpath(RESDIR, "E3.jld2"))
    kappas = E3["kappas"]; ns = E3["ns"]; deltas = E3["deltas"]
    jn = findfirst(==(2), ns); id = findfirst(==(0.2), deltas)
    e = E3["e_inf"][:, jn, id]; r = E3["remainder_inf"][:, jn, id]
    dJ = abs.(E3["J_kappa"][:, jn, id] .- E3["J0_h"][:, jn, id])
    dJ2 = abs.(E3["J_kappa"][:, jn, id] .- E3["J0_h"][:, jn, id] .- kappas .* E3["J1_h"][:, jn, id])
    oe, orr, oJ, oJ2 = (observed_orders(kappas, v) for v in (e, r, dJ, dJ2))
    io = IOBuffer()
    println(io, "\\begin{tabular}{lrrrrrrrr}")
    println(io, "\\toprule")
    println(io, "\$\\kappa\$ & \$\\|e_\\kappa\\|_\\infty\$ & ord. & \$\\|e_\\kappa-\\kappa e_1^h\\|_\\infty\$ & ord. & \$|J_{\\kappa,h}-J_{0,h}|\$ & ord. & \$|J_{\\kappa,h}-J_{0,h}-\\kappa J_{1,h}|\$ & ord.\\\\")
    println(io, "\\midrule")
    for i in eachindex(kappas)
        println(io, kap(kappas[i]), " & ", sci(e[i]), " & ", ord(oe[i]), " & ", sci(r[i]), " & ", ord(orr[i]), " & ", sci(dJ[i]), " & ", ord(oJ[i]), " & ", sci(dJ2[i]), " & ", ord(oJ2[i]), "\\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    write(joinpath(TABDIR, "frozen_expansion.tex"), String(take!(io)))
    println("  wrote frozen_expansion.tex  (J0_h = $(E3["J0_h"][1, jn, id]), J1_h = $(E3["J1_h"][1, jn, id]))")
end

function table_sweep()
    E5 = load_results(joinpath(RESDIR, "E5.jld2"))
    kappas = E5["disk/kappas"]
    dev = E5["disk/dev_norms"]; rem = E5["disk/rem_norms"]; devf = E5["disk/dev_norms_fine"]; remf = E5["disk/rem_norms_fine"]
    exc = E5["disk/excess"]
    o = (v) -> observed_orders(kappas, v)
    io = IOBuffer()
    println(io, "\\begin{tabular}{@{}l@{\\;\\;}rrrrrr@{\\quad}rr@{\\quad}rr@{}}")
    println(io, "\\toprule")
    println(io, " & \\multicolumn{6}{c}{matched discrete reference \$(R_0^h,h_1^h)\$} & \\multicolumn{2}{c}{refined reference} & \\multicolumn{2}{c}{excess}\\\\")
    println(io, "\\cmidrule(lr){2-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}")
    println(io, "\$\\kappa\$ & \$\\|R_\\kappa-R_0^h\\|_{H^0}\$ & ord. & \$\\|R_\\kappa-R_0^h\\|_{H^2}\$ & ord. & \$\\|\\rho_\\kappa^h\\|_{H^0}\$ & ord. & \$\\|R_\\kappa-R_0\\|_{H^0}\$ & \$\\|\\rho_\\kappa\\|_{H^0}\$ & \$\\int d-\\int i_0\$ & ord.\\\\")
    println(io, "\\midrule")
    od0, od2, or0, oex = o(dev[1, :]), o(dev[3, :]), o(rem[1, :]), o(exc)
    for i in eachindex(kappas)
        println(io, kap(kappas[i]), " & ", sci(dev[1, i]), " & ", ord(od0[i]), " & ", sci(dev[3, i]), " & ", ord(od2[i]), " & ", sci(rem[1, i]), " & ", ord(or0[i]),
            " & ", sci(devf[1, i]), " & ", sci(remf[1, i]), " & ", sci(exc[i]), " & ", ord(oex[i]), "\\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    write(joinpath(TABDIR, "conductivity_sweep.tex"), String(take!(io)))
    println("  wrote conductivity_sweep.tex")
end

table_frozen()
table_sweep()
refinement_path = joinpath(RESDIR, "E5_refinement.jld2")
if isfile(refinement_path)
    include(joinpath(@__DIR__, "refinement_table.jl"))
    write_refinement_table(load_results(refinement_path), TABDIR)
    println("  wrote refinement_budget.tex and refinement_budget.csv")
end
coupled_refinement_path = joinpath(RESDIR, "E5_coupled_refinement.jld2")
if isfile(coupled_refinement_path)
    include(joinpath(@__DIR__, "coupled_refinement_table.jl"))
    write_coupled_refinement_table(load_results(coupled_refinement_path), TABDIR)
    println("  wrote coupled_refinement.tex and coupled_refinement.csv")
end
