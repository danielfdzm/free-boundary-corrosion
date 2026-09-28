# Table generation shared by make_tables.jl and run_refinement.jl.
# Entries are measured reference differences, not rigorous error bounds.
function write_refinement_table(record, outdir)
    κ = record["kappa"]
    rows = [
        (raw"Limiting-flow refinement", record["limit_reference_difference"]),
        (raw"$\kappa$ times corrector spatial refinement", κ .* record["corrector_spatial_difference"]),
        (raw"$\kappa$ times corrector time refinement", κ .* record["factor_2/previous_step_difference"]),
        (raw"Matched corrected remainder $\rho_\kappa^h$", record["matched_remainder"]),
        (raw"Refined corrected remainder $\rho_\kappa$", record["refined_remainder"]),
    ]
    scientific(x) = begin
        mantissa, exponent = split(@sprintf("%.2e", x), "e")
        "\$$(mantissa)\\times10^{$(parse(Int, exponent))}\$"
    end
    mkpath(outdir)
    open(joinpath(outdir, "refinement_budget.tex"), "w") do io
        println(io, raw"\begin{tabular}{lrrr}")
        println(io, raw"\toprule")
        println(io, raw"Quantity & $H^0$ & $H^2$ & $H^3$", "\\\\")
        println(io, raw"\midrule")
        for (label, values) in rows
            println(io, label, " & ", join(scientific.(values[[1, 3, 4]]), " & "), "\\\\")
        end
        println(io, raw"\bottomrule")
        println(io, raw"\end{tabular}")
    end
    open(joinpath(outdir, "refinement_budget.csv"), "w") do io
        println(io, "quantity,H0,H1,H2,H3")
        for (label, values) in rows
            println(io, label, ",", join(values, ","))
        end
        println(io, "corrector temporal order,", join(record["time_orders"], ","))
    end
end
