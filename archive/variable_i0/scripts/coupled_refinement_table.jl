# Independent coupled-solution discrepancies and fixed-reference remainders.
function write_coupled_refinement_table(record, outdir)
    rows = [
        (raw"Coupled spatial refinement", record["space_2/production_difference"]),
        (raw"Coupled time refinement ($\Delta t$ to $\Delta t/2$)", record["time_2/previous_step_difference"]),
        (raw"Coupled time refinement ($\Delta t/2$ to $\Delta t/4$)", record["time_4/previous_step_difference"]),
        (raw"Corrected remainder: production solve", record["production/corrected_remainder"]),
        (raw"Corrected remainder: spatially refined solve", record["space_2/corrected_remainder"]),
        (raw"Corrected remainder: $\Delta t/2$ solve", record["time_2/corrected_remainder"]),
        (raw"Corrected remainder: $\Delta t/4$ solve", record["time_4/corrected_remainder"]),
    ]
    scientific(x) = begin
        mantissa, exponent = split(@sprintf("%.2e", x), "e")
        "\$$(mantissa)\\times10^{$(parse(Int, exponent))}\$"
    end
    mkpath(outdir)
    open(joinpath(outdir, "coupled_refinement.tex"), "w") do io
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
    open(joinpath(outdir, "coupled_refinement.csv"), "w") do io
        println(io, "quantity,H0,H1,H2,H3")
        for (label, values) in rows
            println(io, label, ",", join(values, ","))
        end
        println(io, "coupled temporal order,", join(record["time_orders"], ","))
        println(io, "spatially refined discarded Fourier modes,", join(record["space_2/discarded_modes"], ","))
        for name in ("space_2", "time_2", "time_4")
            println(io, "$name difference / production corrected remainder,",
                join(record["$name/production_difference"] ./ record["production/corrected_remainder"], ","))
        end
    end
end
