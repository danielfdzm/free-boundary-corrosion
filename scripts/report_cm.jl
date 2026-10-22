#!/usr/bin/env julia
# Reproducible C^m postprocessing of immutable numerical records.
using FreeBoundaryNumerics, Printf, SHA
include(joinpath(@__DIR__, "cm_diagnostics.jl"))

function main()
    root = normpath(joinpath(@__DIR__, ".."))
    datadir, outdir = joinpath(root, "data"), joinpath(root, "outputs", "cm_diagnostics")
    i = 1
    while i <= length(ARGS)
        if ARGS[i] in ("--help", "-h")
            println("julia --project=. scripts/report_cm.jl [--data DIR] [--out DIR]")
            return
        end
        i < length(ARGS) || error("Missing directory after $(ARGS[i])")
        if ARGS[i] == "--data"
            datadir = abspath(ARGS[i+1])
        elseif ARGS[i] == "--out"
            outdir = abspath(ARGS[i+1])
        else
            error("Unknown argument: $(ARGS[i])")
        end
        i += 2
    end
    mkpath(outdir)
    digest(path) = open(io -> bytes2hex(sha256(io)), path)
    inputs = ["E5.jld2", "E5_coupled_refinement.jld2", "heterogeneous/heterogeneous.jld2"]
    isfile(joinpath(datadir, "E5_bulk_fine.jld2")) && push!(inputs, "E5_bulk_fine.jld2")
    initial_digests = Dict(name => digest(joinpath(datadir, name)) for name in inputs)
    E = load_results(joinpath(datadir, "E5.jld2"))
    C = load_results(joinpath(datadir, "E5_coupled_refinement.jld2"))
    H = load_results(joinpath(datadir, "heterogeneous", "heterogeneous.jld2"))
    disk, coupled, heterogeneous = disk_cm(E), coupled_cm(C), heterogeneous_cm(H)
    disk32, coupled32, heterogeneous32 = disk_cm(E; oversample=32), coupled_cm(C; oversample=32), heterogeneous_cm(H; oversample=32)
    generated = String[]
    function csv(name, header, rows)
        open(joinpath(outdir, name), "w") do io
            println(io, join(header, ','))
            for row in rows
                println(io, join(row, ','))
            end
        end
        push!(generated, name)
    end
    orders(k, matrix) = permutedims(hcat([observed_orders(k, matrix[m, :]) for m in axes(matrix, 1)]...))
    diskorders = map(v -> orders(E["disk/kappas"], v), disk)
    csv("disk_convergence.csv", ["kappa", "m", "deviation_exact", "order_deviation_exact", "remainder_exact", "order_remainder_exact", "deviation_matched", "order_deviation_matched", "remainder_matched", "order_remainder_matched"],
        ([k, m-1, disk.dev_exact[m,j], diskorders.dev_exact[m,j], disk.rem_exact[m,j], diskorders.rem_exact[m,j], disk.dev[m,j], diskorders.dev[m,j], disk.rem[m,j], diskorders.rem[m,j]] for (j,k) in enumerate(E["disk/kappas"]) for m in 1:4))
    Rdelta = cm_norms(E["disk/limit/R"] - E["disk/limit/R_exact"])
    Wdelta = cm_norms(E["disk/limit/W"] - E["disk/limit/W_exact"])
    csv("disk_reference.csv", ["m", "limit_difference", "corrector_difference", "kappa_min_weighted_corrector_difference", "fraction_of_exact_remainder"],
        ([m-1, Rdelta[m], Wdelta[m], E["disk/kappas"][end]*Wdelta[m], (Rdelta[m]+E["disk/kappas"][end]*Wdelta[m])/disk.rem_exact[m,end]] for m in 1:4))
    coupledkeys = ["production/corrected_remainder", "space_2/production_difference", "time_2/production_difference", "time_4/production_difference",
        "time_2/previous_step_difference", "time_4/previous_step_difference", "space_2/corrected_remainder", "time_2/corrected_remainder", "time_4/corrected_remainder",
        "space_2/discarded_modes", "time_2/discarded_modes", "time_4/discarded_modes"]
    csv("coupled_refinement.csv", ["quantity", "m", "value", "fraction_of_production_remainder"],
        ([key, m-1, coupled[key][m], coupled[key][m]/coupled["production/corrected_remainder"][m]] for key in coupledkeys for m in 1:4))
    csv("coupled_time_orders.csv", ["m", "order", "production_to_half_step", "half_to_quarter_step"],
        ([m-1, coupled["time_orders"][m], coupled["time_2/previous_step_difference"][m], coupled["time_4/previous_step_difference"][m]] for m in 1:4))
    csv("heterogeneous_convergence.csv", ["kappa", "m", "deviation", "order_deviation", "remainder", "order_remainder", "matched_deviation", "matched_remainder"],
        ([k, m-1, heterogeneous["errors/deviation"][m,j], heterogeneous["orders/deviation"][m,j], heterogeneous["errors/remainder"][m,j], heterogeneous["orders/remainder"][m,j], heterogeneous["errors/matched_deviation"][m,j], heterogeneous["errors/matched_remainder"][m,j]] for (j,k) in enumerate(H["kappas"]) for m in 1:2))
    heterokeys = sort([key for key in keys(heterogeneous) if startswith(key, "sensitivity/") && !occursin("time_order", key)])
    csv("heterogeneous_refinement.csv", ["quantity", "m", "value", "fraction_of_production_remainder"],
        ([key, m-1, heterogeneous[key][m], heterogeneous[key][m]/heterogeneous["errors/remainder"][m,end]] for key in heterokeys for m in 1:2))
    csv("heterogeneous_time_orders.csv", ["quantity", "C0", "C1"],
        ([key, heterogeneous[key]...] for key in sort(collect(keys(heterogeneous))) if startswith(key, "sensitivity/time_order/")))

    # Separate plotted/convergence quantities from tiny refinement discrepancies.
    samplingrows = Vector{Any}[]
    relative(a, b) = abs(a-b) / max(abs(b), floatmin(Float64))
    for name in keys(disk), m in 1:4
        values, finer = getproperty(disk, name)[m,:], getproperty(disk32, name)[m,:]
        push!(samplingrows, ["disk", string(name), m-1, maximum(relative.(values, finer)), maximum(abs.(values-finer))])
    end
    for key in coupledkeys, m in 1:4
        push!(samplingrows, ["coupled", key, m-1, relative(coupled[key][m], coupled32[key][m]), abs(coupled[key][m]-coupled32[key][m])])
    end
    for key in vcat(["errors/deviation", "errors/remainder", "errors/matched_deviation", "errors/matched_remainder"], heterokeys), m in 1:2
        values, finer = heterogeneous[key], heterogeneous32[key]
        if values isa AbstractMatrix
            push!(samplingrows, ["heterogeneous", key, m-1, maximum(relative.(values[m,:], finer[m,:])), maximum(abs.(values[m,:]-finer[m,:]))])
        else
            push!(samplingrows, ["heterogeneous", key, m-1, relative(values[m], finer[m]), abs(values[m]-finer[m])])
        end
    end
    csv("oversampling.csv", ["dataset", "quantity", "m", "max_relative_change_16_to_32", "max_absolute_change_16_to_32"], samplingrows)

    # Removing |k|>16 here is only a roundoff diagnostic, never the reported norm.
    lowmodes(v) = begin
        coefficients = rfft(v)
        coefficients[18:end] .= 0
        irfft(coefficients, length(v))
    end
    highmodevectors = [
        "production_exact_remainder" => C["production/R"]-C["reference_R"]-C["kappa"]*C["reference_W"],
        "coupled_spatial_difference" => C["space_2/R_common"]-C["production/R"],
        "time_production_to_half" => C["time_2/R_common"]-C["production/R"],
        "time_half_to_quarter" => C["time_4/R_common"]-C["time_2/R_common"],
    ]
    csv("high_mode_diagnostic.csv", ["quantity", "m", "full_norm", "low16_norm", "high16_norm", "high16_fraction_of_full"],
        ([name, m-1, cm_norms(v)[m], cm_norms(lowmodes(v))[m], cm_norms(v-lowmodes(v))[m], cm_norms(v-lowmodes(v))[m]/cm_norms(v)[m]] for (name,v) in highmodevectors for m in 1:4))

    bulkorders = nothing
    if "E5_bulk_fine.jld2" in inputs
        B = load_results(joinpath(datadir, "E5_bulk_fine.jld2"))
        bulk = [maximum(abs, B["$k/phi"]-B["limit/u0"]) for k in B["kappas"]]
        bulkorders = observed_orders(B["kappas"], bulk)
        csv("bulk_c0.csv", ["kappa", "pulled_back_P1_C0", "order_to_next_kappa", "C0_divided_by_kappa"],
            ([k, bulk[j], bulkorders[j], bulk[j]/k] for (j,k) in enumerate(B["kappas"])))
    end

    report = IOBuffer()
    println(report, "C^m diagnostics from stored radii; Julia ", VERSION)
    println(report, "Norm: max over j=0:m of |d_theta^j I_N v| on 16N nodes; real cosine Nyquist convention.")
    println(report, "All convergence errors are at the final time T=0.5; these are numerical checks, not proofs of uniform-in-time rates.")
    println(report, "Disk final-pair conductivity orders (kappa=2^-7 to 2^-8):")
    for m in 1:4
        @printf(report, "  C%d: exact dev %.8f; exact remainder %.8f; matched dev %.8f; matched remainder %.8f; exact remainder %.10e\n",
            m-1, diskorders.dev_exact[m,end-1], diskorders.rem_exact[m,end-1], diskorders.dev[m,end-1], diskorders.rem[m,end-1], disk.rem_exact[m,end])
        @printf(report, "       spatial difference/remainder %.6g; time differences %.6e, %.6e; time order %.6f\n",
            coupled["space_2/production_difference"][m]/disk.rem_exact[m,end], coupled["time_2/previous_step_difference"][m], coupled["time_4/previous_step_difference"][m], coupled["time_orders"][m])
    end
    println(report, "Heterogeneous final-pair conductivity orders (kappa=2^-5 to 2^-6):")
    for m in 1:2
        sensitivities = [heterogeneous[key][m]/heterogeneous["errors/remainder"][m,end] for key in heterokeys if !endswith(key, "fixed_reference_remainder")]
        @printf(report, "  C%d: deviation %.8f; remainder %.8f; remainder %.10e; largest individual refinement/reference sensitivity fraction %.6g; forcing-only relative difference %.6g\n",
            m-1, heterogeneous["orders/deviation"][m,end-1], heterogeneous["orders/remainder"][m,end-1], heterogeneous["errors/remainder"][m,end], maximum(sensitivities), heterogeneous["forcing_only/relative_difference"][m])
    end
    disk_sampling = maximum(row[4] for row in samplingrows if row[1] == "disk" && row[2] in ("dev_exact", "rem_exact"))
    hetero_sampling = maximum(row[4] for row in samplingrows if row[1] == "heterogeneous" && row[2] in ("errors/deviation", "errors/remainder"))
    refinement_sampling = maximum(row[4] for row in samplingrows if row[1] == "coupled" || startswith(row[2], "sensitivity/"))
    @printf(report, "Maximum relative 16N-to-32N change: disk plotted errors %.8e; heterogeneous plotted errors %.8e; refinement quantities %.8e.\n", disk_sampling, hetero_sampling, refinement_sampling)
    println(report, "C3 time-step discrepancies are consistent with amplified high-frequency floating-point noise; they do not establish a temporal order. The full unfiltered norms remain reported.")
    for (name,v) in highmodevectors
        @printf(report, "  %s: full C3 %.8e; modes |k|<=16 C3 %.8e; high modes C3 %.8e\n", name, cm_norms(v)[4], cm_norms(lowmodes(v))[4], cm_norms(v-lowmodes(v))[4])
    end
    !isnothing(bulkorders) && @printf(report, "Bulk pulled-back P1 C0: smallest-conductivity pair order %.8f. No higher spatial derivative rate is inferred from P1 fields.\n", bulkorders[end-1])
    contents = String(take!(report))
    write(joinpath(outdir, "report.txt"), contents)
    push!(generated, "report.txt")
    csv("data_sha256.csv", ["path_relative_to_data", "sha256"], ([name, initial_digests[name]] for name in inputs))
    sources = ["scripts/cm_diagnostics.jl", "scripts/report_cm.jl", "scripts/check_cm.jl", "src/spectral.jl", "Manifest.toml"]
    csv("postprocessing_sha256.csv", ["path_relative_to_project", "sha256"], ([name, digest(joinpath(root, name))] for name in sources))
    artifacts = sort(copy(generated))
    csv("artifacts_sha256.csv", ["path_relative_to_report", "sha256"], ([name, digest(joinpath(outdir, name))] for name in artifacts))
    for name in inputs
        digest(joinpath(datadir, name)) == initial_digests[name] || error("Raw input changed during postprocessing: $name")
    end
    print(contents)
    println("Wrote reproducible reports to ", outdir, "; raw records unchanged.")
end

main()
