#!/usr/bin/env julia
# Bulk fields of fig_bulk.pdf on the doubled E5 mesh; data/E5.jld2 is not read or changed.
#   julia --project=. scripts/run_bulk_fine.jl [--out DIR] [--only NAME] [--factor 2]
# NAME is limit, k0, k3, k6 or k8. Each run is cached under DIR/bulk_fine_runs, so
# separate --only processes may run in parallel before a final call assembles the record.
# With --factor 1 the runs repeat the corresponding E5 production runs.
using FreeBoundaryNumerics, LinearAlgebra, Printf, SHA, Dates
BLAS.set_num_threads(1)

const ROOT = normpath(joinpath(@__DIR__, ".."))
const KS = (0, 3, 6, 8)
digest(path) = open(io -> bytes2hex(sha256(io)), path)

function main()
    outdir = joinpath(ROOT, "outputs", "data")
    only = nothing
    factor = 2
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--out" && i < length(ARGS)
            outdir = abspath(ARGS[i + 1]); i += 2
        elseif ARGS[i] == "--only" && i < length(ARGS)
            only = ARGS[i + 1]; i += 2
        elseif ARGS[i] == "--factor" && i < length(ARGS)
            factor = parse(Int, ARGS[i + 1]); i += 2
        else
            error("Usage: run_bulk_fine.jl [--out DIR] [--only limit|k0|k3|k6|k8] [--factor 2]")
        end
    end
    # The E5 settings with both mesh resolutions multiplied by `factor`.
    nt, nr, dt, T = 512factor, 128factor, 0.0025, 0.5
    P = Params(kappa=1.0)
    names = vcat("limit", ["k$k" for k in KS])
    isnothing(only) || only in names || error("Unknown run $only; choose one of $(join(names, ", "))")
    source_files = ["src/FreeBoundaryNumerics.jl", "src/materials.jl", "src/geometry.jl",
        "src/fem.jl", "src/flows.jl", "src/spectral.jl", "scripts/run_bulk_fine.jl", "Manifest.toml"]
    hashes = Dict(name => digest(joinpath(ROOT, name)) for name in source_files)
    cache = joinpath(outdir, "bulk_fine_runs")
    mkpath(cache)

    function run_case(name)
        path = joinpath(cache, "$(name)_$(nt)x$(nr).jld2")
        if isfile(path)
            d = load_results(path)
            d["source_sha256"] == hashes || error("Cached source mismatch: $path")
            (d["nt"], d["nr"], d["dt"], d["T"]) == (nt, nr, dt, T) || error("Cached parameters mismatch: $path")
            println("Reusing verified ", path)
            return d
        end
        start = string(now(UTC))
        if name == "limit"
            @printf("limit: %d x %d, dt=%.4f, T=%.1f\n", nt, nr, dt, T); flush(stdout)
            lim = evolve_limit(P, M2, ones(nt); dt=dt, T=T)
            # As in E5: the discrete harmonic potential with the limiting trace on the limit mesh.
            u0 = corrector_data(DirichletSolver(AnnulusMesh(nt, nr; B=P.B)), P, M2, lim.R).u0
            d = Dict{String,Any}("kappa" => 0.0, "R" => lim.R, "u0" => u0, "summary" => lim.summary)
        else
            κ = 2.0^-parse(Int, name[2:end])
            @printf("%s: %d x %d, dt=%.4f, T=%.1f, kappa=%g\n", name, nt, nr, dt, T, κ); flush(stdout)
            run = evolve_coupled(withkappa(P, κ), M2, ones(nt); nr=nr, dt=dt, T=T, store_every=typemax(Int))
            d = Dict{String,Any}("kappa" => κ, "R" => run.R, "phi" => run.phi, "summary" => run.summary)
        end
        merge!(d, Dict{String,Any}("nt" => nt, "nr" => nr, "dt" => dt, "T" => T, "source_sha256" => hashes,
            "started_utc" => start, "finished_utc" => string(now(UTC))))
        save_results(path, d)
        @printf("  %s finished in %.1f s\n", name, d["summary"]["elapsed_seconds"]); flush(stdout)
        return d
    end

    if !isnothing(only)
        run_case(only)
        return
    end
    kappas = [2.0^-k for k in KS]
    record = Dict{String,Any}("description" => "E5 disk bulk fields for fig_bulk.pdf on a finer mesh",
        "nt" => nt, "nr" => nr, "dt" => dt, "T" => T, "mesh_factor" => factor, "kappas" => kappas,
        "material" => M2.name, "i_star" => 1.0, "phieq_formula" => "0.30*(x^2-y^2)+0.15*x*y",
        "beta" => P.beta, "A1" => P.A1, "A2" => P.A2, "wall_radius" => P.B,
        "julia_version" => string(VERSION), "blas_threads" => BLAS.get_num_threads(),
        "source_sha256" => hashes)
    lim = run_case("limit")
    record["limit/R"] = lim["R"]; record["limit/u0"] = lim["u0"]; record["limit/summary"] = lim["summary"]
    for k in KS
        d = run_case("k$k")
        κ = d["kappa"]
        record["$κ/R"] = d["R"]; record["$κ/phi"] = d["phi"]; record["$κ/summary"] = d["summary"]
    end
    destination = joinpath(outdir, "E5_bulk_fine.jld2")
    save_results(destination, record)
    println("Wrote ", destination)
end

main()
