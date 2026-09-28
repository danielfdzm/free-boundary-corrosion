#!/usr/bin/env julia
# Independent time refinement of the production E5 disk corrector.
# Archived E5 records are read only; new records go to outputs/refinement/.
using FreeBoundaryNumerics, SHA, Printf, LinearAlgebra

BLAS.set_num_threads(1)

root = normpath(joinpath(@__DIR__, ".."))
datadir = joinpath(root, isdir(joinpath(root, "data")) ? "data" : "results")
outdir = joinpath(root, "outputs", "refinement")
i = 1
while i <= length(ARGS)
    arg = ARGS[i]
    if arg in ("--help", "-h")
        println("julia --project=. scripts/run_refinement.jl [--data DIR] [--out DIR]")
        println("Run the E5 disk corrector at the production spatial resolution with dt/2 and dt/4.")
        exit(0)
    elseif arg in ("--data", "--out")
        i < length(ARGS) || error("$arg requires a directory")
        value = abspath(ARGS[i + 1])
        arg == "--data" ? (global datadir = value) : (global outdir = value)
        global i += 1
    else
        error("Unknown argument: $arg")
    end
    global i += 1
end

source = joinpath(datadir, "E5.jld2")
destination = joinpath(outdir, "E5_refinement.jld2")
normpath(source) != normpath(destination) || error("Output must differ from the archived input")
E5 = load_results(source)
get(E5, "i_star", nothing) == 1.0 || error("Input is not the revised constant-current E5 record")
nt, nr, dt, T = (E5[k] for k in ("nt", "nr", "dt", "T"))
@assert (nt, nr, dt, T) == (512, 128, 0.0025, 0.5) "This check uses the production E5 configuration"
κ = E5["disk/kappas"][end]
norms(v) = [sobolev_norm(v, m) for m in 0:3]
digest(path) = open(io -> bytes2hex(sha256(io)), path)
record = Dict{String,Any}(
    "description" => "E5 disk corrector time refinement at fixed production spatial resolution",
    "nt" => nt, "nr" => nr, "dt" => dt, "T" => T, "kappa" => κ,
    "material" => M2.name, "i_star" => 1.0, "beta" => 0.12, "A1" => 1.0, "A2" => 1.0, "wall_radius" => 2.0,
    "julia_version" => string(VERSION), "blas_threads" => BLAS.get_num_threads(),
    "E5_sha256" => digest(source),
    "source_sha256" => Dict(name => digest(joinpath(root, name)) for name in
        ("src/flows.jl", "src/fem.jl", "src/materials.jl", "src/spectral.jl", "scripts/run_refinement.jl", "Manifest.toml")),
    "limit_reference_difference" => E5["disk/limit/R_matched_minus_fine"],
    "corrector_spatial_difference" => E5["disk/limit/W_matched_minus_fine"],
    "matched_remainder" => E5["disk/rem_norms"][:, end],
    "refined_remainder" => E5["disk/rem_norms_fine"][:, end],
    "dt_factors" => [2, 4],
)
previous = E5["disk/limit/W"]
for factor in record["dt_factors"]
    @printf("E5 corrector time refinement: %d x %d, dt=%.7f, T=%.2f\n", nt, nr, dt / factor, T)
    flush(stdout)
    run = evolve_limit(Params(kappa=1.0), M2, ones(nt); nr=nr, dt=dt / factor, T=T,
        corrector=true, store_every=typemax(Int))
    record["factor_$factor/R"] = run.R
    record["factor_$factor/W"] = run.W
    record["factor_$factor/elapsed_seconds"] = run.summary["elapsed_seconds"]
    record["factor_$factor/previous_step_difference"] = norms(previous .- run.W)
    record["factor_$factor/production_step_difference"] = norms(E5["disk/limit/W"] .- run.W)
    @printf("  completed in %.1f s; kappa-weighted consecutive H0 difference %.6e\n",
        run.summary["elapsed_seconds"], κ * record["factor_$factor/previous_step_difference"][1])
    flush(stdout)
    global previous = run.W
    save_results(destination, record)
end
record["time_orders"] = log2.(record["factor_2/previous_step_difference"] ./ record["factor_4/previous_step_difference"])
save_results(destination, record)
include(joinpath(@__DIR__, "refinement_table.jl"))
write_refinement_table(record, outdir)
println("Wrote ", destination, "; archived E5.jld2 was not modified.")
