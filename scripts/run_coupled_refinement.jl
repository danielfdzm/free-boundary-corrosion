#!/usr/bin/env julia
# Independent spatial and temporal refinement of the smallest-conductivity E5 disk.
# Archived inputs are read only; all runs use the same exact continuum references.
using FreeBoundaryNumerics, SHA, Printf, LinearAlgebra, Dates

BLAS.set_num_threads(1)

root = normpath(joinpath(@__DIR__, ".."))
datadir = joinpath(root, isdir(joinpath(root, "data")) ? "data" : "results")
outdir = joinpath(root, "outputs", "coupled_refinement")
cachedir = nothing
i = 1
while i <= length(ARGS)
    arg = ARGS[i]
    if arg in ("--help", "-h")
        println("julia --project=. scripts/run_coupled_refinement.jl [--data DIR] [--out DIR] [--cache DIR]")
        println("Refine the E5 disk coupled solve at kappa=2^-8: double both mesh resolutions, then halve and quarter dt on the original mesh.")
        exit(0)
    elseif arg in ("--data", "--out", "--cache")
        i < length(ARGS) || error("$arg requires a directory")
        value = abspath(ARGS[i + 1])
        if arg == "--data"
            global datadir = value
        elseif arg == "--out"
            global outdir = value
        else
            global cachedir = value
        end
        global i += 1
    else
        error("Unknown argument: $arg")
    end
    global i += 1
end

source = joinpath(datadir, "E5.jld2")
destination = joinpath(outdir, "E5_coupled_refinement.jld2")
normpath(source) != normpath(destination) || error("Output must differ from the archived input")
E5 = load_results(source)
get(E5, "i_star", nothing) == 1.0 || error("Input is not the revised constant-current E5 record")
nt, nr, dt, T = (E5[k] for k in ("nt", "nr", "dt", "T"))
@assert (nt, nr, dt, T) == (512, 128, 0.0025, 0.5) "This check uses the production E5 configuration"
κ = E5["disk/kappas"][end]
@assert κ == 2.0^-8
norms(v) = [sobolev_norm(v, m) for m in 0:3]
digest(path) = open(io -> bytes2hex(sha256(io)), path)
Rproduction = E5["disk/$κ/R"]
Rreference, Wreference = E5["disk/limit/R_exact"], E5["disk/limit/W_exact"]
record = Dict{String,Any}(
    "description" => "Independent E5 disk coupled-solution spatial and temporal refinement at kappa=2^-8",
    "nt" => nt, "nr" => nr, "dt" => dt, "T" => T, "kappa" => κ,
    "material" => M2.name, "i_star" => 1.0, "beta" => 0.12, "A1" => 1.0, "A2" => 1.0, "wall_radius" => 2.0,
    "julia_version" => string(VERSION), "blas_threads" => BLAS.get_num_threads(),
    "started_utc" => string(now(UTC)), "E5_sha256" => digest(source),
    "source_sha256" => Dict(name => digest(joinpath(root, name)) for name in
        ("src/FreeBoundaryNumerics.jl", "src/flows.jl", "src/fem.jl", "src/geometry.jl",
         "src/materials.jl", "src/spectral.jl", "scripts/run_coupled_refinement.jl", "Manifest.toml")),
    "comparison_grid" => nt,
    "comparison_convention" => "Fourier project all radii onto the production angular grid before taking differences and norms; use the same exact continuum R0 and W for every corrected remainder",
    "reference_R" => Rreference, "reference_W" => Wreference,
    "production/R" => Rproduction,
    "production/corrected_remainder" => norms(Rproduction .- Rreference .- κ .* Wreference),
    "production/matched_remainder" => E5["disk/rem_norms"][:, end],
    "production/elapsed_seconds" => E5["disk/elapsed"][end],
    "run_names" => ["space_2", "time_2", "time_4"],
)

for (name, angular_factor, radial_factor, time_factor) in
        (("space_2", 2, 2, 1), ("time_2", 1, 1, 2), ("time_4", 1, 1, 4))
    @printf("Coupled refinement %s: %d x %d, dt=%.7f, T=%.2f, kappa=%.8f\n",
        name, angular_factor * nt, radial_factor * nr, dt / time_factor, T, κ)
    flush(stdout)
    if isnothing(cachedir)
        run = evolve_coupled(Params(kappa=κ), M2, ones(angular_factor * nt);
            nr=radial_factor * nr, dt=dt / time_factor, T=T,
            store_every=typemax(Int), verbose=true)
    else
        cache_path = joinpath(cachedir, name * ".jld2")
        cached = load_results(cache_path)
        cached["i_star"] == 1.0 || error("Cached run has different material")
        for (relative, expected) in cached["source_sha256"]
            digest(joinpath(root, relative)) == expected || error("Cached source mismatch: $relative")
        end
        summary = cached["summary"]
        (summary["nt"], summary["nr"], summary["dt"], summary["T"], summary["kappa"], summary["material"]) ==
            (angular_factor*nt, radial_factor*nr, dt/time_factor, T, κ, M2.name) || error("Cached parameters mismatch")
        run = (; R=cached["R"], summary=summary)
        record["$name/cached_run_sha256"] = digest(cache_path)
        record["$name/started_utc"] = cached["started_utc"]
        record["$name/finished_utc"] = cached["finished_utc"]
        record["$name/preparation_script_sha256"] = cached["script_sha256"]
    end
    Rcommon = resample(run.R, nt)
    record["$name/R"] = run.R
    record["$name/R_common"] = Rcommon
    record["$name/summary"] = run.summary
    record["$name/production_difference"] = norms(Rcommon .- Rproduction)
    record["$name/corrected_remainder"] = norms(Rcommon .- Rreference .- κ .* Wreference)
    record["$name/discarded_modes"] = norms(run.R .- resample(Rcommon, length(run.R)))
    if name == "time_2"
        record["$name/previous_step_difference"] = record["$name/production_difference"]
    elseif name == "time_4"
        record["$name/previous_step_difference"] = norms(Rcommon .- record["time_2/R_common"])
    end
    @printf("  completed in %.1f s; production H0 difference %.6e; fixed-reference H0 remainder %.6e\n",
        run.summary["elapsed_seconds"], record["$name/production_difference"][1],
        record["$name/corrected_remainder"][1])
    flush(stdout)
    save_results(destination, record)
end
record["time_orders"] = log2.(record["time_2/previous_step_difference"] ./ record["time_4/previous_step_difference"])
record["finished_utc"] = string(now(UTC))
digest(source) == record["E5_sha256"] || error("Archived E5 input changed during the run")
save_results(destination, record)
include(joinpath(@__DIR__, "coupled_refinement_table.jl"))
write_coupled_refinement_table(record, outdir)
println("Wrote ", destination, "; archived E5.jld2 was not modified.")
