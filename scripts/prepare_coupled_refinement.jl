#!/usr/bin/env julia
# Compute independent refinement solves while the production E5 sweep runs.
using FreeBoundaryNumerics, SHA, LinearAlgebra, Printf, Dates
BLAS.set_num_threads(1)
root = normpath(joinpath(@__DIR__, ".."))
outdir = isempty(ARGS) ? joinpath(root, "outputs", "constant_current_refinement_runs") : abspath(ARGS[1])
mkpath(outdir)
source_names = ("src/FreeBoundaryNumerics.jl", "src/flows.jl", "src/fem.jl", "src/geometry.jl", "src/materials.jl", "src/spectral.jl", "Manifest.toml")
digest(path) = open(io -> bytes2hex(sha256(io)), path)
source_hashes = Dict(name => digest(joinpath(root, name)) for name in source_names)
for (name, nt, nr, dt) in (("space_2", 1024, 256, 0.0025), ("time_2", 512, 128, 0.00125), ("time_4", 512, 128, 0.000625))
    @printf("Preparing %s: %d x %d, dt=%.7f\n", name, nt, nr, dt)
    flush(stdout)
    started = string(now(UTC))
    run = evolve_coupled(Params(kappa=2.0^-8), M2, ones(nt); nr=nr, dt=dt, T=0.5, store_every=typemax(Int), verbose=true)
    record = Dict{String,Any}("R"=>run.R, "summary"=>run.summary, "source_sha256"=>source_hashes, "i_star"=>1.0,
        "started_utc"=>started, "finished_utc"=>string(now(UTC)), "julia_version"=>string(VERSION),
        "script_sha256"=>digest(@__FILE__))
    save_results(joinpath(outdir, name*".jld2"), record)
    flush(stdout)
end
