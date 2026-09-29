#!/usr/bin/env julia
# Long-time trajectories of the heterogeneous experiment for its two figures:
# the limiting interface to t = 1 and the coupled interface at kappa = 1/4 to
# t = 2.5. The limiting flow loses smoothness near t = 2.93, when a corner forms
# at the sharpest convex lobe; the coupled trajectory at kappa = 1/4 is no longer
# resolved by 512 Fourier modes after t = 2.6. The sweep record is read only.
using FreeBoundaryNumerics, LinearAlgebra, Printf, SHA, Dates
BLAS.set_num_threads(1)

const ROOT = normpath(joinpath(@__DIR__, ".."))
const P = Params(A2=1.5, A1=0.5, beta=0.12, B=2.0)
const MAT = M2family(1.0)
initial_radius(th) = 1 + 0.12cos(3th) + 0.06sin(2th) + 0.04cos(5th + 0.3)
digest(path) = open(io -> bytes2hex(sha256(io)), path)

function main()
    outdir = joinpath(ROOT, "outputs", "heterogeneous")
    if length(ARGS) == 2 && ARGS[1] == "--out"
        outdir = abspath(ARGS[2])
    elseif !isempty(ARGS)
        error("Usage: run_heterogeneous_long.jl [--out DIR]")
    end
    nt, nr, dt = 512, 128, 0.005
    T_limit, T_coupled, kappa = 1.0, 2.5, 0.25
    source_files = ["src/FreeBoundaryNumerics.jl", "src/materials.jl", "src/geometry.jl",
        "src/fem.jl", "src/flows.jl", "src/spectral.jl", "scripts/run_heterogeneous_long.jl", "Manifest.toml"]
    hashes = Dict(name => digest(joinpath(ROOT, name)) for name in source_files)
    Rinit = initial_radius.(theta_grid(nt))
    record = Dict{String,Any}("description"=>"Long-time trajectories of the heterogeneous experiment",
        "nt"=>nt, "nr"=>nr, "dt"=>dt, "theta"=>theta_grid(nt), "R_initial"=>Rinit,
        "i0_formula"=>"1+0.25*x+0.10*y^2", "A2"=>P.A2, "A1"=>P.A1, "beta"=>P.beta, "B"=>P.B,
        "julia_version"=>string(VERSION), "blas_threads"=>BLAS.get_num_threads(),
        "source_sha256"=>hashes, "started_utc"=>string(now(UTC)))

    @printf("limit: %d nodes, dt=%.3f, T=%.1f\n", nt, dt, T_limit); flush(stdout)
    lim = evolve_limit(P, MAT, Rinit; dt=dt, T=T_limit, store_every=1)
    record["limit"] = Dict{String,Any}("T"=>T_limit, "kappa"=>0.0, "R"=>lim.R, "Rhist"=>lim.Rhist,
        "times"=>lim.saved_times, "history"=>lim.history, "summary"=>lim.summary)

    @printf("coupled: %d x %d, dt=%.3f, T=%.1f, kappa=%g\n", nt, nr, dt, T_coupled, kappa); flush(stdout)
    run = evolve_coupled(withkappa(P, kappa), MAT, Rinit; nr=nr, dt=dt, T=T_coupled, store_every=1)
    record["coupled"] = Dict{String,Any}("T"=>T_coupled, "kappa"=>kappa, "R"=>run.R, "Rhist"=>run.Rhist,
        "times"=>run.saved_times, "history"=>run.history, "summary"=>run.summary)
    @printf("  finished in %.1f s; final maximal curvature %.2f\n",
        run.summary["elapsed_seconds"], run.history["max_curvature"][end])

    record["finished_utc"] = string(now(UTC))
    path = joinpath(outdir, "long_time.jld2")
    save_results(path, record)
    println("Wrote ", path)
end

main()
