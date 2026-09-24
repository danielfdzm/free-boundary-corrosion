#!/usr/bin/env julia
# Verify the archived artifacts and the exact homogeneous benchmark.
using FreeBoundaryNumerics, SHA, Printf

root = normpath(joinpath(@__DIR__, ".."))
for line in eachline(joinpath(root, "checksums.sha256"))
    expected, relative = split(line; limit=2)
    path = joinpath(root, strip(relative))
    actual = open(path) do io
        bytes2hex(sha256(io))
    end
    actual == expected || error("Checksum mismatch: $relative")
end
println("Archived data, figures, previews and tables match their checksums.")

for name in ("E1", "E2", "E3", "E4", "E5", "E6", "E6_phase")
    record = load_results(joinpath(root, "data", name * ".jld2"))
    isempty(record) && error("Empty record: $name")
    println("Opened $name: $(length(record)) stored fields")
end

P = Params(kappa=0.5)
run = evolve_coupled(P, M0, ones(48); nr=12, dt=0.02, T=1.0)
radius_error = maximum(abs, run.R .- (1 - P.beta))
potential_error = maximum(abs, run.phi .- 0.15)
radius_error < 1e-11 || error("Homogeneous radius error: $radius_error")
potential_error < 1e-11 || error("Homogeneous potential error: $potential_error")
run.summary["max_charge_residual"] < 1e-10 || error("Charge residual exceeds tolerance")
@printf("Exact shrinking disk: radius error %.3e, potential error %.3e\n", radius_error, potential_error)
