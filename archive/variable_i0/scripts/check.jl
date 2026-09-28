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

for name in ("E1", "E2", "E3", "E4", "E5", "E5_refinement", "E5_coupled_refinement")
    record = load_results(joinpath(root, "data", name * ".jld2"))
    isempty(record) && error("Empty record: $name")
    println("Opened $name: $(length(record)) stored fields")
end

refinement = load_results(joinpath(root, "data", "E5_refinement.jld2"))
source_digest = open(io -> bytes2hex(sha256(io)), joinpath(root, "data", "E5.jld2"))
refinement["E5_sha256"] == source_digest || error("Refinement supplement uses a different E5 record")
all(isfinite, refinement["time_orders"]) || error("Incomplete time-refinement supplement")
println("Refinement supplement matches its production E5 input and includes both time halvings.")

coupled = load_results(joinpath(root, "data", "E5_coupled_refinement.jld2"))
coupled["E5_sha256"] == source_digest || error("Coupled refinement uses a different E5 record")
all(isfinite, coupled["time_orders"]) || error("Incomplete coupled time refinement")
for (relative, expected) in coupled["source_sha256"]
    actual = open(io -> bytes2hex(sha256(io)), joinpath(root, relative))
    actual == expected || error("Coupled-refinement source differs from the recorded version: $relative")
end
for name in coupled["run_names"]
    summary = coupled["$name/summary"]
    summary["T_final"] == coupled["T"] || error("Incomplete coupled run: $name")
    summary["max_newton_residual"] < 1e-10 || error("Coupled Newton residual exceeds tolerance: $name")
    summary["max_charge_residual"] < 1e-10 || error("Coupled charge residual exceeds tolerance: $name")
    radius = coupled["$name/R_common"]
    difference = radius .- coupled["production/R"]
    remainder = radius .- coupled["reference_R"] .- coupled["kappa"] .* coupled["reference_W"]
    for m in 0:3
        isapprox(sobolev_norm(difference, m), coupled["$name/production_difference"][m + 1]; rtol=1e-12) ||
            error("Inconsistent coupled difference: $name, H$m")
        isapprox(sobolev_norm(remainder, m), coupled["$name/corrected_remainder"][m + 1]; rtol=1e-12) ||
            error("Inconsistent coupled remainder: $name, H$m")
    end
    ratios = coupled["$name/production_difference"] ./ coupled["production/corrected_remainder"]
    all(r -> isfinite(r) && r < 0.1, ratios) || error("Coupled refinement exceeds 10% of the corrected remainder: $name")
end
coupled["space_2/discarded_modes"][4] < 0.02 * coupled["production/corrected_remainder"][4] ||
    error("Discarded spatial modes exceed 2% of the corrected H3 remainder")
println("Coupled refinement matches its production input and source; all three runs and stored norms are consistent.")

P = Params(kappa=0.5)
run = evolve_coupled(P, M0, ones(48); nr=12, dt=0.02, T=1.0)
radius_error = maximum(abs, run.R .- (1 - P.beta))
potential_error = maximum(abs, run.phi .- 0.15)
radius_error < 1e-11 || error("Homogeneous radius error: $radius_error")
potential_error < 1e-11 || error("Homogeneous potential error: $potential_error")
run.summary["max_charge_residual"] < 1e-10 || error("Charge residual exceeds tolerance")
@printf("Exact shrinking disk: radius error %.3e, potential error %.3e\n", radius_error, potential_error)
