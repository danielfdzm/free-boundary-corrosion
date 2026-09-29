#!/usr/bin/env julia
# Verify current artifacts, stored numerical records, and the homogeneous benchmark.
using FreeBoundaryNumerics, SHA, Printf, LinearAlgebra
BLAS.set_num_threads(1)
root = normpath(joinpath(@__DIR__, ".."))
for manifest in ("checksums.sha256", "data/heterogeneous/checksums.sha256")
    for line in eachline(joinpath(root, manifest))
        isempty(strip(line)) && continue
        expected, relative = split(line; limit=2)
        actual = open(io -> bytes2hex(sha256(io)), joinpath(root, strip(relative)))
        actual == expected || error("Checksum mismatch in $manifest: $relative")
    end
end
println("Disk and heterogeneous artifacts match their checksum manifests.")
E5 = load_results(joinpath(root, "data", "E5.jld2"))
E5["i_star"] == 1.0 || error("E5 does not use the revised constant corrosion-current scale")
E5["material"] == M2.name || error("E5 material differs from the current source")
all(M2.i0(x,y) == 1.0 for x in (-1.0,0.0,1.0), y in (-1.0,0.0,1.0)) || error("M2 is not constant-current")
nt, T = E5["nt"], E5["T"]
r = 1.0 - 0.12*T
q = 0.30 .* cos.(2 .* E5["theta"]) .+ 0.075 .* sin.(2 .* E5["theta"])
primitive(s) = -s^2/2 + 4*atan(s^2/4)
maximum(abs, E5["disk/limit/R_exact"] .- r) < 1e-14 || error("Incorrect exact limiting radius")
maximum(abs, E5["disk/limit/W_exact"] .- (primitive(1.0)-primitive(r)).*q) < 1e-14 || error("Incorrect exact corrector")
# Charge is an aggregate of nodal residuals, each solved to 2e-12.
maximum(E5["disk/charge"]) < 1e-9 || error("E5 charge residual exceeds tolerance")
all(>(0), E5["disk/excess"]) || error("Dissolution excess is not positive")
println("E5 uses a constant corrosion-current scale and the exact continuum references.")

coupled = load_results(joinpath(root, "data", "E5_coupled_refinement.jld2"))
source_digest = open(io -> bytes2hex(sha256(io)), joinpath(root, "data", "E5.jld2"))
coupled["E5_sha256"] == source_digest || error("Coupled refinement uses a different E5 record")
coupled["i_star"] == 1.0 || error("Coupled refinement uses a different material")
all(isfinite, coupled["time_orders"]) || error("Incomplete coupled time refinement")
for (relative, expected) in coupled["source_sha256"]
    actual = open(io -> bytes2hex(sha256(io)), joinpath(root, relative))
    actual == expected || error("Coupled-refinement source differs from the recorded version: $relative")
end
for name in coupled["run_names"]
    summary = coupled["$name/summary"]
    summary["T_final"] == coupled["T"] || error("Incomplete coupled run: $name")
    summary["max_newton_residual"] < 1e-10 || error("Coupled Newton residual exceeds tolerance: $name")
    summary["max_charge_residual"] < 1e-9 || error("Coupled charge residual exceeds tolerance: $name")
    radius = coupled["$name/R_common"]
    difference = radius .- coupled["production/R"]
    remainder = radius .- coupled["reference_R"] .- coupled["kappa"] .* coupled["reference_W"]
    for m in 0:3
        isapprox(sobolev_norm(difference, m), coupled["$name/production_difference"][m + 1]; rtol=1e-12) || error("Inconsistent coupled difference: $name, H$m")
        isapprox(sobolev_norm(remainder, m), coupled["$name/corrected_remainder"][m + 1]; rtol=1e-12) || error("Inconsistent coupled remainder: $name, H$m")
    end
end
println("Coupled refinement matches its production input and source; all stored norms are consistent.")

bulk = load_results(joinpath(root, "data", "E5_bulk_fine.jld2"))
(bulk["nt"], bulk["nr"], bulk["dt"], bulk["T"]) == (2nt, 2 * E5["nr"], E5["dt"], T) || error("Bulk record does not double the E5 mesh")
bulk["i_star"] == 1.0 && bulk["material"] == M2.name || error("Bulk record uses a different material")
bulk["kappas"] == [2.0^-k for k in (0, 3, 6, 8)] || error("Incomplete bulk conductivity set")
for (relative, expected) in bulk["source_sha256"]
    actual = open(io -> bytes2hex(sha256(io)), joinpath(root, relative))
    actual == expected || error("Bulk source differs from the recorded version: $relative")
end
bulk["limit/summary"]["T_final"] == T || error("Incomplete bulk limiting flow")
for κ in bulk["kappas"]
    summary = bulk["$κ/summary"]
    summary["T_final"] == T && summary["kappa"] == κ || error("Incomplete bulk run: kappa=$κ")
    summary["max_newton_residual"] < 1e-10 || error("Bulk Newton residual exceeds tolerance: kappa=$κ")
    summary["max_charge_residual"] < 1e-9 || error("Bulk charge residual exceeds tolerance: kappa=$κ")
    length(bulk["$κ/phi"]) == length(bulk["limit/u0"]) == bulk["nt"] * (bulk["nr"] + 1) || error("Invalid bulk field: kappa=$κ")
end
# At kappa = 2^-8 the bulk run repeats the doubled-mesh coupled refinement run.
bulk["$(2.0^-8)/R"] == coupled["space_2/R"] || error("Bulk run differs from the doubled-mesh refinement run")
println("Bulk record doubles the E5 mesh; at kappa=2^-8 it repeats the doubled-mesh refinement run.")

function check_heterogeneous(root)
    d = load_results(joinpath(root, "data", "heterogeneous", "heterogeneous.jld2"))
    nt, nr, dt, T = (d[key] for key in ("nt", "nr", "dt", "T"))
    (nt, nr, dt, T) == (512, 128, 0.005, 0.5) || error("Heterogeneous record is not the production run")
    (d["A2"], d["A1"], d["beta"], d["B"]) == (1.5, 0.5, 0.12, 2.0) || error("Incorrect heterogeneous parameters")
    kappas = d["kappas"]
    kappas == 2.0 .^ (-(2:6)) || error("Incomplete heterogeneous conductivity sweep")
    d["theta"] == theta_grid(nt) || error("Incorrect heterogeneous angular grid")
    initial_radius(th) = 1 + 0.12cos(3th) + 0.06sin(2th) + 0.04cos(5th + 0.3)
    isapprox(d["R_initial"], initial_radius.(d["theta"]); rtol=1e-13) || error("Incorrect heterogeneous initial interface")
    for (relative, expected) in d["source_sha256"]
        actual = open(io -> bytes2hex(sha256(io)), joinpath(root, relative))
        actual == expected || error("Heterogeneous source differs from the recorded version: $relative")
    end

    function check_run(run, label, expected; limit=false)
        (run["nt"], run["nr"], run["dt"], run["kappa"]) == expected || error("Incorrect run parameters: $label")
        run["source_sha256"] == d["source_sha256"] || error("Inconsistent run provenance: $label")
        run["T"] == T == run["summary"]["T_final"] || error("Incomplete heterogeneous run: $label")
        times, radii = run["times"], run["Rhist"]
        first(times) == 0 && last(times) == T && all(>(0), diff(times)) || error("Invalid saved times: $label")
        size(radii) == (run["nt"], length(times)) || error("Invalid radius history: $label")
        all(r -> isfinite(r) && 0 < r < d["B"], radii) || error("Invalid interface geometry: $label")
        radii[:, end] == run["R"] || error("Final radius differs from its history: $label")
        if limit
            size(run["Whist"]) == size(radii) && all(isfinite, run["Whist"]) || error("Invalid corrector history: $label")
            run["Whist"][:, end] == run["W"] || error("Final corrector differs from its history: $label")
        else
            run["summary"]["max_newton_residual"] < 1e-10 || error("Heterogeneous Newton residual exceeds tolerance: $label")
            run["summary"]["max_charge_residual"] < 1e-9 || error("Heterogeneous charge residual exceeds tolerance: $label")
            run["summary"]["min_excess"] > 0 || error("Heterogeneous dissolution excess is not positive: $label")
        end
    end

    hnorms(v) = [sobolev_norm(v, m) for m in 0:1]
    consistent(actual, expected, label) = isapprox(actual, expected; rtol=1e-11, atol=1e-13) || error("Inconsistent heterogeneous $label")
    check_norms(v, key) = consistent(hnorms(v), d[key], key)
    base, reference = d["limit"], d["reference"]
    check_run(base, "limit", (nt, nr, dt, 0.0); limit=true)
    check_run(reference, "reference", (2nt, 2nr, dt/2, 0.0); limit=true)
    Rref, Wref = resample(reference["R"], nt), resample(reference["W"], nt)
    length(d["coupled"]) == length(kappas) || error("Incomplete heterogeneous coupled records")
    for (j, kappa) in enumerate(kappas)
        run = d["coupled"][j]
        check_run(run, "kappa=$kappa", (nt, nr, dt, kappa))
        for (field, difference) in (
            ("deviation", run["R"] - Rref),
            ("remainder", run["R"] - Rref - kappa * Wref),
            ("matched_deviation", run["R"] - base["R"]),
            ("matched_remainder", run["R"] - base["R"] - kappa * base["W"]))
            consistent(hnorms(difference), d["errors/$field"][:, j], "$field at kappa=$kappa")
        end
    end
    for field in ("deviation", "remainder"), m in 1:2
        orders = observed_orders(kappas, d["errors/$field"][m, :])
        consistent(orders[1:end-1], d["orders/$field"][m, 1:end-1], "$field orders in H$(m-1)")
        isnan(d["orders/$field"][m, end]) || error("Missing final-order sentinel: $field")
    end

    kappa = kappas[end]
    for (name, n, r, step) in (("space_2", 2nt, 2nr, dt), ("time_2", nt, nr, dt/2), ("time_4", nt, nr, dt/4))
        lr, cr = d["refinement/limit/$name"], d["refinement/coupled/$name"]
        check_run(lr, "limit/$name", (n, r, step, 0.0); limit=true)
        check_run(cr, "coupled/$name", (n, r, step, kappa))
        check_norms(resample(cr["R"], nt) - d["coupled"][end]["R"], "sensitivity/$name/coupled")
        check_norms(resample(lr["R"], nt) - base["R"], "sensitivity/$name/limit")
        consistent(kappa * hnorms(resample(lr["W"], nt) - base["W"]), d["sensitivity/$name/corrector_weighted"], "$name corrector sensitivity")
        check_norms(resample(cr["R"], nt) - Rref - kappa * Wref, "sensitivity/$name/fixed_reference_remainder")
    end
    for kind in ("limit", "coupled")
        r2, r4 = d["refinement/$kind/time_2"], d["refinement/$kind/time_4"]
        orders = log2.(d["sensitivity/time_2/$kind"] ./ hnorms(r2["R"] - r4["R"]))
        consistent(orders, d["sensitivity/time_order/$kind"], "$kind temporal orders")
    end
    r2, r4 = d["refinement/limit/time_2"], d["refinement/limit/time_4"]
    orders = log2.(d["sensitivity/time_2/corrector_weighted"] ./ (kappa * hnorms(r2["W"] - r4["W"])))
    consistent(orders, d["sensitivity/time_order/corrector"], "corrector temporal orders")
    check_norms(Rref - base["R"], "sensitivity/reference/limit")
    consistent(kappa * hnorms(Wref - base["W"]), d["sensitivity/reference/corrector_weighted"], "reference corrector sensitivity")
    check_norms(base["W"] - d["forcing_only/W"], "forcing_only/difference")
    consistent(hnorms(base["W"] - d["forcing_only/W"]) ./ hnorms(base["W"]), d["forcing_only/relative_difference"], "forcing-only relative difference")
    println("Heterogeneous records match their numerical sources; saved histories, solver diagnostics, and H0/H1 norms are consistent.")
end

check_heterogeneous(root)

function check_long_time(root)
    dir = joinpath(root, "data", "heterogeneous")
    d = load_results(joinpath(dir, "heterogeneous.jld2"))
    L = load_results(joinpath(dir, "long_time.jld2"))
    (L["nt"], L["nr"], L["dt"]) == (d["nt"], d["nr"], d["dt"]) || error("Long-time record uses a different discretization")
    (L["A2"], L["A1"], L["beta"], L["B"]) == (d["A2"], d["A1"], d["beta"], d["B"]) || error("Long-time record uses different parameters")
    L["theta"] == d["theta"] && L["R_initial"] == d["R_initial"] || error("Long-time record uses a different initial interface")
    for (relative, expected) in L["source_sha256"]
        actual = open(io -> bytes2hex(sha256(io)), joinpath(root, relative))
        actual == expected || error("Long-time source differs from the recorded version: $relative")
    end
    L["coupled"]["kappa"] == d["coupled"][1]["kappa"] == 0.25 || error("Long-time coupled run is not at kappa=1/4")
    for (key, T, sweep) in (("limit", 1.0, d["limit"]), ("coupled", 2.5, d["coupled"][1]))
        run = L[key]
        run["T"] == T == run["summary"]["T_final"] || error("Incomplete long-time run: $key")
        times, radii = run["times"], run["Rhist"]
        first(times) == 0 && last(times) == T && all(>(0), diff(times)) || error("Invalid saved times: long-time $key")
        size(radii) == (L["nt"], length(times)) || error("Invalid radius history: long-time $key")
        all(r -> isfinite(r) && 0 < r < L["B"], radii) || error("Invalid interface geometry: long-time $key")
        radii[:, end] == run["R"] || error("Final radius differs from its history: long-time $key")
        # Up to t = 0.5 the long runs repeat the stored sweep runs.
        n = length(sweep["times"])
        times[1:n] == sweep["times"] && maximum(abs, radii[:, 1:n] .- sweep["Rhist"]) <= 1e-12 ||
            error("Long-time $key run departs from the stored sweep run")
    end
    summary = L["coupled"]["summary"]
    summary["max_newton_residual"] < 1e-10 || error("Long-time Newton residual exceeds tolerance")
    summary["max_charge_residual"] < 1e-9 || error("Long-time charge residual exceeds tolerance")
    summary["min_excess"] > 0 || error("Long-time dissolution excess is not positive")
    println("Long-time trajectories match their numerical sources and repeat the stored sweep runs up to t = 0.5.")
end

check_long_time(root)

P = Params(kappa=0.5)
run = evolve_coupled(P, M0, ones(48); nr=12, dt=0.02, T=1.0)
radius_error = maximum(abs, run.R .- (1-P.beta))
potential_error = maximum(abs, run.phi .- 0.15)
radius_error < 1e-11 || error("Homogeneous radius error: $radius_error")
potential_error < 1e-11 || error("Homogeneous potential error: $potential_error")
run.summary["max_charge_residual"] < 1e-9 || error("Charge residual exceeds tolerance")
@printf("Exact shrinking disk: radius error %.3e, potential error %.3e\n", radius_error, potential_error)
