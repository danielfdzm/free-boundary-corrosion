#!/usr/bin/env julia
# Standalone heterogeneous/noncircular experiment; existing E5 files are read only.
using FreeBoundaryNumerics, LinearAlgebra, Printf, SHA, Dates
BLAS.set_num_threads(1)

const ROOT = normpath(joinpath(@__DIR__, ".."))
const P = Params(A2=1.5, A1=0.5, beta=0.12, B=2.0)
const MAT = M2family(1.0)
initial_radius(th) = 1 + 0.12cos(3th) + 0.06sin(2th) + 0.04cos(5th + 0.3)
norms(v) = [sobolev_norm(v, m) for m in 0:1]
digest(path) = open(io -> bytes2hex(sha256(io)), path)

function main()
    nt = 512
    outdir = joinpath(ROOT, "outputs", "heterogeneous")
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--out"
            outdir = abspath(ARGS[i + 1]); i += 2
        elseif ARGS[i] == "--nt"
            nt = parse(Int, ARGS[i + 1]); i += 2
        else
            error("Usage: run_heterogeneous.jl [--out DIR] [--nt 512]")
        end
    end
    nr, dt, T = nt ÷ 4, 0.005, 0.5
    kappas = 2.0 .^ (-(2:6))
    source_files = ["src/FreeBoundaryNumerics.jl", "src/materials.jl", "src/geometry.jl",
        "src/fem.jl", "src/flows.jl", "src/spectral.jl", "scripts/run_heterogeneous.jl", "Manifest.toml"]
    hashes = Dict(name => digest(joinpath(ROOT, name)) for name in source_files)
    cache = joinpath(outdir, "runs")
    mkpath(cache)
    function run_case(kind, name, n, r, step; kappa=0.0)
        path = joinpath(cache, name * ".jld2")
        if isfile(path)
            d = load_results(path)
            d["source_sha256"] == hashes || error("Cached source mismatch: $path")
            (d["nt"], d["nr"], d["dt"], d["T"], d["kappa"]) == (n, r, step, T, kappa) ||
                error("Cached parameters mismatch: $path")
            println("Reusing verified ", name)
            return d
        end
        @printf("%s: %d x %d, dt=%.7f, kappa=%g\n", name, n, r, step, kappa)
        flush(stdout)
        start = string(now(UTC))
        Rinit = initial_radius.(theta_grid(n))
        run = kind == :limit ? evolve_limit(P, MAT, Rinit; nr=r, dt=step, T=T,
                corrector=true, store_every=1) :
            evolve_coupled(withkappa(P, kappa), MAT, Rinit; nr=r, dt=step, T=T,
                store_every=1)
        d = Dict{String,Any}("nt"=>n, "nr"=>r, "dt"=>step, "T"=>T, "kappa"=>kappa,
            "R"=>run.R, "Rhist"=>run.Rhist, "times"=>run.saved_times,
            "summary"=>run.summary, "history"=>run.history, "source_sha256"=>hashes,
            "started_utc"=>start, "finished_utc"=>string(now(UTC)))
        if kind == :limit
            d["W"] = run.W
            d["Whist"] = run.Whist
        end
        save_results(path, d)
        @printf("  finished in %.2f s\n", run.summary["elapsed_seconds"])
        flush(stdout)
        return d
    end

    record = Dict{String,Any}("description"=>"Asymmetric lobed interface, heterogeneous corrosion current, unequal reaction slopes",
        "nt"=>nt, "nr"=>nr, "dt"=>dt, "T"=>T, "kappas"=>kappas,
        "theta"=>theta_grid(nt), "R_initial"=>initial_radius.(theta_grid(nt)),
        "R_initial_formula"=>"1+0.12*cos(3*theta)+0.06*sin(2*theta)+0.04*cos(5*theta+0.3)",
        "i0_formula"=>"1+0.25*x+0.10*y^2", "i0_min_container"=>0.5,
        "phieq_formula"=>"0.30*(x^2-y^2)+0.15*x*y",
        "A2"=>P.A2, "A1"=>P.A1, "beta"=>P.beta, "B"=>P.B,
        "julia_version"=>string(VERSION), "blas_threads"=>BLAS.get_num_threads(),
        "source_sha256"=>hashes, "started_utc"=>string(now(UTC)),
        "reference_convention"=>"Numerical refined reference: 2nt x 2nr, dt/2. Fourier projection onto nt grid for every norm; radial corrector W, not normal displacement.")
    base = run_case(:limit, "limit", nt, nr, dt)
    record["limit"] = base
    coupled = []
    for (j, kappa) in enumerate(kappas)
        push!(coupled, run_case(:coupled, "coupled_k$(j+1)", nt, nr, dt; kappa))
    end
    record["coupled"] = coupled
    for (name, n, r, step) in (("space_2", 2nt, 2nr, dt),
            ("time_2", nt, nr, dt/2), ("time_4", nt, nr, dt/4))
        record["refinement/limit/$name"] = run_case(:limit, "limit_" * name, n, r, step)
        record["refinement/coupled/$name"] = run_case(:coupled, "coupled_" * name, n, r, step; kappa=kappas[end])
    end
    reference = run_case(:limit, "limit_reference", 2nt, 2nr, dt/2)
    record["reference"] = reference
    Rref, Wref = resample(reference["R"], nt), resample(reference["W"], nt)
    dev = hcat([norms(d["R"] - Rref) for d in coupled]...)
    rem = hcat([norms(d["R"] - Rref - k * Wref) for (d, k) in zip(coupled, kappas)]...)
    record["errors/deviation"] = dev
    record["errors/remainder"] = rem
    record["errors/matched_deviation"] = hcat([norms(d["R"] - base["R"]) for d in coupled]...)
    record["errors/matched_remainder"] = hcat([norms(d["R"] - base["R"] - k * base["W"]) for (d, k) in zip(coupled, kappas)]...)
    record["orders/deviation"] = hcat([observed_orders(kappas, dev[m, :]) for m in 1:2]...)'
    record["orders/remainder"] = hcat([observed_orders(kappas, rem[m, :]) for m in 1:2]...)'
    for name in ("space_2", "time_2", "time_4")
        lr, cr = record["refinement/limit/$name"], record["refinement/coupled/$name"]
        record["sensitivity/$name/coupled"] = norms(resample(cr["R"], nt) - coupled[end]["R"])
        record["sensitivity/$name/limit"] = norms(resample(lr["R"], nt) - base["R"])
        record["sensitivity/$name/corrector_weighted"] = kappas[end] * norms(resample(lr["W"], nt) - base["W"])
        record["sensitivity/$name/fixed_reference_remainder"] = norms(resample(cr["R"], nt) - Rref - kappas[end] * Wref)
    end
    for kind in ("limit", "coupled")
        r2, r4 = record["refinement/$kind/time_2"], record["refinement/$kind/time_4"]
        record["sensitivity/time_order/$kind"] = log2.(record["sensitivity/time_2/$kind"] ./ norms(r2["R"] - r4["R"]))
    end
    r2, r4 = record["refinement/limit/time_2"], record["refinement/limit/time_4"]
    record["sensitivity/time_order/corrector"] = log2.(record["sensitivity/time_2/corrector_weighted"] ./ (kappas[end] * norms(r2["W"] - r4["W"])))
    record["sensitivity/reference/limit"] = norms(Rref - base["R"])
    record["sensitivity/reference/corrector_weighted"] = kappas[end] * norms(Wref - base["W"])

    # Diagnostic: integrate only the electrical forcing along the SAME limiting
    # Heun trajectory. Comparing this with W isolates the omitted transport and
    # zeroth-order material/geometry terms, without changing the reference flow.
    D = DirichletSolver(AnnulusMesh(nt, nr; B=P.B))
    th = theta_grid(nt); cs, sn = cos.(th), sin.(th)
    forcing(R) = begin
        cd = corrector_data(D, P, MAT, R)
        -P.beta * P.A2 .* metric_factor(R) .* cd.i0b .* cd.gamma1
    end
    Wforcing = zeros(nt)
    for j in 1:length(base["times"])-1
        R = base["Rhist"][:, j]
        Rp = R + dt * FreeBoundaryNumerics.limit_speed(P, MAT, R, cs, sn)
        Wforcing .+= (dt/2) .* (forcing(R) + forcing(Rp))
    end
    record["forcing_only/W"] = Wforcing
    record["forcing_only/difference"] = norms(base["W"] - Wforcing)
    record["forcing_only/relative_difference"] = norms(base["W"] - Wforcing) ./ norms(base["W"])
    record["finished_utc"] = string(now(UTC))
    save_results(joinpath(outdir, "heterogeneous.jld2"), record)
    open(joinpath(outdir, "convergence.csv"), "w") do io
        println(io, "kappa,deviation_H0,remainder_H0,deviation_H1,remainder_H1,matched_remainder_H0,matched_remainder_H1")
        for j in eachindex(kappas)
            println(io, join([kappas[j], dev[1,j], rem[1,j], dev[2,j], rem[2,j], record["errors/matched_remainder"][1,j], record["errors/matched_remainder"][2,j]], ','))
        end
    end
    open(joinpath(outdir, "refinement.csv"), "w") do io
        println(io, "refinement,norm,coupled,limit,kappa_corrector,fixed_reference_remainder")
        for name in ("space_2", "time_2", "time_4"), m in 1:2
            println(io, join([name, "H$(m-1)", (record["sensitivity/$name/$field"][m] for field in
                ("coupled", "limit", "corrector_weighted", "fixed_reference_remainder"))...], ','))
        end
    end
    println("Refined-reference deviation rates: ", record["orders/deviation"])
    println("Refined-reference remainder rates: ", record["orders/remainder"])
    println("Relative effect of transport/material terms: ", record["forcing_only/relative_difference"])
    println("Wrote ", joinpath(outdir, "heterogeneous.jld2"))
end

main()
