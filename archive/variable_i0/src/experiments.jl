# ---------------------------------------------------------------------------
# Experiments E1–E5 of the paper. Each writes one JLD2 file of plain
# arrays and dictionaries; figures and tables are produced separately.
# ---------------------------------------------------------------------------

function save_results(path::AbstractString, data::Dict{String,Any})
    mkpath(dirname(path))
    jldopen(path, "w") do f
        for (k, v) in data
            f[k] = v
        end
    end
    return path
end

load_results(path::AbstractString) = JLD2.load(path)

_kappas(ks) = [2.0^(-k) for k in ks]

function _norms(dev)
    return [sobolev_norm(dev, m) for m in 0:3]
end

# ===========================================================================
# E1 — verification of the scheme
# ===========================================================================
"""
E1: (a) M0 circle, (b) M0 ellipse against the exact parallel curve,
(c) space–time error map under M2, (d) DtN validation of the flux recovery.
"""
function run_E1(outdir::AbstractString; quick::Bool=false)
    out = Dict{String,Any}()
    P = Params(kappa=0.5)
    println("E1(a): homogeneous circle")
    run = evolve_coupled(P, M0, ones(48); nr=12, dt=0.02, T=1.0)
    out["circle/times"] = run.saved_times
    out["circle/error"] = [maximum(abs, run.Rhist[:, k] .- (1 - P.beta * run.saved_times[k])) for k in eachindex(run.saved_times)]
    out["circle/final_error"] = maximum(abs, run.R .- (1 - P.beta))
    out["circle/summary"] = run.summary
    out["circle/history"] = run.history

    println("E1(b): homogeneous ellipse against the parallel curve")
    A, B = 1.0, 0.6
    K0 = A / B^2
    Tell = 0.8 / (P.beta * K0)          # = 2.4
    dte = 0.0025
    nts = quick ? [64, 128, 256] : [128, 256, 512]
    out["ellipse/A"] = A; out["ellipse/B"] = B; out["ellipse/T"] = Tell; out["ellipse/dt"] = dte
    out["ellipse/nts"] = nts
    for nt in nts
        theta = theta_grid(nt)
        run = evolve_coupled(P, M0, ellipse_radius.(theta, A, B); nr=32, dt=dte, T=Tell, store_every=8)
        errs = similar(run.saved_times); l2 = similar(errs)
        for (k, t) in enumerate(run.saved_times)
            Rex = parallel_ellipse_radius(theta, t, A, B, P.beta)
            errs[k] = maximum(abs, run.Rhist[:, k] .- Rex)
            l2[k] = sobolev_norm(run.Rhist[:, k] .- Rex, 0)
        end
        out["ellipse/$nt/times"] = run.saved_times
        out["ellipse/$nt/error_inf"] = errs
        out["ellipse/$nt/error_L2"] = l2
        out["ellipse/$nt/max_curvature"] = run.history["max_curvature"]
        out["ellipse/$nt/history_times"] = run.history["time"]
        out["ellipse/$nt/summary"] = run.summary
        @printf("   nt=%d  final max error %.3e\n", nt, errs[end])
    end
    # time-step dependence at the finest mesh
    for dtt in (0.005, 0.00125)
        nt = nts[end]; theta = theta_grid(nt)
        run = evolve_coupled(P, M0, ellipse_radius.(theta, A, B); nr=32, dt=dtt, T=Tell, store_every=max(1, round(Int, 0.02 / dtt)))
        Rex = parallel_ellipse_radius(theta, Tell, A, B, P.beta)
        out["ellipse/dt_$dtt/final_error_inf"] = maximum(abs, run.R .- Rex)
    end
    out["ellipse/exact_curvature"] = t -> K0 / (1 - P.beta * K0 * t)

    println("E1(c): space-time error map")
    nts_map = quick ? [32, 64, 128] : [32, 64, 128, 256, 512]
    dts_map = quick ? [0.05, 0.025, 0.0125] : [0.05, 0.025, 0.0125, 0.00625, 0.003125]
    Tm = 0.5
    nt_ref = 2 * nts_map[end]; dt_ref = dts_map[end]
    ref = evolve_coupled(P, M2, ones(nt_ref); nr=nt_ref ÷ 4, dt=dt_ref, T=Tm, store_every=10^9)
    @printf("   reference %d x %d, dt=%g done (%.1fs)\n", nt_ref, nt_ref ÷ 4, dt_ref, ref.summary["elapsed_seconds"])
    E2 = zeros(length(nts_map), length(dts_map)); Einf = similar(E2)
    T2 = zeros(length(nts_map), length(dts_map)); Tinf = similar(T2)     # temporal error against the same-mesh finest step
    for (i, nt) in enumerate(nts_map)
        tref = evolve_coupled(P, M2, ones(nt); nr=nt ÷ 4, dt=dts_map[end] / 2, T=Tm, store_every=10^9)
        for (j, dt) in enumerate(dts_map)
            run = evolve_coupled(P, M2, ones(nt); nr=nt ÷ 4, dt=dt, T=Tm, store_every=10^9)
            d = resample(run.R, nt_ref) .- ref.R
            E2[i, j] = sobolev_norm(d, 0); Einf[i, j] = maximum(abs, d)
            dtm = run.R .- tref.R
            T2[i, j] = sobolev_norm(dtm, 0); Tinf[i, j] = maximum(abs, dtm)
            @printf("   nt=%4d dt=%.6f  L2 err %.3e  temporal %.3e  (newton max %d)\n", nt, dt, E2[i, j], T2[i, j], run.summary["max_newton_iterations"])
        end
    end
    out["map/nts"] = nts_map; out["map/dts"] = dts_map; out["map/L2"] = E2; out["map/Linf"] = Einf
    out["map/temporal_L2"] = T2; out["map/temporal_Linf"] = Tinf
    out["map/nt_ref"] = nt_ref; out["map/dt_ref"] = dt_ref; out["map/T"] = Tm; out["map/kappa"] = P.kappa
    out["map/space_orders"] = observed_orders(1 ./ nts_map, E2[:, end])
    out["map/time_orders"] = [observed_orders(dts_map, T2[i, :]) for i in eachindex(nts_map)]

    println("E1(d): DtN validation")
    nmax = 64
    meshes = quick ? [(128, 32, 1.0), (256, 64, 1.0)] : [(256, 64, 1.0), (512, 128, 1.0), (1024, 128, 1.0), (1024, 128, 2.0)]
    out["dtn/meshes"] = [collect(m) for m in meshes]
    out["dtn/n"] = collect(1:nmax)
    for (nt, nr, grading) in meshes
        mesh = AnnulusMesh(nt, nr; B=P.B, grading=grading)
        D = DirichletSolver(mesh)
        set_geometry!(D, ones(nt))
        theta = mesh.theta
        rel_var = zeros(nmax); rel_naive = zeros(nmax)
        for n in 1:nmax
            g = cos.(n .* theta)
            Λg, u, F = dtn_apply(D, g)
            lam = dtn_eigenvalue(n, 1.0, P.B)
            rel_var[n] = abs(mode_coefficients(Λg, n)[1] - lam) / lam
            naive = dtn_nodal_differencing(D, u)
            rel_naive[n] = abs(mode_coefficients(naive, n)[1] - lam) / lam
        end
        key = "dtn/$(nt)x$(nr)_g$(grading)"
        out["$key/variational"] = rel_var
        out["$key/naive"] = rel_naive
        @printf("   %d x %d (grading %g): rel. error n=1: %.2e, n=8: %.2e, n=64: %.2e (naive n=8: %.2e)\n",
            nt, nr, grading, rel_var[1], rel_var[8], rel_var[64], rel_naive[8])
    end
    delete!(out, "ellipse/exact_curvature")
    out["ellipse/K0"] = K0
    save_results(joinpath(outdir, "E1.jld2"), out)
    return out
end

# ===========================================================================
# E2 — balances, quantitative excess bound, lifetime envelope
# ===========================================================================
function run_E2(outdir::AbstractString; quick::Bool=false)
    out = Dict{String,Any}()
    nt, nr, dt, T = quick ? (128, 32, 0.01, 1.0) : (256, 64, 0.005, 1.0)
    kappas = [1.0, 0.25, 0.0625, 0.015625]
    out["balance/kappas"] = kappas
    for (A2, A1, tag) in ((1.0, 1.0, "sym"), (1.5, 0.5, "asym"))
        ratios = Float64[]
        for κ in kappas
            P = Params(kappa=κ, A2=A2, A1=A1)
            println("E2 balances ($tag): kappa=$κ")
            run = evolve_coupled(P, M2, ones(nt); nr=nr, dt=dt, T=T, store_every=10^9)
            h = run.history
            out["balance/$tag/$κ/time"] = h["time"]
            out["balance/$tag/$κ/excess"] = h["excess"]
            out["balance/$tag/$κ/excess_bound"] = h["excess_bound"]
            out["balance/$tag/$κ/charge"] = h["charge"]
            out["balance/$tag/$κ/relative_area_balance"] = h["relative_area_balance"]
            out["balance/$tag/$κ/summary"] = run.summary
            push!(ratios, h["excess"][end] / h["excess_bound"][end])
            @printf("   excess/bound at T: %.4f   min excess %.3e   max charge %.2e   area balance %.2e\n",
                ratios[end], run.summary["min_excess"], run.summary["max_charge_residual"], run.summary["max_relative_area_balance"])
        end
        out["balance/$tag/ratio_at_T"] = ratios
        out["balance/$tag/predicted_ratio"] = (A2 + A1) / min(A2, A1)
        out["balance/$tag/A2A1"] = [A2, A1]
    end
    out["balance/nt"] = nt; out["balance/nr"] = nr; out["balance/dt"] = dt; out["balance/T"] = T
    # homogeneous polygon area-defect prediction (1-(R(T)/R(0))^2)(1-1/cos(pi/n))
    P0 = Params(kappa=0.5)
    run0 = evolve_coupled(P0, M0, ones(48); nr=12, dt=0.02, T=1.0, store_every=10^9)
    out["balance/homogeneous_area_defect"] = run0.history["area_balance_residual"][end] / polygon_area(ones(48))
    out["balance/homogeneous_area_defect_predicted"] = (1 - (1 - P0.beta)^2) * (1 - 1 / cos(π / 48))

    println("E2 lifetime envelope")
    Pl = Params(kappa=0.5)
    shapes = quick ? [("disk", ones(128), 32, 1.5), ("flower", flower_radius.(theta_grid(128), 0.15, 5), 32, 1.0)] :
        [("disk", ones(256), 64, 1.5), ("ellipse", ellipse_radius.(theta_grid(512), 1.0, 0.6), 48, 1.5),
         ("flower", flower_radius.(theta_grid(512), 0.15, 5), 64, 1.0)]
    for (name, R0, nrs, Ts) in shapes
        snaps = name == "flower" ? [0.0, 0.5, 1.0] : Float64[]
        run = evolve_coupled(Pl, M2, R0; nr=nrs, dt=0.005, T=Ts, store_every=10^9, snapshot_times=snaps)
        for (is, sn) in enumerate(run.snapshots)
            out["lifetime/$name/snapshot/$is/t"] = sn.t
            out["lifetime/$name/snapshot/$is/R"] = sn.R
            out["lifetime/$name/snapshot/$is/phi"] = sn.phi
        end
        out["lifetime/$name/nsnapshots"] = length(run.snapshots)
        out["lifetime/$name/nt"] = length(R0); out["lifetime/$name/nr"] = nrs; out["lifetime/$name/kappa"] = Pl.kappa
        out["lifetime/$name/time"] = run.history["time"]
        out["lifetime/$name/sqrt_area"] = sqrt.(run.history["area_spectral"])
        out["lifetime/$name/area0"] = spectral_area(R0)
        out["lifetime/$name/max_curvature"] = run.history["max_curvature"]
        out["lifetime/$name/summary"] = run.summary
        barrier = sqrt(spectral_area(R0)) .- Pl.beta * M2.imin * sqrt(π) .* run.history["time"]
        @printf("   %s: min(barrier - sqrt area) = %.3e\n", name, minimum(barrier .- sqrt.(run.history["area_spectral"])))
    end
    run = evolve_coupled(Pl, M0, ones(256); nr=64, dt=0.005, T=1.5, store_every=10^9)
    out["lifetime/M0disk/time"] = run.history["time"]
    out["lifetime/M0disk/sqrt_area"] = sqrt.(run.history["area_spectral"])
    out["lifetime/M0disk/area0"] = spectral_area(ones(256))
    out["lifetime/M0disk/barrier_gap"] = maximum(abs, sqrt.(run.history["area_spectral"]) .- (sqrt(π) .- Pl.beta * sqrt(π) .* run.history["time"]))
    out["lifetime/imin"] = M2.imin; out["lifetime/beta"] = Pl.beta
    out["lifetime/T_area_disk"] = sqrt(π) / (Pl.beta * M2.imin * sqrt(π))
    save_results(joinpath(outdir, "E2.jld2"), out)
    return out
end

# ===========================================================================
# E3 — frozen expansion and the transfer factor (flagship)
# ===========================================================================
function run_E3(outdir::AbstractString; quick::Bool=false, grading::Real=1.0)
    out = Dict{String,Any}()
    nt, nr = quick ? (256, 32) : (1024, 128)
    ks = quick ? collect(0:6) : collect(0:10)
    ns = quick ? [1, 2, 4, 8, 16] : [1, 2, 4, 8, 16, 32, 64]
    deltas = [0.05, 0.2]
    kappas = _kappas(ks)
    P = Params(kappa=1.0)
    a = P.A2 + P.A1
    mesh = AnnulusMesh(nt, nr; B=P.B, grading=grading)
    D = DirichletSolver(mesh)
    R = ones(nt)
    theta = mesh.theta
    sz = (length(kappas), length(ns), length(deltas))
    tau0 = zeros(sz); tau2 = zeros(sz); mu = zeros(sz); muh = zeros(sz)
    einf = zeros(sz); rem_inf = zeros(sz); rem0 = zeros(sz); rem2 = zeros(sz)
    Jk = zeros(sz); J0 = zeros(sz); J1 = zeros(sz); dexc = zeros(sz); charge = zeros(sz); newt = zeros(Int, sz)
    lam_h = zeros(length(ns))
    for (id, δ) in enumerate(deltas), (jn, n) in enumerate(ns)
        mat = M1(n, δ)
        cd = corrector_data(D, P, mat, R)
        e1 = cd.gamma1                       # -Lambda_h phi_eq / a  (i0 = 1)
        lam_h[jn] = -a * mode_coefficients(e1, n)[1] / δ
        lam = dtn_eigenvalue(n, 1.0, P.B)
        for (ik, κ) in enumerate(kappas)
            S = RobinSolver(mesh, withkappa(P, κ))
            res = solve_electrical!(S, mat, R)
            g = S.phi[1:nt]
            e = g .- cd.eqb
            tau0[ik, jn, id] = sobolev_norm(e, 0) / sobolev_norm(cd.eqb, 0)
            tau2[ik, jn, id] = sobolev_norm(e, 2) / sobolev_norm(cd.eqb, 2)
            mu[ik, jn, id] = min(1.0, κ * lam / a)
            muh[ik, jn, id] = min(1.0, κ * lam_h[jn] / a)
            einf[ik, jn, id] = linf(e)
            rem_inf[ik, jn, id] = linf(e .- κ .* e1)
            rem0[ik, jn, id] = sobolev_norm(e .- κ .* e1, 0)
            rem2[ik, jn, id] = sobolev_norm(e .- κ .* e1, 2)
            Jk[ik, jn, id] = discrete_energy(S)
            J0[ik, jn, id] = cd.J0; J1[ik, jn, id] = cd.J1
            dexc[ik, jn, id] = res.excess; charge[ik, jn, id] = abs(res.charge); newt[ik, jn, id] = res.iterations
        end
        @printf("E3 delta=%.2f n=%2d: tau(kappa=1)=%.4f tau(kappa=2^-%d)=%.3e  ratio range [%.3f, %.3f]  lam_h/lam=%.5f\n",
            δ, n, tau0[1, jn, id], ks[end], tau0[end, jn, id],
            minimum(tau0[:, jn, id] ./ mu[:, jn, id]), maximum(tau0[:, jn, id] ./ mu[:, jn, id]), lam_h[jn] / lam)
    end
    out["kappas"] = kappas; out["ks"] = ks; out["ns"] = ns; out["deltas"] = deltas
    out["nt"] = nt; out["nr"] = nr; out["grading"] = grading; out["a"] = a
    out["tau0"] = tau0; out["tau2"] = tau2; out["mu"] = mu; out["mu_h"] = muh
    out["e_inf"] = einf; out["remainder_inf"] = rem_inf; out["remainder_H0"] = rem0; out["remainder_H2"] = rem2
    out["J_kappa"] = Jk; out["J0_h"] = J0; out["J1_h"] = J1; out["excess"] = dexc; out["charge"] = charge
    out["newton_iterations"] = newt
    out["lambda_h"] = lam_h; out["lambda"] = [dtn_eigenvalue(n, 1.0, P.B) for n in ns]
    save_results(joinpath(outdir, "E3.jld2"), out)
    return out
end

# ===========================================================================
# E4 — linearization map and absence of parabolic smoothing (optional)
# ===========================================================================
function run_E4(outdir::AbstractString; quick::Bool=false)
    out = Dict{String,Any}()
    nt, nr = quick ? (128, 32) : (512, 128)
    P = Params(kappa=0.5)
    matr = Material("M1r", (x, y) -> 1.0, (x, y) -> (0.0, 0.0), (x, y) -> 0.35 * (x^2 + y^2), 1.0, 0.0, 0.0)
    ms = quick ? collect(1:8) : collect(1:24)
    ks = quick ? collect(0:3) : collect(0:6)
    kappas = _kappas(ks)
    mesh = AnnulusMesh(nt, nr; B=P.B)
    theta = mesh.theta
    ratio = zeros(length(ms), length(kappas))
    eps = 1e-4
    for (ik, κ) in enumerate(kappas)
        S = RobinSolver(mesh, withkappa(P, κ))
        for (im, m) in enumerate(ms)
            wave = cos.(m .* theta)
            lam = dtn_eigenvalue(m, 1.0, P.B)
            sigma = P.beta * P.A2 * 0.7 * κ * lam / (P.A2 + P.A1 + κ * lam)
            plus = solve_electrical!(S, matr, 1 .+ eps .* wave).speed
            minus = solve_electrical!(S, matr, 1 .- eps .* wave).speed
            est = mode_coefficients((plus .- minus) ./ (2eps), m)[1]
            ratio[im, ik] = est / sigma
        end
        @printf("E4 kappa=%.4g: sigma ratio range [%.4f, %.4f]\n", κ, minimum(ratio[:, ik]), maximum(ratio[:, ik]))
    end
    out["sigma/ms"] = ms; out["sigma/kappas"] = kappas; out["sigma/ratio"] = ratio; out["sigma/nt"] = nt; out["sigma/nr"] = nr
    # spectrum transport: perturbed disks under the radial data, kappa = 1/2
    println("E4 spectrum transport")
    modes = quick ? [2, 4] : [2, 4, 8, 16]
    dt, T = quick ? (0.01, 0.5) : (0.005, 1.0)
    base = evolve_coupled(P, matr, ones(nt); nr=nr, dt=dt, T=T, store_every=4)
    out["spectrum/times"] = base.saved_times
    out["spectrum/modes"] = modes
    amp0 = 1e-3
    for m in modes
        run = evolve_coupled(P, matr, 1 .+ amp0 .* cos.(m .* theta); nr=nr, dt=dt, T=T, store_every=4)
        nmodes = nt ÷ 4
        spec = zeros(nmodes, length(run.saved_times))
        for k in eachindex(run.saved_times)
            d = run.Rhist[:, k] .- base.Rhist[:, k]
            for n in 1:nmodes
                spec[n, k] = mode_amplitude(d, n)
            end
        end
        out["spectrum/$m/amplitudes"] = spec
        out["spectrum/$m/mode_amplitude"] = spec[m, :]
        @printf("   m=%d: amplitude %.3e -> %.3e\n", m, spec[m, 1], spec[m, end])
    end
    out["spectrum/amp0"] = amp0; out["spectrum/kappa"] = P.kappa
    save_results(joinpath(outdir, "E4.jld2"), out)
    return out
end

# ===========================================================================
# E5 — the low-conductivity limit and the corrector (flagship)
# ===========================================================================
function run_E5(outdir::AbstractString; quick::Bool=false)
    out = Dict{String,Any}()
    nt, nr, dt, T = quick ? (128, 32, 0.01, 0.5) : (512, 128, 0.0025, 0.5)
    ks = quick ? collect(0:4) : collect(0:8)
    kappas = _kappas(ks)
    P = Params(kappa=1.0)
    theta = theta_grid(nt)
    shapes = quick ? [("disk", ones(nt), ks)] : [("disk", ones(nt), ks), ("flower", flower_radius.(theta, 0.15, 5), [0, 2, 4, 6, 8])]
    out["nt"] = nt; out["nr"] = nr; out["dt"] = dt; out["T"] = T; out["theta"] = theta
    for (name, R0, kss) in shapes
        κs = _kappas(kss)
        println("E5 [$name]: matched limit flow and discrete corrector")
        lim = evolve_limit(P, M2, R0; dt=dt, T=T, nr=nr, corrector=true, store_every=1)
        R0h = lim.R; R1h = lim.W
        out["$name/limit/R"] = R0h; out["$name/limit/W"] = R1h
        out["$name/limit/Rhist"] = lim.Rhist[:, 1:4:end]; out["$name/limit/Whist"] = lim.Whist[:, 1:4:end]
        out["$name/limit/times"] = lim.saved_times[1:4:end]
        println("E5 [$name]: converged references")
        fine = evolve_limit_rk4(P, M2, resample(R0, 4nt); dt=dt / 4, T=T, store_every=10^9)
        R0fine = resample(fine.R, nt)
        limf = evolve_limit(P, M2, resample(R0, 2nt); dt=dt, T=T, nr=2nr, corrector=true, store_every=10^9)
        R1fine = resample(limf.W, nt)
        out["$name/limit/R_fine"] = R0fine; out["$name/limit/W_fine"] = R1fine
        out["$name/limit/R_matched_minus_fine"] = _norms(R0h .- R0fine)
        out["$name/limit/W_matched_minus_fine"] = _norms(R1h .- R1fine)
        out["$name/limit/W_translation_fraction"] = translation_fraction(R1h)
        # geometry data on the limit mesh for the bulk comparison
        meshL = AnnulusMesh(nt, nr; B=P.B)
        DL = DirichletSolver(meshL)
        cdL = corrector_data(DL, P, M2, R0h)
        u0h = cdL.u0
        ML = mass_matrix(meshL, DL.x, DL.y)
        out["$name/limit/u0"] = u0h; out["$name/limit/x"] = copy(DL.x); out["$name/limit/y"] = copy(DL.y)
        out["$name/kappas"] = κs
        nk = length(κs)
        dev_norms = zeros(4, nk); rem_norms = zeros(4, nk); dev_fine = zeros(4, nk); rem_fine = zeros(4, nk)
        rem_wrong = zeros(nk)
        trace_err = zeros(nk); trace_rem = zeros(nk); speed_err = zeros(nk); speed_rem = zeros(nk)
        excess = zeros(nk); intD1 = zeros(nk); bulk_L2 = zeros(nk); bulk_Linf = zeros(nk)
        transl = zeros(nk); newton = zeros(Int, nk); charge = zeros(nk); elapsed = zeros(nk)
        DK = DirichletSolver(AnnulusMesh(nt, nr; B=P.B))
        rays = [1, nt ÷ 4 + 1, nt ÷ 2 + 1, 3nt ÷ 4 + 1]
        for (ik, κ) in enumerate(κs)
            Pk = withkappa(P, κ)
            @printf("E5 [%s]: kappa = 2^-%d\n", name, kss[ik])
            store = (ik == nk) ? 1 : 4
            run = evolve_coupled(Pk, M2, R0; nr=nr, dt=dt, T=T, store_every=store)
            dev = run.R .- R0h
            rem = dev .- κ .* R1h
            dev_norms[:, ik] = _norms(dev); rem_norms[:, ik] = _norms(rem)
            dev_fine[:, ik] = _norms(run.R .- R0fine); rem_fine[:, ik] = _norms(run.R .- R0fine .- κ .* R1fine)
            rem_wrong[ik] = sobolev_norm(dev .+ κ .* R1h, 0)
            transl[ik] = translation_fraction(dev)
            # trace and speed expansions on the kappa-geometry
            cd = corrector_data(DK, P, M2, run.R)
            g = run.phi[1:nt]
            e = g .- cd.eqb
            trace_err[ik] = linf(e); trace_rem[ik] = linf(e .- κ .* cd.gamma1)
            Dk = @. P.beta * cd.i0b * exp(P.A2 * e)
            D0 = @. P.beta * cd.i0b
            D1 = @. P.beta * P.A2 * cd.i0b * cd.gamma1
            speed_err[ik] = linf(Dk .- D0); speed_rem[ik] = linf(Dk .- D0 .- κ .* D1)
            Mi0 = boundary_mass(DK.len, cd.i0b)
            intD1[ik] = P.beta * P.A2 * sum(Mi0 * cd.gamma1)
            excess[ik] = run.history["excess"][end]
            newton[ik] = run.summary["max_newton_iterations"]; charge[ik] = run.summary["max_charge_residual"]
            elapsed[ik] = run.summary["elapsed_seconds"]
            # bulk field: node-wise pullback difference with the limit potential
            diff = run.phi .- u0h
            bulk_L2[ik] = sqrt(dot(diff, ML * diff)); bulk_Linf[ik] = linf(diff)
            prof = zeros(nr + 1, length(rays)); dist = zeros(nr + 1, length(rays))
            for (ir, j) in enumerate(rays), k in 0:nr
                prof[k+1, ir] = diff[k*nt+j] / κ
                dist[k+1, ir] = hypot(run.x[k*nt+j], run.y[k*nt+j]) - run.R[j]
            end
            out["$name/$κ/profiles"] = prof; out["$name/$κ/profile_distances"] = dist
            out["$name/$κ/R"] = run.R; out["$name/$κ/g"] = g; out["$name/$κ/phi"] = run.phi
            out["$name/$κ/x"] = run.x; out["$name/$κ/y"] = run.y
            out["$name/$κ/Rhist"] = run.Rhist; out["$name/$κ/times"] = run.saved_times
            out["$name/$κ/excess_history"] = run.history["excess"]
            out["$name/$κ/gamma1"] = cd.gamma1
            @printf("   |dev|_H0=%.3e  |rem|_H0=%.3e  |rem|_H1=%.3e  trace rem %.3e  speed rem %.3e  excess %.3e  intD1 %.1e  bulk L2 %.3e\n",
                dev_norms[1, ik], rem_norms[1, ik], rem_norms[2, ik], trace_rem[ik], speed_rem[ik], excess[ik], intD1[ik], bulk_L2[ik])
        end
        out["$name/rays"] = rays
        out["$name/dev_norms"] = dev_norms; out["$name/rem_norms"] = rem_norms
        out["$name/dev_norms_fine"] = dev_fine; out["$name/rem_norms_fine"] = rem_fine
        out["$name/rem_wrong_sign"] = rem_wrong
        out["$name/trace_err"] = trace_err; out["$name/trace_rem"] = trace_rem
        out["$name/speed_err"] = speed_err; out["$name/speed_rem"] = speed_rem
        out["$name/excess"] = excess; out["$name/intD1"] = intD1
        out["$name/bulk_L2"] = bulk_L2; out["$name/bulk_Linf"] = bulk_Linf
        out["$name/translation_fraction"] = transl
        out["$name/newton"] = newton; out["$name/charge"] = charge; out["$name/elapsed"] = elapsed
        out["$name/orders/dev"] = [observed_orders(κs, dev_norms[m, :]) for m in 1:4]
        out["$name/orders/rem"] = [observed_orders(κs, rem_norms[m, :]) for m in 1:4]
        out["$name/orders/dev_fine"] = [observed_orders(κs, dev_fine[m, :]) for m in 1:4]
        out["$name/orders/rem_fine"] = [observed_orders(κs, rem_fine[m, :]) for m in 1:4]
        out["$name/orders/trace_err"] = observed_orders(κs, trace_err); out["$name/orders/trace_rem"] = observed_orders(κs, trace_rem)
        out["$name/orders/speed_rem"] = observed_orders(κs, speed_rem); out["$name/orders/excess"] = observed_orders(κs, excess)
        out["$name/orders/bulk_L2"] = observed_orders(κs, bulk_L2)
        out["$name/orders/rem_wrong_sign"] = observed_orders(κs, rem_wrong)
    end
    save_results(joinpath(outdir, "E5.jld2"), out)
    return out
end
