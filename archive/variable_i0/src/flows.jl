# ---------------------------------------------------------------------------
# Time integration: coupled evolution, limiting flow, tangent-linear corrector
# ---------------------------------------------------------------------------

const HISTORY_KEYS = ("time", "charge", "dissolution", "i0_integral", "excess", "excess_bound",
    "newton_iterations", "newton_residual", "min_triangle_area", "phi_min", "phi_max",
    "R_min", "R_max", "area_polygon", "area_spectral", "perimeter", "max_curvature", "min_curvature",
    "area_balance_residual", "relative_area_balance", "graph_stretch")

struct CoupledRun
    theta::Vector{Float64}
    R::Vector{Float64}
    phi::Vector{Float64}
    x::Vector{Float64}
    y::Vector{Float64}
    Rhist::Matrix{Float64}
    saved_times::Vector{Float64}
    history::Dict{String,Vector{Float64}}
    snapshots::Vector{Any}
    summary::Dict{String,Any}
    mesh::AnnulusMesh
    P::Params
end

function _push_history!(hist, t, R, res, phi, area0, integrated_flux)
    K = radial_curvature(R)
    areaP = polygon_area(R)
    push!(hist["time"], t)
    push!(hist["charge"], abs(res.charge))
    push!(hist["dissolution"], res.dissolution)
    push!(hist["i0_integral"], res.i0_integral)
    push!(hist["excess"], res.excess)
    push!(hist["excess_bound"], res.excess_bound)
    push!(hist["newton_iterations"], res.iterations)
    push!(hist["newton_residual"], res.residual)
    push!(hist["min_triangle_area"], res.min_triangle_area)
    push!(hist["phi_min"], minimum(phi)); push!(hist["phi_max"], maximum(phi))
    push!(hist["R_min"], minimum(R)); push!(hist["R_max"], maximum(R))
    push!(hist["area_polygon"], areaP)
    push!(hist["area_spectral"], spectral_area(R))
    push!(hist["perimeter"], spectral_perimeter(R))
    push!(hist["max_curvature"], maximum(K)); push!(hist["min_curvature"], minimum(K))
    push!(hist["area_balance_residual"], area0 - areaP - integrated_flux)
    push!(hist["relative_area_balance"], abs(area0 - areaP - integrated_flux) / area0)
    push!(hist["graph_stretch"], graph_stretch(R))
    return nothing
end

"""
Coupled evolution of the radial graph by Heun's method: the electrical problem
is solved at the beginning and at the predicted end of every step, the radii
are advanced with the averaged speeds, and the fitted mesh is rebuilt at fixed
connectivity. `stop(t, R, history)` may end the run early.
"""
function evolve_coupled(P::Params, mat::Material, R0::AbstractVector; nr::Int, dt::Real, T::Real,
        grading::Real=1.0, store_every::Int=1, snapshot_times=Float64[], stop=nothing, verbose=false)
    nt = length(R0)
    mesh = AnnulusMesh(nt, nr; B=P.B, grading=grading)
    S = RobinSolver(mesh, P)
    R = collect(float.(R0))
    nsteps = round(Int, T / dt)
    abs(nsteps * dt - T) < 1e-10 || error("T must be an integer multiple of dt")
    hist = Dict{String,Vector{Float64}}(k => Float64[] for k in HISTORY_KEYS)
    area0 = polygon_area(R)
    integrated_flux = 0.0
    guess = nothing
    Rhist = Vector{Vector{Float64}}(); saved_times = Float64[]
    snapshots = Any[]
    snap_steps = Set(round(Int, s / dt) for s in snapshot_times)
    started = time()
    res = nothing
    for k in 0:nsteps
        t = k * dt
        res = solve_electrical!(S, mat, R; phi0=guess)
        _push_history!(hist, t, R, res, S.phi, area0, integrated_flux)
        if k % store_every == 0 || k == nsteps
            push!(Rhist, copy(R)); push!(saved_times, t)
        end
        if k in snap_steps
            push!(snapshots, (t=t, R=copy(R), phi=copy(S.phi)))
        end
        verbose && k % 50 == 0 && @printf("  t=%.4f  maxK=%.4f  newton=%d  charge=%.2e\n", t, hist["max_curvature"][end], res.iterations, res.charge)
        if k == nsteps || (stop !== nothing && stop(t, R, hist))
            break
        end
        speed1 = res.speed
        phi1 = copy(S.phi)
        Rp = R .+ dt .* speed1
        res2 = solve_electrical!(S, mat, Rp; phi0=phi1)
        R .+= 0.5 * dt .* (speed1 .+ res2.speed)
        integrated_flux += 0.5 * P.beta * dt * (res.dissolution + res2.dissolution)
        guess = copy(S.phi)
        # keep the final solve consistent with the stored history (the solve at R happens at loop start)
    end
    # make sure phi corresponds to the final R
    if hist["time"][end] != nsteps * dt || S.phi === nothing
        # early stop: recompute nothing; phi already at the last recorded R
    end
    summary = Dict{String,Any}("nt" => nt, "nr" => nr, "dt" => dt, "T" => T, "T_final" => hist["time"][end],
        "max_charge_residual" => maximum(hist["charge"]), "max_newton_iterations" => maximum(hist["newton_iterations"]),
        "max_newton_residual" => maximum(hist["newton_residual"]),
        "max_relative_area_balance" => maximum(hist["relative_area_balance"]),
        "min_excess" => minimum(hist["excess"]), "elapsed_seconds" => time() - started,
        "kappa" => P.kappa, "material" => mat.name)
    x, y = copy(S.x), copy(S.y)
    return CoupledRun(mesh.theta, R, copy(S.phi), x, y, reduce(hcat, Rhist), saved_times, hist, snapshots, summary, mesh, P)
end

# ---------------------------------------------------------------------------
# Limiting flow dR/dt = -beta i0(R e_theta) sqrt(1 + (R_theta/R)^2)
# ---------------------------------------------------------------------------

function limit_speed(P::Params, mat::Material, R::AbstractVector, cs, sn)
    N = metric_factor(R)
    return [-P.beta * mat.i0(R[j] * cs[j], R[j] * sn[j]) * N[j] for j in eachindex(R)]
end

struct LimitRun
    theta::Vector{Float64}
    R::Vector{Float64}
    W::Vector{Float64}
    Rhist::Matrix{Float64}
    Whist::Matrix{Float64}
    saved_times::Vector{Float64}
    history::Dict{String,Vector{Float64}}
    summary::Dict{String,Any}
end

const LIMIT_KEYS = ("time", "R_min", "R_max", "area_spectral", "perimeter", "max_curvature", "min_curvature", "graph_stretch")

function _push_limit_history!(hist, t, R)
    K = radial_curvature(R)
    push!(hist["time"], t)
    push!(hist["R_min"], minimum(R)); push!(hist["R_max"], maximum(R))
    push!(hist["area_spectral"], spectral_area(R)); push!(hist["perimeter"], spectral_perimeter(R))
    push!(hist["max_curvature"], maximum(K)); push!(hist["min_curvature"], minimum(K))
    push!(hist["graph_stretch"], graph_stretch(R))
end

"""
Limiting flow on the same theta-grid and time grid as the coupled scheme
(Heun), optionally together with the discrete corrector W = dR_kappa/dkappa at
kappa = 0. The corrector is the exact tangent-linear model of the coupled
scheme: dW/dt = -beta G'(R)[W] + s(R), G(R) = i0(R e_theta) N(R),
s(R) = -beta N(R) i0 A2 gamma1_h(R), gamma1_h = -M_a^{-1} (K u0_h)|_Gamma,
with the same Fourier derivative, the same Heun stages and the same P1 mesh.
"""
function evolve_limit(P::Params, mat::Material, R0::AbstractVector; dt::Real, T::Real, nr::Int=0,
        corrector::Bool=false, grading::Real=1.0, store_every::Int=1, stop=nothing)
    nt = length(R0)
    theta = theta_grid(nt); cs = cos.(theta); sn = sin.(theta)
    D = corrector ? DirichletSolver(AnnulusMesh(nt, nr; B=P.B, grading=grading)) : nothing
    R = collect(float.(R0)); W = zeros(nt)
    nsteps = round(Int, T / dt)
    abs(nsteps * dt - T) < 1e-10 || error("T must be an integer multiple of dt")
    hist = Dict{String,Vector{Float64}}(k => Float64[] for k in LIMIT_KEYS)
    Rhist = Vector{Vector{Float64}}(); Whist = Vector{Vector{Float64}}(); saved = Float64[]
    function stage(R, W)
        N = metric_factor(R)
        R1 = fourier_derivative(R, 1)
        i0v = [mat.i0(R[j] * cs[j], R[j] * sn[j]) for j in 1:nt]
        k = @. -P.beta * i0v * N
        corrector || return k, zeros(nt)
        cd = corrector_data(D, P, mat, R)
        s = @. -P.beta * N * i0v * P.A2 * cd.gamma1
        gr = [begin gx, gy = mat.gradi0(R[j] * cs[j], R[j] * sn[j]); gx * cs[j] + gy * sn[j] end for j in 1:nt]
        W1 = fourier_derivative(W, 1)
        Gp = @. gr * N * W + (i0v / N) * (R1 * W1 / R^2 - R1^2 * W / R^3)
        l = @. -P.beta * Gp + s
        return k, l
    end
    started = time()
    for k in 0:nsteps
        t = k * dt
        _push_limit_history!(hist, t, R)
        if k % store_every == 0 || k == nsteps
            push!(Rhist, copy(R)); push!(Whist, copy(W)); push!(saved, t)
        end
        if k == nsteps || (stop !== nothing && stop(t, R, hist))
            break
        end
        k1, l1 = stage(R, W)
        Rp = R .+ dt .* k1; Wp = W .+ dt .* l1
        k2, l2 = stage(Rp, Wp)
        R .+= 0.5 * dt .* (k1 .+ k2)
        W .+= 0.5 * dt .* (l1 .+ l2)
    end
    summary = Dict{String,Any}("nt" => nt, "nr" => nr, "dt" => dt, "T" => T, "T_final" => hist["time"][end],
        "corrector" => corrector, "elapsed_seconds" => time() - started, "material" => mat.name)
    return LimitRun(theta, R, W, reduce(hcat, Rhist), reduce(hcat, Whist), saved, hist, summary)
end

"""
Limiting flow by the classical RK4 method with Fourier differentiation, used
for converged references and for the singular-time runs (optional exponential
filter after every step).
"""
function evolve_limit_rk4(P::Params, mat::Material, R0::AbstractVector; dt::Real, T::Real,
        filter::Bool=false, store_every::Int=1, stop=nothing)
    nt = length(R0)
    theta = theta_grid(nt); cs = cos.(theta); sn = sin.(theta)
    R = collect(float.(R0))
    nsteps = round(Int, T / dt)
    abs(nsteps * dt - T) < 1e-10 || error("T must be an integer multiple of dt")
    hist = Dict{String,Vector{Float64}}(k => Float64[] for k in LIMIT_KEYS)
    Rhist = Vector{Vector{Float64}}(); saved = Float64[]
    f(R) = limit_speed(P, mat, R, cs, sn)
    started = time()
    for k in 0:nsteps
        t = k * dt
        _push_limit_history!(hist, t, R)
        if k % store_every == 0 || k == nsteps
            push!(Rhist, copy(R)); push!(saved, t)
        end
        if k == nsteps || (stop !== nothing && stop(t, R, hist))
            break
        end
        k1 = f(R)
        k2 = f(R .+ 0.5dt .* k1)
        k3 = f(R .+ 0.5dt .* k2)
        k4 = f(R .+ dt .* k3)
        R .+= dt / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        filter && exponential_filter!(R)
    end
    summary = Dict{String,Any}("nt" => nt, "dt" => dt, "T" => T, "T_final" => hist["time"][end],
        "filter" => filter, "elapsed_seconds" => time() - started, "material" => mat.name)
    return LimitRun(theta, R, zeros(nt), reduce(hcat, Rhist), zeros(nt, 0), saved, hist, summary)
end
