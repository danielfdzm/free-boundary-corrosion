# Postprocess stored radii without changing the solver or archived Sobolev data.
using FFTW
using FreeBoundaryNumerics: resample, observed_orders

"""
    cm_norms(v; max_order=3, oversample=16)

Return the discrete C^m diagnostics for m = 0:max_order: the maximum, over
derivatives of order 0:m, of the absolute value of the real trigonometric
interpolant on an oversample*length(v) angular grid. These approximate the
continuous supremum norms; increasing `oversample` checks angular sampling.
For even input lengths the Nyquist coefficient represents a cosine and is split
between the positive and negative frequencies, including for odd derivatives.
No high-frequency filtering or differentiation of a piecewise-linear curve is used.
"""
function cm_norms(v::AbstractVector{<:Real}; max_order::Integer=3, oversample::Integer=16)
    n = length(v)
    n >= 2 || throw(ArgumentError("At least two angular samples are required"))
    max_order >= 0 || throw(ArgumentError("max_order must be nonnegative"))
    oversample >= 2 || throw(ArgumentError("oversample must be at least 2 to resolve off-grid Nyquist derivatives"))
    all(isfinite, v) || throw(ArgumentError("Angular samples must be finite"))
    count = oversample * n
    coefficients = rfft(Float64.(v)) ./ n
    # On the source grid the two Nyquist frequencies coincide. In the real
    # interpolant their equal half-coefficients give a*cos((n/2)*theta).
    iseven(n) && (coefficients[end] = real(coefficients[end]) / 2)
    padded = zeros(ComplexF64, count ÷ 2 + 1)
    result = zeros(Float64, max_order + 1)
    accumulated = 0.0
    for derivative in 0:max_order
        fill!(padded, 0)
        for (index, coefficient) in enumerate(coefficients)
            frequency = index - 1
            padded[index] = count * coefficient * (im * frequency)^derivative
        end
        accumulated = max(accumulated, maximum(abs, irfft(padded, count)))
        result[derivative + 1] = accumulated
    end
    return result
end

"C^0 through C^3 disk diagnostics, using both matched and exact references."
function disk_cm(record; oversample::Integer=16)
    norms(v) = cm_norms(v; oversample)
    kappas = record["disk/kappas"]
    radius = [record["disk/$kappa/R"] for kappa in kappas]
    R, W = record["disk/limit/R"], record["disk/limit/W"]
    Rexact, Wexact = record["disk/limit/R_exact"], record["disk/limit/W_exact"]
    return (
        dev=hcat([norms(r - R) for r in radius]...),
        rem=hcat([norms(r - R - kappa * W) for (r, kappa) in zip(radius, kappas)]...),
        dev_exact=hcat([norms(r - Rexact) for r in radius]...),
        rem_exact=hcat([norms(r - Rexact - kappa * Wexact) for (r, kappa) in zip(radius, kappas)]...),
    )
end

"Return a new coupled-refinement dictionary with recomputed C^m diagnostics."
function coupled_cm(record; oversample::Integer=16)
    result = copy(record)
    norms(v) = cm_norms(v; oversample)
    production = record["production/R"]
    n = length(production)
    R, W, kappa = record["reference_R"], record["reference_W"], record["kappa"]
    result["production/corrected_remainder"] = norms(production - R - kappa * W)
    # The matched references are not contained in this refinement record. Do not
    # retain an archived H^m value in a dictionary of C^m diagnostics.
    delete!(result, "production/matched_remainder")
    for name in record["run_names"]
        radius = record["$name/R"]
        common = resample(radius, n)
        result["$name/R_common"] = common
        result["$name/production_difference"] = norms(common - production)
        result["$name/corrected_remainder"] = norms(common - R - kappa * W)
        result["$name/discarded_modes"] = norms(radius - resample(common, length(radius)))
    end
    result["time_2/previous_step_difference"] = copy(result["time_2/production_difference"])
    result["time_4/previous_step_difference"] = norms(result["time_4/R_common"] - result["time_2/R_common"])
    result["time_orders"] = log2.(result["time_2/previous_step_difference"] ./ result["time_4/previous_step_difference"])
    result["diagnostic_norm"] = "max_{0<=j<=m} max_l |d_theta^j I_N v(2pi*l/(oversample*N))|"
    result["diagnostic_oversample"] = oversample
    return result
end

"Return a new heterogeneous dictionary with recomputed C^0/C^1 diagnostics."
function heterogeneous_cm(record; oversample::Integer=16, max_order::Integer=1)
    result = copy(record)
    norms(v) = cm_norms(v; max_order, oversample)
    n, kappas = record["nt"], record["kappas"]
    base, reference = record["limit"], record["reference"]
    R, W = resample(reference["R"], n), resample(reference["W"], n)
    runs = record["coupled"]
    for (name, radius, corrector) in (("", R, W), ("matched_", base["R"], base["W"]))
        result["errors/$(name)deviation"] = hcat([norms(run["R"] - radius) for run in runs]...)
        result["errors/$(name)remainder"] = hcat([norms(run["R"] - radius - kappa * corrector) for (run, kappa) in zip(runs, kappas)]...)
    end
    for name in ("deviation", "remainder")
        result["orders/$name"] = permutedims(hcat([observed_orders(kappas, result["errors/$name"][m, :]) for m in 1:max_order+1]...))
    end
    kappa = kappas[end]
    for name in ("space_2", "time_2", "time_4")
        lr, cr = record["refinement/limit/$name"], record["refinement/coupled/$name"]
        result["sensitivity/$name/coupled"] = norms(resample(cr["R"], n) - runs[end]["R"])
        result["sensitivity/$name/limit"] = norms(resample(lr["R"], n) - base["R"])
        result["sensitivity/$name/corrector_weighted"] = kappa * norms(resample(lr["W"], n) - base["W"])
        result["sensitivity/$name/fixed_reference_remainder"] = norms(resample(cr["R"], n) - R - kappa * W)
    end
    for kind in ("limit", "coupled")
        r2, r4 = record["refinement/$kind/time_2"], record["refinement/$kind/time_4"]
        result["sensitivity/time_order/$kind"] = log2.(result["sensitivity/time_2/$kind"] ./ norms(r2["R"] - r4["R"]))
    end
    r2, r4 = record["refinement/limit/time_2"], record["refinement/limit/time_4"]
    result["sensitivity/time_order/corrector"] = log2.(result["sensitivity/time_2/corrector_weighted"] ./ (kappa * norms(r2["W"] - r4["W"])))
    result["sensitivity/reference/limit"] = norms(R - base["R"])
    result["sensitivity/reference/corrector_weighted"] = kappa * norms(W - base["W"])
    result["forcing_only/difference"] = norms(base["W"] - record["forcing_only/W"])
    result["forcing_only/relative_difference"] = result["forcing_only/difference"] ./ norms(base["W"])
    result["diagnostic_norm"] = "max_{0<=j<=m} max_l |d_theta^j I_N v(2pi*l/(oversample*N))|"
    result["diagnostic_oversample"] = oversample
    return result
end
