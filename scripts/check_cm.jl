#!/usr/bin/env julia
# Analytic checks of the postprocessing norm, including off-grid Nyquist extrema.
using Test, FreeBoundaryNumerics
include(joinpath(@__DIR__, "cm_diagnostics.jl"))

@testset "Real trigonometric C^m diagnostics" begin
    for n in (31, 32)
        theta = 2pi .* (0:n-1) ./ n
        @test cm_norms(fill(-2.5, n)) ≈ fill(2.5, 4)
        @test cm_norms(zeros(n)) == zeros(4)
        @test cm_norms(0.75 .* cos.(4theta); oversample=32) ≈ [0.75 * 4.0^m for m in 0:3] rtol=2e-11
        # Independent analytic derivatives of a mixed sine/cosine polynomial.
        oversample = 16
        fine = 2pi .* (0:oversample*n-1) ./ (oversample*n)
        samples = 0.7 .+ 0.6cos.(3theta) .+ 0.2sin.(5theta)
        expected = [maximum(abs, (j == 0 ? 0.7 : 0.0) .+
            0.6 * 3.0^j .* cos.(3fine .+ j*pi/2) .+
            0.2 * 5.0^j .* sin.(5fine .+ j*pi/2)) for j in 0:3]
        expected = accumulate(max, expected)
        @test cm_norms(samples; oversample) ≈ expected rtol=2e-11
    end
    # For n=32 the first derivative of cos(16theta) vanishes at every source
    # node, but its true supremum is 16, attained between those nodes.
    @test cm_norms([(-1.0)^j for j in 0:31]) ≈ [1.0, 16.0, 256.0, 4096.0] rtol=2e-13
    @test_throws ArgumentError cm_norms(ones(32); oversample=1)
    @test_throws ArgumentError cm_norms([0.0, NaN])
    @test_throws ArgumentError cm_norms(ones(32); max_order=-1)
end

root = normpath(joinpath(@__DIR__, ".."))
datadir = isempty(ARGS) ? joinpath(root, "data") : abspath(ARGS[1])
@testset "Raw records remain unchanged" begin
    E5 = load_results(joinpath(datadir, "E5.jld2"))
    coupled = load_results(joinpath(datadir, "E5_coupled_refinement.jld2"))
    heterogeneous = load_results(joinpath(datadir, "heterogeneous", "heterogeneous.jld2"))
    for (raw, compute) in ((E5, disk_cm), (coupled, coupled_cm), (heterogeneous, heterogeneous_cm))
        before = deepcopy(raw)
        derived = compute(raw)
        @test isequal(raw, before)
        @test derived !== raw
    end
    disk = disk_cm(E5)
    refinement = coupled_cm(coupled)
    @test disk.rem_exact[:, end] ≈ refinement["production/corrected_remainder"] rtol=1e-12
    @test size(disk.rem_exact) == (4, length(E5["disk/kappas"]))
    @test size(heterogeneous_cm(heterogeneous)["errors/remainder"]) == (2, length(heterogeneous["kappas"]))
end
