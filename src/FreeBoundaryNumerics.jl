"""
FreeBoundaryNumerics

Numerical validation of the low-conductivity free boundary problem: P1 finite
elements on a fitted annular mesh, Newton for the nonlinear Robin problem,
Fourier differentiation of the radial graph, Heun time stepping, the limiting
flow, its discrete corrector, and the experiments E1–E6 of the paper.
"""
module FreeBoundaryNumerics

using LinearAlgebra, SparseArrays, FFTW, Printf, Statistics
using JLD2

export Params, Material, M0, M1, M2, M2family, withkappa, kstar, kminus, T_bound, swept_constants,
    excess_constant, dtn_eigenvalue,
    theta_grid, fourier_derivative, sobolev_norm, linf, mode_amplitude, mode_coefficients,
    translation_fraction, resample, exponential_filter!, observed_orders,
    radial_curvature, metric_factor, graph_stretch, polygon_area, spectral_area, polygon_perimeter,
    spectral_perimeter, ellipse_radius, flower_radius, parallel_ellipse_radius, crossing_time,
    extrapolated_corner_time,
    AnnulusMesh, nnodes, node_coordinates, RobinSolver, DirichletSolver, solve_electrical!,
    discrete_energy, harmonic_extension, dtn_apply, corrector_data, dtn_nodal_differencing,
    mass_matrix, set_geometry!, boundary_mass,
    evolve_coupled, evolve_limit, evolve_limit_rk4, CoupledRun, LimitRun,
    run_E1, run_E2, run_E3, run_E4, run_E5, run_E6, run_E6_phase, save_results, load_results

include("materials.jl")
include("spectral.jl")
include("geometry.jl")
include("fem.jl")
include("flows.jl")
include("experiments.jl")

end # module
