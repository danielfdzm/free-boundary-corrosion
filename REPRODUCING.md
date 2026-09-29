# Reproducing the paper's numerics

## Environment and current records

Use the Julia environment pinned in `Project.toml` and `Manifest.toml`.
The recorded environment uses Julia 1.13.0. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/check.jl
```

The check validates both artifact checksum manifests, recorded numerical
source hashes, reference formulas, and stored refinement norms.

The additional noncircular, spatially heterogeneous experiment is documented
separately in [HETEROGENEOUS.md](HETEROGENEOUS.md). Its stored data are under
`data/heterogeneous/`; its two plots show conductivity convergence and a
space-time view of the computed planar interface.

The disk production records are `data/E5.jld2`,
`data/E5_coupled_refinement.jld2`, and `data/E5_bulk_fine.jld2`. They use
constant corrosion current `i_* = 1`. Historical records and their source are preserved under
[`archive/variable_i0/`](archive/variable_i0/); they must not be substituted
for the revised production records.

## Regenerate the paper figures and tables

```bash
bash scripts/make_paper_figures.sh
julia --project=. scripts/make_tables.jl
```

These commands read the versioned records in `data/` and write all four
current PDFs to `outputs/figures/` and tables to `outputs/tables/`.
To regenerate the included PDFs directly, use
`bash scripts/make_paper_figures.sh figures/paper`.

| Figure or table | Input record |
| --- | --- |
| F6: `fig_corrector.pdf` | E5 |
| F7: `fig_bulk.pdf` | E5 bulk fields on the doubled `1024 x 256` mesh |
| `fig_heterogeneous.pdf` | Heterogeneous experiment; limiting interface at `t=1` from its long-time record |
| `fig_evolution3d.pdf` | Long-time heterogeneous coupled trajectory at `kappa=1/4`, `0 <= t <= 2.5`; vertical coordinate is time |
| `conductivity_sweep.tex` | E5 |
| `coupled_refinement.tex` and CSV | E5 coupled refinement |

For individual plots, `make_figures.jl F6 F7` accepts `--data DIR` and
`--out DIR`. Both `plot_heterogeneous.jl` and `plot_evolution3d.jl` accept
the data directory and output directory as positional arguments, defaulting
to `data/heterogeneous` and `outputs/figures`. The table script accepts
`--data DIR` and `--out DIR`.

Recreate all four browser previews from the included PDFs with
`bash scripts/make_previews.sh`.

## Run new simulations

The disk production sequence is:

```bash
julia --project=. scripts/run_experiments.jl E5
julia --project=. scripts/run_coupled_refinement.jl --data outputs/data --out outputs/data
julia --project=. scripts/run_bulk_fine.jl --out outputs/data
julia --project=. scripts/make_figures.jl F6 F7 --data outputs/data
julia --project=. scripts/make_tables.jl --data outputs/data
```

E5 uses the initially unit disk, an insulated outer circle of radius 2,
`beta = 0.12`, `A1 = A2 = 1`, and
`phi_eq(x,y) = 0.30*(x^2-y^2) + 0.15*x*y`. Its production settings are
512 angular nodes, 128 radial intervals, `dt = 0.0025`, and `T = 0.5`.
The conductivity sweep is `kappa = 2^-k` for `k = 0,...,8`.

The bulk fields of `fig_bulk.pdf` repeat the E5 runs with `kappa = 1, 2^-3,
2^-6, 2^-8` and the limiting flow on the doubled `1024 x 256` mesh, with the
same `dt` and `T`. `run_bulk_fine.jl` caches each run under
`outputs/data/bulk_fine_runs/`; to run them in parallel, start one process per
run with `--only limit`, `--only k0`, `--only k3`, `--only k6`, and
`--only k8`, then call the script once more without `--only` to assemble
`E5_bulk_fine.jld2`. With `--factor 1` the runs repeat the E5 production runs
exactly. At `kappa = 2^-8` the doubled-mesh run coincides with the `space_2`
run of the coupled refinement, which `scripts/check.jl` verifies.

For a fresh heterogeneous run and its two plots, use:

```bash
julia --project=. scripts/run_heterogeneous.jl
julia --project=. scripts/run_heterogeneous_long.jl
julia --project=. scripts/plot_heterogeneous.jl outputs/heterogeneous outputs/figures
julia --project=. scripts/plot_evolution3d.jl outputs/heterogeneous outputs/figures
```

The explicit input paths select the fresh simulation. Calling either plot
script without arguments instead reads the versioned record. See
[HETEROGENEOUS.md](HETEROGENEOUS.md) for its parameters and refinements.

A reduced workflow check is available:

```bash
julia --project=. scripts/run_experiments.jl E5 --quick
julia --project=. scripts/make_figures.jl F6 F7 --quick
```

Quick data and plots go to `outputs/data/quick/` and
`outputs/figures/quick/`; the quick `fig_bulk.pdf` uses the fields stored in
the reduced E5 record. The reduced run uses `128 x 32`, `dt = 0.01`,
and conductivities through `2^-4`, at the same final time. It is not the
production experiment. Runtime depends on hardware and Julia compilation.

## Method and reference solutions

The electrolyte is meshed between a positive radial graph and the insulated
circle. Continuous piecewise affine finite elements use three-point Gauss
quadrature on each reactive edge, with the mixed potential interpolated
from its nodal values. Newton's method uses the exact Jacobian and residual
backtracking, with tolerance `2e-12`. The reaction fixes the constant mode.
The interface uses Fourier derivatives and explicit Heun steps; its fitted
mesh is rebuilt at unchanged connectivity.

The disk experiment uses two reference pairs:

1. **Matched discrete:** the limiting flow and its tangent-linear corrector
   use the same angular grid, time grid, finite-element discretization, and
   Heun stages as the coupled solution. This isolates the conductivity
   expansion of the discrete problem.
2. **Exact continuum:** write
   `r(t) = 1 - 0.12*t`,
   `q(theta) = 0.30*cos(2*theta) + 0.075*sin(2*theta)`, and
   `F(r) = -r^2/2 + 4*atan(r^2/4)`. Then
   `R0(t,theta) = r(t)` and
   `R1(t,theta) = q(theta)*(F(1)-F(r(t)))`.
   At `T = 0.5`, `R0 = 0.94`. These explicit references assess the
   conductivity expansion without numerical reference error.

Here `R1` is the radial corrector, opposite in sign to the inward
normal-graph corrector. For compatibility with the plotting schema,
`disk/limit/R_fine` and `disk/limit/W_fine` store the exact references,
also available as `R_exact` and `W_exact`. Fields ending in `_fine` no
longer mean a numerically refined reference. There is no current
`E5_refinement.jld2` corrector-time supplement; the earlier one is archived.

## Independent coupled-solution refinement

At `kappa = 2^-8` and `T = 0.5`, the refinement script compares the production
`512 x 128` coupled solve with three independent runs:

- `1024 x 256` at the production time step;
- `512 x 128` at `dt/2 = 0.00125`;
- `512 x 128` at `dt/4 = 0.000625`.

To use the included E5 record without overwriting it, run:

```bash
julia --project=. scripts/run_coupled_refinement.jl
```

This writes `E5_coupled_refinement.jld2`, `coupled_refinement.tex`, and
`coupled_refinement.csv` to `outputs/coupled_refinement/`. The `--data DIR`
and `--out DIR` options select the input and output directories. There is
no reduced `--quick` version of this production check.

All differences use Fourier projection onto the common 512-point angular
grid and the discrete Fourier Sobolev norms in `src/spectral.jl`. Every
corrected remainder subtracts the same exact continuum limiting radius and
corrector. The record includes final radii, solver summaries, discarded
Fourier-mode norms, consecutive temporal differences, source/data digests,
and execution metadata. The table reports measured discretization
sensitivity, not rigorous error bounds.

The three independent solves can also be prepared while E5 is running:

```bash
julia --project=. scripts/prepare_coupled_refinement.jl
julia --project=. scripts/run_coupled_refinement.jl --data outputs/data --cache outputs/constant_current_refinement_runs
```

The preparation script stores plain final radii and solver summaries with
source and parameter digests. The reporting script verifies these before
reuse. This optional workflow was used for the included revision; it changes
only scheduling, not the meshes, time steps, or reference solutions.

## Provenance and historical material

The revised production data use the constant-current M2 definition in
`src/materials.jl`, and E5 records the value of `i_star` and its reference
convention. The refinement script rejects an E5 record without the revised
`i_star = 1` marker and records its SHA-256 digest. See
[`data/README.md`](data/README.md) for the stored fields and
[`checksums.sha256`](checksums.sha256) for the disk and presentation
artifacts. The heterogeneous records have a separate manifest at
[`data/heterogeneous/checksums.sha256`](data/heterogeneous/checksums.sha256).
Their source digests cover the numerical solver and experiment runner;
figure scripts render those stored results without changing them.

The earlier variable-current distribution, including E1–E4, former E5
records, the corrector-time supplement, and their figures and tables, is
preserved with its source and original documentation in
`archive/variable_i0/`. Current data and figure directories contain the
current experiment artifacts; historical copies remain in that archive.
