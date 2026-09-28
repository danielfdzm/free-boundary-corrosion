# Reproducing the paper's numerics

## Environment

Use Julia **1.13.0**, the version recorded in the paper and `Manifest.toml`.
The environment pins FFTW.jl 1.10.0, JLD2.jl 0.6.7, and CairoMakie.jl 0.15.15,
together with their dependencies. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/check.jl
```

The check verifies the included data and figure checksums, opens every
experiment record, and runs the homogeneous shrinking-disk benchmark. The
first invocation may take longer while Julia compiles dependencies.

## Regenerate the paper figures and tables

```bash
julia --project=. scripts/make_figures.jl
julia --project=. scripts/make_tables.jl
```

These commands read `data/` and write to `outputs/figures/` and
`outputs/tables/`. They retain the exact archived manuscript PDFs in
`figures/paper/` and the included tables in `tables/`.

| Plot selector | Manuscript figure | Input records |
| --- | --- | --- |
| F1 | `fig_hero.pdf` | E2 |
| F4 | `fig_transfer.pdf` | E3 |
| F6 | `fig_corrector.pdf` | E5 |
| F7 | `fig_bulk.pdf` | E5 |

The selectors retain their names from the original plotting script. For one
figure, run `julia --project=. scripts/make_figures.jl F6`. The frozen and sweep
tables are generated from E3 and E5; the refinement-budget table and its CSV
are generated from the supplementary `E5_refinement.jld2` record. The
coupled-refinement table and CSV use `E5_coupled_refinement.jld2`.

All scripts support `--help`. The plot and table scripts accept `--data DIR`
and `--out DIR`; the experiment runner accepts `--out DIR`.

## Run new simulations

A reduced frozen-domain experiment and its transfer plot:

```bash
julia --project=. scripts/run_experiments.jl E3 --quick
julia --project=. scripts/make_figures.jl F4 --quick
```

The `--quick` flag selects reduced configurations, with fresh data and plots
under `outputs/data/quick/` and `outputs/figures/quick/`. These configurations
are for checking the workflow, not for reproducing the production tables.

The complete production sequence is:

```bash
julia --project=. scripts/run_experiments.jl E1 E2 E3 E4 E5
julia --project=. scripts/run_refinement.jl --data outputs/data --out outputs/data
julia --project=. scripts/run_coupled_refinement.jl --data outputs/data --out outputs/data
julia --project=. scripts/make_figures.jl --data outputs/data
julia --project=. scripts/make_tables.jl --data outputs/data
```

With no experiment names, the runner executes all five experiments.
Individual experiments can be run separately; their data files are
independent. Runtime depends on hardware and compilation.

## Method and reference solutions

The electrolyte is meshed between a positive radial graph and the insulated
circle of radius 2. Continuous piecewise affine finite elements use
three-point Gauss quadrature on each reactive edge, with material fields
interpolated from their nodal values. Newton's method uses the exact Jacobian
and residual backtracking, with tolerance `2e-12`. The reaction fixes the
constant mode. The interface uses Fourier derivatives and explicit Heun
steps; its fitted mesh is rebuilt at unchanged connectivity.

The comparisons use three distinct references:

1. **Exact:** homogeneous shrinking disks and inward parallel curves of
   ellipses.
2. **Matched discrete:** the limiting flow and its tangent-linear corrector,
   on the same angular grid, time grid, and finite-element discretization.
3. **Refined:** a spectral RK4 limiting flow at four times the angular
   resolution and one-quarter the time step, and a corrector at twice the
   angular and radial resolution.

The matched references isolate the small-conductivity expansion of the
discrete problem. The refined comparisons assess the common discretization
error.

## Independent corrector time refinement

The supplementary refinement record keeps the production E5 disk mesh
(`512 x 128`) and final time (`T = 0.5`) and solves the limiting flow and
its tangent-linear corrector with Heun steps `dt/2 = 0.00125` and
`dt/4 = 0.000625`. The original step is `dt = 0.0025`. Run:

```bash
julia --project=. scripts/run_refinement.jl
```

This writes `E5_refinement.jld2`, `refinement_budget.tex`, and
`refinement_budget.csv` to `outputs/refinement/`. It reads the archived E5
record without changing it. Use `--data DIR` to choose the source record
and `--out DIR` to choose the destination. This check uses the full
production spatial resolution; it has no reduced `--quick` configuration.

The table records measured differences, not rigorous error bounds. Its
temporal entry is `kappa * ||W_dt - W_dt/2||`, with `kappa = 2^-8`;
the record also includes the next halving and observed temporal orders.
The spatial entry compares the archived correctors on `512 x 128` and
`1024 x 256` meshes at the same original step. The limiting-flow entry
compares the original Heun trajectory with the archived RK4 reference
at four times the angular resolution and one-quarter the step. The
record stores both refined corrector vectors, execution parameters,
Julia version, and SHA-256 digests of the source data and numerical code.
The two conductivity-weighted `H0` differences are `2.02006e-12` and
`5.04985e-13` (observed order `2.00008`). In `H3` they are `2.99e-10`
and `4.73e-10`, so an asymptotic temporal order is not resolved there;
both are much smaller than the `H3` corrected remainder, about `6e-6`.

## Independent coupled-solution refinement

The additional coupled refinement uses the production E5 disk, material M2,
`kappa = 2^-8`, and `T = 0.5`. Starting from the archived `512 x 128` solve
with `dt = 0.0025`, it independently doubles both spatial resolutions to
`1024 x 256` at fixed time step, and halves and quarters the time step on
the original mesh. Run:

```bash
julia --project=. scripts/run_coupled_refinement.jl
```

This writes `E5_coupled_refinement.jld2`, `coupled_refinement.tex`, and
`coupled_refinement.csv` to `outputs/coupled_refinement/`. The `--data DIR`
and `--out DIR` options select the input and output directories. Original
E5 data remain unchanged. There is no reduced `--quick` configuration.

All differences use Fourier projection onto the common production angular
grid of 512 points and the discrete Fourier Sobolev norm implemented in
`src/spectral.jl`. Every corrected remainder subtracts the same archived
refined limiting flow and corrector, so changes in the remainder measure
changes in the coupled solve alone. The record also stores the spatially
refined radius before projection and norms of its discarded Fourier modes.
The tabulated differences are measured discretization sensitivity, not
rigorous error bounds. The second temporal halving assesses the observed
time order separately from the spatial comparison.

The record includes full final radii, solver summaries, all four norm
orders `H0` through `H3`, the common references, UTC start and finish times,
Julia version, thread count, the archived E5 digest, and SHA-256 digests of
the numerical source files, refinement script, and manifest.

In the archived check, the coupled spatial `H0` difference is
`1.74861e-8`, or `4.39%` of the production corrected remainder
`3.97959e-7`. The two consecutive temporal `H0` differences are
`1.30362e-10` and `3.25966e-11`, giving observed order `1.99974`.
The spatially refined corrected remainder is `4.14823e-7`; the half- and
quarter-step values are `3.98057e-7` and `3.98081e-7`, respectively.
Temporal refinement does not resolve an asymptotic order in `H3`: its
consecutive differences are `1.24e-8` and `1.75e-8`, both below `0.3%` of
the corrected `H3` remainder `6.06e-6`. The discarded fine-grid modes have
`H3` norm `5.43e-8`, below `0.9%` of that remainder. All three new runs
completed in about 381 seconds in total in the recorded environment.

## Provenance

The original production `.jld2` records, four publication PDFs, and two table fragments
were copied directly from the numerical distribution used by the manuscript.
The numerical kernels and experiment parameters are preserved. Packaging
changes remove unused plot generators, replace stale numerical references
in comments with result names, and provide paths for fresh outputs. See
[`data/README.md`](data/README.md) for the connection to the paper and
[`checksums.sha256`](checksums.sha256) for the archived artifacts.
The separate `E5_refinement.jld2` and `E5_coupled_refinement.jld2` records
and their tables were added subsequently; the original production records
remain unchanged.

Browser previews are renderings of the included PDFs. To recreate them with
Ghostscript, run `bash scripts/make_previews.sh`.
