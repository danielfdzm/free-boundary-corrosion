# Production data used in the revised paper

The E5 JLD2 records use constant corrosion current `i_* = 1` and mixed potential
`phi_eq(x,y) = 0.30*(x^2-y^2) + 0.15*x*y`. They contain arrays and
dictionaries rather than serialized solver objects. Load a record with:

```julia
using FreeBoundaryNumerics
record = load_results("data/E5.jld2")
record["i_star"]
record["disk/kappas"]
record["disk/dev_norms"]
record["disk/rem_norms"]
```

| Record | Current manuscript content |
| --- | --- |
| `E5.jld2` | Disk conductivity sweep, coupled trajectories, matched discrete limiting flow and corrector, exact continuum references, trace and speed remainders, dissolution excess, and bulk fields; `fig_corrector.pdf`, `fig_bulk.pdf`, and `conductivity_sweep.tex` |
| `E5_coupled_refinement.jld2` | Coupled disk runs at `kappa = 2^-8` with doubled mesh resolution or half and quarter time steps, comparisons against the same exact continuum references, solver summaries, discarded-mode norms, source/data digests, and `coupled_refinement.tex` |
| `heterogeneous/heterogeneous.jld2` | Lobed-interface conductivity sweep with nonconstant `i0`, unequal reaction slopes, a transported radial corrector, refined numerical references, and independent mesh/time refinements; `fig_heterogeneous.pdf`; see [HETEROGENEOUS.md](../HETEROGENEOUS.md) |
| `heterogeneous/long_time.jld2` | Limiting interface to `t = 1` and coupled trajectory at `kappa = 1/4` to `t = 2.5` on the production mesh, with histories, solver summaries, and source digests; `fig_heterogeneous.pdf` panel (a) and `fig_evolution3d.pdf` |

E5 uses `512 x 128`, `dt = 0.0025`, and `T = 0.5`, with conductivities
`2^-k` for `k = 0,...,8`. Its exact limiting radius is `0.94`. The fields
`disk/limit/R_exact` and `disk/limit/W_exact` hold the exact continuum
limiting radius and radial corrector; `R_fine` and `W_fine` are aliases
retained for plotting compatibility. They do not denote numerical
refinements in the revised records. The fields `disk/limit/R` and
`disk/limit/W` hold the matched discrete references.

The disk experiment definition in [`src/experiments.jl`](../src/experiments.jl)
and heterogeneous runners in
[`scripts/run_heterogeneous.jl`](../scripts/run_heterogeneous.jl) and
[`scripts/run_heterogeneous_long.jl`](../scripts/run_heterogeneous_long.jl)
specify the stored fields and parameters. Exact formulas and reproduction
commands are in [`REPRODUCING.md`](../REPRODUCING.md). Artifact checksums are
in [`checksums.sha256`](../checksums.sha256) and
[`heterogeneous/checksums.sha256`](heterogeneous/checksums.sha256).
`scripts/check.jl` validates both manifests and the recorded numerical
source hashes. Figure regeneration reads the stored records and does not
alter them.

The space-time plot uses the long-time coupled trajectory at `kappa=1/4`
from `heterogeneous/long_time.jld2`. Its height coordinate is time; the
simulation itself is planar.

Fresh E5 runs write to `outputs/data/`; reduced `--quick` runs write to
`outputs/data/quick/`. The independent coupled-refinement command defaults
to `outputs/coupled_refinement/`, preserving its E5 input record. Its
`E5_sha256` field identifies that input.

Historical E1–E4 records, the former E5 data and refinement supplements,
and their source and documentation are preserved in
[`archive/variable_i0/`](../archive/variable_i0/). They are not used for the
revised paper. In particular, the earlier `E5_refinement.jld2` time-refinement
record is replaced by the explicit continuum reference formulas.

The archive retains the matching historical source, documentation, and
checksum manifest. The current `data/` directory contains only the records
listed above and their documentation and checksums.
