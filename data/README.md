# Production data used in the paper

These JLD2 files are the production records behind the paper's numerical
section. They contain arrays and dictionaries rather than serialized solver
objects. Load a record with:

```julia
using FreeBoundaryNumerics
record = load_results("data/E5.jld2")
record["disk/kappas"]
record["disk/dev_norms"]
record["disk/rem_norms"]
```

| Record | Content used in the manuscript |
| --- | --- |
| `E1.jld2` | Homogeneous disk and ellipse benchmarks; spatial and temporal refinement; variational boundary flux recovery |
| `E2.jld2` | Charge and area balances; dissolution excess; lifetime barrier; flower snapshots in `fig_hero.pdf` |
| `E3.jld2` | Frozen-mode conductivity sweep; transfer factors and remainders in `fig_transfer.pdf`; `frozen_expansion.tex` |
| `E4.jld2` | Linearized growth-rate and mode-amplitude checks reported in the verification text |
| `E5.jld2` | Coupled and limiting flows, matched and refined correctors, trace and speed remainders, bulk fields; `fig_corrector.pdf`, `fig_bulk.pdf`, `conductivity_sweep.tex` |
| `E6.jld2` | Ellipse lifespan tests, conductivity sweep, refinement checks, curvature histories and measured times in `fig_corner.pdf` |
| `E6_phase.jld2` | The 25-by-11 ellipse/contrast phase diagram in `fig_corner.pdf` |

The experiment definitions in [`src/experiments.jl`](../src/experiments.jl)
specify each stored field and parameter. Production records occupy about
56 MiB. Checksums are in [`checksums.sha256`](../checksums.sha256).

New experiments write to `outputs/data/`, which is ignored by Git. Reduced
`--quick` runs write to `outputs/data/quick/`.
