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
| `E5_refinement.jld2` | Additional E5 disk corrector runs at the production spatial resolution with half and quarter time steps, measured reference discrepancies, source checksums, and `refinement_budget.tex` |
| `E5_coupled_refinement.jld2` | Additional E5 disk coupled runs at `kappa = 2^-8` with doubled mesh resolution or half and quarter time steps, common-reference corrected remainders, solver summaries, source checksums, and `coupled_refinement.tex` |

The experiment definitions in [`src/experiments.jl`](../src/experiments.jl)
specify each stored field and parameter. Production records occupy about
31 MiB. Checksums are in [`checksums.sha256`](../checksums.sha256).

New experiments write to `outputs/data/`, which is ignored by Git. Reduced
`--quick` runs write to `outputs/data/quick/`.
The separate corrector refinement command writes to `outputs/refinement/`;
see [`REPRODUCING.md`](../REPRODUCING.md#independent-corrector-time-refinement).
The corrector supplementary record captures code hashes at run start and at archival.
The only intervening source changes were two comment/docstring notation
edits in `src/fem.jl`; its `provenance_note` records this explicitly.

The coupled-refinement command writes to `outputs/coupled_refinement/`;
see [`REPRODUCING.md`](../REPRODUCING.md#independent-coupled-solution-refinement).
Its data record includes the full refined radii and the fixed archived
references, and records source digests at the start of the run.
