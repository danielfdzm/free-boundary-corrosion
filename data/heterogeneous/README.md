# New heterogeneous-interface records

These records were computed for the asymmetric lobed-interface experiment
described in [HETEROGENEOUS.md](../../HETEROGENEOUS.md). They are independent
of the archived variable-current runs.

- `heterogeneous.jld2`: complete sweep, limiting flow and transported radial
  corrector, refined numerical reference, independent mesh/time checks,
  diagnostics, execution metadata, and numerical-source hashes.
- `long_time.jld2`: the limiting interface to `t = 1` and the coupled
  interface at `kappa = 1/4` to `t = 2.5` on the production mesh, with
  histories, solver diagnostics, execution metadata, and numerical-source
  hashes. Up to `t = 0.5` both runs repeat the corresponding sweep runs.
- `convergence.csv`: H0/H1 deviations and corrected remainders relative to
  the refined pair, together with matched discrete corrected remainders.
- `refinement.csv`: independent changes in the coupled radius, limiting
  radius, and conductivity-weighted corrector; each corrected remainder
  uses the same refined reference.
- `checksums.sha256`: verify from the repository root with
  `shasum -a 256 -c data/heterogeneous/checksums.sha256`.

```julia
using FreeBoundaryNumerics
d = load_results("data/heterogeneous/heterogeneous.jld2")
d["kappas"]
d["errors/deviation"]       # rows H0,H1; conductivity columns
d["errors/remainder"]
d["orders/remainder"]       # consecutive orders; last column NaN
d["sensitivity/space_2/coupled"]
d["sensitivity/space_2/corrector_weighted"]
d["forcing_only/relative_difference"]
```

The production mesh is `512 x 128`, `dt=.005`; the reference uses
`1024 x 256`, `dt=.0025`. The fields `limit` and `reference` are dictionaries
containing their separate `R`, `W`, histories, and metadata. Numerical
references are not exact continuum solutions. A uniform Fourier projection
onto 512 angular points precedes every comparison between meshes.
