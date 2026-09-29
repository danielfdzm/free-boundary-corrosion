# Heterogeneous noncircular experiment

This additional experiment exercises the transported first corrector with
a smooth, asymmetric, nonconvex interface and unequal reaction slopes.
It supplements the exact-reference disk benchmark without changing its
records or figures.

The initial interface is the radial graph

```text
R_init(theta) = 1 + .12*cos(3*theta) + .06*sin(2*theta) + .04*cos(5*theta+.3).
```

The corrosion current, mixed potential, and parameters are

```text
i0(x,y)     = 1 + .25*x + .10*y^2
phi_eq(x,y) = .30*(x^2-y^2) + .15*x*y
A2 = 1.5, A1 = .5, beta = .12, B = 2, T = .5.
```

Thus `i0 >= .5` throughout the insulated circular container. Both the
angular transport coefficient and the material-gradient coefficient of the
radial corrector equation are nonzero. The coupled conductivity sweep uses
`kappa = 2^-k`, `k=2,...,6`, with `512 x 128` fitted meshes and Heun step
`dt=.005`. The solver and its tolerances are unchanged from E5.

## Reference and refinement conventions

The radial limiting speed is `-beta*i0(R*e_theta)*N`, where
`N=sqrt(1+(R_theta/R)^2)`. The radial corrector `W` satisfies

```text
W_t = -beta * { N*(grad(i0) dot e_theta)*W
       + (i0/N)*(R_theta*W_theta/R^2 - R_theta^2*W/R^3) }
      - beta*N*i0*A2*gamma1,
gamma1 = -M_a^(-1)*(K*u0)|_Gamma,     a=(A1+A2)*i0.
```

`evolve_limit` uses precisely this tangent equation at the same Heun stages
as the limiting flow. In particular, `W` is not merely the time integral
of the electrical forcing. Its normal component is `-W/N`, not `-W` as in
the disk benchmark.

The figure uses a numerical reference pair computed on `1024 x 256`
with `dt=.0025`. Every comparison is Fourier projected onto the common
512-point angular grid. The record also contains the matched discrete
reference pair, which isolates the conductivity expansion of the discrete
scheme. Neither pair is an exact continuum reference for this experiment.

At the smallest conductivity the coupled solve, limiting flow, and
corrector are independently recomputed with doubled mesh resolution at
fixed `dt`, and with `dt/2` and `dt/4` at fixed mesh. The CSV separates
the changes in `R_kappa`, `R_0`, and `kappa*W`; corrected remainders always
use the same refined reference. These are measured discretization
sensitivities, not rigorous a posteriori error bounds. The H0/H1 norms use
the Fourier convention documented in `src/spectral.jl`.

An additional diagnostic integrates only the electrical forcing along the
same limiting Heun trajectory. Its difference from the full `W` measures
the effect of the transport and zeroth-order terms without altering the
limiting flow or the electrical discretization.

## Results and measured discretization sensitivity

For the final conductivity pair `1/32,1/64`, the refined-reference H0/H1
deviation orders are `0.97950/0.97059`; the corrected-remainder orders are
`1.97837/1.96322`. At `kappa=1/64` the corrected remainders are
`1.12399e-5` in H0 and `4.70297e-5` in H1. Omitting the transport and
zeroth-order terms changes the corrector by `2.626%` in H0 and `4.502%`
in H1, relative to its full norm.

At `kappa=1/64`, doubling both mesh resolutions at fixed `dt=.005` gives:

| Compared quantity | H0 change | H1 change |
| --- | ---: | ---: |
| Coupled radius | 2.80624e-7 | 1.83551e-6 |
| Limiting radius | 1.22e-15 | 1.89e-13 |
| Conductivity-weighted corrector | 3.00731e-7 | 1.99089e-6 |

At fixed `512 x 128` mesh, changing `dt=.005` to `.00125` changes the
coupled radius by `7.14848e-9` in H0 and `5.56089e-8` in H1. Consecutive
halvings give temporal orders `2.00001` and `2.00004`, respectively; the
limiting radius has the same second-order behavior. All these comparisons
are included in `refinement.csv` and the full JLD2 record.

The largest individual spatial sensitivity is below `4.24%` of the smallest
corrected remainder; the largest individual temporal sensitivity is below
`0.119%`. These percentages compare each of the three quantities above
separately with the corrected remainder in the corresponding norm.

## Long-time trajectories

The figures also use `long_time.jld2`, written by
`scripts/run_heterogeneous_long.jl` with the production mesh and time step:
the limiting interface to `t = 1` and the coupled interface at `kappa = 1/4`
to `t = 2.5`. Up to `t = 0.5` both runs repeat the corresponding sweep runs
exactly. The limiting flow stops being smooth at `T_0 = 2.93`, when its rays
first focus at the sharpest convex lobe and a corner forms there; its maximal
curvature grows from `2.38` at `t = 0` to `16.3` at `t = 2.5`. The coupled
trajectory at `kappa = 1/4` sharpens in the same way (maximal curvature
`16.9` at `t = 2.5`) and is no longer resolved by 512 Fourier modes after
`t = 2.6`. No smooth radial-graph trajectory of this experiment reaches
`t = 5`.

## Included artifacts and reproduction

The versioned records are in `data/heterogeneous/`: the JLD2 sweep record,
the long-time record `long_time.jld2`, `convergence.csv`, `refinement.csv`, and their separate checksum manifest.
The included figures are `figures/paper/fig_heterogeneous.pdf` and
`figures/paper/fig_evolution3d.pdf`, with PNG previews in
`figures/previews/`. The former also shows the limiting interface at
`t = 1`; the latter shows the coupled trajectory at `kappa=1/4` for
`0 <= t <= 2.5` in space-time: its vertical axis represents time, and every
horizontal section is a planar interface.

To regenerate these figures from the included record, run from the
repository root:

```bash
julia --project=. scripts/plot_heterogeneous.jl
julia --project=. scripts/plot_evolution3d.jl
```

Both scripts default to `data/heterogeneous/` as input and
`outputs/figures/` as output; they perform no simulations. They accept an
input directory and output directory as their first two positional
arguments. All four current figures can also be regenerated with
`bash scripts/make_paper_figures.sh`.

For a new simulation and plots of its results:

```bash
julia --project=. scripts/run_heterogeneous.jl
julia --project=. scripts/run_heterogeneous_long.jl
julia --project=. scripts/plot_heterogeneous.jl outputs/heterogeneous outputs/figures
julia --project=. scripts/plot_evolution3d.jl outputs/heterogeneous outputs/figures
```

The runner writes to `outputs/heterogeneous/`, caching individual solves in
its `runs/` subdirectory. It checks their numerical source digests and
parameters before reuse. `--out DIR` changes the output location. A reduced
resolution can be selected with `--nt 256` and a separate output directory.
The long-time runner writes `long_time.jld2` to the same default directory
and also accepts `--out DIR`; it always uses the production resolution.

The stored record includes every final radius, limiting and corrector
history, Newton and geometric diagnostics, execution times, Julia version,
BLAS thread count, and SHA-256 digests of the numerical sources and runner.
`scripts/check.jl` validates these hashes and both artifact checksum
manifests. Presentation scripts use the stored results without changing the
numerical record. All four included browser previews can be regenerated
with `bash scripts/make_previews.sh`.
