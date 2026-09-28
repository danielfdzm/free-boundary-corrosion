#!/usr/bin/env bash
# Regenerate all current paper figures from the versioned numerical records.
set -euo pipefail
repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_dir"
figure_dir="${1:-outputs/figures}"
julia_bin="${JULIA:-julia}"
mkdir -p "$figure_dir"
"$julia_bin" --project=. scripts/make_figures.jl F6 F7 --data data --out "$figure_dir"
"$julia_bin" --project=. scripts/plot_heterogeneous.jl data/heterogeneous "$figure_dir"
"$julia_bin" --project=. scripts/plot_evolution3d.jl data/heterogeneous "$figure_dir"
