#!/usr/bin/env bash
# Render the exact manuscript PDFs for the repository README.
set -euo pipefail
repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$repo_dir/figures/previews"
for name in fig_hero fig_transfer fig_corrector fig_bulk; do
    gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r200 \
        -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
        "-sOutputFile=$repo_dir/figures/previews/$name.png" \
        "$repo_dir/figures/paper/$name.pdf"
done
