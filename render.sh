#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$2/figures" "$2/code"
cp "$repo"/plot/*.csv "$2/figures/"
Rscript "$repo/plot/prepare_local.R" "$1" "$2/figures"
Rscript "$repo/plot/render_figures.R" "$2"
