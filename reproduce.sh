#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")" && pwd)
if [[ -e work || -e outputs ]]; then echo 'Run in an empty analysis directory.' >&2; exit 1; fi
bash "$repo/core/code/run_reproduction.sh"
mkdir -p work/S6_sources outputs/BIOAGE-REL-01/S/S6/code
cp "$repo/supplement/S_status.json" outputs/BIOAGE-REL-01/S/
python3 "$repo/supplement/s_fetch_additional.py"
python3 "$repo/supplement/s6_fetch_more.py"
python3 - <<'PY'
from hashlib import sha256
from pathlib import Path
for path, expected in {
 'work/S6_sources/PH_B.xpt':'fc190c716226078e3c07dd9faf62346e51b5200c98f6bb6b986c22d5fca0e8a4',
 'work/S_sources/NHANES_2001_2002_MORT_2019_PUBLIC.dat':'2788c8c8a02995ea5686882b5fe65b4b9cd415016031525f8a3c1a3646c47f1b',
}.items():
 assert sha256(Path(path).read_bytes()).hexdigest()==expected, path
PY
for step in s1_scores s1_mortality s2_decomposition; do Rscript "$repo/supplement/$step.R"; done
python3 "$repo/supplement/s3_session.py"
Rscript "$repo/supplement/s4_power.R"
Rscript "$repo/supplement/s6_analysis.R"
bash "$repo/render.sh" work/S/s1_results.rds figures-build
