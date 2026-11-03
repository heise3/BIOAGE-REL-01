#!/bin/bash
set -euo pipefail
BIOAGE_PACKAGE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p work/source outputs/BIOAGE-REL-01
cp "$BIOAGE_PACKAGE_DIR"/code/*.py "$BIOAGE_PACKAGE_DIR"/code/*.R work/
cp "$BIOAGE_PACKAGE_DIR"/contract_v1.0.md "$BIOAGE_PACKAGE_DIR"/contract_v1.0.sha256 "$BIOAGE_PACKAGE_DIR"/implementation_lock.md "$BIOAGE_PACKAGE_DIR"/implementation_lock.md.sha256 "$BIOAGE_PACKAGE_DIR"/implementation_audit.md outputs/BIOAGE-REL-01/
python3 - <<'PY'
from pathlib import Path
import json,hashlib
p=Path('outputs/BIOAGE-REL-01');h=hashlib.sha256((p/'contract_v1.0.md').read_bytes()).hexdigest()
(p/'status.json').write_text(json.dumps(dict(project='BIOAGE-REL-01',contract_version='1.0',contract_sha256=h,layer='A1',verdict='NOT_RUN',effects_computed=False)))
PY
python3 work/fetch_sources.py
python3 work/fetch_more.py
python3 work/fetch_kidney.py
python3 - "$BIOAGE_PACKAGE_DIR/source_manifest.csv" <<'PY'
import csv,hashlib,sys
from pathlib import Path
for r in csv.DictReader(open(sys.argv[1])):
 p=Path('work/source')/r['file']
 if r['file']=='nhanes3_datafiles.html':continue
 assert p.exists(),r['file']+' missing'
 assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['file']+' source changed: inspect before proceeding'
PY
python3 work/parse_nhanes.py
Rscript work/a1.R
Rscript work/a2.R
Rscript work/b_analysis.R
Rscript work/c_mortality.R
Rscript work/d_replication.R
Rscript work/qa.R
python3 work/finalize.py
