#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python "${ROOT}/example/create_toy_data.py" --verify-only
python "${ROOT}/example/run_toy.py"
python "${ROOT}/scripts/check_toy_output.py" \
  --observed "${ROOT}/example/output/toy_run" \
  --expected "${ROOT}/example/expected" \
  --data "${ROOT}/example/data"
