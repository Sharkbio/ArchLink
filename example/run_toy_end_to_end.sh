#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"
python example/create_toy_data.py --verify-only
python archlink.py --config example/config.toy.yaml --clustering-mode fast
