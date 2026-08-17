#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
FGS_DIR="${REPO_ROOT}/FragGeneScan-master"
FGS_BIN="${FGS_DIR}/FragGeneScan"

if [[ ! -f "${FGS_DIR}/Makefile" ]]; then
    echo "ERROR: FragGeneScan Makefile not found: ${FGS_DIR}/Makefile" >&2
    exit 1
fi

CC_BIN="${CC:-cc}"

if ! command -v "${CC_BIN}" >/dev/null 2>&1; then
    echo "ERROR: C compiler '${CC_BIN}' is required to build FragGeneScan." >&2
    exit 1
fi

if ! command -v make >/dev/null 2>&1; then
    echo "ERROR: make is required to build FragGeneScan." >&2
    exit 1
fi

echo "Building FragGeneScan in ${FGS_DIR}"
make -C "${FGS_DIR}" clean
make -C "${FGS_DIR}" fgs

if [[ ! -x "${FGS_BIN}" ]]; then
    echo "ERROR: FragGeneScan build completed without an executable: ${FGS_BIN}" >&2
    exit 1
fi

echo "FragGeneScan is ready: ${FGS_BIN}"
