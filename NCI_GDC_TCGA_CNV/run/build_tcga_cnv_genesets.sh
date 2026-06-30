#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
exec "${PYTHON_BIN}" "${REPO_ROOT}/geneset-extractor-dev/NCI_GDC_TCGA_CNV/src/build_tcga_cnv_genesets.py" "$@"
