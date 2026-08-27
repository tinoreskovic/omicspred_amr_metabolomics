#!/usr/bin/env bash
set -euo pipefail
: "${TASK_INDEX:?Set TASK_INDEX to a zero-based trait index}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "$SCRIPT_DIR/06_fit_bayesian_ridge.py" "$TASK_INDEX"
