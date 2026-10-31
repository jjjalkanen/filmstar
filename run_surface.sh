#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -x "$script_dir/venv/bin/python" ]]; then
    python_bin="$script_dir/venv/bin/python"
else
    python_bin="$(command -v python3)"
fi
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/filmstar-matplotlib}"
export PYTHONPATH="$script_dir/src${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
for arg in "$@"; do
    if [[ "$arg" == "--save" || "$arg" == "--write-solution" ]]; then
        exec timeout --signal=TERM --kill-after=5s "${SURFACE_RUN_TIMEOUT:-180s}" \
            "$python_bin" -m filmstar.cli.surfaces "$@"
    fi
done
exec "$python_bin" -m filmstar.cli.surfaces "$@"
