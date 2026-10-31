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
for arg in "$@"; do
    if [[ "$arg" == "--save" ]] && command -v timeout >/dev/null 2>&1; then
        exec timeout --signal=TERM --kill-after=5s "${SPLINE_RUN_TIMEOUT:-180s}" \
            "$python_bin" -m filmstar.cli.splines "$@"
    fi
done
exec "$python_bin" -m filmstar.cli.splines "$@"
