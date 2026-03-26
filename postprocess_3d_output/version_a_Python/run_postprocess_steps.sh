#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_postprocess_steps.sh /path/to/para_file
#
# This script processes steps 16..320 with step=16 using 30 concurrent processes.

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /path/to/para_file" >&2
  exit 1
fi

para_file="$1"
script_dir="$(cd "$(dirname "$0")" && pwd)"
py_script="${script_dir}/postprocess_3d_disp.py"

if [[ ! -f "$para_file" ]]; then
  echo "Parameter file not found: $para_file" >&2
  exit 1
fi

if [[ ! -f "$py_script" ]]; then
  echo "Python script not found: $py_script" >&2
  exit 1
fi

# Read the fixed-order parameter file.
mapfile -t lines < "$para_file"
if [[ ${#lines[@]} -lt 8 ]]; then
  echo "Parameter file must have at least 8 lines" >&2
  exit 1
fi

prefix="${lines[0]}"
step_id="${lines[1]}"
n_cpu_surface="${lines[2]}"
n_cpu_z="${lines[3]}"
n_surface_nodes_per_cpu="${lines[4]}"
n_z_nodes_per_cpu="${lines[5]}"
output_dir="${lines[6]}"
output_prefix="${lines[7]}"

# Dispatch steps in parallel using xargs.
seq 16 16 320 | xargs -n 1 -P 10 -I {} bash -c '
  step=$1
  tmp=$(mktemp)
  {
    echo "'"$prefix"'"
    echo "$step"
    echo "'"$n_cpu_surface"'"
    echo "'"$n_cpu_z"'"
    echo "'"$n_surface_nodes_per_cpu"'"
    echo "'"$n_z_nodes_per_cpu"'"
    echo "'"$output_dir"'"
    echo "'"$output_prefix"'"
  } > "$tmp"
  python "'"$py_script"'" "$tmp"
  rm -f "$tmp"
' _ {}
