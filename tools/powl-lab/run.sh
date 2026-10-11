#!/bin/sh
# Manual opt-in only. Do not invoke from unattended CI.
set -eu
[ "$#" -eq 2 ] || { echo 'Usage: sh run.sh /absolute/system/libm.so.5 /absolute/PR/libm.so.5' >&2; exit 2; }
system_lib=$1
pr_lib=$2
case "$system_lib:$pr_lib" in /*:/*) ;; *) echo 'Use absolute library paths' >&2; exit 2;; esac
[ -f "$system_lib" ] && [ -f "$pr_lib" ] || { echo 'Library not found' >&2; exit 2; }
cd "$(dirname "$0")"
lab_dir=$(pwd -P)
out_dir=$(mktemp -d "$lab_dir/run.XXXXXX")
echo "Logs: $out_dir"
timeout 60 ./build/matrix "$lab_dir/build/libpowl_overlay.so" > "$out_dir/matrix.log" 2>&1
timeout 60 ./build/differential "$system_lib" "$pr_lib" "$lab_dir/build/libpowl_overlay.so" > "$out_dir/differential.log" 2>&1
echo 'Completed. A timeout or failed check terminates this script.'
