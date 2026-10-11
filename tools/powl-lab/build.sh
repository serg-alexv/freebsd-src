#!/bin/sh
# Reconstructed standalone build recipe; not the historical full-libm command.
set -eu
[ "$(uname -s)" = FreeBSD ] || { echo 'Requires native FreeBSD amd64' >&2; exit 2; }
[ "$(uname -m)" = amd64 ] || { echo 'Requires amd64/ld80' >&2; exit 2; }
cd "$(dirname "$0")"
mkdir -p build
cc -shared -fPIC -O2 -I. e_powl.stage3.c -lm -o build/libpowl_overlay.so
cc -O2 -fno-builtin-powl -I/usr/local/include -L/usr/local/lib -Wl,-rpath,/usr/local/lib powl_stage3_probe.c -lmpfr -lgmp -lm -o build/matrix
cc -O2 -fno-builtin-powl -I/usr/local/include -L/usr/local/lib -Wl,-rpath,/usr/local/lib powl_36shots.c -lmpfr -lgmp -lm -o build/differential
echo 'Built standalone overlay and harnesses; no numerical tests executed.'
