# powl ld80 native review lab

Evidence-only companion to FreeBSD PR #2427. No production source is changed by this evidence commit. `sign-fenv.patch` is an experimental follow-up addressing existing signed exceptional returns, not a claim of a PR-introduced regression.

Pinned parent: `e80b8a3e2b8117932ea094695840aed4d32b110a` (PR head). Source base recorded by PR API: `de5fcd16650207f6f2e7c50afd1ebbaeb36612c7`. Original PR: https://github.com/freebsd/freebsd-src/pull/2427 .

## Native observations

FreeBSD 15.1-RELEASE amd64, GENERIC releng/15.1-n283562-96841ea08dcf; VirtualBox; Clang 19.1.7 (llvmorg-19.1.7-0-gcd708029e0b2), MPFR 4.2.2, GMP 6.3.0. Source `e_powl.head.c` was byte-compared against raw PR head. math_accuracy revision: c2b897fb393ed8df33452b225953ccd2003cd389.

| Evidence | What it supports |
| --- | --- |
| stage3-overlay.txt | 16 argument pairs x 4 rounding modes; 64/64 value, sign and standard FENV checks |
| accuracy-pr-t1000-s2427.txt | PR returns +Inf instead of -Inf at (-LDBL_TRUE_MIN,-1) |
| powl_36shots.log | Six repeats x two finite pairs x three libraries: identical results and INEXACT flags at those two pairs |
| accuracy-overlay0/1/2 logs | Fixed-seed worst-pair feedback, 493.223 and 531.180 ulp; local stabilization at threshold 10000 |
| *disasm.txt | Measured x87 loads, sign selection and runtime products in the overlay |

MPFR supplies expected values for the matrix (all inputs are signed powers of two). Expected standard exception flags are explicitly specified in its test table, not derived from MPFR flags. Raw x87 denormal-operand bit 0x2 is separate from INVALID, DIVBYZERO, OVERFLOW, UNDERFLOW and INEXACT.

The two finite errors are inherited at the demonstrated points: system, PR and overlay agree bit-for-bit. This does not prove absence of regressions globally or identify the internal operation losing accuracy. MPFR 256/512 ld80-rounding agreement is a consistency check, not exhaustive proof. Log bit diagnostics are not interchangeable accuracy metrics.

## Supplemental fix

The PR classifier follows parity/sign analysis and already signs its exceptional product. The follow-up signs operands before multiplication in the existing general overflow/underflow returns and in powil(). Signing afterward would produce wrong directed-rounding endpoints. Negative overflow returns -Inf in nearest/downward, but -LDBL_MAX in upward/toward-zero.

Disassembly anchors (specific to measured binaries): powl 0x230f/0x2311 fmulp and 0x2313 faddp; general signed exceptional product 0x247d fmulp; powil 0x2727 fmulp. Existing classifier product 0x242d remains unchanged. Original inlined powil exception return near 0x3368a..0x336a8 has no runtime multiplication. Volatile behavior is verified on this compiler, not asserted as a portable substitute for FENV_ACCESS/strict FP.

## Build and opt-in reproduction

Required: native FreeBSD amd64/ld80, base C compiler, MPFR/GMP development headers/libraries under /usr/local. Install dependencies separately with your normal package workflow.

`sh build.sh` builds only. The recipe is newly reconstructed and has NOT been executed in this publication step. `sh run.sh /lib/libm.so.5 /absolute/path/to/PR/libm.so.5` explicitly runs the matrix and 36-call differential with independent 60-second limits and fresh log paths. No automatic CI or self-hosted runner is configured.

Historical overlay flags: `cc -shared -fPIC -O2 -I/root/powl-headers /root/e_powl.stage3.c -lm -o /root/libpowl_stage3.so`. Harnesses used -fno-builtin-powl and indirect symbol dispatch. Search build retained flags: -O3 -DWORST -DFREEBSD -DFOO=pow -DUSE_LDOUBLE -DNO_FLOAT128 -fno-builtin-powl -fopenmp, random.c, MPFR/libm; OMP_NUM_THREADS=1; seed 2427. Iteration budgets were 1000, 10000 and 10000; feedback pairs are preserved in the three source snapshots. Do not run the default million-budget search unintentionally.

## Binary identities (binaries not distributed)

| Native object | SHA-256 |
| --- | --- |
| system /lib/libm.so.5 | e2e305ba93b76670df4373a7665a0a56daa9e1913715694cdf4d0b5e7becbcfb |
| full PR libm.so.5 | 3f3005add141b7bbe0358f873478a5de74fb5445ef9e41df6e96df96f6fe3717 |
| standalone overlay | 1cf4ce0989f642150f44f53277f4cbdb8493c2112a2b380bf4c2233fe2f1ef20 |

## Explicit provenance gaps

The original full PR-libm compiler/linker command was not retained. A lab-only unrelated header typo was corrected (uint64ptr_t to uintptr_t); its exact pathname is not recovered. The overlay links system helper functions rather than a complete follow-up libm. ATF was not successfully rerun in this lab; PR-author ATF results are distinct.

Earlier original Paul-point NaN/value/FENV observations and baseline true-min exception observations are retained session summaries, not raw transcripts in this public dossier. Do not treat them as independently reconstructible packaged evidence. The raw PR wrong-sign value and overlay correction are included. The terminated million-budget search (SIGTERM/143) is excluded. Guest clock timing is inconsistent and is not a benchmark.

## Licensing and snapshot modifications

Original notices are retained. FreeBSD/msun sources and headers retain their own licenses; this directory does not relicense them. `check_sample2.accuracy*.c`, random.c and random.h retain the upstream GNU GPL notices. A full GPLv3 text is included as COPYING.GPL3 for the unversioned GPL notices; no blanket GPL relicensing of the FreeBSD tree is intended. The check_sample2 snapshots were modified during this lab's feedback experiment by replacing the FreeBSD pow worst-pair entry; these are modified snapshots, not unmodified upstream releases. The changes and observed outputs are preserved separately. Supplemental e_powl.stage3.c changes relative to the pinned source are fully represented by sign-fenv.patch.

No credentials, private email correspondence or VM images are distributed. Draft correspondence is kept out of this public directory. SHA256SUMS covers the published companion files except itself.
