/* Search worst cases of a bivariate function, using a recursive algorithm.

   This program is open-source software distributed under the terms 
   of the GNU General Public License <http://www.fsf.org/copyleft/gpl.html>.

   Compile with:

   gcc -O3 -c random.c
   gcc -DFOO=atan2 -DUSE_xxx -O3 check_sample2.c random.o -lmpfr -lm -fopenmp

   where xxx is FLOAT, DOUBLE, LDOUBLE, or FLOAT128.

   To check functions with integer 2nd argument (pown, rootn, compoundn),
   add -DUSE_INT64.

   For NEWLIB: add -DNEWLIB (to avoid compilation error with __errno).

   You can add -DWORST to use some precomputed values to guide the search.

   You can add -DGLIBC to print the GNU libc release (with -v).

   Command-line options:
   -threshold nnn : set the effort threshold to nnn, the runtime is roughly
                    proportional to nnn
   -seed nnn   sets the random seed to nnn
   -mode k     sets the mode of the heuristic search to k (0 <= k <= 2).
               k=0: considers the maximal error in the range
               k=1: considers the average error in the range
               k=2: considers the estimated maximal error
   -v          verbose
   -rndn       rounding to nearest (default)
   -rndz       rounding towards zero
   -rndu       rounding towards +Inf
   -rndd       rounding towards -Inf
   -nthreads n uses n threads
   -worst xxx,yyy uses xxx,yyy (hexadecimal values) as worst values
   -xmin  xxx  use xxx a minimal value for x
   -xmax  xxx  use xxx a maximal value for x
   -ymin  yyy  use yyy a minimal value for y
   -ymax  yyy  use yyy a maximal value for y

   References and credit:
   * https://www.vinc17.net/research/testlibm/: worst-cases computed by
     Vincent Lefèvre.
   * the idea to sample several intervals instead of only one is due to
     Eric Schneider
*/

#ifndef _GNU_SOURCE
#define _GNU_SOURCE /* to define f128 functions */
#endif

#define RANK /* print the maximal list-rank of best values */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#ifndef NO_FLOAT128
#define MPFR_WANT_FLOAT128
#endif
/* icx does not have _Float128 */
#ifdef __INTEL_CLANG_COMPILER
#define _Float128 __float128
#endif
#include <mpfr.h>
#ifndef ICC
#include <math.h>
#else
#include <mathimf.h>
#endif
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#ifndef NO_OPENMP
#include <omp.h>
#endif
#include <float.h> /* for DBL_MAX */
#include <fenv.h>
#ifdef CRLIBM
#include "crlibm.h"
#endif
#ifdef _MSC_VER
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

/* define GLIBC to print the GNU libc version */
#ifdef GLIBC
#include <gnu/libc-version.h>
#endif

#ifdef FLOAT128
#ifndef __atan2q
extern _Quad __atan2q (_Quad,_Quad);
#endif
#ifndef __hypotq
extern _Quad __hypotq (_Quad,_Quad);
#endif
#ifndef __powq
extern _Quad __powq (_Quad,_Quad);
#endif
#endif

#ifndef atan2pif
extern float atan2pif (float, float);
#endif
#ifndef atan2pi
extern double atan2pi (double, double);
#endif
#ifndef atan2pil
extern long double atan2pil (long double, long double);
#endif
#ifndef pownf
extern float pownf (float, long long);
#endif
#ifndef pown
extern double pown (double, long long);
#endif
#ifndef pownl
extern long double pownl (long double, long long);
#endif
#ifndef rootnf
extern float rootnf (float, long long);
#endif
#ifndef rootn
extern double rootn (double, long long);
#endif
#ifndef rootnl
extern long double rootnl (long double, long long);
#endif
#ifndef compoundnf
extern float compoundnf (float, long long);
#endif
#ifndef compoundn
extern double compoundn (double, long long);
#endif
#ifndef compoundnl
extern long double compoundnl (long double, long long);
#endif

#ifndef NO_FLOAT128
#ifndef atan2pif128
extern _Float128 atan2pif128 (_Float128, _Float128);
#endif
#ifndef pownf128
extern _Float128 pownf128 (_Float128, long long);
#endif
#ifndef rootnf128
extern _Float128 rootnf128 (_Float128, long long);
#endif
#ifndef compoundnf128
extern _Float128 compoundnf128 (_Float128, long long);
#endif
#endif

#define mpfr_rootn mpfr_rootn_si
#define mpfr_compoundn mpfr_compound_si

/* rounding modes */
int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };
mpfr_rnd_t rnd2[] = { MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD };
mpfr_rnd_t rnd = MPFR_RNDN; /* default rounding mode */
int verbose = 0;

/* mode (0,1,2), if -1 set according to omp_get_thread_num() */
int use_mode = -1;

FILE *extra_file = NULL;

#ifdef NEWLIB
/* without this, newlib says: undefined reference to `__errno' */
int errno;
int* __errno () { return &errno; }
#endif

#define CAT1(X,Y) X ## Y
#define CAT2(X,Y) CAT1(X,Y)
#ifndef MPFR_FOO
#define MPFR_FOO CAT2(mpfr_,FOO)
#endif
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define NAME TOSTRING(FOO)

#if !defined(USE_FLOAT) && !defined(USE_DOUBLE) && !defined(USE_LDOUBLE) && \
    !defined(USE_FLOAT128)
#error "Define exactly one of USE_FLOAT, USE_DOUBLE, USE_LDOUBLE, USE_FLOAT128"
#endif

#ifdef USE_FLOAT
#if defined(USE_DOUBLE) || defined(USE_LDOUBLE) || defined(USE_FLOAT128)
#error "Define only one of USE_FLOAT, USE_DOUBLE, USE_LDOUBLE, USE_FLOAT128"
#endif
#define TYPEX float
#define UTYPEX uint32_t
#define FOO2 CAT2(FOO,f)
#define EMAX 128
#define EMIN -149
#define PREC 24
#define mpfr_set_type mpfr_set_flt
#define mpfr_get_type mpfr_get_flt
#define TYPEX_MAX FLT_MAX
#endif

#ifdef USE_DOUBLE
#if defined(USE_FLOAT) || defined(USE_LDOUBLE) || defined(USE_FLOAT128)
#error "Define only one of USE_FLOAT, USE_DOUBLE, USE_LDOUBLE, USE_FLOAT128"
#endif
#define TYPEX double
#define UTYPEX uint64_t
#define FOO2 FOO
#define EMAX 1024
#define EMIN -1074
#define PREC 53
#define mpfr_set_type mpfr_set_d
#define mpfr_get_type mpfr_get_d
#define TYPEX_MAX DBL_MAX
#endif

#ifdef USE_LDOUBLE
#if defined(USE_FLOAT) || defined(USE_DOUBLE) || defined(USE_FLOAT128)
#error "Define only one of USE_FLOAT, USE_DOUBLE, USE_LDOUBLE, USE_FLOAT128"
#endif
#define FOO2 CAT2(FOO,l)
#define TYPEX long double
#define UTYPEX __uint128_t
#define EMAX 16384
#define EMIN -16445
#define PREC 64
#define mpfr_set_type mpfr_set_ld
#define mpfr_get_type mpfr_get_ld
#define TYPEX_MAX LDBL_MAX
#endif

#ifdef USE_FLOAT128
#if defined(USE_FLOAT) || defined(USE_DOUBLE) || defined(USE_LDOUBLE)
#error "Define only one of USE_FLOAT, USE_DOUBLE, USE_LDOUBLE, USE_FLOAT128"
#endif
#if defined(__INTEL_COMPILER) || defined(__INTEL_CLANG_COMPILER)
#define TYPEX _Quad
#define FOO3 CAT2(FOO,q)
#define FOO2 CAT2(__,FOO3)
//extern _Quad FOO2 (_Quad, _Quad);
#define Q(x) (x ## q)
#else
#define TYPEX _Float128
#define FOO2 CAT2(FOO,f128)
#define Q(x) (x ## f128)
#endif
#define TYPEX_MAX Q(0xf.fffffffffffffffffffffffffff8p+16380)
#define UTYPEX __uint128_t
#define EMAX 16384
#define EMIN -16494
#define PREC 113
#define mpfr_set_type mpfr_set_float128
#define mpfr_get_type mpfr_get_float128
#endif

// the 2nd argument is an int64_t
#ifdef USE_INT64
#define TYPEY int64_t
#define UTYPEY uint64_t
#define mpfr_set_typeY mpfr_set_si
#define mpfr_get_typeY mpfr_get_si
#define PRECY 64
#define TYPEY_MAX -(ULONG_MAX-1)
#else
#define TYPEY TYPEX
#define UTYPEY UTYPEX
#define mpfr_set_typeY mpfr_set_type
#define mpfr_get_typeY mpfr_get_type
#define PRECY PREC
#define TYPEY_MAX TYPEX_MAX
#endif


#ifdef ICC /* we use compound to emulate compoundn */
float compoundnf (float x, long long y)
{
  return compoundf (x, (float) y);
}
double compoundn (double x, long long y)
{
  return compound (x, (double) y);
}
long double compoundnl (long double x, long long y)
{
  return compoundl (x, (long double) y);
}
extern _Float128 __compoundq (_Float128, _Float128);
_Float128 __compoundnq (_Float128 x, long long y)
{
  return __compoundq (x, (_Float128) y);
}
#endif

static void
print_type_hex (TYPEX x)
{
#ifdef USE_FLOAT
  printf ("%a", x);
#endif
#ifdef USE_DOUBLE
  printf ("%a", x);
#endif
#ifdef USE_LDOUBLE
  mpfr_t y;
  mpfr_init2 (y, PREC);
  mpfr_set_type (y, x, MPFR_RNDN);
  mpfr_printf ("%Ral", y);
  mpfr_clear (y);
#endif
#ifdef USE_FLOAT128
  mpfr_t y;
  mpfr_init2 (y, PREC);
  mpfr_set_type (y, x, MPFR_RNDN);
  mpfr_printf ("%Ra", y);
  mpfr_clear (y);
#endif
}

static void
print_type_hexY (TYPEY y)
{
#ifndef USE_INT64
  print_type_hex (y);
#else
  printf ("%ld", y);
#endif
}

typedef union { UTYPEX n; TYPEX x; } union_t;
typedef union { UTYPEY n; TYPEY y; } union_Y;

TYPEX
get_type (UTYPEX n)
{
  union_t u;
  u.n = n;
  return u.x;
}

UTYPEX
get_utype (TYPEX x)
{
  union_t u;
  u.x = x;
  return u.n;
}

TYPEY
get_typeY (UTYPEY n)
{
  union_Y u;
  u.n = n;
  return u.y;
}

UTYPEY
get_utypeY (TYPEY y)
{
  union_Y u;
  u.y = y;
  return u.n;
}

#ifdef USE_LDOUBLE
static int
is_valid (TYPEX x)
{
  UTYPEX n = get_utype (x);
  int e = (n >> 64) & 0x7fff; /* exponent */
  uint64_t s = (uint64_t) n;  /* significand */
  if (e == 0) return (s >> 63) == 0;
  else return (n >> 63) & 1;
}
#endif

#ifdef LIBMVEC
/* the libmvec results were obtained with GNU libc revision
   49e2bf58d57758df244eb621d63cedd2ab6d1971 */
/* LIBMVEC should be 128 (sse4.2, default), or 256 (avx2) or 512 (avx512f) */
#if !(LIBMVEC == 128 || LIBMVEC == 256 || LIBMVEC == 512)
#error "LIBMVEC should be 128 (sse4.2) or 256 (avx2) or 512 (avx512f)"
#else
#define	LIBMVEC_N (LIBMVEC/sizeof(TYPEX)/CHAR_BIT)
#endif
#endif

/* return FOO2(x,y) */
static TYPEX
WRAPPER (TYPEX x, TYPEY y)
{
#ifdef LIBMVEC
  TYPEX xx[LIBMVEC_N] = {x,}, yy[LIBMVEC_N] = {y,}, zz[LIBMVEC_N];
  for (int i = 0; i < LIBMVEC_N; i++)
    zz[i] = FOO2 (xx[i], yy[i]);
  return zz[0];
#else
  return FOO2 (x, y);
#endif
}

/* return the distance between FOO2(x,y) and the value computed by MPFR
   with precision 3*PREC, or if CHECK_CR is defined, with the value computed
   by MPFR with precision PREC */
static double
distance (TYPEX x, TYPEY y)
{
  mpfr_t xx, yy, zz;
  TYPEX z;
  double ret = 0;
  int underflow = 0;
  mpfr_exp_t expz = 0;

  if (isnan (x))
    return 0;
#ifndef USE_INT64
  if (isnan (y))
    return 0;
#endif

#ifdef USE_LDOUBLE
#ifndef USE_INT64
  if (!is_valid (x) || !is_valid (y))
#else
    if (!is_valid (x)) // with -DUSE_INT64, y is always valid
#endif
    return 0;
#endif

  mpfr_exp_t emax = mpfr_get_emax ();
  mpfr_set_emax (EMAX);
  z = WRAPPER (x, y);
  mpfr_init2 (xx, PREC);
  mpfr_init2 (yy, PREC+1);
#ifndef CHECK_CR
  mpfr_init2 (zz, 3*PREC);
#else
  /* To check correct rounding, we need to compute zz with the same precision,
     to avoid double rounding issues. */
  mpfr_init2 (zz, PREC);
#endif
  mpfr_set_type (xx, x, MPFR_RNDN);
#ifdef USE_INT64
  MPFR_FOO (zz, xx, y, rnd);
#else
  mpfr_set_typeY (yy, y, MPFR_RNDN);
  MPFR_FOO (zz, xx, yy, rnd);
#endif
  if (isnan (z) && mpfr_nan_p (zz))
    goto end;
  if (isnan (z) || mpfr_nan_p (zz))
  {
    ret = 0x1p1023 + 0x1p1023;
    goto end;
  }
  mpfr_set_emax (mpfr_get_emax_max ());
  if (isinf (z) && mpfr_inf_p (zz) == 0)
    {
      /* libm returns Inf, MPFR not, we assume the libm rounded to the next
         number after the largest finite number in the format */
      mpfr_set_ui_2exp (yy, 1, EMAX, MPFR_RNDN);
      if (z < 0)
        mpfr_neg (yy, yy, MPFR_RNDN);
    }
  else if (!isinf (z) && mpfr_inf_p (zz))
  {
    /* MPFR yields Inf, but libm no: we recompute the MPFR result with
       unbounded exponent */
#ifdef USE_INT64
    MPFR_FOO (zz, xx, y, rnd);
#else
    MPFR_FOO (zz, xx, yy, rnd);
#endif
    expz = mpfr_get_exp (zz);
    underflow = expz <= EMIN;
    mpfr_set_type (yy, z, MPFR_RNDN);
  }
  else
    mpfr_set_type (yy, z, MPFR_RNDN);

  if (mpfr_number_p (zz))
  {
    if (mpfr_zero_p (zz))
      underflow = 1;
    else
    {
      expz = mpfr_get_exp (zz);
      underflow = expz <= EMIN;
    }
  }
  mpfr_sub (zz, zz, yy, MPFR_RNDN);
  mpfr_abs (zz, zz, MPFR_RNDN);
  /* we divide by the ulp of the correctly rounded value (zz)
     which is 2^(expz-PREC) if zz is normalized, and 2^EMIN otherwise */
  if (mpfr_number_p (zz) && !underflow && expz - PREC >= EMIN)
    mpfr_mul_2si (zz, zz, PREC - expz, MPFR_RNDN);
  else
    mpfr_mul_2si (zz, zz, -EMIN, MPFR_RNDN);
  ret = mpfr_get_d (zz, MPFR_RNDU);
 end:
  mpfr_clear (xx);
  mpfr_clear (yy);
  mpfr_clear (zz);
  mpfr_set_emax (emax);
  return ret;
}

uint64_t threshold = 1000000;
double Threshold;

#include "random.h"

#ifdef USE_FLOAT
static inline uint32_t my_random(RandomState_t* rs_ptr, uint32_t n) {
    return bench_rand_u32(rs_ptr, n);
}
#else
#ifdef USE_DOUBLE
static inline uint64_t my_random(RandomState_t* rs_ptr, uint64_t n) {
    return bench_rand_u64(rs_ptr, n);
}
#else /* USE_LDOUBLE or USE_FLOAT128 */
static inline __uint128_t my_random(RandomState_t* rs_ptr, __uint128_t n) {
    return bench_rand_u128(rs_ptr, n);
}
#endif /* USE_DOUBLE */
#endif /* USE_FLOAT */

#ifdef USE_INT64
static inline uint64_t my_randomY(RandomState_t* rs_ptr, uint64_t n) {
    return bench_rand_u64(rs_ptr, n);
}
#else
#define my_randomY my_random
#endif

TYPEX Xmin, Xmax, Xbest;
TYPEY Ymin, Ymax, Ybest;
double Dbest = 0.0;
int Rbest = -1;
int mode_best = 0;

/* return the maximal error on [nxmin,nxmax] x [nymin,nymax],
   and update nxbest,nybest if it improves distance(x,y) */
static double
max_heuristic1 (UTYPEX nxmin, UTYPEX nxmax, UTYPEY nymin, UTYPEY nymax,
                UTYPEX *nxbest, UTYPEY *nybest, RandomState_t* rs_ptr)
{
  TYPEX x;
  TYPEY y;
  UTYPEX i, nx;
  UTYPEY ny;
  double dbest, dmax, d;
  x = get_type (*nxbest);
  y = get_typeY (*nybest);
  dbest = distance (x, y);
  dmax = 0;
  for (i = 0; i < threshold; i++)
    {
      nx = nxmin + my_random (rs_ptr, nxmax - nxmin);
      ny = nymin + my_randomY (rs_ptr, nymax - nymin);
      /* for long double, bit 63 must be set, but if we force it, we might
         espace the range [xmin,xmax] x [ymin,ymax] */
      x = get_type (nx);
      y = get_typeY (ny);
      d = distance (x, y);
      if (d > dmax)
        {
          dmax = d;
          if (d > dbest)
            {
              dbest = d;
              *nxbest = nx;
              *nybest = ny;
            }
        }
    }
  return dmax;
}

/* return the average error on [nxmin,nxmax] x [nymin,nymax],
   given (nxbest,nybest) */
static double
max_heuristic2 (UTYPEX nxmin, UTYPEX nxmax, UTYPEY nymin, UTYPEY nymax,
               UTYPEX *nxbest, UTYPEY *nybest, RandomState_t* rs_ptr)
{
  TYPEX x;
  TYPEY y;
  UTYPEX i, nx;
  UTYPEY ny;
  double dbest, d, s = 0, n = 0;
  x = get_type (*nxbest);
  y = get_typeY (*nybest);
  dbest = distance (x, y);
  for (i = 0; i < threshold; i++)
    {
      nx = nxmin + my_random (rs_ptr, nxmax - nxmin);
      ny = nymin + my_randomY (rs_ptr, nymax - nymin);
      x = get_type (nx);
      y = get_typeY (ny);
      d = distance (x, y);
      if (d != 0)
        {
          s += d;
          n ++;
        }
      if (d > dbest)
        {
          dbest = d;
          *nxbest = nx;
          *nybest = ny;
        }
    }
  if (n != 0)
    s = s / n;
  return s;
}

/* some libraries do not have log(), for example llvm */
static double mylog (double n)
{
  mpfr_t x;
  mpfr_init2 (x, 53);
  mpfr_set_d (x, n, MPFR_RNDN);
  mpfr_log (x, x, MPFR_RNDN);
  double ret = mpfr_get_d (x, MPFR_RNDN);
  mpfr_clear (x);
  return ret;
}

/* return the estimated maximal error on [nxmin,nxmax] x [nymin,nymax],
   taking into account mean and standard deviation,
   and update (nxbest,nybest) */
static double
max_heuristic3 (UTYPEX nxmin, UTYPEX nxmax, UTYPEY nymin, UTYPEY nymax,
               UTYPEX *nxbest, UTYPEY *nybest, RandomState_t* rs_ptr)
{
  TYPEX x;
  TYPEY y;
  UTYPEX i, nx;
  UTYPEY ny;
  double dbest, d, s = 0, v = 0, n = 0;
  x = get_type (*nxbest);
  y = get_typeY (*nybest);
  dbest = distance (x, y);
  for (i = 0; i < threshold; i++)
    {
      nx = nxmin + my_random (rs_ptr, nxmax - nxmin);
      ny = nymin + my_randomY (rs_ptr, nymax - nymin);
      x = get_type (nx);
      y = get_typeY (ny);
      d = distance (x, y);
      if (d != 0)
        {
          s += d;
          v += d * d;
          n ++;
        }
      if (d > dbest)
        {
          dbest = d;
          *nxbest = nx;
          *nybest = ny;
        }
    }
  /* compute mean and standard deviation */
  if (n != 0)
    {
      s = s / n;
      v = v / n - s * s;
      if (v < 0)
        v = 0;
      double sigma = sqrt (v);
      /* we got n non-zero values out of threshold, thus we should get
         n/threshold*(nxmax-nxmin)*(nymin-nymax) */
      n = n * (double) (nxmax - nxmin) * (double) (nymax - nymin)
        / (double) threshold;
      double logn = mylog (n);
      double t = sqrt (2.0 * logn);
      /* Reference: A note on the first moment of extreme order statistics
         from the normal distribution, Max Petzold,
         https://gupea.ub.gu.se/handle/2077/3092 */
      s = s + sigma * (t - (mylog (logn) + 1.3766) / (2.0 * t));
    }
  return s;
}

#ifdef NO_OPENMP
static int omp_get_thread_num (void) { return 0; }
static int omp_get_num_threads (void) { return 1; }
static void omp_set_num_threads (int i) { }
#endif

static int
mode (void)
{
  if (use_mode != -1)
    return use_mode;
  int i = omp_get_thread_num ();
  /* Modes 1,2 seem to be less efficient, thus we use only 1 thread on each. */
  if (i == 1 || i == 2)
    return i;
  return 0;
}

#ifdef STAT
unsigned long eval_heuristic = 0;
unsigned long eval_exhaustive = 0;
#endif

static double
max_heuristic (UTYPEX nxmin, UTYPEX nxmax, UTYPEY nymin, UTYPEY nymax,
               UTYPEX *nxbest, UTYPEY *nybest, RandomState_t* rs_ptr)
{
  if (nxmin == nxmax || nymin == nymax)
    return 0;

  int k = mode ();

#ifdef STAT
#pragma omp atomic update
  eval_heuristic += threshold;
#endif

  if (k == 0)
    return max_heuristic1 (nxmin, nxmax, nymin, nymax, nxbest, nybest, rs_ptr);
  else if (k == 1)
    return max_heuristic2 (nxmin, nxmax, nymin, nymax, nxbest, nybest, rs_ptr);
  else
    return max_heuristic3 (nxmin, nxmax, nymin, nymax, nxbest, nybest, rs_ptr);
}

typedef struct {
  UTYPEX nxmin, nxmax;
  UTYPEY nymin, nymax;
  double d;
#ifdef RANK
  int rank; /* worst rank */
#endif
} chunk_t;

static double
chunk_size (chunk_t c)
{
  return (double) (c.nxmax - c.nxmin) * (double) (c.nymax - c.nymin);
}

static void
chunk_swap (chunk_t *a, chunk_t *b)
{
  UTYPEX tmp;
  tmp = a->nxmin; a->nxmin = b->nxmin; b->nxmin = tmp;
  tmp = a->nxmax; a->nxmax = b->nxmax; b->nxmax = tmp;
  UTYPEY tmpy;
  tmpy = a->nymin; a->nymin = b->nymin; b->nymin = tmpy;
  tmpy = a->nymax; a->nymax = b->nymax; b->nymax = tmpy;
  double ump;
  ump = a->d; a->d = b->d; b->d = ump;
#ifdef RANK
  int rmp;
  rmp = a->rank; a->rank = b->rank; b->rank = rmp;
#endif
}

typedef struct {
  chunk_t *l;
  int size;
} List_struct;
typedef List_struct List_t[1];

/* Idea from Eric Schneider: instead of sampling only one interval at each
   level of the binary splitting tree, we sample up to LIST_ALLOC intervals,
   and keep the most promising ones. We use the default value suggested by
   Eric Schneider (20). */
#define LIST_ALLOC 20

static void
List_init (List_t L)
{
  L->l = (chunk_t*) malloc (LIST_ALLOC * sizeof (chunk_t));
  L->size = 0;
}

static void
List_init2 (List_t L, UTYPEX nxmin, UTYPEX nxmax, UTYPEY nymin, UTYPEY nymax)
{
  L->l = (chunk_t*) malloc (LIST_ALLOC * sizeof (chunk_t));
  L->l[0].nxmin = nxmin;
  L->l[0].nxmax = nxmax;
  L->l[0].nymin = nymin;
  L->l[0].nymax = nymax;
  L->size = 1;
}

#if 0
static void
List_print (List_t L)
{
  for (int i = 0; i < L->size; i++)
    printf ("%.3g ", L->l[i].d);
  printf ("\n");
}

static void
List_check (List_t L)
{
  for (int i = 1; i < L->size; i++)
    if (L->l[i-1].d < L->l[i].d)
      {
        fprintf (stderr, "Error, list not sorted:\n");
        List_print (L);
        exit (1);
      }
}
#endif

static void
List_insert (List_t L, UTYPEX nxmin, UTYPEX nxmax, UTYPEY nymin, UTYPEY nymax,
             double d)
{
  if (nxmin == nxmax || nymin == nymax)
    return;
  int i = L->size;
  if (i < LIST_ALLOC || d > L->l[LIST_ALLOC-1].d)
    {
      /* if list is not full, we insert at position i,
         otherwise we insert at position LIST_ALLOC-1 */
      if (i == LIST_ALLOC)
        i--;
      else
        L->size++;
      L->l[i].nxmin = nxmin;
      L->l[i].nxmax = nxmax;
      L->l[i].nymin = nymin;
      L->l[i].nymax = nymax;
      L->l[i].d = d;
      /* now insertion sort */
      while (i > 0 && d > L->l[i-1].d)
        {
          chunk_swap (L->l + (i-1), L->l + i);
#ifdef RANK
          if (L->l[i].rank < i)
            L->l[i].rank = i; /* updates maximal rank */
#endif
          i--;
        }
#ifdef RANK
      L->l[i].rank = i; /* set initial rank */
#endif
    }
  // List_check (L);
}

static void
List_swap (List_t L1, List_t L2)
{
  chunk_t *tmp;
  tmp = L1->l; L1->l = L2->l; L2->l = tmp;
  int ump;
  ump = L1->size; L1->size = L2->size; L2->size = ump;
}

static void
List_clear (List_t L)
{
  free (L->l);
}

static void
exhaustive_search (chunk_t *c, UTYPEX *nxbest, UTYPEY *nybest, double *dbest,
                   int *rbest)
{
  for (UTYPEX nx = c->nxmin; nx < c->nxmax; nx ++)
    {
      TYPEX x = get_type (nx);
      for (UTYPEY ny = c->nymin; ny < c->nymax; ny ++)
        {
          TYPEY y = get_typeY (ny);
          double d = distance (x, y);
          if (d > *dbest)
            {
              *dbest = d;
              *nxbest = nx;
              *nybest = ny;
#ifdef RANK
              *rbest = c->rank;
#endif
            }
        }
    }
#ifdef STAT
#pragma omp atomic update
  eval_exhaustive += (c->nxmax - c->nxmin) * (c->nymax - c->nymin);
#endif
}

static void
search_aux (UTYPEX nxmin, UTYPEX nxmax, UTYPEY nymin, UTYPEY nymax,
            UTYPEX *nxbest, UTYPEY *nybest, RandomState_t* rs_ptr, List_t L)
{
  double d;
  d = max_heuristic (nxmin, nxmax, nymin, nymax, nxbest, nybest, rs_ptr);
  List_insert (L, nxmin, nxmax, nymin, nymax, d);
}

/* search for nxmin <= nx < nxmax and nymin <= ny < nymax,
   where the worst found so far is (nxbest,nybest,dbest) */
static void
search (UTYPEX nxmin, UTYPEX nxmax, UTYPEY nymin, UTYPEY nymax,
        UTYPEX *nxbest, UTYPEY *nybest, double *dbest, int *rbest,
        RandomState_t* rs_ptr)
{
  int tid = omp_get_thread_num ();
  List_t L;

  List_init2 (L, nxmin, nxmax, nymin, nymax);
  while (1)
    {
      assert (1 <= L->size && L->size <= LIST_ALLOC);
      double width = chunk_size (L->l[0]);
      if (width <= Threshold) /* exhaustive search */
        {
          for (int i = 0; i < L->size; i++)
            exhaustive_search (L->l + i, nxbest, nybest, dbest, rbest);
          break;
        }
      else /* split each chunk in two */
        {
          List_t NewL;
          List_init (NewL);
          for (int i = 0; i < L->size; i++)
            {
              UTYPEX nxmin = L->l[i].nxmin;
              UTYPEX nxmax = L->l[i].nxmax;
              UTYPEX width_x = nxmax - nxmin;
              UTYPEY nymin = L->l[i].nymin;
              UTYPEY nymax = L->l[i].nymax;
              UTYPEY width_y = nymax - nymin;
              if (width_x >= width_y) // cut the x-side in two
              {
                UTYPEX nxmid = nxmin + (nxmax - nxmin) / 2;
                search_aux (nxmin, nxmid, nymin, nymax, nxbest, nybest, rs_ptr, NewL);
                search_aux (nxmid, nxmax, nymin, nymax, nxbest, nybest, rs_ptr, NewL);
              }
              else // cut the y-side in two
              {
                UTYPEY nymid = nymin + (nymax - nymin) / 2;
                search_aux (nxmin, nxmax, nymin, nymid, nxbest, nybest, rs_ptr, NewL);
                search_aux (nxmin, nxmax, nymid, nymax, nxbest, nybest, rs_ptr, NewL);
              }
            }
          List_swap (L, NewL);
#ifdef STAT
          printf ("L: "); List_print (L);
#endif
          List_clear (NewL);
        }
    }
  List_clear (L);
}

#if defined(USE_FLOAT) || defined(USE_DOUBLE) || defined(USE_LDOUBLE)
#define NLIBS 12
#else
#define	NLIBS 3 /* only glibc and icc do provide binary128 */
#endif
#define SIZE (7*NLIBS) /* 7 functions (atan2, hypot, pow, atan2pi, compoundn,
                          pown, rootn) */

#ifdef WORST
/* since CRLIBM is correctly rounded, we should not find any error larger
   than 0.5 ulp for rounding to nearest, or 1 ulp for directed rounding */
#if defined(GLIBC) || defined(CRLIBM)
#define NUMBER 0
#endif
#ifdef ICC
#define NUMBER 1
#endif
#ifdef AMD
#define NUMBER 2
#endif
#ifdef NEWLIB
#define NUMBER 3
#endif
#ifdef OPENLIBM
#define NUMBER 4
#endif
#ifdef MUSL
#define NUMBER 5
#endif
#ifdef APPLE
#define NUMBER 6
#endif
#ifdef LLVM
#ifndef USE_FLOAT128
#define NUMBER 7
#else
#define NUMBER 2
#endif
#endif
#ifdef LIBMVEC
#define NUMBER 8
#endif
#ifdef MSVC
#define NUMBER 9
#endif
#ifdef FREEBSD
#define NUMBER 10
#endif
#ifdef ARM
#define NUMBER 11
#endif
#ifndef NUMBER
#error "NUMBER undefined"
#endif
#endif

static void
setround (int rnd)
{
#ifndef CRLIBM
  fesetround (rnd);
#endif
}

static void
doit (unsigned int seed)
{
  UTYPEX nxbest = 0;
  UTYPEY nybest = 0;
  double dbest = 0;
  int rbest = -1;

  const int tid = omp_get_thread_num();
  RandomState_t randomState;
  /* set the random seed of the current thread */
  bench_srand(&randomState, seed + tid);

  /* set the rounding mode of the current thread */
  setround (rnd1[rnd]);

  /* for thread 0, get the "worst" value so far as initial point */
  if (tid == 0)
    {
      dbest = Dbest;
      nxbest = get_utype (Xbest);
      nybest = get_utypeY (Ybest);
    }

  search (get_utype (Xmin), get_utype (Xmax) + 1,
          get_utypeY (Ymin), get_utypeY (Ymax) + 1,
          &nxbest, &nybest, &dbest, &rbest, &randomState);
#pragma omp critical
  if (dbest > Dbest)
    {
      Dbest = dbest;
#ifdef RANK
      Rbest = rbest;
#endif
      mode_best = mode ();
      Xbest = get_type (nxbest);
      Ybest = get_typeY (nybest);
    }
  mpfr_free_cache (); /* free cache of current thread */
}

#ifdef WORST

#define SIZE_EXTRA 50

typedef struct { TYPEX x; TYPEY y; } pair_t;

#ifdef USE_FLOAT
  pair_t worst[SIZE] = {
#ifndef USE_INT64 // to avoid compiler warning with too large y values
    /* atan2(x,y) */
    {0x1.960176p-63,0x1.f88dp+80},      /* GNU libc       0.5 */
    {-0x1.58a7ecp-118,0x1.58a7bep-123}, /* icx            0.549965 */
    {0x1.fffe24p+59,0x1.000adcp+73},    /* AMD LibM       0.583286 */
    {-0x1.f9cf48p+49,0x1.f60598p+51},   /* Newlib         1.51236 */
    {0x1.a10104p+123,0x1.99f182p+125},  /* OpenLibm       1.54065 */
    {0x1.a10104p+123,0x1.99f182p+125},  /* Musl           1.54065 */
    {-0x1.ce62cep-116,0x1.cbf9bp-113},  /* Darwin 20.4.0  0.721704 */
    {0x1.960176p-63,0x1.f88dp+80},      /* llvm 0.5 */
#if LIBMVEC==128 /* sse4.2 */
    {0x1.b9895cp+92,0x1.c6fe64p+92},    /* libmvec 3.84566 */
#elif LIBMVEC==256 /* avx2 */
    {-0x1.bddbe6p+68,0x1.a76414p+69},   /* libmvec 3.45800 */
#else /* avx512 */
    {0x1.001998p-125,0x1.ffe0e4p-123},  /* libmvec 3.86596 */
#endif
    {0x1.fffe24p+59,0x1.000adcp+73},    /* msvc */
    {0x1.a10104p+123,0x1.99f182p+125},  /* freebsd 1.54065 */
    {-0x1.9be82ap-101,0x1.92cfb2p-101}, /* ArmPL 2.92727 */

    /* hypot(x,y) */
    {-0x1.003222p-20,-0x1.6a2d58p-32},  /* GNU libc       0.500000001392678 */
    {-0x1.003222p-20,-0x1.6a2d58p-32},  /* icx            0.500000001392678 */
    {-0x1.003222p-20,-0x1.6a2d58p-32},  /* AMD LibM       0.500000001392678 */
    {-0x1.6b05c4p-127,0x1.6b3146p-126}, /* Newlib         1.20805 */
    {-0x1.6b05c4p-127,0x1.6b3146p-126}, /* OpenLibm       1.20806 */
    {0x1.26b188p-127,-0x1.a4f2fp-128},  /* Musl           0.926707 */
    {-0x1.003222p-20,-0x1.6a2d58p-32},  /* Darwin 20.4.0  0.500001 */
    {-0x1.65e2cp+124,0x1.ac1p+112},     /* llvm           */
    /* https://reviews.llvm.org/D114726 */
#if LIBMVEC==128 /* sse4.2 */
    {-0x1.6af304p-2,-0x1.01d9d4p-2},    /* libmvec 1.4287 */
#elif LIBMVEC==256 /* avx2 */
    {-0x1.6a11aap+20,-0x1.a9216ep+12},  /* libmvec 1.20676 */
#else /* avx512 */
    {0x1.6be8b4p-62,-0x1.b8703cp-66},   /* libmvec 1.29410 */
#endif
    {-0x1.003222p-20,-0x1.6a2d58p-32},  /* msvc */
    {-0x1.6b05c4p-127,0x1.6b3146p-126}, /* freebsd 1.20806 */
    {-0x1.003222p-20,-0x1.6a2d58p-32},  /* ArmPL */

    /* pow(x,y) */
    {0x1.025736p+0,0x1.309f94p+13},     /* GNU libc       0.816737 */
    {0x1.fe7782p-1,-0x1.c361cap+14},    /* icc            0.514374 */
    {0x1.10fff4p+0,0x1.58fd76p+10},     /* AMD LibM       1.55284 */
    {0x1.d6411cp-102,0x1.793482p+0},    /* Newlib         1.00000 */
    /* the following is with -D__OBSOLETE_MATH_DEFAULT=1 */
    // {0x1.025736p+0,0x1.309f94p+13},     /* Newlib         0.816736 */
    {0x1.343e4ep+0,0x1.af3c4p+8},       /* OpenLibm       0.969482 */
    {0x1.025736p+0,0x1.309f94p+13},     /* Musl           0.816736 */
    {0x1.95e42cp-27,0x1.369972p-14},    /* Darwin         0.513872 */
    {0x1.ffff38p-1,0x1.bd5894p+23},     /* llvm 0.500001 */
#if LIBMVEC==128 /* sse4.2 */
    {0x1.b85816p+89,-0x1.90194ap-32},   /* libmvec 0.62007 */
#elif LIBMVEC==256 /* avx2 */
    {0x1.b85816p+89,-0x1.90194ap-32},   /* libmvec 0.62007 */
#else /* avx512 */
    {0x1.ff7f02p-1,0x1.55026cp+16},     /* libmvec 0.961324 */
#endif
    {0x1.e7051p+41,-0x1.c6e076p+1},     /* msvc 1.00 */
    {0x1.343e4ep+0,0x1.af3c4p+8},       /* freebsd 0.969482 */
    {0x1.025736p+0,0x1.309f94p+13},     /* ArmPL */

    /* atan2pi(x,y) */
    {-0x1.e64632p-74,-0x1.1e706cp-95},  /* GNU libc 0.5 */
    {0x1.896284p-52,0x1.f4dfccp+63},    /* icx 0.524994 */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */
#endif // USE_INT64

    /* compoundn(x,y) */
    {0x1.d0004cp-32,5263617841566842879}, /* glibc Inf */
    {-0x1.bfbc42p-19,14034238},         /* icx 0.500001 */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */

    /* pown(x,y) */
    {0x1.ffff7cp-1,21325791},           /* glibc 1.00846 */
    {0,0},                              /* icx */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */

    /* rootn(x,y) */
    {0x1.05986cp-123,15},               /* glibc 1.83243 */
    {0,0},                              /* icx */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */
  };
TYPEX extra[SIZE_EXTRA][2] = {
  {0x1.a10104p+123,0x1.99f182p+125},
  /* the following are hard-to-round cases of hypot with many identical bits after
     the round bit */
  {0x1.900004p+34,0x1.400002p+23}, /* 45 identical bits */
  {0x1.05555p+34,0x1.bffffep+23},  /* 44 identical bits */
  {0x1.e5fffap+34,0x1.affffep+23}, /* 45 identical bits */
  {0x1.260002p+34,0x1.500002p+23}, /* 45 identical bits */
  {0x1.fffffap+34,0x1.fffffep+23}, /* 45 identical bits */
  {0x1.8ffffap+34,0x1.3ffffep+23}, /* 45 identical bits */
  {0x1.87fffcp+35,0x1.bffffep+23}, /* 47 identical bits */
  /* the following exhibit errors in a previous llvm-libc implementation,
     see https://reviews.llvm.org/D118157 */
  {0x1.faf49ep+25,0x1.480002p+23},
  {0x1.014bf2p+26,0x1.780012p+23},
  {0x1.097d9ep+27,0x1.ac0024p+23},
  {0x1.0a776ap+27,0x1.d80024p+23},
  {0x1.055176p+28,0x1.400048p+23},
  {0x1.216c96p+29,0x1.7c009p+23},
  {0x1.0bec92p+30,0x1.b4012p+23},
  {0x1.463d9ap+31,0x1.ec024p+23},
  {0x1.0d6d6ap+32,0x1.7c048p+23},
  {0x1.2ebd96p+33,0x1.d809p+23},
  {0x1.ffffecp-1,-0x1.000002p+27}, /* bug in Openlibm 0.7.5 for powf */
  {-0x1.8p-49,0x1.8p+1},           /* bug in LLVM 18.1.8 powf RNDN */
  {0x1.9a1p-20,0x1.ep+0},          /* bug in LLVM 18.1.8 powf RNDZ */
  {0x1.8p-48,0x1.8p+1},            /* bug in LLVM 18.1.8 powf RNDU */
  {0x1.9a1p-20,0x1.ep+0},          /* bug in LLVM 18.1.8 powf RNDD */
  {0x1.d55902p-1,-0x1.fe037ep+9},  /* Newlib < 4.5.0 powf 169 RNDN */
  };
#endif

#ifdef USE_DOUBLE
  pair_t worst[SIZE] = {
#ifndef USE_INT64

    /* atan2(x,y) */
    {0x1.ed6060626eecfp-429,0x1.f42ebb62994dcp-426},   /* glibc 0.523171 */
    {0x1.b77ade79a36d5p-326,0x1.ff6a37b72b52bp-319},   /* icx 0.547508 */
    {-0x1.ec2003fe89e36p-433,0x1.ec173ed5fd8f5p-427},  /* AMD 1.50098 */
    {-0x1.358bb5eb25bdcp+813,0x1.2f86b82481a0ap+815},  /* Newlib 1.54908 */
    {-0x1.358bb5eb25bdcp+813,0x1.2f86b82481a0ap+815},  /* OpenLibm 1.54908 */
    {-0x1.358bb5eb25bdcp+813,0x1.2f86b82481a0ap+815},  /* Musl 1.54908 */
    {-0x1.6a539153430d8p-416,0x1.d2b5b9dc716d8p-415},  /* Apple 0.746154 */
    {-0x1.bdb2c358b0b81p-93,0x1.c717a5969be93p-86},    /* llvm 0.500065 */
#if LIBMVEC==128 /* sse4.2 */
    /* see https://sourceware.org/bugzilla/show_bug.cgi?id=28765:
       error was larger than 4 ulps */
    {0x1.abe93a1719613p-948,0x1.aab7bbb5ca811p-948},   /* libmvec 2.47316 */
#elif LIBMVEC==256 /* avx2 */
    {0x1.a83f842ef3f73p-633,0x1.a799d8a6677ep-633},    /* libmvec 3.46942 */
#else /* avx512 */
    {0x1.499c920038ab4p+559,0x1.4939bd8e01601p+559},   /* libmvec 3.41468 */
#endif
    {0x1.fffffffffffffp-627,0x1.0000000000000p+448},   /* msvc 1.00 */
    {-0x1.358bb5eb25bdcp+813,0x1.2f86b82481a0ap+815},  /* freebsd 1.54909 */
    {0x1.d5de7a294d493p-935,0x1.d030d7b608be1p-935},   /* ArmPL 2.22671 */

    /* hypot(x,y) */
    {0x0.603e52daf0bfdp-1022,-0x0.a622d0a9a433bp-1022},  /* glibc 0.791665 */
    {0x0.19deaac345ffap-1022,0x0.92c8727c389b6p-1022},   /* icc 0.7500000255880764 */
    {-0x1.a9a537a197f59p+427,0x1.254f84c458eb4p+428},    /* AMD 1.02908 */
    {0x1.6a0a41410b1abp-1004,-0x0.a24afe71b539fp-1022},  /* Newlib 1.20710 */
    {0x1.6a0a41410b1abp-1004,-0x0.a24afe71b539fp-1022},  /* OpenLibm 1.20710 */
    {0x1.00014d4b1c6b9p-1015,-0x1.000105ba9bf4p-1015},   /* Musl 1.03031 */
    {0x1.6a0a41410b1abp-1004,-0x1.4495fce36a73ep-1023},  /* Apple 1.20710 */
    /* llvm 14.0.6: tested up to target=1000000000 */
    {0x1.a308e1455f447p+0,0x1.9d931a83ef879p+0},         /* llvm */
#if LIBMVEC==128 /* sse4.2 */
    {0x1.205e5d5fee071p-9,0x1.a71193d4eb838p-9},         /* libmvec 2.66580 */
#elif LIBMVEC==256 /* avx2 */
    {0x1.6d080c1f5339ep+25,0x1.149ee0ad66632p+13},       /* libmvec 1.38804 */
#else /* avx512 */
    {-0x1.72b48b14296a7p-510,-0x1.3dcd53d99b107p-518},   /* libmvec 1.50331 */
#endif
    {-0x1.6a5a0ce661358p+890,-0x1.0151c108425b1p+890},   /* msvc 1.21768 */
    {0x1.6a0a41410b1abp-1004,-0x1.4495fce36a73ep-1023},  /* freebsd 1.2071 */
    {-0x0.5a22c27a3893p-1022,0x0.9cfea180c00dap-1022},   /* ArmPL 0.948812 */

    /* pow(x,y) */
    {0x1.010e2e7ee71aep+0,0x1.44bf0047427f6p+17},    /* glibc 0.522909 */
    {0x1.fffff9c61ce4p-1,0x1.c4e304ed4c734p+31},     /* icc 1.72811 */
    {0x1.156c99097c8cdp-753,0x1.6d3023561e86dp+0},   /* AMD 0.999324 */
    {0x1.a26d96e778dc8p+163,0x1.90d5de9f964f6p-9},   /* Newlib 0.894018 */
    // {0x1.000002c5e2e99p+0,0x1.c9eee35374af6p+31},    /* OpenLibm 635.33 */
    {0x1.a26d96e778dc8p+163,0x1.90d5de9f964f6p-9}, // OpenLibm patched 0.894019
    {0x1.010e2e7ec0c83p+0,0x1.44bf00479249dp+17},    /* Musl 0.524293 */
    {0x1.111616f835fb1p-72,0x1.c6cfa07925d49p+3},    /* Apple 0.756470 */
    {0x1.b4234939ab74ap+997,0x1.59ccddd614389p+62},  /* llvm Inf */
#if LIBMVEC==128 /* sse4.2 */
    {0x1.e174ee53813b7p+859,-0x1.d929d0607bf52p-12}, /* libmvec 1.00093 */
#elif LIBMVEC==256 /* avx2 */
    {0x1.bb393b102aa6p+246,-0x1.9c23caed44f1fp-10},  /* libmvec 0.999950 */
#else /* avx512 */
    {0x1.65f5c9d0c7bc9p-828,0x1.eba10d43b8f54p-12},  /* libmvec 0.999925 */
#endif
    {0x1.b8fa241e68c0ep-434,0x1.4ec0af8e9feb1p-18},  /* msvc 0.505 */
    {0x1.000002c5e2e99p+0,0x1.c9eee35374af6p+31},    /* freebsd 635.33 */
    {0x1.010e2e7ee71aep+0,0x1.44bf0047427f6p+17},    /* ArmPL 0.522909 */

    /* atan2pi(x,y) */
    {-0x1.07ee414995c6dp-99,0x1.9e8b07666db9ep-138},   /* GNU libc 1.48772 */
    {-0x1.c6bad820b26a9p-471,0x1.21801e79bf379p-449},  /* icx 0.527941 */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */

#else // USE_INT64

    /* compoundn(x,y) */
    {0x1.62e33ec26d2bbp-53,4611733747624586752}, /* glibc 711.264 */
    {-0x1.b7153a586c46p-58,6356776891531843}, /* icx 0.518808 */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */

    /* pown(x,y) */
    {0x1.4509b1baf557ep-1,531},         /* glibc 0.507087 */
    {0,0},                              /* icx */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */

    /* rootn(x,y) */
    {0x1.e23e265264585p-949,41},        /* glibc 1.90258 */
    {0,0},                              /* icx */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */

#endif // USE_INT64
  };
TYPEX extra[SIZE_EXTRA][2] = {
  {0x1.9c6d80093ebf8p+446,0x1.add8969f3a07ep+421},    /* as_hypot 0.5000000000000002 */
  {-0x0.5a22c27a3893p-1022,0x0.9cfea180c00dap-1022},  /* llvm with patch */
  /* the following is for glibc hypot (patch series 5463 without fma) */
  {0x0.603e52daf0bfdp-1022,-0x0.a622d0a9a433bp-1022}, /* 0.791665 */
  /* the following is for glibc hypot (patch series 5463 with fma) */
  {-0x0.5a22c27a3893p-1022,0x0.9cfea180c00dap-1022}, /* 0.948812 */
  {0x1.171436ab898ccp-796,-0x1.db3136be02ca6p-63},   /* bug crlibm pow_rn */
  /* the following is a hard-to-round case from Carlos Borges's article
     "Algorithm 1014: An Improved Algorithm for hypot(x,y)" */
  {0x1.a308e1455f447p+0,0x1.9d931a83ef879p+0},
  {0x0.89c127b351f2bp-1022,0x0.4fd9b86bb4fe8p-1022},
  /* the following gave a huge error for pow in Newlib:
     https://sourceware.org/pipermail/newlib/2023/020164.html */
  {-0x1.647ff80007ff8p-576,-0x1.3d018267f12fp+48},
  /* the following y,x pair, kindly communicated by Nikita Astafev, should give an error of
     1.76622 under CUDA. It is a worst case of x/y. */
  {0x1.2632db402a593p+52,0x1.255022bb76d28p+52},
  /* the following yields an error of 2.13e+85 ulps with pow, LLVM 21.1.8 and
     rndz (result is -0x1.5d4592cfe84d6p-791 instead of 0) */
  {0x0.fdcd1e1b3c0ccp-1022,0x1.1fa923f758348p+312},
  /* the following yields an error of inf ulps with pow, LLVM 21.1.8 and
     rndd (result is 0x1.fffffffffffffp+1023 instead of -inf) */
  {-0x0.0000000000003p-1022,-0x1.fcc4000000002p+51},
  };
#endif

#ifdef USE_LDOUBLE
  pair_t worst[SIZE] = {
#ifndef USE_INT64 // to avoid compiler warning with too large y values
/* atan2(x,y) */
{-0x7.9301460b8463cbp+15368l,0xf.25cd5eb1280b4d1p+15372l},/* glibc 0.750106 */
{-0x5.c0c9cc5a59632f88p+16340l,0x5.db7810fba1ce4908p+16348l},/* icc 0.500209 */
{0,0},                                                    /* AMD Libm */
{0,0},                                                    /* Newlib */
{0x3.d34c9d81dcd29354p+5568l,0xf.3afc4f6c9f5c4a2p+5568l}, /* OpenLibm 1.68651 */
{-0x7.9301460b8463cbp+15368l,0xf.25cd5eb1280b4d1p+15372l},/* Musl 0.750106 */
    {0,0},                                                /* Apple */
    {0,0},                                                /* llvm */
    {0,0},                                                /* libmvec */
    {0,0},                                                /* msvc */
{0x3.d34c9d81dcd29354p+5568l,0xf.3afc4f6c9f5c4a2p+5568l}, /* freebsd 1.68651 */
    {0,0},                                                /* ArmPL */

    /* hypot(x,y) */
{-0x2.97b86706043d619p+7240l,0x1.8256bdd12d2e163ep+7240l},/* glibc 0.583306 */
{-0x3.00bad8a56d87a0cp-16384l,-0xe.6d794db04791398p-16388l}, /* icc 0.75001 */
{0,0},                                                    /* AMD Libm */
{0x6.0d9f0aa73c9f1438p+11112l,0x1.82fdacf00c626b3ap+11116l}, /* Newlib inf */
{0x1.73f339f61eda21dp-16384l,0x2.e45f9f9500877e2p-16384l},/* OpenLibm 0.980552 */
{0x2.00007da75fd5903cp-8960l,0x2.d42207352184bff4p-8960l},/* Musl 1.07731 */
    {0,0},                                                /* Apple */
    {0,0},                                                /* llvm */
    {0,0},                                                /* libmvec */
    {0,0},                                                /* msvc */
{-0x8.d1ab82969da23p-16392l,0x2.e05f59c2fee13bfp-16384l}, /* freebsd 1.26817 */
    {0,0},                                                /* ArmPL */

    /* pow(x,y) */
{0x2.21dda4bcec55b158p-3616l,0x7.ef1ef5fbe3df50dp-16l},   /* glibc 0.913900 */
{0xc.b80572af668bb57p+152l,-0x6.8a6d3d7b442f3c18p+4l},    /* icc 0.500654 */
{0,0},                                                    /* AMD 0.913476 */
{0,0},                                                    /* Newlib */
{0x1.98p-3072l,0xa.ab43b6dba9a6383p+16364l},              /* OpenLibm Inf */
{-0x8p-16448l,-0x1p+0l},                                  /* Musl Inf */
    {0,0},                                                /* Apple */
    {0,0},                                                /* llvm */
    {0,0},                                                /* libmvec */
    {0,0},                                                /* msvc */
/* the following FreeBSD value is with OMP_NUM_THREADS=1,
   with several threads FreeBSD's powl returns nonsense like OpenLibm:
   https://github.com/JuliaMath/openlibm/issues/222 */
    {0x1.98p-3072l,0xa.ab43b6dba9a6383p+16364l},          /* freebsd inf */
    {0,0},                                                /* ArmPL */

    /* atan2pi(x,y) */
    {0x4.8010e21a5d13ad38p+212l,0x5.bb4879ec6325337p+10940l}, /* GNU libc 1.59019 */
    {0,0},                              /* icx */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */
#endif // USE_INT64

    /* compoundn(x,y) */
    {0xd.e923ca0615c7327p-64l,9140000946556511015}, /* glibc 3.72199 */
    {-0x1.8519710bc613bacep-60l,7189822741325238},  /* icx 0.501219 */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */

    /* pown(x,y) */
    {0x3.0e19903093cf4a18p-4l,3},       /* glibc 1.23520 */
    {0,0},                              /* icx */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */

    /* rootn(x,y) */
    {0x1.14f740a2240daa3cp+14776l,165}, /* glibc 2.05886 */
    {0,0},                              /* icx */
    {0,0},                              /* AMD LibM */
    {0,0},                              /* Newlib */
    {0,0},                              /* OpenLibm */
    {0,0},                              /* Musl */
    {0,0},                              /* Darwin */
    {0,0},                              /* llvm */
    {0,0},                              /* libmvec */
    {0,0},                              /* msvc */
    {0,0},                              /* freebsd */
    {0,0},                              /* ArmPL */
  };
TYPEX extra[SIZE_EXTRA][2] = {
  {-0x4p-16320l,0x4.0000000000000008p-16384l},
  {0x4p-16320l,0x4.0000000000000008p-16384l},
  {0x4p-16384l,0x4.0000000000000008p-16320l},
  {0x4p-16384l,-0x4.0000000000000008p-16320l},
  {0x4p-16320l,-0x4.0000000000000008p-16384l},
  {-0x4p-16384l,0x4.0000000000000008p-16320l},
  {-0x4p-16384l,-0x4.0000000000000008p-16320l},
  {-0x4p-16320l,-0x4.0000000000000008p-16384l},
  {-0x5.e7da736ff042237p-16320l,-0x5.1c6dd26cc70b3cc8p-16384l},
  {0x5.480b3781eecccadp-16384l,0x4.0687a05a152d8b18p-16320l},
  {0x6.c4e6069efd8b27b8p-16320l,-0x5.21e077ba3d981c7p-16384l},
  {-0x6.2d56393a2d535118p-16384l,-0x5.974a475e5cf1823p-16320l},
  /* the following is for musl hypot (non deterministic) */
  {0x1.868c729944440d8p-16388l,0x5.a9300a1cf661ee58p-16384l}, /* 1.20301 */
  {0x1.3c563460396af68p-16388l,0x5.ba04a504ccdc064p-16384l}, /* 1.19331 */
  {-0x9.4027deaa9f77748p-16388l,-0x5.a936e8226632238p-16384l}, /* 1.19401 */
  {-0x1.f8b424b1875c998p-16388l,-0x5.b1d67270d564d3p-16384l}, /* 1.19713 */
  /* the following is for glibc hypot (patch series 5463 with/without fma) */
  {-0x2.97b86706043d619p+7240l,0x1.8256bdd12d2e163ep+7240l}, /* 0.58330 */
  };
#endif

#ifdef USE_FLOAT128
  pair_t worst[SIZE] = {
#ifndef USE_INT64 // to avoid compiler warning with too large y values
    /* atan2(x,y) */
    {Q(0x1.41df5aa214612c7e019fa6ade88p-13316),
     Q(0x5.e53b26a270a29eb9f77ef8ef7af8p-13316)},  /* glibc 1.88155 */
    {Q(-0x1.fb41ff205f5ade930a9fcbba8ea8p-16384),
     Q(0x2.23f098fd6b8799dbeb03219bfa08p-10520)},  /* icc 0.500500 */
    {Q(-0x3.c2c372668762d218p-1876),
     Q(0x1.e161b93343b169856p+6196)},   /* llvm 1.03846e+34 */
    /* hypot(x,y) */
    {Q(-0x1.80e7403e1b344c4a78edeced92e4p-16384),
     Q(-0x2.986c750d01c32e4c807c12ad685p-16384)},  /* glibc 0.791623 */
    {Q(0x8.79ec30b61f9b839fe507bbdf414p-11908),
     Q(0xb.94f6832f64d0729ebd68035ed7a8p-11908)},  /* icc 0.500324 */
    {0,0},                                         /* llvm */
    /* pow(x,y) */
    {Q(0x1.364dcbbad0512d7bacaae2a8d56bp+0),
     Q(-0xe.68759219434c37725fdf30d17d2p+12)},     /* glibc 30.2577 */
    {Q(0x4p-16496),
     Q(0x3.ffffff39c102f0aa11bb2c8a91dp-128)},     /* icc 1.395601741238769 */
    {0,0},                                         /* llvm */
    /* atan2pi(x,y) */
    {Q(-0x2.c6994d7f40fae88117550b428404p+16040),
     Q(0xe.181a2c7696ba5cf23a6eab680b2p+16040)},   /* glibc 2.80451 */
    {Q(0x0p0),Q(0x0p0)},                           /* icx */
    {0,0},                                         /* llvm */
#endif // USE_INT64

    /* compoundn(x,y) */
    {Q(0x3.cdb067015dda5a46da054ddef174p-4),4082},   /* glibc 3.79611 */
    {Q(0x5.02ff9de4fcf096b4dd848385e0ap-52),8242432133810478537}, /* icx 0.500413 */
    {0,0},                                         /* llvm */

    /* pown(x,y) */
    {Q(0x1.364dcbbad0512d7bacaae2a8d56bp+0),-59015}, /* glibc 15.6418 */
    {Q(0x0p0),0},                                    /* icx */
    {0,0},                                         /* llvm */

    /* rootn(x,y) */
    {Q(0x8.89f8c91f0d8322325a5e27930618p-15372),1875}, /* glibc 1.96185 */
    {Q(0x0p0),0},                                      /* icx */
    {0,0},                                         /* llvm */
  };
TYPEX extra[SIZE_EXTRA][2] = {
  /* the following is for glibc hypot (patch series 5463 without fma) */
  {Q(0x2.2d5faf4036d6e68566f01054612p-8192),
   Q(0x3.5738e8e2505f5d1fc2973716f05p-8192)},      /* 0.748943 */
  };
#endif
#endif

static void
print_error (double d)
{
  mpfr_t x;
  mpfr_init2 (x, 53);
  mpfr_set_d (x, d, MPFR_RNDU);
  if (d < 0.999)
    mpfr_printf ("%.3RUf", x);
  else if (d < 9.99)
    mpfr_printf ("%.2RUf", x);
  else if (d < 99.9)
    mpfr_printf ("%.1RUf", x);
  else
    mpfr_printf ("%.2RUe", x);
  mpfr_clear (x);
}

static void
init_Threshold (void)
{
  UTYPEX boundx = get_utype (Xmax) - get_utype (Xmin);
  UTYPEY boundy = get_utypeY (Ymax) - get_utypeY (Ymin);
  double w = (double) boundx * (double) boundy;
  double s = 0;              /* number of evaluations so far */
  while (s < w)
    {
      w = w / 2.0; // cut in 2
      s += 2 * (double) threshold; // use threshold evaluations in each half
    }
  Threshold = w;
}

static int
is_smaller (TYPEX x1, TYPEX x2)
{
  TYPEX zero = 0;
  if (x1 >= zero && x2 >= zero)
    return x1 <= x2;
  else if (x1 >= zero && x2 < zero)
    return 1; // negative numbers have larger encoding
  else if (x1 <= zero && x2 >= zero)
    return 0;
  else
    return -x1 <= -x2;
}

static int
in_range (TYPEX x, TYPEY y) {
  // disable for now: there is an issue with long double (e.g. OpenLibm)
  return 1;
  int in_range_x = is_smaller (Xmin, x) && is_smaller (x, Xmax);
  int in_range_y = get_utypeY (Ymin) <= get_utypeY (y) && get_utypeY (y) < get_utypeY (Ymax);
  return in_range_x && in_range_y;
}

int
main (int argc, char *argv[])
{
  unsigned int seed = 0;
  int nthreads = 0;
  double worst_input_x = 0, worst_input_y = 0;

#ifdef CRLIBM
  /* We should call the crlibm_init() function, which sets the rounding
     precision to double. */
  crlibm_init ();
#endif

  /* we should have get_utype (Xmin) = UTYPE_MIN and
     get_utype (Xmax) = UTYPE_MAX */
  Xmin = 0;
  Xmax = -TYPEX_MAX;
  Ymin = 0;
  Ymax = -TYPEY_MAX;

  while (argc >= 2 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-v") == 0)
        {
          verbose ++;
          argc --;
          argv ++;
        }
      else if (argc >= 3 && strcmp (argv[1], "-threshold") == 0)
        {
          threshold = strtoul (argv[2], NULL, 10);
          argv += 2;
          argc -= 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-seed") == 0)
        {
          /* Note that we cannot set the seed when hardware RNG is used. */
          seed = strtoul (argv[2], NULL, 10);
          argv += 2;
          argc -= 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-mode") == 0)
        {
          use_mode = strtoul (argv[2], NULL, 10);
          assert (0 <= use_mode && use_mode <= 2);
          argv += 2;
          argc -= 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-extra") == 0)
        {
          extra_file = fopen (argv[2], "r");
          if (extra_file == NULL) {
            fprintf (stderr, "Error while opening %s\n", argv[2]);
            exit (1);
          }
          argv += 2;
          argc -= 2;
#ifndef WORST
          fprintf (stderr, "Error, -extra requires -DWORST\n");
          exit (1);
#endif
        }
      else if (strcmp (argv[1], "-rndn") == 0)
        {
          rnd = (mpfr_rnd_t) 0;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "-rndz") == 0)
        {
          rnd = (mpfr_rnd_t) 1;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "-rndu") == 0)
        {
          rnd = (mpfr_rnd_t) 2;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "-rndd") == 0)
        {
          rnd = (mpfr_rnd_t) 3;
          argc --;
          argv ++;
        }
      else if (argc >= 3 && strcmp (argv[1], "-nthreads") == 0)
        {
          nthreads = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-worst") == 0)
        {
          sscanf (argv[2], "%la,%la", &worst_input_x, &worst_input_y);
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-xmin") == 0)
        {
          long double t;
          sscanf (argv[2], "%La", &t);
          Xmin = t;
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-xmax") == 0)
        {
          long double t;
          sscanf (argv[2], "%La", &t);
          Xmax = t;
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-ymin") == 0)
        {
          long double t;
          sscanf (argv[2], "%La", &t);
          Ymin = t;
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-ymax") == 0)
        {
          long double t;
          sscanf (argv[2], "%La", &t);
          Ymax = t;
          argc -= 2;
          argv += 2;
        }
      else
        {
          fprintf (stderr, "Unknown option %s\n", argv[1]);
          exit (1);
        }
    }
  assert (threshold > 0);
  /* divide threshold by LIST_ALLOC so that the total number of evaluations
     does not vary with LIST_ALLOC */
  threshold = 1 + (threshold - 1) / LIST_ALLOC;
  init_Threshold ();

#ifdef GLIBC
  if (verbose)
    {
      printf("GNU libc version: %s\n", gnu_get_libc_version ());
      printf("GNU libc release: %s\n", gnu_get_libc_release ());
    }
#endif
#ifdef _MSC_VER
  if (verbose)
    printf ("Using Microsoft math library %d\n", _MSC_VER);
#endif

  assert(sizeof(RandomState_t) == sizeofRandomState);
  if (seed == 0)
    seed = getpid ();

  if (verbose)
    printf ("Using seed %lu\n", (unsigned long) seed);
  if (verbose)
  {
    printf ("Using xmin=");
    print_type_hex (Xmin);
    printf (" xmax=");
    print_type_hex (Xmax);
    printf (" ymin=");
    print_type_hex (Ymin);
    printf (" ymax=");
    print_type_hexY (Ymax);
    printf ("\n");
  }

#ifdef WORST
  TYPEX x;
  TYPEY y;
  double d, Dbest0;
  setround (rnd1[rnd]);
  /* first check the given -worst value, if any */
  if (worst_input_x != 0 || worst_input_y != 0)
  {
    x = worst_input_x;
    y = worst_input_y;
    d = distance (x, y);
      if (d > Dbest)
        {
          Dbest = d;
          Xbest = x;
          Ybest = y;
        }
  }
  /* first check the 'worst' values for the given library, with the given
     rounding mode */
  for (int i = NUMBER; i < SIZE; i+=NLIBS)
    {
      /* FIXME: when USE_INT64 is defined, we test extra 0 values at the
         end, since SIZE is overestimated */
      x = worst[i].x;
      y = worst[i].y;
      d = distance (x, y);
      if (d > Dbest && in_range (x, y))
        {
          Dbest = d;
          Xbest = x;
          Ybest = y;
        }
    }
  Dbest0 = Dbest;

  /* then check the 'worst' values for the other libraries */
  for (int i = 0; i < SIZE; i++)
    {
      if ((i % NLIBS) == NUMBER)
        continue;
      x = worst[i].x;
      y = worst[i].y;
      d = distance (x, y);
      if (d > Dbest && in_range (x, y))
        {
          Dbest = d;
          Xbest = x;
          Ybest = y;
        }
    }

  /* then check the 'extra' values if any */
  for (int i = 0; i < SIZE_EXTRA; i++)
    {
      x = extra[i][0];
      y = extra[i][1];
      if (x == 0 && y == 0)
        break;
      d = distance (x, y);
      if (d > Dbest && in_range (x, y))
        {
          Dbest = d;
          Xbest = x;
          Ybest = y;
        }
    }

  /* then read extra file if given */
  if (extra_file != NULL)
  {
    char *buf = NULL;
    size_t buflength = 0;
    ssize_t nn;
    int extras = 0;
    while ((nn = getline(&buf, &buflength, extra_file)) >= 0) {
      if (nn > 0 && buf[0] == '#')
        continue;
#ifdef USE_FLOAT
      if (sscanf (buf, "%a,%a", &x, &y) == 2)
#elif defined(USE_DOUBLE)
      if (sscanf (buf, "%la,%la", &x, &y) == 2)
#elif defined(USE_LDOUBLE)
        if (sscanf (buf, "%La,%La", &x, &y) == 2)
#else
      fprintf (stderr, "extra_file not implemented for this type\n");
      exit (1);
#endif
      {
        d = distance (x, y);
        if (d > Dbest)
        {
          Dbest = d;
          Xbest = x;
          Ybest = y;
        }
        extras ++;
      }
    }
    printf ("tested %d numbers from extra file\n", extras);
  }
#endif

  /* Apparently Visual Studio does not properly set the number of threads. */
#pragma omp parallel
  if (nthreads <= 0)
    nthreads = omp_get_num_threads ();

  if (verbose)
    printf ("Using %d thread(s)\n", nthreads);
  omp_set_num_threads (nthreads);

#pragma omp parallel
  doit(seed);

#ifdef WORST
  if (Dbest > Dbest0)
    printf ("NEW ");
#endif
  printf ("%s %d %d ", NAME, mode_best, Rbest);
  print_type_hex (Xbest);
  printf (",");
  print_type_hexY (Ybest);
  printf (" [");
  print_error (Dbest);
  setround (FE_UPWARD);
  printf ("] %.6g %.16g\n", Dbest, Dbest);
  /* reset to the current rounding mode */
  setround (rnd1[rnd]);
#ifdef CRLIBM /* force verification for CRLIBM */
  if (1)
#else
  if (verbose)
#endif
    {
      mpfr_t xx, yy;
      TYPEY y;
      TYPEX x,z;
      z = WRAPPER (Xbest, Ybest);
      printf ("libm gives ");
      print_type_hex (z);
      printf ("\n");
      mpfr_set_emin (EMIN + 1);
      mpfr_set_emax (EMAX);
      mpfr_init2 (xx, PREC);
      mpfr_init2 (yy, PRECY);
      mpfr_set_type (xx, Xbest, MPFR_RNDN);
#ifdef USE_INT64
      y = get_typeY (Ybest);
      int ret = MPFR_FOO (xx, xx, y, rnd);
#else
      mpfr_set_typeY (yy, Ybest, MPFR_RNDN);
      int ret = MPFR_FOO (xx, xx, yy, rnd);
#endif
      mpfr_subnormalize (xx, ret, rnd2[rnd]);
      x = mpfr_get_type (xx, MPFR_RNDN);
      printf ("mpfr gives ");
      print_type_hex (x);
      printf ("\n");
#ifdef CRLIBM
      if (x != z)
        printf ("BUG in MPFR or CRLIBM\n");
#endif
      fflush (stdout);
      mpfr_clear (xx);
      mpfr_clear (yy);
    }
  fflush (stdout);
  mpfr_free_cache ();
#ifdef STAT
  printf ("eval_heuristic=%lu eval_exhaustive=%lu\n",
          eval_heuristic, eval_exhaustive);
#endif

  if (extra_file != NULL)
    fclose (extra_file);

  return 0;
}
