/*  This program is open-source software distributed under the terms of
 *  the GNU General Public License <http://www.fsf.org/copyleft/gpl.html>.
 */

#include "random.h"

const size_t sizeofRandomState = sizeof(RandomState_t);

#ifndef LEGACY_RAND

/* Adapted from the public domain code https://prng.di.unimi.it/splitmix64.c
   written in 2015 by Sebastiano Vigna. */
static uint64_t splitmix64(uint64_t* x) {
    uint64_t z = (*x += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

/* Initializes the internal state.  The seed can be any value. */
void bench_srand(RandomState_t* rs_ptr, unsigned int seed) {
    uint64_t x = (uint64_t)seed;
    rs_ptr->state[0] = splitmix64(&x);
    rs_ptr->state[1] = splitmix64(&x);
    rs_ptr->state[2] = splitmix64(&x);
    rs_ptr->state[3] = splitmix64(&x);
}

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t bench_rand(RandomState_t* rs_ptr) {
    /* The next 8 lines of code were adapted from the public domain code
       https://prng.di.unimi.it/xoshiro256plus.c
       written in 2018 by David Blackman and Sebastiano Vigna. */
    const uint64_t result = rs_ptr->state[0] + rs_ptr->state[3];
    const uint64_t t = rs_ptr->state[1] << 17;
    rs_ptr->state[2] ^= rs_ptr->state[0];
    rs_ptr->state[3] ^= rs_ptr->state[1];
    rs_ptr->state[1] ^= rs_ptr->state[2];
    rs_ptr->state[0] ^= rs_ptr->state[3];
    rs_ptr->state[2] ^= t;
    rs_ptr->state[3] = rotl(rs_ptr->state[3], 45);
    return result;
}

uint32_t bench_rand_u32(RandomState_t* rs_ptr, uint32_t n) {
    uint32_t ret = (uint32_t)(bench_rand(rs_ptr) >> 32);
    return ret % n;
}

uint64_t bench_rand_u64(RandomState_t* rs_ptr, uint64_t n) {
    return bench_rand(rs_ptr) % n;
}

#ifndef MSVC
__uint128_t bench_rand_u128(RandomState_t* rs_ptr, __uint128_t n) {
    __uint128_t ret = bench_rand(rs_ptr);
    if (n >> 64) {
        ret = (ret << 64) | bench_rand(rs_ptr);
    }
    return ret % n;
}
#endif

#else  /* LEGACY_RAND */

#include <assert.h>
#include <stdlib.h>

#if defined(__aarch64__) && defined(__ARM_FEATURE_RNG)
  #define ARM_USE_HARDWARE_RNG
  #include <arm_acle.h>
#endif

void bench_srand(RandomState_t* rs_ptr, unsigned int seed) {
    *rs_ptr = seed;
}

static uint32_t my_rand_wrapper(RandomState_t* rs_ptr)
{
  uint32_t val = 0;

#ifdef ARM_USE_HARDWARE_RNG
  /* Utilize hardware RNG instructions on the ARM processors that support them.
   *
   * During profiling, it was discovered that almost all of the processing time
   * was spent servicing the rand calls due to continual exhaustion of the
   * available system entropy. Using hardware RNG resolved this issue.
   *
   * TODO: The rndr instruction returns a 64 bit unsigned integer, of which we
   *       throw out 32 bits in order to be compatible with the existing code
   *       as it was written with rand() in mind, which returns a 32 bit
   *       integer. Ideally, we would refactor the code so that we don't throw
   *       out the bits just to request another 32 bits to build a full 64 bit
   *       value when USE_DOUBLE is defined. */
  uint64_t temp = 0;
  int res = __rndr(&temp);
  assert(res == 0);
  val = (uint32_t)temp;
  (void)rs_ptr;  /* suppresses "maybe unused" warning */
#else
  val = rand_r(rs_ptr);
#endif

  return val;
}

#if RAND_MAX == 2147483647

uint32_t bench_rand_u32(RandomState_t* rs_ptr, uint32_t n) {
  uint32_t ret = my_rand_wrapper(rs_ptr);
  if (n > RAND_MAX)
    ret = (ret << 31) | my_rand_wrapper(rs_ptr);
  return ret % n;
}

uint64_t bench_rand_u64(RandomState_t* rs_ptr, uint64_t n) {
  uint64_t ret = my_rand_wrapper(rs_ptr);
  if (n > RAND_MAX)
    {
      ret = (ret << 31) | my_rand_wrapper(rs_ptr);
      /* now ret <= 2^62-1 */
      if (n >> 62)
        ret = (ret << 31) | my_rand_wrapper(rs_ptr);
    }
  return ret % n;
}

#ifndef MSVC
__uint128_t bench_rand_u128(RandomState_t* rs_ptr, __uint128_t n) {
  __uint128_t ret = my_rand_wrapper(rs_ptr);
  if (n > RAND_MAX)
    {
      ret = (ret << 31) | my_rand_wrapper(rs_ptr);
      /* now ret <= 2^62-1 */
      if (n >> 62)
        {
          ret = (ret << 31) | my_rand_wrapper(rs_ptr);
          /* now ret <= 2^93-1 */
          if (n >> 93)
            {
              ret = (ret << 31) | my_rand_wrapper(rs_ptr);
              /* now ret <= 2^93-1 */
              if (n >> 124)
                ret = (ret << 31) | my_rand_wrapper(rs_ptr);
            }
        }
    }
  return ret % n;
}
#endif

#elif RAND_MAX == 0x7fff /* Microsoft Visual C/C++ */

uint32_t bench_rand_u32(RandomState_t* rs_ptr, uint32_t n) {
  uint32_t ret = my_rand_wrapper(rs_ptr);
  if (n > RAND_MAX)
    {
      ret = (ret << 15) | my_rand_wrapper(rs_ptr);
      /* now ret <= 2^30-1 */
      if (n >> 30)
        ret = (ret << 15) | my_rand_wrapper(rs_ptr);
    }
  return ret % n;
}

uint64_t bench_rand_u64(RandomState_t* rs_ptr, uint64_t n) {
  uint64_t ret = my_rand_wrapper(rs_ptr);
  if (n > RAND_MAX)
    {
      ret = (ret << 15) | my_rand_wrapper(rs_ptr);
      /* now ret <= 2^30-1 */
      if (n >> 30)
        {
          ret = (ret << 15) | my_rand_wrapper(rs_ptr);
          /* now ret <= 2^45-1 */
          if (n >> 45)
            {
              ret = (ret << 15) | my_rand_wrapper(rs_ptr);
              /* now ret <= 2^60-1 */
              if (n >> 60)
                ret = (ret << 15) | my_rand_wrapper(rs_ptr);
            }
        }
    }
  return ret % n;
}

#ifndef MSVC
__uint128_t bench_rand_u128(RandomState_t* rs_ptr, __uint128_t n) {
    assert(false);  /* Not yet implemented for 128-bit type. */
    return 0U;
}
#endif

#else
#error "Unexpected value of RAND_MAX"
#endif

#endif  /* LEGACY_RAND */
