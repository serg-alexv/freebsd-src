#ifndef BENCH_RANDOM_H_INCLUDED
#define BENCH_RANDOM_H_INCLUDED

#include <stddef.h>
#include <stdint.h>

/* Type that holds internal state for the random number generator */
#ifndef LEGACY_RAND
typedef struct {
    uint64_t state[4];
} RandomState_t;
#else
typedef unsigned int RandomState_t;
#endif

/* The size, in bytes, of a RandomState_t */
extern const size_t sizeofRandomState;

void bench_srand(RandomState_t* rs_ptr, unsigned int seed);

uint32_t    bench_rand_u32 (RandomState_t* rs_ptr, uint32_t n);
uint64_t    bench_rand_u64 (RandomState_t* rs_ptr, uint64_t n);
#ifndef MSVC
__uint128_t bench_rand_u128(RandomState_t* rs_ptr, __uint128_t n);
#endif

#endif  /* BENCH_RANDOM_H_INCLUDED */
