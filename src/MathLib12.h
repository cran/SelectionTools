#ifndef MATHLIB12_H
#define MATHLIB12_H

#include <stdint.h>

#define ST_RNG_MAXTHREADS 256

int st_rand(void);
int st_rand_stream(int stream);
int st_rand_index(int n);
int st_rand_index_stream(int stream, int n);
uint64_t st_rand_index64(uint64_t n);
uint64_t st_rand_index64_stream(int stream, uint64_t n);
double st_rand_uniform(void);
double st_rand_uniform_stream(int stream);
double st_rand_uniform_pos(void);
double st_rand_uniform_pos_stream(int stream);
void st_rng_cleanup(void);
/* Called only from the serial initialization path. */
void st_rng_set_seed(uint32_t seed);
int st_rng_is_seeded(void);

#endif
