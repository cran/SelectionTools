/* MathLib12: random number generation
 *
 * Integer ran_array algorithm described by D. E. Knuth,
 * The Art of Computer Programming, Vol. 2, Section 3.6.
 * This is an independent implementation for SelectionTools.
 *
 * Each top-level OpenMP thread has its own heap-allocated state.
 * No substantial automatic arrays are used. Nested/invalid streams and
 * allocation failures use one protected static fallback state.
 */

#include "MathLib12.h"
#include <stdlib.h>
#include <time.h>

#if defined(_OPENMP) && !defined(ST)
#define ST_RNG_MT
#include <omp.h>
#endif

#define ST_RNG_K 100
#define ST_RNG_L 37
#define ST_RNG_T 70
#define ST_RNG_Q 1009
#define ST_RNG_M (UINT32_C(1) << 30)
#define ST_RNG_MASK (ST_RNG_M - UINT32_C(1))

typedef struct {
   uint32_t x[ST_RNG_K];
   uint32_t b[ST_RNG_Q];
   int p;
} st_rng_state;

static st_rng_state *st_rng_state_ptr[ST_RNG_MAXTHREADS] = {0};
static st_rng_state st_rng_fallback;
static int st_rng_fallback_ready = 0;
static uint32_t st_rng_base_seed = 0;
static int st_rng_base_ready = 0;

static uint32_t st_rng_diff(uint32_t a, uint32_t b) {
   return (a - b) & ST_RNG_MASK;
}

static void st_rng_array(st_rng_state *s, uint32_t *a, int n) {
   int i, j;
   for (j = 0; j < ST_RNG_K; j++) a[j] = s->x[j];
   for (; j < n; j++) a[j] = st_rng_diff(a[j-ST_RNG_K],a[j-ST_RNG_L]);
   for (i = 0; i < ST_RNG_L; i++, j++)
      s->x[i] = st_rng_diff(a[j-ST_RNG_K],a[j-ST_RNG_L]);
   for (; i < ST_RNG_K; i++, j++)
      s->x[i] = st_rng_diff(a[j-ST_RNG_K],s->x[i-ST_RNG_L]);
}

static void st_rng_seed_state(st_rng_state *s, uint32_t seed) {
   uint32_t *a = s->b;
   uint32_t ss = (seed + UINT32_C(2)) & (ST_RNG_M-UINT32_C(2));
   int j, t;
   for (j = 0; j < ST_RNG_K; j++) {
      a[j] = ss;
      ss <<= 1;
      if (ss >= ST_RNG_M) ss -= ST_RNG_M-UINT32_C(2);
   }
   a[1]++;
   ss = seed & ST_RNG_MASK;
   t = ST_RNG_T-1;
   while (t) {
      for (j = ST_RNG_K-1; j > 0; j--) {
         a[j+j] = a[j];
         a[j+j-1] = 0;
      }
      for (j = 2*ST_RNG_K-2; j >= ST_RNG_K; j--) {
         a[j-(ST_RNG_K-ST_RNG_L)] =
            st_rng_diff(a[j-(ST_RNG_K-ST_RNG_L)],a[j]);
         a[j-ST_RNG_K] = st_rng_diff(a[j-ST_RNG_K],a[j]);
      }
      if (ss & UINT32_C(1)) {
         for (j = ST_RNG_K; j > 0; j--) a[j] = a[j-1];
         a[0] = a[ST_RNG_K];
         a[ST_RNG_L] = st_rng_diff(a[ST_RNG_L],a[ST_RNG_K]);
      }
      if (ss) ss >>= 1;
      else t--;
   }
   for (j = 0; j < ST_RNG_L; j++) s->x[j+ST_RNG_K-ST_RNG_L] = a[j];
   for (; j < ST_RNG_K; j++) s->x[j-ST_RNG_L] = a[j];
   for (j = 0; j < 10; j++) st_rng_array(s,a,2*ST_RNG_K-1);
   s->p = ST_RNG_K;
}

static uint32_t st_rng_seed_for_stream(int stream) {
   uint32_t seed;
#ifdef ST_RNG_MT
#pragma omp critical(st_rng_seed_lock)
#endif
   {
      if (!st_rng_base_ready) {
         st_rng_base_seed = (uint32_t)time(NULL) & ST_RNG_MASK;
         st_rng_base_ready = 1;
      }
      seed = (st_rng_base_seed + (uint32_t)(stream+1)) & ST_RNG_MASK;
   }
   return seed;
}

static void st_rng_init_fallback(void) {
#ifdef ST_RNG_MT
#pragma omp critical(st_rng_fallback_init_lock)
#endif
   {
      if (!st_rng_fallback_ready) {
         st_rng_seed_state(&st_rng_fallback,st_rng_seed_for_stream(ST_RNG_MAXTHREADS));
         st_rng_fallback_ready = 1;
      }
   }
}

static st_rng_state *st_rng_get(int stream, int *shared) {
   st_rng_state *s;
   *shared = 0;
#ifdef ST_RNG_MT
   if (omp_get_level() > 1) {
      st_rng_init_fallback();
      *shared = 1;
      return &st_rng_fallback;
   }
#endif
   if (stream < 0 || stream >= ST_RNG_MAXTHREADS) {
      st_rng_init_fallback();
      *shared = 1;
      return &st_rng_fallback;
   }
   s = st_rng_state_ptr[stream];
   if (s) return s;
#ifdef ST_RNG_MT
#pragma omp critical(st_rng_stream_init_lock)
#endif
   {
      s = st_rng_state_ptr[stream];
      if (!s) {
         s = (st_rng_state*)malloc(sizeof(*s));
         if (s) {
            st_rng_seed_state(s,st_rng_seed_for_stream(stream));
            st_rng_state_ptr[stream] = s;
         }
      }
   }
   if (s) return s;
   st_rng_init_fallback();
   *shared = 1;
   return &st_rng_fallback;
}

static int st_rng_next(st_rng_state *s) {
   if (s->p >= ST_RNG_K) {
      st_rng_array(s,s->b,ST_RNG_Q);
      s->p = 0;
   }
   return (int)s->b[s->p++];
}

int st_rand_stream(int stream) {
   st_rng_state *s;
   int shared;
   s = st_rng_get(stream,&shared);
#ifdef ST_RNG_MT
   if (shared) {
      int r;
#pragma omp critical(st_rng_fallback_use_lock)
      r = st_rng_next(s);
      return r;
   }
#else
   (void)shared;
#endif
   return st_rng_next(s);
}

int st_rand(void) {
   int stream = 0;
#ifdef ST_RNG_MT
   stream = omp_get_thread_num();
#endif
   return st_rand_stream(stream);
}

uint64_t st_rand_index64_stream(int stream, uint64_t n) {
   uint64_t r, threshold;
   if (n <= UINT64_C(1)) return 0;
   if (n <= (uint64_t)ST_RNG_M) {
      uint32_t u = (uint32_t)n, x, limit;
      limit = ST_RNG_M-(ST_RNG_M%u);
      do x = (uint32_t)st_rand_stream(stream); while (x >= limit);
      return (uint64_t)(x%u);
   }
   threshold = (UINT64_C(0)-n)%n;
   do {
      r = ((uint64_t)(uint32_t)st_rand_stream(stream) << 34) |
          ((uint64_t)(uint32_t)st_rand_stream(stream) << 4) |
          ((uint64_t)(uint32_t)st_rand_stream(stream) & UINT64_C(15));
   } while (r < threshold);
   return r%n;
}

uint64_t st_rand_index64(uint64_t n) {
   int stream = 0;
#ifdef ST_RNG_MT
   stream = omp_get_thread_num();
#endif
   return st_rand_index64_stream(stream,n);
}

int st_rand_index_stream(int stream, int n) {
   if (n <= 1) return 0;
   return (int)st_rand_index64_stream(stream,(uint64_t)(unsigned int)n);
}

int st_rand_index(int n) {
   int stream = 0;
#ifdef ST_RNG_MT
   stream = omp_get_thread_num();
#endif
   return st_rand_index_stream(stream,n);
}

double st_rand_uniform_stream(int stream) {
   return (double)(uint32_t)st_rand_stream(stream)/(double)ST_RNG_M;
}

double st_rand_uniform(void) {
   int stream = 0;
#ifdef ST_RNG_MT
   stream = omp_get_thread_num();
#endif
   return st_rand_uniform_stream(stream);
}

double st_rand_uniform_pos_stream(int stream) {
   double u;
   do u = st_rand_uniform_stream(stream); while (u <= 0.0);
   return u;
}

double st_rand_uniform_pos(void) {
   int stream = 0;
#ifdef ST_RNG_MT
   stream = omp_get_thread_num();
#endif
   return st_rand_uniform_pos_stream(stream);
}

void st_rng_cleanup(void) {
   int i;
#ifdef ST_RNG_MT
#pragma omp critical(st_rng_cleanup_lock)
#endif
   {
      for (i = 0; i < ST_RNG_MAXTHREADS; i++) {
         free(st_rng_state_ptr[i]);
         st_rng_state_ptr[i] = NULL;
      }
      st_rng_fallback_ready = 0;
      st_rng_base_ready = 0;
      st_rng_base_seed = 0;
   }
}
