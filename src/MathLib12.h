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

void st_rng_set_seed(uint32_t seed);
int st_rng_is_seeded(void);

typedef double st_optim_function(int, double *, void *);
typedef void st_optim_gradient(int, double *, double *, void *);
typedef void st_optim_monitor(int, int, double, double, const double *, int,
                              void *);

typedef struct
{
    int status, iterations, function_calls, gradient_calls;
    double value, projected_gradient;
    const char *message;
} st_lbfgsb_result;

void st_lbfgsb_minimize(int n, int history, double *x, const double *lower,
                        const double *upper, const int *bounds,
                        st_optim_function *function,
                        st_optim_gradient *gradient, void *data, int maxit,
                        double factr, double pgtol, st_optim_monitor *monitor,
                        void *monitor_data, st_lbfgsb_result *result);

#endif
