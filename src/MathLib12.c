

#include "MathLib12.h"
#include "SelectionToolsAllele.h"
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>

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

typedef struct
{
    uint32_t x[ST_RNG_K];
    uint32_t b[ST_RNG_Q];
    int p;
} st_rng_state;

static st_rng_state *st_rng_state_ptr[ST_RNG_MAXTHREADS] = {0};
static st_rng_state st_rng_fallback;
static int st_rng_fallback_ready = 0;
static uint32_t st_rng_base_seed = 0;
static int st_rng_base_ready = 0;

static uint32_t st_rng_diff(uint32_t a, uint32_t b)
{
    return (a - b) & ST_RNG_MASK;
}

static void st_rng_array(st_rng_state *s, uint32_t *a, int n)
{
    int i, j;
    for (j = 0; j < ST_RNG_K; j++)
        a[j] = s->x[j];
    for (; j < n; j++)
        a[j] = st_rng_diff(a[j - ST_RNG_K], a[j - ST_RNG_L]);
    for (i = 0; i < ST_RNG_L; i++, j++)
        s->x[i] = st_rng_diff(a[j - ST_RNG_K], a[j - ST_RNG_L]);
    for (; i < ST_RNG_K; i++, j++)
        s->x[i] = st_rng_diff(a[j - ST_RNG_K], s->x[i - ST_RNG_L]);
}

static void st_rng_seed_state(st_rng_state *s, uint32_t seed)
{
    uint32_t *a = s->b;
    uint32_t ss = (seed + UINT32_C(2)) & (ST_RNG_M - UINT32_C(2));
    int j, t;
    for (j = 0; j < ST_RNG_K; j++)
    {
        a[j] = ss;
        ss <<= 1;
        if (ss >= ST_RNG_M)
            ss -= ST_RNG_M - UINT32_C(2);
    }
    a[1]++;
    ss = seed & ST_RNG_MASK;
    t = ST_RNG_T - 1;
    while (t)
    {
        for (j = ST_RNG_K - 1; j > 0; j--)
        {
            a[j + j] = a[j];
            a[j + j - 1] = 0;
        }
        for (j = 2 * ST_RNG_K - 2; j >= ST_RNG_K; j--)
        {
            a[j - (ST_RNG_K - ST_RNG_L)] =
                st_rng_diff(a[j - (ST_RNG_K - ST_RNG_L)], a[j]);
            a[j - ST_RNG_K] = st_rng_diff(a[j - ST_RNG_K], a[j]);
        }
        if (ss & UINT32_C(1))
        {
            for (j = ST_RNG_K; j > 0; j--)
                a[j] = a[j - 1];
            a[0] = a[ST_RNG_K];
            a[ST_RNG_L] = st_rng_diff(a[ST_RNG_L], a[ST_RNG_K]);
        }
        if (ss)
            ss >>= 1;
        else
            t--;
    }
    for (j = 0; j < ST_RNG_L; j++)
        s->x[j + ST_RNG_K - ST_RNG_L] = a[j];
    for (; j < ST_RNG_K; j++)
        s->x[j - ST_RNG_L] = a[j];
    for (j = 0; j < 10; j++)
        st_rng_array(s, a, 2 * ST_RNG_K - 1);
    s->p = ST_RNG_K;
}

static uint32_t st_rng_seed_for_stream(int stream)
{
    uint32_t seed;
#ifdef ST_RNG_MT
#pragma omp critical(st_rng_seed_lock)
#endif
    {
        if (!st_rng_base_ready)
        {
            st_rng_base_seed = (uint32_t)time(NULL) & ST_RNG_MASK;
            st_rng_base_ready = 1;
        }
        seed = (st_rng_base_seed + (uint32_t)(stream + 1)) & ST_RNG_MASK;
    }
    return seed;
}

static void st_rng_init_fallback(void)
{
#ifdef ST_RNG_MT
#pragma omp critical(st_rng_fallback_init_lock)
#endif
    {
        if (!st_rng_fallback_ready)
        {
            st_rng_seed_state(&st_rng_fallback,
                              st_rng_seed_for_stream(ST_RNG_MAXTHREADS));
            st_rng_fallback_ready = 1;
        }
    }
}

static st_rng_state *st_rng_get(int stream, int *shared)
{
    st_rng_state *s;
    *shared = 0;
#ifdef ST_RNG_MT
    if (omp_get_level() > 1)
    {
        st_rng_init_fallback();
        *shared = 1;
        return &st_rng_fallback;
    }
#endif
    if (stream < 0 || stream >= ST_RNG_MAXTHREADS)
    {
        st_rng_init_fallback();
        *shared = 1;
        return &st_rng_fallback;
    }
    s = st_rng_state_ptr[stream];
    if (s)
        return s;
#ifdef ST_RNG_MT
#pragma omp critical(st_rng_stream_init_lock)
#endif
    {
        s = st_rng_state_ptr[stream];
        if (!s)
        {
            s = (st_rng_state *)malloc(sizeof(*s));
            if (s)
            {
                st_rng_seed_state(s, st_rng_seed_for_stream(stream));
                st_rng_state_ptr[stream] = s;
            }
        }
    }
    if (s)
        return s;
    st_rng_init_fallback();
    *shared = 1;
    return &st_rng_fallback;
}

static int st_rng_next(st_rng_state *s)
{
    if (s->p >= ST_RNG_K)
    {
        st_rng_array(s, s->b, ST_RNG_Q);
        s->p = 0;
    }
    return (int)s->b[s->p++];
}

int st_rand_stream(int stream)
{
    st_rng_state *s;
    int shared;
    s = st_rng_get(stream, &shared);
#ifdef ST_RNG_MT
    if (shared)
    {
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

int st_rand(void)
{
    int stream = 0;
#ifdef ST_RNG_MT
    stream = omp_get_thread_num();
#endif
    return st_rand_stream(stream);
}

uint64_t st_rand_index64_stream(int stream, uint64_t n)
{
    uint64_t r, threshold;
    if (n <= UINT64_C(1))
        return 0;
    if (n <= (uint64_t)ST_RNG_M)
    {
        uint32_t u = (uint32_t)n, x, limit;
        limit = ST_RNG_M - (ST_RNG_M % u);
        do
            x = (uint32_t)st_rand_stream(stream);
        while (x >= limit);
        return (uint64_t)(x % u);
    }
    threshold = (UINT64_C(0) - n) % n;
    do
    {
        r = ((uint64_t)(uint32_t)st_rand_stream(stream) << 34) |
            ((uint64_t)(uint32_t)st_rand_stream(stream) << 4) |
            ((uint64_t)(uint32_t)st_rand_stream(stream) & UINT64_C(15));
    } while (r < threshold);
    return r % n;
}

uint64_t st_rand_index64(uint64_t n)
{
    int stream = 0;
#ifdef ST_RNG_MT
    stream = omp_get_thread_num();
#endif
    return st_rand_index64_stream(stream, n);
}

int st_rand_index_stream(int stream, int n)
{
    if (n <= 1)
        return 0;
    return (int)st_rand_index64_stream(stream, (uint64_t)(unsigned int)n);
}

int st_rand_index(int n)
{
    int stream = 0;
#ifdef ST_RNG_MT
    stream = omp_get_thread_num();
#endif
    return st_rand_index_stream(stream, n);
}

double st_rand_uniform_stream(int stream)
{
    return (double)(uint32_t)st_rand_stream(stream) / (double)ST_RNG_M;
}

double st_rand_uniform(void)
{
    int stream = 0;
#ifdef ST_RNG_MT
    stream = omp_get_thread_num();
#endif
    return st_rand_uniform_stream(stream);
}

double st_rand_uniform_pos_stream(int stream)
{
    double u;
    do
        u = st_rand_uniform_stream(stream);
    while (u <= 0.0);
    return u;
}

double st_rand_uniform_pos(void)
{
    int stream = 0;
#ifdef ST_RNG_MT
    stream = omp_get_thread_num();
#endif
    return st_rand_uniform_pos_stream(stream);
}

void st_rng_cleanup(void)
{
    int i;
#ifdef ST_RNG_MT
#pragma omp critical(st_rng_cleanup_lock)
#endif
    {
        for (i = 0; i < ST_RNG_MAXTHREADS; i++)
        {
            free(st_rng_state_ptr[i]);
            st_rng_state_ptr[i] = NULL;
        }
        st_rng_fallback_ready = 0;
        st_rng_base_ready = 0;
        st_rng_base_seed = 0;
    }
}

void st_rng_set_seed(uint32_t seed)
{
    st_rng_cleanup();
    st_rng_base_seed = seed & ST_RNG_MASK;
    st_rng_base_ready = 1;
}

typedef struct
{
    double time;
    int index;
} st_lbfgsb_breakpoint;

typedef struct
{
    int n, capacity, count;
    double theta;
    double *s, *y, *w, *sy;
} st_lbfgsb_memory;

static double st_lbfgsb_dot(const double *a, const double *b, int n)
{
    double answer = 0;
    for (int i = 0; i < n; i++)
        answer += a[i] * b[i];
    return answer;
}

static int st_lbfgsb_has_lower(int bound)
{
    return bound == 1 || bound == 2;
}

static int st_lbfgsb_has_upper(int bound)
{
    return bound == 2 || bound == 3;
}

static double st_lbfgsb_project(double value, double lower, double upper,
                                int bound)
{
    if (st_lbfgsb_has_lower(bound) && value < lower)
        value = lower;
    if (st_lbfgsb_has_upper(bound) && value > upper)
        value = upper;
    return value;
}

static double st_lbfgsb_projected(int n, const double *x, const double *g,
                                  const double *lower, const double *upper,
                                  const int *bounds)
{
    double largest = 0;
    for (int i = 0; i < n; i++)
    {
        double value = g[i];
        if (value > 0 && st_lbfgsb_has_lower(bounds[i]))
            value = fmin(value, x[i] - lower[i]);
        else if (value < 0 && st_lbfgsb_has_upper(bounds[i]))
            value = fmax(value, x[i] - upper[i]);
        largest = fmax(largest, fabs(value));
    }
    return largest;
}

static void st_lbfgsb_product(const st_lbfgsb_memory *memory,
                              const double *input, double *output)
{
    int n = memory->n;
    for (int i = 0; i < n; i++)
        output[i] = memory->theta * input[i];
    for (int j = 0; j < 2 * memory->count; j++)
    {
        const double *column = memory->w + (size_t)n * j;
        double coefficient = st_lbfgsb_dot(column, input, n);
        if (j & 1)
            coefficient = -coefficient;
        for (int i = 0; i < n; i++)
            output[i] += coefficient * column[i];
    }
}

static int st_lbfgsb_build(st_lbfgsb_memory *memory)
{
    int n = memory->n, count = memory->count;
    if (!count)
    {
        memory->theta = 1;
        return 1;
    }
    const double *last = memory->y + (size_t)n * (count - 1);
    memory->theta = st_lbfgsb_dot(last, last, n) / memory->sy[count - 1];
    if (!(memory->theta > 0) || !isfinite(memory->theta))
        return 0;
    for (int j = 0; j < count; j++)
    {
        const double *step = memory->s + (size_t)n * j;
        const double *difference = memory->y + (size_t)n * j;
        double *positive = memory->w + (size_t)n * (2 * j);
        double *negative = positive + n;
        for (int i = 0; i < n; i++)
            negative[i] = memory->theta * step[i];
        for (int k = 0; k < 2 * j; k++)
        {
            const double *column = memory->w + (size_t)n * k;
            double coefficient = st_lbfgsb_dot(column, step, n);
            if (k & 1)
                coefficient = -coefficient;
            for (int i = 0; i < n; i++)
                negative[i] += coefficient * column[i];
        }
        double curvature = st_lbfgsb_dot(step, negative, n);
        if (!(curvature > 0) || !isfinite(curvature))
            return 0;
        double scale_positive = 1 / sqrt(memory->sy[j]);
        double scale_negative = 1 / sqrt(curvature);
        for (int i = 0; i < n; i++)
        {
            positive[i] = difference[i] * scale_positive;
            negative[i] *= scale_negative;
            if (!isfinite(positive[i]) || !isfinite(negative[i]))
                return 0;
        }
    }
    return 1;
}

static int st_lbfgsb_break_compare(const void *left, const void *right)
{
    const st_lbfgsb_breakpoint *a = left, *b = right;
    if (a->time < b->time)
        return -1;
    if (a->time > b->time)
        return 1;
    return (a->index > b->index) - (a->index < b->index);
}

static int st_lbfgsb_cauchy(const st_lbfgsb_memory *memory, const double *x,
                            const double *g, const double *lower,
                            const double *upper, const int *bounds,
                            double *point, int *free_variable,
                            st_lbfgsb_breakpoint *breakpoints, double *wd,
                            double *wz)
{
    int n = memory->n, rank = 2 * memory->count, count = 0;
    double norm2 = 0;
    for (int j = 0; j < rank; j++)
        wd[j] = wz[j] = 0;
    for (int i = 0; i < n; i++)
    {
        int low = st_lbfgsb_has_lower(bounds[i]);
        int high = st_lbfgsb_has_upper(bounds[i]);
        int active = (low && high && lower[i] == upper[i]) ||
                     (low && x[i] <= lower[i] && g[i] >= 0) ||
                     (high && x[i] >= upper[i] && g[i] <= 0);
        free_variable[i] = !active;
        point[i] = x[i];
        if (active || g[i] == 0)
            continue;
        double direction = -g[i];
        norm2 += direction * direction;
        for (int j = 0; j < rank; j++)
            wd[j] += memory->w[i + (size_t)n * j] * direction;
        double breakpoint = INFINITY;
        if (direction > 0 && high)
            breakpoint = (upper[i] - x[i]) / direction;
        else if (direction < 0 && low)
            breakpoint = (lower[i] - x[i]) / direction;
        if (isfinite(breakpoint))
        {
            breakpoints[count].time = fmax(0, breakpoint);
            breakpoints[count++].index = i;
        }
    }
    if (!isfinite(norm2))
        return 0;
    qsort(breakpoints, count, sizeof(*breakpoints), st_lbfgsb_break_compare);
    double time = 0;
    int next = 0;
    while (norm2 > 0)
    {
        double first = (memory->theta * time - 1) * norm2;
        double second = memory->theta * norm2;
        for (int j = 0; j < rank; j++)
        {
            double sign = (j & 1) ? -1 : 1;
            first += sign * wd[j] * wz[j];
            second += sign * wd[j] * wd[j];
        }
        if (!isfinite(first) || !isfinite(second) || !(second > 0))
            return 0;
        if (first >= 0)
            break;
        double minimum = -first / second;
        double edge = next < count ? breakpoints[next].time : INFINITY;
        double distance = edge - time;
        if (minimum < distance || next == count)
        {
            time += minimum;
            break;
        }
        for (int j = 0; j < rank; j++)
            wz[j] += distance * wd[j];
        time = edge;
        do
        {
            int i = breakpoints[next++].index;
            double direction = -g[i];
            point[i] = direction > 0 ? upper[i] : lower[i];
            free_variable[i] = 0;
            norm2 = fmax(0, norm2 - direction * direction);
            for (int j = 0; j < rank; j++)
                wd[j] -= memory->w[i + (size_t)n * j] * direction;
        } while (next < count && breakpoints[next].time == edge);
    }
    for (int i = 0; i < n; i++)
        if (free_variable[i])
        {
            point[i] = st_lbfgsb_project(x[i] - time * g[i], lower[i], upper[i],
                                         bounds[i]);
            if ((st_lbfgsb_has_lower(bounds[i]) && point[i] == lower[i]) ||
                (st_lbfgsb_has_upper(bounds[i]) && point[i] == upper[i]))
                free_variable[i] = 0;
            if (!isfinite(point[i]))
                return 0;
        }
    return 1;
}

static int st_lbfgsb_small_solve(double *matrix, double *rhs, int n)
{
    for (int j = 0; j < n; j++)
    {
        int pivot = j;
        for (int i = j + 1; i < n; i++)
            if (fabs(matrix[i + (size_t)n * j]) >
                fabs(matrix[pivot + (size_t)n * j]))
                pivot = i;
        double diagonal = matrix[pivot + (size_t)n * j];
        if (!isfinite(diagonal) || fabs(diagonal) < DBL_MIN)
            return 0;
        if (pivot != j)
        {
            for (int k = j; k < n; k++)
            {
                double value = matrix[j + (size_t)n * k];
                matrix[j + (size_t)n * k] = matrix[pivot + (size_t)n * k];
                matrix[pivot + (size_t)n * k] = value;
            }
            double value = rhs[j];
            rhs[j] = rhs[pivot];
            rhs[pivot] = value;
        }
        for (int i = j + 1; i < n; i++)
        {
            double factor = matrix[i + (size_t)n * j] / diagonal;
            for (int k = j + 1; k < n; k++)
                matrix[i + (size_t)n * k] -= factor * matrix[j + (size_t)n * k];
            rhs[i] -= factor * rhs[j];
        }
    }
    for (int j = n - 1; j >= 0; j--)
    {
        double value = rhs[j];
        for (int k = j + 1; k < n; k++)
            value -= matrix[j + (size_t)n * k] * rhs[k];
        rhs[j] = value / matrix[j + (size_t)n * j];
        if (!isfinite(rhs[j]))
            return 0;
    }
    return 1;
}

static int st_lbfgsb_subspace(const st_lbfgsb_memory *memory, const double *x,
                              const double *g, const double *lower,
                              const double *upper, const int *bounds,
                              const int *free_variable, double *point,
                              double *step, double *residual, double *work,
                              double *matrix, double *rhs)
{
    int n = memory->n, rank = 2 * memory->count, free_count = 0;
    for (int i = 0; i < n; i++)
    {
        free_count += free_variable[i];
        step[i] = point[i] - x[i];
    }
    if (!free_count)
        return 1;
    st_lbfgsb_product(memory, step, residual);
    for (int i = 0; i < n; i++)
        residual[i] += g[i];
    if (free_count == n)
    {
        memcpy(work, g, (size_t)n * sizeof(double));
        for (int j = memory->count - 1; j >= 0; j--)
        {
            const double *s = memory->s + (size_t)n * j;
            const double *y = memory->y + (size_t)n * j;
            rhs[j] = st_lbfgsb_dot(s, work, n) / memory->sy[j];
            for (int i = 0; i < n; i++)
                work[i] -= rhs[j] * y[i];
        }
        for (int i = 0; i < n; i++)
            work[i] /= memory->theta;
        for (int j = 0; j < memory->count; j++)
        {
            const double *s = memory->s + (size_t)n * j;
            const double *y = memory->y + (size_t)n * j;
            double coefficient =
                rhs[j] - st_lbfgsb_dot(y, work, n) / memory->sy[j];
            for (int i = 0; i < n; i++)
                work[i] += coefficient * s[i];
        }
        for (int i = 0; i < n; i++)
            step[i] = x[i] - point[i] - work[i];
    }
    else
    {
        for (int j = 0; j < rank; j++)
        {
            const double *column = memory->w + (size_t)n * j;
            rhs[j] = 0;
            for (int i = 0; i < n; i++)
                if (free_variable[i])
                    rhs[j] += column[i] * residual[i] / memory->theta;
            for (int k = 0; k <= j; k++)
            {
                const double *other = memory->w + (size_t)n * k;
                double value = 0;
                for (int i = 0; i < n; i++)
                    if (free_variable[i])
                        value += column[i] * other[i] / memory->theta;
                if (k == j)
                    value += (j & 1) ? -1 : 1;
                matrix[k + (size_t)rank * j] = value;
                matrix[j + (size_t)rank * k] = value;
            }
        }
        if (!st_lbfgsb_small_solve(matrix, rhs, rank))
            return 0;
        for (int i = 0; i < n; i++)
        {
            step[i] = 0;
            if (!free_variable[i])
                continue;
            double value = -residual[i];
            for (int j = 0; j < rank; j++)
                value += memory->w[i + (size_t)n * j] * rhs[j];
            step[i] = value / memory->theta;
            if (!isfinite(step[i]))
                return 0;
        }
    }
    double fraction = 1;
    int hit = -1;
    for (int i = 0; i < n; i++)
    {
        double bound_fraction = INFINITY;
        if (step[i] < 0 && st_lbfgsb_has_lower(bounds[i]))
            bound_fraction = (lower[i] - point[i]) / step[i];
        else if (step[i] > 0 && st_lbfgsb_has_upper(bounds[i]))
            bound_fraction = (upper[i] - point[i]) / step[i];
        if (bound_fraction < fraction)
        {
            fraction = fmax(0, bound_fraction);
            hit = i;
        }
    }
    for (int i = 0; i < n; i++)
        point[i] = st_lbfgsb_project(point[i] + fraction * step[i], lower[i],
                                     upper[i], bounds[i]);
    if (hit >= 0)
        point[hit] = step[hit] > 0 ? upper[hit] : lower[hit];
    return 1;
}

typedef struct
{
    double step, value, slope;
} st_lbfgsb_line_point;

static double st_lbfgsb_cubic(st_lbfgsb_line_point left,
                              st_lbfgsb_line_point right, double *leading)
{
    double width = right.step - left.step;
    double linear = width * left.slope;
    double quadratic =
        3 * (right.value - left.value) - width * (2 * left.slope + right.slope);
    double cubic =
        2 * (left.value - right.value) + width * (left.slope + right.slope);
    *leading = cubic;
    double scale = fmax(fabs(linear), fmax(fabs(quadratic), fabs(cubic)));
    if (!(scale > 0) || !isfinite(scale) || width == 0)
        return NAN;
    double a = 3 * cubic / scale, b = 2 * quadratic / scale, c = linear / scale;
    if (fabs(a) <= DBL_EPSILON * fabs(b))
        return b > 0 ? left.step - width * c / b : NAN;
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0 || !isfinite(discriminant))
        return NAN;
    double q = -.5 * (b + copysign(sqrt(discriminant), b));
    double roots[2] = {q / a, q != 0 ? c / q : NAN};
    for (int j = 0; j < 2; j++)
        if (isfinite(roots[j]) && b + 2 * a * roots[j] > 0)
            return left.step + width * roots[j];
    return NAN;
}

static double st_lbfgsb_trial(st_lbfgsb_line_point best,
                              st_lbfgsb_line_point trial,
                              st_lbfgsb_line_point other, int bracketed)
{
    double leading, cubic = st_lbfgsb_cubic(best, trial, &leading);
    double width = trial.step - best.step;
    double secant =
        trial.step - trial.slope * width / (trial.slope - best.slope);
    double answer;
    if (trial.value > best.value)
    {
        double denominator =
            2 * (trial.value - best.value - width * best.slope);
        double quadratic = best.step - best.slope * width * width / denominator;
        if (!isfinite(cubic))
            cubic = .5 * (best.step + trial.step);
        if (!isfinite(quadratic))
            quadratic = .5 * (best.step + trial.step);
        answer = fabs(cubic - best.step) < fabs(quadratic - best.step)
                     ? cubic
                     : .5 * (cubic + quadratic);
    }
    else if (signbit(trial.slope) != signbit(best.slope))
    {
        if (!isfinite(cubic))
            cubic = secant;
        answer = fabs(cubic - trial.step) >= fabs(secant - trial.step) ? cubic
                                                                       : secant;
    }
    else if (fabs(trial.slope) < fabs(best.slope))
    {
        int beyond =
            isfinite(cubic) && leading > 0 && (cubic - trial.step) * width > 0;
        answer = beyond && fabs(cubic - trial.step) < fabs(secant - trial.step)
                     ? cubic
                     : secant;
        if (bracketed)
        {
            double safeguard = trial.step + .66 * (other.step - trial.step);
            answer =
                width > 0 ? fmin(answer, safeguard) : fmax(answer, safeguard);
        }
    }
    else
        answer = bracketed ? st_lbfgsb_cubic(trial, other, &leading)
                           : trial.step + 4 * width;
    return answer;
}

static int st_lbfgsb_evaluate(int n, double *point, double *g, double *value,
                              st_optim_function *function,
                              st_optim_gradient *gradient, void *data,
                              st_lbfgsb_result *result)
{
    *value = function(n, point, data);
    result->function_calls++;
    if (!isfinite(*value))
    {
        result->status = -1;
        result->message = "L-BFGS-B needs finite values of 'fn'";
        return 0;
    }
    gradient(n, point, g, data);
    result->gradient_calls++;
    for (int i = 0; i < n; i++)
        if (!isfinite(g[i]))
        {
            result->status = -1;
            result->message = "L-BFGS-B needs finite gradient values";
            return 0;
        }
    return 1;
}

static int st_lbfgsb_line(int n, const double *x, const double *g,
                          const double *direction, const double *lower,
                          const double *upper, const int *bounds, double start,
                          double maximum, double *point, double *new_gradient,
                          double *best_gradient, double *new_value,
                          st_optim_function *function,
                          st_optim_gradient *gradient, void *data,
                          st_lbfgsb_result *result)
{
    const double sufficient = 1e-3, curvature = .9;
    double initial = result->value;
    double initial_slope = st_lbfgsb_dot(g, direction, n);
    double shift = sufficient * initial_slope;
    st_lbfgsb_line_point best = {0, initial, initial_slope};
    st_lbfgsb_line_point other = best;
    double step = fmin(start, maximum);
    double previous_width = maximum, older_width = 2 * maximum;
    int bracketed = 0, auxiliary = 1;
    if (!(initial_slope < 0) || !(step > 0) || !isfinite(step))
        return 0;
    memcpy(best_gradient, g, (size_t)n * sizeof(double));
    for (int evaluation = 0; evaluation < 20; evaluation++)
    {
        for (int i = 0; i < n; i++)
            point[i] = st_lbfgsb_project(x[i] + step * direction[i], lower[i],
                                         upper[i], bounds[i]);
        if (!st_lbfgsb_evaluate(n, point, new_gradient, new_value, function,
                                gradient, data, result))
            return -1;
        st_lbfgsb_line_point trial = {
            step, *new_value, st_lbfgsb_dot(new_gradient, direction, n)};
        int decrease = trial.value <= initial + shift * step;
        if (decrease && (fabs(trial.slope) <= -curvature * initial_slope ||
                         (step == maximum && trial.slope < 0)))
            return 1;
        if (auxiliary && decrease && trial.slope >= 0)
            auxiliary = 0;
        st_lbfgsb_line_point a = best, b = trial, c = other;
        if (auxiliary)
        {
            a.value -= shift * a.step;
            b.value -= shift * b.step;
            c.value -= shift * c.step;
            a.slope -= shift;
            b.slope -= shift;
            c.slope -= shift;
        }
        double proposal = st_lbfgsb_trial(a, b, c, bracketed);
        double movement = trial.step - best.step;
        if (b.value > a.value)
        {
            other = trial;
            bracketed = 1;
        }
        else
        {
            if (b.slope * (a.step - b.step) < 0)
            {
                other = best;
                bracketed = 1;
            }
            best = trial;
            memcpy(best_gradient, new_gradient, (size_t)n * sizeof(double));
        }
        if (bracketed)
        {
            double low = fmin(best.step, other.step);
            double high = fmax(best.step, other.step);
            double width = high - low;

            if (best.step > 0 && best.value <= initial + shift * best.step &&
                width <= .1 * high)
            {
                for (int i = 0; i < n; i++)
                    point[i] =
                        st_lbfgsb_project(x[i] + best.step * direction[i],
                                          lower[i], upper[i], bounds[i]);
                memcpy(new_gradient, best_gradient, (size_t)n * sizeof(double));
                *new_value = best.value;
                return 1;
            }
            if (!(proposal > low && proposal < high) ||
                width >= .66 * older_width)
                proposal = .5 * (best.step + other.step);
            older_width = previous_width;
            previous_width = width;
            if (width <= 4 * DBL_EPSILON * fmax(fabs(low), fabs(high)))
                return 0;
        }
        else
        {
            double from = trial.step + 1.1 * movement;
            double to = trial.step + 4 * movement;
            if (!isfinite(proposal))
                proposal = to;
            proposal = fmax(fmin(from, to), fmin(fmax(from, to), proposal));
        }
        proposal = fmax(0, fmin(maximum, proposal));
        if (!(proposal > 0) || proposal == step || !isfinite(proposal))
            return 0;
        step = proposal;
    }
    return 0;
}

void st_lbfgsb_minimize(int n, int history, double *x, const double *lower,
                        const double *upper, const int *bounds,
                        st_optim_function *function,
                        st_optim_gradient *gradient, void *data, int maxit,
                        double factr, double pgtol, st_optim_monitor *monitor,
                        void *monitor_data, st_lbfgsb_result *result)
{
    memset(result, 0, sizeof(*result));
    result->status = 52;
    result->message = "ERROR: INVALID L-BFGS-B CONTROLS";
    result->value = INFINITY;
    result->projected_gradient = INFINITY;
    if (n < 0 || history < 1 || history > (INT_MAX - 7) / 4 || maxit < 0 ||
        !isfinite(factr) || factr < 0 || !isfinite(pgtol) || pgtol < 0)
        return;
    for (int i = 0; i < n; i++)
    {
        if (bounds[i] < 0 || bounds[i] > 3 || !isfinite(x[i]) ||
            (st_lbfgsb_has_lower(bounds[i]) && isnan(lower[i])) ||
            (st_lbfgsb_has_upper(bounds[i]) && isnan(upper[i])) ||
            (bounds[i] == 2 && lower[i] > upper[i]))
        {
            result->message = "ERROR: INVALID L-BFGS-B BOUNDS OR START";
            return;
        }
        x[i] = st_lbfgsb_project(x[i], lower[i], upper[i], bounds[i]);
        if (!isfinite(x[i]))
        {
            result->message = "ERROR: INVALID L-BFGS-B BOUNDS OR START";
            return;
        }
    }
    size_t rank = (size_t)2 * history;
    if (rank > SIZE_MAX / rank || (size_t)n > SIZE_MAX / (2 * rank + 7) ||
        rank * rank > SIZE_MAX - (size_t)n * (2 * rank + 7) ||
        rank * rank + (size_t)n * (2 * rank + 7) >
            SIZE_MAX - (size_t)7 * history)
    {
        result->message = "ERROR: L-BFGS-B WORKSPACE SIZE OVERFLOW";
        return;
    }
    size_t size =
        rank * rank + (size_t)n * (2 * rank + 7) + (size_t)7 * history;
    double *workspace = st_work_alloc(size ? size : 1, sizeof(double));
    int *free_variable = st_work_alloc(n ? n : 1, sizeof(int));
    st_lbfgsb_breakpoint *breakpoints =
        st_work_alloc(n ? n : 1, sizeof(*breakpoints));
    double *next = workspace;
    st_lbfgsb_memory memory = {n, history, 0, 1, NULL, NULL, NULL, NULL};
    memory.s = next;
    next += (size_t)n * history;
    memory.y = next;
    next += (size_t)n * history;
    memory.w = next;
    next += (size_t)n * rank;
    memory.sy = next;
    next += history;
    double *g = next;
    next += n;
    double *point = next;
    next += n;
    double *new_gradient = next;
    next += n;
    double *direction = next;
    next += n;
    double *step = next;
    next += n;
    double *residual = next;
    next += n;
    double *work = next;
    next += n;
    double *wd = next;
    next += rank;
    double *wz = next;
    next += rank;
    double *rhs = next;
    next += rank;
    double *matrix = next;
    result->status = 0;
    result->message = "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL";
    if (!st_lbfgsb_evaluate(n, x, g, &result->value, function, gradient, data,
                            result))
        goto finish;
    result->projected_gradient =
        st_lbfgsb_projected(n, x, g, lower, upper, bounds);
    if (monitor)
        monitor(0, result->function_calls, result->value,
                result->projected_gradient, x, n, monitor_data);
    while (result->projected_gradient > pgtol)
    {
        if (result->iterations > maxit || result->iterations == INT_MAX)
        {
            result->status = 1;
            result->message = "STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT";
            break;
        }
        int accepted = 0;
        double new_value = result->value;
        for (int attempt = 0; attempt < 2 && !accepted; attempt++)
        {
            if (!st_lbfgsb_build(&memory) ||
                !st_lbfgsb_cauchy(&memory, x, g, lower, upper, bounds, point,
                                  free_variable, breakpoints, wd, wz) ||
                !st_lbfgsb_subspace(&memory, x, g, lower, upper, bounds,
                                    free_variable, point, step, residual, work,
                                    matrix, rhs))
            {
                if (!memory.count)
                    break;
                memory.count = 0;
                continue;
            }
            double norm2 = 0, maximum = 1e10;
            int bounded = 0, boxed = 1;
            for (int i = 0; i < n; i++)
            {
                direction[i] = point[i] - x[i];
                norm2 += direction[i] * direction[i];
                bounded |= bounds[i] != 0;
                boxed &= bounds[i] == 2;
                if (direction[i] < 0 && st_lbfgsb_has_lower(bounds[i]))
                    maximum = fmin(maximum, (lower[i] - x[i]) / direction[i]);
                else if (direction[i] > 0 && st_lbfgsb_has_upper(bounds[i]))
                    maximum = fmin(maximum, (upper[i] - x[i]) / direction[i]);
            }
            if (!result->iterations && bounded)
                maximum = fmin(1, maximum);
            double start = 1;
            if (!result->iterations && !boxed && norm2 > 0)
                start = fmin(1, 1 / sqrt(norm2));
            int status =
                st_lbfgsb_line(n, x, g, direction, lower, upper, bounds, start,
                               maximum, point, new_gradient, work, &new_value,
                               function, gradient, data, result);
            if (status < 0)
                goto finish;
            if (status > 0)
                accepted = 1;
            else if (memory.count)
                memory.count = 0;
            else
                break;
        }
        if (!accepted)
        {
            result->status = 52;
            result->message = "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH";
            break;
        }
        for (int i = 0; i < n; i++)
        {
            step[i] = point[i] - x[i];
            residual[i] = new_gradient[i] - g[i];
        }
        double sy = st_lbfgsb_dot(step, residual, n);
        double yy = st_lbfgsb_dot(residual, residual, n);
        double slope = -st_lbfgsb_dot(g, step, n);
        if (isfinite(sy) && isfinite(yy) && yy > 0 &&
            sy > DBL_EPSILON * fmax(0, slope))
        {
            if (memory.count == history)
            {
                memmove(memory.s, memory.s + n,
                        (size_t)n * (history - 1) * sizeof(double));
                memmove(memory.y, memory.y + n,
                        (size_t)n * (history - 1) * sizeof(double));
                memmove(memory.sy, memory.sy + 1,
                        (size_t)(history - 1) * sizeof(double));
                memory.count--;
            }
            int j = memory.count++;
            memcpy(memory.s + (size_t)n * j, step, (size_t)n * sizeof(double));
            memcpy(memory.y + (size_t)n * j, residual,
                   (size_t)n * sizeof(double));
            memory.sy[j] = sy;
        }
        double previous = result->value;
        result->value = new_value;
        memcpy(x, point, (size_t)n * sizeof(double));
        memcpy(g, new_gradient, (size_t)n * sizeof(double));
        result->iterations++;
        result->projected_gradient =
            st_lbfgsb_projected(n, x, g, lower, upper, bounds);
        if (monitor)
            monitor(result->iterations, result->function_calls, result->value,
                    result->projected_gradient, x, n, monitor_data);
        if (result->projected_gradient <= pgtol)
            break;
        if (previous - result->value <=
            factr * DBL_EPSILON *
                fmax(1, fmax(fabs(previous), fabs(result->value))))
        {
            result->message = "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH";
            break;
        }
    }
finish:
    st_work_detach(breakpoints);
    free(breakpoints);
    st_work_detach(free_variable);
    free(free_variable);
    st_work_detach(workspace);
    free(workspace);
}
