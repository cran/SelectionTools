

#define AUTHOR "$Author: frisch-m $"
#define DATE "$Date: 2019/06/10 00:00:00 $"
#define REVISION "$Revision: 19.2 $"

#include <stdio.h>
#include "SelectionToolsAllele.h"
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include <stdint.h>
#include <limits.h>
#include <float.h>
#include <wchar.h>
#include <wctype.h>

#include <stdint.h>

#ifdef _OPENMP
#ifndef ST
#define MT
#include <omp.h>
#endif
#endif

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Random.h>

#define MAXTHREADS 256

#include "MathLib12.h"

#if MAXTHREADS != ST_RNG_MAXTHREADS
#error "MAXTHREADS and ST_RNG_MAXTHREADS differ"
#endif

#ifndef MT
#define omp_get_thread_num() 0
#define omp_set_num_threads(thr)                                               \
    {                                                                          \
    }
#define omp_get_num_threads() 0
#define omp_set_dynamic(thr)                                                   \
    {                                                                          \
    }
#endif

void info(int level, char *info, ...);
void remove_genotype_population(char **PopNames);
void remove_evaluate_population(char **PopNames);
void flush_population(char **PopName, int mode);
void effmap_remove_all(void);

int rng_initialized = 0;
static int rng_seeded_mode = 0;
#ifdef MT
static int rng_previous_dynamic = 0;
#endif

int st_rng_is_seeded(void)
{
    return rng_seeded_mode;
}

#ifdef use_gsl

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *rng[MAXTHREADS];

const gsl_rng_type **rng_default;
const gsl_rng_type **all_types;

char rng_def_nme[] = "taus";
char *pp_rng_def_nme = rng_def_nme;
char **p_rng_def_nme = &pp_rng_def_nme;

void rng_init(void);

void rng_info(void)
{
    info(0, "Random number generator type: '%s'\n", gsl_rng_name(rng[0]));
}

void rng_choose(char **name)
{
    for (rng_default = all_types; *rng_default != 0; rng_default++)
    {
        if (0 == strcmp((*rng_default)->name, *name))
            break;
    }
    if (*rng_default == 0)
        rng_choose(p_rng_def_nme);
    rng_init();
    rng_info();
}

void rng_list(void)
{
    const gsl_rng_type **t;

    info(0, "Available generators:\n");
    for (t = all_types; *t != 0; t++)
    {
        if (0 == strcmp(gsl_rng_name(rng[0]), (*t)->name))
        {
            info(0, "--> %s <--\n", (*t)->name);
        }
        else
        {
            info(0, "    %s\n", (*t)->name);
        }
    }
}

#endif

#ifndef use_gsl

#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

typedef struct
{
    long seed;
    long iyy;
    long ivv[NTAB];
} gsl_rng;

gsl_rng *rng[MAXTHREADS];

static inline long pm_step(long x)
{

    const long k = x / IQ;
    long xnew = IA * (x - k * IQ) - IR * k;
    if (xnew < 0)
        xnew += IM;
    return xnew;
}

double gsl_rng_uniform(gsl_rng *rng)
{

    int j;
    long *idum = &(rng->seed);
    long *iy = &(rng->iyy);
    long *iv = &(rng->ivv[0]);

    if ((*idum) <= 0 || (*iy) == 0)
    {
        long x = *idum;

        if (x == 0)
            x = 1;
        if (x < 0)
            x = -x;

        for (j = NTAB + 7; j >= 0; --j)
        {
            x = pm_step(x);
            if (j < NTAB)
                iv[j] = x;
        }

        *iy = iv[0];
        *idum = x;
    }

    *idum = pm_step(*idum);

    j = (int)((*iy) / NDIV);
    if (j < 0)
        j = 0;
    if (j >= NTAB)
        j = NTAB - 1;

    *iy = iv[j];
    iv[j] = *idum;

    {
        const double u = AM * (double)(*iy);
        return (u > RNMX) ? RNMX : u;
    }
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

double gammln(double xx)
{

    return lgamma(xx);
}

#define rPI 3.141592654

double gsl_rng_uniform_pos(gsl_rng *rng);

unsigned int gsl_ran_poisson(gsl_rng *rng, double xm)
{

    if (!(xm >= 0.0) || !isfinite(xm))
        return 0U;
    if (xm == 0.0)
        return 0U;

    if (xm < 12.0)
    {
        const double L = exp(-xm);
        unsigned int k = 0U;
        double p = 1.0;

        do
        {
            ++k;
            p *= gsl_rng_uniform_pos(rng);
        } while (p > L);

        return k - 1U;
    }

    {
        const double mu = xm;
        const double s = sqrt(mu);
        const double b = 0.931 + 2.53 * s;
        const double a = -0.059 + 0.02483 * b;
        const double inv_alpha = 1.1239 + 1.1328 / (b - 3.4);
        const double v_r = 0.9277 - 3.6224 / (b - 2.0);

        for (;;)
        {
            const double U = gsl_rng_uniform_pos(rng) - 0.5;
            const double V = gsl_rng_uniform_pos(rng);

            const double us = 0.5 - fabs(U);
            double k = floor((2.0 * a / us + b) * U + mu + 0.43);

            if (us >= 0.07 && V <= v_r)
                return (unsigned int)k;

            if (k < 0.0)
                continue;

            {
                const double logV = log(V * inv_alpha / (a / (us * us) + b));
                const double logP = k * log(mu) - gammln(k + 1.0) - mu;
                if (logV <= logP)
                    return (unsigned int)k;
            }
        }
    }
}

#undef rPI

double gsl_rng_uniform_pos(gsl_rng *rng)
{

    double r;
    do
    {
        r = gsl_rng_uniform(rng);
    } while (r <= 0.0 || r >= 1.0);
    return r;
}

unsigned long int gsl_rng_uniform_int(gsl_rng *rng, unsigned long int n)
{
    const unsigned long int range = 2147483646UL;
    unsigned long int x, limit, q, r, h;

    if (n <= 1UL)
        return 0UL;

    if (n <= range)
    {
        q = range / n;
        limit = q * n;

        do
        {
            (void)gsl_rng_uniform(rng);
            x = (unsigned long int)(rng->iyy - 1L);
        } while (x >= limit);

        return x / q;
    }

    q = n / range;
    r = n % range;

    if (r == 0UL)
    {
        h = gsl_rng_uniform_int(rng, q);
        x = gsl_rng_uniform_int(rng, range);
        return h * range + x;
    }

    do
    {
        h = gsl_rng_uniform_int(rng, q + 1UL);
        x = gsl_rng_uniform_int(rng, range);
    } while ((h == q) && (x >= r));

    return h * range + x;
}

gsl_rng *gsl_rng_alloc(void)
{

    gsl_rng *r = (gsl_rng *)malloc(sizeof(gsl_rng));
    if (!r)
        return NULL;
    r->seed = -1;
    r->iyy = 0;
    for (int i = 0; i < (int)(sizeof(r->ivv) / sizeof(r->ivv[0])); ++i)
        r->ivv[i] = 0;
    return r;
}

void gsl_rng_set(gsl_rng *rng, unsigned long int s)
{

    if (s == 0UL)
        s = 1UL;
    rng->seed = -(long)s;
    rng->iyy = 0;
}

void gsl_rng_free(gsl_rng *r)
{
    free(r);
}

void rng_choose(char **name_rng)
{
    (void)name_rng;
    info(0,
         "Random number generator: Park-Miller (clean-room) + shuffle table\n");
}

void rng_info(void)
{
    info(0, "Random number generator type: Park-Miller (clean-room) + shuffle "
            "table\n");
}

void rng_list(void)
{
    info(0, "Available generators: Park-Miller (clean-room) + shuffle table\n");
}

#define p_rng_default

#endif

void rng_init(void)
{

    if (rng_seeded_mode)
    {
        st_rng_cleanup();
        rng_seeded_mode = 0;
#ifdef MT
        omp_set_dynamic(rng_previous_dynamic);
#endif
    }

    int i, j;
    unsigned long int seed;
    gsl_rng *tmp_rng[MAXTHREADS] = {NULL};

#ifdef use_gsl
    if (0 == rng_initialized)
    {
        all_types = gsl_rng_types_setup();

        for (rng_default = all_types; *rng_default != 0; rng_default++)
        {
            if (0 == strcmp((*rng_default)->name, rng_def_nme))
                break;
        }
    }
    const gsl_rng_type *p_rng_default = *rng_default;
#endif

    for (i = 0; i < MAXTHREADS; i++)
    {
        seed = (unsigned long)st_rand() + 1UL;
        tmp_rng[i] = gsl_rng_alloc(p_rng_default);

        if (NULL == tmp_rng[i])
            goto fail;

        gsl_rng_set(tmp_rng[i], seed);
    }

    if (1 == rng_initialized)
    {
        for (i = 0; i < MAXTHREADS; i++)
            gsl_rng_free(rng[i]);
        info(1, "Re-initializing random number generator \n");
    }
    else
    {
        info(1, "Initializing random number generator \n");
    }

    for (i = 0; i < MAXTHREADS; i++)
        rng[i] = tmp_rng[i];

    rng_initialized = 1;
    return;

fail:

    for (j = 0; j < i; j++)
        gsl_rng_free(tmp_rng[j]);

    info(-2, "Memory allocation failed in rng_init\n");
    return;
}

void rng_init_seed(int *seeds, int *nseeds, int *retval)
{
    gsl_rng *tmp_rng[MAXTHREADS] = {NULL};
    int i, j;
    *retval = -2;
    if (*nseeds != MAXTHREADS + 1)
        return;
    if (seeds[0] < 0 || seeds[0] > 1073741823)
        return;
    for (i = 1; i <= MAXTHREADS; i++)
        if (seeds[i] < 1 || seeds[i] > 1073741824)
            return;
#ifdef use_gsl
    const gsl_rng_type *p_rng_default = *rng_default;
#endif
    for (i = 0; i < MAXTHREADS; i++)
    {
        tmp_rng[i] = gsl_rng_alloc(p_rng_default);
        if (NULL == tmp_rng[i])
        {
            for (j = 0; j < i; j++)
                gsl_rng_free(tmp_rng[j]);
            return;
        }
        gsl_rng_set(tmp_rng[i], (unsigned long)seeds[i + 1]);
    }
    st_rng_set_seed((uint32_t)seeds[0]);
    for (i = 0; i < MAXTHREADS; i++)
    {
        if (rng_initialized)
            gsl_rng_free(rng[i]);
        rng[i] = tmp_rng[i];
    }
#ifdef MT
    if (!rng_seeded_mode)
        rng_previous_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
#endif
    rng_seeded_mode = 1;
    rng_initialized = 1;
    *retval = 0;
}

#define MAXNAME 256
#define MAXLINE 100
#define CPOS double
#define POPULATION long int
#define MAP_DISPERSE_FACTOR_DEFAULT 100.0

#define CCTIME __TIME__
#define CCDATE __DATE__

#define yCalloc(n, t) (t *)malloc((size_t)(n) * sizeof(t))

static void *yRealloc_impl(void *ptr, size_t count, size_t size)
{
    void *new_ptr;

    if ((count == 0) || (size == 0))
        goto fail;
    if (count > SIZE_MAX / size)
        goto fail;

    new_ptr = realloc(ptr, count * size);

    if (new_ptr == NULL)
        goto fail;

    return new_ptr;

fail:
    free(ptr);
    return NULL;
}

#define yRealloc(p, n, t)                                                      \
    (t *)yRealloc_impl((void *)(p), (size_t)(n), sizeof(t))

#define yFree(p) (free((void *)(p)), (p) = NULL)

typedef struct
{
    CPOS Pos;
    int chrom;
    char *name;
    char *class;
    int counter;
} TMappoint;

typedef struct
{
    double ChromLen;
    int NrMappoints;
    TMappoint *Mappoint;
} TChromDescr;

typedef struct
{
    CPOS Pos;
    unsigned short All;
} TLoc;

typedef struct st_sim_storage st_sim_storage;

typedef struct
{
    int NoLociAlloc;
    int NoLociUsed;
    TLoc *Loc;
    st_sim_storage *Owner;
    size_t Offset;
} THom;

typedef struct
{
    THom *Hom;
} TChrom;

struct st_sim_storage
{
    size_t NoChroms;
    size_t NoHoms;
    size_t NoSlots;
    size_t NoLociAlloc;
    size_t NoLociReserved;
    TChrom *Chrom;
    THom *Hom;
    TLoc *Loc;
};

typedef struct
{
    TChrom *Chrom;
    st_sim_storage *Storage;
    float *GValue;
    float *PValue;
    int *genotype;
    char *info;
    char *name;
    char *pedigree;
} TInd;

typedef struct
{
    char Name[MAXNAME];
    POPULATION NoIndsAlloc;
    POPULATION NoInds;
    TInd *Ind;
} TPop;

typedef struct
{
    int sizeRows, sizeCols;
    int nPop, nMarker;
    int *nAll, *nInd;
    int *nAll_off, *nInd_off;
    int *data;
} TPopMatrix;

struct TIntLinkList
{
    int number;
    struct TIntLinkList *child, *parent;
    struct TIntLinkList *next, *prev;
};

typedef struct TIntLinkList IntNode;

typedef struct SPopListEl
{
    TPop *Pop;
    struct SPopListEl *Prev;
    struct SPopListEl *Next;
} TPopListEl;

typedef struct
{
    TMappoint *MP;
    int NrAll;
    int *All;
} TGen;

typedef struct
{
    float Eff;
    int NrLoci;
    TGen *Loc;
} TEff;

typedef struct
{
    float GValue_avg;
    float weight;
    int NrEff;
    int Pos;
    TEff *Eff;
    char Name[3 * MAXNAME];
    int MPs;
    TMappoint **MPset;
} TEffects;

typedef struct
{
    int betrag;
    int *member;
} TSet;

int NoChroms = 0;
int NoHoms = 0;
int NoLociInit = 10;
IntNode *startNode = NULL;
double COFreq = 1;

TPopListEl *PopListEl = NULL;
int bg, missing, doubleHetero;
int dimLA, dimLB, locusA, locusB, nrGeno;
int **table, *genotype;

TPopListEl *PopListStart = NULL;
TChromDescr *PopDescr = NULL;
int cMP = 0;
int cMPhom = 0;
int INFOLEVEL = 0;

TEffects *Effects = NULL;
int NrEffects = 0;
TSet *AllSet;
unsigned long int MDP = 0;

static uint64_t st_identity_serial = 0;
static const char **st_cross_left = NULL, **st_cross_right = NULL;
static uint64_t st_cross_first = 0;
static size_t st_cross_count = 0;
void ps_identity_observe(const char *name)
{
    if (!name)
        return;
    for (const char *dot = strchr(name, '.'); dot; dot = strchr(dot + 1, '.'))
    {
        const char *p = dot + 1;
        uint64_t value = 0;
        int digits = 0, overflow = 0;
        while (*p >= '0' && *p <= '9')
        {
            unsigned int digit = (unsigned int)(*p - '0');
            if (value > (UINT64_MAX - digit) / 10)
                overflow = 1;
            if (!overflow)
                value = value * 10 + digit;
            p++;
            digits = 1;
        }
        if (digits && !overflow && (*p == '(' || !*p) &&
            value > st_identity_serial)
            st_identity_serial = value;
    }
}
static int st_sim_reserve(POPULATION n, uint64_t *first)
{
    if (n < 0 || (uint64_t)n > UINT64_MAX - st_identity_serial)
        return -1;
    *first = st_identity_serial + 1;
    st_identity_serial += (uint64_t)n;
    return 0;
}
typedef struct st_sim_pool_block st_sim_pool_block;
typedef struct st_sim_pool_span st_sim_pool_span;

struct st_sim_pool_span
{
    st_sim_pool_block *block;
    st_sim_pool_span *prev;
    st_sim_pool_span *next;
    st_sim_pool_span *free_prev;
    st_sim_pool_span *free_next;
    size_t capacity;
    size_t requested;
    int used;
};

struct st_sim_pool_block
{
    st_sim_pool_block *next;
    st_sim_pool_span *first;
    size_t bytes;
};

typedef union
{
    long double ld;
    long long ll;
    void *pointer;
    size_t size;
    TLoc locus;
    THom homologue;
    TChrom chromosome;
    TInd individual;
    TPop population;
    TPopListEl list;
    st_sim_pool_span span;
    st_sim_pool_block block;
} st_sim_pool_alignment;

typedef struct
{
    char first;
    st_sim_pool_alignment aligned;
} st_sim_pool_alignment_probe;

typedef struct
{
    size_t allocation_requests;
    size_t resize_requests;
    size_t free_requests;
    size_t backing_requests;
    size_t backing_allocations;
    size_t backing_frees;
    size_t failed_requests;
    size_t reused_allocations;
    size_t inplace_resizes;
    size_t moving_resizes;
    size_t splits;
    size_t coalesces;
    size_t disposal_rejections;
    size_t live_allocations;
    size_t live_bytes;
    size_t live_capacity_bytes;
    size_t reserved_bytes;
    size_t peak_live_allocations;
    size_t peak_live_bytes;
    size_t peak_live_capacity_bytes;
    size_t peak_reserved_bytes;
} st_sim_pool_statistics;

#define ST_SIM_POOL_BINS (sizeof(size_t) * CHAR_BIT)
#define ST_SIM_POOL_ALIGN offsetof(st_sim_pool_alignment_probe, aligned)

size_t add_block_size = 100000000;

static st_sim_pool_block *st_sim_pool_blocks;
static st_sim_pool_span *st_sim_pool_bins[ST_SIM_POOL_BINS];
static st_sim_pool_statistics st_sim_pool_stats;

#ifdef ST_SIM_POOL_DIAGNOSTICS
static size_t st_sim_pool_fail_after = SIZE_MAX;
#define ST_SIM_POOL_COUNT(field)                                                \
    do                                                                        \
    {                                                                         \
        if (st_sim_pool_stats.field != SIZE_MAX)                               \
            st_sim_pool_stats.field++;                                         \
    } while (0)
#else
#define ST_SIM_POOL_COUNT(field) ((void)0)
#endif

static int st_sim_pool_round(size_t value, size_t quantum, size_t *rounded)
{
    size_t remainder;

    if (quantum == 0)
        return 0;
    remainder = value % quantum;
    if (remainder != 0)
    {
        size_t extra = quantum - remainder;

        if (value > SIZE_MAX - extra)
            return 0;
        value += extra;
    }
    *rounded = value;
    return 1;
}

static size_t st_sim_pool_header_size(void)
{
    size_t result = 0;

    st_sim_pool_round(sizeof(st_sim_pool_span), ST_SIM_POOL_ALIGN, &result);
    return result;
}

static size_t st_sim_pool_bin_index(size_t capacity)
{
    size_t index = 0;

    while (capacity > 1)
    {
        capacity >>= 1;
        index++;
    }
    return index;
}

static void st_sim_pool_bin_insert(st_sim_pool_span *span)
{
    size_t index = st_sim_pool_bin_index(span->capacity);

    span->free_prev = NULL;
    span->free_next = st_sim_pool_bins[index];
    if (span->free_next)
        span->free_next->free_prev = span;
    st_sim_pool_bins[index] = span;
}

static void st_sim_pool_bin_remove(st_sim_pool_span *span)
{
    size_t index = st_sim_pool_bin_index(span->capacity);

    if (span->free_prev)
        span->free_prev->free_next = span->free_next;
    else
        st_sim_pool_bins[index] = span->free_next;
    if (span->free_next)
        span->free_next->free_prev = span->free_prev;
    span->free_prev = NULL;
    span->free_next = NULL;
}

static void st_sim_pool_update_peaks(void)
{
    if (st_sim_pool_stats.live_allocations >
        st_sim_pool_stats.peak_live_allocations)
        st_sim_pool_stats.peak_live_allocations =
            st_sim_pool_stats.live_allocations;
    if (st_sim_pool_stats.live_bytes > st_sim_pool_stats.peak_live_bytes)
        st_sim_pool_stats.peak_live_bytes = st_sim_pool_stats.live_bytes;
    if (st_sim_pool_stats.live_capacity_bytes >
        st_sim_pool_stats.peak_live_capacity_bytes)
        st_sim_pool_stats.peak_live_capacity_bytes =
            st_sim_pool_stats.live_capacity_bytes;
    if (st_sim_pool_stats.reserved_bytes >
        st_sim_pool_stats.peak_reserved_bytes)
        st_sim_pool_stats.peak_reserved_bytes = st_sim_pool_stats.reserved_bytes;
}

static void st_sim_pool_span_init(st_sim_pool_span *span,
                                  st_sim_pool_block *block, size_t capacity,
                                  st_sim_pool_span *prev,
                                  st_sim_pool_span *next)
{
    span->block = block;
    span->prev = prev;
    span->next = next;
    span->free_prev = NULL;
    span->free_next = NULL;
    span->capacity = capacity;
    span->requested = 0;
    span->used = 0;
}

static st_sim_pool_span *st_sim_pool_add_block(size_t capacity)
{
    size_t first, need, bytes, available;
    size_t header = st_sim_pool_header_size();
    st_sim_pool_block *block;
    st_sim_pool_span *span;

    if (!st_sim_pool_round(sizeof(st_sim_pool_block), ST_SIM_POOL_ALIGN,
                           &first) ||
        first > SIZE_MAX - header ||
        capacity > SIZE_MAX - first - header)
        return NULL;
    need = first + header + capacity;
    if (!st_sim_pool_round(need, add_block_size, &bytes) ||
        bytes > (size_t)PTRDIFF_MAX ||
        bytes > SIZE_MAX - st_sim_pool_stats.reserved_bytes)
        return NULL;
    ST_SIM_POOL_COUNT(backing_requests);
    block = (st_sim_pool_block *)malloc(bytes);
    if (!block)
        return NULL;
    available = bytes - first - header;
    available -= available % ST_SIM_POOL_ALIGN;
    span = (st_sim_pool_span *)((unsigned char *)block + first);
    st_sim_pool_span_init(span, block, available, NULL, NULL);
    block->first = span;
    block->bytes = bytes;
    block->next = st_sim_pool_blocks;
    st_sim_pool_blocks = block;
    st_sim_pool_stats.reserved_bytes += bytes;
    ST_SIM_POOL_COUNT(backing_allocations);
    st_sim_pool_update_peaks();
    st_sim_pool_bin_insert(span);
    return span;
}

static st_sim_pool_span *st_sim_pool_split(st_sim_pool_span *span,
                                          size_t capacity)
{
    size_t header = st_sim_pool_header_size();
    st_sim_pool_span *tail;

    if (capacity > span->capacity || span->capacity - capacity < header ||
        span->capacity - capacity - header < ST_SIM_POOL_ALIGN)
        return NULL;
    tail = (st_sim_pool_span *)((unsigned char *)span + header + capacity);
    st_sim_pool_span_init(tail, span->block,
                          span->capacity - capacity - header, span,
                          span->next);
    if (tail->next)
        tail->next->prev = tail;
    span->next = tail;
    span->capacity = capacity;
    ST_SIM_POOL_COUNT(splits);
    return tail;
}

static st_sim_pool_span *st_sim_pool_merge_next(st_sim_pool_span *span)
{
    st_sim_pool_span *next = span->next;

    if (next && !next->used)
    {
        st_sim_pool_bin_remove(next);
        span->capacity += st_sim_pool_header_size() + next->capacity;
        span->next = next->next;
        if (span->next)
            span->next->prev = span;
        ST_SIM_POOL_COUNT(coalesces);
    }
    return span;
}

static void *st_sim_pool_alloc_locked(size_t bytes)
{
    st_sim_pool_span *span = NULL, *tail;
    size_t capacity, index;
    int reused = 1;

    if (bytes == 0)
        return NULL;
    if (!st_sim_pool_round(bytes, ST_SIM_POOL_ALIGN, &capacity) ||
        capacity > (size_t)PTRDIFF_MAX ||
        st_sim_pool_stats.live_allocations == SIZE_MAX ||
        bytes > SIZE_MAX - st_sim_pool_stats.live_bytes)
        goto fail;
#ifdef ST_SIM_POOL_DIAGNOSTICS
    if (st_sim_pool_fail_after == 0)
        goto fail;
#endif
    index = st_sim_pool_bin_index(capacity);
    for (; index < ST_SIM_POOL_BINS; index++)
    {
        for (span = st_sim_pool_bins[index]; span; span = span->free_next)
            if (span->capacity >= capacity)
                break;
        if (span)
            break;
    }
    if (!span)
    {
        span = st_sim_pool_add_block(capacity);
        reused = 0;
        if (!span)
            goto fail;
    }
    st_sim_pool_bin_remove(span);
    tail = st_sim_pool_split(span, capacity);
    if (tail)
        st_sim_pool_bin_insert(tail);
    span->used = 1;
    span->requested = bytes;
    st_sim_pool_stats.live_allocations++;
    st_sim_pool_stats.live_bytes += bytes;
    st_sim_pool_stats.live_capacity_bytes += span->capacity;
    if (reused)
        ST_SIM_POOL_COUNT(reused_allocations);
#ifdef ST_SIM_POOL_DIAGNOSTICS
    if (st_sim_pool_fail_after != SIZE_MAX)
        st_sim_pool_fail_after--;
#endif
    st_sim_pool_update_peaks();
    return (unsigned char *)span + st_sim_pool_header_size();

fail:
    ST_SIM_POOL_COUNT(failed_requests);
    return NULL;
}

static st_sim_pool_span *st_sim_pool_span_from_pointer(const void *pointer)
{
    return (st_sim_pool_span *)((unsigned char *)pointer -
                                st_sim_pool_header_size());
}

static void st_sim_pool_free_locked(void *pointer)
{
    st_sim_pool_span *span;

    if (!pointer)
        return;
    span = st_sim_pool_span_from_pointer(pointer);
    st_sim_pool_stats.live_allocations--;
    st_sim_pool_stats.live_bytes -= span->requested;
    st_sim_pool_stats.live_capacity_bytes -= span->capacity;
    span->requested = 0;
    span->used = 0;
    span = st_sim_pool_merge_next(span);
    if (span->prev && !span->prev->used)
    {
        st_sim_pool_span *prev = span->prev;

        st_sim_pool_bin_remove(prev);
        prev->capacity += st_sim_pool_header_size() + span->capacity;
        prev->next = span->next;
        if (prev->next)
            prev->next->prev = prev;
        span = prev;
        ST_SIM_POOL_COUNT(coalesces);
    }
    st_sim_pool_bin_insert(span);
}

static void *st_sim_pool_resize_locked(void *pointer, size_t bytes)
{
    st_sim_pool_span *span, *tail;
    size_t capacity, old_capacity, old_requested;
    void *replacement;

    if (!pointer)
        return st_sim_pool_alloc_locked(bytes);
    if (bytes == 0)
    {
        st_sim_pool_free_locked(pointer);
        return NULL;
    }
    span = st_sim_pool_span_from_pointer(pointer);
    old_capacity = span->capacity;
    old_requested = span->requested;
    if (!st_sim_pool_round(bytes, ST_SIM_POOL_ALIGN, &capacity) ||
        capacity > (size_t)PTRDIFF_MAX ||
        (bytes > old_requested &&
         bytes - old_requested > SIZE_MAX - st_sim_pool_stats.live_bytes))
        goto fail;
    if (capacity > span->capacity && span->next && !span->next->used &&
        capacity - span->capacity <=
            st_sim_pool_header_size() + span->next->capacity)
        st_sim_pool_merge_next(span);
    if (capacity <= span->capacity)
    {
        tail = st_sim_pool_split(span, capacity);
        if (tail)
        {
            st_sim_pool_merge_next(tail);
            st_sim_pool_bin_insert(tail);
        }
        span->requested = bytes;
        st_sim_pool_stats.live_bytes -= old_requested;
        st_sim_pool_stats.live_bytes += bytes;
        st_sim_pool_stats.live_capacity_bytes -= old_capacity;
        st_sim_pool_stats.live_capacity_bytes += span->capacity;
        ST_SIM_POOL_COUNT(inplace_resizes);
        st_sim_pool_update_peaks();
        return pointer;
    }
    replacement = st_sim_pool_alloc_locked(bytes);
    if (!replacement)
        return NULL;
    memcpy(replacement, pointer, old_capacity < bytes ? old_capacity : bytes);
    st_sim_pool_free_locked(pointer);
    ST_SIM_POOL_COUNT(moving_resizes);
    return replacement;

fail:
    ST_SIM_POOL_COUNT(failed_requests);
    return NULL;
}

static void *st_sim_pool_alloc(size_t bytes)
{
    void *result;

#ifdef MT
#pragma omp critical(st_sim_pool_lock)
#endif
    {
        ST_SIM_POOL_COUNT(allocation_requests);
        result = st_sim_pool_alloc_locked(bytes);
    }
    return result;
}

static void *st_sim_pool_resize(void *pointer, size_t bytes)
{
    void *result;

#ifdef MT
#pragma omp critical(st_sim_pool_lock)
#endif
    {
        ST_SIM_POOL_COUNT(resize_requests);
        result = st_sim_pool_resize_locked(pointer, bytes);
    }
    return result;
}

static void st_sim_pool_free(void *pointer)
{
    if (!pointer)
        return;
#ifdef MT
#pragma omp critical(st_sim_pool_lock)
#endif
    {
        ST_SIM_POOL_COUNT(free_requests);
        st_sim_pool_free_locked(pointer);
    }
}

static inline size_t st_sim_pool_capacity(const void *pointer)
{
    return pointer ? st_sim_pool_span_from_pointer(pointer)->capacity : 0;
}

void st_sim_pool_dispose(void)
{
#ifdef MT
#pragma omp critical(st_sim_pool_lock)
#endif
    {
        if (st_sim_pool_stats.live_allocations != 0)
        {
            ST_SIM_POOL_COUNT(disposal_rejections);
        }
        else
        {
            st_sim_pool_block *block = st_sim_pool_blocks;

            while (block)
            {
                st_sim_pool_block *next = block->next;

                free(block);
                ST_SIM_POOL_COUNT(backing_frees);
                block = next;
            }
            st_sim_pool_blocks = NULL;
            memset(st_sim_pool_bins, 0, sizeof(st_sim_pool_bins));
            st_sim_pool_stats.reserved_bytes = 0;
        }
    }
}

#ifdef ST_SIM_POOL_DIAGNOSTICS
void st_sim_pool_get_stats(st_sim_pool_statistics *result)
{
#ifdef MT
#pragma omp critical(st_sim_pool_lock)
#endif
    { *result = st_sim_pool_stats; }
}

void st_sim_pool_test_fail_after(size_t successful_allocations)
{
#ifdef MT
#pragma omp critical(st_sim_pool_lock)
#endif
    { st_sim_pool_fail_after = successful_allocations; }
}

size_t st_sim_pool_validate(void)
{
    size_t errors = 0;

#ifdef MT
#pragma omp critical(st_sim_pool_lock)
#endif
    {
        size_t reserved = 0, live = 0, requested = 0, capacity = 0;
        size_t free_spans = 0, indexed_free_spans = 0;
        size_t header = st_sim_pool_header_size();
        st_sim_pool_block *block;

        for (block = st_sim_pool_blocks; block; block = block->next)
        {
            st_sim_pool_span *span, *prev = NULL;
            size_t first, offset, iterations = 0;

            st_sim_pool_round(sizeof(st_sim_pool_block), ST_SIM_POOL_ALIGN,
                              &first);
            reserved += block->bytes;
            offset = first;
            for (span = block->first; span; span = span->next)
            {
                if (++iterations > block->bytes / header)
                {
                    errors++;
                    break;
                }
                if ((unsigned char *)span != (unsigned char *)block + offset ||
                    span->block != block || span->prev != prev ||
                    span->capacity % ST_SIM_POOL_ALIGN != 0 ||
                    span->capacity > block->bytes - offset - header ||
                    span->requested > span->capacity ||
                    (span->used != 0 && span->used != 1))
                    errors++;
                if (span->used)
                {
                    live++;
                    requested += span->requested;
                    capacity += span->capacity;
                    if (span->free_prev || span->free_next)
                        errors++;
                }
                else
                {
                    free_spans++;
                    if (span->requested != 0 ||
                        (span->next && !span->next->used))
                        errors++;
                }
                offset += header + span->capacity;
                prev = span;
            }
            if (offset > block->bytes ||
                block->bytes - offset >= ST_SIM_POOL_ALIGN)
                errors++;
        }
        for (size_t i = 0; i < ST_SIM_POOL_BINS; i++)
        {
            st_sim_pool_span *span, *prev = NULL;
            size_t iterations = 0;

            for (span = st_sim_pool_bins[i]; span; span = span->free_next)
            {
                if (++iterations > free_spans)
                {
                    errors++;
                    break;
                }
                indexed_free_spans++;
                if (span->used || st_sim_pool_bin_index(span->capacity) != i ||
                    span->free_prev != prev)
                    errors++;
                prev = span;
            }
        }
        if (reserved != st_sim_pool_stats.reserved_bytes ||
            live != st_sim_pool_stats.live_allocations ||
            requested != st_sim_pool_stats.live_bytes ||
            capacity != st_sim_pool_stats.live_capacity_bytes ||
            free_spans != indexed_free_spans)
            errors++;
    }
    return errors;
}
#endif

#undef ST_SIM_POOL_COUNT

static void *st_sim_alloc_array(size_t count, size_t size)
{
    if (!count || !size || count > SIZE_MAX / size)
        return NULL;
    return st_sim_pool_alloc(count * size);
}

static void *st_sim_resize_array(void *ptr, size_t count, size_t size)
{
    if (!count || !size || count > SIZE_MAX / size)
        return NULL;
    return st_sim_pool_resize(ptr, count * size);
}

static char *st_sim_text_dup(const char *text)
{
    size_t bytes;
    char *copy;

    if (!text)
        return NULL;
    bytes = strlen(text);
    if (bytes >= INT_MAX || bytes == SIZE_MAX)
        return NULL;
    bytes++;
    copy = st_sim_pool_alloc(bytes);
    if (copy)
        memcpy(copy, text, bytes);
    return copy;
}

#define st_sim_alloc(n, t) \
    ((t *)st_sim_alloc_array((size_t)(n), sizeof(t)))
#define st_sim_free(p) (st_sim_pool_free((void *)(p)), (p) = NULL)

typedef struct
{
    char pad;
    TChrom value;
} st_sim_chrom_alignment;

typedef struct
{
    char pad;
    THom value;
} st_sim_hom_alignment;

static int st_sim_layout_append(size_t *bytes, size_t alignment,
                                size_t count, size_t size, size_t *offset)
{
    size_t padding, payload;

    if (!alignment || !count || !size || count > SIZE_MAX / size)
        return -1;
    padding = *bytes % alignment;
    if (padding)
        padding = alignment - padding;
    payload = count * size;
    if (*bytes > SIZE_MAX - padding ||
        *bytes + padding > SIZE_MAX - payload)
        return -1;
    *offset = *bytes + padding;
    *bytes = *offset + payload;
    return 0;
}

static int st_sim_loci_capacity(int current, size_t needed, int *capacity)
{
    size_t result;

    if (current < 0 || needed > INT_MAX || NoLociInit <= 0)
        return -1;
    if (needed <= (size_t)current)
    {
        *capacity = current;
        return 0;
    }
    result = (size_t)current;
    if (result < (size_t)NoLociInit)
        result = (size_t)NoLociInit;
    while (result < needed)
    {
        if (result > (size_t)INT_MAX / 2)
            result = INT_MAX;
        else
            result *= 2;
    }
    *capacity = (int)result;
    return 0;
}

static void st_sim_storage_views(st_sim_storage *storage)
{
    size_t i;

    for (i = 0; i < storage->NoSlots; i++)
        storage->Hom[i].Loc = storage->Loc + storage->Hom[i].Offset;
}

static int st_sim_storage_grow(st_sim_storage *storage, size_t needed)
{
    size_t capacity, maximum = SIZE_MAX / sizeof(TLoc);
    TLoc *loci;

    if (needed <= storage->NoLociAlloc)
        return 0;
    if (needed > maximum)
        return -1;
    capacity = storage->NoLociAlloc;
    if (capacity <= maximum / 2)
        capacity *= 2;
    else
        capacity = maximum;
    if (capacity < needed)
        capacity = needed;
    loci = (TLoc *)st_sim_resize_array(storage->Loc, capacity, sizeof(TLoc));
    if (!loci)
        return -1;
    storage->Loc = loci;
    storage->NoLociAlloc = st_sim_pool_capacity(loci) / sizeof(TLoc);
    st_sim_storage_views(storage);
    return 0;
}

static int st_sim_hom_reserve(THom *hom, size_t needed)
{
    st_sim_storage *storage = hom->Owner;
    size_t index, delta, total, i;
    int capacity;

    if (!storage || st_sim_loci_capacity(hom->NoLociAlloc, needed, &capacity))
        return -1;
    if (capacity == hom->NoLociAlloc)
        return 0;
    index = (size_t)(hom - storage->Hom);
    if (index >= storage->NoSlots)
        return -1;
    delta = (size_t)capacity - (size_t)hom->NoLociAlloc;
    if (storage->NoLociReserved > SIZE_MAX - delta)
        return -1;
    total = storage->NoLociReserved + delta;
    if (st_sim_storage_grow(storage, total))
        return -1;
    for (i = storage->NoSlots; i > index + 1;)
    {
        THom *tail = &storage->Hom[--i];
        size_t offset = tail->Offset + delta;

        if (tail->NoLociUsed)
            memmove(storage->Loc + offset, storage->Loc + tail->Offset,
                    (size_t)tail->NoLociUsed * sizeof(TLoc));
        tail->Offset = offset;
        tail->Loc = storage->Loc + offset;
    }
    hom->NoLociAlloc = capacity;
    storage->NoLociReserved = total;
    return 0;
}

static int st_sim_prepare_capacity(const st_sim_storage *storage,
                                   const st_sim_storage *left,
                                   const st_sim_storage *right, size_t index,
                                   int mating, int *capacity)
{
    size_t needed;

    if (mating)
    {
        size_t start = (index / storage->NoHoms) * storage->NoHoms;
        const st_sim_storage *parent =
            index % storage->NoHoms ? right : left;

        if (parent->Hom[start].NoLociUsed < 0 ||
            parent->Hom[start + 1].NoLociUsed < 0)
            return -1;
        needed = (size_t)parent->Hom[start].NoLociUsed +
                 (size_t)parent->Hom[start + 1].NoLociUsed;
        if (NoLociInit <= 0 || needed > (size_t)INT_MAX - NoLociInit)
            return -1;
        needed += (size_t)NoLociInit;
    }
    else
    {
        if (left->Hom[index].NoLociUsed < 0)
            return -1;
        needed = (size_t)left->Hom[index].NoLociUsed;
    }
    return st_sim_loci_capacity(storage->Hom[index].NoLociAlloc, needed,
                                capacity);
}

static int st_sim_prepare_individual(TInd *destination, const TInd *left,
                                     const TInd *right, int mating)
{
    st_sim_storage *storage = destination->Storage;
    const st_sim_storage *a = left->Storage;
    const st_sim_storage *b = right ? right->Storage : a;
    size_t total = 0, offset, i;

    if (!storage || !a || !b ||
        storage->NoChroms != a->NoChroms || storage->NoHoms != a->NoHoms ||
        storage->NoChroms != b->NoChroms || storage->NoHoms != b->NoHoms ||
        (mating && storage->NoHoms != 2))
        return -1;
    for (i = 0; i < storage->NoSlots; i++)
    {
        int capacity;

        if (st_sim_prepare_capacity(storage, a, b, i, mating, &capacity) ||
            total > SIZE_MAX - (size_t)capacity)
            return -1;
        total += (size_t)capacity;
    }
    if (total == storage->NoLociReserved)
        return 0;
    if (st_sim_storage_grow(storage, total))
        return -1;
    offset = total;
    for (i = storage->NoSlots; i > 0;)
    {
        THom *hom = &storage->Hom[--i];
        int capacity;

        if (st_sim_prepare_capacity(storage, a, b, i, mating, &capacity))
            return -1;
        offset -= (size_t)capacity;
        if (hom->NoLociUsed && offset != hom->Offset)
            memmove(storage->Loc + offset, storage->Loc + hom->Offset,
                    (size_t)hom->NoLociUsed * sizeof(TLoc));
        hom->NoLociAlloc = capacity;
        hom->Offset = offset;
        hom->Loc = storage->Loc + offset;
    }
    storage->NoLociReserved = total;
    return 0;
}

static int st_sim_compact_capacity(const THom *hom, int *capacity)
{
    size_t result, minimum;

    if (NoLociInit <= 0 || hom->NoLociUsed < 0 ||
        hom->NoLociUsed > hom->NoLociAlloc)
        return -1;
    minimum = (size_t)NoLociInit;
    result = hom->NoLociUsed ? (size_t)hom->NoLociUsed : 1;
    if (result % minimum)
        result += minimum - result % minimum;
    if (result > (size_t)hom->NoLociAlloc)
        result = (size_t)hom->NoLociAlloc;
    *capacity = (int)result;
    return 0;
}

static int st_sim_storage_compact(st_sim_storage *storage)
{
    size_t total = 0, offset = 0, i;
    TLoc *loci;

    if (!storage)
        return -1;
    for (i = 0; i < storage->NoSlots; i++)
    {
        int capacity;

        if (st_sim_compact_capacity(&storage->Hom[i], &capacity) ||
            total > SIZE_MAX - (size_t)capacity)
            return -1;
        total += (size_t)capacity;
    }
    for (i = 0; i < storage->NoSlots; i++)
    {
        THom *hom = &storage->Hom[i];
        int capacity;

        if (st_sim_compact_capacity(hom, &capacity))
            return -1;
        if (hom->NoLociUsed && offset != hom->Offset)
            memmove(storage->Loc + offset, storage->Loc + hom->Offset,
                    (size_t)hom->NoLociUsed * sizeof(TLoc));
        hom->NoLociAlloc = capacity;
        hom->Offset = offset;
        offset += (size_t)capacity;
    }
    storage->NoLociReserved = total;
    loci = (TLoc *)st_sim_resize_array(storage->Loc, total, sizeof(TLoc));
    if (loci)
        storage->Loc = loci;
    storage->NoLociAlloc = st_sim_pool_capacity(storage->Loc) / sizeof(TLoc);
    st_sim_storage_views(storage);
    return 0;
}

static int st_sim_text_alias(const char *text, const char *buffer)
{
    uintptr_t source, destination;

    if (!text || !buffer)
        return 0;
    source = (uintptr_t)(const void *)text;
    destination = (uintptr_t)(const void *)buffer;
    return source >= destination &&
           source - destination < st_sim_pool_capacity(buffer);
}

static int st_sim_identity_set(TInd *ind, const char *name,
                               const char *pedigree)
{
    const char *source_pedigree = pedigree ? pedigree : name;
    char *new_name = ind->name, *new_pedigree = ind->pedigree;
    size_t name_bytes, pedigree_bytes;

    if (!name)
        return -1;
    name_bytes = strlen(name);
    pedigree_bytes = strlen(source_pedigree);
    if (name_bytes >= INT_MAX || pedigree_bytes >= INT_MAX)
        return -1;
    name_bytes++;
    pedigree_bytes++;
    if (!new_name || st_sim_pool_capacity(new_name) < name_bytes)
        new_name = (char *)st_sim_pool_alloc(name_bytes);
    if (!new_name)
        goto fail;
    if (!new_pedigree || st_sim_pool_capacity(new_pedigree) < pedigree_bytes)
        new_pedigree = (char *)st_sim_pool_alloc(pedigree_bytes);
    if (!new_pedigree)
        goto fail;
    if (st_sim_text_alias(source_pedigree, new_name) &&
        st_sim_text_alias(name, new_pedigree))
    {
        new_name = (char *)st_sim_pool_alloc(name_bytes);
        if (!new_name)
            goto fail;
    }
    if (st_sim_text_alias(source_pedigree, new_name))
    {
        memmove(new_pedigree, source_pedigree, pedigree_bytes);
        memmove(new_name, name, name_bytes);
    }
    else
    {
        memmove(new_name, name, name_bytes);
        memmove(new_pedigree, source_pedigree, pedigree_bytes);
    }
    if (ind->name && ind->name != new_name)
        st_sim_pool_free(ind->name);
    if (ind->pedigree && ind->pedigree != new_pedigree)
        st_sim_pool_free(ind->pedigree);
    ind->name = new_name;
    ind->pedigree = new_pedigree;
    ps_identity_observe(new_name);
    return 0;

fail:

    if (new_name && new_name != ind->name)
        st_sim_pool_free(new_name);
    if (new_pedigree && new_pedigree != ind->pedigree)
        st_sim_pool_free(new_pedigree);
    return -1;
}
static int st_sim_identity_copy(TInd *dst, const TInd *src)
{
    if (!src->name)
        return -1;
    return st_sim_identity_set(dst, src->name, src->pedigree);
}
static int st_sim_identity_new(TInd *ind, const char *population,
                               uint64_t serial, const TInd *left,
                               const TInd *right, const char *operation)
{
    const char *a = left ? left->pedigree : NULL;
    const char *b = right ? right->pedigree : NULL;
    char *name = ind->name, *pedigree = ind->pedigree;
    char serial_text[32], *cursor;
    size_t la, lb, lo, population_len, pedigree_len, serial_len;
    size_t name_bytes, pedigree_bytes;
    int serial_chars;

    if (!population || (right && !left) || (left && !right && !operation))
        return -1;
    if (left && !a)
        a = left->name;
    if (right && !b)
        b = right->name;
    if ((left && !a) || (right && !b))
        return -1;
    la = a ? strlen(a) : 0;
    lb = b ? strlen(b) : 0;
    lo = operation ? strlen(operation) : 0;
    if (la > (size_t)INT_MAX - 8 || lb > (size_t)INT_MAX - 8 - la ||
        lo > (size_t)INT_MAX - 8 - la - lb)
        return -1;
    pedigree_len = right ? la + lb + 5 : left ? la + lo + 2 : 0;
    serial_chars = snprintf(serial_text, sizeof(serial_text), ".%llu",
                            (unsigned long long)serial);
    if (serial_chars < 0 || (size_t)serial_chars >= sizeof(serial_text))
        return -1;
    serial_len = (size_t)serial_chars;
    population_len = strlen(population);
    if (population_len > (size_t)INT_MAX - serial_len - 2 ||
        pedigree_len > (size_t)INT_MAX - serial_len - 2 - population_len)
        return -1;
    name_bytes = population_len + serial_len + pedigree_len + 1;
    pedigree_bytes = left ? pedigree_len + 1 : name_bytes;
    if (!name || st_sim_pool_capacity(name) < name_bytes)
        name = (char *)st_sim_pool_alloc(name_bytes);
    if (!name)
        goto fail;
    if (!pedigree || st_sim_pool_capacity(pedigree) < pedigree_bytes ||
        (left && (st_sim_text_alias(a, pedigree) ||
                  st_sim_text_alias(b, pedigree) ||
                  st_sim_text_alias(population, pedigree) ||
                  st_sim_text_alias(operation, pedigree))))
        pedigree = (char *)st_sim_pool_alloc(pedigree_bytes);
    if (!pedigree)
        goto fail;
    if (left)
    {
        cursor = pedigree;
        *cursor++ = '(';
        memcpy(cursor, a, la);
        cursor += la;
        if (right)
        {
            memcpy(cursor, "_x_", 3);
            cursor += 3;
            memcpy(cursor, b, lb);
            cursor += lb;
        }
        else
        {
            memcpy(cursor, operation, lo);
            cursor += lo;
        }
        *cursor++ = ')';
        *cursor = '\0';
    }
    memmove(name, population, population_len);
    memcpy(name + population_len, serial_text, serial_len);
    if (left)
        memcpy(name + population_len + serial_len, pedigree, pedigree_len + 1);
    else
    {
        name[population_len + serial_len] = '\0';
        memcpy(pedigree, name, name_bytes);
    }
    if (ind->name && ind->name != name)
        st_sim_pool_free(ind->name);
    if (ind->pedigree && ind->pedigree != pedigree)
        st_sim_pool_free(ind->pedigree);
    if (ind->info)
        st_sim_pool_free(ind->info);
    ind->name = name;
    ind->pedigree = pedigree;
    ind->info = NULL;
    if (st_cross_left && right && serial >= st_cross_first &&
        serial - st_cross_first < st_cross_count)
    {
        size_t i = (size_t)(serial - st_cross_first);
        st_cross_left[i] = left->name;
        st_cross_right[i] = right->name;
    }
    return 0;

fail:

    if (name && name != ind->name)
        st_sim_pool_free(name);
    if (pedigree && pedigree != ind->pedigree)
        st_sim_pool_free(pedigree);
    return -1;
}
static int st_sim_name_cmp(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}
static int st_sim_export_identity(TPop *pop, char ***names, char ***pedigree)
{
    size_t n = (size_t)pop->NoInds;
    char **order = calloc(n ? n : 1, sizeof(*order));
    if (!order)
        return -1;
    for (size_t i = 0; i < n; i++)
    {
        if (!pop->Ind[i].name)
        {
            free(order);
            return -1;
        }
        order[i] = pop->Ind[i].name;
    }
    qsort(order, n, sizeof(*order), st_sim_name_cmp);
    for (size_t i = 1; i < n; i++)
        if (!strcmp(order[i - 1], order[i]))
        {
            free(order);
            return -2;
        }
    free(order);
    *names = calloc(n ? n : 1, sizeof(**names));
    *pedigree = calloc(n ? n : 1, sizeof(**pedigree));
    if (!*names || !*pedigree)
        return -1;
    for (size_t i = 0; i < n; i++)
    {
        (*names)[i] = st_text_dup(pop->Ind[i].name);
        (*pedigree)[i] = st_text_dup(pop->Ind[i].pedigree);
        if (!(*names)[i] || !(*pedigree)[i])
            return -1;
    }
    return 0;
}
int ps_population_exists(const char *name)
{
    for (TPopListEl *el = PopListStart; el; el = el->Next)
        if (!strcmp(el->Pop->Name, name))
            return 1;
    return 0;
}
SEXP ps_population_pedigree(const char *name)
{
    TPop *pop = NULL;
    for (TPopListEl *el = PopListStart; el; el = el->Next)
        if (!strcmp(el->Pop->Name, name))
        {
            pop = el->Pop;
            break;
        }
    if (!pop)
        Rf_error("Simulation population '%s' does not exist", name);
    if (pop->NoInds > INT_MAX)
        Rf_error("Population exceeds R data-frame limits");
    int n = (int)pop->NoInds;
    SEXP ans = PROTECT(allocVector(VECSXP, 2)),
         nm = PROTECT(allocVector(STRSXP, 2)),
         ids = PROTECT(allocVector(STRSXP, n)),
         ped = PROTECT(allocVector(STRSXP, n)),
         rn = PROTECT(allocVector(INTSXP, n ? 2 : 0));
    for (int i = 0; i < n; i++)
    {
        SET_STRING_ELT(ids, i,
                       pop->Ind[i].name ? mkCharCE(pop->Ind[i].name, CE_UTF8)
                                        : NA_STRING);
        SET_STRING_ELT(ped, i,
                       pop->Ind[i].pedigree
                           ? mkCharCE(pop->Ind[i].pedigree, CE_UTF8)
                           : NA_STRING);
    }
    SET_STRING_ELT(nm, 0, mkChar("name"));
    SET_STRING_ELT(nm, 1, mkChar("pedigree"));
    SET_VECTOR_ELT(ans, 0, ids);
    SET_VECTOR_ELT(ans, 1, ped);
    setAttrib(ans, R_NamesSymbol, nm);
    if (n)
    {
        INTEGER(rn)[0] = NA_INTEGER;
        INTEGER(rn)[1] = -n;
    }
    setAttrib(ans, R_RowNamesSymbol, rn);
    setAttrib(ans, R_ClassSymbol, mkString("data.frame"));
    UNPROTECT(5);
    return ans;
}

static int st_console_output_allowed(void)
{
#ifdef MT
    if (omp_in_parallel())
        return 0;
#endif

    return 1;
}

static int st_parallel_info_error = 0;

static void st_parallel_info_error_clear(void)
{
    st_parallel_info_error = 0;
}

static void st_parallel_info_error_set(void)
{
#ifdef MT
#pragma omp atomic write
#endif

    st_parallel_info_error = 1;
}

static int st_parallel_info_error_get(void)
{
    return st_parallel_info_error;
}

void info(int level, char *info, ...)
{
    if (INFOLEVEL >= level)
    {

        if (!st_console_output_allowed())
        {
            if (-2 >= level)
                st_parallel_info_error_set();

            return;
        }

        va_list ap;
        char Info[200];

        va_start(ap, info);
        vsnprintf(Info, sizeof(Info), info, ap);
        va_end(ap);

        if (0 <= level)
            Rprintf("M: %s", Info);
        else if (-1 == level)
            Rprintf("W: %s", Info);
        else if (-2 == level)
            Rprintf("E: %s", Info);
    }
}

void r_info(int *level, char **message)
{
    info(*level, "%s\n", *message);
}

void plabsim_info_cctime_get(char **message)
{
    snprintf(*message, strlen(CCTIME) + strlen(CCDATE) + 4, "%s - %s", CCTIME,
             CCDATE);
}

void print_version(void)
{
    info(0, "%s (%s): Rev. %s \n", AUTHOR, DATE, REVISION);
    info(0, "compiled: %s - %s\n", CCTIME, CCDATE);
}

void print_version_2(char **message)
{
    char *co_date_complete = DATE;
    char co_date[20];
    strncpy(co_date, co_date_complete + 7, 10);
    co_date[10] = '\0';
    info(0, "plabsim c %s\n", co_date);
    info(0, "plabsim r %s\n", *message);
}


#define FILE st_io_stream
#define fopen st_io_open
#define fclose st_io_close
#define fprintf st_io_printf
#define fwrite st_io_write
#define fgetc st_io_getc
#define fgets st_io_gets
#define feof st_io_eof
#define ferror st_io_error
#define fscanf st_io_scan
#define remove st_io_remove

FILE *fileopen(char *path, char *mode)
{
    FILE *fp;

    if ((fp = fopen(path, mode)) == NULL)
        info(-2, "Can't open the file %s\n", path);

    return fp;
}

void fileclose(FILE *fp)
{
    if (fclose(fp) != 0)
        info(-2, "Can't close the file stream\n");
}

void set_info_level(int *level)
{
    INFOLEVEL = *level;
    if (INFOLEVEL < -2)
        INFOLEVEL = -2;
    if (INFOLEVEL > 10)
        INFOLEVEL = 10;
    info(2, "Info level set to %d\n", INFOLEVEL);
}

void get_info_level(int *i)
{
    *i = INFOLEVEL;
}

void find_list_element(char *Name, TPopListEl **PopListEl)
{
    *PopListEl = PopListStart;

    while ((*PopListEl != NULL) && (strcmp((*PopListEl)->Pop->Name, Name) != 0))
    {
        (*PopListEl) = (*PopListEl)->Next;
    }
}

void remove_individuum(TInd *DelInd)
{
    if (DelInd->Storage)
    {
        st_sim_pool_free(DelInd->Storage->Loc);
        st_sim_pool_free(DelInd->Storage);
    }
    DelInd->Storage = NULL;
    DelInd->Chrom = NULL;

    if (NULL != DelInd->GValue)
        st_sim_free(DelInd->GValue);
    if (NULL != DelInd->PValue)
        st_sim_free(DelInd->PValue);
    if (NULL != DelInd->genotype)
        st_sim_free(DelInd->genotype);
    if (NULL != DelInd->info)
        st_sim_free(DelInd->info);
    st_sim_pool_free(DelInd->name);
    DelInd->name = NULL;
    st_sim_pool_free(DelInd->pedigree);
    DelInd->pedigree = NULL;
}

int new_individuum(TInd *NewInd)
{
    st_sim_storage *storage;
    size_t bytes = sizeof(st_sim_storage), chrom_offset, hom_offset;
    size_t slots, loci, i;

    NewInd->Chrom = NULL;
    NewInd->Storage = NULL;
    NewInd->GValue = NULL;
    NewInd->PValue = NULL;
    NewInd->genotype = NULL;
    NewInd->info = NULL;
    NewInd->name = NULL;
    NewInd->pedigree = NULL;

    if (NoChroms <= 0 || NoHoms <= 0 || NoLociInit <= 0 ||
        (size_t)NoChroms > SIZE_MAX / (size_t)NoHoms)
        return -1;
    slots = (size_t)NoChroms * (size_t)NoHoms;
    if (slots > SIZE_MAX / (size_t)NoLociInit)
        return -1;
    loci = slots * (size_t)NoLociInit;
    if (loci > SIZE_MAX / sizeof(TLoc) ||
        st_sim_layout_append(&bytes, offsetof(st_sim_chrom_alignment, value),
                             (size_t)NoChroms, sizeof(TChrom), &chrom_offset) ||
        st_sim_layout_append(&bytes, offsetof(st_sim_hom_alignment, value),
                             slots, sizeof(THom), &hom_offset))
        return -1;
    storage = (st_sim_storage *)st_sim_pool_alloc(bytes);
    if (!storage)
        return -1;
    storage->Loc = (TLoc *)st_sim_alloc_array(loci, sizeof(TLoc));
    if (!storage->Loc)
    {
        st_sim_pool_free(storage);
        return -1;
    }
    storage->NoChroms = (size_t)NoChroms;
    storage->NoHoms = (size_t)NoHoms;
    storage->NoSlots = slots;
    storage->NoLociAlloc =
        st_sim_pool_capacity(storage->Loc) / sizeof(TLoc);
    storage->NoLociReserved = loci;
    storage->Chrom = (TChrom *)((unsigned char *)storage + chrom_offset);
    storage->Hom = (THom *)((unsigned char *)storage + hom_offset);
    for (i = 0; i < storage->NoChroms; i++)
        storage->Chrom[i].Hom = storage->Hom + i * storage->NoHoms;
    for (i = 0; i < slots; i++)
    {
        THom *hom = &storage->Hom[i];

        hom->NoLociAlloc = NoLociInit;
        hom->NoLociUsed = 0;
        hom->Offset = i * (size_t)NoLociInit;
        hom->Owner = storage;
        hom->Loc = storage->Loc + hom->Offset;
    }
    NewInd->Storage = storage;
    NewInd->Chrom = storage->Chrom;
    return 0;
}

TInd *get_individual(TPop *Pop, POPULATION nrInd)
{
    return &(Pop->Ind[nrInd]);
}

static void st_clear_population_genotype(TPop *Pop)
{
    POPULATION ind;

    if (NULL == Pop)
        return;

    for (ind = 0; ind < Pop->NoIndsAlloc; ind++)
    {
        st_sim_pool_free(Pop->Ind[ind].genotype);
        Pop->Ind[ind].genotype = NULL;
    }
}

static void st_clear_population_pvalue(TPop *Pop)
{
    POPULATION ind;

    if (NULL == Pop)
        return;

    for (ind = 0; ind < Pop->NoIndsAlloc; ind++)
    {
        st_sim_pool_free(Pop->Ind[ind].PValue);
        Pop->Ind[ind].PValue = NULL;
    }
}

static void st_clear_population_gvalue(TPop *Pop)
{
    POPULATION ind;

    if (NULL == Pop)
        return;

    for (ind = 0; ind < Pop->NoIndsAlloc; ind++)
    {
        st_sim_pool_free(Pop->Ind[ind].GValue);
        Pop->Ind[ind].GValue = NULL;
    }

    st_clear_population_pvalue(Pop);
}

static void st_clear_all_population_gvalue(void)
{
    TPopListEl *el = PopListStart;

    while (NULL != el)
    {
        st_clear_population_gvalue(el->Pop);
        el = el->Next;
    }
}

static int st_population_all_genotyped(TPop *Pop)
{
    POPULATION ind;

    if (NULL == Pop)
        return 0;

    if ((0 < Pop->NoInds) && (NULL == Pop->Ind))
        return 0;

    for (ind = 0; ind < Pop->NoInds; ind++)
    {
        if (NULL == Pop->Ind[ind].genotype)
            return 0;
    }

    return 1;
}

static int st_population_all_gvalue(TPop *Pop)
{
    POPULATION ind;

    if (NULL == Pop)
        return 0;

    if ((0 < Pop->NoInds) && (NULL == Pop->Ind))
        return 0;

    for (ind = 0; ind < Pop->NoInds; ind++)
    {
        if (NULL == Pop->Ind[ind].GValue)
            return 0;
    }

    return 1;
}

static int st_population_all_pvalue(TPop *Pop)
{
    POPULATION ind;

    if (NULL == Pop)
        return 0;

    if ((0 < Pop->NoInds) && (NULL == Pop->Ind))
        return 0;

    for (ind = 0; ind < Pop->NoInds; ind++)
    {
        if (NULL == Pop->Ind[ind].PValue)
            return 0;
    }

    return 1;
}

static int st_copy_value_array(float **target, const float *source)
{
    float *tmp;
    size_t n;

    if (NULL == target)
        return 0;

    if (NULL == source)
    {
        st_sim_pool_free(*target);
        *target = NULL;
        return 1;
    }

    n = (size_t)NrEffects + 1;

    if (n > SIZE_MAX / sizeof(*tmp))
        return 0;

    tmp = (float *)st_sim_alloc_array(n, sizeof(*tmp));

    if (NULL == tmp)
        return 0;

    memcpy(tmp, source, n * sizeof(*tmp));
    st_sim_pool_free(*target);
    *target = tmp;

    return 1;
}

static int st_effect_position(char **EffName, int *eff);
static int st_selection_criterion(char **criterion);
static float *st_individual_value_array(TInd *Ind, int criterion);
static int st_population_selection_ready(TPop *Pop, int criterion, int eff);

void set_nolociinit(int *NewLociInit)
{
    if (0 < *NewLociInit)
        NoLociInit = *NewLociInit;
}

void remove_map(void)
{
    TPopListEl *PopListEl;
    POPULATION ind;
    int chrom, i;
    int mapRemoved = 0;

    PopListEl = PopListStart;

    while (NULL != PopListEl)
    {
        for (ind = 0; ind < PopListEl->Pop->NoIndsAlloc; ind++)
        {
            st_sim_pool_free(PopListEl->Pop->Ind[ind].GValue);
            PopListEl->Pop->Ind[ind].GValue = NULL;

            st_sim_pool_free(PopListEl->Pop->Ind[ind].PValue);
            PopListEl->Pop->Ind[ind].PValue = NULL;

            st_sim_pool_free(PopListEl->Pop->Ind[ind].genotype);
            PopListEl->Pop->Ind[ind].genotype = NULL;
        }

        PopListEl = PopListEl->Next;
    }

    effmap_remove_all();

    if (NULL != PopDescr)
    {
        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            if (0 < PopDescr[chrom].NrMappoints)
            {
                mapRemoved = 1;

                for (i = 0; i < PopDescr[chrom].NrMappoints; i++)
                {
                    free(PopDescr[chrom].Mappoint[i].name);

                    free(PopDescr[chrom].Mappoint[i].class);
                }

                free(PopDescr[chrom].Mappoint);

                PopDescr[chrom].Mappoint = NULL;
                PopDescr[chrom].NrMappoints = 0;
            }
        }
    }

    cMP = 0;
    cMPhom = 0;

    if (mapRemoved)
        info(1, "Linkage map removed\n");
}

static int set_genome_par_impl(int *_NoChroms, int *_NoHoms, double *_ChromLen)
{
    time_t ltime;
    char tme[64];
    TChromDescr *tmpPopDescr;
    int newNoChroms = *_NoChroms;
    int newNoHoms = *_NoHoms;
    int chrom;

    time(&ltime);
    strcpy(tme, ctime(&ltime));

    if (newNoChroms <= 0)
    {
        info(-2, "Number of chromosomes must be positive\n");
        return 0;
    }

    if (newNoHoms <= 0)
    {
        info(-2, "Number of homologues must be positive\n");
        return 0;
    }

    for (chrom = 0; chrom < newNoChroms; chrom++)
    {
        if (!isfinite(_ChromLen[chrom]) || (_ChromLen[chrom] <= 0))
        {
            info(-2, "Chromosome length must be finite and positive\n");
            return 0;
        }
    }

    if (NULL != PopListStart)
    {
        info(-2, "First remove all populations\n");
        return 0;
    }

    if (NrEffects > 0)
    {
        info(-2, "First remove all effects\n");
        return 0;
    }

    if ((PopDescr != NULL) && (cMP > 0))
    {
        info(-2, "First remove previously loaded linkage map\n");
        return 0;
    }

    if ((size_t)newNoChroms > SIZE_MAX / sizeof(TChromDescr))
    {
        info(-2, "Memory allocation failed in set_genome_par\n");
        return 0;
    }

    if (NULL == PopDescr)
        tmpPopDescr =
            (TChromDescr *)malloc((size_t)newNoChroms * sizeof(TChromDescr));
    else
        tmpPopDescr = (TChromDescr *)realloc(PopDescr, (size_t)newNoChroms *
                                                           sizeof(TChromDescr));

    if (NULL == tmpPopDescr)
    {
        info(-2, "Memory allocation failed in set_genome_par\n");
        return 0;
    }

    PopDescr = tmpPopDescr;
    NoChroms = newNoChroms;
    NoHoms = newNoHoms;

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        PopDescr[chrom].ChromLen = _ChromLen[chrom];

        PopDescr[chrom].NrMappoints = 0;
        PopDescr[chrom].Mappoint = NULL;
    }

    info(0, "%d chromosomes defined\n", NoChroms);

    return 1;
}

void set_genome_par(int *_NoChroms, int *_NoHoms, double *_ChromLen)
{
    (void)set_genome_par_impl(_NoChroms, _NoHoms, _ChromLen);
}

void optimize_population(char **PopNames)
{
    TPopListEl *PopListEl;
    char *PopName = strtok(*PopNames, "[' ','\n']");

    while (NULL != PopName)
    {
        TPop *Pop;
        POPULATION ind;

        find_list_element(PopName, &PopListEl);
        if (PopListEl == NULL)
        {
            info(-2, "Population not defined: %s\n", PopName);
            return;
        }
        Pop = PopListEl->Pop;
        for (ind = 0; ind < Pop->NoIndsAlloc; ind++)
        {
            if (st_sim_storage_compact(Pop->Ind[ind].Storage))
            {
                info(-2, "Cannot compact population segment storage\n");
                return;
            }
        }
        if (Pop->NoInds != Pop->NoIndsAlloc)
        {
            info(-1, "NoIndsAlloc is greater than NoInds\n");
            info(-1, "use resize to remove unused individuals\n");
        }
        PopName = strtok(NULL, "[' ','\n']");
    }
}

void int_resize_population(char **PopName, POPULATION newSize)
{
    POPULATION ind;
    TPopListEl *PopListEl;
    TPop *Pop = NULL;
    TInd *tmpInd;

    find_list_element(*PopName, &PopListEl);

    if (PopListEl == NULL)
    {
        info(-2, "Population not defined: %s\n", *PopName);
        return;
    }

    if (newSize < 0)
    {
        info(-2, "Population size must be a positive number \n");
        return;
    }

    Pop = PopListEl->Pop;

    if (newSize > Pop->NoIndsAlloc)
    {
        tmpInd = (TInd *)st_sim_resize_array(
            Pop->Ind, (size_t)newSize, sizeof(TInd));

        if (NULL == tmpInd)
        {
            info(-2, "Memory allocation failed in int_resize_population\n");
            return;
        }

        Pop->Ind = tmpInd;

        for (ind = Pop->NoIndsAlloc; ind < newSize; ind++)
            if (0 != new_individuum(&(Pop->Ind[ind])))
                goto fail;

        Pop->NoIndsAlloc = newSize;
    }

    if (newSize < Pop->NoIndsAlloc)
    {
        if (Pop->NoInds > newSize)
            Pop->NoInds = newSize;

        for (ind = newSize; ind < Pop->NoIndsAlloc; ind++)
            remove_individuum(&(Pop->Ind[ind]));

        if (0 == newSize)
        {
            st_sim_free(Pop->Ind);
        }
        else
        {
            tmpInd = (TInd *)st_sim_resize_array(
                Pop->Ind, (size_t)newSize, sizeof(TInd));

            if (NULL != tmpInd)
                Pop->Ind = tmpInd;
        }

        Pop->NoIndsAlloc = newSize;
    }

    info(1, "Population %s resized to %ld individuals\n", *PopName, newSize);
    return;

fail:

    while (ind > Pop->NoIndsAlloc)
    {
        ind--;
        remove_individuum(&(Pop->Ind[ind]));
    }

    if (0 == Pop->NoIndsAlloc)
    {
        st_sim_free(Pop->Ind);
    }
    else
    {
        tmpInd = (TInd *)st_sim_resize_array(
            Pop->Ind, (size_t)Pop->NoIndsAlloc, sizeof(TInd));

        if (NULL != tmpInd)
            Pop->Ind = tmpInd;
    }

    info(-2, "Memory allocation failed in int_resize_population\n");
    return;
}

void resize_population(char **PopName, char **cnewSize)
{
    POPULATION newSize = atol(*cnewSize);
    POPULATION oldSize;
    TPopListEl *PopListEl;

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        int_resize_population(PopName, newSize);
        return;
    }

    oldSize = PopListEl->Pop->NoIndsAlloc;

    int_resize_population(PopName, newSize);

    if ((oldSize != newSize) && (PopListEl->Pop->NoIndsAlloc == newSize))
    {
        st_clear_population_genotype(PopListEl->Pop);
        st_clear_population_gvalue(PopListEl->Pop);
    }
}

TPop *new_population(char **Name, POPULATION *NoInds)
{
    POPULATION ind = 0;
    TPop *NewPop = NULL;
    TPopListEl *NewPopListEl = NULL;

    if (2 > strlen(*Name))
    {
        info(-2, "Population name is too short\n");
        return NULL;
    }

    if (strlen(*Name) >= MAXNAME)
    {
        info(-2, "Population name is too long\n");
        return NULL;
    }

    if (0 > *NoInds)
    {
        info(-2, "Incorrect number of individuals\n");
        return NULL;
    }

    find_list_element(*Name, &NewPopListEl);

    if (NewPopListEl != NULL)
        NewPop = NewPopListEl->Pop;

    if (NULL != NewPop)
    {
        info(-2, "Population already exists\n");
        return NULL;
    }

    NewPop = (TPop *)st_sim_alloc(1, TPop);

    if (NULL == NewPop)
        goto fail_alloc;

    NewPopListEl = (TPopListEl *)st_sim_alloc(1, TPopListEl);

    if (NULL == NewPopListEl)
        goto fail_alloc;

    strcpy(NewPop->Name, *Name);

    NewPop->NoIndsAlloc = *NoInds;
    NewPop->NoInds = 0;
    NewPopListEl->Pop = NewPop;
    NewPopListEl->Prev = NULL;
    NewPopListEl->Next = PopListStart;

    if (NewPopListEl->Next != NULL)
        NewPopListEl->Next->Prev = NewPopListEl;

    PopListStart = NewPopListEl;

    if (NewPop->NoIndsAlloc > 0)
    {
        NewPop->Ind = (TInd *)st_sim_alloc(NewPop->NoIndsAlloc, TInd);

        if (NULL == NewPop->Ind)
            goto fail;

        for (ind = 0; ind < NewPop->NoIndsAlloc; ind++)
        {
            if (0 != new_individuum(&(NewPop->Ind[ind])))
                goto fail;
        }
    }
    else
    {
        NewPop->Ind = NULL;
    }

    info(1, "Population %s with %ld individuals generated\n", *Name, *NoInds);
    return (NewPop);

fail_alloc:

    st_sim_free(NewPopListEl);
    st_sim_free(NewPop);

    info(-2, "Memory allocation failed in new_population\n");
    return NULL;

fail:

    while (ind > 0)
    {
        ind--;
        remove_individuum(&(NewPop->Ind[ind]));
    }

    st_sim_free(NewPop->Ind);

    PopListStart = NewPopListEl->Next;

    if (NULL != PopListStart)
        PopListStart->Prev = NULL;

    st_sim_free(NewPop);
    st_sim_free(NewPopListEl);

    info(-2, "Memory allocation failed in new_population\n");
    return NULL;
}

static void st_remove_population_entry(TPopListEl *entry)
{
    for (POPULATION ind = 0; ind < entry->Pop->NoIndsAlloc; ind++)
        remove_individuum(&(entry->Pop->Ind[ind]));

    if (entry->Prev)
        entry->Prev->Next = entry->Next;
    else
        PopListStart = entry->Next;
    if (entry->Next)
        entry->Next->Prev = entry->Prev;

    if (entry->Pop->Ind)
        st_sim_free(entry->Pop->Ind);
    st_sim_free(entry->Pop);
    st_sim_free(entry);
}

void st_simulation_unload(void)
{
    while (PopListStart)
        st_remove_population_entry(PopListStart);
    st_sim_pool_dispose();
}

void remove_population(char **PopNames)
{
    char *token;
    TPopListEl *DelPopListEl;
    char *_PopNames;

    _PopNames = (char *)yCalloc(strlen(*PopNames) + 2, char);

    if (NULL == _PopNames)
    {
        info(-2, "Memory allocation failed in remove_population\n");
        return;
    }

    strcpy(_PopNames, *PopNames);
    token = strtok(_PopNames, "[' ','\n']");

    while (NULL != token)
    {
        find_list_element(token, &DelPopListEl);
        if (DelPopListEl != NULL)
        {
            st_remove_population_entry(DelPopListEl);
            info(1, "Population %s deleted\n", token);
        }

        token = strtok(NULL, "[' ','\n']");
    }
    yFree(_PopNames);
}

TPop *exist_population(char *PopName, POPULATION NoInds, int empty)
{
    TPopListEl *PopListEl;

    find_list_element(PopName, &PopListEl);

    if (PopListEl == NULL)
    {
        return new_population(&PopName, &NoInds);
    }

    if (PopListEl->Pop->NoIndsAlloc != NoInds)
    {
        int_resize_population(&PopName, NoInds);

        if (PopListEl->Pop->NoIndsAlloc != NoInds)
            return NULL;
    }

    if (1 == empty)
    {
        remove_genotype_population(&PopName);
        remove_evaluate_population(&PopName);
    }

    return PopListEl->Pop;
}

int insertLocus(unsigned short All, CPOS Pos, THom *Hom)
{
    if ((0 == Hom->NoLociUsed) || ((Hom->Loc[Hom->NoLociUsed - 1].Pos < Pos) &&
                                   (Hom->Loc[Hom->NoLociUsed - 1].All != All)))
    {
        TLoc *Loc;

        if (Hom->NoLociUsed >= Hom->NoLociAlloc &&
            st_sim_hom_reserve(Hom, (size_t)Hom->NoLociUsed + 1))
        {
            info(-2, "Cannot grow population segment storage\n");
            return -1;
        }
        Loc = &Hom->Loc[Hom->NoLociUsed];
        Loc->Pos = Pos;
        Loc->All = All;
        Hom->NoLociUsed++;
        return 0;
    }
    return 1;
}

void init_population(char **Name, char **cCount, int *_Ind, int *_Chrom,
                     int *_Hom, double *_Pos, int *_All)
{
    int min_chrom = 0, max_chrom = 0, min_hom = 0, max_hom = 0;
    POPULATION row, Count = atol(*cCount);
    POPULATION max_Ind = 0, min_Ind = 0;
    TPopListEl *PopListEl;
    TPop *Pop;
    THom *Hom;

    if (Count < 1)
    {
        info(-2, "Population data contain no loci\n");
        return;
    }

    find_list_element(*Name, &PopListEl);

    if (PopListEl != NULL)
    {
        info(-2, "First remove population %s\n", *Name);
        return;
    }

    for (row = 0; row < Count; row++)
    {
        max_chrom =
            (*(_Chrom + row) - 1 > max_chrom) ? *(_Chrom + row) - 1 : max_chrom;
        min_chrom =
            (*(_Chrom + row) - 1 < min_chrom) ? *(_Chrom + row) - 1 : min_chrom;
        max_hom = (*(_Hom + row) - 1 > max_hom) ? *(_Hom + row) - 1 : max_hom;
        min_hom = (*(_Hom + row) - 1 < min_hom) ? *(_Hom + row) - 1 : min_hom;
        max_Ind = (*(_Ind + row) - 1 > max_Ind) ? *(_Ind + row) - 1 : max_Ind;
        min_Ind = (*(_Ind + row) - 1 < min_Ind) ? *(_Ind + row) - 1 : min_Ind;
    }

    if ((min_chrom < 0) || (max_chrom >= NoChroms))
    {
        info(-2, "Error in the chromosome column\n");
        return;
    }

    if ((min_hom < 0) || (max_hom >= NoHoms))
    {
        info(-2, "Error in the homologue column\n");
        return;
    }

    if (min_Ind < 0)
    {
        info(-2, "Error in the individual column\n");
        return;
    }

    Pop = exist_population(*Name, max_Ind + 1, 1);
    if (NULL == Pop)
        return;

    Pop->NoInds = max_Ind + 1;
    uint64_t first;
    if (st_sim_reserve(Pop->NoInds, &first))
    {
        remove_population(Name);
        return;
    }
    for (POPULATION i = 0; i < Pop->NoInds; i++)
        if (st_sim_identity_new(&Pop->Ind[i], Pop->Name, first + (uint64_t)i,
                                NULL, NULL, NULL))
        {
            remove_population(Name);
            return;
        }

    for (row = 0; row < Count; row++)
    {
        TInd *Ind = get_individual(Pop, *(_Ind + row) - 1);
        Hom = &(Ind->Chrom[*(_Chrom + row) - 1].Hom[*(_Hom + row) - 1]);

        if (-1 == insertLocus(_All[row], _Pos[row], Hom))
        {
            remove_population(Name);
            return;
        }

        if (0 == _All[row])
            info(-1, "Don't use 0 to code alleles\n");
    }
}

void get_allele_number(char **Name, int *Count)
{
    int chrom, hom;
    POPULATION ind;
    TPopListEl *PopListEl;
    TPop *Pop;
    TInd *Ind;
    TChrom *Chrom;
    THom *Hom;

    find_list_element(*Name, &PopListEl);

    if (PopListEl != NULL)
    {
        *Count = 0;
        Pop = PopListEl->Pop;
        for (ind = 0; ind < Pop->NoInds; ind++)
        {
            Ind = get_individual(PopListEl->Pop, ind);
            for (chrom = 0; chrom < NoChroms; chrom++)
            {
                Chrom = &(Ind->Chrom[chrom]);
                for (hom = 0; hom < NoHoms; hom++)
                {
                    Hom = &(Chrom->Hom[hom]);
                    *Count += Hom->NoLociUsed;
                }
            }
        }
    }
}

void get_population(char **PopName, int *_Ind, int *_Chrom, int *_Hom,
                    double *_Pos, int *_All)
{
    int chrom, hom, loc;
    POPULATION ind;
    POPULATION row = 0;
    TPopListEl *PopListEl;
    TPop *Pop;
    TInd *Ind;
    TChrom *Chrom;
    THom *Hom;
    TLoc *Loc;

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined: %s\n", *PopName);
        return;
    }

    Pop = PopListEl->Pop;

    for (ind = 0; ind < Pop->NoInds; ind++)
    {
        Ind = get_individual(PopListEl->Pop, ind);
        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            Chrom = &(Ind->Chrom[chrom]);
            for (hom = 0; hom < NoHoms; hom++)
            {
                Hom = &(Chrom->Hom[hom]);
                for (loc = 0; loc < Hom->NoLociUsed; loc++)
                {
                    Loc = &(Hom->Loc[loc]);
                    _Ind[row] = ind + 1;
                    _Chrom[row] = chrom + 1;
                    _Hom[row] = hom + 1;
                    _Pos[row] = (double)Loc->Pos;
                    _All[row++] = Loc->All;
                }
            }
        }
    }
}

void list_populations_parameter(int *lPopName, int *NoPops)
{
    TPopListEl *PopListEl;
    PopListEl = PopListStart;

    while (PopListEl != NULL)
    {
        (*NoPops)++;
        (*lPopName) += strlen(PopListEl->Pop->Name) + 1;
        PopListEl = PopListEl->Next;
    }
    (*lPopName) += 10;
}

void list_populations(char **all_pops, int *inds)
{
    TPopListEl *PopListEl;
    PopListEl = PopListStart;

    strcpy(*all_pops, "");

    while (PopListEl != NULL)
    {
        *inds = PopListEl->Pop->NoInds;
        inds = inds + 1;
        strcat(*all_pops, PopListEl->Pop->Name);
        strcat(*all_pops, " ");
        PopListEl = PopListEl->Next;
    }
}

void set_co_freq(double *cofreq)
{
    if (0 <= *cofreq)
    {
        COFreq = *cofreq;
        info(1, "Crossoverfrequency set to %f\n", *cofreq);
    }
    else
        info(-1, "Recombination Frequency must be a positive number\n");
}

int compare_doubles(const void *a, const void *b)
{
    const double *da = (const double *)a;
    const double *db = (const double *)b;
    return (*da > *db) - (*da < *db);
}

static int meiosis_impl(THom *PgHom, THom *PHom, double Length,
                        double *significant, int r)
{
    unsigned int NuCO = gsl_ran_poisson(rng[r], COFreq * Length);
    size_t needed;
    unsigned int i, aHom;
    int aPos[] = {0, 0};
    CPOS pos;

    if (PgHom->NoLociUsed < 0 || PHom[0].NoLociUsed < 0 ||
        PHom[1].NoLociUsed < 0 || NuCO > (unsigned int)INT_MAX - 2)
        return -1;
    needed = (size_t)NuCO + 2;
    if ((size_t)PHom[0].NoLociUsed > (size_t)INT_MAX - needed)
        return -1;
    needed += (size_t)PHom[0].NoLociUsed;
    if ((size_t)PHom[1].NoLociUsed > (size_t)INT_MAX - needed)
        return -1;
    needed += (size_t)PHom[1].NoLociUsed;
    if ((size_t)PgHom->NoLociUsed > (size_t)INT_MAX - needed)
        return -1;
    needed += (size_t)PgHom->NoLociUsed;
    if (st_sim_hom_reserve(PgHom, needed))
        return -1;

    CPOS posCO[(size_t)NuCO + 2];

    if (-1 < *significant)
    {
        for (i = 1; i <= NuCO; i++)
        {
            posCO[i] =
                fround(gsl_rng_uniform_pos(rng[r]) * Length, *significant);
            posCO[i] = (*(posCO + i) > Length) ? Length : *(posCO + i);
        }
    }
    else
    {
        for (i = 1; i < NuCO + 1; i++)
        {
            posCO[i] = gsl_rng_uniform_pos(rng[r]) * Length;
        }
    }

    posCO[0] = 0;
    posCO[NuCO + 1] = Length;

    switch (NuCO)
    {

    case 0:

    case 1:
        break;

    case 2:
        if ((posCO[1] > posCO[2]))
        {
            pos = posCO[1];
            posCO[1] = posCO[2];
            posCO[2] = pos;
        }
        break;

    case 3:
        if (!((posCO[1] < posCO[2]) && (posCO[2] < posCO[3])))
            qsort(posCO + 1, NuCO, sizeof(double), compare_doubles);
        break;

    default:

        qsort(posCO + 1, NuCO, sizeof(double), compare_doubles);
    }

    aHom = (int)gsl_rng_uniform_int(rng[r], 2);

    for (i = 0; i <= NuCO; i++)
    {
        while ((aPos[aHom] + 1 < PHom[aHom].NoLociUsed) &&
               (PHom[aHom].Loc[aPos[aHom]].Pos < posCO[i]) &&
               (PHom[aHom].Loc[aPos[aHom] + 1].Pos <= posCO[i]))
        {
            aPos[aHom]++;
        }
        while ((aPos[aHom] < PHom[aHom].NoLociUsed) &&
               (PHom[aHom].Loc[aPos[aHom]].Pos < posCO[i + 1]))
        {
            pos = PHom[aHom].Loc[aPos[aHom]].Pos;
            pos = (pos > posCO[i]) ? pos : posCO[i];
            if (-1 == insertLocus(PHom[aHom].Loc[aPos[aHom]].All, pos, PgHom))
                return -1;

            aPos[aHom]++;
        }
        if (aPos[aHom] > 0)
            aPos[aHom]--;
        aHom = (aHom == 1) ? 0 : 1;
    }
    return 0;
}

void meiosis(THom *PgHom, THom *PHom, double Length, double *significant, int r)
{
    (void)meiosis_impl(PgHom, PHom, Length, significant, r);
}

int intcmp(const void *a, const void *b)
{
    const int *ia = (const int *)a;
    const int *ib = (const int *)b;
    return *ia - *ib;
}

int unsignedshortcmp(const void *a, const void *b)
{
    const unsigned short *ia = (const unsigned short *)a;
    const unsigned short *ib = (const unsigned short *)b;
    return *ia - *ib;
}

void cross(char **NamePg, char **NameP1, char **NameP2, double *Self,
           char **cNoPg, double *significant, int *mode, int *crossing_scheme)
{
    POPULATION ind, NoPg = atol(*cNoPg);

    TPopListEl *PopListEl;

    TPop *PgPop = NULL, *P1Pop = NULL, *P2Pop = NULL;
    int parallel_error = 0;

    if (2 != NoHoms)
    {
        info(-2, "Implemented for diploids only\n");
        return;
    }

    find_list_element(*NamePg, &PopListEl);
    if (PopListEl != NULL)
        PgPop = PopListEl->Pop;

    find_list_element(*NameP1, &PopListEl);
    if (PopListEl != NULL)
        P1Pop = PopListEl->Pop;

    find_list_element(*NameP2, &PopListEl);
    if (PopListEl != NULL)
        P2Pop = PopListEl->Pop;

    if (-1 == NoPg)
    {
        if ((NULL == PgPop) || (PgPop->NoIndsAlloc < 1))
        {
            info(-2, "Progeny population must exist when NoPg is -1\n");
            return;
        }

        NoPg = PgPop->NoIndsAlloc;
    }
    else if (1 > NoPg)
    {
        info(-2, "Number of progenies must be positive\n");
        return;
    }

    if ((P1Pop == NULL) || (P2Pop == NULL))
    {
        info(-2, "Parent population not defined\n");
        return;
    }

    if ((P1Pop == PgPop) || (P2Pop == PgPop))
    {
        info(-2, "Progeny population incorrect\n");
        return;
    }

    if ((0 == P1Pop->NoInds) || (0 == P2Pop->NoInds))
    {
        info(-2, "Parental population of size 0\n");
        return;
    }

    PgPop = exist_population(*NamePg, NoPg, 1);
    if (NULL == PgPop)
        return;

    PgPop->NoInds = NoPg;
    uint64_t first;
    if (st_sim_reserve(NoPg, &first))
    {
        remove_population(NamePg);
        return;
    }

    st_parallel_info_error_clear();

#ifdef MT
#pragma omp parallel reduction(| : parallel_error)
#endif

    {
        int chrom;
        double r;
        POPULATION ind1, ind2;
        TInd *P1Ind, *P2Ind, *PgInd;
        int local_error = 0;

        int i = omp_get_thread_num();

#ifdef MT
        omp_set_schedule(omp_sched_static, 0);
#pragma omp for schedule(runtime)
#endif

        for (ind = 0; ind < PgPop->NoInds; ind++)
        {

            if ((0 == *crossing_scheme) || (1 == *crossing_scheme))
            {
                P1Ind = get_individual(
                    P1Pop, ind1 = (POPULATION)gsl_rng_uniform_int(
                               rng[i], (unsigned long int)P1Pop->NoInds));
            }
            else
            {
                P1Ind = get_individual(P1Pop, ind1 = ind % P1Pop->NoInds);
            }

            if ((0 == *crossing_scheme) || (2 == *crossing_scheme))
            {
                P2Ind = get_individual(
                    P2Pop, ind2 = (POPULATION)gsl_rng_uniform_int(
                               rng[i], (unsigned long int)P2Pop->NoInds));
            }
            else
            {
                P2Ind = get_individual(P2Pop, ind2 = ind % P2Pop->NoInds);
            }

            if ((0 == *mode) && (P1Ind == P2Ind) && (1 < P1Pop->NoInds))
            {
                while (P1Ind == P2Ind)
                {
                    P2Ind = get_individual(
                        P2Pop, ind2 = (POPULATION)gsl_rng_uniform_int(
                                   rng[i], (unsigned long int)P2Pop->NoInds));
                }
            }

            r = gsl_rng_uniform_pos(rng[i]);
            PgInd = get_individual(PgPop, ind);

            if (*Self / 2 > r)
            {
                P2Ind = P1Ind;
                info(2, "  %s (%ld) selfed -> %s (%ld)\n", *NameP1, ind1,
                     *NamePg, ind);
            }
            else
            {
                if (*Self > r)
                {
                    P1Ind = P2Ind;
                    info(2, "  %s (%ld) selfed -> %s (%ld)\n", *NameP2, ind2,
                         *NamePg, ind);
                }
                else
                {
                    info(2, "  %s (%ld) x %s (%ld) -> %s (%ld)\n", *NameP1,
                         ind1, *NameP2, ind2, *NamePg, ind);
                }
            }

            if (st_sim_prepare_individual(PgInd, P1Ind, P2Ind, 1) ||
                st_sim_identity_new(PgInd, PgPop->Name, first + (uint64_t)ind,
                                    P1Ind, P2Ind, NULL))
            {
                parallel_error = 1;
                continue;
            }

            if (2 == NoHoms)
            {
                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    PgInd->Chrom[chrom].Hom[0].NoLociUsed = 0;
                    PgInd->Chrom[chrom].Hom[1].NoLociUsed = 0;

                    if (-1 == meiosis_impl(PgInd->Chrom[chrom].Hom,
                                           (*P1Ind).Chrom[chrom].Hom,
                                           PopDescr[chrom].ChromLen,
                                           significant, i))
                        local_error = 1;

                    if (!local_error &&
                        (-1 == meiosis_impl(&(PgInd->Chrom[chrom].Hom[1]),
                                            (*P2Ind).Chrom[chrom].Hom,
                                            PopDescr[chrom].ChromLen,
                                            significant, i)))
                        local_error = 1;

                    if (local_error)
                    {
                        parallel_error = 1;
                        break;
                    }
                }
            }
        }
    }

    if (parallel_error || st_parallel_info_error_get())
    {
        st_parallel_info_error_clear();
        remove_population(NamePg);
        info(-2, "Error during parallel crossing\n");
        return;
    }

    st_parallel_info_error_clear();
}

void ps_set_num_threads(int *activate)
{
    if (*activate > MAXTHREADS)
        *activate = MAXTHREADS;
    if (*activate < 1)
        *activate = 1;
    omp_set_num_threads(*activate);
}

void ssd_mating(char **NamePg, char **NameP1, char **cNoPg, int *maxcycles,
                double *significant)
{
    POPULATION ind, NoPg = atol(*cNoPg), nPg = atol(*cNoPg), prog;
    int chrom;

    TPopListEl *PopListEl;
    TPop *PgPop = NULL, *P1Pop = NULL;
    TInd *PgInd, *P1Ind;
    THom *tmpHomA = NULL, *tmpHomB;
    int cycle, nocycles = *maxcycles, hom, loc;
    int destination_prepared = 0;

    if (2 != NoHoms)
    {
        info(-2, "Implemented only for diploids\n");
        return;
    }

    if (1 > NoPg)
    {
        info(-2, "Number of progenies must be a positive number\n");
        return;
    }

    if (1 > nocycles)
    {
        info(-2, "Number of mating cycles must be positive\n");
        return;
    }

    find_list_element(*NamePg, &PopListEl);

    if (PopListEl != NULL)
        PgPop = PopListEl->Pop;

    find_list_element(*NameP1, &PopListEl);

    if (PopListEl != NULL)
        P1Pop = PopListEl->Pop;

    if (P1Pop == NULL)
    {
        info(-2, "Population not defined: %s\n", *NameP1);
        return;
    }

    if (P1Pop == PgPop)
    {
        info(-2,
             "Progeny population name (%s) is identical to parent population "
             "name (%s)\n",
             *NamePg, *NameP1);
        return;
    }

    if (0 == P1Pop->NoInds)
    {
        info(-2, "Parental population of size 0\n");
        return;
    }

    if (P1Pop->NoInds > LONG_MAX / NoPg)
    {
        info(-2, "Number of progenies is too large\n");
        return;
    }

    NoPg = P1Pop->NoInds * NoPg;

    tmpHomA = (THom *)malloc((size_t)NoHoms * sizeof(THom));

    if (NULL == tmpHomA)
    {
        info(-2, "Memory allocation failed in ssd_mating\n");
        return;
    }

    for (hom = 0; hom < NoHoms; hom++)
    {
        tmpHomA[hom].Loc = NULL;
        tmpHomA[hom].Owner = NULL;
        tmpHomA[hom].Offset = 0;
    }

    for (hom = 0; hom < NoHoms; hom++)
    {
        tmpHomA[hom].NoLociAlloc = NoLociInit;
        tmpHomA[hom].NoLociUsed = 0;
        tmpHomA[hom].Loc = (TLoc *)malloc((size_t)NoLociInit * sizeof(TLoc));

        if (NULL == tmpHomA[hom].Loc)
        {
            info(-2, "Memory allocation failed in ssd_mating\n");
            goto cleanup;
        }
    }

    PgPop = exist_population(*NamePg, NoPg, 1);
    if (NULL == PgPop)
        goto cleanup;

    destination_prepared = 1;
    PgPop->NoInds = NoPg;
    uint64_t first;
    if (st_sim_reserve(NoPg, &first))
        goto fail_population;

    for (ind = 0; ind < P1Pop->NoInds; ind++)
    {
        P1Ind = get_individual(P1Pop, ind);
        for (prog = 0; prog < nPg; prog++)
        {
            PgInd = get_individual(PgPop, ind * nPg + prog);
            if (st_sim_prepare_individual(PgInd, P1Ind, P1Ind, 1))
                goto fail_population;
            if (st_sim_identity_new(PgInd, PgPop->Name,
                                    first + (uint64_t)(ind * nPg + prog), P1Ind,
                                    NULL, "SSD"))
                goto fail_population;
            for (chrom = 0; chrom < NoChroms; chrom++)
            {
                int diff = 0;
                cycle = nocycles;
                tmpHomB = P1Ind->Chrom[chrom].Hom;

                do
                {
                    diff = 0;

                    for (hom = 0; hom < NoHoms; hom++)
                    {
                        if (tmpHomA[hom].NoLociAlloc < tmpHomB[hom].NoLociUsed)
                        {
                            int newNoLociAlloc = tmpHomB[hom].NoLociUsed;
                            TLoc *tmpLoc;

                            tmpLoc = (TLoc *)realloc(tmpHomA[hom].Loc,
                                                     (size_t)newNoLociAlloc *
                                                         sizeof(TLoc));

                            if (NULL == tmpLoc)
                            {
                                info(
                                    -2,
                                    "Memory allocation failed in ssd_mating\n");
                                goto fail_population;
                            }

                            tmpHomA[hom].Loc = tmpLoc;
                            tmpHomA[hom].NoLociAlloc = newNoLociAlloc;
                        }

                        tmpHomA[hom].NoLociUsed = tmpHomB[hom].NoLociUsed;

                        for (loc = 0; loc < tmpHomA[hom].NoLociUsed; loc++)
                        {
                            tmpHomA[hom].Loc[loc].All =
                                tmpHomB[hom].Loc[loc].All;
                            tmpHomA[hom].Loc[loc].Pos =
                                tmpHomB[hom].Loc[loc].Pos;
                        }
                    }

                    PgInd->Chrom[chrom].Hom[0].NoLociUsed = 0;
                    PgInd->Chrom[chrom].Hom[1].NoLociUsed = 0;

                    if (-1 == meiosis_impl(PgInd->Chrom[chrom].Hom, tmpHomA,
                                           PopDescr[chrom].ChromLen,
                                           significant, 0))
                        goto fail_population;

                    if (-1 == meiosis_impl(&(PgInd->Chrom[chrom].Hom[1]),
                                           tmpHomA, PopDescr[chrom].ChromLen,
                                           significant, 0))
                        goto fail_population;

                    if (PgInd->Chrom[chrom].Hom[0].NoLociUsed !=
                        PgInd->Chrom[chrom].Hom[1].NoLociUsed)
                    {
                        diff++;
                    }
                    else
                    {
                        for (loc = 0;
                             loc < PgInd->Chrom[chrom].Hom[0].NoLociUsed; loc++)
                        {
                            if ((PgInd->Chrom[chrom].Hom[0].Loc[loc].All !=
                                 PgInd->Chrom[chrom].Hom[1].Loc[loc].All) ||
                                (PgInd->Chrom[chrom].Hom[0].Loc[loc].Pos !=
                                 PgInd->Chrom[chrom].Hom[1].Loc[loc].Pos))
                            {
                                diff++;
                                break;
                            }
                        }
                    }

                    tmpHomB = PgInd->Chrom[chrom].Hom;
                } while ((0 != --cycle) && (diff));
            }
        }
    }

    destination_prepared = 0;
    goto cleanup;

fail_population:

    if (destination_prepared)
    {
        remove_population(NamePg);
        PgPop = NULL;
        destination_prepared = 0;
    }

cleanup:

    if (NULL != tmpHomA)
    {
        for (hom = 0; hom < NoHoms; hom++)
            free(tmpHomA[hom].Loc);

        free(tmpHomA);
    }
}

void dh(char **NamePg, char **NameP1, char **cNoPg, double *significant)
{
    POPULATION ind, NoPg = atol(*cNoPg);

    int chrom;

    TPopListEl *PopListEl;

    TPop *PgPop = NULL, *P1Pop = NULL;

    TInd *PgInd, *P1Ind;
    int destination_prepared = 0;

    if (2 != NoHoms)
    {
        info(-2, "Implemented for diploids only\n");
        return;
    }

    find_list_element(*NamePg, &PopListEl);
    if (PopListEl != NULL)
        PgPop = PopListEl->Pop;

    find_list_element(*NameP1, &PopListEl);
    if (PopListEl != NULL)
        P1Pop = PopListEl->Pop;

    if (-1 == NoPg)
    {
        if ((NULL == PgPop) || (PgPop->NoIndsAlloc < 1))
        {
            info(-2, "Progeny population must exist when NoPg is -1\n");
            return;
        }

        NoPg = PgPop->NoIndsAlloc;
    }
    else if (1 > NoPg)
    {
        info(-2, "Number of progenies must be positive\n");
        return;
    }

    if (P1Pop == NULL)
    {
        info(-2, "Parent population not defined\n");
        return;
    }

    if (P1Pop == PgPop)
    {
        info(-2, "Progeny population incorrect\n");
        return;
    }

    if (0 == P1Pop->NoInds)
    {
        info(-2, "Parental population of size 0\n");
        return;
    }

    PgPop = exist_population(*NamePg, NoPg, 1);
    if (NULL == PgPop)
        return;

    destination_prepared = 1;
    PgPop->NoInds = NoPg;
    uint64_t first;
    if (st_sim_reserve(NoPg, &first))
        goto fail_population;

    for (ind = 0; ind < PgPop->NoInds; ind++)
    {
        P1Ind = get_individual(P1Pop,
                               (POPULATION)gsl_rng_uniform_int(
                                   rng[0], (unsigned long int)P1Pop->NoInds));
        PgInd = get_individual(PgPop, ind);
        if (st_sim_prepare_individual(PgInd, P1Ind, P1Ind, 1) ||
            st_sim_identity_new(PgInd, PgPop->Name, first + (uint64_t)ind,
                                P1Ind, NULL, "DH"))
            goto fail_population;

        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            int loc, hom;

            PgInd->Chrom[chrom].Hom[0].NoLociUsed = 0;
            PgInd->Chrom[chrom].Hom[1].NoLociUsed = 0;

            if (-1 == meiosis_impl(PgInd->Chrom[chrom].Hom,
                                   P1Ind->Chrom[chrom].Hom,
                                   PopDescr[chrom].ChromLen, significant, 0))
                goto fail_population;

            for (hom = 1; hom < NoHoms; hom++)
            {
                THom *P1Hom = &(PgInd->Chrom[chrom]).Hom[hom];
                THom *P2Hom = &(PgInd->Chrom[chrom]).Hom[0];
                if (st_sim_hom_reserve(P1Hom, (size_t)P2Hom->NoLociUsed))
                    goto fail_population;
                PgInd->Chrom[chrom].Hom[hom].NoLociUsed = 0;

                for (loc = 0; loc < P2Hom->NoLociUsed; loc++)
                {
                    if (-1 == insertLocus(P2Hom->Loc[loc].All,
                                          P2Hom->Loc[loc].Pos, P1Hom))
                        goto fail_population;
                }
            }
        }
    }

    return;

fail_population:

    if (destination_prepared)
        remove_population(NamePg);

    return;
}

static int parse_int_token(const char *token, int *value)
{
    char *end;
    long x;

    if ((NULL == token) || ('\0' == *token))
        return 0;

    errno = 0;
    x = strtol(token, &end, 10);

    if ((ERANGE == errno) || (end == token) || ('\0' != *end) ||
        (x < INT_MIN) || (x > INT_MAX))
        return 0;

    *value = (int)x;

    return 1;
}

static int parse_double_token(const char *token, CPOS *value)
{
    char *end;
    double x;

    if ((NULL == token) || ('\0' == *token))
        return 0;

    errno = 0;
    x = strtod(token, &end);

    if ((ERANGE == errno) || (end == token) || ('\0' != *end) || (!isfinite(x)))
        return 0;

    *value = (CPOS)x;

    return 1;
}

static int parse_chrom_pos_token(const char *token, int *chrom, CPOS *pos)
{
    char *endChrom;
    char *endPos;
    long c;
    double p;

    if ((NULL == token) || ('\0' == *token))
        return 0;

    errno = 0;
    c = strtol(token, &endChrom, 10);

    if ((ERANGE == errno) || (endChrom == token) || ('/' != *endChrom) ||
        (c < INT_MIN) || (c > INT_MAX))
        return 0;

    errno = 0;
    p = strtod(endChrom + 1, &endPos);

    if ((ERANGE == errno) || (endPos == endChrom + 1) || ('\0' != *endPos) ||
        (!isfinite(p)))
        return 0;

    *chrom = (int)c;
    *pos = (CPOS)p;

    return 1;
}

static void free_mappoints(TMappoint *points, int n)
{
    int i;

    if (NULL != points)
    {
        for (i = 0; i < n; i++)
        {
            free(points[i].name);
            free(points[i].class);
        }

        free(points);
    }
}

static int append_mappoint(TMappoint **points, int *n, int chrom, CPOS pos,
                           const char *name, const char *class)
{
    TMappoint *tmp;
    TMappoint *point;

    if (*n == INT_MAX)
        return 0;

    if ((size_t)(*n + 1) > SIZE_MAX / sizeof(TMappoint))
        return 0;

    tmp = (TMappoint *)realloc(*points, (size_t)(*n + 1) * sizeof(TMappoint));

    if (NULL == tmp)
        return 0;

    *points = tmp;

    point = &(tmp[*n]);

    point->chrom = chrom;
    point->Pos = pos;
    point->counter = *n;
    point->name = NULL;
    point->class = NULL;

    point->name = (char *)malloc(strlen(name) + 1);

    if (NULL == point->name)
        return 0;

    point->class = (char *)malloc(strlen(class) + 1);

    if (NULL == point->class)
    {
        free(point->name);
        point->name = NULL;
        return 0;
    }

    strcpy(point->name, name);
    strcpy(point->class, class);

    (*n)++;

    return 1;
}

static int mnamecmp(const void *i, const void *j)
{
    const TMappoint *m0 = (const TMappoint *)i;
    const TMappoint *m1 = (const TMappoint *)j;

    return strcmp(m0->name, m1->name);
}

int mccmp(const void *i, const void *j)
{
    const TMappoint *m0 = (const TMappoint *)i;
    const TMappoint *m1 = (const TMappoint *)j;

    if (m0->chrom > m1->chrom)
        return 1;
    if (m0->chrom < m1->chrom)
        return -1;

    if (m0->Pos > m1->Pos)
        return 1;
    if (m0->Pos < m1->Pos)
        return -1;

    if (m0->counter > m1->counter)
        return 1;
    if (m0->counter < m1->counter)
        return -1;

    return 0;
}

static int validate_map_genome(void)
{
    int chrom;

    if ((NULL == PopDescr) || (NoChroms < 1) || (NoHoms < 1))
    {
        info(-2, "No valid genome defined\n");
        return 0;
    }

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        if (!isfinite(PopDescr[chrom].ChromLen) ||
            (PopDescr[chrom].ChromLen <= 0))
        {
            info(-2, "Invalid chromosome length\n");
            return 0;
        }
    }

    return 1;
}

static int disperse_mappoints_impl(TMappoint *points, int n, double factor,
                                   int *modified, int strict)
{
    int start = 0;
    int i;

    *modified = 0;

    if (!isfinite(factor) || (strict ? factor < 1 : factor <= 0))
    {
        info(-2, "disperse.factor must be finite and >= 1\n");

        return 0;
    }

    while (start < n)
    {
        int end = start;
        int size;
        int chrom;
        int j;
        int backward = 0;
        CPOS pos;
        CPOS distance;
        CPOS delta;

        while ((end + 1 < n) &&
               (points[end + 1].chrom == points[start].chrom) &&
               (points[end + 1].Pos == points[start].Pos))
            end++;

        size = end - start + 1;

        if (size > 1)
        {
            chrom = points[start].chrom;
            pos = points[start].Pos;

            if ((end + 1 < n) && (points[end + 1].chrom == chrom))
            {
                distance = points[end + 1].Pos - pos;
            }

            else if (pos == PopDescr[chrom - 1].ChromLen)
            {
                backward = 1;

                if ((start > 0) && (points[start - 1].chrom == chrom))
                {
                    distance = pos - points[start - 1].Pos;
                }

                else
                {
                    distance = pos;
                }
            }

            else
            {
                distance = PopDescr[chrom - 1].ChromLen - pos;
            }

            delta = distance / (factor * (double)size);

            if (strict && (!(delta > 0) || !isfinite(delta)))
            {
                info(-2,
                     "Unable to disperse loci at chromosome "
                     "%d position %f\n",
                     chrom, pos);

                return 0;
            }

            if (backward)
            {
                for (j = 0; j < size; j++)
                {
                    points[start + j].Pos -= (size - 1 - j) * delta;
                }
            }

            else
            {
                for (j = 0; j < size; j++)
                {
                    points[start + j].Pos += j * delta;
                }
            }

            *modified += size - 1;
        }

        start = end + 1;
    }

    for (i = 0; strict && i < n; i++)
    {
        int chrom = points[i].chrom;

        if (!isfinite(points[i].Pos) || (points[i].Pos < 0) ||
            (points[i].Pos > PopDescr[chrom - 1].ChromLen))
        {
            info(-2, "Invalid dispersed locus position on chromosome %d\n",
                 chrom);

            return 0;
        }

        if ((i > 0) && (points[i - 1].chrom == chrom) &&
            !(points[i].Pos > points[i - 1].Pos))
        {
            info(-2,
                 "Unable to obtain distinct locus positions on "
                 "chromosome %d\n",
                 chrom);

            return 0;
        }
    }

    return 1;
}

static int disperse_mappoints(TMappoint *points, int n, double factor,
                              int *modified)
{
    return disperse_mappoints_impl(points, n, factor, modified, 1);
}

static int prepare_mappoints(TMappoint *points, int n, int disperse,
                             double factor, int *modified)
{
    int i;

    *modified = 0;

    if (!validate_map_genome())
        return 0;

    if ((NULL == points) || (n < 1))
    {
        info(-2, "Linkage map contains no loci\n");
        return 0;
    }

    if (n > INT_MAX / NoHoms)
    {
        info(-2, "Linkage map contains too many loci\n");
        return 0;
    }

    if ((0 != disperse) && (1 != disperse))
    {
        info(-2, "Invalid disperse argument\n");
        return 0;
    }

    for (i = 0; i < n; i++)
    {
        int chrom = points[i].chrom;

        if ((chrom < 1) || (chrom > NoChroms))
        {
            info(-2, "Wrong chromosome (%d)\n", chrom);

            return 0;
        }

        if (!isfinite(points[i].Pos) || (points[i].Pos < 0) ||
            (points[i].Pos > PopDescr[chrom - 1].ChromLen))
        {
            info(-2,
                 "Chromosome %d: Locusposition outside of the chromosome "
                 "(0 <= %f <= %f)\n",
                 chrom, points[i].Pos, PopDescr[chrom - 1].ChromLen);

            return 0;
        }

        if ((NULL == points[i].name) || (NULL == points[i].class) ||
            ('\0' == points[i].name[0]) || ('\0' == points[i].class[0]))
        {
            info(-2, "Empty locus name or class\n");
            return 0;
        }

        if ((strlen(points[i].name) >= MAXNAME) ||
            (strlen(points[i].class) >= MAXNAME))
        {
            info(-2, "Locus name or class is too long\n");
            return 0;
        }

        if ((NULL != strpbrk(points[i].name, " \t\r\n")) ||
            (NULL != strpbrk(points[i].class, " \t\r\n")))
        {
            info(-2, "Whitespace not allowed in locus names or classes\n");

            return 0;
        }
    }

    qsort(points, n, sizeof(TMappoint), mnamecmp);

    for (i = 1; i < n; i++)
    {
        if (0 == strcmp(points[i - 1].name, points[i].name))
        {
            info(-2, "Duplicate locus name: %s\n", points[i].name);

            return 0;
        }
    }

    qsort(points, n, sizeof(TMappoint), mccmp);

    if (disperse)
    {
        if (!disperse_mappoints(points, n, factor, modified))
            return 0;
    }

    return 1;
}

static int read_mappoint_file(char *fName, TMappoint **points, int *nMP)
{
    FILE *fp = NULL;
    TMappoint *newMPs = NULL;
    char line[4 * MAXNAME + 128];
    char *cChrom;
    char *cPos;
    char *lociname;
    char *lociclass;
    char *extra;
    CPOS pos;
    int lineNo = 0;
    int chrom;

    *points = NULL;
    *nMP = 0;

    if ((NULL == fName) || ('\0' == *fName))
    {
        info(-2, "Invalid linkage map file\n");
        return 0;
    }

    fp = fileopen(fName, "r");

    if (NULL == fp)
        goto fail;

    while (NULL != fgets(line, sizeof(line), fp))
    {
        int ch;

        lineNo++;

        if (NULL == strchr(line, '\n'))
        {
            ch = fgetc(fp);

            if (EOF != ch)
            {
                while (('\n' != ch) && (EOF != ch))
                    ch = fgetc(fp);

                info(-2, "Map file line %d is too long\n", lineNo);

                goto fail;
            }
        }

        cChrom = strtok(line, " \t\r\n");
        cPos = strtok(NULL, " \t\r\n");
        lociname = strtok(NULL, " \t\r\n");
        lociclass = strtok(NULL, " \t\r\n");
        extra = strtok(NULL, " \t\r\n");

        if ((NULL == cChrom) || (NULL == cPos) || (NULL == lociname) ||
            (NULL == lociclass) || (NULL != extra))
        {
            info(-2, "Error in map file line %d\n", lineNo);

            goto fail;
        }

        if (!parse_int_token(cChrom, &chrom))
        {
            info(-2, "Invalid chromosome in map file line %d\n", lineNo);

            goto fail;
        }

        if (!parse_double_token(cPos, &pos))
        {
            info(-2, "Invalid locus position in map file line %d\n", lineNo);

            goto fail;
        }

        if (!append_mappoint(&newMPs, nMP, chrom, pos, lociname, lociclass))
            goto fail_alloc;
    }

    if (ferror(fp))
    {
        info(-2, "Error reading map file %s\n", fName);

        goto fail;
    }

    if (0 != fclose(fp))
    {
        fp = NULL;

        info(-2, "Can't close the file stream\n");

        goto fail;
    }

    fp = NULL;

    if (0 == *nMP)
    {
        info(-2, "Linkage map contains no loci\n");
        goto fail;
    }

    *points = newMPs;

    return 1;

fail_alloc:

    info(-2, "Memory allocation failed while reading linkage map\n");

fail:

    if (NULL != fp)
        fileclose(fp);

    free_mappoints(newMPs, *nMP);

    *nMP = 0;

    return 0;
}

static int write_mappoint_file(char *fName, const TMappoint *points, int n)
{
    FILE *fp;
    int i;

    fp = fileopen(fName, "w");

    if (NULL == fp)
        return 0;

    for (i = 0; i < n; i++)
    {
        if (0 > fprintf(fp, "%d\t%.17g\t%s\t%s\n", points[i].chrom,
                        points[i].Pos, points[i].name, points[i].class))
        {
            info(-2, "Error writing map file %s\n", fName);

            goto fail;
        }
    }

    if (0 != fclose(fp))
    {
        info(-2, "Can't close the file stream\n");

        remove(fName);

        return 0;
    }

    return 1;

fail:

    fclose(fp);
    remove(fName);

    return 0;
}

static int install_mappoints(TMappoint *points, int n)
{
    TMappoint **newMap = NULL;
    int *count = NULL;
    int chrom, i, j;
    int counter = 0;

    if ((NULL == PopDescr) || (NoChroms < 1) || (NoHoms < 1) || (n < 1))
        return 0;

    if (n > INT_MAX / NoHoms)
        return 0;

    if ((size_t)NoChroms > SIZE_MAX / sizeof(TMappoint *))
        return 0;

    newMap = (TMappoint **)malloc((size_t)NoChroms * sizeof(TMappoint *));

    if (NULL == newMap)
        return 0;

    for (chrom = 0; chrom < NoChroms; chrom++)
        newMap[chrom] = NULL;

    if ((size_t)NoChroms > SIZE_MAX / sizeof(int))
        goto fail;

    count = (int *)malloc((size_t)NoChroms * sizeof(int));

    if (NULL == count)
        goto fail;

    for (chrom = 0; chrom < NoChroms; chrom++)
        count[chrom] = 0;

    for (i = 0; i < n; i++)
    {
        if ((points[i].chrom < 1) || (points[i].chrom > NoChroms))
            goto fail;

        count[points[i].chrom - 1]++;
    }

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        if (0 < count[chrom])
        {
            if ((size_t)count[chrom] > SIZE_MAX / sizeof(TMappoint))
                goto fail;

            newMap[chrom] =
                (TMappoint *)malloc((size_t)count[chrom] * sizeof(TMappoint));

            if (NULL == newMap[chrom])
                goto fail;
        }
    }

    i = 0;

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        for (j = 0; j < count[chrom]; j++)
        {
            if ((i >= n) || (points[i].chrom != chrom + 1))
                goto fail;

            newMap[chrom][j] = points[i];
            newMap[chrom][j].counter = counter++;

            i++;
        }
    }

    if (i != n)
        goto fail;

    for (i = 0; i < n; i++)
    {
        points[i].name = NULL;
        points[i].class = NULL;
    }

    remove_map();

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        PopDescr[chrom].Mappoint = newMap[chrom];

        PopDescr[chrom].NrMappoints = count[chrom];
    }

    cMP = n;
    cMPhom = NoHoms * n;

    free(newMap);
    free(count);

    return 1;

fail:

    if (NULL != newMap)
    {
        for (chrom = 0; chrom < NoChroms; chrom++)
            free(newMap[chrom]);
    }

    free(newMap);
    free(count);

    return 0;
}

void linkage_map_load(char **fName, int *disperse, char **fNameDisperse,
                      double *factor, int *retval)
{
    TMappoint *points = NULL;
    int nMP = 0;
    int modified = 0;

    if (NULL == retval)
    {
        info(-2, "Invalid linkage map arguments\n");
        return;
    }

    *retval = -2;

    if ((NULL == fName) || (NULL == *fName) || (NULL == disperse) ||
        (NULL == factor))
    {
        info(-2, "Invalid linkage map arguments\n");
        return;
    }

    if (!validate_map_genome())
        return;

    if (!read_mappoint_file(*fName, &points, &nMP))
        return;

    if (!prepare_mappoints(points, nMP, *disperse, *factor, &modified))
        goto fail;

    if (*disperse && (NULL != fNameDisperse) && (NULL != *fNameDisperse) &&
        ('\0' != (*fNameDisperse)[0]))
    {
        if (!write_mappoint_file(*fNameDisperse, points, nMP))
            goto fail;
    }

    if (!install_mappoints(points, nMP))
    {
        info(-2, "Memory allocation failed in linkage_map_load\n");

        goto fail;
    }

    free_mappoints(points, nMP);

    if (*disperse)
    {
        if (0 < modified)
            info(0, "%d loci positions modified\n", modified);
    }

    if (*disperse && (NULL != fNameDisperse) && (NULL != *fNameDisperse) &&
        ('\0' != (*fNameDisperse)[0]))
    {
        info(0, "Dispersed linkage map saved to file: %s\n", *fNameDisperse);
    }

    if (!strncmp(*fName,"@SelectionTools:",strlen("@SelectionTools:")))
        info(0,"Linkage map loaded\n");
    else
        info(0, "Map file %s loaded\n", *fName);

    info(0, "%d loci defined\n", nMP);

    *retval = 0;

    return;

fail:

    free_mappoints(points, nMP);

    return;
}

void define_map(char **fName, int *retval)
{
    int disperse = 1;
    double factor = MAP_DISPERSE_FACTOR_DEFAULT;

    linkage_map_load(fName, &disperse, NULL, &factor, retval);
}

void define_linkage_map(int *dim, int *chrom, double *pos, char **names,
                        char **class)
{
    TMappoint *newMPs = NULL;
    int nMP = 0;
    int modified = 0;
    int i;

    if (!validate_map_genome())
        return;

    if ((NULL == dim) || (*dim < 1))
    {
        info(-2, "Linkage map contains no loci\n");
        return;
    }

    if ((NULL == chrom) || (NULL == pos) || (NULL == names) || (NULL == class))
    {
        info(-2, "Invalid linkage map data\n");
        return;
    }

    if (*dim > INT_MAX / NoHoms)
    {
        info(-2, "Linkage map contains too many loci\n");
        return;
    }

    for (i = 0; i < *dim; i++)
    {
        if ((NULL == names[i]) || (NULL == class[i]))
        {
            info(-2, "Invalid linkage map data\n");
            goto fail;
        }

        if (!append_mappoint(&newMPs, &nMP, chrom[i], pos[i], names[i],
                             class[i]))
            goto fail_alloc;
    }

    if (!prepare_mappoints(newMPs, nMP, 1, MAP_DISPERSE_FACTOR_DEFAULT,
                           &modified))
        goto fail;

    if (!install_mappoints(newMPs, nMP))
        goto fail_alloc;

    free_mappoints(newMPs, nMP);

    if (0 < modified)
        info(0, "%d loci positions modified\n", modified);

    info(0, "%d loci defined\n", nMP);

    return;

fail_alloc:

    info(-2, "Memory allocation failed in define_linkage_map\n");

fail:

    free_mappoints(newMPs, nMP);

    return;
}

void generate_map_file(char **fName, char **Description, int *retval)
{
    TMappoint *newMPs = NULL;
    char *description = NULL;
    char *token;
    char name[MAXNAME];
    char classname[MAXNAME];
    int status = 0;
    int number = 1;
    int cMPs = 0;
    int maxMP;
    int i;
    FILE *fp = NULL;
    int outputOpened = 0;

    *retval = -2;

    if ((NULL == PopDescr) || (NoChroms < 1) || (NoHoms < 1))
    {
        info(-2, "No valid genome defined\n");
        return;
    }

    for (i = 0; i < NoChroms; i++)
    {
        if (!isfinite(PopDescr[i].ChromLen) || (PopDescr[i].ChromLen <= 0))
        {
            info(-2, "Invalid chromosome length\n");
            return;
        }
    }

    maxMP = INT_MAX / NoHoms;

    description = (char *)malloc(strlen(*Description) + 1);

    if (NULL == description)
        goto fail_alloc;

    strcpy(description, *Description);

    token = strtok(description, "[]', \t\r\n");

    while (NULL != token)
    {
        switch (status)
        {
        case 0:

            if ('*' == token[0])
            {
                if (!parse_int_token(token + 1, &number) || (number < 1))
                {
                    info(-2, "Invalid number of loci: %s\n", token);

                    goto fail;
                }

                status = 1;
            }

            else
            {
                if (strlen(token) >= MAXNAME)
                {
                    info(-2, "Locus name too long: %s\n", token);

                    goto fail;
                }

                strcpy(name, token);
                number = 1;
                status = 2;
            }

            break;

        case 1:

            if (strlen(token) >= MAXNAME)
            {
                info(-2, "Locus name too long: %s\n", token);

                goto fail;
            }

            strcpy(name, token);
            status = 2;
            break;

        case 2:

            if (strlen(token) >= MAXNAME)
            {
                info(-2, "Locus class too long: %s\n", token);

                goto fail;
            }

            strcpy(classname, token);
            status = 3;
            break;

        case 3:

            if (number > maxMP - cMPs)
            {
                info(-2, "Linkage map contains too many loci\n");

                goto fail;
            }

            if (0 == strcmp(token, "random"))
            {
                double genomelength = 0;
                int lNR = 0;

                for (i = 0; i < NoChroms; i++)
                    genomelength += PopDescr[i].ChromLen;

                if (!isfinite(genomelength) || (genomelength <= 0))
                {
                    info(-2, "Invalid genome length\n");
                    goto fail;
                }

                for (i = 0; i < number; i++)
                {
                    char tmpstring[MAXNAME];
                    int nchars;
                    int chrom = -1;
                    CPOS pos = -gsl_rng_uniform_pos(rng[0]) * genomelength;

                    while ((pos < 0) && (chrom + 1 < NoChroms))
                        pos += PopDescr[++chrom].ChromLen;

                    if ((chrom < 0) || (chrom >= NoChroms) || (pos < 0) ||
                        (pos > PopDescr[chrom].ChromLen))
                    {
                        info(-2, "Error generating random map position\n");

                        goto fail;
                    }

                    nchars = snprintf(tmpstring, sizeof(tmpstring), "%s%d",
                                      name, lNR++);

                    if ((nchars < 0) || (nchars >= MAXNAME))
                    {
                        info(-2, "Generated locus name is too long\n");

                        goto fail;
                    }

                    if (!append_mappoint(&newMPs, &cMPs, chrom + 1, pos,
                                         tmpstring, classname))
                        goto fail_alloc;
                }

                number = 1;
                status = 0;
            }

            else if (0 == strcmp(token, "terminal"))
            {
                int chrom, j;
                int lNR = 0;

                if (0 != number % NoChroms)
                {
                    info(-2, "Number of terminal loci must be a multiple "
                             "of the number of chromosomes\n");

                    goto fail;
                }

                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    int nPerChrom = number / NoChroms;

                    CPOS delta = (nPerChrom > 1) ? PopDescr[chrom].ChromLen /
                                                       (nPerChrom - 1)
                                                 : 0;

                    for (j = 0; j < nPerChrom; j++)
                    {
                        char tmpstring[MAXNAME];
                        int nchars;
                        CPOS pos = delta * j;

                        nchars = snprintf(tmpstring, sizeof(tmpstring), "%s%d",
                                          name, lNR++);

                        if ((nchars < 0) || (nchars >= MAXNAME))
                        {
                            info(-2, "Generated locus name is too long\n");

                            goto fail;
                        }

                        if (!append_mappoint(&newMPs, &cMPs, chrom + 1, pos,
                                             tmpstring, classname))
                            goto fail_alloc;
                    }
                }

                number = 1;
                status = 0;
            }

            else if (0 == strcmp(token, "nonterminal"))
            {
                int chrom, j;
                int lNR = 0;

                if (0 != number % NoChroms)
                {
                    info(-2, "Number of nonterminal loci must be a multiple "
                             "of the number of chromosomes\n");

                    goto fail;
                }

                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    int nPerChrom = number / NoChroms;

                    CPOS delta = PopDescr[chrom].ChromLen / (nPerChrom + 1);

                    for (j = 1; j <= nPerChrom; j++)
                    {
                        char tmpstring[MAXNAME];
                        int nchars;
                        CPOS pos = delta * j;

                        nchars = snprintf(tmpstring, sizeof(tmpstring), "%s%d",
                                          name, lNR++);

                        if ((nchars < 0) || (nchars >= MAXNAME))
                        {
                            info(-2, "Generated locus name is too long\n");

                            goto fail;
                        }

                        if (!append_mappoint(&newMPs, &cMPs, chrom + 1, pos,
                                             tmpstring, classname))
                            goto fail_alloc;
                    }
                }

                number = 1;
                status = 0;
            }

            else
            {
                int chrom;
                CPOS pos;

                if ((1 != number) ||
                    !parse_chrom_pos_token(token, &chrom, &pos))
                {
                    info(-2, "Invalid map specification: %s\n", token);

                    goto fail;
                }

                if ((chrom < 1) || (chrom > NoChroms))
                {
                    info(-2, "Wrong chromosome (%d)\n", chrom);

                    goto fail;
                }

                if ((pos < 0) || (pos > PopDescr[chrom - 1].ChromLen))
                {
                    info(-2,
                         "Chromosome %d: Locusposition outside of "
                         "the chromosome (0 <= %f <= %f)\n",
                         chrom, pos, PopDescr[chrom - 1].ChromLen);

                    goto fail;
                }

                if (!append_mappoint(&newMPs, &cMPs, chrom, pos, name,
                                     classname))
                    goto fail_alloc;

                number = 1;
                status = 0;
            }

            break;
        }

        token = strtok(NULL, "[]', \t\r\n");
    }

    if (0 != status)
    {
        info(-2, "Incomplete map description\n");
        goto fail;
    }

    if (0 == cMPs)
    {
        info(-2, "Linkage map contains no loci\n");
        goto fail;
    }

    qsort(newMPs, cMPs, sizeof(TMappoint), mnamecmp);

    for (i = 1; i < cMPs; i++)
    {
        if (0 == strcmp(newMPs[i - 1].name, newMPs[i].name))
        {
            info(-2, "Duplicate locus name: %s\n", newMPs[i].name);

            goto fail;
        }
    }

    qsort(newMPs, cMPs, sizeof(TMappoint), mccmp);

    fp = fileopen(*fName, "w");

    if (NULL == fp)
        goto fail;

    outputOpened = 1;

    for (i = 0; i < cMPs; i++)
    {
        if (0 > fprintf(fp, "%d\t%.17g\t%s\t%s\n", newMPs[i].chrom,
                        newMPs[i].Pos, newMPs[i].name, newMPs[i].class))
        {
            info(-2, "Error writing map file %s\n", *fName);

            goto fail;
        }
    }

    if (0 != fclose(fp))
    {
        fp = NULL;

        info(-2, "Can't close the file stream\n");

        goto fail;
    }

    fp = NULL;

    free_mappoints(newMPs, cMPs);
    free(description);

    *retval = 0;
    return;

fail_alloc:

    info(-2, "Memory allocation failed in generate_map_file\n");

fail:

    if (NULL != fp)
    {
        fclose(fp);
        fp = NULL;
    }

    if (outputOpened)
        remove(*fName);

    free_mappoints(newMPs, cMPs);
    free(description);

    return;
}

void generate_effect_file(char **fName, char **Description)
{

    int status = 0;
    char *token, classname[MAXNAME * 3], all0[MAXNAME * 3], all1[MAXNAME * 3];

    FILE *fp = fileopen(*fName, "w");

    if (NULL == fp)
        return;

    fprintf(fp, "# Automatisch erzeugt\n");
    fprintf(fp, "# %s\n", *Description);

    GetRNGstate();

    token = strtok(*Description, "[' ',',']");

    while (NULL != token)
    {
        switch (status)
        {

        case 0:
        {
            int tmpi;

            if (1 == sscanf(token, "%d", &tmpi))
            {
                fprintf(fp, "%d\n", tmpi);
            }

            else
            {
                if (strlen(token) >= sizeof(classname))
                {
                    info(-2, "Class name too long: %s\n", token);
                    goto fail;
                }

                strcpy(classname, token);
                strcpy(all0, "");
                strcpy(all1, "");
                status = 2;
            }

            break;
        }

        case 2:
            if (strlen(token) >= sizeof(all0))
            {
                info(-2, "Allele name too long: %s\n", token);
                goto fail;
            }

            strcpy(all0, token);
            status = 4;
            break;

        case 3:
        {
            float mean, var;
            int chrom, i;

            if ((2 == sscanf(token, "%f/%f", &mean, &var)) && (0.0 <= var))
            {
                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    for (i = 0; i < PopDescr[chrom].NrMappoints; i++)
                    {
                        if (NULL != strstr(PopDescr[chrom].Mappoint[i].class,
                                           classname))
                        {
                            if (0 == strcmp(all1, ""))
                                fprintf(fp, "%s\t%s\t%f\n",
                                        PopDescr[chrom].Mappoint[i].name, all0,
                                        norm_rand() * sqrt(var) + mean);
                            else
                                fprintf(fp, "%s\t%s\t%s\t%s\t%f\n",
                                        PopDescr[chrom].Mappoint[i].name, all0,
                                        PopDescr[chrom].Mappoint[i].name, all1,
                                        norm_rand() * sqrt(var) + mean);
                        }
                    }
                }

                status = 0;
            }

            else
            {
                info(-2, "Error\n");
                goto fail;
            }
            break;
        }

        case 4:

            if (0 == strcmp(token, "normal"))
                status = 3;
            else if (0 == strcmp(token, "uniform"))
                status = 6;
            else
            {
                if (strlen(token) >= sizeof(all1))
                {
                    info(-2, "Allele name too long: %s\n", token);
                    goto fail;
                }

                strcpy(all1, token);
                status = 5;
            }

            break;

        case 5:

            if (0 == strcmp(token, "normal"))
                status = 3;
            else if (0 == strcmp(token, "uniform"))
                status = 6;
            else
            {
                info(-2, "Fehler Verteilungsart erwartet: %s\n", token);
                goto fail;
            }

            break;

        case 6:
        {
            float eff;
            int chrom, i;

            if (1 == sscanf(token, "%f", &eff))
            {
                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    for (i = 0; i < PopDescr[chrom].NrMappoints; i++)
                    {
                        if (NULL != strstr(PopDescr[chrom].Mappoint[i].class,
                                           classname))
                        {
                            if (0 == strcmp(all1, ""))
                                fprintf(fp, "%s\t%s\t%f\n",
                                        PopDescr[chrom].Mappoint[i].name, all0,
                                        eff);
                            else
                                fprintf(fp, "%s\t%s\t%s\t%s\t%f\n",
                                        PopDescr[chrom].Mappoint[i].name, all0,
                                        PopDescr[chrom].Mappoint[i].name, all1,
                                        eff);
                        }
                    }
                }

                status = 0;
            }

            else
            {
                info(-2, "Error\n");
                goto fail;
            }
            break;
        }
        }

        token = strtok(NULL, "[' ',',']");
    }

fail:

    PutRNGstate();
    fileclose(fp);
}

unsigned short get_allele(THom *Hom, CPOS Pos)
{
    int start, end, pos;

    if ((NULL == Hom) || (NULL == Hom->Loc) || (1 > Hom->NoLociUsed))
        return 0;

    start = 0;
    end = Hom->NoLociUsed - 1;

    while (start <= end)
    {
        pos = start + (end - start) / 2;

        if (Pos < Hom->Loc[pos].Pos)
        {
            end = pos - 1;
        }
        else
        {
            if ((pos + 1 < Hom->NoLociUsed) && (Hom->Loc[pos + 1].Pos <= Pos))
            {
                start = pos + 1;
            }
            else
            {
                return Hom->Loc[pos].All;
            }
        }
    }

    return Hom->Loc[0].All;
}

IntNode *insertIntLinkList(IntNode *Genotype, int *Phaenotyp, int counter,
                           int *nrGt)
{
    IntNode *nodepos = Genotype;

    if (NULL == Genotype)
    {
        Genotype = (IntNode *)yCalloc(1, IntNode);
        Genotype->child = NULL;
        Genotype->parent = NULL;
        Genotype->prev = NULL;
        Genotype->next = NULL;

        if (0 == counter)
        {
            Genotype->number = 1;
            (*nrGt)++;
        }
        else
        {
            Genotype->number = Phaenotyp[0];
            Genotype->child = insertIntLinkList(Genotype->child, &Phaenotyp[1],
                                                counter - 1, nrGt);
            Genotype->child->parent = Genotype;
        }

        return Genotype;
    }

    if (0 == counter)
    {
        (Genotype->number)++;
        return Genotype;
    }
    do
    {
        if (nodepos->number != Phaenotyp[0])
        {
            if (nodepos->number > Phaenotyp[0])
            {

                IntNode *tmp;
                tmp = (IntNode *)yCalloc(1, IntNode);
                tmp->child = NULL;
                tmp->parent = nodepos->parent;
                tmp->prev = NULL;
                tmp->next = nodepos;
                nodepos->prev = tmp;

                if (nodepos->parent != NULL)
                    nodepos->parent->child = tmp;

                tmp->number = Phaenotyp[0];

                Genotype = nodepos = tmp;
            }
            else
            {
                if ((NULL != nodepos->next) &&
                    (nodepos->next->number <= Phaenotyp[0]))
                {
                    nodepos = nodepos->next;
                }
                else
                {
                    IntNode *tmp;
                    tmp = (IntNode *)yCalloc(1, IntNode);
                    tmp->child = NULL;
                    tmp->parent = nodepos->parent;
                    tmp->prev = nodepos;
                    tmp->next = nodepos->next;
                    tmp->number = Phaenotyp[0];
                    nodepos->next = tmp;

                    if (tmp->next != NULL)
                    {
                        tmp->next->prev = tmp;
                    }

                    nodepos = tmp;
                }
            }
        }
    }

    while (nodepos->number < Phaenotyp[0]);

    nodepos->child =
        insertIntLinkList(nodepos->child, &Phaenotyp[1], counter - 1, nrGt);

    nodepos->child->parent = nodepos;

    return Genotype;
}

IntNode *go_down(IntNode *Node)
{
    while (Node->child != NULL)
        Node = Node->child;
    return Node;
}

void print_path(IntNode *Node, int *_nrLoci, int *_nrGt, int *m)
{
    int i = *_nrLoci;

    do
    {
        m[(i--) * (*_nrGt)] = Node->number;
    } while (NULL != (Node = Node->parent));
}

IntNode *find_next_leaf(IntNode *Node)
{
    IntNode *tmp;

    while (NULL != Node->parent)
    {
        tmp = Node;
        Node = Node->parent;
        yFree(tmp);

        if (Node->next != NULL)
        {
            tmp = Node->next;
            yFree(Node);
            return (go_down(tmp));
        }
    }

    yFree(Node);
    return (NULL);
}

SEXP genotype_evaluate(SEXP _PopName, SEXP _Loci, SEXP _mode)
{
    TPopListEl *PopListEl = NULL;
    int mode = 0;
    char *PopName = NULL;

    if (0 == cMP)
    {
        info(-2, "No linkage map defined\n");
        return R_NilValue;
    }

    if (!isInteger(_mode) || (length(_mode) != 1))
    {
        info(-2, "Value of mode is not an integer\n");
        return R_NilValue;
    }
    else
    {
        mode = INTEGER(_mode)[0];
    }

    if (!isString(_PopName) || (length(_PopName) != 1))
    {
        info(-2, "Name of population is not a single string\n");
        return R_NilValue;
    }

    PopName = (char *)yCalloc(strlen(CHAR(STRING_ELT(_PopName, 0))) + 2, char);

    strcpy(PopName, CHAR(STRING_ELT(_PopName, 0)));

    if (!isString(_Loci))
    {
        info(-2, "Name of Loci is not a string\n");
        goto fail;
    }

    find_list_element(PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population in not defined: %s\n", PopName);
        goto fail;
    }

    else
    {

        typedef struct
        {
            int nr_marker;
            TMappoint *Mappoint;
            int *order;
        }

        TEvChrom;

        TEvChrom *EvChrom = (TEvChrom *)yCalloc(NoChroms, TEvChrom);
        int chrom, i, j, l;

        int nrLoci = 0;
        int *pnrLoci = &nrLoci;
        int nrGt = 0;
        int *pnrGt = &nrGt;

        SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;
        int *locus_idx = (int *)yCalloc(length(_Loci), int);

        POPULATION ind;
        CPOS Pos;

        startNode = NULL;

        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            EvChrom[chrom].nr_marker = 0;
            EvChrom[chrom].Mappoint = NULL;
            EvChrom[chrom].order = NULL;
        }

        for (l = 0; l < length(_Loci); l++)
        {
            locus_idx[l] = 0;

            for (chrom = 0; ((chrom < NoChroms) && (0 == locus_idx[l]));
                 chrom++)
            {
                for (j = 0;
                     ((j < PopDescr[chrom].NrMappoints) && (0 == locus_idx[l]));
                     j++)
                {
                    if (0 == strcmp(PopDescr[chrom].Mappoint[j].name,
                                    CHAR(STRING_ELT(_Loci, l))))
                    {
                        locus_idx[l] = 1;
                        EvChrom[chrom].nr_marker++;

                        if (EvChrom[chrom].nr_marker == 1)
                        {
                            EvChrom[chrom].Mappoint =
                                (TMappoint *)yCalloc(1, TMappoint);

                            EvChrom[chrom].order = (int *)yCalloc(1, int);
                        }

                        else
                        {
                            EvChrom[chrom].Mappoint = (TMappoint *)yRealloc(
                                EvChrom[chrom].Mappoint,
                                EvChrom[chrom].nr_marker, TMappoint);

                            EvChrom[chrom].order =
                                (int *)yRealloc(EvChrom[chrom].order,
                                                EvChrom[chrom].nr_marker, int);
                        }

                        EvChrom[chrom].Mappoint[EvChrom[chrom].nr_marker - 1] =
                            PopDescr[chrom].Mappoint[j];

                        EvChrom[chrom].order[EvChrom[chrom].nr_marker - 1] =
                            nrLoci;

                        nrLoci += NoHoms;
                    }
                }
            }
        }

        if (0 == nrLoci)
        {
            info(-2, "Loci not found\n");
            goto fail_data;
        }

        if (0 == PopListEl->Pop->NoInds)
        {
            info(-2, "Empty population: %s\n", PopName);
            goto fail_data;
        }

        {
            TChrom *Chrom;
            TInd *Ind;
            IntNode *Node;
            int *alleles = (int *)yCalloc(NoHoms, int);
            int *IndGT = (int *)yCalloc(nrLoci + 1, int);
            int hom, tmp_counter = 0;

            for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
            {
                Ind = get_individual(PopListEl->Pop, ind);

                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    for (l = 0; l < EvChrom[chrom].nr_marker; l++)
                    {
                        int offset = EvChrom[chrom].order[l];

                        if (NULL == Ind->genotype)
                        {
                            Chrom = &(Ind->Chrom[chrom]);
                            Pos = EvChrom[chrom].Mappoint[l].Pos;

                            for (hom = 0; hom < NoHoms; hom++)
                            {
                                alleles[hom] =
                                    get_allele(&(Chrom->Hom[hom]), Pos);
                            }
                        }

                        else
                        {
                            for (hom = 0; hom < NoHoms; hom++)
                            {
                                alleles[hom] =
                                    Ind->genotype
                                        [EvChrom[chrom].Mappoint[l].counter *
                                             NoHoms +
                                         hom];
                            }
                        }

                        switch (mode)
                        {

                        case 0:
                            break;

                        case 1:

                            if (EvChrom[chrom].nr_marker == 1)
                                R_isort(&alleles[0], NoHoms);

                            break;

                        default:
                            R_isort(&alleles[0], NoHoms);

                            break;
                        }

                        for (hom = 0; hom < NoHoms; hom++)
                        {
                            IndGT[offset + hom] = alleles[hom];
                        }
                    }
                }

                startNode = insertIntLinkList(startNode, IndGT, nrLoci, pnrGt);
            }

            {
                int *_matrix = (int *)yCalloc((nrLoci + 1) * nrGt, int);
                int offset;
                char tmp[1024];
                int N = 0;

                Node = startNode;

                if (NULL != (Node = go_down(Node)))
                {
                    do
                    {
                        print_path(Node, pnrLoci, pnrGt,
                                   &_matrix[tmp_counter++]);
                    }

                    while (NULL != (Node = find_next_leaf(Node)));

                    startNode = NULL;
                }

                PROTECT(matrix = allocMatrix(INTSXP, nrGt, nrLoci + 2));
                PROTECT(matrixDim = allocVector(INTSXP, 2));
                INTEGER(matrixDim)[0] = nrGt;
                INTEGER(matrixDim)[1] = nrLoci + 2;
                setAttrib(matrix, R_DimSymbol, matrixDim);

                PROTECT(matrixDimNames = allocVector(VECSXP, 2));
                PROTECT(xAxis = allocVector(STRSXP, nrLoci + 2));
                PROTECT(yAxis = allocVector(STRSXP, nrGt));

                for (i = offset = 0; i < nrLoci + 1; i++)
                {
                    for (j = 0; j < nrGt; j++, offset++)
                    {
                        INTEGER(matrix)[offset] = _matrix[offset];

                        if (i == nrLoci)
                        {
                            N += _matrix[offset];
                        }
                    }
                }

                for (j = 0; j < nrGt; j++, offset++)
                {
                    INTEGER(matrix)[offset] = N;
                }

                for (i = 0, offset = 0; i < length(_Loci); i++)
                {
                    if (locus_idx[i])
                    {
                        for (hom = 1; hom <= NoHoms; hom++)
                        {
                            snprintf(tmp, sizeof(tmp), "%s-%d ",
                                     CHAR(STRING_ELT(_Loci, i)), hom);

                            SET_STRING_ELT(xAxis, offset++, mkChar(tmp));
                        }
                    }
                }

                SET_STRING_ELT(xAxis, offset++, mkChar("count"));

                SET_STRING_ELT(xAxis, offset++, mkChar("frequency"));

                for (i = 0; i < nrGt; i++)
                {
                    snprintf(tmp, sizeof(tmp), "%d ", i + 1);

                    SET_STRING_ELT(yAxis, i, mkChar(tmp));
                }

                SET_VECTOR_ELT(matrixDimNames, 1, xAxis);
                SET_VECTOR_ELT(matrixDimNames, 0, yAxis);
                setAttrib(matrix, R_DimNamesSymbol, matrixDimNames);

                UNPROTECT(5);

                yFree(PopName);

                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    yFree(EvChrom[chrom].Mappoint);
                    yFree(EvChrom[chrom].order);
                }

                yFree(EvChrom);
                yFree(locus_idx);
                yFree(alleles);
                yFree(IndGT);
                yFree(_matrix);

                return (matrix);
            }
        }

    fail_data:

        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            yFree(EvChrom[chrom].Mappoint);
            yFree(EvChrom[chrom].order);
        }

        yFree(EvChrom);
        yFree(locus_idx);
    }

fail:

    yFree(PopName);
    return (R_NilValue);
}

SEXP allele_evaluate(SEXP _PopName, SEXP _Loci, SEXP _missing, SEXP _bg)

{
    char *PopName = NULL;
    int bg = 0;
    int missing = 0;
    TPopListEl *PopListEl = NULL;

    if (0 == cMP)
    {
        info(-2, "No linkage map defined\n");
        return R_NilValue;
    }

    if (!isString(_PopName) || (length(_PopName) != 1))
    {
        info(-2, "Name of population is not a single string\n");
        return R_NilValue;
    }

    PopName = (char *)yCalloc(strlen(CHAR(STRING_ELT(_PopName, 0))) + 2, char);

    strcpy(PopName, CHAR(STRING_ELT(_PopName, 0)));

    if (!isString(_Loci))
    {
        info(-2, "Name of Loci is not a string\n");
        goto fail;
    }

    if (!isInteger(_bg) || (length(_bg) != 1))
    {
        info(-2, "Value of background allele is not an integer\n");
        goto fail;
    }
    else
        bg = INTEGER(_bg)[0];

    if (!isInteger(_missing) || (length(_missing) != 1))
    {
        info(-2, "Value of missing is not an integer\n");
        goto fail;
    }
    else
        missing = INTEGER(_missing)[0];

    find_list_element(PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined: %s\n", PopName);
        goto fail;
    }

    else
    {

        int *locus_idx = (int *)yCalloc(length(_Loci), int);

        int nrGt = 0;
        int *pnrGt = &nrGt;

        typedef struct
        {
            int nr_marker;
            TMappoint *Mappoint;
            int *order;
        }

        TEvChrom;

        TEvChrom *EvChrom = (TEvChrom *)yCalloc(NoChroms, TEvChrom);
        int j, l, chrom;
        int nrLoci = 0;
        int *pnrLoci = &nrLoci;

        startNode = NULL;

        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            EvChrom[chrom].nr_marker = 0;
            EvChrom[chrom].Mappoint = NULL;
            EvChrom[chrom].order = NULL;
        }

        for (l = 0; l < length(_Loci); l++)
        {
            locus_idx[l] = 0;

            for (chrom = 0; ((chrom < NoChroms) && (0 == locus_idx[l]));
                 chrom++)
            {
                for (j = 0;
                     ((j < PopDescr[chrom].NrMappoints) && (0 == locus_idx[l]));
                     j++)
                {
                    if (0 == strcmp(PopDescr[chrom].Mappoint[j].name,
                                    CHAR(STRING_ELT(_Loci, l))))
                    {
                        locus_idx[l] = 1;
                        EvChrom[chrom].nr_marker++;

                        if (1 == EvChrom[chrom].nr_marker)
                        {
                            EvChrom[chrom].Mappoint =
                                (TMappoint *)yCalloc(1, TMappoint);
                            EvChrom[chrom].order = (int *)yCalloc(1, int);
                        }

                        else
                        {
                            EvChrom[chrom].Mappoint = (TMappoint *)yRealloc(
                                EvChrom[chrom].Mappoint,
                                EvChrom[chrom].nr_marker, TMappoint);

                            EvChrom[chrom].order =
                                (int *)yRealloc(EvChrom[chrom].order,
                                                EvChrom[chrom].nr_marker, int);
                        }

                        EvChrom[chrom].Mappoint[EvChrom[chrom].nr_marker - 1] =
                            PopDescr[chrom].Mappoint[j];

                        EvChrom[chrom].order[EvChrom[chrom].nr_marker - 1] =
                            nrLoci;

                        nrLoci++;
                    }
                }
            }
        }

        if (0 == nrLoci)
        {
            info(-2, "Loci not found\n");
            goto fail_data;
        }

        else if (0 == PopListEl->Pop->NoInds)
        {
            info(-2, "Empty population %s\n", PopName);
            goto fail_data;
        }

        else
        {
            CPOS Pos;
            int chrom, hom, l, temp, i;
            int *all = (int *)yCalloc(nrLoci, int);
            POPULATION ind;
            TInd *Ind;
            TChrom *Chrom;

            for (i = 0; i < nrLoci; i++)
            {
                all[i] = 0;
            }

            for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
            {
                Ind = get_individual(PopListEl->Pop, ind);

                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    for (l = 0; l < EvChrom[chrom].nr_marker; l++)
                    {
                        temp = EvChrom[chrom].order[l];

                        if (NULL == Ind->genotype)
                        {
                            Chrom = &(Ind->Chrom[chrom]);
                            Pos = EvChrom[chrom].Mappoint[l].Pos;

                            for (hom = 0; hom < NoHoms; hom++)
                            {
                                if (0 == (all[temp] = get_allele(
                                              &(Chrom->Hom[hom]), Pos)))
                                    info(-2, "Allele 0 not allowed\n");

                                if ((1 == missing) && (all[temp] == bg))
                                {
                                    all[temp] = 0;
                                }

                                else
                                {
                                    startNode = insertIntLinkList(
                                        startNode, all, nrLoci, pnrGt);
                                }
                            }
                        }

                        else
                        {
                            for (hom = 0; hom < NoHoms; hom++)
                            {
                                if (0 == (all[temp] =
                                              Ind->genotype[EvChrom[chrom]
                                                                    .Mappoint[l]
                                                                    .counter *
                                                                NoHoms +
                                                            hom]))
                                    info(-2, "Allele 0 not allowed\n");

                                if ((1 == missing) && (all[temp] == bg))
                                {
                                    all[temp] = 0;
                                }

                                else
                                {
                                    startNode = insertIntLinkList(
                                        startNode, all, nrLoci, pnrGt);
                                }
                            }
                        }

                        all[temp] = 0;
                    }
                }
            }

            for (chrom = 0; chrom < NoChroms; chrom++)
            {
                yFree(EvChrom[chrom].Mappoint);
                yFree(EvChrom[chrom].order);
            }

            {
                int *_matrix = (int *)yCalloc((nrLoci + 1) * nrGt, int);
                int offset;
                char tmp[1024];
                int N = 0;
                IntNode *Node;
                int tmp_counter = 0;
                SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;

                Node = startNode;

                if ((NULL != Node) && (NULL != (Node = go_down(Node))))
                {
                    do
                    {
                        print_path(Node, pnrLoci, pnrGt,
                                   &_matrix[tmp_counter++]);
                    }

                    while (NULL != (Node = find_next_leaf(Node)));

                    startNode = NULL;
                }

                PROTECT(matrix = allocMatrix(INTSXP, nrGt, nrLoci + 2));
                PROTECT(matrixDim = allocVector(INTSXP, 2));
                INTEGER(matrixDim)[0] = nrGt;
                INTEGER(matrixDim)[1] = nrLoci + 2;
                setAttrib(matrix, R_DimSymbol, matrixDim);
                PROTECT(matrixDimNames = allocVector(VECSXP, 2));
                PROTECT(xAxis = allocVector(STRSXP, nrLoci + 2));
                PROTECT(yAxis = allocVector(STRSXP, nrGt));

                for (i = offset = 0; i < nrLoci + 1; i++)
                {
                    for (j = 0; j < nrGt; j++, offset++)
                    {
                        INTEGER(matrix)[offset] = _matrix[offset];

                        if (i == nrLoci)
                        {
                            N += _matrix[offset];
                        }
                    }
                }

                for (j = 0; j < nrGt; j++, offset++)
                {
                    INTEGER(matrix)[offset] = N;
                }

                for (i = 0, offset = 0; i < length(_Loci); i++)
                {
                    if (locus_idx[i])
                    {
                        snprintf(tmp, sizeof(tmp), "%s",
                                 CHAR(STRING_ELT(_Loci, i)));

                        SET_STRING_ELT(xAxis, offset++, mkChar(tmp));
                    }
                }

                SET_STRING_ELT(xAxis, offset++, mkChar("count"));
                SET_STRING_ELT(xAxis, offset++, mkChar("frequency"));

                for (i = 0; i < nrGt; i++)
                {
                    snprintf(tmp, sizeof(tmp), "%d ", i + 1);

                    SET_STRING_ELT(yAxis, i, mkChar(tmp));
                }

                SET_VECTOR_ELT(matrixDimNames, 1, xAxis);
                SET_VECTOR_ELT(matrixDimNames, 0, yAxis);
                setAttrib(matrix, R_DimNamesSymbol, matrixDimNames);

                UNPROTECT(5);

                yFree(PopName);
                yFree(locus_idx);
                yFree(EvChrom);
                yFree(all);
                yFree(_matrix);

                return (matrix);
            }
        }

    fail_data:

        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            yFree(EvChrom[chrom].Mappoint);
            yFree(EvChrom[chrom].order);
        }

        yFree(locus_idx);
        yFree(EvChrom);
    }

fail:

    yFree(PopName);
    return R_NilValue;
}

void evaluate_locus(char **PopName, char **Loci, int *_matrix, int *noLoci,
                    int *bg)

{
    TPopListEl *PopListEl;
    char *pch;
    int chrom, i = 0, j;
    POPULATION noHomo, noHetero, noMissing, ind;

    if (2 != NoHoms)
    {
        info(-2, "This function is only defined for diploid organisms!\n");
        return;
    }

    if (0 == cMP)
    {
        info(-2, "No linkage map defined\n");
        return;
    }

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined: %s\n", *PopName);
        return;
    }

    pch = strtok(*Loci, " ");

    while (pch != NULL)
    {
        noHomo = 0;
        noHetero = 0;
        noMissing = 0;

        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            if (PopDescr[chrom].NrMappoints > 0)
            {
                for (j = 0; j < PopDescr[chrom].NrMappoints; j++)
                {
                    if (0 == strcmp(PopDescr[chrom].Mappoint[j].name, pch))
                    {
                        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
                        {
                            TInd *Ind = get_individual(PopListEl->Pop, ind);

                            if (NULL == Ind->genotype)
                            {
                                info(-2,
                                     "first genotype_population (chrom=%d, "
                                     "j=%d, ind=%ld, pch=%s)",
                                     chrom, j, ind, pch);
                                return;
                            }

                            if ((Ind->genotype[2 * PopDescr[chrom]
                                                       .Mappoint[j]
                                                       .counter] == *bg) ||
                                (Ind->genotype
                                     [2 * PopDescr[chrom].Mappoint[j].counter +
                                      1] == *bg))
                            {
                                noMissing++;
                            }

                            else if (Ind->genotype[2 * PopDescr[chrom]
                                                           .Mappoint[j]
                                                           .counter] !=
                                     Ind->genotype[2 * PopDescr[chrom]
                                                           .Mappoint[j]
                                                           .counter +
                                                   1])
                            {
                                noHetero++;
                            }

                            else
                            {
                                noHomo++;
                            }
                        }
                    }
                }
            }
        }

        _matrix[i] = noHomo;
        _matrix[i + (*noLoci)] = noHetero;
        _matrix[i + 2 * (*noLoci)] = noMissing;
        i++;
        pch = strtok(NULL, " ");
    }
}

SEXP evaluate_haplotype(SEXP _PopName, SEXP _Loci, SEXP _missing, SEXP _bg)
{
    int bg = 0, missing = 0;
    char *locus;
    int *values = NULL;
    int output_rows = 0, output_columns = 0;

    if (!isString(_PopName) || (length(_PopName) != 1))
    {
        info(-2, "Name of population is not a single string\n");
        return R_NilValue;
    }

    if (!isString(_Loci) || (length(_Loci) != 1))
    {
        info(-2, "Name of Loci is not a single string\n");
        return R_NilValue;
    }

    if (!isInteger(_bg) || (length(_bg) != 1))
    {
        info(-2, "Value of background allele is not an integer\n");
        return R_NilValue;
    }
    else
    {
        bg = INTEGER(_bg)[0];
    }

    if (!isInteger(_missing) || (length(_missing) != 1))
    {
        info(-2, "Value of missing variable is not an integer\n");
        return R_NilValue;
    }
    else
    {
        missing = INTEGER(_missing)[0];
    }

    {
        TPopListEl *PopListEl;
        POPULATION ind;
        TMappoint *Mappoint = NULL;

        int nLoci = 0;
        int nGt = 0;
        int *pGt = &nGt;
        int *hap = NULL;
        int free_startNode = 0;

        char *PopName =
            (char *)yCalloc(strlen(CHAR(STRING_ELT(_PopName, 0))) + 2, char);
        char *Loci =
            (char *)yCalloc(strlen(CHAR(STRING_ELT(_Loci, 0))) + 2, char);

        strcpy(PopName, CHAR(STRING_ELT(_PopName, 0)));
        strcpy(Loci, CHAR(STRING_ELT(_Loci, 0)));

        find_list_element(PopName, &PopListEl);

        if (!PopListEl)
        {
            info(-2, "Population not defined: %s\n", PopName);
            goto fail;
        }

        if (1 > PopListEl->Pop->NoInds)
        {
            info(-2, "Empty population %s\n", PopName);
            goto fail;
        }

        {
            int lastLoci;
            int chrom, loc;

            locus = strtok(Loci, " ");

            while (locus != NULL)
            {
                lastLoci = nLoci;

                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    if (PopDescr[chrom].NrMappoints > 0)
                    {
                        for (loc = 0; loc < PopDescr[chrom].NrMappoints; loc++)
                        {
                            if (0 == strcmp(PopDescr[chrom].Mappoint[loc].name,
                                            locus))
                            {
                                if (0 == nLoci)
                                {
                                    Mappoint =
                                        (TMappoint *)yCalloc(1, TMappoint);
                                }

                                else
                                {
                                    Mappoint = (TMappoint *)yRealloc(
                                        Mappoint, nLoci + 1, TMappoint);
                                }

                                Mappoint[nLoci] =
                                    (PopDescr[chrom].Mappoint[loc]);
                                nLoci++;
                                break;
                            }
                        }
                    }
                }

                if (lastLoci == nLoci)
                {
                    info(-2, "Cant find locus %s\n ", locus);
                    goto fail;
                }

                locus = strtok(NULL, " ");
            }
        }

        {
            TInd *Ind;
            int loc, hom, nrmissing;

            hap = (int *)yCalloc(nLoci, int);

            startNode = NULL;

            for (loc = 0; loc < nLoci; loc++)
            {
                hap[loc] = 0;
            }

            startNode = insertIntLinkList(startNode, hap, nLoci, pGt);
            free_startNode = 1;

            for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
            {
                Ind = get_individual(PopListEl->Pop, ind);

                if (NULL == Ind->genotype)
                {
                    info(-2, "First genotype_population\n");
                    goto fail;
                }

                for (hom = 0; hom < NoHoms; hom++)
                {
                    nrmissing = 0;

                    for (loc = 0; loc < nLoci; loc++)
                    {
                        hap[loc] =
                            Ind->genotype[Mappoint[loc].counter * NoHoms + hom];

                        if ((1 == missing) && (bg == hap[loc]))
                        {
                            nrmissing++;
                        }
                    }

                    if (0 == nrmissing)
                    {
                        startNode =
                            insertIntLinkList(startNode, hap, nLoci, pGt);
                    }
                }
            }

            yFree(hap);
        }

        {
            int y_dim = nGt;
            int x_dim = nLoci + 1;
            IntNode *Node = startNode;
            int tmpcounter = 0;

            if (y_dim < 0 || x_dim < 1 || (size_t)y_dim > SIZE_MAX / sizeof(int) / (size_t)x_dim)
            {
                info(-2,"Haplotype result dimensions overflow\n");
                goto fail;
            }
            values = calloc(y_dim ? (size_t)y_dim * x_dim : 1,sizeof(int));
            if (!values)
            {
                info(-2,"Cannot allocate haplotype counts\n");
                goto fail;
            }
            output_rows = y_dim;
            output_columns = x_dim;
            if ((NULL != Node) && (NULL != (Node = go_down(Node))))
            {
                int *path = (int *)yCalloc(nGt * x_dim, int);
                int ttt = 0;
                int iii;

                do
                {
                    print_path(Node, &nLoci, &nGt, &path[ttt++]);
                }

                while (NULL != (Node = find_next_leaf(Node)));

                startNode = NULL;

                for (iii = 0; iii < x_dim; iii++)
                {
                    for (ttt = 0; ttt < y_dim; ttt++, tmpcounter++)
                    {
                        values[tmpcounter] = path[tmpcounter];
                    }
                }

                yFree(path);
            }


        }

        if (2 <= nLoci)
        {
            int x;
            int **table = NULL;
            int dimLA = 0, dimLB = 0;
            int lx, ly, xPos, yPos;

            for (x = 0; x < nGt; x++)
            {
                xPos = 0;
                yPos = 0;

                if ((0 == dimLA) && (0 == dimLB))
                {
                    dimLA = 1;
                    dimLB = 1;
                    table = yCalloc(dimLA + 1, int *);
                    table[0] = (int *)yCalloc(dimLB + 1, int);
                    table[1] = (int *)yCalloc(dimLB + 1, int);
                    table[0][0] = values[x + 2 * nGt];
                    table[1][0] = values[x];
                    table[0][1] = values[x + nGt];
                    table[1][1] = values[x + 2 * nGt];
                }

                else
                {
                    for (lx = 1; lx <= dimLB; lx++)
                    {
                        if (table[0][lx] == values[x + nGt])
                        {
                            xPos = lx;
                        }
                    }

                    if (0 == xPos)
                    {
                        dimLB++;

                        for (ly = 0; ly <= dimLA; ly++)
                        {
                            table[ly] =
                                (int *)yRealloc(table[ly], dimLB + 1, int);
                            table[ly][dimLB] = 0;
                        }

                        lx = 1;

                        while (table[0][lx] < values[x + nGt])
                        {
                            lx++;

                            if (lx >= dimLB)
                                break;
                        }

                        xPos = lx;

                        if (lx != dimLB)
                        {
                            for (ly = 0; ly <= dimLA; ly++)
                            {
                                for (lx = dimLB; lx > xPos; lx--)
                                {
                                    table[ly][lx] = table[ly][lx - 1];
                                }
                            }

                            for (ly = 0; ly <= dimLA; ly++)
                            {
                                table[ly][xPos] = 0;
                            }
                        }

                        table[0][xPos] = values[x + nGt];
                    }

                    for (ly = 1; ly <= dimLA; ly++)
                    {
                        if (table[ly][0] == values[x])
                        {
                            yPos = ly;
                        }
                    }

                    if (0 == yPos)
                    {
                        dimLA++;
                        table = (int **)yRealloc(table, dimLA + 1, int *);
                        table[dimLA] = (int *)yCalloc(dimLB + 1, int);

                        table[dimLA][0] = values[x];

                        for (lx = 1; lx <= dimLB; lx++)
                        {
                            table[dimLA][lx] = 0;
                        }

                        yPos = dimLA;
                    }

                    if ((xPos != 0) && (yPos != 0))
                    {
                        table[0][0] += values[x + 2 * nGt];
                        table[yPos][xPos] = values[x + 2 * nGt];
                    }
                }
            }

            for (ly = 0; ly <= dimLA; ly++)
            {
                yFree(table[ly]);
            }

            yFree(table);
        }

        yFree(Mappoint);
        yFree(PopName);
        yFree(Loci);

        st_work_adopt(values,free);
        int rows = output_rows ? output_rows-1 : 0, np = 0;
        SEXP result = PROTECT(allocVector(VECSXP,output_columns));np++;
        SEXP names = PROTECT(allocVector(STRSXP,output_columns));np++;
        SEXP rownames = PROTECT(allocVector(INTSXP,rows ? 2 : 0));np++;
        SEXP cls = PROTECT(mkString("data.frame"));np++;
        for (int j = 0; j < output_columns; j++)
        {
            SEXP column = allocVector(INTSXP,rows);
            SET_VECTOR_ELT(result,j,column);
            for (int i = 0; i < rows; i++) INTEGER(column)[i] = values[i+1+(size_t)output_rows*j];
        }
        char *labels = st_work_alloc(strlen(CHAR(STRING_ELT(_Loci,0)))+1,1);
        strcpy(labels,CHAR(STRING_ELT(_Loci,0)));
        char *label = strtok(labels," ");int at = 0;
        while (label)
        {
            if (at >= output_columns-1) error("Invalid haplotype locus layout");
            SET_STRING_ELT(names,at++,mkChar(label));
            label = strtok(NULL," ");
        }
        if (at != output_columns-1) error("Invalid haplotype locus layout");
        SET_STRING_ELT(names,at,mkChar("count"));
        if (rows) { INTEGER(rownames)[0] = NA_INTEGER;INTEGER(rownames)[1] = -rows; }
        setAttrib(result,R_NamesSymbol,names);
        setAttrib(result,R_RowNamesSymbol,rownames);
        setAttrib(result,R_ClassSymbol,cls);
        UNPROTECT(np);return result;


    fail:

        free(values);

        if (free_startNode)
        {
            IntNode *Node = go_down(startNode);

            while (NULL != (Node = find_next_leaf(Node)))
            {
            }

            startNode = NULL;
        }

        yFree(hap);
        yFree(Mappoint);
        yFree(PopName);
        yFree(Loci);

        return R_NilValue;
    }
}

void haplotype_init(char **PopName, int *_missing, int *_doubleHetero, int *_bg)

{
    PopListEl = NULL;
    find_list_element(*PopName, &PopListEl);

    if (!PopListEl)
    {
        info(-2, "Population not defined: %s\n", *PopName);
        return;
    }

    if (1 > PopListEl->Pop->NoInds)
    {
        info(-2, "Population %s is empty\n", *PopName);
        return;
    }

    bg = *_bg;

    missing = *_missing;

    doubleHetero = *_doubleHetero;
}

void haplotype_next(int *_locusA, int *_locusB, int *_noTI, int *_noI,
                    int *_nodH, int *_dLocA, int *_dLocB)

{
    int nGt = 0;
    int *pGt = &nGt;
    TInd *Ind;
    POPULATION ind;
    int nLoci = 2, hap[2], dH = 0;

    *_noTI = 0;
    *_noI = 0;
    *_nodH = 0;
    *_dLocA = 0;
    *_dLocB = 0;

    if ((NULL == PopListEl) || (NULL == PopListEl->Pop))
    {
        info(-2, "First initialize haplotype evaluation\n");
        return;
    }

    if (2 != NoHoms)
    {
        info(-2, "Haplotype evaluation is implemented for diploids only\n");
        return;
    }

    if ((0 > *_locusA) || (cMP <= *_locusA) || (0 > *_locusB) ||
        (cMP <= *_locusB))
    {
        info(-2, "Locus index out of range\n");
        return;
    }

    if (NULL != table)
    {
        int ly;

        for (ly = 0; ly <= dimLA; ly++)
        {
            yFree(table[ly]);
        }

        yFree(table);
    }

    yFree(genotype);
    nrGeno = 0;
    dimLA = 0;
    dimLB = 0;

    locusA = *_locusA;
    locusA *= 2;
    locusB = *_locusB;
    locusB *= 2;

    startNode = NULL;
    *_noTI = PopListEl->Pop->NoInds;

    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
        Ind = get_individual(PopListEl->Pop, ind);

        if (NULL == Ind->genotype)
        {
            IntNode *Node;

            info(-2, "first genotype_population\n");

            if (NULL != startNode)
            {
                Node = go_down(startNode);

                while (NULL != (Node = find_next_leaf(Node)))
                {
                }

                startNode = NULL;
            }

            yFree(genotype);
            nrGeno = 0;
            *_noI = 0;

            return;
        }

        if (doubleHetero &&
            (Ind->genotype[locusA] != Ind->genotype[locusA + 1]) &&
            (Ind->genotype[locusB] != Ind->genotype[locusB + 1]))
        {
            dH++;
        }

        else
        {
            if (missing || ((Ind->genotype[locusA] != bg) &&
                            (Ind->genotype[locusA + 1] != bg) &&
                            (Ind->genotype[locusB] != bg) &&
                            (Ind->genotype[locusB + 1] != bg)))
            {
                hap[0] = Ind->genotype[locusA];
                hap[1] = Ind->genotype[locusB];
                startNode = insertIntLinkList(startNode, hap, 2, pGt);
                hap[0] = Ind->genotype[locusA + 1];
                hap[1] = Ind->genotype[locusB + 1];
                (*_noI)++;
                nrGeno += 4;

                if (4 == nrGeno)
                    genotype = (int *)yCalloc(nrGeno, int);
                else
                    genotype = (int *)yRealloc(genotype, nrGeno, int);

                genotype[nrGeno - 4] = Ind->genotype[locusA];
                genotype[nrGeno - 3] = Ind->genotype[locusB];
                genotype[nrGeno - 2] = Ind->genotype[locusA + 1];
                genotype[nrGeno - 1] = Ind->genotype[locusB + 1];

                startNode = insertIntLinkList(startNode, hap, 2, pGt);
            }
        }
    }

    {
        IntNode *Node = startNode;

        if ((NULL != Node) && (NULL != (Node = go_down(Node))))
        {
            int ttt = 0;
            int x, xPos = 0, yPos = 0;
            int *path = (int *)yCalloc(nGt * 3, int);
            int lx, ly;

            for (x = 0; x < nGt * 3; x++)
            {
                path[x] = 0;
            }

            do
            {
                print_path(Node, &nLoci, &nGt, &path[ttt++]);
            } while (NULL != (Node = find_next_leaf(Node)));

            startNode = NULL;

            for (x = 0; x < nGt; x++)
            {
                if ((0 == dimLA) && (0 == dimLB))
                {
                    dimLA = 1;
                    dimLB = 1;
                    table = yCalloc(dimLA + 1, int *);
                    table[0] = yCalloc(dimLB + 1, int);
                    table[1] = yCalloc(dimLB + 1, int);
                    table[0][0] = path[x + 2 * nGt];
                    table[1][0] = path[x];
                    table[0][1] = path[x + nGt];
                    table[1][1] = path[x + 2 * nGt];
                }
                else
                {
                    for (lx = 1, xPos = 0; lx <= dimLB; lx++)
                    {
                        if (table[0][lx] == path[x + nGt])
                        {
                            xPos = lx;
                        }
                    }
                    if (0 == xPos)
                    {
                        dimLB++;

                        for (ly = 0; ly <= dimLA; ly++)
                        {
                            table[ly] =
                                (int *)yRealloc(table[ly], dimLB + 1, int);
                            table[ly][dimLB] = 0;
                        }

                        lx = 1;

                        while (table[0][lx] < path[x + nGt])
                        {
                            lx++;

                            if (lx >= dimLB)
                                break;
                        }

                        xPos = lx;

                        if (lx != dimLB)
                        {
                            for (ly = 0; ly <= dimLA; ly++)
                            {
                                for (lx = dimLB; lx > xPos; lx--)
                                {
                                    table[ly][lx] = table[ly][lx - 1];
                                }
                            }

                            for (ly = 0; ly <= dimLA; ly++)
                            {
                                table[ly][xPos] = 0;
                            }
                        }

                        table[0][xPos] = path[x + nGt];
                    }

                    for (ly = 1, yPos = 0; ly <= dimLA; ly++)
                    {
                        if (table[ly][0] == path[x])
                        {
                            yPos = ly;
                        }
                    }

                    if (0 == yPos)
                    {
                        dimLA++;
                        table = (int **)yRealloc(table, dimLA + 1, int *);
                        table[dimLA] = (int *)yCalloc(dimLB + 1, int);

                        table[dimLA][0] = path[x];

                        for (lx = 1; lx <= dimLB; lx++)
                        {
                            table[dimLA][lx] = 0;
                        }

                        yPos = dimLA;
                    }

                    if ((xPos != 0) && (yPos != 0))
                    {
                        table[0][0] += path[x + 2 * nGt];
                        table[yPos][xPos] = path[x + 2 * nGt];
                    }
                }
            }

            yFree(path);
        }
    }

    *_nodH = dH;
    *_dLocA = dimLA;
    *_dLocB = dimLB;

    {
        int maxAllelA = 0;
        int maxAllelB = 0;

        int *transAllelA;
        int *transAllelB;

        int i;

        for (i = 1; i <= dimLA; i++)
        {
            if (maxAllelA < table[i][0])
            {
                maxAllelA = table[i][0];
            }
        }

        for (i = 1; i <= dimLB; i++)
        {
            if (maxAllelB < table[0][i])
            {
                maxAllelB = table[0][i];
            }
        }

        transAllelA = (int *)yCalloc(maxAllelA + 1, int);
        transAllelB = (int *)yCalloc(maxAllelB + 1, int);

        for (i = 0; i <= maxAllelA; i++)
        {
            transAllelA[i] = 0;
        }

        for (i = 0; i <= maxAllelB; i++)
        {
            transAllelB[i] = 0;
        }

        for (i = 1; i <= dimLA; i++)
        {
            transAllelA[table[i][0]] = i - 1;
        }

        for (i = 1; i <= dimLB; i++)
        {
            transAllelB[table[0][i]] = (i - 1) * dimLA;
        }

        for (i = 0; i < nrGeno; i += 2)
        {
            genotype[i] = transAllelA[genotype[i]];
            genotype[i + 1] = transAllelB[genotype[i + 1]];
        }

        yFree(transAllelA);
        yFree(transAllelB);
    }
}

SEXP haplotype_return_table(void)
{
    SEXP out = PROTECT(allocVector(VECSXP,dimLB)),
         names = PROTECT(allocVector(STRSXP,dimLB)),
         rows = PROTECT(allocVector(dimLA ? STRSXP : INTSXP,dimLA)),
         cls = PROTECT(mkString("data.frame"));
    for (int x = 0; x < dimLB; x++)
    {
        char label[200];snprintf(label,sizeof(label),"locusB.%d ",table[0][x+1]);
        SET_STRING_ELT(names,x,mkChar(label));
        SEXP column = allocVector(INTSXP,dimLA);SET_VECTOR_ELT(out,x,column);
        for (int y = 0; y < dimLA; y++) INTEGER(column)[y] = table[y+1][x+1];
    }
    for (int y = 0; y < dimLA; y++)
    {
        char label[200];snprintf(label,sizeof(label),"locusA.%d ",table[y+1][0]);
        SET_STRING_ELT(rows,y,mkChar(label));
    }
    setAttrib(out,R_NamesSymbol,names);setAttrib(out,R_RowNamesSymbol,rows);setAttrib(out,R_ClassSymbol,cls);
    UNPROTECT(4);return out;
}

SEXP haplotype_return_genotype(void)
{
    int n = nrGeno/2;
    SEXP out = PROTECT(allocVector(VECSXP,2)),
         names = PROTECT(allocVector(STRSXP,2)),
         rows = PROTECT(allocVector(n ? STRSXP : INTSXP,n)),
         cls = PROTECT(mkString("data.frame"));
    SET_STRING_ELT(names,0,mkChar("LocusA"));SET_STRING_ELT(names,1,mkChar("LocusB"));
    for (int j = 0; j < 2; j++)
    {
        SEXP column = allocVector(INTSXP,n);SET_VECTOR_ELT(out,j,column);
        for (int i = 0; i < n; i++) INTEGER(column)[i] = genotype[2*i+j];
    }
    for (int i = 0; i < n; i++)
    {
        char label[200];snprintf(label,sizeof(label),"Ind%d.%d",i/2+1,i%2+1);
        SET_STRING_ELT(rows,i,mkChar(label));
    }
    setAttrib(out,R_NamesSymbol,names);setAttrib(out,R_RowNamesSymbol,rows);setAttrib(out,R_ClassSymbol,cls);
    UNPROTECT(4);return out;
}

void haplotype_calc_d(double *_dd, double *_rr, double *_dprime)

{
    int x, y, xx, yy, sumxx, sumyy, sumTotal;
    double D, DD = 0, RA = 1, RB = 1, dPrime = 0;
    double p, q;

    RA = 1;
    RB = 1;

    for (x = 1, sumTotal = 0; x <= dimLB; x++)
    {
        for (y = 1; y <= dimLA; y++)
        {
            sumTotal += table[y][x];
        }
    }

    for (x = 1; x <= dimLB; x++)
    {

        for (yy = 1, sumyy = 0; yy <= dimLA; yy++)
        {
            sumyy += table[yy][x];
        }

        for (y = 1; y <= dimLA; y++)
        {

            for (xx = 1, sumxx = 0; xx <= dimLB; xx++)
            {
                sumxx += table[y][xx];
            }

            D = (((double)table[y][x] -
                  ((double)sumxx * (double)sumyy) / ((double)sumTotal)) /
                 ((double)sumTotal));

            DD += D * D;

            if (1 == x)
            {
                RA -= (double)(((double)sumxx / ((double)sumTotal)) *
                               ((double)sumxx / ((double)sumTotal)));
            }

            if (1 == y)
            {
                RB -= (double)(((double)sumyy / ((double)sumTotal)) *
                               ((double)sumyy / ((double)sumTotal)));
            }

            p = (double)sumxx / ((double)sumTotal);
            q = (double)sumyy / ((double)sumTotal);

            if (D < 0)
            {

                if ((p * q) < ((1 - p) * (1 - q)))
                {
                    dPrime += (p * q * D * (-1) / (p * q));
                }
                else
                {
                    dPrime += (p * q * D * (-1) / ((1 - p) * (1 - q)));
                }
            }
            else
            {

                if (((1 - p) * q) < (p * (1 - q)))
                {
                    dPrime += (p * q * D / ((1 - p) * q));
                }
                else
                {
                    dPrime += (p * q * D / (p * (1 - q)));
                }
            }
        }
    }

    *_dd = DD;
    *_rr = DD / (RA * RB);
    *_dprime = dPrime;
}

void haplotype_free(void)

{
    int ly;

    if (NULL != table)
    {
        for (ly = 0; ly <= dimLA; ly++)
        {
            yFree(table[ly]);
        }

        yFree(table);
    }

    yFree(genotype);
    nrGeno = 0;
    dimLA = 0;
    dimLB = 0;
}

void concat_population(char **NameP1, char **NameP2)
{
    TInd tmpInd;
    TInd *tmpP1Ind;
    POPULATION ind = 0, j, oldP1Alloc, newP1Alloc;
    TPopListEl *PopListEl;
    TPop *P1Pop = NULL;
    TPop *P2Pop = NULL;

    find_list_element(*NameP2, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined: %s!", *NameP2);
        return;
    }

    P2Pop = PopListEl->Pop;

    if (0 == strcmp(*NameP1, *NameP2))
    {
        info(-2, "Source and destination population must be different\n");
        return;
    }

    find_list_element(*NameP1, &PopListEl);

    if (NULL == PopListEl)
    {
        P1Pop = new_population(NameP1, &ind);
        if (NULL == P1Pop)
            return;
    }

    else
    {
        P1Pop = PopListEl->Pop;
    }

    ind = P1Pop->NoInds, j = 0, oldP1Alloc = P1Pop->NoIndsAlloc;
    if (P2Pop->NoIndsAlloc > LONG_MAX - oldP1Alloc)
    {
        info(-2, "Resulting population is too large\n");
        return;
    }
    newP1Alloc = oldP1Alloc + P2Pop->NoIndsAlloc;

    if (NULL == P1Pop->Ind)
    {
        if (0 < newP1Alloc)
        {
            tmpP1Ind = (TInd *)st_sim_alloc_array(
                (size_t)newP1Alloc, sizeof(*tmpP1Ind));

            if (NULL == tmpP1Ind)
            {
                info(-2, "Memory allocation failed in concat_population\n");
                return;
            }

            P1Pop->Ind = tmpP1Ind;
        }
    }

    else if (oldP1Alloc != newP1Alloc)
    {
        tmpP1Ind = (TInd *)st_sim_resize_array(
            P1Pop->Ind, (size_t)newP1Alloc, sizeof(*tmpP1Ind));

        if (NULL == tmpP1Ind)
        {
            info(-2, "Memory allocation failed in concat_population\n");
            return;
        }

        P1Pop->Ind = tmpP1Ind;
    }

    remove_evaluate_population(NameP2);
    remove_genotype_population(NameP2);

    if (NULL != PopListEl)
    {
        remove_evaluate_population(NameP1);
        remove_genotype_population(NameP1);
    }

    P1Pop->NoInds += P2Pop->NoInds;
    P1Pop->NoIndsAlloc = newP1Alloc;

    for (; j < P2Pop->NoIndsAlloc; ind++, j++)
    {
        if (ind < oldP1Alloc)
        {
            tmpInd = P1Pop->Ind[ind];
            P1Pop->Ind[ind] = P2Pop->Ind[j];
            P1Pop->Ind[P1Pop->NoIndsAlloc - j - 1] = tmpInd;
        }
        else
        {
            P1Pop->Ind[ind] = P2Pop->Ind[j];
        }
    }

    P2Pop->NoIndsAlloc = 0;
    remove_population(NameP2);

    info(1, "Population %s concatenated to population %s\n", *NameP2, *NameP1);
}

void sample_population(char **NameP1, char **NameP2, char **cNoI, int *rep)
{
    int chrom, hom, loc, r;
    int copy_genotype, copy_gvalue, copy_pvalue;
    int genotype_alloc_error = 0, gvalue_alloc_error = 0;
    int pvalue_alloc_error = 0, info_alloc_error = 0;
    POPULATION ind, NoI = atol(*cNoI), *samples = NULL;
    POPULATION *tmp_samples;
    TInd *P1Ind, *P2Ind;
    TPopListEl *PopListEl;

    TPop *P1Pop = NULL, *P2Pop = NULL;

    find_list_element(*NameP2, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined: %s \n", *NameP2);
        return;
    }

    P2Pop = PopListEl->Pop;

    if (0 == P2Pop->NoInds)
    {
        info(-2, "Population >%s< is empty\n", *NameP2);
        return;
    }

    if (0 == strcmp(*NameP1, *NameP2))
    {
        info(-2, "Source and destination population must be different\n");
        return;
    }

    if (-1 == NoI)
        NoI = P2Pop->NoInds;

    if (0 > NoI)
    {
        info(-2, "Incorrect number of individuals\n");
        return;
    }

    if ((0 == *rep) && (P2Pop->NoInds < NoI))
    {
        info(-2, "sampling only with repetition possible\n");
        return;
    }

    copy_genotype = st_population_all_genotyped(P2Pop);
    copy_gvalue = st_population_all_gvalue(P2Pop);
    copy_pvalue = st_population_all_pvalue(P2Pop);

    P1Pop = exist_population(*NameP1, NoI, 1);
    if (NULL == P1Pop)
        return;

    if (NoI != P1Pop->NoIndsAlloc)
        return;

    P1Pop->NoInds = NoI;

    for (ind = 0; ind < P1Pop->NoInds; ind++)
    {

        r = (int)gsl_rng_uniform_int(rng[0], (unsigned long int)P2Pop->NoInds);

        if (0 == *rep)
        {
            int j = 0;
            POPULATION ind_i;

            do
            {
                j = 0;

                for (ind_i = 0; ind_i < ind; ind_i++)
                {
                    if (r == *(samples + ind_i))
                    {
                        j = 1;
                        break;
                    }
                }

                if (0 == j)
                {
                    if (0 == ind)
                        tmp_samples =
                            (POPULATION *)malloc(sizeof(*tmp_samples));
                    else
                        tmp_samples = (POPULATION *)realloc(
                            samples, (size_t)(ind + 1) * sizeof(*tmp_samples));

                    if (NULL == tmp_samples)
                    {
                        info(-2,
                             "Memory allocation failed in sample_population\n");
                        goto fail_population;
                    }

                    samples = tmp_samples;
                    *(samples + ind) = r;
                }

                else
                    r = (int)gsl_rng_uniform_int(
                        rng[0], (unsigned long int)P2Pop->NoInds);
            }

            while (1 == j);
        }

        P1Ind = get_individual(P1Pop, ind);
        P2Ind = get_individual(P2Pop, r);
        if (st_sim_prepare_individual(P1Ind, P2Ind, NULL, 0) ||
            st_sim_identity_copy(P1Ind, P2Ind))
            goto fail_population;

        if (copy_genotype)
        {
            int i;

            P1Ind->genotype = (int *)st_sim_alloc(cMPhom, int);

            if (NULL == P1Ind->genotype)
            {
                genotype_alloc_error = 1;
                copy_genotype = 0;
            }
            else
            {
                for (i = 0; i < cMPhom; i++)
                    P1Ind->genotype[i] = P2Ind->genotype[i];
            }
        }

        else
        {
            P1Ind->genotype = NULL;
        }

        if (!st_copy_value_array(&(P1Ind->GValue),
                                 copy_gvalue ? P2Ind->GValue : NULL))
        {
            gvalue_alloc_error = 1;
            copy_gvalue = 0;
        }

        if (!st_copy_value_array(&(P1Ind->PValue),
                                 copy_pvalue ? P2Ind->PValue : NULL))
        {
            pvalue_alloc_error = 1;
            copy_pvalue = 0;
        }

        if (NULL != P1Ind->info)
            st_sim_free(P1Ind->info);

        if (NULL != P2Ind->info)
        {
            P1Ind->info = (char *)st_sim_pool_alloc(MAXNAME);

            if (NULL == P1Ind->info)
                info_alloc_error = 1;
            else
                strcpy(P1Ind->info, P2Ind->info);
        }

        else
        {
            P1Ind->info = NULL;
        }

        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            for (hom = 0; hom < NoHoms; hom++)
            {
                THom *P1Hom = &(P1Ind->Chrom[chrom]).Hom[hom];
                THom *P2Hom = &(P2Ind->Chrom[chrom]).Hom[hom];
                P1Hom->NoLociUsed = 0;

                for (loc = 0; loc < P2Hom->NoLociUsed; loc++)
                {
                    if (-1 == insertLocus(P2Hom->Loc[loc].All,
                                          P2Hom->Loc[loc].Pos, P1Hom))
                        goto fail_population;
                }
            }
        }
    }

    if (genotype_alloc_error)
        st_clear_population_genotype(P1Pop);

    if (gvalue_alloc_error)
        st_clear_population_gvalue(P1Pop);

    if (pvalue_alloc_error)
        st_clear_population_pvalue(P1Pop);

    free(samples);

    if (genotype_alloc_error || gvalue_alloc_error || pvalue_alloc_error ||
        info_alloc_error)
        info(-2, "Memory allocation failed in sample_population\n");

    return;

fail_population:

    free(samples);

    if (NULL != P1Pop)
        remove_population(NameP1);

    return;
}

void copy_population(char **NameP1, char **NameP2, char **cstart, char **cn)

{
    int chrom, hom, loc;
    int copy_genotype, copy_gvalue, copy_pvalue;
    int genotype_alloc_error = 0, gvalue_alloc_error = 0;
    int pvalue_alloc_error = 0, info_alloc_error = 0;
    POPULATION ind, start = atol(*cstart), n = atol(*cn);
    TInd *P1Ind, *P2Ind;
    TPop *P1Pop = NULL, *P2Pop = NULL;
    TPopListEl *PopListEl;

    find_list_element(*NameP2, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined: %s !\n", *NameP2);
        return;
    }

    P2Pop = PopListEl->Pop;

    if (0 == P2Pop->NoInds)
    {
        info(-2, "Population >%s< is empty\n", *NameP2);
        return;
    }

    if (0 == strcmp(*NameP1, *NameP2))
    {
        info(-2, "Source and destination population must be different\n");
        return;
    }

    start--;

    if (0 > start)
        start = 0;

    if (P2Pop->NoInds <= start)
    {
        info(-2, "Start index outside population %s\n", *NameP2);
        return;
    }

    if (0 > n)
        n = P2Pop->NoInds;

    if (P2Pop->NoInds < start + n)
        n = P2Pop->NoInds - start;

    copy_genotype = st_population_all_genotyped(P2Pop);
    copy_gvalue = st_population_all_gvalue(P2Pop);
    copy_pvalue = st_population_all_pvalue(P2Pop);

    P1Pop = exist_population(*NameP1, n, 1);
    if (NULL == P1Pop)
        return;

    if (n != P1Pop->NoIndsAlloc)
        return;

    P1Pop->NoInds = n;

    for (ind = 0; ind < n; ind++)
    {
        P1Ind = get_individual(P1Pop, ind);
        P2Ind = get_individual(P2Pop, start + ind);
        if (st_sim_prepare_individual(P1Ind, P2Ind, NULL, 0) ||
            st_sim_identity_copy(P1Ind, P2Ind))
            goto fail_population;

        if (copy_genotype)
        {
            int i;

            P1Ind->genotype = (int *)st_sim_alloc(cMPhom, int);

            if (NULL == P1Ind->genotype)
            {
                genotype_alloc_error = 1;
                copy_genotype = 0;
            }
            else
            {
                for (i = 0; i < cMPhom; i++)
                    P1Ind->genotype[i] = P2Ind->genotype[i];
            }
        }

        else
        {
            P1Ind->genotype = NULL;
        }

        if (!st_copy_value_array(&(P1Ind->GValue),
                                 copy_gvalue ? P2Ind->GValue : NULL))
        {
            gvalue_alloc_error = 1;
            copy_gvalue = 0;
        }

        if (!st_copy_value_array(&(P1Ind->PValue),
                                 copy_pvalue ? P2Ind->PValue : NULL))
        {
            pvalue_alloc_error = 1;
            copy_pvalue = 0;
        }

        if (NULL != P1Ind->info)
            st_sim_free(P1Ind->info);

        if (NULL != P2Ind->info)
        {
            P1Ind->info = (char *)st_sim_pool_alloc(MAXNAME);

            if (NULL == P1Ind->info)
                info_alloc_error = 1;
            else
                strcpy(P1Ind->info, P2Ind->info);
        }

        else
        {
            P1Ind->info = NULL;
        }

        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            for (hom = 0; hom < NoHoms; hom++)
            {
                THom *P1Hom = &(P1Ind->Chrom[chrom]).Hom[hom];
                THom *P2Hom = &(P2Ind->Chrom[chrom]).Hom[hom];
                P1Hom->NoLociUsed = 0;

                for (loc = 0; loc < P2Hom->NoLociUsed; loc++)
                {
                    if (-1 == insertLocus(P2Hom->Loc[loc].All,
                                          P2Hom->Loc[loc].Pos, P1Hom))
                        goto fail_population;
                }
            }
        }
    }

    if (genotype_alloc_error)
        st_clear_population_genotype(P1Pop);

    if (gvalue_alloc_error)
        st_clear_population_gvalue(P1Pop);

    if (pvalue_alloc_error)
        st_clear_population_pvalue(P1Pop);

    if (genotype_alloc_error || gvalue_alloc_error || pvalue_alloc_error ||
        info_alloc_error)
    {
        info(-2, "Memory allocation failed in copy_population\n");
        return;
    }

    info(1,
         "Copied %ld individuals (Index %ld to %ld) from population %s to "
         "population %s.\n",
         n, start + 1, start + n, *NameP2, *NameP1);
    return;

fail_population:

    if (NULL != P1Pop)
        remove_population(NameP1);

    return;
}

void append_population(char **NameP1, char **NameP2)
{
    int chrom, hom, loc;
    int copy_genotype = 0, copy_gvalue = 0, copy_pvalue = 0;
    int genotype_alloc_error = 0, gvalue_alloc_error = 0;
    int pvalue_alloc_error = 0, info_alloc_error = 0;
    int source_genotyped, source_gvalue, source_pvalue;
    int destination_existed = 0;
    POPULATION P1Inds = 0, P1Alloc = 0, NoIndsNew, ind;
    TInd *P1Ind, *P2Ind;
    TPopListEl *PopListEl;
    TPop *P1Pop = NULL, *P2Pop = NULL;

    find_list_element(*NameP2, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined: %s \n", *NameP2);
        return;
    }

    P2Pop = PopListEl->Pop;

    if (0 == P2Pop->NoInds)
    {
        info(-2, "Population >%s< is empty\n", *NameP2);
        return;
    }

    if (0 == strcmp(*NameP1, *NameP2))
    {
        info(-2, "Source and destination population must be different\n");
        return;
    }

    source_genotyped = st_population_all_genotyped(P2Pop);
    source_gvalue = st_population_all_gvalue(P2Pop);
    source_pvalue = st_population_all_pvalue(P2Pop);

    find_list_element(*NameP1, &PopListEl);

    if (NULL != PopListEl)
    {
        destination_existed = 1;
        P1Inds = PopListEl->Pop->NoInds;
        P1Alloc = PopListEl->Pop->NoIndsAlloc;
    }

    if (P2Pop->NoInds > LONG_MAX - P1Inds)
    {
        info(-2, "Resulting population is too large\n");
        return;
    }

    NoIndsNew = P1Inds + P2Pop->NoInds;

    P1Pop = exist_population(*NameP1, NoIndsNew, 0);
    if (NULL == P1Pop)
        return;

    if (NoIndsNew != P1Pop->NoIndsAlloc)
    {
        info(-2, "Population %s could not be resized\n", *NameP1);
        return;
    }

    if (0 == P1Inds)
    {
        copy_gvalue = source_gvalue;
        copy_pvalue = source_pvalue;
        copy_genotype = source_genotyped;
    }
    else
    {
        copy_gvalue = st_population_all_gvalue(P1Pop) && source_gvalue;

        copy_pvalue = st_population_all_pvalue(P1Pop) && source_pvalue;

        copy_genotype = st_population_all_genotyped(P1Pop) && source_genotyped;
    }

    for (ind = 0; ind < P2Pop->NoInds; ind++)
    {
        P1Ind = get_individual(P1Pop, ind + P1Inds);
        P2Ind = get_individual(P2Pop, ind);
        if (st_sim_prepare_individual(P1Ind, P2Ind, NULL, 0) ||
            st_sim_identity_copy(P1Ind, P2Ind))
            goto fail_population;

        if (NULL != P1Ind->genotype)
            st_sim_free(P1Ind->genotype);

        if (NULL != P1Ind->GValue)
            st_sim_free(P1Ind->GValue);

        if (NULL != P1Ind->PValue)
            st_sim_free(P1Ind->PValue);

        if (NULL != P1Ind->info)
            st_sim_free(P1Ind->info);

        if (copy_genotype)
        {
            int i;

            P1Ind->genotype = (int *)st_sim_alloc(cMPhom, int);

            if (NULL == P1Ind->genotype)
            {
                genotype_alloc_error = 1;
                copy_genotype = 0;
            }
            else
            {
                for (i = 0; i < cMPhom; i++)
                    P1Ind->genotype[i] = P2Ind->genotype[i];
            }
        }

        if (!st_copy_value_array(&(P1Ind->GValue),
                                 copy_gvalue ? P2Ind->GValue : NULL))
        {
            gvalue_alloc_error = 1;
            copy_gvalue = 0;
        }

        if (!st_copy_value_array(&(P1Ind->PValue),
                                 copy_pvalue ? P2Ind->PValue : NULL))
        {
            pvalue_alloc_error = 1;
            copy_pvalue = 0;
        }

        if (NULL != P2Ind->info)
        {
            P1Ind->info = (char *)st_sim_pool_alloc(MAXNAME);

            if (NULL == P1Ind->info)
                info_alloc_error = 1;
            else
                strcpy(P1Ind->info, P2Ind->info);
        }

        else
        {
            P1Ind->info = NULL;
        }

        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            for (hom = 0; hom < NoHoms; hom++)
            {
                THom *P1Hom = &(P1Ind->Chrom[chrom]).Hom[hom];
                THom *P2Hom = &(P2Ind->Chrom[chrom]).Hom[hom];
                P1Hom->NoLociUsed = 0;

                for (loc = 0; loc < P2Hom->NoLociUsed; loc++)
                {
                    if (-1 == insertLocus(P2Hom->Loc[loc].All,
                                          P2Hom->Loc[loc].Pos, P1Hom))
                        goto fail_population;
                }
            }
        }
    }

    P1Pop->NoInds = NoIndsNew;

    if (!copy_genotype || genotype_alloc_error)
        st_clear_population_genotype(P1Pop);

    if (!copy_gvalue || gvalue_alloc_error)
        st_clear_population_gvalue(P1Pop);

    if (!copy_pvalue || pvalue_alloc_error)
        st_clear_population_pvalue(P1Pop);

    if (genotype_alloc_error || gvalue_alloc_error || pvalue_alloc_error ||
        info_alloc_error)
    {
        info(-2, "Memory allocation failed in append_population\n");
        return;
    }

    info(1, "Population %s (%ld individuals) appended to population %s.\n",
         *NameP2, P2Pop->NoInds, *NameP1);
    return;

fail_population:

    for (ind = P1Inds; (ind < NoIndsNew) && (ind < P1Pop->NoIndsAlloc); ind++)
    {
        P1Ind = get_individual(P1Pop, ind);

        st_sim_pool_free(P1Ind->genotype);
        P1Ind->genotype = NULL;
        st_sim_pool_free(P1Ind->GValue);
        P1Ind->GValue = NULL;
        st_sim_pool_free(P1Ind->PValue);
        P1Ind->PValue = NULL;
        st_sim_pool_free(P1Ind->info);
        P1Ind->info = NULL;
        st_sim_pool_free(P1Ind->name);
        P1Ind->name = NULL;
        st_sim_pool_free(P1Ind->pedigree);
        P1Ind->pedigree = NULL;

        for (chrom = 0; chrom < NoChroms; chrom++)
            for (hom = 0; hom < NoHoms; hom++)
                P1Ind->Chrom[chrom].Hom[hom].NoLociUsed = 0;
    }

    P1Pop->NoInds = P1Inds;

    if (!destination_existed)
        remove_population(NameP1);
    else if (P1Pop->NoIndsAlloc != P1Alloc)
        int_resize_population(NameP1, P1Alloc);

    return;
}

void select_all_n_best(char **NamePop, char **indexbest, int *n,
                       double *minscore, double *maxscore,
                       char **criterion_name, char **EffName, int *retval)
{
    TPopListEl *PopListEl;
    TPop *Pop;
    TInd *Ind0;
    TInd *Ind1;
    POPULATION ind = 0;
    int criterion;
    int eff;
    int _n = *n;
    float *value0;
    float *value1;

    *retval = -2;

    criterion = st_selection_criterion(criterion_name);

    if (0 > criterion)
        return;

    if (!st_effect_position(EffName, &eff))
        return;

    find_list_element(*NamePop, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined: %s\n", *NamePop);
        return;
    }

    Pop = PopListEl->Pop;

    if (0 == Pop->NoInds)
    {
        info(-2, "Population %s is empty\n", *NamePop);
        return;
    }

    if (1 > _n)
    {
        info(-2, "Number of individuals to select must be at least 1\n");
        return;
    }

    if (!st_population_selection_ready(Pop, criterion, eff))
    {
        info(-2, (0 == criterion) ? "Population not evaluated\n"
                                  : "Population not phenotyped\n");
        return;
    }

    Ind0 = get_individual(Pop, ind);
    value0 = st_individual_value_array(Ind0, criterion);
    *maxscore = value0[eff];

    for (ind = 1; ind < Pop->NoInds; ind++)
    {
        Ind1 = get_individual(Pop, ind);
        value1 = st_individual_value_array(Ind1, criterion);

        if (0.0f != (value0[eff] - value1[eff]))
        {
            if (0 == --_n)
                break;
        }

        Ind0 = Ind1;
        value0 = value1;
    }

    *minscore = value0[eff];

    snprintf(*indexbest, strlen(*indexbest), "%ld", ind);

    info(1, "%ld individuals selected, (score %f .. %f)\n", ind, *minscore,
         *maxscore);

    *retval = 0;
}

void divide_population(char **NameP1, char **NameP2, char **cNoI)
{
    POPULATION NoI = atol(*cNoI);
    TPopListEl *PopListEl;
    TPop *P1Pop, *P2Pop = NULL;

    find_list_element(*NameP2, &PopListEl);

    if (PopListEl != NULL)
    {
        P2Pop = PopListEl->Pop;
        P1Pop = exist_population(*NameP1, NoI, 1);
        if (NULL == P1Pop)
            return;

        if (NoI > P2Pop->NoInds)
        {
            info(-1,
                 "population.divide - population size adjusted %ld -> %ld \n",
                 NoI, P2Pop->NoInds);

            NoI = P2Pop->NoInds;
        }

        {
            TInd tmpInd;
            POPULATION ind, ind_j;

            P1Pop->NoInds = NoI;

            for (ind = 0; ind < P1Pop->NoInds; ind++)
            {
                tmpInd = P1Pop->Ind[ind];
                P1Pop->Ind[ind] = P2Pop->Ind[ind];
                P2Pop->Ind[ind] = tmpInd;
            }

            for (ind = NoI, ind_j = 0; ind < P2Pop->NoInds; ind++, ind_j++)
            {
                tmpInd = P2Pop->Ind[ind];
                P2Pop->Ind[ind] = P2Pop->Ind[ind_j];
                P2Pop->Ind[ind_j] = tmpInd;
            }

            P2Pop->NoInds -= NoI;
        }
    }
    else
    {
        info(-2, "Population %s is not defined.", *NameP2);
        return;
    }
}

void rename_population(char **OldName, char **NewName)

{
    TPopListEl *PopListEl1, *PopListEl2;
    TPop *Pop = NULL;

    find_list_element(*OldName, &PopListEl1);

    if (PopListEl1 == NULL)
    {
        info(-2, "Population not defined: %s\n", *OldName);
        return;
    }

    find_list_element(*NewName, &PopListEl2);

    if (PopListEl2 != NULL)
    {
        info(-2, "Populationname >%s< is already used\n", *NewName);
        return;
    }

    if (2 > strlen(*NewName))
    {
        info(-2, "New populationname is too short\n");
        return;
    }

    if (strlen(*NewName) >= MAXNAME)
    {
        info(-2, "New populationname is too long\n");
        return;
    }

    Pop = PopListEl1->Pop;

    strcpy(Pop->Name, *NewName);
}

void swap_population_name(char **NameP1, char **NameP2)

{
    TPopListEl *PopListEl;
    find_list_element(*NameP1, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined: %s\n", *NameP1);
        return;
    }

    find_list_element(*NameP2, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population not defined:  %s\n", *NameP2);
        return;
    }

    if (0 == strcmp(*NameP1, *NameP2))
    {
        info(-2, "You have to specify two different populations\n");
        return;
    }

    {

        {
            TPopListEl *PopListEl;
            char tmpName[MAXNAME];
            TPop *P1Pop = NULL, *P2Pop = NULL;

            find_list_element(*NameP1, &PopListEl);

            if (PopListEl != NULL)
                P1Pop = PopListEl->Pop;

            if (NULL == P1Pop)
            {
                info(-2, "Population not defined: %s\n", *NameP1);
                return;
            }

            find_list_element(*NameP2, &PopListEl);

            if (PopListEl != NULL)
                P2Pop = PopListEl->Pop;

            if (NULL == P2Pop)
            {
                info(-2, "Population not defined:  %s\n", *NameP2);
                return;
            }

            strcpy(tmpName, P1Pop->Name);

            strcpy(P1Pop->Name, P2Pop->Name);

            strcpy(P2Pop->Name, tmpName);
        }
    }
}

void get_population_size(char **PopNames, char **cNoInds, char **cNoIndsAlloc)

{
    TPopListEl *PopListEl;
    POPULATION NoInds = 0, NoIndsAlloc = 0;
    char *PopName;

    PopName = strtok(*PopNames, " ");

    while (NULL != PopName)
    {
        find_list_element(PopName, &PopListEl);

        if (NULL == PopListEl)
        {
            info(-2, "Population not defined: %s\n", PopName);
            return;
        }

        NoInds += PopListEl->Pop->NoInds;

        NoIndsAlloc += PopListEl->Pop->NoIndsAlloc;

        PopName = strtok(NULL, " ");
    }

    snprintf(*cNoInds, strlen(*cNoInds), "%ld", NoInds);
    snprintf(*cNoIndsAlloc, strlen(*cNoIndsAlloc), "%ld", NoIndsAlloc);
}

void get_genome_par(int *_NoChroms, int *_NoHoms)

{
    *_NoChroms = NoChroms;
    *_NoHoms = NoHoms;
}

void get_genome_par_b(double *_ChromLen)

{
    int chrom;

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        _ChromLen[chrom] = PopDescr[chrom].ChromLen;
    }
}

void get_map_a(int *_NoMappoints, int *_maxname, char **_classlgth,
               char **_namelgth)
{
    int chrom;
    int point;
    unsigned long int class_lgth = 0;
    unsigned long int name_lgth = 0;

    *_NoMappoints = 0;
    *_maxname = MAXNAME;

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        *_NoMappoints += PopDescr[chrom].NrMappoints;

        for (point = 0; point < PopDescr[chrom].NrMappoints; point++)
        {
            name_lgth += (1 + strlen(PopDescr[chrom].Mappoint[point].name));
            class_lgth += (1 + strlen(PopDescr[chrom].Mappoint[point].class));
        }
    }

    snprintf(*_classlgth, strlen(*_classlgth), "%lu", class_lgth);

    snprintf(*_namelgth, strlen(*_namelgth), "%lu", name_lgth);
}

void get_map_b(double *_chrom, double *_pos, char **_name, char **_class)
{
    int chrom;
    int point;
    long counter = 0;
    char *p = *_name;
    char *q = *_class;

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        for (point = 0; point < PopDescr[chrom].NrMappoints; point++)
        {
            size_t len;

            _chrom[counter] = chrom + 1;
            _pos[counter] = PopDescr[chrom].Mappoint[point].Pos;

            len = strlen(PopDescr[chrom].Mappoint[point].name);
            memcpy(p, PopDescr[chrom].Mappoint[point].name, len);
            p += len;
            *p++ = ' ';

            len = strlen(PopDescr[chrom].Mappoint[point].class);
            memcpy(q, PopDescr[chrom].Mappoint[point].class, len);
            q += len;
            *q++ = ' ';

            counter++;
        }
    }
}

TMappoint *find_mappoint(char *pch)
{
    int chrom, i;

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        if (PopDescr[chrom].NrMappoints > 0)
        {
            for (i = 0; i < PopDescr[chrom].NrMappoints; i++)
            {
                if (0 == strcmp(PopDescr[chrom].Mappoint[i].name, pch))
                    return &(PopDescr[chrom].Mappoint[i]);
            }
        }
    }

    return NULL;
}

int insertmember(TSet *Set, int newmember)
{
    int i;
    int *tmp;
    size_t count;

    for (i = 0; i < Set->betrag; i++)
    {
        if (Set->member[i] == newmember)
            return 1;
    }

    if ((Set->betrag < 0) || (Set->betrag == INT_MAX))
        return 0;

    count = (size_t)Set->betrag + 1;

    if (count > SIZE_MAX / sizeof(*Set->member))
        return 0;

    tmp = (int *)realloc(Set->member, count * sizeof(*Set->member));

    if (NULL == tmp)
        return 0;

    Set->member = tmp;
    Set->member[Set->betrag] = newmember;
    Set->betrag++;

    return 1;
}

TEffects *find_Effect(char *Name)
{
    int i;

    for (i = 0; i < NrEffects; i++)
    {
        if (0 == strcmp(Name, Effects[i].Name))
            return &(Effects[i]);
    }

    return NULL;
}

void list_eff_parameter(int *lEffName, int *nrEffects)
{
    int i;
    *nrEffects = NrEffects;

    for (i = 0; i < NrEffects; i++)
    {
        *lEffName += strlen(Effects[i].Name) + 1;
    }

    (*lEffName) += 10;
}

void list_eff(char **all_eff, double *weigth_eff)
{
    int i;
    strcpy(*all_eff, "");

    for (i = 0; i < NrEffects; i++)
    {
        strcat(*all_eff, Effects[i].Name);
        strcat(*all_eff, " ");
        weigth_eff[i] = Effects[i].weight;
    }
}

static int st_effect_position(char **EffName, int *eff)
{
    TEffects *Effect;

    if (NULL == eff)
        return 0;

    *eff = 0;

    if ((NULL == EffName) || (NULL == *EffName) || ('\0' == **EffName))
        return 1;

    Effect = find_Effect(*EffName);

    if (NULL == Effect)
    {
        info(-2, "Effect %s not defined\n", *EffName);
        return 0;
    }

    *eff = Effect->Pos;

    return 1;
}

static int st_selection_criterion(char **criterion)
{
    if ((NULL == criterion) || (NULL == *criterion) ||
        (0 == strcmp(*criterion, "G")))
        return 0;

    if (0 == strcmp(*criterion, "P"))
        return 1;

    info(-2, "Selection criterion must be G or P\n");
    return -1;
}

static float *st_individual_value_array(TInd *Ind, int criterion)
{
    if (NULL == Ind)
        return NULL;

    if (0 == criterion)
        return Ind->GValue;

    return Ind->PValue;
}

static int st_population_selection_ready(TPop *Pop, int criterion, int eff)
{
    POPULATION ind;

    if ((NULL == Pop) || (0 > eff) || (eff > NrEffects))
        return 0;

    for (ind = 0; ind < Pop->NoInds; ind++)
    {
        float *value = st_individual_value_array(&(Pop->Ind[ind]), criterion);

        if ((NULL == value) || !isfinite((double)value[eff]))
            return 0;
    }

    return 1;
}

void get_population_value(char **PopName, double *value, char **EffName,
                          char **criterion_name, int *retval)
{
    int criterion;
    int eff;
    POPULATION ind;
    TPopListEl *PopListEl;

    *retval = -2;

    criterion = st_selection_criterion(criterion_name);

    if (0 > criterion)
        return;

    if (!st_effect_position(EffName, &eff))
        return;

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population %s not defined\n", *PopName);
        return;
    }

    if (!st_population_selection_ready(PopListEl->Pop, criterion, eff))
    {
        info(-2, (0 == criterion) ? "Population not evaluated\n"
                                  : "Population not phenotyped\n");
        return;
    }

    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
        TInd *Ind = get_individual(PopListEl->Pop, ind);
        float *v = st_individual_value_array(Ind, criterion);

        value[ind] = (double)v[eff];
    }

    *retval = 0;
}

void set_population_gvalue(char **PopName, double *gvalue, int *ngvalue,
                           int *retval)
{
    POPULATION ind;
    TPopListEl *PopListEl;
    float **new_value = NULL;
    int need_allocate = 0;

    *retval = -2;

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population %s not defined\n", *PopName);
        return;
    }

    if (0 == PopListEl->Pop->NoInds)
    {
        info(-2, "Population %s is empty\n", *PopName);
        return;
    }

    if ((0 > *ngvalue) || ((POPULATION)*ngvalue != PopListEl->Pop->NoInds))
    {
        info(-2, "Number of genotypic values does not match population size\n");
        return;
    }

    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
        if (!isfinite(gvalue[ind]) || (fabs(gvalue[ind]) > FLT_MAX))
        {
            info(-2, "Genotypic values must be finite float values\n");
            return;
        }
    }

    if (!st_population_all_gvalue(PopListEl->Pop))
    {
        if (0 != NrEffects)
        {
            info(-2, "First evaluate population %s\n", *PopName);
            return;
        }

        need_allocate = 1;
        new_value = (float **)calloc((size_t)PopListEl->Pop->NoInds,
                                     sizeof(*new_value));

        if (NULL == new_value)
            goto fail_alloc;

        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
            new_value[ind] = (float *)st_sim_pool_alloc(sizeof(**new_value));

            if (NULL == new_value[ind])
                goto fail_alloc;

            new_value[ind][0] = (float)gvalue[ind];
        }
    }

    if (need_allocate)
    {
        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
            TInd *Ind = get_individual(PopListEl->Pop, ind);

            st_sim_pool_free(Ind->GValue);
            Ind->GValue = new_value[ind];
            new_value[ind] = NULL;
        }
    }
    else
    {
        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
            TInd *Ind = get_individual(PopListEl->Pop, ind);

            Ind->GValue[0] = (float)gvalue[ind];
        }
    }

    st_clear_population_pvalue(PopListEl->Pop);
    free(new_value);

    *retval = 0;
    return;

fail_alloc:

    if (NULL != new_value)
    {
        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
            st_sim_pool_free(new_value[ind]);

        free(new_value);
    }

    info(-2, "Memory allocation failed in set_population_gvalue\n");
}

void phenotype_population(char **PopName, double *hsq, int *nhsq, int *retval)
{
    POPULATION ind;
    TPopListEl *PopListEl;
    float **new_value = NULL;
    double *mean = NULL;
    double *variance = NULL;
    double *sd = NULL;
    int i;

    *retval = -2;

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population %s not defined\n", *PopName);
        return;
    }

    if (2 > PopListEl->Pop->NoInds)
    {
        info(-2,
             "At least two individuals are required to generate phenotypes\n");
        return;
    }

    if (0 == NrEffects)
    {
        info(-2, "No effects defined\n");
        return;
    }

    if (*nhsq != NrEffects)
    {
        info(-2, "One heritability must be supplied for each effect\n");
        return;
    }

    if (!st_population_all_gvalue(PopListEl->Pop))
    {
        info(-2, "First evaluate population %s\n", *PopName);
        return;
    }

    for (i = 0; i < NrEffects; i++)
    {
        if (!isfinite(hsq[i]) || (0.0 >= hsq[i]) || (1.0 < hsq[i]))
        {
            info(-2, "Heritabilities must be in (0,1]\n");
            return;
        }
    }

    mean = (double *)calloc((size_t)NrEffects, sizeof(*mean));
    variance = (double *)calloc((size_t)NrEffects, sizeof(*variance));
    sd = (double *)calloc((size_t)NrEffects, sizeof(*sd));
    new_value =
        (float **)calloc((size_t)PopListEl->Pop->NoInds, sizeof(*new_value));

    if ((NULL == mean) || (NULL == variance) || (NULL == sd) ||
        (NULL == new_value))
        goto fail_alloc;

    for (i = 1; i < NrEffects + 1; i++)
    {
        double ss = 0.0;

        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
            TInd *Ind = get_individual(PopListEl->Pop, ind);

            if (!isfinite((double)Ind->GValue[i]))
            {
                info(-2, "Genotypic values must be finite\n");
                goto fail;
            }

            mean[i - 1] += (double)Ind->GValue[i];
        }

        mean[i - 1] /= (double)PopListEl->Pop->NoInds;

        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
            TInd *Ind = get_individual(PopListEl->Pop, ind);
            double d = (double)Ind->GValue[i] - mean[i - 1];

            ss += d * d;
        }

        variance[i - 1] = ss / (double)(PopListEl->Pop->NoInds - 1);

        if (!isfinite(variance[i - 1]) || (0.0 > variance[i - 1]))
        {
            info(-2, "Invalid genotypic variance\n");
            goto fail;
        }

        sd[i - 1] = sqrt(variance[i - 1] / hsq[i - 1] - variance[i - 1]);

        if (!isfinite(sd[i - 1]))
        {
            info(-2, "Invalid residual standard deviation\n");
            goto fail;
        }
    }

    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
        if ((size_t)NrEffects + 1 > SIZE_MAX / sizeof(**new_value))
            goto fail_alloc;

        new_value[ind] = (float *)st_sim_alloc_array(
            (size_t)NrEffects + 1, sizeof(**new_value));

        if (NULL == new_value[ind])
            goto fail_alloc;
    }

    GetRNGstate();

    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
        TInd *Ind = get_individual(PopListEl->Pop, ind);
        double index = 0.0;

        for (i = 1; i < NrEffects + 1; i++)
        {
            double p = (double)Ind->GValue[i] + norm_rand() * sd[i - 1];

            if (!isfinite(p) || (fabs(p) > FLT_MAX))
            {
                PutRNGstate();
                info(-2, "Generated phenotypic value is not finite\n");
                goto fail;
            }

            new_value[ind][i] = (float)p;
            index += p * (double)Effects[i - 1].weight;
        }

        if (!isfinite(index) || (fabs(index) > FLT_MAX))
        {
            PutRNGstate();
            info(-2, "Phenotypic selection index is not finite\n");
            goto fail;
        }

        new_value[ind][0] = (float)index;
    }

    PutRNGstate();

    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
        TInd *Ind = get_individual(PopListEl->Pop, ind);

        st_sim_pool_free(Ind->PValue);
        Ind->PValue = new_value[ind];
        new_value[ind] = NULL;
    }

    free(new_value);
    free(mean);
    free(variance);
    free(sd);

    *retval = 0;
    return;

fail_alloc:

    info(-2, "Memory allocation failed in phenotype_population\n");

fail:

    if (NULL != new_value)
    {
        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
            st_sim_pool_free(new_value[ind]);
    }

    free(new_value);
    free(mean);
    free(variance);
    free(sd);
}

void set_population_info(char **PopName, char **_ind, char **_info)
{
    POPULATION ind = atol(*_ind) - 1;
    TPopListEl *PopListEl;
    TInd *Ind;

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population %s not defined\n", *PopName);
        return;
    }

    Ind = (TInd *)get_individual(PopListEl->Pop, ind);

    if (strlen(*_info) >= MAXNAME)
    {
        info(-2, "Info string is too long\n");
        return;
    }

    if (NULL == Ind->info)
    {
        Ind->info = st_sim_alloc(MAXNAME, char);
        if (NULL == Ind->info)
        {
            info(-2, "Memory allocation failed in set_population_info\n");
            return;
        }
    }

    strcpy(Ind->info, *_info);
}

void get_population_info(char **PopName, char **_ind, char **_info)
{
    POPULATION ind = atol(*_ind) - 1;
    TPopListEl *PopListEl;
    TInd *Ind;

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Population %s not defined\n", *PopName);
        return;
    }

    Ind = (TInd *)get_individual(PopListEl->Pop, ind);

    if (NULL == Ind->info)
        strcpy(*_info, "none");
    else
        strcpy(*_info, Ind->info);
}

void set_eff_weight(char **Name, double *weight)
{
    TEffects *Effect = find_Effect(*Name);

    if (NULL == Effect)
    {
        info(-2, "Can't find effect %s\n", *Name);
        return;
    }

    Effect->weight = (float)*weight;
    st_clear_all_population_gvalue();
}

void effect_weight_set_all(double *weight)
{
    int i;

    for (i = 0; i < NrEffects; i++)
    {
        Effects[i].weight = *weight;
    }

    st_clear_all_population_gvalue();
}

void Insert_ALL_Loc(TEff *Eff, TMappoint *MP, int All)
{
    int i;

    for (i = 0; i < Eff->NrLoci; i++)
    {
        if (MP == Eff->Loc[i].MP)
        {
            if (Eff->Loc[i].NrAll++ < NoHoms)
            {
                Eff->Loc[i].All[Eff->Loc[i].NrAll - 1] = All;
            }

            else
            {
                Eff->Loc[i].NrAll--;
                info(-2, "Error in the effect-file - too many alleles\n");
            }

            return;
        }
    }

    if (0 == Eff->NrLoci++)
        Eff->Loc = (TGen *)yCalloc(1, TGen);
    else
        Eff->Loc = (TGen *)yRealloc(Eff->Loc, Eff->NrLoci, TGen);

    Eff->Loc[Eff->NrLoci - 1].MP = MP;
    Eff->Loc[Eff->NrLoci - 1].NrAll = 1;
    Eff->Loc[Eff->NrLoci - 1].All = (int *)yCalloc(NoHoms, int);
    Eff->Loc[Eff->NrLoci - 1].All[0] = All;

    for (i = 1; i < NoHoms; i++)
        Eff->Loc[Eff->NrLoci - 1].All[i] = -1;
}

void load_effmap_mod(char **effectName, char **effectSpec, double *effectWeight)

{
    char line[10 * MAXLINE], *token;
    int i, j, k, tmpi, status, maxLoc = 0;
    float tmpf;
    TMappoint *tmpMP = NULL;
    TEffects *Effect;

    if (NULL != find_Effect(*effectName))
    {
        info(-2, "This effect is already loaded (same name)\n");
        return;
    }

    if (strlen(*effectSpec) >= sizeof(line))
    {
        info(-2, "Effect specification too long\n");
        return;
    }

    if (0 == NrEffects++)
        Effects = (TEffects *)yCalloc(1, TEffects);
    else
        Effects = (TEffects *)yRealloc(Effects, NrEffects, TEffects);

    Effect = &Effects[NrEffects - 1];
    Effect->GValue_avg = 0;
    Effect->weight = 1;
    Effect->NrEff = 0;
    Effect->Pos = NrEffects;
    Effect->Eff = NULL;
    Effect->MPset = NULL;
    strcpy(Effect->Name, *effectName);
    Effect->weight = (float)*effectWeight;
    Effect->GValue_avg = 0;

    for (i = 0; i < 1; i++)
    {
        strcpy(line, *effectSpec);
        status = 0;
        token = strtok(line, "[ \t]");

        while ((NULL != token) && (status != -1) && (status != 1))
        {
            switch (status)
            {

            case 0:

                if (0 == strcmp(line, "\n"))
                    status = -1;
                else
                {
                    status = 2;

                    if (0 == (Effect->NrEff)++)
                        Effect->Eff = (TEff *)yCalloc(1, TEff);
                    else
                        Effect->Eff =
                            (TEff *)yRealloc(Effect->Eff, Effect->NrEff, TEff);

                    Effect->Eff[Effect->NrEff - 1].Eff = 0;
                    Effect->Eff[Effect->NrEff - 1].NrLoci = 0;
                    Effect->Eff[Effect->NrEff - 1].Loc = NULL;
                }

                break;

            case 2:

                if (NULL != (tmpMP = find_mappoint(token)))
                {
                    status = 3;
                }

                else
                {
                    status = -1;
                    info(-1, "Marker %s not found.\nline %d\n", token, i);
                }

                break;

            case 3:

                if (1 == sscanf(token, "%d", &tmpi))
                {
                    status = 4;
                    Insert_ALL_Loc(&Effect->Eff[Effect->NrEff - 1], tmpMP,
                                   tmpi);
                }

                else
                {
                    status = -1;
                    info(-1, "Integer expected\nline %d\n", i);
                }

                break;

            case 4:

                if (1 == sscanf(token, "%f", &tmpf))
                {
                    status = 1;
                    Effect->Eff[Effect->NrEff - 1].Eff = tmpf;

                    if (maxLoc < Effect->Eff[Effect->NrEff - 1].NrLoci)
                    {
                        maxLoc = Effect->Eff[Effect->NrEff - 1].NrLoci;
                    }
                }

                else
                {
                    status = 2;
                }

                break;
            }

            if (2 != status)
                token = strtok(NULL, "[ \t]");
        }
    }

    if (1 != status)
    {
        if (-1 != status)
            info(-1, "Incomplete effect specification\n");
        goto fail;
    }

    Effect->MPs = 0;

    for (i = 0; i < Effect->NrEff; i++)
    {
        for (j = 0; j < Effect->Eff[i].NrLoci; j++)
        {
            for (k = 0; k < Effect->MPs; k++)
            {
                if (*(Effect->MPset + k) == (Effect->Eff[i].Loc[j].MP))
                    break;
            }

            if (k == Effect->MPs)
            {
                Effect->MPs++;

                if (1 == Effect->MPs)
                    Effect->MPset = (TMappoint **)yCalloc(1, TMappoint *);
                else
                    Effect->MPset = (TMappoint **)yRealloc(
                        Effect->MPset, Effect->MPs, TMappoint *);

                Effect->MPset[Effect->MPs - 1] = Effect->Eff[i].Loc[j].MP;
            }
        }
    }

    st_clear_all_population_gvalue();
    return;

fail:

    for (i = 0; i < Effect->NrEff; i++)
    {
        for (j = 0; j < Effect->Eff[i].NrLoci; j++)
            yFree(Effect->Eff[i].Loc[j].All);

        yFree(Effect->Eff[i].Loc);
    }

    yFree(Effect->Eff);
    yFree(Effect->MPset);

    NrEffects--;

    if (0 == NrEffects)
        yFree(Effects);
    else
    {
        TEffects *tmp = realloc(Effects, NrEffects * sizeof(*tmp));
        if (NULL != tmp)
            Effects = tmp;
    }
}

void load_effmap(char **file, char **name)
{
    FILE *fp;
    char line[10 * MAXLINE], *token;
    int i, j, k, tmpi, status, maxLoc = 0;
    float tmpf;
    TMappoint *tmpMP = NULL;
    TEffects *Effect;

    if (NULL != find_Effect(*name))
    {
        info(-2, "This effect is already loaded (same name)\n");
        return;
    }

    fp = fileopen(*file, "r");

    if (NULL == fp)
    {
        info(-2, "Effec file can not be opened\n");
        return;
    }

    if (0 == NrEffects++)
        Effects = (TEffects *)yCalloc(1, TEffects);
    else
        Effects = (TEffects *)yRealloc(Effects, NrEffects, TEffects);

    Effect = &Effects[NrEffects - 1];

    Effect->GValue_avg = 0;

    Effect->weight = 1;

    Effect->NrEff = 0;

    Effect->Pos = NrEffects;

    Effect->Eff = NULL;

    Effect->MPset = NULL;

    strcpy(Effect->Name, *name);

    for (i = 0;; i++)
    {
        if (fgets(line, 10 * MAXLINE, fp) == NULL)
            break;

        if ((unsigned)(j = strcspn(line, "#")) != strlen(line))
        {
            if (j == 0)
                strcpy(line, "");
            else
            {
                line[j++] = '\n';
                line[j] = '\0';
            }
        }

        status = 0;
        token = strtok(line, "[ \t]");

        while ((NULL != token) && (status != -1) && (status != 1))
        {
            switch (status)
            {

            case 0:

                if (0 == strcmp(token, "Weigth:"))
                {
                    status = 5;
                }

                else if (1 == sscanf(token, "%f", &Effect->GValue_avg))
                {
                    status = 1;
                }

                else
                {
                    if (0 == strcmp(line, "\n"))
                        status = -1;
                    else
                    {
                        status = 2;

                        if (0 == (Effect->NrEff)++)
                            Effect->Eff = (TEff *)yCalloc(1, TEff);
                        else
                            Effect->Eff = (TEff *)yRealloc(Effect->Eff,
                                                           Effect->NrEff, TEff);

                        Effect->Eff[(Effect->NrEff) - 1].Eff = 0;

                        Effect->Eff[(Effect->NrEff) - 1].NrLoci = 0;
                        Effect->Eff[(Effect->NrEff) - 1].Loc = NULL;
                    }
                }

                break;

            case 2:

                if (NULL != (tmpMP = find_mappoint(token)))
                {
                    status = 3;
                }

                else
                {
                    status = -1;
                    info(-1, "Marker %s not found.\nline %d\n", token, i);
                }

                break;

            case 3:

                if (1 == sscanf(token, "%d", &tmpi))
                {
                    status = 4;
                    Insert_ALL_Loc(&Effect->Eff[Effect->NrEff - 1], tmpMP,
                                   tmpi);
                }

                else
                {
                    status = -1;
                    info(-1, "Integer expected\nline %d\n", i);
                }

                break;

            case 4:

                if (1 == sscanf(token, "%f", &tmpf))
                {
                    status = 1;
                    Effect->Eff[Effect->NrEff - 1].Eff = tmpf;

                    if (maxLoc < Effect->Eff[Effect->NrEff - 1].NrLoci)
                    {
                        maxLoc = Effect->Eff[Effect->NrEff - 1].NrLoci;
                    }
                }

                else
                {
                    status = 2;
                }

                break;

            case 5:

                if (1 == sscanf(token, "%f", &Effect->weight))
                {
                    status = 1;
                }

                else
                {
                    status = -1;
                    info(-1, "Error: weight");
                }

                break;
            }

            if (2 != status)
                token = strtok(NULL, "[ \t]");
        }
    }

    fileclose(fp);

    Effect->MPs = 0;

    for (i = 0; i < Effect->NrEff; i++)
    {
        for (j = 0; j < Effect->Eff[i].NrLoci; j++)
        {
            for (k = 0; k < Effect->MPs; k++)
            {
                if (*(Effect->MPset + k) == (Effect->Eff[i].Loc[j].MP))
                    break;
            }

            if (k == Effect->MPs)
            {
                Effect->MPs++;

                if (1 == Effect->MPs)
                    Effect->MPset = (TMappoint **)yCalloc(1, TMappoint *);
                else
                    Effect->MPset = (TMappoint **)yRealloc(
                        Effect->MPset, Effect->MPs, TMappoint *);

                Effect->MPset[Effect->MPs - 1] = Effect->Eff[i].Loc[j].MP;
            }
        }
    }

    {
        int l, m;

        for (i = 0; i < maxLoc; i++)
        {
            for (k = 0; k < NoHoms; k++)
            {
                int oldEff = Effect->NrEff;
                int n;

                for (n = k + 1; n < NoHoms; n++)
                {
                    for (j = 0; j < oldEff; j++)
                    {
                        if ((i < Effect->Eff[j].NrLoci) &&
                            (k < Effect->Eff[j].Loc[i].NrAll) &&
                            (Effect->Eff[j].Loc[i].All[k] !=
                             Effect->Eff[j].Loc[i].All[n]))
                        {
                            Effect->NrEff++;
                            Effect->Eff = (TEff *)yRealloc(Effect->Eff,
                                                           Effect->NrEff, TEff);
                            Effect->Eff[Effect->NrEff - 1].Eff =
                                Effect->Eff[j].Eff;
                            Effect->Eff[Effect->NrEff - 1].NrLoci =
                                Effect->Eff[j].NrLoci;
                            Effect->Eff[Effect->NrEff - 1].Loc =
                                (TGen *)yCalloc(Effect->Eff[j].NrLoci, TGen);

                            for (l = 0; l < Effect->Eff[j].NrLoci; l++)
                            {
                                Effect->Eff[Effect->NrEff - 1].Loc[l].MP =
                                    Effect->Eff[j].Loc[l].MP;
                                Effect->Eff[Effect->NrEff - 1].Loc[l].NrAll =
                                    Effect->Eff[j].Loc[l].NrAll;
                                Effect->Eff[Effect->NrEff - 1].Loc[l].All =
                                    (int *)yCalloc(NoHoms, int);

                                for (m = 0; m < NoHoms; m++)
                                {
                                    if ((l == i) && (m == k))
                                    {
                                        Effect->Eff[Effect->NrEff - 1]
                                            .Loc[l]
                                            .All[m] =
                                            Effect->Eff[j].Loc[l].All[n];
                                    }

                                    else if ((l == i) && (m == n))
                                    {
                                        Effect->Eff[Effect->NrEff - 1]
                                            .Loc[l]
                                            .All[m] =
                                            Effect->Eff[j].Loc[l].All[k];
                                    }

                                    else
                                    {
                                        Effect->Eff[Effect->NrEff - 1]
                                            .Loc[l]
                                            .All[m] =
                                            Effect->Eff[j].Loc[l].All[m];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    st_clear_all_population_gvalue();
}

void permutate_effmap(void)
{
}

void effmap_remove(char **effectName)
{
    if (NULL == find_Effect(*effectName))
    {
        info(-2, "Couldn't find effect: %s\n", *effectName);
        return;
    }

    int i, j, k;

    for (i = NrEffects - 1; i >= 0; i--)
    {
        if (0 == strcmp(*effectName, Effects[i].Name))
        {

            for (j = 0; j < Effects[i].NrEff; j++)
            {
                for (k = 0; k < Effects[i].Eff[j].NrLoci; k++)
                    if (NULL != Effects[i].Eff[j].Loc[k].All)
                        yFree(Effects[i].Eff[j].Loc[k].All);

                yFree(Effects[i].Eff[j].Loc);
            }

            if (NULL != Effects[i].Eff)
                yFree(Effects[i].Eff);

            if (NULL != Effects[i].MPset)
                yFree(Effects[i].MPset);

            for (j = i + 1; j < NrEffects; j++)
            {
                Effects[j - 1] = Effects[j];
                Effects[j - 1].Pos = j;
            }

            NrEffects--;

            if (0 == NrEffects)
                yFree(Effects);
            else
                Effects = (TEffects *)yRealloc(Effects, NrEffects, TEffects);

            st_clear_all_population_gvalue();
            info(1, "Effect %s removed\n", *effectName);
            break;
        }
    }
}

void effmap_remove_all(void)
{
    int i, j, k;

    for (i = 0; i < NrEffects; i++)
    {
        for (j = 0; j < Effects[i].NrEff; j++)
        {
            for (k = 0; k < Effects[i].Eff[j].NrLoci; k++)
            {
                if (NULL != Effects[i].Eff[j].Loc[k].All)
                {
                    yFree(Effects[i].Eff[j].Loc[k].All);
                }
            }

            yFree(Effects[i].Eff[j].Loc);
        }

        if (NULL != Effects[i].Eff)
            yFree(Effects[i].Eff);
        if (NULL != Effects[i].MPset)
            yFree(Effects[i].MPset);
    }

    if (NULL != Effects)
        yFree(Effects);

    if (0 < NrEffects)
        info(1, "%d effects removed\n", NrEffects);

    NrEffects = 0;
    st_clear_all_population_gvalue();
}

void remove_genotype_population(char **PopNames)
{
    char *PopName;

    PopName = strtok(*PopNames, "[' ','\n']");

    while (NULL != PopName)
    {
        TPopListEl *PopListEl;

        find_list_element(PopName, &PopListEl);

        if (NULL == PopListEl)
            info(-1, "Population not defined: %s\n", PopName);
        else
        {
            POPULATION ind;

            for (ind = 0; ind < PopListEl->Pop->NoIndsAlloc; ind++)
            {
                if (NULL != PopListEl->Pop->Ind[ind].genotype)
                {
                    st_sim_pool_free(PopListEl->Pop->Ind[ind].genotype);

                    PopListEl->Pop->Ind[ind].genotype = NULL;
                }
            }
            st_clear_population_pvalue(PopListEl->Pop);

            info(1, "Genotype removed for population %s\n", PopName);
        }

        PopName = strtok(NULL, "[' ','\n']");
    }
}

void genotype_population(char **PopNames)
{
    char *PopName;

    TPopListEl *PopListEl;
    int alloc_error;

    if (0 == cMP)
    {
        info(-2, "No linkage map defined\n");
        return;
    }

    PopName = strtok(*PopNames, "[' ','\n']");

    while (NULL != PopName)
    {

        find_list_element(PopName, &PopListEl);

        if (NULL == PopListEl)
            info(-1, "Population not defined: %s\n", PopName);
        else if (0 == PopListEl->Pop->NoInds)
            info(-1, "Population %s is empty\n", PopName);
        else
        {
            if (st_population_all_genotyped(PopListEl->Pop))
            {
                info(1,
                     "Population %s already genotyped - genotyping skipped\n",
                     PopName);
            }
            else
            {

                alloc_error = 0;

#ifdef MT
#pragma omp parallel
#endif

                {

                    POPULATION ind;
                    int chrom, hom;

                    TInd *Ind;
                    TChrom *Chrom;
                    int i, m;

#ifdef MT
#pragma omp for
#endif

                    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
                    {
                        int *tmp;

                        Ind = get_individual(PopListEl->Pop, ind);

                        tmp = (int *)st_sim_resize_array(
                            Ind->genotype, (size_t)cMPhom, sizeof(*tmp));

                        if (NULL == tmp)
                        {
#ifdef MT
#pragma omp atomic write
#endif
                            alloc_error = 1;
                            continue;
                        }

                        Ind->genotype = tmp;

                        m = 0;

                        for (chrom = 0; chrom < NoChroms; chrom++)
                        {

                            for (i = 0; i < PopDescr[chrom].NrMappoints;
                                 i++, m++)
                            {

                                Chrom = &(Ind->Chrom[chrom]);

                                for (hom = 0; hom < NoHoms; hom++)
                                {
                                    Ind->genotype[m * NoHoms + hom] =
                                        get_allele(
                                            &(Chrom->Hom[hom]),
                                            PopDescr[chrom].Mappoint[i].Pos);
                                }
                            }
                        }
                    }
                }

                if (alloc_error)
                {
                    st_clear_population_genotype(PopListEl->Pop);
                    info(-2,
                         "Memory allocation failed in genotype_population\n");
                    return;
                }

                info(1, "Population %s genotyped\n", PopName);
            }
        }

        PopName = strtok(NULL, "[' ','\n']");
    }
}

void evaluate_population(char **PopNames, char **EffName)
{
    char *PopName;
    int eff = -1;
    int alloc_error = 0;

    TPopListEl *PopListEl;

    if (0 == NrEffects)
    {
        info(-2, "No effects defined\n");
        return;
    }

    PopName = strtok(*PopNames, "[' ','\n']");

    while (NULL != PopName)
    {
        find_list_element(PopName, &PopListEl);

        if (NULL == PopListEl)
            info(-1, "Population not defined: %s\n", PopName);
        else
        {
            if (0 == PopListEl->Pop->NoInds)
            {
                info(-2, "Population %s is empty\n", PopName);
                return;
            }

            if (!st_population_all_genotyped(PopListEl->Pop))
            {
                info(-2, "First genotype population %s\n", PopName);
                return;
            }

            if (NULL != EffName)
            {
                TEffects *Effect;

                if (NULL == (Effect = find_Effect(*EffName)))
                {
                    info(-2, "Can't find the effect %s\n", *EffName);
                    return;
                }

                eff = Effect->Pos;
            }
            st_clear_population_pvalue(PopListEl->Pop);

            alloc_error = 0;

#ifdef MT
#pragma omp parallel
#endif

            {

                POPULATION ind;
                TInd *Ind;
                int i, j, hom, loc, status;

#ifdef MT
#pragma omp for
#endif

                for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
                {
                    float *tmp;

                    Ind = (TInd *)get_individual(PopListEl->Pop, ind);

                    tmp = (float *)st_sim_resize_array(
                        Ind->GValue, (size_t)NrEffects + 1, sizeof(*tmp));

                    if (NULL == tmp)
                    {
#ifdef MT
#pragma omp atomic write
#endif
                        alloc_error = 1;
                        continue;
                    }

                    Ind->GValue = tmp;

                    Ind->GValue[0] = 0;

                    for (i = 1; i < NrEffects + 1; i++)
                    {
                        Ind->GValue[i] = Effects[i - 1].GValue_avg;

                        for (j = 0; j < Effects[i - 1].NrEff; j++)
                        {
                            status = 0;
                            loc = 0;

                            while ((0 == status) &&
                                   (loc < Effects[i - 1].Eff[j].NrLoci))
                            {
                                hom = 0;

                                while ((0 == status) && (hom < NoHoms))
                                {
                                    if ((-1 < Effects[i - 1]
                                                  .Eff[j]
                                                  .Loc[loc]
                                                  .All[hom]) &&
                                        (Effects[i - 1]
                                             .Eff[j]
                                             .Loc[loc]
                                             .All[hom] !=
                                         Ind->genotype[Effects[i - 1]
                                                               .Eff[j]
                                                               .Loc[loc]
                                                               .MP->counter *
                                                           NoHoms +
                                                       hom]))
                                    {
                                        status = 1;
                                    }

                                    hom++;
                                }

                                loc++;
                            }

                            if (0 == status)
                                Ind->GValue[i] += Effects[i - 1].Eff[j].Eff;
                        }

                        Ind->GValue[0] +=
                            Ind->GValue[i] * Effects[i - 1].weight;
                    }

                    if (eff > -1)
                        Ind->GValue[0] = Ind->GValue[eff];
                }
            }

            if (alloc_error)
            {
                st_clear_population_gvalue(PopListEl->Pop);
                info(-2, "Memory allocation failed in evaluate_population\n");
                return;
            }

            info(1, "Population %s evaluated; selection index: %s\n", PopName,
                 (-1 == eff) ? "weigthed sum" : *EffName);
        }

        PopName = strtok(NULL, "[' ','\n']");
    }
}

void remove_evaluate_population(char **PopNames)
{
    char *PopName;

    PopName = strtok(*PopNames, "[' ','\n']");

    while (NULL != PopName)
    {
        TPopListEl *PopListEl;

        find_list_element(PopName, &PopListEl);

        if (NULL == PopListEl)
            info(-1, "Population not defined: %s\n", PopName);
        else
        {
            st_clear_population_gvalue(PopListEl->Pop);

            info(1, "Evaluation removed for population %s\n", PopName);
        }

        PopName = strtok(NULL, "[' ','\n']");
    }
}

static int st_sort_criterion = 0;
static int st_sort_eff = 0;
static int st_sort_decreasing = 1;

static int st_indcmp(const void *a, const void *b)
{
    const TInd *ia = (const TInd *)a;
    const TInd *ib = (const TInd *)b;
    const float *va = (0 == st_sort_criterion) ? ia->GValue : ia->PValue;
    const float *vb = (0 == st_sort_criterion) ? ib->GValue : ib->PValue;

    if (va[st_sort_eff] > vb[st_sort_eff])
        return st_sort_decreasing ? -1 : 1;

    if (va[st_sort_eff] < vb[st_sort_eff])
        return st_sort_decreasing ? 1 : -1;

    return 0;
}

void sort_population(char **PopNames, int *decreasing, char **criterion_name,
                     char **EffName, int *retval)
{
    TPopListEl *PopListEl;
    char *PopName;
    int criterion;
    int eff;

    *retval = -2;

    criterion = st_selection_criterion(criterion_name);

    if (0 > criterion)
        return;

    if (!st_effect_position(EffName, &eff))
        return;

    PopName = strtok(*PopNames, " ");

    while (NULL != PopName)
    {
        find_list_element(PopName, &PopListEl);

        if (NULL == PopListEl)
        {
            info(-2, "Population not defined: %s\n", PopName);
            return;
        }

        if (0 == PopListEl->Pop->NoInds)
        {
            info(-2, "Population %s is empty\n", PopName);
            return;
        }

        if (!st_population_selection_ready(PopListEl->Pop, criterion, eff))
        {
            info(-2,
                 (0 == criterion) ? "First evaluate population %s\n"
                                  : "First phenotype population %s\n",
                 PopName);
            return;
        }

        st_sort_criterion = criterion;
        st_sort_eff = eff;
        st_sort_decreasing = *decreasing ? 1 : 0;

        qsort(PopListEl->Pop->Ind, PopListEl->Pop->NoInds, sizeof(TInd),
              st_indcmp);

        PopName = strtok(NULL, " ");
    }

    *retval = 0;
}

int getPos(int all, TSet Set)
{
    int i;

    for (i = 0; i < Set.betrag; i++)
    {
        if (Set.member[i] == all)
        {
            return i;
        }
    }

    info(-2, "Fehler in getPos()\n");
    return 0;
}

typedef struct
{
    TSet *sets;
    int markers;
    char *names;
} st_population_export_workspace;
static void st_population_export_destroy(void *payload)
{
    st_population_export_workspace *workspace = payload;
    if (workspace->sets)
        for (int i = 0; i < workspace->markers; i++) free(workspace->sets[i].member);
    free(workspace->sets);free(workspace->names);free(workspace);
}

static SEXP st_legacy_argument(SEXP args, const char *key);

static SEXP st_mab_io_keep(SEXP value, int *np)
{
    PROTECT(value);
    ++*np;
    return value;
}

static void st_population_mab_export(SEXP names, SEXP args, int background,
                                      int nalleles, int nindividuals)
{
    const char **prefix = st_work_alloc(nalleles, sizeof(*prefix)),
               **label = st_work_alloc(nalleles, sizeof(*label));
    TSet *sets = AllSet;
    int row = 0;
    for (int chrom = 0, marker = 0; chrom < NoChroms; chrom++)
        for (int m = 0; m < PopDescr[chrom].NrMappoints; m++, marker++)
            for (int a = 0; a < sets[marker].betrag; a++)
            {
                char buffer[2 * MAXNAME];
                snprintf(buffer, sizeof(buffer), "%s.%d ", PopDescr[chrom].Mappoint[m].name, sets[marker].member[a]);
                const char *dot = strchr(buffer, '.');
                if (!dot) dot = strchr(buffer, '_');
                size_t n = dot ? (size_t)(dot - buffer) : strlen(buffer);
                char *part = st_work_alloc(n + 1, 1);
                memcpy(part, buffer, n); part[n] = '\0'; prefix[row] = part;
                const char *rest = dot ? dot + 1 : "";
                char *copy = st_work_alloc(strlen(rest) + 1, 1);
                strcpy(copy, rest); label[row++] = copy;
            }
    int nr = 0;
    for (int i = 0; i < nalleles;)
    {
        i += i + 1 < nalleles && !strcmp(prefix[i], prefix[i + 1]) ? 2 : 1;
        nr++;
    }
    TInd **individuals = st_work_alloc(nindividuals, sizeof(*individuals));
    const char *text = CHAR(STRING_ELT(names, 0));
    char *copy = st_work_alloc(strlen(text) + 1, 1); strcpy(copy, text);
    char *population = strtok(copy, " "); int at = 0;
    while (population)
    {
        TPopListEl *entry = NULL; find_list_element(population, &entry);
        for (POPULATION i = 0; entry && i < entry->Pop->NoInds; i++)
            individuals[at++] = get_individual(entry->Pop, i);
        population = strtok(NULL, " ");
    }
    int **columns = st_work_alloc(nindividuals, sizeof(*columns));
    for (int j = 0; j < nindividuals; j++)
    {
        columns[j] = st_work_alloc(nalleles, sizeof(int));
        memset(columns[j], 0, (size_t)nalleles * sizeof(int));
        for (int m = 0, offset = 0; m < cMP; offset += sets[m++].betrag)
            for (int h = 0; h < NoHoms; h++)
            {
                int allele = individuals[j]->genotype[m * NoHoms + h];
                if (allele != background) columns[j][offset + getPos(allele, sets[m])] = 1;
            }
    }
    int np = 0;
    SEXP utils = st_mab_io_keep(R_FindNamespace(st_mab_io_keep(mkString("utils"), &np)), &np),
         path_function = st_legacy_argument(args, "path.function"),
         directory = st_legacy_argument(args, "directory"),
         name = st_legacy_argument(args, "name"),
         file_prefix = st_legacy_argument(args, "prefix");
    const char *first = NULL, *second = NULL;
    for (int j = 0; j < (nindividuals ? nindividuals : 2); j++)
    {
        const int *presence = nindividuals ? columns[j] : NULL;
        int before = np;
        SEXP matrix = nalleles ? st_mab_io_keep(allocMatrix(STRSXP, nr, 3), &np) : R_NilValue;
        int r = 0;
        for (int i = 0; i < nalleles;)
        {
            int paired = i + 1 < nalleles && !strcmp(prefix[i], prefix[i + 1]);
            if (paired)
            {
                if (presence[i] && !presence[i + 1]) first = second = label[i];
                else if (presence[i] && presence[i + 1]) first = label[i], second = label[i + 1];
                else if (!presence[i] && presence[i + 1]) first = second = label[i + 1];
                if (!first) error("object 'mgt' not found");
            }
            else first = second = label[i];
            SET_STRING_ELT(matrix, r, mkChar(prefix[i]));
            SET_STRING_ELT(matrix, r + nr, mkChar(first));
            SET_STRING_ELT(matrix, r + 2 * nr, mkChar(second));
            i += paired ? 2 : 1; r++;
        }
        char digits[32]; snprintf(digits, sizeof(digits), "%03d", nindividuals ? j + 1 : j ? 0 : 1);
        SEXP suffix = st_mab_io_keep(mkString(digits), &np),
             empty = st_mab_io_keep(mkString(""), &np),
             call = st_mab_io_keep(lang4(findFun(install("paste"), R_BaseEnv), name, suffix, empty), &np);
        SET_TAG(CDDDR(call), install("sep"));
        SEXP iname = st_mab_io_keep(eval(call, R_BaseEnv), &np);
        call = st_mab_io_keep(lang5(findFun(install("paste"), R_BaseEnv), file_prefix, suffix,
                                    st_mab_io_keep(mkString(".mda"), &np), empty), &np);
        SET_TAG(CDR(CDDDR(call)), install("sep"));
        SEXP fname = st_mab_io_keep(eval(call, R_BaseEnv), &np);
        call = st_mab_io_keep(lang3(path_function, directory, fname), &np);
        SEXP path = st_mab_io_keep(eval(call, R_BaseEnv), &np),
             no = st_mab_io_keep(ScalarLogical(FALSE), &np),
             yes = st_mab_io_keep(ScalarLogical(TRUE), &np);
        call = st_mab_io_keep(lang4(findFun(install("write"), R_BaseEnv), iname, path, no), &np);
        SET_TAG(CDDR(call), install("file")); SET_TAG(CDDDR(call), install("append"));
        eval(call, R_BaseEnv);
        SEXP tail = st_mab_io_keep(allocList(6), &np);
        call = st_mab_io_keep(lcons(findFun(install("write.table"), utils), tail), &np);
        SEXP values[] = {matrix, path, no, no, no, yes};
        const char *keys[] = {"", "file", "quote", "row.names", "col.names", "append"};
        for (int i = 0; i < 6; i++, tail = CDR(tail))
        {
            SETCAR(tail, values[i]); if (i) SET_TAG(tail, install(keys[i]));
        }
        eval(call, R_BaseEnv);
        UNPROTECT(np - before); np = before;
    }
    UNPROTECT(np);
}

static SEXP st_population_export(SEXP _PopNames, SEXP _bg, SEXP mab, int *success)
{
    int i_bg = 0;
    int nprotect = 0;
    int chrom, hom, marker, cmarker, allele;
    int NoAll = 0, NoInds = 0;
    R_xlen_t col = 0, row = 0;
    POPULATION ind;
    TPopListEl *PopListEl;
    char *PopNames = NULL;
    char *PopName = NULL;
    char tmp[2 * MAXNAME];
    int *bg = &i_bg;
    size_t pop_names_length;
    SEXP result = R_NilValue;
    SEXP matrix, xAxis, yAxis, cls;
    st_population_export_workspace *workspace = NULL;

    if (!isString(_PopNames) || length(_PopNames) != 1)
    {
        info(-2, "Name of populations is not a single string.\n");
        return R_NilValue;
    }

    if (!isInteger(_bg) || length(_bg) != 1)
    {
        info(-2, "Value of background allele is not an integer.\n");
        return R_NilValue;
    }

    workspace = calloc(1,sizeof(*workspace));
    if (!workspace) error("Cannot allocate population export workspace");
    st_work_adopt(workspace,st_population_export_destroy);
    i_bg = INTEGER(_bg)[0];

    pop_names_length = strlen(CHAR(STRING_ELT(_PopNames, 0)));

    if (pop_names_length > SIZE_MAX - 2)
    {
        info(-2, "Invalid population name in return_population\n");
        goto cleanup;
    }

    PopNames = (char *)malloc((pop_names_length + 2) * sizeof(*PopNames));

    if (PopNames == NULL)
    {
        info(-2, "Memory allocation failed in return_population\n");
        goto cleanup;
    }

    workspace->names = PopNames;
    strcpy(PopNames, CHAR(STRING_ELT(_PopNames, 0)));

    if ((cMP < 0) || ((size_t)cMP > SIZE_MAX / sizeof(*AllSet)))
    {
        info(-2, "Invalid number of markers in return_population\n");
        goto cleanup;
    }

    if (cMP > 0)
    {
        AllSet = calloc((size_t)cMP,sizeof(*AllSet));
        workspace->sets = AllSet;workspace->markers = cMP;

        if (AllSet == NULL)
        {
            info(-2, "Memory allocation failed in return_population\n");
            goto cleanup;
        }
    }

    for (marker = 0; marker < cMP; marker++)
    {
        AllSet[marker].member = NULL;
        AllSet[marker].betrag = 0;
    }

    PopName = strtok(PopNames, " ");

    while (PopName != NULL)
    {
        find_list_element(PopName, &PopListEl);

        if (PopListEl == NULL)
        {
            info(-2, "Population not defined: %s\n", PopName);
            goto cleanup;
        }

        if (PopListEl->Pop->NoInds < 1)
        {
            info(-2, "Empty population %s\n", PopName);
            goto cleanup;
        }

        if ((unsigned long int)PopListEl->Pop->NoInds >
            (unsigned long int)(INT_MAX - NoInds))
        {
            info(-2, "Too many individuals in return_population\n");
            goto cleanup;
        }

        NoInds += (int)PopListEl->Pop->NoInds;

        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
            TInd *Ind = get_individual(PopListEl->Pop, ind);

            if (Ind->genotype == NULL)
            {
                info(-2, "First genotype_population\n");
                goto cleanup;
            }

            for (marker = 0; marker < cMP; marker++)
            {
                for (hom = 0; hom < NoHoms; hom++)
                {
                    if (*bg != Ind->genotype[NoHoms * marker + hom])
                    {
                        if (!insertmember(AllSet + marker,
                                          Ind->genotype[NoHoms * marker + hom]))
                        {
                            info(-2, "Memory allocation failed in "
                                     "return_population\n");
                            goto cleanup;
                        }
                    }
                }
            }
        }

        PopName = strtok(NULL, " ");
    }

    for (marker = 0; marker < cMP; marker++)
    {
        if (AllSet[marker].betrag > 1)
        {
            qsort(AllSet[marker].member, AllSet[marker].betrag, sizeof(int),
                  intcmp);
        }

        if (AllSet[marker].betrag > INT_MAX - NoAll)
        {
            info(-2, "Too many alleles in return_population\n");
            goto cleanup;
        }

        NoAll += AllSet[marker].betrag;
    }

    if (!isNull(mab))
    {
        st_population_mab_export(_PopNames, mab, i_bg, NoAll, NoInds);
        *success = 1;
        goto cleanup;
    }

    PROTECT(matrix = allocVector(VECSXP,NoInds));nprotect++;
    PROTECT(xAxis = allocVector(STRSXP,NoInds));nprotect++;
    PROTECT(yAxis = allocVector(NoAll ? STRSXP : INTSXP,NoAll));nprotect++;
    PROTECT(cls = mkString("data.frame"));nprotect++;
    for (int j = 0; j < NoInds; j++)
    {
        SEXP column = allocVector(INTSXP,NoAll);SET_VECTOR_ELT(matrix,j,column);
        memset(INTEGER(column),0,(size_t)NoAll*sizeof(int));
    }
    strcpy(PopNames, CHAR(STRING_ELT(_PopNames, 0)));

    PopName = strtok(PopNames, " ");

    while (PopName != NULL)
    {
        find_list_element(PopName, &PopListEl);

        if (PopListEl == NULL)
        {
            info(-2, "Population not defined: %s\n", PopName);
            goto cleanup;
        }

        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
            TInd *Ind = get_individual(PopListEl->Pop, ind);

            snprintf(tmp, sizeof(tmp), "%s.%ld ", PopListEl->Pop->Name,
                     (long)(ind + 1));

            SET_STRING_ELT(xAxis, col, mkChar(tmp));
            int offset = 0;
            int *column = INTEGER(VECTOR_ELT(matrix,col++));

            for (marker = 0; marker < cMP; marker++)
            {
                if (Ind->genotype == NULL)
                {
                    info(-2, "First genotype_population\n");
                    goto cleanup;
                }

                for (hom = 0; hom < NoHoms; hom++)
                {
                    if (*bg != Ind->genotype[NoHoms * marker + hom])
                    {
                        column
                        [offset + getPos(Ind->genotype[marker * NoHoms + hom],
                                         AllSet[marker])] = 1;
                    }
                }

                offset += AllSet[marker].betrag;
            }
        }

        PopName = strtok(NULL, " ");
    }

    for (chrom = 0, marker = 0, row = 0; chrom < NoChroms; chrom++)
    {
        for (cmarker = 0; cmarker < PopDescr[chrom].NrMappoints;
             cmarker++, marker++)
        {
            for (allele = 0; allele < AllSet[marker].betrag; allele++)
            {
                snprintf(tmp, sizeof(tmp), "%s.%d ",
                         PopDescr[chrom].Mappoint[cmarker].name,
                         AllSet[marker].member[allele]);

                SET_STRING_ELT(yAxis, row++, mkChar(tmp));
            }
        }
    }

    setAttrib(matrix,R_NamesSymbol,xAxis);
    setAttrib(matrix,R_RowNamesSymbol,yAxis);
    setAttrib(matrix,R_ClassSymbol,cls);
    result = matrix;

cleanup:

    if (nprotect > 0)
        UNPROTECT(nprotect);

    return result;
}

SEXP return_population(SEXP names, SEXP background)
{
    return st_population_export(names, background, R_NilValue, NULL);
}


void evaluate_genome(char **PopName, int *all, double *ratio, double *length,
                     double *total, int *blocks)
{
    TPopListEl *PopListEl;
    POPULATION ind;
    int chrom, hom, loc;
    double gLength = 0;
    TInd *Ind;

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Can't find Population %s\n", *PopName);
        return;
    }

    for (chrom = 0; chrom < NoChroms; chrom++)
        gLength += PopDescr[chrom].ChromLen * NoHoms;

    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
        Ind = get_individual(PopListEl->Pop, ind);
        for (chrom = 0; chrom < NoChroms; chrom++)
        {
            for (hom = 0; hom < NoHoms; hom++)
            {
                for (loc = 0; loc < Ind->Chrom[chrom].Hom[hom].NoLociUsed;
                     loc++)
                {
                    if (*all == Ind->Chrom[chrom].Hom[hom].Loc[loc].All)
                    {
                        if (loc + 1 < Ind->Chrom[chrom].Hom[hom].NoLociUsed)
                        {
                            length[ind] +=
                                Ind->Chrom[chrom].Hom[hom].Loc[loc + 1].Pos -
                                Ind->Chrom[chrom].Hom[hom].Loc[loc].Pos;
                        }
                        else
                        {
                            length[ind] +=
                                PopDescr[chrom].ChromLen -
                                Ind->Chrom[chrom].Hom[hom].Loc[loc].Pos;
                        }
                        blocks[ind]++;
                    }
                }
            }
        }
        ratio[ind] = length[ind] / gLength;
        total[ind] = gLength;
    }
}

void evaluate_genome_region(char **PopName, int *all, int *chrom, double *b,
                            double *e, double *ratio, double *length,
                            double *total, int *blocks)
{
    TPopListEl *PopListEl;
    POPULATION ind;
    int hom, loc;
    double gLength = 0;
    double begin, end;
    TInd *Ind;

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Can't find Population %s\n", *PopName);
        return;
    }

    if ((*chrom < 0) || (*chrom >= NoChroms))
    {
        info(-2, "Chromosome out of range: %d\n", *chrom + 1);
        return;
    }

    gLength += (*e - *b) * NoHoms;

    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
        Ind = get_individual(PopListEl->Pop, ind);

        for (hom = 0; hom < NoHoms; hom++)
        {
            for (loc = 0; loc < Ind->Chrom[*chrom].Hom[hom].NoLociUsed; loc++)
            {
                if (*all == Ind->Chrom[*chrom].Hom[hom].Loc[loc].All)
                {
                    begin = Ind->Chrom[*chrom].Hom[hom].Loc[loc].Pos;

                    if (loc + 1 < Ind->Chrom[*chrom].Hom[hom].NoLociUsed)
                    {
                        end = Ind->Chrom[*chrom].Hom[hom].Loc[loc + 1].Pos;
                    }
                    else
                    {
                        end = PopDescr[*chrom].ChromLen;
                    }

                    if ((end >= *b) && (begin <= *e))
                    {
                        if (begin < *b)
                            begin = *b;
                        if (end > *e)
                            end = *e;
                        length[ind] += end - begin;
                        blocks[ind]++;
                    }
                }
            }
        }
        ratio[ind] = length[ind] / gLength;
        total[ind] = gLength;
    }
}

void reset_mdp(void)
{
    MDP = 0;
    info(2, "MDP counter reset\n");
}

void set_mdp(char **cMDP)
{
    MDP = atol(*cMDP);
}

void evaluate_mdp(char **progeny, char **parent1, char **parent2,
                  char **effectfile)
{
    TEffects *Effect = NULL;
    unsigned long int oldMDP = MDP;
    int i;

    if (0 == NrEffects)
    {
        info(-2, "No Effectfiles loaded\n");
        return;
    }

    if (NULL == effectfile)
    {
        info(-2, "No effect specified\n");
        return;
    }

    if ((NULL != effectfile) && (NULL == (Effect = find_Effect(*effectfile))))
    {
        info(-2, "Can't find the effect %s\n", *effectfile);
        return;
    }

    {
        TPopListEl *PopListEl;
        POPULATION ind;
        int *counter = NULL;
        int *all = NULL;
        int hom;
        TInd *Ind;
        int error = 0;

        if (Effect->MPs > 0)
        {
            counter = (int *)malloc((size_t)Effect->MPs * sizeof(*counter));
            all = (int *)malloc((size_t)Effect->MPs * sizeof(*all));

            if ((NULL == counter) || (NULL == all))
            {
                info(-2, "Memory allocation failed in evaluate_mdp\n");
                error = 1;
                goto cleanup;
            }
        }

        for (i = 0; i < Effect->MPs; i++)
        {
            counter[i] = 0;
            all[i] = -1;
        }

        find_list_element(*parent1, &PopListEl);

        if (NULL == PopListEl)
        {
            info(-2, "Population not defined: %s\n", *parent1);
            error = 1;
            goto cleanup;
        }

        if (!st_population_all_genotyped(PopListEl->Pop))
        {
            info(-2, "First genotype population %s\n", *parent1);
            error = 1;
            goto cleanup;
        }

        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
            Ind = get_individual(PopListEl->Pop, ind);
            for (i = 0; i < Effect->MPs; i++)
            {
                if (counter[i] < 2)
                {
                    for (hom = 0; hom < NoHoms; hom++)
                    {
                        if (all[i] !=
                            Ind->genotype[NoHoms * Effect->MPset[i]->counter +
                                          hom])
                        {
                            all[i] =
                                Ind->genotype[NoHoms *
                                                  Effect->MPset[i]->counter +
                                              hom];
                            counter[i]++;
                        }
                    }
                }
            }
        }

        for (i = 0; i < Effect->MPs; i++)
        {
            if (1 == counter[i])
            {
                counter[i] = 0;
                all[i] = -1;
            }
        }

        find_list_element(*parent2, &PopListEl);

        if (NULL == PopListEl)
        {
            info(-2, "Population not defined: %s\n", *parent2);
            error = 1;
            goto cleanup;
        }

        if (!st_population_all_genotyped(PopListEl->Pop))
        {
            info(-2, "First genotype population %s\n", *parent2);
            error = 1;
            goto cleanup;
        }

        for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
            Ind = get_individual(PopListEl->Pop, ind);
            for (i = 0; i < Effect->MPs; i++)
            {
                if (counter[i] < 2)
                {
                    for (hom = 0; hom < NoHoms; hom++)
                    {
                        if (all[i] !=
                            Ind->genotype[NoHoms * Effect->MPset[i]->counter +
                                          hom])
                        {
                            all[i] =
                                Ind->genotype[NoHoms *
                                                  Effect->MPset[i]->counter +
                                              hom];
                            counter[i]++;
                        }
                    }
                }
            }
        }

        find_list_element(*progeny, &PopListEl);

        if (NULL == PopListEl)
        {
            info(-2, "Population not defined: %s\n", *progeny);
            error = 1;
            goto cleanup;
        }

        for (i = 0; i < Effect->MPs; i++)
        {
            if (counter[i] == 2)
            {
                MDP += PopListEl->Pop->NoInds;
            }
        }

    cleanup:

        free(counter);
        free(all);

        if (error)
            return;
    }

    info(1, "%lu marker data points required in this step (total %lu)\n",
         MDP - oldMDP, MDP);
}

void get_mdp(char **cmdp)
{
    snprintf(*cmdp, strlen(*cmdp), "%lu", MDP);
}

typedef struct
{
    int oldPos;
    TMappoint *Mappoint;
} TGenPop;

TGenPop *PGenPop;

int genpopcmp(const void *i, const void *j)
{
    int chrom0, chrom1;
    CPOS pos0, pos1;
    chrom0 = ((TGenPop *)i)->Mappoint->chrom;
    chrom1 = ((TGenPop *)j)->Mappoint->chrom;
    pos0 = ((TGenPop *)i)->Mappoint->Pos;
    pos1 = ((TGenPop *)j)->Mappoint->Pos;

    if (chrom0 > chrom1)
        return 1;

    if (chrom0 < chrom1)
        return -1;

    if (pos0 > pos1)
        return 1;

    if (pos0 < pos1)
        return -1;

    return 0;
}

void generate_population(int *no_pop, int *no_ind, int *no_mar, int *no_all,
                         int *no_pop_u, int *pop_u, int *no_ind_u, int *ind_u,
                         int *no_mar_u, int *mar_u, double *data, double *freq,
                         char **pop_list, char **mar_list, int *al_list,
                         int *bg, int *bc)
{
    int i, j, k, l;
    double *ptr_data, *ptr_dataC, *dataC = NULL;
    double count;
    int rowpos_d;
    size_t colpos, dataC_count;
    int nrow_d, ncol_d;
    int *coffset = NULL;
    int *roff_d = NULL;
    int *roff_f = NULL;
    int *_ind = NULL, *_chrom = NULL, *_hom = NULL, *_all = NULL;
    int *candidates = NULL;
    int *tmp_int;
    double *_pos = NULL, *tmp_double;

    (void)*freq;

    PGenPop = NULL;

    if ((*no_pop < 1) || (*no_mar < 1) || (*no_pop_u < 1) || (*no_mar_u < 1) ||
        (NoHoms < 1))
        goto invalid_or_allocation_fail;

    coffset = (int *)malloc((size_t)*no_pop * sizeof(*coffset));
    roff_d = (int *)malloc((size_t)*no_mar * sizeof(*roff_d));
    roff_f = (int *)malloc((size_t)*no_mar_u * sizeof(*roff_f));

    if ((NULL == coffset) || (NULL == roff_d) || (NULL == roff_f))
        goto invalid_or_allocation_fail;

    coffset[0] = 0;
    for (k = 1; k < *no_pop; k++)
    {
        coffset[k] = coffset[k - 1] + no_ind[k - 1];
    }

    roff_d[0] = 0;
    for (i = 1; i < *no_mar; i++)
    {
        roff_d[i] = roff_d[i - 1] + no_all[i - 1];
    }

    roff_f[0] = 0;
    for (i = 1; i < *no_mar_u; i++)
    {
        roff_f[i] = roff_f[i - 1] + no_all[mar_u[i - 1] - 1];
    }

    nrow_d = roff_d[*no_mar - 1] + no_all[*no_mar - 1];
    ncol_d = coffset[*no_pop - 1] + no_ind[*no_pop - 1];

    {
        int sum_ind_u = 0;

        if ((nrow_d < 1) || (ncol_d < 1) ||
            ((size_t)nrow_d > SIZE_MAX / (size_t)ncol_d))
            goto invalid_or_allocation_fail;

        dataC_count = (size_t)nrow_d * (size_t)ncol_d;
        if (dataC_count > SIZE_MAX / sizeof(*dataC))
            goto invalid_or_allocation_fail;

        dataC = (double *)malloc(dataC_count * sizeof(*dataC));
        if (NULL == dataC)
            goto invalid_or_allocation_fail;

        for (k = 0; k < *no_pop_u; k++)
        {
            for (l = 0; l < no_ind_u[k]; l++)
            {
                colpos = ((size_t)coffset[pop_u[k] - 1] +
                          (size_t)(ind_u[l + sum_ind_u] - 1)) *
                         (size_t)nrow_d;
                for (i = 0; i < *no_mar_u; i++)
                {
                    count = 0;
                    for (j = 0; j < no_all[mar_u[i] - 1]; j++)
                    {
                        rowpos_d = (roff_d[mar_u[i] - 1] + j);
                        ptr_data = &data[colpos + rowpos_d];
                        if (!ISNA(*ptr_data) && (*ptr_data > 0))
                        {
                            *ptr_data = 1;
                            count++;
                        }
                    }
                    for (j = 0; j < no_all[mar_u[i] - 1]; j++)
                    {
                        rowpos_d = (roff_d[mar_u[i] - 1] + j);
                        ptr_data = &data[colpos + rowpos_d];
                        ptr_dataC = &dataC[colpos + rowpos_d];
                        if (count > 0)
                            (*ptr_dataC) = (*ptr_data) / (double)count;
                        else
                            (*ptr_dataC) = (*ptr_data);
                    }
                }
            }
            sum_ind_u += no_ind_u[k];
        }

        {
            int chrom, hom, tokencounter = 0;
            char *token;

            if (0 == cMP)
            {
                info(-2, "No linkage map defined\n");
                goto fail;
            }

            PGenPop = (TGenPop *)malloc((size_t)*no_mar * sizeof(*PGenPop));
            if (NULL == PGenPop)
                goto invalid_or_allocation_fail;
            token = strtok(*mar_list, " ");

            while (NULL != token)
            {
                PGenPop[tokencounter].Mappoint = NULL;
                for (chrom = 0; chrom < NoChroms; chrom++)
                {
                    for (i = 0; i < PopDescr[chrom].NrMappoints; i++)
                    {
                        if (0 ==
                            strcmp(PopDescr[chrom].Mappoint[i].name, token))
                        {
                            PGenPop[tokencounter].oldPos = tokencounter;
                            PGenPop[tokencounter].Mappoint =
                                &(PopDescr[chrom].Mappoint[i]);
                        }
                    }
                }

                if (NULL == PGenPop[tokencounter].Mappoint)
                {
                    info(-2, "Can't find marker %s!", token);
                    goto fail;
                }

                tokencounter++;

                token = strtok(NULL, " ");
            }

            qsort(PGenPop, *no_mar, sizeof(TGenPop), genpopcmp);

            {
                size_t offset = 0;
                token = strtok(*pop_list, " ");

                for (k = 0; k < *no_pop_u; k++)
                {
                    int nrow = 0;
                    int ind;

                    _ind = (int *)malloc(sizeof(*_ind));
                    _chrom = (int *)malloc(sizeof(*_chrom));
                    _hom = (int *)malloc(sizeof(*_hom));
                    _pos = (double *)malloc(sizeof(*_pos));
                    _all = (int *)malloc(sizeof(*_all));

                    if ((NULL == _ind) || (NULL == _chrom) || (NULL == _hom) ||
                        (NULL == _pos) || (NULL == _all))
                        goto invalid_or_allocation_fail;

                    for (ind = 0; ind < no_ind_u[k]; ind++)
                    {
                        for (chrom = 0; chrom < NoChroms; chrom++)
                        {
                            int tmpcounter = 0;
                            CPOS lastpos = 0;
                            for (i = 0; i < *no_mar; i++)
                            {
                                int nr_cand = NoHoms, c, act_cand = 0;

                                candidates = (int *)malloc((size_t)NoHoms *
                                                           sizeof(*candidates));
                                if (NULL == candidates)
                                    goto invalid_or_allocation_fail;

                                if (PGenPop[i].Mappoint->chrom - 1 == chrom)
                                {
                                    for (c = 0; c < NoHoms; c++)
                                    {
                                        candidates[c] = *bg;
                                    }
                                    for (j = 0; j < no_all[PGenPop[i].oldPos];
                                         j++)
                                    {
                                        while (
                                            dataC[offset +
                                                  (size_t)ind * (size_t)nrow_d +
                                                  roff_d[PGenPop[i].oldPos] +
                                                  j] > 0)
                                        {
                                            if (nr_cand == act_cand)
                                            {
                                                ++nr_cand;
                                                tmp_int = (int *)realloc(
                                                    candidates,
                                                    (size_t)nr_cand *
                                                        sizeof(*candidates));
                                                if (NULL == tmp_int)
                                                    goto invalid_or_allocation_fail;
                                                candidates = tmp_int;
                                            }
                                            candidates[act_cand++] =
                                                (int)al_list
                                                    [roff_d[PGenPop[i].oldPos] +
                                                     j];
                                            dataC[offset +
                                                  (size_t)ind * (size_t)nrow_d +
                                                  roff_d[PGenPop[i].oldPos] +
                                                  j] -= (float)1 / NoHoms;
                                        }
                                    }

                                    qsort(candidates, nr_cand, sizeof(int),
                                          intcmp);

                                    for (hom = 0; hom < NoHoms; hom++)
                                    {
                                        nrow++;

                                        tmp_int = (int *)realloc(
                                            _ind, (size_t)nrow * sizeof(*_ind));
                                        if (NULL == tmp_int)
                                            goto invalid_or_allocation_fail;
                                        _ind = tmp_int;

                                        tmp_int = (int *)realloc(
                                            _chrom,
                                            (size_t)nrow * sizeof(*_chrom));
                                        if (NULL == tmp_int)
                                            goto invalid_or_allocation_fail;
                                        _chrom = tmp_int;

                                        tmp_int = (int *)realloc(
                                            _hom, (size_t)nrow * sizeof(*_hom));
                                        if (NULL == tmp_int)
                                            goto invalid_or_allocation_fail;
                                        _hom = tmp_int;

                                        tmp_double = (double *)realloc(
                                            _pos, (size_t)nrow * sizeof(*_pos));
                                        if (NULL == tmp_double)
                                            goto invalid_or_allocation_fail;
                                        _pos = tmp_double;

                                        tmp_int = (int *)realloc(
                                            _all, (size_t)nrow * sizeof(*_all));
                                        if (NULL == tmp_int)
                                            goto invalid_or_allocation_fail;
                                        _all = tmp_int;

                                        *(_ind + nrow - 1) = ind + 1;
                                        *(_chrom + nrow - 1) = chrom + 1;
                                        *(_hom + nrow - 1) = hom + 1;

                                        if (0 == tmpcounter)
                                        {
                                            *(_pos + nrow - 1) = 0;
                                        }
                                        else
                                        {
                                            *(_pos + nrow - 1) =
                                                (PGenPop[i].Mappoint->Pos +
                                                 lastpos) /
                                                2;
                                        }

                                        {
                                            int sel_all;

                                            if (1 == *bc)
                                            {
                                                sel_all = hom;
                                            }
                                            else
                                            {
                                                sel_all =
                                                    (int)gsl_rng_uniform_int(
                                                        rng[0],
                                                        (unsigned long int)
                                                            nr_cand);
                                            }
                                            nr_cand--;
                                            *(_all + nrow - 1) =
                                                candidates[sel_all];
                                            candidates[sel_all] =
                                                candidates[nr_cand];
                                        }
                                    }

                                    tmpcounter++;
                                    lastpos = PGenPop[i].Mappoint->Pos;

                                    if (nr_cand > 0)
                                    {
                                        info(-1,
                                             "More than two alleles. Marker %s "
                                             "Pop: %s Ind: %d\n",
                                             PGenPop[i].Mappoint->name, token,
                                             ind + 1);
                                    }
                                }

                                free(candidates);
                                candidates = NULL;
                            }

                            if (0 == tmpcounter)
                            {
                                if (0 == ind)
                                {
                                    info(-1,
                                         "No marker defined for chromosome "
                                         "%d.\n",
                                         chrom + 1);
                                }

                                for (hom = 0; hom < NoHoms; hom++)
                                {
                                    nrow++;

                                    tmp_int = (int *)realloc(
                                        _ind, (size_t)nrow * sizeof(*_ind));
                                    if (NULL == tmp_int)
                                        goto invalid_or_allocation_fail;
                                    _ind = tmp_int;

                                    tmp_int = (int *)realloc(
                                        _chrom, (size_t)nrow * sizeof(*_chrom));
                                    if (NULL == tmp_int)
                                        goto invalid_or_allocation_fail;
                                    _chrom = tmp_int;

                                    tmp_int = (int *)realloc(
                                        _hom, (size_t)nrow * sizeof(*_hom));
                                    if (NULL == tmp_int)
                                        goto invalid_or_allocation_fail;
                                    _hom = tmp_int;

                                    tmp_double = (double *)realloc(
                                        _pos, (size_t)nrow * sizeof(*_pos));
                                    if (NULL == tmp_double)
                                        goto invalid_or_allocation_fail;
                                    _pos = tmp_double;

                                    tmp_int = (int *)realloc(
                                        _all, (size_t)nrow * sizeof(*_all));
                                    if (NULL == tmp_int)
                                        goto invalid_or_allocation_fail;
                                    _all = tmp_int;

                                    *(_ind + nrow - 1) = ind + 1;
                                    *(_chrom + nrow - 1) = chrom + 1;
                                    *(_hom + nrow - 1) = hom + 1;
                                    *(_pos + nrow - 1) = 0;
                                    *(_all + nrow - 1) = *bg;
                                }
                            }
                        }
                    }

                    offset += (size_t)ind * (size_t)nrow_d;

                    {
                        char cnrow[50];
                        char *prow = cnrow;

                        snprintf(cnrow, sizeof(cnrow), "%d", nrow);

                        {
                            TPopListEl *generated = NULL;

                            find_list_element(token, &generated);
                            if (NULL != generated)
                            {
                                info(-2, "First remove population %s\n", token);
                                goto fail;
                            }

                            init_population(&token, &prow, &(_ind[0]),
                                            &(_chrom[0]), &(_hom[0]),
                                            &(_pos[0]), &(_all[0]));

                            find_list_element(token, &generated);
                            if (NULL == generated)
                                goto fail;
                        }
                    }

                    token = strtok(NULL, " ");
                    nrow = 0;
                    free(_ind);
                    free(_chrom);
                    free(_hom);
                    free(_pos);
                    free(_all);
                    _ind = NULL;
                    _chrom = NULL;
                    _hom = NULL;
                    _pos = NULL;
                    _all = NULL;
                }
            }

            free(PGenPop);
            PGenPop = NULL;
        }

        free(dataC);
    }

    free(coffset);
    free(roff_d);
    free(roff_f);
    return;

invalid_or_allocation_fail:

    info(-2, "Invalid dimensions or memory allocation failed in "
             "generate_population\n");

fail:

    free(candidates);
    free(_ind);
    free(_chrom);
    free(_hom);
    free(_pos);
    free(_all);
    free(PGenPop);
    PGenPop = NULL;
    free(dataC);
    free(coffset);
    free(roff_d);
    free(roff_f);
    return;
}

SEXP splitdt(SEXP in_string)
{
    const char *ma;
    const char *a;
    int i;
    int n = length(in_string);
    SEXP element;
    SEXP out_string;

    PROTECT(out_string = allocMatrix(STRSXP, n, 2));
    for (i = 0; i < n; i++)
    {
        element = STRING_ELT(in_string, i);

        if (NA_STRING == element)
        {
            SET_STRING_ELT(out_string, i, NA_STRING);
            SET_STRING_ELT(out_string, n + i, NA_STRING);
            continue;
        }

        ma = CHAR(element);

        a = strchr(ma, '.');
        if (NULL == a)
            a = strchr(ma, '_');

        if (NULL == a)
        {
            SET_STRING_ELT(out_string, i, element);
            SET_STRING_ELT(out_string, n + i, mkChar(""));
        }
        else
        {
            SET_STRING_ELT(out_string, i, mkCharLen(ma, (int)(a - ma)));
            SET_STRING_ELT(out_string, n + i, mkChar(a + 1));
        }
    }
    UNPROTECT(1);
    return (out_string);
}

void pq_genotype(int *nomar, int *noind, int *mgt, int *k)
{
    int i, j;
    int *p_m = mgt;
    int *p_k = k;

    for (j = 0; j < *nomar; j++)
        for (i = 0; i < *noind; i++)
        {
            if ((*p_m == 0) && (*(1 + p_m) == 1))
                *p_k = 0;
            else if ((*p_m == 1) && (*(1 + p_m) == 1))
                *p_k = 1;
            else if ((*p_m == 1) && (*(1 + p_m) == 0))
                *p_k = 2;

            p_k++;

            p_m++;

            p_m++;
        }
}

void evaluate_ld(char **PopName, int *all, int *chr, double *pos, double *ld)

{
    TPopListEl *PopListEl;
    POPULATION ind;
    int chrom, hom, loc;
    double gLength = 0;
    TInd *Ind;
    double upper, lower;

    (void)gLength;

    if ((*chr < 1) || (*chr > NoChroms))
    {
        info(-2, "Chromosome out of range: %d\n", *chr);
        return;
    }

    find_list_element(*PopName, &PopListEl);

    if (NULL == PopListEl)
    {
        info(-2, "Can't find Population %s\n", *PopName);
        return;
    }

    for (chrom = 0; chrom < NoChroms; chrom++)
        gLength += PopDescr[chrom].ChromLen * NoHoms;

    for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
        Ind = get_individual(PopListEl->Pop, ind);

        chrom = *chr - 1;
        {
            for (hom = 0; hom < NoHoms; hom++)
            {
                for (loc = 0; loc < Ind->Chrom[chrom].Hom[hom].NoLociUsed;
                     loc++)
                {
                    if (*all == Ind->Chrom[chrom].Hom[hom].Loc[loc].All)
                    {
                        if (loc + 1 < Ind->Chrom[chrom].Hom[hom].NoLociUsed)
                        {
                            upper = Ind->Chrom[chrom].Hom[hom].Loc[loc + 1].Pos;
                        }
                        else
                        {
                            upper = PopDescr[chrom].ChromLen;
                        }
                        lower = Ind->Chrom[chrom].Hom[hom].Loc[loc].Pos;
                        if ((lower <= *pos) && (*pos < upper))
                        {
                            if (ld[ind] == 0)
                            {
                                ld[ind] = upper - lower;
                            }
                            else
                            {
                                ld[ind] = ld[ind] / 2 + (upper - lower) / 2;
                            }
                        }
                    }
                }
            }
        }
    }
}

#include <R_ext/Random.h>

static int AllowZeroFrequencies = 0;

void set_allow_zeros(int *allow)
{
    AllowZeroFrequencies = *allow;
}

void get_allow_zeros(int *allow)
{
    *allow = AllowZeroFrequencies;
}

SEXP splitdot(SEXP in_string)
{

    const char *ma, *a;
    int i;
    int n = length(in_string);
    SEXP element;
    SEXP out_string;

    PROTECT(out_string = allocMatrix(STRSXP, n, 2));
    for (i = 0; i < n; i++)
    {
        element = STRING_ELT(in_string, i);

        if (element == NA_STRING)
        {
            SET_STRING_ELT(out_string, i, NA_STRING);
            SET_STRING_ELT(out_string, n + i, NA_STRING);
            continue;
        }

        ma = CHAR(element);

        a = strchr(ma, '.');
        if (NULL == a)
            a = strchr(ma, '_');

        if (NULL == a)
        {
            SET_STRING_ELT(out_string, i, element);
            SET_STRING_ELT(out_string, n + i, mkChar(""));
        }
        else
        {
            SET_STRING_ELT(out_string, i, mkCharLen(ma, (int)(a - ma)));
            SET_STRING_ELT(out_string, n + i, mkChar(a + 1));
        }
    }
    UNPROTECT(1);
    return (out_string);
}

void allele_freq(int *no_pop, int *no_ind, int *no_mar, int *no_all,
                 int *no_ind_u, int *ind_u, double *data, double *freq)
{
    size_t *coff = NULL;
    size_t *roff = NULL;
    int i, j, k, l;
    size_t rowpos, colpos, nrow, ncol_d, data_count;
    double *ptr_data, *ptr_freq, *ptr_dataC;
    double *dataC = NULL;
    double count;

    if ((*no_pop < 1) || (*no_mar < 1))
        goto fail;

    if ((size_t)*no_pop > SIZE_MAX / sizeof(*coff))
        goto fail;
    if ((size_t)*no_mar > SIZE_MAX / sizeof(*roff))
        goto fail;

    coff = (size_t *)malloc((size_t)(*no_pop) * sizeof(*coff));
    if (NULL == coff)
        goto fail;

    roff = (size_t *)malloc((size_t)(*no_mar) * sizeof(*roff));
    if (NULL == roff)
        goto fail;

    roff[0] = coff[0] = 0;

    for (k = 1; k < *no_pop; k++)
    {
        if ((no_ind[k - 1] < 0) ||
            (coff[k - 1] > SIZE_MAX - (size_t)no_ind[k - 1]))
            goto fail;
        coff[k] = coff[k - 1] + (size_t)no_ind[k - 1];
    }

    for (i = 1; i < *no_mar; i++)
    {
        if ((no_all[i - 1] < 0) ||
            (roff[i - 1] > SIZE_MAX - (size_t)no_all[i - 1]))
            goto fail;
        roff[i] = roff[i - 1] + (size_t)no_all[i - 1];
    }

    if ((no_all[*no_mar - 1] < 0) ||
        (roff[*no_mar - 1] > SIZE_MAX - (size_t)no_all[*no_mar - 1]))
        goto fail;

    nrow = roff[*no_mar - 1] + (size_t)no_all[*no_mar - 1];

    if ((no_ind[*no_pop - 1] < 0) ||
        (coff[*no_pop - 1] > SIZE_MAX - (size_t)no_ind[*no_pop - 1]))
        goto fail;

    ncol_d = coff[*no_pop - 1] + (size_t)no_ind[*no_pop - 1];

    if ((nrow != 0) && (ncol_d > SIZE_MAX / nrow))
        goto fail;
    data_count = nrow * ncol_d;

    if (data_count > SIZE_MAX / sizeof(*dataC))
        goto fail;

    if ((nrow != 0) && ((size_t)*no_pop > SIZE_MAX / nrow))
        goto fail;

    if (data_count > 0)
    {
        dataC = (double *)malloc(data_count * sizeof(*dataC));
        if (NULL == dataC)
            goto fail;
    }

    for (k = 0; k < *no_pop; k++)
    {
        for (l = 0; l < no_ind_u[k]; l++)
        {
            colpos =
                (coff[k] + (size_t)(ind_u[coff[k] + (size_t)l] - 1)) * nrow;
            for (i = 0; i < *no_mar; i++)
            {
                count = 0;
                for (j = 0; j < no_all[i]; j++)
                {
                    rowpos = roff[i] + (size_t)j;
                    ptr_data = &data[colpos + rowpos];
                    if (*ptr_data > 0)
                    {
                        *ptr_data = 1;
                        count += 1;
                    }
                }
                for (j = 0; j < no_all[i]; j++)
                {
                    rowpos = roff[i] + (size_t)j;
                    ptr_data = &data[colpos + rowpos];
                    ptr_dataC = &dataC[colpos + rowpos];

                    if (AllowZeroFrequencies)
                    {
                        *ptr_dataC = *ptr_data;
                    }
                    else
                    {
                        if (count > 0)
                        {
                            if (*ptr_data > 0)
                                (*ptr_dataC) = (*ptr_data) / (double)count;
                            else
                                (*ptr_dataC) = 0;
                        }
                        else
                        {
                            *ptr_dataC = -1;
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < *no_mar; i++)
        for (j = 0; j < no_all[i]; j++)
        {
            rowpos = roff[i] + (size_t)j;
            for (k = 0; k < *no_pop; k++)
            {
                ptr_freq = &freq[(size_t)k * nrow + rowpos];
                *ptr_freq = -1;
                count = 0;
                for (l = 0; l < no_ind_u[k]; l++)
                {
                    colpos =
                        (coff[k] + (size_t)(ind_u[coff[k] + (size_t)l] - 1)) *
                        nrow;
                    ptr_dataC = &dataC[colpos + rowpos];
                    if (*ptr_dataC != -1)
                    {
                        count += 1;
                        if (*ptr_freq == -1)
                            (*ptr_freq) = 0;
                        (*ptr_freq) += (*ptr_dataC);
                    }
                }
                if (count > 0)
                {
                    if (*ptr_freq > 0)
                        *ptr_freq = *ptr_freq / (double)count;
                    else
                        *ptr_freq = 0;
                }
                else
                    *ptr_freq = -1;
            }
        }

    free(dataC);
    free(roff);
    free(coff);
    return;

fail:

    free(dataC);
    free(roff);
    free(coff);

    info(-2, "Invalid dimensions or memory allocation failed in allele_freq\n");
    return;
}

void correct_missing(int *no_pop, int *no_ind, int *no_mar, int *no_all,
                     int *no_ind_u, int *ind_u, double *data, double *freq)
{

    (void)freq;

    size_t *coff = NULL;
    size_t *roff = NULL;
    int i, j, k, l;
    size_t rowpos, colpos, nrow, ncol_d, data_count;
    double *ptr_data, *ptr_dataC;
    double *dataC = NULL;

    double count;

    if ((*no_pop < 1) || (*no_mar < 1))
        goto fail;

    if ((size_t)*no_pop > SIZE_MAX / sizeof(*coff))
        goto fail;
    if ((size_t)*no_mar > SIZE_MAX / sizeof(*roff))
        goto fail;

    coff = (size_t *)malloc((size_t)(*no_pop) * sizeof(*coff));
    if (NULL == coff)
        goto fail;

    roff = (size_t *)malloc((size_t)(*no_mar) * sizeof(*roff));
    if (NULL == roff)
        goto fail;

    roff[0] = coff[0] = 0;

    for (k = 1; k < *no_pop; k++)
    {
        if ((no_ind[k - 1] < 0) ||
            (coff[k - 1] > SIZE_MAX - (size_t)no_ind[k - 1]))
            goto fail;
        coff[k] = coff[k - 1] + (size_t)no_ind[k - 1];
    }

    for (i = 1; i < *no_mar; i++)
    {
        if ((no_all[i - 1] < 0) ||
            (roff[i - 1] > SIZE_MAX - (size_t)no_all[i - 1]))
            goto fail;
        roff[i] = roff[i - 1] + (size_t)no_all[i - 1];
    }

    if ((no_all[*no_mar - 1] < 0) ||
        (roff[*no_mar - 1] > SIZE_MAX - (size_t)no_all[*no_mar - 1]))
        goto fail;

    nrow = roff[*no_mar - 1] + (size_t)no_all[*no_mar - 1];

    if ((no_ind[*no_pop - 1] < 0) ||
        (coff[*no_pop - 1] > SIZE_MAX - (size_t)no_ind[*no_pop - 1]))
        goto fail;

    ncol_d = coff[*no_pop - 1] + (size_t)no_ind[*no_pop - 1];

    if ((nrow != 0) && (ncol_d > SIZE_MAX / nrow))
        goto fail;
    data_count = nrow * ncol_d;

    if (data_count > SIZE_MAX / sizeof(*dataC))
        goto fail;

    if (data_count > 0)
    {
        dataC = (double *)malloc(data_count * sizeof(*dataC));
        if (NULL == dataC)
            goto fail;
    }

    for (k = 0; k < *no_pop; k++)
    {
        for (l = 0; l < no_ind_u[k]; l++)
        {
            colpos = (coff[k] + (size_t)(ind_u[coff[k] + l] - 1)) * nrow;
            for (i = 0; i < *no_mar; i++)
            {
                count = 0;
                for (j = 0; j < no_all[i]; j++)
                {
                    rowpos = roff[i] + (size_t)j;
                    ptr_data = &data[colpos + rowpos];
                    if (*ptr_data > 0)
                    {
                        *ptr_data = 1;
                        count += 1;
                    }
                }
                for (j = 0; j < no_all[i]; j++)
                {
                    rowpos = roff[i] + (size_t)j;
                    ptr_data = &data[colpos + rowpos];
                    ptr_dataC = &dataC[colpos + rowpos];

                    if (AllowZeroFrequencies)
                    {
                        *ptr_dataC = *ptr_data;
                    }
                    else
                    {
                        if (count > 0)
                        {
                            if (*ptr_data > 0)
                                (*ptr_dataC) = (*ptr_data) / (double)count;
                            else
                                (*ptr_dataC) = 0;
                        }
                        else
                        {
                            *ptr_dataC = -1;
                        }
                    }
                }
            }
        }
    }

    for (k = 0; k < *no_pop; k++)
    {
        for (l = 0; l < no_ind_u[k]; l++)
        {
            colpos = (coff[k] + (size_t)(ind_u[coff[k] + l] - 1)) * nrow;
            for (i = 0; i < *no_mar; i++)
            {
                for (j = 0; j < no_all[i]; j++)
                {
                    rowpos = roff[i] + (size_t)j;
                    ptr_data = &data[colpos + rowpos];
                    ptr_dataC = &dataC[colpos + rowpos];
                    (*ptr_data) = (*ptr_dataC);
                }
            }
        }
    }

    free(dataC);
    free(roff);
    free(coff);
    return;

fail:

    free(dataC);
    free(roff);
    free(coff);

    info(-2,
         "Invalid dimensions or memory allocation failed in correct_missing\n");
    return;
}

void rd(int *no_pop, int *no_mar, int *no_all, int *no_mar_u, int *mar_u,
        double *freq, double *dist)
{

    long i, j, k1, k2, m, observed, nrow, rowpos;
    double t1, t2;
    double *ptr_dist, *ptr_freq1, *ptr_freq2;
    int *roff;

    if (*no_mar < 1)
    {
        info(-2, "Invalid number of markers in rd\n");
        return;
    }

    if ((size_t)*no_mar > SIZE_MAX / sizeof(*roff))
    {
        info(-2, "Memory allocation failed in rd\n");
        return;
    }

    roff = (int *)malloc((size_t)*no_mar * sizeof(*roff));

    if (NULL == roff)
    {
        info(-2, "Memory allocation failed in rd\n");
        return;
    }

    roff[0] = 0;
    for (i = 1; i < *no_mar; i++)
        roff[i] = roff[i - 1] + no_all[i - 1];
    nrow = roff[*no_mar - 1] + no_all[*no_mar - 1];

    ptr_dist = dist;

    for (k1 = 0; k1 < (*no_pop - 1); k1++)
        for (k2 = k1 + 1; k2 < (*no_pop); k2++)
        {
            t2 = 0;
            m = 0;
            for (i = 0; i < *no_mar_u; i++)
            {
                observed = 1;
                t1 = 0;
                for (j = 0; j < no_all[mar_u[i] - 1]; j++)
                {
                    rowpos = (roff[mar_u[i] - 1] + j);
                    ptr_freq1 = &freq[k1 * nrow + rowpos];
                    ptr_freq2 = &freq[k2 * nrow + rowpos];
                    if ((*ptr_freq1 != -1) && (*ptr_freq2 != -1))
                        t1 += pow((*ptr_freq1 - *ptr_freq2), 2);
                    else
                        observed = 0;
                }
                if (observed)
                {
                    t2 += sqrt(t1 / 2);
                    m += 1;
                }
            }
            (*ptr_dist) = t2 / (double)m;
            ptr_dist++;
        }

    free(roff);
}

void mrd(int *no_pop, int *no_mar, int *no_all, int *no_mar_u, int *mar_u,
         double *freq, double *dist)
{
    long i, j, k1, k2, m, observed, nrow, rowpos;
    double t1;
    double *ptr_dist, *ptr_freq1, *ptr_freq2;
    int *roff;

    if (*no_mar < 1)
    {
        info(-2, "Invalid number of markers in mrd\n");
        return;
    }

    if ((size_t)*no_mar > SIZE_MAX / sizeof(*roff))
    {
        info(-2, "Memory allocation failed in mrd\n");
        return;
    }

    roff = (int *)malloc((size_t)*no_mar * sizeof(*roff));

    if (NULL == roff)
    {
        info(-2, "Memory allocation failed in mrd\n");
        return;
    }

    roff[0] = 0;
    for (i = 1; i < *no_mar; i++)
        roff[i] = roff[i - 1] + no_all[i - 1];
    nrow = roff[*no_mar - 1] + no_all[*no_mar - 1];

    ptr_dist = dist;

    for (k1 = 0; k1 < (*no_pop - 1); k1++)
        for (k2 = k1 + 1; k2 < (*no_pop); k2++)
        {
            t1 = 0;
            m = 0;
            for (i = 0; i < *no_mar_u; i++)
            {
                observed = 1;
                for (j = 0; j < no_all[mar_u[i] - 1]; j++)
                {
                    rowpos = (roff[mar_u[i] - 1] + j);
                    ptr_freq1 = &freq[k1 * nrow + rowpos];
                    ptr_freq2 = &freq[k2 * nrow + rowpos];
                    if ((*ptr_freq1 != -1) && (*ptr_freq2 != -1))
                        t1 += pow((*ptr_freq1 - *ptr_freq2), 2);
                    else
                        observed = 0;
                }
                if (observed)
                    m += 1;
            }
            (*ptr_dist) = sqrt(t1 / 2 / (double)m);
            ptr_dist++;
        }

    free(roff);
}

void euc(int *no_pop, int *no_mar, int *no_all, int *no_mar_u, int *mar_u,
         double *freq, double *dist)
{
    long i, j, k1, k2, m, observed, nrow, rowpos;
    double t1;
    double *ptr_dist, *ptr_freq1, *ptr_freq2;
    int *roff;

    if (*no_mar < 1)
    {
        info(-2, "Invalid number of markers in euc\n");
        return;
    }

    if ((size_t)*no_mar > SIZE_MAX / sizeof(*roff))
    {
        info(-2, "Memory allocation failed in euc\n");
        return;
    }

    roff = (int *)malloc((size_t)*no_mar * sizeof(*roff));

    if (NULL == roff)
    {
        info(-2, "Memory allocation failed in euc\n");
        return;
    }

    (void)m;

    roff[0] = 0;
    for (i = 1; i < *no_mar; i++)
        roff[i] = roff[i - 1] + no_all[i - 1];
    nrow = roff[*no_mar - 1] + no_all[*no_mar - 1];

    ptr_dist = dist;

    for (k1 = 0; k1 < (*no_pop - 1); k1++)
        for (k2 = k1 + 1; k2 < (*no_pop); k2++)
        {
            t1 = 0;
            m = 0;
            for (i = 0; i < *no_mar_u; i++)
            {
                observed = 1;
                for (j = 0; j < no_all[mar_u[i] - 1]; j++)
                {
                    rowpos = (roff[mar_u[i] - 1] + j);
                    ptr_freq1 = &freq[k1 * nrow + rowpos];
                    ptr_freq2 = &freq[k2 * nrow + rowpos];
                    if ((*ptr_freq1 != -1) && (*ptr_freq2 != -1))
                        t1 += pow((*ptr_freq1 - *ptr_freq2), 2);
                    else
                        observed = 0;
                }
                if (observed)
                    m += 1;
            }
            (*ptr_dist) = sqrt(t1);
            ptr_dist++;
        }

    free(roff);
}

void jac(int *no_pop, int *no_mar, int *no_all, int *no_mar_u, int *mar_u,
         double *freq, double *dist)
{
    long i, j, k1, k2, nrow, rowpos;
    double zaehler, nenner, fa, fb;
    double *ptr_dist, *ptr_freq1, *ptr_freq2;
    int *roff;

    if (*no_mar < 1)
    {
        info(-2, "Invalid number of markers in jac\n");
        return;
    }

    if ((size_t)*no_mar > SIZE_MAX / sizeof(*roff))
    {
        info(-2, "Memory allocation failed in jac\n");
        return;
    }

    roff = (int *)malloc((size_t)*no_mar * sizeof(*roff));

    if (NULL == roff)
    {
        info(-2, "Memory allocation failed in jac\n");
        return;
    }

    roff[0] = 0;
    for (i = 1; i < *no_mar; i++)
        roff[i] = roff[i - 1] + no_all[i - 1];
    nrow = roff[*no_mar - 1] + no_all[*no_mar - 1];

    ptr_dist = dist;

    for (k1 = 0; k1 < (*no_pop - 1); k1++)
        for (k2 = k1 + 1; k2 < (*no_pop); k2++)
        {
            zaehler = nenner = 0;
            for (i = 0; i < *no_mar_u; i++)
            {
                for (j = 0; j < no_all[mar_u[i] - 1]; j++)
                {
                    rowpos = (roff[mar_u[i] - 1] + j);
                    ptr_freq1 = &freq[k1 * nrow + rowpos];
                    ptr_freq2 = &freq[k2 * nrow + rowpos];
                    if ((*ptr_freq1 != -1) && (*ptr_freq2 != -1))
                    {
                        fa = (double)(*ptr_freq1 > 0 ? 1 : 0);
                        fb = (double)(*ptr_freq2 > 0 ? 1 : 0);
                        zaehler += fa * fb;
                        nenner += fa + fb - fa * fb;
                    }
                }
            }
            if (nenner > 0)
                (*ptr_dist) = zaehler / nenner;
            else
                (*ptr_dist) = 0;
            ptr_dist++;
        }

    free(roff);
}

void dic(int *no_pop, int *no_mar, int *no_all, int *no_mar_u, int *mar_u,
         double *freq, double *dist)
{
    long i, j, k1, k2, nrow, rowpos;
    double zaehler, nenner, fa, fb;
    double *ptr_dist, *ptr_freq1, *ptr_freq2;
    int *roff;

    if (*no_mar < 1)
    {
        info(-2, "Invalid number of markers in dic\n");
        return;
    }

    if ((size_t)*no_mar > SIZE_MAX / sizeof(*roff))
    {
        info(-2, "Memory allocation failed in dic\n");
        return;
    }

    roff = (int *)malloc((size_t)*no_mar * sizeof(*roff));

    if (NULL == roff)
    {
        info(-2, "Memory allocation failed in dic\n");
        return;
    }

    roff[0] = 0;
    for (i = 1; i < *no_mar; i++)
        roff[i] = roff[i - 1] + no_all[i - 1];
    nrow = roff[*no_mar - 1] + no_all[*no_mar - 1];

    ptr_dist = dist;

    for (k1 = 0; k1 < (*no_pop - 1); k1++)
        for (k2 = k1 + 1; k2 < (*no_pop); k2++)
        {
            zaehler = nenner = 0;
            for (i = 0; i < *no_mar_u; i++)
            {
                for (j = 0; j < no_all[mar_u[i] - 1]; j++)
                {
                    rowpos = (roff[mar_u[i] - 1] + j);
                    ptr_freq1 = &freq[k1 * nrow + rowpos];
                    ptr_freq2 = &freq[k2 * nrow + rowpos];
                    if ((*ptr_freq1 != -1) && (*ptr_freq2 != -1))
                    {
                        fa = (double)(*ptr_freq1 > 0 ? 1 : 0);
                        fb = (double)(*ptr_freq2 > 0 ? 1 : 0);
                        zaehler += fa * fb;
                        nenner += fa + fb;
                    }
                }
            }
            if (nenner > 0)
                (*ptr_dist) = 2 * zaehler / nenner;
            else
                (*ptr_dist) = 0;
            ptr_dist++;
        }

    free(roff);
}

void sma(int *no_pop, int *no_mar, int *no_all, int *no_mar_u, int *mar_u,
         double *freq, double *dist)
{
    long i, j, k1, k2, nrow, rowpos;
    double zaehler, nenner, fa, fb;
    double *ptr_dist, *ptr_freq1, *ptr_freq2;
    int *roff;

    if (*no_mar < 1)
    {
        info(-2, "Invalid number of markers in sma\n");
        return;
    }

    if ((size_t)*no_mar > SIZE_MAX / sizeof(*roff))
    {
        info(-2, "Memory allocation failed in sma\n");
        return;
    }

    roff = (int *)malloc((size_t)*no_mar * sizeof(*roff));

    if (NULL == roff)
    {
        info(-2, "Memory allocation failed in sma\n");
        return;
    }

    roff[0] = 0;
    for (i = 1; i < *no_mar; i++)
        roff[i] = roff[i - 1] + no_all[i - 1];
    nrow = roff[*no_mar - 1] + no_all[*no_mar - 1];

    ptr_dist = dist;

    for (k1 = 0; k1 < (*no_pop - 1); k1++)
        for (k2 = k1 + 1; k2 < (*no_pop); k2++)
        {
            zaehler = nenner = 0;
            for (i = 0; i < *no_mar_u; i++)
            {
                for (j = 0; j < no_all[mar_u[i] - 1]; j++)
                {
                    rowpos = (roff[mar_u[i] - 1] + j);
                    ptr_freq1 = &freq[k1 * nrow + rowpos];
                    ptr_freq2 = &freq[k2 * nrow + rowpos];
                    if ((*ptr_freq1 != -1) && (*ptr_freq2 != -1))
                    {
                        fa = (double)(*ptr_freq1 > 0 ? 1 : 0);
                        fb = (double)(*ptr_freq2 > 0 ? 1 : 0);
                        zaehler += fa * fb + (1 - fa) * (1 - fb);
                        nenner += 1;
                    }
                }
            }
            if (nenner > 0)
                (*ptr_dist) = zaehler / nenner;
            else
                (*ptr_dist) = 0;
            ptr_dist++;
        }

    free(roff);
}

void gen_dist(int *measure, int *no_pop, int *no_ind, int *no_mar, int *no_all,
              int *no_pop_u, int *pop_u, int *no_ind_u, int *ind_u,
              int *no_mar_u, int *mar_u, double *data, double *freq,
              double *dist, int *s, double *sdev, char **par)
{

    void (*distance)(int *, int *, int *, int *, int *, double *, double *) =
        NULL;
    double *dist_j = NULL, *dist_b = NULL;
    double *sum_dist_sq = NULL, *sum_dist = NULL;
    int *mar_u_j = NULL, *ind_u_boot = NULL, *mar_u_boot = NULL;

    (void)no_pop_u;
    (void)pop_u;

    if ((*no_pop < 1) || (*no_mar < 1))
        goto fail;

    switch (*measure)
    {
    case 1:
        distance = rd;
        break;
    case 2:
        distance = mrd;
        break;
    case 3:
        distance = euc;
        break;
    case 11:
        distance = dic;
        break;
    case 12:
        distance = jac;
        break;
    case 13:
        distance = sma;
        break;
    default:
        distance = euc;
    }

    allele_freq(no_pop, no_ind, no_mar, no_all, no_ind_u, ind_u, data, freq);
    (*distance)(no_pop, no_mar, no_all, no_mar_u, mar_u, freq, dist);

    if (*s == 1)
    {

        int i, j, k, n_dist;
        size_t n_dist_size, n_dist_a, n_dist_b;
        double f = (double)(*no_mar - 1) / (double)(*no_mar);
        int no_mar_u[1];
        *no_mar_u = *no_mar - 1;

        if ((*no_pop < 2) || (*no_mar_u < 1))
            goto fail;

        n_dist_a = (size_t)*no_pop;
        n_dist_b = (size_t)(*no_pop - 1);
        if (0 == n_dist_a % 2)
            n_dist_a /= 2;
        else
            n_dist_b /= 2;
        if (n_dist_a > SIZE_MAX / n_dist_b)
            goto fail;
        n_dist_size = n_dist_a * n_dist_b;
        if (n_dist_size > INT_MAX)
            goto fail;
        n_dist = (int)n_dist_size;

        {
            dist_j = (double *)yRealloc_impl(NULL, n_dist_size, sizeof(double));
            if (NULL == dist_j)
                goto fail;
            sum_dist_sq =
                (double *)yRealloc_impl(NULL, n_dist_size, sizeof(double));
            if (NULL == sum_dist_sq)
                goto fail;
            sum_dist =
                (double *)yRealloc_impl(NULL, n_dist_size, sizeof(double));
            if (NULL == sum_dist)
                goto fail;
            mar_u_j =
                (int *)yRealloc_impl(NULL, (size_t)*no_mar_u, sizeof(int));
            if (NULL == mar_u_j)
                goto fail;

            for (j = 0; j < n_dist; j++)
                sum_dist_sq[j] = sum_dist[j] = 0;

            for (i = 0; i < *no_mar; i++)
            {
                for (j = k = 0; j < *no_mar; j++)
                    if (i != j)
                    {
                        mar_u_j[k] = j + 1;
                        k++;
                    }
                for (j = 0; j < n_dist; j++)
                    dist_j[j] = 0;

                (*distance)(no_pop, no_mar, no_all, no_mar_u, mar_u_j, freq,
                            dist_j);

                for (j = 0; j < n_dist; j++)
                {
                    sum_dist_sq[j] += pow(dist_j[j], 2);
                    sum_dist[j] += dist_j[j];
                }
            }

            for (j = 0; j < n_dist; j++)
                sdev[j] = sqrt(f * (sum_dist_sq[j] -
                                    pow(sum_dist[j], 2) / (double)(*no_mar)));

            yFree(dist_j);
            yFree(sum_dist_sq);
            yFree(sum_dist);
            yFree(mar_u_j);
        }
    }

    if (*s > 1)
    {

        int i, j, k, l, lfdnr, n_dist, sumall = 0;
        size_t n_dist_size, n_dist_a, n_dist_b, sumind = 0;
        double f = (double)1 / (double)(*s);

        if ((*no_pop < 2) || (*no_mar_u < 1))
            goto fail;

        for (i = 0; i < *no_pop; i++)
        {
            if ((no_ind[i] < 0) || ((size_t)no_ind[i] > SIZE_MAX - sumind))
                goto fail;
            sumind += (size_t)no_ind[i];
        }

        if ((0 == sumind) || (sumind > INT_MAX))
            goto fail;

        n_dist_a = (size_t)*no_pop;
        n_dist_b = (size_t)(*no_pop - 1);
        if (0 == n_dist_a % 2)
            n_dist_a /= 2;
        else
            n_dist_b /= 2;
        if (n_dist_a > SIZE_MAX / n_dist_b)
            goto fail;
        n_dist_size = n_dist_a * n_dist_b;
        if (n_dist_size > INT_MAX)
            goto fail;
        n_dist = (int)n_dist_size;

        (void)sumall;

        {

            ind_u_boot = (int *)yRealloc_impl(NULL, sumind, sizeof(int));
            if (NULL == ind_u_boot)
                goto fail;
            mar_u_boot =
                (int *)yRealloc_impl(NULL, (size_t)*no_mar_u, sizeof(int));
            if (NULL == mar_u_boot)
                goto fail;
            dist_b = (double *)yRealloc_impl(NULL, n_dist_size, sizeof(double));
            if (NULL == dist_b)
                goto fail;
            sum_dist_sq =
                (double *)yRealloc_impl(NULL, n_dist_size, sizeof(double));
            if (NULL == sum_dist_sq)
                goto fail;
            sum_dist =
                (double *)yRealloc_impl(NULL, n_dist_size, sizeof(double));
            if (NULL == sum_dist)
                goto fail;

            for (j = 0; j < n_dist; j++)
                sum_dist_sq[j] = sum_dist[j] = 0;

            for (i = 0; i < *s; i++)
            {
                if (strstr(*par, "m"))
                {
                    sumall = 0;
                    for (k = 0; k < *no_mar_u; k++)
                    {
                        mar_u_boot[k] = 1 + st_rand_index(*no_mar_u);
                        sumall += no_all[mar_u_boot[k] - 1];
                    }
                }
                else
                {
                    sumall = 0;
                    for (k = 0; k < *no_mar_u; k++)
                    {
                        mar_u_boot[k] = mar_u[k];
                        sumall += no_all[mar_u_boot[k] - 1];
                    }
                }

                if (strstr(*par, "i"))
                {
                    lfdnr = 0;
                    for (k = 0; k < *no_pop; k++)
                    {
                        for (l = 0; l < no_ind[k]; l++)
                        {
                            ind_u_boot[lfdnr] = 1 + st_rand_index(no_ind[k]);
                            lfdnr++;
                        }
                    }
                }
                else
                {
                    lfdnr = 0;
                    for (k = 0; k < *no_pop; k++)
                    {
                        for (l = 0; l < no_ind[k]; l++)
                        {
                            ind_u_boot[lfdnr] = ind_u[lfdnr];
                            lfdnr++;
                        }
                    }
                }

                for (j = 0; j < n_dist; j++)
                    dist_b[j] = 0;

                allele_freq(no_pop, no_ind, no_mar, no_all, no_ind_u,
                            ind_u_boot, data, freq);
                (*distance)(no_pop, no_mar, no_all, no_mar_u, mar_u_boot, freq,
                            dist_b);

                for (j = 0; j < n_dist; j++)
                {
                    sum_dist_sq[j] += pow(dist_b[j], 2);
                    sum_dist[j] += dist_b[j];
                }
            }

            for (j = 0; j < n_dist; j++)
                sdev[j] = sqrt(
                    f * (sum_dist_sq[j] - pow(sum_dist[j], 2) / (double)(*s)));

            yFree(ind_u_boot);
            yFree(mar_u_boot);
            yFree(dist_b);
            yFree(sum_dist_sq);
            yFree(sum_dist);
        }
    }

    return;

fail:

    yFree(dist_j);
    yFree(dist_b);
    yFree(sum_dist_sq);
    yFree(sum_dist);
    yFree(mar_u_j);
    yFree(ind_u_boot);
    yFree(mar_u_boot);

    info(-2, "Invalid dimensions or memory allocation failed in gen_dist\n");
    return;
}

static int ps_get_simpop_population_name_valid(const char *name)
{
    const unsigned char *p = (const unsigned char *)name;
    size_t n;

    if ((NULL == name) || ('\0' == name[0]))
        return 0;

    n = strlen(name);
    if ((n < 2) || (n >= MAXNAME))
        return 0;

    for (; '\0' != *p; p++)
        if (!(('A' <= *p && *p <= 'Z') || ('a' <= *p && *p <= 'z') ||
              ('0' <= *p && *p <= '9')))
            return 0;

    return 1;
}

static int ps_get_simpop_allele_pair(unsigned short a1, unsigned short a2,
                                     unsigned short bg, GS_ALLELE_TYPE *out1,
                                     GS_ALLELE_TYPE *out2)
{
    if ((bg == a1) || (bg == a2) || (USHRT_MAX == a1) || (USHRT_MAX == a2))
    {
        *out1 = GS_ALLELE_MISSING;
        *out2 = GS_ALLELE_MISSING;
        return 1;
    }

    if ((0 == a1) || (0 == a2) || (a1 > GS_ALLELE_MAX) || (a2 > GS_ALLELE_MAX))
        return 0;

    *out1 = (GS_ALLELE_TYPE)a1;
    *out2 = (GS_ALLELE_TYPE)a2;

    return 1;
}

void ps_get_simpop(char **pop_name, gs_UNSIGNED *NoInd, gs_UNSIGNED *NoMar,
                   char ***IndName, char ***Pedigree,
                   char (**MarName)[gs_MAXNME], GS_ALLELE_TYPE *(*All1),
                   GS_ALLELE_TYPE *(*All2), gs_UNSIGNED *NLoc,
                   gs_mappoint_type **map_S, char **map_str_S,
                   size_t *map_str_N, gs_mappoint_type ***map,
                   gs_UNSIGNED **mdl, int *bg, int *Sretval)
{

    gs_UNSIGNED m;
    size_t data_count;
    size_t map_pool_size = 0;
    size_t map_offset = 0;
    int c, l;
    int allele_error = 0;
    unsigned short bg_allele;

    if ((*bg < 0) || (*bg > USHRT_MAX))
    {
        info(-2, "SIM missing allele indicator is invalid\n");
        *Sretval = -2;
        return;
    }
    bg_allele = (unsigned short)*bg;

    TPopListEl *PopListEl;
    TPop *Pop = NULL;

    find_list_element(*pop_name, &PopListEl);

    if (PopListEl == NULL)
    {
        info(-2, "SIM Population is not defined\n");
        *Sretval = -2;
        return;
    }

    Pop = PopListEl->Pop;

    if (!ps_get_simpop_population_name_valid(Pop->Name))
    {
        info(-2, "SIM population name must contain 2 to 255 characters "
                 "from A-Z, a-z, and 0-9\n");
        *Sretval = -2;
        return;
    }

    if (1 > Pop->NoInds)
    {
        info(-2, "SIM Population is empty\n ");
        *Sretval = -2;
        return;
    }

    if (2 != NoHoms)
    {
        info(-2, "SIM transfer to GS data set implemented for diploids only\n");
        *Sretval = -2;
        return;
    }

    *NoInd = Pop->NoInds;

    int identity_status = st_sim_export_identity(Pop, IndName, Pedigree);
    if (identity_status)
    {
        info(-2, identity_status == -2 ? "Duplicate individual IDs cannot be "
                                         "exported to a marker dataset\n"
                                       : "Cannot export population identity\n");
        goto fail;
    }

    *NoMar = cMP;
    *MarName = (char(*)[gs_MAXNME])yRealloc_impl(*MarName, *NoMar,
                                                 sizeof(char[gs_MAXNME]));
    if (NULL == *MarName)
        goto fail;

    if (NULL == PopDescr)
    {
        info(-2, "SIM No linkage map loaded\n");
        goto fail;
    }

    m = 0;
    for (c = 0; c < NoChroms; c++)
    {
        for (l = 0; l < PopDescr[c].NrMappoints; l++)
        {
            if (!st_identifier_valid(PopDescr[c].Mappoint[l].name) ||
                strlen(PopDescr[c].Mappoint[l].name) >= gs_MAXNME)
            {
                info(-2,
                     "SIM marker name invalid for SelectionTools transfer\n");
                goto fail;
            }
            strcpy((*MarName)[m], PopDescr[c].Mappoint[l].name);
            m++;
        }
    }

    if ((size_t)*NoMar > SIZE_MAX / (size_t)*NoInd)
        goto fail;

    data_count = (size_t)*NoMar * (size_t)*NoInd;

    *All1 = (GS_ALLELE_TYPE *)yRealloc(*All1, data_count, GS_ALLELE_TYPE);
    if (NULL == *All1)
        goto fail;
    *All2 = (GS_ALLELE_TYPE *)yRealloc(*All2, data_count, GS_ALLELE_TYPE);
    if (NULL == *All2)
        goto fail;

#ifdef MT
#pragma omp parallel
#endif

    {

        POPULATION ind;
        int chrom;

        TInd *Ind;
        TChrom *Chrom;
        int i, m;

        size_t idx;

#ifdef MT
#pragma omp for
#endif

        for (ind = 0; ind < Pop->NoInds; ind++)
        {
            Ind = get_individual(PopListEl->Pop, ind);

            for (chrom = 0, m = 0; chrom < NoChroms; chrom++)
            {
                Chrom = &(Ind->Chrom[chrom]);

                for (i = 0; i < PopDescr[chrom].NrMappoints; i++, m++)
                {

                    idx = (size_t)m * (size_t)*NoInd + (size_t)ind;

                    {
                        unsigned short a1 = get_allele(
                            &(Chrom->Hom[0]), PopDescr[chrom].Mappoint[i].Pos);
                        unsigned short a2 = get_allele(
                            &(Chrom->Hom[1]), PopDescr[chrom].Mappoint[i].Pos);

                        if (!ps_get_simpop_allele_pair(a1, a2, bg_allele,
                                                       &((*All1)[idx]),
                                                       &((*All2)[idx])))
                        {
#ifdef MT
#pragma omp atomic write
#endif
                            allele_error = 1;
                            (*All1)[idx] = GS_ALLELE_MISSING;
                            (*All2)[idx] = GS_ALLELE_MISSING;
                        }
                    }
                }
            }
        }
    }

    if (allele_error)
    {
        info(-2, "SIM allele code cannot be represented in GS marker data\n");
        goto fail;
    }

    *NLoc = cMP;

    for (c = 0; c < NoChroms; c++)
    {
        for (l = 0; l < PopDescr[c].NrMappoints; l++)
        {
            const char *name = PopDescr[c].Mappoint[l].name;
            const char *class_name = PopDescr[c].Mappoint[l].class;
            size_t name_len;
            size_t class_len;

            if ((NULL == name) || (NULL == class_name))
            {
                info(-2, "SIM Marker name or class invalid\n");
                goto fail;
            }

            name_len = strlen(name);
            class_len = strlen(class_name);

            if (name_len >= gs_MAXSTR)
            {
                info(-2, "SIM Marker name too long\n");
                goto fail;
            }

            if (class_len >= gs_MAXSTR)
            {
                info(-2, "SIM Marker class too long\n");
                goto fail;
            }

            if ((name_len + 1 > SIZE_MAX - map_pool_size) ||
                (class_len + 1 > SIZE_MAX - (map_pool_size + name_len + 1)))
            {
                info(-2, "SIM linkage map string storage too large\n");
                goto fail;
            }

            map_pool_size += name_len + 1;
            map_pool_size += class_len + 1;
        }
    }

    *map_S = (gs_mappoint_type *)yRealloc(*map_S, *NLoc, gs_mappoint_type);
    if (NULL == *map_S)
        goto fail;
    *map_str_S = (char *)yRealloc(*map_str_S, map_pool_size, char);
    if (NULL == *map_str_S)
        goto fail;
    *map_str_N = map_pool_size;
    *map = (gs_mappoint_type **)yRealloc(*map, *NLoc, gs_mappoint_type *);
    if (NULL == *map)
        goto fail;
    *mdl = (gs_UNSIGNED *)yRealloc(*mdl, *NLoc, gs_UNSIGNED);
    if (NULL == *mdl)
        goto fail;

    for (c = 0, m = 0; c < NoChroms; c++)
    {
        for (l = 0; l < PopDescr[c].NrMappoints; l++, m++)
        {
            const char *name = PopDescr[c].Mappoint[l].name;
            const char *class_name = PopDescr[c].Mappoint[l].class;
            size_t name_len = strlen(name) + 1;
            size_t class_len = strlen(class_name) + 1;

            (*map_S)[m].chrom = (int)PopDescr[c].Mappoint[l].chrom;
            (*map_S)[m].pos = (double)PopDescr[c].Mappoint[l].Pos * 100;

            memcpy(*map_str_S + map_offset, name, name_len);
            (*map_S)[m].name = *map_str_S + map_offset;
            map_offset += name_len;

            memcpy(*map_str_S + map_offset, class_name, class_len);
            (*map_S)[m].class = *map_str_S + map_offset;
            map_offset += class_len;
        }
    }

    if (map_offset != *map_str_N)
        goto fail;

    for (m = 0; m < *NLoc; m++)
    {
        (*map)[m] = &((*map_S)[m]);
        (*mdl)[m] = m;
        (*map_S)[m].nr = m;
    }

    *Sretval = 0;
    return;

fail:

    st_text_array_free(*IndName, *NoInd);
    *IndName = NULL;
    st_text_array_free(*Pedigree, *NoInd);
    *Pedigree = NULL;
    yFree(*MarName);
    yFree(*All1);
    yFree(*All2);
    yFree(*map_S);
    yFree(*map_str_S);
    *map_str_N = 0;
    yFree(*map);
    yFree(*mdl);

    *NoInd = 0;
    *NoMar = 0;
    *NLoc = 0;
    *Sretval = -2;
    return;
}

void ps_get_simpop2(char **pop_name, gs_UNSIGNED *NoInd, gs_UNSIGNED *NoMar,
                    char ***IndName, char ***Pedigree,
                    char (**MarName)[gs_MAXNME], GS_ALLELE_TYPE *(*All1),
                    GS_ALLELE_TYPE *(*All2), gs_UNSIGNED *NLoc,
                    gs_mappoint_type **map_S, gs_mappoint_type ***map,
                    gs_UNSIGNED **mdl, int *bg, int *Sretval)

{
    (void)NLoc;
    (void)map_S;
    (void)map;
    (void)mdl;

    gs_UNSIGNED m;
    size_t data_count;
    int c, l;
    int allele_error = 0;
    unsigned short bg_allele;

    if ((*bg < 0) || (*bg > USHRT_MAX))
    {
        info(-2, "SIM missing allele indicator is invalid\n");
        *Sretval = -2;
        return;
    }
    bg_allele = (unsigned short)*bg;

    TPopListEl *PopListEl;
    TPop *Pop = NULL;

    find_list_element(*pop_name, &PopListEl);

    if (PopListEl == NULL)
    {
        info(-2, "SIM Population is not defined\n");
        *Sretval = -2;
        return;
    }

    Pop = PopListEl->Pop;

    if (!ps_get_simpop_population_name_valid(Pop->Name))
    {
        info(-2, "SIM population name must contain 2 to 255 characters "
                 "from A-Z, a-z, and 0-9\n");
        *Sretval = -2;
        return;
    }

    if (1 > Pop->NoInds)
    {
        info(-2, "SIM Population is empty\n ");
        *Sretval = -2;
        return;
    }

    if (2 != NoHoms)
    {
        info(-2, "SIM transfer to GS data set implemented for diploids only\n");
        *Sretval = -2;
        return;
    }

    *NoInd = Pop->NoInds;

    int identity_status = st_sim_export_identity(Pop, IndName, Pedigree);
    if (identity_status)
    {
        info(-2, identity_status == -2 ? "Duplicate individual IDs cannot be "
                                         "exported to a marker dataset\n"
                                       : "Cannot export population identity\n");
        goto fail;
    }

    *NoMar = cMP;
    *MarName = (char(*)[gs_MAXNME])yRealloc_impl(*MarName, *NoMar,
                                                 sizeof(char[gs_MAXNME]));
    if (NULL == *MarName)
        goto fail;

    if (NULL == PopDescr)
    {
        info(-2, "SIM No linkage map loaded\n");
        goto fail;
    }

    m = 0;
    for (c = 0; c < NoChroms; c++)
    {
        for (l = 0; l < PopDescr[c].NrMappoints; l++)
        {
            if (!st_identifier_valid(PopDescr[c].Mappoint[l].name) ||
                strlen(PopDescr[c].Mappoint[l].name) >= gs_MAXNME)
            {
                info(-2,
                     "SIM marker name invalid for SelectionTools transfer\n");
                goto fail;
            }
            strcpy((*MarName)[m], PopDescr[c].Mappoint[l].name);
            m++;
        }
    }

    if ((size_t)*NoMar > SIZE_MAX / (size_t)*NoInd)
        goto fail;

    data_count = (size_t)*NoMar * (size_t)*NoInd;

    *All1 = (GS_ALLELE_TYPE *)yRealloc(*All1, data_count, GS_ALLELE_TYPE);
    if (NULL == *All1)
        goto fail;
    *All2 = (GS_ALLELE_TYPE *)yRealloc(*All2, data_count, GS_ALLELE_TYPE);
    if (NULL == *All2)
        goto fail;

#ifdef MT
#pragma omp parallel
#endif

    {

        POPULATION ind;
        int chrom;

        TInd *Ind;
        TChrom *Chrom;
        int i, m;

        gs_UNSIGNED idx;

#ifdef MT
#pragma omp for
#endif

        for (ind = 0; ind < Pop->NoInds; ind++)
        {
            Ind = get_individual(PopListEl->Pop, ind);

            for (chrom = 0, m = 0; chrom < NoChroms; chrom++)
            {
                Chrom = &(Ind->Chrom[chrom]);

                for (i = 0; i < PopDescr[chrom].NrMappoints; i++, m++)
                {

                    idx = m * *NoInd + ind;

                    {
                        unsigned short a1 = get_allele(
                            &(Chrom->Hom[0]), PopDescr[chrom].Mappoint[i].Pos);
                        unsigned short a2 = get_allele(
                            &(Chrom->Hom[1]), PopDescr[chrom].Mappoint[i].Pos);

                        if (!ps_get_simpop_allele_pair(a1, a2, bg_allele,
                                                       &((*All1)[idx]),
                                                       &((*All2)[idx])))
                        {
#ifdef MT
#pragma omp atomic write
#endif
                            allele_error = 1;
                            (*All1)[idx] = GS_ALLELE_MISSING;
                            (*All2)[idx] = GS_ALLELE_MISSING;
                        }
                    }
                }
            }
        }
    }

    if (allele_error)
    {
        info(-2, "SIM allele code cannot be represented in GS marker data\n");
        goto fail;
    }

    *Sretval = 0;
    return;

fail:

    st_text_array_free(*IndName, *NoInd);
    *IndName = NULL;
    st_text_array_free(*Pedigree, *NoInd);
    *Pedigree = NULL;
    yFree(*MarName);
    yFree(*All1);
    yFree(*All2);

    *NoInd = 0;
    *NoMar = 0;
    *Sretval = -2;
    return;
}

static int st_copy_sim_mappoints(gs_UNSIGNED NLoc, gs_mappoint_type **map,
                                 TMappoint **points, int *nMP)
{
    gs_UNSIGNED lne;

    *points = NULL;
    *nMP = 0;

    if ((0 == NLoc) || (NLoc > (gs_UNSIGNED)INT_MAX) || (NULL == map))
    {
        info(-2, "Invalid GS linkage map\n");
        return 0;
    }

    for (lne = 0; lne < NLoc; lne++)
    {
        if (NULL == map[lne])
        {
            info(-2, "Invalid GS linkage map pointer\n");
            goto fail;
        }

        if ((NULL == map[lne]->name) || (NULL == map[lne]->class))
        {
            info(-2, "Invalid GS linkage map name or class\n");
            goto fail;
        }

        if (!append_mappoint(points, nMP, map[lne]->chrom,
                             map[lne]->pos / 100.0, map[lne]->name,
                             map[lne]->class))
        {
            info(-2, "Memory allocation failed in st_copy_sim_mappoints\n");
            goto fail;
        }
    }

    return 1;

fail:

    free_mappoints(*points, *nMP);
    *points = NULL;
    *nMP = 0;

    return 0;
}

static int st_sim_allele(GS_ALLELE_TYPE allele, unsigned short *sim_allele)
{
    if (GS_ALLELE_MISSING == allele)
    {
        *sim_allele = USHRT_MAX;
        return 1;
    }

    if ((0 == allele) || (allele > GS_ALLELE_MAX))
        return 0;

    *sim_allele = (unsigned short)allele;

    return 1;
}

void st_set_sim_gp_provide_data(int *n, double **chrlen, int *retval,
                                char **set_name);

void st_set_sim_gp(int *retval, char **set_name)
{
    int n = 0;
    int Sretval = -2;
    int hom = 2;
    double *chrlen = NULL;

    *retval = -2;

    st_set_sim_gp_provide_data(&n, &chrlen, &Sretval, set_name);

    if ((0 != Sretval) || (n < 1) || (NULL == chrlen))
    {
        free(chrlen);
        return;
    }

    if (!set_genome_par_impl(&n, &hom, chrlen))
    {
        free(chrlen);
        return;
    }

    free(chrlen);

    *retval = 0;
    return;
}

void st_set_sim_mp_provide_data(gs_UNSIGNED *NLoc, gs_mappoint_type ***map,
                                int *retval, char **set_name);

void st_set_sim_mp(int *retval, char **set_name)
{
    gs_UNSIGNED NLoc = 0;
    gs_mappoint_type **map = NULL;
    TMappoint *points = NULL;
    int nMP = 0;
    int modified = 0;
    int Sretval = -2;

    *retval = -2;

    if (cMP != 0)
    {
        info(-2, "First remove previously loaded linkage map\n");
        return;
    }

    if (NULL == PopDescr)
    {
        info(-2, "No genome parameters defined\n");
        return;
    }

    st_set_sim_mp_provide_data(&NLoc, &map, &Sretval, set_name);

    if (0 != Sretval)
        return;

    if (!st_copy_sim_mappoints(NLoc, map, &points, &nMP))
        return;

    if (!prepare_mappoints(points, nMP, 1, MAP_DISPERSE_FACTOR_DEFAULT,
                           &modified))
        goto fail;

    if (!install_mappoints(points, nMP))
    {
        info(-2, "Memory allocation failed in st_set_sim_mp\n");
        goto fail;
    }

    free_mappoints(points, nMP);

    if (0 < modified)
        info(0, "%d loci positions modified\n", modified);

    info(0, "%d loci defined\n", nMP);

    *retval = 0;
    return;

fail:

    free_mappoints(points, nMP);
}

void st_set_sim_ef_provide_data(gs_UNSIGNED *NoEff, char (**EffNme)[gs_MAXNME],
                                int **EffAll, int **CmpAll, double **u,
                                int *retval, char **set_name);

static void st_sim_ef_free(TEffects *Effect)
{
    int i, j;

    if (NULL == Effect)
        return;

    if (NULL != Effect->Eff)
    {
        for (i = 0; i < Effect->NrEff; i++)
        {
            if (NULL != Effect->Eff[i].Loc)
            {
                for (j = 0; j < Effect->Eff[i].NrLoci; j++)
                    free(Effect->Eff[i].Loc[j].All);

                free(Effect->Eff[i].Loc);
            }
        }

        free(Effect->Eff);
    }

    free(Effect->MPset);

    Effect->Eff = NULL;
    Effect->MPset = NULL;
    Effect->NrEff = 0;
    Effect->MPs = 0;
}

static int st_sim_ef_finalize(TEffects *Effect)
{
    int i, j, k, l, m, n;
    int maxLoc = 0;

    Effect->MPs = 0;
    free(Effect->MPset);
    Effect->MPset = NULL;

    for (i = 0; i < Effect->NrEff; i++)
    {
        if (maxLoc < Effect->Eff[i].NrLoci)
            maxLoc = Effect->Eff[i].NrLoci;

        for (j = 0; j < Effect->Eff[i].NrLoci; j++)
        {
            TMappoint **tmp;

            for (k = 0; k < Effect->MPs; k++)
                if (Effect->MPset[k] == Effect->Eff[i].Loc[j].MP)
                    break;

            if (k != Effect->MPs)
                continue;

            if (INT_MAX == Effect->MPs)
                return 0;

            tmp = realloc(Effect->MPset,
                          (size_t)(Effect->MPs + 1) * sizeof(*Effect->MPset));

            if (NULL == tmp)
                return 0;

            Effect->MPset = tmp;
            Effect->MPset[Effect->MPs] = Effect->Eff[i].Loc[j].MP;
            Effect->MPs++;
        }
    }

    for (i = 0; i < maxLoc; i++)
    {
        for (k = 0; k < NoHoms; k++)
        {
            int oldEff = Effect->NrEff;

            for (n = k + 1; n < NoHoms; n++)
            {
                for (j = 0; j < oldEff; j++)
                {
                    TEff *tmp;
                    TEff *newEff;

                    if ((i >= Effect->Eff[j].NrLoci) ||
                        (k >= Effect->Eff[j].Loc[i].NrAll) ||
                        (Effect->Eff[j].Loc[i].All[k] ==
                         Effect->Eff[j].Loc[i].All[n]))
                        continue;

                    if (INT_MAX == Effect->NrEff)
                        return 0;

                    tmp = realloc(Effect->Eff, (size_t)(Effect->NrEff + 1) *
                                                   sizeof(*Effect->Eff));

                    if (NULL == tmp)
                        return 0;

                    Effect->Eff = tmp;
                    newEff = &(Effect->Eff[Effect->NrEff]);
                    newEff->Eff = Effect->Eff[j].Eff;
                    newEff->NrLoci = Effect->Eff[j].NrLoci;
                    newEff->Loc = NULL;
                    Effect->NrEff++;

                    if (0 < newEff->NrLoci)
                    {
                        if ((size_t)newEff->NrLoci >
                            SIZE_MAX / sizeof(*newEff->Loc))
                            return 0;

                        newEff->Loc = malloc((size_t)newEff->NrLoci *
                                             sizeof(*newEff->Loc));

                        if (NULL == newEff->Loc)
                            return 0;

                        for (l = 0; l < newEff->NrLoci; l++)
                        {
                            newEff->Loc[l].MP = NULL;
                            newEff->Loc[l].NrAll = 0;
                            newEff->Loc[l].All = NULL;
                        }
                    }

                    for (l = 0; l < newEff->NrLoci; l++)
                    {
                        newEff->Loc[l].MP = Effect->Eff[j].Loc[l].MP;
                        newEff->Loc[l].NrAll = Effect->Eff[j].Loc[l].NrAll;

                        if ((NoHoms < 1) ||
                            ((size_t)NoHoms >
                             SIZE_MAX / sizeof(*newEff->Loc[l].All)))
                            return 0;

                        newEff->Loc[l].All = malloc(
                            (size_t)NoHoms * sizeof(*newEff->Loc[l].All));

                        if (NULL == newEff->Loc[l].All)
                            return 0;

                        for (m = 0; m < NoHoms; m++)
                        {
                            if ((l == i) && (m == k))
                                newEff->Loc[l].All[m] =
                                    Effect->Eff[j].Loc[l].All[n];
                            else if ((l == i) && (m == n))
                                newEff->Loc[l].All[m] =
                                    Effect->Eff[j].Loc[l].All[k];
                            else
                                newEff->Loc[l].All[m] =
                                    Effect->Eff[j].Loc[l].All[m];
                        }
                    }
                }
            }
        }
    }

    return 1;
}

void st_set_sim_ef(int *retval, char **effect_name, char **set_name)
{
    gs_UNSIGNED NoEff = 0;
    char(*EffNme)[gs_MAXNME] = NULL;
    int *EffAll = NULL;
    int *CmpAll = NULL;
    double *u = NULL;
    double gvalue_avg;
    int Sretval = -2;
    TEffects Effect;
    TEffects *tmpEffects;
    gs_UNSIGNED i;
    gs_UNSIGNED nmarker;

    *retval = -2;

    st_set_sim_ef_provide_data(&NoEff, &EffNme, &EffAll, &CmpAll, &u, &Sretval,
                               set_name);

    if (0 != Sretval)
        return;

    if (0 == NoEff)
    {
        *retval = 0;
        return;
    }

    if ((NULL == effect_name) || (NULL == *effect_name) ||
        ('\0' == (*effect_name)[0]) ||
        (strlen(*effect_name) >= sizeof(Effect.Name)))
    {
        info(-2, "Invalid simulation effect name\n");
        return;
    }

    if (NULL != find_Effect(*effect_name))
    {
        info(-2, "This effect is already loaded (same name)\n");
        return;
    }

    if (2 != NoHoms)
    {
        info(-2, "GS effect transfer is implemented for diploids only\n");
        return;
    }

    if ((NULL == PopDescr) || (cMP < 1))
    {
        info(-2, "No simulation linkage map defined\n");
        return;
    }

    if ((NULL == EffNme) || (NULL == EffAll) || (NULL == CmpAll) || (NULL == u))
    {
        info(-2, "Invalid marker effect data\n");
        return;
    }

    memset(&Effect, 0, sizeof(Effect));
    gvalue_avg = u[0];
    Effect.GValue_avg = 0;
    Effect.weight = 1;
    strcpy(Effect.Name, *effect_name);

    nmarker = NoEff - 1;

    if (nmarker > (gs_UNSIGNED)(INT_MAX / 2))
    {
        info(-2, "Too many marker effects\n");
        return;
    }

    if (0 < nmarker)
    {
        int baseEff = 2 * (int)nmarker;

        if ((size_t)baseEff > SIZE_MAX / sizeof(*Effect.Eff))
        {
            info(-2, "Memory allocation failed in st_set_sim_ef\n");
            return;
        }

        Effect.Eff = malloc((size_t)baseEff * sizeof(*Effect.Eff));

        if (NULL == Effect.Eff)
        {
            info(-2, "Memory allocation failed in st_set_sim_ef\n");
            return;
        }

        Effect.NrEff = 0;
    }

    for (i = 1; i < NoEff; i++)
    {
        TMappoint *tmpMP;
        int e;

        tmpMP = find_mappoint(EffNme[i]);

        if (NULL == tmpMP)
        {
            info(-1, "Marker %s not found; effect ignored\n", EffNme[i]);
            continue;
        }

        e = Effect.NrEff;

        Effect.Eff[e].Eff = 0;
        Effect.Eff[e].NrLoci = 0;
        Effect.Eff[e].Loc = NULL;

        Effect.Eff[e + 1].Eff = 0;
        Effect.Eff[e + 1].NrLoci = 0;
        Effect.Eff[e + 1].Loc = NULL;

        Insert_ALL_Loc(&(Effect.Eff[e]), tmpMP, EffAll[i]);
        Effect.Eff[e].Eff = (float)(u[i] / 2.0);

        Insert_ALL_Loc(&(Effect.Eff[e + 1]), tmpMP, CmpAll[i]);
        Effect.Eff[e + 1].Eff = (float)(-u[i] / 2.0);

        gvalue_avg += u[i];

        Effect.NrEff += 2;
    }

    if (0 == Effect.NrEff)
    {
        free(Effect.Eff);
        Effect.Eff = NULL;
    }

    Effect.GValue_avg = (float)gvalue_avg;

    if (!st_sim_ef_finalize(&Effect))
    {
        info(-2, "Memory allocation failed in st_set_sim_ef\n");
        goto fail;
    }

    if (INT_MAX == NrEffects)
    {
        info(-2, "Too many simulation effects\n");
        goto fail;
    }

    if ((size_t)(NrEffects + 1) > SIZE_MAX / sizeof(*Effects))
    {
        info(-2, "Memory allocation failed in st_set_sim_ef\n");
        goto fail;
    }

    tmpEffects = realloc(Effects, (size_t)(NrEffects + 1) * sizeof(*Effects));

    if (NULL == tmpEffects)
    {
        info(-2, "Memory allocation failed in st_set_sim_ef\n");
        goto fail;
    }

    Effects = tmpEffects;
    Effect.Pos = NrEffects + 1;
    Effects[NrEffects] = Effect;
    NrEffects++;

    {
        TPopListEl *el = PopListStart;

        while (NULL != el)
        {
            st_clear_population_gvalue(el->Pop);
            el = el->Next;
        }
    }

    *retval = 0;
    return;

fail:

    st_sim_ef_free(&Effect);
}

typedef struct
{
    const char *name;
    gs_UNSIGNED row;
} st_sim_marker_index_type;

static int st_sim_marker_index_cmp(const void *a, const void *b)
{
    const st_sim_marker_index_type *ia = (const st_sim_marker_index_type *)a;
    const st_sim_marker_index_type *ib = (const st_sim_marker_index_type *)b;

    return strcmp(ia->name, ib->name);
}

static gs_UNSIGNED
st_sim_marker_index_find(const st_sim_marker_index_type *index, gs_UNSIGNED n,
                         const char *name)
{
    st_sim_marker_index_type key;
    const st_sim_marker_index_type *hit;

    key.name = name;
    key.row = 0;

    hit = bsearch(&key, index, (size_t)n, sizeof(*index),
                  st_sim_marker_index_cmp);

    if (NULL == hit)
        return ULONG_MAX;

    return hit->row;
}

void st_set_sim_pp_provide_data(gs_UNSIGNED *NoInd, gs_UNSIGNED *NoMar,
                                char ***IndName, char ***Pedigree,
                                char (**MarName)[gs_MAXNME],
                                GS_ALLELE_TYPE **All1, GS_ALLELE_TYPE **All2,
                                int *retval, char **set_name);

void st_set_sim_pp(int *retval, char **pop_name, char **set_name)
{
    char **IndName = NULL, **Pedigree = NULL;
    gs_UNSIGNED NoInd = 0;
    gs_UNSIGNED NoMar = 0;
    gs_UNSIGNED i, j;
    gs_UNSIGNED *marker_row = NULL;
    char(*MarName)[gs_MAXNME] = NULL;
    GS_ALLELE_TYPE *All1 = NULL;
    GS_ALLELE_TYPE *All2 = NULL;
    st_sim_marker_index_type *marker_index = NULL;
    TMappoint **points = NULL;
    int Sretval = -2;
    double *PB = NULL;
    int *_Ind = NULL;
    int *_Chrom = NULL;
    int *_Hom = NULL;
    double *_Pos = NULL;
    int *_All = NULL;
    size_t idx;
    size_t loci_count;
    size_t row;
    size_t Cnt;
    TPopListEl *checkPopListEl = NULL;
    TPop *Pop = NULL;
    int population_created = 0;
    int chrom;
    int loc;
    char d[64];
    char *cCountC = d;
    char **cCount = &cCountC;

    *retval = -2;

    if (!ps_get_simpop_population_name_valid(*pop_name))
    {
        info(-2, "Population name must contain 2 to 255 characters "
                 "from A-Z, a-z, and 0-9\n");
        return;
    }

    find_list_element(*pop_name, &checkPopListEl);

    if (NULL != checkPopListEl)
    {
        info(-2, "First remove population %s\n", *pop_name);
        return;
    }

    st_set_sim_pp_provide_data(&NoInd, &NoMar, &IndName, &Pedigree, &MarName,
                               &All1, &All2, &Sretval, set_name);

    if (0 != Sretval)
        return;

    if ((2 != NoHoms) || (NULL == PopDescr) || (cMP < 1) || (0 == NoInd) ||
        (0 == NoMar) || (NULL == MarName) || (NULL == All1) || (NULL == All2))
    {
        info(-2, "Simulation map or GS marker data are not available\n");
        return;
    }

    if (NoInd > (gs_UNSIGNED)INT_MAX)
    {
        info(-2, "Too many individuals in GS data set\n");
        return;
    }

    if ((size_t)NoMar > SIZE_MAX / sizeof(*marker_index))
    {
        info(-2, "GS marker data are too large\n");
        return;
    }

    marker_index = malloc((size_t)NoMar * sizeof(*marker_index));

    if (NULL == marker_index)
        goto fail_alloc;

    for (i = 0; i < NoMar; i++)
    {
        if (NULL == memchr(MarName[i], '\0', gs_MAXNME))
        {
            info(-2, "Invalid marker name in GS data set\n");
            goto fail;
        }

        marker_index[i].name = MarName[i];
        marker_index[i].row = i;
    }

    qsort(marker_index, (size_t)NoMar, sizeof(*marker_index),
          st_sim_marker_index_cmp);

    for (i = 1; i < NoMar; i++)
    {
        if (0 == strcmp(marker_index[i - 1].name, marker_index[i].name))
        {
            info(-2, "Duplicate marker name in GS data set: %s\n",
                 marker_index[i].name);
            goto fail;
        }
    }

    if ((size_t)cMP > SIZE_MAX / sizeof(*points))
    {
        info(-2, "Simulation linkage map is too large\n");
        goto fail;
    }

    points = calloc((size_t)cMP, sizeof(*points));
    marker_row = malloc((size_t)cMP * sizeof(*marker_row));
    PB = malloc((size_t)cMP * sizeof(*PB));

    if ((NULL == points) || (NULL == marker_row) || (NULL == PB))
        goto fail_alloc;

    for (chrom = 0; chrom < NoChroms; chrom++)
    {
        for (loc = 0; loc < PopDescr[chrom].NrMappoints; loc++)
        {
            TMappoint *point = &(PopDescr[chrom].Mappoint[loc]);

            if ((point->counter < 0) || (point->counter >= cMP) ||
                (NULL != points[point->counter]))
            {
                info(-2, "Invalid installed simulation linkage map\n");
                goto fail;
            }

            points[point->counter] = point;
        }
    }

    for (i = 0; i < (gs_UNSIGNED)cMP; i++)
    {
        if (NULL == points[i])
        {
            info(-2, "Invalid installed simulation linkage map\n");
            goto fail;
        }

        marker_row[i] =
            st_sim_marker_index_find(marker_index, NoMar, points[i]->name);

        if (ULONG_MAX == marker_row[i])
        {
            info(-2, "Marker %s not found in GS data set\n", points[i]->name);
            goto fail;
        }

        if ((0 == i) || (points[i]->chrom != points[i - 1]->chrom))
        {
            PB[i] = 0.0;
        }
        else
        {
            PB[i] = points[i - 1]->Pos +
                    (points[i]->Pos - points[i - 1]->Pos) / 2.0;
        }
    }

    if ((size_t)cMP > SIZE_MAX / (size_t)NoInd)
    {
        info(-2, "GS population data are too large\n");
        goto fail;
    }

    loci_count = (size_t)cMP * (size_t)NoInd;

    if (loci_count > SIZE_MAX / 2)
    {
        info(-2, "GS population data are too large\n");
        goto fail;
    }

    Cnt = 2 * loci_count;

    if ((Cnt > (size_t)LONG_MAX) || (Cnt > SIZE_MAX / sizeof(double)))
    {
        info(-2, "GS population data are too large\n");
        goto fail;
    }

    _Ind = malloc(Cnt * sizeof(int));
    _Chrom = malloc(Cnt * sizeof(int));
    _Hom = malloc(Cnt * sizeof(int));
    _Pos = malloc(Cnt * sizeof(double));
    _All = malloc(Cnt * sizeof(int));

    if ((NULL == _Ind) || (NULL == _Chrom) || (NULL == _Hom) ||
        (NULL == _Pos) || (NULL == _All))
        goto fail_alloc;

    row = 0;

    for (j = 0; j < NoInd; j++)
    {
        for (i = 0; i < (gs_UNSIGNED)cMP; i++)
        {
            unsigned short a1;
            unsigned short a2;

            idx = (size_t)marker_row[i] * (size_t)NoInd + (size_t)j;

            if (!st_sim_allele(All1[idx], &a1) ||
                !st_sim_allele(All2[idx], &a2))
            {
                info(-2, "Allele code cannot be represented in simulation\n");
                goto fail;
            }

            _Ind[row] = 1 + (int)j;
            _Chrom[row] = points[i]->chrom;
            _Hom[row] = 1;
            _Pos[row] = PB[i];
            _All[row] = (int)a1;
            row++;

            _Ind[row] = 1 + (int)j;
            _Chrom[row] = points[i]->chrom;
            _Hom[row] = 2;
            _Pos[row] = PB[i];
            _All[row] = (int)a2;
            row++;
        }
    }

    snprintf(d, sizeof(d), "%ld", (long)Cnt);

    init_population(pop_name, cCount, _Ind, _Chrom, _Hom, _Pos, _All);

    find_list_element(*pop_name, &checkPopListEl);

    if (NULL == checkPopListEl)
    {
        info(-2, "Population import from GS data failed\n");
        goto fail;
    }

    population_created = 1;
    Pop = checkPopListEl->Pop;

    if (Pop->NoInds != (POPULATION)NoInd)
    {
        info(-2, "Population import from GS data failed\n");
        goto fail_import;
    }

    for (gs_UNSIGNED i = 0; i < NoInd; i++)
        if (st_sim_identity_set(&Pop->Ind[i], IndName[i],
                                Pedigree ? Pedigree[i] : IndName[i]))
            goto fail_import;

    for (j = 0; j < NoInd; j++)
    {
        TInd *Ind = get_individual(Pop, (POPULATION)j);

        for (i = 0; i < (gs_UNSIGNED)cMP; i++)
        {
            unsigned short expected1;
            unsigned short expected2;

            chrom = points[i]->chrom - 1;

            idx = (size_t)marker_row[i] * (size_t)NoInd + (size_t)j;

            if (!st_sim_allele(All1[idx], &expected1) ||
                !st_sim_allele(All2[idx], &expected2))
                goto fail_import;

            if ((get_allele(&(Ind->Chrom[chrom].Hom[0]), points[i]->Pos) !=
                 expected1) ||
                (get_allele(&(Ind->Chrom[chrom].Hom[1]), points[i]->Pos) !=
                 expected2))
            {
                info(-2,
                     "Population import from GS data failed verification\n");
                goto fail_import;
            }
        }
    }

    free(_Ind);
    free(_Chrom);
    free(_Hom);
    free(_Pos);
    free(_All);
    free(PB);
    free(points);
    free(marker_row);
    free(marker_index);

    *retval = 0;
    return;

fail_alloc:

    info(-2, "Memory allocation failed in st_set_sim_pp\n");
    goto fail;

fail_import:

    if (population_created)
        remove_population(pop_name);

fail:

    free(_Ind);
    free(_Chrom);
    free(_Hom);
    free(_Pos);
    free(_All);
    free(PB);
    free(points);
    free(marker_row);
    free(marker_index);
}

#include "SelectionTools02a.h"

SEXP st_simulation_identity(SEXP population, SEXP metadata, SEXP serial)
{
    if (TYPEOF(population) != STRSXP || XLENGTH(population) != 1 ||
        STRING_ELT(population, 0) == NA_STRING)
        Rf_error("population must be one name");
    const char *name = CHAR(STRING_ELT(population, 0));
    TPop *pop = NULL;
    for (TPopListEl *el = PopListStart; el; el = el->Next)
        if (!strcmp(el->Pop->Name, name))
        {
            pop = el->Pop;
            break;
        }
    if (!pop)
        Rf_error("Simulation population does not exist");
    if (!isNull(metadata))
    {
        if (TYPEOF(metadata) != VECSXP || LENGTH(metadata) != 2 ||
            TYPEOF(VECTOR_ELT(metadata, 0)) != STRSXP ||
            TYPEOF(VECTOR_ELT(metadata, 1)) != STRSXP ||
            XLENGTH(VECTOR_ELT(metadata, 0)) != pop->NoInds ||
            XLENGTH(VECTOR_ELT(metadata, 1)) != pop->NoInds)
            Rf_error("Invalid simulation identity metadata");
        for (POPULATION i = 0; i < pop->NoInds; i++)
        {
            SEXP n = STRING_ELT(VECTOR_ELT(metadata, 0), i),
                 p = STRING_ELT(VECTOR_ELT(metadata, 1), i);
            if (n == NA_STRING || p == NA_STRING ||
                !st_identifier_valid(CHAR(n)) || !*CHAR(p))
                Rf_error("Invalid simulation identifier or pedigree");
        }
        char **names = calloc((size_t)pop->NoInds ? (size_t)pop->NoInds : 1,
                              sizeof(*names));
        char **ped =
            calloc((size_t)pop->NoInds ? (size_t)pop->NoInds : 1, sizeof(*ped));
        if (!names || !ped)
        {
            free(names);
            free(ped);
            Rf_error("Cannot stage identity metadata");
        }
        for (POPULATION i = 0; i < pop->NoInds; i++)
        {
            names[i] =
                st_sim_text_dup(CHAR(STRING_ELT(VECTOR_ELT(metadata, 0), i)));
            ped[i] = st_sim_text_dup(CHAR(STRING_ELT(VECTOR_ELT(metadata, 1), i)));
            if (!names[i] || !ped[i])
            {
                for (POPULATION j = 0; j < pop->NoInds; j++)
                {
                    st_sim_pool_free(names[j]);
                    st_sim_pool_free(ped[j]);
                }
                free(names);
                free(ped);
                Rf_error("Cannot stage identity metadata");
            }
        }
        for (POPULATION i = 0; i < pop->NoInds; i++)
        {
            st_sim_pool_free(pop->Ind[i].name);
            st_sim_pool_free(pop->Ind[i].pedigree);
            pop->Ind[i].name = names[i];
            pop->Ind[i].pedigree = ped[i];
            ps_identity_observe(names[i]);
        }
        free(names);
        free(ped);
    }
    if (!isNull(serial))
    {
        if (TYPEOF(serial) != STRSXP || XLENGTH(serial) != 1 ||
            STRING_ELT(serial, 0) == NA_STRING)
            Rf_error("Invalid identity serial state");
        const char *text = CHAR(STRING_ELT(serial, 0));
        uint64_t value = 0;
        if (!*text)
            Rf_error("Invalid identity serial state");
        for (const char *p = text; *p; p++)
        {
            if (*p < '0' || *p > '9' ||
                value > (UINT64_MAX - (uint64_t)(*p - '0')) / 10)
                Rf_error("Invalid identity serial state");
            value = value * 10 + (uint64_t)(*p - '0');
        }
        if (value > st_identity_serial)
            st_identity_serial = value;
    }
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%llu",
             (unsigned long long)st_identity_serial);
    return mkString(buffer);
}

SEXP st_population_cross(SEXP offspring, SEXP left, SEXP right, SEXP number)
{
    if (TYPEOF(offspring) != STRSXP || XLENGTH(offspring) != 1 ||
        STRING_ELT(offspring, 0) == NA_STRING || TYPEOF(left) != STRSXP ||
        XLENGTH(left) != 1 || STRING_ELT(left, 0) == NA_STRING ||
        TYPEOF(right) != STRSXP || XLENGTH(right) != 1 ||
        STRING_ELT(right, 0) == NA_STRING)
        Rf_error("Cross population names must be scalar strings");
    const char *a = CHAR(STRING_ELT(left, 0)), *b = CHAR(STRING_ELT(right, 0));
    if (!strcmp(a, b) || !strcmp(a, "name") || !strcmp(b, "name"))
        Rf_error("Crossing-design output requires distinct parental pool names "
                 "other than 'name'");
    double value = asReal(number);
    if (!R_FINITE(value) || value < 1 || value > INT_MAX ||
        value != floor(value))
        Rf_error("Number of hybrids must be a positive supported integer");
    int n = (int)value;
    const char **first = (const char **)st_work_alloc(n, sizeof(*first)),
               **second = (const char **)st_work_alloc(n, sizeof(*second));
    memset(first, 0, (size_t)n * sizeof(*first));
    memset(second, 0, (size_t)n * sizeof(*second));
    st_cross_left = first;
    st_cross_right = second;
    st_cross_first = st_identity_serial + 1;
    st_cross_count = (size_t)n;
    char count[32];
    snprintf(count, sizeof(count), "%d", n);
    char *pg = (char *)CHAR(STRING_ELT(offspring, 0)), *p1 = (char *)a,
         *p2 = (char *)b, *np = count;
    double self = 0, significant = -1;
    int mode = 1, scheme = 0;
    cross(&pg, &p1, &p2, &self, &np, &significant, &mode, &scheme);
    st_cross_left = st_cross_right = NULL;
    st_cross_count = 0;
    TPop *pop = NULL;
    for (TPopListEl *el = PopListStart; el; el = el->Next)
        if (!strcmp(el->Pop->Name, pg))
        {
            pop = el->Pop;
            break;
        }
    if (!pop || pop->NoInds != n)
        Rf_error("Cross did not produce the requested offspring population");
    for (int i = 0; i < n; i++)
        if (!first[i] || !second[i])
            Rf_error("Crossing-design capture was incomplete");
    SEXP out = PROTECT(allocVector(VECSXP, 3)),
         names = PROTECT(allocVector(STRSXP, 3)),
         ids = PROTECT(allocVector(STRSXP, n)),
         l = PROTECT(allocVector(STRSXP, n)),
         r = PROTECT(allocVector(STRSXP, n)),
         rn = PROTECT(allocVector(INTSXP, 2));
    for (int i = 0; i < n; i++)
    {
        SET_STRING_ELT(ids, i, mkCharCE(pop->Ind[i].name, CE_UTF8));
        SET_STRING_ELT(l, i, mkCharCE(first[i], CE_UTF8));
        SET_STRING_ELT(r, i, mkCharCE(second[i], CE_UTF8));
    }
    SET_VECTOR_ELT(out, 0, ids);
    SET_VECTOR_ELT(out, 1, l);
    SET_VECTOR_ELT(out, 2, r);
    SET_STRING_ELT(names, 0, mkChar("name"));
    SET_STRING_ELT(names, 1, STRING_ELT(left, 0));
    SET_STRING_ELT(names, 2, STRING_ELT(right, 0));
    setAttrib(out, R_NamesSymbol, names);
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -n;
    setAttrib(out, R_RowNamesSymbol, rn);
    setAttrib(out, R_ClassSymbol, mkString("data.frame"));
    UNPROTECT(6);
    return out;
}

static SEXP st_sim_list(int n, const char **keys, int *np)
{
    SEXP x = PROTECT(allocVector(VECSXP, n)),
         names = PROTECT(allocVector(STRSXP, n));
    *np += 2;
    for (int i = 0; i < n; i++)
        SET_STRING_ELT(names, i, mkChar(keys[i]));
    setAttrib(x, R_NamesSymbol, names);
    return x;
}
static void st_sim_frame(SEXP x, int n, int *np)
{
    SEXP rows = PROTECT(allocVector(INTSXP, n ? 2 : 0)),
         cls = PROTECT(mkString("data.frame"));
    *np += 2;
    if (n)
    {
        INTEGER(rows)[0] = NA_INTEGER;
        INTEGER(rows)[1] = -n;
    }
    setAttrib(x, R_RowNamesSymbol, rows);
    setAttrib(x, R_ClassSymbol, cls);
}
static SEXP st_sim_map_snapshot(int *np)
{
    const char *keys[] = {"chrom", "pos", "name", "class"};
    SEXP out = st_sim_list(4, keys, np);
    for (int j = 0; j < 4; j++)
    {
        SEXP col = PROTECT(allocVector(j < 2 ? REALSXP : STRSXP, cMP));
        (*np)++;
        SET_VECTOR_ELT(out, j, col);
    }
    int i = 0;
    for (int c = 0; c < NoChroms; c++)
        for (int j = 0; j < PopDescr[c].NrMappoints; j++, i++)
        {
            TMappoint *m = PopDescr[c].Mappoint + j;
            REAL(VECTOR_ELT(out, 0))[i] = c + 1;
            REAL(VECTOR_ELT(out, 1))[i] = m->Pos;
            SET_STRING_ELT(VECTOR_ELT(out, 2), i, mkCharCE(m->name, CE_UTF8));
            SET_STRING_ELT(VECTOR_ELT(out, 3), i, mkCharCE(m->class, CE_UTF8));
        }
    st_sim_frame(out, cMP, np);
    return out;
}
SEXP st_population_marker_map(SEXP number)
{
    if ((TYPEOF(number) != REALSXP && TYPEOF(number) != INTSXP) ||
        XLENGTH(number) != 1 || !validate_map_genome())
        error("n.markers must cover every chromosome of a defined genome");
    double value = asReal(number), total = 0;
    if (!R_FINITE(value) || value != floor(value) || value < NoChroms ||
        value > INT_MAX / NoHoms)
        error("n.markers must cover every chromosome within supported "
              "dimensions");
    int n = (int)value;
    for (int c = 0; c < NoChroms; c++)
        total += PopDescr[c].ChromLen;
    if (!R_FINITE(total) || total <= 0)
        error("Invalid genome length");
    TMappoint *points = calloc((size_t)n, sizeof(*points));
    if (!points)
        error("Cannot allocate marker map");

    for (int i = 0; i < n; i++)
    {
        char label[32];
        snprintf(label, sizeof(label), "M%d", i + 1);
        points[i].name = st_text_dup(label);
        points[i].class = st_text_dup("m");
        if (!points[i].name || !points[i].class)
        {
            free_mappoints(points, n);
            error("Cannot allocate marker labels");
        }
    }
    GetRNGstate();
    for (int i = 0; i < n; i++)
    {
        int c = i;
        if (i >= NoChroms)
        {
            double draw = unif_rand() * total, limit = PopDescr[0].ChromLen;
            c = 0;
            while (c + 1 < NoChroms && draw >= limit)
                limit += PopDescr[++c].ChromLen;
        }
        points[i].chrom = c + 1;
        points[i].Pos = unif_rand() * PopDescr[c].ChromLen;
    }
    PutRNGstate();
    qsort(points, n, sizeof(*points), mccmp);
    if (!install_mappoints(points, n))
    {
        free_mappoints(points, n);
        error("Cannot install marker map");
    }
    free_mappoints(points, n);
    int np = 0;
    SEXP out = st_sim_map_snapshot(&np);
    UNPROTECT(np);
    return out;
}

typedef struct
{
    TMappoint *points;
    int n;
} st_map_extension;
static void st_map_extension_free(void *pointer)
{
    st_map_extension *stage = pointer;
    if (stage)
    {
        free_mappoints(stage->points, stage->n);
        free(stage);
    }
}
static SEXP st_sim_points_snapshot(const TMappoint *points, int n, int *np)
{
    const char *keys[] = {"chrom", "pos", "name", "class"};
    SEXP out = st_sim_list(4, keys, np);
    for (int j = 0; j < 4; j++)
    {
        SEXP col = PROTECT(allocVector(j < 2 ? REALSXP : STRSXP, n));
        (*np)++;
        SET_VECTOR_ELT(out, j, col);
    }
    for (int i = 0; i < n; i++)
    {
        REAL(VECTOR_ELT(out, 0))[i] = points[i].chrom;
        REAL(VECTOR_ELT(out, 1))[i] = points[i].Pos;
        SET_STRING_ELT(VECTOR_ELT(out, 2), i,
                       mkCharCE(points[i].name, CE_UTF8));
        SET_STRING_ELT(VECTOR_ELT(out, 3), i,
                       mkCharCE(points[i].class, CE_UTF8));
    }
    st_sim_frame(out, n, np);
    return out;
}
SEXP st_extend_marker_map(SEXP number)
{
    if ((TYPEOF(number) != REALSXP && TYPEOF(number) != INTSXP) ||
        XLENGTH(number) != 1 || !validate_map_genome() || cMP < 1)
        error("Map extension requires an existing map and one total marker "
              "count");
    double value = asReal(number), total = 0;
    if (!R_FINITE(value) || value != floor(value) || value < cMP ||
        value > INT_MAX / NoHoms ||
        value > (double)(SIZE_MAX / sizeof(TMappoint)))
        error("n.markers must be at least the existing map size within "
              "supported dimensions");
    int n = (int)value, old = cMP, np = 0;
    if (n == old)
    {
        SEXP out = st_sim_map_snapshot(&np);
        UNPROTECT(np);
        return out;
    }
    for (int c = 0; c < NoChroms; c++)
        total += PopDescr[c].ChromLen;
    if (!R_FINITE(total) || total <= 0)
        error("Invalid genome length");
    st_map_extension *stage = calloc(1, sizeof(*stage));
    if (!stage)
        error("Cannot stage marker map");
    st_work_adopt(stage, st_map_extension_free);
    stage->n = n;
    stage->points = calloc((size_t)n, sizeof(*stage->points));
    if (!stage->points)
        error("Cannot allocate extended marker map");
    TMappoint *points = stage->points,
              *names = st_work_alloc(old, sizeof(*names));
    int at = 0;
    for (int c = 0; c < NoChroms; c++)
        for (int j = 0; j < PopDescr[c].NrMappoints; j++)
        {
            if (at >= old)
                error("Inconsistent existing map dimensions");
            TMappoint *source = PopDescr[c].Mappoint + j, *dest = points + at;
            *dest = *source;
            dest->chrom = c + 1;
            dest->counter = at;
            dest->name = NULL;
            dest->class = NULL;
            dest->name = st_text_dup(source->name);
            dest->class = st_text_dup(source->class);
            if (!dest->name || !dest->class)
                error("Cannot copy existing marker labels");
            names[at++] = *dest;
        }
    if (at != old)
        error("Inconsistent existing map dimensions");
    qsort(names, old, sizeof(*names), mnamecmp);
    for (int i = 1; i < old; i++)
        if (!strcmp(names[i - 1].name, names[i].name))
            error("Existing marker names must be unique");

    uint64_t serial = 0;
    for (int i = old; i < n; i++)
    {
        char label[32];
        TMappoint key = {0};
        do
        {
            snprintf(label, sizeof(label), "M%llu",
                     (unsigned long long)++serial);
            key.name = label;
        } while (bsearch(&key, names, old, sizeof(*names), mnamecmp));
        points[i].counter = i;
        points[i].name = st_text_dup(label);
        points[i].class = st_text_dup("m");
        if (!points[i].name || !points[i].class)
            error("Cannot allocate marker labels");
    }
    GetRNGstate();
    for (int i = old; i < n; i++)
    {
        double draw = unif_rand() * total, limit = PopDescr[0].ChromLen;
        int c = 0;
        while (c + 1 < NoChroms && draw >= limit)
            limit += PopDescr[++c].ChromLen;
        points[i].chrom = c + 1;
        points[i].Pos = unif_rand() * PopDescr[c].ChromLen;
    }
    PutRNGstate();
    qsort(points, n, sizeof(*points), mccmp);
    SEXP out = st_sim_points_snapshot(points, n, &np);
    if (!install_mappoints(points, n))
        error("Cannot install extended marker map");
    UNPROTECT(np);
    return out;
}

#include <R_ext/Lapack.h>

static SEXP st_legacy_argument(SEXP args, const char *key)
{
    SEXP names = getAttrib(args, R_NamesSymbol);
    if (TYPEOF(args) != VECSXP || TYPEOF(names) != STRSXP)
        error("Legacy interface arguments must be named");
    for (R_xlen_t i = 0; i < XLENGTH(args); i++)
        if (!strcmp(CHAR(STRING_ELT(names, i)), key))
            return VECTOR_ELT(args, i);
    return R_NilValue;
}

static double st_legacy_number(SEXP x, R_xlen_t i)
{
    if (TYPEOF(x) == REALSXP)
        return REAL(x)[i];
    if (TYPEOF(x) == INTSXP || TYPEOF(x) == LGLSXP)
        return INTEGER(x)[i] == NA_INTEGER ? NA_REAL : INTEGER(x)[i];
    error("Numeric input required");
    return NA_REAL;
}

static char *st_legacy_text(const char *source)
{
    char *out = st_work_alloc(strlen(source) + 1, 1);
    strcpy(out, source);
    return out;
}

static double st_legacy_mean(const double *x, size_t n)
{
    if (!n)
        return R_NaN;
    long double total = 0, correction = 0;
    for (size_t i = 0; i < n; i++)
        total += x[i];
    double mean = (double)(total / n);
    if (R_FINITE(mean))
    {
        for (size_t i = 0; i < n; i++)
            correction += (long double)x[i] - mean;
        mean += (double)(correction / n);
    }
    return mean;
}

static int st_legacy_double_compare(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void st_legacy_summary_values(double *x, size_t n, double out[7])
{
    out[0] = (double)n;
    out[1] = st_legacy_mean(x, n);
    long double ss = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (ISNAN(x[i]))
            error("missing values and NaN's not allowed if 'na.rm' is FALSE");
        long double d = (long double)x[i] - out[1];
        ss += d * d;
    }
    out[2] = n > 1 ? sqrt((double)(ss / (n - 1))) : NA_REAL;
    qsort(x, n, sizeof(double), st_legacy_double_compare);
    out[3] = n ? x[0] : R_PosInf;
    out[6] = n ? x[n - 1] : R_NegInf;
    for (int j = 0; j < 2; j++)
    {
        double index = n ? (n - 1) * (j ? 0.5 : 0.1) : 0;
        size_t lo = (size_t)floor(index), hi = (size_t)ceil(index);
        double fraction = index - lo;
        out[4 + j] = !n ? NA_REAL
                     : lo == hi || x[lo] == x[hi]
                         ? x[lo]
                         : (1 - fraction) * x[lo] + fraction * x[hi];
    }
}

static SEXP st_legacy_population_summary(SEXP args)
{
    SEXP populations = st_legacy_argument(args, "pops");
    const char *kind = CHAR(asChar(st_legacy_argument(args, "kind")));
    if (TYPEOF(populations) != STRSXP || XLENGTH(populations) > INT_MAX)
        error("pops must be character population names");
    int count = (int)XLENGTH(populations);
    int allele = asInteger(st_legacy_argument(args, "allele"));
    SEXP chromosome = st_legacy_argument(args, "chromosome");
    int chrom = isNull(chromosome) ? -1 : asInteger(chromosome) - 1;
    SEXP begin_arg = st_legacy_argument(args, "begin");
    SEXP end_arg = st_legacy_argument(args, "end");
    double begin = isNull(begin_arg) ? 0 : asReal(begin_arg);
    if (chrom >= NoChroms || chrom < -1)
        error("Chromosome out of range");
    double end = isNull(end_arg) ? (chrom >= 0 ? PopDescr[chrom].ChromLen : 0)
                                 : asReal(end_arg);
    SEXP out = PROTECT(allocMatrix(REALSXP, count, 7));
    for (int p = 0; p < count; p++)
    {
        if (STRING_ELT(populations, p) == NA_STRING)
            error("Population names must not be missing");
        char *name = st_legacy_text(CHAR(STRING_ELT(populations, p)));
        TPopListEl *el = NULL;
        find_list_element(name, &el);
        if (!el)
            error("Population %s not defined", name);
        size_t n = (size_t)el->Pop->NoInds;
        void *mark = st_work_mark();
        double *value = st_work_alloc(n ? n : 1, sizeof(double));
        memset(value, 0, n * sizeof(double));
        if (!strcmp(kind, "gvalue"))
        {
            char *evaluation = st_legacy_text(name);
            genotype_population(&evaluation);
            evaluation = st_legacy_text(name);
            evaluate_population(&evaluation, NULL);
            SEXP effect_arg = st_legacy_argument(args, "effect");
            char *effect =
                isNull(effect_arg) ? "" : (char *)CHAR(asChar(effect_arg));
            char *criterion = "G";
            int status = -2;
            get_population_value(&name, value, &effect, &criterion, &status);
        }
        else if (!strcmp(kind, "drag"))
        {
            int chromosome1 = chrom + 1;
            double position = asReal(st_legacy_argument(args, "pos"));
            evaluate_ld(&name, &allele, &chromosome1, &position, value);
        }
        else
        {
            double *length = st_work_alloc(n ? n : 1, sizeof(double));
            double *total = st_work_alloc(n ? n : 1, sizeof(double));
            int *blocks = st_work_alloc(n ? n : 1, sizeof(int));
            memset(length, 0, n * sizeof(double));
            memset(total, 0, n * sizeof(double));
            memset(blocks, 0, n * sizeof(int));
            if (chrom < 0)
                evaluate_genome(&name, &allele, value, length, total, blocks);
            else
                evaluate_genome_region(&name, &allele, &chrom, &begin, &end,
                                       value, length, total, blocks);
            if (!strcmp(kind, "segments"))
                for (size_t i = 0; i < n; i++)
                    value[i] = blocks[i];
        }
        double summary[7];
        st_legacy_summary_values(value, n, summary);
        for (int j = 0; j < 7; j++)
            REAL(out)[p + (size_t)count * j] = summary[j];
        st_work_reset(mark);
    }
    const char *labels[] = {"Obs", "Mean", "SDev", "Min", "Q10", "Med", "Max"};
    SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
    SEXP cols = PROTECT(allocVector(STRSXP, 7));
    for (int i = 0; i < 7; i++)
        SET_STRING_ELT(cols, i, mkChar(labels[i]));
    SET_VECTOR_ELT(dimnames, 0, populations);
    SET_VECTOR_ELT(dimnames, 1, cols);
    setAttrib(out, R_DimNamesSymbol, dimnames);
    UNPROTECT(3);
    return out;
}

static void st_legacy_disperse(TMappoint *points, int n, double factor)
{
    int changed = 0;
    if (!disperse_mappoints_impl(points, n, factor, &changed, 0))
        error("disperse.factor must be positive and finite");
    info(0, "%d loci positions modified\n", changed);
}

static SEXP st_legacy_map(SEXP args, int define_genome)
{
    SEXP map = st_legacy_argument(args, "map");
    if (!inherits(map, "data.frame") || LENGTH(map) < 3)
        error(
            "map must be a data frame with chromosome, position and name columns");
    SEXP chrom = VECTOR_ELT(map, 0), position = VECTOR_ELT(map, 1);
    SEXP names = VECTOR_ELT(map, 2);
    SEXP classes = LENGTH(map) >= 4 ? VECTOR_ELT(map, 3) : R_NilValue;
    if (!isNumeric(chrom) || !isNumeric(position) || TYPEOF(names) != STRSXP ||
        (!isNull(classes) && TYPEOF(classes) != STRSXP) || !XLENGTH(chrom) ||
        XLENGTH(chrom) > INT_MAX / (NoHoms ? NoHoms : 2))
        error("Invalid linkage-map columns");
    int n = (int)XLENGTH(chrom);
    if (XLENGTH(position) != n || XLENGTH(names) != n ||
        (!isNull(classes) && XLENGTH(classes) != n))
        error("Linkage-map columns have different lengths");
    st_map_extension *stage = calloc(1, sizeof(*stage));
    if (!stage)
        error("Cannot allocate linkage-map staging");
    st_work_adopt(stage, st_map_extension_free);
    stage->n = n;
    stage->points = calloc(n, sizeof(TMappoint));
    if (!stage->points)
        error("Cannot allocate linkage-map points");
    TMappoint *points = stage->points;
    int *identifiers = st_work_alloc(n, sizeof(int));
    for (int i = 0; i < n; i++)
    {
        double chromosome = st_legacy_number(chrom, i);
        double pos = st_legacy_number(position, i);
        if (!R_FINITE(chromosome) || chromosome != floor(chromosome) ||
            chromosome < 1 || chromosome > INT_MAX || !R_FINITE(pos) ||
            pos < 0 || STRING_ELT(names, i) == NA_STRING)
            error("Invalid linkage-map value");
        points[i].chrom = identifiers[i] = (int)chromosome;
        points[i].Pos = pos;
        points[i].counter = i;
        points[i].name = st_text_dup(CHAR(STRING_ELT(names, i)));
        if (isNull(classes))
        {
            const char *name = CHAR(STRING_ELT(names, i));
            size_t name_length = strlen(name);
            points[i].class = malloc(name_length + 6);
            if (points[i].class)
            {
                memcpy(points[i].class, "Class", 5);
                memcpy(points[i].class + 5, name, name_length + 1);
            }
        }
        else
        {
            if (STRING_ELT(classes, i) == NA_STRING)
                error("Missing linkage-map class");
            points[i].class = st_text_dup(CHAR(STRING_ELT(classes, i)));
        }
        if (!points[i].name || !points[i].class)
            error("Cannot allocate linkage-map labels");
    }
    if (define_genome)
    {
        R_isort(identifiers, n);
        int nc = 0;
        for (int i = 0; i < n; i++)
            if (!i || identifiers[i] != identifiers[i - 1])
                identifiers[nc++] = identifiers[i];
        for (int i = 0; i < n; i++)
            for (int c = 0; c < nc; c++)
                if (points[i].chrom == identifiers[c])
                {
                    points[i].chrom = c + 1;
                    break;
                }
        qsort(points, n, sizeof(TMappoint), mccmp);
        SEXP lengths_arg = st_legacy_argument(args, "lengths");
        if (!isNull(lengths_arg) &&
            (!isNumeric(lengths_arg) || XLENGTH(lengths_arg) != nc))
            error(
                "Vector provided for chromosome length must match number of chromosomes");
        double *lengths = st_work_alloc(nc, sizeof(double));
        memset(lengths, 0, nc * sizeof(double));
        for (int i = 0; i < n; i++)
            lengths[points[i].chrom - 1] =
                fmax(lengths[points[i].chrom - 1], points[i].Pos);
        int increased = 0;
        for (int c = 0; c < nc; c++)
        {
            if (!isNull(lengths_arg))
                lengths[c] = st_legacy_number(lengths_arg, c);
            if (lengths[c] == 0)
                lengths[c] = 0.01, increased = 1;
            if (!R_FINITE(lengths[c]) || lengths[c] <= 0)
            {
                info(-2, "Chromosome length must be finite and positive\n");
                return R_NilValue;
            }
        }
        if (increased)
            warning("Chromosome length cannot be 0; increased to 0.01");
        while (PopListStart)
        {
            info(1, "Population %s deleted\n", PopListStart->Pop->Name);
            st_remove_population_entry(PopListStart);
        }
        info(1, "All populations removed\n");
        remove_map();
        int homs = 2;
        if (!set_genome_par_impl(&nc, &homs, lengths))
            return R_NilValue;
    }
    for (int i = 0; i < n; i++)
        if (points[i].chrom > NoChroms)
            error("Chromosome is outside the defined genome");
    int disperse =
        !define_genome && asLogical(st_legacy_argument(args, "disperse")) == 1;
    if (disperse)
        st_legacy_disperse(points, n,
                           asReal(st_legacy_argument(args, "factor")));
    for (TPopListEl *el = PopListStart; el; el = el->Next)
    {
        st_clear_population_genotype(el->Pop);
        st_clear_population_gvalue(el->Pop);
        st_clear_population_pvalue(el->Pop);
    }
    effmap_remove_all();
    SEXP out = R_NilValue;
    int np = 0;
    if (disperse)
    {
        out = PROTECT(duplicate(map));
        np++;
        SEXP pos = PROTECT(allocVector(REALSXP, n));
        np++;
        for (int i = 0; i < n; i++)
            REAL(pos)[i] = points[i].Pos;
        SET_VECTOR_ELT(out, 1, pos);
    }
    if (!install_mappoints(points, n))
        error("Cannot install linkage map");
    UNPROTECT(np);
    return out;
}

static const char *st_legacy_prefix(const char *name)
{
    const char *end = strchr(name, '.');
    if (!end)
        end = strchr(name, '_');
    size_t n = end ? (size_t)(end - name) : strlen(name);
    char *out = st_work_alloc(n + 1, 1);
    memcpy(out, name, n);
    out[n] = '\0';
    return out;
}

typedef struct
{
    int nr, nc, npop, nmar, frame;
    int *no_all, *no_ind, *individuals, *markers;
    const char **marker, **population;
    SEXP source, rownames, colnames;
    double *values;
    const char **row_labels;
} st_legacy_data;

static int st_legacy_name_compare(const void *left, const void *right)
{
    return strcmp(*(const char *const *)left, *(const char *const *)right);
}

static int st_legacy_primer_position(const char *text, char target)
{
    if (!text) return NA_INTEGER;
    int character = 1;
    for (const unsigned char *p = (const unsigned char *)text; *p; p++)
    {
        if (*p == (unsigned char)target) return character;
        if ((*p & 0xc0) != 0x80) character++;
    }
    return -1;
}

static const char *st_legacy_primer_substring(const char *text, int first, int last)
{
    if (!text || first == NA_INTEGER || last == NA_INTEGER) return "NA";
    if (first < 1) first = 1;
    const char *start = text, *end = text;
    int at = 1;
    while (*start && at < first)
    {
        start++;
        while (((unsigned char)*start & 0xc0) == 0x80) start++;
        at++;
    }
    end = start;
    while (*end && at <= last)
    {
        end++;
        while (((unsigned char)*end & 0xc0) == 0x80) end++;
        at++;
    }
    size_t n = (size_t)(end - start);
    char *result = st_work_alloc(n + 1, 1);
    memcpy(result, start, n); result[n] = '\0';
    return result;
}

static const char **st_legacy_primer_labels(SEXP data)
{
    int frame = inherits(data, "data.frame"), nr = 0;
    SEXP names = R_NilValue;
    if (frame)
    {
        nr = LENGTH(data) ? LENGTH(VECTOR_ELT(data, 0)) : 0;
        names = getAttrib(data, R_RowNamesSymbol);
    }
    else
    {
        SEXP dim = getAttrib(data, R_DimSymbol), dn = getAttrib(data, R_DimNamesSymbol);
        if (TYPEOF(dim) == INTSXP && LENGTH(dim) >= 1) nr = INTEGER(dim)[0];
        if (TYPEOF(dn) == VECSXP && LENGTH(dn)) names = VECTOR_ELT(dn, 0);
    }
    int n = isNull(names) ? 0 : frame ? nr : LENGTH(names);
    const char **labels = st_work_alloc(n ? n : 1, sizeof(*labels));
    for (int i = 0; i < n; i++)
    {
        if (TYPEOF(names) == STRSXP)
            labels[i] = STRING_ELT(names, i) == NA_STRING ? NULL : translateCharUTF8(STRING_ELT(names, i));
        else
        {
            char number[32]; snprintf(number, sizeof(number), "%d", i + 1);
            labels[i] = st_legacy_text(number);
        }
    }
    for (int stage = 0; stage < 2; stage++)
    {
        int *positions = st_work_alloc(n, sizeof(int)), *lengths = st_work_alloc(n, sizeof(int)),
            *good = st_work_alloc(n, sizeof(int)), ngood = 0;
        for (int i = 0; i < n; i++)
        {
            positions[i] = st_legacy_primer_position(labels[i], stage ? '+' : '.');
            lengths[i] = positions[i] == NA_INTEGER ? NA_INTEGER : positions[i] < 0 ? -1 : 1;
            if (positions[i] > 0 || positions[i] == NA_INTEGER) good[ngood++] = positions[i];
        }
        if (n && !ngood) error("invalid substring arguments");
        if (ngood && n % ngood) warning("longer object length is not a multiple of shorter object length");
        const char **next = st_work_alloc(n ? n : 1, sizeof(*next));
        for (int i = 0; i < n; i++)
        {
            int position = good[i % ngood], last = position == NA_INTEGER ? NA_INTEGER : position - 1,
                first = position == NA_INTEGER || lengths[i] == NA_INTEGER ? NA_INTEGER : position + lengths[i];
            const char *left = st_legacy_primer_substring(labels[i], 1, last),
                       *right = st_legacy_primer_substring(labels[i], first, 1000000);
            size_t length = strlen(left) + strlen(right) + (stage ? 2 : 1);
            char *joined = st_work_alloc(length, 1);
            snprintf(joined, length, stage ? "%s.%s" : "%s%s", left, right);
            next[i] = joined;
        }
        if (!n && stage) next[0] = ".", n = 1;
        labels = next;
    }
    if (n != nr)
        error(frame ? "invalid 'row.names' length" : "length of 'dimnames' [1] not equal to array extent");
    if (frame)
    {
        int duplicated = 0, unique = 0, *indices = st_work_alloc(n, sizeof(int));
        for (int i = 0; i < n; i++)
        {
            int duplicate = 0;
            for (int j = 0; j < i; j++) duplicate |= !strcmp(labels[i], labels[j]);
            if (!duplicate) continue;
            duplicated++;
            int seen = 0;
            for (int j = 0; j < unique; j++) seen |= !strcmp(labels[i], labels[indices[j]]);
            if (!seen) indices[unique++] = i;
        }
        if (duplicated)
        {
            SEXP values = PROTECT(allocVector(STRSXP, unique));
            for (int i = 0; i < unique; i++) SET_STRING_ELT(values, i, mkCharCE(labels[indices[i]], CE_UTF8));
            SEXP call = PROTECT(lang2(findFun(install("sQuote"), R_BaseEnv), values)),
                 quoted = PROTECT(eval(call, R_BaseEnv));
            const char *prefix = duplicated == 1 ? "non-unique value when setting 'row.names': "
                                                 : "non-unique values when setting 'row.names': ";
            size_t size = strlen(prefix) + 1 + (size_t)unique * 2;
            for (int i = 0; i < unique; i++) size += strlen(CHAR(STRING_ELT(quoted, i)));
            char *message = st_work_alloc(size, 1); strcpy(message, prefix);
            for (int i = 0; i < unique; i++)
            {
                if (i) strcat(message, ", ");
                strcat(message, CHAR(STRING_ELT(quoted, i)));
            }
            warning("%s", message);
            UNPROTECT(3);
            error("duplicate 'row.names' are not allowed");
        }
    }
    return labels;
}

static st_legacy_data st_legacy_data_read_labels(SEXP data, int missing, const char **row_labels)
{
    SEXP dims = getAttrib(data, R_DimSymbol);
    int nr, nc, frame = inherits(data, "data.frame");
    SEXP rownames, colnames;
    if (frame)
    {
        nc = LENGTH(data);
        nr = nc ? LENGTH(VECTOR_ELT(data, 0)) : 0;
        rownames = getAttrib(data, R_RowNamesSymbol);
        colnames = getAttrib(data, R_NamesSymbol);
    }
    else
    {
        if (TYPEOF(dims) != INTSXP || LENGTH(dims) != 2)
            error("dta must be a numeric matrix or data frame");
        nr = INTEGER(dims)[0];
        nc = INTEGER(dims)[1];
        SEXP dn = getAttrib(data, R_DimNamesSymbol);
        if (TYPEOF(dn) != VECSXP || LENGTH(dn) != 2)
            error("dta must have row and column names");
        rownames = VECTOR_ELT(dn, 0);
        colnames = VECTOR_ELT(dn, 1);
    }
    if (nr < 1 || nc < 1 || TYPEOF(colnames) != STRSXP ||
        XLENGTH(colnames) != nc)
        error("dta must have marker rows and named population columns");
    int *no_all = st_work_alloc(nr, sizeof(int));
    int *no_ind = st_work_alloc(nc, sizeof(int));
    const char **marker = st_work_alloc(nr, sizeof(char *));
    const char **population = st_work_alloc(nc, sizeof(char *));
    int nmar = 0, npop = 0;
    for (int i = 0; i < nr; i++)
    {
        char number[32];
        snprintf(number, sizeof(number), "%d", i + 1);
        const char *name =
            row_labels ? row_labels[i] : TYPEOF(rownames) == STRSXP ? CHAR(STRING_ELT(rownames, i)) : number;
        const char *prefix = st_legacy_prefix(name);
        if (!i || strcmp(prefix, marker[nmar - 1]))
        {
            marker[nmar] = prefix;
            no_all[nmar++] = 0;
        }
        no_all[nmar - 1]++;
    }
    for (int i = 0; i < nc; i++)
    {
        const char *prefix = st_legacy_prefix(CHAR(STRING_ELT(colnames, i)));
        if (!i || strcmp(prefix, population[npop - 1]))
        {
            population[npop] = prefix;
            no_ind[npop++] = 0;
        }
        no_ind[npop - 1]++;
    }
    if (missing != 2)
    {
        const char **ordered =
            st_work_alloc(nmar > npop ? nmar : npop, sizeof(*ordered));
        memcpy(ordered, marker, nmar * sizeof(*ordered));
        qsort(ordered, nmar, sizeof(*ordered), st_legacy_name_compare);
        for (int i = 1; i < nmar; i++)
            if (!strcmp(ordered[i - 1], ordered[i]))
                error("Rows for each marker must be consecutive");
        memcpy(ordered, population, npop * sizeof(*ordered));
        qsort(ordered, npop, sizeof(*ordered), st_legacy_name_compare);
        for (int i = 1; i < npop; i++)
            if (!strcmp(ordered[i - 1], ordered[i]))
                error("Columns for each population must be consecutive");
    }
    size_t cells = (size_t)nr * nc;
    double *values = missing < 0 ? NULL : st_work_alloc(cells, sizeof(double));
    for (int j = 0; missing >= 0 && j < nc; j++)
        for (int i = 0; i < nr; i++)
        {
            double v = frame ? st_legacy_number(VECTOR_ELT(data, j), i)
                             : st_legacy_number(data, i + (size_t)nr * j);
            values[i + (size_t)nr * j] = missing == 1 && ISNAN(v) ? -1 : v;
        }
    int *individuals = st_work_alloc(nc, sizeof(int));
    int *markers = st_work_alloc(nmar, sizeof(int));
    for (int p = 0, offset = 0; p < npop; p++)
        for (int i = 0; i < no_ind[p]; i++)
            individuals[offset++] = i + 1;
    for (int i = 0; i < nmar; i++)
        markers[i] = i + 1;
    st_legacy_data result = {nr,         nc,     npop,        nmar,     frame,
                             no_all,     no_ind, individuals, markers,  marker,
                             population, data,   rownames,    colnames, values, row_labels};
    return result;
}

static st_legacy_data st_legacy_data_read(SEXP data, int missing)
{
    return st_legacy_data_read_labels(data, missing, NULL);
}

static SEXP st_legacy_pcoa(SEXP args)
{
    st_legacy_data d = st_legacy_data_read(st_legacy_argument(args, "data"), 1);
    int nr = d.nr, npop = d.npop, nmar = d.nmar;
    int *no_all = d.no_all, *no_ind = d.no_ind;
    int *individuals = d.individuals, *markers = d.markers;
    const char **population = d.population;
    double *values = d.values;
    int k = asInteger(st_legacy_argument(args, "k"));
    if (k < 1 || k >= npop)
        error("'k' must be in {1, 2, ..  n - 1}");
    size_t pairs = (size_t)npop * (npop - 1) / 2;
    double *freq = st_work_alloc((size_t)nr * npop, sizeof(double));
    double *distance = st_work_alloc(pairs, sizeof(double));
    allele_freq(&npop, no_ind, &nmar, no_all, no_ind, individuals, values,
                freq);
    mrd(&npop, &nmar, no_all, &nmar, markers, freq, distance);
    double *matrix = st_work_alloc((size_t)npop * npop, sizeof(double));
    double *means = st_work_alloc(npop, sizeof(double));
    memset(matrix, 0, (size_t)npop * npop * sizeof(double));
    for (int i = 0, at = 0; i < npop - 1; i++)
        for (int j = i + 1; j < npop; j++, at++)
        {
            if (!R_FINITE(distance[at]))
                error("NA/NaN/Inf in principal-coordinate distance matrix");
            matrix[i + (size_t)npop * j] = matrix[j + (size_t)npop * i] =
                distance[at] * distance[at];
        }
    for (int i = 0; i < npop; i++)
    {
        double sum = 0;
        for (int j = 0; j < npop; j++)
            sum += matrix[i + (size_t)npop * j];
        means[i] = sum / npop;
        for (int j = 0; j < npop; j++)
            matrix[i + (size_t)npop * j] -= means[i];
    }
    for (int j = 0; j < npop; j++)
    {
        double sum = 0;
        for (int i = 0; i < npop; i++)
            sum += matrix[i + (size_t)npop * j];
        means[j] = sum / npop;
        for (int i = 0; i < npop; i++)
            matrix[i + (size_t)npop * j] =
                -(matrix[i + (size_t)npop * j] - means[j]) / 2;
    }
    double *eigenvalues = st_work_alloc(npop, sizeof(double)), query;
    double *eigenvectors = st_work_alloc((size_t)npop * npop, sizeof(double));
    int *support = st_work_alloc((size_t)npop * 2, sizeof(int));
    int lwork = -1, liwork = -1, iquery, info = 0, found = 0;
    int first = 1, last = npop;
    double lower = 0, upper_bound = 0, abstol = 0;
    char vectors = 'V', range = 'A', upper = 'L';
    F77_CALL(dsyevr)(&vectors, &range, &upper, &npop, matrix, &npop, &lower,
                     &upper_bound, &first, &last, &abstol, &found, eigenvalues,
                     eigenvectors, &npop, support, &query, &lwork, &iquery,
                     &liwork, &info FCONE FCONE FCONE);
    if (info || !R_FINITE(query) || query > INT_MAX || iquery < 1)
        error("Cannot size principal-coordinate eigensolver workspace");
    lwork = (int)query;
    liwork = iquery;
    double *work = st_work_alloc(lwork, sizeof(double));
    int *iwork = st_work_alloc(liwork, sizeof(int));
    F77_CALL(dsyevr)(&vectors, &range, &upper, &npop, matrix, &npop, &lower,
                     &upper_bound, &first, &last, &abstol, &found, eigenvalues,
                     eigenvectors, &npop, support, work, &lwork, iwork, &liwork,
                     &info FCONE FCONE FCONE);
    if (info || found != npop)
        error("Principal-coordinate eigensolver failed");
    int kept = 0;
    for (int j = 0; j < k; j++)
        if (eigenvalues[npop - j - 1] > 0)
            kept++;
    if (kept < k)
        warning("only %d of the first %d eigenvalues are > 0", kept, k);
    SEXP out = PROTECT(allocMatrix(REALSXP, npop, kept));
    for (int j = 0, target = 0; j < k; j++)
    {
        int source = npop - j - 1;
        if (eigenvalues[source] <= 0)
            continue;
        double scale = sqrt(eigenvalues[source]);
        for (int i = 0; i < npop; i++)
            REAL(out)
            [i + (size_t)npop * target] =
                eigenvectors[i + (size_t)npop * source] * scale;
        target++;
    }
    SEXP labels = PROTECT(allocMatrix(STRSXP, npop, 1));
    SEXP dn = PROTECT(allocVector(VECSXP, 2));
    for (int i = 0; i < npop; i++)
        SET_STRING_ELT(labels, i, mkChar(population[i]));
    SET_VECTOR_ELT(dn, 0, labels);
    setAttrib(out, R_DimNamesSymbol, dn);
    UNPROTECT(3);
    return out;
}

static SEXP st_legacy_rng_init(void)
{
    int count = MAXTHREADS + 1, status = -2;
    int *seeds = st_work_alloc(count, sizeof(int));
    GetRNGstate();
    for (int i = 0; i < count; i++)
        seeds[i] = (int)floor(runif(0, 1) * 1073741824.0) + (i != 0);
    PutRNGstate();
    rng_init_seed(seeds, &count, &status);
    int np = 0;
    const char *keys[] = {"", "", "retval"};
    SEXP out = st_sim_list(3, keys, &np);
    SEXP seed_out = PROTECT(allocVector(INTSXP, count));
    np++;
    memcpy(INTEGER(seed_out), seeds, count * sizeof(int));
    SET_VECTOR_ELT(out, 0, seed_out);
    SET_VECTOR_ELT(out, 1, ScalarInteger(count));
    SET_VECTOR_ELT(out, 2, ScalarInteger(status));
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_plant(SEXP args)
{
    SEXP genotype = st_legacy_argument(args, "genotype");
    int constant = asLogical(st_legacy_argument(args, "constant"));
    double number = asReal(st_legacy_argument(args, "n"));
    if (!R_FINITE(number) || number < 1 || number != floor(number) ||
        number > INT_MAX || NoHoms < 1 || NoChroms < 1 ||
        number * NoHoms * NoChroms > INT_MAX)
        error("NoInd and genome dimensions must be positive integers in range");
    if (!isNumeric(genotype) || (constant && XLENGTH(genotype) != 1))
        error("Genotypes must contain numeric allele codes");
    int n = (int)number, rows = n * NoHoms * NoChroms, np = 0;
    const char *keys[] = {"ind", "chrom", "hom", "pos", "all"};
    SEXP out = st_sim_list(5, keys, &np);
    for (int j = 0; j < 5; j++)
        SET_VECTOR_ELT(out, j, allocVector(REALSXP, rows));
    for (int i = 0, at = 0; i < n; i++)
        for (int c = 0; c < NoChroms; c++)
            for (int h = 0; h < NoHoms; h++, at++)
            {
                REAL(VECTOR_ELT(out, 0))[at] = i + 1;
                REAL(VECTOR_ELT(out, 1))[at] = c + 1;
                REAL(VECTOR_ELT(out, 2))[at] = h + 1;
                REAL(VECTOR_ELT(out, 3))[at] = 0;
                REAL(VECTOR_ELT(out, 4))
                [at] = constant                ? asReal(genotype)
                       : h < XLENGTH(genotype) ? st_legacy_number(genotype, h)
                                               : NA_REAL;
            }
    st_sim_frame(out, rows, &np);
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_sdev(SEXP args)
{
    SEXP data = st_legacy_argument(args, "data");
    int remove_missing = asLogical(st_legacy_argument(args, "na.rm"));
    int frame = inherits(data, "data.frame"), matrix = isMatrix(data);
    SEXP dims = getAttrib(data, R_DimSymbol), names = R_NilValue;
    int columns = frame ? LENGTH(data) : matrix ? INTEGER(dims)[1] : 1;
    if (frame)
        names = getAttrib(data, R_NamesSymbol);
    else if (matrix)
    {
        SEXP dn = getAttrib(data, R_DimNamesSymbol);
        if (!isNull(dn))
            names = VECTOR_ELT(dn, 1);
    }
    SEXP out = PROTECT(allocVector(REALSXP, columns));
    for (int j = 0; j < columns; j++)
    {
        SEXP input = frame ? VECTOR_ELT(data, j) : data;
        if (!isNumeric(input) && TYPEOF(input) != LGLSXP)
            error("x must be numeric");
        R_xlen_t n = matrix ? INTEGER(dims)[0] : XLENGTH(input);
        size_t offset = matrix ? (size_t)j * n : 0, count = 0;
        int missing = 0;
        long double sum = 0;
        for (R_xlen_t i = 0; i < n; i++)
        {
            double x = st_legacy_number(input, offset + i);
            if (ISNAN(x))
                missing = 1;
            else
                sum += x, count++;
        }
        if ((missing && !remove_missing) || count < 2)
        {
            REAL(out)[j] = NA_REAL;
            continue;
        }
        double mean = (double)(sum / count);
        long double correction = 0;
        for (R_xlen_t i = 0; i < n; i++)
        {
            double x = st_legacy_number(input, offset + i);
            if (!ISNAN(x))
                correction += x - mean;
        }
        mean += (double)(correction / count);
        long double squares = 0;
        for (R_xlen_t i = 0; i < n; i++)
        {
            double x = st_legacy_number(input, offset + i);
            if (!ISNAN(x))
            {
                long double difference = x - mean;
                squares += difference * difference;
            }
        }
        REAL(out)[j] = sqrt((double)(squares / (count - 1)));
    }
    if (!isNull(names))
        setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(1);
    return out;
}

static const char *st_legacy_data_rowname(const st_legacy_data *d, int i)
{
    if (d->row_labels) return d->row_labels[i];
    if (TYPEOF(d->rownames) == STRSXP)
        return CHAR(STRING_ELT(d->rownames, i));
    char number[32];
    snprintf(number, sizeof(number), "%d", i + 1);
    return st_legacy_text(number);
}

static SEXP st_legacy_data_parameters(const st_legacy_data *d)
{
    const char *keys[] = {"colhead",  "rowhead",  "no.mar",
                          "no.all",   "no.pop",   "no.ind",
                          "ind.list", "mar.list", "pop.list"};
    int np = 0;
    SEXP out = st_sim_list(9, keys, &np);
    for (int type = 0; type < 2; type++)
    {
        int n = type ? d->nr : d->nc;
        SEXP head = PROTECT(allocMatrix(STRSXP, n, 2));
        np++;
        for (int i = 0; i < n; i++)
        {
            const char *name = type ? st_legacy_data_rowname(d, i)
                                    : CHAR(STRING_ELT(d->colnames, i));
            const char *split = strchr(name, '.');
            if (!split)
                split = strchr(name, '_');
            SET_STRING_ELT(
                head, i, split ? mkCharLen(name, split - name) : mkChar(name));
            SET_STRING_ELT(head, i + n, mkChar(split ? split + 1 : ""));
        }
        SET_VECTOR_ELT(out, type, head);
    }
    SET_VECTOR_ELT(out, 2, ScalarInteger(d->nmar));
    SET_VECTOR_ELT(out, 4, ScalarInteger(d->npop));
    SEXP alleles = PROTECT(allocMatrix(REALSXP, 1, d->nmar));
    SEXP individuals = PROTECT(allocMatrix(REALSXP, 1, d->npop));
    SEXP index = PROTECT(allocVector(INTSXP, d->nc));
    SEXP markers = PROTECT(allocMatrix(STRSXP, d->nmar, 1));
    SEXP populations = PROTECT(allocMatrix(STRSXP, d->npop, 1));
    np += 5;
    for (int i = 0; i < d->nmar; i++)
    {
        REAL(alleles)[i] = d->no_all[i];
        SET_STRING_ELT(markers, i, mkChar(d->marker[i]));
    }
    for (int i = 0; i < d->npop; i++)
    {
        REAL(individuals)[i] = d->no_ind[i];
        SET_STRING_ELT(populations, i, mkChar(d->population[i]));
    }
    memcpy(INTEGER(index), d->individuals, d->nc * sizeof(int));
    SET_VECTOR_ELT(out, 3, alleles);
    SET_VECTOR_ELT(out, 5, individuals);
    SET_VECTOR_ELT(out, 6, index);
    SET_VECTOR_ELT(out, 7, markers);
    SET_VECTOR_ELT(out, 8, populations);
    UNPROTECT(np);
    return out;
}

static int st_legacy_default(SEXP x)
{
    return isNull(x) || (XLENGTH(x) == 1 && asReal(x) == 0);
}

static int *st_legacy_indices(SEXP x, const int *fallback, int n, int upper,
                              const char *name)
{
    int *out = st_work_alloc(n ? n : 1, sizeof(int));
    if (st_legacy_default(x))
    {
        memcpy(out, fallback, n * sizeof(int));
        return out;
    }
    if (!isNumeric(x) || XLENGTH(x) != n)
        error("Incorrect length of %s", name);
    for (int i = 0; i < n; i++)
    {
        double v = st_legacy_number(x, i);
        if (!R_FINITE(v) || v < 1 || v != floor(v) || v > upper)
            error("Invalid %s index", name);
        out[i] = (int)v;
    }
    return out;
}

static char *st_legacy_join(const char **strings, int n)
{
    size_t length = 1;
    for (int i = 0; i < n; i++)
    {
        size_t size = strlen(strings[i]);
        if (size > SIZE_MAX - length - 1)
            error("Joined identifier list is too long");
        length += size + 1;
    }
    char *out = st_work_alloc(length, 1), *p = out;
    for (int i = 0; i < n; i++)
    {
        if (i)
            *p++ = ' ';
        size_t size = strlen(strings[i]);
        memcpy(p, strings[i], size);
        p += size;
    }
    *p = '\0';
    return out;
}

static SEXP st_legacy_generate(st_legacy_data *d, SEXP args)
{
    int *populations = st_work_alloc(d->npop, sizeof(int));
    int *alleles = st_work_alloc(d->nr, sizeof(int));
    for (int i = 0; i < d->npop; i++)
        populations[i] = i + 1;
    for (int i = 0; i < d->nr; i++)
    {
        const char *name = st_legacy_data_rowname(d, i),
                   *suffix = strchr(name, '.');
        if (!suffix)
            suffix = strchr(name, '_');
        char *end = NULL;
        double allele = suffix ? strtod(suffix + 1, &end) : NA_REAL;
        if (!suffix || end == suffix + 1 || *end || !R_FINITE(allele) ||
            allele == 0 || allele < INT_MIN || allele > INT_MAX)
        {
            info(
                -2,
                "Error in allele encoding: only objects of type integer allowed\n");
            return R_NilValue;
        }
        alleles[i] = (int)allele;
    }
    int background = asInteger(st_legacy_argument(args, "background"));
    int backcross = asLogical(st_legacy_argument(args, "backcross"));
    char *pop_list = st_legacy_join(d->population, d->npop);
    char *mar_list = st_legacy_join(d->marker, d->nmar);
    size_t nf = (size_t)d->nr * d->npop;
    double *freq = st_work_alloc(nf, sizeof(double));
    for (size_t i = 0; i < nf; i++)
        freq[i] = 16;
    generate_population(&d->npop, d->no_ind, &d->nmar, d->no_all, &d->npop,
                        populations, d->no_ind, d->individuals, &d->nmar,
                        d->markers, d->values, freq, &pop_list, &mar_list,
                        alleles, &background, &backcross);
    const char *keys[] = {"", "",    "",     "", "", "", "", "", "",
                          "", "DTA", "FREQ", "", "", "", "", ""};
    int np = 0;
    SEXP out = st_sim_list(17, keys, &np);
    const int *arrays[] = {&d->npop, d->no_ind,   &d->nmar,  d->no_all,
                           &d->npop, populations, d->no_ind, d->individuals,
                           &d->nmar, d->markers};
    int lengths[] = {1,       d->npop, 1,     d->nmar, 1,
                     d->npop, d->npop, d->nc, 1,       d->nmar};
    for (int j = 0; j < 10; j++)
    {
        SEXP col = PROTECT(allocVector(INTSXP, lengths[j]));
        np++;
        memcpy(INTEGER(col), arrays[j], lengths[j] * sizeof(int));
        SET_VECTOR_ELT(out, j, col);
    }
    SEXP values = PROTECT(allocVector(REALSXP, (R_xlen_t)d->nr * d->nc));
    SEXP frequencies = PROTECT(allocVector(REALSXP, nf));
    SEXP allele_out = PROTECT(allocVector(INTSXP, d->nr));
    np += 3;
    memcpy(REAL(values), d->values, (size_t)d->nr * d->nc * sizeof(double));
    memcpy(REAL(frequencies), freq, nf * sizeof(double));
    memcpy(INTEGER(allele_out), alleles, d->nr * sizeof(int));
    SET_VECTOR_ELT(out, 10, values);
    SET_VECTOR_ELT(out, 11, frequencies);
    SET_VECTOR_ELT(out, 12, mkString(pop_list));
    SET_VECTOR_ELT(out, 13, mkString(mar_list));
    SET_VECTOR_ELT(out, 14, allele_out);
    SET_VECTOR_ELT(out, 15, ScalarInteger(background));
    SET_VECTOR_ELT(out, 16, ScalarInteger(backcross));
    UNPROTECT(np);
    return out;
}

static int st_legacy_marker_name_order(const void *left, const void *right)
{
    return strcoll(*(const char *const *)left, *(const char *const *)right);
}

static SEXP st_legacy_data_diagnostics(st_legacy_data *d, int irregular)
{
    if (irregular)
    {
        double unused = 0;
        correct_missing(&d->npop, d->no_ind, &d->nmar, d->no_all, d->no_ind,
                        d->individuals, d->values, &unused);
    }
    size_t bound = (size_t)d->nmar * d->nc;
    if (bound > INT_MAX)
        error("Genetic-data diagnostic table is too large");
    const char **names = st_work_alloc(bound ? bound : 1, sizeof(*names));
    int *individual = st_work_alloc(bound ? bound : 1, sizeof(*individual));
    int rows = 0;
    for (int j = 0; j < d->nc; j++)
    {
        int first = rows, offset = 0;
        for (int m = 0; m < d->nmar; m++)
        {
            int flagged = 0;
            for (int a = 0; a < d->no_all[m]; a++, offset++)
            {
                double value = d->values[offset + (size_t)d->nr * j];
                if (irregular ? !ISNAN(value) && value != -1 && value != 0 &&
                                    value != 1 && value != .5
                              : ISNAN(value))
                    flagged = 1;
            }
            if (flagged)
                names[rows++] = d->marker[m];
        }
        qsort(names + first, rows - first, sizeof(*names),
              st_legacy_marker_name_order);
        int unique = first;
        for (int i = first; i < rows; i++)
            if (unique == first || strcmp(names[i], names[unique - 1]))
                names[unique++] = names[i];
        rows = unique;
        for (int i = first; i < rows; i++)
            individual[i] = j;
    }
    int np = 0;
    const char *keys[] = {"ind", "marker"};
    SEXP out = st_sim_list(2, keys, &np);
    SEXP ids = PROTECT(allocVector(STRSXP, rows));
    SEXP markers = PROTECT(allocVector(STRSXP, rows));
    SEXP rn = PROTECT(allocVector(INTSXP, rows));
    np += 3;
    for (int i = 0; i < rows; i++)
    {
        SET_STRING_ELT(ids, i, STRING_ELT(d->colnames, individual[i]));
        SET_STRING_ELT(markers, i, mkChar(names[i]));
        INTEGER(rn)[i] = i + 2;
    }
    SET_VECTOR_ELT(out, 0, ids);
    SET_VECTOR_ELT(out, 1, markers);
    st_sim_frame(out, rows, &np);
    setAttrib(out, R_RowNamesSymbol, rn);
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_distance_matrix(const st_legacy_data *d, SEXP values)
{
    int n = d->npop, np = 0;
    R_xlen_t pairs = (R_xlen_t)n * (n - 1) / 2, count = XLENGTH(values);
    if (pairs && TYPEOF(values) == RAWSXP)
        error(
            "incompatible types (from double to raw) in subassignment type fix");
    if (pairs && !count)
        error("replacement has length zero");
    if (pairs && pairs % count)
    {
        warning(
            "number of items to replace is not a multiple of replacement length");
        warning(
            "number of items to replace is not a multiple of replacement length");
    }
    if (pairs && TYPEOF(values) == VECSXP)
        error("'dimnames' applied to non-array");
    SEXPTYPE type = TYPEOF(values) == STRSXP    ? STRSXP
                    : TYPEOF(values) == VECSXP  ? VECSXP
                    : TYPEOF(values) == CPLXSXP ? CPLXSXP
                                                : REALSXP;
    if (!pairs)
        type = REALSXP;
    SEXP out = PROTECT(allocMatrix(type, n, n));
    SEXP labels = PROTECT(allocMatrix(STRSXP, n, 1));
    SEXP names = PROTECT(allocVector(VECSXP, 2));
    np += 3;
    SEXP zero = R_NilValue;
    if (type == VECSXP)
    {
        zero = PROTECT(ScalarReal(0));
        np++;
    }
    for (int i = 0; i < n; i++)
    {
        SET_STRING_ELT(labels, i, mkChar(d->population[i]));
        R_xlen_t at = i + (R_xlen_t)n * i;
        if (type == STRSXP)
            SET_STRING_ELT(out, at, mkChar("0"));
        else if (type == VECSXP)
            SET_VECTOR_ELT(out, at, zero);
        else if (type == CPLXSXP)
        {
            COMPLEX(out)[at].r = 0;
            COMPLEX(out)[at].i = 0;
        }
        else
            REAL(out)[at] = 0;
    }
    R_xlen_t at = 0;
    for (int j = 0; j < n - 1; j++)
        for (int i = j + 1; i < n; i++, at++)
        {
            R_xlen_t source = at % count, lower = i + (R_xlen_t)n * j,
                     upper = j + (R_xlen_t)n * i;
            if (type == STRSXP)
            {
                SET_STRING_ELT(out, lower, STRING_ELT(values, source));
                SET_STRING_ELT(out, upper, STRING_ELT(values, source));
            }
            else if (type == VECSXP)
            {
                SET_VECTOR_ELT(out, lower, VECTOR_ELT(values, source));
                SET_VECTOR_ELT(out, upper, VECTOR_ELT(values, source));
            }
            else if (type == CPLXSXP)
                COMPLEX(out)
                [lower] = COMPLEX(out)[upper] = COMPLEX(values)[source];
            else
                REAL(out)
                [lower] = REAL(out)[upper] =
                    TYPEOF(values) == RAWSXP ? RAW(values)[source]
                                             : st_legacy_number(values, source);
        }
    SET_VECTOR_ELT(names, 0, labels);
    SET_VECTOR_ELT(names, 1, labels);
    setAttrib(out, R_DimNamesSymbol, names);
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_distance_result(const st_legacy_data *d,
                                      const double *distance,
                                      const double *deviation, int sd,
                                      int matrix)
{
    size_t pairs = (size_t)d->npop * (d->npop - 1) / 2;
    if (matrix)
    {
        SEXP out = PROTECT(allocVector(REALSXP, pairs));
        SEXP labels = PROTECT(allocMatrix(STRSXP, d->npop, 1));
        memcpy(REAL(out), distance, pairs * sizeof(*distance));
        for (int i = 0; i < d->npop; i++)
            SET_STRING_ELT(labels, i, mkChar(d->population[i]));
        setAttrib(out, install("Size"), ScalarInteger(d->npop));
        setAttrib(out, install("Labels"), labels);
        setAttrib(out, install("Diag"), ScalarLogical(FALSE));
        setAttrib(out, install("Upper"), ScalarLogical(TRUE));
        setAttrib(out, install("method"), mkString("euclidean"));
        setAttrib(out, R_ClassSymbol, mkString("dist"));
        UNPROTECT(2);
        return out;
    }
    const char *keys[] = {"OTU1", "OTU2", "Measure", "SDev"};
    int np = 0;
    SEXP out = st_sim_list(sd ? 4 : 3, keys, &np);
    for (int j = 0; j < (sd ? 4 : 3); j++)
        SET_VECTOR_ELT(out, j, allocVector(j < 2 ? STRSXP : REALSXP, pairs));
    for (int i = 0, at = 0; i < d->npop - 1; i++)
        for (int j = i + 1; j < d->npop; j++, at++)
        {
            SET_STRING_ELT(VECTOR_ELT(out, 0), at, mkChar(d->population[i]));
            SET_STRING_ELT(VECTOR_ELT(out, 1), at, mkChar(d->population[j]));
            REAL(VECTOR_ELT(out, 2))[at] = distance[at];
            if (sd)
                REAL(VECTOR_ELT(out, 3))[at] = deviation[at];
        }
    st_sim_frame(out, (int)pairs, &np);
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_distance_values(st_legacy_data *d, int measure, int sd,
                                      const char *parameters, int *no_ind,
                                      int *ind, int nmar, int *mar, int matrix,
                                      int auxiliary)
{
    int *populations = st_work_alloc(d->npop, sizeof(int));
    for (int i = 0; i < d->npop; i++)
        populations[i] = i + 1;
    size_t pairs = (size_t)d->npop * (d->npop - 1) / 2;
    if (pairs > INT_MAX)
        error("Too many pairwise genetic distances");
    double *freq = st_work_alloc((size_t)d->nr * d->npop, sizeof(double));
    double *distance = st_work_alloc(pairs, sizeof(double));
    double *deviation = st_work_alloc(pairs, sizeof(double));
    for (size_t i = 0; i < pairs; i++)
        distance[i] = deviation[i] = i + 1;
    char *par = st_legacy_text(parameters);
    gen_dist(&measure, &d->npop, d->no_ind, &d->nmar, d->no_all, &d->npop,
             populations, no_ind, ind, &nmar, mar, d->values, freq, distance,
             &sd, deviation, &par);
    SEXP result =
        PROTECT(st_legacy_distance_result(d, distance, deviation, sd, matrix));
    if (!auxiliary)
    {
        UNPROTECT(1);
        return result;
    }
    SEXP table =
        matrix
            ? PROTECT(st_legacy_distance_result(d, distance, deviation, sd, 0))
            : PROTECT(result);
    int np = 2;
    const char *keys[] = {"result", "table"};
    SEXP out = st_sim_list(2, keys, &np);
    SET_VECTOR_ELT(out, 0, result);
    SET_VECTOR_ELT(out, 1, table);
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_data_operation(SEXP args)
{
    const char *operation = CHAR(asChar(st_legacy_argument(args, "operation")));
    int generate = !strcmp(operation, "generate");
    SEXP source = st_legacy_argument(args, "data");
    const char **row_labels = asLogical(st_legacy_argument(args, "primers")) == TRUE
        ? st_legacy_primer_labels(source) : NULL;
    st_legacy_data d = st_legacy_data_read_labels(
        source,
        !strcmp(operation, "parameters") || !strcmp(operation, "matrix") ? -1
        : !strcmp(operation, "missing")                                  ? 2
                                        : !generate, row_labels);
    if (!strcmp(operation, "matrix"))
        return st_legacy_distance_matrix(&d,
                                         st_legacy_argument(args, "distance"));
    if (!strcmp(operation, "missing") || !strcmp(operation, "irregular"))
        return st_legacy_data_diagnostics(&d, !strcmp(operation, "irregular"));
    if (generate)
        return st_legacy_generate(&d, args);
    if (!strcmp(operation, "parameters"))
        return st_legacy_data_parameters(&d);
    int frequencies = !strcmp(operation, "frequencies");
    int correct = !strcmp(operation, "correct");
    if (frequencies || correct)
    {
        double *freq = st_work_alloc((size_t)d.nr * d.npop, sizeof(double));
        if (correct)
            correct_missing(&d.npop, d.no_ind, &d.nmar, d.no_all, d.no_ind,
                            d.individuals, d.values, freq);
        else
            allele_freq(&d.npop, d.no_ind, &d.nmar, d.no_all, d.no_ind,
                        d.individuals, d.values, freq);
        int columns = correct ? d.nc : d.npop, np = 0;
        SEXP out = PROTECT(allocVector(VECSXP, columns));
        SEXP names = PROTECT(allocVector(STRSXP, columns));
        SEXP rows = PROTECT(allocVector(STRSXP, d.nr));
        np += 3;
        for (int i = 0; i < d.nr; i++)
            SET_STRING_ELT(rows, i, mkChar(st_legacy_data_rowname(&d, i)));
        for (int j = 0; j < columns; j++)
        {
            SEXP col = PROTECT(allocVector(REALSXP, d.nr));
            np++;
            for (int i = 0; i < d.nr; i++)
            {
                double value = correct ? d.values[i + (size_t)d.nr * j]
                                       : freq[i + (size_t)d.nr * j];
                REAL(col)
                [i] = value == -1              ? NA_REAL
                      : correct && value == .5 ? 1
                                               : value;
            }
            SET_VECTOR_ELT(out, j, col);
            SET_STRING_ELT(names, j,
                           correct ? STRING_ELT(d.colnames, j)
                                   : mkChar(d.population[j]));
        }
        st_sim_frame(out, d.nr, &np);
        setAttrib(out, R_RowNamesSymbol, rows);
        if (!correct || d.frame)
            setAttrib(out, R_NamesSymbol, names);
        UNPROTECT(np);
        return out;
    }
    if (strcmp(operation, "distance"))
        error("Unknown genetic-data operation");
    const char *measure_name =
        CHAR(asChar(st_legacy_argument(args, "measure")));
    const char *measure_names[] = {"rd", "mrd", "euc", "dic", "jac", "sma"};
    int codes[] = {1, 2, 3, 11, 12, 13}, measure = 0;
    for (int i = 0; i < 6; i++)
        if (!strcmp(measure_name, measure_names[i]))
            measure = codes[i];
    if (!measure || d.npop < 2)
        error("A supported distance measure and two populations are required");
    const char *parameters = CHAR(asChar(st_legacy_argument(args, "par")));
    int sd = !strncmp(parameters, "sdev/jack", 9);
    if (!strncmp(parameters, "sdev/boot", 9))
        sd = asInteger(st_legacy_argument(args, "boot"));
    if (sd < 0)
        error("boot must be nonnegative");
    int *no_ind = st_legacy_indices(st_legacy_argument(args, "no.ind.u"),
                                    d.no_ind, d.npop, d.nc, "no.ind.u");
    int *ind = st_legacy_indices(st_legacy_argument(args, "ind.u"),
                                 d.individuals, d.nc, d.nc, "ind.u");
    for (int p = 0, offset = 0; p < d.npop; offset += d.no_ind[p++])
    {
        if (no_ind[p] > d.no_ind[p])
            error("Too many individuals selected from population");
        for (int i = 0; i < no_ind[p]; i++)
            if (ind[offset + i] > d.no_ind[p])
                error("Individual index exceeds population size");
    }
    SEXP mar_count = st_legacy_argument(args, "no.mar.u");
    int nmar = st_legacy_default(mar_count) ? d.nmar : asInteger(mar_count);
    if (nmar < 1 || nmar > d.nmar)
        error("Invalid number of selected markers");
    int *mar = st_legacy_indices(st_legacy_argument(args, "mar.u"), d.markers,
                                 nmar, d.nmar, "mar.u");
    return st_legacy_distance_values(&d, measure, sd, parameters, no_ind, ind,
                                     nmar, mar, 0, 0);
}

extern void st_legacy_marker_presence(const char *, int *, int *, int *, int **,
                                      const char ***, const char ***,
                                      double **);

static const char *st_legacy_population_label(const char *name)
{
    size_t length = strlen(name), used = 0;
    char *out = st_work_alloc(length + 2, 1);
    mbstate_t state = {0};
    const char *cursor = name;
    int first = 1;
    while (*cursor)
    {
        wchar_t character;
        size_t width = mbrtowc(&character, cursor, strlen(cursor), &state);
        if (width == (size_t)-1 || width == (size_t)-2)
            error("invalid multibyte string in NTSys individual name");
        if (first && !iswalpha(character))
            out[used++] = 'X';
        first = 0;
        if (!(iswalnum(character) || character == '_'))
            break;
        memcpy(out + used, cursor, width);
        used += width;
        cursor += width;
    }
    out[used] = '\0';
    return out;
}

static SEXP st_legacy_marker_distance(SEXP args)
{
    st_legacy_data d = {0};
    const char **individual_names;
    st_legacy_marker_presence(
        CHAR(asChar(st_legacy_argument(args, "data.set"))), &d.nr, &d.nc,
        &d.nmar, &d.no_all, &d.marker, &individual_names, &d.values);
    if (!d.nr || !d.nc)
        error("dta must have marker rows and named population columns");
    d.no_ind = st_work_alloc(d.nc, sizeof(*d.no_ind));
    d.individuals = st_work_alloc(d.nc, sizeof(*d.individuals));
    d.population = st_work_alloc(d.nc, sizeof(*d.population));
    d.markers = st_work_alloc(d.nmar, sizeof(*d.markers));
    for (int i = 0; i < d.nc; i++)
    {
        const char *name = st_legacy_population_label(individual_names[i]);
        if (!i || strcmp(name, d.population[d.npop - 1]))
        {
            d.population[d.npop] = name;
            d.no_ind[d.npop++] = 0;
        }
        d.individuals[i] = ++d.no_ind[d.npop - 1];
    }
    const char **ordered = st_work_alloc(d.npop, sizeof(*ordered));
    memcpy(ordered, d.population, d.npop * sizeof(*ordered));
    qsort(ordered, d.npop, sizeof(*ordered), st_legacy_name_compare);
    for (int i = 1; i < d.npop; i++)
        if (!strcmp(ordered[i - 1], ordered[i]))
            error("Columns for each population must be consecutive");
    for (int i = 0; i < d.nmar; i++)
        d.markers[i] = i + 1;

    SEXP measure_arg = st_legacy_argument(args, "measure"),
         format_arg = st_legacy_argument(args, "format");
    if (isNull(measure_arg) || !XLENGTH(measure_arg))
        error("missing value where TRUE/FALSE needed");
    if (XLENGTH(measure_arg) > 1)
        error("'length = %lld' in coercion to 'logical(1)'",
              (long long)XLENGTH(measure_arg));
    SEXP measure_text = asChar(measure_arg);
    if (measure_text == NA_STRING)
        error("missing value where TRUE/FALSE needed");
    const char *name = CHAR(measure_text);
    int measure = !strcmp(name, "rd")    ? 1
                  : !strcmp(name, "mrd") ? 2
                  : !strcmp(name, "euc") ? 3
                                         : 0;
    if (!measure)
        info(-2, "Distance measure not defined\n");
    else
    {
        int allow;
        get_allow_zeros(&allow);
        if (d.nr == d.nmar)
        {
            if (!allow)
            {
                info(-1, "Dominant marker data?\n");
                info(-1, "Consider 'gd.allow.zero.frequencies(1)'\n");
            }
            if (measure != 3)
                info(-1, "Consider the Euclidean distance for dominant data\n");
            if (d.nc > d.npop)
                info(-1, "analysis on population level for dominant data\n");
        }
        if (d.nr > d.nmar && allow == 1)
        {
            info(-1, "Co-dominant marker data?\n");
            info(-1, "Consider 'gd.allow.zero.frequencies(0)'\n");
        }
        if (d.npop < 2)
            error(
                "A supported distance measure and two populations are required");
    }
    if (isNull(format_arg) || !XLENGTH(format_arg))
        error("argument is of length zero");
    if (XLENGTH(format_arg) > 1)
        error("the condition has length > 1");
    SEXP format_text = asChar(format_arg);
    if (format_text == NA_STRING)
        error("missing value where TRUE/FALSE needed");
    int matrix = !strcmp(CHAR(format_text), "m");
    int auxiliary = asLogical(st_legacy_argument(args, "auxiliary")) == TRUE;
    if (!measure)
    {
        if (matrix)
            error("attempt to set an attribute on NULL");
        if (!auxiliary)
            return R_NilValue;
        const char *keys[] = {"result", "table"};
        int np = 0;
        SEXP out = st_sim_list(2, keys, &np);
        UNPROTECT(np);
        return out;
    }
    return st_legacy_distance_values(&d, measure, 0, "dist", d.no_ind,
                                     d.individuals, d.nmar, d.markers, matrix,
                                     auxiliary);
}

static SEXP st_legacy_evaluate(SEXP args, int genotype);
static SEXP st_legacy_mab(SEXP args);

static SEXP st_legacy_genome_evaluate(SEXP args);
static SEXP st_legacy_select_segments(SEXP args);
static SEXP st_legacy_single_cross(SEXP args);
static SEXP st_legacy_cross_plan_execute(SEXP args);
static SEXP st_legacy_mdp(SEXP args);

static SEXP st_legacy_target_export(SEXP args);
static SEXP st_legacy_select_genotype(SEXP args);
static SEXP st_legacy_locus_effect(SEXP args);

extern void gs_set_performance_data_01_GV(int *, double *, int *, char **);

static SEXP st_legacy_performance_copy(SEXP args)
{
    SEXP population = st_legacy_argument(args, "population"),
         dataset = st_legacy_argument(args, "data.set");
    if (TYPEOF(population) != STRSXP || !XLENGTH(population) ||
        TYPEOF(dataset) != STRSXP || !XLENGTH(dataset))
        error("Invalid performance transfer arguments");
    char *name = (char *)CHAR(STRING_ELT(population, 0)),
         *data_set = (char *)CHAR(STRING_ELT(dataset, 0));
    TPopListEl *entry = NULL;
    find_list_element(name, &entry);
    if (!entry)
    {
        info(-2, "Population not defined: %s\n", name);
        error("invalid 'times' argument");
    }
    if (entry->Pop->NoInds > INT_MAX)
    {
        warning("NAs introduced by coercion to integer range");
        error("invalid 'times' argument");
    }
    int n = (int)entry->Pop->NoInds, status = -2;
    double *values = st_work_alloc(n, sizeof(*values));
    memset(values, 0, (size_t)n * sizeof(*values));
    char *effect = "", *criterion = "G";
    get_population_value(&name, values, &effect, &criterion, &status);
    gs_set_performance_data_01_GV(&n, values, &status, &data_set);
    return ScalarInteger(status);
}

typedef struct
{
    const char *name, *source, *column;
    double scale;
    int regional, quantile, marker, parameter;
} st_mab_report;

static R_xlen_t st_mab_report_length(SEXP value)
{
    return isNull(value) ? 0 : XLENGTH(value);
}

static int st_mab_report_flag(SEXP value)
{
    if (!st_mab_report_length(value))
        error("argument is of length zero");
    if (st_mab_report_length(value) > 1)
        error("the condition has length > 1");
    int flag = asLogical(value);
    if (flag == NA_LOGICAL)
    {
        if (TYPEOF(value) == STRSXP && STRING_ELT(value, 0) != NA_STRING)
            error("argument is not interpretable as logical");
        error("missing value where TRUE/FALSE needed");
    }
    return flag;
}

static const char *st_mab_report_real(double value)
{
    if (ISNA(value))
        return NULL;
    if (ISNAN(value))
        return "NaN";
    if (!R_FINITE(value))
        return value < 0 ? "-Inf" : "Inf";
    if (value == 0)
        value = 0;
    char scientific[80], output[768];
    snprintf(scientific, sizeof(scientific), "%.*e", DBL_DIG - 1, value);
    char *exponent = strchr(scientific, 'e'), *point = strchr(scientific, '.');
    int power = (int)strtol(exponent + 1, NULL, 10);
    char *last = exponent - 1;
    while (point && last > point && *last == '0')
        last--;
    if (point && last == point)
        last--;
    memmove(last + 1, exponent, strlen(exponent) + 1);
    exponent = last + 1;
    int sign = scientific[0] == '-', significant = 0;
    for (char *at = scientific + sign; at < exponent; at++)
        significant += *at >= '0' && *at <= '9';
    int decimals = significant - 1 - power;
    if (decimals < 0)
        decimals = 0;
    int fixed_width =
        sign + (power >= 0 ? power + 1 : 1) + decimals + (decimals > 0);
    SEXP option = GetOption1(install("scipen"));
    int penalty = isNull(option) ? 0 : asInteger(option);
    if (penalty == NA_INTEGER)
        penalty = 0;
    if (fixed_width <= (int)strlen(scientific) + penalty)
    {
        snprintf(output, sizeof(output), "%*.*f", fixed_width, decimals, value);
        char *dot = strchr(output, '.');
        if (dot)
        {
            char *end = output + strlen(output) - 1;
            while (end > dot && *end == '0')
                *end-- = '\0';
            if (end == dot)
                *end = '\0';
        }
    }
    else
        strcpy(output, scientific);
    const char *decimal = ".";
    option = GetOption1(install("OutDec"));
    if (TYPEOF(option) == STRSXP && LENGTH(option) == 1 &&
        STRING_ELT(option, 0) != NA_STRING)
        decimal = CHAR(STRING_ELT(option, 0));
    char *dot = strchr(output, '.');
    if (dot && strcmp(decimal, "."))
    {
        size_t left = (size_t)(dot - output), mark = strlen(decimal),
               right = strlen(dot + 1);
        char *result = st_work_alloc(left + mark + right + 1, 1);
        memcpy(result, output, left);
        memcpy(result + left, decimal, mark);
        memcpy(result + left + mark, dot + 1, right + 1);
        return result;
    }
    return st_legacy_text(output);
}

static const char *st_mab_report_text(SEXP value, R_xlen_t index, int factor)
{
    if (TYPEOF(value) == STRSXP)
        return STRING_ELT(value, index) == NA_STRING
                   ? NULL
                   : translateCharUTF8(STRING_ELT(value, index));
    if (TYPEOF(value) == INTSXP && factor && inherits(value, "factor"))
    {
        int code = INTEGER(value)[index];
        SEXP levels = getAttrib(value, R_LevelsSymbol);
        return code == NA_INTEGER
                   ? NULL
                   : translateCharUTF8(STRING_ELT(levels, code - 1));
    }
    if (TYPEOF(value) == LGLSXP)
        return LOGICAL(value)[index] == NA_LOGICAL ? NULL
               : LOGICAL(value)[index]             ? "TRUE"
                                                   : "FALSE";
    if (TYPEOF(value) == INTSXP)
    {
        if (INTEGER(value)[index] == NA_INTEGER)
            return NULL;
        char buffer[32];
        snprintf(buffer, sizeof(buffer), "%d", INTEGER(value)[index]);
        return st_legacy_text(buffer);
    }
    if (TYPEOF(value) == REALSXP)
        return st_mab_report_real(REAL(value)[index]);
    if (TYPEOF(value) == CPLXSXP)
    {
        Rcomplex z = COMPLEX(value)[index];
        if (ISNA(z.r) || ISNA(z.i))
            return NULL;
        const char *real = st_mab_report_real(z.r),
                   *imaginary = st_mab_report_real(fabs(z.i));
        size_t size = strlen(real) + strlen(imaginary) + 3;
        char *text = st_work_alloc(size, 1);
        snprintf(text, size, "%s%c%si", real, z.i < 0 ? '-' : '+', imaginary);
        return text;
    }
    if (TYPEOF(value) == RAWSXP)
    {
        char buffer[3];
        snprintf(buffer, sizeof(buffer), "%02x", RAW(value)[index]);
        return st_legacy_text(buffer);
    }
    error("Unsupported simulation report metadata");
    return NULL;
}

static SEXP st_mab_report_field(SEXP value, int collapse)
{
    if (!collapse)
    {
        if (st_mab_report_length(value) != 1)
            error("number of items to replace is not a multiple of "
                  "replacement length");
        const char *text = st_mab_report_text(value, 0, 0);
        return text ? mkCharCE(text, CE_UTF8) : NA_STRING;
    }
    R_xlen_t n = st_mab_report_length(value);
    const char **parts = st_work_alloc(n, sizeof(*parts));
    size_t size = n ? n : 1;
    for (R_xlen_t i = 0; i < n; i++)
    {
        parts[i] = st_mab_report_text(value, i, 1);
        if (!parts[i])
            parts[i] = "NA";
        if (strlen(parts[i]) > SIZE_MAX - size)
            error("Simulation report text is too long");
        size += strlen(parts[i]);
    }
    char *text = st_work_alloc(size, 1), *at = text;
    for (R_xlen_t i = 0; i < n; i++)
    {
        if (i)
            *at++ = '/';
        size_t length = strlen(parts[i]);
        memcpy(at, parts[i], length);
        at += length;
    }
    *at = '\0';
    return mkCharCE(text, CE_UTF8);
}

static const char *st_mab_report_name(const char *prefix, int index,
                                      const char *suffix)
{
    size_t size = strlen(prefix) + strlen(suffix) + 32;
    char *text = st_work_alloc(size, 1);
    snprintf(text, size, "%s%d%s", prefix, index, suffix);
    return text;
}

static SEXP st_mab_report_column(SEXP source, const char *key, R_xlen_t *offset,
                                 R_xlen_t *length)
{
    *offset = *length = 0;
    if (isNull(source))
        return R_NilValue;
    SEXP names = R_NilValue, dim = getAttrib(source, R_DimSymbol);
    int matrix = TYPEOF(dim) == INTSXP && LENGTH(dim) == 2;
    if (matrix)
    {
        SEXP dn = getAttrib(source, R_DimNamesSymbol);
        if (TYPEOF(dn) == VECSXP && LENGTH(dn) == 2)
            names = VECTOR_ELT(dn, 1);
    }
    else if (TYPEOF(source) == VECSXP)
        names = getAttrib(source, R_NamesSymbol);
    int found = -1, partial = -1, matches = 0;
    if (TYPEOF(names) == STRSXP)
        for (int i = 0; i < LENGTH(names); i++)
        {
            const char *name = CHAR(STRING_ELT(names, i));
            if (!strcmp(name, key))
            {
                found = i;
                break;
            }
            if (!strncmp(name, key, strlen(key)))
            {
                partial = i;
                matches++;
            }
        }
    if (found < 0 && matches == 1)
        found = partial;
    if (found < 0)
        return R_NilValue;
    if (matrix)
    {
        *length = INTEGER(dim)[0];
        *offset = *length * found;
        return source;
    }
    SEXP column = VECTOR_ELT(source, found);
    *length = XLENGTH(column);
    return column;
}

static void st_mab_report_dimnames(SEXP matrix, SEXP rows, SEXP columns)
{
    SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, rows);
    SET_VECTOR_ELT(dimnames, 1, columns);
    setAttrib(matrix, R_DimNamesSymbol, dimnames);
    UNPROTECT(1);
}

static SEXP st_mab_report_parameters(SEXP data, SEXP columns, int extended)
{
    const char *short_keys[] = {"repetitions",  "donor.parent",
                                "gen.type",     "population.size",
                                "sel.strategy", "no.selected"};
    const char *short_labels[] = {"Reps",     "DP",      "Pop type",
                                  "Pop size", "Sel str", "n Sel"};
    const char *long_keys[] = {
        "repetitions",  "recurrent.parent", "donor.parent",   "linkage.map",
        "target.loci",  "flanking.loci",    "gen.type",       "population.size",
        "sel.strategy", "no.selected",      "no.preselected", "success.factor",
        "reg.chr",      "reg.begin",        "reg.end"};
    const char *long_labels[] = {"Reps",    "RP",      "DP",       "Map",
                                 "Target",  "Flm",     "Pop type", "Pop size",
                                 "Sel str", "n Sel",   "n PreSel", "Succ f",
                                 "Reg chr", "Reg beg", "Reg end"};
    const char **keys = extended ? long_keys : short_keys,
               **labels = extended ? long_labels : short_labels;
    int n = LENGTH(data),
        rows = !extended ? 6
               : isNull(st_legacy_argument(VECTOR_ELT(data, 0), "reg.chr"))
                   ? 12
                   : 15;
    SEXP out = PROTECT(allocMatrix(STRSXP, rows, n)),
         rownames = PROTECT(allocVector(STRSXP, rows));
    for (int i = 0; i < rows; i++)
        SET_STRING_ELT(rownames, i, mkChar(labels[i]));
    for (int j = 0; j < n; j++)
        for (int i = 0; i < rows; i++)
            SET_STRING_ELT(out, i + (size_t)rows * j,
                           st_mab_report_field(
                               st_legacy_argument(VECTOR_ELT(data, j), keys[i]),
                               i >= (extended ? 4 : 2)));
    st_mab_report_dimnames(out, rownames, columns);
    UNPROTECT(2);
    return out;
}

static SEXP st_mab_report_numeric(SEXP data, SEXP columns, int rows,
                                  const st_mab_report *report, SEXP digits)
{
    int n = LENGTH(data);
    SEXP *sources = st_work_alloc(n, sizeof(*sources));
    R_xlen_t *offset = st_work_alloc(n, sizeof(*offset)),
             *length = st_work_alloc(n, sizeof(*length));
    if (report->marker && rows == 1)
        error("length of 'dimnames' [1] not equal to array extent");
    for (int j = 0; j < n; j++)
    {
        sources[j] = st_mab_report_column(
            st_legacy_argument(VECTOR_ELT(data, j), report->source),
            report->column, offset + j, length + j);
        if ((!report->quantile && !length[j]) ||
            (report->quantile && length[j] > 1))
            error("number of items to replace is not a multiple of "
                  "replacement length");
        if (length[j] > rows)
            error("subscript out of bounds");
    }
    if (report->quantile)
        error("argument is missing, with no default");
    R_xlen_t nd = st_mab_report_length(digits);
    double regional_digits = 0;
    if (report->regional)
    {
        if (nd)
        {
            if (TYPEOF(digits) == STRSXP)
            {
                const char *text = CHAR(STRING_ELT(digits, nd - 1));
                char *end;
                regional_digits = R_strtod(text, &end);
                if (*end)
                    error("object '%s' not found", text);
            }
            else
                regional_digits = st_legacy_number(digits, nd - 1);
        }
        nd = 1;
    }
    else
    {
        if (!nd)
            error("invalid second argument of length 0");
        if ((TYPEOF(digits) != REALSXP && TYPEOF(digits) != INTSXP &&
             TYPEOF(digits) != LGLSXP) ||
            inherits(digits, "factor"))
            error("non-numeric argument to mathematical function");
    }
    R_xlen_t size = (R_xlen_t)rows * n, count = size > nd ? size : nd;
    SEXP out = PROTECT(count == size ? allocMatrix(REALSXP, rows, n)
                                     : allocVector(REALSXP, count));
    for (R_xlen_t i = 0; i < count; i++)
    {
        R_xlen_t cell = i % size, row = cell % rows;
        int scenario = (int)(cell / rows);
        double value = NA_REAL;
        if (row < length[scenario])
        {
            SEXP source = sources[scenario];
            if (TYPEOF(source) != REALSXP && TYPEOF(source) != INTSXP &&
                TYPEOF(source) != LGLSXP)
                error("non-numeric argument to mathematical function");
            value = st_legacy_number(source, offset[scenario] + row);
        }
        double precision = report->regional ? regional_digits
                                            : st_legacy_number(digits, i % nd);
        REAL(out)[i] = fround(report->scale * value, precision);
    }
    if (count == size)
    {
        SEXP rownames = PROTECT(allocVector(STRSXP, rows));
        SET_STRING_ELT(rownames, 0, mkChar(report->marker ? "Total" : "P1"));
        for (int i = 1; i < rows; i++)
        {
            char label[32];
            snprintf(label, sizeof(label), "G%d", i - 1);
            SET_STRING_ELT(rownames, i, mkChar(label));
        }
        st_mab_report_dimnames(out, rownames, columns);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return out;
}

static SEXP st_legacy_mab_tabulate(SEXP args)
{
    SEXP data = st_legacy_argument(args, "data"),
         scenarios = st_legacy_argument(args, "scenarios"),
         digits = st_legacy_argument(args, "digits");
    if (TYPEOF(data) != VECSXP || !LENGTH(data))
        error("Simulation report data must contain scenarios");
    int n = LENGTH(data), rows = 1;
    for (int i = 0; i < n; i++)
    {
        R_xlen_t size = st_mab_report_length(
            st_legacy_argument(VECTOR_ELT(data, i), "gen.type"));
        if (size >= INT_MAX)
            error("Too many simulation generations");
        if (size + 1 > rows)
            rows = (int)size + 1;
    }
    SEXP last = VECTOR_ELT(data, n - 1);
    int regions =
            (int)st_mab_report_length(st_legacy_argument(last, "reg.chr")),
        targets =
            (int)st_mab_report_length(st_legacy_argument(last, "target.loci"));
    if (regions > (INT_MAX - 16) / 6 || targets > (INT_MAX - 16) / 6)
        error("Too many simulation report regions");
    st_mab_report *reports = st_work_alloc(
        16 + (size_t)4 * regions + (size_t)2 * targets, sizeof(*reports));
    int count = 0;
    SEXP par = st_legacy_argument(args, "par.summary"),
         ext = st_legacy_argument(args, "ext.summary");
    if (TYPEOF(ext) != LGLSXP && TYPEOF(ext) != INTSXP &&
        TYPEOF(ext) != REALSXP && TYPEOF(ext) != CPLXSXP)
        error("invalid argument type");
    if (TYPEOF(par) != LGLSXP && TYPEOF(par) != INTSXP &&
        TYPEOF(par) != REALSXP && TYPEOF(par) != CPLXSXP && !isNull(par))
        error(
            "operations are possible only for numeric, logical or complex types");
    if (!st_mab_report_length(par))
        error("argument is of length zero");
    if (st_mab_report_length(par) > 1)
        error("the condition has length > 1");
    int extended = st_mab_report_flag(ext),
        ordinary = extended ? 0 : st_mab_report_flag(par);
    if (ordinary && !extended)
        reports[count++] =
            (st_mab_report){.name = "par.summary", .parameter = 1};
    if (extended)
        reports[count++] =
            (st_mab_report){.name = "ext.summary", .parameter = 2};
    double proportion =
        asLogical(st_legacy_argument(args, "percent")) ? 100 : 1;
    const char *flags[] = {
        "TargetAlleles",     "RPG.gw.mean", "RPG.gw.Q10",
        "RPG.reg.mean",      "RPG.reg.Q10", "DonorSegments.gw",
        "DonorSegments.reg", "LinkageDrag", "MarkerAnalyses"};
    const char *sources[] = {
        "TargetAlleles",     "RPG.gw",      "RPG.gw",
        "RPG.reg",           "RPG.reg",     "DonorSegments.gw",
        "DonorSegments.reg", "LinkageDrag", "SM"};
    const char *names[] = {
        "TargetAlleles.mean",    "RPG.gw.mean", "RPG.gw.Q10", NULL,     NULL,
        "DonorSegments.gw.mean", NULL,          NULL,         "SM.mean"};
    for (int f = 0; f < 9; f++)
    {
        if (!st_mab_report_flag(st_legacy_argument(args, flags[f])))
            continue;
        int regional = f == 3 || f == 4 || f == 6 || f == 7;
        int repetitions = !regional ? 1 : f == 7 ? targets : regions;
        if (!repetitions && (f == 3 || f == 7))
            repetitions = 2;
        for (int j = 0; j < repetitions; j++)
        {
            int index =
                (f == 3 && !regions) || (f == 7 && !targets) ? 1 - j : j + 1;
            reports[count++] = (st_mab_report){
                .name = regional ? st_mab_report_name(sources[f], index,
                                                      f == 4 ? ".Q10" : ".mean")
                                 : names[f],
                .source = regional ? st_mab_report_name(sources[f], index, "")
                                   : sources[f],
                .column = f == 2 || f == 4 ? "Q10" : "Mean",
                .scale = f >= 1 && f <= 4 ? proportion
                         : f == 7         ? 100
                                          : 1,
                .regional = regional,
                .quantile = f == 4,
                .marker = f == 8};
        }
        if (f == 8)
            reports[count++] = (st_mab_report){.name = "HT.mean",
                                               .source = "HT",
                                               .column = "Mean",
                                               .scale = 1,
                                               .marker = 1};
    }
    if (!count)
        return R_NilValue;
    SEXP out = PROTECT(allocVector(VECSXP, count)),
         output_names = PROTECT(allocVector(STRSXP, count)),
         columns = PROTECT(allocVector(STRSXP, n));
    for (int i = 0; i < n; i++)
    {
        const char *left = st_mab_report_text(scenarios, i, 1),
                   *right = st_mab_report_text(scenarios, i + n, 1);
        if (!left)
            left = "NA";
        if (!right)
            right = "NA";
        size_t length = strlen(left) + strlen(right) + 2;
        char *label = st_work_alloc(length, 1);
        snprintf(label, length, "%s-%s", left, right);
        SET_STRING_ELT(columns, i, mkCharCE(label, CE_UTF8));
    }
    for (int i = 0; i < count; i++)
    {
        st_mab_report *report = reports + i;
        SET_VECTOR_ELT(
            out, i,
            report->parameter
                ? st_mab_report_parameters(data, columns,
                                           report->parameter == 2)
                : st_mab_report_numeric(data, columns, rows, report, digits));
        SET_STRING_ELT(output_names, i, mkChar(report->name));
    }
    setAttrib(out, R_NamesSymbol, output_names);
    UNPROTECT(3);
    return out;
}

static int st_legacy_population_match(const char *const *names, R_xlen_t n,
                                      SEXP input)
{
    if (TYPEOF(input) == ENVSXP || TYPEOF(input) == CLOSXP ||
        TYPEOF(input) == BUILTINSXP || TYPEOF(input) == SPECIALSXP)
        error("comparison (==) is possible only for atomic and list types");
    if (inherits(input, "data.frame"))
    {
        int found = 0;
        for (R_xlen_t j = 0; j < XLENGTH(input); j++)
        {
            int value =
                st_legacy_population_match(names, n, VECTOR_ELT(input, j));
            if (value == NA_LOGICAL)
                found = NA_LOGICAL;
            else if (value && found != NA_LOGICAL)
                found = TRUE;
        }
        return found;
    }
    SEXP supplied = input;
    int np = 0;
    if (TYPEOF(input) == SYMSXP)
    {
        if (!n)
            return FALSE;
        const char *name = translateCharUTF8(PRINTNAME(input));
        for (R_xlen_t i = 0; i < n; i++)
            if (!strcmp(names[i], name))
                return TRUE;
        return FALSE;
    }
    int type = TYPEOF(input),
        factor = type != STRSXP && inherits(input, "factor"),
        native = type == STRSXP || type == INTSXP || type == REALSXP ||
                 type == LGLSXP || type == RAWSXP || type == CPLXSXP;
    SEXP levels = factor ? getAttrib(input, R_LevelsSymbol) : R_NilValue;
    if (factor)
    {
        native = type == INTSXP && TYPEOF(levels) == STRSXP;
        for (R_xlen_t i = 0; native && i < XLENGTH(input); i++)
        {
            int code = INTEGER(input)[i];
            if (code != NA_INTEGER && (code < 1 || code > XLENGTH(levels)))
                native = 0;
        }
    }
    if (!isNull(input) && type != STRSXP && !native)
    {
        supplied =
            PROTECT(inherits(input, "factor") ? asCharacterFactor(input)
                                              : coerceVector(input, STRSXP));
        np++;
    }
    R_xlen_t m = st_mab_report_length(supplied),
             length = n && m ? (n > m ? n : m) : 0;
    const char **converted = NULL;
    if (native && !factor && type != STRSXP)
    {
        converted = st_work_alloc((size_t)m, sizeof(*converted));
        for (R_xlen_t i = 0; i < m; i++)
            converted[i] = st_mab_report_text(input, i, 0);
    }
    if (length && length % (n < m ? n : m))
        warning(
            "longer object length is not a multiple of shorter object length");
    SEXP dimensions = getAttrib(input, R_DimSymbol);
    if (!isNull(dimensions) && length && length != m)
        error("dims [product %lld] do not match the length of object [%lld]",
              (long long)m, (long long)length);
    int found = FALSE;
    for (R_xlen_t i = 0; i < length; i++)
    {
        R_xlen_t at = i % m;
        const char *name;
        if (native && factor)
        {
            int code = INTEGER(input)[at];
            SEXP label = code == NA_INTEGER ? NA_STRING
                                            : STRING_ELT(levels, code - 1);
            name = label == NA_STRING ? NULL : translateCharUTF8(label);
        }
        else if (native && type != STRSXP)
            name = converted[at];
        else
        {
            SEXP label = STRING_ELT(supplied, at);
            name = label == NA_STRING ? NULL : translateCharUTF8(label);
        }
        if (!name)
            found = NA_LOGICAL;
        else if (found != NA_LOGICAL &&
                 !strcmp(names[i % n], name))
            found = TRUE;
    }
    UNPROTECT(np);
    return found;
}

static SEXP st_legacy_population_exists(SEXP args)
{
    SEXP input = st_legacy_argument(args, "name");
    size_t bytes = 1, populations = 0;
    for (TPopListEl *entry = PopListStart; entry; entry = entry->Next)
    {
        size_t length = strlen(entry->Pop->Name);
        if (length >= SIZE_MAX - bytes)
            error("Population names are too long");
        bytes += length + 1;
        populations++;
    }
    char *buffer = st_work_alloc(bytes, 1), *at = buffer;
    size_t tokens = 0;
    for (TPopListEl *entry = PopListStart; entry; entry = entry->Next)
    {
        const char *name = entry->Pop->Name;
        for (; *name; name++)
        {
            *at++ = *name;
            tokens += *name == ' ';
        }
        *at++ = ' ';
        tokens++;
    }
    *at = '\0';
    if (populations && tokens % populations)
        error("arguments imply differing number of rows: %llu, %llu",
              (unsigned long long)tokens, (unsigned long long)populations);
    const char **names = st_work_alloc(tokens, sizeof(*names));
    size_t index = 0;
    char *begin = buffer;
    for (char *position = buffer; *position; position++)
        if (*position == ' ')
        {
            names[index++] = begin;
            *position = '\0';
            begin = position + 1;
        }
    return ScalarLogical(
        st_legacy_population_match(names, (R_xlen_t)tokens, input));
}

static SEXP st_legacy_define_description(SEXP args, int map)
{
    SEXP description = st_legacy_argument(args,"description");
    int invalid = TYPEOF(description) != STRSXP || !XLENGTH(description);
    if (!invalid)
        for (R_xlen_t i = 0; i < XLENGTH(description); i++)
            invalid |= STRING_ELT(description,i) == NA_STRING;
    if (invalid)
    {
        if (!map) error("Invalid description");
        info(-2,"Invalid map description\n");
        int np = 0;const char *key = "retval";
        SEXP result = st_sim_list(1,&key,&np);
        SET_VECTOR_ELT(result,0,ScalarInteger(-2));
        UNPROTECT(np);return result;
    }
    size_t size = 1;
    for (R_xlen_t i = 0; i < XLENGTH(description); i++)
    {
        if (STRING_ELT(description,i) == NA_STRING) error("Invalid description");
        size += strlen(CHAR(STRING_ELT(description,i))) + 2;
    }
    char *text = st_work_alloc(size,1);text[0] = '\0';
    for (R_xlen_t i = 0; i < (map ? XLENGTH(description) : 1); i++)
    {
        if (i) strcat(text,", ");
        strcat(text,CHAR(STRING_ELT(description,i)));
    }
    st_io_begin(); char *path = (char *)st_io_memory("definition");
    int status = -2, loaded = 0;
    if (map)
    {
        generate_map_file(&path,&text,&status);
        if (!status)
        {
            int disperse = 1; double factor = 100; char *output = "";
            linkage_map_load(&path,&disperse,&output,&factor,&status);loaded = 1;
        }
    }
    else
    {
        char *name = (char *)CHAR(asChar(st_legacy_argument(args,"name")));
        generate_effect_file(&path,&text); load_effmap(&path,&name);
    }
    st_io_end();
    if (!map) return ScalarInteger(0);
    int np = 0;
    SEXP out, names;PROTECT(out = allocVector(VECSXP,loaded ? 5 : 3));np++;
    PROTECT(names = allocVector(STRSXP,loaded ? 5 : 3));np++;
    for (int i = 0; i < LENGTH(out); i++) SET_STRING_ELT(names,i,mkChar(i == LENGTH(out)-1 ? "retval" : ""));
    setAttrib(out,R_NamesSymbol,names);
    SET_VECTOR_ELT(out,0,mkString(""));
    if (loaded)
    {
        SET_VECTOR_ELT(out,1,ScalarInteger(1));SET_VECTOR_ELT(out,2,mkString(""));SET_VECTOR_ELT(out,3,ScalarReal(100));
    }
    else SET_VECTOR_ELT(out,1,mkString(text));
    SET_VECTOR_ELT(out,LENGTH(out)-1,ScalarInteger(status));
    UNPROTECT(np);return out;
}

static SEXP st_legacy_output_frame(int rows, int columns, const char **keys,
                                    const int *types, int *np)
{
    SEXP out = st_sim_list(columns,keys,np);
    for (int j = 0; j < columns; j++)
    {
        SEXP column = PROTECT(allocVector(types[j],rows)); (*np)++;
        SET_VECTOR_ELT(out,j,column);
        if (types[j] == REALSXP) memset(REAL(column),0,(size_t)rows*sizeof(double));
        else if (types[j] == INTSXP) memset(INTEGER(column),0,(size_t)rows*sizeof(int));
    }
    st_sim_frame(out,rows,np);
    return out;
}
static SEXP st_legacy_locus_output(SEXP args)
{
    SEXP populations = st_legacy_argument(args,"population"), loci = st_legacy_argument(args,"loci");
    if (TYPEOF(populations) != STRSXP || !XLENGTH(populations) || TYPEOF(loci) != STRSXP)
        error("Population and loci must be character strings");
    size_t capacity = 0;
    for (R_xlen_t i = 0; i < XLENGTH(loci); i++)
    {
        if (STRING_ELT(loci,i) == NA_STRING) error("Missing locus name");
        const char *s = CHAR(STRING_ELT(loci,i));
        if (strlen(s) > SIZE_MAX-capacity-1) error("Locus name size overflow");
        capacity += strlen(s)+1;
    }
    const char **labels = st_work_alloc(capacity ? capacity : 1,sizeof(*labels));
    int n = 0;
    for (R_xlen_t i = 0; i < XLENGTH(loci); i++)
    {
        char *text = st_legacy_text(CHAR(STRING_ELT(loci,i))), *start = text;
        while (*start)
        {
            char *end = strchr(start,' ');
            if (n == INT_MAX) error("Too many loci");
            labels[n++] = start;
            if (!end) break;
            *end = '\0';start = end+1;
        }
    }
    st_sim_marker_index_type *order = st_work_alloc(n ? n : 1,sizeof(*order));
    for (int i = 0; i < n; i++) { order[i].name = labels[i];order[i].row = i; }
    qsort(order,n,sizeof(*order),st_sim_marker_index_cmp);
    for (int i = 1; i < n; i++)
        if (!strcmp(order[i-1].name,order[i].name)) error("duplicate 'row.names' are not allowed");
    int *counts = st_work_alloc((size_t)(n ? n : 1)*3,sizeof(int));
    memset(counts,0,(size_t)n*3*sizeof(int));
    char *population = st_legacy_text(CHAR(STRING_ELT(populations,0))),
         *input = XLENGTH(loci) ? st_legacy_text(CHAR(STRING_ELT(loci,0))) : st_legacy_text("");
    int bg = asInteger(st_legacy_argument(args,"background"));
    genotype_population(&population);
    population = st_legacy_text(CHAR(STRING_ELT(populations,0)));
    evaluate_locus(&population,&input,counts,&n,&bg);
    int np = 0;
    const char *keys[] = {"homozygote","heterozygote","missing"};
    const int types[] = {INTSXP,INTSXP,INTSXP};
    SEXP out = st_legacy_output_frame(n,3,keys,types,&np);
    SEXP rows = PROTECT(allocVector(STRSXP,n));np++;
    for (int i = 0; i < n; i++) SET_STRING_ELT(rows,i,mkChar(labels[i]));
    setAttrib(out,R_RowNamesSymbol,rows);
    for (int j = 0; j < 3; j++) memcpy(INTEGER(VECTOR_ELT(out,j)),counts+(size_t)n*j,(size_t)n*sizeof(int));
    UNPROTECT(np);return out;
}
static SEXP st_legacy_select_best(SEXP args)
{
    SEXP populations = st_legacy_argument(args,"population");
    char *population = st_legacy_text(CHAR(asChar(populations)));
    SEXP effect_arg = st_legacy_argument(args,"effect");
    char *effect = isNull(effect_arg) ? "" : st_legacy_text(CHAR(asChar(effect_arg))),
         *criterion = st_legacy_text(CHAR(asChar(st_legacy_argument(args,"criterion"))));
    int n = asInteger(st_legacy_argument(args,"n")), status = -2;
    if (n == NA_INTEGER) error("Missing selection count");
    double minimum = 0, maximum = 0;
    char text[64];memset(text,'a',sizeof(text)-1);text[sizeof(text)-1] = '\0';char *index = text;
    select_all_n_best(&population,&index,&n,&minimum,&maximum,&criterion,&effect,&status);
    if (status) return R_NilValue;
    char *end;double chosen = strtod(text,&end);
    int np = 0;const char *keys[] = {"n","minscore","maxscore"};const int types[] = {INTSXP,REALSXP,REALSXP};
    SEXP out = st_legacy_output_frame(1,3,keys,types,&np);
    INTEGER(VECTOR_ELT(out,0))[0] = end == text || *end || chosen <= INT_MIN || chosen > INT_MAX ? NA_INTEGER : (int)chosen;
    REAL(VECTOR_ELT(out,1))[0] = minimum;REAL(VECTOR_ELT(out,2))[0] = maximum;
    UNPROTECT(np);return out;
}

static SEXP st_legacy_output(SEXP args)
{
    const char *operation = CHAR(asChar(st_legacy_argument(args,"operation")));
    int np = 0; SEXP out = R_NilValue;
    if (!strcmp(operation,"locus"))
        return st_legacy_locus_output(args);
    if (!strcmp(operation,"select.best"))
        return st_legacy_select_best(args);
    if (!strcmp(operation,"map"))
    {
        int rows = 0;
        for (int c = 0; c < NoChroms; c++)
        {
            if (PopDescr[c].NrMappoints > INT_MAX-rows) error("Too many map points");
            rows += PopDescr[c].NrMappoints;
        }
        const char *keys[] = {"chrom","pos","name","class"};
        const int types[] = {REALSXP,REALSXP,STRSXP,STRSXP};
        out = st_legacy_output_frame(rows,4,keys,types,&np);
        int at = 0;
        for (int c = 0; c < NoChroms; c++)
            for (int j = 0; j < PopDescr[c].NrMappoints; j++,at++)
            {
                TMappoint *p = PopDescr[c].Mappoint+j;
                REAL(VECTOR_ELT(out,0))[at] = c+1;
                REAL(VECTOR_ELT(out,1))[at] = p->Pos;
                SET_STRING_ELT(VECTOR_ELT(out,2),at,mkChar(p->name));
                SET_STRING_ELT(VECTOR_ELT(out,3),at,mkChar(p->class));
            }
    }
    else if (!strcmp(operation,"genome"))
    {
        const char *keys[] = {"no.chrom","no.hom","chrom.len"};
        out = st_sim_list(3,keys,&np);
        SET_VECTOR_ELT(out,0,ScalarInteger(NoChroms));
        SET_VECTOR_ELT(out,1,ScalarInteger(NoHoms));
        SEXP lengths = PROTECT(allocVector(REALSXP,NoChroms)); np++;
        SET_VECTOR_ELT(out,2,lengths);
        for (int c = 0; c < NoChroms; c++) REAL(lengths)[c] = PopDescr[c].ChromLen;
    }
    else if (!strcmp(operation,"effects") || !strcmp(operation,"populations"))
    {
        int effects = !strcmp(operation,"effects"), rows = effects ? NrEffects : 0;
        if (!effects)
            for (TPopListEl *p = PopListStart; p; p = p->Next)
            {
                if (rows == INT_MAX) error("Too many populations");
                rows++;
            }
        const char *keys[] = {effects ? "effect" : "PopName",effects ? "weight" : "count"};
        const int types[] = {STRSXP,effects ? REALSXP : INTSXP};
        out = st_legacy_output_frame(rows,2,keys,types,&np);
        if (effects)
            for (int i = 0; i < rows; i++)
            {
                SET_STRING_ELT(VECTOR_ELT(out,0),i,mkChar(Effects[i].Name));
                REAL(VECTOR_ELT(out,1))[i] = Effects[i].weight;
            }
        else
        {
            int i = 0;
            for (TPopListEl *p = PopListStart; p; p = p->Next,i++)
            {
                SET_STRING_ELT(VECTOR_ELT(out,0),i,mkChar(p->Pop->Name));
                INTEGER(VECTOR_ELT(out,1))[i] = p->Pop->NoInds > INT_MAX ? NA_INTEGER : (int)p->Pop->NoInds;
            }
        }
    }
    else if (!strcmp(operation,"mdp"))
    {
        char buffer[64]; memset(buffer,'a',sizeof(buffer)-1);buffer[sizeof(buffer)-1] = '\0';
        char *p = buffer;get_mdp(&p);char *end;double v = strtod(buffer,&end);
        out = PROTECT(ScalarInteger(end == buffer || *end || v <= INT_MIN || v > INT_MAX ? NA_INTEGER : (int)v));np++;
    }
    else
    {
        SEXP name_arg = st_legacy_argument(args,"population");
        if (TYPEOF(name_arg) != STRSXP || !XLENGTH(name_arg)) error("Population name is required");
        if (!strcmp(operation,"size"))
        {
            POPULATION n = 0, allocated = 0;
            for (R_xlen_t j = 0; j < XLENGTH(name_arg); j++)
            {
                char *text = st_legacy_text(CHAR(STRING_ELT(name_arg,j)));
                for (char *name = strtok(text," "); name; name = strtok(NULL," "))
                {
                    TPopListEl *p = NULL;find_list_element(name,&p);
                    if (!p) { info(-2,"Population not defined: %s\n",name);return R_NilValue; }
                    n += p->Pop->NoInds; allocated += p->Pop->NoIndsAlloc;
                }
            }
            const char *keys[] = {"NoInds","NoIndsAlloc"};const int types[] = {INTSXP,INTSXP};
            out = st_legacy_output_frame(1,2,keys,types,&np);
            INTEGER(VECTOR_ELT(out,0))[0] = n > INT_MAX ? NA_INTEGER : (int)n;
            INTEGER(VECTOR_ELT(out,1))[0] = allocated > INT_MAX ? NA_INTEGER : (int)allocated;
        }
        else
        {
            char *name = st_legacy_text(CHAR(STRING_ELT(name_arg,0)));
            TPopListEl *population = NULL;find_list_element(name,&population);
            if (!population)
            {
                info(-2,"Population not defined: %s\n",name);
                if (!strcmp(operation,"population"))
                {
                    const char *keys[] = {"ind","chrom","hom","pos","all"};
                    const int types[] = {INTSXP,INTSXP,INTSXP,REALSXP,INTSXP};
                    out = st_legacy_output_frame(0,5,keys,types,&np);
                    UNPROTECT(np);return out;
                }
                error("invalid 'times' argument");
            }
            if (!strcmp(operation,"population"))
            {
                int rows = 0;
                for (POPULATION i = 0; i < population->Pop->NoInds; i++)
                    for (int c = 0; c < NoChroms; c++)
                        for (int h = 0; h < NoHoms; h++)
                        {
                            int n = get_individual(population->Pop,i)->Chrom[c].Hom[h].NoLociUsed;
                            if (n > INT_MAX-rows) error("Population output exceeds integer row range");
                            rows += n;
                        }
                const char *keys[] = {"ind","chrom","hom","pos","all"};
                const int types[] = {INTSXP,INTSXP,INTSXP,REALSXP,INTSXP};
                out = st_legacy_output_frame(rows,5,keys,types,&np);
                get_population(&name,INTEGER(VECTOR_ELT(out,0)),INTEGER(VECTOR_ELT(out,1)),INTEGER(VECTOR_ELT(out,2)),REAL(VECTOR_ELT(out,3)),INTEGER(VECTOR_ELT(out,4)));
            }
            else if (!strcmp(operation,"gvalue") || !strcmp(operation,"pvalue") || !strcmp(operation,"ld"))
            {
                if (population->Pop->NoInds > INT_MAX) error("Population output exceeds integer row range");
                int rows = population->Pop->NoInds;
                const char *key = !strcmp(operation,"ld") ? "length" : operation;
                const int type = REALSXP;out = st_legacy_output_frame(rows,1,&key,&type,&np);
                if (!strcmp(operation,"ld"))
                {
                    int allele = asInteger(st_legacy_argument(args,"allele")), chrom = asInteger(st_legacy_argument(args,"chromosome"));
                    double pos = asReal(st_legacy_argument(args,"position"));
                    evaluate_ld(&name,&allele,&chrom,&pos,REAL(VECTOR_ELT(out,0)));
                }
                else
                {
                    SEXP effect_arg = st_legacy_argument(args,"effect");
                    char *effect = isNull(effect_arg) ? "" : st_legacy_text(CHAR(asChar(effect_arg)));
                    char *criterion = !strcmp(operation,"gvalue") ? "G" : "P";int status = -2;
                    get_population_value(&name,REAL(VECTOR_ELT(out,0)),&effect,&criterion,&status);
                }
            }
            else error("Unknown simulation output operation");
        }
    }
    UNPROTECT(np);return out;
}

static SEXP st_legacy_info_cat(SEXP args)
{
    SEXP environment = st_legacy_argument(args, "environment");
    if (TYPEOF(environment) != ENVSXP)
        error("Information output needs its calling environment");
    int np = 0, level;
    get_info_level(&level);
    SEXP requested = PROTECT(eval(install("lev"), environment));
    np++;
    int direct = !isObject(requested) &&
        (TYPEOF(requested) == INTSXP || TYPEOF(requested) == REALSXP ||
         TYPEOF(requested) == LGLSXP) && XLENGTH(requested) == 1;
    if (direct)
    {
        double value = TYPEOF(requested) == REALSXP ? REAL(requested)[0]
                       : TYPEOF(requested) == INTSXP ? INTEGER(requested)[0]
                                                     : LOGICAL(requested)[0];
        if (ISNAN(value) || (TYPEOF(requested) != REALSXP && value == NA_INTEGER))
            error("missing value where TRUE/FALSE needed");
        if (value < level)
        {
            UNPROTECT(np);
            return R_NilValue;
        }
    }
    SEXP format = PROTECT(mkString("M: %s\r"));
    np++;
    SEXP text = PROTECT(lang3(install("sprintf"), format, install("msg")));
    np++;
    SEXP newline = PROTECT(ScalarLogical(FALSE));
    np++;
    SEXP message = PROTECT(lang3(install("message"), text, newline));
    np++;
    SET_TAG(CDDR(message), install("appendLF"));
    if (direct)
        eval(message, environment);
    else
    {
        SEXP operand = PROTECT(ScalarInteger(level));
        np++;
        SEXP comparison = PROTECT(lang3(install(">="), install("lev"), operand));
        np++;
        SEXP condition = PROTECT(lang3(install("if"), comparison, message));
        np++;
        eval(condition, environment);
    }
    UNPROTECT(np);
    return R_NilValue;
}

static int st_legacy_haplotype_index(SEXP value, int argument)
{
    if (TYPEOF(value) != INTSXP || !XLENGTH(value))
        error("Haplotype locus index must be an integer");
    for (R_xlen_t i = 0; i < XLENGTH(value); i++)
        if (INTEGER(value)[i] == NA_INTEGER)
            error("NAs in foreign function call (arg %d)", argument);
    return INTEGER(value)[0];
}

static SEXP st_legacy_haplotype_next(SEXP args)
{
    int a = st_legacy_haplotype_index(st_legacy_argument(args, "locusA"), 1),
        b = st_legacy_haplotype_index(st_legacy_argument(args, "locusB"), 2);
    SEXP out = PROTECT(allocVector(INTSXP, 5));
    haplotype_next(&a, &b, INTEGER(out), INTEGER(out) + 1, INTEGER(out) + 2,
                   INTEGER(out) + 3, INTEGER(out) + 4);
    UNPROTECT(1);
    return out;
}

static SEXP st_legacy_haplotype_d(void)
{
    SEXP out = PROTECT(allocVector(REALSXP, 3));
    haplotype_calc_d(REAL(out), REAL(out) + 2, REAL(out) + 1);
    UNPROTECT(1);
    return out;
}

static SEXP st_legacy_population_info(SEXP args)
{
    SEXP population = st_legacy_argument(args, "population"),
         individual = st_legacy_argument(args, "individual");
    if (TYPEOF(population) != STRSXP || !XLENGTH(population) ||
        TYPEOF(individual) != STRSXP || !XLENGTH(individual))
        error("Population information needs a population name and individual");
    char *name = (char *)CHAR(STRING_ELT(population, 0)),
         *index = (char *)CHAR(STRING_ELT(individual, 0));
    char buffer[257];
    memset(buffer, 'a', 256);
    buffer[256] = '\0';
    char *text = buffer;
    get_population_info(&name, &index, &text);
    return mkString(buffer);
}

static SEXP st_legacy_remove_all(void)
{
    if (PopListStart)
    {
        size_t length = 1;
        for (TPopListEl *p = PopListStart; p; p = p->Next)
        {
            size_t n = strlen(p->Pop->Name);
            if (n >= SIZE_MAX - length)
                error("Population name list is too large");
            length += n + 1;
        }
        char *names = st_work_alloc(length, sizeof(char));
        size_t at = 0;
        for (TPopListEl *p = PopListStart; p; p = p->Next)
        {
            if (at) names[at++] = ' ';
            size_t n = strlen(p->Pop->Name);
            memcpy(names + at, p->Pop->Name, n);
            at += n;
        }
        names[at] = '\0';
        remove_population(&names);
    }
    int level = 1;
    char *message = "All populations removed";
    r_info(&level, &message);
    return R_NilValue;
}

SEXP st_legacy_sim(SEXP args)
{
    const char *action = CHAR(asChar(st_legacy_argument(args, "action")));
    if (!strcmp(action, "output"))
        return st_legacy_output(args);
    if (!strcmp(action, "info.cat"))
        return st_legacy_info_cat(args);
    if (!strcmp(action, "remove.all"))
        return st_legacy_remove_all();
    if (!strcmp(action, "hap.next"))
        return st_legacy_haplotype_next(args);
    if (!strcmp(action, "hap.d"))
        return st_legacy_haplotype_d();
    if (!strcmp(action, "population.info"))
        return st_legacy_population_info(args);
    if (!strcmp(action, "allow.zeros"))
    {
        int allow;
        get_allow_zeros(&allow);
        return ScalarInteger(allow);
    }
    if (!strcmp(action, "define.map"))
        return st_legacy_define_description(args,1);
    if (!strcmp(action, "define.effects"))
        return st_legacy_define_description(args,0);
    if (!strcmp(action, "marker.distance"))
        return st_legacy_marker_distance(args);
    if (!strcmp(action, "population.exists"))
        return st_legacy_population_exists(args);
    if (!strcmp(action, "performance.copy"))
        return st_legacy_performance_copy(args);
    if (!strcmp(action, "target.classify"))
        return st_legacy_target_export(args);
    if (!strcmp(action, "select.genotype"))
        return st_legacy_select_genotype(args);
    if (!strcmp(action, "effect.locus"))
        return st_legacy_locus_effect(args);
    if (!strcmp(action, "genome.evaluate"))
        return st_legacy_genome_evaluate(args);
    if (!strcmp(action, "select.segments"))
        return st_legacy_select_segments(args);
    if (!strcmp(action, "cross.plan.execute"))
        return st_legacy_cross_plan_execute(args);
    if (!strcmp(action, "single.cross"))
        return st_legacy_single_cross(args);
    if (!strcmp(action, "mdp"))
        return st_legacy_mdp(args);
    if (!strcmp(action, "rng.init"))
        return st_legacy_rng_init();
    if (!strcmp(action, "plant"))
        return st_legacy_plant(args);
    if (!strcmp(action, "sdev"))
        return st_legacy_sdev(args);
    if (!strcmp(action, "data"))
        return st_legacy_data_operation(args);
    if (!strcmp(action, "mab"))
        return st_legacy_mab(args);
    if (!strcmp(action, "mab.individuals"))
    {
        SEXP names = st_legacy_argument(args, "population");
        if (TYPEOF(names) == STRSXP && LENGTH(names) == 1)
        {
            char *text = (char *)st_legacy_text(CHAR(STRING_ELT(names, 0)));
            genotype_population(&text);
        }
        int success = 0;
        st_population_export(names, st_legacy_argument(args, "background"), args, &success);
        if (!success) error("argument of length 0");
        return R_NilValue;
    }
    if (!strcmp(action, "mab.tabulate"))
        return st_legacy_mab_tabulate(args);
    if (!strcmp(action, "allele"))
        return st_legacy_evaluate(args, 0);
    if (!strcmp(action, "genotype"))
        return st_legacy_evaluate(args, 1);
    if (!strcmp(action, "summary"))
        return st_legacy_population_summary(args);
    if (!strcmp(action, "map"))
        return st_legacy_map(args, 0);
    if (!strcmp(action, "genome"))
        return st_legacy_map(args, 1);
    if (!strcmp(action, "pcoa"))
        return st_legacy_pcoa(args);
    error("Unknown legacy simulation operation");
    return R_NilValue;
}

typedef struct
{
    int *value;
    int width;
    size_t count;
} st_legacy_genotype_row;

typedef struct st_legacy_allele_row
{
    int locus, allele, count;
    struct st_legacy_allele_row *next;
} st_legacy_allele_row;

static int st_legacy_genotype_compare(const void *a, const void *b)
{
    const st_legacy_genotype_row *x = a, *y = b;
    for (int j = 0; j < x->width; j++)
        if (x->value[j] != y->value[j])
            return (x->value[j] > y->value[j]) ? 1 : -1;
    return 0;
}

static int st_legacy_allele_compare(const void *a, const void *b)
{
    const st_legacy_allele_row *x = *(st_legacy_allele_row *const *)a;
    const st_legacy_allele_row *y = *(st_legacy_allele_row *const *)b;
    if (!x->allele || !y->allele)
        return (x->allele > y->allele) - (x->allele < y->allele);
    if (x->locus != y->locus)
        return (x->locus < y->locus) - (x->locus > y->locus);
    return (x->allele > y->allele) - (x->allele < y->allele);
}

static int st_legacy_individual_allele(TInd *ind, TMappoint *point, int hom)
{
    if (ind->genotype)
        return ind->genotype[point->counter * NoHoms + hom];
    return get_allele(&ind->Chrom[point->chrom - 1].Hom[hom], point->Pos);
}

static SEXP st_legacy_evaluate(SEXP args, int genotype)
{
    SEXP pop_arg = st_legacy_argument(args, "population");
    SEXP loci_arg = st_legacy_argument(args, "loci");
    int bg = asInteger(st_legacy_argument(args, "background"));
    if (!cMP)
    {
        info(-2, "No linkage map defined\n");
        return R_NilValue;
    }
    if (TYPEOF(pop_arg) != STRSXP || XLENGTH(pop_arg) != 1 ||
        STRING_ELT(pop_arg, 0) == NA_STRING || TYPEOF(loci_arg) != STRSXP)
        error("Population and locus names must be character strings");
    TPopListEl *el = NULL;
    find_list_element((char *)CHAR(STRING_ELT(pop_arg, 0)), &el);
    if (!el || !el->Pop->NoInds)
    {
        info(-2, "Population not defined or empty\n");
        return R_NilValue;
    }
    int missing = asInteger(st_legacy_argument(args, "missing")) == 1;
    int drop_empty = asLogical(st_legacy_argument(args, "drop.empty")) == 1;
    size_t ni = (size_t)el->Pop->NoInds;
    if (XLENGTH(loci_arg) > INT_MAX)
        error("Too many loci");
    TMappoint **points = st_work_alloc(
        XLENGTH(loci_arg) ? XLENGTH(loci_arg) : 1, sizeof(TMappoint *));
    int nloci = 0;
    for (R_xlen_t i = 0; i < XLENGTH(loci_arg); i++)
    {
        if (STRING_ELT(loci_arg, i) == NA_STRING)
            error("Locus names must not be missing");
        TMappoint *point = find_mappoint((char *)CHAR(STRING_ELT(loci_arg, i)));
        if (!point)
            continue;
        if (!genotype && missing && drop_empty)
        {
            int observed = 0;
            for (size_t j = 0; j < ni && !observed; j++)
            {
                TInd *ind = get_individual(el->Pop, j);
                int complete = NoHoms == 2;
                for (int h = 0; h < NoHoms; h++)
                    if (st_legacy_individual_allele(ind, point, h) == bg)
                        complete = 0;
                observed |= complete;
            }
            if (!observed)
                continue;
        }
        points[nloci++] = point;
    }
    if (!nloci)
    {
        info(-2, "Loci not found\n");
        return R_NilValue;
    }
    int nprotect = 0;
    if (genotype)
    {
        if (NoHoms < 1 || nloci > (INT_MAX - 2) / NoHoms)
            error("Too many genotype columns");
        int width = nloci * NoHoms;
        int mode = asInteger(st_legacy_argument(args, "mode"));
        int omit = asInteger(st_legacy_argument(args, "omit")) == 1;
        int *per_chrom = st_work_alloc(NoChroms, sizeof(int));
        memset(per_chrom, 0, NoChroms * sizeof(int));
        for (int j = 0; j < nloci; j++)
            per_chrom[points[j]->chrom - 1]++;
        int *values = st_work_alloc(ni * (size_t)width, sizeof(int));
        st_legacy_genotype_row *rows = st_work_alloc(ni, sizeof(*rows));
        size_t kept = 0;
        for (size_t i = 0; i < ni; i++)
        {
            TInd *ind = get_individual(el->Pop, i);
            int *row = values + kept * width, complete = 1;
            for (int j = 0; j < nloci; j++)
            {
                for (int h = 0; h < NoHoms; h++)
                {
                    int value = st_legacy_individual_allele(ind, points[j], h);
                    row[j * NoHoms + h] = value;
                    if (value == bg)
                        complete = 0;
                }
                if (mode != 0 &&
                    (mode != 1 || per_chrom[points[j]->chrom - 1] == 1))
                    R_isort(row + j * NoHoms, NoHoms);
            }
            if (omit && !complete)
                continue;
            rows[kept++] = (st_legacy_genotype_row){row, width, 1};
        }
        qsort(rows, kept, sizeof(*rows), st_legacy_genotype_compare);
        size_t nr = 0;
        for (size_t i = 0; i < kept; i++)
        {
            if (nr && !st_legacy_genotype_compare(&rows[nr - 1], &rows[i]))
                rows[nr - 1].count++;
            else
                rows[nr++] = rows[i];
        }
        if (nr > INT_MAX)
            error("Too many observed genotypes");
        SEXP selected_arg = st_legacy_argument(args, "alleles");
        size_t selected = SIZE_MAX;
        if (XLENGTH(selected_arg) && st_legacy_number(selected_arg, 0) != 0)
        {
            if (XLENGTH(selected_arg) > width + 2)
                error("Allele pattern exceeds the genotype result width");
            for (size_t i = 0; i < nr; i++)
            {
                int matches = 1;
                for (R_xlen_t j = 0; j < XLENGTH(selected_arg); j++)
                {
                    double value = j < width    ? rows[i].value[j]
                                   : j == width ? (double)rows[i].count
                                                : (double)rows[i].count / kept;
                    if (value != st_legacy_number(selected_arg, j))
                        matches = 0;
                }
                if (matches)
                {
                    selected = i;
                    break;
                }
            }
        }
        SEXP names = PROTECT(allocVector(STRSXP, width + 2));
        nprotect++;
        char text[1024];
        for (int j = 0; j < nloci; j++)
            for (int h = 0; h < NoHoms; h++)
            {
                snprintf(text, sizeof(text), "%s-%d ", points[j]->name, h + 1);
                SET_STRING_ELT(names, j * NoHoms + h, mkChar(text));
            }
        SET_STRING_ELT(names, width, mkChar("count"));
        SET_STRING_ELT(names, width + 1, mkChar("frequency"));
        SEXP out;
        if (selected != SIZE_MAX)
        {
            const char *keys[] = {"result"};
            out = st_sim_list(1, keys, &nprotect);
            SEXP value = PROTECT(allocVector(REALSXP, width + 2));
            nprotect++;
            for (int j = 0; j < width; j++)
                REAL(value)[j] = rows[selected].value[j];
            REAL(value)[width] = rows[selected].count;
            REAL(value)[width + 1] = (double)rows[selected].count / kept;
            SET_VECTOR_ELT(out, 0, value);
            st_sim_frame(out, width + 2, &nprotect);
            setAttrib(out, R_RowNamesSymbol, names);
        }
        else
        {
            out = PROTECT(allocVector(VECSXP, width + 2));
            nprotect++;
            setAttrib(out, R_NamesSymbol, names);
            for (int j = 0; j < width + 2; j++)
            {
                SEXP col = PROTECT(allocVector(REALSXP, nr));
                for (size_t i = 0; i < nr; i++)
                    REAL(col)
                    [i] = j < width    ? rows[i].value[j]
                          : j == width ? (double)rows[i].count
                                       : (double)rows[i].count / kept;
                SET_VECTOR_ELT(out, j, col);
                UNPROTECT(1);
            }
            st_sim_frame(out, (int)nr, &nprotect);
            if (!omit)
            {
                SEXP rn = PROTECT(allocVector(STRSXP, nr));
                nprotect++;
                for (size_t i = 0; i < nr; i++)
                {
                    snprintf(text, sizeof(text), "%zu ", i + 1);
                    SET_STRING_ELT(rn, i, mkChar(text));
                }
                setAttrib(out, R_RowNamesSymbol, rn);
            }
        }
        UNPROTECT(nprotect);
        return out;
    }
    int frequencies = asLogical(st_legacy_argument(args, "frequencies")) == 1;
    int *counts = st_work_alloc((size_t)USHRT_MAX + 1, sizeof(int));
    size_t nr = 0;
    int zero_count = 0;
    st_legacy_allele_row *first = NULL, **tail = &first;
    for (int j = 0; j < nloci; j++)
    {
        memset(counts, 0, ((size_t)USHRT_MAX + 1) * sizeof(int));
        for (size_t i = 0; i < ni; i++)
        {
            TInd *ind = get_individual(el->Pop, i);
            for (int h = 0; h < NoHoms; h++)
            {
                int value = st_legacy_individual_allele(ind, points[j], h);
                if (value < 0 || value > USHRT_MAX)
                    error("Allele code exceeds the simulation range");
                if (!missing || value != bg)
                {
                    if (counts[value] == INT_MAX)
                        error("Allele count exceeds the integer range");
                    counts[value]++;
                }
            }
        }
        if (counts[0] > INT_MAX - zero_count)
            error("Allele count exceeds the integer range");
        zero_count += counts[0];
        for (int value = 1; value <= USHRT_MAX; value++)
            if (counts[value])
            {
                st_legacy_allele_row *row = st_work_alloc(1, sizeof(*row));
                *row = (st_legacy_allele_row){j, value, counts[value], NULL};
                *tail = row;
                tail = &row->next;
                nr++;
            }
    }
    if (zero_count)
    {
        st_legacy_allele_row *row = st_work_alloc(1, sizeof(*row));
        *row = (st_legacy_allele_row){0, 0, zero_count, NULL};
        *tail = row;
        nr++;
    }
    if (nr > INT_MAX)
        error("Too many observed alleles");
    st_legacy_allele_row **rows = st_work_alloc(nr ? nr : 1, sizeof(*rows));
    size_t at = 0;
    double total = 0;
    for (st_legacy_allele_row *row = first; row; row = row->next)
    {
        rows[at++] = row;
        total += row->count;
    }
    qsort(rows, nr, sizeof(*rows), st_legacy_allele_compare);
    int nc = nloci + 1 + frequencies;
    SEXP out = PROTECT(allocVector(VECSXP, nc));
    SEXP names = PROTECT(allocVector(STRSXP, nc));
    SEXP rn = PROTECT(allocVector(STRSXP, nr));
    nprotect += 3;
    for (int j = 0; j < nc; j++)
    {
        SEXP col = PROTECT(allocVector(j > nloci ? REALSXP : INTSXP, nr));
        for (size_t i = 0; i < nr; i++)
        {
            if (j < nloci)
                INTEGER(col)[i] = rows[i]->locus == j ? rows[i]->allele : 0;
            else if (j == nloci)
                INTEGER(col)[i] = rows[i]->count;
            else
                REAL(col)[i] = rows[i]->count / total;
        }
        SET_VECTOR_ELT(out, j, col);
        SET_STRING_ELT(names, j,
                       mkChar(j < nloci    ? points[j]->name
                              : j == nloci ? "count"
                                           : "freq"));
        UNPROTECT(1);
    }
    for (size_t i = 0; i < nr; i++)
    {
        char text[32];
        snprintf(text, sizeof(text), "%zu ", i + 1);
        SET_STRING_ELT(rn, i, mkChar(text));
    }
    setAttrib(out, R_NamesSymbol, names);
    st_sim_frame(out, (int)nr, &nprotect);
    setAttrib(out, R_RowNamesSymbol, rn);
    UNPROTECT(nprotect);
    return out;
}

typedef struct
{
    const char *name, *first, *second;
} st_mab_parent_marker;

typedef struct
{
    char population[40], first[40], second[40], third[40], selected[40],
        total[40];
    const char *generation, *strategy;
    int size, selected_n, preselected_n;
} st_mab_generation;

static int st_mab_parent_compare(const void *a, const void *b)
{
    return strcmp(((const st_mab_parent_marker *)a)->name,
                  ((const st_mab_parent_marker *)b)->name);
}

static st_mab_parent_marker *st_mab_parent_read(SEXP data, int *n)
{
    if (TYPEOF(data) != VECSXP || LENGTH(data) < 2 || LENGTH(data) > 3)
        error("Parent marker input must contain two or three columns");
    *n = LENGTH(VECTOR_ELT(data, 0));
    for (int j = 0; j < LENGTH(data); j++)
        if (TYPEOF(VECTOR_ELT(data, j)) != STRSXP ||
            LENGTH(VECTOR_ELT(data, j)) != *n)
            error(
                "Parent marker columns must be character vectors of equal length");
    st_mab_parent_marker *rows = st_work_alloc(*n ? *n : 1, sizeof(*rows));
    for (int i = 0; i < *n; i++)
    {
        for (int j = 0; j < LENGTH(data); j++)
            if (STRING_ELT(VECTOR_ELT(data, j), i) == NA_STRING)
                error("Parent marker values must not be missing");
        rows[i].name = CHAR(STRING_ELT(VECTOR_ELT(data, 0), i));
        rows[i].first = CHAR(STRING_ELT(VECTOR_ELT(data, 1), i));
        rows[i].second =
            CHAR(STRING_ELT(VECTOR_ELT(data, LENGTH(data) - 1), i));
    }
    qsort(rows, *n, sizeof(*rows), st_mab_parent_compare);
    for (int i = 1; i < *n; i++)
        if (!strcmp(rows[i - 1].name, rows[i].name))
            error("Parent marker names must be unique");
    return rows;
}

static st_mab_parent_marker *st_mab_parent_find(st_mab_parent_marker *rows,
                                                int n, const char *name)
{
    st_mab_parent_marker key = {name, NULL, NULL};
    return bsearch(&key, rows, n, sizeof(*rows), st_mab_parent_compare);
}

static int st_mab_contains(SEXP strings, const char *name)
{
    if (isNull(strings))
        return 0;
    if (TYPEOF(strings) != STRSXP)
        error("Locus lists must contain character names");
    for (R_xlen_t i = 0; i < XLENGTH(strings); i++)
        if (STRING_ELT(strings, i) != NA_STRING &&
            !strcmp(CHAR(STRING_ELT(strings, i)), name))
            return 1;
    return 0;
}

static void st_mab_uniform_effect(const char *name, int allele)
{
    if (find_Effect((char *)name))
        error("MAB effect already exists: %s", name);
    size_t count = 0;
    for (int c = 0; c < NoChroms; c++)
        for (int i = 0; i < PopDescr[c].NrMappoints; i++)
            if (strstr(PopDescr[c].Mappoint[i].class, name))
                count++;
    if (count > INT_MAX)
        error("Too many MAB effect loci");
    TEffects *resized = realloc(Effects, (NrEffects + 1) * sizeof(*Effects));
    if (!resized)
        error("Cannot allocate MAB effects");
    Effects = resized;
    TEffects *effect = &Effects[NrEffects];
    memset(effect, 0, sizeof(*effect));
    effect->weight = 1;
    effect->Pos = NrEffects + 1;
    snprintf(effect->Name, sizeof(effect->Name), "%s", name);
    effect->Eff = calloc(count ? count : 1, sizeof(TEff));
    effect->MPset = calloc(count ? count : 1, sizeof(TMappoint *));
    NrEffects++;
    if (!effect->Eff || !effect->MPset)
        error("Cannot allocate MAB locus effects");
    for (int c = 0; c < NoChroms; c++)
        for (int i = 0; i < PopDescr[c].NrMappoints; i++)
        {
            TMappoint *point = &PopDescr[c].Mappoint[i];
            if (!strstr(point->class, name))
                continue;
            TEff *entry = &effect->Eff[effect->NrEff++];
            entry->Eff = 1;
            entry->Loc = calloc(1, sizeof(TGen));
            if (!entry->Loc)
                error("Cannot allocate MAB genotype effect");
            entry->NrLoci = 1;
            entry->Loc[0].MP = point;
            entry->Loc[0].NrAll = 1;
            entry->Loc[0].All = malloc(NoHoms * sizeof(int));
            if (!entry->Loc[0].All)
                error("Cannot allocate MAB effect alleles");
            entry->Loc[0].All[0] = allele;
            for (int h = 1; h < NoHoms; h++)
                entry->Loc[0].All[h] = -1;
            effect->MPset[effect->MPs++] = point;
        }
    if (!st_sim_ef_finalize(effect))
        error("Cannot expand MAB allele effects");
    GetRNGstate();
    PutRNGstate();
    st_clear_all_population_gvalue();
}

static void st_mab_init_homo(char *name, int allele)
{
    size_t n = (size_t)NoChroms * NoHoms;
    void *mark = st_work_mark();
    int *ind = st_work_alloc(n, sizeof(int)),
        *chrom = st_work_alloc(n, sizeof(int));
    int *hom = st_work_alloc(n, sizeof(int)),
        *all = st_work_alloc(n, sizeof(int));
    double *pos = st_work_alloc(n, sizeof(double));
    for (int c = 0, at = 0; c < NoChroms; c++)
        for (int h = 0; h < NoHoms; h++, at++)
        {
            ind[at] = 1;
            chrom[at] = c + 1;
            hom[at] = h + 1;
            pos[at] = 0;
            all[at] = allele;
        }
    char text[32], *number = text;
    snprintf(text, sizeof(text), "%zu", n);
    init_population(&name, &number, ind, chrom, hom, pos, all);
    st_work_reset(mark);
}

static POPULATION st_mab_population_size(const char *name)
{
    TPopListEl *el = NULL;
    find_list_element((char *)name, &el);
    if (!el)
        error("MAB population not available: %s", name);
    return el->Pop->NoInds;
}

static void st_mab_sample(char *target, char *source, int n, int replace)
{
    char text[32], *number = text;
    snprintf(text, sizeof(text), "%d", n);
    sample_population(&target, &source, &number, &replace);
}

static void st_mab_select_parameters(char *target, char *source,
                                     const char *effect_name, int n,
                                     int all_best, int segments, int allele,
                                     int decreasing)
{
    if (st_mab_population_size(source) == 0)
    {
        info(-2, "Population %s is of size 0\n", source);
        return;
    }
    if (segments && (POPULATION)n > st_mab_population_size(source))
    {
        info(-2, "Population %s has too few individuals for selection\n",
             source);
        return;
    }
    char *effect = (char *)effect_name, *criterion = "G", *none = "";
    char *name = st_legacy_text(source);
    if (segments)
    {
        size_t count = st_mab_population_size(source);
        double *ratio = st_work_alloc(count, sizeof(double));
        double *length = st_work_alloc(count, sizeof(double));
        double *total = st_work_alloc(count, sizeof(double));
        double *value = st_work_alloc(count, sizeof(double));
        int *blocks = st_work_alloc(count, sizeof(int));
        memset(ratio, 0, count * sizeof(double));
        memset(length, 0, count * sizeof(double));
        memset(total, 0, count * sizeof(double));
        memset(blocks, 0, count * sizeof(int));
        int nvalue = (int)count, status;
        evaluate_genome(&name, &allele, ratio, length, total, blocks);
        for (size_t i = 0; i < count; i++)
            value[i] = blocks[i];
        genotype_population(&name);
        name = st_legacy_text(source);
        evaluate_population(&name, NULL);
        name = st_legacy_text(source);
        set_population_gvalue(&name, value, &nvalue, &status);
        if (status)
            error("Cannot set MAB segment-selection values");
    }
    else
    {
        genotype_population(&name);
        name = st_legacy_text(source);
        evaluate_population(&name, &effect);
    }
    int status = -2;
    name = st_legacy_text(source);
    sort_population(&name, &decreasing, &criterion, &none, &status);
    if (status)
        error("MAB population sorting failed");
    char text[64], *number = text;
    if (all_best)
    {
        memset(text, 'a', 50);
        text[50] = '\0';
        int one = 1;
        double minimum = 0, maximum = 0;
        select_all_n_best(&source, &number, &one, &minimum, &maximum,
                          &criterion, &none, &status);
        if (status)
            error("MAB best-class selection failed");
    }
    else
        snprintf(text, sizeof(text), "%d", n);
    divide_population(&target, &source, &number);
}

static void st_mab_select(char *target, char *source, const char *effect_name,
                          int n, int all_best, int segments)
{
    st_mab_select_parameters(target, source, effect_name, n, all_best, segments,
                             1, !segments);
}

static void st_mab_fill(char *target, char *source, int required)
{
    int selected = (int)st_mab_population_size(target);
    if (selected < required)
    {
        char *temporary = "tmp-fillup";
        st_mab_select(temporary, source, "fm", required - selected, 0, 0);
        concat_population(&target, &temporary);
    }
}

static void st_mab_statistics(const char *population, const char *kind,
                              int allele, int chromosome, double begin,
                              double end, double position,
                              const char *effect_name, double summary[7])
{
    size_t n = st_mab_population_size(population);
    void *mark = st_work_mark();
    char *name = st_legacy_text(population);
    double *value = st_work_alloc(n ? n : 1, sizeof(double));
    memset(value, 0, n * sizeof(double));
    if (!strcmp(kind, "gvalue"))
    {
        genotype_population(&name);
        name = st_legacy_text(population);
        evaluate_population(&name, NULL);
        char *effect = (char *)effect_name, *criterion = "G";
        int status = -2;
        get_population_value(&name, value, &effect, &criterion, &status);
    }
    else if (!strcmp(kind, "drag"))
    {
        int chrom = chromosome + 1;
        evaluate_ld(&name, &allele, &chrom, &position, value);
    }
    else
    {
        double *length = st_work_alloc(n ? n : 1, sizeof(double));
        double *total = st_work_alloc(n ? n : 1, sizeof(double));
        int *blocks = st_work_alloc(n ? n : 1, sizeof(int));
        memset(length, 0, n * sizeof(double));
        memset(total, 0, n * sizeof(double));
        memset(blocks, 0, n * sizeof(int));
        if (chromosome < 0)
            evaluate_genome(&name, &allele, value, length, total, blocks);
        else
            evaluate_genome_region(&name, &allele, &chromosome, &begin, &end,
                                   value, length, total, blocks);
        if (!strcmp(kind, "segments"))
            for (size_t i = 0; i < n; i++)
                value[i] = blocks[i];
    }
    st_legacy_summary_values(value, n, summary);
    st_work_reset(mark);
}

static SEXP st_mab_summary(st_mab_generation *generations, int n,
                           const char *kind, int allele, int chromosome,
                           double begin, double end, double position,
                           const char *effect, int generation_labels)
{
    SEXP out = PROTECT(allocMatrix(REALSXP, n + 1, 7));
    SEXP rows = PROTECT(allocVector(STRSXP, n + 1));
    SEXP cols = PROTECT(allocVector(STRSXP, 7));
    SEXP dn = PROTECT(allocVector(VECSXP, 2));
    const char *names[] = {"Obs", "Mean", "SDev", "Min", "Q10", "Med", "Max"};
    for (int j = 0; j < 7; j++)
        SET_STRING_ELT(cols, j, mkChar(names[j]));
    for (int i = 0; i <= n; i++)
    {
        const char *population = i ? generations[i - 1].total : "P1";
        double values[7];
        st_mab_statistics(population, kind, allele, chromosome, begin, end,
                          position, effect, values);
        for (int j = 0; j < 7; j++)
            REAL(out)[i + (size_t)(n + 1) * j] = values[j];
        if (generation_labels && i)
        {
            char label[128];
            snprintf(label, sizeof(label), "G%d-%s", i - 1,
                     generations[i - 1].generation);
            SET_STRING_ELT(rows, i, mkChar(label));
        }
        else
            SET_STRING_ELT(rows, i, mkChar(population));
    }
    SET_VECTOR_ELT(dn, 0, rows);
    SET_VECTOR_ELT(dn, 1, cols);
    setAttrib(out, R_DimNamesSymbol, dn);
    UNPROTECT(4);
    return out;
}

static SEXP st_mab_marker_means(const double *totals,
                                st_mab_generation *generation, int n,
                                int repetitions)
{
    SEXP out = PROTECT(allocMatrix(REALSXP, n + 1, 1));
    SEXP rows = PROTECT(allocVector(STRSXP, n + 1));
    SEXP cols = PROTECT(mkString("Mean"));
    SEXP dn = PROTECT(allocVector(VECSXP, 2));
    SET_STRING_ELT(rows, 0, mkChar("Total"));
    long double total = 0;
    for (int i = 0; i < n; i++)
    {
        REAL(out)[i + 1] = totals[i] / repetitions;
        total += REAL(out)[i + 1];
        char name[128];
        snprintf(name, sizeof(name), "G%d-%s", i, generation[i].generation);
        SET_STRING_ELT(rows, i + 1, mkChar(name));
    }
    REAL(out)[0] = (double)total;
    SET_VECTOR_ELT(dn, 0, rows);
    SET_VECTOR_ELT(dn, 1, cols);
    setAttrib(out, R_DimNamesSymbol, dn);
    UNPROTECT(4);
    return out;
}

static void st_mab_output_add(SEXP output, SEXP names, int *index,
                              const char *key, SEXP value)
{
    SET_VECTOR_ELT(output, *index, value);
    SET_STRING_ELT(names, (*index)++, mkChar(key));
}

static SEXP st_legacy_mab(SEXP args)
{
    SEXP generation_arg = st_legacy_argument(args, "gen.type");
    SEXP sizes_arg = st_legacy_argument(args, "population.size");
    SEXP strategies_arg = st_legacy_argument(args, "sel.strategy");
    SEXP selected_arg = st_legacy_argument(args, "no.selected");
    SEXP preselected_arg = st_legacy_argument(args, "no.preselected");
    SEXP region_chr = st_legacy_argument(args, "reg.chr");
    SEXP region_begin = st_legacy_argument(args, "reg.begin");
    SEXP region_end = st_legacy_argument(args, "reg.end");
    int ng = LENGTH(generation_arg), nr = LENGTH(region_chr);
    if (LENGTH(sizes_arg) != ng || LENGTH(strategies_arg) != ng ||
        (!isNull(selected_arg) && LENGTH(selected_arg) != ng) ||
        (!isNull(preselected_arg) && LENGTH(preselected_arg) != ng))
    {
        info(-2, "Definition of generations incorrect\n");
        return R_NilValue;
    }
    if (LENGTH(region_begin) != nr || LENGTH(region_end) != nr)
    {
        info(-2, "Definition of target regions incorrect\n");
        return R_NilValue;
    }
    double reps = asReal(st_legacy_argument(args, "repetitions"));
    if (!R_FINITE(reps) || reps < 1 || reps > INT_MAX || reps != floor(reps))
        error("repetitions must be a positive integer");
    int repetitions = (int)reps;
    double success = asReal(st_legacy_argument(args, "success.factor"));
    if (!R_FINITE(success) || success < 0 || success > 1)
        error("success.factor must lie between zero and one");
    st_mab_generation *generations =
        st_work_alloc(ng ? ng : 1, sizeof(*generations));
    for (int i = 0; i < ng; i++)
    {
        st_mab_generation *g = &generations[i];
        if (TYPEOF(generation_arg) != STRSXP ||
            TYPEOF(strategies_arg) != STRSXP ||
            STRING_ELT(generation_arg, i) == NA_STRING ||
            STRING_ELT(strategies_arg, i) == NA_STRING)
            error(
                "Generation types and selection strategies must be character strings");
        g->generation = CHAR(STRING_ELT(generation_arg, i));
        g->strategy = CHAR(STRING_ELT(strategies_arg, i));
        if (!strcmp(g->strategy, "0"))
            g->strategy = "n";
        else if (!strcmp(g->strategy, "1"))
            g->strategy = "t";
        else if (!strcmp(g->strategy, "2"))
            g->strategy = "tb";
        else if (!strcmp(g->strategy, "3"))
            g->strategy = "tfb";
        double size = st_legacy_number(sizes_arg, i);
        double take =
            isNull(selected_arg) ? 1 : st_legacy_number(selected_arg, i);
        double pre =
            isNull(preselected_arg) ? -1 : st_legacy_number(preselected_arg, i);
        if (!R_FINITE(size) || size < 1 || size > INT_MAX ||
            size != floor(size) || !R_FINITE(take) || take < 1 ||
            take > INT_MAX || take != floor(take) || !R_FINITE(pre) ||
            pre < -1 || pre > INT_MAX || pre != floor(pre))
            error(
                "MAB population and selection sizes must be integers in range");
        g->size = (int)size;
        g->selected_n = (int)take;
        g->preselected_n = (int)pre;
        snprintf(g->population, sizeof(g->population), "GEN%d", i + 1);
        snprintf(g->first, sizeof(g->first), "GEN%ds1", i + 1);
        snprintf(g->second, sizeof(g->second), "GEN%ds2", i + 1);
        snprintf(g->third, sizeof(g->third), "GEN%ds3", i + 1);
        snprintf(g->selected, sizeof(g->selected), "GEN%dse", i + 1);
        snprintf(g->total, sizeof(g->total), "GEN%dst", i + 1);
    }
    const char *donor_name =
        CHAR(asChar(st_legacy_argument(args, "donor.parent")));
    const char *recurrent_name =
        CHAR(asChar(st_legacy_argument(args, "recurrent.parent")));
    int actual_parents = strcmp(donor_name, "H") && strcmp(recurrent_name, "H");
    if ((strcmp(donor_name, "H") != 0) != (strcmp(recurrent_name, "H") != 0))
    {
        info(-2, "Marker data must be supplied for both parents\n");
        return R_NilValue;
    }
    int nleft = 0, nright = 0;
    st_mab_parent_marker *left = NULL, *right = NULL;
    if (actual_parents)
    {
        left = st_mab_parent_read(st_legacy_argument(args, "parent1"), &nleft);
        right =
            st_mab_parent_read(st_legacy_argument(args, "parent2"), &nright);
    }
    const char *missing =
        CHAR(asChar(st_legacy_argument(args, "missing.allele")));
    SEXP map = st_legacy_argument(args, "map");
    if (TYPEOF(map) != VECSXP || LENGTH(map) != 3 ||
        TYPEOF(VECTOR_ELT(map, 2)) != STRSXP)
        error(
            "MAB linkage map must contain chromosome, position and marker columns");
    int map_n = LENGTH(VECTOR_ELT(map, 0));
    if (map_n < 1 || LENGTH(VECTOR_ELT(map, 1)) != map_n ||
        LENGTH(VECTOR_ELT(map, 2)) != map_n)
        error("Invalid MAB linkage map dimensions");
    st_map_extension *stage = calloc(1, sizeof(*stage));
    if (!stage)
        error("Cannot allocate MAB map staging");
    st_work_adopt(stage, st_map_extension_free);
    stage->n = map_n;
    stage->points = calloc(map_n, sizeof(TMappoint));
    if (!stage->points)
        error("Cannot allocate MAB linkage map");
    TMappoint *points = stage->points;
    int nmap = 0, nchrom = 0;
    for (int i = 0; i < map_n; i++)
    {
        if (STRING_ELT(VECTOR_ELT(map, 2), i) == NA_STRING)
            error("Missing MAB marker name");
        const char *name = CHAR(STRING_ELT(VECTOR_ELT(map, 2), i));
        if (actual_parents)
        {
            st_mab_parent_marker *a = st_mab_parent_find(left, nleft, name);
            st_mab_parent_marker *b = st_mab_parent_find(right, nright, name);
            if (!a || !b || !strcmp(a->first, missing) ||
                !strcmp(a->second, missing) || !strcmp(b->first, missing) ||
                !strcmp(b->second, missing))
                continue;
        }
        double chromosome = st_legacy_number(VECTOR_ELT(map, 0), i);
        double position = st_legacy_number(VECTOR_ELT(map, 1), i) / 100;
        if (!R_FINITE(chromosome) || chromosome < 1 || chromosome > INT_MAX ||
            chromosome != floor(chromosome) || !R_FINITE(position) ||
            position < 0)
            error("Invalid MAB map chromosome or position");
        TMappoint *point = &points[nmap++];
        point->chrom = (int)chromosome;
        point->Pos = position;
        point->counter = i;
        point->name = st_text_dup(name);
        point->class = malloc(4);
        if (!point->name || !point->class)
            error("Cannot allocate MAB marker names");
        strcpy(point->class, "m");
        if (point->chrom > nchrom)
            nchrom = point->chrom;
    }
    if (!nmap)
        error("No common nonmissing MAB markers");
    qsort(points, nmap, sizeof(*points), mccmp);
    SEXP target_loci = st_legacy_argument(args, "target.loci");
    SEXP flank_loci = st_legacy_argument(args, "flanking.loci");
    SEXP recipient_loci = st_legacy_argument(args, "recipient.loci");
    SEXP locus_lists[] = {target_loci, flank_loci, recipient_loci};
    const char *classes[] = {"t", "fm", "r"};
    for (int group = 0; group < 3; group++)
    {
        int found = 0;
        for (int i = 0; i < nmap; i++)
            if (st_mab_contains(locus_lists[group], points[i].name))
                strcpy(points[i].class, classes[group]), found = 1;
        if (LENGTH(locus_lists[group]) && !found)
        {
            info(-2,
                 "Map position or marker data for requested locus missing\n");
            return R_NilValue;
        }
    }
    double *lengths = st_work_alloc(nchrom, sizeof(double));
    int *chrom_present = st_work_alloc(nchrom, sizeof(int));
    memset(lengths, 0, nchrom * sizeof(double));
    memset(chrom_present, 0, nchrom * sizeof(int));
    for (int i = 0; i < nmap; i++)
    {
        int c = points[i].chrom - 1;
        lengths[c] = fmax(lengths[c], points[i].Pos);
        chrom_present[c] = 1;
    }
    for (int c = 0; c < nchrom; c++)
    {
        if (!chrom_present[c])
            error(
                "MAB map chromosomes must be numbered consecutively from one");
        lengths[c] += 0.001;
    }
    int homs = 2;
    if (!set_genome_par_impl(&nchrom, &homs, lengths))
        return R_NilValue;
    int modified = 0;
    if (!prepare_mappoints(points, nmap, 1, 100, &modified) ||
        !install_mappoints(points, nmap))
        error("Cannot install MAB linkage map");
    st_mab_uniform_effect("t", 1);
    st_mab_uniform_effect("fm", 8);
    st_mab_uniform_effect("r", 8);
    st_mab_uniform_effect("m", 8);
    char *p1 = "P1", *p2 = "P2";
    st_mab_init_homo(p1, 1);
    st_mab_init_homo(p2, 8);
    if (actual_parents)
    {
        char *both = "P1 P2";
        remove_population(&both);
        size_t count = (size_t)nmap * 2;
        int *ind = st_work_alloc(count, sizeof(int));
        int *chrom = st_work_alloc(count, sizeof(int));
        int *hom = st_work_alloc(count, sizeof(int));
        int *donor = st_work_alloc(count, sizeof(int));
        int *recipient = st_work_alloc(count, sizeof(int));
        double *position = st_work_alloc(count, sizeof(double));
        int *recode = st_work_alloc(nmap, sizeof(int));
        TMappoint *export = st_work_alloc(nmap, sizeof(TMappoint));
        int mi = 0;
        for (int c = 0; c < NoChroms; c++)
            for (int j = 0; j < PopDescr[c].NrMappoints; j++, mi++)
            {
                TMappoint *point = &PopDescr[c].Mappoint[j];
                st_mab_parent_marker *a =
                    st_mab_parent_find(left, nleft, point->name);
                st_mab_parent_marker *b =
                    st_mab_parent_find(right, nright, point->name);
                int shared = (!strcmp(a->first, b->first) ||
                              !strcmp(a->first, b->second)) +
                             (!strcmp(a->second, b->first) ||
                              !strcmp(a->second, b->second));
                if (st_mab_contains(target_loci, point->name))
                {
                    if (shared == 2)
                    {
                        info(-2, "Monomorphic target locus '%s'\n",
                             point->name);
                        return R_NilValue;
                    }
                    if (shared == 0 && strcmp(a->first, a->second))
                    {
                        shared = 1;
                        info(0, "Heterozygous target locus '%s'\n",
                             point->name);
                    }
                }
                recode[mi] = shared;
                export[mi] = *point;
                export[mi].counter = mi;
                for (int h = 0; h < 2; h++)
                {
                    size_t at = (size_t)mi * 2 + h;
                    ind[at] = 1;
                    chrom[at] = c + 1;
                    hom[at] = h + 1;
                    position[at] =
                        j ? (point->Pos + PopDescr[c].Mappoint[j - 1].Pos) / 2
                          : 0;
                    donor[at] = shared == 0 ? 1 : shared == 2 ? 8 : h ? 8 : 1;
                    recipient[at] = 8;
                }
            }
        char text[32], *number = text;
        snprintf(text, sizeof(text), "%zu", count);
        init_population(&p1, &number, ind, chrom, hom, position, donor);
        init_population(&p2, &number, ind, chrom, hom, position, recipient);
        SEXP recode_path = st_legacy_argument(args, "recode.path");
        if (!isNull(recode_path))
        {
            const char *path = CHAR(asChar(recode_path));
            FILE *file = fopen(path, "w");
            if (!file)
                error("Cannot open recoded parent output");
            SEXP labels = st_legacy_argument(args, "parent.labels");
            fprintf(file, "%s %s\n", CHAR(STRING_ELT(labels, 0)),
                    CHAR(STRING_ELT(labels, 1)));
            qsort(export, nmap, sizeof(*export), mnamecmp);
            for (int i = 0; i < nmap; i++)
            {
                int shared = recode[export[i].counter];
                fprintf(file, "%s.1  %d 0\n%s.8  %d 1\n", export[i].name,
                        shared != 2, export[i].name, shared != 0);
            }
            int failed = ferror(file);
            if (fclose(file) || failed)
                error("Cannot write recoded parent output");
        }
    }
    double *sm = st_work_alloc(ng ? ng : 1, sizeof(double));
    double *ht = st_work_alloc(ng ? ng : 1, sizeof(double));
    memset(sm, 0, ng * sizeof(double));
    memset(ht, 0, ng * sizeof(double));
    for (int i = 0; i < ng; i++)
    {
        char *name = generations[i].total;
        remove_population(&name);
    }
    for (int repetition = 0; repetition < repetitions && ng; repetition++)
    {
        reset_mdp();
        for (int i = 0; i < ng; i++)
        {
            void *mark = st_work_mark();
            st_mab_generation *g = &generations[i];
            char *population = g->population;
            char *previous = i ? generations[i - 1].selected : p1;
            char *first_parent =
                !i || !strcmp(g->generation, "f1") ? p1 : previous;
            char *second_parent =
                !strcmp(g->generation, "s") && i ? previous : p2;
            char count[32], *number = count;
            double significant = -1;
            if (!strcmp(g->generation, "dh") && i)
            {
                snprintf(count, sizeof(count), "%d", generations[i - 1].size);
                dh(&population, &previous, &number, &significant);
            }
            else
            {
                if (i && strcmp(g->generation, "f1") &&
                    strcmp(g->generation, "bc") && strcmp(g->generation, "s"))
                    error("Unsupported MAB generation type: %s", g->generation);
                snprintf(count, sizeof(count), "%d", g->size);
                double self = i && !strcmp(g->generation, "s") ? 1 : 0;
                int mode = 1, scheme = 0;
                cross(&population, &first_parent, &second_parent, &self,
                      &number, &significant, &mode, &scheme);
            }
            const char *strategy = g->strategy;
            if (!strcmp(strategy, "n"))
                st_mab_sample(g->selected, population, g->selected_n, 0);
            else
            {
                st_mab_select(g->first, population, "t", 1, 1, 0);
                if (!strcmp(strategy, "t"))
                    st_mab_sample(g->selected, g->first, g->selected_n, 1);
                else if (!strcmp(strategy, "tb") || !strcmp(strategy, "tr") ||
                         !strcmp(strategy, "ts"))
                {
                    ht[i] += st_mab_population_size(g->first);
                    st_mab_select(g->selected, g->first,
                                  !strcmp(strategy, "tr") ? "r" : "m",
                                  g->selected_n, 0, !strcmp(strategy, "ts"));
                }
                else if (!strcmp(strategy, "tf") || !strcmp(strategy, "tfb") ||
                         !strcmp(strategy, "tfrb") || !strcmp(strategy, "tfs"))
                {
                    reset_mdp();
                    char *candidate = g->first, *flanking = "fm";
                    genotype_population(&previous);
                    genotype_population(&p2);
                    evaluate_mdp(&candidate, &previous, &p2, &flanking);
                    sm[i] += MDP;
                    reset_mdp();
                    if (!strcmp(strategy, "tf"))
                        st_mab_select(g->selected, g->first, "fm",
                                      g->selected_n, 0, 0);
                    else
                    {
                        st_mab_select(g->second, g->first, "fm", 1, 1, 0);
                        if (g->preselected_n >= 0)
                            st_mab_fill(g->second, g->first, g->preselected_n);
                        ht[i] += st_mab_population_size(g->second);
                        char *candidate_final = g->second;
                        if (!strcmp(strategy, "tfrb"))
                        {
                            st_mab_select(g->third, g->second, "r", 1, 1, 0);
                            candidate_final = g->third;
                        }
                        st_mab_select(g->selected, candidate_final, "m",
                                      g->selected_n, 0,
                                      !strcmp(strategy, "tfs"));
                    }
                }
                else if (!strcmp(strategy, "trb") || !strcmp(strategy, "trfb"))
                {
                    ht[i] += st_mab_population_size(g->first);
                    st_mab_select(g->second, g->first, "r", 1, 1, 0);
                    char *candidate = g->second;
                    if (!strcmp(strategy, "trfb"))
                    {
                        st_mab_select(g->third, g->second, "fm", 1, 1, 0);
                        if (g->preselected_n >= 0)
                            st_mab_fill(g->third, g->second, g->preselected_n);
                        candidate = g->third;
                    }
                    st_mab_select(g->selected, candidate, "m", g->selected_n, 0,
                                  0);
                }
                else
                    error("Selection strategy '%s' not implemented", strategy);
            }
            if (success < 1)
            {
                GetRNGstate();
                double draw = rbinom(g->selected_n, success);
                PutRNGstate();
                int survivors = draw < 1 ? 1 : (int)draw;
                char *temporary = "temp", *chosen = g->selected;
                st_mab_sample(temporary, chosen, survivors, 0);
                swap_population_name(&temporary, &chosen);
            }
            char *total = g->total, *chosen = g->selected;
            append_population(&total, &chosen);
            st_work_reset(mark);
        }
        R_CheckUserInterrupt();
    }
    const char *metadata[] = {
        "simulation.name",  "simulation.run",  "repetitions",
        "recurrent.parent", "donor.parent",    "linkage.map",
        "target.loci",      "flanking.loci",   "recipient.loci",
        "gen.type",         "population.size", "sel.strategy",
        "no.selected",      "no.preselected",  "reg.chr",
        "reg.begin",        "reg.end",         "success.factor"};
    int metadata_count = 0;
    for (int i = 0; i < 18; i++)
        if (i == 12 || !isNull(st_legacy_argument(args, metadata[i])))
            metadata_count++;
    int nt = LENGTH(target_loci), target_summary = !isNull(target_loci);
    int outputs =
        metadata_count + target_summary + 2 + 2 * nr + nt + (ng ? 2 : 0);
    SEXP out = PROTECT(allocVector(VECSXP, outputs));
    SEXP names = PROTECT(allocVector(STRSXP, outputs));
    int index = 0;
    for (int i = 0; i < 18; i++)
    {
        SEXP value = st_legacy_argument(args, metadata[i]);
        if (i == 11)
        {
            SEXP canonical = PROTECT(allocVector(STRSXP, ng));
            for (int j = 0; j < ng; j++)
                SET_STRING_ELT(canonical, j, mkChar(generations[j].strategy));
            st_mab_output_add(out, names, &index, metadata[i], canonical);
            UNPROTECT(1);
        }
        else if (i == 12 && isNull(value))
        {
            SEXP counts = PROTECT(allocVector(REALSXP, ng));
            for (int j = 0; j < ng; j++)
                REAL(counts)[j] = generations[j].selected_n;
            st_mab_output_add(out, names, &index, metadata[i], counts);
            UNPROTECT(1);
        }
        else if (!isNull(value))
            st_mab_output_add(out, names, &index, metadata[i], value);
    }
    if (target_summary)
        st_mab_output_add(
            out, names, &index, "TargetAlleles",
            st_mab_summary(generations, ng, "gvalue", 0, -1, 0, 0, 0, "t", 1));
    st_mab_output_add(
        out, names, &index, "RPG.gw",
        st_mab_summary(generations, ng, "ratio", 8, -1, 0, 0, 0, "", 1));
    double *begins = st_work_alloc(nr ? nr : 1, sizeof(double));
    double *ends = st_work_alloc(nr ? nr : 1, sizeof(double));
    int *chromosomes = st_work_alloc(nr ? nr : 1, sizeof(int));
    for (int i = 0; i < nr; i++)
    {
        chromosomes[i] = (int)st_legacy_number(region_chr, i) - 1;
        if (chromosomes[i] < 0 || chromosomes[i] >= NoChroms)
            error("MAB region chromosome is out of range");
        begins[i] = st_legacy_number(region_begin, i);
        ends[i] = st_legacy_number(region_end, i);
        if (ends[i] == 0)
            ends[i] = PopDescr[chromosomes[i]].ChromLen;
    }
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nr; j++)
            begins[j] *= 100;
        char key[64];
        snprintf(key, sizeof(key), "RPG.reg%d", i + 1);
        st_mab_output_add(out, names, &index, key,
                          st_mab_summary(generations, ng, "ratio", 8,
                                         chromosomes[i], begins[i], ends[i], 0,
                                         "", 0));
    }
    st_mab_output_add(
        out, names, &index, "DonorSegments.gw",
        st_mab_summary(generations, ng, "segments", 1, -1, 0, 0, 0, "", 1));
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nr; j++)
            begins[j] *= 100;
        char key[64];
        snprintf(key, sizeof(key), "DonorSegments.reg%d", i + 1);
        st_mab_output_add(out, names, &index, key,
                          st_mab_summary(generations, ng, "segments", 1,
                                         chromosomes[i], begins[i], ends[i], 0,
                                         "", 0));
    }
    for (int i = 0; i < nt; i++)
    {
        TMappoint *point =
            find_mappoint((char *)CHAR(STRING_ELT(target_loci, i)));
        if (!point)
            error("MAB target locus not found");
        char key[64];
        snprintf(key, sizeof(key), "LinkageDrag%d", i + 1);
        st_mab_output_add(out, names, &index, key,
                          st_mab_summary(generations, ng, "drag", 1,
                                         point->chrom - 1, 0, 0, point->Pos, "",
                                         1));
    }
    if (ng)
    {
        st_mab_output_add(
            out, names, &index, "SM",
            st_mab_marker_means(sm, generations, ng, repetitions));
        st_mab_output_add(
            out, names, &index, "HT",
            st_mab_marker_means(ht, generations, ng, repetitions));
    }
    if (index != outputs)
        error("MAB result assembly mismatch");
    setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

static SEXP st_legacy_genome_evaluate(SEXP args)
{
    const char *population =
        CHAR(asChar(st_legacy_argument(args, "population")));
    TPopListEl *el = NULL;
    find_list_element((char *)population, &el);
    if (!el || el->Pop->NoInds > INT_MAX)
        error("Simulation population is unavailable or too large to export");
    int n = (int)el->Pop->NoInds,
        allele = asInteger(st_legacy_argument(args, "allele"));
    int verbose = asLogical(st_legacy_argument(args, "verbose")), np = 0;
    const char *keys[] = {"ratio", "length", "total", "noBlocks"};
    SEXP out = st_sim_list(verbose ? 4 : 1, keys, &np);
    for (int j = 0; j < (verbose ? 4 : 1); j++)
    {
        SEXP col = PROTECT(allocVector(j == 3 ? INTSXP : REALSXP, n));
        np++;
        SET_VECTOR_ELT(out, j, col);
        if (j == 3)
            memset(INTEGER(col), 0, n * sizeof(int));
        else
            memset(REAL(col), 0, n * sizeof(double));
    }
    double *ratio = REAL(VECTOR_ELT(out, 0));
    double *length = verbose ? REAL(VECTOR_ELT(out, 1))
                             : st_work_alloc(n ? n : 1, sizeof(double));
    double *total = verbose ? REAL(VECTOR_ELT(out, 2))
                            : st_work_alloc(n ? n : 1, sizeof(double));
    int *blocks = verbose ? INTEGER(VECTOR_ELT(out, 3))
                          : st_work_alloc(n ? n : 1, sizeof(int));
    if (!verbose)
    {
        memset(length, 0, n * sizeof(double));
        memset(total, 0, n * sizeof(double));
        memset(blocks, 0, n * sizeof(int));
    }
    char *name = st_legacy_text(population);
    SEXP chromosome = st_legacy_argument(args, "chromosome");
    if (isNull(chromosome))
        evaluate_genome(&name, &allele, ratio, length, total, blocks);
    else
    {
        int chrom = asInteger(chromosome) - 1;
        if (chrom < 0 || chrom >= NoChroms)
            error("Invalid chromosome number");
        SEXP begin_arg = st_legacy_argument(args, "begin"),
             end_arg = st_legacy_argument(args, "end");
        double begin = isNull(begin_arg) ? 0 : asReal(begin_arg);
        double end =
            isNull(end_arg) ? PopDescr[chrom].ChromLen : asReal(end_arg);
        evaluate_genome_region(&name, &allele, &chrom, &begin, &end, ratio,
                               length, total, blocks);
    }
    st_sim_frame(out, n, &np);
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_select_segments(SEXP args)
{
    char *target =
        st_legacy_text(CHAR(asChar(st_legacy_argument(args, "target"))));
    char *source =
        st_legacy_text(CHAR(asChar(st_legacy_argument(args, "source"))));
    double count = asReal(st_legacy_argument(args, "n"));
    int allele = asInteger(st_legacy_argument(args, "allele"));
    int decreasing = asLogical(st_legacy_argument(args, "decreasing"));
    if (!R_FINITE(count) || count < 1 || count > INT_MAX ||
        count != floor(count))
        error("n must be a positive integer");
    POPULATION size = st_mab_population_size(source);
    if (!size || count > size)
    {
        info(-2, "Population %s has too few individuals for selection\n",
             source);
        return R_NilValue;
    }
    st_mab_select_parameters(target, source, "", (int)count, 0, 1, allele,
                             decreasing);
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    char number[64];
    snprintf(number, sizeof(number), "%.0f", count);
    SET_VECTOR_ELT(out, 0, mkString(target));
    SET_VECTOR_ELT(out, 1, mkString(source));
    SET_VECTOR_ELT(out, 2, mkString(number));
    UNPROTECT(1);
    return out;
}

typedef struct
{
    TPop *founders;
    int families, progeny, total, dh, created, success;
    POPULATION *left, *right;
    int *family;
    char (*targets)[MAXNAME];
    TInd f1;
} st_cross_plan_execution;

static void st_cross_plan_execution_cleanup(void *payload)
{
    st_cross_plan_execution *context = payload;
    remove_individuum(&context->f1);
    if (!context->success)
        for (int i = 0; i < context->created; i++)
        {
            TPopListEl *entry = NULL;
            find_list_element(context->targets[i], &entry);
            if (entry) st_remove_population_entry(entry);
        }
}

static int st_cross_plan_gametes(TInd *offspring, TInd *left, TInd *right)
{
    double significant = -1;
    if (st_sim_prepare_individual(offspring, left, right, 1)) return -1;
    for (int chrom = 0; chrom < NoChroms; chrom++)
    {
        THom *hom = offspring->Chrom[chrom].Hom;
        hom[0].NoLociUsed = hom[1].NoLociUsed = 0;
        if (meiosis_impl(&hom[0], left->Chrom[chrom].Hom,
                         PopDescr[chrom].ChromLen, &significant, 0) ||
            meiosis_impl(&hom[1], right->Chrom[chrom].Hom,
                         PopDescr[chrom].ChromLen, &significant, 0))
            return -1;
    }
    return 0;
}

static int st_cross_plan_doubled_haploid(TInd *offspring, TInd *f1)
{
    double significant = -1;
    if (st_sim_prepare_individual(offspring, f1, f1, 1)) return -1;
    for (int chrom = 0; chrom < NoChroms; chrom++)
    {
        THom *hom = offspring->Chrom[chrom].Hom;
        hom[0].NoLociUsed = hom[1].NoLociUsed = 0;
        if (meiosis_impl(&hom[0], f1->Chrom[chrom].Hom,
                         PopDescr[chrom].ChromLen, &significant, 0) ||
            st_sim_hom_reserve(&hom[1], (size_t)hom[0].NoLociUsed))
            return -1;
        for (int loc = 0; loc < hom[0].NoLociUsed; loc++)
            if (insertLocus(hom[0].Loc[loc].All, hom[0].Loc[loc].Pos,
                             &hom[1]) == -1)
                return -1;
    }
    return 0;
}

static SEXP st_cross_plan_execution_body(void *payload)
{
    st_cross_plan_execution *context = payload;
    const char *summary_names[] = {"family", "population", "P1Name",
        "P2Name", "n.progeny", "scheme"};
    const char *design_names[] = {"family", "population", "name", "P1Name",
        "P2Name", "F1", "pedigree", "scheme"};
    int np = 0, row = 0;
    SEXP summary = st_sim_list(6, summary_names, &np),
         design = st_sim_list(8, design_names, &np);
    for (int j = 0; j < 6; j++)
        SET_VECTOR_ELT(summary, j, allocVector(j == 0 || j == 4 ? INTSXP :
            STRSXP, context->families));
    for (int j = 0; j < 8; j++)
        SET_VECTOR_ELT(design, j, allocVector(j == 0 ? INTSXP : STRSXP,
            context->total));
    st_sim_frame(summary, context->families, &np);
    st_sim_frame(design, context->total, &np);
    SEXP scheme = PROTECT(mkChar(context->dh ? "DH" : "F1")); np++;
    for (int i = 0; i < context->families; i++)
    {
        R_CheckUserInterrupt();
        TInd *left = &context->founders->Ind[context->left[i]],
             *right = &context->founders->Ind[context->right[i]];
        uint64_t first, f1_serial = 0;
        if ((context->dh && st_sim_reserve(1, &f1_serial)) ||
            st_sim_reserve((POPULATION)context->progeny, &first))
            error("Simulation individual identities are exhausted");
        if (context->dh)
        {
            char f1_name[MAXNAME];
            size_t length = strlen(context->targets[i]);
            memcpy(f1_name, context->targets[i], length);
            memcpy(f1_name + length, "F1", 3);
            if (new_individuum(&context->f1) ||
                st_cross_plan_gametes(&context->f1, left, right) ||
                st_sim_identity_new(&context->f1, f1_name, f1_serial, left,
                    right, NULL))
                error("The planned F1 could not be generated");
        }
        char *name = context->targets[i];
        POPULATION count = (POPULATION)context->progeny;
        TPop *offspring = new_population(&name, &count);
        if (!offspring)
            error("A planned progeny population could not be allocated");
        context->created = i + 1;
        offspring->NoInds = count;
        for (int k = 0; k < context->progeny; k++, row++)
        {
            if (!(k % 128)) R_CheckUserInterrupt();
            TInd *individual = &offspring->Ind[k];
            int failed = context->dh
                             ? st_cross_plan_doubled_haploid(individual,
                                 &context->f1)
                             : st_cross_plan_gametes(individual, left, right);
            if (failed ||
                st_sim_identity_new(individual, offspring->Name,
                                    first + (uint64_t)k,
                                    context->dh ? &context->f1 : left,
                                    context->dh ? NULL : right,
                                    context->dh ? "DH" : NULL))
                error("The planned offspring could not be generated");
            INTEGER(VECTOR_ELT(design, 0))[row] = context->family[i];
            SET_STRING_ELT(VECTOR_ELT(design, 1), row, mkChar(offspring->Name));
            SET_STRING_ELT(VECTOR_ELT(design, 2), row,
                mkCharCE(individual->name, CE_UTF8));
            SET_STRING_ELT(VECTOR_ELT(design, 3), row, mkCharCE(left->name,
                CE_UTF8));
            SET_STRING_ELT(VECTOR_ELT(design, 4), row, mkCharCE(right->name,
                CE_UTF8));
            SET_STRING_ELT(VECTOR_ELT(design, 5), row,
                           context->dh ? mkCharCE(context->f1.name, CE_UTF8)
                               : NA_STRING);
            SET_STRING_ELT(VECTOR_ELT(design, 6), row,
                mkCharCE(individual->pedigree, CE_UTF8));
            SET_STRING_ELT(VECTOR_ELT(design, 7), row, scheme);
        }
        INTEGER(VECTOR_ELT(summary, 0))[i] = context->family[i];
        SET_STRING_ELT(VECTOR_ELT(summary, 1), i, mkChar(offspring->Name));
        SET_STRING_ELT(VECTOR_ELT(summary, 2), i, mkCharCE(left->name,
            CE_UTF8));
        SET_STRING_ELT(VECTOR_ELT(summary, 3), i, mkCharCE(right->name,
            CE_UTF8));
        INTEGER(VECTOR_ELT(summary, 4))[i] = context->progeny;
        SET_STRING_ELT(VECTOR_ELT(summary, 5), i, scheme);
        remove_individuum(&context->f1);
    }
    setAttrib(summary, install("design"), design);
    context->success = 1;
    UNPROTECT(np);
    return summary;
}

static const char *st_cross_plan_execution_string(SEXP input, const char *name)
{
    if (TYPEOF(input) != STRSXP || XLENGTH(input) != 1 ||
        STRING_ELT(input, 0) == NA_STRING || !*CHAR(STRING_ELT(input, 0)))
        error("%s must be one nonmissing name", name);
    return CHAR(STRING_ELT(input, 0));
}

static void st_cross_plan_founder_check(TInd *individual)
{
    if (!individual->Storage || individual->Storage->NoChroms !=
        (size_t)NoChroms ||
        individual->Storage->NoHoms != 2 || !individual->name)
        error("Selected founders must have the current diploid simulation"
            " genome");
    for (int chrom = 0; chrom < NoChroms; chrom++)
    {
        THom *hom = individual->Chrom[chrom].Hom;
        for (int h = 0; h < 2; h++)
        {
            if (!hom[h].Loc || hom[h].NoLociUsed < 1)
                error("Selected founders have incomplete simulation"
                    " haplotypes");
            for (int loc = 0; loc < hom[h].NoLociUsed; loc++)
                if (!R_FINITE(hom[h].Loc[loc].Pos) || hom[h].Loc[loc].Pos < 0 ||
                    hom[h].Loc[loc].Pos > PopDescr[chrom].ChromLen ||
                    (loc && hom[h].Loc[loc - 1].Pos >= hom[h].Loc[loc].Pos))
                    error("Selected founders have invalid simulation"
                        " haplotypes");
        }
        for (int marker = 0; marker < PopDescr[chrom].NrMappoints; marker++)
        {
            double pos = PopDescr[chrom].Mappoint[marker].Pos;
            unsigned short a = get_allele(&hom[0], pos),
                           b = get_allele(&hom[1], pos);
            if (!a || !b || a == USHRT_MAX || b == USHRT_MAX)
                error("Cross-plan execution requires complete selected"
                    " founders at every marker");
            if (a != b)
                error("Cross-plan execution currently requires homozygous"
                    " selected founders");
        }
    }
}

static SEXP st_legacy_cross_plan_execute(SEXP args)
{
    SEXP plan = st_legacy_argument(args, "plan"),
         crosses = st_legacy_argument(plan, "crosses"),
         left = st_legacy_argument(crosses, "P1Name"),
         right = st_legacy_argument(crosses, "P2Name"),
         sizes = st_legacy_argument(crosses, "n.progeny"),
         family = st_legacy_argument(crosses, "family");
    const char *population =
        st_cross_plan_execution_string(st_legacy_argument(args,
        "population"), "population"),
               *prefix =
                   st_cross_plan_execution_string(st_legacy_argument(args,
                   "name"), "name"),
               *scheme =
                   st_cross_plan_execution_string(st_legacy_argument(args,
                   "scheme"), "scheme");
    if (strcmp(scheme, "DH") && strcmp(scheme, "F1"))
        error("scheme must be DH or F1");
    size_t prefix_length = strlen(prefix);
    for (size_t i = 0; i < prefix_length; i++)
        if (!((prefix[i] >= 'A' && prefix[i] <= 'Z') ||
              (prefix[i] >= 'a' && prefix[i] <= 'z') ||
              (prefix[i] >= '0' && prefix[i] <= '9')))
            error("name must contain only ASCII letters and digits");
    if (TYPEOF(left) != STRSXP || TYPEOF(right) != STRSXP ||
        XLENGTH(left) < 1 || XLENGTH(left) > INT_MAX ||
        XLENGTH(left) != XLENGTH(right) || XLENGTH(left) != XLENGTH(sizes) ||
        XLENGTH(left) != XLENGTH(family) ||
        (TYPEOF(sizes) != INTSXP && TYPEOF(sizes) != REALSXP) ||
        (TYPEOF(family) != INTSXP && TYPEOF(family) != REALSXP))
        error("plan$crosses requires equal-length family, P1Name, P2Name and"
            " n.progeny columns");
    if (inherits(sizes, "factor") || inherits(family, "factor"))
        error("family and n.progeny must be numeric without factor levels");
    st_cross_plan_execution context;
    memset(&context, 0, sizeof(context));
    context.families = (int)XLENGTH(left);
    context.dh = !strcmp(scheme, "DH");
    char suffix[32];
    int suffix_length = snprintf(suffix, sizeof(suffix), "%03d",
        context.families);
    if (prefix_length + (size_t)suffix_length + 2 >= MAXNAME)
        error("name is too long for the planned family populations");
    TPopListEl *entry = NULL;
    find_list_element((char *)population, &entry);
    if (!entry || entry->Pop->NoInds < 2)
        error("population must name an existing founder population with at"
            " least two individuals");
    context.founders = entry->Pop;
    if (NoHoms != 2 || NoChroms < 1 || !PopDescr || cMP < 1)
        error("Cross-plan execution requires an installed diploid simulation"
            " map");
    if (!rng_initialized || !rng[0])
        error("The simulation random number generator must be initialized"
            " before cross-plan execution");
    if (!R_FINITE(COFreq) || COFreq < 0)
        error("The simulation crossover frequency must be finite and"
            " nonnegative");
    int *seen = st_work_alloc((size_t)cMP, sizeof(*seen));
    memset(seen, 0, (size_t)cMP * sizeof(*seen));
    int map_markers = 0;
    for (int chrom = 0; chrom < NoChroms; chrom++)
    {
        if (!R_FINITE(PopDescr[chrom].ChromLen) ||
            PopDescr[chrom].ChromLen <= 0 ||
            PopDescr[chrom].NrMappoints < 1 || !PopDescr[chrom].Mappoint)
            error("Cross-plan execution requires a valid simulation map for"
                " every chromosome");
        for (int marker = 0; marker < PopDescr[chrom].NrMappoints; marker++)
        {
            TMappoint *point = &PopDescr[chrom].Mappoint[marker];
            if (!R_FINITE(point->Pos) || point->Pos < 0 ||
                point->Pos > PopDescr[chrom].ChromLen ||
                point->chrom != chrom + 1 ||
                point->counter < 0 || point->counter >= cMP ||
                    seen[point->counter] ||
                (marker && PopDescr[chrom].Mappoint[marker - 1].Pos >
                    point->Pos))
                error("Cross-plan execution requires a valid installed"
                    " simulation marker map");
            seen[point->counter] = 1;
            if (map_markers == INT_MAX) error("The simulation marker map is"
                " too large");
            map_markers++;
        }
    }
    if (map_markers != cMP) error("The simulation marker map is incomplete");
    context.left = st_work_alloc((size_t)context.families,
        sizeof(*context.left));
    context.right = st_work_alloc((size_t)context.families,
        sizeof(*context.right));
    context.family = st_work_alloc((size_t)context.families,
        sizeof(*context.family));
    context.targets = st_work_alloc((size_t)context.families,
        sizeof(*context.targets));
    for (int i = 0; i < context.families; i++)
    {
        if (STRING_ELT(left, i) == NA_STRING || STRING_ELT(right, i) ==
            NA_STRING)
            error("Planned parent IDs must be nonmissing");
        const char *a = CHAR(STRING_ELT(left, i)), *b =
            CHAR(STRING_ELT(right, i));
        double count = st_legacy_number(sizes, i), f =
            st_legacy_number(family, i);
        if (!R_FINITE(count) || count < 1 || count > INT_MAX || count !=
            floor(count) ||
            (i && count != context.progeny))
            error("All planned families require the same positive whole"
                " n.progeny");
        if (!R_FINITE(f) || f < 1 || f > INT_MAX || f != floor(f))
            error("Planned family IDs must be positive whole numbers");
        context.progeny = (int)count;
        context.family[i] = (int)f;
        context.left[i] = context.right[i] = -1;
        for (POPULATION j = 0; j < context.founders->NoInds; j++)
        {
            const char *id = context.founders->Ind[j].name;
            if (id && !strcmp(id, a))
            {
                if (context.left[i] >= 0) error("A selected founder ID is"
                    " ambiguous");
                context.left[i] = j;
            }
            if (id && !strcmp(id, b))
            {
                if (context.right[i] >= 0) error("A selected founder ID is"
                    " ambiguous");
                context.right[i] = j;
            }
        }
        if (context.left[i] < 0 || context.right[i] < 0)
            error("Every planned parent ID must identify a founder in"
                " population");
        if (context.left[i] == context.right[i]) error("Self crosses are not"
            " supported");
        for (int j = 0; j < i; j++)
            if (context.family[i] == context.family[j] ||
                ((context.left[i] == context.left[j] && context.right[i] ==
                    context.right[j]) ||
                 (context.left[i] == context.right[j] && context.right[i] ==
                     context.left[j])))
                error("Planned families and unordered parent pairs must be"
                    " distinct");
        st_cross_plan_founder_check(&context.founders->Ind[context.left[i]]);
        st_cross_plan_founder_check(&context.founders->Ind[context.right[i]]);
        snprintf(context.targets[i], MAXNAME, "%s%03d", prefix, i + 1);
        if (ps_population_exists(context.targets[i]))
            error("Planned target population '%s' already exists",
                context.targets[i]);
    }
    if (context.families > INT_MAX / context.progeny)
        error("The total planned progeny exceeds the supported design size");
    context.total = context.families * context.progeny;
    uint64_t identities = (uint64_t)context.total + (context.dh ?
        (uint64_t)context.families : 0);
    if (identities > UINT64_MAX - st_identity_serial)
        error("The plan exceeds the available simulation individual"
            " identities");
    return R_ExecWithCleanup(st_cross_plan_execution_body, &context,
                             st_cross_plan_execution_cleanup, &context);
}


static SEXP st_legacy_single_cross(SEXP args)
{
    char *target =
        st_legacy_text(CHAR(asChar(st_legacy_argument(args, "target"))));
    char *source =
        st_legacy_text(CHAR(asChar(st_legacy_argument(args, "source"))));
    double size_arg = asReal(st_legacy_argument(args, "size"));
    if (!R_FINITE(size_arg) || size_arg < 1 || size_arg != floor(size_arg) ||
        size_arg > INT_MAX)
        error("classSize must be a positive integer");
    POPULATION size = st_mab_population_size(source),
               class_size = (POPULATION)size_arg;
    POPULATION groups = size / class_size;
    if (groups < 2)
        error("At least two classes are required for single.cross");
    remove_population(&target);
    char *left = "single-cross-a", *right = "single-cross-b",
         *temporary = "single-cross-tmp";
    for (POPULATION i = 1; i < groups; i++)
    {
        char start_text[64], count_text[64], *start = start_text,
                                             *count = count_text;
        snprintf(start_text, sizeof(start_text), "%llu",
                 (unsigned long long)(i * class_size + 1));
        snprintf(count_text, sizeof(count_text), "-1");
        copy_population(&right, &source, &start, &count);
        POPULATION right_size = st_mab_population_size(right);
        for (POPULATION j = 1; j <= class_size; j++)
        {
            snprintf(start_text, sizeof(start_text), "%llu",
                     (unsigned long long)((i - 1) * class_size + j));
            snprintf(count_text, sizeof(count_text), "1");
            copy_population(&left, &source, &start, &count);
            snprintf(count_text, sizeof(count_text), "%llu",
                     (unsigned long long)right_size);
            double self = 0, significant = -1;
            int mode = 1, scheme = 3;
            cross(&temporary, &left, &right, &self, &count, &significant, &mode,
                  &scheme);
            concat_population(&target, &temporary);
        }
        remove_population(&right);
        R_CheckUserInterrupt();
    }
    remove_population(&left);
    SEXP out = PROTECT(allocVector(VECSXP, 1));
    SET_VECTOR_ELT(out, 0, mkString(left));
    UNPROTECT(1);
    return out;
}

static SEXP st_legacy_mdp(SEXP args)
{
    const char *pg = CHAR(asChar(st_legacy_argument(args, "progeny")));
    const char *p1 = CHAR(asChar(st_legacy_argument(args, "parent1")));
    const char *p2 = CHAR(asChar(st_legacy_argument(args, "parent2")));
    SEXP effect_arg = st_legacy_argument(args, "effect");
    char *progeny = st_legacy_text(pg), *parent1 = st_legacy_text(p1),
         *parent2 = st_legacy_text(p2);
    genotype_population(&parent1);
    genotype_population(&parent2);
    if (isNull(effect_arg))
    {
        for (int i = 0; i < NrEffects; i++)
            if (Effects[i].NrEff)
            {
                char *effect = Effects[i].Name;
                parent1 = st_legacy_text(p1);
                parent2 = st_legacy_text(p2);
                genotype_population(&parent1);
                genotype_population(&parent2);
                evaluate_mdp(&progeny, &parent1, &parent2, &effect);
            }
        return R_NilValue;
    }
    char *effect = st_legacy_text(CHAR(asChar(effect_arg)));
    evaluate_mdp(&progeny, &parent1, &parent2, &effect);
    SEXP out = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out, 0, mkString(progeny));
    SET_VECTOR_ELT(out, 1, mkString(parent1));
    SET_VECTOR_ELT(out, 2, mkString(parent2));
    SET_VECTOR_ELT(out, 3, mkString(effect));
    UNPROTECT(1);
    return out;
}

typedef struct
{
    int mask, count;
    const char *label[2];
} st_legacy_target;

static st_legacy_target st_legacy_target_classify(SEXP criteria)
{
    st_legacy_target target = {0, 0, {NULL, NULL}};
    if (TYPEOF(criteria) != VECSXP || LENGTH(criteria) < 3)
        error(
            "criteria must be a data frame with locus and two allele columns");
    SEXP a = VECTOR_ELT(criteria, 1), b = VECTOR_ELT(criteria, 2);
    int n = LENGTH(VECTOR_ELT(criteria, 0));
    if (LENGTH(a) != n || LENGTH(b) != n)
        error("Selection criterion columns must have equal lengths");
    if (n == 3)
    {
        target.mask = 7;
        target.count = 1;
        target.label[0] = "All";
        return target;
    }
    if (n < 1 || n > 3)
    {
        info(-2, "More than three possible selection criteria\n");
        return target;
    }
    double first = st_legacy_number(a, 0) + st_legacy_number(b, 0);
    double second =
        n == 2 ? st_legacy_number(a, 1) + st_legacy_number(b, 1) : first;
    if (first == second)
    {
        if (first == 2)
            target.mask = 1, target.label[0] = "H1";
        else if (first == 3)
            target.mask = 2, target.label[0] = "MIX";
        else if (first == 4)
            target.mask = 4, target.label[0] = "H2";
        if (target.mask)
            target.count = 1;
    }
    else
    {
        double total = first + second;
        if (total == 5)
            target.mask = 3, target.label[0] = "H1", target.label[1] = "MIX";
        else if (total == 6)
            target.mask = 5, target.label[0] = "H2", target.label[1] = "H1";
        else if (total == 7)
            target.mask = 3, target.label[0] = "MIX", target.label[1] = "H1";
        if (target.mask)
            target.count = 2;
    }
    if (!target.mask)
        info(-2, "Faulty selection criterion\n");
    return target;
}

static SEXP st_legacy_target_export(SEXP args)
{
    st_legacy_target target =
        st_legacy_target_classify(st_legacy_argument(args, "criteria"));
    if (!target.count)
        return R_NilValue;
    SEXP out = PROTECT(allocVector(STRSXP, target.count));
    for (int i = 0; i < target.count; i++)
        SET_STRING_ELT(out, i, mkChar(target.label[i]));
    UNPROTECT(1);
    return out;
}

static void st_legacy_sort(char *population, int decreasing)
{
    char *criterion = "G", *none = "";
    int status = -2;
    sort_population(&population, &decreasing, &criterion, &none, &status);
    if (status)
        error("Cannot sort genotype-selection population");
}

static void st_legacy_best_number(char *population, int classes,
                                  char number[64])
{
    memset(number, 'a', 50);
    number[50] = '\0';
    char *count = number, *criterion = "G", *none = "";
    double minimum = 0, maximum = 0;
    int status = -2;
    select_all_n_best(&population, &count, &classes, &minimum, &maximum,
                      &criterion, &none, &status);
    if (status)
        error("Cannot determine genotype-selection class size");
}

static void st_legacy_best_divide(char *target, char *source, int decreasing,
                                  int classes)
{
    st_legacy_sort(source, decreasing);
    char number[64], *count = number;
    st_legacy_best_number(source, classes, number);
    divide_population(&target, &source, &count);
}

static SEXP st_legacy_select_genotype(SEXP args)
{
    st_legacy_target criterion =
        st_legacy_target_classify(st_legacy_argument(args, "criteria"));
    if (!criterion.mask)
        error("No valid genotype selection criterion");
    char *source =
        st_legacy_text(CHAR(asChar(st_legacy_argument(args, "source"))));
    char *target =
        st_legacy_text(CHAR(asChar(st_legacy_argument(args, "target"))));
    char *destination =
        st_legacy_text(CHAR(asChar(st_legacy_argument(args, "destination"))));
    char *effect =
        st_legacy_text(CHAR(asChar(st_legacy_argument(args, "effect"))));
    genotype_population(&source);
    evaluate_population(&source, &effect);
    if (criterion.mask == 7)
        st_legacy_best_divide(target, source, 1, 3);
    else if (criterion.mask == 3)
        st_legacy_best_divide(target, source, 1, 2);
    else if (criterion.mask == 5)
    {
        char *temporary = "tmp";
        st_legacy_best_divide(target, source, 1, 1);
        st_legacy_best_divide(temporary, source, 0, 1);
        concat_population(&target, &temporary);
    }
    else if (criterion.mask == 1)
        st_legacy_best_divide(target, source, 1, 1);
    else if (criterion.mask == 4)
        st_legacy_best_divide(target, source, 0, 1);
    else
    {
        char *temporary = "tempPop", number[64], *count = number;
        st_legacy_best_divide(temporary, source, 1, 2);
        st_legacy_sort(temporary, 0);
        st_legacy_best_number(source, 1, number);
        divide_population(&target, &temporary, &count);
        concat_population(&source, &temporary);
    }
    char *start = "1", *count = "-1";
    copy_population(&destination, &target, &start, &count);
    return R_NilValue;
}

static SEXP st_legacy_locus_effect(SEXP args)
{
    const char *name = CHAR(asChar(st_legacy_argument(args, "name")));
    const char *locus = CHAR(asChar(st_legacy_argument(args, "locus")));
    int allele = 1;
    double value = 1;
    st_define_effects_native(name, &locus, &allele, &value, 1);
    return R_NilValue;
}
