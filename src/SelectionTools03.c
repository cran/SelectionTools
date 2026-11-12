

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <ctype.h>
#ifdef _OPENMP
#ifndef ST
#define MT
#include <omp.h>
#endif
#endif

#include <stdint.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include "SelectionToolsAllele.h"

extern gs_varset_type *gs_fdta(char *set_name);
extern void gs_info(gs_varset_type *gv, int level, char *message);
extern int gs_string_cmp(const void *a, const void *b);

static int gs_size_mul(size_t a, size_t b, size_t *out)
{
    if ((0 != a) && (b > ((size_t)-1) / a))
        return -1;
    *out = a * b;
    return 0;
}

static int gs_map_pool_valid(const gs_mappoint_type *map_S, gs_UNSIGNED nloc,
                             const char *pool, size_t pool_n)
{
    size_t offset = 0;

    if (0 == nloc)
        return ((NULL == map_S) && (NULL == pool) && (0 == pool_n));
    if ((NULL == map_S) || (NULL == pool) || (0 == pool_n))
        return 0;

    for (gs_UNSIGNED k = 0; k < nloc; k++)
    {
        const char *end;
        size_t len;

        if ((offset >= pool_n) || (map_S[k].name != pool + offset))
            return 0;
        end = memchr(pool + offset, '\0', pool_n - offset);
        if (NULL == end)
            return 0;
        len = (size_t)(end - (pool + offset));
        if (len >= gs_MAXSTR)
            return 0;
        offset += len + 1;

        if ((offset >= pool_n) || (map_S[k].class != pool + offset))
            return 0;
        end = memchr(pool + offset, '\0', pool_n - offset);
        if (NULL == end)
            return 0;
        len = (size_t)(end - (pool + offset));
        if (len >= gs_MAXSTR)
            return 0;
        offset += len + 1;
    }

    return (offset == pool_n);
}

static int gs_marker_genotypes_valid(const gs_varset_type *gv)
{
    size_t n;

    if ((NULL == gv) || !st_allele_present(gv) || (0 == gv->NoInd) ||
        (0 == gv->NoMar) ||
        (0 != gs_size_mul((size_t)gv->NoInd, (size_t)gv->NoMar, &n)))
        return 0;

    for (size_t idx = 0; idx < n; idx++)
    {
        GS_ALLELE_TYPE a = st_allele_get(gv, idx, 0);
        GS_ALLELE_TYPE b = st_allele_get(gv, idx, 1);
        int ma = (GS_ALLELE_MISSING == a);
        int mb = (GS_ALLELE_MISSING == b);
        if (ma != mb)
            return 0;
        if (!ma && ((0 == a) || (0 == b)))
            return 0;
    }
    return 1;
}

static gs_UNSIGNED *gs_map_storage_index(gs_varset_type *gv,
                                         int require_complete,
                                         int require_unique_names)
{
    gs_UNSIGNED *storage = NULL;
    unsigned char *seen = NULL;
    char **names = NULL;
    gs_UNSIGNED k, lne, m;
    int ok = 0;

    if ((0 == gv->NoMar) || (NULL == gv->MarName))
    {
        gs_info(gv, -2, "No marker data available");
        return NULL;
    }
    if ((0 == gv->NLoc) || (NULL == gv->map_S) || (NULL == gv->map_str_S) ||
        (0 == gv->map_str_N) || (NULL == gv->map) || (NULL == gv->mdl))
    {
        gs_info(gv, -2, "No linkage map defined");
        return NULL;
    }
    if (gv->NLoc > gv->NoMar)
    {
        gs_info(gv, -2, "More linkage-map loci than marker-data rows");
        return NULL;
    }
    if (require_complete && (gv->NLoc != gv->NoMar))
    {
        gs_info(gv, -2, "Linkage map is not matching the marker data");
        return NULL;
    }
    if (gv->NoMar > (gs_UNSIGNED)(((size_t)-1) / sizeof(*storage)))
    {
        gs_info(gv, -2, "Memory allocation error");
        return NULL;
    }

    storage = malloc((size_t)gv->NoMar * sizeof(*storage));
    seen = calloc((size_t)gv->NoMar, sizeof(*seen));
    if ((NULL == storage) || (NULL == seen))
    {
        gs_info(gv, -2, "Memory allocation error");
        goto cleanup;
    }
    for (m = 0; m < gv->NoMar; m++)
        storage[m] = gs_UNSIGNED_MAX;

    if (!gs_map_pool_valid(gv->map_S, gv->NLoc, gv->map_str_S, gv->map_str_N))
    {
        gs_info(gv, -2, "Invalid linkage map name or class storage");
        goto cleanup;
    }

    for (k = 0; k < gv->NLoc; k++)
    {
        m = gv->map_S[k].nr;
        if ((m >= gv->NoMar) || (gs_UNSIGNED_MAX != storage[m]))
        {
            gs_info(gv, -2, "Invalid linkage-map data-row index");
            goto cleanup;
        }
        storage[m] = k;
    }

    for (lne = 0; lne < gv->NLoc; lne++)
    {
        if (NULL == gv->map[lne])
        {
            gs_info(gv, -2, "Invalid linkage map pointer");
            goto cleanup;
        }
        m = gv->mdl[lne];
        if ((m >= gv->NoMar) || (gs_UNSIGNED_MAX == storage[m]))
        {
            gs_info(gv, -2, "Invalid map to marker-data index");
            goto cleanup;
        }
        if (gv->map[lne] != &(gv->map_S[storage[m]]))
        {
            gs_info(gv, -2, "Linkage map pointer does not match map storage");
            goto cleanup;
        }
        if (gv->map[lne]->nr != m)
        {
            gs_info(gv, -2, "Inconsistent map to marker-data index");
            goto cleanup;
        }
        if (NULL == memchr(gv->MarName[m], '\0', gs_MAXNME))
        {
            gs_info(gv, -2, "Invalid marker name");
            goto cleanup;
        }
        if (0 != strcmp(gv->map[lne]->name, gv->MarName[m]))
        {
            gs_info(gv, -2, "Linkage-map marker is not matching marker data");
            goto cleanup;
        }
        if (seen[m])
        {
            gs_info(gv, -2,
                    "Marker-data row used more than once by linkage map");
            goto cleanup;
        }
        seen[m] = 1;
    }

    if (require_complete)
        for (m = 0; m < gv->NoMar; m++)
            if (!seen[m])
            {
                gs_info(gv, -2,
                        "Marker data are not completely covered by map");
                goto cleanup;
            }

    if (require_unique_names)
    {
        if (gv->NoMar > (gs_UNSIGNED)(((size_t)-1) / sizeof(*names)))
        {
            gs_info(gv, -2, "Memory allocation error");
            goto cleanup;
        }
        names = malloc((size_t)gv->NoMar * sizeof(*names));
        if (NULL == names)
        {
            gs_info(gv, -2, "Memory allocation error");
            goto cleanup;
        }
        for (m = 0; m < gv->NoMar; m++)
            names[m] = gv->MarName[m];
        qsort(names, gv->NoMar, sizeof(*names), gs_string_cmp);
        for (m = 1; m < gv->NoMar; m++)
            if (0 == strcmp(names[m - 1], names[m]))
            {
                snprintf(gv->msg, sizeof(gv->msg),
                         "Duplicate marker name %s in marker data", names[m]);
                gs_info(gv, -2, gv->msg);
                goto cleanup;
            }
    }

    ok = 1;
cleanup:
    free(names);
    free(seen);
    if (!ok)
    {
        free(storage);
        storage = NULL;
    }
    return storage;
}

typedef struct
{
    const char *name;
    size_t index;
} st_phase_r_name_index;

static int st_phase_r_name_index_cmp(const void *a, const void *b)
{
    const st_phase_r_name_index *x = (const st_phase_r_name_index *)a;
    const st_phase_r_name_index *y = (const st_phase_r_name_index *)b;
    int cmp = strcmp(x->name, y->name);
    if (0 != cmp)
        return cmp;
    return (x->index > y->index) - (x->index < y->index);
}

extern gs_varset_type *gv;
extern gs_varset_type **gs_data;
extern gs_UNSIGNED gs_ndta;

static gs_varset_type *st_phase_step1_find_dataset(const char *name)
{
    if ((NULL == name) || ('\0' == name[0]))
        return NULL;
    if ((NULL != gv) && (0 == strcmp(name, "default")))
        return gv;
    for (gs_UNSIGNED i = 0u; i < gs_ndta; i++)
        if ((NULL != gs_data) && (NULL != gs_data[i]) &&
            (0 == strcmp(name, gs_data[i]->name)))
            return gs_data[i];
    return NULL;
}

static int st_phase_step1_dataset_name(SEXP x, const char **name)
{
    if ((NULL == name) || (STRSXP != TYPEOF(x)) || (1 != XLENGTH(x)) ||
        (NA_STRING == STRING_ELT(x, 0)))
        return 0;
    *name = CHAR(STRING_ELT(x, 0));
    return ('\0' != (*name)[0]) && (strlen(*name) < gs_MAXSTR);
}

static int st_phase_match_names(char *const *from, size_t n_from,
                                char *const *to, size_t n_to,
                                size_t **match_out)
{
    st_phase_r_name_index *index = NULL;
    size_t *match = NULL;
    unsigned char *seen = NULL;
    int status = -1;

    if ((NULL == from) || (NULL == to) || (NULL == match_out) ||
        (0u == n_from) || (n_from != n_to))
        return -1;
    *match_out = NULL;

    index = malloc(n_to * sizeof(*index));
    match = malloc(n_from * sizeof(*match));
    seen = calloc(n_to, sizeof(*seen));
    if ((NULL == index) || (NULL == match) || (NULL == seen))
        goto cleanup;

    for (size_t i = 0u; i < n_to; i++)
    {
        if (NULL == to[i])
            goto cleanup;
        index[i].name = to[i];
        index[i].index = i;
    }
    qsort(index, n_to, sizeof(*index), st_phase_r_name_index_cmp);
    for (size_t i = 1u; i < n_to; i++)
        if (0 == strcmp(index[i - 1u].name, index[i].name))
            goto cleanup;

    for (size_t i = 0u; i < n_from; i++)
    {
        size_t lo = 0u, hi = n_to;
        const char *want;
        if (NULL == from[i])
            goto cleanup;
        want = from[i];
        while (lo < hi)
        {
            size_t mid = lo + (hi - lo) / 2u;
            int cmp = strcmp(index[mid].name, want);
            if (cmp < 0)
                lo = mid + 1u;
            else
                hi = mid;
        }
        if ((lo >= n_to) || (0 != strcmp(index[lo].name, want)))
            goto cleanup;
        match[i] = index[lo].index;
        if (seen[match[i]])
            goto cleanup;
        seen[match[i]] = 1u;
    }

    *match_out = match;
    match = NULL;
    status = 0;

cleanup:
    free(index);
    free(match);
    free(seen);
    return status;
}

static int st_switch_marker_alleles(const gs_varset_type *gv, size_t marker,
                                    GS_ALLELE_TYPE allele[2], size_t *n_allele)
{
    size_t n = 0u;

    if ((NULL == gv) || (NULL == allele) || (NULL == n_allele) ||
        (marker >= (size_t)gv->NoMar) || !st_allele_present(gv))
        return -1;

    for (size_t i = 0u; i < (size_t)gv->NoInd; i++)
    {
        size_t idx = marker * (size_t)gv->NoInd + i;
        GS_ALLELE_TYPE x[2] = {st_allele_get(gv, idx, 0),
                               st_allele_get(gv, idx, 1)};
        for (size_t h = 0u; h < 2u; h++)
        {
            size_t k;
            if (GS_ALLELE_MISSING == x[h])
                continue;
            if (0u == x[h])
                return -1;
            for (k = 0u; k < n; k++)
                if (allele[k] == x[h])
                    break;
            if (k < n)
                continue;
            if (n >= 2u)
                return 1;
            allele[n++] = x[h];
        }
    }

    if ((2u == n) && (allele[1] < allele[0]))
    {
        GS_ALLELE_TYPE z = allele[0];
        allele[0] = allele[1];
        allele[1] = z;
    }
    *n_allele = n;
    return 0;
}

static int st_phase_allele_rank(GS_ALLELE_TYPE x,
                                const GS_ALLELE_TYPE allele[2], size_t n)
{
    if (GS_ALLELE_MISSING == x)
        return -1;
    for (size_t k = 0u; k < n; k++)
        if (x == allele[k])
            return (int)k;
    return -2;
}

SEXP st_destroy_phase(SEXP data_set)
{
    const char *name = NULL;
    gs_varset_type *gv0 = NULL;
    size_t ncell = 0u, nhet = 0u, nswap = 0u;
    SEXP ans, names;

    if (!st_phase_step1_dataset_name(data_set, &name))
        error("data.set must name an existing SelectionTools data set");
    gv0 = st_phase_step1_find_dataset(name);
    if (NULL == gv0)
        error("SelectionTools data set not found: %s", name);
    if (!gs_marker_genotypes_valid(gv0) ||
        (0 != gs_size_mul((size_t)gv0->NoMar, (size_t)gv0->NoInd, &ncell)))
        error("Invalid marker genotypes in data.set '%s'", name);

    GetRNGstate();
    for (size_t idx = 0u; idx < ncell; idx++)
        if ((GS_ALLELE_MISSING != st_allele_get(gv0, idx, 0)) &&
            (st_allele_get(gv0, idx, 0) != st_allele_get(gv0, idx, 1)))
        {
            nhet++;
            if (unif_rand() < 0.5)
            {
                GS_ALLELE_TYPE a = st_allele_get(gv0, idx, 0),
                               b = st_allele_get(gv0, idx, 1);
                if (st_allele_write_known_pair(gv0, idx, b, a))
                    error("Invalid phase-destruction allele dictionary");
                nswap++;
            }
        }
    PutRNGstate();
    gv0->phased = 0;

    PROTECT(ans = allocVector(VECSXP, 2));
    PROTECT(names = allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("heterozygotes"));
    SET_STRING_ELT(names, 1, mkChar("swapped"));
    SET_VECTOR_ELT(ans, 0, ScalarReal((double)nhet));
    SET_VECTOR_ELT(ans, 1, ScalarReal((double)nswap));
    setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(2);
    return ans;
}

SEXP st_switch_error(SEXP true_data_set, SEXP data_set)
{
    const char *true_name = NULL, *name = NULL;
    gs_varset_type *truth = NULL, *test = NULL;
    size_t *ind_match = NULL, *marker_match = NULL;
    gs_UNSIGNED *truth_storage = NULL, *test_storage = NULL;
    GS_ALLELE_TYPE(*truth_code)[2] = NULL, (*test_code)[2] = NULL;
    unsigned char *truth_n = NULL, *test_n = NULL;
    size_t nhet = 0u, nscored = 0u, nmissing = 0u, nmismatch = 0u;
    size_t nlink = 0u, nerror = 0u, M;
    SEXP ans, names;

    if (!st_phase_step1_dataset_name(true_data_set, &true_name))
        error("true.data.set must name an existing SelectionTools data set");
    if (!st_phase_step1_dataset_name(data_set, &name))
        error("data.set must name an existing SelectionTools data set");
    truth = st_phase_step1_find_dataset(true_name);
    test = st_phase_step1_find_dataset(name);
    if (NULL == truth)
        error("SelectionTools data set not found: %s", true_name);
    if (NULL == test)
        error("SelectionTools data set not found: %s", name);
    if (!gs_marker_genotypes_valid(truth) || !gs_marker_genotypes_valid(test))
        error("Invalid marker genotypes in switch-error data set");
    if ((truth->NoInd != test->NoInd) || (truth->NoMar != test->NoMar))
        error("Switch-error data sets must contain the same numbers of "
              "individuals and markers");

    truth_storage = gs_map_storage_index(truth, 1, 1);
    test_storage = gs_map_storage_index(test, 1, 1);
    if ((NULL == truth_storage) || (NULL == test_storage))
    {
        free(truth_storage);
        free(test_storage);
        error("Switch-error data sets require complete linkage maps with "
              "unique marker names");
    }
    free(truth_storage);
    truth_storage = NULL;
    free(test_storage);
    test_storage = NULL;

    char **truth_marker_names =
        (char **)st_work_alloc((size_t)truth->NoMar, sizeof(char *));
    char **test_marker_names =
        (char **)st_work_alloc((size_t)test->NoMar, sizeof(char *));
    for (gs_UNSIGNED m = 0; m < truth->NoMar; m++)
        truth_marker_names[m] = truth->MarName[m];
    for (gs_UNSIGNED m = 0; m < test->NoMar; m++)
        test_marker_names[m] = test->MarName[m];
    if ((0 != st_phase_match_names(truth->IndName, (size_t)truth->NoInd,
                                   test->IndName, (size_t)test->NoInd,
                                   &ind_match)) ||
        (0 != st_phase_match_names(truth_marker_names, (size_t)truth->NoMar,
                                   test_marker_names, (size_t)test->NoMar,
                                   &marker_match)))
    {
        free(ind_match);
        free(marker_match);
        error("Switch-error data sets must contain identical unique "
              "individual and marker names");
    }

    M = (size_t)truth->NoMar;
    truth_code = malloc(M * sizeof(*truth_code));
    test_code = malloc(M * sizeof(*test_code));
    truth_n = malloc(M * sizeof(*truth_n));
    test_n = malloc(M * sizeof(*test_n));
    if ((NULL == truth_code) || (NULL == test_code) || (NULL == truth_n) ||
        (NULL == test_n))
    {
        free(ind_match);
        free(marker_match);
        free(truth_code);
        free(test_code);
        free(truth_n);
        free(test_n);
        error("Memory allocation error in st.switch.error");
    }

    for (size_t m = 0u; m < M; m++)
    {
        size_t nt = 0u, nx = 0u, mt = marker_match[m];
        int rt = st_switch_marker_alleles(truth, m, truth_code[m], &nt);
        int rx = st_switch_marker_alleles(test, mt, test_code[m], &nx);
        if ((1 == rt) || (1 == rx))
        {
            free(ind_match);
            free(marker_match);
            free(truth_code);
            free(test_code);
            free(truth_n);
            free(test_n);
            error("st.switch.error currently requires biallelic markers");
        }
        if ((0 != rt) || (0 != rx))
        {
            free(ind_match);
            free(marker_match);
            free(truth_code);
            free(test_code);
            free(truth_n);
            free(test_n);
            error("Invalid allele coding in switch-error data set");
        }
        truth_n[m] = (unsigned char)nt;
        test_n[m] = (unsigned char)nx;
    }

    for (size_t i = 0u; i < (size_t)truth->NoInd; i++)
    {
        int previous_chrom = INT_MIN;
        int have_previous = 0, previous_orientation = 0;
        size_t j = ind_match[i];

        for (size_t l = 0u; l < (size_t)truth->NLoc; l++)
        {
            size_t m = (size_t)truth->mdl[l];
            size_t mt = marker_match[m];
            int chrom = truth->map[l]->chrom;
            size_t idx = m * (size_t)truth->NoInd + i;
            size_t idt = mt * (size_t)test->NoInd + j;
            GS_ALLELE_TYPE ta1 = st_allele_get(truth, idx, 0),
                           ta2 = st_allele_get(truth, idx, 1);
            GS_ALLELE_TYPE xa1 = st_allele_get(test, idt, 0),
                           xa2 = st_allele_get(test, idt, 1);
            int rt1, rt2, rx1, rx2, orientation;

            if (chrom != previous_chrom)
            {
                previous_chrom = chrom;
                have_previous = 0;
            }
            if ((GS_ALLELE_MISSING == ta1) || (GS_ALLELE_MISSING == ta2) ||
                (ta1 == ta2))
                continue;
            nhet++;

            if ((GS_ALLELE_MISSING == xa1) || (GS_ALLELE_MISSING == xa2))
            {
                nmissing++;
                continue;
            }

            rt1 = st_phase_allele_rank(ta1, truth_code[m], truth_n[m]);
            rt2 = st_phase_allele_rank(ta2, truth_code[m], truth_n[m]);
            rx1 = st_phase_allele_rank(xa1, test_code[m], test_n[m]);
            rx2 = st_phase_allele_rank(xa2, test_code[m], test_n[m]);
            if ((2u != truth_n[m]) || (2u != test_n[m]) || (rt1 < 0) ||
                (rt2 < 0) || (rx1 < 0) || (rx2 < 0) || (rt1 == rt2) ||
                (rx1 == rx2))
            {
                nmismatch++;
                continue;
            }

            orientation = (rt1 == rx1) ? 0 : 1;
            nscored++;
            if (have_previous)
            {
                nlink++;
                if (orientation != previous_orientation)
                    nerror++;
            }
            previous_orientation = orientation;
            have_previous = 1;
        }
    }

    free(ind_match);
    ind_match = NULL;
    free(marker_match);
    marker_match = NULL;
    free(truth_code);
    truth_code = NULL;
    free(test_code);
    test_code = NULL;
    free(truth_n);
    truth_n = NULL;
    free(test_n);
    test_n = NULL;

    double value[] = {
        (double)nerror, (double)nlink,
        nlink ? (double)nerror / (double)nlink : NA_REAL,
        nlink ? 100.0 * (double)nerror / (double)nlink : NA_REAL,
        (double)nhet, (double)nscored,
        nhet ? (double)nscored / (double)nhet : NA_REAL,
        (double)nmissing, (double)nmismatch};
    const char *labels[] = {"errors", "links", "switch.error", "switch.error.percent",
                            "heterozygotes", "scored.heterozygotes", "coverage",
                            "missing.heterozygotes", "genotype.mismatches"};
    PROTECT(ans = allocVector(VECSXP,1));
    PROTECT(names = allocVector(STRSXP,9));
    SEXP result = PROTECT(allocVector(REALSXP,9)),
         classes = PROTECT(allocVector(STRSXP,2));
    SET_VECTOR_ELT(ans,0,result);
    memcpy(REAL(result),value,sizeof(value));
    for (int i = 0; i < 9; i++) SET_STRING_ELT(names,i,mkChar(labels[i]));
    SET_STRING_ELT(classes,0,mkChar("st.switch.error"));
    SET_STRING_ELT(classes,1,mkChar("data.frame"));
    setAttrib(ans,R_NamesSymbol,mkString("result"));
    setAttrib(ans,R_RowNamesSymbol,names);
    setAttrib(ans,R_ClassSymbol,classes);
    UNPROTECT(4);
    return ans;
}

#include <R_ext/Rdynload.h>
#include <setjmp.h>
#include <stdatomic.h>
#include <inttypes.h>
#include <stdarg.h>

#ifndef MT
static int omp_get_thread_num(void)
{
    return 0;
}
static int omp_get_num_threads(void)
{
    return 1;
}
#endif

typedef struct
{
    jmp_buf jump;
    int stage;
} st_guard;
static _Thread_local st_guard *st_active_guard;
static _Atomic int st_worker_error, st_max_team;
static const char *st_phase_dataset;

enum
{
    ST_PHASE_MEMORY = 1,
    ST_PHASE_SIZE,
    ST_PHASE_ALLELES,
    ST_PHASE_DONORS,
    ST_PHASE_SAMPLE_COUNT,
    ST_PHASE_REFERENCE,
    ST_PHASE_PBWT
};

static _Noreturn void st_phase_error(const char *format, ...)
{
    char message[gs_MAXSTR];
    va_list args;
    va_start(args, format);
    vsnprintf(message, sizeof(message), format, args);
    va_end(args);
    if (st_phase_dataset)
        Rf_error("Data set '%s': %s", st_phase_dataset, message);
    Rf_error("st.phase: %s", message);
}
static _Noreturn void st_report_error(int code, int stage)
{
    const char *message;
    switch (code)
    {
    case ST_PHASE_MEMORY:
        message = "Memory allocation error";
        break;
    case ST_PHASE_SIZE:
        message = "Memory allocation size overflow";
        break;
    case ST_PHASE_ALLELES:
        message = "Selected markers must be monomorphic or biallelic";
        break;
    case ST_PHASE_DONORS:
        message = "No donor haplotypes available";
        break;
    case ST_PHASE_SAMPLE_COUNT:
        message = "Invalid sample count for parameter estimation";
        break;
    case ST_PHASE_REFERENCE:
        message = "Phasing attempted to change a reference haplotype";
        break;
    case ST_PHASE_PBWT:
        message = "Invalid donor-search workspace dimensions";
        break;
    default:
        message = "Internal phasing error";
    }
    const char *stages[] = {"phasing",
                            "parameter estimation",
                            "common-marker phasing",
                            "haplotype coding",
                            "common-marker donor search",
                            "rare-marker donor search",
                            "rare-marker phasing and imputation"};
    if (stage < 0 || stage > 6)
        stage = 0;
    st_phase_error("%s during %s", message, stages[stage]);
}
static _Noreturn void st_failure(int code)
{
    if (st_active_guard)
    {
        int expected = 0;
        atomic_compare_exchange_strong(&st_worker_error, &expected,
                                       code + 256 * st_active_guard->stage);
        longjmp(st_active_guard->jump, 1);
    }

    st_report_error(code, 0);
}
static void st_worker_run(void (*work)(void *), void *data, int stage)
{
    if (atomic_load(&st_worker_error))
        return;
    int team = omp_get_num_threads(), previous = atomic_load(&st_max_team);
    while (previous < team &&
           !atomic_compare_exchange_weak(&st_max_team, &previous, team))
    {
    }
    st_guard guard;
    guard.stage = stage;
    st_active_guard = &guard;
    if (setjmp(guard.jump) == 0)
    {
        work(data);
    }
    st_active_guard = NULL;
}
static void st_workers_check(void)
{
    int code = atomic_load(&st_worker_error);
    if (code)
        st_report_error(code % 256, code / 256);
    R_CheckUserInterrupt();
}

typedef union st_block st_block;
union st_block
{
    struct
    {
        st_block *previous, *next;
        size_t size;
    } info;
    max_align_t alignment;
};
static st_block *st_blocks;
static int st_busy;
static size_t st_block_count;
static void *st_memory(size_t size)
{
    if (size > SIZE_MAX - sizeof(st_block))
        st_failure(ST_PHASE_SIZE);
    st_block *b = malloc(sizeof(*b) + size);
    if (!b)
        st_failure(ST_PHASE_MEMORY);
#ifdef MT
#pragma omp critical(st_allocation_list)
#endif
    {
        b->info.previous = NULL;
        b->info.next = st_blocks;
        b->info.size = size;
        if (st_blocks)
            st_blocks->info.previous = b;
        st_blocks = b;
        ++st_block_count;
    }
    return b + 1;
}
static void st_release(void *p)
{
    if (!p)
        return;
    st_block *b = (st_block *)p - 1;
#ifdef MT
#pragma omp critical(st_allocation_list)
#endif
    {
        if (b->info.previous)
            b->info.previous->info.next = b->info.next;
        else
            st_blocks = b->info.next;
        if (b->info.next)
            b->info.next->info.previous = b->info.previous;
        --st_block_count;
    }
    free(b);
}
static void *st_zero(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size)
        st_failure(ST_PHASE_SIZE);
    void *p = st_memory(n * size);
    memset(p, 0, n * size);
    return p;
}
static void *st_resize(void *p, size_t n)
{
    if (!p)
        return st_memory(n);
    size_t old = ((st_block *)p - 1)->info.size;
    void *out = st_memory(n);
    memcpy(out, p, old < n ? old : n);
    st_release(p);
    return out;
}
#define malloc st_memory
#define calloc st_zero
#define realloc st_resize
#define free st_release
#define ST_PHASE_CHECK_INTERRUPT() R_CheckUserInterrupt()

static int er_thread_limit = 1;
static int er_set_threads(int n)
{
    if (n < 1)
        return 0;
#ifndef MT
    if (n != 1)
        return 0;
#endif
    er_thread_limit = n;
    return 1;
}

typedef struct
{
    uint64_t state;
} j_random;
static void jr_seed(j_random *r, uint64_t seed)
{
    r->state = (seed ^ UINT64_C(0x5deece66d)) & ((UINT64_C(1) << 48) - 1);
}
static uint32_t jr_bits(j_random *r, unsigned bits)
{
    r->state =
        (r->state * UINT64_C(0x5deece66d) + 11) & ((UINT64_C(1) << 48) - 1);
    return (uint32_t)(r->state >> (48 - bits));
}
static uint32_t jr_int(j_random *r, uint32_t bound)
{
    if ((bound & (bound - 1)) == 0)
        return (uint32_t)(((uint64_t)bound * jr_bits(r, 31)) >> 31);
    uint32_t b, v;
    do
    {
        b = jr_bits(r, 31);
        v = b % bound;
    } while ((uint64_t)b - v + (bound - 1) > INT32_MAX);
    return v;
}

typedef struct
{
    int hap, start, last, index;
} cs_segment;
typedef struct
{
    int M, H, S, K, min_steps, words;
    int *starts, *donors;
    uint64_t *haps;
} cs_input;
typedef struct
{
    int n, *last;
    cs_segment *heap;
    uint64_t *comp;
} cs_work;

static void *cs_alloc(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size)
        st_failure(ST_PHASE_SIZE);
    void *p = calloc(n ? n : 1, size);
    if (!p)
        st_failure(ST_PHASE_MEMORY);
    return p;
}

static cs_work cs_create(const cs_input *in)
{
    cs_work w = {0, cs_alloc(in->H, sizeof(int)),
                 cs_alloc(in->K, sizeof(cs_segment)),
                 cs_alloc((size_t)in->K * in->words, sizeof(uint64_t))};
    return w;
}

static void cs_destroy(cs_work *w)
{
    free(w->last);
    free(w->heap);
    free(w->comp);
}

static void cs_offer(cs_work *w, cs_segment x)
{
    int k = w->n++;
    while (k > 0)
    {
        int parent = (k - 1) >> 1;
        if (x.last >= w->heap[parent].last)
            break;
        w->heap[k] = w->heap[parent];
        k = parent;
    }
    w->heap[k] = x;
}

static cs_segment cs_poll(cs_work *w)
{
    cs_segment out = w->heap[0], x = w->heap[--w->n];
    int k = 0, half = w->n >> 1;
    while (k < half)
    {
        int child = 2 * k + 1, right = child + 1;
        if (right < w->n && w->heap[child].last > w->heap[right].last)
            child = right;
        if (x.last <= w->heap[child].last)
            break;
        w->heap[k] = w->heap[child];
        k = child;
    }
    if (w->n)
        w->heap[k] = x;
    return out;
}

static void cs_copy(const cs_input *in, cs_work *w, cs_segment x, int end,
                    FILE *trace)
{
    if (trace)
        fprintf(trace, "SEG %d %d %d %d\n", x.index, x.hap, x.start, end);
    uint64_t *dst = w->comp + (size_t)x.index * in->words;
    const uint64_t *src = in->haps + (size_t)x.hap * in->words;
    for (int start = x.start; start < end;)
    {
        int word = start >> 6, lo = start & 63;
        int stop = end < (word + 1) * 64 ? end : (word + 1) * 64;
        int hi = stop - word * 64;
        uint64_t mask = UINT64_MAX << lo;
        if (hi < 64)
            mask &= (UINT64_C(1) << hi) - 1;
        dst[word] = (dst[word] & ~mask) | (src[word] & mask);
        start = stop;
    }
}

static void cs_trace_queue(FILE *trace, const cs_work *w, int step, int hap)
{
    if (!trace)
        return;
    fprintf(trace, "Q %d %d %d", step, hap, w->n);
    for (int i = 0; i < w->n; ++i)
    {
        cs_segment x = w->heap[i];
        fprintf(trace, " %d %d %d %d %d", x.hap, x.start, x.last, x.index,
                w->last[x.hap]);
    }
    fputc('\n', trace);
}

static void cs_add(const cs_input *in, cs_work *w, int hap, int step,
                   FILE *trace)
{
    if (w->last[hap] == -103)
    {
        while (w->n && w->heap[0].last != w->last[w->heap[0].hap])
        {
            cs_segment head = cs_poll(w);
            head.last = w->last[head.hap];
            cs_offer(w, head);
        }
        if (w->n == in->K || (w->n && step - w->heap[0].last >= in->min_steps))
        {
            cs_segment head = cs_poll(w);
            int next = in->starts[(head.last + step) >> 1];
            w->last[head.hap] = -103;
            cs_copy(in, w, head, next, trace);
            head.hap = hap;
            head.start = next;
            head.last = step;
            cs_offer(w, head);
        }
        else
        {
            cs_segment head = {hap, 0, step, w->n};
            cs_offer(w, head);
        }
    }
    w->last[hap] = step;
    cs_trace_queue(trace, w, step, hap);
}

static int cs_run(const cs_input *in, cs_work *w, int sample, uint64_t seed,
                  FILE *trace)
{
    w->n = 0;
    for (int h = 0; h < in->H; ++h)
        w->last[h] = -103;
    for (int step = 0; step < in->S; ++step)
    {
        for (int j = 0; j < 2; ++j)
        {
            int hap = in->donors[2 * step + j];
            if (hap >= 0)
                cs_add(in, w, hap, step, trace);
        }
    }
    if (!w->n)
    {
        int count = in->H - 2 < in->K ? in->H - 2 : in->K;
        if (count <= 0)
            st_failure(ST_PHASE_DONORS);
        j_random rng;
        jr_seed(&rng, seed + (uint64_t)sample);
        int index = 0;
        for (int j = 0; j < count; ++j)
        {
            int h;
            do
            {
                h = (int)jr_int(&rng, (uint32_t)in->H);
                if (trace)
                    fprintf(trace, "DRAW %d\n", h);
            } while ((h >> 1) == sample);
            if (w->last[h] == -103)
            {
                cs_segment head = {h, 0, 0, index++};
                cs_offer(w, head);
                w->last[h] = 0;
            }
            cs_trace_queue(trace, w, -1, h);
        }
    }
    int count = w->n;
    while (w->n)
        cs_copy(in, w, cs_poll(w), in->M, trace);
    if (trace)
    {
        fprintf(trace, "FINAL %d\n", count);
        for (int c = 0; c < count; ++c)
        {
            fprintf(trace, "ALLELES %d", c);
            for (int m = 0; m < in->M; ++m)
                fprintf(trace, " %d",
                        (int)((w->comp[(size_t)c * in->words + (m >> 6)] >>
                               (m & 63)) &
                              1));
            fputc('\n', trace);
        }
    }
    return count;
}

typedef struct
{
    int count;
    double mismatch, distance, switches;
} pe_result;

typedef struct
{
    pe_result initial;
    int markers, states;
    float mismatch;
    const float *distance, *recomb;
    const unsigned char *discord;
} pe_input;

typedef struct
{
    float *fwd, *bwd, *saved;
} pe_work;

static uint32_t pe_float_bits(float x)
{
    uint32_t b;
    memcpy(&b, &x, sizeof b);
    return b;
}
static uint64_t pe_double_bits(double x)
{
    uint64_t b;
    memcpy(&b, &x, sizeof b);
    return b;
}
static pe_work pe_create(int markers, int states)
{
    pe_work w;
    w.fwd = calloc((size_t)states, sizeof(float));
    w.bwd = calloc((size_t)states, sizeof(float));
    w.saved = calloc((size_t)markers * states, sizeof(float));
    if (!w.fwd || !w.bwd || !w.saved)
    {
        st_failure(ST_PHASE_MEMORY);
    }
    return w;
}
static void pe_destroy(pe_work *w)
{
    free(w->fwd);
    free(w->bwd);
    free(w->saved);
}
static void pe_row(FILE *trace, const char *tag, int marker, const float *a,
                   int n)
{
    if (!trace)
        return;
    fprintf(trace, "%s %d", tag, marker);
    for (int k = 0; k < n; ++k)
        fprintf(trace, " %u", pe_float_bits(a[k]));
    fputc('\n', trace);
}
static void pe_totals(FILE *trace, const char *tag, int marker, pe_result r)
{
    if (trace)
        fprintf(trace, "%s %d %d %llu %llu %llu\n", tag, marker, r.count,
                (unsigned long long)pe_double_bits(r.mismatch),
                (unsigned long long)pe_double_bits(r.distance),
                (unsigned long long)pe_double_bits(r.switches));
}

static pe_result pe_run(const pe_input *in, pe_work *w, FILE *trace)
{
    pe_result r = in->initial;
    int M = in->markers, K = in->states;
    if (K <= 1)
        return r;
    float emission[2] = {1.0f - in->mismatch, in->mismatch};
    for (int k = 0; k < K; ++k)
        w->bwd[k] = w->saved[(size_t)(M - 1) * K + k] = 1.0f;
    pe_row(trace, "B", M - 1, w->bwd, K);
    for (int m = M - 2; m >= 0; --m)
    {
        float sum = 0.0f;
        const unsigned char *d = in->discord + (size_t)(m + 1) * K;
        for (int k = 0; k < K; ++k)
        {
            w->bwd[k] *= emission[d[k]];
            sum += w->bwd[k];
        }
        float shift = in->recomb[m + 1] / K;
        float scale = (1.0f - in->recomb[m + 1]) / sum;
        for (int k = 0; k < K; ++k)
            w->bwd[k] = scale * w->bwd[k] + shift;
        memcpy(w->saved + (size_t)m * K, w->bwd, (size_t)K * sizeof(float));
        pe_row(trace, "B", m, w->bwd, K);
    }
    float hfactor = K / (K - 1.0f);
    for (int k = 0; k < K; ++k)
        w->fwd[k] = 1.0f / K;
    float last_sum = 1.0f;
    for (int m = 0; m < M; ++m)
    {
        float ps = in->recomb[m], shift = ps / K;
        float scale = (1.0f - ps) / last_sum;
        float no_switch_scale = ((1.0f - ps) + shift) / last_sum;
        float joint = 0.0f, state_sum = 0.0f, fwd_sum = 0.0f,
              mismatch_sum = 0.0f;
        const float *b = w->saved + (size_t)m * K;
        const unsigned char *d = in->discord + (size_t)m * K;
        for (int k = 0; k < K; ++k)
        {
            float em = emission[d[k]];
            joint += b[k] * em * no_switch_scale * w->fwd[k];
            w->fwd[k] = em * (scale * w->fwd[k] + shift);
            fwd_sum += w->fwd[k];
            float prob = w->fwd[k] * b[k];
            state_sum += prob;
            if (d[k] > 0)
                mismatch_sum += prob;
        }
        ++r.count;
        r.mismatch += mismatch_sum / state_sum;
        double switch_prob = hfactor * (1.0f - joint / state_sum);
        if (switch_prob > 0.0)
        {
            r.distance += in->distance[m];
            r.switches += switch_prob;
        }
        last_sum = fwd_sum;
        pe_row(trace, "F", m, w->fwd, K);
        pe_totals(trace, "S", m, r);
    }
    return r;
}

static int pe_compare_mismatch(const void *aa, const void *bb)
{
    const pe_result *a = aa, *b = bb;
    if (a->mismatch != b->mismatch)
        return a->mismatch < b->mismatch ? -1 : 1;
    return (a->count > b->count) - (a->count < b->count);
}
static int pe_compare_switch(const void *aa, const void *bb)
{
    const pe_result *a = aa, *b = bb;
    if (a->distance != b->distance)
        return a->distance < b->distance ? -1 : 1;
    return (a->switches > b->switches) - (a->switches < b->switches);
}
static void pe_estimates(const pe_result *rows, int n, float *mismatch,
                         float *intensity)
{
    pe_result *tmp = calloc((size_t)(n ? n : 1), sizeof(*tmp));
    if (!tmp)
        st_failure(ST_PHASE_MEMORY);
    int count = 0;
    for (int i = 0; i < n; ++i)
    {
        if (rows[i].count > 0 && rows[i].mismatch > 0 &&
            isfinite(rows[i].mismatch))
            tmp[count++] = rows[i];
    }
    qsort(tmp, (size_t)count, sizeof(*tmp), pe_compare_mismatch);
    int64_t markers = 0;
    double sum = 0.0;
    for (int i = 0; i < count; ++i)
    {
        markers += tmp[i].count;
        sum += tmp[i].mismatch;
    }
    *mismatch = markers == 0 ? NAN : (float)(sum / markers);
    count = 0;
    for (int i = 0; i < n; ++i)
    {
        if (rows[i].distance > 0 && rows[i].switches > 0 &&
            isfinite(rows[i].distance) && isfinite(rows[i].switches))
            tmp[count++] = rows[i];
    }
    qsort(tmp, (size_t)count, sizeof(*tmp), pe_compare_switch);
    double distance = 0.0, switches = 0.0;
    for (int i = 0; i < count; ++i)
    {
        switches += tmp[i].switches;
        distance += tmp[i].distance;
    }
    *intensity = distance == 0.0 ? NAN : (float)(switches / distance);
    free(tmp);
}

static int pe_samples(int *indices, int samples, j_random *rng)
{
    for (int i = 0; i < samples; ++i)
        indices[i] = i;
    if (samples <= 500)
        return samples;
    for (int j = 0; j < 500; ++j)
    {
        int k = j + (int)jr_int(rng, (uint32_t)(samples - j));
        int tmp = indices[j];
        indices[j] = indices[k];
        indices[k] = tmp;
    }
    return 500;
}

typedef struct
{
    cs_input states;
    int samples, initialize;
    uint64_t seed;
    int *donors;
    float *distance;
    float mismatch, intensity;
} er_input;
typedef struct
{
    float mismatch, intensity;
    int rounds;
} er_result;
static uint64_t er_hash(const unsigned char *x, size_t n)
{
    uint64_t h = UINT64_C(1469598103934665603);
    for (size_t j = 0; j < n; ++j)
        h = (h ^ x[j]) * UINT64_C(1099511628211);
    return h;
}
typedef struct
{
    cs_input ci;
    cs_work cw;
    pe_work work;
    unsigned char *discord;
} er_worker;
typedef struct
{
    const er_input *in;
    const int *selected;
    float mismatch;
    const float *recomb;
    pe_result *rows;
    FILE *trace;
} er_context;
static void er_process_sample(const er_context *ctx, er_worker *ws, int j)
{
    cs_input ci = ws->ci;
    cs_work *cw = &ws->cw;
    pe_work *work = &ws->work;
    unsigned char *discord = ws->discord;
    int sample = ctx->selected[j];
    for (int s = 0; s < ci.S; ++s)
        for (int h = 0; h < 2; ++h)
            ci.donors[2 * s + h] =
                ctx->in->donors[(size_t)s * ci.H + 2 * sample + h];
    int K = cs_run(&ci, cw, sample, ctx->in->seed, NULL);
    for (int h = 0; h < 2; ++h)
        for (int m = 0; m < ci.M; ++m)
        {
            int allele =
                (ci.haps[(size_t)(2 * sample + h) * ci.words + (m >> 6)] >>
                 (m & 63)) &
                1;
            for (int k = 0; k < K; ++k)
            {
                int ref =
                    (cw->comp[(size_t)k * ci.words + (m >> 6)] >> (m & 63)) & 1;
                discord[((size_t)h * ci.M + m) * K + k] =
                    (unsigned char)(allele != ref);
            }
        }
    if (ctx->trace)
        fprintf(ctx->trace, "STATES %d %d %llu\n", sample, K,
                (unsigned long long)er_hash(discord, (size_t)2 * ci.M * K));
    pe_input pi = {{0, 0, 0, 0},      ci.M,        K,      ctx->mismatch,
                   ctx->in->distance, ctx->recomb, discord};
    pe_result r = pe_run(&pi, work, NULL);
    pi.initial = r;
    pi.discord = discord + (size_t)ci.M * K;
    ctx->rows[j] = pe_run(&pi, work, NULL);
    pe_totals(ctx->trace, "SAMPLE", sample, ctx->rows[j]);
}

typedef struct
{
    const er_context *ctx;
    er_worker *ws;
    int j;
} er_task;
static void er_task_run(void *data)
{
    er_task *t = data;
    er_process_sample(t->ctx, t->ws, t->j);
}

static er_result er_run(const er_input *in, FILE *trace)
{
    cs_input ci = in->states;
    int workers = trace ? 1 : er_thread_limit;
    if (workers > in->samples)
        workers = in->samples;
    er_worker *pool = cs_alloc(workers, sizeof(er_worker));
    for (int w = 0; w < workers; ++w)
    {
        pool[w].ci = ci;
        pool[w].ci.donors = cs_alloc((size_t)2 * ci.S, sizeof(int));
        pool[w].cw = cs_create(&pool[w].ci);
        pool[w].work = pe_create(ci.M, ci.K);
        pool[w].discord = cs_alloc((size_t)2 * ci.M * ci.K, 1);
    }
    int *selected = cs_alloc(in->samples, sizeof(int));
    pe_result *rows = cs_alloc(in->samples, sizeof(pe_result));
    float *recomb = cs_alloc(ci.M, sizeof(float));
    j_random rng;
    jr_seed(&rng, in->seed);
    er_result out = {in->mismatch, in->intensity, 0};
    int limit = in->initialize ? 15 : 1;
    for (int round = 0; round < limit; ++round)
    {
        float old_intensity = out.intensity;
        for (int m = 0; m < ci.M; ++m)
            recomb[m] = (float)-expm1(-(double)out.intensity * in->distance[m]);
        int n = pe_samples(selected, in->samples, &rng);
        if (n < 1 || n > 500)
            st_failure(ST_PHASE_SAMPLE_COUNT);
        if (trace)
        {
            fprintf(trace, "SELECT %d", round);
            for (int j = 0; j < n; ++j)
                fprintf(trace, " %d", selected[j]);
            fputc('\n', trace);
            fprintf(trace, "RECOMB %d", round);
            for (int m = 0; m < ci.M; ++m)
                fprintf(trace, " %u", pe_float_bits(recomb[m]));
            fputc('\n', trace);
        }
        er_context context = {in, selected, out.mismatch, recomb, rows, trace};

#ifdef MT
#pragma omp parallel for num_threads(workers) schedule(static) if (workers > 1)
#endif
        for (int j = 0; j < n; ++j)
        {
            int tid = 0;
#ifdef MT
            tid = omp_get_thread_num();
#endif
            er_task task = {&context, pool + tid, j};
            st_worker_run(er_task_run, &task, 1);
        }
        st_workers_check();
        float pm, ri;
        pe_estimates(rows, n, &pm, &ri);
        if (isfinite(pm) && pm > out.mismatch)
            out.mismatch = pm;
        if (isfinite(ri) && ri > 0)
            out.intensity = ri;
        ++out.rounds;
        if (trace)
            fprintf(trace, "UPDATE %u %u %u %u\n", pe_float_bits(pm),
                    pe_float_bits(ri), pe_float_bits(out.mismatch),
                    pe_float_bits(out.intensity));
        if (in->initialize &&
            fabsf(out.intensity - old_intensity) <= 0.1 * (double)old_intensity)
            break;
    }
    if (trace)
        fprintf(trace, "FINAL %d %u %u\n", out.rounds,
                pe_float_bits(out.mismatch), pe_float_bits(out.intensity));
    for (int w = 0; w < workers; ++w)
    {
        pe_destroy(&pool[w].work);
        cs_destroy(&pool[w].cw);
        free(pool[w].ci.donors);
        free(pool[w].discord);
    }
    free(pool);
    free(selected);
    free(rows);
    free(recomb);
    return out;
}

static int st_phase_threads = 1;

#define it_find_donors st_serial_find_donors

typedef struct
{
    int H;
    int *p, *offset, *a, *d;
} pbwt_workspace;
static int pbwt_alloc(pbwt_workspace *w, int H)
{
    memset(w, 0, sizeof(*w));
    if (H < 1 || (size_t)H > SIZE_MAX / sizeof(int))
        return -1;
    w->H = H;
    w->p = calloc((size_t)H, sizeof(int));
    w->offset = calloc((size_t)H, sizeof(int));
    w->a = calloc((size_t)H, sizeof(int));
    w->d = calloc((size_t)H, sizeof(int));
    if (!w->p || !w->offset || !w->a || !w->d)
    {
        free(w->p);
        free(w->offset);
        free(w->a);
        free(w->d);
        memset(w, 0, sizeof(*w));
        return -1;
    }
    return 0;
}
static void pbwt_free(pbwt_workspace *w)
{
    free(w->p);
    free(w->offset);
    free(w->a);
    free(w->d);
    memset(w, 0, sizeof(*w));
}
static void pbwt_update(pbwt_workspace *w, const int *rec, int K, int marker,
                        int backward, int *a, int *d)
{
    int H = w->H;
    memset(w->offset, 0, (size_t)K * sizeof(int));
    for (int h = 0; h < H; ++h)
        ++w->offset[rec[h]];
    int n = 0;
    for (int k = 0; k < K; ++k)
    {
        int count = w->offset[k];
        w->offset[k] = n;
        n += count;
        w->p[k] = marker + (backward ? -1 : 1);
    }
    for (int i = 0; i < H; ++i)
    {
        int h = a[i], allele = rec[h];
        for (int k = 0; k < K; ++k)
        {
            if (backward ? d[i] < w->p[k] : d[i] > w->p[k])
                w->p[k] = d[i];
        }
        int at = w->offset[allele]++;
        w->a[at] = h;
        w->d[at] = w->p[allele];
        w->p[allele] = backward ? INT_MAX : INT_MIN;
    }
    memcpy(a, w->a, (size_t)H * sizeof(int));
    memcpy(d, w->d, (size_t)H * sizeof(int));
}

typedef struct
{
    int sample, start, end;
} ibs_segment;
typedef struct
{
    ibs_segment *v;
    size_t n, cap;
} ibs_list;
typedef struct
{
    int M, N;
    ibs_list *sample;
} ibs_result;
typedef struct
{
    int *v, n, hom;
} ibs_cluster;
typedef struct
{
    ibs_cluster *v;
    size_t n, cap;
} ibs_clusters;
static void *ibs_calloc(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size)
        st_failure(ST_PHASE_SIZE);
    void *x = calloc(n ? n : 1, size);
    if (!x)
        st_failure(ST_PHASE_MEMORY);
    return x;
}
static void *ibs_grow(void *old, size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size)
        st_failure(ST_PHASE_SIZE);
    void *x = realloc(old, n * size);
    if (!x)
        st_failure(ST_PHASE_MEMORY);
    return x;
}
static void ibs_append(ibs_list *l, ibs_segment x)
{
    if (l->n == l->cap)
    {
        l->cap = l->cap ? 2 * l->cap : 8;
        l->v = ibs_grow(l->v, l->cap, sizeof(*l->v));
    }
    l->v[l->n++] = x;
}
static void ibs_cluster_append(ibs_clusters *l, ibs_cluster x)
{
    if (l->n == l->cap)
    {
        l->cap = l->cap ? 2 * l->cap : 8;
        l->v = ibs_grow(l->v, l->cap, sizeof(*l->v));
    }
    l->v[l->n++] = x;
}
static int ibs_gt(const int *g, int N, int m, int s)
{
    size_t z = ((size_t)m * N + s) * 2;
    return g[z] < 0 || g[z + 1] < 0 ? -1 : g[z] + g[z + 1];
}
static int ibs_next(const double *pos, int M, int start, int *use)
{
    double min = pos[start] + (double)0.02f;
    int next = start + 1, n = 0;
    while (next < M && n < 50)
    {
        if (use[next])
        {
            if (pos[next] < min)
                use[next] = 0;
            else
            {
                ++n;
                min = pos[next] + (double)0.02f;
            }
        }
        ++next;
    }
    return next;
}
static void ibs_trace_list(FILE *f, const char *stage, int s, const ibs_list *l)
{
    if (!f)
        return;
    fprintf(f, "%s %d", stage, s);
    for (size_t j = 0; j < l->n; ++j)
        fprintf(f, " %d %d %d", l->v[j].sample, l->v[j].start, l->v[j].end);
    fputc('\n', f);
}
static int ibs_cmp(const void *av, const void *bv)
{
    const ibs_segment *a = av, *b = bv;
    if (a->sample != b->sample)
        return (a->sample > b->sample) - (a->sample < b->sample);
    if (a->start != b->start)
        return (a->start > b->start) - (a->start < b->start);
    return (a->end > b->end) - (a->end < b->end);
}
static void ibs_merge(ibs_list *l, const double *pos)
{
    if (l->n < 2)
        return;
    size_t out = 0;
    ibs_segment prev = l->v[0];
    for (size_t j = 1; j < l->n; ++j)
    {
        ibs_segment next = l->v[j];
        if (prev.sample == next.sample &&
            pos[next.start] - pos[prev.end] <= 4.0)
        {
            prev.end = next.end;
        }
        else
        {
            l->v[out++] = prev;
            prev = next;
        }
    }
    l->v[out++] = prev;
    l->n = out;
}
static int ibs_compatible(const int *g, int N, int m, int s, int t)
{
    size_t a = ((size_t)m * N + s) * 2, b = ((size_t)m * N + t) * 2;
    int a1 = g[a], a2 = g[a + 1], b1 = g[b], b2 = g[b + 1];
    return ((a1 < 0 || b1 < 0 || a1 == b1) && (a2 < 0 || b2 < 0 || a2 == b2)) ||
           ((a1 < 0 || b2 < 0 || a1 == b2) && (a2 < 0 || b1 < 0 || a2 == b1));
}
static ibs_result ibs_build(const int *g, const double *pos, const float *maf,
                            int M, int N, FILE *trace)
{
    ibs_result out = {M, N, ibs_calloc(N, sizeof(ibs_list))};
    int *use = ibs_calloc(M, sizeof(int)), *starts = ibs_calloc(M, sizeof(int)),
        *markers = ibs_calloc(M, sizeof(int));
    float miss_product = 0.1f * (float)(2 * N);
    int maxmiss = (int)ceil((double)miss_product);
    for (int m = 0; m < M; ++m)
    {
        int missing = 0;
        for (int h = 0; h < 2 * N; ++h)
            missing += g[(size_t)m * 2 * N + h] < 0;
        use[m] = maf[m] >= 0.1f && missing <= maxmiss;
    }
    int W = 0, last = 0, next = ibs_next(pos, M, 0, use);
    while (next < M)
    {
        starts[W++] = last;
        last = next;
        next = ibs_next(pos, M, next, use);
    }
    if (trace)
    {
        fprintf(trace, "USE");
        for (int m = 0; m < M; ++m)
            fprintf(trace, " %d", use[m]);
        fprintf(trace, "\nSTARTS");
        for (int w = 0; w < W; ++w)
            fprintf(trace, " %d", starts[w]);
        fputc('\n', trace);
    }
    unsigned char *pairs = ibs_calloc((size_t)N * N, 1);
    for (int w = 0; w < W; ++w)
    {
        int start = starts[w], end = w + 1 < W ? starts[w + 1] : M, L = 0;
        for (int m = start; m < end; ++m)
            if (use[m])
                markers[L++] = m;
        float p = 0.1f * (float)L;
        int allowed = (int)floor((double)p);
        ibs_cluster root = {ibs_calloc(N, sizeof(int)), 0, 1};
        for (int s = 0; s < N; ++s)
        {
            int miss = 0;
            for (int j = 0; j < L; ++j)
                miss += ibs_gt(g, N, markers[j], s) < 0;
            if (miss <= allowed)
                root.v[root.n++] = s;
        }
        ibs_clusters current = {0};
        ibs_cluster_append(&current, root);
        for (int j = 0; j < L; ++j)
        {
            int m = markers[j];
            ibs_clusters children = {0};
            for (size_t c = 0; c < current.n; ++c)
            {
                ibs_cluster parent = current.v[c];
                int present[3] = {0, 0, 0};
                for (int i = 0; i < parent.n; ++i)
                {
                    int code = ibs_gt(g, N, m, parent.v[i]);
                    if (code >= 0)
                        present[code] = 1;
                }
                for (int code = 0; code < 3; ++code)
                    if (present[code])
                    {
                        ibs_cluster child = {ibs_calloc(parent.n, sizeof(int)),
                                             0, parent.hom && code != 1};
                        for (int i = 0; i < parent.n; ++i)
                        {
                            int s = parent.v[i], v = ibs_gt(g, N, m, s);
                            if (v < 0 || v == code)
                                child.v[child.n++] = s;
                        }
                        if (child.n > 1)
                            ibs_cluster_append(&children, child);
                        else
                            free(child.v);
                    }
                free(parent.v);
            }
            free(current.v);
            current = children;
        }
        memset(pairs, 0, (size_t)N * N);
        for (size_t c = 0; c < current.n; ++c)
        {
            ibs_cluster x = current.v[c];
            if (!x.hom)
                for (int i = 0; i < x.n; ++i)
                    for (int j = 0; j < x.n; ++j)
                        pairs[(size_t)x.v[i] * N + x.v[j]] = 1;
            free(x.v);
        }
        free(current.v);
        for (int s = 0; s < N; ++s)
        {
            if (trace)
                fprintf(trace, "SET %d %d", w, s);
            for (int t = 0; t < N; ++t)
                if (pairs[(size_t)s * N + t])
                {
                    if (trace)
                        fprintf(trace, " %d", t);
                    if (s != t)
                        ibs_append(out.sample + s,
                                   (ibs_segment){t, start, end - 1});
                }
            if (trace)
                fputc('\n', trace);
        }
    }
    for (int s = 0; s < N; ++s)
    {
        ibs_list *l = out.sample + s;
        ibs_trace_list(trace, "RAW", s, l);
        if (l->n)
            qsort(l->v, l->n, sizeof(ibs_segment), ibs_cmp);
        ibs_merge(l, pos);
        ibs_trace_list(trace, "MERGE1", s, l);
        for (size_t j = 0; j < l->n; ++j)
        {
            ibs_segment *v = l->v + j;
            while (v->start > 0 &&
                   ibs_compatible(g, N, v->start - 1, s, v->sample))
                --v->start;
            while (v->end + 1 < M &&
                   ibs_compatible(g, N, v->end + 1, s, v->sample))
                ++v->end;
        }
        ibs_trace_list(trace, "EXTEND", s, l);
        ibs_merge(l, pos);
        ibs_trace_list(trace, "MERGE2", s, l);
        size_t n = 0;
        for (size_t j = 0; j < l->n; ++j)
            if (pos[l->v[j].end] - pos[l->v[j].start] >= 2.0)
                l->v[n++] = l->v[j];
        l->n = n;
        ibs_trace_list(trace, "FINAL", s, l);
    }
    free(use);
    free(starts);
    free(markers);
    free(pairs);
    return out;
}
static int ibs_query(const ibs_result *x, int s, int t, int start, int end)
{
    if (s == t)
        return 1;
    const ibs_list *l = x->sample + s;
    for (size_t j = 0; j < l->n; ++j)
    {
        ibs_segment v = l->v[j];
        if (v.sample == t && start <= v.end && v.start <= end)
            return 1;
    }
    return 0;
}
static void ibs_free(ibs_result *x)
{
    for (int s = 0; s < x->N; ++s)
        free(x->sample[s].v);
    free(x->sample);
    x->sample = NULL;
}

static void choose_donors(int H, int T, int candidates, int step, int backward,
                          uint64_t seed, const int *a, const int *d,
                          const ibs_result *ibs, int start, int end,
                          int *selected)
{
    j_random r;
    jr_seed(&r, seed);
    for (int i = 0; i < H; ++i)
        if (a[i] < T)
        {
            int hap = a[i], u = i, v = i + 1, left = d[u], right = d[v];
            while (v - u < candidates &&
                   (backward ? (step <= left || step <= right)
                             : (left <= step || right <= step)))
            {
                if (backward)
                {
                    if (left <= right)
                    {
                        ++v;
                        if (d[v] < right)
                            right = d[v];
                    }
                    else
                    {
                        --u;
                        if (d[u] < left)
                            left = d[u];
                    }
                }
                else
                {
                    if (right <= left)
                    {
                        ++v;
                        if (d[v] > right)
                            right = d[v];
                    }
                    else
                    {
                        --u;
                        if (d[u] > left)
                            left = d[u];
                    }
                }
            }
            int n = v - u;
            selected[hap] = -1;
            if (n > 1)
            {
                int index = u + (int)jr_int(&r, (uint32_t)n);
                for (int j = 0; j < n; ++j, ++index)
                {
                    if (index == v)
                        index = u;
                    if (index != i &&
                        !ibs_query(ibs, hap / 2, a[index] / 2, start, end))
                    {
                        selected[hap] = a[index];
                        break;
                    }
                }
            }
        }
}

typedef struct
{
    int n;
    int *end, *type;
} hi_clusters;
static void *hi_alloc(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size)
        st_failure(ST_PHASE_SIZE);
    void *x = calloc(n ? n : 1, size);
    if (!x)
        st_failure(ST_PHASE_MEMORY);
    return x;
}
static hi_clusters hi_cluster(int M, const int *a, const int *b,
                              const double *pos, int U, const int *unph, int D,
                              const int *missing)
{
    hi_clusters out = {0, hi_alloc(M, sizeof(int)), hi_alloc(M, sizeof(int))};
    double max_end = pos[0] + (double)0.005f;
    int prev_special = 0, last_end = 0, mi = 0, ui = 0, prev_type = 2;
    int next_m = mi < D ? missing[mi++] : -1;
    int next_u = ui < U ? unph[ui++] : -1;
    for (int m = 0; m < M; ++m)
    {
        int type = m == next_m ? 0 : a[m] == b[m] ? 2 : m == next_u ? 4 : 3;
        if (type == 0)
            next_m = mi < D ? missing[mi++] : -1;
        else if (type == 4)
            next_u = ui < U ? unph[ui++] : -1;
        int special = type != 2;
        if (special || prev_special || pos[m] > max_end || m - last_end == 255)
        {
            if (m > 0)
            {
                out.end[out.n] = m;
                out.type[out.n++] = prev_type;
                max_end = pos[m] + (double)0.005f;
                last_end = m;
            }
            prev_type = type;
        }
        prev_special = special;
    }
    out.end[out.n] = M;
    out.type[out.n++] = prev_type;
    return out;
}
static void hi_mask_group(hi_clusters *cl, const int *clusters,
                          const int *markers, int n, const int *bp)
{
    if (n == 2)
        cl->type[clusters[0]] = 1;
    else if (n == 3 && bp[markers[1]] - bp[markers[0]] <= 3000)
    {
        cl->type[clusters[0]] = 1;
        cl->type[clusters[1]] = 1;
    }
}
static void hi_mask(hi_clusters *cl, const int *bp)
{
    int *clusters = hi_alloc(cl->n, sizeof(int)),
        *markers = hi_alloc(cl->n, sizeof(int));
    int n = 0, start = 0;
    for (int c = 0; c < cl->n; ++c)
    {
        if (cl->type[c] == 3)
        {
            hi_mask_group(cl, clusters, markers, n, bp);
            n = 0;
        }
        else if (cl->type[c] == 4)
        {
            clusters[n] = c;
            markers[n++] = start;
        }
        start = cl->end[c];
    }
    hi_mask_group(cl, clusters, markers, n, bp);
    free(clusters);
    free(markers);
}
static void hi_recomb(const hi_clusters *cl, const float *p, float *out)
{
    out[0] = 0.0f;
    for (int c = 1; c < cl->n; ++c)
    {
        float no = 1.0f;
        for (int m = cl->end[c - 1]; m < cl->end[c]; ++m)
            no *= 1.0f - p[m];
        out[c] = 1.0f - no;
    }
}

static int hi_equal(const uint64_t *a, const uint64_t *b, int start, int end)
{
    while (start < end)
    {
        int w = start >> 6, stop = end < (w + 1) * 64 ? end : (w + 1) * 64;
        uint64_t mask = UINT64_MAX << (start & 63);
        int hi = stop - w * 64;
        if (hi < 64)
            mask &= (UINT64_C(1) << hi) - 1;
        if ((a[w] ^ b[w]) & mask)
            return 0;
        start = stop;
    }
    return 1;
}
static void hi_mismatches(const hi_clusters *cl, int K, int words,
                          const uint64_t *a, const uint64_t *b,
                          const uint64_t *refs, unsigned char *mismatch,
                          int *missing)
{
    int start = 0, miss = 0;
    size_t plane = (size_t)cl->n * K;
    memset(mismatch, 0, plane * 3);
    for (int c = 0; c < cl->n; ++c)
    {
        int end = cl->end[c];
        if (cl->type[c] == 0 || cl->type[c] == 1)
        {
            for (int k = 0; k < K; ++k)
                missing[(size_t)miss * K + k] =
                    (int)((refs[(size_t)k * words + (start >> 6)] >>
                           (start & 63)) &
                          1);
            ++miss;
        }
        else if ((start >> 6) == ((end - 1) >> 6))
        {

            int word = start >> 6, hi = end - word * 64;
            uint64_t mask = UINT64_MAX << (start & 63);
            if (hi < 64)
                mask &= (UINT64_C(1) << hi) - 1;
            uint64_t av = a[word] & mask, bv = b[word] & mask;
            int same = av == bv;
            for (int k = 0; k < K; ++k)
            {
                uint64_t rv = refs[(size_t)k * words + word] & mask;
                size_t z = (size_t)c * K + k;
                mismatch[plane + z] = av != rv;
                mismatch[2 * plane + z] = bv != rv;
                if (same)
                    mismatch[z] = av != rv;
            }
        }
        else
        {
            int same = hi_equal(a, b, start, end);
            for (int k = 0; k < K; ++k)
            {
                const uint64_t *ref = refs + (size_t)k * words;
                size_t z = (size_t)c * K + k;
                int x = !hi_equal(a, ref, start, end);
                mismatch[plane + z] = x;
                mismatch[2 * plane + z] =
                    same ? x : !hi_equal(b, ref, start, end);
                if (same)
                    mismatch[z] = x;
            }
        }
        start = end;
    }
}
static void hi_free(hi_clusters *cl)
{
    free(cl->end);
    free(cl->type);
}

typedef struct
{
    int M, C, K, burnin;
    float mismatch, threshold;
    int *ends, *types, *a, *b, *refs;
    float *recomb;
    unsigned char *mm;
} ba_input;
typedef struct
{
    int *types, *a, *b, *unph;
    unsigned char *mm, **rows;
    float *fwd, *bwd, sums[3], *bmiss1, *bmiss2, *bhet1, *bhet2;
    int nmiss, nunph, swap, swaps;
} ba_work;
static void *ba_alloc(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size)
        st_failure(ST_PHASE_SIZE);
    void *x = calloc(n ? n : 1, size);
    if (!x)
        st_failure(ST_PHASE_MEMORY);
    return x;
}
static uint32_t ba_bits(float x)
{
    uint32_t b;
    memcpy(&b, &x, sizeof b);
    return b;
}
static void ba_vector(FILE *f, const float *v, int n)
{
    for (int k = 0; k < n; ++k)
        fprintf(f, " %" PRIu32, ba_bits(v[k]));
}
static float ba_forward(float *fwd, float sum, float rec, const float em[2],
                        const unsigned char *mm, int K)
{
    float shift = rec / K, scale = (1.0f - rec) / sum;
    sum = 0.0f;
    for (int k = 0; k < K; ++k)
    {
        fwd[k] = em[mm[k]] * (scale * fwd[k] + shift);
        sum += fwd[k];
    }
    return sum;
}
static void ba_backward(float *bwd, float rec, const float em[2],
                        const unsigned char *mm, int K)
{
    float sum = 0.0f;
    for (int k = 0; k < K; ++k)
    {
        bwd[k] *= em[mm[k]];
        sum += bwd[k];
    }
    float shift = rec / K, scale = (1.0f - rec) / sum;
    for (int k = 0; k < K; ++k)
        bwd[k] = scale * bwd[k] + shift;
}
static ba_work ba_create(const ba_input *in)
{
    ba_work w = {0};
    w.types = ba_alloc(in->C, sizeof(int));
    w.a = ba_alloc(in->M, sizeof(int));
    w.b = ba_alloc(in->M, sizeof(int));
    w.unph = ba_alloc(in->C, sizeof(int));
    w.mm = ba_alloc((size_t)3 * in->C * in->K, 1);
    w.rows = ba_alloc((size_t)3 * in->C, sizeof(unsigned char *));
    w.fwd = ba_alloc((size_t)3 * in->K, sizeof(float));
    w.bwd = ba_alloc((size_t)3 * in->K, sizeof(float));
    for (int c = 0; c < in->C; ++c)
    {
        w.nmiss += in->types[c] < 2;
        w.nunph += in->types[c] == 4;
    }
    w.bmiss1 = ba_alloc((size_t)w.nmiss * in->K, sizeof(float));
    w.bmiss2 = ba_alloc((size_t)w.nmiss * in->K, sizeof(float));
    w.bhet1 = ba_alloc((size_t)w.nunph * in->K, sizeof(float));
    w.bhet2 = ba_alloc((size_t)w.nunph * in->K, sizeof(float));
    return w;
}
static void ba_reset(const ba_input *in, ba_work *w)
{
    memcpy(w->types, in->types, (size_t)in->C * sizeof(int));
    memcpy(w->a, in->a, (size_t)in->M * sizeof(int));
    memcpy(w->b, in->b, (size_t)in->M * sizeof(int));
    memcpy(w->mm, in->mm, (size_t)3 * in->C * in->K);
    int u = 0;
    for (int c = 0; c < in->C; ++c)
        if (in->types[c] == 4)
            w->unph[u++] = c;
    for (int r = 0; r < 3 * in->C; ++r)
        w->rows[r] = w->mm + (size_t)r * in->K;
    w->swap = w->swaps = 0;
}
static void ba_step(const ba_input *in, ba_work *w, int c, int backward,
                    FILE *trace)
{
    int start = c ? in->ends[c - 1] : 0;
    float e = (in->ends[c] - start) * in->mismatch;
    if (e >= 0.5f)
        e = 0.5f;
    float em[2] = {1.0f - e, e};
    for (int t = 0; t < 3; ++t)
    {
        if (backward)
            ba_backward(w->bwd + t * in->K, in->recomb[c], em,
                        w->rows[t * in->C + c], in->K);
        else
            w->sums[t] =
                ba_forward(w->fwd + t * in->K, w->sums[t], in->recomb[c], em,
                           w->rows[t * in->C + c], in->K);
    }
    if (trace)
    {
        fprintf(trace, "%c %d", backward ? 'B' : 'F', backward ? c - 1 : c);
        if (!backward)
            ba_vector(trace, w->sums, 3);
        ba_vector(trace, backward ? w->bwd : w->fwd, 3 * in->K);
        fputc('\n', trace);
    }
}
static void ba_phase_het(const ba_input *in, ba_work *w, int c, int u,
                         FILE *trace)
{
    float *b1 = w->bhet1 + (size_t)u * in->K,
          *b2 = w->bhet2 + (size_t)u * in->K;
    float p11 = 0, p12 = 0, p21 = 0, p22 = 0;
    for (int k = 0; k < in->K; ++k)
    {
        p11 += w->fwd[in->K + k] * b1[k];
        p12 += w->fwd[in->K + k] * b2[k];
        p21 += w->fwd[2 * in->K + k] * b1[k];
        p22 += w->fwd[2 * in->K + k] * b2[k];
    }
    float num = p11 * p22, den = p12 * p21;
    int last = w->swap;
    w->swap = num < den;
    if (last != w->swap)
        ++w->swaps;
    if (!in->burnin &&
        (num >= den * in->threshold || (w->swap && den >= num * in->threshold)))
        w->types[c] = 3;
    if (trace)
        fprintf(trace, "H %d %" PRIu32 " %" PRIu32 " %d %d %d\n", c,
                ba_bits(num), ba_bits(den), w->swap, w->swaps, w->types[c]);
}
static void ba_impute(const ba_input *in, ba_work *w, int c, int i, FILE *trace)
{
    float *p1 = w->bmiss1 + (size_t)i * in->K,
          *p2 = w->bmiss2 + (size_t)i * in->K;
    if (w->swap)
    {
        float *tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
    float a1[2] = {0, 0}, a2[2] = {0, 0};
    for (int k = 0; k < in->K; ++k)
    {
        p1[k] *= w->fwd[in->K + k];
        p2[k] *= w->fwd[2 * in->K + k];
    }
    for (int k = 0; k < in->K; ++k)
    {
        int al = in->refs[(size_t)i * in->K + k];
        a1[al] += p1[k];
        a2[al] += p2[k];
    }
    int m = c ? in->ends[c - 1] : 0;
    if (w->types[c] == 0)
    {
        w->a[m] = a1[1] > a1[0];
        w->b[m] = a2[1] > a2[0];
    }
    else
    {
        int x = w->a[m], y = w->b[m];
        float no = a1[x] * a2[y], yes = a1[y] * a2[x];
        if (yes > no)
        {
            w->a[m] = y;
            w->b[m] = x;
            if (yes >= in->threshold * no)
                w->types[c] = 3;
        }
        else if (no >= in->threshold * yes)
            w->types[c] = 3;
    }
    if (trace)
    {
        fprintf(trace, "I %d", c);
        ba_vector(trace, a1, 2);
        ba_vector(trace, a2, 2);
        fprintf(trace, " %d %d %d\n", w->a[m], w->b[m], w->types[c]);
    }
}
static void ba_run(const ba_input *in, ba_work *w, FILE *trace)
{
    int K = in->K, C = in->C, miss = w->nmiss - 1, unph = w->nunph - 1;
    size_t bytes = (size_t)K * sizeof(float);
    for (int k = 0; k < 3 * K; ++k)
        w->bwd[k] = 1.0f / K;
    if (w->types[C - 1] < 2)
    {
        memcpy(w->bmiss1 + (size_t)miss * K, w->bwd, bytes);
        memcpy(w->bmiss2 + (size_t)miss * K, w->bwd, bytes);
        --miss;
    }
    for (int c = C - 2; c >= 0; --c)
    {
        ba_step(in, w, c + 1, 1, trace);
        if (w->types[c] < 2)
        {
            memcpy(w->bmiss1 + (size_t)miss * K, w->bwd + K, bytes);
            memcpy(w->bmiss2 + (size_t)miss * K, w->bwd + 2 * K, bytes);
            --miss;
        }
        if (w->types[c + 1] == 4)
        {
            memcpy(w->bhet1 + (size_t)unph * K, w->bwd + K, bytes);
            memcpy(w->bhet2 + (size_t)unph * K, w->bwd + 2 * K, bytes);
            memcpy(w->bwd + K, w->bwd, bytes);
            memcpy(w->bwd + 2 * K, w->bwd, bytes);
            --unph;
        }
    }
    for (int k = 0; k < 3 * K; ++k)
        w->fwd[k] = 1.0f / K;
    for (int t = 0; t < 3; ++t)
        w->sums[t] = 1.0f;
    miss = unph = 0;
    for (int c = 0; c < C; ++c)
    {
        if (w->types[c] == 4)
        {
            ba_phase_het(in, w, c, unph++, trace);
            if (w->swap)
            {
                int end = unph < w->nunph ? w->unph[unph] : C;
                for (int j = c; j < end; ++j)
                {
                    unsigned char *tmp = w->rows[C + j];
                    w->rows[C + j] = w->rows[2 * C + j];
                    w->rows[2 * C + j] = tmp;
                }
                for (int m = c ? in->ends[c - 1] : 0; m < in->ends[end - 1];
                     ++m)
                {
                    int tmp = w->a[m];
                    w->a[m] = w->b[m];
                    w->b[m] = tmp;
                }
            }
            memcpy(w->fwd + K, w->fwd, bytes);
            memcpy(w->fwd + 2 * K, w->fwd, bytes);
            w->sums[1] = w->sums[2] = w->sums[0];
        }
        ba_step(in, w, c, 0, trace);
        if (w->types[c] < 2)
            ba_impute(in, w, c, miss++, trace);
    }
}
static void ba_destroy(ba_work *w)
{
    free(w->types);
    free(w->a);
    free(w->b);
    free(w->unph);
    free(w->mm);
    free(w->rows);
    free(w->fwd);
    free(w->bwd);
    free(w->bmiss1);
    free(w->bmiss2);
    free(w->bhet1);
    free(w->bhet2);
}

typedef struct
{
    int M, N, it, burnin, iterations, threads, max_states;
    int T;
    uint64_t seed;
    float step, buffer, mismatch, threshold;
    double *pos;
    int *bp, *gt;
    uint64_t *haps;
    hi_clusters *clusters;
    float *recomb;
} it_input;
typedef struct
{
    int S, per_batch, batches, overlap, backoff, candidates;
    int *ends, *codes, *sizes, *donors;
} it_donors;
static void it_ints(FILE *f, const char *tag, int index, const int *v, int n)
{
    fprintf(f, "%s", tag);
    if (index >= 0)
        fprintf(f, " %d", index);
    for (int j = 0; j < n; ++j)
        fprintf(f, " %d", v[j]);
    fputc('\n', f);
}
static it_donors it_find_donors(const it_input *in, const ibs_result *ibs,
                                FILE *trace)
{
    int M = in->M, H = 2 * in->N, words = (M + 63) / 64;
    it_donors out = {0};
    out.ends = ba_alloc(M, sizeof(int));
    for (int end = 0; end < M;)
    {
        double boundary = in->pos[end] + in->step;
        ++end;
        while (end < M && in->pos[end] < boundary)
            ++end;
        out.ends[out.S++] = end;
    }
    out.per_batch = (out.S + in->threads - 1) / in->threads;
    out.batches = (out.S + out.per_batch - 1) / out.per_batch;
    out.overlap = (int)rint((double)(in->buffer / in->step));
    out.backoff = (int)rint((double)(0.3f / in->step));
    int candidates = 100;
    if (in->it >= in->burnin)
    {
        double remaining = in->burnin + in->iterations - in->it;
        double proportion = remaining / in->iterations;
        candidates = (int)floor(proportion * 90 + 0.5);
        if (candidates < 5)
            candidates = 5;
    }
    if (in->it >= in->burnin + in->iterations)
        candidates = 10;
    out.candidates = candidates < H ? candidates : H;
    out.codes = ba_alloc((size_t)out.S * H, sizeof(int));
    out.sizes = ba_alloc(out.S, sizeof(int));
    out.donors = ba_alloc((size_t)out.S * H, sizeof(int));
    uint32_t *hash = ba_alloc(H, sizeof(uint32_t));
    for (int s = 0, start = 0; s < out.S; start = out.ends[s++])
    {
        int *code = out.codes + (size_t)s * H;
        for (int h = 0; h < H; ++h)
        {
            uint32_t value = 0;
            for (int m = start; m < out.ends[s]; ++m)
                if ((in->haps[(size_t)h * words + (m >> 6)] >> (m & 63)) & 1)
                    value ^= UINT32_C(1) << (m % 32);
            hash[h] = value;
            int previous = 0;
            while (previous < h && hash[previous] != value)
                ++previous;
            code[h] = previous < h ? code[previous] : out.sizes[s]++;
        }
    }
    free(hash);
    if (trace)
    {
        fprintf(trace, "SCHED %d %d %d %d %d %d\n", out.S, out.per_batch,
                out.batches, out.overlap, out.backoff, out.candidates);
        it_ints(trace, "ENDS", -1, out.ends, out.S);
        for (int s = 0; s < out.S; ++s)
            it_ints(trace, "CODE", s, out.codes + (size_t)s * H, H);
    }
    pbwt_workspace w;
    if (pbwt_alloc(&w, H))
        st_failure(ST_PHASE_PBWT);
    int *a = ba_alloc(H, sizeof(int)), *d = ba_alloc(H + 1, sizeof(int));
    int backward = (in->it & 1) == 0;
    for (int batch = 0; batch < out.batches; ++batch)
    {
        int start = batch * out.per_batch, end = start + out.per_batch;
        if (end > out.S)
            end = out.S;
        int first = backward ? end + out.overlap : start - out.overlap;
        if (first < 0)
            first = 0;
        if (first > out.S)
            first = out.S;
        int initial = backward ? first - 1 : first;
        for (int h = 0; h < H; ++h)
            a[h] = h;
        for (int h = 0; h <= H; ++h)
            d[h] = initial;
        for (int s = backward ? first - 1 : first;
             backward ? s >= start : s < end; s += backward ? -1 : 1)
        {
            pbwt_update(&w, out.codes + (size_t)s * H, out.sizes[s], s,
                        backward, a, d);
            if (s >= start && s < end)
            {
                d[0] = d[H] = s + (backward ? -2 : 2);
                choose_donors(H, 2 * in->T, out.candidates, s, backward,
                              in->seed + (uint64_t)s, a, d, ibs,
                              s ? out.ends[s - 1] : 0, out.ends[s] - 1,
                              out.donors + (size_t)s * H);
            }
        }
    }
    pbwt_free(&w);
    free(a);
    free(d);
    if (trace)
        for (int s = 0; s < out.S; ++s)
            it_ints(trace, "DONOR", s, out.donors + (size_t)s * H, H);
    return out;
}

#undef it_find_donors

typedef struct
{
    const uint64_t *haps;
    const int *ends;
    int H, words, begin, end;
    int *codes, *sizes;
} st_code_task;
static uint32_t st_hash_slot(uint32_t x)
{
    x ^= x >> 16;
    x *= UINT32_C(0x7feb352d);
    x ^= x >> 15;
    x *= UINT32_C(0x846ca68b);
    return x ^ (x >> 16);
}
static void st_codes_worker(void *data)
{
    st_code_task *t = data;
    size_t capacity = 1;
    while (capacity < (size_t)2 * t->H)
        capacity <<= 1;
    uint32_t *keys = ba_alloc(capacity, sizeof(uint32_t));
    int *first = ba_alloc(capacity, sizeof(int));
    for (int s = t->begin; s < t->end; ++s)
    {
        memset(first, 0, capacity * sizeof(int));
        int start = s ? t->ends[s - 1] : 0;
        int *code = t->codes + (size_t)s * t->H;
        for (int h = 0; h < t->H; ++h)
        {
            uint32_t value = 0;
            for (int m = start; m < t->ends[s]; ++m)
                if ((t->haps[(size_t)h * t->words + (m >> 6)] >> (m & 63)) & 1)
                    value ^= UINT32_C(1) << (m % 32);
            size_t slot = st_hash_slot(value) & (capacity - 1);
            while (first[slot] && keys[slot] != value)
                slot = (slot + 1) & (capacity - 1);
            if (first[slot])
                code[h] = code[first[slot] - 1];
            else
            {
                first[slot] = h + 1;
                keys[slot] = value;
                code[h] = t->sizes[s]++;
            }
        }
    }
    free(keys);
    free(first);
}
static void st_build_codes(const uint64_t *haps, int H, int words, int S,
                           const int *ends, int *codes, int *sizes, int workers)
{
    if (workers > S)
        workers = S;
    if (workers < 1)
        workers = 1;
#ifdef MT
#pragma omp parallel num_threads(workers) if (workers > 1)
#endif
    {
        int tid = omp_get_thread_num(), team = omp_get_num_threads();
        st_code_task t = {haps,  ends,           H,
                          words, S * tid / team, S * (tid + 1) / team,
                          codes, sizes};
        st_worker_run(st_codes_worker, &t, 3);
    }
    st_workers_check();
}

typedef struct
{
    const it_input *in;
    const ibs_result *ibs;
    it_donors *out;
    const int *a, *d, *steps;
    int begin, end, backward, stride;
} st_donor_task;
static void st_donors_worker(void *data)
{
    st_donor_task *t = data;
    int H = 2 * t->in->N;
    for (int j = t->begin; j < t->end; j += t->stride)
    {
        int s = t->steps[j];
        choose_donors(H, 2 * t->in->T, t->out->candidates, s, t->backward,
                      t->in->seed + (uint64_t)s, t->a + (size_t)j * H,
                      t->d + (size_t)j * (H + 1), t->ibs,
                      s ? t->out->ends[s - 1] : 0, t->out->ends[s] - 1,
                      t->out->donors + (size_t)s * H);
    }
}
static void st_donors_block(const it_input *in, const ibs_result *ibs,
                            it_donors *out, int *a, int *d, int *steps,
                            int count, int backward, int workers)
{
    if (workers > count)
        workers = count;
#ifdef MT
#pragma omp parallel num_threads(workers) if (workers > 1)
#endif
    {
        int tid = omp_get_thread_num(), team = omp_get_num_threads();
        st_donor_task t = {in,    ibs, out,   a,        d,
                           steps, tid, count, backward, team};
        st_worker_run(st_donors_worker, &t, 4);
    }
    st_workers_check();
}
static it_donors it_find_donors(const it_input *in, const ibs_result *ibs,
                                FILE *trace)
{
    if (trace)
        return st_serial_find_donors(in, ibs, trace);
    int M = in->M, H = 2 * in->N, words = (M + 63) / 64,
        workers = st_phase_threads;
    it_donors out = {0};
    out.ends = ba_alloc(M, sizeof(int));
    for (int end = 0; end < M;)
    {
        double boundary = in->pos[end] + in->step;
        ++end;
        while (end < M && in->pos[end] < boundary)
            ++end;
        out.ends[out.S++] = end;
    }

    out.per_batch = (out.S + in->threads - 1) / in->threads;
    out.batches = (out.S + out.per_batch - 1) / out.per_batch;
    out.overlap = (int)rint((double)(in->buffer / in->step));
    out.backoff = (int)rint((double)(0.3f / in->step));
    int candidates = 100;
    if (in->it >= in->burnin)
    {
        double remaining = in->burnin + in->iterations - in->it;
        double proportion = remaining / in->iterations;
        candidates = (int)floor(proportion * 90 + 0.5);
        if (candidates < 5)
            candidates = 5;
    }
    if (in->it >= in->burnin + in->iterations)
        candidates = 10;
    out.candidates = candidates < H ? candidates : H;
    out.codes = ba_alloc((size_t)out.S * H, sizeof(int));
    out.sizes = ba_alloc(out.S, sizeof(int));
    out.donors = ba_alloc((size_t)out.S * H, sizeof(int));
    st_build_codes(in->haps, H, words, out.S, out.ends, out.codes, out.sizes,
                   workers);
    pbwt_workspace w;
    if (pbwt_alloc(&w, H))
        st_failure(ST_PHASE_PBWT);
    int *a = ba_alloc(H, sizeof(int)), *d = ba_alloc(H + 1, sizeof(int));
    int backward = (in->it & 1) == 0;
    enum
    {
        BLOCK = 128
    };
    int *snap_a = ba_alloc((size_t)BLOCK * H, sizeof(int));
    int *snap_d = ba_alloc((size_t)BLOCK * (H + 1), sizeof(int));
    int steps[BLOCK], count = 0;
    for (int batch = 0; batch < out.batches; ++batch)
    {
        int start = batch * out.per_batch, end = start + out.per_batch;
        if (end > out.S)
            end = out.S;
        int first = backward ? end + out.overlap : start - out.overlap;
        if (first < 0)
            first = 0;
        if (first > out.S)
            first = out.S;
        int initial = backward ? first - 1 : first;
        for (int h = 0; h < H; ++h)
            a[h] = h;
        for (int h = 0; h <= H; ++h)
            d[h] = initial;
        for (int s = backward ? first - 1 : first;
             backward ? s >= start : s < end; s += backward ? -1 : 1)
        {
            pbwt_update(&w, out.codes + (size_t)s * H, out.sizes[s], s,
                        backward, a, d);
            if (s >= start && s < end)
            {
                d[0] = d[H] = s + (backward ? -2 : 2);
                if (workers == 1)
                {
                    choose_donors(H, 2 * in->T, out.candidates, s, backward,
                                  in->seed + (uint64_t)s, a, d, ibs,
                                  s ? out.ends[s - 1] : 0, out.ends[s] - 1,
                                  out.donors + (size_t)s * H);
                }
                else
                {
                    memcpy(snap_a + (size_t)count * H, a,
                           (size_t)H * sizeof(int));
                    memcpy(snap_d + (size_t)count * (H + 1), d,
                           (size_t)(H + 1) * sizeof(int));
                    steps[count++] = s;
                    if (count == BLOCK)
                    {
                        st_donors_block(in, ibs, &out, snap_a, snap_d, steps,
                                        count, backward, workers);
                        count = 0;
                    }
                }
            }
        }
    }
    if (count)
        st_donors_block(in, ibs, &out, snap_a, snap_d, steps, count, backward,
                        workers);
    pbwt_free(&w);
    free(a);
    free(d);
    free(snap_a);
    free(snap_d);
    return out;
}

typedef struct
{
    it_input input;
    float *distance, intensity, initial_lr;
    int em;
} sq_input;
typedef struct
{
    uint64_t checksum;
    int iterations;
    float mismatch, intensity;
} sq_result;
typedef struct
{
    uint64_t checksum, unphased, swaps;
} sq_counts;
typedef struct
{
    it_input *in;
    it_donors donor;
    cs_input ci;
    uint64_t *next;
    FILE *trace;
    int begin, end, stride;
    sq_counts *counts;
} sq_task;

static void sq_task_run(void *data)
{
    sq_task *task = data;
    it_input *in = task->in;
    it_donors donor = task->donor;
    cs_input ci = task->ci;
    uint64_t *next = task->next;
    FILE *trace = task->trace;
    int M = in->M, H = 2 * in->N, words = (M + 63) / 64;
    uint64_t checksum = 0, unphased_count = 0, swap_count = 0;
    ci.donors = ba_alloc((size_t)2 * donor.S, sizeof(int));
    cs_work cw = cs_create(&ci);
    for (int sample = task->begin; sample < task->end; sample += task->stride)
    {
        hi_clusters cl = {in->clusters[sample].n, in->clusters[sample].end,
                          ba_alloc(in->clusters[sample].n, sizeof(int))};
        memcpy(cl.type, in->clusters[sample].type, (size_t)cl.n * sizeof(int));
        if (in->threshold < 50.0f)
            hi_mask(&cl, in->bp);
        int missing = 0, unph = 0;
        for (int c = 0; c < cl.n; ++c)
        {
            missing += cl.type[c] < 2;
            unph += cl.type[c] == 4;
        }
        const uint64_t *h1 = in->haps + (size_t)2 * sample * words,
                       *h2 = h1 + words;
        int *a = ba_alloc(M, sizeof(int)), *b = ba_alloc(M, sizeof(int));
        for (int m = 0; m < M; ++m)
        {
            a[m] = (h1[m >> 6] >> (m & 63)) & 1;
            b[m] = (h2[m >> 6] >> (m & 63)) & 1;
        }
        if (missing || unph)
        {
            for (int s = 0; s < donor.S; ++s)
                for (int j = 0; j < 2; ++j)
                    ci.donors[2 * s + j] =
                        donor.donors[(size_t)s * H + 2 * sample + j];
            int K = cs_run(&ci, &cw, sample, in->seed, NULL);
            ba_input bi = {M,
                           cl.n,
                           K,
                           in->it < in->burnin,
                           in->mismatch,
                           in->threshold,
                           cl.end,
                           cl.type,
                           a,
                           b,
                           NULL,
                           NULL,
                           NULL};
            bi.refs = ba_alloc((size_t)missing * K, sizeof(int));
            bi.recomb = ba_alloc(cl.n, sizeof(float));
            bi.mm = ba_alloc((size_t)3 * cl.n * K, 1);
            hi_recomb(&cl, in->recomb, bi.recomb);
            hi_mismatches(&cl, K, words, h1, h2, cw.comp, bi.mm, bi.refs);
            ba_work bw = ba_create(&bi);
            ba_reset(&bi, &bw);
            ba_run(&bi, &bw, NULL);
            memcpy(a, bw.a, (size_t)M * sizeof(int));
            memcpy(b, bw.b, (size_t)M * sizeof(int));
            memcpy(cl.type, bw.types, (size_t)cl.n * sizeof(int));
            unphased_count += (uint64_t)unph;
            swap_count += (uint64_t)bw.swaps;
            ba_destroy(&bw);
            free(bi.refs);
            free(bi.recomb);
            free(bi.mm);
        }
        if (trace)
        {
            it_ints(trace, "TYPES", sample, cl.type, cl.n);
            it_ints(trace, "A", sample, a, M);
            it_ints(trace, "B", sample, b, M);
        }
        for (int m = 0; m < M; ++m)
            checksum += (uint64_t)(sample + 1) * (m + 1) * (a[m] + 3 * b[m]);
        for (int c = 0; c < cl.n; ++c)
            checksum += (uint64_t)(sample + 1) * (c + 1) * cl.type[c];
        for (int m = 0; m < M; ++m)
        {
            next[(size_t)2 * sample * words + (m >> 6)] |= (uint64_t)a[m]
                                                           << (m & 63);
            next[(size_t)(2 * sample + 1) * words + (m >> 6)] |= (uint64_t)b[m]
                                                                 << (m & 63);
        }
        memcpy(in->clusters[sample].type, cl.type, (size_t)cl.n * sizeof(int));
        free(a);
        free(b);
        free(cl.type);
    }
    cs_destroy(&cw);
    free(ci.donors);
    *task->counts = (sq_counts){checksum, unphased_count, swap_count};
}
static uint64_t sq_phase(it_input *in, it_donors donor, uint64_t *unphased,
                         uint64_t *swaps, FILE *trace)
{
    int M = in->M, H = 2 * in->N, words = (M + 63) / 64;
    cs_input ci = {M,     H,    donor.S, in->max_states, 200,
                   words, NULL, NULL,    in->haps};
    int min_steps = (int)ceil((double)(1.0f / in->step));
    if (min_steps > ci.min_steps)
        ci.min_steps = min_steps;
    ci.starts = ba_alloc(donor.S, sizeof(int));
    for (int s = 1; s < donor.S; ++s)
        ci.starts[s] = donor.ends[s - 1];
    uint64_t checksum = 0, unphased_count = 0, swap_count = 0;
    uint64_t *next = ba_alloc((size_t)H * words, sizeof(uint64_t));

    memcpy(next + (size_t)2 * in->T * words,
           in->haps + (size_t)2 * in->T * words,
           (size_t)2 * (in->N - in->T) * words * sizeof(uint64_t));
    int workers = trace ? 1 : st_phase_threads;
    if (workers > in->T)
        workers = in->T;
    sq_counts *counts = ba_alloc(workers, sizeof(sq_counts));

#ifdef MT
#pragma omp parallel num_threads(workers) if (workers > 1)
#endif
    {
        int tid = 0, team = 1;
#ifdef MT
        tid = omp_get_thread_num();
        team = omp_get_num_threads();
#endif
        sq_task task = {in,  donor, ci,   next,        trace,
                        tid, in->T, team, counts + tid};
        st_worker_run(sq_task_run, &task, 2);
    }
    st_workers_check();
    for (int w = 0; w < workers; ++w)
    {
        checksum += counts[w].checksum;
        unphased_count += counts[w].unphased;
        swap_count += counts[w].swaps;
    }
    free(counts);
    *unphased = unphased_count;
    *swaps = swap_count;
    free(ci.starts);
    memcpy(in->haps, next, (size_t)H * words * sizeof(uint64_t));
    free(next);
    return checksum;
}

static float sq_threshold(const it_input *in, float initial)
{
    if (in->it < in->burnin)
        return INFINITY;
    int last = in->iterations - 1;
    if (in->it == last + in->burnin)
        return 1.0f;
    double exponent = (double)(last - (in->it - in->burnin)) / last;
    return (float)(4.0 * pow((double)initial / 4.0, exponent));
}
static void sq_check_reference(const it_input *in)
{
    int H = 2 * in->N, words = (in->M + 63) / 64;
    for (int h = 2 * in->T; h < H; ++h)
        for (int m = 0; m < in->M; ++m)
            if ((int)((in->haps[(size_t)h * words + (m >> 6)] >> (m & 63)) &
                      1) != in->gt[(size_t)m * H + h])
                st_failure(ST_PHASE_REFERENCE);
}
static sq_result sq_run(sq_input *seq, const ibs_result *ibs, FILE *trace)
{
    it_input *in = &seq->input;
    uint64_t base_seed = in->seed - (uint64_t)in->it;
    sq_result result = {0, 0, in->mismatch, seq->intensity};
    while (in->it < in->burnin + in->iterations)
    {
        in->seed = base_seed + (uint64_t)in->it;
        in->threshold = sq_threshold(in, seq->initial_lr);
        if (trace)
            fprintf(trace, "ITER %d %u\n", in->it,
                    pe_float_bits(in->threshold));
        it_donors donor = it_find_donors(in, ibs, trace);
        if (seq->em && (in->it == 0 || in->it < in->burnin))
        {
            int min_steps = (int)ceil((double)(1.0f / in->step));
            if (min_steps < 200)
                min_steps = 200;
            int *starts = ba_alloc(donor.S, sizeof(int));
            for (int s = 1; s < donor.S; ++s)
                starts[s] = donor.ends[s - 1];
            er_input est = {{in->M, 2 * in->N, donor.S, in->max_states,
                             min_steps, (in->M + 63) / 64, starts, NULL,
                             in->haps},
                            in->T,
                            in->it == 0,
                            in->seed,
                            donor.donors,
                            seq->distance,
                            in->mismatch,
                            seq->intensity};
            er_result e = er_run(&est, NULL);
            in->mismatch = e.mismatch;
            seq->intensity = e.intensity;
            for (int m = 0; m < in->M; ++m)
                in->recomb[m] =
                    (float)-expm1(-(double)e.intensity * seq->distance[m]);
            free(starts);
        }
        if (trace)
        {
            fprintf(trace, "PARAM %u %u\nRECOMB", pe_float_bits(in->mismatch),
                    pe_float_bits(seq->intensity));
            for (int m = 0; m < in->M; ++m)
                fprintf(trace, " %u", pe_float_bits(in->recomb[m]));
            fputc('\n', trace);
        }
        uint64_t unphased, swaps;
        result.checksum = sq_phase(in, donor, &unphased, &swaps, trace);
        sq_check_reference(in);
        if (trace)
            fprintf(trace, "COUNTS %llu %llu\n", (unsigned long long)unphased,
                    (unsigned long long)swaps);
        free(donor.ends);
        free(donor.codes);
        free(donor.sizes);
        free(donor.donors);
        ++result.iterations;
        ++in->it;
        double rate = unphased ? (double)swaps / unphased : NAN;
        if (in->it < in->burnin && rate <= 0.01)
            in->it = in->burnin;
    }
    result.mismatch = in->mismatch;
    result.intensity = seq->intensity;
    return result;
}

static uint64_t rng;
static uint32_t next_bits(unsigned bits)
{
    rng = (rng * UINT64_C(0x5deece66d) + 11) & ((UINT64_C(1) << 48) - 1);
    return (uint32_t)(rng >> (48 - bits));
}
static uint32_t next_int(uint32_t bound)
{
    if ((bound & (bound - 1)) == 0)
        return (uint32_t)(((uint64_t)bound * next_bits(31)) >> 31);
    uint32_t bits, val;
    do
    {
        bits = next_bits(31);
        val = bits % bound;
    } while ((uint64_t)bits - val + (bound - 1) > INT32_MAX);
    return val;
}
static int vote(int h, int a1, int a2, const int *order, const int *rank,
                const int *allele, const int *unph, int H)
{
    int r = rank[h], total = 0;
    for (int side = 0; side < 2; ++side)
    {
        int j = r + (side ? 1 : -1);
        if (j < 0 || j >= H)
            continue;
        int d = order[j];
        if (!unph[d / 2])
            total += (allele[d] == a1) - (allele[d] == a2);
    }
    return total;
}
static int fill(int h, const int *order, const int *rank, const int *a,
                const int *unph, int H)
{
    int r = rank[h], left = -1, right = -1;
    if (r > 0 && !unph[order[r - 1] / 2])
        left = a[order[r - 1]];
    if (r + 1 < H && !unph[order[r + 1] / 2])
        right = a[order[r + 1]];
    if (left >= 0 && (left == right || right < 0))
        return left;
    if (left < 0 && right >= 0)
        return right;
    return -1;
}
static void partial(int *a, int *unph, const int *order, const int *rank, int N)
{
    int threshold = 2, change = 1;
    while (threshold > 0 || change)
    {
        change = 0;
        for (int s = 0; s < N; ++s)
        {
            int h = 2 * s;
            if (unph[s])
            {
                int count =
                    vote(h, a[h], a[h + 1], order, rank, a, unph, 2 * N) +
                    vote(h + 1, a[h + 1], a[h], order, rank, a, unph, 2 * N);
                if (count >= threshold)
                {
                    unph[s] = 0;
                    change = 1;
                }
                else if (count <= -threshold)
                {
                    int t = a[h];
                    a[h] = a[h + 1];
                    a[h + 1] = t;
                    unph[s] = 0;
                    change = 1;
                }
            }
            else
            {
                for (int side = 0; side < 2; ++side)
                    if (a[h + side] < 0)
                    {
                        a[h + side] =
                            fill(h + side, order, rank, a, unph, 2 * N);
                        change |= a[h + side] >= 0;
                    }
            }
        }
        if (!change)
            --threshold;
    }
}

typedef struct
{
    int M, N, burn, iterations, K, em;
    uint64_t seed;
    float rare, scale, ne, buffer, lr;
    double *pos;
    int *bp, *gt;
    int T;
} np_input;
static int np_cmp(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
static void np_free(sq_input *seq)
{
    it_input *in = &seq->input;
    for (int s = 0; s < in->N; ++s)
        hi_free(in->clusters + s);
    free(in->clusters);
    free(in->gt);
    free(in->pos);
    free(in->bp);
    free(in->haps);
    free(in->recomb);
    free(seq->distance);
}

static sq_input fp_select(const np_input *raw, const int *na, int **selected)
{
    int M = raw->M, N = raw->N, H = 2 * N;
    int *ix = ba_alloc(M, sizeof(int));
    int S = 0;
    int max = (int)floor((double)((float)N * raw->rare));
    if (max < 3)
        max = 3;
    for (int m = 0; m < M; ++m)
    {
        int *cnt = ba_alloc(na[m], sizeof(int));
        for (int s = 0; s < N; ++s)
        {
            int a = raw->gt[(size_t)m * H + 2 * s],
                b = raw->gt[(size_t)m * H + 2 * s + 1];
            if (a >= 0)
                ++cnt[a];
            if (b >= 0 && b != a)
                ++cnt[b];
        }
        int high = 0;
        for (int a = 0; a < na[m]; ++a)
            high += cnt[a] > max;
        free(cnt);
        if (high > 1)
            ix[S++] = m;
    }
    if (S < 2 || (float)S > 0.75f * (float)M)
    {
        S = M;
        for (int m = 0; m < M; ++m)
            ix[m] = m;
    }
    for (int j = 0; j < S; ++j)
        if (na[ix[j]] > 2)
        {
            st_failure(ST_PHASE_ALLELES);
        }
    sq_input out = {0};
    it_input *in = &out.input;
    in->M = S;
    in->N = N;
    in->T = raw->T;
    in->burnin = raw->burn;
    in->iterations = raw->iterations;
    in->threads = 1;
    in->max_states = raw->K;
    in->buffer = raw->buffer;
    j_random jr;
    jr_seed(&jr, raw->seed);
    uint64_t high = (uint64_t)(int64_t)(int32_t)jr_bits(&jr, 32);
    int64_t low = (int32_t)jr_bits(&jr, 32);
    in->seed = (high << 32) + (uint64_t)low;
    in->pos = ba_alloc(S, sizeof(double));
    in->bp = ba_alloc(S, sizeof(int));
    in->gt = ba_alloc((size_t)S * H, sizeof(int));
    in->haps = ba_alloc((size_t)H * ((S + 63) / 64), sizeof(uint64_t));
    in->clusters = ba_alloc(N, sizeof(hi_clusters));
    in->recomb = ba_alloc(S, sizeof(float));
    out.distance = ba_alloc(S, sizeof(float));
    out.initial_lr = raw->lr;
    out.em = raw->em;
    double *diff = ba_alloc(S - 1, sizeof(double));
    for (int m = 0; m < S; ++m)
    {
        in->pos[m] = raw->pos[ix[m]];
        in->bp[m] = raw->bp[ix[m]];
        memcpy(in->gt + (size_t)m * H, raw->gt + (size_t)ix[m] * H,
               (size_t)H * sizeof(int));
        if (m)
        {
            out.distance[m] = (float)(in->pos[m] - in->pos[m - 1]);
            diff[m - 1] = in->pos[m] - in->pos[m - 1];
        }
    }
    qsort(diff, S - 1, sizeof(double), np_cmp);
    int n = S - 1;
    float median = 0.5f * (float)(diff[(n - 1) / 2] + diff[n / 2]);
    in->step = raw->scale * median;
    free(diff);
    *selected = ix;
    return out;
}

typedef struct
{
    int H, T, S, per_batch, overlap, backoff, candidates, backward;
    uint64_t seed;
    int *sizes, *codes, *list_start, *lengths, **lists;
    const ibs_result *ibs;
    const int *starts, *ends;
} rp_input;
static int rp_min(int a, int b)
{
    return a < b ? a : b;
}
static int rp_max(int a, int b)
{
    return a > b ? a : b;
}
static int rp_cmp(const void *a, const void *b)
{
    int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}
static int rp_excluded(const rp_input *in, int step, int hap, int other)
{
    return ibs_query(in->ibs, hap / 2, other / 2, in->starts[step],
                     in->ends[step] - 1);
}
static void rp_neighbors(const rp_input *in, int step, const int *a, int *inv,
                         int *prev, int *next, int *scratch)
{
    for (int i = 0; i < in->H; ++i)
    {
        inv[a[i]] = i;
        prev[i] = INT_MIN;
        next[i] = INT_MAX;
    }
    for (int list = in->list_start[step]; list < in->list_start[step + 1];
         ++list)
    {
        int n = in->lengths[list];
        for (int j = 0; j < n; ++j)
            scratch[j] = inv[in->lists[list][j]];
        qsort(scratch, (size_t)n, sizeof(int), rp_cmp);
        for (int j = 1; j < n; ++j)
        {
            int x = scratch[j - 1], y = scratch[j];
            if (x > prev[y])
                prev[y] = x;
            if (y < next[x])
                next[x] = y;
        }
    }
}
static int rp_best(const rp_input *in, int step, int i, const int *a,
                   const int *d, const int *prev, const int *next)
{
    int best_prev = -1, best_next = -1;
    int left = in->backward ? in->S - 1 : 0, right = left;
    int edge = (i + 1) < in->H ? (in->backward ? rp_max(d[i], d[i + 1])
                                               : rp_min(d[i], d[i + 1]))
                               : d[i];
    int limit = in->backward ? rp_max(edge - in->backoff, step)
                             : rp_min(edge + in->backoff, step);
    int pi = prev[i];
    while (pi > INT_MIN && rp_excluded(in, step, a[i], a[pi]))
        pi = prev[pi];
    if (pi > INT_MIN)
    {
        int u = i;
        while (u - 1 != pi && (in->backward ? d[u] >= limit : d[u] <= limit))
        {
            left = in->backward ? rp_min(left, d[u]) : rp_max(left, d[u]);
            --u;
        }
        if (u - 1 == pi && (in->backward ? d[u] >= limit : d[u] <= limit))
        {
            left = in->backward ? rp_min(left, d[u]) : rp_max(left, d[u]);
            best_prev = pi;
        }
    }
    int ni = next[i];
    while (ni < INT_MAX && rp_excluded(in, step, a[i], a[ni]))
        ni = next[ni];
    if (ni < INT_MAX)
    {
        int v = i;
        while (v + 1 != ni &&
               (in->backward ? d[v + 1] >= limit : d[v + 1] <= limit))
        {
            ++v;
            right = in->backward ? rp_min(right, d[v]) : rp_max(right, d[v]);
        }
        if (v + 1 == ni &&
            (in->backward ? d[v + 1] >= limit : d[v + 1] <= limit))
        {
            ++v;
            right = in->backward ? rp_min(right, d[v]) : rp_max(right, d[v]);
            best_next = ni;
        }
    }

    if ((in->backward ? left > right : left < right) && best_prev != -1)
        return best_prev;
    return best_next;
}
static void rp_select(const rp_input *in, int step, const int *a, int *d,
                      const int *prev, const int *next, int *out)
{
    j_random rng;
    jr_seed(&rng, in->seed + (uint64_t)step);
    d[in->H] = step + (in->backward ? -1 : 1);
    for (int i = 0; i < in->H; ++i)
        if (a[i] < in->T)
        {
            int best = rp_best(in, step, i, a, d, prev, next);
            if (best >= 0)
            {
                out[a[i]] = a[best];
                continue;
            }
            int u = i, v = i + 1, left = d[u], right = d[v];
            while (v - u < in->candidates &&
                   (in->backward ? (step <= left || step <= right)
                                 : (left <= step || right <= step)))
            {
                if (in->backward)
                {
                    if (left <= right)
                    {
                        ++v;
                        right = rp_min(d[v], right);
                    }
                    else
                    {
                        --u;
                        left = rp_min(d[u], left);
                    }
                }
                else
                {
                    if (right <= left)
                    {
                        ++v;
                        right = rp_max(d[v], right);
                    }
                    else
                    {
                        --u;
                        left = rp_max(d[u], left);
                    }
                }
            }
            int n = v - u, match = -1;
            if (n > 1)
            {
                int index = u + (int)jr_int(&rng, (uint32_t)n);
                for (int j = 0; j < n && match == -1; ++j)
                {
                    if (!rp_excluded(in, step, a[i], a[index]))
                        match = a[index];
                    if (++index == v)
                        index = u;
                }
            }
            out[a[i]] = match;
        }
}
static void rp_print(FILE *f, const char *name, int step, const int *x, int n)
{
    if (!f)
        return;
    fprintf(f, "%s %d", name, step);
    for (int j = 0; j < n; ++j)
        fprintf(f, " %d", x[j]);
    fputc('\n', f);
}
static void rp_run(const rp_input *in, int *selected, FILE *trace)
{
    int H = in->H;
    int *a = malloc((size_t)H * sizeof(int)),
        *d = malloc((size_t)(H + 1) * sizeof(int));
    int *inv = malloc((size_t)H * sizeof(int)),
        *prev = malloc((size_t)H * sizeof(int)),
        *next = malloc((size_t)H * sizeof(int)),
        *tmp = malloc((size_t)H * sizeof(int));
    pbwt_workspace w;
    if (!a || !d || !inv || !prev || !next || !tmp || pbwt_alloc(&w, H))
        st_failure(ST_PHASE_PBWT);
    for (int start = 0; start < in->S; start += in->per_batch)
    {
        int end = rp_min(start + in->per_batch, in->S);
        int buffer = in->backward ? rp_min(end + in->overlap, in->S)
                                  : rp_max(0, start - in->overlap);
        for (int i = 0; i < H; ++i)
            a[i] = i;
        for (int i = 0; i <= H; ++i)
            d[i] = buffer - (in->backward ? 1 : 0);
        int first = in->backward ? buffer - 1 : buffer,
            stop = in->backward ? start - 1 : end,
            delta = in->backward ? -1 : 1;
        for (int s = first; s != stop; s += delta)
        {
            pbwt_update(&w, in->codes + (size_t)s * H, in->sizes[s], s,
                        in->backward, a, d);
            if (s < start || s >= end)
                continue;
            rp_neighbors(in, s, a, inv, prev, next, tmp);
            int *out = selected + (size_t)s * in->T;
            rp_select(in, s, a, d, prev, next, out);
            rp_print(trace, "A", s, a, H);
            rp_print(trace, "D", s, d, H + 1);
            rp_print(trace, "PREV", s, prev, H);
            rp_print(trace, "NEXT", s, next, H);
            rp_print(trace, "DONOR", s, out, in->T);
        }
    }
    pbwt_free(&w);
    free(a);
    free(d);
    free(inv);
    free(prev);
    free(next);
    free(tmp);
}

typedef struct
{
    int M, H, S, K, min_steps, target;
    uint64_t seed;
    const int *starts, *donors, *alleles;
} ls_input;
typedef struct
{
    int hap, end, next;
} ls_node;
typedef struct
{
    cs_work queue;
    int used, *first, *last, *cursor;
    ls_node *nodes;
    int *ids;
    unsigned char *mismatch;
} ls_work;
static ls_work ls_create(const ls_input *in)
{
    ls_work w = {0};
    w.queue.last = cs_alloc(in->H, sizeof(int));
    w.queue.heap = cs_alloc(in->K, sizeof(cs_segment));
    w.first = cs_alloc(in->K, sizeof(int));
    w.last = cs_alloc(in->K, sizeof(int));
    w.cursor = cs_alloc(in->K, sizeof(int));
    w.nodes = cs_alloc((size_t)2 * in->S + in->K, sizeof(ls_node));
    w.ids = cs_alloc((size_t)in->M * in->K, sizeof(int));
    w.mismatch = cs_alloc((size_t)in->M * in->K, sizeof(unsigned char));
    return w;
}
static void ls_destroy(ls_work *w)
{
    cs_destroy(&w->queue);
    free(w->ids);
    free(w->mismatch);
    free(w->first);
    free(w->last);
    free(w->cursor);
    free(w->nodes);
}
static void ls_append(ls_work *w, cs_segment segment, int end)
{
    int n = w->used++, index = segment.index;
    w->nodes[n] = (ls_node){segment.hap, end, -1};
    if (w->first[index] < 0)
        w->first[index] = n;
    else
        w->nodes[w->last[index]].next = n;
    w->last[index] = n;
}
static void ls_add(const ls_input *in, ls_work *w, int hap, int step)
{
    if (hap < 0)
        return;
    cs_work *q = &w->queue;
    if (q->last[hap] == -103)
    {
        while (q->n && q->heap[0].last != q->last[q->heap[0].hap])
        {
            cs_segment head = cs_poll(q);
            head.last = q->last[head.hap];
            cs_offer(q, head);
        }
        if (q->n == in->K || (q->n && step - q->heap[0].last >= in->min_steps))
        {
            cs_segment head = cs_poll(q);
            int next = in->starts[(head.last + step) >> 1];
            q->last[head.hap] = -103;
            ls_append(w, head, next);
            head.hap = hap;
            head.start = next;
            head.last = step;
            cs_offer(q, head);
        }
        else
        {
            cs_segment head = {hap, 0, step, q->n};
            cs_offer(q, head);
        }
    }
    q->last[hap] = step;
}
static int ls_run(const ls_input *in, ls_work *w)
{
    cs_work *q = &w->queue;
    q->n = 0;
    w->used = 0;
    for (int j = 0; j < in->K; ++j)
        w->first[j] = w->last[j] = -1;
    for (int h = 0; h < in->H; ++h)
        q->last[h] = -103;
    for (int step = 0; step < in->S; ++step)
    {
        ls_add(in, w, in->donors[2 * step], step);
        ls_add(in, w, in->donors[2 * step + 1], step);
    }
    if (!q->n)
    {
        int count = in->H - 2 < in->K ? in->H - 2 : in->K;
        if (count <= 0)
            return 0;
        j_random rng;
        jr_seed(&rng, in->seed + (uint64_t)in->target);
        for (int j = 0; j < count; ++j)
        {
            int h;
            do
            {
                h = (int)jr_int(&rng, (uint32_t)in->H);
            } while ((h >> 1) == (in->target >> 1));

            cs_segment head = {h, 0, 0, j};
            cs_offer(q, head);
        }
    }
    int count = q->n;
    while (q->n)
        ls_append(w, cs_poll(q), in->M);
    for (int j = 0; j < count; ++j)
        w->cursor[j] = w->first[j];
    for (int m = 0; m < in->M; ++m)
    {
        const int *alleles = in->alleles + (size_t)m * in->H;
        for (int j = 0; j < count; ++j)
        {
            size_t index = (size_t)m * in->K + j;
            int cursor = w->cursor[j];

            if (m == w->nodes[cursor].end)
                cursor = w->nodes[cursor].next;
            w->cursor[j] = cursor;
            w->ids[index] = w->nodes[cursor].hap;
            w->mismatch[index] = alleles[w->ids[index]] != alleles[in->target];
        }
    }
    return count;
}

typedef struct
{
    int markers, states;
    float emission[2];
    const float *recomb;
    const unsigned char *mismatch;
} sp_input;
static void sp_forward(const sp_input *in, float *restrict probs)
{
    int K = in->states;
    float last_sum = 0.0f;
    for (int j = 0; j < K; ++j)
    {
        probs[j] = in->emission[in->mismatch[j]];
        last_sum += probs[j];
    }
    for (int m = 1; m < in->markers; ++m)
    {
        size_t row = (size_t)m * K, previous = row - K;
        float p_rec = in->recomb[m];
        float shift = p_rec / K;
        float scale = (1.0f - p_rec) / last_sum;
        last_sum = 0.0f;
        for (int j = 0; j < K; ++j)
        {
            float emission = in->emission[in->mismatch[row + j]];
            probs[row + j] = emission * (scale * probs[previous + j] + shift);
            last_sum += probs[row + j];
        }
    }
}
static void sp_backward(const sp_input *in, float *restrict probs,
                        float *restrict backward)
{
    int K = in->states;
    for (int j = 0; j < K; ++j)
        backward[j] = 1.0f / K;

    for (int m = in->markers - 2; m >= 0; --m)
    {
        size_t row = (size_t)m * K, next = row + K;
        float sum = 0.0f;
        for (int j = 0; j < K; ++j)
        {
            backward[j] *= in->emission[in->mismatch[next + j]];
            sum += backward[j];
        }
        float p_rec = in->recomb[m + 1];
        float scale = (1.0f - p_rec) / sum;
        float shift = p_rec / K;
        sum = 0.0f;
        for (int j = 0; j < K; ++j)
        {
            backward[j] = scale * backward[j] + shift;
            probs[row + j] *= backward[j];
            sum += probs[row + j];
        }
        for (int j = 0; j < K; ++j)
            probs[row + j] /= sum;
    }
}

typedef struct
{
    int a, b, reference;
    float left, right;
} rr_donor;
typedef struct
{
    int alleles, a, b, count[2];
    float weight;
    int *rare;
    rr_donor *donors[2];
} rr_input;
static void rr_probabilities(const rr_input *in, int bit, int impute,
                             float *out)
{
    memset(out, 0, (size_t)in->alleles * sizeof(*out));
    for (int j = 0; j < in->count[bit]; ++j)
    {
        rr_donor d = in->donors[bit][j];
        if (d.a < 0 || d.b < 0)
            continue;
        float probability = in->weight * d.left + (1.0f - in->weight) * d.right;
        if (d.a == d.b || (impute && d.reference))
        {
            out[d.a] += probability;
        }
        else if (!impute)
        {
            int x = in->rare[in->a] && (in->a == d.a || in->a == d.b);
            int y = in->rare[in->b] && (in->b == d.a || in->b == d.b);
            if (x != y)
                out[x ? in->a : in->b] += probability;
        }
        else
        {
            int x = in->rare[d.a], y = in->rare[d.b];

            double wa = x != y ? (x ? 0.55 : 0.45) : 0.5;
            double wb = x != y ? (x ? 0.45 : 0.55) : 0.5;
            out[d.a] = (float)((double)out[d.a] + wa * probability);
            out[d.b] = (float)((double)out[d.b] + wb * probability);
        }
    }
}
static int rr_max(const float *p, int n)
{
    int best = 0;
    for (int j = 1; j < n; ++j)
        if (p[j] > p[best])
            best = j;
    return best;
}
static void rr_run(const rr_input *in, j_random *rng, float *p, float *q,
                   int *a, int *b)
{
    *a = in->a;
    *b = in->b;
    if (*a >= 0 && *b >= 0 && *a == *b)
        return;
    int impute = *a < 0 || *b < 0;
    rr_probabilities(in, 0, impute, p);
    rr_probabilities(in, 1, impute, q);
    if (impute)
    {
        *a = rr_max(p, in->alleles);
        *b = rr_max(q, in->alleles);
    }
    else
    {
        float x = p[*a] * q[*b], y = p[*b] * q[*a];
        if (x < y || (x == y && jr_bits(rng, 1)))
        {
            int tmp = *a;
            *a = *b;
            *b = tmp;
        }
    }
}

typedef struct
{
    int n, rare, *samples;
} fp_carrier;
typedef struct
{
    const np_input *raw;
    const int *na;
    int *selected, *offset, *previous, *full_to_common;
    float *weight;
    fp_carrier *carriers;
    sq_input seq;
    ibs_result ibs;
    uint64_t phase_seed;
} fp_data;
static fp_data fp_prepare(const np_input *raw, const int *na,
                          const np_input *spliced)
{
    fp_data f = {0};
    f.raw = spliced;
    f.na = na;
    f.seq = fp_select(raw, na, &f.selected);
    it_input *in = &f.seq.input;
    int M = raw->M, N = raw->N, H = 2 * N, C = in->M;
    for (int m = 0; m < C; ++m)
        memcpy(in->gt + (size_t)m * H, spliced->gt + (size_t)f.selected[m] * H,
               (size_t)H * sizeof(int));
    f.phase_seed = in->seed;
    float *maf = ba_alloc(C, sizeof(float));
    for (int m = 0; m < C; ++m)
    {
        int count[2] = {0, 0};
        for (int h = 0; h < H; ++h)
        {
            int a = in->gt[(size_t)m * H + h];
            if (a >= 0)
                ++count[a];
        }
        int n = count[0] + count[1];
        maf[m] = n ? (float)((double)rp_min(count[0], count[1]) / n) : 0;
    }

    if (in->T == N)
        f.ibs = ibs_build(in->gt, in->pos, maf, C, N, NULL);
    else
    {
        int *target = ba_alloc((size_t)C * 2 * in->T, sizeof(int));
        for (int m = 0; m < C; ++m)
            memcpy(target + (size_t)m * 2 * in->T, in->gt + (size_t)m * H,
                   (size_t)2 * in->T * sizeof(int));
        f.ibs = ibs_build(target, in->pos, maf, C, in->T, NULL);
        free(target);
    }
    free(maf);
    f.offset = ba_alloc(M + 1, sizeof(int));
    f.previous = ba_alloc(M, sizeof(int));
    f.weight = ba_alloc(M, sizeof(float));
    f.full_to_common = ba_alloc(M, sizeof(int));
    for (int m = 0; m < M; ++m)
    {
        f.offset[m + 1] = f.offset[m] + na[m];
        f.full_to_common[m] = -1;
    }
    for (int m = 0; m < C; ++m)
        f.full_to_common[f.selected[m]] = m;
    f.carriers = ba_alloc(f.offset[M], sizeof(fp_carrier));
    int limit = (int)floor((double)((float)N * raw->rare));
    if (limit < 3)
        limit = 3;
    for (int m = 0; m < M; ++m)
    {
        fp_carrier *c = f.carriers + f.offset[m];
        for (int a = 0; a < na[m]; ++a)
            c[a].samples = ba_alloc(limit, sizeof(int));
        for (int s = 0; s < N; ++s)
        {
            int a = raw->gt[(size_t)m * H + 2 * s],
                b = raw->gt[(size_t)m * H + 2 * s + 1];
            if (a >= 0)
            {
                if (c[a].n < limit)
                    c[a].samples[c[a].n] = s;
                ++c[a].n;
            }
            if (b >= 0 && b != a)
            {
                if (c[b].n < limit)
                    c[b].samples[c[b].n] = s;
                ++c[b].n;
            }
        }
        for (int a = 0; a < na[m]; ++a)
            c[a].rare = (C < M && c[a].n <= limit);
    }
    for (int m = 0, left = 0; m < M; ++m)
    {
        while (left + 1 < C && f.selected[left + 1] <= m)
            ++left;
        f.previous[m] = left;
        f.weight[m] = (m <= f.selected[0] || left == C - 1)
                          ? 1.0f
                          : (float)((in->pos[left + 1] - raw->pos[m]) /
                                    (in->pos[left + 1] - in->pos[left]));
    }
    return f;
}

static void fp_codes(const it_input *in, rp_input *rp, int **starts, int **ends)
{
    int M = in->M, H = 2 * in->N, words = (M + 63) / 64;
    *starts = ba_alloc(M, sizeof(int));
    *ends = ba_alloc(M, sizeof(int));
    for (int end = 0; end < M;)
    {
        (*starts)[rp->S] = end;
        double boundary = in->pos[end] + in->step;
        ++end;
        while (end < M && in->pos[end] < boundary)
            ++end;
        (*ends)[rp->S++] = end;
    }
    rp->H = H;
    rp->T = 2 * in->T;
    rp->per_batch = (rp->S + in->threads - 1) / in->threads;
    rp->overlap = (int)rint((double)(in->buffer / in->step));
    rp->backoff = (int)rint((double)(0.3f / in->step));
    rp->candidates = rp_min(10, H);
    rp->seed = in->seed;
    rp->starts = *starts;
    rp->ends = *ends;
    rp->codes = ba_alloc((size_t)rp->S * H, sizeof(int));
    rp->sizes = ba_alloc(rp->S, sizeof(int));
    st_build_codes(in->haps, H, words, rp->S, *ends, rp->codes, rp->sizes,
                   st_phase_threads);
}
typedef struct
{
    rp_input input;
    int *donors;
} fp_direction_task;
static void fp_direction_run(void *data)
{
    fp_direction_task *t = data;
    rp_run(&t->input, t->donors, NULL);
}
typedef struct
{
    fp_data *f;
    int *out, *phase;
    rp_input *rp;
    int *starts, **donors;
    int begin, end;
} fp_sample_task;
static void fp_samples_run(void *data)
{
    fp_sample_task *t = data;
    fp_data *f = t->f;
    int *out = t->out, *phase = t->phase;
    rp_input rp = *t->rp;
    int *starts = t->starts, **donors = t->donors;
    it_input *in = &f->seq.input;
    const np_input *raw = f->raw;
    int M = raw->M, C = in->M, N = raw->N, H = 2 * N;
    int Kmax = in->max_states / 2,
        min_steps = rp_max(200, (int)ceil((double)(1.0f / in->step)));
    int *pair = ba_alloc((size_t)2 * rp.S, sizeof(int));
    ls_input li = {C, H,        rp.S,   Kmax, min_steps,
                   0, in->seed, starts, pair, phase};
    ls_work lw = ls_create(&li);
    int *ids[2] = {ba_alloc((size_t)C * Kmax, sizeof(int)),
                   ba_alloc((size_t)C * Kmax, sizeof(int))};
    float *prob[2] = {ba_alloc((size_t)C * Kmax, sizeof(float)),
                      ba_alloc((size_t)C * Kmax, sizeof(float))};
    unsigned char *mm = ba_alloc((size_t)C * Kmax, 1);
    float *back = ba_alloc(Kmax, sizeof(float));
    rr_donor *rd[2] = {ba_alloc(Kmax, sizeof(rr_donor)),
                       ba_alloc(Kmax, sizeof(rr_donor))};
    int max_alleles = 0;
    for (int m = 0; m < M; ++m)
        max_alleles = rp_max(max_alleles, f->na[m]);
    int *rare = ba_alloc(max_alleles, sizeof(int));
    float *p = ba_alloc(max_alleles, sizeof(float)),
          *q = ba_alloc(max_alleles, sizeof(float));
    for (int sample = t->begin; sample < t->end; ++sample)
    {
        int K[2];
        for (int bit = 0; bit < 2; ++bit)
        {
            li.target = 2 * sample + bit;
            for (int s = 0; s < rp.S; ++s)
            {
                pair[2 * s] = donors[0][(size_t)s * rp.T + li.target];
                pair[2 * s + 1] = donors[1][(size_t)s * rp.T + li.target];
            }
            K[bit] = ls_run(&li, &lw);
            memcpy(ids[bit], lw.ids, (size_t)C * Kmax * sizeof(int));
            for (int m = 0; m < C; ++m)
                memcpy(mm + (size_t)m * K[bit], lw.mismatch + (size_t)m * Kmax,
                       K[bit]);
            sp_input si = {
                C, K[bit], {1.0f - in->mismatch, in->mismatch}, in->recomb, mm};
            sp_forward(&si, prob[bit]);
            sp_backward(&si, prob[bit], back);
        }
        j_random rng;
        jr_seed(&rng, in->seed + (uint64_t)sample);
        for (int m = 0; m < M; ++m)
            if (f->full_to_common[m] < 0)
            {
                int left = f->previous[m], right = rp_min(left + 1, C - 1);
                for (int a = 0; a < f->na[m]; ++a)
                    rare[a] = f->carriers[f->offset[m] + a].rare;
                for (int bit = 0; bit < 2; ++bit)
                    for (int j = 0; j < K[bit]; ++j)
                    {
                        int h = ids[bit][(size_t)left * Kmax + j];
                        rd[bit][j] = (rr_donor){
                            raw->gt[(size_t)m * H + h],
                            raw->gt[(size_t)m * H + (h ^ 1)], h >= 2 * in->T,
                            prob[bit][(size_t)left * K[bit] + j],
                            prob[bit][(size_t)right * K[bit] + j]};
                    }
                rr_input ri = {f->na[m],
                               raw->gt[(size_t)m * H + 2 * sample],
                               raw->gt[(size_t)m * H + 2 * sample + 1],
                               {K[0], K[1]},
                               f->weight[m],
                               rare,
                               {rd[0], rd[1]}};
                int a, b;
                rr_run(&ri, &rng, p, q, &a, &b);
                out[(size_t)m * H + 2 * sample] = a;
                out[(size_t)m * H + 2 * sample + 1] = b;
            }
    }
    ls_destroy(&lw);
    free(pair);
    free(ids[0]);
    free(ids[1]);
    free(prob[0]);
    free(prob[1]);
    free(mm);
    free(back);
    free(rd[0]);
    free(rd[1]);
    free(rare);
    free(p);
    free(q);
}
static void fp_stage2(fp_data *f, int *out, FILE *trace)
{
    it_input *in = &f->seq.input;
    const np_input *raw = f->raw;
    int M = raw->M, C = in->M, N = raw->N, H = 2 * N, words = (C + 63) / 64;
    int *phase = ba_alloc((size_t)C * H, sizeof(int));
    for (int m = 0; m < M; ++m)
        memcpy(out + (size_t)m * H + 2 * in->T,
               raw->gt + (size_t)m * H + 2 * in->T,
               (size_t)2 * (N - in->T) * sizeof(int));
    for (int m = 0; m < C; ++m)
        for (int h = 0; h < H; ++h)
        {
            int a = (in->haps[(size_t)h * words + (m >> 6)] >> (m & 63)) & 1;
            phase[(size_t)m * H + h] = a;
            out[(size_t)f->selected[m] * H + h] = a;
        }
    if (C == M)
    {
        free(phase);
        return;
    }
    rp_input rp = {0};
    int *starts, *ends;
    fp_codes(in, &rp, &starts, &ends);
    rp.ibs = &f->ibs;
    rp.list_start = ba_alloc(rp.S + 1, sizeof(int));
    int capacity = f->offset[M];
    rp.lists = ba_alloc(capacity, sizeof(int *));
    rp.lengths = ba_alloc(capacity, sizeof(int));
    int lists = 0;
    for (int s = 0; s < rp.S; ++s)
    {
        rp.list_start[s] = lists;
        int begin = s ? f->selected[starts[s]] : 0,
            end = s + 1 < rp.S ? f->selected[starts[s + 1]] : M;
        for (int m = begin; m < end; ++m)
            for (int a = 0; a < f->na[m]; ++a)
            {
                fp_carrier c = f->carriers[f->offset[m] + a];
                if (c.rare && c.n > 1)
                {
                    rp.lengths[lists] = 2 * c.n;
                    int *v = ba_alloc(2 * c.n, sizeof(int));
                    for (int j = 0; j < c.n; ++j)
                    {
                        v[2 * j] = 2 * c.samples[j];
                        v[2 * j + 1] = 2 * c.samples[j] + 1;
                    }
                    rp.lists[lists++] = v;
                }
            }
    }
    rp.list_start[rp.S] = lists;
    int *donors[2] = {ba_alloc((size_t)rp.S * rp.T, sizeof(int)),
                      ba_alloc((size_t)rp.S * rp.T, sizeof(int))};
#ifdef MT
    int directions = trace ? 1 : (st_phase_threads < 2 ? 1 : 2);
#pragma omp parallel for num_threads(directions) if (directions > 1)
#endif
    for (int b = 0; b < 2; ++b)
    {
        fp_direction_task t = {rp, donors[b]};
        t.input.backward = b;
        st_worker_run(fp_direction_run, &t, 5);
    }
    st_workers_check();
    if (trace)
        for (int b = 0; b < 2; ++b)
            for (int s = 0; s < rp.S; ++s)
                it_ints(trace, b ? "BWD" : "FWD", s,
                        donors[b] + (size_t)s * rp.T, rp.T);
    int workers = trace ? 1 : st_phase_threads;
    if (workers > in->T)
        workers = in->T;
#ifdef MT
#pragma omp parallel num_threads(workers) if (workers > 1)
#endif
    {
        int tid = omp_get_thread_num(), team = omp_get_num_threads();
        fp_sample_task t = {f,
                            out,
                            phase,
                            &rp,
                            starts,
                            donors,
                            in->T * tid / team,
                            in->T * (tid + 1) / team};
        st_worker_run(fp_samples_run, &t, 6);
    }
    st_workers_check();
    free(phase);
    free(starts);
    free(ends);
    free(rp.codes);
    free(rp.sizes);
    free(rp.list_start);
    for (int j = 0; j < lists; ++j)
        free(rp.lists[j]);
    free(rp.lists);
    free(rp.lengths);
    free(donors[0]);
    free(donors[1]);
}
static void fp_free(fp_data *f)
{
    for (int a = 0; a < f->offset[f->raw->M]; ++a)
        free(f->carriers[a].samples);
    free(f->carriers);
    free(f->offset);
    free(f->previous);
    free(f->weight);
    free(f->full_to_common);
    free(f->selected);
    ibs_free(&f->ibs);
    np_free(&f->seq);
}

static int *np_phase_overlap(int M, int N, const int *g, uint64_t seed,
                             int overlap, int T)
{
    int H = 2 * N;
    size_t Z = (size_t)M * H;
    int *rev = calloc(Z, sizeof(int)), *fwd = calloc(Z, sizeof(int));
    int *a = calloc(H, sizeof(int)), *order = calloc(H, sizeof(int)),
        *next = calloc(H, sizeof(int));
    int *rank = calloc(H, sizeof(int)), *unph = calloc(N, sizeof(int)),
        *missing = calloc(N, sizeof(int));
    int *last = calloc(N, sizeof(int));
    if (!rev || !fwd || !a || !order || !next || !rank || !unph || !missing ||
        !last)
        st_failure(ST_PHASE_MEMORY);
    rng = ((uint64_t)seed ^ UINT64_C(0x5deece66d)) & ((UINT64_C(1) << 48) - 1);
    for (int pass = 0; pass < 2; ++pass)
    {
        for (int h = 0; h < H; ++h)
            order[h] = h;
        for (int s = 0; s < N; ++s)
            last[s] = -1;
        for (int t = 0; t < M; ++t)
        {
            int m = pass ? t : M - 1 - t;
            if (t)
            {
                int zeros = 0;
                for (int h = 0; h < H; ++h)
                    zeros += a[h] == 0;
                int offset[2] = {0, zeros};
                for (int j = 0; j < H; ++j)
                {
                    int h = order[j];
                    next[offset[a[h]]++] = h;
                }
                memcpy(order, next, (size_t)H * sizeof(int));
            }
            for (int h = 0; h < H; ++h)
                rank[order[h]] = h;
            memcpy(a, g + (size_t)m * H, (size_t)H * sizeof(int));
            int cdf[2] = {0, 0};
            for (int s = 0; s < N; ++s)
            {
                int h = 2 * s;
                missing[s] = a[h] < 0 || a[h + 1] < 0;
                unph[s] = s < T && !missing[s] && a[h] != a[h + 1];
                if (a[h] >= 0)
                    ++cdf[a[h]];
                if (a[h + 1] >= 0)
                    ++cdf[a[h + 1]];
            }
            cdf[1] += cdf[0];
            if (m >= overlap)
                partial(a, unph, order, rank, N);
            if (m >= overlap)
                for (int s = 0; s < T; ++s)
                {
                    int h = 2 * s, prev = last[s];
                    if (unph[s])
                    {
                        int swap = 0;
                        if (!pass)
                            swap = (int)next_bits(1);
                        else if (prev >= 0)
                        {
                            size_t q = (size_t)prev * H + h,
                                   z = (size_t)m * H + h;
                            int sameRev =
                                (rev[q] < rev[q + 1]) == (rev[z] < rev[z + 1]);
                            int sameFwd =
                                (fwd[q] < fwd[q + 1]) == (a[h] < a[h + 1]);
                            swap = sameRev != sameFwd;
                        }
                        if (swap)
                        {
                            int x = a[h];
                            a[h] = a[h + 1];
                            a[h + 1] = x;
                        }
                        unph[s] = 0;
                    }
                    else
                    {
                        for (int side = 0; side < 2; ++side)
                            if (a[h + side] < 0)
                            {
                                if (!pass)
                                    a[h + side] =
                                        cdf[1] ? (next_int((uint32_t)cdf[1]) >=
                                                  (uint32_t)cdf[0])
                                               : 0;
                                else
                                {
                                    int aligned = side;
                                    if (prev >= 0)
                                    {
                                        size_t q = (size_t)prev * H + h;
                                        if ((rev[q] < rev[q + 1]) !=
                                            (fwd[q] < fwd[q + 1]))
                                            aligned ^= 1;
                                    }
                                    a[h + side] =
                                        rev[(size_t)m * H + h + aligned];
                                }
                            }
                    }
                }
            memcpy((pass ? fwd : rev) + (size_t)m * H, a,
                   (size_t)H * sizeof(int));
            if (pass)
            {
                for (int s = 0; s < N; ++s)
                    if (!missing[s] && a[2 * s] != a[2 * s + 1])
                        last[s] = m;
            }
        }
    }
    free(rev);
    free(a);
    free(order);
    free(next);
    free(rank);
    free(unph);
    free(missing);
    free(last);
    return fwd;
}

static int gm_search_double(const double *a, int n, double key)
{
    int lo = 0, hi = n - 1;
    while (lo <= hi)
    {
        int mid = (lo + hi) >> 1;
        if (a[mid] < key)
            lo = mid + 1;
        else if (a[mid] > key)
            hi = mid - 1;
        else
        {
            int64_t abits, kbits;
            memcpy(&abits, a + mid, sizeof(abits));
            memcpy(&kbits, &key, sizeof(kbits));
            if (abits == kbits)
                return mid;
            if (abits < kbits)
                lo = mid + 1;
            else
                hi = mid - 1;
        }
    }
    return -lo - 1;
}
static double gm_position(const int *bp, const double *cm, int n, int x)
{
    int lo = 0, hi = n - 1;
    while (lo <= hi)
    {
        int mid = (lo + hi) >> 1;
        if (bp[mid] < x)
            lo = mid + 1;
        else if (bp[mid] > x)
            hi = mid - 1;
        else
            return cm[mid];
    }
    int a = lo - 1, b = lo;
    if (a == n - 1)
    {
        int i = gm_search_double(cm, n, cm[n - 1] - 5.0);
        if (i < 0)
            i = -i - 2;
        a = i > 0 ? i : 0;
        b = n - 1;
    }
    else if (b == 0)
    {
        int i = gm_search_double(cm, n, cm[0] + 5.0);
        if (i < 0)
            i = -i - 1;
        a = 0;
        b = i < n - 1 ? i : n - 1;
    }
    return cm[a] +
           ((double)(x - bp[a]) / (double)(bp[b] - bp[a])) * (cm[b] - cm[a]);
}

static int gm_prepare(const int *bp, const double *cm, int n,
                      const int *positions, int m, double *out)
{
    if (m < 2 || positions[0] == positions[m - 1] || n == 1 || n < 0)
        return 0;
    for (int i = 0; i < m; ++i)
        if (positions[i] < 0 || (i && positions[i] < positions[i - 1]))
            return 0;
    for (int i = 0; i < n; ++i)
        if (bp[i] < 0 || !isfinite(cm[i]) ||
            (i && (bp[i] <= bp[i - 1] || cm[i] < cm[i - 1])))
            return 0;
    if (n && cm[0] == cm[n - 1])
        return 0;
    double first =
        n ? gm_position(bp, cm, n, positions[0]) : positions[0] * 1e-6;
    double last =
        n ? gm_position(bp, cm, n, positions[m - 1]) : positions[m - 1] * 1e-6;
    double minimum =
        fmax(fabs(last - first) / (positions[m - 1] - positions[0]), 1e-8);
    double previous = first;
    out[0] = first;
    for (int i = 1; i < m; ++i)
    {
        double current =
            n ? gm_position(bp, cm, n, positions[i]) : positions[i] * 1e-6;
        out[i] = out[i - 1] + fmax(current - previous, minimum);
        previous = current;
    }
    return 1;
}

static int st_pbwt_to(const double *pos, int m, double key)
{
    int i = gm_search_double(pos, m, key);
    return i < 0 ? -i - 1 : i + 1;
}
static int *np_phase_windows(int M, int N, const int *gt, const double *pos,
                             uint64_t seed, int overlap, int T)
{
    double advance = fmax(2.0, pos[M - 1] - pos[0]);
    int end = st_pbwt_to(pos, M, pos[0] + advance);
    if (end == M)
        return np_phase_overlap(M, N, gt, seed, overlap, T);
    int haps = 2 * N, start = 0, previous_end = 0, part = 0;
    int *out = ba_alloc((size_t)M * haps, sizeof(int));
    int *hets = ba_alloc(M, sizeof(int));
    for (;;)
    {
        int fixed = overlap > start ? overlap - start : 0;
        int *phase = np_phase_overlap(end - start, N, gt + (size_t)start * haps,
                                      seed + part, fixed, T);
        int copy_start = (start + previous_end) >> 1;
        for (int s = 0; s < N; ++s)
        {
            int swap = 0;
            if (start > 0 && s < T)
            {
                int n = 0, seen = 0;
                for (int m = 0; m < M; ++m)
                {
                    int a = gt[(size_t)m * haps + 2 * s],
                        b = gt[(size_t)m * haps + 2 * s + 1];
                    if (a >= 0 && b >= 0 && a != b)
                    {
                        if (seen && m >= overlap)
                            hets[n++] = m;
                        seen = 1;
                    }
                }
                if (n)
                {
                    int i = 0;
                    while (i < n && hets[i] < copy_start)
                        ++i;
                    if (i == n || (hets[i] >= previous_end && i > 0))
                        --i;
                    int het = hets[i];
                    if (start <= het && het < previous_end)
                    {
                        size_t a = (size_t)het * haps + 2 * s;
                        size_t b = (size_t)(het - start) * haps + 2 * s;
                        swap = out[a] == phase[b + 1] && out[a + 1] == phase[b];
                    }
                }
            }
            for (int m = copy_start; m < end; ++m)
            {
                size_t a = (size_t)m * haps + 2 * s;
                size_t b = (size_t)(m - start) * haps + 2 * s;
                out[a] = phase[b + swap];
                out[a + 1] = phase[b + (swap ^ 1)];
            }
        }
        free(phase);
        if (end == M)
            break;
        int next_start = gm_search_double(pos, M, pos[end] - 0.5);
        start = next_start < 0 ? -next_start - 1 : next_start;
        previous_end = end;
        end = st_pbwt_to(pos, M, pos[end] + advance);
        ++part;
    }
    free(hets);
    return out;
}

static void np_initialize_overlap(sq_input *seq, float ne, FILE *trace,
                                  int overlap)
{
    it_input *in = &seq->input;
    int M = in->M, N = in->N, H = 2 * N, words = (M + 63) / 64;
    int *phase =
        np_phase_windows(M, N, in->gt, in->pos, in->seed, overlap, in->T);
    memset(in->haps, 0, (size_t)H * words * sizeof(uint64_t));
    int *a = ba_alloc(M, sizeof(int)), *b = ba_alloc(M, sizeof(int)),
        *unph = ba_alloc(M, sizeof(int)), *miss = ba_alloc(M, sizeof(int));
    for (int s = 0; s < N; ++s)
    {
        int U = 0, D = 0, seen = 0;
        for (int m = 0; m < M; ++m)
        {
            a[m] = phase[(size_t)m * H + 2 * s];
            b[m] = phase[(size_t)m * H + 2 * s + 1];
            int ga = in->gt[(size_t)m * H + 2 * s],
                gb = in->gt[(size_t)m * H + 2 * s + 1];
            if (ga < 0 || gb < 0)
                miss[D++] = m;
            else if (s < in->T && ga != gb)
            {
                if (seen && m >= overlap)
                    unph[U++] = m;
                seen = 1;
            }
            in->haps[(size_t)(2 * s) * words + (m >> 6)] |= (uint64_t)a[m]
                                                            << (m & 63);
            in->haps[(size_t)(2 * s + 1) * words + (m >> 6)] |= (uint64_t)b[m]
                                                                << (m & 63);
        }
        hi_free(in->clusters + s);
        in->clusters[s] = hi_cluster(M, a, b, in->pos, U, unph, D, miss);
        if (trace)
        {
            it_ints(trace, "INITENDS", s, in->clusters[s].end,
                    in->clusters[s].n);
            it_ints(trace, "INITTYPES", s, in->clusters[s].type,
                    in->clusters[s].n);
            it_ints(trace, "INITA", s, a, M);
            it_ints(trace, "INITB", s, b, M);
        }
    }
    double theta = 1.0 / (log(H) + 0.5);
    in->mismatch = (float)(theta / (2 * (theta + H)));
    seq->intensity = (0.04f * ne) / H;
    in->it = 0;
    in->threshold = sq_threshold(in, seq->initial_lr);
    for (int m = 0; m < M; ++m)
        in->recomb[m] =
            (float)-expm1(-(double)seq->intensity * seq->distance[m]);
    if (trace)
    {
        fprintf(trace, "INITPARAM %u %u\n", pe_float_bits(in->mismatch),
                pe_float_bits(seq->intensity));
    }
    free(phase);
    free(a);
    free(b);
    free(unph);
    free(miss);
}

typedef struct
{
    double window_cm;
    double overlap_cm;
    int max_markers;
} st_window_options;

typedef struct
{
    int start, end;
    int emit_start, emit_end, overlap_start, fixed_prefix;
} st_window;

static int st_map_base(const int *bp, const double *cm, int n, double genetic,
                       int *out)
{
    double value;
    if (!n)
        value = genetic * 1e6;
    else
    {
        int index = gm_search_double(cm, n, genetic);
        if (index >= 0)
        {
            *out = bp[index];
            return 1;
        }
        int a = -index - 2, b = -index - 1;
        if (a == n - 1)
        {
            int i = gm_search_double(cm, n, cm[n - 1] - 5.0);
            if (i < 0)
                i = -i - 2;
            a = i > 0 ? i : 0;
            b = n - 1;
        }
        else if (b == 0)
        {
            int i = gm_search_double(cm, n, cm[0] + 5.0);
            if (i < 0)
                i = -i - 1;
            a = 0;
            b = i < n - 1 ? i : n - 1;
        }
        if (cm[b] <= cm[a])
            return 0;
        value = bp[a] + ((genetic - cm[a]) / (cm[b] - cm[a])) * (bp[b] - bp[a]);
    }
    if (!isfinite(value) || value >= INT_MAX || value < INT_MIN)
        return 0;
    *out = (int)floor(value + 0.5);
    return 1;
}

static int st_plan_windows(const int *map_bp, const double *map_cm, int map_n,
                           const int *positions, int m,
                           st_window_options options, st_window **plan,
                           int *count)
{
    *plan = NULL;
    *count = 0;
    if (!positions || m < 2 || map_n < 0 || map_n == 1 ||
        (map_n && (!map_bp || !map_cm)))
        return 0;
    int whole = isinf(options.window_cm) && options.window_cm > 0;
    if ((!whole && (!isfinite(options.window_cm) || options.window_cm <= 0)) ||
        !isfinite(options.overlap_cm) || options.overlap_cm < 0 ||
        options.max_markers < 4)
        return 0;
    float window = (float)options.window_cm,
          overlap = (float)options.overlap_cm;
    if (!whole && (!isfinite(window) || !isfinite(overlap) || window <= 0 ||
                   window < 1.1 * overlap))
        return 0;
    for (int i = 0; i < m; ++i)
        if (positions[i] < 0 || (i && positions[i] < positions[i - 1]))
            return 0;
    if (positions[0] == positions[m - 1])
        return 0;
    for (int i = 0; i < map_n; ++i)
        if (map_bp[i] < 0 || !isfinite(map_cm[i]) ||
            (i && (map_bp[i] <= map_bp[i - 1] || map_cm[i] < map_cm[i - 1])))
            return 0;
    if (map_n && map_cm[0] == map_cm[map_n - 1])
        return 0;
    st_window *windows = calloc((size_t)m, sizeof(*windows));
    if (!windows)
        return 0;
    int start = 0, next = 0, used = 0;
    while (next < m)
    {
        int carried = next - start, end = next;
        if (whole)
            end = m;
        else
        {
            double cm =
                map_n ? gm_position(map_bp, map_cm, map_n, positions[next])
                      : positions[next] * 1e-6;
            cm += carried ? (window - overlap) : window;
            int end_bp;
            if (!st_map_base(map_bp, map_cm, map_n, cm, &end_bp))
                goto invalid;
            while (end < m && positions[end] < end_bp &&
                   end - start < options.max_markers)
                ++end;
        }
        int n = end - start;
        if (end <= next || n < 2 || positions[start] == positions[end - 1])
            goto invalid;
        int keep = n;
        if (end < m)
        {
            double cm =
                map_n ? gm_position(map_bp, map_cm, map_n, positions[end - 1])
                      : positions[end - 1] * 1e-6;
            int key;
            if (!st_map_base(map_bp, map_cm, map_n, cm - overlap, &key))
                goto invalid;
            int lo = n - (options.max_markers >> 2), hi = n - 1;
            if (lo < 0)
                lo = 0;
            while (lo <= hi)
            {
                int mid = (lo + hi) >> 1;
                if (positions[start + mid] < key)
                    lo = mid + 1;
                else if (positions[start + mid] > key)
                    hi = mid - 1;
                else
                {
                    hi = mid;
                    break;
                }
            }
            keep = hi > 0 ? hi : 0;
            while (keep > 0 &&
                   positions[start + keep - 1] == positions[start + keep])
                --keep;
        }
        int emit_start = carried >> 1, emit_end = (n + keep) >> 1;
        if (emit_end <= emit_start ||
            (used && windows[used - 1].start + windows[used - 1].emit_end !=
                         start + emit_start))
            goto invalid;
        windows[used++] =
            (st_window){start, end, emit_start, emit_end, keep, emit_start};
        start += keep;
        next = end;
    }
    *plan = windows;
    *count = used;
    return 1;
invalid:
    free(windows);
    return 0;
}

static int st_phase_windowed(const np_input *chromosome, const int *alleles,
                             const int *map_bp, const double *map_cm, int map_n,
                             st_window_options options, double initial_error,
                             j_random *master, int chrom, int *window_number,
                             int *output, FILE *trace)
{
    if (!chromosome || !alleles || !master || !window_number || !output ||
        chrom < 1)
        return 0;
    if (!isnan(initial_error) &&
        (!isfinite(initial_error) || initial_error < 0 || initial_error > 1))
        return 0;
    const np_input *c = chromosome;
    if (c->M < 2 || c->N < 2 || c->N > 5000 || c->T < 1 || c->T > c->N ||
        c->K < 2 || c->burn < 0 || c->iterations < 1 || !c->gt || !c->bp ||
        !isfinite(c->rare) || c->rare < 0 || c->rare > 1 ||
        !isfinite(c->scale) || c->scale <= 0 || !isfinite(c->ne) ||
        c->ne <= 0 || !isfinite(c->buffer) || c->buffer <= 0 ||
        !isfinite(c->lr) || c->lr <= 0)
        return 0;
    int haps = 2 * c->N;
    for (int m = 0; m < c->M; ++m)
    {
        if (alleles[m] < 1 || alleles[m] > 100)
            return 0;
        for (int h = 0; h < haps; ++h)
            if (c->gt[(size_t)m * haps + h] < -1 ||
                c->gt[(size_t)m * haps + h] >= alleles[m] ||
                (h >= 2 * c->T && c->gt[(size_t)m * haps + h] < 0))
                return 0;
    }

    for (int s = 0; s < c->N; ++s)
    {
        int observed = 0;
        for (int m = 0; m < c->M && !observed; ++m)
            observed = c->gt[(size_t)m * haps + 2 * s] >= 0 ||
                       c->gt[(size_t)m * haps + 2 * s + 1] >= 0;
        if (!observed)
            return 0;
    }
    st_window *plan = NULL;
    int count = 0;
    if (!st_plan_windows(map_bp, map_cm, map_n, c->bp, c->M, options, &plan,
                         &count))
        return 0;
    int *previous = NULL, previous_start = 0, previous_end = 0, ok = 1;
    for (int j = 0; j < count; ++j)
    {
#ifdef ST_PHASE_CHECK_INTERRUPT
        ST_PHASE_CHECK_INTERRUPT();
#endif
        st_window w = plan[j];
        np_input raw = *c;
        raw.M = w.end - w.start;
        raw.bp = c->bp + w.start;
        raw.gt = c->gt + (size_t)w.start * haps;
        raw.pos = ba_alloc(raw.M, sizeof(double));
        if (!gm_prepare(map_bp, map_cm, map_n, raw.bp, raw.M, raw.pos))
        {
            free(raw.pos);
            ok = 0;
            break;
        }
        np_input spliced = raw;
        spliced.gt = ba_alloc((size_t)raw.M * haps, sizeof(int));
        memcpy(spliced.gt, raw.gt, (size_t)raw.M * haps * sizeof(int));
        if (w.fixed_prefix)
        {
            if (!previous || previous_end - previous_start != w.fixed_prefix)
            {
                free(raw.pos);
                free(spliced.gt);
                ok = 0;
                break;
            }
            memcpy(spliced.gt, previous + (size_t)previous_start * haps,
                   (size_t)w.fixed_prefix * haps * sizeof(int));
        }
        fp_data f = fp_prepare(&raw, alleles + w.start, &spliced);
        if (!isfinite(f.seq.input.step) || f.seq.input.step <= 0 ||
            (double)(1.0f / f.seq.input.step) > INT_MAX - 1024.0)
        {
            fp_free(&f);
            free(raw.pos);
            free(spliced.gt);
            ok = 0;
            break;
        }
        int common_fixed = 0;
        while (common_fixed < f.seq.input.M &&
               f.selected[common_fixed] < w.fixed_prefix)
            ++common_fixed;

        uint64_t high = (uint64_t)(int64_t)(int32_t)jr_bits(master, 32);
        int64_t low = (int32_t)jr_bits(master, 32);
        f.phase_seed = (high << 32) + (uint64_t)low;
        it_input *in = &f.seq.input;
        in->seed = f.phase_seed;
        np_initialize_overlap(&f.seq, raw.ne, NULL, common_fixed);
        sq_check_reference(in);
        if (!isnan(initial_error))
            in->mismatch = (float)initial_error;
        sq_result result = sq_run(&f.seq, &f.ibs, NULL);
        in->seed = f.phase_seed + (uint64_t)in->it;
        int *out = ba_alloc((size_t)raw.M * haps, sizeof(int));
        fp_stage2(&f, out, NULL);
        for (int m = 0; m < raw.M; ++m)
            for (int h = 2 * c->T; h < haps; ++h)
                if (out[(size_t)m * haps + h] != raw.gt[(size_t)m * haps + h])
                    st_failure(ST_PHASE_REFERENCE);
        memcpy(output + (size_t)(w.start + w.emit_start) * haps,
               out + (size_t)w.emit_start * haps,
               (size_t)(w.emit_end - w.emit_start) * haps * sizeof(int));
        ++*window_number;
        if (trace)
        {
            fprintf(trace, "WINDOW %d %d %d %d %d %d\n", *window_number, chrom,
                    w.emit_start, w.emit_end, w.overlap_start, w.fixed_prefix);
            it_ints(trace, "SELECTED", -1, f.selected, f.seq.input.M);
            fprintf(trace, "WINDOWSEED %" PRIu64 "\n", f.phase_seed);
            fprintf(trace, "FINAL %d %" PRIu64 " %u %u\n", result.iterations,
                    in->seed, pe_float_bits(in->mismatch),
                    pe_float_bits(f.seq.intensity));
            for (int m = 0; m < raw.M; ++m)
                it_ints(trace, "CALL", m, out + (size_t)m * haps, haps);
        }
        free(previous);
        previous = out;
        previous_start = w.overlap_start;
        previous_end = w.emit_end;
        fp_free(&f);
        free(raw.pos);
        free(spliced.gt);
    }
    free(previous);
    free(plan);
    return ok;
}

#undef malloc
#undef calloc
#undef realloc
#undef free
#undef ST_PHASE_CHECK_INTERRUPT

static void st_cleanup(void *unused)
{
    (void)unused;
    while (st_blocks)
        st_release(st_blocks + 1);
    er_set_threads(1);
    st_phase_threads = 1;
    atomic_store(&st_worker_error, 0);
    atomic_store(&st_max_team, 1);
    st_active_guard = NULL;
    st_phase_dataset = NULL;
    st_busy = 0;
}

extern gs_varset_type *gs_fdta(char *set_name);

enum
{
    STP_WINDOW,
    STP_OVERLAP,
    STP_MARKERS,
    STP_BURNIN,
    STP_ITERATIONS,
    STP_STATES,
    STP_NE,
    STP_ERR,
    STP_EM,
    STP_LR,
    STP_SCALE,
    STP_RARE,
    STP_BUFFER,
    STP_SEED,
    STP_OPTIONS
};

static const char *stp_option_names[] = {
    "window",     "overlap",    "window.markers",
    "burnin",     "iterations", "phase.states",
    "ne",         "err",        "em",
    "initial.lr", "step.scale", "rare",
    "buffer",     "seed"};

static SEXP stp_list(int n, const char *const *labels)
{
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; i++)
        SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
    Rf_setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(2);
    return result;
}

static int stp_field_index(SEXP x, const char *name)
{
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(x) != VECSXP || TYPEOF(names) != STRSXP ||
        XLENGTH(names) != XLENGTH(x))
        st_phase_error("Invalid STvcf fields");
    for (int i = 0; i < Rf_length(x); i++)
        if (STRING_ELT(names, i) != NA_STRING &&
            !strcmp(CHAR(STRING_ELT(names, i)), name))
            return i;
    st_phase_error("Missing STvcf field: %s", name);
    return -1;
}

static SEXP stp_field(SEXP x, const char *name)
{
    return VECTOR_ELT(x, stp_field_index(x, name));
}

static double stp_number_at(SEXP x, R_xlen_t i)
{
    if (TYPEOF(x) == REALSXP)
        return REAL(x)[i];
    return INTEGER(x)[i] == NA_INTEGER ? NA_REAL : INTEGER(x)[i];
}

static int stp_logical(SEXP x, const char *name)
{
    if (TYPEOF(x) != LGLSXP || XLENGTH(x) != 1 || LOGICAL(x)[0] == NA_LOGICAL)
        st_phase_error("%s must be TRUE or FALSE", name);
    return LOGICAL(x)[0];
}

static int stp_automatic(SEXP x)
{
    return TYPEOF(x) == STRSXP && XLENGTH(x) == 1 &&
           STRING_ELT(x, 0) != NA_STRING &&
           !strcmp(CHAR(STRING_ELT(x, 0)), "auto");
}

typedef struct
{
    double value[STP_OPTIONS];
    double *window, *overlap;
    R_xlen_t nwindow, noverlap;
    int automatic, error_set, seed_set;
} stp_controls;

typedef struct
{
    const char *name;
    int index;
} stp_name_entry;

typedef struct
{
    int markers, samples, nlevels;
    int *chrom;
    double *pos;
    const char **id, **levels;
} stp_metadata;

typedef struct
{
    int count;
    const char **label;
} stp_alleles;

static int stp_name_compare(const void *left, const void *right)
{
    const stp_name_entry *x = left, *y = right;
    return strcmp(x->name, y->name);
}

static int stp_names_unique(const char *const *names, size_t n)
{
    stp_name_entry *ordered = st_zero(n, sizeof(*ordered));
    int valid = 1;
    for (size_t i = 0; i < n; i++)
    {
        if (!names[i] || !names[i][0])
        {
            valid = 0;
            break;
        }
        ordered[i].name = names[i];
    }
    if (valid)
    {
        qsort(ordered, n, sizeof(*ordered), stp_name_compare);
        for (size_t i = 1; i < n; i++)
            if (!strcmp(ordered[i - 1].name, ordered[i].name))
            {
                valid = 0;
                break;
            }
    }
    st_release(ordered);
    return valid;
}

static char *stp_input_text(SEXP x)
{
    const char *text =
        Rf_getCharCE(x) == CE_BYTES ? CHAR(x) : Rf_translateCharUTF8(x);
    size_t length = strlen(text);
    char *out = st_zero(length + 1, 1);
    memcpy(out, text, length);
    return out;
}

static int stp_input_names_unique(SEXP x)
{
    if (TYPEOF(x) != STRSXP)
        return 0;
    size_t n = (size_t)XLENGTH(x);
    const char **names = st_zero(n, sizeof(*names));
    for (size_t i = 0; i < n; i++)
        if (STRING_ELT(x, (R_xlen_t)i) != NA_STRING)
            names[i] = stp_input_text(STRING_ELT(x, (R_xlen_t)i));
    int valid = stp_names_unique(names, n);
    for (size_t i = 0; i < n; i++)
        st_release((void *)names[i]);
    st_release(names);
    return valid;
}

static void stp_names(SEXP x, const char *role, const char *field)
{
    if (!stp_input_names_unique(x))
        st_phase_error("Invalid or duplicated %s %s", role, field);
}

static void stp_number(SEXP x, int index, double lower, double upper, int whole,
                       int vector, int infinity)
{
    if ((TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP) || XLENGTH(x) < 1 ||
        (!vector && XLENGTH(x) != 1))
        st_phase_error("Invalid option: %s", stp_option_names[index]);
    for (R_xlen_t i = 0; i < XLENGTH(x); i++)
    {
        double value = stp_number_at(x, i);
        if (ISNAN(value) || (!infinity && !R_FINITE(value)) || value < lower ||
            value > upper || (whole && value != trunc(value)))
            st_phase_error("Invalid option: %s", stp_option_names[index]);
    }
}

static double *stp_number_copy(SEXP x)
{
    R_xlen_t n = XLENGTH(x);
    double *out = st_zero((size_t)n, sizeof(*out));
    for (R_xlen_t i = 0; i < n; i++)
        out[i] = stp_number_at(x, i);
    return out;
}

static void stp_options(SEXP input, stp_controls *out, SEXP *seed_input)
{
    SEXP names = Rf_getAttrib(input, R_NamesSymbol);
    if (TYPEOF(input) != VECSXP ||
        (XLENGTH(input) &&
         (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(input) ||
          !stp_input_names_unique(names))))
        st_phase_error("phasing.options must be a uniquely named list");

    const double defaults[] = {0, 2, 4000000, 3, 12,   280, 100000,
                               0, 1, 100000,  3, .002, 1,   0};
    memset(out, 0, sizeof(*out));
    memcpy(out->value, defaults, sizeof(defaults));
    *seed_input = R_NilValue;
    SEXP values[STP_OPTIONS];
    int supplied[STP_OPTIONS] = {0};
    for (int i = 0; i < STP_OPTIONS; i++)
        values[i] = R_NilValue;
    for (R_xlen_t i = 0; i < XLENGTH(input); i++)
    {
        const char *name = CHAR(STRING_ELT(names, i));
        size_t size = strlen(name);
        char *normalized = st_zero(size + 1, 1);
        for (size_t j = 0; j <= size; j++)
            normalized[j] = name[j] == '-' ? '.' : name[j];
        if (!strcmp(normalized, "nthreads"))
            st_phase_error(
                "Set threads with st.set.num.threads(); nthreads is not a "
                "phasing.options entry");
        if (!strcmp(normalized, "ref"))
            st_phase_error("Use the top-level ref argument");
        if (!strcmp(normalized, "impute"))
            st_phase_error(
                "Use the top-level impute argument for missing calls at typed "
                "markers");
        int index = -1;
        for (int j = 0; j < STP_OPTIONS; j++)
            if (!strcmp(normalized, stp_option_names[j]))
                index = j;
        if (index < 0)
            st_phase_error(
                "Unsupported phasing option: %s. See help('st.phase').",
                normalized);
        if (supplied[index])
            st_phase_error("Duplicate phasing option after name normalization");
        supplied[index] = 1;
        values[index] = VECTOR_ELT(input, i);
        st_release(normalized);
    }

    out->automatic = !supplied[STP_WINDOW] || stp_automatic(values[STP_WINDOW]);
    if (!out->automatic)
    {
        stp_number(values[STP_WINDOW], STP_WINDOW, DBL_MIN, R_PosInf, 0, 1, 1);
        out->nwindow = XLENGTH(values[STP_WINDOW]);
        out->window = stp_number_copy(values[STP_WINDOW]);
    }
    if (!out->automatic ||
        (supplied[STP_OVERLAP] && values[STP_OVERLAP] != R_NilValue))
    {
        if (supplied[STP_OVERLAP])
        {
            stp_number(values[STP_OVERLAP], STP_OVERLAP, 0, R_PosInf, 0, 1, 0);
            out->noverlap = XLENGTH(values[STP_OVERLAP]);
            out->overlap = stp_number_copy(values[STP_OVERLAP]);
        }
        else
        {
            out->noverlap = 1;
            out->overlap = st_zero(1, sizeof(*out->overlap));
            out->overlap[0] = defaults[STP_OVERLAP];
        }
    }
    if (!out->automatic)
    {
        R_xlen_t nw = out->nwindow, no = out->noverlap;
        if (nw > 1 && no > 1 && nw != no)
            st_phase_error("Window/overlap lengths do not match");
        for (R_xlen_t i = 0; i < (nw > no ? nw : no); i++)
        {
            double w = out->window[i % nw], o = out->overlap[i % no];
            if (R_FINITE(w) && w < 1.1 * o)
                st_phase_error("window must be at least 1.1 times overlap");
        }
    }
    const int indices[] = {STP_BURNIN, STP_ITERATIONS, STP_STATES, STP_MARKERS,
                           STP_NE,     STP_SCALE,      STP_BUFFER, STP_LR,
                           STP_RARE,   STP_ERR,        STP_SEED};
    const double lower[] = {1,       1, 2,       4, DBL_MIN, DBL_MIN,
                            DBL_MIN, 1, DBL_MIN, 0, -INT_MAX};
    const double upper[] = {INT_MAX, INT_MAX, 1000000, INT_MAX, 3.4e38, 3.4e38,
                            3.4e38,  3.4e38,  .5,      1,       INT_MAX};
    for (size_t i = 0; i < sizeof(indices) / sizeof(indices[0]); i++)
    {
        int k = indices[i];
        if (!supplied[k] ||
            ((k == STP_ERR || k == STP_SEED) && values[k] == R_NilValue))
            continue;
        stp_number(values[k], k, lower[i], upper[i], i < 4 || k == STP_SEED, 0,
                   0);
        out->value[k] = stp_number_at(values[k], 0);
    }
    if (supplied[STP_EM])
        out->value[STP_EM] = stp_logical(values[STP_EM], "em");
    out->error_set = values[STP_ERR] != R_NilValue;
    out->seed_set = values[STP_SEED] != R_NilValue;
    *seed_input = values[STP_SEED];
}

static stp_alleles stp_allele_labels(SEXP ref, SEXP alt, const char *role)
{
    const char *r = stp_input_text(ref), *a = stp_input_text(alt);
    if (!r[0] || !a[0] || !strcmp(r, "."))
        st_phase_error("Invalid %s allele label", role);
    for (const char *p = r; *p; p++)
        if (*p == ',' || isspace((unsigned char)*p))
            st_phase_error("Invalid %s allele label", role);
    int n = !strcmp(a, ".") ? 1 : 2;
    for (const char *p = a; *p; p++)
    {
        if (isspace((unsigned char)*p))
            st_phase_error("Invalid %s allele label", role);
        if (*p == ',')
            n++;
    }
    if (n > 100)
        st_phase_error("Too many %s alleles", role);
    stp_alleles out = {n, st_zero((size_t)n, sizeof(*out.label))};
    out.label[0] = r;
    const char *start = a;
    for (int j = 1; j < n; j++)
    {
        const char *end = strchr(start, ',');
        size_t length = end ? (size_t)(end - start) : strlen(start);
        if (!length || length > INT_MAX || (length == 1 && *start == '.'))
            st_phase_error("Invalid %s allele label", role);
        char *label = st_zero(length + 1, 1);
        memcpy(label, start, length);
        out.label[j] = label;
        start += length + (end != NULL);
    }
    if (!stp_names_unique(out.label, (size_t)n))
        st_phase_error("Invalid or duplicated %s alleles", role);
    st_release((void *)a);
    return out;
}

static stp_alleles *stp_validate(SEXP x, const char *role,
                                 stp_metadata *metadata)
{
    if (!Rf_inherits(x, "STvcf"))
        st_phase_error("%s must be an STvcf or a loaded dataset name", role);
    SEXP id = stp_field(x, "ID"), ind = stp_field(x, "individual");
    SEXP chr = stp_field(x, "CHROM"), pos = stp_field(x, "POS");
    SEXP levels = stp_field(x, "CHROM.levels");
    SEXP a = stp_field(x, "allele1"), b = stp_field(x, "allele2");
    SEXP da = Rf_getAttrib(a, R_DimSymbol), db = Rf_getAttrib(b, R_DimSymbol);
    int m = Rf_length(id), n = Rf_length(ind);
    if (m < 1 || n < 1 || TYPEOF(chr) != INTSXP ||
        (TYPEOF(pos) != INTSXP && TYPEOF(pos) != REALSXP) ||
        XLENGTH(chr) != m || XLENGTH(pos) != m || TYPEOF(a) != RAWSXP ||
        TYPEOF(b) != RAWSXP || TYPEOF(da) != INTSXP || TYPEOF(db) != INTSXP ||
        XLENGTH(da) != 2 || XLENGTH(db) != 2 || INTEGER(da)[0] != m ||
        INTEGER(da)[1] != n || INTEGER(db)[0] != m || INTEGER(db)[1] != n ||
        XLENGTH(a) != (R_xlen_t)m * n || XLENGTH(b) != XLENGTH(a))
        st_phase_error("Invalid %s dimensions, map or allele matrices", role);
    stp_names(id, role, "ID");
    stp_names(ind, role, "individual");
    stp_names(levels, role, "CHROM.levels");
    metadata->markers = m;
    metadata->samples = n;
    metadata->nlevels = Rf_length(levels);
    metadata->chrom = st_zero((size_t)m, sizeof(*metadata->chrom));
    metadata->pos = st_zero((size_t)m, sizeof(*metadata->pos));
    metadata->id = st_zero((size_t)m, sizeof(*metadata->id));
    metadata->levels =
        st_zero((size_t)metadata->nlevels, sizeof(*metadata->levels));
    for (int k = 0; k < metadata->nlevels; k++)
        metadata->levels[k] = stp_input_text(STRING_ELT(levels, k));
    for (int j = 0; j < m; j++)
    {
        int chromosome = INTEGER(chr)[j];
        double position = stp_number_at(pos, j);
        if (chromosome < 1 || chromosome > XLENGTH(levels) ||
            !R_FINITE(position) || position < 0)
            st_phase_error("Invalid %s dimensions, map or allele matrices",
                           role);
        metadata->chrom[j] = chromosome;
        metadata->pos[j] = position;
        metadata->id[j] = stp_input_text(STRING_ELT(id, j));
    }
    SEXP ref = stp_field(x, "REF"), alt = stp_field(x, "ALT");
    if (TYPEOF(ref) != STRSXP || TYPEOF(alt) != STRSXP || XLENGTH(ref) != m ||
        XLENGTH(alt) != m)
        st_phase_error("Invalid %s allele definitions", role);
    stp_alleles *definitions = st_zero((size_t)m, sizeof(*definitions));
    for (int j = 0; j < m; j++)
    {
        if (STRING_ELT(ref, j) == NA_STRING || STRING_ELT(alt, j) == NA_STRING)
            st_phase_error("Invalid %s allele definitions", role);
        definitions[j] =
            stp_allele_labels(STRING_ELT(ref, j), STRING_ELT(alt, j), role);
        int count = definitions[j].count;
        for (int i = 0; i < n; i++)
        {
            R_xlen_t q = j + (R_xlen_t)m * i;
            if ((RAW(a)[q] != 255 && RAW(a)[q] >= count) ||
                (RAW(b)[q] != 255 && RAW(b)[q] >= count))
                st_phase_error("%s allele index outside REF/ALT at %s", role,
                               CHAR(STRING_ELT(id, j)));
        }
        if (!(j % 1024))
            R_CheckUserInterrupt();
    }
    return definitions;
}

static int stp_double_compare(const void *left, const void *right)
{
    double x = *(const double *)left, y = *(const double *)right;
    return (x > y) - (x < y);
}

static double stp_expanded(const double *value, R_xlen_t length, int index,
                           int chromosomes)
{
    if (length != 1 && length != chromosomes)
        st_phase_error(
            "Window/overlap must be scalar or have one value per chromosome");
    return value[length == 1 ? 0 : index];
}

typedef struct
{
    gs_varset_type *target, *reference;
    int markers, targets, samples, reference_markers;
    const unsigned char *external_a, *external_b;
    int *reference_index, *counts;
    unsigned char *target_code, *reference_code, *decode;
    GS_ALLELE_TYPE *staged1, *staged2;
    R_xlen_t filled;
    int impute;
} stp_native_view;

static void stp_native_metadata(gs_varset_type *x, stp_metadata *out)
{
    if (!x || !x->NoMar || !x->NoInd || x->NoMar > INT_MAX || x->NoInd > 5000 ||
        !st_allele_present(x) || !x->MarName || !x->IndName || !x->NLoc ||
        !x->map || !x->mdl)
        st_phase_error("No valid marker data and linkage map available");
    int m = (int)x->NoMar, n = (int)x->NoInd;
    out->markers = m;
    out->samples = n;
    out->chrom = st_zero((size_t)m, sizeof(*out->chrom));
    out->pos = st_zero((size_t)m, sizeof(*out->pos));
    out->id = st_zero((size_t)m, sizeof(*out->id));
    out->levels = st_zero((size_t)m, sizeof(*out->levels));
    int *seen = st_zero((size_t)m, sizeof(int));
    int *levels = st_zero((size_t)m, sizeof(int)), used = 0;
    for (gs_UNSIGNED l = 0; l < x->NLoc; l++)
    {
        gs_UNSIGNED marker = x->mdl[l];
        if (marker == gs_UNSIGNED_MAX)
            continue;
        if (marker >= x->NoMar || !x->map[l] || seen[marker]++)
            st_phase_error("Invalid or duplicated marker linkage-map entry");
        gs_mappoint_type *point = x->map[l];
        if (point->chrom < 1 || !R_FINITE(point->pos) || point->pos < 0 ||
            strcmp(point->name, x->MarName[marker]))
            st_phase_error("Invalid marker linkage-map entry");
        out->chrom[marker] = point->chrom;
        out->pos[marker] = point->pos;
    }
    for (int j = 0; j < m; j++)
    {
        if (!seen[j])
            st_phase_error("Marker has no linkage-map entry");
        int chr = out->chrom[j], c = 0;
        while (c < used && levels[c] != chr)
            c++;
        if (c == used)
        {
            char *label = st_zero(64, 1);
            levels[used++] = chr;
            snprintf(label, 64, "%d", chr);
            out->levels[c] = label;
        }
        out->chrom[j] = c + 1;
        out->id[j] = x->MarName[j];
    }
    out->nlevels = used;
    if (!stp_names_unique(out->id, (size_t)m))
        st_phase_error("Invalid or duplicated native ID");
    if (!stp_names_unique((const char *const *)x->IndName, (size_t)n))
        st_phase_error("Invalid or duplicated native individual");
    st_release(seen);
    st_release(levels);
}

static int stp_native_codes(gs_varset_type *x, int marker,
                            unsigned char *encode, unsigned char *decode)
{
    unsigned char used[255] = {0};
    memset(encode, 255, 256);
    for (gs_UNSIGNED i = 0; i < x->NoInd; i++)
    {
        size_t at = (size_t)marker * x->NoInd + i;
        int a = st_allele_get(x, at, 0), b = st_allele_get(x, at, 1);
        if (a != 255)
            used[a] = 1;
        if (b != 255)
            used[b] = 1;
    }
    int count = 0;
    for (int a = 0; a < 255; a++)
        if (used[a])
        {
            if (count == 100)
                st_phase_error("Too many observed phasing alleles");
            encode[a] = (unsigned char)count;
            decode[count++] = (unsigned char)a;
        }
    if (!count)
    {
        decode[0] = 1;
        count = 1;
    }
    return count;
}

static int stp_native_read(const stp_native_view *v, int marker, int sample,
                           int hom)
{
    int a, b, code;
    if (sample < v->targets)
    {
        size_t at = (size_t)marker * v->target->NoInd + sample;
        a = st_allele_get(v->target, at, 0);
        b = st_allele_get(v->target, at, 1);
        code = hom ? b : a;
        if (a == 255 || b == 255)
            return 255;
        return v->target_code[(size_t)marker * 256 + code];
    }
    int row = v->reference_index[marker], i = sample - v->targets;
    if (v->reference)
    {
        size_t at = (size_t)row * v->reference->NoInd + i;
        code = hom ? st_allele_get(v->reference, at, 1)
                   : st_allele_get(v->reference, at, 0);
    }
    else
    {
        size_t at = row + (size_t)v->reference_markers * i;
        code = (hom ? v->external_b : v->external_a)[at];
    }
    return v->reference_code[(size_t)marker * 256 + code];
}

static void stp_native_prepare(stp_native_view *v, const stp_metadata *metadata,
                               const stp_metadata *refmeta,
                               const stp_alleles *refdefs, int phased)
{
    int m = v->markers;
    v->target_code = st_zero((size_t)m, 256);
    v->reference_code = st_zero((size_t)m, 256);
    v->decode = st_zero((size_t)m, 100);
    v->counts = st_zero((size_t)m, sizeof(int));
    v->reference_index = st_zero((size_t)m, sizeof(int));
    stp_name_entry *reference_names = NULL;
    if (refmeta)
    {
        reference_names =
            st_zero((size_t)refmeta->markers, sizeof(*reference_names));
        for (int r = 0; r < refmeta->markers; r++)
        {
            reference_names[r].name = refmeta->id[r];
            reference_names[r].index = r;
        }
        qsort(reference_names, (size_t)refmeta->markers,
              sizeof(*reference_names), stp_name_compare);
    }
    for (int j = 0; j < m; j++)
    {
        if (!(j % 1024))
            R_CheckUserInterrupt();
        unsigned char *enc = v->target_code + (size_t)j * 256;
        unsigned char *dec = v->decode + (size_t)j * 100;
        int count = stp_native_codes(v->target, j, enc, dec),
            original_count = count;
        if (refmeta)
        {
            stp_name_entry key = {metadata->id[j], 0};
            const stp_name_entry *matched =
                bsearch(&key, reference_names, (size_t)refmeta->markers,
                        sizeof(*reference_names), stp_name_compare);
            int r = matched ? matched->index : -1;
            if (r < 0)
                st_phase_error("Reference panel lacks target marker: %s",
                               v->target->MarName[j]);
            v->reference_index[j] = r;
            if (strcmp(metadata->levels[metadata->chrom[j] - 1],
                       refmeta->levels[refmeta->chrom[r] - 1]))
                st_phase_error(
                    "Reference chromosome labels disagree with target markers");
            if (fabs(metadata->pos[j] - refmeta->pos[r]) > 1e-6)
                st_phase_error(
                    "Reference genetic positions disagree with target markers");
            unsigned char *renc = v->reference_code + (size_t)j * 256;
            memset(renc, 255, 256);
            int reference_count = 0;
            if (v->reference)
            {
                unsigned char rdecode[100];
                reference_count =
                    stp_native_codes(v->reference, r, renc, rdecode);
                if (reference_count > count)
                    count = reference_count;
            }
            else
            {
                const stp_alleles *rd = &refdefs[r];
                reference_count = rd->count;
                const char *labels[100];
                char target_labels[100][16];
                for (int k = 0; k < count; k++)
                {
                    snprintf(target_labels[k], sizeof(target_labels[k]), "%d",
                             k + 1);
                    labels[k] = target_labels[k];
                }
                int shared = 0;
                for (int k = 0; k < reference_count; k++)
                {
                    int code = -1;
                    for (int l = 0; l < count; l++)
                        if (!strcmp(labels[l], rd->label[k]))
                            code = l;
                    if (code >= 0)
                        shared = 1;
                    else
                    {
                        if (count == 100)
                            st_phase_error("Too many combined phasing alleles");
                        code = count++;
                        labels[code] = rd->label[k];
                    }
                    renc[k] = (unsigned char)code;
                }
                if (!shared)
                    st_phase_error(
                        "No shared allele definition at reference marker: %s",
                        v->target->MarName[j]);
            }
            unsigned char used[255] = {0};
            int maximum = -1;
            for (int k = 0; k < original_count; k++)
            {
                used[dec[k]] = 1;
                if (dec[k] > maximum)
                    maximum = dec[k];
            }
            for (int k = original_count; k < count; k++)
            {
                int code = maximum < 1 ? 1 : maximum + 1;
                if (code >= 255)
                {
                    code = 1;
                    while (code < 255 && used[code])
                        code++;
                }
                if (code >= 255)
                    st_phase_error("No available internal allele code");
                dec[k] = (unsigned char)code;
                used[code] = 1;
                if (code > maximum)
                    maximum = code;
            }
            for (int i = v->targets; i < v->samples; i++)
            {
                int a = stp_native_read(v, j, i, 0),
                    b = stp_native_read(v, j, i, 1);
                if (a == 255 || b == 255)
                    st_phase_error("Reference calls must be complete");
                if (a != b && !phased)
                    st_phase_error("Heterozygous reference genotypes must be "
                                   "declared phased");
            }
        }
        v->counts[j] = count;
    }
    st_release(reference_names);
    if ((size_t)m > SIZE_MAX / (size_t)v->targets)
        st_phase_error("Phasing staging size overflow");
    size_t cells = (size_t)m * v->targets;
    v->staged1 = st_zero(cells, sizeof(*v->staged1));
    v->staged2 = st_zero(cells, sizeof(*v->staged2));
}

typedef struct
{
    int chromosome, markers, common, windows;
    double span, median, window, overlap;
    char reason[128];
} stp_window_result;

static stp_window_result *stp_windows(const stp_metadata *input,
                                      const stp_controls *options,
                                      int chromosomes,
                                      const stp_native_view *native)
{
    int automatic = options->automatic;
    stp_window_result *out = st_zero((size_t)chromosomes, sizeof(*out));
    const int *chr = input->chrom;
    const double *pos = input->pos;
    int m = input->markers, n = native->samples;
    double *common =
        automatic ? (double *)st_zero((size_t)m, sizeof(double)) : NULL;
    double *gaps =
        automatic ? (double *)st_zero((size_t)m, sizeof(double)) : NULL;
    float carrier_product = (float)n * (float)options->value[STP_RARE];
    int limit = (int)floor(carrier_product);
    if (limit < 3)
        limit = 3;
    int start = 0;
    for (int c = 0; c < chromosomes; c++)
    {
        int end = start + 1;
        while (end < m && chr[end] == chr[start])
            end++;
        out[c].chromosome = chr[start];
        double window, overlap;
        const char *reason = "explicit";
        char reason_buffer[128];
        if (automatic)
        {
            int count = 0;
            double span = pos[end - 1] - pos[start];
            if (end - start < 2 || span <= 0)
                st_phase_error("Window selection requires an ordered, nonzero "
                               "genetic map on "
                               "chromosome %d",
                               chr[start]);
            for (int j = start; j < end; j++)
            {
                if (j > start && pos[j] < pos[j - 1])
                    st_phase_error(
                        "Window selection requires an ordered genetic map");
                int carriers[100] = {0}, common_alleles = 0;
                for (int i = 0; i < n; i++)
                {
                    int x = stp_native_read(native, j, i, 0);
                    int y = stp_native_read(native, j, i, 1);
                    if (x != 255)
                        carriers[x]++;
                    if (y != 255 && y != x)
                        carriers[y]++;
                }
                for (int k = 0; k < 100; k++)
                    common_alleles += carriers[k] > limit;
                if (common_alleles >= 2)
                    common[count++] = pos[j];
                if (!(j % 1024))
                    R_CheckUserInterrupt();
            }
            window = R_PosInf;
            overlap = 2;
            reason = "fewer_than_two_common_marker_blocks";
            if (count >= 800)
            {
                double full_span = 0, overlap_span = 0;
                for (int j = 400; j < count; j++)
                    full_span = fmax(full_span, common[j] - common[j - 400]);
                for (int j = 40; j < count; j++)
                    overlap_span =
                        fmax(overlap_span, common[j] - common[j - 40]);
                window = fmax(40, full_span * (1 + 1e-6));
                overlap = fmax(2, overlap_span * (1 + 1e-6));
                window = fmax(window, 2 * overlap);
                reason = "common_marker_spacing";
                if (window >= span)
                {
                    window = R_PosInf;
                    reason = "required_context_spans_chromosome";
                }
            }
            if (!R_FINITE(window))
                overlap = 2;
            double median = NA_REAL;
            if (count > 1)
            {
                for (int j = 1; j < count; j++)
                    gaps[j - 1] = common[j] - common[j - 1];
                qsort(gaps, count - 1, sizeof(double), stp_double_compare);
                int ngaps = count - 1;
                median = ngaps % 2
                             ? gaps[ngaps / 2]
                             : (gaps[ngaps / 2 - 1] + gaps[ngaps / 2]) / 2;
            }
            out[c].markers = end - start;
            out[c].common = count;
            out[c].span = span;
            out[c].median = median;
            if (options->noverlap)
            {
                overlap = stp_expanded(options->overlap, options->noverlap, c,
                                       chromosomes);
                window = fmax(window, 1.1 * overlap);
                snprintf(reason_buffer, sizeof(reason_buffer),
                         "%s;explicit_overlap", reason);
                reason = reason_buffer;
            }
            if (!R_FINITE(window) && end - start > options->value[STP_MARKERS])
                st_phase_error(
                    "Automatic whole-chromosome choice exceeds window.markers; "
                    "supply explicit windows");
        }
        else
        {
            window =
                stp_expanded(options->window, options->nwindow, c, chromosomes);
            overlap = stp_expanded(options->overlap, options->noverlap, c,
                                   chromosomes);
        }
        out[c].window = window;
        out[c].overlap = overlap;
        snprintf(out[c].reason, sizeof(out[c].reason), "%s", reason);
        start = end;
    }
    st_release(common);
    st_release(gaps);
    return out;
}

static SEXP stp_windows_result(const stp_window_result *windows, int automatic,
                               int chromosomes)
{
    const char *auto_names[] = {"chromosome",
                                "markers",
                                "common.markers",
                                "span.cm",
                                "median.common.gap.cm",
                                "minimum.common",
                                "window",
                                "overlap",
                                "reason"};
    const char *explicit_names[] = {"chromosome", "window", "overlap",
                                    "reason"};
    int fields = automatic ? 9 : 4, wi = automatic ? 6 : 1, oi = wi + 1,
        ri = wi + 2;
    SEXP out =
        PROTECT(stp_list(fields, automatic ? auto_names : explicit_names));
    for (int k = 0; k < fields; k++)
    {
        SEXPTYPE type = (k == 0 || (automatic && (k == 1 || k == 2 || k == 5)))
                            ? INTSXP
                            : REALSXP;
        if (k == ri)
            type = STRSXP;
        SET_VECTOR_ELT(out, k, Rf_allocVector(type, chromosomes));
    }
    SEXP rownames = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(rownames)[0] = NA_INTEGER;
    INTEGER(rownames)[1] = -chromosomes;
    Rf_setAttrib(out, R_RowNamesSymbol, rownames);
    SEXP cls = PROTECT(Rf_mkString("data.frame"));
    Rf_setAttrib(out, R_ClassSymbol, cls);
    for (int c = 0; c < chromosomes; c++)
    {
        const stp_window_result *w = &windows[c];
        INTEGER(VECTOR_ELT(out, 0))[c] = w->chromosome;
        if (automatic)
        {
            INTEGER(VECTOR_ELT(out, 1))[c] = w->markers;
            INTEGER(VECTOR_ELT(out, 2))[c] = w->common;
            REAL(VECTOR_ELT(out, 3))[c] = w->span;
            REAL(VECTOR_ELT(out, 4))[c] = w->median;
            INTEGER(VECTOR_ELT(out, 5))[c] = 400;
        }
        REAL(VECTOR_ELT(out, wi))[c] = w->window;
        REAL(VECTOR_ELT(out, oi))[c] = w->overlap;
        SET_STRING_ELT(VECTOR_ELT(out, ri), c, Rf_mkChar(w->reason));
    }
    UNPROTECT(3);
    return out;
}

typedef struct
{
    SEXP dataset, reference, requested, impute, commit;
} st_phase_call;

static SEXP st_phase_body(void *data)
{
    st_phase_call *c = data;
    const char *name = NULL;
    st_busy = 1;
    if (!st_phase_step1_dataset_name(c->dataset, &name))
        st_phase_error("Invalid source data set name");
    st_phase_dataset = name;
    stp_native_view v = {0};
    v.target = st_phase_step1_find_dataset(name);
    v.impute = stp_logical(c->impute, "impute");
    int commit = stp_logical(c->commit, "commit");
    if (!v.target)
        st_phase_error("Unknown source data set");
    stp_metadata metadata = {0}, refmeta = {0};
    stp_native_metadata(v.target, &metadata);
    v.markers = (int)v.target->NoMar;
    v.targets = (int)v.target->NoInd;
    v.samples = v.targets;
    stp_alleles *refdefs = NULL;
    int reference_phased = 0;
    int with_reference = c->reference != R_NilValue;
    if (with_reference)
    {
        if (TYPEOF(c->reference) == STRSXP && XLENGTH(c->reference) == 1 &&
            STRING_ELT(c->reference, 0) != NA_STRING)
        {
            v.reference =
                st_phase_step1_find_dataset(CHAR(STRING_ELT(c->reference, 0)));
            if (!v.reference)
                st_phase_error("Unknown reference data set");
            if (commit && v.reference == v.target)
                st_phase_error("Reference and target must be different data "
                               "sets when commit=TRUE");
            stp_native_metadata(v.reference, &refmeta);
            reference_phased = v.reference->phased;
        }
        else
        {
            refdefs = stp_validate(c->reference, "reference", &refmeta);
            v.external_a = RAW(stp_field(c->reference, "allele1"));
            v.external_b = RAW(stp_field(c->reference, "allele2"));
            reference_phased = Rf_asLogical(Rf_getAttrib(
                                   c->reference, Rf_install("phased"))) == TRUE;
        }
        v.reference_markers = refmeta.markers;
        int rn = refmeta.samples;
        if (rn > 5000 - v.samples)
            st_phase_error("Phasing supports at most 5000 target and reference "
                           "samples combined");
        v.samples += rn;
    }
    if (v.markers < 2 || v.samples < 2 || v.samples > 5000)
        st_phase_error("Phasing requires at least 2 markers and 2 to 5000 "
                       "samples combined");
    stp_native_prepare(&v, &metadata, with_reference ? &refmeta : NULL, refdefs,
                       reference_phased);
    const int *chr = metadata.chrom;
    const double *cm = metadata.pos;
    int chromosomes = 0;
    for (int j = 0; j < v.markers; j++)
    {
        if (j && chr[j] < chr[j - 1])
            st_phase_error("Phasing requires markers ordered by chromosome");
        chromosomes += !j || chr[j] != chr[j - 1];
    }
    stp_controls options;
    SEXP seed_input = R_NilValue;
    stp_options(c->requested, &options, &seed_input);
    stp_window_result *windows =
        stp_windows(&metadata, &options, chromosomes, &v);
    int threads = 1;
#ifdef MT
    threads = omp_get_max_threads();
#endif
    if (threads < 1)
        st_phase_error("Invalid SelectionTools thread setting");
    int actual_seed = (int)options.value[STP_SEED];
    if (!options.seed_set)
    {
        GetRNGstate();
        actual_seed = 1 + (int)R_unif_index((double)INT_MAX);
        PutRNGstate();
    }
    int *bp = st_zero((size_t)v.markers, sizeof(int));
    for (int j = 0; j < v.markers; j++)
    {
        double coordinate = floor(cm[j] * 1000000 + .5) + 1;
        if (j && chr[j] == chr[j - 1])
            coordinate = fmax(coordinate, (double)bp[j - 1] + 1);
        if (!R_FINITE(coordinate) || coordinate < 1 || coordinate > INT_MAX)
            st_phase_error("Map cannot be represented by integer positions");
        bp[j] = (int)coordinate;
    }
    j_random master;
    jr_seed(&master, (uint64_t)(int64_t)actual_seed);
    atomic_store(&st_max_team, 1);
    int total_windows = 0, column = 0;
    const int indices[] = {STP_BURNIN, STP_ITERATIONS, STP_STATES, STP_RARE,
                           STP_SCALE,  STP_NE,         STP_BUFFER, STP_LR};
    for (int start = 0; start < v.markers; column++)
    {
        R_CheckUserInterrupt();
        int chromosome = chr[start], end = start + 1;
        while (end < v.markers && chr[end] == chromosome)
            end++;
        int m = end - start;
        double o[14] = {0};
        o[0] = windows[column].window;
        o[1] = windows[column].overlap;
        for (int k = 0; k < 8; k++)
            o[k + 2] = options.value[indices[k]];
        o[10] = options.value[STP_EM];
        o[11] = options.error_set ? options.value[STP_ERR] : R_NaN;
        o[12] = options.value[STP_MARKERS];
        o[13] = threads;
        er_set_threads(threads);
        st_phase_threads = threads;
        np_input raw = {0};
        raw.M = m;
        raw.N = v.samples;
        raw.T = v.targets;
        raw.seed = (uint64_t)(int64_t)actual_seed;
        raw.burn = (int)o[2];
        raw.iterations = (int)o[3];
        raw.K = (int)o[4];
        raw.rare = (float)o[5];
        raw.scale = (float)o[6];
        raw.ne = (float)o[7];
        raw.buffer = (float)o[8];
        raw.lr = (float)o[9];
        raw.em = (int)o[10];
        raw.bp = bp + start;
        raw.gt = st_zero((size_t)m * 2 * v.samples, sizeof(int));
        int *output = st_zero((size_t)m * 2 * v.samples, sizeof(int));
        for (int j = 0; j < m; j++)
            for (int i = 0; i < v.samples; i++)
            {
                size_t at = (size_t)j * 2 * v.samples + 2 * i;
                int a = stp_native_read(&v, start + j, i, 0),
                    b = stp_native_read(&v, start + j, i, 1);
                raw.gt[at] = a == 255 ? -1 : a;
                raw.gt[at + 1] = b == 255 ? -1 : b;
            }
        int before = total_windows;
        if (!st_phase_windowed(&raw, v.counts + start, raw.bp, cm + start, m,
                               (st_window_options){o[0], o[1], (int)o[12]},
                               o[11], &master, chromosome, &total_windows,
                               output, NULL))
            st_phase_error("Invalid data or options for chromosome %d",
                           chromosome);
        windows[column].windows = total_windows - before;
        for (int j = 0; j < m; j++)
            for (int i = 0; i < v.samples; i++)
            {
                int marker = start + j;
                size_t at = (size_t)j * 2 * v.samples + 2 * i;
                int a = stp_native_read(&v, marker, i, 0),
                    b = stp_native_read(&v, marker, i, 1);
                int x = output[at] < 0 ? 255 : output[at],
                    y = output[at + 1] < 0 ? 255 : output[at + 1];
                if ((x != 255 && x >= v.counts[marker]) ||
                    (y != 255 && y >= v.counts[marker]))
                    st_phase_error("Invalid phased allele");
                if (i >= v.targets)
                {
                    if (x != a || y != b)
                        st_phase_error("Phasing changed a reference haplotype");
                    continue;
                }
                int missing = a == 255 || b == 255;
                if (!missing && !((x == a && y == b) || (x == b && y == a)))
                    st_phase_error(
                        "Phasing attempted to change an observed genotype");
                if (missing && !v.impute)
                {
                    x = a;
                    y = b;
                }
                v.filled += missing && x != 255 && y != 255;
                size_t target = (size_t)marker * v.targets + i;
                unsigned char *decode = v.decode + (size_t)marker * 100;
                v.staged1[target] = x == 255 ? GS_ALLELE_MISSING : decode[x];
                v.staged2[target] = y == 255 ? GS_ALLELE_MISSING : decode[y];
            }
        st_release(raw.gt);
        st_release(output);
        start = end;
    }
    if (commit)
    {
        st_allele_replace(v.target, v.staged1, v.staged2);
        v.target->phased = 1;
        return R_NilValue;
    }

    SEXP per_chromosome = PROTECT(Rf_allocVector(INTSXP, chromosomes));
    SEXP chromosome_names = PROTECT(Rf_allocVector(STRSXP, chromosomes));
    int np = 2;
    for (int i = 0; i < chromosomes; i++)
    {
        INTEGER(per_chromosome)[i] = windows[i].windows;
        char label[32];
        snprintf(label, sizeof(label), "%d", windows[i].chromosome);
        SET_STRING_ELT(chromosome_names, i, Rf_mkChar(label));
    }
    Rf_setAttrib(per_chromosome, R_NamesSymbol, chromosome_names);
    const char *keys[] = {
        "actual.seed",        "windows",           "maximum.worker.team",
        "chromosome.windows", "filled.genotypes",  "committed",
        "window.settings",    "requested.threads", "reference"};
    SEXP diagnostics = PROTECT(stp_list(9, keys));
    np++;
    SET_VECTOR_ELT(diagnostics, 0,
                   options.seed_set ? seed_input
                                    : Rf_ScalarInteger(actual_seed));
    SET_VECTOR_ELT(diagnostics, 1, Rf_ScalarInteger(total_windows));
    SET_VECTOR_ELT(diagnostics, 2, Rf_ScalarInteger(atomic_load(&st_max_team)));
    SET_VECTOR_ELT(diagnostics, 3, per_chromosome);
    SET_VECTOR_ELT(diagnostics, 4, Rf_ScalarReal((double)v.filled));
    SET_VECTOR_ELT(diagnostics, 5, Rf_ScalarLogical(FALSE));
    SET_VECTOR_ELT(diagnostics, 6,
                   stp_windows_result(windows, options.automatic, chromosomes));
    SET_VECTOR_ELT(diagnostics, 7, Rf_ScalarInteger(threads));
    if (with_reference)
    {
        const char *rkeys[] = {"samples", "haplotypes", "unchanged",
                               "ignored.markers"};
        SEXP info = PROTECT(stp_list(4, rkeys));
        np++;
        SET_VECTOR_ELT(info, 0, Rf_ScalarInteger(v.samples - v.targets));
        SET_VECTOR_ELT(info, 1, Rf_ScalarInteger(2 * (v.samples - v.targets)));
        SET_VECTOR_ELT(info, 2, Rf_ScalarLogical(TRUE));
        SET_VECTOR_ELT(info, 3,
                       Rf_ScalarInteger(v.reference_markers - v.markers));
        SET_VECTOR_ELT(diagnostics, 8, info);
    }
    UNPROTECT(np);
    return diagnostics;
}

SEXP st_phase(SEXP dataset, SEXP reference, SEXP requested, SEXP impute,
              SEXP commit)
{
    if (st_busy)
        Rf_error("st.phase is already running");
    st_phase_call c = {dataset, reference, requested, impute, commit};
    return R_ExecWithCleanup(st_phase_body, &c, st_cleanup, NULL);
}
