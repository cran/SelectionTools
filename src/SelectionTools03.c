/* st.phase.v1 */

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

#include "SelectionToolsAllele.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include "SelectionToolsDataSet.h"

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

static int gs_map_pool_valid(const gs_mappoint_type *map_S,
                             gs_UNSIGNED nloc, const char *pool,
                             size_t pool_n)
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

  if ((NULL == gv) || (NULL == gv->All1) || (NULL == gv->All2) ||
      (0 == gv->NoInd) || (0 == gv->NoMar) ||
      (0 != gs_size_mul((size_t)gv->NoInd, (size_t)gv->NoMar, &n)))
    return 0;

  for (size_t idx = 0; idx < n; idx++)
  {
    GS_ALLELE_TYPE a = gv->All1[idx];
    GS_ALLELE_TYPE b = gv->All2[idx];
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

  if (!gs_map_pool_valid(gv->map_S, gv->NLoc, gv->map_str_S,
                         gv->map_str_N))
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
        gs_info(gv, -2, "Marker data are not completely covered by map");
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

#define ST_PHASE_STATE_MISSING ((uint8_t)255u)

typedef struct
{
  size_t n_ind;
  size_t n_marker;
  size_t n_chrom;

  gs_UNSIGNED *source_marker;
  int *chrom_label;
  double *pos_cm;

  size_t *chrom_start;
  size_t *chrom_length;

  size_t *allele_offset;
  uint8_t *n_alleles;
  GS_ALLELE_TYPE *allele_value;
  size_t n_allele_value;

  uint8_t *geno1;
  uint8_t *geno2;
} st_phase_data;

static int st_phase_size_mul(size_t a, size_t b, size_t *out)
{
  if ((0 != a) && (b > SIZE_MAX / a))
    return -1;

  *out = a * b;
  return 0;
}

static void st_phase_data_init(st_phase_data *d)
{
  if (NULL != d)
    memset(d, 0, sizeof(*d));
}

static void st_phase_data_free(st_phase_data *d)
{
  if (NULL == d)
    return;

  free(d->source_marker);
  free(d->chrom_label);
  free(d->pos_cm);
  free(d->chrom_start);
  free(d->chrom_length);
  free(d->allele_offset);
  free(d->n_alleles);
  free(d->allele_value);
  free(d->geno1);
  free(d->geno2);

  st_phase_data_init(d);
}

static int st_phase_marker_alleles(const gs_varset_type *gv,
                                   gs_UNSIGNED marker,
                                   unsigned char present[256],
                                   uint8_t *n_alleles)
{
  unsigned int n = 0;

  memset(present, 0, 256u * sizeof(*present));

  for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
  {
    size_t idx = (size_t)marker * (size_t)gv->NoInd + (size_t)i;
    GS_ALLELE_TYPE a = gv->All1[idx];
    GS_ALLELE_TYPE b = gv->All2[idx];

    if (GS_ALLELE_MISSING == a)
      continue;

    if ((0 == a) || (0 == b) || (GS_ALLELE_MISSING == b))
      return -1;

    if (!present[(unsigned int)a])
    {
      present[(unsigned int)a] = 1;
      n++;
    }
    if (!present[(unsigned int)b])
    {
      present[(unsigned int)b] = 1;
      n++;
    }
  }

  if (n > 254u)
    return -1;

  *n_alleles = (uint8_t)n;
  return 0;
}

static int st_phase_build_data(gs_varset_type *gv, st_phase_data *d)
{
  gs_UNSIGNED *map_storage = NULL;
  unsigned char present[256];
  size_t n_genotype = 0;
  size_t allele_total = 0;
  size_t chrom_count = 0;
  int previous_chrom = INT_MIN;
  double previous_pos = 0.0;

  if ((NULL == gv) || (NULL == d))
    return -1;

  st_phase_data_init(d);

  if (!gs_marker_genotypes_valid(gv))
  {
    gs_info(gv, -2, "Invalid marker genotype representation");
    return -1;
  }

  map_storage = gs_map_storage_index(gv, 1, 1);
  if (NULL == map_storage)
    return -1;
  free(map_storage);
  map_storage = NULL;

  if ((0 == gv->NLoc) || (gv->NLoc != gv->NoMar) || (NULL == gv->map) ||
      (NULL == gv->mdl))
  {
    gs_info(gv, -2, "Linkage map is not matching marker data");
    return -1;
  }

  if (0 !=
      st_phase_size_mul((size_t)gv->NoInd, (size_t)gv->NLoc, &n_genotype))
  {
    gs_info(gv, -2, "Memory allocation error");
    return -1;
  }

  d->n_ind = (size_t)gv->NoInd;
  d->n_marker = (size_t)gv->NLoc;

  d->source_marker = malloc(d->n_marker * sizeof(*d->source_marker));
  d->chrom_label = malloc(d->n_marker * sizeof(*d->chrom_label));
  d->pos_cm = malloc(d->n_marker * sizeof(*d->pos_cm));
  d->allele_offset =
      malloc((d->n_marker + 1u) * sizeof(*d->allele_offset));
  d->n_alleles = malloc(d->n_marker * sizeof(*d->n_alleles));
  d->geno1 = malloc(n_genotype * sizeof(*d->geno1));
  d->geno2 = malloc(n_genotype * sizeof(*d->geno2));

  if ((NULL == d->source_marker) || (NULL == d->chrom_label) ||
      (NULL == d->pos_cm) || (NULL == d->allele_offset) ||
      (NULL == d->n_alleles) || (NULL == d->geno1) || (NULL == d->geno2))
  {
    gs_info(gv, -2, "Memory allocation error");
    goto fail;
  }

  d->allele_offset[0] = 0;

  for (size_t l = 0; l < d->n_marker; l++)
  {
    gs_UNSIGNED marker = gv->mdl[l];
    gs_mappoint_type *mp = gv->map[l];
    uint8_t n_alleles = 0;

    if ((NULL == mp) || (marker >= gv->NoMar) || (mp->nr != marker) ||
        (mp->chrom < 1) || !isfinite(mp->pos) || (mp->pos < 0.0))
    {
      gs_info(gv, -2, "Invalid linkage map for phasing");
      goto fail;
    }

    if (0 == l)
    {
      chrom_count = 1;
    }
    else if (mp->chrom < previous_chrom)
    {
      gs_info(gv, -2, "Linkage map is not ordered by chromosome");
      goto fail;
    }
    else if (mp->chrom == previous_chrom)
    {
      if (mp->pos < previous_pos)
      {
        gs_info(gv, -2, "Linkage map is not ordered by position");
        goto fail;
      }
    }
    else
    {
      chrom_count++;
    }

    if (0 != st_phase_marker_alleles(gv, marker, present, &n_alleles))
    {
      gs_info(gv, -2, "Invalid marker alleles for phasing");
      goto fail;
    }

    if (allele_total > SIZE_MAX - (size_t)n_alleles)
    {
      gs_info(gv, -2, "Memory allocation error");
      goto fail;
    }

    d->source_marker[l] = marker;
    d->chrom_label[l] = mp->chrom;
    d->pos_cm[l] = mp->pos;
    d->n_alleles[l] = n_alleles;
    allele_total += (size_t)n_alleles;
    d->allele_offset[l + 1u] = allele_total;

    previous_chrom = mp->chrom;
    previous_pos = mp->pos;
  }

  d->n_chrom = chrom_count;
  d->n_allele_value = allele_total;

  d->chrom_start = malloc(d->n_chrom * sizeof(*d->chrom_start));
  d->chrom_length = malloc(d->n_chrom * sizeof(*d->chrom_length));

  if (allele_total > 0)
    d->allele_value = malloc(allele_total * sizeof(*d->allele_value));

  if ((NULL == d->chrom_start) || (NULL == d->chrom_length) ||
      ((allele_total > 0) && (NULL == d->allele_value)))
  {
    gs_info(gv, -2, "Memory allocation error");
    goto fail;
  }

  {
    size_t c = 0;
    size_t start = 0;

    d->chrom_start[0] = 0;

    for (size_t l = 1; l < d->n_marker; l++)
      if (d->chrom_label[l] != d->chrom_label[l - 1u])
      {
        d->chrom_length[c] = l - start;
        c++;
        start = l;
        d->chrom_start[c] = l;
      }

    d->chrom_length[c] = d->n_marker - start;
  }

  for (size_t l = 0; l < d->n_marker; l++)
  {
    gs_UNSIGNED marker = d->source_marker[l];
    uint8_t state_by_value[256];
    uint8_t state = 0;
    size_t offset = d->allele_offset[l];

    memset(state_by_value, ST_PHASE_STATE_MISSING, sizeof(state_by_value));

    if (0 !=
        st_phase_marker_alleles(gv, marker, present, &d->n_alleles[l]))
    {
      gs_info(gv, -2, "Invalid marker alleles for phasing");
      goto fail;
    }

    for (unsigned int value = 1; value < 255u; value++)
      if (present[value])
      {
        state_by_value[value] = state;
        d->allele_value[offset + (size_t)state] = (GS_ALLELE_TYPE)value;
        state++;
      }

    if (state != d->n_alleles[l])
    {
      gs_info(gv, -2, "Internal allele-state construction error");
      goto fail;
    }

    for (size_t i = 0; i < d->n_ind; i++)
    {
      size_t source_idx = (size_t)marker * d->n_ind + i;
      size_t target_idx = l * d->n_ind + i;
      GS_ALLELE_TYPE a = gv->All1[source_idx];
      GS_ALLELE_TYPE b = gv->All2[source_idx];
      uint8_t x, y;

      if (GS_ALLELE_MISSING == a)
      {
        d->geno1[target_idx] = ST_PHASE_STATE_MISSING;
        d->geno2[target_idx] = ST_PHASE_STATE_MISSING;
        continue;
      }

      x = state_by_value[(unsigned int)a];
      y = state_by_value[(unsigned int)b];

      if ((ST_PHASE_STATE_MISSING == x) || (ST_PHASE_STATE_MISSING == y))
      {
        gs_info(gv, -2, "Internal allele-state lookup error");
        goto fail;
      }

      if (x <= y)
      {
        d->geno1[target_idx] = x;
        d->geno2[target_idx] = y;
      }
      else
      {
        d->geno1[target_idx] = y;
        d->geno2[target_idx] = x;
      }
    }
  }

  return 0;

fail:
  st_phase_data_free(d);
  return -1;
}

typedef struct
{
  size_t n_hap;
  size_t position;
  size_t *order;
  size_t *divergence;
  size_t *next_order;
  size_t *next_divergence;
} st_phase_pbwt;

static void st_phase_pbwt_init_empty(st_phase_pbwt *p)
{
  if (NULL != p)
    memset(p, 0, sizeof(*p));
}

static void st_phase_pbwt_free(st_phase_pbwt *p)
{
  if (NULL == p)
    return;

  free(p->order);
  free(p->divergence);
  free(p->next_order);
  free(p->next_divergence);
  st_phase_pbwt_init_empty(p);
}

static int st_phase_pbwt_init(st_phase_pbwt *p, size_t n_hap)
{
  if ((NULL == p) || (0 == n_hap) || (n_hap > SIZE_MAX / sizeof(size_t)))
    return -1;

  st_phase_pbwt_init_empty(p);

  p->order = malloc(n_hap * sizeof(*p->order));
  p->divergence = calloc(n_hap, sizeof(*p->divergence));
  p->next_order = malloc(n_hap * sizeof(*p->next_order));
  p->next_divergence = malloc(n_hap * sizeof(*p->next_divergence));

  if ((NULL == p->order) || (NULL == p->divergence) ||
      (NULL == p->next_order) || (NULL == p->next_divergence))
  {
    st_phase_pbwt_free(p);
    return -1;
  }

  p->n_hap = n_hap;
  p->position = 0;

  for (size_t i = 0; i < n_hap; i++)
    p->order[i] = i;

  return 0;
}

static int st_phase_pbwt_update(st_phase_pbwt *p, const uint8_t *allele)
{
  size_t n_zero = 0;
  size_t iz, io;
  size_t p_zero, p_one;
  size_t *tmp;

  if ((NULL == p) || (NULL == allele) || (0 == p->n_hap) ||
      (NULL == p->order) || (NULL == p->divergence) ||
      (NULL == p->next_order) || (NULL == p->next_divergence) ||
      (SIZE_MAX == p->position))
    return -1;

  for (size_t h = 0; h < p->n_hap; h++)
  {
    if (allele[h] > 1u)
      return -1;
    if (0u == allele[h])
      n_zero++;
  }

  iz = 0;
  io = n_zero;
  p_zero = p->position + 1u;
  p_one = p->position + 1u;

  for (size_t i = 0; i < p->n_hap; i++)
  {
    size_t h = p->order[i];
    size_t d = p->divergence[i];

    if (h >= p->n_hap)
      return -1;

    if (d > p_zero)
      p_zero = d;
    if (d > p_one)
      p_one = d;

    if (0u == allele[h])
    {
      p->next_order[iz] = h;
      p->next_divergence[iz] = p_zero;
      iz++;
      p_zero = 0;
    }
    else
    {
      p->next_order[io] = h;
      p->next_divergence[io] = p_one;
      io++;
      p_one = 0;
    }
  }

  if ((iz != n_zero) || (io != p->n_hap))
    return -1;

  tmp = p->order;
  p->order = p->next_order;
  p->next_order = tmp;

  tmp = p->divergence;
  p->divergence = p->next_divergence;
  p->next_divergence = tmp;

  p->position++;
  return 0;
}

typedef struct
{
  size_t hap1;
  size_t hap2;
  size_t start;
  size_t end;
} st_phase_ibs_match;

typedef int (*st_phase_ibs_emit_fn)(void *context,
                                    const st_phase_ibs_match *match);

typedef struct
{
  size_t hap;
  size_t ibs_length;
} st_phase_candidate;

typedef struct
{
  size_t n_marker;
  size_t n_hap;
  size_t max_states;
  size_t n_cell;
  size_t n_slot;
  size_t *count;
  st_phase_candidate *state;
} st_phase_candidate_set;

static int st_phase_candidate_better(const st_phase_candidate *a,
                                     const st_phase_candidate *b)
{
  if (a->ibs_length > b->ibs_length)
    return 1;
  if (a->ibs_length < b->ibs_length)
    return 0;
  return a->hap < b->hap;
}

typedef struct
{
  st_phase_candidate_set *set;
} st_phase_candidate_context;

static int st_phase_hmm_error_valid(double error_rate)
{
  return isfinite(error_rate) && (error_rate >= 0.0) &&
         (error_rate <= 1.0);
}

static int st_phase_hmm_hap_emission(uint8_t observed, uint8_t donor,
                                     uint8_t n_alleles, double error_rate,
                                     double *probability)
{
  double p;

  if ((NULL == probability) || !st_phase_hmm_error_valid(error_rate) ||
      (0u == n_alleles) || (ST_PHASE_STATE_MISSING == n_alleles) ||
      (donor >= n_alleles))
    return -1;

  if (ST_PHASE_STATE_MISSING == observed)
    p = 1.0;
  else if (observed >= n_alleles)
    return -1;
  else if (1u == n_alleles)
    p = 1.0;
  else if (observed == donor)
    p = 1.0 - error_rate;
  else
    p = error_rate / (double)(n_alleles - 1u);

  *probability = p;
  return 0;
}

static int st_phase_hmm_break_probability(double distance_cm,
                                          double break_rate_per_morgan,
                                          double *probability)
{
  double x, p;

  if ((NULL == probability) || !isfinite(distance_cm) ||
      !isfinite(break_rate_per_morgan) || (distance_cm < 0.0) ||
      (break_rate_per_morgan < 0.0))
    return -1;

  x = (distance_cm / 100.0) * break_rate_per_morgan;

  if (!isfinite(x))
    p = 1.0;
  else
    p = -expm1(-x);

  if (p < 0.0)
    p = 0.0;
  if (p > 1.0)
    p = 1.0;

  *probability = p;
  return 0;
}

typedef struct
{
  size_t n_marker;
  size_t max_states;
  size_t n_slot;

  double *forward;

  double *log_scale;

  size_t *viterbi_hap;

  double log_likelihood;
  double viterbi_log_probability;
  int possible;
} st_phase_hmm_result;

static double st_phase_hmm_logadd(double a, double b)
{
  double hi, lo;

  if (-INFINITY == a)
    return b;
  if (-INFINITY == b)
    return a;

  if (a >= b)
  {
    hi = a;
    lo = b;
  }
  else
  {
    hi = b;
    lo = a;
  }

  return hi + log1p(exp(lo - hi));
}

static double st_phase_hmm_log_probability(double p)
{
  return (p > 0.0) ? log(p) : -INFINITY;
}

typedef struct
{
  size_t n_marker;
  size_t max_states;
  size_t n_slot;

  double *backward;

  double *posterior;

  double log_likelihood;
  int possible;
} st_phase_hmm_posterior_result;

static int st_phase_hmm_posterior_log_row(const double *log_forward,
                                          const double *log_backward,
                                          size_t n, double *posterior)
{
  double log_norm = -INFINITY;

  if ((NULL == log_forward) || (NULL == log_backward) ||
      (NULL == posterior) || (0u == n))
    return -1;

  for (size_t j = 0; j < n; j++)
    if ((log_forward[j] > -INFINITY) && (log_backward[j] > -INFINITY))
      log_norm =
          st_phase_hmm_logadd(log_norm, log_forward[j] + log_backward[j]);

  if (!(log_norm > -INFINITY))
    return -1;

  for (size_t j = 0; j < n; j++)
    posterior[j] =
        ((log_forward[j] > -INFINITY) && (log_backward[j] > -INFINITY))
            ? exp(log_forward[j] + log_backward[j] - log_norm)
            : 0.0;

  return 0;
}

typedef struct
{
  size_t hap1;
  size_t hap2;
  size_t slot1;
  size_t slot2;
} st_phase_pair_state;

typedef struct
{
  size_t n_marker;
  size_t n_hap;
  size_t target1;
  size_t target2;
  size_t max_side_states;
  size_t max_pair_states;
  size_t n_slot;

  size_t *left_count;
  size_t *right_count;
  size_t *count;

  st_phase_pair_state *state;
} st_phase_pair_state_set;

static void st_phase_pair_state_set_init(st_phase_pair_state_set *s)
{
  if (NULL != s)
  {
    memset(s, 0, sizeof(*s));
    s->target1 = SIZE_MAX;
    s->target2 = SIZE_MAX;
  }
}

static void st_phase_pair_state_set_free(st_phase_pair_state_set *s)
{
  if (NULL == s)
    return;

  free(s->left_count);
  free(s->right_count);
  free(s->count);
  free(s->state);
  st_phase_pair_state_set_init(s);
}

static const st_phase_pair_state *
st_phase_pair_state_get(const st_phase_pair_state_set *s, size_t marker,
                        size_t *count, size_t *left_count,
                        size_t *right_count)
{
  if (NULL != count)
    *count = 0;
  if (NULL != left_count)
    *left_count = 0;
  if (NULL != right_count)
    *right_count = 0;

  if ((NULL == s) || (marker >= s->n_marker) || (NULL == s->left_count) ||
      (NULL == s->right_count) || (NULL == s->count) ||
      (0u == s->max_pair_states) || (NULL == s->state))
    return NULL;

  if (NULL != count)
    *count = s->count[marker];
  if (NULL != left_count)
    *left_count = s->left_count[marker];
  if (NULL != right_count)
    *right_count = s->right_count[marker];

  return s->state + marker * s->max_pair_states;
}

static int st_phase_pair_state_row_valid(const st_phase_pair_state_set *s,
                                         size_t marker)
{
  size_t n, n1, n2, product, expected_pair, expected_slot;
  const st_phase_pair_state *v;

  if ((NULL == s) || (0u == s->n_marker) || (0u == s->n_hap) ||
      (marker >= s->n_marker) || (s->target1 >= s->n_hap) ||
      (s->target2 >= s->n_hap) || (0u == s->max_side_states) ||
      (NULL == s->left_count) || (NULL == s->right_count) ||
      (NULL == s->count) || (NULL == s->state) ||
      (0 != st_phase_size_mul(s->max_side_states, s->max_side_states,
                              &expected_pair)) ||
      (expected_pair != s->max_pair_states) ||
      (0 != st_phase_size_mul(s->n_marker, s->max_pair_states,
                              &expected_slot)) ||
      (expected_slot != s->n_slot))
    return 0;

  v = st_phase_pair_state_get(s, marker, &n, &n1, &n2);
  if ((NULL == v) || (0u == n1) || (0u == n2) ||
      (n1 > s->max_side_states) || (n2 > s->max_side_states) ||
      (0 != st_phase_size_mul(n1, n2, &product)) || (n != product) ||
      (n > s->max_pair_states))
    return 0;

  for (size_t i = 0; i < n1; i++)
    for (size_t j = 0; j < n2; j++)
    {
      size_t z = i * n2 + j;

      if ((v[z].slot1 != i) || (v[z].slot2 != j) ||
          (SIZE_MAX == v[z].hap1) || (SIZE_MAX == v[z].hap2) ||
          (v[z].hap1 >= s->n_hap) || (v[z].hap2 >= s->n_hap))
        return 0;

      if ((j > 0u) && (v[z].hap1 != v[i * n2].hap1))
        return 0;
      if ((i > 0u) && (v[z].hap2 != v[j].hap2))
        return 0;
    }

  for (size_t i = 0; i < n1; i++)
    for (size_t z = 0; z < i; z++)
      if (v[i * n2].hap1 == v[z * n2].hap1)
        return 0;

  for (size_t j = 0; j < n2; j++)
    for (size_t z = 0; z < j; z++)
      if (v[j].hap2 == v[z].hap2)
        return 0;

  return 1;
}

static int st_phase_hmm_pair_emission_components(
    uint8_t observed1, uint8_t observed2, uint8_t donor1, uint8_t donor2,
    uint8_t n_alleles, double error_rate, double *phase12, double *phase21,
    double *total)
{
  double a1, a2, b1 = 0.0, b2 = 0.0, p12, p21 = 0.0, sum;
  int m1 = (ST_PHASE_STATE_MISSING == observed1);
  int m2 = (ST_PHASE_STATE_MISSING == observed2);

  if ((NULL == phase12) || (NULL == phase21) || (NULL == total) ||
      !st_phase_hmm_error_valid(error_rate) || (0u == n_alleles) ||
      (ST_PHASE_STATE_MISSING == n_alleles) || (donor1 >= n_alleles) ||
      (donor2 >= n_alleles) || (m1 != m2))
    return -1;

  if (m1)
  {
    *phase12 = 1.0;
    *phase21 = 0.0;
    *total = 1.0;
    return 0;
  }

  if ((observed1 >= n_alleles) || (observed2 >= n_alleles))
    return -1;

  if ((0 != st_phase_hmm_hap_emission(observed1, donor1, n_alleles,
                                      error_rate, &a1)) ||
      (0 != st_phase_hmm_hap_emission(observed2, donor2, n_alleles,
                                      error_rate, &a2)))
    return -1;

  p12 = a1 * a2;

  if (observed1 != observed2)
  {
    if ((0 != st_phase_hmm_hap_emission(observed2, donor1, n_alleles,
                                        error_rate, &b1)) ||
        (0 != st_phase_hmm_hap_emission(observed1, donor2, n_alleles,
                                        error_rate, &b2)))
      return -1;
    p21 = b1 * b2;
  }

  sum = p12 + p21;
  *phase12 = p12;
  *phase21 = p21;
  *total = sum;
  return 0;
}

static int st_phase_hmm_pair_emission_row(
    const st_phase_pair_state_set *pair, size_t marker,
    const uint8_t *donor_allele, size_t n_hap, uint8_t observed1,
    uint8_t observed2, uint8_t n_alleles, double error_rate,
    double *phase12, double *phase21, double *total)
{
  size_t n, n1, n2;
  const st_phase_pair_state *row;
  int m1 = (ST_PHASE_STATE_MISSING == observed1);
  int m2 = (ST_PHASE_STATE_MISSING == observed2);

  if ((NULL == donor_allele) || (NULL == phase12) || (NULL == phase21) ||
      (NULL == total) || (n_hap != ((NULL != pair) ? pair->n_hap : 0u)) ||
      !st_phase_hmm_error_valid(error_rate) || (0u == n_alleles) ||
      (ST_PHASE_STATE_MISSING == n_alleles) || (m1 != m2) ||
      (!m1 && ((observed1 >= n_alleles) || (observed2 >= n_alleles))) ||
      !st_phase_pair_state_row_valid(pair, marker))
    return -1;

  row = st_phase_pair_state_get(pair, marker, &n, &n1, &n2);
  (void)n1;
  (void)n2;

  for (size_t z = 0; z < n; z++)
    if ((row[z].hap1 >= n_hap) || (row[z].hap2 >= n_hap) ||
        (donor_allele[row[z].hap1] >= n_alleles) ||
        (donor_allele[row[z].hap2] >= n_alleles))
      return -1;

  for (size_t z = 0; z < n; z++)
  {
    double p12, p21, p;

    if (0 != st_phase_hmm_pair_emission_components(
                 observed1, observed2, donor_allele[row[z].hap1],
                 donor_allele[row[z].hap2], n_alleles, error_rate, &p12,
                 &p21, &p))
      return -1;

    phase12[z] = p12;
    phase21[z] = p21;
    total[z] = p;
  }

  return 0;
}

typedef struct
{
  size_t n_marker;
  size_t max_pair_states;
  size_t n_slot;

  double *forward;
  double *backward;
  double *posterior;

  double *log_scale;

  size_t *viterbi_hap1;
  size_t *viterbi_hap2;

  double log_likelihood;
  double viterbi_log_probability;
  int possible;
} st_phase_pair_hmm_result;

static void st_phase_pair_hmm_result_init(st_phase_pair_hmm_result *r)
{
  if (NULL != r)
  {
    memset(r, 0, sizeof(*r));
    r->log_likelihood = -INFINITY;
    r->viterbi_log_probability = -INFINITY;
  }
}

static void st_phase_pair_hmm_result_free(st_phase_pair_hmm_result *r)
{
  if (NULL == r)
    return;

  free(r->forward);
  free(r->backward);
  free(r->posterior);
  free(r->log_scale);
  free(r->viterbi_hap1);
  free(r->viterbi_hap2);
  st_phase_pair_hmm_result_init(r);
}

static int st_phase_pair_hmm_set_valid(const st_phase_pair_state_set *pair)
{
  size_t expected_pair, expected_slot;

  if ((NULL == pair) || (0u == pair->n_marker) || (0u == pair->n_hap) ||
      (0u == pair->max_side_states) || (0u == pair->max_pair_states) ||
      (NULL == pair->left_count) || (NULL == pair->right_count) ||
      (NULL == pair->count) || (NULL == pair->state) ||
      (0 != st_phase_size_mul(pair->max_side_states, pair->max_side_states,
                              &expected_pair)) ||
      (expected_pair != pair->max_pair_states) ||
      (0 != st_phase_size_mul(pair->n_marker, pair->max_pair_states,
                              &expected_slot)) ||
      (expected_slot != pair->n_slot))
    return 0;

  for (size_t marker = 0; marker < pair->n_marker; marker++)
    if (!st_phase_pair_state_row_valid(pair, marker))
      return 0;

  return 1;
}

static void st_phase_pair_hmm_compensated_add(double value, double *sum,
                                              double *compensation)
{
  double y, t;

  y = value - *compensation;
  t = *sum + y;
  *compensation = (t - *sum) - y;
  *sum = t;
}

typedef enum
{
  ST_PHASE_ORIENTATION_MISSING = 0,
  ST_PHASE_ORIENTATION_HOMOZYGOUS = 1,
  ST_PHASE_ORIENTATION_HETEROZYGOUS = 2
} st_phase_orientation_kind;

typedef struct
{
  size_t n_marker;
  st_phase_orientation_kind *kind;
  double *phase12;
  double *phase21;
} st_phase_orientation_result;

static void
st_phase_orientation_result_init(st_phase_orientation_result *r)
{
  if (NULL != r)
    memset(r, 0, sizeof(*r));
}

static void
st_phase_orientation_result_free(st_phase_orientation_result *r)
{
  if (NULL == r)
    return;

  free(r->kind);
  free(r->phase12);
  free(r->phase21);
  st_phase_orientation_result_init(r);
}

static int st_phase_orientation_probability_valid(double x)
{
  return isfinite(x) && (x >= 0.0) && (x <= 1.0);
}

static int st_phase_orientation_posterior(
    const st_phase_pair_state_set *pair,
    const st_phase_pair_hmm_result *hmm, const uint8_t *observed1,
    const uint8_t *observed2, const double *phase12_component,
    const double *phase21_component, st_phase_orientation_result *out)
{
  st_phase_orientation_result tmp;
  size_t M, P;

  st_phase_orientation_result_init(&tmp);

  if ((NULL == out) || (NULL == observed1) || (NULL == observed2) ||
      (NULL == phase12_component) || (NULL == phase21_component) ||
      !st_phase_pair_hmm_set_valid(pair) || (NULL == hmm) ||
      (hmm->n_marker != pair->n_marker) ||
      (hmm->max_pair_states != pair->max_pair_states) ||
      (hmm->n_slot != pair->n_slot) || (NULL == hmm->posterior) ||
      ((0 != hmm->possible) && (1 != hmm->possible)))
    return -1;

  M = pair->n_marker;
  P = pair->max_pair_states;

  tmp.kind = malloc(M * sizeof(*tmp.kind));
  tmp.phase12 = calloc(M, sizeof(*tmp.phase12));
  tmp.phase21 = calloc(M, sizeof(*tmp.phase21));
  if ((NULL == tmp.kind) || (NULL == tmp.phase12) || (NULL == tmp.phase21))
    goto fail;
  tmp.n_marker = M;

  for (size_t marker = 0; marker < M; marker++)
  {
    size_t n;
    const st_phase_pair_state *row =
        st_phase_pair_state_get(pair, marker, &n, NULL, NULL);
    int m1 = (ST_PHASE_STATE_MISSING == observed1[marker]);
    int m2 = (ST_PHASE_STATE_MISSING == observed2[marker]);
    double post_sum = 0.0, post_comp = 0.0;

    if ((NULL == row) || (0u == n) || (n > P) || (m1 != m2) ||
        (!m1 && (observed2[marker] < observed1[marker])))
      goto fail;
    (void)row;

    if (m1)
      tmp.kind[marker] = ST_PHASE_ORIENTATION_MISSING;
    else if (observed1[marker] == observed2[marker])
      tmp.kind[marker] = ST_PHASE_ORIENTATION_HOMOZYGOUS;
    else
      tmp.kind[marker] = ST_PHASE_ORIENTATION_HETEROZYGOUS;

    for (size_t z = 0; z < n; z++)
    {
      size_t slot = marker * P + z;
      double post = hmm->posterior[slot];
      double e12 = phase12_component[slot];
      double e21 = phase21_component[slot];

      if (!st_phase_orientation_probability_valid(post) ||
          !st_phase_orientation_probability_valid(e12) ||
          !st_phase_orientation_probability_valid(e21))
        goto fail;

      st_phase_pair_hmm_compensated_add(post, &post_sum, &post_comp);

      if (ST_PHASE_ORIENTATION_HETEROZYGOUS == tmp.kind[marker])
      {
        double denom = e12 + e21;

        if (denom > 1.0 + 16.0 * DBL_EPSILON)
          goto fail;

        if (!(denom > 0.0) && (post > 0.0))
          goto fail;
      }
      else if (ST_PHASE_ORIENTATION_HOMOZYGOUS == tmp.kind[marker])
      {
        if (e21 != 0.0)
          goto fail;
      }
    }

    if (ST_PHASE_ORIENTATION_HETEROZYGOUS == tmp.kind[marker])
    {
      double s12 = 0.0, c12 = 0.0, s21 = 0.0, c21 = 0.0;

      for (size_t z = 0; z < n; z++)
      {
        size_t slot = marker * P + z;
        double post = hmm->posterior[slot];
        double e12 = phase12_component[slot];
        double e21 = phase21_component[slot];
        double denom = e12 + e21;

        if (denom > 0.0)
        {
          st_phase_pair_hmm_compensated_add(post * (e12 / denom), &s12,
                                            &c12);
          st_phase_pair_hmm_compensated_add(post * (e21 / denom), &s21,
                                            &c21);
        }
      }

      if (hmm->possible)
      {
        double support_sum = s12 + s21;
        double tol = 64.0 * DBL_EPSILON * (double)((n > 1u) ? n : 1u);

        if ((fabs(post_sum - 1.0) > tol) || !(support_sum > 0.0) ||
            (fabs(support_sum - 1.0) > tol))
          goto fail;

        tmp.phase12[marker] = s12 / support_sum;
        tmp.phase21[marker] = s21 / support_sum;
      }
      else
      {
        if (post_sum != 0.0)
          goto fail;
        tmp.phase12[marker] = 0.0;
        tmp.phase21[marker] = 0.0;
      }
    }
    else if (ST_PHASE_ORIENTATION_HOMOZYGOUS == tmp.kind[marker])
    {
      if (hmm->possible)
      {
        double tol = 64.0 * DBL_EPSILON * (double)((n > 1u) ? n : 1u);
        if (fabs(post_sum - 1.0) > tol)
          goto fail;
        tmp.phase12[marker] = 1.0;
        tmp.phase21[marker] = 0.0;
      }
      else
      {
        if (post_sum != 0.0)
          goto fail;
        tmp.phase12[marker] = 0.0;
        tmp.phase21[marker] = 0.0;
      }
    }
    else
    {

      if (hmm->possible)
      {
        double tol = 64.0 * DBL_EPSILON * (double)((n > 1u) ? n : 1u);
        if (fabs(post_sum - 1.0) > tol)
          goto fail;
      }
      else if (post_sum != 0.0)
        goto fail;
      tmp.phase12[marker] = 0.0;
      tmp.phase21[marker] = 0.0;
    }

    for (size_t z = n; z < P; z++)
    {
      size_t slot = marker * P + z;
      if ((hmm->posterior[slot] != 0.0) ||
          (phase12_component[slot] != 0.0) ||
          (phase21_component[slot] != 0.0))
        goto fail;
    }
  }

  st_phase_orientation_result_free(out);
  *out = tmp;
  return 0;

fail:
  st_phase_orientation_result_free(&tmp);
  return -1;
}

typedef enum
{
  ST_PHASE_HARD_MISSING = 0,
  ST_PHASE_HARD_HOMOZYGOUS = 1,
  ST_PHASE_HARD_12 = 2,
  ST_PHASE_HARD_21 = 3,
  ST_PHASE_HARD_UNRESOLVED = 4
} st_phase_hard_call;

typedef struct
{
  size_t n_marker;
  st_phase_hard_call *call;
  double *confidence;
} st_phase_hard_call_result;

static void st_phase_hard_call_result_init(st_phase_hard_call_result *r)
{
  if (NULL != r)
    memset(r, 0, sizeof(*r));
}

static void st_phase_hard_call_result_free(st_phase_hard_call_result *r)
{
  if (NULL == r)
    return;

  free(r->call);
  free(r->confidence);
  st_phase_hard_call_result_init(r);
}

static int st_phase_hard_threshold_valid(double threshold)
{
  return isfinite(threshold) && (threshold >= 0.5) && (threshold <= 1.0);
}

static int st_phase_step0_probability_close(double x, double y,
                                            double scale);

static int st_phase_hard_call_from_orientation(
    const st_phase_orientation_result *orientation, double threshold,
    st_phase_hard_call_result *out)
{
  st_phase_hard_call_result tmp;
  size_t M;

  st_phase_hard_call_result_init(&tmp);

  if ((NULL == out) || (NULL == orientation) ||
      (0u == orientation->n_marker) || (NULL == orientation->kind) ||
      (NULL == orientation->phase12) || (NULL == orientation->phase21) ||
      !st_phase_hard_threshold_valid(threshold))
    return -1;

  M = orientation->n_marker;
  tmp.call = malloc(M * sizeof(*tmp.call));
  tmp.confidence = malloc(M * sizeof(*tmp.confidence));
  if ((NULL == tmp.call) || (NULL == tmp.confidence))
    goto fail;
  tmp.n_marker = M;

  for (size_t marker = 0; marker < M; marker++)
  {
    st_phase_orientation_kind kind = orientation->kind[marker];
    double p12 = orientation->phase12[marker];
    double p21 = orientation->phase21[marker];

    if (!st_phase_orientation_probability_valid(p12) ||
        !st_phase_orientation_probability_valid(p21))
      goto fail;

    switch (kind)
    {
    case ST_PHASE_ORIENTATION_MISSING:
      if ((0.0 != p12) || (0.0 != p21))
        goto fail;
      tmp.call[marker] = ST_PHASE_HARD_MISSING;
      tmp.confidence[marker] = 0.0;
      break;

    case ST_PHASE_ORIENTATION_HOMOZYGOUS:
      if ((1.0 == p12) && (0.0 == p21))
      {
        tmp.call[marker] = ST_PHASE_HARD_HOMOZYGOUS;
        tmp.confidence[marker] = 1.0;
      }
      else if ((0.0 == p12) && (0.0 == p21))
      {
        tmp.call[marker] = ST_PHASE_HARD_UNRESOLVED;
        tmp.confidence[marker] = 0.0;
      }
      else
        goto fail;
      break;

    case ST_PHASE_ORIENTATION_HETEROZYGOUS:
      if ((0.0 == p12) && (0.0 == p21))
      {
        tmp.call[marker] = ST_PHASE_HARD_UNRESOLVED;
        tmp.confidence[marker] = 0.0;
      }
      else
      {
        double sum = p12 + p21;
        double tol = 64.0 * DBL_EPSILON;
        double q12, q21, best;

        if (!(sum > 0.0) || (fabs(sum - 1.0) > tol))
          goto fail;

        q12 = p12 / sum;
        q21 = p21 / sum;
        best = (q12 > q21) ? q12 : q21;
        tmp.confidence[marker] = best;

        if ((q12 > 0.5) && (q12 >= threshold) &&
            !st_phase_step0_probability_close(q12, q21, 1.0))
          tmp.call[marker] = ST_PHASE_HARD_12;
        else if ((q21 > 0.5) && (q21 >= threshold) &&
                 !st_phase_step0_probability_close(q12, q21, 1.0))
          tmp.call[marker] = ST_PHASE_HARD_21;
        else
          tmp.call[marker] = ST_PHASE_HARD_UNRESOLVED;
      }
      break;

    default:
      goto fail;
    }
  }

  st_phase_hard_call_result_free(out);
  *out = tmp;
  return 0;

fail:
  st_phase_hard_call_result_free(&tmp);
  return -1;
}

typedef struct
{
  size_t individual;
  size_t chromosome;
  int chrom_label;
  size_t marker_start;
  size_t n_marker;
  uint8_t *observed1;
  uint8_t *observed2;
} st_phase_diploid_input;

typedef struct
{
  size_t individual;
  size_t chromosome;
  int chrom_label;
  size_t marker_start;
  size_t n_marker;
  gs_UNSIGNED *source_marker;
  GS_ALLELE_TYPE *allele1;
  GS_ALLELE_TYPE *allele2;
  st_phase_hard_call *call;
  double *confidence;
  unsigned char *write_phase;
} st_phase_proposed_output;

#define st_phase_adapter_malloc malloc

static void st_phase_diploid_input_init(st_phase_diploid_input *x)
{
  if (NULL != x)
    memset(x, 0, sizeof(*x));
}

static void st_phase_diploid_input_free(st_phase_diploid_input *x)
{
  if (NULL == x)
    return;

  free(x->observed1);
  free(x->observed2);
  st_phase_diploid_input_init(x);
}

static void st_phase_proposed_output_init(st_phase_proposed_output *x)
{
  if (NULL != x)
    memset(x, 0, sizeof(*x));
}

static void st_phase_proposed_output_free(st_phase_proposed_output *x)
{
  if (NULL == x)
    return;

  free(x->source_marker);
  free(x->allele1);
  free(x->allele2);
  free(x->call);
  free(x->confidence);
  free(x->write_phase);
  st_phase_proposed_output_init(x);
}

static int st_phase_data_adapter_valid(const st_phase_data *d)
{
  size_t n_genotype;

  if ((NULL == d) || (0u == d->n_ind) || (0u == d->n_marker) ||
      (0u == d->n_chrom) || (NULL == d->source_marker) ||
      (NULL == d->chrom_label) || (NULL == d->pos_cm) ||
      (NULL == d->chrom_start) || (NULL == d->chrom_length) ||
      (NULL == d->allele_offset) || (NULL == d->n_alleles) ||
      (NULL == d->geno1) || (NULL == d->geno2) ||
      ((d->n_allele_value > 0u) && (NULL == d->allele_value)) ||
      (0 != st_phase_size_mul(d->n_ind, d->n_marker, &n_genotype)))
    return 0;

  (void)n_genotype;

  if ((0u != d->allele_offset[0]) ||
      (d->allele_offset[d->n_marker] != d->n_allele_value) ||
      (0u != d->chrom_start[0]))
    return 0;

  {
    size_t expected_start = 0u;

    for (size_t c = 0; c < d->n_chrom; c++)
    {
      size_t start = d->chrom_start[c];
      size_t len = d->chrom_length[c];

      if ((start != expected_start) || (0u == len) ||
          (start > d->n_marker) || (len > d->n_marker - start))
        return 0;

      for (size_t m = start; m < start + len; m++)
      {
        if (d->chrom_label[m] != d->chrom_label[start])
          return 0;
        if ((m > start) && (d->pos_cm[m] < d->pos_cm[m - 1u]))
          return 0;
        if (!isfinite(d->pos_cm[m]) || (d->pos_cm[m] < 0.0))
          return 0;
      }

      if ((c > 0u) && (d->chrom_label[start] <=
                       d->chrom_label[d->chrom_start[c - 1u]]))
        return 0;

      expected_start = start + len;
    }

    if (expected_start != d->n_marker)
      return 0;
  }

  for (size_t m = 0; m < d->n_marker; m++)
  {
    size_t begin = d->allele_offset[m];
    size_t end = d->allele_offset[m + 1u];
    size_t n = (size_t)d->n_alleles[m];

    if ((end < begin) || (end > d->n_allele_value) || (end - begin != n) ||
        (n > 254u))
      return 0;

    for (size_t s = 0; s < n; s++)
    {
      GS_ALLELE_TYPE value = d->allele_value[begin + s];
      if ((0 == value) || (GS_ALLELE_MISSING == value) ||
          ((s > 0u) && !(d->allele_value[begin + s - 1u] < value)))
        return 0;
    }

    for (size_t i = 0; i < d->n_ind; i++)
    {
      size_t slot = m * d->n_ind + i;
      uint8_t a = d->geno1[slot];
      uint8_t b = d->geno2[slot];
      int ma = (ST_PHASE_STATE_MISSING == a);
      int mb = (ST_PHASE_STATE_MISSING == b);

      if (ma != mb)
        return 0;
      if (!ma &&
          ((a >= d->n_alleles[m]) || (b >= d->n_alleles[m]) || (a > b)))
        return 0;
    }
  }

  return 1;
}

static int st_phase_prepare_diploid_input(const st_phase_data *d,
                                          size_t individual,
                                          size_t chromosome,
                                          st_phase_diploid_input *out)
{
  st_phase_diploid_input tmp;
  size_t start, len;

  st_phase_diploid_input_init(&tmp);

  if ((NULL == out) || !st_phase_data_adapter_valid(d) ||
      (individual >= d->n_ind) || (chromosome >= d->n_chrom))
    return -1;

  start = d->chrom_start[chromosome];
  len = d->chrom_length[chromosome];

  tmp.observed1 = st_phase_adapter_malloc(len * sizeof(*tmp.observed1));
  tmp.observed2 = st_phase_adapter_malloc(len * sizeof(*tmp.observed2));
  if ((NULL == tmp.observed1) || (NULL == tmp.observed2))
    goto fail;

  tmp.individual = individual;
  tmp.chromosome = chromosome;
  tmp.chrom_label = d->chrom_label[start];
  tmp.marker_start = start;
  tmp.n_marker = len;

  for (size_t j = 0; j < len; j++)
  {
    size_t slot = (start + j) * d->n_ind + individual;
    tmp.observed1[j] = d->geno1[slot];
    tmp.observed2[j] = d->geno2[slot];
  }

  st_phase_diploid_input_free(out);
  *out = tmp;
  return 0;

fail:
  st_phase_diploid_input_free(&tmp);
  return -1;
}

static int
st_phase_diploid_input_matches(const st_phase_data *d,
                               const st_phase_diploid_input *input)
{
  size_t start, len;

  if ((NULL == input) || !st_phase_data_adapter_valid(d) ||
      (NULL == input->observed1) || (NULL == input->observed2) ||
      (input->individual >= d->n_ind) || (input->chromosome >= d->n_chrom))
    return 0;

  start = d->chrom_start[input->chromosome];
  len = d->chrom_length[input->chromosome];

  if ((input->marker_start != start) || (input->n_marker != len) ||
      (input->chrom_label != d->chrom_label[start]))
    return 0;

  for (size_t j = 0; j < len; j++)
  {
    size_t slot = (start + j) * d->n_ind + input->individual;
    if ((input->observed1[j] != d->geno1[slot]) ||
        (input->observed2[j] != d->geno2[slot]))
      return 0;
  }

  return 1;
}

static int st_phase_build_proposed_output(
    const st_phase_data *d, const st_phase_diploid_input *input,
    const st_phase_hard_call_result *hard, st_phase_proposed_output *out)
{
  st_phase_proposed_output tmp;
  size_t M;

  st_phase_proposed_output_init(&tmp);

  if ((NULL == out) || (NULL == hard) ||
      !st_phase_diploid_input_matches(d, input) ||
      (hard->n_marker != input->n_marker) || (0u == hard->n_marker) ||
      (NULL == hard->call) || (NULL == hard->confidence))
    return -1;

  M = input->n_marker;
  tmp.source_marker =
      st_phase_adapter_malloc(M * sizeof(*tmp.source_marker));
  tmp.allele1 = st_phase_adapter_malloc(M * sizeof(*tmp.allele1));
  tmp.allele2 = st_phase_adapter_malloc(M * sizeof(*tmp.allele2));
  tmp.call = st_phase_adapter_malloc(M * sizeof(*tmp.call));
  tmp.confidence = st_phase_adapter_malloc(M * sizeof(*tmp.confidence));
  tmp.write_phase = st_phase_adapter_malloc(M * sizeof(*tmp.write_phase));
  if ((NULL == tmp.source_marker) || (NULL == tmp.allele1) ||
      (NULL == tmp.allele2) || (NULL == tmp.call) ||
      (NULL == tmp.confidence) || (NULL == tmp.write_phase))
    goto fail;

  tmp.individual = input->individual;
  tmp.chromosome = input->chromosome;
  tmp.chrom_label = input->chrom_label;
  tmp.marker_start = input->marker_start;
  tmp.n_marker = M;

  for (size_t j = 0; j < M; j++)
  {
    size_t marker = input->marker_start + j;
    size_t offset = d->allele_offset[marker];
    uint8_t x = input->observed1[j];
    uint8_t y = input->observed2[j];
    st_phase_hard_call call = hard->call[j];
    double confidence = hard->confidence[j];
    int missing = (ST_PHASE_STATE_MISSING == x);
    int homo = (!missing && (x == y));

    if (!isfinite(confidence) || (confidence < 0.0) ||
        (confidence > 1.0) || (missing != (ST_PHASE_STATE_MISSING == y)))
      goto fail;

    tmp.source_marker[j] = d->source_marker[marker];
    tmp.call[j] = call;
    tmp.confidence[j] = confidence;
    tmp.write_phase[j] = 0u;

    if (missing)
    {
      if ((ST_PHASE_HARD_MISSING != call) || (0.0 != confidence))
        goto fail;
      tmp.allele1[j] = GS_ALLELE_MISSING;
      tmp.allele2[j] = GS_ALLELE_MISSING;
      continue;
    }

    if ((x >= d->n_alleles[marker]) || (y >= d->n_alleles[marker]))
      goto fail;

    {
      GS_ALLELE_TYPE a = d->allele_value[offset + (size_t)x];
      GS_ALLELE_TYPE b = d->allele_value[offset + (size_t)y];

      if (homo)
      {
        if (ST_PHASE_HARD_HOMOZYGOUS == call)
        {
          if (1.0 != confidence)
            goto fail;
        }
        else if (ST_PHASE_HARD_UNRESOLVED == call)
        {
          if (0.0 != confidence)
            goto fail;
        }
        else
          goto fail;

        tmp.allele1[j] = a;
        tmp.allele2[j] = a;
      }
      else
      {
        switch (call)
        {
        case ST_PHASE_HARD_12:
          if (!(confidence >= 0.5))
            goto fail;
          tmp.allele1[j] = a;
          tmp.allele2[j] = b;
          tmp.write_phase[j] = 1u;
          break;

        case ST_PHASE_HARD_21:
          if (!(confidence >= 0.5))
            goto fail;
          tmp.allele1[j] = b;
          tmp.allele2[j] = a;
          tmp.write_phase[j] = 1u;
          break;

        case ST_PHASE_HARD_UNRESOLVED:
          if (!((0.0 == confidence) || (confidence >= 0.5)))
            goto fail;
          tmp.allele1[j] = a;
          tmp.allele2[j] = b;
          break;

        default:
          goto fail;
        }
      }
    }
  }

  st_phase_proposed_output_free(out);
  *out = tmp;
  return 0;

fail:
  st_phase_proposed_output_free(&tmp);
  return -1;
}

typedef enum
{
  ST_PHASE_STAGE_OK = 0,
  ST_PHASE_STAGE_INVALID = -1,
  ST_PHASE_STAGE_STALE = -2,
  ST_PHASE_STAGE_OVERLAP = -3,
  ST_PHASE_STAGE_CONFLICT = -4
} st_phase_stage_status;

typedef struct
{
  size_t n_ind;
  size_t n_marker;
  size_t n_genotype;
  size_t n_written;
  size_t n_covered;

  GS_ALLELE_TYPE *base1;
  GS_ALLELE_TYPE *base2;

  GS_ALLELE_TYPE *staged1;
  GS_ALLELE_TYPE *staged2;

  unsigned char *claimed;

  unsigned char *covered;
} st_phase_dataset_stage;

#define st_phase_stage_malloc malloc

static void st_phase_dataset_stage_init(st_phase_dataset_stage *x)
{
  if (NULL != x)
    memset(x, 0, sizeof(*x));
}

static void st_phase_dataset_stage_free(st_phase_dataset_stage *x)
{
  if (NULL == x)
    return;

  free(x->base1);
  free(x->base2);
  free(x->staged1);
  free(x->staged2);
  free(x->claimed);
  free(x->covered);
  st_phase_dataset_stage_init(x);
}

static int st_phase_dataset_stage_structurally_valid(
    const st_phase_dataset_stage *stage)
{
  size_t n;

  if ((NULL == stage) || (0u == stage->n_ind) || (0u == stage->n_marker) ||
      (0u == stage->n_genotype) || (NULL == stage->base1) ||
      (NULL == stage->base2) || (NULL == stage->staged1) ||
      (NULL == stage->staged2) || (NULL == stage->claimed) ||
      (NULL == stage->covered) ||
      (0 != st_phase_size_mul(stage->n_ind, stage->n_marker, &n)) ||
      (n != stage->n_genotype) || (stage->n_written > n) ||
      (stage->n_covered > n))
    return 0;

  {
    size_t count = 0u, covered_count = 0u;
    for (size_t k = 0; k < n; k++)
    {
      if ((stage->claimed[k] > 1u) || (stage->covered[k] > 1u) ||
          (stage->claimed[k] && !stage->covered[k]))
        return 0;
      if (stage->covered[k])
        covered_count++;
      if (stage->claimed[k])
      {
        count++;
        if ((stage->base1[k] == stage->base2[k]) ||
            (GS_ALLELE_MISSING == stage->base1[k]) ||
            (GS_ALLELE_MISSING == stage->base2[k]))
          return 0;
        if (!(((stage->staged1[k] == stage->base1[k]) &&
               (stage->staged2[k] == stage->base2[k])) ||
              ((stage->staged1[k] == stage->base2[k]) &&
               (stage->staged2[k] == stage->base1[k]))))
          return 0;
      }
      else if ((stage->staged1[k] != stage->base1[k]) ||
               (stage->staged2[k] != stage->base2[k]))
        return 0;
    }
    if ((count != stage->n_written) || (covered_count != stage->n_covered))
      return 0;
  }

  return 1;
}

static int
st_phase_stage_workspace_matches_source(const gs_varset_type *gv,
                                        const st_phase_data *d)
{
  size_t n;

  if ((NULL == gv) || !st_phase_data_adapter_valid(d) ||
      (NULL == gv->All1) || (NULL == gv->All2) || (NULL == gv->map) ||
      (NULL == gv->mdl) || ((size_t)gv->NoInd != d->n_ind) ||
      ((size_t)gv->NoMar != d->n_marker) ||
      ((size_t)gv->NLoc != d->n_marker) ||
      (0 != st_phase_size_mul(d->n_ind, d->n_marker, &n)))
    return 0;

  (void)n;

  for (size_t m = 0; m < d->n_marker; m++)
  {
    gs_UNSIGNED source = d->source_marker[m];
    size_t offset = d->allele_offset[m];

    if ((source >= gv->NoMar) || (gv->mdl[m] != source) ||
        (NULL == gv->map[m]) || (gv->map[m]->nr != source) ||
        (gv->map[m]->chrom != d->chrom_label[m]) ||
        (gv->map[m]->pos != d->pos_cm[m]) || (d->n_alleles[m] > 254u))
      return 0;

    for (size_t i = 0; i < d->n_ind; i++)
    {
      size_t raw = (size_t)source * d->n_ind + i;
      size_t compact = m * d->n_ind + i;
      GS_ALLELE_TYPE a = gv->All1[raw];
      GS_ALLELE_TYPE b = gv->All2[raw];
      uint8_t x = d->geno1[compact];
      uint8_t y = d->geno2[compact];

      if (ST_PHASE_STATE_MISSING == x)
      {
        if ((ST_PHASE_STATE_MISSING != y) || (GS_ALLELE_MISSING != a) ||
            (GS_ALLELE_MISSING != b))
          return 0;
      }
      else
      {
        GS_ALLELE_TYPE dx, dy;

        if ((ST_PHASE_STATE_MISSING == y) || (x >= d->n_alleles[m]) ||
            (y >= d->n_alleles[m]))
          return 0;
        dx = d->allele_value[offset + (size_t)x];
        dy = d->allele_value[offset + (size_t)y];
        if (!(((a == dx) && (b == dy)) || ((a == dy) && (b == dx))))
          return 0;
      }
    }
  }

  return 1;
}

static int st_phase_dataset_stage_build(const gs_varset_type *gv,
                                        const st_phase_data *d,
                                        st_phase_dataset_stage *out)
{
  st_phase_dataset_stage tmp;
  size_t n, bytes;

  st_phase_dataset_stage_init(&tmp);

  if ((NULL == out) || !st_phase_stage_workspace_matches_source(gv, d) ||
      (0 != st_phase_size_mul(d->n_ind, d->n_marker, &n)) ||
      (0 != st_phase_size_mul(n, sizeof(GS_ALLELE_TYPE), &bytes)))
    return ST_PHASE_STAGE_INVALID;

  tmp.base1 = st_phase_stage_malloc(bytes);
  tmp.base2 = st_phase_stage_malloc(bytes);
  tmp.staged1 = st_phase_stage_malloc(bytes);
  tmp.staged2 = st_phase_stage_malloc(bytes);
  tmp.claimed = st_phase_stage_malloc(n * sizeof(*tmp.claimed));
  tmp.covered = st_phase_stage_malloc(n * sizeof(*tmp.covered));
  if ((NULL == tmp.base1) || (NULL == tmp.base2) ||
      (NULL == tmp.staged1) || (NULL == tmp.staged2) ||
      (NULL == tmp.claimed) || (NULL == tmp.covered))
    goto fail;

  tmp.n_ind = d->n_ind;
  tmp.n_marker = d->n_marker;
  tmp.n_genotype = n;
  tmp.n_written = 0u;
  tmp.n_covered = 0u;
  memcpy(tmp.base1, gv->All1, bytes);
  memcpy(tmp.base2, gv->All2, bytes);
  memcpy(tmp.staged1, gv->All1, bytes);
  memcpy(tmp.staged2, gv->All2, bytes);
  memset(tmp.claimed, 0, n * sizeof(*tmp.claimed));
  memset(tmp.covered, 0, n * sizeof(*tmp.covered));

  if (!st_phase_dataset_stage_structurally_valid(&tmp))
    goto fail;

  st_phase_dataset_stage_free(out);
  *out = tmp;
  return ST_PHASE_STAGE_OK;

fail:
  st_phase_dataset_stage_free(&tmp);
  return ST_PHASE_STAGE_INVALID;
}

static int
st_phase_dataset_stage_source_current(const gs_varset_type *gv,
                                      const st_phase_dataset_stage *stage)
{
  size_t bytes;

  if ((NULL == gv) || !st_phase_dataset_stage_structurally_valid(stage) ||
      (NULL == gv->All1) || (NULL == gv->All2) ||
      ((size_t)gv->NoInd != stage->n_ind) ||
      ((size_t)gv->NoMar != stage->n_marker) ||
      (0 != st_phase_size_mul(stage->n_genotype, sizeof(GS_ALLELE_TYPE),
                              &bytes)))
    return 0;

  return ((0 == memcmp(stage->base1, gv->All1, bytes)) &&
          (0 == memcmp(stage->base2, gv->All2, bytes)));
}

static int
st_phase_stage_proposal_valid(const st_phase_data *d,
                              const st_phase_dataset_stage *stage,
                              const st_phase_proposed_output *p)
{
  size_t start, len;

  if (!st_phase_data_adapter_valid(d) ||
      !st_phase_dataset_stage_structurally_valid(stage) || (NULL == p) ||
      (0u == p->n_marker) || (NULL == p->source_marker) ||
      (NULL == p->allele1) || (NULL == p->allele2) || (NULL == p->call) ||
      (NULL == p->confidence) || (NULL == p->write_phase) ||
      (p->individual >= d->n_ind) || (p->individual >= stage->n_ind) ||
      (p->chromosome >= d->n_chrom) || (d->n_ind != stage->n_ind) ||
      (d->n_marker != stage->n_marker))
    return 0;

  start = d->chrom_start[p->chromosome];
  len = d->chrom_length[p->chromosome];
  if ((p->marker_start != start) || (p->n_marker != len) ||
      (p->chrom_label != d->chrom_label[start]))
    return 0;

  for (size_t j = 0; j < len; j++)
  {
    size_t marker = start + j;
    gs_UNSIGNED source = d->source_marker[marker];
    size_t slot;
    GS_ALLELE_TYPE a, b, x, y, lo, hi;
    st_phase_hard_call call = p->call[j];
    double conf = p->confidence[j];
    unsigned char write = p->write_phase[j];

    if ((p->source_marker[j] != source) || (source >= stage->n_marker) ||
        (write > 1u) || !isfinite(conf) || (conf < 0.0) || (conf > 1.0))
      return 0;

    slot = (size_t)source * stage->n_ind + p->individual;
    a = stage->base1[slot];
    b = stage->base2[slot];
    x = p->allele1[j];
    y = p->allele2[j];

    if ((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b))
    {
      if (!((GS_ALLELE_MISSING == a) && (GS_ALLELE_MISSING == b) &&
            (GS_ALLELE_MISSING == x) && (GS_ALLELE_MISSING == y) &&
            (ST_PHASE_HARD_MISSING == call) && (0.0 == conf) &&
            (0u == write)))
        return 0;
      continue;
    }

    if ((0 == a) || (0 == b) ||
        !(((a == x) && (b == y)) || ((a == y) && (b == x))))
      return 0;

    if (a == b)
    {
      if ((x != a) || (y != a) || (0u != write))
        return 0;
      if (ST_PHASE_HARD_HOMOZYGOUS == call)
      {
        if (1.0 != conf)
          return 0;
      }
      else if (ST_PHASE_HARD_UNRESOLVED == call)
      {
        if (0.0 != conf)
          return 0;
      }
      else
        return 0;
      continue;
    }

    lo = (a < b) ? a : b;
    hi = (a < b) ? b : a;
    switch (call)
    {
    case ST_PHASE_HARD_12:
      if ((1u != write) || !(conf >= 0.5) || (x != lo) || (y != hi))
        return 0;
      break;

    case ST_PHASE_HARD_21:
      if ((1u != write) || !(conf >= 0.5) || (x != hi) || (y != lo))
        return 0;
      break;

    case ST_PHASE_HARD_UNRESOLVED:
      if ((0u != write) || !((0.0 == conf) || (conf >= 0.5)) ||
          (x != lo) || (y != hi))
        return 0;
      break;

    default:
      return 0;
    }
  }

  return 1;
}

static int st_phase_dataset_stage_merge(const gs_varset_type *gv,
                                        const st_phase_data *d,
                                        const st_phase_proposed_output *p,
                                        st_phase_dataset_stage *stage)
{
  size_t new_writes = 0u;
  int saw_overlap = 0;
  int saw_conflict = 0;

  if (!st_phase_dataset_stage_structurally_valid(stage) ||
      !st_phase_data_adapter_valid(d) || (d->n_ind != stage->n_ind) ||
      (d->n_marker != stage->n_marker))
    return ST_PHASE_STAGE_INVALID;

  if (!st_phase_dataset_stage_source_current(gv, stage) ||
      !st_phase_stage_workspace_matches_source(gv, d))
    return ST_PHASE_STAGE_STALE;

  if (!st_phase_stage_proposal_valid(d, stage, p))
    return ST_PHASE_STAGE_INVALID;

  for (size_t j = 0; j < p->n_marker; j++)
  {
    size_t slot =
        (size_t)p->source_marker[j] * stage->n_ind + p->individual;

    if (stage->covered[slot])
      saw_overlap = 1;

    if (p->write_phase[j])
    {
      if (stage->claimed[slot] &&
          ((stage->staged1[slot] != p->allele1[j]) ||
           (stage->staged2[slot] != p->allele2[j])))
        saw_conflict = 1;
      if (!stage->claimed[slot])
        new_writes++;
    }
  }

  if (saw_conflict)
    return ST_PHASE_STAGE_CONFLICT;
  if (saw_overlap)
    return ST_PHASE_STAGE_OVERLAP;

  if ((new_writes > stage->n_genotype - stage->n_written) ||
      (p->n_marker > stage->n_genotype - stage->n_covered))
    return ST_PHASE_STAGE_INVALID;

  for (size_t j = 0; j < p->n_marker; j++)
  {
    size_t slot =
        (size_t)p->source_marker[j] * stage->n_ind + p->individual;
    stage->covered[slot] = 1u;
    if (p->write_phase[j])
    {
      stage->staged1[slot] = p->allele1[j];
      stage->staged2[slot] = p->allele2[j];
      stage->claimed[slot] = 1u;
    }
  }

  stage->n_written += new_writes;
  stage->n_covered += p->n_marker;
  return ST_PHASE_STAGE_OK;
}

typedef enum
{
  ST_PHASE_AUDIT_OK = 0,
  ST_PHASE_AUDIT_INVALID = -1,
  ST_PHASE_AUDIT_STALE = -2,
  ST_PHASE_AUDIT_INCOMPLETE = -3,
  ST_PHASE_AUDIT_CORRUPT = -4
} st_phase_audit_status;

typedef struct
{
  size_t n_ind;
  size_t n_marker;
  size_t n_genotype;

  size_t n_covered;
  size_t n_uncovered;
  size_t n_written;
  size_t n_covered_unwritten;

  size_t n_missing;
  size_t n_homozygous;
  size_t n_heterozygous;

  size_t n_covered_missing;
  size_t n_covered_homozygous;
  size_t n_covered_heterozygous;
  size_t n_covered_heterozygous_unresolved;

  size_t n_written_same_source_order;
  size_t n_written_swapped_source_order;

  int complete;
  int source_bytes_current;
  int workspace_current;
} st_phase_stage_audit;

static void st_phase_stage_audit_init(st_phase_stage_audit *x)
{
  if (NULL != x)
    memset(x, 0, sizeof(*x));
}

typedef enum
{
  ST_PHASE_COMMIT_OK = 0,
  ST_PHASE_COMMIT_INVALID = -1,
  ST_PHASE_COMMIT_STALE = -2,
  ST_PHASE_COMMIT_INCOMPLETE = -3,
  ST_PHASE_COMMIT_CORRUPT = -4,
  ST_PHASE_COMMIT_ROLLED_BACK = -5
} st_phase_commit_status;

#define st_phase_commit_test_should_fail(point) 0

static int st_phase_commit_from_audit_status(int audit_status)
{
  switch (audit_status)
  {
  case ST_PHASE_AUDIT_OK:
    return ST_PHASE_COMMIT_OK;
  case ST_PHASE_AUDIT_INVALID:
    return ST_PHASE_COMMIT_INVALID;
  case ST_PHASE_AUDIT_STALE:
    return ST_PHASE_COMMIT_STALE;
  case ST_PHASE_AUDIT_INCOMPLETE:
    return ST_PHASE_COMMIT_INCOMPLETE;
  case ST_PHASE_AUDIT_CORRUPT:
    return ST_PHASE_COMMIT_CORRUPT;
  default:
    return ST_PHASE_COMMIT_INVALID;
  }
}

typedef struct
{
  size_t individual;
  size_t chromosome;
  size_t n_marker;
  size_t n_ref_hap;
  int possible;
  double log_likelihood;
  double viterbi_log_probability;

  st_phase_orientation_result orientation;
  st_phase_hard_call_result hard;
  st_phase_proposed_output proposal;
} st_phase_chromosome_result;

static void st_phase_chromosome_result_init(st_phase_chromosome_result *r)
{
  if (NULL != r)
  {
    memset(r, 0, sizeof(*r));
    r->log_likelihood = -INFINITY;
    r->viterbi_log_probability = -INFINITY;
    st_phase_orientation_result_init(&r->orientation);
    st_phase_hard_call_result_init(&r->hard);
    st_phase_proposed_output_init(&r->proposal);
  }
}

static void st_phase_chromosome_result_free(st_phase_chromosome_result *r)
{
  if (NULL == r)
    return;

  st_phase_orientation_result_free(&r->orientation);
  st_phase_hard_call_result_free(&r->hard);
  st_phase_proposed_output_free(&r->proposal);
  st_phase_chromosome_result_init(r);
}

#define st_phase_driver_malloc malloc
#define st_phase_driver_calloc calloc

typedef struct
{
  size_t n_marker;
  size_t n_ref_hap;
  size_t max_states;
  size_t n_slot;
  size_t *count;
  st_phase_candidate *state;
} st_phase_unphased_candidate_set;

static void
st_phase_unphased_candidate_set_init(st_phase_unphased_candidate_set *s)
{
  if (NULL != s)
    memset(s, 0, sizeof(*s));
}

static void
st_phase_unphased_candidate_set_free(st_phase_unphased_candidate_set *s)
{
  if (NULL == s)
    return;

  free(s->count);
  free(s->state);
  st_phase_unphased_candidate_set_init(s);
}

#define st_phase_unphased_malloc malloc
#define st_phase_unphased_calloc calloc

static st_phase_candidate *
st_phase_unphased_candidate_row(st_phase_unphased_candidate_set *s,
                                size_t marker)
{
  return s->state + marker * s->max_states;
}

static const st_phase_candidate *st_phase_unphased_candidate_row_const(
    const st_phase_unphased_candidate_set *s, size_t marker)
{
  return s->state + marker * s->max_states;
}

static int
st_phase_unphased_candidate_insert(st_phase_unphased_candidate_set *s,
                                   size_t marker, size_t donor,
                                   size_t score)
{
  size_t n, pos;
  st_phase_candidate x, *v;

  if ((NULL == s) || (NULL == s->count) || (marker >= s->n_marker) ||
      (donor >= s->n_ref_hap))
    return -1;

  if (0u == s->max_states)
    return 0;

  n = s->count[marker];
  if ((n > s->max_states) || (NULL == s->state))
    return -1;

  v = st_phase_unphased_candidate_row(s, marker);
  x.hap = donor;
  x.ibs_length = score;

  for (pos = 0u; pos < n; pos++)
    if (v[pos].hap == donor)
    {
      if (v[pos].ibs_length >= score)
        return 0;
      v[pos] = x;
      while ((pos > 0u) &&
             st_phase_candidate_better(&v[pos], &v[pos - 1u]))
      {
        st_phase_candidate tmp = v[pos - 1u];
        v[pos - 1u] = v[pos];
        v[pos] = tmp;
        pos--;
      }
      return 0;
    }

  if (n < s->max_states)
  {
    v[n] = x;
    s->count[marker] = n + 1u;
    pos = n;
  }
  else
  {
    if (!st_phase_candidate_better(&x, &v[n - 1u]))
      return 0;
    v[n - 1u] = x;
    pos = n - 1u;
  }

  while ((pos > 0u) && st_phase_candidate_better(&v[pos], &v[pos - 1u]))
  {
    st_phase_candidate tmp = v[pos - 1u];
    v[pos - 1u] = v[pos];
    v[pos] = tmp;
    pos--;
  }

  return 0;
}

static int st_phase_unphased_marker_compatible(uint8_t observed1,
                                               uint8_t observed2,
                                               uint8_t donor)
{
  if (ST_PHASE_STATE_MISSING == observed1)
    return 1;
  return (donor == observed1) || (donor == observed2);
}

static int st_phase_unphased_candidate_input_valid(
    const uint8_t *observed1, const uint8_t *observed2,
    const uint8_t *n_alleles, const uint8_t *reference_hap,
    size_t n_marker, size_t n_ref_hap)
{
  if ((NULL == observed1) || (NULL == observed2) || (NULL == n_alleles) ||
      (NULL == reference_hap) || (0u == n_marker) || (0u == n_ref_hap))
    return 0;

  if (n_ref_hap > SIZE_MAX / n_marker)
    return 0;

  for (size_t marker = 0u; marker < n_marker; marker++)
  {
    uint8_t A = n_alleles[marker];
    uint8_t a = observed1[marker];
    uint8_t b = observed2[marker];
    int ma = (ST_PHASE_STATE_MISSING == a);
    int mb = (ST_PHASE_STATE_MISSING == b);

    if ((0u == A) || (A > 2u) || (ST_PHASE_STATE_MISSING == A) ||
        (ma != mb))
      return 0;

    if (!ma && ((a >= A) || (b >= A)))
      return 0;

    for (size_t h = 0u; h < n_ref_hap; h++)
      if (reference_hap[marker * n_ref_hap + h] >= A)
        return 0;
  }

  return 1;
}

static int st_phase_unphased_candidate_build(
    const uint8_t *observed1, const uint8_t *observed2,
    const uint8_t *n_alleles, const uint8_t *reference_hap,
    size_t n_marker, size_t n_ref_hap, size_t max_states,
    st_phase_unphased_candidate_set *out)
{
  st_phase_unphased_candidate_set tmp;
  size_t capacity, n_slot;

  st_phase_unphased_candidate_set_init(&tmp);

  if ((NULL == out) || (0u == max_states) ||
      !st_phase_unphased_candidate_input_valid(observed1, observed2,
                                               n_alleles, reference_hap,
                                               n_marker, n_ref_hap))
    return -1;

  capacity = (max_states < n_ref_hap) ? max_states : n_ref_hap;
  if ((0u == capacity) ||
      (0 != st_phase_size_mul(n_marker, capacity, &n_slot)) ||
      (n_marker > SIZE_MAX / sizeof(*tmp.count)) ||
      (n_slot > SIZE_MAX / sizeof(*tmp.state)))
    return -1;

  tmp.n_marker = n_marker;
  tmp.n_ref_hap = n_ref_hap;
  tmp.max_states = capacity;
  tmp.n_slot = n_slot;

  tmp.count = st_phase_unphased_calloc(n_marker, sizeof(*tmp.count));
  if (NULL == tmp.count)
    goto fail;

  tmp.state = st_phase_unphased_malloc(n_slot * sizeof(*tmp.state));
  if (NULL == tmp.state)
    goto fail;

  for (size_t z = 0u; z < n_slot; z++)
  {
    tmp.state[z].hap = SIZE_MAX;
    tmp.state[z].ibs_length = 0u;
  }

  for (size_t donor = 0u; donor < n_ref_hap; donor++)
  {
    size_t marker = 0u;

    while (marker < n_marker)
    {
      if (!st_phase_unphased_marker_compatible(
              observed1[marker], observed2[marker],
              reference_hap[marker * n_ref_hap + donor]))
      {
        marker++;
        continue;
      }

      {
        size_t start = marker;
        size_t length;

        do
          marker++;
        while ((marker < n_marker) &&
               st_phase_unphased_marker_compatible(
                   observed1[marker], observed2[marker],
                   reference_hap[marker * n_ref_hap + donor]));

        length = marker - start;
        for (size_t j = start; j < marker; j++)
          if (0 !=
              st_phase_unphased_candidate_insert(&tmp, j, donor, length))
            goto fail;
      }
    }
  }

  for (size_t marker = 0u; marker < n_marker; marker++)
  {
    for (size_t donor = 0u;
         (tmp.count[marker] < capacity) && (donor < n_ref_hap); donor++)
    {
      int present = 0;
      st_phase_candidate *row =
          st_phase_unphased_candidate_row(&tmp, marker);

      for (size_t z = 0u; z < tmp.count[marker]; z++)
        if (row[z].hap == donor)
        {
          present = 1;
          break;
        }

      if (!present && (0 != st_phase_unphased_candidate_insert(
                                &tmp, marker, donor, 0u)))
        goto fail;
    }

    if (tmp.count[marker] != capacity)
      goto fail;
  }

  st_phase_unphased_candidate_set_free(out);
  *out = tmp;
  return 0;

fail:
  st_phase_unphased_candidate_set_free(&tmp);
  return -1;
}

static const st_phase_candidate *
st_phase_unphased_candidate_get(const st_phase_unphased_candidate_set *s,
                                size_t marker, size_t *count)
{
  if (NULL != count)
    *count = 0u;

  if ((NULL == s) || (marker >= s->n_marker) || (0u == s->max_states) ||
      (NULL == s->count) || (NULL == s->state) ||
      (s->count[marker] > s->max_states))
    return NULL;

  if (NULL != count)
    *count = s->count[marker];

  return st_phase_unphased_candidate_row_const(s, marker);
}

static int st_phase_unphased_candidate_set_valid(
    const st_phase_unphased_candidate_set *s)
{
  size_t expected_slot;

  if ((NULL == s) || (0u == s->n_marker) || (0u == s->n_ref_hap) ||
      (0u == s->max_states) || (s->max_states > s->n_ref_hap) ||
      (NULL == s->count) || (NULL == s->state) ||
      (0 !=
       st_phase_size_mul(s->n_marker, s->max_states, &expected_slot)) ||
      (expected_slot != s->n_slot))
    return 0;

  for (size_t marker = 0u; marker < s->n_marker; marker++)
  {
    size_t n = s->count[marker];
    const st_phase_candidate *row;

    if ((0u == n) || (n != s->max_states))
      return 0;

    row = st_phase_unphased_candidate_row_const(s, marker);
    if (NULL == row)
      return 0;

    for (size_t z = 0u; z < n; z++)
    {
      if ((SIZE_MAX == row[z].hap) || (row[z].hap >= s->n_ref_hap) ||
          (row[z].ibs_length > s->n_marker))
        return 0;

      if ((z > 0u) && st_phase_candidate_better(&row[z], &row[z - 1u]))
        return 0;

      for (size_t j = 0u; j < z; j++)
        if (row[j].hap == row[z].hap)
          return 0;
    }
  }

  return 1;
}

static int st_phase_candidate_row_insert_local(st_phase_candidate *row,
                                               size_t *n, size_t capacity,
                                               size_t donor, size_t score)
{
  st_phase_candidate x;
  size_t pos;
  if ((NULL == row) || (NULL == n) || (0u == capacity) || (*n > capacity))
    return -1;
  x.hap = donor;
  x.ibs_length = score;

  if (*n == capacity &&
      !st_phase_candidate_better(&x, &row[capacity - 1u]))
    return 0;
  for (pos = 0u; pos < *n; pos++)
    if (row[pos].hap == donor)
    {
      if (row[pos].ibs_length >= score)
        return 0;
      row[pos] = x;
      while ((pos > 0u) &&
             st_phase_candidate_better(&row[pos], &row[pos - 1u]))
      {
        st_phase_candidate z = row[pos - 1u];
        row[pos - 1u] = row[pos];
        row[pos] = z;
        pos--;
      }
      return 0;
    }
  if (*n < capacity)
  {
    row[*n] = x;
    pos = *n;
    (*n)++;
  }
  else
  {
    if (!st_phase_candidate_better(&x, &row[capacity - 1u]))
      return 0;
    row[capacity - 1u] = x;
    pos = capacity - 1u;
  }
  while ((pos > 0u) &&
         st_phase_candidate_better(&row[pos], &row[pos - 1u]))
  {
    st_phase_candidate z = row[pos - 1u];
    row[pos - 1u] = row[pos];
    row[pos] = z;
    pos--;
  }
  return 0;
}

static int st_phase_unphased_candidate_input_valid_cohort(
    const uint8_t *observed1, const uint8_t *observed2,
    const uint8_t *n_alleles, const uint8_t *reference_hap,
    size_t n_marker, size_t n_ref_hap)
{
  if (!observed1 || !observed2 || !n_alleles || !reference_hap ||
      0u == n_marker || 0u == n_ref_hap || n_ref_hap > SIZE_MAX / n_marker)
    return 0;
  for (size_t m = 0u; m < n_marker; m++)
  {
    uint8_t A = n_alleles[m], a = observed1[m], b = observed2[m];
    int ma = (ST_PHASE_STATE_MISSING == a),
        mb = (ST_PHASE_STATE_MISSING == b);
    if (A > 2u || ST_PHASE_STATE_MISSING == A || ma != mb)
      return 0;
    if (0u == A)
    {
      if (!ma)
        return 0;

      for (size_t h = 0u; h < n_ref_hap; h++)
        if (reference_hap[m * n_ref_hap + h] != 0u)
          return 0;
      continue;
    }
    if (!ma && (a >= A || b >= A))
      return 0;
    for (size_t h = 0u; h < n_ref_hap; h++)
      if (reference_hap[m * n_ref_hap + h] >= A)
        return 0;
  }
  return 1;
}

static int st_phase_unphased_candidate_build_cohort_balanced(
    const uint8_t *observed1, const uint8_t *observed2,
    const uint8_t *n_alleles, const uint8_t *reference_hap,
    size_t n_marker, size_t n_ref_hap, size_t max_states,
    st_phase_unphased_candidate_set *out)
{
  st_phase_unphased_candidate_set tmp;
  st_phase_candidate *by0 = NULL, *by1 = NULL;
  size_t *n0 = NULL, *n1 = NULL;
  size_t capacity, n_slot;
  uint8_t *donor_major = NULL;
  int status = -1;
  st_phase_unphased_candidate_set_init(&tmp);
  if ((NULL == out) || (0u == max_states) ||
      !st_phase_unphased_candidate_input_valid_cohort(
          observed1, observed2, n_alleles, reference_hap, n_marker,
          n_ref_hap))
    return -1;
  capacity = (max_states < n_ref_hap) ? max_states : n_ref_hap;
  if ((0u == capacity) ||
      (0 != st_phase_size_mul(n_marker, capacity, &n_slot)) ||
      (n_marker > SIZE_MAX / sizeof(*tmp.count)) ||
      (n_slot > SIZE_MAX / sizeof(*tmp.state)))
    return -1;
  tmp.n_marker = n_marker;
  tmp.n_ref_hap = n_ref_hap;
  tmp.max_states = capacity;
  tmp.n_slot = n_slot;
  tmp.count = st_phase_unphased_calloc(n_marker, sizeof(*tmp.count));
  tmp.state = st_phase_unphased_malloc(n_slot * sizeof(*tmp.state));
  n0 = st_phase_unphased_calloc(n_marker, sizeof(*n0));
  n1 = st_phase_unphased_calloc(n_marker, sizeof(*n1));
  by0 = st_phase_unphased_malloc(n_slot * sizeof(*by0));
  by1 = st_phase_unphased_malloc(n_slot * sizeof(*by1));
  if ((NULL == tmp.count) || (NULL == tmp.state) || (NULL == n0) ||
      (NULL == n1) || (NULL == by0) || (NULL == by1))
    goto cleanup;
  for (size_t z = 0u; z < n_slot; z++)
  {
    tmp.state[z].hap = by0[z].hap = by1[z].hap = SIZE_MAX;
    tmp.state[z].ibs_length = by0[z].ibs_length = by1[z].ibs_length = 0u;
  }

  donor_major = st_phase_unphased_malloc(n_marker * n_ref_hap);
  if (!donor_major)
    goto cleanup;
  for (size_t bm = 0u; bm < n_marker; )
  {
    size_t me = n_marker - bm < 32u ? n_marker : bm + 32u;
    for (size_t br = 0u; br < n_ref_hap; )
    {
      size_t re = n_ref_hap - br < 32u ? n_ref_hap : br + 32u;
      for (size_t r = br; r < re; r++)
        for (size_t m = bm; m < me; m++)
          donor_major[r * n_marker + m] = reference_hap[m * n_ref_hap + r];
      br = re;
    }
    bm = me;
  }
  for (size_t donor = 0u; donor < n_ref_hap; donor++)
  {
    const uint8_t *donor_hap = donor_major + donor * n_marker;
    size_t marker = 0u;
    while (marker < n_marker)
    {
      if (!st_phase_unphased_marker_compatible(
              observed1[marker], observed2[marker],
              donor_hap[marker]))
      {
        marker++;
        continue;
      }
      size_t start = marker;
      do
        marker++;
      while ((marker < n_marker) &&
             st_phase_unphased_marker_compatible(
                 observed1[marker], observed2[marker],
                 donor_hap[marker]));
      size_t len = marker - start;
      for (size_t m = start; m < marker; m++)
      {
        uint8_t a = donor_hap[m];
        st_phase_candidate *row =
            (0u == a) ? by0 + m * capacity : by1 + m * capacity;
        size_t *nr = (0u == a) ? n0 + m : n1 + m;
        if (0 != st_phase_candidate_row_insert_local(row, nr, capacity,
                                                     donor, len))
          goto cleanup;
      }
    }
  }

  for (size_t m = 0u; m < n_marker; m++)
    for (size_t donor = 0u; donor < n_ref_hap; donor++)
    {

      if (n0[m] >= capacity &&
          (n1[m] >= capacity || n_alleles[m] <= 1u))
        break;
      uint8_t a = reference_hap[m * n_ref_hap + donor];
      st_phase_candidate *row =
          (0u == a) ? by0 + m * capacity : by1 + m * capacity;
      size_t *nr = (0u == a) ? n0 + m : n1 + m;
      if (*nr < capacity && 0 != st_phase_candidate_row_insert_local(
                                     row, nr, capacity, donor, 0u))
        goto cleanup;
    }

  for (size_t m = 0u; m < n_marker; m++)
  {
    st_phase_candidate *dst = tmp.state + m * capacity;
    st_phase_candidate *r0 = by0 + m * capacity, *r1 = by1 + m * capacity;
    size_t want0 = 0u, want1 = 0u, take0, take1, n = 0u;
    uint8_t A = n_alleles[m], a = observed1[m], b = observed2[m];
    if (1u == A)
    {
      want0 = capacity;
      want1 = 0u;
    }
    else if ((ST_PHASE_STATE_MISSING == a) || (a != b))
    {
      want0 = (capacity + 1u) / 2u;
      want1 = capacity / 2u;
    }
    else if (0u == a)
    {
      want0 = capacity;
      want1 = 0u;
    }
    else
    {
      want0 = 0u;
      want1 = capacity;
    }
    take0 = (want0 < n0[m]) ? want0 : n0[m];
    take1 = (want1 < n1[m]) ? want1 : n1[m];
    while (take0 + take1 < capacity)
    {
      int can0 = (take0 < n0[m]), can1 = (take1 < n1[m]);
      if (!can0 && !can1)
        goto cleanup;
      if (can0 &&
          (!can1 || st_phase_candidate_better(&r0[take0], &r1[take1])))
        take0++;
      else
        take1++;
    }
    for (size_t z = 0u; z < take0; z++)
      dst[n++] = r0[z];
    for (size_t z = 0u; z < take1; z++)
      dst[n++] = r1[z];
    if (n != capacity)
      goto cleanup;
    for (size_t x = 1u; x < n; x++)
    {
      st_phase_candidate v = dst[x];
      size_t y = x;
      while ((y > 0u) && st_phase_candidate_better(&v, &dst[y - 1u]))
      {
        dst[y] = dst[y - 1u];
        y--;
      }
      dst[y] = v;
    }
    tmp.count[m] = capacity;
  }
  if (!st_phase_unphased_candidate_set_valid(&tmp))
    goto cleanup;
  st_phase_unphased_candidate_set_free(out);
  *out = tmp;
  st_phase_unphased_candidate_set_init(&tmp);
  status = 0;
cleanup:
  free(donor_major);
  free(n0);
  free(n1);
  free(by0);
  free(by1);
  st_phase_unphased_candidate_set_free(&tmp);
  return status;
}

#define st_phase_bounded_malloc malloc
#define st_phase_bounded_calloc calloc

static int st_phase_pair_state_build_unphased(
    const st_phase_unphased_candidate_set *candidate,
    st_phase_pair_state_set *out)
{
  st_phase_pair_state_set tmp;
  size_t max_pair, n_slot;

  st_phase_pair_state_set_init(&tmp);

  if ((NULL == out) || !st_phase_unphased_candidate_set_valid(candidate) ||
      (0 != st_phase_size_mul(candidate->max_states, candidate->max_states,
                              &max_pair)) ||
      (0u == max_pair) ||
      (0 != st_phase_size_mul(candidate->n_marker, max_pair, &n_slot)) ||
      (candidate->n_marker > SIZE_MAX / sizeof(size_t)) ||
      (n_slot > SIZE_MAX / sizeof(st_phase_pair_state)))
    return -1;

  tmp.n_marker = candidate->n_marker;
  tmp.n_hap = candidate->n_ref_hap;
  tmp.target1 = 0u;
  tmp.target2 = 0u;
  tmp.max_side_states = candidate->max_states;
  tmp.max_pair_states = max_pair;
  tmp.n_slot = n_slot;

  tmp.left_count =
      st_phase_bounded_malloc(tmp.n_marker * sizeof(*tmp.left_count));
  tmp.right_count =
      st_phase_bounded_malloc(tmp.n_marker * sizeof(*tmp.right_count));
  tmp.count = st_phase_bounded_malloc(tmp.n_marker * sizeof(*tmp.count));
  tmp.state = st_phase_bounded_malloc(tmp.n_slot * sizeof(*tmp.state));

  if ((NULL == tmp.left_count) || (NULL == tmp.right_count) ||
      (NULL == tmp.count) || (NULL == tmp.state))
    goto fail;

  for (size_t z = 0u; z < tmp.n_slot; z++)
  {
    tmp.state[z].hap1 = SIZE_MAX;
    tmp.state[z].hap2 = SIZE_MAX;
    tmp.state[z].slot1 = SIZE_MAX;
    tmp.state[z].slot2 = SIZE_MAX;
  }

  for (size_t marker = 0u; marker < tmp.n_marker; marker++)
  {
    size_t n = 0u;
    const st_phase_candidate *row =
        st_phase_unphased_candidate_get(candidate, marker, &n);
    st_phase_pair_state *pair_row =
        tmp.state + marker * tmp.max_pair_states;
    size_t pair_n;

    if ((NULL == row) || (0u == n) || (n != candidate->max_states) ||
        (0 != st_phase_size_mul(n, n, &pair_n)) ||
        (pair_n > tmp.max_pair_states))
      goto fail;

    tmp.left_count[marker] = n;
    tmp.right_count[marker] = n;
    tmp.count[marker] = pair_n;

    for (size_t i = 0u; i < n; i++)
      for (size_t j = 0u; j < n; j++)
      {
        size_t z = i * n + j;
        pair_row[z].hap1 = row[i].hap;
        pair_row[z].hap2 = row[j].hap;
        pair_row[z].slot1 = i;
        pair_row[z].slot2 = j;
      }
  }

  for (size_t marker = 0u; marker < tmp.n_marker; marker++)
    if (!st_phase_pair_state_row_valid(&tmp, marker))
      goto fail;

  st_phase_pair_state_set_free(out);
  *out = tmp;
  return 0;

fail:
  st_phase_pair_state_set_free(&tmp);
  return -1;
}

#define st_phase_factor_malloc malloc
#define st_phase_factor_calloc calloc

static void st_phase_factor_best_update(double score, size_t identity,
                                        double *best_score,
                                        size_t *best_identity)
{
  if ((score > *best_score) ||
      ((score == *best_score) &&
       ((SIZE_MAX == *best_identity) || (identity < *best_identity))))
  {
    *best_score = score;
    *best_identity = identity;
  }
}

static void st_phase_factor_top2_update(double score, size_t identity,
                                        double *first_score,
                                        size_t *first_identity,
                                        double *second_score,
                                        size_t *second_identity)
{
  if ((score > *first_score) ||
      ((score == *first_score) &&
       ((SIZE_MAX == *first_identity) || (identity < *first_identity))))
  {
    if (*first_identity != identity)
    {
      *second_score = *first_score;
      *second_identity = *first_identity;
    }
    *first_score = score;
    *first_identity = identity;
  }
  else if ((identity != *first_identity) &&
           ((score > *second_score) ||
            ((score == *second_score) && ((SIZE_MAX == *second_identity) ||
                                          (identity < *second_identity)))))
  {
    *second_score = score;
    *second_identity = identity;
  }
}

static int st_phase_factor_interval_maps(
    const st_phase_pair_state_set *pair, size_t marker,
    size_t *prev_to_cur1, size_t *prev_to_cur2, size_t *cur_to_prev1,
    size_t *cur_to_prev2, size_t *np1_out, size_t *np2_out,
    size_t *nc1_out, size_t *nc2_out)
{
  size_t np, np1, np2, nc, nc1, nc2;
  const st_phase_pair_state *previous, *current;

  if ((NULL == pair) || (0u == marker) || (marker >= pair->n_marker) ||
      (NULL == prev_to_cur1) || (NULL == prev_to_cur2) ||
      (NULL == cur_to_prev1) || (NULL == cur_to_prev2) ||
      (NULL == np1_out) || (NULL == np2_out) || (NULL == nc1_out) ||
      (NULL == nc2_out) ||
      !st_phase_pair_state_row_valid(pair, marker - 1u) ||
      !st_phase_pair_state_row_valid(pair, marker))
    return -1;

  previous = st_phase_pair_state_get(pair, marker - 1u, &np, &np1, &np2);
  current = st_phase_pair_state_get(pair, marker, &nc, &nc1, &nc2);
  if ((NULL == previous) || (NULL == current) || (np != np1 * np2) ||
      (nc != nc1 * nc2) || (0u == np1) || (0u == np2) || (0u == nc1) ||
      (0u == nc2))
    return -1;

  for (size_t a = 0u; a < np1; a++)
    prev_to_cur1[a] = SIZE_MAX;
  for (size_t b = 0u; b < np2; b++)
    prev_to_cur2[b] = SIZE_MAX;
  for (size_t c = 0u; c < nc1; c++)
    cur_to_prev1[c] = SIZE_MAX;
  for (size_t d = 0u; d < nc2; d++)
    cur_to_prev2[d] = SIZE_MAX;

  for (size_t a = 0u; a < np1; a++)
    for (size_t c = 0u; c < nc1; c++)
      if (previous[a * np2].hap1 == current[c * nc2].hap1)
      {
        prev_to_cur1[a] = c;
        cur_to_prev1[c] = a;
      }

  for (size_t b = 0u; b < np2; b++)
    for (size_t d = 0u; d < nc2; d++)
      if (previous[b].hap2 == current[d].hap2)
      {
        prev_to_cur2[b] = d;
        cur_to_prev2[d] = b;
      }

  *np1_out = np1;
  *np2_out = np2;
  *nc1_out = nc1;
  *nc2_out = nc2;
  return 0;
}

static double st_phase_factor_log_base(double q, size_t n_current,
                                       int retained)
{
  if (!retained)
    return -log((double)n_current);
  return st_phase_hmm_log_probability(q / (double)n_current);
}

static double st_phase_factor_log_same(double q, size_t n_current)
{
  return st_phase_hmm_log_probability((1.0 - q) + q / (double)n_current);
}

static int st_phase_factor_forward_transition(
    const double *log_prev, size_t np1, size_t np2, size_t nc1, size_t nc2,
    const size_t *prev_to_cur1, const size_t *prev_to_cur2,
    const size_t *cur_to_prev1, const size_t *cur_to_prev2, double q1,
    double q2, double *side1, double *side2, double *log_mass)
{
  double log_const = -INFINITY;
  double log_stay1 = st_phase_hmm_log_probability(1.0 - q1);
  double log_stay2 = st_phase_hmm_log_probability(1.0 - q2);

  if ((NULL == log_prev) || (NULL == prev_to_cur1) ||
      (NULL == prev_to_cur2) || (NULL == cur_to_prev1) ||
      (NULL == cur_to_prev2) || (NULL == side1) || (NULL == side2) ||
      (NULL == log_mass) || (0u == np1) || (0u == np2) || (0u == nc1) ||
      (0u == nc2) || !isfinite(q1) || !isfinite(q2) || (q1 < 0.0) ||
      (q1 > 1.0) || (q2 < 0.0) || (q2 > 1.0))
    return -1;

  const double retained_base1 = st_phase_factor_log_base(q1, nc1, 1);
  const double absent_base1 = st_phase_factor_log_base(q1, nc1, 0);
  const double retained_base2 = st_phase_factor_log_base(q2, nc2, 1);
  const double absent_base2 = st_phase_factor_log_base(q2, nc2, 0);

  for (size_t c = 0u; c < nc1; c++)
    side1[c] = -INFINITY;
  for (size_t d = 0u; d < nc2; d++)
    side2[d] = -INFINITY;

  for (size_t a = 0u; a < np1; a++)
  {
    double lb1 =
        (SIZE_MAX != prev_to_cur1[a] ? retained_base1 : absent_base1);
    for (size_t b = 0u; b < np2; b++)
    {
      double lf = log_prev[a * np2 + b];
      double lb2 =
          (SIZE_MAX != prev_to_cur2[b] ? retained_base2 : absent_base2);
      if ((lf > -INFINITY) && (lb1 > -INFINITY) && (lb2 > -INFINITY))
        log_const = st_phase_hmm_logadd(log_const, lf + lb1 + lb2);
    }
  }

  if (log_stay1 > -INFINITY)
    for (size_t c = 0u; c < nc1; c++)
    {
      size_t a = cur_to_prev1[c];
      double v = -INFINITY;
      if (SIZE_MAX == a)
        continue;
      for (size_t b = 0u; b < np2; b++)
      {
        double lf = log_prev[a * np2 + b];
        double lb2 =
            (SIZE_MAX != prev_to_cur2[b] ? retained_base2 : absent_base2);
        if ((lf > -INFINITY) && (lb2 > -INFINITY))
          v = st_phase_hmm_logadd(v, lf + lb2);
      }
      if (v > -INFINITY)
        side1[c] = v + log_stay1;
    }

  if (log_stay2 > -INFINITY)
    for (size_t d = 0u; d < nc2; d++)
    {
      size_t b = cur_to_prev2[d];
      double v = -INFINITY;
      if (SIZE_MAX == b)
        continue;
      for (size_t a = 0u; a < np1; a++)
      {
        double lf = log_prev[a * np2 + b];
        double lb1 =
            (SIZE_MAX != prev_to_cur1[a] ? retained_base1 : absent_base1);
        if ((lf > -INFINITY) && (lb1 > -INFINITY))
          v = st_phase_hmm_logadd(v, lf + lb1);
      }
      if (v > -INFINITY)
        side2[d] = v + log_stay2;
    }

  for (size_t c = 0u; c < nc1; c++)
  {

    double prefix = st_phase_hmm_logadd(log_const, side1[c]);
    for (size_t d = 0u; d < nc2; d++)
    {
      double v = prefix;
      size_t a = cur_to_prev1[c], b = cur_to_prev2[d];
      v = st_phase_hmm_logadd(v, side2[d]);
      if ((SIZE_MAX != a) && (SIZE_MAX != b) && (log_stay1 > -INFINITY) &&
          (log_stay2 > -INFINITY))
      {
        double lf = log_prev[a * np2 + b];
        if (lf > -INFINITY)
          v = st_phase_hmm_logadd(v, lf + log_stay1 + log_stay2);
      }
      log_mass[c * nc2 + d] = v;
    }
  }

  return 0;
}

static int st_phase_factor_viterbi_side2(
    const double *vit_prev, size_t np1, size_t np2, size_t nc2,
    const size_t *prev_to_cur2, const size_t *cur_to_prev2, double q2,
    double *stage_score, size_t *stage_prev_b)
{
  double log_same = st_phase_factor_log_same(q2, nc2);
  double log_other = st_phase_hmm_log_probability(q2 / (double)nc2);
  double log_absent = -log((double)nc2);

  if ((NULL == vit_prev) || (NULL == prev_to_cur2) ||
      (NULL == cur_to_prev2) || (NULL == stage_score) ||
      (NULL == stage_prev_b) || (0u == np1) || (0u == np2) || (0u == nc2))
    return -1;

  for (size_t a = 0u; a < np1; a++)
  {
    double r1 = -INFINITY, r2 = -INFINITY, ab = -INFINITY;
    size_t rb1 = SIZE_MAX, rb2 = SIZE_MAX, abb = SIZE_MAX;

    for (size_t b = 0u; b < np2; b++)
    {
      double raw = vit_prev[a * np2 + b];
      if (SIZE_MAX != prev_to_cur2[b])
        st_phase_factor_top2_update(raw, b, &r1, &rb1, &r2, &rb2);
      else
        st_phase_factor_best_update(raw, b, &ab, &abb);
    }

    for (size_t d = 0u; d < nc2; d++)
    {
      double best = -INFINITY;
      size_t best_b = SIZE_MAX;
      size_t same = cur_to_prev2[d];

      if ((SIZE_MAX != same) && (vit_prev[a * np2 + same] > -INFINITY) &&
          (log_same > -INFINITY))
        st_phase_factor_best_update(vit_prev[a * np2 + same] + log_same,
                                    same, &best, &best_b);

      {
        size_t other = (rb1 != same) ? rb1 : rb2;
        double raw = (rb1 != same) ? r1 : r2;
        if ((SIZE_MAX != other) && (raw > -INFINITY) &&
            (log_other > -INFINITY))
          st_phase_factor_best_update(raw + log_other, other, &best,
                                      &best_b);
      }

      if ((SIZE_MAX != abb) && (ab > -INFINITY) &&
          (log_absent > -INFINITY))
        st_phase_factor_best_update(ab + log_absent, abb, &best, &best_b);

      stage_score[a * nc2 + d] = best;
      stage_prev_b[a * nc2 + d] = best_b;
    }
  }

  return 0;
}

static int st_phase_factor_viterbi_side1(
    const double *stage_score, const size_t *stage_prev_b, size_t np1,
    size_t np2, size_t nc1, size_t nc2, const size_t *prev_to_cur1,
    const size_t *cur_to_prev1, double q1, double *vit_cur,
    size_t *back_row)
{
  double log_same = st_phase_factor_log_same(q1, nc1);
  double log_other = st_phase_hmm_log_probability(q1 / (double)nc1);
  double log_absent = -log((double)nc1);

  if ((NULL == stage_score) || (NULL == stage_prev_b) ||
      (NULL == prev_to_cur1) || (NULL == cur_to_prev1) ||
      (NULL == vit_cur) || (NULL == back_row) || (0u == np1) ||
      (0u == np2) || (0u == nc1) || (0u == nc2))
    return -1;

  for (size_t d = 0u; d < nc2; d++)
  {
    double r1 = -INFINITY, r2 = -INFINITY, ab = -INFINITY;
    size_t ra1 = SIZE_MAX, ra2 = SIZE_MAX, aba = SIZE_MAX;

    for (size_t a = 0u; a < np1; a++)
    {
      double raw = stage_score[a * nc2 + d];
      if (SIZE_MAX != prev_to_cur1[a])
        st_phase_factor_top2_update(raw, a, &r1, &ra1, &r2, &ra2);
      else
        st_phase_factor_best_update(raw, a, &ab, &aba);
    }

    for (size_t c = 0u; c < nc1; c++)
    {
      double best = -INFINITY;
      size_t best_flat = SIZE_MAX;
      size_t same = cur_to_prev1[c];
      size_t z = c * nc2 + d;

#define ST_PHASE_FACTOR_VIT_CANDIDATE(A, SCORE, LOGP)                     \
  do                                                                      \
  {                                                                       \
    size_t st_phase_factor_a_ = (A);                                      \
    double st_phase_factor_raw_ = (SCORE);                                \
    double st_phase_factor_lp_ = (LOGP);                                  \
    if ((SIZE_MAX != st_phase_factor_a_) &&                               \
        (st_phase_factor_raw_ > -INFINITY) &&                             \
        (st_phase_factor_lp_ > -INFINITY))                                \
    {                                                                     \
      size_t st_phase_factor_b_ =                                         \
          stage_prev_b[st_phase_factor_a_ * nc2 + d];                     \
      if (SIZE_MAX != st_phase_factor_b_)                                 \
      {                                                                   \
        size_t st_phase_factor_flat_ =                                    \
            st_phase_factor_a_ * np2 + st_phase_factor_b_;                \
        double st_phase_factor_v_ =                                       \
            st_phase_factor_raw_ + st_phase_factor_lp_;                   \
        if ((st_phase_factor_v_ > best) ||                                \
            ((st_phase_factor_v_ == best) &&                              \
             ((SIZE_MAX == best_flat) ||                                  \
              (st_phase_factor_flat_ < best_flat))))                      \
        {                                                                 \
          best = st_phase_factor_v_;                                      \
          best_flat = st_phase_factor_flat_;                              \
        }                                                                 \
      }                                                                   \
    }                                                                     \
  } while (0)

      if (SIZE_MAX != same)
        ST_PHASE_FACTOR_VIT_CANDIDATE(same, stage_score[same * nc2 + d],
                                      log_same);

      {
        size_t other = (ra1 != same) ? ra1 : ra2;
        double raw = (ra1 != same) ? r1 : r2;
        ST_PHASE_FACTOR_VIT_CANDIDATE(other, raw, log_other);
      }

      ST_PHASE_FACTOR_VIT_CANDIDATE(aba, ab, log_absent);
#undef ST_PHASE_FACTOR_VIT_CANDIDATE

      vit_cur[z] = best;
      back_row[z] = best_flat;
    }
  }

  return 0;
}

static int st_phase_factor_backward_transition(
    const double *log_next, const double *emission_current, size_t np1,
    size_t np2, size_t nc1, size_t nc2, const size_t *prev_to_cur1,
    const size_t *prev_to_cur2, double q1, double q2, double *row_log,
    double *col_log, double *x_log, double *log_previous)
{
  double total = -INFINITY;
  double log_stay1 = st_phase_hmm_log_probability(1.0 - q1);
  double log_stay2 = st_phase_hmm_log_probability(1.0 - q2);

  if ((NULL == log_next) || (NULL == emission_current) ||
      (NULL == prev_to_cur1) || (NULL == prev_to_cur2) ||
      (NULL == row_log) || (NULL == col_log) || (NULL == x_log) ||
      (NULL == log_previous) || (0u == np1) || (0u == np2) ||
      (0u == nc1) || (0u == nc2))
    return -1;

  const double retained_base1 = st_phase_factor_log_base(q1, nc1, 1);
  const double absent_base1 = st_phase_factor_log_base(q1, nc1, 0);
  const double retained_base2 = st_phase_factor_log_base(q2, nc2, 1);
  const double absent_base2 = st_phase_factor_log_base(q2, nc2, 0);

  for (size_t c = 0u; c < nc1; c++)
    row_log[c] = -INFINITY;
  for (size_t d = 0u; d < nc2; d++)
    col_log[d] = -INFINITY;

  for (size_t c = 0u; c < nc1; c++)
    for (size_t d = 0u; d < nc2; d++)
    {
      size_t z = c * nc2 + d;
      double le = st_phase_hmm_log_probability(emission_current[z]);
      double x = ((le > -INFINITY) && (log_next[z] > -INFINITY))
                     ? le + log_next[z]
                     : -INFINITY;
      x_log[z] = x;
      total = st_phase_hmm_logadd(total, x);
      row_log[c] = st_phase_hmm_logadd(row_log[c], x);
      col_log[d] = st_phase_hmm_logadd(col_log[d], x);
    }

  for (size_t a = 0u; a < np1; a++)
  {
    int retained1 = (SIZE_MAX != prev_to_cur1[a]);
    double lb1 = (retained1 ? retained_base1 : absent_base1);
    size_t c = prev_to_cur1[a];
    double prefix[2];
    int prefix_ready[2] = {0, 0};

    for (size_t b = 0u; b < np2; b++)
    {
      int retained2 = (SIZE_MAX != prev_to_cur2[b]);
      double lb2 = (retained2 ? retained_base2 : absent_base2);
      size_t d = prev_to_cur2[b];
      double v = -INFINITY;

      if (!prefix_ready[retained2])
      {
        if ((total > -INFINITY) && (lb1 > -INFINITY) && (lb2 > -INFINITY))
          v = lb1 + lb2 + total;

        if (retained1 && (log_stay1 > -INFINITY) &&
            (row_log[c] > -INFINITY) && (lb2 > -INFINITY))
          v = st_phase_hmm_logadd(v, log_stay1 + lb2 + row_log[c]);
        prefix[retained2] = v;
        prefix_ready[retained2] = 1;
      }
      else
        v = prefix[retained2];

      if (retained2 && (log_stay2 > -INFINITY) &&
          (col_log[d] > -INFINITY) && (lb1 > -INFINITY))
        v = st_phase_hmm_logadd(v, lb1 + log_stay2 + col_log[d]);

      if (retained1 && retained2 && (log_stay1 > -INFINITY) &&
          (log_stay2 > -INFINITY) && (x_log[c * nc2 + d] > -INFINITY))
        v = st_phase_hmm_logadd(v, log_stay1 + log_stay2 +
                                       x_log[c * nc2 + d]);

      log_previous[a * np2 + b] = v;
    }
  }

  return 0;
}

static int st_phase_pair_hmm_forward_backward_viterbi_factorized(
    const st_phase_pair_state_set *pair, const double *initial_weight,
    const double *emission, const double *q1, const double *q2,
    st_phase_pair_hmm_result *out)
{
  st_phase_pair_hmm_result tmp;
  const size_t M = (NULL != pair) ? pair->n_marker : 0u;
  const size_t P = (NULL != pair) ? pair->max_pair_states : 0u;
  const size_t K = (NULL != pair) ? pair->max_side_states : 0u;
  size_t local_slot = 0u, n_first = 0u;
  double max_initial = 0.0, scaled_initial_sum = 0.0;
  double *log_forward = NULL, *log_prev = NULL, *log_cur = NULL;
  double *vit_prev = NULL, *vit_cur = NULL;
  double *log_next = NULL, *log_back_cur = NULL;
  double *factor_mass = NULL, *side1 = NULL, *side2 = NULL;
  double *stage_score = NULL, *row_log = NULL, *col_log = NULL,
         *x_log = NULL;
  size_t *stage_prev_b = NULL, *back = NULL;
  size_t *prev_to_cur1 = NULL, *prev_to_cur2 = NULL;
  size_t *cur_to_prev1 = NULL, *cur_to_prev2 = NULL;
  double log_likelihood = 0.0, log_likelihood_comp = 0.0;
  double viterbi_offset = 0.0, viterbi_offset_comp = 0.0;

  if ((NULL == out) || (NULL == initial_weight) || (NULL == emission) ||
      !st_phase_pair_hmm_set_valid(pair) ||
      ((M > 1u) && ((NULL == q1) || (NULL == q2))) ||
      (0 != st_phase_size_mul(M, P, &local_slot)))
    return -1;

  {
    const st_phase_pair_state *first =
        st_phase_pair_state_get(pair, 0u, &n_first, NULL, NULL);
    if ((NULL == first) || (0u == n_first) || (n_first > P))
      return -1;
    (void)first;
  }

  for (size_t j = 0u; j < n_first; j++)
  {
    double w = initial_weight[j];
    if (!isfinite(w) || (w < 0.0))
      return -1;
    if (w > max_initial)
      max_initial = w;
  }
  if (!(max_initial > 0.0))
    return -1;
  for (size_t j = 0u; j < n_first; j++)
    scaled_initial_sum += initial_weight[j] / max_initial;
  if (!(scaled_initial_sum > 0.0) || !isfinite(scaled_initial_sum))
    return -1;

  for (size_t marker = 0u; marker < M; marker++)
  {
    size_t n;
    const st_phase_pair_state *row =
        st_phase_pair_state_get(pair, marker, &n, NULL, NULL);
    if ((NULL == row) || (0u == n) || (n > P))
      return -1;
    (void)row;
    for (size_t j = 0u; j < n; j++)
    {
      double e = emission[marker * P + j];
      if (!isfinite(e) || (e < 0.0) || (e > 1.0))
        return -1;
    }
  }

  for (size_t marker = 1u; marker < M; marker++)
  {
    double a = q1[marker - 1u], b = q2[marker - 1u];
    if (!isfinite(a) || !isfinite(b) || (a < 0.0) || (a > 1.0) ||
        (b < 0.0) || (b > 1.0))
      return -1;
  }

  st_phase_pair_hmm_result_init(&tmp);
  tmp.n_marker = M;
  tmp.max_pair_states = P;
  tmp.n_slot = local_slot;

  if ((local_slot > SIZE_MAX / sizeof(double)) ||
      (local_slot > SIZE_MAX / sizeof(size_t)) ||
      (M > SIZE_MAX / sizeof(double)) || (M > SIZE_MAX / sizeof(size_t)) ||
      (P > SIZE_MAX / sizeof(double)) || (P > SIZE_MAX / sizeof(size_t)) ||
      (K > SIZE_MAX / sizeof(double)) || (K > SIZE_MAX / sizeof(size_t)))
    return -1;

  tmp.forward = st_phase_factor_calloc(local_slot, sizeof(*tmp.forward));
  tmp.backward = st_phase_factor_calloc(local_slot, sizeof(*tmp.backward));
  tmp.posterior =
      st_phase_factor_calloc(local_slot, sizeof(*tmp.posterior));
  tmp.log_scale = st_phase_factor_malloc(M * sizeof(*tmp.log_scale));
  tmp.viterbi_hap1 = st_phase_factor_malloc(M * sizeof(*tmp.viterbi_hap1));
  tmp.viterbi_hap2 = st_phase_factor_malloc(M * sizeof(*tmp.viterbi_hap2));
  log_forward = st_phase_factor_malloc(local_slot * sizeof(*log_forward));
  log_prev = st_phase_factor_malloc(P * sizeof(*log_prev));
  log_cur = st_phase_factor_malloc(P * sizeof(*log_cur));
  vit_prev = st_phase_factor_malloc(P * sizeof(*vit_prev));
  vit_cur = st_phase_factor_malloc(P * sizeof(*vit_cur));
  log_next = st_phase_factor_malloc(P * sizeof(*log_next));
  log_back_cur = st_phase_factor_malloc(P * sizeof(*log_back_cur));
  factor_mass = st_phase_factor_malloc(P * sizeof(*factor_mass));
  stage_score = st_phase_factor_malloc(P * sizeof(*stage_score));
  stage_prev_b = st_phase_factor_malloc(P * sizeof(*stage_prev_b));
  x_log = st_phase_factor_malloc(P * sizeof(*x_log));
  back = st_phase_factor_malloc(local_slot * sizeof(*back));
  side1 = st_phase_factor_malloc(K * sizeof(*side1));
  side2 = st_phase_factor_malloc(K * sizeof(*side2));
  row_log = st_phase_factor_malloc(K * sizeof(*row_log));
  col_log = st_phase_factor_malloc(K * sizeof(*col_log));
  prev_to_cur1 = st_phase_factor_malloc(K * sizeof(*prev_to_cur1));
  prev_to_cur2 = st_phase_factor_malloc(K * sizeof(*prev_to_cur2));
  cur_to_prev1 = st_phase_factor_malloc(K * sizeof(*cur_to_prev1));
  cur_to_prev2 = st_phase_factor_malloc(K * sizeof(*cur_to_prev2));

  if ((NULL == tmp.forward) || (NULL == tmp.backward) ||
      (NULL == tmp.posterior) || (NULL == tmp.log_scale) ||
      (NULL == tmp.viterbi_hap1) || (NULL == tmp.viterbi_hap2) ||
      (NULL == log_forward) || (NULL == log_prev) || (NULL == log_cur) ||
      (NULL == vit_prev) || (NULL == vit_cur) || (NULL == log_next) ||
      (NULL == log_back_cur) || (NULL == factor_mass) ||
      (NULL == stage_score) || (NULL == stage_prev_b) || (NULL == x_log) ||
      (NULL == back) || (NULL == side1) || (NULL == side2) ||
      (NULL == row_log) || (NULL == col_log) || (NULL == prev_to_cur1) ||
      (NULL == prev_to_cur2) || (NULL == cur_to_prev1) ||
      (NULL == cur_to_prev2))
    goto fail;

  for (size_t z = 0u; z < local_slot; z++)
  {
    log_forward[z] = -INFINITY;
    back[z] = SIZE_MAX;
  }
  for (size_t marker = 0u; marker < M; marker++)
  {
    tmp.log_scale[marker] = -INFINITY;
    tmp.viterbi_hap1[marker] = SIZE_MAX;
    tmp.viterbi_hap2[marker] = SIZE_MAX;
  }

  {
    double log_scaled_sum = log(scaled_initial_sum);
    double log_norm = -INFINITY, vit_norm = -INFINITY;

    for (size_t j = 0u; j < P; j++)
    {
      log_prev[j] = -INFINITY;
      vit_prev[j] = -INFINITY;
    }
    for (size_t j = 0u; j < n_first; j++)
    {
      double lw =
          (initial_weight[j] > 0.0)
              ? log(initial_weight[j]) - log(max_initial) - log_scaled_sum
              : -INFINITY;
      double le = st_phase_hmm_log_probability(emission[j]);
      double raw =
          ((lw > -INFINITY) && (le > -INFINITY)) ? lw + le : -INFINITY;
      log_prev[j] = raw;
      vit_prev[j] = raw;
      log_norm = st_phase_hmm_logadd(log_norm, raw);
      if (raw > vit_norm)
        vit_norm = raw;
    }

    tmp.log_scale[0] = log_norm;
    if (!(log_norm > -INFINITY))
      goto impossible;
    for (size_t j = 0u; j < n_first; j++)
    {
      if (log_prev[j] > -INFINITY)
      {
        log_prev[j] -= log_norm;
        log_forward[j] = log_prev[j];
        tmp.forward[j] = exp(log_prev[j]);
      }
      if (vit_prev[j] > -INFINITY)
        vit_prev[j] -= vit_norm;
    }
    log_likelihood = log_norm;
    viterbi_offset = vit_norm;
  }

  for (size_t marker = 1u; marker < M; marker++)
  {
    size_t np1, np2, nc1, nc2, n_previous, n_current;
    const st_phase_pair_state *previous = st_phase_pair_state_get(
        pair, marker - 1u, &n_previous, NULL, NULL);
    const st_phase_pair_state *current =
        st_phase_pair_state_get(pair, marker, &n_current, NULL, NULL);
    double log_norm = -INFINITY, vit_norm = -INFINITY;

    if ((NULL == previous) || (NULL == current) ||
        (0 != st_phase_factor_interval_maps(
                  pair, marker, prev_to_cur1, prev_to_cur2, cur_to_prev1,
                  cur_to_prev2, &np1, &np2, &nc1, &nc2)) ||
        (n_previous != np1 * np2) || (n_current != nc1 * nc2))
      goto fail;

    for (size_t j = 0u; j < P; j++)
    {
      log_cur[j] = -INFINITY;
      vit_cur[j] = -INFINITY;
    }

    if (0 != st_phase_factor_forward_transition(
                 log_prev, np1, np2, nc1, nc2, prev_to_cur1, prev_to_cur2,
                 cur_to_prev1, cur_to_prev2, q1[marker - 1u],
                 q2[marker - 1u], side1, side2, factor_mass))
      goto fail;

    if ((0 != st_phase_factor_viterbi_side2(
                  vit_prev, np1, np2, nc2, prev_to_cur2, cur_to_prev2,
                  q2[marker - 1u], stage_score, stage_prev_b)) ||
        (0 != st_phase_factor_viterbi_side1(stage_score, stage_prev_b, np1,
                                            np2, nc1, nc2, prev_to_cur1,
                                            cur_to_prev1, q1[marker - 1u],
                                            vit_cur, back + marker * P)))
      goto fail;

    for (size_t j = 0u; j < n_current; j++)
    {
      double le = st_phase_hmm_log_probability(emission[marker * P + j]);
      log_cur[j] = factor_mass[j];
      if ((log_cur[j] > -INFINITY) && (le > -INFINITY))
        log_cur[j] += le;
      else
        log_cur[j] = -INFINITY;
      if ((vit_cur[j] > -INFINITY) && (le > -INFINITY))
        vit_cur[j] += le;
      else
        vit_cur[j] = -INFINITY;
      log_norm = st_phase_hmm_logadd(log_norm, log_cur[j]);
      if (vit_cur[j] > vit_norm)
        vit_norm = vit_cur[j];
    }

    tmp.log_scale[marker] = log_norm;
    if (!(log_norm > -INFINITY))
      goto impossible;

    for (size_t j = 0u; j < n_current; j++)
    {
      if (log_cur[j] > -INFINITY)
      {
        log_cur[j] -= log_norm;
        log_forward[marker * P + j] = log_cur[j];
        tmp.forward[marker * P + j] = exp(log_cur[j]);
      }
      if (vit_cur[j] > -INFINITY)
        vit_cur[j] -= vit_norm;
    }

    st_phase_pair_hmm_compensated_add(log_norm, &log_likelihood,
                                      &log_likelihood_comp);
    st_phase_pair_hmm_compensated_add(vit_norm, &viterbi_offset,
                                      &viterbi_offset_comp);
    {
      double *swap = log_prev;
      log_prev = log_cur;
      log_cur = swap;
      swap = vit_prev;
      vit_prev = vit_cur;
      vit_cur = swap;
    }
  }

  tmp.log_likelihood = log_likelihood;
  tmp.viterbi_log_probability = viterbi_offset;
  tmp.possible = 1;

  {
    size_t n_last, slot = SIZE_MAX;
    const st_phase_pair_state *last =
        st_phase_pair_state_get(pair, M - 1u, &n_last, NULL, NULL);
    double best = -INFINITY;
    if ((NULL == last) || (0u == n_last) || (n_last > P))
      goto fail;
    for (size_t j = 0u; j < n_last; j++)
      if ((vit_prev[j] > best) ||
          ((vit_prev[j] == best) && ((SIZE_MAX == slot) || (j < slot))))
      {
        best = vit_prev[j];
        slot = j;
      }
    if (SIZE_MAX == slot)
      goto fail;
    for (size_t marker = M; marker > 0u; marker--)
    {
      size_t m = marker - 1u, n;
      const st_phase_pair_state *row =
          st_phase_pair_state_get(pair, m, &n, NULL, NULL);
      if ((NULL == row) || (slot >= n))
        goto fail;
      tmp.viterbi_hap1[m] = row[slot].hap1;
      tmp.viterbi_hap2[m] = row[slot].hap2;
      if (m > 0u)
      {
        slot = back[m * P + slot];
        if (SIZE_MAX == slot)
          goto fail;
      }
    }
  }

  {
    size_t n_last;
    const st_phase_pair_state *last =
        st_phase_pair_state_get(pair, M - 1u, &n_last, NULL, NULL);
    double log_terminal;
    if ((NULL == last) || (0u == n_last) || (n_last > P))
      goto fail;
    (void)last;
    log_terminal = -log((double)n_last);
    for (size_t j = 0u; j < P; j++)
      log_next[j] = -INFINITY;
    for (size_t j = 0u; j < n_last; j++)
    {
      log_next[j] = log_terminal;
      tmp.backward[(M - 1u) * P + j] = 1.0 / (double)n_last;
    }
    if (0 != st_phase_hmm_posterior_log_row(log_forward + (M - 1u) * P,
                                            log_next, n_last,
                                            tmp.posterior + (M - 1u) * P))
      goto fail;
  }

  for (size_t marker = M - 1u; marker > 0u; marker--)
  {
    size_t np1, np2, nc1, nc2, n_previous, n_current;
    const st_phase_pair_state *previous = st_phase_pair_state_get(
        pair, marker - 1u, &n_previous, NULL, NULL);
    const st_phase_pair_state *current =
        st_phase_pair_state_get(pair, marker, &n_current, NULL, NULL);
    double log_norm = -INFINITY;

    if ((NULL == previous) || (NULL == current) ||
        (0 != st_phase_factor_interval_maps(
                  pair, marker, prev_to_cur1, prev_to_cur2, cur_to_prev1,
                  cur_to_prev2, &np1, &np2, &nc1, &nc2)) ||
        (n_previous != np1 * np2) || (n_current != nc1 * nc2))
      goto fail;
    (void)previous;
    (void)current;

    for (size_t i = 0u; i < P; i++)
      log_back_cur[i] = -INFINITY;
    if (0 != st_phase_factor_backward_transition(
                 log_next, emission + marker * P, np1, np2, nc1, nc2,
                 prev_to_cur1, prev_to_cur2, q1[marker - 1u],
                 q2[marker - 1u], row_log, col_log, x_log, log_back_cur))
      goto fail;

    for (size_t i = 0u; i < n_previous; i++)
      log_norm = st_phase_hmm_logadd(log_norm, log_back_cur[i]);
    if (!(log_norm > -INFINITY))
      goto fail;

    for (size_t i = 0u; i < n_previous; i++)
    {
      if (log_back_cur[i] > -INFINITY)
        log_back_cur[i] -= log_norm;
      tmp.backward[(marker - 1u) * P + i] =
          (log_back_cur[i] > -INFINITY) ? exp(log_back_cur[i]) : 0.0;
    }

    if (0 != st_phase_hmm_posterior_log_row(
                 log_forward + (marker - 1u) * P, log_back_cur, n_previous,
                 tmp.posterior + (marker - 1u) * P))
      goto fail;
    {
      double *swap = log_next;
      log_next = log_back_cur;
      log_back_cur = swap;
    }
  }

  free(log_forward);
  free(log_prev);
  free(log_cur);
  free(vit_prev);
  free(vit_cur);
  free(log_next);
  free(log_back_cur);
  free(factor_mass);
  free(stage_score);
  free(stage_prev_b);
  free(x_log);
  free(back);
  free(side1);
  free(side2);
  free(row_log);
  free(col_log);
  free(prev_to_cur1);
  free(prev_to_cur2);
  free(cur_to_prev1);
  free(cur_to_prev2);
  st_phase_pair_hmm_result_free(out);
  *out = tmp;
  return 0;

impossible:
  tmp.log_likelihood = -INFINITY;
  tmp.viterbi_log_probability = -INFINITY;
  tmp.possible = 0;
  free(log_forward);
  free(log_prev);
  free(log_cur);
  free(vit_prev);
  free(vit_cur);
  free(log_next);
  free(log_back_cur);
  free(factor_mass);
  free(stage_score);
  free(stage_prev_b);
  free(x_log);
  free(back);
  free(side1);
  free(side2);
  free(row_log);
  free(col_log);
  free(prev_to_cur1);
  free(prev_to_cur2);
  free(cur_to_prev1);
  free(cur_to_prev2);
  st_phase_pair_hmm_result_free(out);
  *out = tmp;
  return 0;

fail:
  free(log_forward);
  free(log_prev);
  free(log_cur);
  free(vit_prev);
  free(vit_cur);
  free(log_next);
  free(log_back_cur);
  free(factor_mass);
  free(stage_score);
  free(stage_prev_b);
  free(x_log);
  free(back);
  free(side1);
  free(side2);
  free(row_log);
  free(col_log);
  free(prev_to_cur1);
  free(prev_to_cur2);
  free(cur_to_prev1);
  free(cur_to_prev2);
  st_phase_pair_hmm_result_free(&tmp);
  return -1;
}

static int st_phase_pair_hmm_viterbi_factorized_path(
    const st_phase_pair_state_set *pair, const double *initial_weight,
    const double *emission, const double *q1, const double *q2,
    size_t *viterbi_hap1, size_t *viterbi_hap2, int *possible)
{
  const size_t M = (NULL != pair) ? pair->n_marker : 0u;
  const size_t P = (NULL != pair) ? pair->max_pair_states : 0u;
  const size_t K = (NULL != pair) ? pair->max_side_states : 0u;
  size_t n_first = 0u, ncell = 0u;
  double max_initial = 0.0, scaled_initial_sum = 0.0;
  double *vit_prev = NULL, *vit_cur = NULL, *stage_score = NULL;
  size_t *stage_prev_b = NULL, *back = NULL, *p2c1 = NULL, *p2c2 = NULL,
         *c2p1 = NULL, *c2p2 = NULL;
  int status = -1;
  if ((NULL == possible) || (NULL == viterbi_hap1) ||
      (NULL == viterbi_hap2) || (NULL == initial_weight) ||
      (NULL == emission) || !st_phase_pair_hmm_set_valid(pair) ||
      ((M > 1u) && ((NULL == q1) || (NULL == q2))) ||
      (0 != st_phase_size_mul(M, P, &ncell)))
    return -1;
  *possible = 0;
  for (size_t m = 0u; m < M; m++)
  {
    viterbi_hap1[m] = SIZE_MAX;
    viterbi_hap2[m] = SIZE_MAX;
  }
  {
    const st_phase_pair_state *r =
        st_phase_pair_state_get(pair, 0u, &n_first, NULL, NULL);
    if ((NULL == r) || (0u == n_first) || (n_first > P))
      return -1;
  }
  for (size_t j = 0u; j < n_first; j++)
  {
    double w = initial_weight[j];
    if (!isfinite(w) || w < 0.0)
      return -1;
    if (w > max_initial)
      max_initial = w;
  }
  if (!(max_initial > 0.0))
    return -1;
  for (size_t j = 0u; j < n_first; j++)
    scaled_initial_sum += initial_weight[j] / max_initial;
  if (!(scaled_initial_sum > 0.0) || !isfinite(scaled_initial_sum))
    return -1;
  for (size_t m = 0u; m < M; m++)
  {
    size_t n = 0u;
    const st_phase_pair_state *r =
        st_phase_pair_state_get(pair, m, &n, NULL, NULL);
    if ((NULL == r) || (0u == n) || (n > P))
      return -1;
    for (size_t j = 0u; j < n; j++)
    {
      double e = emission[m * P + j];
      if (!isfinite(e) || e < 0.0 || e > 1.0)
        return -1;
    }
  }
  for (size_t m = 1u; m < M; m++)
    if (!isfinite(q1[m - 1u]) || !isfinite(q2[m - 1u]) ||
        q1[m - 1u] < 0.0 || q1[m - 1u] > 1.0 || q2[m - 1u] < 0.0 ||
        q2[m - 1u] > 1.0)
      return -1;
  vit_prev = st_phase_factor_malloc(P * sizeof(*vit_prev));
  vit_cur = st_phase_factor_malloc(P * sizeof(*vit_cur));
  stage_score = st_phase_factor_malloc(P * sizeof(*stage_score));
  stage_prev_b = st_phase_factor_malloc(P * sizeof(*stage_prev_b));
  back = st_phase_factor_malloc(ncell * sizeof(*back));
  p2c1 = st_phase_factor_malloc(K * sizeof(*p2c1));
  p2c2 = st_phase_factor_malloc(K * sizeof(*p2c2));
  c2p1 = st_phase_factor_malloc(K * sizeof(*c2p1));
  c2p2 = st_phase_factor_malloc(K * sizeof(*c2p2));
  if (!vit_prev || !vit_cur || !stage_score || !stage_prev_b || !back ||
      !p2c1 || !p2c2 || !c2p1 || !c2p2)
    goto cleanup;
  for (size_t z = 0u; z < ncell; z++)
    back[z] = SIZE_MAX;
  {
    double lsum = log(scaled_initial_sum), norm = -INFINITY;
    for (size_t j = 0u; j < P; j++)
      vit_prev[j] = -INFINITY;
    for (size_t j = 0u; j < n_first; j++)
    {
      double lw = (initial_weight[j] > 0.0)
                      ? log(initial_weight[j]) - log(max_initial) - lsum
                      : -INFINITY;
      double le = st_phase_hmm_log_probability(emission[j]);
      double x =
          ((lw > -INFINITY) && (le > -INFINITY)) ? lw + le : -INFINITY;
      vit_prev[j] = x;
      if (x > norm)
        norm = x;
    }
    if (!(norm > -INFINITY))
    {
      status = 0;
      goto cleanup;
    }
    for (size_t j = 0u; j < n_first; j++)
      if (vit_prev[j] > -INFINITY)
        vit_prev[j] -= norm;
  }
  for (size_t m = 1u; m < M; m++)
  {
    size_t np1, np2, nc1, nc2, np, nc;
    const st_phase_pair_state *pr = st_phase_pair_state_get(
                                  pair, m - 1u, &np, NULL, NULL),
                              *cr = st_phase_pair_state_get(pair, m, &nc,
                                                            NULL, NULL);
    double norm = -INFINITY;
    if (!pr || !cr ||
        0 != st_phase_factor_interval_maps(pair, m, p2c1, p2c2, c2p1, c2p2,
                                           &np1, &np2, &nc1, &nc2) ||
        np != np1 * np2 || nc != nc1 * nc2)
      goto cleanup;
    for (size_t j = 0u; j < P; j++)
      vit_cur[j] = -INFINITY;
    if (0 != st_phase_factor_viterbi_side2(vit_prev, np1, np2, nc2, p2c2,
                                           c2p2, q2[m - 1u], stage_score,
                                           stage_prev_b) ||
        0 != st_phase_factor_viterbi_side1(
                 stage_score, stage_prev_b, np1, np2, nc1, nc2, p2c1, c2p1,
                 q1[m - 1u], vit_cur, back + m * P))
      goto cleanup;
    for (size_t j = 0u; j < nc; j++)
    {
      double le = st_phase_hmm_log_probability(emission[m * P + j]);
      if ((vit_cur[j] > -INFINITY) && (le > -INFINITY))
        vit_cur[j] += le;
      else
        vit_cur[j] = -INFINITY;
      if (vit_cur[j] > norm)
        norm = vit_cur[j];
    }
    if (!(norm > -INFINITY))
    {
      status = 0;
      goto cleanup;
    }
    for (size_t j = 0u; j < nc; j++)
      if (vit_cur[j] > -INFINITY)
        vit_cur[j] -= norm;
    {
      double *z = vit_prev;
      vit_prev = vit_cur;
      vit_cur = z;
    }
  }
  {
    size_t nl = 0u, slot = SIZE_MAX;
    double best = -INFINITY;
    const st_phase_pair_state *r =
        st_phase_pair_state_get(pair, M - 1u, &nl, NULL, NULL);
    if (!r || 0u == nl || nl > P)
      goto cleanup;
    for (size_t j = 0u; j < nl; j++)
      if ((vit_prev[j] > best) ||
          ((vit_prev[j] == best) && ((SIZE_MAX == slot) || (j < slot))))
      {
        best = vit_prev[j];
        slot = j;
      }
    if (SIZE_MAX == slot)
      goto cleanup;
    for (size_t mm = M; mm > 0u; mm--)
    {
      size_t m = mm - 1u, n = 0u;
      const st_phase_pair_state *row =
          st_phase_pair_state_get(pair, m, &n, NULL, NULL);
      if (!row || slot >= n)
        goto cleanup;
      viterbi_hap1[m] = row[slot].hap1;
      viterbi_hap2[m] = row[slot].hap2;
      if (m > 0u)
      {
        slot = back[m * P + slot];
        if (SIZE_MAX == slot)
          goto cleanup;
      }
    }
  }
  *possible = 1;
  status = 0;
cleanup:
  free(vit_prev);
  free(vit_cur);
  free(stage_score);
  free(stage_prev_b);
  free(back);
  free(p2c1);
  free(p2c2);
  free(c2p1);
  free(c2p2);
  return status;
}

typedef struct
{
  size_t individual;
  size_t chromosome;
} st_phase_dataset_request;

typedef enum
{
  ST_PHASE_DATASET_RUN_OK = 0,
  ST_PHASE_DATASET_RUN_INVALID = -1,
  ST_PHASE_DATASET_RUN_PHASING_FAILED = -2,
  ST_PHASE_DATASET_RUN_MERGE_FAILED = -3,
  ST_PHASE_DATASET_RUN_AUDIT_FAILED = -4
} st_phase_dataset_run_status;

typedef struct
{
  size_t n_request;
  size_t n_possible;
  size_t n_impossible;
  size_t n_requested_markers;
  int audit_status;
  st_phase_dataset_stage stage;
  st_phase_stage_audit audit;
} st_phase_dataset_run_result;

#define st_phase_dataset_run_calloc calloc

static int st_phase_dataset_requests_valid(
    const st_phase_data *d, const st_phase_dataset_request *request,
    size_t n_request, size_t *requested_markers)
{
  unsigned char *seen = NULL;
  size_t total = 0u, max_request;
  int valid = 0;

  if (!st_phase_data_adapter_valid(d) || (NULL == request) ||
      (0u == n_request) ||
      (0 != st_phase_size_mul(d->n_ind, d->n_chrom, &max_request)) ||
      (0u == max_request) || (n_request > max_request) ||
      (NULL == requested_markers))
    return 0;

  seen = st_phase_dataset_run_calloc(max_request, sizeof(*seen));
  if (NULL == seen)
    return 0;

  for (size_t r = 0u; r < n_request; r++)
  {
    size_t len, slot;

    if ((request[r].individual >= d->n_ind) ||
        (request[r].chromosome >= d->n_chrom))
      goto cleanup;

    slot = request[r].individual * d->n_chrom + request[r].chromosome;
    if (seen[slot])
      goto cleanup;
    seen[slot] = 1u;

    len = d->chrom_length[request[r].chromosome];
    if ((0u == len) || (total > SIZE_MAX - len))
      goto cleanup;
    total += len;
  }

  *requested_markers = total;
  valid = 1;

cleanup:
  free(seen);
  return valid;
}

typedef enum
{
  ST_PHASE_PIPELINE_OK = 0,
  ST_PHASE_PIPELINE_INVALID = -1,
  ST_PHASE_PIPELINE_DATA_FAILED = -2,
  ST_PHASE_PIPELINE_RUN_FAILED = -3,
  ST_PHASE_PIPELINE_COMMIT_FAILED = -4
} st_phase_pipeline_status;

typedef struct
{
  int pipeline_status;
  int dataset_run_status;
  int audit_status;
  int commit_requested;
  int commit_attempted;
  int commit_status;
  int committed;

  size_t n_ind;
  size_t n_marker;
  size_t n_chrom;
  size_t n_genotype;

  size_t n_request;
  size_t n_possible;
  size_t n_impossible;
  size_t n_requested_markers;

  size_t n_covered;
  size_t n_uncovered;
  size_t n_written;
  size_t n_covered_unwritten;
  size_t n_changed_genotypes;

  size_t n_missing;
  size_t n_homozygous;
  size_t n_heterozygous;
  size_t n_covered_missing;
  size_t n_covered_homozygous;
  size_t n_covered_heterozygous;
  size_t n_covered_heterozygous_unresolved;
  size_t n_written_same_source_order;
  size_t n_written_swapped_source_order;

  int audit_complete;
  int audit_source_bytes_current;
  int audit_workspace_current;
} st_phase_pipeline_diagnostics;

static void
st_phase_pipeline_diagnostics_init(st_phase_pipeline_diagnostics *x)
{
  if (NULL != x)
  {
    memset(x, 0, sizeof(*x));
    x->pipeline_status = ST_PHASE_PIPELINE_INVALID;
    x->dataset_run_status = ST_PHASE_DATASET_RUN_INVALID;
    x->audit_status = ST_PHASE_AUDIT_INVALID;
    x->commit_status = ST_PHASE_COMMIT_INVALID;
  }
}

static size_t st_phase_dataset_stage_changed_genotypes(
    const st_phase_dataset_stage *stage)
{
  size_t changed = 0u;

  if ((NULL == stage) || (NULL == stage->base1) ||
      (NULL == stage->base2) || (NULL == stage->staged1) ||
      (NULL == stage->staged2))
    return 0u;

  for (size_t z = 0u; z < stage->n_genotype; z++)
    if ((stage->base1[z] != stage->staged1[z]) ||
        (stage->base2[z] != stage->staged2[z]))
      changed++;

  return changed;
}

static int st_phase_pipeline_parameters_valid(
    const gs_varset_type *gv, const st_phase_dataset_request *request,
    size_t n_request, const uint8_t *reference_hap,
    size_t reference_n_marker, size_t n_ref_hap, size_t max_states,
    const double *initial_pair_weight, size_t initial_pair_count,
    double error_rate, double break_rate1_per_morgan,
    double break_rate2_per_morgan, double confidence_threshold,
    int commit_requested)
{
  size_t ref_cells, pair_count;
  double prior_sum = 0.0;

  if ((NULL == gv) || (NULL == request) || (0u == n_request) ||
      (NULL == reference_hap) || (NULL == initial_pair_weight) ||
      (0u == reference_n_marker) || (0u == n_ref_hap) ||
      (0u == max_states) ||
      (commit_requested != 0 && commit_requested != 1) ||
      !st_phase_hmm_error_valid(error_rate) ||
      !isfinite(break_rate1_per_morgan) ||
      !isfinite(break_rate2_per_morgan) ||
      (break_rate1_per_morgan < 0.0) || (break_rate2_per_morgan < 0.0) ||
      !st_phase_hard_threshold_valid(confidence_threshold) ||
      (reference_n_marker != (size_t)gv->NoMar) ||
      (0 !=
       st_phase_size_mul(reference_n_marker, n_ref_hap, &ref_cells)) ||
      (0u == ref_cells) ||
      (0 != st_phase_size_mul(n_ref_hap, n_ref_hap, &pair_count)) ||
      (initial_pair_count != pair_count))
    return 0;

  for (size_t z = 0u; z < ref_cells; z++)
    if (reference_hap[z] > 1u)
      return 0;

  for (size_t z = 0u; z < pair_count; z++)
  {
    double w = initial_pair_weight[z];
    if (!isfinite(w) || (w < 0.0))
      return 0;
    prior_sum += w;
    if (!isfinite(prior_sum))
      return 0;
  }

  return prior_sum > 0.0;
}

#define st_phase_relative_malloc malloc
#define st_phase_relative_calloc calloc

typedef struct
{
  size_t n_link;
  size_t *left_marker;
  size_t *right_marker;
  double *same;
  double *switch_probability;
  int possible;
} st_phase_relative_orientation_result;

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

static int st_phase_r_scalar_size(SEXP x, size_t *out)
{
  double v;
  if ((NULL == out) || (XLENGTH(x) != 1))
    return 0;
  if (INTSXP == TYPEOF(x))
  {
    int z = INTEGER(x)[0];
    if ((NA_INTEGER == z) || (z <= 0))
      return 0;
    *out = (size_t)z;
    return 1;
  }
  if (REALSXP != TYPEOF(x))
    return 0;
  v = REAL(x)[0];
  if (!isfinite(v) || (v < 1.0) || (floor(v) != v) ||
      (v > (double)SIZE_MAX))
    return 0;
  *out = (size_t)v;
  return 1;
}

static int st_phase_r_scalar_double(SEXP x, double *out)
{
  double v;
  if ((NULL == out) || (XLENGTH(x) != 1))
    return 0;
  if (REALSXP == TYPEOF(x))
    v = REAL(x)[0];
  else if (INTSXP == TYPEOF(x))
  {
    int z = INTEGER(x)[0];
    if (NA_INTEGER == z)
      return 0;
    v = (double)z;
  }
  else
    return 0;
  if (!isfinite(v))
    return 0;
  *out = v;
  return 1;
}

static int st_phase_r_scalar_logical(SEXP x, int *out)
{
  int v;
  if ((NULL == out) || (LGLSXP != TYPEOF(x)) || (XLENGTH(x) != 1))
    return 0;
  v = LOGICAL(x)[0];
  if (NA_LOGICAL == v)
    return 0;
  *out = (0 != v);
  return 1;
}

static int st_phase_r_matrix_dims(SEXP x, size_t *nr, size_t *nc)
{
  SEXP dim = getAttrib(x, R_DimSymbol);
  if ((NULL == nr) || (NULL == nc) || (INTSXP != TYPEOF(dim)) ||
      (XLENGTH(dim) != 2))
    return 0;
  if ((INTEGER(dim)[0] <= 0) || (INTEGER(dim)[1] <= 0))
    return 0;
  *nr = (size_t)INTEGER(dim)[0];
  *nc = (size_t)INTEGER(dim)[1];
  return 1;
}

static int st_phase_r_build_prior(SEXP prior, size_t R, double **out)
{
  size_t nr, nc, ncell;
  double *x = NULL;
  double sum = 0.0;

  if ((NULL == out) || !st_phase_r_matrix_dims(prior, &nr, &nc) ||
      (nr != R) || (nc != R) ||
      ((REALSXP != TYPEOF(prior)) && (INTSXP != TYPEOF(prior))) ||
      (0 != st_phase_size_mul(R, R, &ncell)))
    return -1;
  x = malloc(ncell * sizeof(*x));
  if (NULL == x)
    return -1;

  for (size_t i = 0; i < R; i++)
    for (size_t j = 0; j < R; j++)
    {
      R_xlen_t idx = (R_xlen_t)i + (R_xlen_t)R * (R_xlen_t)j;
      double v;
      if (REALSXP == TYPEOF(prior))
        v = REAL(prior)[idx];
      else
      {
        int z = INTEGER(prior)[idx];
        if (NA_INTEGER == z)
        {
          free(x);
          return -1;
        }
        v = (double)z;
      }
      if (!isfinite(v) || (v < 0.0))
      {
        free(x);
        return -1;
      }
      x[i * R + j] = v;
      sum += v;
      if (!isfinite(sum))
      {
        free(x);
        return -1;
      }
    }
  if (!(sum > 0.0))
  {
    free(x);
    return -1;
  }
  *out = x;
  return 0;
}

static int st_phase_r_individuals(const gs_varset_type *gv, SEXP x,
                                  size_t **index_out, size_t *n_out)
{
  size_t *index = NULL;
  unsigned char *seen = NULL;
  size_t n;
  int status = -1;

  if ((NULL == gv) || (NULL == index_out) || (NULL == n_out) ||
      (0 == gv->NoInd))
    return -1;

  if (R_NilValue == x)
  {
    n = (size_t)gv->NoInd;
    index = malloc(n * sizeof(*index));
    if (NULL == index)
      return -1;
    for (size_t i = 0; i < n; i++)
      index[i] = i;
    *index_out = index;
    *n_out = n;
    return 0;
  }

  n = (size_t)XLENGTH(x);
  if (0u == n)
    return -1;
  index = malloc(n * sizeof(*index));
  seen = calloc((size_t)gv->NoInd, sizeof(*seen));
  if ((NULL == index) || (NULL == seen))
    goto cleanup;

  if (STRSXP == TYPEOF(x))
  {
    st_phase_r_name_index *names =
        malloc((size_t)gv->NoInd * sizeof(*names));
    if (NULL == names)
      goto cleanup;
    for (size_t i = 0; i < (size_t)gv->NoInd; i++)
    {
      names[i].name = gv->IndName[i];
      names[i].index = i;
    }
    qsort(names, (size_t)gv->NoInd, sizeof(*names),
          st_phase_r_name_index_cmp);
    for (size_t k = 1; k < (size_t)gv->NoInd; k++)
      if (0 == strcmp(names[k - 1].name, names[k].name))
      {
        free(names);
        goto cleanup;
      }
    for (size_t k = 0; k < n; k++)
    {
      SEXP s = STRING_ELT(x, (R_xlen_t)k);
      const char *want;
      size_t lo = 0, hi = (size_t)gv->NoInd;
      if (NA_STRING == s)
      {
        free(names);
        goto cleanup;
      }
      want = CHAR(s);
      while (lo < hi)
      {
        size_t mid = lo + (hi - lo) / 2;
        int cmp = strcmp(names[mid].name, want);
        if (cmp < 0)
          lo = mid + 1;
        else
          hi = mid;
      }
      if ((lo >= (size_t)gv->NoInd) || (0 != strcmp(names[lo].name, want)))
      {
        free(names);
        goto cleanup;
      }
      index[k] = names[lo].index;
      if (seen[index[k]])
      {
        free(names);
        goto cleanup;
      }
      seen[index[k]] = 1u;
    }
    free(names);
  }
  else if ((INTSXP == TYPEOF(x)) || (REALSXP == TYPEOF(x)))
  {
    for (size_t k = 0; k < n; k++)
    {
      double v;
      if (INTSXP == TYPEOF(x))
      {
        int z = INTEGER(x)[k];
        if (NA_INTEGER == z)
          goto cleanup;
        v = (double)z;
      }
      else
        v = REAL(x)[k];
      if (!isfinite(v) || (floor(v) != v) || (v < 1.0) ||
          (v > (double)gv->NoInd))
        goto cleanup;
      index[k] = (size_t)v - 1u;
      if (seen[index[k]])
        goto cleanup;
      seen[index[k]] = 1u;
    }
  }
  else
    goto cleanup;

  *index_out = index;
  *n_out = n;
  index = NULL;
  status = 0;
cleanup:
  free(index);
  free(seen);
  return status;
}

static int st_phase_r_chromosomes(const st_phase_data *d, SEXP x,
                                  size_t **index_out, size_t *n_out)
{
  size_t *index = NULL;
  unsigned char *seen = NULL;
  size_t n;
  int status = -1;

  if ((NULL == d) || (NULL == index_out) || (NULL == n_out) ||
      (0u == d->n_chrom))
    return -1;
  if (R_NilValue == x)
  {
    n = d->n_chrom;
    index = malloc(n * sizeof(*index));
    if (NULL == index)
      return -1;
    for (size_t c = 0; c < n; c++)
      index[c] = c;
    *index_out = index;
    *n_out = n;
    return 0;
  }
  if ((INTSXP != TYPEOF(x)) && (REALSXP != TYPEOF(x)))
    return -1;
  n = (size_t)XLENGTH(x);
  if (0u == n)
    return -1;
  index = malloc(n * sizeof(*index));
  seen = calloc(d->n_chrom, sizeof(*seen));
  if ((NULL == index) || (NULL == seen))
    goto cleanup;

  for (size_t k = 0; k < n; k++)
  {
    double v;
    int label;
    size_t found = SIZE_MAX;
    if (INTSXP == TYPEOF(x))
    {
      int z = INTEGER(x)[k];
      if (NA_INTEGER == z)
        goto cleanup;
      v = (double)z;
    }
    else
      v = REAL(x)[k];
    if (!isfinite(v) || (floor(v) != v) || (v < (double)INT_MIN) ||
        (v > (double)INT_MAX))
      goto cleanup;
    label = (int)v;
    for (size_t c = 0; c < d->n_chrom; c++)
      if (d->chrom_label[d->chrom_start[c]] == label)
      {
        found = c;
        break;
      }
    if ((SIZE_MAX == found) || seen[found])
      goto cleanup;
    seen[found] = 1u;
    index[k] = found;
  }
  *index_out = index;
  *n_out = n;
  index = NULL;
  status = 0;
cleanup:
  free(index);
  free(seen);
  return status;
}

static const char *st_phase_r_pipeline_status_name(int x)
{
  switch (x)
  {
  case ST_PHASE_PIPELINE_OK:
    return "ok";
  case ST_PHASE_PIPELINE_INVALID:
    return "invalid";
  case ST_PHASE_PIPELINE_DATA_FAILED:
    return "data_failed";
  case ST_PHASE_PIPELINE_RUN_FAILED:
    return "run_failed";
  case ST_PHASE_PIPELINE_COMMIT_FAILED:
    return "commit_failed";
  default:
    return "unknown";
  }
}
static const char *st_phase_r_run_status_name(int x)
{
  switch (x)
  {
  case ST_PHASE_DATASET_RUN_OK:
    return "ok";
  case ST_PHASE_DATASET_RUN_INVALID:
    return "invalid";
  case ST_PHASE_DATASET_RUN_PHASING_FAILED:
    return "phasing_failed";
  case ST_PHASE_DATASET_RUN_MERGE_FAILED:
    return "merge_failed";
  case ST_PHASE_DATASET_RUN_AUDIT_FAILED:
    return "audit_failed";
  default:
    return "unknown";
  }
}
static const char *st_phase_r_audit_status_name(int x)
{
  switch (x)
  {
  case ST_PHASE_AUDIT_OK:
    return "ok";
  case ST_PHASE_AUDIT_INVALID:
    return "invalid";
  case ST_PHASE_AUDIT_STALE:
    return "stale";
  case ST_PHASE_AUDIT_INCOMPLETE:
    return "incomplete";
  case ST_PHASE_AUDIT_CORRUPT:
    return "corrupt";
  default:
    return "unknown";
  }
}
static const char *st_phase_r_commit_status_name(int x)
{
  switch (x)
  {
  case ST_PHASE_COMMIT_OK:
    return "ok";
  case ST_PHASE_COMMIT_INVALID:
    return "invalid";
  case ST_PHASE_COMMIT_STALE:
    return "stale";
  case ST_PHASE_COMMIT_INCOMPLETE:
    return "incomplete";
  case ST_PHASE_COMMIT_CORRUPT:
    return "corrupt";
  case ST_PHASE_COMMIT_ROLLED_BACK:
    return "rolled_back";
  default:
    return "unknown";
  }
}

static SEXP st_phase_r_count(size_t x)
{
  return ScalarReal((double)x);
}

static SEXP st_phase_r_diagnostics(const st_phase_pipeline_diagnostics *d)
{
  static const char *names[] = {"success",
                                "pipeline.status",
                                "pipeline.status.name",
                                "dataset.run.status",
                                "dataset.run.status.name",
                                "audit.status",
                                "audit.status.name",
                                "commit.requested",
                                "commit.attempted",
                                "commit.status",
                                "commit.status.name",
                                "committed",
                                "n.ind",
                                "n.marker",
                                "n.chrom",
                                "n.genotype",
                                "n.request",
                                "n.possible",
                                "n.impossible",
                                "n.requested.markers",
                                "n.covered",
                                "n.uncovered",
                                "n.written",
                                "n.covered.unwritten",
                                "n.changed.genotypes",
                                "n.missing",
                                "n.homozygous",
                                "n.heterozygous",
                                "n.covered.missing",
                                "n.covered.homozygous",
                                "n.covered.heterozygous",
                                "n.covered.heterozygous.unresolved",
                                "n.written.same.source.order",
                                "n.written.swapped.source.order",
                                "audit.complete",
                                "audit.source.bytes.current",
                                "audit.workspace.current"};
  enum
  {
    N = (int)(sizeof(names) / sizeof(names[0]))
  };
  SEXP ans, nms;
  int k = 0;
  PROTECT(ans = allocVector(VECSXP, N));
  PROTECT(nms = allocVector(STRSXP, N));
#define ADD(value)                                                        \
  do                                                                      \
  {                                                                       \
    SET_VECTOR_ELT(ans, k, (value));                                      \
    SET_STRING_ELT(nms, k, mkChar(names[k]));                             \
    k++;                                                                  \
  } while (0)
  ADD(ScalarLogical(ST_PHASE_PIPELINE_OK == d->pipeline_status));
  ADD(ScalarInteger(d->pipeline_status));
  ADD(mkString(st_phase_r_pipeline_status_name(d->pipeline_status)));
  ADD(ScalarInteger(d->dataset_run_status));
  ADD(mkString(st_phase_r_run_status_name(d->dataset_run_status)));
  ADD(ScalarInteger(d->audit_status));
  ADD(mkString(st_phase_r_audit_status_name(d->audit_status)));
  ADD(ScalarLogical(d->commit_requested));
  ADD(ScalarLogical(d->commit_attempted));
  ADD(ScalarInteger(d->commit_status));
  ADD(mkString(st_phase_r_commit_status_name(d->commit_status)));
  ADD(ScalarLogical(d->committed));
  ADD(st_phase_r_count(d->n_ind));
  ADD(st_phase_r_count(d->n_marker));
  ADD(st_phase_r_count(d->n_chrom));
  ADD(st_phase_r_count(d->n_genotype));
  ADD(st_phase_r_count(d->n_request));
  ADD(st_phase_r_count(d->n_possible));
  ADD(st_phase_r_count(d->n_impossible));
  ADD(st_phase_r_count(d->n_requested_markers));
  ADD(st_phase_r_count(d->n_covered));
  ADD(st_phase_r_count(d->n_uncovered));
  ADD(st_phase_r_count(d->n_written));
  ADD(st_phase_r_count(d->n_covered_unwritten));
  ADD(st_phase_r_count(d->n_changed_genotypes));
  ADD(st_phase_r_count(d->n_missing));
  ADD(st_phase_r_count(d->n_homozygous));
  ADD(st_phase_r_count(d->n_heterozygous));
  ADD(st_phase_r_count(d->n_covered_missing));
  ADD(st_phase_r_count(d->n_covered_homozygous));
  ADD(st_phase_r_count(d->n_covered_heterozygous));
  ADD(st_phase_r_count(d->n_covered_heterozygous_unresolved));
  ADD(st_phase_r_count(d->n_written_same_source_order));
  ADD(st_phase_r_count(d->n_written_swapped_source_order));
  ADD(ScalarLogical(d->audit_complete));
  ADD(ScalarLogical(d->audit_source_bytes_current));
  ADD(ScalarLogical(d->audit_workspace_current));
#undef ADD
  setAttrib(ans, R_NamesSymbol, nms);
  UNPROTECT(2);
  return ans;
}

typedef struct
{
  size_t n_ind;
  size_t n_marker;
  size_t n_hap;
  size_t n_heterozygous;
  size_t n_orientation_informative;
  size_t n_orientation_tie;
  size_t n_missing_proxy;
  uint8_t *hap;
} st_phase_cohort_initial;

static void st_phase_cohort_initial_init(st_phase_cohort_initial *x)
{
  if (NULL != x)
    memset(x, 0, sizeof(*x));
}

static void st_phase_cohort_initial_free(st_phase_cohort_initial *x)
{
  if (NULL == x)
    return;
  free(x->hap);
  st_phase_cohort_initial_init(x);
}

#define st_phase_cohort_malloc malloc

static int st_phase_cohort_genotype_valid(uint8_t a, uint8_t b)
{
  if ((ST_PHASE_STATE_MISSING == a) || (ST_PHASE_STATE_MISSING == b))
    return ((ST_PHASE_STATE_MISSING == a) &&
            (ST_PHASE_STATE_MISSING == b));
  return ((a <= 1u) && (b <= 1u));
}

static int st_phase_cohort_genotype_contains(uint8_t a, uint8_t b,
                                             uint8_t allele)
{
  if (!st_phase_cohort_genotype_valid(a, b) ||
      (ST_PHASE_STATE_MISSING == a) || (allele > 1u))
    return 0;
  return ((a == allele) || (b == allele));
}

static int st_phase_cohort_pbwt_prefix_support(
    const st_phase_pbwt *p, const uint8_t *g1, const uint8_t *g2,
    size_t n_ind, uint8_t allele, size_t *best)
{
  size_t H;
  size_t last1 = SIZE_MAX, last2 = SIZE_MAX;
  size_t maxd1 = 0u, maxd2 = 0u;

  if ((NULL == p) || (NULL == g1) || (NULL == g2) || (NULL == best) ||
      (0u == n_ind) || (allele > 1u) ||
      (0 != st_phase_size_mul(n_ind, 2u, &H)) || (p->n_hap != H) ||
      (NULL == p->order) || (NULL == p->divergence))
    return -1;

  for (size_t h = 0; h < H; h++)
    best[h] = 0u;

  for (size_t pos = 0; pos < H; pos++)
  {
    size_t h = p->order[pos];
    size_t target_ind;
    size_t candidate = SIZE_MAX;
    size_t candidate_maxd = 0u;

    if (h >= H)
      return -1;
    target_ind = h / 2u;

    if (SIZE_MAX != last1)
    {
      if (p->divergence[pos] > maxd1)
        maxd1 = p->divergence[pos];
      if (SIZE_MAX != last2 && p->divergence[pos] > maxd2)
        maxd2 = p->divergence[pos];

      if ((last1 / 2u) != target_ind)
      {
        candidate = last1;
        candidate_maxd = maxd1;
      }
      else if ((SIZE_MAX != last2) && ((last2 / 2u) != target_ind))
      {
        candidate = last2;
        candidate_maxd = maxd2;
      }

      if (SIZE_MAX != candidate)
      {
        if (candidate_maxd > p->position)
          return -1;
        best[h] = p->position - candidate_maxd;
      }
    }

    if (st_phase_cohort_genotype_contains(g1[target_ind], g2[target_ind],
                                          allele))
    {
      last2 = last1;
      maxd2 = maxd1;
      last1 = h;
      maxd1 = 0u;
    }
  }

  last1 = SIZE_MAX;
  last2 = SIZE_MAX;
  maxd1 = 0u;
  maxd2 = 0u;
  for (size_t pos = H; pos-- > 0u;)
  {
    size_t h = p->order[pos];
    size_t target_ind;
    size_t candidate = SIZE_MAX;
    size_t candidate_maxd = 0u;

    if (h >= H)
      return -1;
    target_ind = h / 2u;

    if (SIZE_MAX != last1)
    {
      size_t d = p->divergence[pos + 1u];
      if (d > maxd1)
        maxd1 = d;
      if (SIZE_MAX != last2 && d > maxd2)
        maxd2 = d;

      if ((last1 / 2u) != target_ind)
      {
        candidate = last1;
        candidate_maxd = maxd1;
      }
      else if ((SIZE_MAX != last2) && ((last2 / 2u) != target_ind))
      {
        candidate = last2;
        candidate_maxd = maxd2;
      }

      if (SIZE_MAX != candidate)
      {
        size_t length;
        if (candidate_maxd > p->position)
          return -1;
        length = p->position - candidate_maxd;
        if (length > best[h])
          best[h] = length;
      }
    }

    if (st_phase_cohort_genotype_contains(g1[target_ind], g2[target_ind],
                                          allele))
    {
      last2 = last1;
      maxd2 = maxd1;
      last1 = h;
      maxd1 = 0u;
    }
  }

  return 0;
}

static uint8_t st_phase_cohort_missing_proxy(const uint8_t *g1,
                                             const uint8_t *g2,
                                             size_t n_ind)
{
  size_t count0 = 0u, count1 = 0u;

  for (size_t i = 0; i < n_ind; i++)
    if ((ST_PHASE_STATE_MISSING != g1[i]) &&
        (ST_PHASE_STATE_MISSING != g2[i]))
    {
      count0 += (0u == g1[i]) + (0u == g2[i]);
      count1 += (1u == g1[i]) + (1u == g2[i]);
    }

  return (count1 > count0) ? 1u : 0u;
}

static int st_phase_cohort_initialize_binary(const uint8_t *geno1,
                                             const uint8_t *geno2,
                                             size_t n_marker, size_t n_ind,
                                             st_phase_cohort_initial *out)
{
  st_phase_cohort_initial tmp;
  st_phase_pbwt pbwt;
  uint8_t *marker_hap = NULL;
  size_t *best0 = NULL, *best1 = NULL;
  size_t H, ncell;

  st_phase_cohort_initial_init(&tmp);
  st_phase_pbwt_init_empty(&pbwt);

  if ((NULL == out) || (NULL == geno1) || (NULL == geno2) ||
      (0u == n_marker) || (n_ind < 2u) ||
      (0 != st_phase_size_mul(n_ind, 2u, &H)) ||
      (0 != st_phase_size_mul(n_marker, n_ind, &ncell)) ||
      (H > SIZE_MAX / n_marker))
    return -1;

  for (size_t z = 0; z < ncell; z++)
    if (!st_phase_cohort_genotype_valid(geno1[z], geno2[z]))
      return -1;

  if ((n_marker * H) > SIZE_MAX / sizeof(*tmp.hap))
    return -1;

  tmp.hap = st_phase_cohort_malloc(n_marker * H * sizeof(*tmp.hap));
  marker_hap = st_phase_cohort_malloc(H * sizeof(*marker_hap));
  best0 = st_phase_cohort_malloc(H * sizeof(*best0));
  best1 = st_phase_cohort_malloc(H * sizeof(*best1));
  if ((NULL == tmp.hap) || (NULL == marker_hap) || (NULL == best0) ||
      (NULL == best1))
    goto fail;

  tmp.n_ind = n_ind;
  tmp.n_marker = n_marker;
  tmp.n_hap = H;

  if (0 != st_phase_pbwt_init(&pbwt, H))
    goto fail;

  for (size_t marker = 0; marker < n_marker; marker++)
  {
    const uint8_t *g1 = geno1 + marker * n_ind;
    const uint8_t *g2 = geno2 + marker * n_ind;
    uint8_t missing_proxy = st_phase_cohort_missing_proxy(g1, g2, n_ind);

    if ((0 != st_phase_cohort_pbwt_prefix_support(&pbwt, g1, g2, n_ind, 0u,
                                                  best0)) ||
        (0 != st_phase_cohort_pbwt_prefix_support(&pbwt, g1, g2, n_ind, 1u,
                                                  best1)))
      goto fail;

    for (size_t i = 0; i < n_ind; i++)
    {
      size_t h0 = 2u * i, h1 = h0 + 1u;
      uint8_t a = g1[i], b = g2[i];

      if (ST_PHASE_STATE_MISSING == a)
      {
        marker_hap[h0] = missing_proxy;
        marker_hap[h1] = missing_proxy;
        tmp.n_missing_proxy++;
      }
      else if (a == b)
      {
        marker_hap[h0] = a;
        marker_hap[h1] = a;
      }
      else
      {
        size_t score01, score10;

        tmp.n_heterozygous++;
        score01 = best0[h0] + best1[h1];
        score10 = best1[h0] + best0[h1];

        if (score10 > score01)
        {
          marker_hap[h0] = 1u;
          marker_hap[h1] = 0u;
          tmp.n_orientation_informative++;
        }
        else
        {
          marker_hap[h0] = 0u;
          marker_hap[h1] = 1u;
          if (score01 > score10)
            tmp.n_orientation_informative++;
          else
            tmp.n_orientation_tie++;
        }
      }
    }

    memcpy(tmp.hap + marker * H, marker_hap, H * sizeof(*marker_hap));
    if (0 != st_phase_pbwt_update(&pbwt, marker_hap))
      goto fail;
  }

  st_phase_pbwt_free(&pbwt);
  free(marker_hap);
  free(best0);
  free(best1);
  st_phase_cohort_initial_free(out);
  *out = tmp;
  return 0;

fail:
  st_phase_pbwt_free(&pbwt);
  free(marker_hap);
  free(best0);
  free(best1);
  st_phase_cohort_initial_free(&tmp);
  return -1;
}

typedef struct
{
  gs_UNSIGNED source_marker;
  int chrom_label;
  double pos_cm;
} st_phase_cohort_snapshot_marker;

typedef struct
{
  size_t start;
  size_t length;
  int chrom_label;
  size_t n_heterozygous;
  size_t n_orientation_informative;
  size_t n_orientation_tie;
  size_t n_missing_proxy;
} st_phase_cohort_snapshot_chrom;

typedef struct
{
  size_t n_ind;
  size_t n_marker;
  size_t n_chrom;
  size_t n_hap;

  size_t n_heterozygous;
  size_t n_orientation_informative;
  size_t n_orientation_tie;
  size_t n_missing_proxy;

  uint8_t *hap;

  uint8_t *observed_missing;

  st_phase_cohort_snapshot_marker *marker;
  st_phase_cohort_snapshot_chrom *chrom;
} st_phase_cohort_snapshot;

static void st_phase_cohort_snapshot_init(st_phase_cohort_snapshot *s)
{
  if (NULL != s)
    memset(s, 0, sizeof(*s));
}

static void st_phase_cohort_snapshot_free(st_phase_cohort_snapshot *s)
{
  if (NULL == s)
    return;

  free(s->hap);
  free(s->observed_missing);
  free(s->marker);
  free(s->chrom);
  st_phase_cohort_snapshot_init(s);
}

static int st_phase_cohort_snapshot_source_valid(const st_phase_data *d)
{
  size_t expected_start = 0u;

  if ((NULL == d) || (d->n_ind < 2u) || (0u == d->n_marker) ||
      (0u == d->n_chrom) || (NULL == d->source_marker) ||
      (NULL == d->chrom_label) || (NULL == d->pos_cm) ||
      (NULL == d->chrom_start) || (NULL == d->chrom_length) ||
      (NULL == d->n_alleles) || (NULL == d->geno1) || (NULL == d->geno2))
    return 0;

  for (size_t c = 0; c < d->n_chrom; c++)
  {
    size_t start = d->chrom_start[c];
    size_t length = d->chrom_length[c];
    int label;

    if ((0u == length) || (start != expected_start) ||
        (start >= d->n_marker) || (length > d->n_marker - start))
      return 0;

    label = d->chrom_label[start];
    if (label < 1)
      return 0;

    for (size_t m = start; m < start + length; m++)
    {
      if ((d->chrom_label[m] != label) || !isfinite(d->pos_cm[m]) ||
          (d->pos_cm[m] < 0.0) ||
          ((m > start) && (d->pos_cm[m] < d->pos_cm[m - 1u])))
        return 0;
    }

    expected_start = start + length;
  }

  if (expected_start != d->n_marker)
    return 0;

  for (size_t m = 0; m < d->n_marker; m++)
  {
    if (d->n_alleles[m] > 2u)
      return 0;

    for (size_t i = 0; i < d->n_ind; i++)
    {
      size_t z = m * d->n_ind + i;
      uint8_t a = d->geno1[z], b = d->geno2[z];

      if (!st_phase_cohort_genotype_valid(a, b))
        return 0;
    }
  }

  return 1;
}

static int
st_phase_cohort_snapshot_valid(const st_phase_data *d,
                               const st_phase_cohort_snapshot *s)
{
  size_t H, n_genotype, n_hap_cell;
  size_t total_het = 0u, total_info = 0u, total_tie = 0u,
         total_missing = 0u;
  size_t expected_start = 0u;

  if (!st_phase_cohort_snapshot_source_valid(d) || (NULL == s) ||
      (0 != st_phase_size_mul(d->n_ind, 2u, &H)) ||
      (0 != st_phase_size_mul(d->n_marker, d->n_ind, &n_genotype)) ||
      (0 != st_phase_size_mul(d->n_marker, H, &n_hap_cell)) ||
      (s->n_ind != d->n_ind) || (s->n_marker != d->n_marker) ||
      (s->n_chrom != d->n_chrom) || (s->n_hap != H) || (NULL == s->hap) ||
      (NULL == s->observed_missing) || (NULL == s->marker) ||
      (NULL == s->chrom))
    return 0;

  (void)n_genotype;
  (void)n_hap_cell;

  for (size_t m = 0; m < d->n_marker; m++)
  {
    if ((s->marker[m].source_marker != d->source_marker[m]) ||
        (s->marker[m].chrom_label != d->chrom_label[m]) ||
        (s->marker[m].pos_cm != d->pos_cm[m]))
      return 0;

    for (size_t i = 0; i < d->n_ind; i++)
    {
      size_t z = m * d->n_ind + i;
      size_t h = m * H + 2u * i;
      uint8_t a = d->geno1[z], b = d->geno2[z];
      uint8_t x = s->hap[h], y = s->hap[h + 1u];
      int missing = (ST_PHASE_STATE_MISSING == a);

      if ((x > 1u) || (y > 1u) ||
          (s->observed_missing[z] != (uint8_t)(missing ? 1u : 0u)))
        return 0;

      if (missing)
        total_missing++;
      else
      {
        if (!(((x == a) && (y == b)) || ((x == b) && (y == a))))
          return 0;
        if (a != b)
          total_het++;
      }
    }
  }

  for (size_t c = 0; c < s->n_chrom; c++)
  {
    const st_phase_cohort_snapshot_chrom *sc = s->chrom + c;
    size_t start = d->chrom_start[c], length = d->chrom_length[c];
    size_t chrom_het = 0u, chrom_missing = 0u;

    for (size_t m = start; m < start + length; m++)
      for (size_t i = 0; i < d->n_ind; i++)
      {
        size_t z = m * d->n_ind + i;
        if (ST_PHASE_STATE_MISSING == d->geno1[z])
          chrom_missing++;
        else if (d->geno1[z] != d->geno2[z])
          chrom_het++;
      }

    if ((sc->start != start) || (sc->length != length) ||
        (sc->start != expected_start) ||
        (sc->chrom_label != d->chrom_label[start]) ||
        (sc->n_heterozygous != chrom_het) ||
        (sc->n_missing_proxy != chrom_missing) ||
        (sc->n_orientation_informative > sc->n_heterozygous) ||
        (sc->n_orientation_tie > sc->n_heterozygous) ||
        (sc->n_orientation_informative + sc->n_orientation_tie !=
         sc->n_heterozygous))
      return 0;

    total_info += sc->n_orientation_informative;
    total_tie += sc->n_orientation_tie;
    expected_start += length;
  }

  if ((expected_start != s->n_marker) ||
      (total_het != s->n_heterozygous) ||
      (total_info != s->n_orientation_informative) ||
      (total_tie != s->n_orientation_tie) ||
      (total_missing != s->n_missing_proxy) ||
      (s->n_orientation_informative + s->n_orientation_tie !=
       s->n_heterozygous))
    return 0;

  {
    size_t chrom_missing = 0u;
    for (size_t c = 0; c < s->n_chrom; c++)
      chrom_missing += s->chrom[c].n_missing_proxy;
    if (chrom_missing != s->n_missing_proxy)
      return 0;
  }

  return 1;
}

static int st_phase_cohort_snapshot_build(const st_phase_data *d,
                                          st_phase_cohort_snapshot *out)
{
  st_phase_cohort_snapshot tmp;
  size_t H, n_genotype, n_hap_cell;

  st_phase_cohort_snapshot_init(&tmp);

  if ((NULL == out) || !st_phase_cohort_snapshot_source_valid(d) ||
      (0 != st_phase_size_mul(d->n_ind, 2u, &H)) ||
      (0 != st_phase_size_mul(d->n_marker, d->n_ind, &n_genotype)) ||
      (0 != st_phase_size_mul(d->n_marker, H, &n_hap_cell)) ||
      (n_hap_cell > SIZE_MAX / sizeof(*tmp.hap)) ||
      (n_genotype > SIZE_MAX / sizeof(*tmp.observed_missing)) ||
      (d->n_marker > SIZE_MAX / sizeof(*tmp.marker)) ||
      (d->n_chrom > SIZE_MAX / sizeof(*tmp.chrom)))
    return -1;

  tmp.hap = st_phase_cohort_malloc(n_hap_cell * sizeof(*tmp.hap));
  tmp.observed_missing =
      st_phase_cohort_malloc(n_genotype * sizeof(*tmp.observed_missing));
  tmp.marker = st_phase_cohort_malloc(d->n_marker * sizeof(*tmp.marker));
  tmp.chrom = st_phase_cohort_malloc(d->n_chrom * sizeof(*tmp.chrom));

  if ((NULL == tmp.hap) || (NULL == tmp.observed_missing) ||
      (NULL == tmp.marker) || (NULL == tmp.chrom))
    goto fail;

  tmp.n_ind = d->n_ind;
  tmp.n_marker = d->n_marker;
  tmp.n_chrom = d->n_chrom;
  tmp.n_hap = H;

  for (size_t m = 0; m < d->n_marker; m++)
  {
    tmp.marker[m].source_marker = d->source_marker[m];
    tmp.marker[m].chrom_label = d->chrom_label[m];
    tmp.marker[m].pos_cm = d->pos_cm[m];

    for (size_t i = 0; i < d->n_ind; i++)
    {
      size_t z = m * d->n_ind + i;
      tmp.observed_missing[z] =
          (ST_PHASE_STATE_MISSING == d->geno1[z]) ? 1u : 0u;
    }
  }

  for (size_t c = 0; c < d->n_chrom; c++)
  {
    st_phase_cohort_initial init;
    size_t start = d->chrom_start[c];
    size_t length = d->chrom_length[c];
    st_phase_cohort_snapshot_chrom *sc = tmp.chrom + c;

    st_phase_cohort_initial_init(&init);

    if (0 != st_phase_cohort_initialize_binary(d->geno1 + start * d->n_ind,
                                               d->geno2 + start * d->n_ind,
                                               length, d->n_ind, &init))
    {
      st_phase_cohort_initial_free(&init);
      goto fail;
    }

    if ((init.n_ind != d->n_ind) || (init.n_marker != length) ||
        (init.n_hap != H) || (NULL == init.hap))
    {
      st_phase_cohort_initial_free(&init);
      goto fail;
    }

    memcpy(tmp.hap + start * H, init.hap, length * H * sizeof(*tmp.hap));

    sc->start = start;
    sc->length = length;
    sc->chrom_label = d->chrom_label[start];
    sc->n_heterozygous = init.n_heterozygous;
    sc->n_orientation_informative = init.n_orientation_informative;
    sc->n_orientation_tie = init.n_orientation_tie;
    sc->n_missing_proxy = init.n_missing_proxy;

    tmp.n_heterozygous += init.n_heterozygous;
    tmp.n_orientation_informative += init.n_orientation_informative;
    tmp.n_orientation_tie += init.n_orientation_tie;
    tmp.n_missing_proxy += init.n_missing_proxy;

    st_phase_cohort_initial_free(&init);
  }

  if (!st_phase_cohort_snapshot_valid(d, &tmp))
    goto fail;

  st_phase_cohort_snapshot_free(out);
  *out = tmp;
  return 0;

fail:
  st_phase_cohort_snapshot_free(&tmp);
  return -1;
}

#define st_phase_cohort_rephase_malloc malloc

typedef struct
{
  const st_phase_cohort_snapshot *snapshot;
  size_t marker_start;
  size_t n_marker;
  size_t n_hap;
  size_t hap_index;
} st_phase_cohort_reference_descriptor;

static int st_phase_cohort_reference_descriptor_cmp(const void *va,
                                                    const void *vb)
{
  const st_phase_cohort_reference_descriptor *a =
      (const st_phase_cohort_reference_descriptor *)va;
  const st_phase_cohort_reference_descriptor *b =
      (const st_phase_cohort_reference_descriptor *)vb;
  size_t n = (a->n_marker < b->n_marker) ? a->n_marker : b->n_marker;

  for (size_t m = 0u; m < n; m++)
  {
    uint8_t x =
        a->snapshot->hap[(a->marker_start + m) * a->n_hap + a->hap_index];
    uint8_t y =
        b->snapshot->hap[(b->marker_start + m) * b->n_hap + b->hap_index];
    if (x < y)
      return -1;
    if (x > y)
      return 1;
  }
  if (a->n_marker < b->n_marker)
    return -1;
  if (a->n_marker > b->n_marker)
    return 1;
  return 0;
}

static int st_phase_cohort_build_loo_reference(
    const st_phase_cohort_snapshot *frozen, size_t target,
    size_t chromosome, uint8_t **reference_out, size_t *n_ref_out)
{
  st_phase_cohort_reference_descriptor *desc = NULL;
  uint8_t *reference = NULL;
  size_t start, M, H, n_donor = 0u, R, ncell;
  int status = -1;

  if ((NULL == frozen) || (NULL == reference_out) || (NULL == n_ref_out) ||
      (target >= frozen->n_ind) || (chromosome >= frozen->n_chrom) ||
      (NULL == frozen->hap) || (NULL == frozen->observed_missing) ||
      (NULL == frozen->chrom) || (frozen->n_ind < 2u))
    return -1;

  *reference_out = NULL;
  *n_ref_out = 0u;
  start = frozen->chrom[chromosome].start;
  M = frozen->chrom[chromosome].length;
  H = frozen->n_hap;
  if ((0u == M) || (H != 2u * frozen->n_ind))
    return -1;

  for (size_t i = 0u; i < frozen->n_ind; i++)
    if (i != target)
    {
      int complete = 1;
      for (size_t m = 0u; m < M; m++)
        if (frozen->observed_missing[(start + m) * frozen->n_ind + i])
        {
          complete = 0;
          break;
        }
      if (complete)
        n_donor++;
    }

  if (0u == n_donor)
    return 0;
  if (n_donor > SIZE_MAX / 2u)
    return -1;
  R = 2u * n_donor;
  if ((R > SIZE_MAX / sizeof(*desc)) ||
      (0 != st_phase_size_mul(M, R, &ncell)) ||
      (ncell > SIZE_MAX / sizeof(*reference)))
    return -1;

  desc = st_phase_cohort_rephase_malloc(R * sizeof(*desc));
  reference = st_phase_cohort_rephase_malloc(ncell * sizeof(*reference));
  if ((NULL == desc) || (NULL == reference))
    goto cleanup;

  {
    size_t r = 0u;
    for (size_t i = 0u; i < frozen->n_ind; i++)
      if (i != target)
      {
        int complete = 1;
        for (size_t m = 0u; m < M; m++)
          if (frozen->observed_missing[(start + m) * frozen->n_ind + i])
          {
            complete = 0;
            break;
          }
        if (!complete)
          continue;
        for (size_t side = 0u; side < 2u; side++)
        {
          desc[r].snapshot = frozen;
          desc[r].marker_start = start;
          desc[r].n_marker = M;
          desc[r].n_hap = H;
          desc[r].hap_index = 2u * i + side;
          r++;
        }
      }
    if (r != R)
      goto cleanup;
  }

  qsort(desc, R, sizeof(*desc), st_phase_cohort_reference_descriptor_cmp);
  for (size_t m = 0u; m < M; m++)
    for (size_t r = 0u; r < R; r++)
    {
      uint8_t x = frozen->hap[(start + m) * H + desc[r].hap_index];
      if (x > 1u)
        goto cleanup;
      reference[m * R + r] = x;
    }

  *reference_out = reference;
  *n_ref_out = R;
  reference = NULL;
  status = 0;
cleanup:
  free(desc);
  free(reference);
  return status;
}

typedef struct
{
  int possible;
  st_phase_relative_orientation_result relative;
} st_phase_cohort_relative_chromosome_result;

typedef struct
{
  st_phase_cohort_snapshot snapshot;
  size_t n_target_chromosome;
  size_t n_possible;
  size_t n_impossible;
  size_t n_unavailable;
  size_t n_relative_link;
  size_t n_posterior_same;
  size_t n_posterior_switch;
  size_t n_posterior_tie;
  size_t n_changed_relative_link;
  size_t n_changed_heterozygous_marker;
} st_phase_cohort_rephase_once_result;

static int
st_phase_cohort_snapshot_clone(const st_phase_cohort_snapshot *src,
                               st_phase_cohort_snapshot *out)
{
  st_phase_cohort_snapshot tmp;
  size_t hap_n, geno_n;

  st_phase_cohort_snapshot_init(&tmp);
  if ((NULL == src) || (NULL == out) || (0u == src->n_ind) ||
      (0u == src->n_marker) || (0u == src->n_chrom) ||
      (src->n_hap != 2u * src->n_ind) ||
      (0 != st_phase_size_mul(src->n_marker, src->n_hap, &hap_n)) ||
      (0 != st_phase_size_mul(src->n_marker, src->n_ind, &geno_n)) ||
      (hap_n > SIZE_MAX / sizeof(*tmp.hap)) ||
      (geno_n > SIZE_MAX / sizeof(*tmp.observed_missing)) ||
      (src->n_marker > SIZE_MAX / sizeof(*tmp.marker)) ||
      (src->n_chrom > SIZE_MAX / sizeof(*tmp.chrom)) ||
      (NULL == src->hap) || (NULL == src->observed_missing) ||
      (NULL == src->marker) || (NULL == src->chrom))
    return -1;

  tmp.hap = st_phase_cohort_rephase_malloc(hap_n * sizeof(*tmp.hap));
  tmp.observed_missing = st_phase_cohort_rephase_malloc(
      geno_n * sizeof(*tmp.observed_missing));
  tmp.marker =
      st_phase_cohort_rephase_malloc(src->n_marker * sizeof(*tmp.marker));
  tmp.chrom =
      st_phase_cohort_rephase_malloc(src->n_chrom * sizeof(*tmp.chrom));
  if ((NULL == tmp.hap) || (NULL == tmp.observed_missing) ||
      (NULL == tmp.marker) || (NULL == tmp.chrom))
    goto fail;

  tmp.n_ind = src->n_ind;
  tmp.n_marker = src->n_marker;
  tmp.n_chrom = src->n_chrom;
  tmp.n_hap = src->n_hap;
  tmp.n_heterozygous = src->n_heterozygous;
  tmp.n_orientation_informative = src->n_orientation_informative;
  tmp.n_orientation_tie = src->n_orientation_tie;
  tmp.n_missing_proxy = src->n_missing_proxy;

  memcpy(tmp.hap, src->hap, hap_n * sizeof(*tmp.hap));
  memcpy(tmp.observed_missing, src->observed_missing,
         geno_n * sizeof(*tmp.observed_missing));
  memcpy(tmp.marker, src->marker, src->n_marker * sizeof(*tmp.marker));
  memcpy(tmp.chrom, src->chrom, src->n_chrom * sizeof(*tmp.chrom));

  st_phase_cohort_snapshot_free(out);
  *out = tmp;
  return 0;

fail:
  st_phase_cohort_snapshot_free(&tmp);
  return -1;
}

static int st_phase_cohort_relative_relation_from_snapshot(
    const st_phase_data *d, const st_phase_cohort_snapshot *s,
    size_t individual, size_t left, size_t right, int *same_out)
{
  size_t H, hleft, hright;
  uint8_t aL, bL, aR, bR, xL, xR;

  if ((NULL == d) || (NULL == s) || (NULL == same_out) ||
      (individual >= d->n_ind) || (left >= d->n_marker) ||
      (right >= d->n_marker) || (left >= right) ||
      (s->n_hap != 2u * d->n_ind))
    return -1;

  H = s->n_hap;
  aL = d->geno1[left * d->n_ind + individual];
  bL = d->geno2[left * d->n_ind + individual];
  aR = d->geno1[right * d->n_ind + individual];
  bR = d->geno2[right * d->n_ind + individual];
  if ((ST_PHASE_STATE_MISSING == aL) || (ST_PHASE_STATE_MISSING == aR) ||
      !(aL < bL) || !(aR < bR))
    return -1;

  hleft = left * H + 2u * individual;
  hright = right * H + 2u * individual;
  xL = s->hap[hleft];
  xR = s->hap[hright];
  if (!((xL == aL) || (xL == bL)) || !((xR == aR) || (xR == bR)))
    return -1;

  *same_out = ((xL == aL) == (xR == aR)) ? 1 : 0;
  return 0;
}

#define st_phase_relative_factor_malloc malloc
#define st_phase_relative_factor_calloc calloc

typedef struct
{
  size_t iteration;
  size_t n_target_chromosome;
  size_t n_possible;
  size_t n_impossible;
  size_t n_unavailable;
  size_t n_relative_link;
  size_t n_posterior_same;
  size_t n_posterior_switch;
  size_t n_posterior_tie;
  size_t n_changed_relative_link;
  size_t n_changed_heterozygous_marker;
  double change_fraction;
} st_phase_cohort_burnin_iteration;

typedef struct
{
  st_phase_cohort_snapshot snapshot;
  size_t max_iterations;
  size_t min_iterations;
  double convergence_fraction;
  size_t n_iteration_completed;
  int converged;
  st_phase_cohort_burnin_iteration *iteration;
} st_phase_cohort_burnin_result;

#define st_phase_cohort_burnin_malloc malloc

typedef struct
{
  int possible;
  size_t n_marker;
  uint8_t *orientation;
  uint8_t *missing1;
  uint8_t *missing2;
} st_phase_cohort_viterbi_chromosome_result;
static void st_phase_cohort_viterbi_chromosome_result_init(
    st_phase_cohort_viterbi_chromosome_result *r)
{
  if (r)
    memset(r, 0, sizeof(*r));
}
static void st_phase_cohort_viterbi_chromosome_result_free(
    st_phase_cohort_viterbi_chromosome_result *r)
{
  if (!r)
    return;
  free(r->orientation);
  free(r->missing1);
  free(r->missing2);
  st_phase_cohort_viterbi_chromosome_result_init(r);
}

#define ST_PHASE_COHORT_SCAFFOLD_FREE ((uint8_t)0u)
#define ST_PHASE_COHORT_SCAFFOLD_FIXED_SAME ((uint8_t)1u)
#define ST_PHASE_COHORT_SCAFFOLD_FIXED_SWITCH ((uint8_t)2u)

typedef struct
{
  size_t n_ind;
  size_t n_marker;
  size_t n_chrom;
  size_t n_fixed_same;
  size_t n_fixed_switch;

  uint8_t *relation;
} st_phase_cohort_scaffold;

#define st_phase_cohort_scaffold_malloc malloc

static int st_phase_cohort_observed_heterozygous(const st_phase_data *d,
                                                 size_t marker,
                                                 size_t individual)
{
  uint8_t a, b;
  if ((NULL == d) || (marker >= d->n_marker) || (individual >= d->n_ind))
    return 0;
  a = d->geno1[marker * d->n_ind + individual];
  b = d->geno2[marker * d->n_ind + individual];
  return (ST_PHASE_STATE_MISSING != a) && (ST_PHASE_STATE_MISSING != b) &&
         (a < b);
}

typedef struct
{
  st_phase_cohort_snapshot snapshot;
  size_t n_target_chromosome;
  size_t n_possible;
  size_t n_impossible;
  size_t n_unavailable;
  size_t n_relative_link;
  size_t n_posterior_same;
  size_t n_posterior_switch;
  size_t n_posterior_tie;
  size_t n_unconstrained_changed_relative_link;
  size_t n_unconstrained_changed_heterozygous_marker;
  size_t n_fixed_same;
  size_t n_fixed_switch;
  size_t n_fixed_overrode_unconstrained;
  size_t n_final_changed_relative_link;
  size_t n_final_changed_heterozygous_marker;
} st_phase_cohort_rephase_scaffold_result;

typedef struct
{
  st_phase_cohort_scaffold scaffold;
  double confidence_threshold;
  size_t n_valid_link;
  size_t n_fixed_before;
  size_t n_free_before;
  size_t n_free_evaluable;
  size_t n_free_not_evaluable;
  size_t n_free_tie;
  size_t n_free_below_threshold;
  size_t n_new_fixed_same;
  size_t n_new_fixed_switch;
  size_t n_fixed_after;
  size_t n_free_after;
} st_phase_cohort_fix_confident_result;

#define st_phase_cohort_fix_malloc malloc

typedef struct
{
  size_t iteration;
  size_t n_target_chromosome;
  size_t n_possible;
  size_t n_impossible;
  size_t n_unavailable;
  size_t n_relative_link;
  size_t n_posterior_same;
  size_t n_posterior_switch;
  size_t n_posterior_tie;
  size_t n_fixed_overrode_unconstrained;
  size_t n_changed_relative_link;
  size_t n_changed_heterozygous_marker;
  double change_fraction;
  size_t n_valid_link;
  size_t n_fixed_before;
  size_t n_free_before;
  size_t n_free_evaluable;
  size_t n_free_not_evaluable;
  size_t n_free_tie;
  size_t n_free_below_threshold;
  size_t n_new_fixed_same;
  size_t n_new_fixed_switch;
  size_t n_fixed_after;
  size_t n_free_after;
} st_phase_cohort_progressive_iteration;

typedef struct
{
  st_phase_cohort_snapshot snapshot;
  st_phase_cohort_scaffold scaffold;
  size_t max_iterations;
  size_t min_iterations;
  double convergence_fraction;
  double confidence_threshold;
  size_t max_states;
  double error_rate;
  double break_rate1_per_morgan;
  double break_rate2_per_morgan;
  size_t n_iteration_completed;
  int stopped_all_fixed;
  int stopped_converged;
  st_phase_cohort_progressive_iteration *iteration;
} st_phase_cohort_progressive_result;

#define st_phase_cohort_progressive_malloc malloc

typedef struct
{
  st_phase_cohort_snapshot snapshot;
  st_phase_cohort_scaffold scaffold;

  size_t initial_n_heterozygous;
  size_t initial_n_orientation_informative;
  size_t initial_n_orientation_tie;
  size_t initial_n_missing_proxy;

  size_t max_states;
  double error_rate;
  double break_rate1_per_morgan;
  double break_rate2_per_morgan;

  size_t burnin_max_iterations;
  size_t burnin_min_iterations;
  double burnin_convergence_fraction;
  size_t n_burnin_iteration_completed;
  int burnin_converged;
  st_phase_cohort_burnin_iteration *burnin_iteration;

  size_t progressive_max_iterations;
  size_t progressive_min_iterations;
  double progressive_convergence_fraction;
  double confidence_threshold;
  size_t n_progressive_iteration_completed;
  int progressive_stopped_all_fixed;
  int progressive_stopped_converged;
  st_phase_cohort_progressive_iteration *progressive_iteration;
} st_phase_cohort_end_to_end_result;

typedef struct
{
  st_phase_dataset_stage stage;
  st_phase_stage_audit audit;
  int audit_status;

  size_t n_ind;
  size_t n_marker;
  size_t n_chrom;
  size_t n_target_chromosome;
  size_t n_phase_ready_chromosome;
  size_t n_partial_fixed_chromosome;
  size_t n_unresolved_chromosome;
  size_t n_no_heterozygote_chromosome;

  size_t n_valid_link;
  size_t n_fixed_link;
  size_t n_free_link;
  size_t n_fixed_same;
  size_t n_fixed_switch;

  size_t n_observed_heterozygous;
  size_t n_confident_heterozygous;
  size_t n_unresolved_heterozygous;
  size_t n_confident_same_source_order;
  size_t n_confident_swapped_source_order;

  unsigned char *chromosome_phase_ready;
} st_phase_cohort_stage_result;

#define st_phase_cohort_stage_malloc malloc

typedef enum
{
  ST_PHASE_COHORT_PIPELINE_OK = 0,
  ST_PHASE_COHORT_PIPELINE_INVALID = -1,
  ST_PHASE_COHORT_PIPELINE_DATA_FAILED = -2,
  ST_PHASE_COHORT_PIPELINE_PHASING_FAILED = -3,
  ST_PHASE_COHORT_PIPELINE_STAGING_FAILED = -4,
  ST_PHASE_COHORT_PIPELINE_COMMIT_FAILED = -5
} st_phase_cohort_pipeline_status;

typedef struct
{
  int pipeline_status;
  int phasing_succeeded;
  int staging_succeeded;
  int audit_status;
  int commit_requested;
  int commit_attempted;
  int commit_status;
  int committed;

  size_t n_ind;
  size_t n_marker;
  size_t n_chrom;
  size_t n_genotype;

  size_t initial_n_heterozygous;
  size_t initial_n_orientation_informative;
  size_t initial_n_orientation_tie;
  size_t initial_n_missing_proxy;

  size_t burnin_iterations_completed;
  int burnin_converged;
  size_t burnin_last_relative_link;
  size_t burnin_last_changed_relative_link;
  size_t burnin_last_changed_heterozygous_marker;
  double burnin_last_change_fraction;

  size_t progressive_iterations_completed;
  int progressive_stopped_all_fixed;
  int progressive_stopped_converged;
  size_t progressive_last_relative_link;
  size_t progressive_last_fixed_before;
  size_t progressive_last_new_fixed_same;
  size_t progressive_last_new_fixed_switch;
  size_t progressive_last_fixed_after;
  size_t progressive_last_free_after;
  size_t progressive_last_changed_relative_link;
  double progressive_last_change_fraction;

  size_t n_target_chromosome;
  size_t n_phase_ready_chromosome;
  size_t n_partial_fixed_chromosome;
  size_t n_unresolved_chromosome;
  size_t n_no_heterozygote_chromosome;

  size_t n_valid_link;
  size_t n_fixed_link;
  size_t n_free_link;
  size_t n_fixed_same;
  size_t n_fixed_switch;

  size_t n_observed_heterozygous;
  size_t n_confident_heterozygous;
  size_t n_unresolved_heterozygous;

  size_t n_covered;
  size_t n_uncovered;
  size_t n_written;
  size_t n_covered_unwritten;
  size_t n_changed_genotypes;
  size_t n_missing;
  size_t n_homozygous;
  size_t n_heterozygous;
  size_t n_covered_heterozygous_unresolved;
  size_t n_written_same_source_order;
  size_t n_written_swapped_source_order;

  int audit_complete;
  int audit_source_bytes_current;
  int audit_workspace_current;
} st_phase_cohort_pipeline_diagnostics;

typedef struct
{
  size_t n_marker;
  size_t n_missing;
  size_t n_imputed;
  size_t n_unresolved;
  size_t n_imputed_hom0;
  size_t n_imputed_hom1;
  size_t n_imputed_heterozygous;
  size_t n_imputed_heterozygous_phased;
  size_t n_imputed_heterozygous_unphased;
  unsigned char *is_missing;
  unsigned char *imputed;
  unsigned char *phased;
  uint8_t *allele1;
  uint8_t *allele2;
  double *genotype_probability;
  double *phase_probability;
} st_phase_step0_imputation_result;

static void
st_phase_step0_imputation_result_init(st_phase_step0_imputation_result *x)
{
  if (NULL != x)
    memset(x, 0, sizeof(*x));
}

static void
st_phase_step0_imputation_result_free(st_phase_step0_imputation_result *x)
{
  if (NULL == x)
    return;
  free(x->is_missing);
  free(x->imputed);
  free(x->phased);
  free(x->allele1);
  free(x->allele2);
  free(x->genotype_probability);
  free(x->phase_probability);
  st_phase_step0_imputation_result_init(x);
}

static int st_phase_step0_imputation_threshold_valid(double x)
{
  return isfinite(x) && (x >= 0.0) && (x <= 1.0);
}

static int st_phase_step0_probability_close(double x, double y,
                                            double scale)
{
  double tol =
      1024.0 * DBL_EPSILON * (1.0 + fabs(x) + fabs(y) + fabs(scale));
  return fabs(x - y) <= tol;
}

typedef struct
{
  st_phase_diploid_input input;
  st_phase_unphased_candidate_set candidate;
  st_phase_pair_state_set pair;
  st_phase_pair_hmm_result hmm;
  double *phase12;
  double *phase21;
  double *emission;
  double *q1;
  double *q2;
  double *initial;
  size_t M;
  size_t P;
} st_phase_step0_hmm_work;

static void st_phase_step0_hmm_work_init(st_phase_step0_hmm_work *w)
{
  if (NULL == w)
    return;
  memset(w, 0, sizeof(*w));
  st_phase_diploid_input_init(&w->input);
  st_phase_unphased_candidate_set_init(&w->candidate);
  st_phase_pair_state_set_init(&w->pair);
  st_phase_pair_hmm_result_init(&w->hmm);
}

static void st_phase_step0_hmm_work_free(st_phase_step0_hmm_work *w)
{
  if (NULL == w)
    return;
  free(w->phase12);
  free(w->phase21);
  free(w->emission);
  free(w->q1);
  free(w->q2);
  free(w->initial);
  st_phase_pair_hmm_result_free(&w->hmm);
  st_phase_pair_state_set_free(&w->pair);
  st_phase_unphased_candidate_set_free(&w->candidate);
  st_phase_diploid_input_free(&w->input);
  st_phase_step0_hmm_work_init(w);
}

static int st_phase_step0_hmm_build(
    const st_phase_data *d, size_t individual, size_t chromosome,
    const uint8_t *reference_hap, size_t n_ref_hap, size_t max_states,
    const double *initial_pair_weight, size_t initial_pair_count,
    double error_rate, double break_rate1_per_morgan,
    double break_rate2_per_morgan, st_phase_step0_hmm_work *out)
{
  st_phase_step0_hmm_work tmp;
  size_t pair_count = 0u;
  int status = -1;

  st_phase_step0_hmm_work_init(&tmp);
  if ((NULL == out) || (NULL == reference_hap) || (0u == n_ref_hap) ||
      (0u == max_states) || !st_phase_data_adapter_valid(d) ||
      !st_phase_hmm_error_valid(error_rate) ||
      !isfinite(break_rate1_per_morgan) ||
      !isfinite(break_rate2_per_morgan) ||
      (break_rate1_per_morgan < 0.0) || (break_rate2_per_morgan < 0.0) ||
      (0 != st_phase_size_mul(n_ref_hap, n_ref_hap, &pair_count)))
    goto cleanup;

  if (NULL != initial_pair_weight)
  {
    double sum = 0.0;
    if (initial_pair_count != pair_count)
      goto cleanup;
    for (size_t z = 0u; z < pair_count; z++)
    {
      double x = initial_pair_weight[z];
      if (!isfinite(x) || (x < 0.0))
        goto cleanup;
      sum += x;
      if (!isfinite(sum))
        goto cleanup;
    }
    if (!(sum > 0.0))
      goto cleanup;
  }
  else if (0u != initial_pair_count)
    goto cleanup;

  if (0 != st_phase_prepare_diploid_input(d, individual, chromosome,
                                          &tmp.input))
    goto cleanup;
  tmp.M = tmp.input.n_marker;
  if (0u == tmp.M)
    goto cleanup;

  if (0 != st_phase_unphased_candidate_build(
               tmp.input.observed1, tmp.input.observed2,
               d->n_alleles + tmp.input.marker_start, reference_hap, tmp.M,
               n_ref_hap, max_states, &tmp.candidate))
    goto cleanup;
  if (0 != st_phase_pair_state_build_unphased(&tmp.candidate, &tmp.pair))
    goto cleanup;
  tmp.P = tmp.pair.max_pair_states;
  if ((0u == tmp.P) || (tmp.pair.n_slot != tmp.M * tmp.P))
    goto cleanup;

  tmp.phase12 = calloc(tmp.pair.n_slot, sizeof(*tmp.phase12));
  tmp.phase21 = calloc(tmp.pair.n_slot, sizeof(*tmp.phase21));
  tmp.emission = calloc(tmp.pair.n_slot, sizeof(*tmp.emission));
  if ((NULL == tmp.phase12) || (NULL == tmp.phase21) ||
      (NULL == tmp.emission))
    goto cleanup;

  for (size_t marker = 0u; marker < tmp.M; marker++)
  {
    size_t g = tmp.input.marker_start + marker;
    if (0 != st_phase_hmm_pair_emission_row(
                 &tmp.pair, marker, reference_hap + marker * n_ref_hap,
                 n_ref_hap, tmp.input.observed1[marker],
                 tmp.input.observed2[marker], d->n_alleles[g], error_rate,
                 tmp.phase12 + marker * tmp.P,
                 tmp.phase21 + marker * tmp.P,
                 tmp.emission + marker * tmp.P))
      goto cleanup;
  }

  if (tmp.M > 1u)
  {
    tmp.q1 = malloc((tmp.M - 1u) * sizeof(*tmp.q1));
    tmp.q2 = malloc((tmp.M - 1u) * sizeof(*tmp.q2));
    if ((NULL == tmp.q1) || (NULL == tmp.q2))
      goto cleanup;
    for (size_t marker = 1u; marker < tmp.M; marker++)
    {
      size_t g = tmp.input.marker_start + marker;
      double distance = d->pos_cm[g] - d->pos_cm[g - 1u];
      if ((0 != st_phase_hmm_break_probability(distance,
                                               break_rate1_per_morgan,
                                               &tmp.q1[marker - 1u])) ||
          (0 != st_phase_hmm_break_probability(distance,
                                               break_rate2_per_morgan,
                                               &tmp.q2[marker - 1u])))
        goto cleanup;
    }
  }

  tmp.initial = calloc(tmp.P, sizeof(*tmp.initial));
  if (NULL == tmp.initial)
    goto cleanup;
  {
    size_t n_first = 0u;
    const st_phase_pair_state *first =
        st_phase_pair_state_get(&tmp.pair, 0u, &n_first, NULL, NULL);
    if ((NULL == first) || (0u == n_first) || (n_first > tmp.P))
      goto cleanup;
    for (size_t z = 0u; z < n_first; z++)
    {
      if ((first[z].hap1 >= n_ref_hap) || (first[z].hap2 >= n_ref_hap))
        goto cleanup;
      tmp.initial[z] =
          (NULL == initial_pair_weight)
              ? 1.0
              : initial_pair_weight[first[z].hap1 * n_ref_hap +
                                    first[z].hap2];
    }
  }

  if (0 !=
      st_phase_pair_hmm_forward_backward_viterbi_factorized(
          &tmp.pair, tmp.initial, tmp.emission, tmp.q1, tmp.q2, &tmp.hmm))
    goto cleanup;

  st_phase_step0_hmm_work_free(out);
  *out = tmp;
  st_phase_step0_hmm_work_init(&tmp);
  status = 0;
cleanup:
  st_phase_step0_hmm_work_free(&tmp);
  return status;
}

static int st_phase_step0_complete_reference_calls_viterbi(
    const st_phase_step0_hmm_work *w, st_phase_hard_call_result *hard)
{
  st_phase_hard_call current = ST_PHASE_HARD_12;

  if ((NULL == w) || (NULL == hard) || (hard->n_marker != w->M) ||
      (NULL == hard->call) || (NULL == hard->confidence) ||
      !st_phase_pair_hmm_set_valid(&w->pair) ||
      (w->hmm.n_marker != w->M) || (NULL == w->hmm.viterbi_hap1) ||
      (NULL == w->hmm.viterbi_hap2) || (NULL == w->phase12) ||
      (NULL == w->phase21))
    return -1;

  for (size_t m = 0u; m < w->M; m++)
  {
    uint8_t a = w->input.observed1[m], b = w->input.observed2[m];
    if ((ST_PHASE_STATE_MISSING == a) || (a == b))
      continue;

    if ((ST_PHASE_HARD_12 == hard->call[m]) ||
        (ST_PHASE_HARD_21 == hard->call[m]))
    {
      current = hard->call[m];
      continue;
    }

    if (ST_PHASE_HARD_UNRESOLVED != hard->call[m])
      return -1;
    if (!w->hmm.possible || SIZE_MAX == w->hmm.viterbi_hap1[m] ||
        SIZE_MAX == w->hmm.viterbi_hap2[m])
      continue;

    {
      size_t n = 0u, slot = SIZE_MAX;
      const st_phase_pair_state *row =
          st_phase_pair_state_get(&w->pair, m, &n, NULL, NULL);
      if ((NULL == row) || (0u == n) || (n > w->P))
        return -1;
      for (size_t z = 0u; z < n; z++)
        if ((row[z].hap1 == w->hmm.viterbi_hap1[m]) &&
            (row[z].hap2 == w->hmm.viterbi_hap2[m]))
        {
          slot = z;
          break;
        }
      if (SIZE_MAX == slot)
        return -1;
      {
        double p12 = w->phase12[m * w->P + slot];
        double p21 = w->phase21[m * w->P + slot];
        if (!isfinite(p12) || !isfinite(p21) || (p12 < 0.0) || (p21 < 0.0))
          return -1;
        if (p12 > p21)
          current = ST_PHASE_HARD_12;
        else if (p21 > p12)
          current = ST_PHASE_HARD_21;
        hard->call[m] = current;
        if (hard->confidence[m] < 0.5)
          hard->confidence[m] = 0.5;
      }
    }
  }
  return 0;
}

static int st_phase_step0_decode_missing(
    const st_phase_data *d, const uint8_t *reference_hap, size_t n_ref_hap,
    const st_phase_step0_hmm_work *w, double confidence_threshold,
    double imputation_threshold, st_phase_step0_imputation_result *out)
{
  st_phase_step0_imputation_result tmp;
  size_t M, P;

  st_phase_step0_imputation_result_init(&tmp);
  if ((NULL == out) || (NULL == reference_hap) || (0u == n_ref_hap) ||
      (NULL == w) || !st_phase_data_adapter_valid(d) ||
      !st_phase_hard_threshold_valid(confidence_threshold) ||
      !st_phase_step0_imputation_threshold_valid(imputation_threshold) ||
      !st_phase_pair_hmm_set_valid(&w->pair) ||
      (w->hmm.n_marker != w->pair.n_marker) ||
      (w->hmm.max_pair_states != w->pair.max_pair_states) ||
      (w->hmm.n_slot != w->pair.n_slot) || (NULL == w->hmm.posterior) ||
      (w->input.n_marker != w->pair.n_marker))
    return -1;

  M = w->pair.n_marker;
  P = w->pair.max_pair_states;
  tmp.n_marker = M;
  tmp.is_missing = calloc(M, sizeof(*tmp.is_missing));
  tmp.imputed = calloc(M, sizeof(*tmp.imputed));
  tmp.phased = calloc(M, sizeof(*tmp.phased));
  tmp.allele1 = malloc(M * sizeof(*tmp.allele1));
  tmp.allele2 = malloc(M * sizeof(*tmp.allele2));
  tmp.genotype_probability = calloc(M, sizeof(*tmp.genotype_probability));
  tmp.phase_probability = calloc(M, sizeof(*tmp.phase_probability));
  if ((NULL == tmp.is_missing) || (NULL == tmp.imputed) ||
      (NULL == tmp.phased) || (NULL == tmp.allele1) ||
      (NULL == tmp.allele2) || (NULL == tmp.genotype_probability) ||
      (NULL == tmp.phase_probability))
    goto fail;
  for (size_t marker = 0u; marker < M; marker++)
  {
    tmp.allele1[marker] = ST_PHASE_STATE_MISSING;
    tmp.allele2[marker] = ST_PHASE_STATE_MISSING;
  }

  for (size_t marker = 0u; marker < M; marker++)
  {
    uint8_t o1 = w->input.observed1[marker],
            o2 = w->input.observed2[marker];
    size_t g = w->input.marker_start + marker, n = 0u;
    const st_phase_pair_state *row;
    double p00 = 0.0, p01 = 0.0, p10 = 0.0, p11 = 0.0, total = 0.0;
    double gp[3], best = -1.0, second = -1.0;
    size_t best_g = SIZE_MAX;

    if ((ST_PHASE_STATE_MISSING == o1) != (ST_PHASE_STATE_MISSING == o2))
      goto fail;
    if (ST_PHASE_STATE_MISSING != o1)
      continue;
    tmp.is_missing[marker] = 1u;
    tmp.n_missing++;

    if ((0u == d->n_alleles[g]) || (d->n_alleles[g] > 2u) ||
        !w->hmm.possible)
    {
      tmp.n_unresolved++;
      continue;
    }

    row = st_phase_pair_state_get(&w->pair, marker, &n, NULL, NULL);
    if ((NULL == row) || (0u == n) || (n > P))
      goto fail;
    for (size_t z = 0u; z < n; z++)
    {
      double post = w->hmm.posterior[marker * P + z];
      uint8_t a, b;
      if (!isfinite(post) || (post < 0.0) || (post > 1.0) ||
          (row[z].hap1 >= n_ref_hap) || (row[z].hap2 >= n_ref_hap))
        goto fail;
      a = reference_hap[marker * n_ref_hap + row[z].hap1];
      b = reference_hap[marker * n_ref_hap + row[z].hap2];
      if ((a >= d->n_alleles[g]) || (b >= d->n_alleles[g]))
        goto fail;
      if ((0u == a) && (0u == b))
        p00 += post;
      else if ((0u == a) && (1u == b))
        p01 += post;
      else if ((1u == a) && (0u == b))
        p10 += post;
      else if ((1u == a) && (1u == b))
        p11 += post;
      else
        goto fail;
      total += post;
    }
    if (!isfinite(total) ||
        !st_phase_step0_probability_close(total, 1.0, (double)n))
      goto fail;
    if (total <= 0.0)
    {
      tmp.n_unresolved++;
      continue;
    }
    p00 /= total;
    p01 /= total;
    p10 /= total;
    p11 /= total;
    gp[0] = p00;
    gp[1] = p01 + p10;
    gp[2] = p11;
    for (size_t q = 0u; q < 3u; q++)
    {
      double v = gp[q];
      if ((v > best) || ((v == best) && (q < best_g)))
      {
        second = best;
        best = v;
        best_g = q;
      }
      else if (v > second)
        second = v;
    }
    if ((SIZE_MAX == best_g) || !isfinite(best) ||
        (best < imputation_threshold) ||
        st_phase_step0_probability_close(best, second, 1.0))
    {
      tmp.n_unresolved++;
      continue;
    }

    tmp.imputed[marker] = 1u;
    tmp.genotype_probability[marker] = best;
    tmp.n_imputed++;
    if (0u == best_g)
    {
      tmp.allele1[marker] = 0u;
      tmp.allele2[marker] = 0u;
      tmp.phase_probability[marker] = 1.0;
      tmp.n_imputed_hom0++;
    }
    else if (2u == best_g)
    {
      if (d->n_alleles[g] < 2u)
        goto fail;
      tmp.allele1[marker] = 1u;
      tmp.allele2[marker] = 1u;
      tmp.phase_probability[marker] = 1.0;
      tmp.n_imputed_hom1++;
    }
    else
    {
      double het = gp[1], q01, q10, phase_best;
      if ((d->n_alleles[g] < 2u) || !(het > 0.0))
        goto fail;
      q01 = p01 / het;
      q10 = p10 / het;
      phase_best = (q01 > q10) ? q01 : q10;
      tmp.n_imputed_heterozygous++;
      tmp.phase_probability[marker] = phase_best;
      if ((phase_best >= confidence_threshold) &&
          !st_phase_step0_probability_close(q01, q10, 1.0))
      {
        tmp.phased[marker] = 1u;
        if (q10 > q01)
        {
          tmp.allele1[marker] = 1u;
          tmp.allele2[marker] = 0u;
        }
        else
        {
          tmp.allele1[marker] = 0u;
          tmp.allele2[marker] = 1u;
        }
        tmp.n_imputed_heterozygous_phased++;
      }
      else
      {
        tmp.allele1[marker] = 0u;
        tmp.allele2[marker] = 1u;
        tmp.n_imputed_heterozygous_unphased++;
      }
    }
  }

  if ((tmp.n_imputed + tmp.n_unresolved != tmp.n_missing) ||
      (tmp.n_imputed_hom0 + tmp.n_imputed_hom1 +
           tmp.n_imputed_heterozygous !=
       tmp.n_imputed) ||
      (tmp.n_imputed_heterozygous_phased +
           tmp.n_imputed_heterozygous_unphased !=
       tmp.n_imputed_heterozygous))
    goto fail;

  st_phase_step0_imputation_result_free(out);
  *out = tmp;
  st_phase_step0_imputation_result_init(&tmp);
  return 0;
fail:
  st_phase_step0_imputation_result_free(&tmp);
  return -1;
}

static int st_phase_step0_stage_structurally_valid(
    const st_phase_data *d, const st_phase_dataset_stage *stage)
{
  size_t n = 0u, claimed = 0u, covered = 0u;
  if (!st_phase_data_adapter_valid(d) || (NULL == stage) ||
      (stage->n_ind != d->n_ind) || (stage->n_marker != d->n_marker) ||
      (0 != st_phase_size_mul(d->n_ind, d->n_marker, &n)) ||
      (stage->n_genotype != n) || (NULL == stage->base1) ||
      (NULL == stage->base2) || (NULL == stage->staged1) ||
      (NULL == stage->staged2) || (NULL == stage->claimed) ||
      (NULL == stage->covered))
    return 0;
  for (size_t m = 0u; m < d->n_marker; m++)
    for (size_t i = 0u; i < d->n_ind; i++)
    {
      size_t raw = (size_t)d->source_marker[m] * d->n_ind + i;
      GS_ALLELE_TYPE a = stage->base1[raw], b = stage->base2[raw];
      GS_ALLELE_TYPE x = stage->staged1[raw], y = stage->staged2[raw];
      unsigned char c = stage->claimed[raw], v = stage->covered[raw];
      if ((c > 1u) || (v > 1u) || (c && !v))
        return 0;
      if (c)
        claimed++;
      if (v)
        covered++;
      if ((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b))
      {
        int fx = 0, fy = 0;
        if (!((GS_ALLELE_MISSING == a) && (GS_ALLELE_MISSING == b)))
          return 0;
        if (!c)
        {
          if ((x != a) || (y != b))
            return 0;
          continue;
        }
        if ((GS_ALLELE_MISSING == x) || (GS_ALLELE_MISSING == y) ||
            (0u == x) || (0u == y) || (0u == d->n_alleles[m]) ||
            (d->n_alleles[m] > 2u))
          return 0;
        for (size_t q = 0u; q < d->n_alleles[m]; q++)
        {
          GS_ALLELE_TYPE z = d->allele_value[d->allele_offset[m] + q];
          if (z == x)
            fx = 1;
          if (z == y)
            fy = 1;
        }
        if (!fx || !fy)
          return 0;
      }
      else if (a == b)
      {
        if (c || (x != a) || (y != b))
          return 0;
      }
      else if (c)
      {
        if (!(((x == a) && (y == b)) || ((x == b) && (y == a))))
          return 0;
      }
      else if ((x != a) || (y != b))
        return 0;
    }
  return (claimed == stage->n_written) && (covered == stage->n_covered);
}

static int st_phase_step0_stage_apply_imputation(
    const st_phase_data *d, size_t individual, size_t chromosome,
    const st_phase_step0_imputation_result *imp,
    st_phase_dataset_stage *stage)
{
  size_t start, len, new_writes = 0u;
  if (!st_phase_data_adapter_valid(d) ||
      !st_phase_step0_stage_structurally_valid(d, stage) ||
      (NULL == imp) || (individual >= d->n_ind) ||
      (chromosome >= d->n_chrom) || (stage->n_ind != d->n_ind) ||
      (stage->n_marker != d->n_marker))
    return -1;
  start = d->chrom_start[chromosome];
  len = d->chrom_length[chromosome];
  if ((imp->n_marker != len) || (NULL == imp->is_missing) ||
      (NULL == imp->imputed) || (NULL == imp->phased) ||
      (NULL == imp->allele1) || (NULL == imp->allele2) ||
      (NULL == imp->genotype_probability) ||
      (NULL == imp->phase_probability))
    return -1;

  for (size_t j = 0u; j < len; j++)
    if (imp->imputed[j])
    {
      size_t m = start + j,
             raw = (size_t)d->source_marker[m] * d->n_ind + individual;
      size_t off = d->allele_offset[m];
      uint8_t x = imp->allele1[j], y = imp->allele2[j];
      if (!imp->is_missing[j] || (x >= d->n_alleles[m]) ||
          (y >= d->n_alleles[m]) || !stage->covered[raw] ||
          stage->claimed[raw] ||
          (GS_ALLELE_MISSING != stage->base1[raw]) ||
          (GS_ALLELE_MISSING != stage->base2[raw]) ||
          (GS_ALLELE_MISSING != stage->staged1[raw]) ||
          (GS_ALLELE_MISSING != stage->staged2[raw]) ||
          !isfinite(imp->genotype_probability[j]) ||
          !isfinite(imp->phase_probability[j]) ||
          (imp->genotype_probability[j] < 0.0) ||
          (imp->genotype_probability[j] > 1.0) ||
          (imp->phase_probability[j] < 0.0) ||
          (imp->phase_probability[j] > 1.0) ||
          (0u == d->allele_value[off + (size_t)x]) ||
          (0u == d->allele_value[off + (size_t)y]) ||
          (GS_ALLELE_MISSING == d->allele_value[off + (size_t)x]) ||
          (GS_ALLELE_MISSING == d->allele_value[off + (size_t)y]))
        return -1;
      new_writes++;
    }
    else if (imp->is_missing[j] &&
             ((ST_PHASE_STATE_MISSING != imp->allele1[j]) ||
              (ST_PHASE_STATE_MISSING != imp->allele2[j])))
      return -1;

  if (new_writes > stage->n_genotype - stage->n_written)
    return -1;
  for (size_t j = 0u; j < len; j++)
    if (imp->imputed[j])
    {
      size_t m = start + j,
             raw = (size_t)d->source_marker[m] * d->n_ind + individual;
      size_t off = d->allele_offset[m];
      stage->staged1[raw] = d->allele_value[off + (size_t)imp->allele1[j]];
      stage->staged2[raw] = d->allele_value[off + (size_t)imp->allele2[j]];
      stage->claimed[raw] = 1u;
    }
  stage->n_written += new_writes;
  return 0;
}

typedef struct
{
  st_phase_stage_audit audit;
  size_t n_imputed;
  size_t n_imputed_homozygous;
  size_t n_imputed_heterozygous;
  size_t n_source_missing_unresolved;
  size_t n_observed_phase_written;
} st_phase_step0_stage_audit;

static void st_phase_step0_stage_audit_init(st_phase_step0_stage_audit *x)
{
  if (NULL != x)
  {
    memset(x, 0, sizeof(*x));
    st_phase_stage_audit_init(&x->audit);
  }
}

static int st_phase_step0_stage_audit_readonly(
    const gs_varset_type *gv, const st_phase_data *d,
    const st_phase_dataset_stage *stage, st_phase_step0_stage_audit *out)
{
  st_phase_step0_stage_audit tmp;
  size_t n = 0u, covered = 0u, written = 0u;
  st_phase_step0_stage_audit_init(&tmp);
  if ((NULL == gv) || (NULL == d) || (NULL == stage) || (NULL == out) ||
      !st_phase_data_adapter_valid(d) || (stage->n_ind != d->n_ind) ||
      (stage->n_marker != d->n_marker) ||
      (0 != st_phase_size_mul(d->n_ind, d->n_marker, &n)) ||
      (stage->n_genotype != n) || (NULL == stage->base1) ||
      (NULL == stage->base2) || (NULL == stage->staged1) ||
      (NULL == stage->staged2) || (NULL == stage->claimed) ||
      (NULL == stage->covered))
    return ST_PHASE_AUDIT_INVALID;

  tmp.audit.n_ind = d->n_ind;
  tmp.audit.n_marker = d->n_marker;
  tmp.audit.n_genotype = n;
  for (size_t m = 0u; m < d->n_marker; m++)
    for (size_t i = 0u; i < d->n_ind; i++)
    {
      size_t raw = (size_t)d->source_marker[m] * d->n_ind + i;
      GS_ALLELE_TYPE a = stage->base1[raw], b = stage->base2[raw];
      GS_ALLELE_TYPE x = stage->staged1[raw], y = stage->staged2[raw];
      unsigned char claim = stage->claimed[raw], cov = stage->covered[raw];
      if ((claim > 1u) || (cov > 1u) || (claim && !cov))
        return ST_PHASE_AUDIT_CORRUPT;
      if (cov)
      {
        covered++;
        tmp.audit.n_covered++;
      }
      if (claim)
      {
        written++;
        tmp.audit.n_written++;
      }

      if ((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b))
      {
        int foundx = 0, foundy = 0;
        if (!((GS_ALLELE_MISSING == a) && (GS_ALLELE_MISSING == b)))
          return ST_PHASE_AUDIT_CORRUPT;
        tmp.audit.n_missing++;
        if (cov)
          tmp.audit.n_covered_missing++;
        if (!claim)
        {
          if ((GS_ALLELE_MISSING != x) || (GS_ALLELE_MISSING != y))
            return ST_PHASE_AUDIT_CORRUPT;
          tmp.n_source_missing_unresolved++;
          continue;
        }
        if ((GS_ALLELE_MISSING == x) || (GS_ALLELE_MISSING == y) ||
            (0u == x) || (0u == y) || (0u == d->n_alleles[m]) ||
            (d->n_alleles[m] > 2u))
          return ST_PHASE_AUDIT_CORRUPT;
        for (size_t q = 0u; q < d->n_alleles[m]; q++)
        {
          GS_ALLELE_TYPE v = d->allele_value[d->allele_offset[m] + q];
          if (v == x)
            foundx = 1;
          if (v == y)
            foundy = 1;
        }
        if (!foundx || !foundy)
          return ST_PHASE_AUDIT_CORRUPT;
        tmp.n_imputed++;
        if (x == y)
          tmp.n_imputed_homozygous++;
        else
          tmp.n_imputed_heterozygous++;
        continue;
      }

      if ((0u == a) || (0u == b) || (GS_ALLELE_MISSING == x) ||
          (GS_ALLELE_MISSING == y) || (0u == x) || (0u == y))
        return ST_PHASE_AUDIT_CORRUPT;
      if (a == b)
      {
        if (claim || (x != a) || (y != a))
          return ST_PHASE_AUDIT_CORRUPT;
        tmp.audit.n_homozygous++;
        if (cov)
          tmp.audit.n_covered_homozygous++;
      }
      else
      {
        tmp.audit.n_heterozygous++;
        if (cov)
          tmp.audit.n_covered_heterozygous++;
        if (claim)
        {
          tmp.n_observed_phase_written++;
          if ((x == a) && (y == b))
            tmp.audit.n_written_same_source_order++;
          else if ((x == b) && (y == a))
            tmp.audit.n_written_swapped_source_order++;
          else
            return ST_PHASE_AUDIT_CORRUPT;
        }
        else
        {
          if ((x != a) || (y != b))
            return ST_PHASE_AUDIT_CORRUPT;
          if (cov)
            tmp.audit.n_covered_heterozygous_unresolved++;
        }
      }
    }

  if ((covered != stage->n_covered) || (written != stage->n_written) ||
      (tmp.audit.n_missing + tmp.audit.n_homozygous +
           tmp.audit.n_heterozygous !=
       n) ||
      (tmp.audit.n_covered_missing + tmp.audit.n_covered_homozygous +
           tmp.audit.n_covered_heterozygous !=
       tmp.audit.n_covered) ||
      (tmp.n_imputed_homozygous + tmp.n_imputed_heterozygous !=
       tmp.n_imputed) ||
      (tmp.n_observed_phase_written !=
       tmp.audit.n_written_same_source_order +
           tmp.audit.n_written_swapped_source_order) ||
      (tmp.audit.n_written !=
       tmp.n_observed_phase_written + tmp.n_imputed))
    return ST_PHASE_AUDIT_CORRUPT;

  tmp.audit.n_uncovered = n - tmp.audit.n_covered;
  tmp.audit.n_covered_unwritten =
      tmp.audit.n_covered - tmp.audit.n_written;
  tmp.audit.complete = (tmp.audit.n_covered == n) ? 1 : 0;
  tmp.audit.source_bytes_current = 0;
  {
    size_t bytes = 0u;
    if ((NULL != gv->All1) && (NULL != gv->All2) &&
        ((size_t)gv->NoInd == stage->n_ind) &&
        ((size_t)gv->NoMar == stage->n_marker) &&
        (0 == st_phase_size_mul(stage->n_genotype, sizeof(GS_ALLELE_TYPE),
                                &bytes)) &&
        (0 == memcmp(stage->base1, gv->All1, bytes)) &&
        (0 == memcmp(stage->base2, gv->All2, bytes)))
      tmp.audit.source_bytes_current = 1;
  }
  tmp.audit.workspace_current =
      st_phase_stage_workspace_matches_source(gv, d) ? 1 : 0;
  *out = tmp;
  if (!tmp.audit.source_bytes_current || !tmp.audit.workspace_current)
    return ST_PHASE_AUDIT_STALE;
  if (!tmp.audit.complete)
    return ST_PHASE_AUDIT_INCOMPLETE;
  return ST_PHASE_AUDIT_OK;
}

static int
st_phase_step0_dataset_stage_commit(gs_varset_type *gv,
                                    const st_phase_data *d,
                                    const st_phase_dataset_stage *stage)
{
  st_phase_step0_stage_audit audit;
  size_t bytes;
  int astatus;
  if ((NULL == gv) || (NULL == d) || (NULL == stage) ||
      (0 != st_phase_size_mul(stage->n_genotype, sizeof(GS_ALLELE_TYPE),
                              &bytes)))
    return ST_PHASE_COMMIT_INVALID;
  astatus = st_phase_step0_stage_audit_readonly(gv, d, stage, &audit);
  if (ST_PHASE_AUDIT_OK != astatus)
    return st_phase_commit_from_audit_status(astatus);
  if (st_phase_commit_test_should_fail(1))
    return ST_PHASE_COMMIT_ROLLED_BACK;
  memcpy(gv->All1, stage->staged1, bytes);
  if (st_phase_commit_test_should_fail(2))
    goto rollback;
  memcpy(gv->All2, stage->staged2, bytes);
  if (st_phase_commit_test_should_fail(3))
    goto rollback;
  if ((0 != memcmp(gv->All1, stage->staged1, bytes)) ||
      (0 != memcmp(gv->All2, stage->staged2, bytes)))
    goto rollback;
  return ST_PHASE_COMMIT_OK;
rollback:
  memcpy(gv->All1, stage->base1, bytes);
  memcpy(gv->All2, stage->base2, bytes);
  return ST_PHASE_COMMIT_ROLLED_BACK;
}

typedef struct
{
  st_phase_chromosome_result phase;
  st_phase_step0_imputation_result imputation;
} st_phase_step0_chromosome_result;

static void
st_phase_step0_chromosome_result_init(st_phase_step0_chromosome_result *x)
{
  if (NULL != x)
  {
    memset(x, 0, sizeof(*x));
    st_phase_chromosome_result_init(&x->phase);
    st_phase_step0_imputation_result_init(&x->imputation);
  }
}
static void
st_phase_step0_chromosome_result_free(st_phase_step0_chromosome_result *x)
{
  if (NULL == x)
    return;
  st_phase_chromosome_result_free(&x->phase);
  st_phase_step0_imputation_result_free(&x->imputation);
  st_phase_step0_chromosome_result_init(x);
}

static int st_phase_step0_run_chromosome_reference(
    const st_phase_data *d, size_t individual, size_t chromosome,
    const uint8_t *reference_hap, size_t n_ref_hap, size_t max_states,
    const double *initial_pair_weight, size_t initial_pair_count,
    double error_rate, double break_rate1_per_morgan,
    double break_rate2_per_morgan, double confidence_threshold,
    int impute_requested, double imputation_threshold,
    st_phase_step0_chromosome_result *out)
{
  st_phase_step0_chromosome_result tmp;
  st_phase_step0_hmm_work w;
  st_phase_orientation_result orientation;
  st_phase_hard_call_result hard;
  st_phase_proposed_output proposal;
  int status = -1;
  st_phase_step0_chromosome_result_init(&tmp);
  st_phase_step0_hmm_work_init(&w);
  st_phase_orientation_result_init(&orientation);
  st_phase_hard_call_result_init(&hard);
  st_phase_proposed_output_init(&proposal);
  if ((NULL == out) ||
      ((0 != impute_requested) && (1 != impute_requested)) ||
      !st_phase_hard_threshold_valid(confidence_threshold) ||
      !st_phase_step0_imputation_threshold_valid(imputation_threshold))
    goto cleanup;
  if (0 != st_phase_step0_hmm_build(
               d, individual, chromosome, reference_hap, n_ref_hap,
               max_states, initial_pair_weight, initial_pair_count,
               error_rate, break_rate1_per_morgan, break_rate2_per_morgan,
               &w))
    goto cleanup;
  if (0 != st_phase_orientation_posterior(
               &w.pair, &w.hmm, w.input.observed1, w.input.observed2,
               w.phase12, w.phase21, &orientation))
    goto cleanup;
  if (0 != st_phase_hard_call_from_orientation(
               &orientation, confidence_threshold, &hard))
    goto cleanup;
  if (0 != st_phase_step0_complete_reference_calls_viterbi(&w, &hard))
    goto cleanup;
  if (0 != st_phase_build_proposed_output(d, &w.input, &hard, &proposal))
    goto cleanup;
  if (impute_requested &&
      (0 != st_phase_step0_decode_missing(
                d, reference_hap, n_ref_hap, &w, confidence_threshold,
                imputation_threshold, &tmp.imputation)))
    goto cleanup;

  tmp.phase.individual = individual;
  tmp.phase.chromosome = chromosome;
  tmp.phase.n_marker = w.M;
  tmp.phase.n_ref_hap = n_ref_hap;
  tmp.phase.possible = w.hmm.possible;
  tmp.phase.log_likelihood = w.hmm.log_likelihood;
  tmp.phase.viterbi_log_probability = w.hmm.viterbi_log_probability;
  tmp.phase.orientation = orientation;
  st_phase_orientation_result_init(&orientation);
  tmp.phase.hard = hard;
  st_phase_hard_call_result_init(&hard);
  tmp.phase.proposal = proposal;
  st_phase_proposed_output_init(&proposal);
  st_phase_step0_chromosome_result_free(out);
  *out = tmp;
  st_phase_step0_chromosome_result_init(&tmp);
  status = 0;
cleanup:
  st_phase_proposed_output_free(&proposal);
  st_phase_hard_call_result_free(&hard);
  st_phase_orientation_result_free(&orientation);
  st_phase_step0_hmm_work_free(&w);
  st_phase_step0_chromosome_result_free(&tmp);
  return status;
}

typedef struct
{
  size_t n_request, n_possible, n_impossible, n_requested_markers;
  size_t n_missing_requested, n_imputed, n_imputation_unresolved;
  size_t n_imputed_hom0, n_imputed_hom1, n_imputed_heterozygous;
  size_t n_imputed_heterozygous_phased, n_imputed_heterozygous_unphased;
  int audit_status;
  st_phase_dataset_stage stage;
  st_phase_step0_stage_audit audit;
} st_phase_step0_dataset_result;

static void
st_phase_step0_dataset_result_init(st_phase_step0_dataset_result *x)
{
  if (NULL != x)
  {
    memset(x, 0, sizeof(*x));
    x->audit_status = ST_PHASE_AUDIT_INVALID;
    st_phase_dataset_stage_init(&x->stage);
    st_phase_step0_stage_audit_init(&x->audit);
  }
}
static void
st_phase_step0_dataset_result_free(st_phase_step0_dataset_result *x)
{
  if (NULL == x)
    return;
  st_phase_dataset_stage_free(&x->stage);
  st_phase_step0_dataset_result_init(x);
}

static int st_phase_step0_run_dataset_reference(
    const gs_varset_type *gv, const st_phase_data *d,
    const st_phase_dataset_request *request, size_t n_request,
    const uint8_t *reference_hap, size_t reference_n_marker,
    size_t n_ref_hap, size_t max_states, const double *initial_pair_weight,
    size_t initial_pair_count, double error_rate,
    double break_rate1_per_morgan, double break_rate2_per_morgan,
    double confidence_threshold, int impute_requested,
    double imputation_threshold, st_phase_step0_dataset_result *out)
{
  st_phase_step0_dataset_result tmp;
  st_phase_step0_imputation_result *imp_by_request = NULL;
  size_t requested_markers = 0u, ref_cells = 0u, pair_count = 0u;
  int status = ST_PHASE_DATASET_RUN_INVALID, abort_merge = 0;
  st_phase_step0_dataset_result_init(&tmp);
  if ((NULL == gv) || (NULL == out) || (NULL == reference_hap) ||
      (NULL == initial_pair_weight) || (0u == n_ref_hap) ||
      (0u == max_states) || !st_phase_data_adapter_valid(d) ||
      (reference_n_marker != d->n_marker) ||
      ((0 != impute_requested) && (1 != impute_requested)) ||
      !st_phase_step0_imputation_threshold_valid(imputation_threshold) ||
      (0 !=
       st_phase_size_mul(reference_n_marker, n_ref_hap, &ref_cells)) ||
      (0 != st_phase_size_mul(n_ref_hap, n_ref_hap, &pair_count)) ||
      (initial_pair_count != pair_count) ||
      !st_phase_dataset_requests_valid(d, request, n_request,
                                       &requested_markers))
    goto cleanup;
  for (size_t z = 0u; z < ref_cells; z++)
    if (reference_hap[z] > 1u)
      goto cleanup;
  if (impute_requested)
  {
    if (n_request > SIZE_MAX / sizeof(*imp_by_request))
      goto cleanup;
    imp_by_request = calloc(n_request, sizeof(*imp_by_request));
    if (NULL == imp_by_request)
      goto cleanup;
    for (size_t r = 0u; r < n_request; r++)
      st_phase_step0_imputation_result_init(imp_by_request + r);
  }
  if (ST_PHASE_STAGE_OK != st_phase_dataset_stage_build(gv, d, &tmp.stage))
    goto cleanup;
  tmp.n_request = n_request;
  tmp.n_requested_markers = requested_markers;

#ifdef MT
#pragma omp parallel for schedule(static, 1) ordered
#endif
  for (size_t r = 0u; r < n_request; r++)
  {
    st_phase_step0_chromosome_result chr;
    size_t c = request[r].chromosome, start = d->chrom_start[c],
           len = d->chrom_length[c];
    size_t ref_offset = 0u;
    int run_failed = 0;
    st_phase_step0_chromosome_result_init(&chr);
    if ((start > reference_n_marker) ||
        (len > reference_n_marker - start) ||
        (0 != st_phase_size_mul(start, n_ref_hap, &ref_offset)))
      run_failed = 1;
    else if (0 != st_phase_step0_run_chromosome_reference(
                      d, request[r].individual, c,
                      reference_hap + ref_offset, n_ref_hap, max_states,
                      initial_pair_weight, initial_pair_count, error_rate,
                      break_rate1_per_morgan, break_rate2_per_morgan,
                      confidence_threshold, impute_requested,
                      imputation_threshold, &chr))
      run_failed = 1;
#ifdef MT
#pragma omp ordered
#endif
    {
      if (!abort_merge)
      {
        if (run_failed)
        {
          status = ST_PHASE_DATASET_RUN_PHASING_FAILED;
          abort_merge = 1;
        }
        else
        {
          if (chr.phase.possible)
            tmp.n_possible++;
          else
            tmp.n_impossible++;
          if (ST_PHASE_STAGE_OK !=
              st_phase_dataset_stage_merge(gv, d, &chr.phase.proposal,
                                           &tmp.stage))
          {
            status = ST_PHASE_DATASET_RUN_MERGE_FAILED;
            abort_merge = 1;
          }
          else if (impute_requested)
          {
            imp_by_request[r] = chr.imputation;
            st_phase_step0_imputation_result_init(&chr.imputation);
            tmp.n_missing_requested += imp_by_request[r].n_missing;
            tmp.n_imputed += imp_by_request[r].n_imputed;
            tmp.n_imputation_unresolved += imp_by_request[r].n_unresolved;
            tmp.n_imputed_hom0 += imp_by_request[r].n_imputed_hom0;
            tmp.n_imputed_hom1 += imp_by_request[r].n_imputed_hom1;
            tmp.n_imputed_heterozygous +=
                imp_by_request[r].n_imputed_heterozygous;
            tmp.n_imputed_heterozygous_phased +=
                imp_by_request[r].n_imputed_heterozygous_phased;
            tmp.n_imputed_heterozygous_unphased +=
                imp_by_request[r].n_imputed_heterozygous_unphased;
          }
        }
      }
    }
    st_phase_step0_chromosome_result_free(&chr);
  }
  if (abort_merge)
    goto cleanup;
  if (impute_requested)
    for (size_t r = 0u; r < n_request; r++)
      if (0 != st_phase_step0_stage_apply_imputation(
                   d, request[r].individual, request[r].chromosome,
                   imp_by_request + r, &tmp.stage))
      {
        status = ST_PHASE_DATASET_RUN_MERGE_FAILED;
        goto cleanup;
      }

  tmp.audit_status =
      st_phase_step0_stage_audit_readonly(gv, d, &tmp.stage, &tmp.audit);
  if ((ST_PHASE_AUDIT_OK != tmp.audit_status) &&
      (ST_PHASE_AUDIT_INCOMPLETE != tmp.audit_status))
  {
    status = ST_PHASE_DATASET_RUN_AUDIT_FAILED;
    goto cleanup;
  }
  if ((tmp.n_possible + tmp.n_impossible != tmp.n_request) ||
      (tmp.stage.n_covered != tmp.n_requested_markers) ||
      (tmp.audit.audit.n_covered != tmp.n_requested_markers) ||
      (tmp.n_imputed + tmp.n_imputation_unresolved !=
       tmp.n_missing_requested))
  {
    status = ST_PHASE_DATASET_RUN_AUDIT_FAILED;
    goto cleanup;
  }
  st_phase_step0_dataset_result_free(out);
  *out = tmp;
  st_phase_step0_dataset_result_init(&tmp);
  status = ST_PHASE_DATASET_RUN_OK;
cleanup:
  if (NULL != imp_by_request)
  {
    for (size_t r = 0u; r < n_request; r++)
      st_phase_step0_imputation_result_free(imp_by_request + r);
    free(imp_by_request);
  }
  st_phase_step0_dataset_result_free(&tmp);
  return status;
}

typedef struct
{
  st_phase_pipeline_diagnostics base;
  int impute_requested;
  double imputation_threshold;
  size_t n_missing_requested;
  size_t n_imputed;
  size_t n_imputation_unresolved;
  size_t n_imputed_hom0;
  size_t n_imputed_hom1;
  size_t n_imputed_heterozygous;
  size_t n_imputed_heterozygous_phased;
  size_t n_imputed_heterozygous_unphased;
} st_phase_step0_pipeline_diagnostics;

static void st_phase_step0_pipeline_diagnostics_init(
    st_phase_step0_pipeline_diagnostics *x)
{
  if (NULL != x)
  {
    memset(x, 0, sizeof(*x));
    st_phase_pipeline_diagnostics_init(&x->base);
  }
}

int st_phase_dataset_pipeline_binary_reference_bounded_factorized_step0(
    gs_varset_type *gv, const st_phase_dataset_request *request,
    size_t n_request, const uint8_t *reference_hap,
    size_t reference_n_marker, size_t n_ref_hap, size_t max_states,
    const double *initial_pair_weight, size_t initial_pair_count,
    double error_rate, double break_rate1_per_morgan,
    double break_rate2_per_morgan, double confidence_threshold,
    int impute_requested, double imputation_threshold,
    int commit_requested, st_phase_step0_pipeline_diagnostics *diagnostics)
{
  st_phase_step0_pipeline_diagnostics diag;
  st_phase_data d;
  st_phase_step0_dataset_result result;
  int status = ST_PHASE_PIPELINE_INVALID,
      run_status = ST_PHASE_DATASET_RUN_INVALID;
  st_phase_step0_pipeline_diagnostics_init(&diag);
  st_phase_data_init(&d);
  st_phase_step0_dataset_result_init(&result);
  diag.base.commit_requested = commit_requested;
  diag.impute_requested = impute_requested;
  diag.imputation_threshold = imputation_threshold;
  if (NULL == diagnostics)
    goto cleanup;
  if (!st_phase_pipeline_parameters_valid(
          gv, request, n_request, reference_hap, reference_n_marker,
          n_ref_hap, max_states, initial_pair_weight, initial_pair_count,
          error_rate, break_rate1_per_morgan, break_rate2_per_morgan,
          confidence_threshold, commit_requested) ||
      ((0 != impute_requested) && (1 != impute_requested)) ||
      !st_phase_step0_imputation_threshold_valid(imputation_threshold))
    goto cleanup;
  if (0 != st_phase_build_data(gv, &d))
  {
    status = ST_PHASE_PIPELINE_DATA_FAILED;
    goto cleanup;
  }
  diag.base.n_ind = d.n_ind;
  diag.base.n_marker = d.n_marker;
  diag.base.n_chrom = d.n_chrom;
  (void)st_phase_size_mul(d.n_ind, d.n_marker, &diag.base.n_genotype);
  run_status = st_phase_step0_run_dataset_reference(
      gv, &d, request, n_request, reference_hap, reference_n_marker,
      n_ref_hap, max_states, initial_pair_weight, initial_pair_count,
      error_rate, break_rate1_per_morgan, break_rate2_per_morgan,
      confidence_threshold, impute_requested, imputation_threshold,
      &result);
  diag.base.dataset_run_status = run_status;
  if (ST_PHASE_DATASET_RUN_OK != run_status)
  {
    status = ST_PHASE_PIPELINE_RUN_FAILED;
    goto cleanup;
  }
  diag.base.n_request = result.n_request;
  diag.base.n_possible = result.n_possible;
  diag.base.n_impossible = result.n_impossible;
  diag.base.n_requested_markers = result.n_requested_markers;
  diag.base.audit_status = result.audit_status;
  diag.base.n_covered = result.audit.audit.n_covered;
  diag.base.n_uncovered = result.audit.audit.n_uncovered;
  diag.base.n_written = result.audit.audit.n_written;
  diag.base.n_covered_unwritten = result.audit.audit.n_covered_unwritten;
  diag.base.n_changed_genotypes =
      st_phase_dataset_stage_changed_genotypes(&result.stage);
  diag.base.n_missing = result.audit.audit.n_missing;
  diag.base.n_homozygous = result.audit.audit.n_homozygous;
  diag.base.n_heterozygous = result.audit.audit.n_heterozygous;
  diag.base.n_covered_missing = result.audit.audit.n_covered_missing;
  diag.base.n_covered_homozygous = result.audit.audit.n_covered_homozygous;
  diag.base.n_covered_heterozygous =
      result.audit.audit.n_covered_heterozygous;
  diag.base.n_covered_heterozygous_unresolved =
      result.audit.audit.n_covered_heterozygous_unresolved;
  diag.base.n_written_same_source_order =
      result.audit.audit.n_written_same_source_order;
  diag.base.n_written_swapped_source_order =
      result.audit.audit.n_written_swapped_source_order;
  diag.base.audit_complete = result.audit.audit.complete;
  diag.base.audit_source_bytes_current =
      result.audit.audit.source_bytes_current;
  diag.base.audit_workspace_current = result.audit.audit.workspace_current;
  diag.n_missing_requested = result.n_missing_requested;
  diag.n_imputed = result.n_imputed;
  diag.n_imputation_unresolved = result.n_imputation_unresolved;
  diag.n_imputed_hom0 = result.n_imputed_hom0;
  diag.n_imputed_hom1 = result.n_imputed_hom1;
  diag.n_imputed_heterozygous = result.n_imputed_heterozygous;
  diag.n_imputed_heterozygous_phased =
      result.n_imputed_heterozygous_phased;
  diag.n_imputed_heterozygous_unphased =
      result.n_imputed_heterozygous_unphased;
  if (commit_requested)
  {
    diag.base.commit_attempted = 1;
    diag.base.commit_status =
        st_phase_step0_dataset_stage_commit(gv, &d, &result.stage);
    if (ST_PHASE_COMMIT_OK != diag.base.commit_status)
    {
      status = ST_PHASE_PIPELINE_COMMIT_FAILED;
      goto cleanup;
    }
    diag.base.committed = 1;
  }
  status = ST_PHASE_PIPELINE_OK;
cleanup:
  diag.base.pipeline_status = status;
  if (NULL != diagnostics)
    *diagnostics = diag;
  st_phase_step0_dataset_result_free(&result);
  st_phase_data_free(&d);
  return status;
}

typedef struct
{
  st_phase_cohort_pipeline_diagnostics base;
  int impute_requested;
  double imputation_threshold;
  size_t n_missing_requested, n_imputed, n_imputation_unresolved;
  size_t n_imputed_hom0, n_imputed_hom1, n_imputed_heterozygous;
  size_t n_imputed_heterozygous_phased, n_imputed_heterozygous_unphased;
} st_phase_step0_cohort_pipeline_diagnostics;

static SEXP
st_phase_step0_r_diagnostics(const st_phase_step0_pipeline_diagnostics *d)
{
  static const char *extra_names[] = {"impute",
                                      "imputation.threshold",
                                      "n.imputation.missing.requested",
                                      "n.imputed",
                                      "n.imputation.unresolved",
                                      "n.imputed.hom0",
                                      "n.imputed.hom1",
                                      "n.imputed.heterozygous",
                                      "n.imputed.heterozygous.phased",
                                      "n.imputed.heterozygous.unphased"};
  enum
  {
    NEXTRA = (int)(sizeof(extra_names) / sizeof(extra_names[0]))
  };
  SEXP base = PROTECT(st_phase_r_diagnostics(&d->base));
  SEXP base_names = getAttrib(base, R_NamesSymbol);
  R_xlen_t nb = XLENGTH(base);
  SEXP ans = PROTECT(allocVector(VECSXP, nb + NEXTRA));
  SEXP nms = PROTECT(allocVector(STRSXP, nb + NEXTRA));
  for (R_xlen_t k = 0; k < nb; k++)
  {
    SET_VECTOR_ELT(ans, k, VECTOR_ELT(base, k));
    SET_STRING_ELT(nms, k, STRING_ELT(base_names, k));
  }
#define ST_PHASE_STEP0_R_ADD(slot, name, value)                           \
  do                                                                      \
  {                                                                       \
    SET_VECTOR_ELT(ans, nb + (slot), (value));                            \
    SET_STRING_ELT(nms, nb + (slot), mkChar(name));                       \
  } while (0)
  ST_PHASE_STEP0_R_ADD(0, extra_names[0],
                       ScalarLogical(d->impute_requested));
  ST_PHASE_STEP0_R_ADD(1, extra_names[1],
                       ScalarReal(d->imputation_threshold));
  ST_PHASE_STEP0_R_ADD(2, extra_names[2],
                       st_phase_r_count(d->n_missing_requested));
  ST_PHASE_STEP0_R_ADD(3, extra_names[3], st_phase_r_count(d->n_imputed));
  ST_PHASE_STEP0_R_ADD(4, extra_names[4],
                       st_phase_r_count(d->n_imputation_unresolved));
  ST_PHASE_STEP0_R_ADD(5, extra_names[5],
                       st_phase_r_count(d->n_imputed_hom0));
  ST_PHASE_STEP0_R_ADD(6, extra_names[6],
                       st_phase_r_count(d->n_imputed_hom1));
  ST_PHASE_STEP0_R_ADD(7, extra_names[7],
                       st_phase_r_count(d->n_imputed_heterozygous));
  ST_PHASE_STEP0_R_ADD(8, extra_names[8],
                       st_phase_r_count(d->n_imputed_heterozygous_phased));
  ST_PHASE_STEP0_R_ADD(
      9, extra_names[9],
      st_phase_r_count(d->n_imputed_heterozygous_unphased));
#undef ST_PHASE_STEP0_R_ADD
  setAttrib(ans, R_NamesSymbol, nms);
  UNPROTECT(3);
  return ans;
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

static int st_phase_match_names(char (*from)[gs_MAXNME], size_t n_from,
                                char (*to)[gs_MAXNME], size_t n_to,
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
    if (NULL == memchr(to[i], '\0', gs_MAXNME))
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
    if (NULL == memchr(from[i], '\0', gs_MAXNME))
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

static int st_switch_marker_alleles(const gs_varset_type *gv,
                                    size_t marker,
                                    GS_ALLELE_TYPE allele[2],
                                    size_t *n_allele)
{
  size_t n = 0u;

  if ((NULL == gv) || (NULL == allele) || (NULL == n_allele) ||
      (marker >= (size_t)gv->NoMar) || (NULL == gv->All1) ||
      (NULL == gv->All2))
    return -1;

  for (size_t i = 0u; i < (size_t)gv->NoInd; i++)
  {
    size_t idx = marker * (size_t)gv->NoInd + i;
    GS_ALLELE_TYPE x[2] = {gv->All1[idx], gv->All2[idx]};
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
    if ((GS_ALLELE_MISSING != gv0->All1[idx]) &&
        (gv0->All1[idx] != gv0->All2[idx]))
    {
      nhet++;
      if (unif_rand() < 0.5)
      {
        GS_ALLELE_TYPE z = gv0->All1[idx];
        gv0->All1[idx] = gv0->All2[idx];
        gv0->All2[idx] = z;
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
  if (!gs_marker_genotypes_valid(truth) ||
      !gs_marker_genotypes_valid(test))
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

  if ((0 != st_phase_match_names(truth->IndName, (size_t)truth->NoInd,
                                 test->IndName, (size_t)test->NoInd,
                                 &ind_match)) ||
      (0 != st_phase_match_names(truth->MarName, (size_t)truth->NoMar,
                                 test->MarName, (size_t)test->NoMar,
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
      GS_ALLELE_TYPE ta1 = truth->All1[idx], ta2 = truth->All2[idx];
      GS_ALLELE_TYPE xa1 = test->All1[idt], xa2 = test->All2[idt];
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

  PROTECT(ans = allocVector(VECSXP, 9));
  PROTECT(names = allocVector(STRSXP, 9));
#define ST_SW_SET(I, N, V)                                                \
  do                                                                      \
  {                                                                       \
    SET_STRING_ELT(names, (I), mkChar((N)));                              \
    SET_VECTOR_ELT(ans, (I), (V));                                        \
  } while (0)
  ST_SW_SET(0, "errors", ScalarReal((double)nerror));
  ST_SW_SET(1, "links", ScalarReal((double)nlink));
  ST_SW_SET(2, "switch.error",
            ScalarReal((0u == nlink) ? NA_REAL
                                     : (double)nerror / (double)nlink));
  ST_SW_SET(3, "switch.error.percent",
            ScalarReal((0u == nlink)
                           ? NA_REAL
                           : 100.0 * (double)nerror / (double)nlink));
  ST_SW_SET(4, "heterozygotes", ScalarReal((double)nhet));
  ST_SW_SET(5, "scored.heterozygotes", ScalarReal((double)nscored));
  ST_SW_SET(
      6, "coverage",
      ScalarReal((0u == nhet) ? NA_REAL : (double)nscored / (double)nhet));
  ST_SW_SET(7, "missing.heterozygotes", ScalarReal((double)nmissing));
  ST_SW_SET(8, "genotype.mismatches", ScalarReal((double)nmismatch));
#undef ST_SW_SET
  setAttrib(ans, R_NamesSymbol, names);
  UNPROTECT(2);
  return ans;
}

static int st_phase_step1_build_reference_dataset(
    const gs_varset_type *target_gv, const st_phase_data *target,
    const gs_varset_type *reference_gv, SEXP reference_individuals,
    uint8_t **reference_out, size_t *n_reference_haplotype_out,
    size_t *n_reference_individual_out, int *heterozygous_phased_out)
{
  st_phase_r_name_index *marker_index = NULL;
  size_t *selected = NULL;
  size_t *reference_marker = NULL;
  uint8_t *compact = NULL;
  size_t N = 0u, H = 0u, ncell = 0u, M;
  int has_heterozygous = 0;
  int status = -1;

  if ((NULL == target_gv) || (NULL == target) || (NULL == reference_gv) ||
      (NULL == reference_out) || (NULL == n_reference_haplotype_out) ||
      (NULL == n_reference_individual_out) ||
      (NULL == heterozygous_phased_out) ||
      (NULL == reference_gv->MarName) || (NULL == reference_gv->IndName) ||
      (NULL == reference_gv->All1) || (NULL == reference_gv->All2) ||
      (0u == reference_gv->NoInd) || (0u == reference_gv->NoMar))
    return -1;

  *reference_out = NULL;
  *n_reference_haplotype_out = 0u;
  *n_reference_individual_out = 0u;
  *heterozygous_phased_out = 0;

  M = target->n_marker;
  if ((size_t)reference_gv->NoMar != M)
    return -1;

  if (0 != st_phase_r_individuals(reference_gv, reference_individuals,
                                  &selected, &N))
    goto cleanup;
  if (0u == N)
    goto cleanup;

  marker_index = malloc(M * sizeof(*marker_index));
  reference_marker = malloc(M * sizeof(*reference_marker));
  if ((NULL == marker_index) || (NULL == reference_marker))
    goto cleanup;

  for (size_t r = 0u; r < M; r++)
  {
    const char *name = reference_gv->MarName[r];
    if ((NULL == name) || ('\0' == name[0]))
      goto cleanup;
    marker_index[r].name = name;
    marker_index[r].index = r;
  }
  qsort(marker_index, M, sizeof(*marker_index), st_phase_r_name_index_cmp);
  for (size_t r = 1u; r < M; r++)
    if (0 == strcmp(marker_index[r - 1u].name, marker_index[r].name))
      goto cleanup;

  for (size_t m = 0u; m < M; m++)
  {
    const char *want = target_gv->MarName[target->source_marker[m]];
    size_t lo = 0u, hi = M;
    if ((NULL == want) || ('\0' == want[0]))
      goto cleanup;
    while (lo < hi)
    {
      size_t mid = lo + (hi - lo) / 2u;
      int cmp = strcmp(marker_index[mid].name, want);
      if (cmp < 0)
        lo = mid + 1u;
      else
        hi = mid;
    }
    if ((lo >= M) || (0 != strcmp(marker_index[lo].name, want)))
      goto cleanup;
    reference_marker[m] = marker_index[lo].index;
  }

  for (size_t m = 0u; m < M; m++)
  {
    const size_t off = target->allele_offset[m];
    const uint8_t A = target->n_alleles[m];
    if ((A < 1u) || (A > 2u))
      goto cleanup;
    for (size_t j = 0u; j < N; j++)
    {
      const size_t raw =
          reference_marker[m] * (size_t)reference_gv->NoInd + selected[j];
      const GS_ALLELE_TYPE a = reference_gv->All1[raw];
      const GS_ALLELE_TYPE b = reference_gv->All2[raw];
      int found_a = 0, found_b = 0;
      if ((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b))
        goto cleanup;
      for (uint8_t q = 0u; q < A; q++)
      {
        if (target->allele_value[off + (size_t)q] == a)
          found_a = 1;
        if (target->allele_value[off + (size_t)q] == b)
          found_b = 1;
      }
      if (!found_a || !found_b)
        goto cleanup;
      if (a != b)
        has_heterozygous = 1;
    }
  }

  if (has_heterozygous && !reference_gv->phased)
  {
    status = -2;
    goto cleanup;
  }

  if (has_heterozygous)
  {
    if (N > SIZE_MAX / 2u)
      goto cleanup;
    H = 2u * N;
  }
  else
    H = N;
  if (0 != st_phase_size_mul(M, H, &ncell))
    goto cleanup;
  compact = malloc(ncell * sizeof(*compact));
  if (NULL == compact)
    goto cleanup;

  for (size_t m = 0u; m < M; m++)
  {
    const size_t off = target->allele_offset[m];
    const uint8_t A = target->n_alleles[m];
    for (size_t j = 0u; j < N; j++)
    {
      const size_t raw =
          reference_marker[m] * (size_t)reference_gv->NoInd + selected[j];
      const GS_ALLELE_TYPE a = reference_gv->All1[raw];
      const GS_ALLELE_TYPE b = reference_gv->All2[raw];
      uint8_t sa = ST_PHASE_STATE_MISSING, sb = ST_PHASE_STATE_MISSING;
      for (uint8_t q = 0u; q < A; q++)
      {
        if (target->allele_value[off + (size_t)q] == a)
          sa = q;
        if (target->allele_value[off + (size_t)q] == b)
          sb = q;
      }
      if ((ST_PHASE_STATE_MISSING == sa) || (ST_PHASE_STATE_MISSING == sb))
        goto cleanup;
      if (has_heterozygous)
      {
        compact[m * H + 2u * j] = sa;
        compact[m * H + 2u * j + 1u] = sb;
      }
      else
        compact[m * H + j] = sa;
    }
  }

  *reference_out = compact;
  *n_reference_haplotype_out = H;
  *n_reference_individual_out = N;
  *heterozygous_phased_out = has_heterozygous;
  compact = NULL;
  status = 0;

cleanup:
  free(marker_index);
  free(reference_marker);
  free(selected);
  free(compact);
  return status;
}

typedef enum
{
  ST_PHASE_POP_DH = 0,
  ST_PHASE_POP_F2,
  ST_PHASE_POP_F3,
  ST_PHASE_POP_F4,
  ST_PHASE_POP_S2,
  ST_PHASE_POP_S3,
  ST_PHASE_POP_S4,
  ST_PHASE_POP_SSD,
  ST_PHASE_POP_FINFTY
} st_phase_population_type;

typedef struct
{
  st_phase_population_type population_type;
  double expected_breaks_per_morgan;
  int breaks_estimated;
  double break_rate1_per_morgan;
  double break_rate2_per_morgan;
  double error_rate;
  size_t max_states;
} st_phase_step1_configuration;

static int st_phase_step1_ascii_equal_nocase(const char *a, const char *b)
{
  unsigned char x, y;
  if ((NULL == a) || (NULL == b))
    return 0;
  while (('\0' != *a) && ('\0' != *b))
  {
    x = (unsigned char)*a++;
    y = (unsigned char)*b++;
    if ((x >= 'a') && (x <= 'z'))
      x = (unsigned char)(x - 'a' + 'A');
    if ((y >= 'a') && (y <= 'z'))
      y = (unsigned char)(y - 'a' + 'A');
    if (x != y)
      return 0;
  }
  return ('\0' == *a) && ('\0' == *b);
}

static const char *
st_phase_step1_population_name(st_phase_population_type x)
{
  switch (x)
  {
  case ST_PHASE_POP_DH:
    return "DH";
  case ST_PHASE_POP_F2:
    return "F2";
  case ST_PHASE_POP_F3:
    return "F3";
  case ST_PHASE_POP_F4:
    return "F4";
  case ST_PHASE_POP_S2:
    return "S2";
  case ST_PHASE_POP_S3:
    return "S3";
  case ST_PHASE_POP_S4:
    return "S4";
  case ST_PHASE_POP_SSD:
    return "SSD";
  case ST_PHASE_POP_FINFTY:
    return "Finfty";
  default:
    return "unknown";
  }
}

static int st_phase_step1_population_parse(SEXP x,
                                           st_phase_population_type *out)
{
  const char *z;
  if ((NULL == out) || (R_NilValue == x) || (STRSXP != TYPEOF(x)) ||
      (1 != XLENGTH(x)) || (NA_STRING == STRING_ELT(x, 0)))
    return 0;
  z = CHAR(STRING_ELT(x, 0));
  if (st_phase_step1_ascii_equal_nocase(z, "DH"))
    *out = ST_PHASE_POP_DH;
  else if (st_phase_step1_ascii_equal_nocase(z, "F2"))
    *out = ST_PHASE_POP_F2;
  else if (st_phase_step1_ascii_equal_nocase(z, "F3"))
    *out = ST_PHASE_POP_F3;
  else if (st_phase_step1_ascii_equal_nocase(z, "F4"))
    *out = ST_PHASE_POP_F4;
  else if (st_phase_step1_ascii_equal_nocase(z, "S2"))
    *out = ST_PHASE_POP_S2;
  else if (st_phase_step1_ascii_equal_nocase(z, "S3"))
    *out = ST_PHASE_POP_S3;
  else if (st_phase_step1_ascii_equal_nocase(z, "S4"))
    *out = ST_PHASE_POP_S4;
  else if (st_phase_step1_ascii_equal_nocase(z, "SSD"))
    *out = ST_PHASE_POP_SSD;
  else if (st_phase_step1_ascii_equal_nocase(z, "FINFTY"))
    *out = ST_PHASE_POP_FINFTY;
  else
    return 0;
  return 1;
}

static double st_phase_step1_population_breaks(st_phase_population_type x)
{
  switch (x)
  {
  case ST_PHASE_POP_DH:
    return 1.0;
  case ST_PHASE_POP_F2:
    return 1.0;
  case ST_PHASE_POP_F3:
    return 1.5;
  case ST_PHASE_POP_F4:
    return 2.0;
  case ST_PHASE_POP_S2:
    return 1.5;
  case ST_PHASE_POP_S3:
    return 1.75;
  case ST_PHASE_POP_S4:
    return 1.875;
  case ST_PHASE_POP_SSD:
    return 2.0;
  case ST_PHASE_POP_FINFTY:
    return NAN;
  default:
    return NAN;
  }
}

static int st_phase_step1_size_default(SEXP x, size_t value, size_t *out)
{
  if (NULL == out)
    return 0;
  if (R_NilValue == x)
  {
    *out = value;
    return 1;
  }
  return st_phase_r_scalar_size(x, out);
}

static int st_phase_step1_double_default(SEXP x, double value, double *out)
{
  if (NULL == out)
    return 0;
  if (R_NilValue == x)
  {
    *out = value;
    return 1;
  }
  return st_phase_r_scalar_double(x, out);
}

static int st_phase_step1_default_prior(size_t R, double **out)
{
  size_t ncell;
  double *x;
  if ((NULL == out) || (0u == R) || (0 != st_phase_size_mul(R, R, &ncell)))
    return -1;
  x = calloc(ncell, sizeof(*x));
  if (NULL == x)
    return -1;
  for (size_t i = 0u; i < R; i++)
    for (size_t j = i; j < R; j++)
      x[i * R + j] = 1.0;
  *out = x;
  return 0;
}

static int st_phase_step1_double_cmp(const void *a, const void *b)
{
  const double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

static double st_phase_step1_median(double *x, size_t n)
{
  if ((NULL == x) || (0u == n))
    return NAN;
  qsort(x, n, sizeof(*x), st_phase_step1_double_cmp);
  if (n & 1u)
    return x[n / 2u];
  return 0.5 * (x[n / 2u - 1u] + x[n / 2u]);
}

static double st_phase_step1_quantile7(const double *sorted, size_t n,
                                       double p)
{
  double h, g;
  size_t j;
  if ((NULL == sorted) || (0u == n) || !isfinite(p) || (p < 0.0) ||
      (p > 1.0))
    return NAN;
  if (1u == n)
    return sorted[0];
  h = (double)(n - 1u) * p;
  j = (size_t)floor(h);
  g = h - (double)j;
  if (j >= n - 1u)
    return sorted[n - 1u];
  return sorted[j] + g * (sorted[j + 1u] - sorted[j]);
}

static size_t st_phase_step1_sample_index(size_t k, size_t n_sample,
                                          size_t n_ind)
{
  double x;
  if ((n_sample <= 1u) || (n_ind <= 1u))
    return 0u;
  x = (double)k * (double)(n_ind - 1u) / (double)(n_sample - 1u);
  return (size_t)floor(x + 0.5);
}

static int st_phase_step1_marker_common(const st_phase_data *d,
                                        size_t marker, size_t n_sample)
{
  size_t n_observed = 0u, n_alt = 0u, min_observed;
  if ((NULL == d) || (marker >= d->n_marker) ||
      (2u != d->n_alleles[marker]) || (0u == n_sample))
    return 0;
  for (size_t k = 0u; k < n_sample; k++)
  {
    const size_t i = st_phase_step1_sample_index(k, n_sample, d->n_ind);
    const size_t z = marker * d->n_ind + i;
    const uint8_t a = d->geno1[z], b = d->geno2[z];
    if ((ST_PHASE_STATE_MISSING == a) || (ST_PHASE_STATE_MISSING == b))
      continue;
    n_observed++;
    n_alt += (size_t)a + (size_t)b;
  }
  min_observed = (size_t)ceil(0.5 * (double)n_sample);
  if (min_observed < 10u)
    min_observed = 10u;
  if (n_observed < min_observed)
    return 0;
  {
    const double af = (double)n_alt / (2.0 * (double)n_observed);
    return (af >= 0.05) && (af <= 0.95);
  }
}

static double st_phase_step1_marker_pair_r2(const st_phase_data *d,
                                            size_t m1, size_t m2,
                                            size_t n_sample)
{
  size_t n = 0u;
  double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
  if ((NULL == d) || (m1 >= d->n_marker) || (m2 >= d->n_marker))
    return NAN;
  for (size_t k = 0u; k < n_sample; k++)
  {
    const size_t i = st_phase_step1_sample_index(k, n_sample, d->n_ind);
    const size_t z1 = m1 * d->n_ind + i, z2 = m2 * d->n_ind + i;
    const uint8_t a1 = d->geno1[z1], b1 = d->geno2[z1];
    const uint8_t a2 = d->geno1[z2], b2 = d->geno2[z2];
    double x, y;
    if ((ST_PHASE_STATE_MISSING == a1) || (ST_PHASE_STATE_MISSING == b1) ||
        (ST_PHASE_STATE_MISSING == a2) || (ST_PHASE_STATE_MISSING == b2))
      continue;
    x = (double)a1 + (double)b1;
    y = (double)a2 + (double)b2;
    n++;
    sx += x;
    sy += y;
    sxx += x * x;
    syy += y * y;
    sxy += x * y;
  }
  if (n < 2u)
    return NAN;
  {
    const double vx = (double)n * sxx - sx * sx;
    const double vy = (double)n * syy - sy * sy;
    const double cv = (double)n * sxy - sx * sy;
    double r2;
    if (!(vx > 0.0) || !(vy > 0.0))
      return NAN;
    r2 = (cv * cv) / (vx * vy);
    if (!isfinite(r2))
      return NAN;
    if (r2 < 0.0)
      r2 = 0.0;
    if (r2 > 1.0)
      r2 = 1.0;
    return r2;
  }
}

static size_t st_phase_step1_collect_ld_pairs(
    const st_phase_data *d, const unsigned char *common, size_t n_sample,
    double lo_fraction, double hi_fraction, double *distance, double *r2,
    size_t capacity)
{
  size_t n_pair = 0u;
  size_t max_per_chromosome;
  if ((NULL == d) || (NULL == common) || (NULL == distance) ||
      (NULL == r2) || (0u == capacity) || (0u == d->n_chrom))
    return 0u;
  max_per_chromosome = (capacity + d->n_chrom - 1u) / d->n_chrom;
  if (max_per_chromosome < 200u)
    max_per_chromosome = 200u;

  for (size_t c = 0u; (c < d->n_chrom) && (n_pair < capacity); c++)
  {
    const size_t start = d->chrom_start[c], length = d->chrom_length[c];
    const double span =
        (d->pos_cm[start + length - 1u] - d->pos_cm[start]) / 100.0;
    size_t *marker = NULL, n_common = 0u, n_chromosome = 0u, n_order;
    if ((length < 2u) || !isfinite(span) || !(span > 0.0))
      continue;
    marker = malloc(length * sizeof(*marker));
    if (NULL == marker)
      return n_pair;
    for (size_t q = 0u; q < length; q++)
      if (common[start + q])
        marker[n_common++] = start + q;
    if (n_common < 2u)
    {
      free(marker);
      continue;
    }
    n_order = (n_common < 300u) ? n_common : 300u;
    for (size_t kk = 0u; (kk < n_order) && (n_pair < capacity) &&
                         (n_chromosome < max_per_chromosome);
         kk++)
    {
      size_t ii;
      if (n_order <= 1u)
        ii = 0u;
      else
        ii = (size_t)floor((double)kk * (double)(n_common - 1u) /
                               (double)(n_order - 1u) +
                           0.5);
      if (ii >= n_common - 1u)
        continue;
      {
        size_t jj_max = ii + 64u;
        if (jj_max >= n_common)
          jj_max = n_common - 1u;
        for (size_t jj = ii + 1u; (jj <= jj_max) && (n_pair < capacity) &&
                                  (n_chromosome < max_per_chromosome);
             jj++)
        {
          const double dist =
              (d->pos_cm[marker[jj]] - d->pos_cm[marker[ii]]) / 100.0;
          double rr;
          if (!isfinite(dist) || (dist < lo_fraction * span) ||
              (dist > hi_fraction * span))
            continue;
          rr = st_phase_step1_marker_pair_r2(d, marker[ii], marker[jj],
                                             n_sample);
          if (!isfinite(rr))
            continue;
          distance[n_pair] = dist;
          r2[n_pair] = rr;
          n_pair++;
          n_chromosome++;
        }
      }
    }
    free(marker);
  }
  return n_pair;
}

static double st_phase_step1_fininfty_objective(double log_lambda,
                                                const double *distance,
                                                const double *r2, size_t n,
                                                double sampling_r2,
                                                double n_haplotype)
{
  const double lambda = exp(log_lambda);
  double value = 0.0;
  for (size_t i = 0u; i < n; i++)
  {
    const double pred =
        sampling_r2 +
        (1.0 - sampling_r2) / (1.0 + n_haplotype * lambda * distance[i]);
    const double a = log(fmax(r2[i], 1e-5));
    const double b = log(fmax(pred, 1e-5));
    const double d = a - b;
    value += d * d;
  }
  return value;
}

static double st_phase_step1_fininfty_fit(const double *distance,
                                          const double *r2, size_t n,
                                          double sampling_r2,
                                          double n_haplotype)
{
  const double ratio = 0.6180339887498948482;
  double a = log(0.05), b = log(1e6);
  double c = b - ratio * (b - a), d = a + ratio * (b - a);
  double fc = st_phase_step1_fininfty_objective(c, distance, r2, n,
                                                sampling_r2, n_haplotype);
  double fd = st_phase_step1_fininfty_objective(d, distance, r2, n,
                                                sampling_r2, n_haplotype);
  for (int iter = 0; iter < 80; iter++)
  {
    if (fc <= fd)
    {
      b = d;
      d = c;
      fd = fc;
      c = b - ratio * (b - a);
      fc = st_phase_step1_fininfty_objective(c, distance, r2, n,
                                             sampling_r2, n_haplotype);
    }
    else
    {
      a = c;
      c = d;
      fc = fd;
      d = a + ratio * (b - a);
      fd = st_phase_step1_fininfty_objective(d, distance, r2, n,
                                             sampling_r2, n_haplotype);
    }
  }
  return exp(0.5 * (a + b));
}

static double
st_phase_step1_estimate_fininfty_breaks(const st_phase_data *d)
{
  enum
  {
    MAX_PAIR = 20000
  };
  unsigned char *common = NULL;
  double *distance = NULL, *r2 = NULL, *sorted = NULL, *tmp1 = NULL,
         *tmp2 = NULL;
  double bd[10], br2[10], breaks[11], fallback, result;
  size_t n_sample, n_pair = 0u, n_break = 0u, n_bin = 0u;
  if ((NULL == d) || (0u == d->n_ind))
    return NAN;
  fallback = 200000.0 / (double)d->n_ind;
  if ((d->n_ind < 10u) || (d->n_marker < 20u))
    return fallback;
  n_sample = (d->n_ind < 500u) ? d->n_ind : 500u;
  common = calloc(d->n_marker, sizeof(*common));
  distance = malloc(MAX_PAIR * sizeof(*distance));
  r2 = malloc(MAX_PAIR * sizeof(*r2));
  sorted = malloc(MAX_PAIR * sizeof(*sorted));
  tmp1 = malloc(MAX_PAIR * sizeof(*tmp1));
  tmp2 = malloc(MAX_PAIR * sizeof(*tmp2));
  if ((NULL == common) || (NULL == distance) || (NULL == r2) ||
      (NULL == sorted) || (NULL == tmp1) || (NULL == tmp2))
  {
    result = fallback;
    goto cleanup;
  }
  for (size_t m = 0u; m < d->n_marker; m++)
    common[m] =
        (unsigned char)st_phase_step1_marker_common(d, m, n_sample);
  n_pair = st_phase_step1_collect_ld_pairs(d, common, n_sample, 0.005,
                                           0.05, distance, r2, MAX_PAIR);
  if (n_pair < 100u)
    n_pair = st_phase_step1_collect_ld_pairs(d, common, n_sample, 0.002,
                                             0.10, distance, r2, MAX_PAIR);
  if (n_pair < 50u)
  {
    result = fallback;
    goto cleanup;
  }
  memcpy(sorted, distance, n_pair * sizeof(*sorted));
  qsort(sorted, n_pair, sizeof(*sorted), st_phase_step1_double_cmp);
  for (int q = 0; q <= 10; q++)
  {
    const double x =
        st_phase_step1_quantile7(sorted, n_pair, (double)q / 10.0);
    if ((0u == n_break) || (x != breaks[n_break - 1u]))
      breaks[n_break++] = x;
  }
  if (n_break < 4u)
  {
    result = fallback;
    goto cleanup;
  }
  for (size_t q = 0u; q + 1u < n_break; q++)
  {
    size_t n = 0u;
    for (size_t j = 0u; j < n_pair; j++)
      if ((distance[j] >= breaks[q]) && (distance[j] <= breaks[q + 1u]))
      {
        tmp1[n] = distance[j];
        tmp2[n] = r2[j];
        n++;
      }
    if (n >= 10u)
    {
      bd[n_bin] = st_phase_step1_median(tmp1, n);
      br2[n_bin] = st_phase_step1_median(tmp2, n);
      n_bin++;
    }
  }
  if (n_bin < 3u)
  {
    result = fallback;
    goto cleanup;
  }
  result = st_phase_step1_fininfty_fit(
      bd, br2, n_bin, 1.0 / fmax(2.0 * (double)n_sample, 2.0),
      2.0 * (double)d->n_ind);
  if (!isfinite(result) || (result < 0.0))
    result = fallback;
cleanup:
  free(common);
  free(distance);
  free(r2);
  free(sorted);
  free(tmp1);
  free(tmp2);
  return result;
}

static int st_phase_step1_resolve_population(
    SEXP population_type_, int reference_mode,
    size_t n_reference_individual, int heterozygous_reference,
    st_phase_step1_configuration *cfg, const char **message)
{
  if ((NULL == cfg) || (NULL == message))
    return 0;
  *message = NULL;
  if (R_NilValue != population_type_)
  {
    if (!st_phase_step1_population_parse(population_type_,
                                         &cfg->population_type))
    {
      *message = "population.type must be one of DH, F2, F3, F4, S2, "
                 "S3, S4, SSD, or Finfty";
      return 0;
    }
  }
  else if (reference_mode && (2u == n_reference_individual) &&
           !heterozygous_reference)
    cfg->population_type = ST_PHASE_POP_F2;
  else
    cfg->population_type = ST_PHASE_POP_FINFTY;
  return 1;
}

static int st_phase_step1_resolve_breaks(
    SEXP population_type_, SEXP expected_breaks_, SEXP break_rate1,
    SEXP break_rate2, const st_phase_data *d,
    st_phase_step1_configuration *cfg, const char **message)
{
  const int population_supplied = (R_NilValue != population_type_);
  const int expected_supplied = (R_NilValue != expected_breaks_);
  const int b1_supplied = (R_NilValue != break_rate1);
  const int b2_supplied = (R_NilValue != break_rate2);
  double x;
  if ((NULL == cfg) || (NULL == message))
    return 0;
  *message = NULL;

  if (b1_supplied || b2_supplied)
  {
    if (!(b1_supplied && b2_supplied))
    {
      *message = "break.rate1.per.morgan and break.rate2.per.morgan must "
                 "be supplied together";
      return 0;
    }
    if (expected_supplied)
    {
      *message = "expected.breaks.per.morgan cannot be combined with the "
                 "legacy break-rate pair";
      return 0;
    }
    if (population_supplied)
    {
      *message = "population.type cannot be combined with the legacy "
                 "break-rate pair";
      return 0;
    }
    if (!st_phase_r_scalar_double(break_rate1,
                                  &cfg->break_rate1_per_morgan) ||
        !st_phase_r_scalar_double(break_rate2,
                                  &cfg->break_rate2_per_morgan) ||
        (cfg->break_rate1_per_morgan < 0.0) ||
        (cfg->break_rate2_per_morgan < 0.0))
    {
      *message = "break rates must be finite and nonnegative";
      return 0;
    }
    cfg->expected_breaks_per_morgan = NA_REAL;
    cfg->breaks_estimated = 0;
    return 1;
  }

  if (expected_supplied)
  {
    if (!st_phase_r_scalar_double(expected_breaks_, &x) || (x < 0.0))
    {
      *message = "expected.breaks.per.morgan must be one finite "
                 "nonnegative number";
      return 0;
    }
    cfg->expected_breaks_per_morgan = x;
    cfg->breaks_estimated = 0;
  }
  else if (ST_PHASE_POP_FINFTY == cfg->population_type)
  {
    if (NULL == d)
    {
      *message = "Internal error: Finfty estimation requires phasing data";
      return 0;
    }
    x = st_phase_step1_estimate_fininfty_breaks(d);
    if (!isfinite(x) || (x < 0.0))
    {
      *message = "Cannot estimate Finfty copying breaks";
      return 0;
    }
    cfg->expected_breaks_per_morgan = x;
    cfg->breaks_estimated = 1;
  }
  else
  {
    x = st_phase_step1_population_breaks(cfg->population_type);
    if (!isfinite(x) || (x < 0.0))
    {
      *message = "Internal error resolving population-specific breaks";
      return 0;
    }
    cfg->expected_breaks_per_morgan = x;
    cfg->breaks_estimated = 0;
  }
  cfg->break_rate1_per_morgan = cfg->expected_breaks_per_morgan;
  cfg->break_rate2_per_morgan = cfg->expected_breaks_per_morgan;
  return 1;
}

static const char *st_phase_step1_cohort_status_name(int x)
{
  switch (x)
  {
  case ST_PHASE_COHORT_PIPELINE_OK:
    return "ok";
  case ST_PHASE_COHORT_PIPELINE_INVALID:
    return "invalid";
  case ST_PHASE_COHORT_PIPELINE_DATA_FAILED:
    return "data_failed";
  case ST_PHASE_COHORT_PIPELINE_PHASING_FAILED:
    return "phasing_failed";
  case ST_PHASE_COHORT_PIPELINE_STAGING_FAILED:
    return "staging_failed";
  case ST_PHASE_COHORT_PIPELINE_COMMIT_FAILED:
    return "commit_failed";
  default:
    return "unknown";
  }
}

static SEXP st_phase_step1_reference_diagnostics(
    const st_phase_step0_pipeline_diagnostics *d,
    const char *reference_name, size_t n_reference_individual,
    size_t n_reference_haplotype, int heterozygous_phased)
{
  SEXP base = PROTECT(st_phase_step0_r_diagnostics(d));
  SEXP base_names = getAttrib(base, R_NamesSymbol);
  R_xlen_t nb = XLENGTH(base);
  SEXP ans = PROTECT(allocVector(VECSXP, nb + 5));
  SEXP nms = PROTECT(allocVector(STRSXP, nb + 5));
  for (R_xlen_t k = 0; k < nb; k++)
  {
    SET_VECTOR_ELT(ans, k, VECTOR_ELT(base, k));
    SET_STRING_ELT(nms, k, STRING_ELT(base_names, k));
  }
  SET_VECTOR_ELT(ans, nb, mkString("reference"));
  SET_STRING_ELT(nms, nb, mkChar("mode"));
  SET_VECTOR_ELT(ans, nb + 1, mkString(reference_name));
  SET_STRING_ELT(nms, nb + 1, mkChar("reference.data.set"));
  SET_VECTOR_ELT(ans, nb + 2, st_phase_r_count(n_reference_individual));
  SET_STRING_ELT(nms, nb + 2, mkChar("n.reference.individuals"));
  SET_VECTOR_ELT(ans, nb + 3, st_phase_r_count(n_reference_haplotype));
  SET_STRING_ELT(nms, nb + 3, mkChar("n.reference.haplotypes"));
  SET_VECTOR_ELT(ans, nb + 4, ScalarLogical(heterozygous_phased));
  SET_STRING_ELT(nms, nb + 4,
                 mkChar("heterozygous.reference.interpreted.as.phased"));
  setAttrib(ans, R_NamesSymbol, nms);
  UNPROTECT(3);
  return ans;
}

typedef struct
{
  int pipeline_status;
  int phasing_succeeded;
  int staging_succeeded;
  int audit_status;
  int commit_requested;
  int commit_attempted;
  int commit_status;
  int committed;
  int impute_requested;
  double imputation_threshold;
  size_t n_ind, n_marker, n_chrom, n_genotype;
  size_t stage1_markers, stage2_markers, stage1_minor_allele_count;
  size_t burnin_states, local_states, max_states;
  double provisional_threshold;
  size_t burnin_iterations_completed;
  size_t burnin_last_changed_link;
  size_t burnin_relative_link;
  double burnin_last_change_fraction;
  size_t progressive_iterations_completed;
  size_t progressive_last_changed_link;
  size_t progressive_relative_link;
  double progressive_last_change_fraction;
  size_t n_confident_missing_haplotype;
  size_t n_target_chromosome, n_phase_ready_chromosome;
  size_t n_partial_fixed_chromosome, n_unresolved_chromosome;
  size_t n_no_heterozygote_chromosome;
  size_t n_valid_link, n_confident_link, n_unresolved_link;
  size_t n_fixed_same, n_fixed_switch;
  size_t n_observed_heterozygous, n_confident_heterozygous;
  size_t n_unresolved_heterozygous;
  size_t n_covered, n_uncovered, n_written, n_covered_unwritten;
  size_t n_changed_genotypes;
  size_t n_missing, n_homozygous, n_heterozygous;
  size_t n_covered_heterozygous_unresolved;
  size_t n_written_same_source_order, n_written_swapped_source_order;
  int audit_complete, audit_source_bytes_current, audit_workspace_current;
  size_t n_missing_requested, n_imputed, n_imputation_unresolved;
  size_t n_imputed_hom0, n_imputed_hom1, n_imputed_heterozygous;
  size_t n_imputed_heterozygous_phased, n_imputed_heterozygous_unphased;
} st_phase_two_stage_pipeline_diagnostics;

static void st_phase_two_stage_pipeline_diagnostics_init(
    st_phase_two_stage_pipeline_diagnostics *x)
{
  if (NULL != x)
  {
    memset(x, 0, sizeof(*x));
    x->pipeline_status = ST_PHASE_COHORT_PIPELINE_INVALID;
    x->audit_status = ST_PHASE_AUDIT_INVALID;
    x->commit_status = ST_PHASE_COMMIT_INVALID;
  }
}

static int st_phase_dataset_pipeline_cohort_two_stage(
    gs_varset_type *gv, size_t burnin_max, size_t burnin_min,
    double burnin_convergence, size_t progressive_max,
    size_t progressive_min, double progressive_convergence,
    size_t burnin_states, size_t local_states, size_t max_states,
    size_t stage1_minor_allele_count, double provisional_threshold,
    double error_rate, double break_rate1_per_morgan,
    double break_rate2_per_morgan, double confidence_threshold,
    int impute_requested, double imputation_threshold,
    int commit_requested,
    st_phase_two_stage_pipeline_diagnostics *diagnostics);

static SEXP st_phase_step1_cohort_diagnostics(
    const st_phase_two_stage_pipeline_diagnostics *d)
{
  static const char *names[] = {"success",
                                "mode",
                                "pipeline.status",
                                "pipeline.status.name",
                                "phasing.succeeded",
                                "staging.succeeded",
                                "audit.status",
                                "audit.status.name",
                                "commit.requested",
                                "commit.attempted",
                                "commit.status",
                                "commit.status.name",
                                "committed",
                                "n.ind",
                                "n.marker",
                                "n.chrom",
                                "n.genotype",
                                "stage1.markers",
                                "stage2.markers",
                                "stage1.minor.allele.count",
                                "burnin.states",
                                "local.states",
                                "max.states",
                                "provisional.threshold",
                                "burnin.iterations.completed",
                                "burnin.last.changed.relative.link",
                                "burnin.relative.links",
                                "burnin.last.change.fraction",
                                "progressive.iterations.completed",
                                "progressive.last.changed.relative.link",
                                "progressive.relative.links",
                                "progressive.last.change.fraction",
                                "n.confident.missing.haplotypes",
                                "n.target.chromosome",
                                "n.phase.ready.chromosome",
                                "n.partial.fixed.chromosome",
                                "n.unresolved.chromosome",
                                "n.no.heterozygote.chromosome",
                                "n.valid.link",
                                "n.confident.link",
                                "n.unresolved.link",
                                "n.fixed.same",
                                "n.fixed.switch",
                                "n.observed.heterozygous",
                                "n.confident.heterozygous",
                                "n.unresolved.heterozygous",
                                "n.covered",
                                "n.uncovered",
                                "n.written",
                                "n.covered.unwritten",
                                "n.changed.genotypes",
                                "n.missing",
                                "n.homozygous",
                                "n.heterozygous",
                                "n.covered.heterozygous.unresolved",
                                "n.written.same.source.order",
                                "n.written.swapped.source.order",
                                "audit.complete",
                                "audit.source.bytes.current",
                                "audit.workspace.current",
                                "impute",
                                "imputation.threshold",
                                "n.imputation.missing.requested",
                                "n.imputed",
                                "n.imputation.unresolved",
                                "n.imputed.hom0",
                                "n.imputed.hom1",
                                "n.imputed.heterozygous",
                                "n.imputed.heterozygous.phased",
                                "n.imputed.heterozygous.unphased"};
  enum
  {
    N = (int)(sizeof(names) / sizeof(names[0]))
  };
  SEXP ans = PROTECT(allocVector(VECSXP, N));
  SEXP nms = PROTECT(allocVector(STRSXP, N));
  int k = 0;
#define ADD(value)                                                        \
  do                                                                      \
  {                                                                       \
    SET_VECTOR_ELT(ans, k, (value));                                      \
    SET_STRING_ELT(nms, k, mkChar(names[k]));                             \
    k++;                                                                  \
  } while (0)
  ADD(ScalarLogical(ST_PHASE_COHORT_PIPELINE_OK == d->pipeline_status));
  ADD(mkString("cohort"));
  ADD(ScalarInteger(d->pipeline_status));
  ADD(mkString(st_phase_step1_cohort_status_name(d->pipeline_status)));
  ADD(ScalarLogical(d->phasing_succeeded));
  ADD(ScalarLogical(d->staging_succeeded));
  ADD(ScalarInteger(d->audit_status));
  ADD(mkString(st_phase_r_audit_status_name(d->audit_status)));
  ADD(ScalarLogical(d->commit_requested));
  ADD(ScalarLogical(d->commit_attempted));
  ADD(ScalarInteger(d->commit_status));
  ADD(mkString(st_phase_r_commit_status_name(d->commit_status)));
  ADD(ScalarLogical(d->committed));
  ADD(st_phase_r_count(d->n_ind));
  ADD(st_phase_r_count(d->n_marker));
  ADD(st_phase_r_count(d->n_chrom));
  ADD(st_phase_r_count(d->n_genotype));
  ADD(st_phase_r_count(d->stage1_markers));
  ADD(st_phase_r_count(d->stage2_markers));
  ADD(st_phase_r_count(d->stage1_minor_allele_count));
  ADD(st_phase_r_count(d->burnin_states));
  ADD(st_phase_r_count(d->local_states));
  ADD(st_phase_r_count(d->max_states));
  ADD(ScalarReal(d->provisional_threshold));
  ADD(st_phase_r_count(d->burnin_iterations_completed));
  ADD(st_phase_r_count(d->burnin_last_changed_link));
  ADD(st_phase_r_count(d->burnin_relative_link));
  ADD(ScalarReal(d->burnin_last_change_fraction));
  ADD(st_phase_r_count(d->progressive_iterations_completed));
  ADD(st_phase_r_count(d->progressive_last_changed_link));
  ADD(st_phase_r_count(d->progressive_relative_link));
  ADD(ScalarReal(d->progressive_last_change_fraction));
  ADD(st_phase_r_count(d->n_confident_missing_haplotype));
  ADD(st_phase_r_count(d->n_target_chromosome));
  ADD(st_phase_r_count(d->n_phase_ready_chromosome));
  ADD(st_phase_r_count(d->n_partial_fixed_chromosome));
  ADD(st_phase_r_count(d->n_unresolved_chromosome));
  ADD(st_phase_r_count(d->n_no_heterozygote_chromosome));
  ADD(st_phase_r_count(d->n_valid_link));
  ADD(st_phase_r_count(d->n_confident_link));
  ADD(st_phase_r_count(d->n_unresolved_link));
  ADD(st_phase_r_count(d->n_fixed_same));
  ADD(st_phase_r_count(d->n_fixed_switch));
  ADD(st_phase_r_count(d->n_observed_heterozygous));
  ADD(st_phase_r_count(d->n_confident_heterozygous));
  ADD(st_phase_r_count(d->n_unresolved_heterozygous));
  ADD(st_phase_r_count(d->n_covered));
  ADD(st_phase_r_count(d->n_uncovered));
  ADD(st_phase_r_count(d->n_written));
  ADD(st_phase_r_count(d->n_covered_unwritten));
  ADD(st_phase_r_count(d->n_changed_genotypes));
  ADD(st_phase_r_count(d->n_missing));
  ADD(st_phase_r_count(d->n_homozygous));
  ADD(st_phase_r_count(d->n_heterozygous));
  ADD(st_phase_r_count(d->n_covered_heterozygous_unresolved));
  ADD(st_phase_r_count(d->n_written_same_source_order));
  ADD(st_phase_r_count(d->n_written_swapped_source_order));
  ADD(ScalarLogical(d->audit_complete));
  ADD(ScalarLogical(d->audit_source_bytes_current));
  ADD(ScalarLogical(d->audit_workspace_current));
  ADD(ScalarLogical(d->impute_requested));
  ADD(ScalarReal(d->imputation_threshold));
  ADD(st_phase_r_count(d->n_missing_requested));
  ADD(st_phase_r_count(d->n_imputed));
  ADD(st_phase_r_count(d->n_imputation_unresolved));
  ADD(st_phase_r_count(d->n_imputed_hom0));
  ADD(st_phase_r_count(d->n_imputed_hom1));
  ADD(st_phase_r_count(d->n_imputed_heterozygous));
  ADD(st_phase_r_count(d->n_imputed_heterozygous_phased));
  ADD(st_phase_r_count(d->n_imputed_heterozygous_unphased));
#undef ADD
  if (k != N)
    error("Internal st.phase cohort diagnostics length mismatch");
  setAttrib(ans, R_NamesSymbol, nms);
  UNPROTECT(2);
  return ans;
}

static SEXP st_phase_step1_append_configuration(
    SEXP base, const st_phase_step1_configuration *cfg,
    int include_max_states)
{
  SEXP base_names, ans, nms, cls;
  R_xlen_t nb, nadd = include_max_states ? 7 : 6, k;
  if (NULL == cfg)
    error("Internal st.phase configuration error");
  PROTECT(base);
  base_names = getAttrib(base, R_NamesSymbol);
  nb = XLENGTH(base);
  PROTECT(ans = allocVector(VECSXP, nb + nadd));
  PROTECT(nms = allocVector(STRSXP, nb + nadd));
  for (k = 0; k < nb; k++)
  {
    SET_VECTOR_ELT(ans, k, VECTOR_ELT(base, k));
    SET_STRING_ELT(nms, k, STRING_ELT(base_names, k));
  }
#define ST_PHASE_CFG_ADD(name, value)                                     \
  do                                                                      \
  {                                                                       \
    SET_VECTOR_ELT(ans, k, (value));                                      \
    SET_STRING_ELT(nms, k, mkChar(name));                                 \
    k++;                                                                  \
  } while (0)
  ST_PHASE_CFG_ADD(
      "population.type",
      mkString(st_phase_step1_population_name(cfg->population_type)));
  ST_PHASE_CFG_ADD("expected.breaks.per.morgan",
                   ScalarReal(cfg->expected_breaks_per_morgan));
  ST_PHASE_CFG_ADD("breaks.estimated",
                   ScalarLogical(cfg->breaks_estimated));
  ST_PHASE_CFG_ADD("break.rate1.per.morgan",
                   ScalarReal(cfg->break_rate1_per_morgan));
  ST_PHASE_CFG_ADD("break.rate2.per.morgan",
                   ScalarReal(cfg->break_rate2_per_morgan));
  ST_PHASE_CFG_ADD("error.rate", ScalarReal(cfg->error_rate));
  if (include_max_states)
    ST_PHASE_CFG_ADD("max.states", st_phase_r_count(cfg->max_states));
#undef ST_PHASE_CFG_ADD
  if (k != nb + nadd)
    error("Internal st.phase configuration length mismatch");
  setAttrib(ans, R_NamesSymbol, nms);
  PROTECT(cls = allocVector(STRSXP, 2));
  SET_STRING_ELT(cls, 0, mkChar("st.phase.result"));
  SET_STRING_ELT(cls, 1, mkChar("list"));
  setAttrib(ans, R_ClassSymbol, cls);
  UNPROTECT(4);
  return ans;
}

SEXP st_phase_03_step1(SEXP data_set, SEXP reference_data_set,
                       SEXP reference_individuals, SEXP prior,
                       SEXP population_type_, SEXP expected_breaks_,
                       SEXP max_states, SEXP individuals, SEXP chromosomes,
                       SEXP error_rate, SEXP break_rate1, SEXP break_rate2,
                       SEXP confidence_threshold, SEXP impute,
                       SEXP imputation_threshold, SEXP burnin_max,
                       SEXP burnin_min, SEXP burnin_convergence,
                       SEXP progressive_max, SEXP progressive_min,
                       SEXP progressive_convergence, SEXP burnin_states_,
                       SEXP local_states_, SEXP stage1_minor_allele_count_,
                       SEXP provisional_threshold_, SEXP commit)
{
  const char *set_name = NULL, *reference_name = NULL, *cfg_error = NULL;
  gs_varset_type *target_gv = NULL;
  st_phase_step1_configuration cfg;
  double threshold = 0.0, imp_threshold = 0.0;
  int do_impute = 0, do_commit = 0;

  memset(&cfg, 0, sizeof(cfg));
  if (!st_phase_step1_dataset_name(data_set, &set_name))
    error("data.set must name an existing SelectionTools data set");
  target_gv = st_phase_step1_find_dataset(set_name);
  if (NULL == target_gv)
    error("SelectionTools data set not found: %s", set_name);
  if (!st_phase_step1_double_default(error_rate, 0.001, &cfg.error_rate) ||
      !st_phase_hmm_error_valid(cfg.error_rate))
    error("error.rate must be between 0 and 1");
  if (!st_phase_step1_double_default(confidence_threshold, 0.50,
                                     &threshold) ||
      (threshold < 0.5) || (threshold > 1.0))
    error("confidence.threshold must be between 0.5 and 1");
  if (!st_phase_r_scalar_logical(impute, &do_impute))
    error("impute must be TRUE or FALSE");
  if (!st_phase_step1_double_default(imputation_threshold, 0.90,
                                     &imp_threshold) ||
      !st_phase_step0_imputation_threshold_valid(imp_threshold))
    error("imputation.threshold must be between 0 and 1");
  if (!st_phase_r_scalar_logical(commit, &do_commit))
    error("commit must be TRUE or FALSE");

  if (R_NilValue != reference_data_set)
  {
    gs_varset_type *reference_gv = NULL;
    st_phase_data d;
    st_phase_step0_pipeline_diagnostics diag;
    uint8_t *ref = NULL;
    double *w = NULL;
    size_t *ind = NULL, *chrom = NULL;
    st_phase_dataset_request *request = NULL;
    size_t R = 0u, n_reference_individual = 0u, n_ind = 0u, n_chrom = 0u,
           n_request = 0u;
    size_t default_states;
    int heterozygous_reference = 0;
    SEXP base;

    if (!st_phase_step1_dataset_name(reference_data_set, &reference_name))
      error("reference.data.set must be NULL or name an existing "
            "SelectionTools data set");
    reference_gv = st_phase_step1_find_dataset(reference_name);
    if (NULL == reference_gv)
      error("Reference SelectionTools data set not found: %s",
            reference_name);
    if ((R_NilValue != burnin_max) || (R_NilValue != burnin_min) ||
        (R_NilValue != burnin_convergence) ||
        (R_NilValue != progressive_max) ||
        (R_NilValue != progressive_min) ||
        (R_NilValue != progressive_convergence) ||
        (R_NilValue != burnin_states_) || (R_NilValue != local_states_) ||
        (R_NilValue != stage1_minor_allele_count_) ||
        (R_NilValue != provisional_threshold_))
      error("Reference-data-set mode does not accept cohort-specific "
            "arguments");

    st_phase_data_init(&d);
    st_phase_step0_pipeline_diagnostics_init(&diag);
    if (0 != st_phase_build_data(target_gv, &d))
      error("Cannot build phasing workspace for data.set '%s'", set_name);
    {
      int reference_status = st_phase_step1_build_reference_dataset(
          target_gv, &d, reference_gv, reference_individuals, &ref,
          &R, &n_reference_individual, &heterozygous_reference);
      if (0 != reference_status)
      {
        st_phase_data_free(&d);
        if (-2 == reference_status)
          error("reference.data.set contains heterozygous genotypes but "
                "is not marked phased");
        error("reference.data.set/reference.individuals must define "
              "completely observed reference individuals with the same "
              "marker set and compatible allele codes");
      }
    }
    if (!st_phase_step1_resolve_population(
            population_type_, 1, n_reference_individual,
            heterozygous_reference, &cfg, &cfg_error) ||
        !st_phase_step1_resolve_breaks(population_type_, expected_breaks_,
                                       break_rate1, break_rate2, &d, &cfg,
                                       &cfg_error))
    {
      free(ref);
      st_phase_data_free(&d);
      error("%s", (NULL != cfg_error)
                      ? cfg_error
                      : "Invalid st.phase population/break configuration");
    }
    default_states = (R < 60u) ? R : 60u;
    if (!st_phase_step1_size_default(max_states, default_states,
                                     &cfg.max_states))
    {
      free(ref);
      st_phase_data_free(&d);
      error("max.states must be one positive whole number");
    }
    if (R_NilValue == prior)
    {
      if (0 != st_phase_step1_default_prior(R, &w))
      {
        free(ref);
        st_phase_data_free(&d);
        error("Cannot construct default initial.pair.weights");
      }
    }
    else if (0 != st_phase_r_build_prior(prior, R, &w))
    {
      free(ref);
      st_phase_data_free(&d);
      error("initial.pair.weights must be NULL or a nonnegative R x R "
            "matrix with positive total weight, where R is the number of "
            "reference haplotypes");
    }
    if (0 != st_phase_r_individuals(target_gv, individuals, &ind, &n_ind))
    {
      free(w);
      free(ref);
      st_phase_data_free(&d);
      error("individuals are invalid, duplicated, or not found");
    }
    if (0 != st_phase_r_chromosomes(&d, chromosomes, &chrom, &n_chrom))
    {
      free(ind);
      free(w);
      free(ref);
      st_phase_data_free(&d);
      error("chromosomes are invalid, duplicated, or not found");
    }
    if ((0 != st_phase_size_mul(n_ind, n_chrom, &n_request)) ||
        (0u == n_request) || (n_request > SIZE_MAX / sizeof(*request)))
    {
      free(chrom);
      free(ind);
      free(w);
      free(ref);
      st_phase_data_free(&d);
      error("Too many phasing requests");
    }
    request = malloc(n_request * sizeof(*request));
    if (NULL == request)
    {
      free(chrom);
      free(ind);
      free(w);
      free(ref);
      st_phase_data_free(&d);
      error("Memory allocation error");
    }
    {
      size_t z = 0u;
      for (size_t i = 0u; i < n_ind; i++)
        for (size_t c = 0u; c < n_chrom; c++)
        {
          request[z].individual = ind[i];
          request[z].chromosome = chrom[c];
          z++;
        }
    }
    {
      int phase_status =
        st_phase_dataset_pipeline_binary_reference_bounded_factorized_step0(
            target_gv, request, n_request, ref, d.n_marker, R,
            cfg.max_states, w, R * R, cfg.error_rate,
            cfg.break_rate1_per_morgan, cfg.break_rate2_per_morgan,
            threshold, do_impute, imp_threshold, do_commit, &diag);
      if ((ST_PHASE_PIPELINE_OK == phase_status) && do_commit &&
          diag.base.committed && (n_ind == (size_t)target_gv->NoInd) &&
          (n_chrom == d.n_chrom) &&
          (0u == diag.base.n_covered_heterozygous_unresolved))
        target_gv->phased = 1;
    }
    free(request);
    free(chrom);
    free(ind);
    free(w);
    free(ref);
    st_phase_data_free(&d);
    base = st_phase_step1_reference_diagnostics(&diag, reference_name,
                                                n_reference_individual, R,
                                                heterozygous_reference);
    return st_phase_step1_append_configuration(base, &cfg, 1);
  }
  else
  {
    st_phase_two_stage_pipeline_diagnostics diag;
    st_phase_data estimate_data;
    size_t bmax = 0u, bmin = 0u, pmax = 0u, pmin = 0u;
    size_t burn_states = 0u, local_states = 0u, mac_cutoff = 0u;
    double bconv = 0.0, pconv = 0.0, provisional = 0.0;
    int need_estimate_data = 0;
    SEXP base;

    if ((R_NilValue != reference_individuals) || (R_NilValue != prior) ||
        (R_NilValue != individuals) || (R_NilValue != chromosomes))
      error("reference.individuals, initial.pair.weights, individuals and "
            "chromosomes must be NULL in reference-free cohort mode");
    if (!st_phase_step1_size_default(max_states, 60u, &cfg.max_states))
      error("max.states must be one positive whole number");
    if (!st_phase_step1_size_default(burnin_max, 4u, &bmax) ||
        !st_phase_step1_size_default(burnin_min, 2u, &bmin) ||
        !st_phase_step1_double_default(burnin_convergence, 0.001,
                                       &bconv) ||
        !st_phase_step1_size_default(progressive_max, 2u, &pmax) ||
        !st_phase_step1_size_default(progressive_min, 1u, &pmin) ||
        !st_phase_step1_double_default(progressive_convergence, 0.001,
                                       &pconv) ||
        !st_phase_step1_size_default(burnin_states_, 12u, &burn_states) ||
        !st_phase_step1_size_default(local_states_, 8u, &local_states) ||
        !st_phase_step1_size_default(stage1_minor_allele_count_, 4u,
                                     &mac_cutoff) ||
        !st_phase_step1_double_default(provisional_threshold_, 0.80,
                                       &provisional) ||
        (bmin > bmax) || (pmin > pmax) || (bconv < 0.0) || (bconv > 1.0) ||
        (pconv < 0.0) || (pconv > 1.0) || (provisional < 0.5) ||
        (provisional > 1.0))
      error("Reference-free cohort mode has invalid burn-in, progressive "
            "or two-stage arguments");
    if (!st_phase_step1_resolve_population(population_type_, 0, 0u, 0,
                                           &cfg, &cfg_error))
      error("%s", (NULL != cfg_error)
                      ? cfg_error
                      : "Invalid st.phase population.type");

    need_estimate_data = (R_NilValue == expected_breaks_) &&
                         (R_NilValue == break_rate1) &&
                         (R_NilValue == break_rate2) &&
                         (ST_PHASE_POP_FINFTY == cfg.population_type);
    st_phase_data_init(&estimate_data);
    if (need_estimate_data)
    {
      if (0 != st_phase_build_data(target_gv, &estimate_data))
        error("Cannot build phasing workspace for Finfty estimation");
    }
    if (!st_phase_step1_resolve_breaks(
            population_type_, expected_breaks_, break_rate1, break_rate2,
            need_estimate_data ? &estimate_data : NULL, &cfg, &cfg_error))
    {
      st_phase_data_free(&estimate_data);
      error("%s", (NULL != cfg_error)
                      ? cfg_error
                      : "Invalid st.phase population/break configuration");
    }
    st_phase_data_free(&estimate_data);

    st_phase_two_stage_pipeline_diagnostics_init(&diag);
    {
      int phase_status = st_phase_dataset_pipeline_cohort_two_stage(
        target_gv, bmax, bmin, bconv, pmax, pmin, pconv, burn_states,
        local_states, cfg.max_states, mac_cutoff, provisional,
        cfg.error_rate, cfg.break_rate1_per_morgan,
        cfg.break_rate2_per_morgan, threshold, do_impute, imp_threshold,
        do_commit, &diag);
      if ((ST_PHASE_COHORT_PIPELINE_OK == phase_status) && do_commit &&
          diag.committed && (0u == diag.n_unresolved_heterozygous) &&
          (0u == diag.n_covered_heterozygous_unresolved))
        target_gv->phased = 1;
    }
    base = st_phase_step1_cohort_diagnostics(&diag);
    return st_phase_step1_append_configuration(base, &cfg, 0);
  }
}

static void st_phase_factor_normalize_max(double *x, size_t n)
{
  double best = -INFINITY;
  if (NULL == x)
    return;
  for (size_t j = 0u; j < n; j++)
    if (x[j] > best)
      best = x[j];
  if (best > -INFINITY)
    for (size_t j = 0u; j < n; j++)
      if (x[j] > -INFINITY)
        x[j] -= best;
}

typedef struct
{
  int possible;
  size_t n_marker;
  size_t n_link;
  uint8_t *relation;
  double *margin;
} st_phase_cohort_maxlink_result;

typedef struct
{
  int possible;
  size_t n_marker, n_link;
  uint8_t *relation;
  double *log_ratio;
} st_phase_cohort_haplink_result;
static void
st_phase_cohort_haplink_result_init(st_phase_cohort_haplink_result *r)
{
  if (r)
    memset(r, 0, sizeof(*r));
}
static void
st_phase_cohort_haplink_result_free(st_phase_cohort_haplink_result *r)
{
  if (!r)
    return;
  free(r->relation);
  free(r->log_ratio);
  st_phase_cohort_haplink_result_init(r);
}

static double st_phase_logsum_array(const double *x, size_t n)
{
  double s = -INFINITY;
  for (size_t j = 0u; j < n; j++)
    s = st_phase_hmm_logadd(s, x[j]);
  return s;
}

static int st_phase_haploid_full_transition(const double *prev, size_t R,
                                            double q, double *cur)
{
  double ls, stay, sw;
  if (!prev || !cur || 0u == R || !isfinite(q) || q < 0.0 || q > 1.0)
    return -1;
  ls = st_phase_logsum_array(prev, R);
  stay = st_phase_hmm_log_probability(1.0 - q);
  sw = st_phase_hmm_log_probability(q / (double)R);
  for (size_t r = 0u; r < R; r++)
  {
    double a = (prev[r] > -INFINITY && stay > -INFINITY) ? prev[r] + stay
                                                         : -INFINITY;
    double b = (ls > -INFINITY && sw > -INFINITY) ? ls + sw : -INFINITY;
    cur[r] = st_phase_hmm_logadd(a, b);
  }
  return 0;
}

static int st_phase_haploid_endpoint_emission(uint8_t allele,
                                              const uint8_t *refrow,
                                              size_t R, uint8_t A,
                                              double e, double *loge)
{
  if (!refrow || !loge || ST_PHASE_STATE_MISSING == A ||
      !st_phase_hmm_error_valid(e))
    return -1;
  if (0u == A)
  {
    if (ST_PHASE_STATE_MISSING != allele)
      return -1;
    for (size_t r = 0u; r < R; r++)
      loge[r] = 0.0;
    return 0;
  }
  if (ST_PHASE_STATE_MISSING != allele && allele >= A)
    return -1;
  for (size_t r = 0u; r < R; r++)
  {
    double p;
    uint8_t donor = refrow[r];
    if (ST_PHASE_STATE_MISSING == allele)
      p = 1.0;
    else if (ST_PHASE_STATE_MISSING == donor)
      p = 1.0 / (double)A;
    else if (0 != st_phase_hmm_hap_emission(allele, donor, A, e, &p))
      return -1;
    loge[r] = st_phase_hmm_log_probability(p);
  }
  return 0;
}

static int st_phase_haploid_base_emission(uint8_t a, uint8_t b,
                                          const uint8_t *refrow, size_t R,
                                          uint8_t A, double e,
                                          double *loge)
{
  uint8_t obs = ((ST_PHASE_STATE_MISSING != a) && (a == b))
                    ? a
                    : ST_PHASE_STATE_MISSING;
  return st_phase_haploid_endpoint_emission(obs, refrow, R, A, e, loge);
}

static double st_phase_haploid_link_pattern_score(
    size_t L, size_t RR, uint8_t alleleL, uint8_t alleleR,
    const st_phase_diploid_input *input, const uint8_t *n_alleles,
    const uint8_t *ref, size_t R, double e, const double *q,
    const double *forward, const double *backward, double *v, double *w,
    double *loge)
{
  size_t M = input->n_marker;
  (void)M;
  if (0 != st_phase_haploid_endpoint_emission(alleleL, ref + L * R, R,
                                              n_alleles[L], e, loge))
    return -INFINITY;
  for (size_t r = 0u; r < R; r++)
    v[r] = (forward[L * R + r] > -INFINITY && loge[r] > -INFINITY)
               ? forward[L * R + r] + loge[r]
               : -INFINITY;
  for (size_t m = L + 1u; m <= RR; m++)
  {
    if (0 != st_phase_haploid_full_transition(v, R, q[m - 1u], w))
      return -INFINITY;
    if (m < RR)
    {
      if (0 != st_phase_haploid_base_emission(
                   input->observed1[m], input->observed2[m], ref + m * R,
                   R, n_alleles[m], e, loge))
        return -INFINITY;
    }
    else
    {
      if (0 != st_phase_haploid_endpoint_emission(alleleR, ref + m * R, R,
                                                  n_alleles[m], e, loge))
        return -INFINITY;
    }
    for (size_t r = 0u; r < R; r++)
      w[r] = (w[r] > -INFINITY && loge[r] > -INFINITY) ? w[r] + loge[r]
                                                       : -INFINITY;
    {
      double *z = v;
      v = w;
      w = z;
    }
  }
  {
    double s = -INFINITY;
    for (size_t r = 0u; r < R; r++)
      if (v[r] > -INFINITY && backward[RR * R + r] > -INFINITY)
        s = st_phase_hmm_logadd(s, v[r] + backward[RR * R + r]);
    return s;
  }
}

static int st_phase_cohort_haploid_links_chromosome(
    const st_phase_data *d, size_t individual, size_t chromosome,
    const uint8_t *ref, size_t R, double e, double br,
    st_phase_cohort_haplink_result *out)
{
  st_phase_cohort_haplink_result tmp;
  st_phase_diploid_input input;
  double *q = NULL, *f = NULL, *b = NULL, *v = NULL, *w = NULL,
         *loge = NULL, *x = NULL;
  size_t M;
  int status = -1;
  st_phase_cohort_haplink_result_init(&tmp);
  st_phase_diploid_input_init(&input);
  if (!out || !ref || 0u == R || !st_phase_data_adapter_valid(d) ||
      !st_phase_hmm_error_valid(e) || !isfinite(br) || br < 0.0)
    goto cleanup;
  if (0 !=
      st_phase_prepare_diploid_input(d, individual, chromosome, &input))
    goto cleanup;
  M = input.n_marker;
  if (0u == M)
    goto cleanup;
  if (M > SIZE_MAX / R || M * R > SIZE_MAX / sizeof(double))
    goto cleanup;
  f = malloc(M * R * sizeof(*f));
  b = malloc(M * R * sizeof(*b));
  v = malloc(R * sizeof(*v));
  w = malloc(R * sizeof(*w));
  loge = malloc(R * sizeof(*loge));
  x = malloc(R * sizeof(*x));
  tmp.relation = malloc(M * sizeof(*tmp.relation));
  tmp.log_ratio = malloc(M * sizeof(*tmp.log_ratio));
  if (!f || !b || !v || !w || !loge || !x || !tmp.relation ||
      !tmp.log_ratio)
    goto cleanup;
  for (size_t m = 0u; m < M; m++)
  {
    tmp.relation[m] = 255u;
    tmp.log_ratio[m] = NAN;
    for (size_t r = 0u; r < R; r++)
    {
      uint8_t x = ref[m * R + r];
      if (ST_PHASE_STATE_MISSING != x &&
          x >= d->n_alleles[input.marker_start + m])
        goto cleanup;
    }
  }
  if (M > 1u)
  {
    q = malloc((M - 1u) * sizeof(*q));
    if (!q)
      goto cleanup;
    for (size_t m = 1u; m < M; m++)
    {
      size_t g = input.marker_start + m;
      double dist = d->pos_cm[g] - d->pos_cm[g - 1u];
      if (0 != st_phase_hmm_break_probability(dist, br, &q[m - 1u]))
        goto cleanup;
    }
  }

  if (0 != st_phase_haploid_base_emission(
               input.observed1[0], input.observed2[0], ref, R,
               d->n_alleles[input.marker_start], e, loge))
    goto cleanup;
  for (size_t r = 0u; r < R; r++)
    f[r] = (-log((double)R)) + loge[r];
  st_phase_factor_normalize_max(f, R);
  for (size_t m = 1u; m < M; m++)
  {
    if (0 != st_phase_haploid_full_transition(f + (m - 1u) * R, R,
                                              q[m - 1u], v))
      goto cleanup;
    if (0 != st_phase_haploid_base_emission(
                 input.observed1[m], input.observed2[m], ref + m * R, R,
                 d->n_alleles[input.marker_start + m], e, loge))
      goto cleanup;
    for (size_t r = 0u; r < R; r++)
      f[m * R + r] = (v[r] > -INFINITY && loge[r] > -INFINITY)
                         ? v[r] + loge[r]
                         : -INFINITY;
    st_phase_factor_normalize_max(f + m * R, R);
  }

  for (size_t r = 0u; r < R; r++)
    b[(M - 1u) * R + r] = 0.0;
  for (size_t m = M - 1u; m > 0u; m--)
  {
    if (0 != st_phase_haploid_base_emission(
                 input.observed1[m], input.observed2[m], ref + m * R, R,
                 d->n_alleles[input.marker_start + m], e, loge))
      goto cleanup;
    for (size_t r = 0u; r < R; r++)
      x[r] = (b[m * R + r] > -INFINITY && loge[r] > -INFINITY)
                 ? b[m * R + r] + loge[r]
                 : -INFINITY;

    if (0 != st_phase_haploid_full_transition(x, R, q[m - 1u],
                                              b + (m - 1u) * R))
      goto cleanup;
    st_phase_factor_normalize_max(b + (m - 1u) * R, R);
  }
  {
    size_t L = SIZE_MAX;
    for (size_t RR = 0u; RR < M; RR++)
    {
      uint8_t aR = input.observed1[RR], bR = input.observed2[RR];
      if (ST_PHASE_STATE_MISSING == aR || ST_PHASE_STATE_MISSING == bR ||
          aR == bR)
        continue;
      if (SIZE_MAX == L)
      {
        L = RR;
        continue;
      }
      uint8_t aL = input.observed1[L], bL = input.observed2[L];
      double aa = st_phase_haploid_link_pattern_score(
          L, RR, aL, aR, &input, d->n_alleles + input.marker_start, ref, R,
          e, q, f, b, v, w, loge);
      double bb = st_phase_haploid_link_pattern_score(
          L, RR, bL, bR, &input, d->n_alleles + input.marker_start, ref, R,
          e, q, f, b, v, w, loge);
      double ab = st_phase_haploid_link_pattern_score(
          L, RR, aL, bR, &input, d->n_alleles + input.marker_start, ref, R,
          e, q, f, b, v, w, loge);
      double ba = st_phase_haploid_link_pattern_score(
          L, RR, bL, aR, &input, d->n_alleles + input.marker_start, ref, R,
          e, q, f, b, v, w, loge);
      double same = aa + bb, sw = ab + ba;
      if (!isfinite(same) && !isfinite(sw))
        goto cleanup;
      tmp.relation[RR] = (same >= sw) ? 1u : 0u;
      tmp.log_ratio[RR] = fabs(same - sw);
      tmp.n_link++;
      L = RR;
    }
  }
  tmp.possible = 1;
  tmp.n_marker = M;
  st_phase_cohort_haplink_result_free(out);
  *out = tmp;
  st_phase_cohort_haplink_result_init(&tmp);
  status = 0;
cleanup:
  free(q);
  free(f);
  free(b);
  free(v);
  free(w);
  free(loge);
  free(x);
  st_phase_diploid_input_free(&input);
  st_phase_cohort_haplink_result_free(&tmp);
  return status;
}

static int st_phase_cohort_build_loo_reference_masked(
    const st_phase_cohort_snapshot *frozen, size_t target,
    size_t chromosome, uint8_t **proxy_out, uint8_t **masked_out,
    size_t *R_out)
{
  st_phase_cohort_reference_descriptor *desc = NULL;
  uint8_t *proxy = NULL, *masked = NULL;
  size_t start, M, H, R, ncell;
  int status = -1;
  if (!frozen || !proxy_out || !masked_out || !R_out ||
      target >= frozen->n_ind || chromosome >= frozen->n_chrom ||
      frozen->n_ind < 2u || !frozen->hap || !frozen->observed_missing)
    return -1;
  start = frozen->chrom[chromosome].start;
  M = frozen->chrom[chromosome].length;
  H = frozen->n_hap;
  if (0u == M || H != 2u * frozen->n_ind ||
      (frozen->n_ind - 1u) > SIZE_MAX / 2u)
    return -1;
  R = 2u * (frozen->n_ind - 1u);
  if (R > SIZE_MAX / sizeof(*desc) ||
      0 != st_phase_size_mul(M, R, &ncell) ||
      ncell > SIZE_MAX / sizeof(*proxy))
    return -1;
  desc = malloc(R * sizeof(*desc));
  proxy = malloc(ncell);
  masked = malloc(ncell);
  if (!desc || !proxy || !masked)
    goto cleanup;
  {
    size_t r = 0u;
    for (size_t i = 0u; i < frozen->n_ind; i++)
      if (i != target)
        for (size_t side = 0u; side < 2u; side++)
        {
          desc[r].snapshot = frozen;
          desc[r].marker_start = start;
          desc[r].n_marker = M;
          desc[r].n_hap = H;
          desc[r].hap_index = 2u * i + side;
          r++;
        }
    if (r != R)
      goto cleanup;
  }
  qsort(desc, R, sizeof(*desc), st_phase_cohort_reference_descriptor_cmp);
  for (size_t m = 0u; m < M; m++)
    for (size_t r = 0u; r < R; r++)
    {
      size_t hi = desc[r].hap_index, di = hi / 2u;
      uint8_t x = frozen->hap[(start + m) * H + hi];
      if (x > 1u)
        goto cleanup;
      proxy[m * R + r] = x;
      masked[m * R + r] =
          frozen->observed_missing[(start + m) * frozen->n_ind + di]
              ? ST_PHASE_STATE_MISSING
              : x;
    }
  *proxy_out = proxy;
  *masked_out = masked;
  *R_out = R;
  proxy = NULL;
  masked = NULL;
  status = 0;
cleanup:
  free(desc);
  free(proxy);
  free(masked);
  return status;
}

static int st_phase_pair_emission_row_neutral(
    const st_phase_pair_state_set *pair, size_t marker,
    const uint8_t *donor_allele, size_t n_hap, uint8_t observed1,
    uint8_t observed2, uint8_t A, double e, double *phase12,
    double *phase21, double *total)
{
  size_t n;
  const st_phase_pair_state *row;
  int m1 = (ST_PHASE_STATE_MISSING == observed1),
      m2 = (ST_PHASE_STATE_MISSING == observed2);
  if (!pair || !donor_allele || !phase12 || !phase21 || !total ||
      n_hap != pair->n_hap || !st_phase_hmm_error_valid(e) ||
      ST_PHASE_STATE_MISSING == A || m1 != m2 ||
      (!m1 && (0u == A || observed1 >= A || observed2 >= A)) ||
      !st_phase_pair_state_row_valid(pair, marker))
    return -1;
  row = st_phase_pair_state_get(pair, marker, &n, NULL, NULL);
  if (!row)
    return -1;
  if (0u == A)
  {
    if (!m1)
      return -1;
    for (size_t z = 0u; z < n; z++)
    {
      phase12[z] = 0.5;
      phase21[z] = 0.5;
      total[z] = 1.0;
    }
    return 0;
  }
  for (size_t z = 0; z < n; z++)
  {
    uint8_t d1 = donor_allele[row[z].hap1], d2 = donor_allele[row[z].hap2];
    double a1, a2, b1, b2, p12, p21, p;
    if (ST_PHASE_STATE_MISSING != d1 && d1 >= A)
      return -1;
    if (ST_PHASE_STATE_MISSING != d2 && d2 >= A)
      return -1;
    if (m1)
    {
      phase12[z] = 0.5;
      phase21[z] = 0.5;
      total[z] = 1.0;
      continue;
    }
#define HEMIT(obs, don, ret)                                              \
  do                                                                      \
  {                                                                       \
    if (ST_PHASE_STATE_MISSING == (don))                                  \
      (ret) = 1.0 / (double)A;                                            \
    else if (0 != st_phase_hmm_hap_emission((obs), (don), A, e, &(ret)))  \
      return -1;                                                          \
  } while (0)
    HEMIT(observed1, d1, a1);
    HEMIT(observed2, d2, a2);
    if (observed1 == observed2)
    {
      p12 = a1 * a2;
      p21 = 0.0;
      p = p12;
    }
    else
    {
      HEMIT(observed2, d1, b1);
      HEMIT(observed1, d2, b2);
      p12 = a1 * a2;
      p21 = b1 * b2;
      p = p12 + p21;
    }
#undef HEMIT
    phase12[z] = p12;
    phase21[z] = p21;
    total[z] = p;
  }
  return 0;
}

static int st_phase_cohort_viterbi_chromosome_neutral(
    const st_phase_data *d, size_t individual, size_t chromosome,
    const uint8_t *rank_ref, const uint8_t *emit_ref, size_t R, size_t K,
    double e, double br1, double br2,
    st_phase_cohort_viterbi_chromosome_result *out)
{
  st_phase_cohort_viterbi_chromosome_result tmp;
  st_phase_diploid_input input;
  st_phase_unphased_candidate_set cand;
  st_phase_pair_state_set pair;
  double *p12 = NULL, *p21 = NULL, *em = NULL, *q1 = NULL, *q2 = NULL,
         *initial = NULL;
  size_t *h1 = NULL, *h2 = NULL;
  size_t M, P;
  int possible = 0, status = -1;
  st_phase_cohort_viterbi_chromosome_result_init(&tmp);
  st_phase_diploid_input_init(&input);
  st_phase_unphased_candidate_set_init(&cand);
  st_phase_pair_state_set_init(&pair);
  if (!d || !rank_ref || !emit_ref || 0u == R || 0u == K ||
      !st_phase_hmm_error_valid(e) || !isfinite(br1) || !isfinite(br2) ||
      br1 < 0.0 || br2 < 0.0)
    goto cleanup;
  if (0 !=
      st_phase_prepare_diploid_input(d, individual, chromosome, &input))
    goto cleanup;
  M = input.n_marker;
  if (0u == M)
    goto cleanup;
  if (0 != st_phase_unphased_candidate_build_cohort_balanced(
               input.observed1, input.observed2,
               d->n_alleles + input.marker_start, rank_ref, M, R, K,
               &cand) ||
      0 != st_phase_pair_state_build_unphased(&cand, &pair))
    goto cleanup;
  P = pair.max_pair_states;
  if (0u == P || pair.n_slot != M * P)
    goto cleanup;
  p12 = malloc(pair.n_slot * sizeof(*p12));
  p21 = malloc(pair.n_slot * sizeof(*p21));
  em = malloc(pair.n_slot * sizeof(*em));
  h1 = malloc(M * sizeof(*h1));
  h2 = malloc(M * sizeof(*h2));
  tmp.orientation = malloc(M);
  tmp.missing1 = malloc(M);
  tmp.missing2 = malloc(M);
  if (!p12 || !p21 || !em || !h1 || !h2 || !tmp.orientation ||
      !tmp.missing1 || !tmp.missing2)
    goto cleanup;
  memset(tmp.orientation, 255, M);
  memset(tmp.missing1, 255, M);
  memset(tmp.missing2, 255, M);
  for (size_t m = 0; m < M; m++)
  {
    size_t g = input.marker_start + m;
    if (0 != st_phase_pair_emission_row_neutral(
                 &pair, m, emit_ref + m * R, R, input.observed1[m],
                 input.observed2[m], d->n_alleles[g], e, p12 + m * P,
                 p21 + m * P, em + m * P))
      goto cleanup;
  }
  if (M > 1u)
  {
    q1 = malloc((M - 1u) * sizeof(*q1));
    q2 = malloc((M - 1u) * sizeof(*q2));
    if (!q1 || !q2)
      goto cleanup;
    for (size_t m = 1; m < M; m++)
    {
      size_t g = input.marker_start + m;
      double dist = d->pos_cm[g] - d->pos_cm[g - 1u];
      if (0 != st_phase_hmm_break_probability(dist, br1, &q1[m - 1u]) ||
          0 != st_phase_hmm_break_probability(dist, br2, &q2[m - 1u]))
        goto cleanup;
    }
  }
  initial = malloc(P * sizeof(*initial));
  if (!initial)
    goto cleanup;
  for (size_t z = 0; z < P; z++)
    initial[z] = 0.0;
  {
    size_t n = 0;
    const st_phase_pair_state *r =
        st_phase_pair_state_get(&pair, 0, &n, NULL, NULL);
    if (!r || 0u == n || n > P)
      goto cleanup;
    for (size_t z = 0; z < n; z++)
      initial[z] = 1.0;
  }
  if (0 != st_phase_pair_hmm_viterbi_factorized_path(
               &pair, initial, em, q1, q2, h1, h2, &possible))
    goto cleanup;
  tmp.possible = possible;
  tmp.n_marker = M;
  if (possible)
    for (size_t m = 0; m < M; m++)
    {
      uint8_t a = input.observed1[m], b = input.observed2[m];
      if (ST_PHASE_STATE_MISSING == a || a == b)
        continue;
      size_t n = 0, slot = SIZE_MAX;
      const st_phase_pair_state *r =
          st_phase_pair_state_get(&pair, m, &n, NULL, NULL);
      if (!r)
        goto cleanup;
      for (size_t z = 0; z < n; z++)
        if (r[z].hap1 == h1[m] && r[z].hap2 == h2[m])
        {
          slot = z;
          break;
        }
      if (SIZE_MAX == slot)
        goto cleanup;
      double x = p12[m * P + slot], y = p21[m * P + slot],
             tol = 256.0 * DBL_EPSILON;
      if (x > y + tol)
        tmp.orientation[m] = 0u;
      else if (y > x + tol)
        tmp.orientation[m] = 1u;
    }
  if (possible)
    for (size_t m = 0; m < M; m++)
      if (ST_PHASE_STATE_MISSING == input.observed1[m])
      {
        uint8_t x = emit_ref[m * R + h1[m]], y = emit_ref[m * R + h2[m]];
        if (ST_PHASE_STATE_MISSING == x)
          x = rank_ref[m * R + h1[m]];
        if (ST_PHASE_STATE_MISSING == y)
          y = rank_ref[m * R + h2[m]];
        if (x < d->n_alleles[input.marker_start + m] &&
            y < d->n_alleles[input.marker_start + m])
        {
          tmp.missing1[m] = x;
          tmp.missing2[m] = y;
        }
      }
  st_phase_cohort_viterbi_chromosome_result_free(out);
  *out = tmp;
  st_phase_cohort_viterbi_chromosome_result_init(&tmp);
  status = 0;
cleanup:
  free(p12);
  free(p21);
  free(em);
  free(q1);
  free(q2);
  free(initial);
  free(h1);
  free(h2);
  st_phase_pair_state_set_free(&pair);
  st_phase_unphased_candidate_set_free(&cand);
  st_phase_diploid_input_free(&input);
  st_phase_cohort_viterbi_chromosome_result_free(&tmp);
  return status;
}

static int st_phase_cohort_rephase_once_viterbi_neutral(
    const st_phase_data *d, const st_phase_cohort_snapshot *frozen,
    size_t K, double e, double br, st_phase_cohort_snapshot *out,
    size_t *changed_out)
{
  st_phase_cohort_snapshot tmp;
  size_t changed = 0u;
  int failed = 0;
  st_phase_cohort_snapshot_init(&tmp);
  if (!d || !out || !changed_out || 0u == K ||
      !st_phase_cohort_snapshot_valid(d, frozen))
    return -1;
  if (0 != st_phase_cohort_snapshot_clone(frozen, &tmp))
    return -1;

#ifdef MT
#pragma omp parallel for schedule(static) reduction(+ : changed)
#endif
  for (size_t i = 0u; i < d->n_ind; i++)
  {
    int bad = 0;
    for (size_t c = 0u; c < d->n_chrom && !bad; c++)
    {
      uint8_t *proxy = NULL, *masked = NULL;
      size_t R = 0u, start = d->chrom_start[c], prev = SIZE_MAX;
      uint8_t pori = 255u;
      st_phase_cohort_viterbi_chromosome_result ph;
      st_phase_cohort_viterbi_chromosome_result_init(&ph);
      if (0 != st_phase_cohort_build_loo_reference_masked(
                   frozen, i, c, &proxy, &masked, &R) ||
          0u == R ||
          0 != st_phase_cohort_viterbi_chromosome_neutral(
                   d, i, c, proxy, masked, R, K, e, br, br, &ph))
      {
        free(proxy);
        free(masked);
        st_phase_cohort_viterbi_chromosome_result_free(&ph);
        bad = 1;
        break;
      }
      free(proxy);
      free(masked);
      if (!ph.possible)
      {
        st_phase_cohort_viterbi_chromosome_result_free(&ph);
        continue;
      }
      for (size_t local = 0u; local < ph.n_marker; local++)
      {
        size_t m = start + local, z = m * d->n_ind + i;
        if (ST_PHASE_STATE_MISSING == d->geno1[z] && ph.missing1 &&
            ph.missing2 && ph.missing1[local] <= 1u &&
            ph.missing2[local] <= 1u)
        {
          size_t hh = m * tmp.n_hap + 2u * i;
          tmp.hap[hh] = ph.missing1[local];
          tmp.hap[hh + 1u] = ph.missing2[local];
        }
      }
      for (size_t local = 0u; local < ph.n_marker; local++)
      {
        size_t m = start + local, z = m * d->n_ind + i;
        uint8_t a = d->geno1[z], b = d->geno2[z],
                ori = ph.orientation[local];
        if (ST_PHASE_STATE_MISSING == a || a == b)
          continue;
        if (SIZE_MAX == prev)
        {
          prev = m;
          pori = ori;
          continue;
        }
        {
          int oldsame, newsame;
          size_t hr = m * tmp.n_hap + 2u * i;
          uint8_t lf = tmp.hap[prev * tmp.n_hap + 2u * i];
          if (0 != st_phase_cohort_relative_relation_from_snapshot(
                       d, frozen, i, prev, m, &oldsame))
          {
            bad = 1;
            break;
          }
          newsame = (pori <= 1u && ori <= 1u) ? (pori == ori) : oldsame;
          if (newsame != oldsame)
            changed++;
          {
            int left12 = (lf == d->geno1[prev * d->n_ind + i]);
            if (!newsame)
              left12 = !left12;
            tmp.hap[hr] = left12 ? a : b;
            tmp.hap[hr + 1u] = left12 ? b : a;
          }
          prev = m;
          pori = ori;
        }
      }
      st_phase_cohort_viterbi_chromosome_result_free(&ph);
    }
    if (bad)
    {
#ifdef MT
#pragma omp atomic write
#endif
      failed = 1;
    }
  }
  if (failed || !st_phase_cohort_snapshot_valid(d, &tmp))
    goto fail;
  st_phase_cohort_snapshot_free(out);
  *out = tmp;
  st_phase_cohort_snapshot_init(&tmp);
  *changed_out = changed;
  return 0;
fail:
  st_phase_cohort_snapshot_free(&tmp);
  return -1;
}

typedef struct
{
  size_t donor, count;
  uint64_t length_sum;
} st_phase_cohort_donor_support;
static int st_phase_cohort_donor_support_better(
    const st_phase_cohort_donor_support *a,
    const st_phase_cohort_donor_support *b)
{
  if (a->count != b->count)
    return a->count > b->count;
  if (a->length_sum != b->length_sum)
    return a->length_sum > b->length_sum;
  return a->donor < b->donor;
}

static int st_phase_cohort_build_raw_subset_reference_dual(
    const st_phase_diploid_input *input, const uint8_t *n_alleles,
    const uint8_t *rank_ref, const uint8_t *hmm_ref, size_t R,
    size_t local_pool, size_t max_subset, uint8_t **subset_out,
    size_t *S_out)
{
  st_phase_unphased_candidate_set cand;
  st_phase_cohort_donor_support *sup = NULL;
  uint8_t *sub = NULL;
  size_t S, M;
  int status = -1;
  st_phase_unphased_candidate_set_init(&cand);
  if (!input || !n_alleles || !rank_ref || !hmm_ref || !subset_out ||
      !S_out || 0u == R || 0u == local_pool || 0u == max_subset)
    return -1;
  *subset_out = NULL;
  *S_out = 0u;
  M = input->n_marker;
  S = (max_subset < R) ? max_subset : R;
  if (0 != st_phase_unphased_candidate_build_cohort_balanced(
               input->observed1, input->observed2, n_alleles, rank_ref, M,
               R, local_pool, &cand))
    goto cleanup;
  sup = calloc(R, sizeof(*sup));
  if (!sup)
    goto cleanup;
  for (size_t r = 0; r < R; r++)
    sup[r].donor = r;
  for (size_t m = 0; m < M; m++)
  {
    size_t n = cand.count[m];
    const st_phase_candidate *row =
        st_phase_unphased_candidate_row_const(&cand, m);
    if (!row || 0u == n)
      goto cleanup;
    for (size_t z = 0; z < n; z++)
    {
      size_t r = row[z].hap;
      if (r >= R)
        goto cleanup;
      sup[r].count++;
      sup[r].length_sum += (uint64_t)row[z].ibs_length;
    }
  }
  for (size_t x = 1; x < R; x++)
  {
    st_phase_cohort_donor_support v = sup[x];
    size_t y = x;
    while (y > 0 && st_phase_cohort_donor_support_better(&v, &sup[y - 1]))
    {
      sup[y] = sup[y - 1];
      y--;
    }
    sup[y] = v;
  }
  if (M > SIZE_MAX / S || M * S > SIZE_MAX / sizeof(*sub))
    goto cleanup;
  sub = malloc(M * S);
  if (!sub)
    goto cleanup;
  for (size_t m = 0; m < M; m++)
    for (size_t z = 0; z < S; z++)
    {
      uint8_t a = hmm_ref[m * R + sup[z].donor];
      if (ST_PHASE_STATE_MISSING != a && a >= n_alleles[m])
        goto cleanup;
      sub[m * S + z] = a;
    }
  *subset_out = sub;
  *S_out = S;
  sub = NULL;
  status = 0;
cleanup:
  free(sub);
  free(sup);
  st_phase_unphased_candidate_set_free(&cand);
  return status;
}

typedef struct
{
  size_t n_marker;
  uint8_t *allele;
  double *confidence;
} st_phase_haploid_missing_posterior;

static void st_phase_haploid_missing_posterior_init(
    st_phase_haploid_missing_posterior *x)
{
  if (x)
    memset(x, 0, sizeof(*x));
}

static void st_phase_haploid_missing_posterior_free(
    st_phase_haploid_missing_posterior *x)
{
  if (!x)
    return;
  free(x->allele);
  free(x->confidence);
  st_phase_haploid_missing_posterior_init(x);
}

static int st_phase_haploid_emission_probability_posterior(
    uint8_t observed, uint8_t donor, uint8_t A, double e, double *p)
{
  if (!p || ST_PHASE_STATE_MISSING == A || !st_phase_hmm_error_valid(e))
    return -1;
  if (0u == A)
  {
    if (ST_PHASE_STATE_MISSING != observed)
      return -1;
    *p = 1.0;
    return 0;
  }
  if (ST_PHASE_STATE_MISSING != observed && observed >= A)
    return -1;
  if (ST_PHASE_STATE_MISSING != donor && donor >= A)
    return -1;
  if (ST_PHASE_STATE_MISSING == observed)
  {
    *p = 1.0;
    return 0;
  }
  if (ST_PHASE_STATE_MISSING == donor)
  {
    *p = 1.0 / (double)A;
    return 0;
  }
  return st_phase_hmm_hap_emission(observed, donor, A, e, p);
}

static int st_phase_haploid_missing_posterior_fb(
    const st_phase_data *d, size_t marker_start, size_t M,
    const uint8_t *observed, const uint8_t *ref, size_t S, double e,
    double br, double threshold, st_phase_haploid_missing_posterior *out)
{
  st_phase_haploid_missing_posterior tmp;
  double *f = NULL, *b = NULL, *emit = NULL, *z = NULL;
  double *q = NULL;
  int status = -1;

  st_phase_haploid_missing_posterior_init(&tmp);
  if (!d || !observed || !ref || !out || 0u == M || 0u == S ||
      marker_start + M > d->n_marker || !st_phase_hmm_error_valid(e) ||
      !isfinite(br) || br < 0.0 || !isfinite(threshold) ||
      threshold < 0.5 || threshold > 1.0)
    return -1;
  if (M > SIZE_MAX / S || M * S > SIZE_MAX / sizeof(*f))
    return -1;

  f = malloc(M * S * sizeof(*f));
  b = malloc(M * S * sizeof(*b));
  emit = malloc(S * sizeof(*emit));
  z = malloc(S * sizeof(*z));
  tmp.allele = malloc(M * sizeof(*tmp.allele));
  tmp.confidence = malloc(M * sizeof(*tmp.confidence));
  if (!f || !b || !emit || !z || !tmp.allele || !tmp.confidence)
    goto cleanup;

  for (size_t m = 0u; m < M; m++)
  {
    tmp.allele[m] = ST_PHASE_STATE_MISSING;
    tmp.confidence[m] = NAN;
  }

  if (M > 1u)
  {
    q = malloc((M - 1u) * sizeof(*q));
    if (!q)
      goto cleanup;
    for (size_t m = 1u; m < M; m++)
    {
      size_t g = marker_start + m;
      double dist = d->pos_cm[g] - d->pos_cm[g - 1u];
      if (0 != st_phase_hmm_break_probability(dist, br, &q[m - 1u]))
        goto cleanup;
    }
  }

  {
    uint8_t A = d->n_alleles[marker_start];
    double sum = 0.0;
    for (size_t r = 0u; r < S; r++)
    {
      double p;
      if (0 != st_phase_haploid_emission_probability_posterior(
                   observed[0], ref[r], A, e, &p))
        goto cleanup;
      f[r] = p / (double)S;
      sum += f[r];
    }
    if (!(sum > 0.0) || !isfinite(sum))
      goto cleanup;
    for (size_t r = 0u; r < S; r++)
      f[r] /= sum;
  }

  for (size_t m = 1u; m < M; m++)
  {
    size_t g = marker_start + m;
    uint8_t A = d->n_alleles[g];
    double sum_prev = 0.0, sum = 0.0;
    for (size_t r = 0u; r < S; r++)
      sum_prev += f[(m - 1u) * S + r];
    if (!(sum_prev > 0.0) || !isfinite(sum_prev))
      goto cleanup;
    for (size_t r = 0u; r < S; r++)
    {
      double ep, prior = (1.0 - q[m - 1u]) * f[(m - 1u) * S + r] +
                         (q[m - 1u] / (double)S) * sum_prev;
      if (0 != st_phase_haploid_emission_probability_posterior(
                   observed[m], ref[m * S + r], A, e, &ep))
        goto cleanup;
      f[m * S + r] = prior * ep;
      sum += f[m * S + r];
    }
    if (!(sum > 0.0) || !isfinite(sum))
      goto cleanup;
    for (size_t r = 0u; r < S; r++)
      f[m * S + r] /= sum;
  }

  for (size_t r = 0u; r < S; r++)
    b[(M - 1u) * S + r] = 1.0 / (double)S;
  for (size_t m = M - 1u; m > 0u; m--)
  {
    size_t g = marker_start + m;
    uint8_t A = d->n_alleles[g];
    double sumz = 0.0, sumb = 0.0;
    for (size_t s = 0u; s < S; s++)
    {
      double ep;
      if (0 != st_phase_haploid_emission_probability_posterior(
                   observed[m], ref[m * S + s], A, e, &ep))
        goto cleanup;
      z[s] = ep * b[m * S + s];
      sumz += z[s];
    }
    if (!(sumz > 0.0) || !isfinite(sumz))
      goto cleanup;
    for (size_t r = 0u; r < S; r++)
    {
      b[(m - 1u) * S + r] =
          (1.0 - q[m - 1u]) * z[r] + (q[m - 1u] / (double)S) * sumz;
      sumb += b[(m - 1u) * S + r];
    }
    if (!(sumb > 0.0) || !isfinite(sumb))
      goto cleanup;
    for (size_t r = 0u; r < S; r++)
      b[(m - 1u) * S + r] /= sumb;
  }

  for (size_t m = 0u; m < M; m++)
    if (ST_PHASE_STATE_MISSING == observed[m])
    {
      uint8_t A = d->n_alleles[marker_start + m];
      double mass[2] = {0.0, 0.0}, den = 0.0;
      if (A > 2u)
        goto cleanup;
      if (0u == A)
      {
        tmp.confidence[m] = NAN;
        tmp.allele[m] = ST_PHASE_STATE_MISSING;
        continue;
      }
      for (size_t r = 0u; r < S; r++)
      {
        double w = f[m * S + r] * b[m * S + r];
        uint8_t a = ref[m * S + r];
        if (!(w >= 0.0) || !isfinite(w))
          goto cleanup;
        den += w;
        if (ST_PHASE_STATE_MISSING == a)
        {
          for (uint8_t x = 0u; x < A; x++)
            mass[x] += w / (double)A;
        }
        else
        {
          if (a >= A)
            goto cleanup;
          mass[a] += w;
        }
      }
      if (!(den > 0.0) || !isfinite(den))
        goto cleanup;
      for (uint8_t x = 0u; x < A; x++)
        mass[x] /= den;
      {
        uint8_t best = 0u;
        double bestp = mass[0];
        if (A > 1u && mass[1] > bestp)
        {
          best = 1u;
          bestp = mass[1];
        }
        tmp.confidence[m] = bestp;
        if (bestp >= threshold &&
            (1u == A || fabs(mass[0] - mass[1]) > 256.0 * DBL_EPSILON))
          tmp.allele[m] = best;
      }
    }

  tmp.n_marker = M;
  st_phase_haploid_missing_posterior_free(out);
  *out = tmp;
  st_phase_haploid_missing_posterior_init(&tmp);
  status = 0;

cleanup:
  free(f);
  free(b);
  free(emit);
  free(z);
  free(q);
  st_phase_haploid_missing_posterior_free(&tmp);
  return status;
}

static int st_phase_cohort_build_loo_reference_posterior(
    const st_phase_cohort_snapshot *frozen, const uint8_t *hap_confident,
    size_t target, size_t chromosome, uint8_t **proxy_out,
    uint8_t **emit_out, size_t *R_out)
{
  st_phase_cohort_reference_descriptor *desc = NULL;
  uint8_t *proxy = NULL, *emit = NULL;
  size_t start, M, H, R, ncell;
  int status = -1;

  if (!frozen || !hap_confident || !proxy_out || !emit_out || !R_out ||
      target >= frozen->n_ind || chromosome >= frozen->n_chrom ||
      frozen->n_ind < 2u || !frozen->hap || !frozen->observed_missing)
    return -1;
  start = frozen->chrom[chromosome].start;
  M = frozen->chrom[chromosome].length;
  H = frozen->n_hap;
  if (0u == M || H != 2u * frozen->n_ind ||
      (frozen->n_ind - 1u) > SIZE_MAX / 2u)
    return -1;
  R = 2u * (frozen->n_ind - 1u);
  if (R > SIZE_MAX / sizeof(*desc) ||
      0 != st_phase_size_mul(M, R, &ncell) ||
      ncell > SIZE_MAX / sizeof(*proxy))
    return -1;

  desc = malloc(R * sizeof(*desc));
  proxy = malloc(ncell);
  emit = malloc(ncell);
  if (!desc || !proxy || !emit)
    goto cleanup;

  {
    size_t r = 0u;
    for (size_t i = 0u; i < frozen->n_ind; i++)
      if (i != target)
        for (size_t side = 0u; side < 2u; side++)
        {
          desc[r].snapshot = frozen;
          desc[r].marker_start = start;
          desc[r].n_marker = M;
          desc[r].n_hap = H;
          desc[r].hap_index = 2u * i + side;
          r++;
        }
    if (r != R)
      goto cleanup;
  }
  qsort(desc, R, sizeof(*desc), st_phase_cohort_reference_descriptor_cmp);

  for (size_t m = 0u; m < M; m++)
    for (size_t r = 0u; r < R; r++)
    {
      size_t hi = desc[r].hap_index;
      size_t di = hi / 2u;
      size_t gm = start + m;
      uint8_t x = frozen->hap[gm * H + hi];
      if (x > 1u)
        goto cleanup;
      proxy[m * R + r] = x;
      if (!frozen->observed_missing[gm * frozen->n_ind + di])
        emit[m * R + r] = x;
      else
        emit[m * R + r] =
            hap_confident[gm * H + hi] ? x : ST_PHASE_STATE_MISSING;
    }

  *proxy_out = proxy;
  *emit_out = emit;
  *R_out = R;
  proxy = NULL;
  emit = NULL;
  status = 0;
cleanup:
  free(desc);
  free(proxy);
  free(emit);
  return status;
}

static void
st_phase_cohort_confidence_initialize(const st_phase_cohort_snapshot *s,
                                      uint8_t *conf)
{
  size_t H = s->n_hap;
  for (size_t m = 0u; m < s->n_marker; m++)
    for (size_t i = 0u; i < s->n_ind; i++)
    {
      uint8_t v = s->observed_missing[m * s->n_ind + i] ? 0u : 1u;
      conf[m * H + 2u * i] = v;
      conf[m * H + 2u * i + 1u] = v;
    }
}

static size_t st_phase_cohort_count_confident_missing_haps(
    const st_phase_cohort_snapshot *s, const uint8_t *conf)
{
  size_t n = 0u, H = s->n_hap;
  for (size_t m = 0u; m < s->n_marker; m++)
    for (size_t i = 0u; i < s->n_ind; i++)
      if (s->observed_missing[m * s->n_ind + i])
      {
        n += conf[m * H + 2u * i] ? 1u : 0u;
        n += conf[m * H + 2u * i + 1u] ? 1u : 0u;
      }
  return n;
}

static int st_phase_cohort_rephase_once_haploid_subset_posterior(
    const st_phase_data *d, const st_phase_cohort_snapshot *frozen,
    const uint8_t *frozen_conf, size_t local_pool, size_t max_subset,
    double e, double br, double posterior_threshold,
    double phase_threshold, st_phase_cohort_snapshot *out,
    uint8_t *out_conf, uint8_t *link_conf, size_t *changed_link_out,
    size_t *confident_missing_hap_out)
{
  st_phase_cohort_snapshot tmp;
  size_t changed = 0u, confident_missing = 0u;
  int failed = 0;

  st_phase_cohort_snapshot_init(&tmp);
  if (!d || !frozen || !frozen_conf || !out || !out_conf ||
      !changed_link_out || !confident_missing_hap_out ||
      0u == local_pool || 0u == max_subset ||
      !st_phase_cohort_snapshot_valid(d, frozen) ||
      !st_phase_hmm_error_valid(e) || !isfinite(br) || br < 0.0 ||
      !isfinite(posterior_threshold) || posterior_threshold < 0.5 ||
      posterior_threshold > 1.0 || !isfinite(phase_threshold) ||
      phase_threshold < 0.5 || phase_threshold > 1.0 || !link_conf)
    return -1;

  if (0 != st_phase_cohort_snapshot_clone(frozen, &tmp))
    return -1;
  st_phase_cohort_confidence_initialize(frozen, out_conf);
  memset(link_conf, 0, d->n_marker * d->n_ind * sizeof(*link_conf));

#ifdef MT
#pragma omp parallel for schedule(static) reduction(+ : changed)
#endif
  for (size_t i = 0u; i < d->n_ind; i++)
  {
    int bad = 0;
    for (size_t c = 0u; c < d->n_chrom && !bad; c++)
    {
      uint8_t *proxy = NULL, *emit = NULL, *sub = NULL, *obs = NULL;
      size_t R = 0u, S = 0u, start = d->chrom_start[c],
             M = d->chrom_length[c], end = start + M;
      size_t left = SIZE_MAX;
      st_phase_diploid_input input;
      st_phase_cohort_haplink_result hr;
      st_phase_haploid_missing_posterior mp;
      st_phase_diploid_input_init(&input);
      st_phase_cohort_haplink_result_init(&hr);
      st_phase_haploid_missing_posterior_init(&mp);

      if (0 != st_phase_cohort_build_loo_reference_posterior(
                   frozen, frozen_conf, i, c, &proxy, &emit, &R) ||
          0u == R ||
          0 != st_phase_prepare_diploid_input(d, i, c, &input) ||
          0 != st_phase_cohort_build_raw_subset_reference_dual(
                   &input, d->n_alleles + input.marker_start, proxy, emit,
                   R, local_pool, max_subset, &sub, &S) ||
          0u == S ||
          0 != st_phase_cohort_haploid_links_chromosome(d, i, c, sub, S, e,
                                                        br, &hr))
      {
        free(proxy);
        free(emit);
        free(sub);
        free(obs);
        st_phase_diploid_input_free(&input);
        st_phase_cohort_haplink_result_free(&hr);
        st_phase_haploid_missing_posterior_free(&mp);
        bad = 1;
        break;
      }
      free(proxy);
      free(emit);
      proxy = NULL;
      emit = NULL;

      for (size_t m = start; m < end; m++)
        if (st_phase_cohort_observed_heterozygous(d, m, i))
        {
          if (SIZE_MAX == left)
          {
            left = m;
            continue;
          }
          {
            size_t local = m - start, h = m * tmp.n_hap + 2u * i;
            int oldsame, left12;
            uint8_t a = d->geno1[m * d->n_ind + i],
                    b = d->geno2[m * d->n_ind + i];
            if (255u == hr.relation[local] ||
                0 != st_phase_cohort_relative_relation_from_snapshot(
                         d, frozen, i, left, m, &oldsame))
            {
              bad = 1;
              break;
            }
            if ((int)hr.relation[local] != oldsame)
              changed++;
            {
              double lr = hr.log_ratio[local];
              double p =
                  (isfinite(lr) && lr >= 0.0)
                      ? ((lr > 40.0) ? 1.0 : (1.0 / (1.0 + exp(-lr))))
                      : 0.5;
              if (p > 0.5 && p >= phase_threshold)
                link_conf[m * d->n_ind + i] = 1u;
            }
            left12 = (tmp.hap[left * tmp.n_hap + 2u * i] ==
                      d->geno1[left * d->n_ind + i]);
            if (!hr.relation[local])
              left12 = !left12;
            tmp.hap[h] = left12 ? a : b;
            tmp.hap[h + 1u] = left12 ? b : a;
            left = m;
          }
        }
      if (bad)
      {
        free(sub);
        st_phase_diploid_input_free(&input);
        st_phase_cohort_haplink_result_free(&hr);
        st_phase_haploid_missing_posterior_free(&mp);
        break;
      }

      obs = malloc(M * sizeof(*obs));
      if (!obs)
      {
        free(sub);
        st_phase_diploid_input_free(&input);
        st_phase_cohort_haplink_result_free(&hr);
        st_phase_haploid_missing_posterior_free(&mp);
        bad = 1;
        break;
      }
      for (size_t side = 0u; side < 2u; side++)
      {
        for (size_t lm = 0u; lm < M; lm++)
        {
          size_t g = start + lm, z = g * d->n_ind + i;
          obs[lm] = (ST_PHASE_STATE_MISSING == d->geno1[z])
                        ? ST_PHASE_STATE_MISSING
                        : tmp.hap[g * tmp.n_hap + 2u * i + side];
        }
        st_phase_haploid_missing_posterior_free(&mp);
        if (0 !=
            st_phase_haploid_missing_posterior_fb(
                d, start, M, obs, sub, S, e, br, posterior_threshold, &mp))
        {
          bad = 1;
          break;
        }
        for (size_t lm = 0u; lm < M; lm++)
        {
          size_t g = start + lm, z = g * d->n_ind + i,
                 hi = g * tmp.n_hap + 2u * i + side;
          if (ST_PHASE_STATE_MISSING == d->geno1[z])
          {
            if (ST_PHASE_STATE_MISSING != mp.allele[lm])
            {
              tmp.hap[hi] = mp.allele[lm];
              out_conf[hi] = 1u;
            }
            else
              out_conf[hi] = 0u;
          }
        }
      }
      free(obs);
      free(sub);
      st_phase_diploid_input_free(&input);
      st_phase_cohort_haplink_result_free(&hr);
      st_phase_haploid_missing_posterior_free(&mp);
    }
    if (bad)
    {
#ifdef MT
#pragma omp atomic write
#endif
      failed = 1;
    }
  }

  if (failed || !st_phase_cohort_snapshot_valid(d, &tmp))
    goto fail;
  confident_missing =
      st_phase_cohort_count_confident_missing_haps(&tmp, out_conf);
  st_phase_cohort_snapshot_free(out);
  *out = tmp;
  st_phase_cohort_snapshot_init(&tmp);
  *changed_link_out = changed;
  *confident_missing_hap_out = confident_missing;
  return 0;
fail:
  st_phase_cohort_snapshot_free(&tmp);
  return -1;
}

static int st_phase_data_subset_marker_mask(const st_phase_data *src,
                                            const uint8_t *keep,
                                            st_phase_data *out,
                                            size_t **full_index_out)
{
  st_phase_data tmp;
  size_t *map = NULL;
  size_t M = 0u, C = 0u, alleles = 0u, ncell = 0u, k = 0u;
  int last_chrom = INT_MIN;
  st_phase_data_init(&tmp);
  if (!src || !keep || !out || !full_index_out ||
      !st_phase_data_adapter_valid(src))
    return -1;
  *full_index_out = NULL;
  for (size_t m = 0u; m < src->n_marker; m++)
    if (keep[m])
    {
      M++;
      if (src->chrom_label[m] != last_chrom)
      {
        C++;
        last_chrom = src->chrom_label[m];
      }
      if (alleles > SIZE_MAX - (size_t)src->n_alleles[m])
        goto fail;
      alleles += (size_t)src->n_alleles[m];
    }
  if (0u == M || 0u == C || 0 != st_phase_size_mul(M, src->n_ind, &ncell))
    goto fail;
  tmp.n_ind = src->n_ind;
  tmp.n_marker = M;
  tmp.n_chrom = C;
  tmp.n_allele_value = alleles;
  tmp.source_marker = malloc(M * sizeof(*tmp.source_marker));
  tmp.chrom_label = malloc(M * sizeof(*tmp.chrom_label));
  tmp.pos_cm = malloc(M * sizeof(*tmp.pos_cm));
  tmp.chrom_start = malloc(C * sizeof(*tmp.chrom_start));
  tmp.chrom_length = malloc(C * sizeof(*tmp.chrom_length));
  tmp.allele_offset = malloc((M + 1u) * sizeof(*tmp.allele_offset));
  tmp.n_alleles = malloc(M * sizeof(*tmp.n_alleles));
  tmp.allele_value =
      (alleles ? malloc(alleles * sizeof(*tmp.allele_value)) : NULL);
  tmp.geno1 = malloc(ncell * sizeof(*tmp.geno1));
  tmp.geno2 = malloc(ncell * sizeof(*tmp.geno2));
  map = malloc(M * sizeof(*map));
  if (!tmp.source_marker || !tmp.chrom_label || !tmp.pos_cm ||
      !tmp.chrom_start || !tmp.chrom_length || !tmp.allele_offset ||
      !tmp.n_alleles || (alleles && !tmp.allele_value) || !tmp.geno1 ||
      !tmp.geno2 || !map)
    goto fail;
  tmp.allele_offset[0] = 0u;
  size_t ao = 0u;
  for (size_t m = 0u; m < src->n_marker; m++)
    if (keep[m])
    {
      size_t A = src->n_alleles[m], so = src->allele_offset[m];
      map[k] = m;
      tmp.source_marker[k] = src->source_marker[m];
      tmp.chrom_label[k] = src->chrom_label[m];
      tmp.pos_cm[k] = src->pos_cm[m];
      tmp.n_alleles[k] = src->n_alleles[m];
      for (size_t a = 0u; a < A; a++)
        tmp.allele_value[ao + a] = src->allele_value[so + a];
      ao += A;
      tmp.allele_offset[k + 1u] = ao;
      memcpy(tmp.geno1 + k * src->n_ind, src->geno1 + m * src->n_ind,
             src->n_ind * sizeof(*tmp.geno1));
      memcpy(tmp.geno2 + k * src->n_ind, src->geno2 + m * src->n_ind,
             src->n_ind * sizeof(*tmp.geno2));
      k++;
    }
  {
    size_t c = 0u, start = 0u;
    tmp.chrom_start[0] = 0u;
    for (size_t m = 1u; m < M; m++)
      if (tmp.chrom_label[m] != tmp.chrom_label[m - 1u])
      {
        tmp.chrom_length[c] = m - start;
        c++;
        start = m;
        tmp.chrom_start[c] = m;
      }
    tmp.chrom_length[c] = M - start;
  }
  if (k != M || ao != alleles || !st_phase_data_adapter_valid(&tmp))
    goto fail;
  st_phase_data_free(out);
  *out = tmp;
  st_phase_data_init(&tmp);
  *full_index_out = map;
  map = NULL;
  return 0;
fail:
  free(map);
  st_phase_data_free(&tmp);
  return -1;
}

static int st_phase_two_stage_highfreq_mask(const st_phase_data *d,
                                            size_t cutoff,
                                            uint8_t **mask_out,
                                            size_t *n_high_out,
                                            size_t *n_low_out)
{
  uint8_t *mask = NULL;
  size_t hi = 0u, lo = 0u;
  if (!d || !mask_out || !n_high_out || !n_low_out || 0u == cutoff ||
      !st_phase_data_adapter_valid(d))
    return -1;
  mask = calloc(d->n_marker, sizeof(*mask));
  if (!mask)
    return -1;
  for (size_t m = 0u; m < d->n_marker; m++)
  {
    size_t c0 = 0u, c1 = 0u;
    uint8_t A = d->n_alleles[m];
    if (A > 2u)
    {
      free(mask);
      return -1;
    }
    for (size_t i = 0u; i < d->n_ind; i++)
    {
      size_t z = m * d->n_ind + i;
      uint8_t a = d->geno1[z], b = d->geno2[z];
      if (ST_PHASE_STATE_MISSING == a)
        continue;
      c0 += (0u == a) + (0u == b);
      c1 += (1u == a) + (1u == b);
    }
    if (A == 2u && ((c0 < c1 ? c0 : c1) >= cutoff))
    {
      mask[m] = 1u;
      hi++;
    }
    else if (A == 2u && c0 > 0u && c1 > 0u)
      lo++;
  }
  *mask_out = mask;
  *n_high_out = hi;
  *n_low_out = lo;
  return hi ? 0 : -1;
}

static int st_phase_two_stage_scaffold_reference(
    const st_phase_cohort_snapshot *s, size_t target, size_t chromosome,
    uint8_t **ref_out, size_t **donor_ind_out, size_t *R_out)
{
  uint8_t *ref = NULL;
  size_t *di = NULL;
  size_t M, H, R, start, r = 0u;
  if (!s || !ref_out || !donor_ind_out || !R_out || target >= s->n_ind ||
      chromosome >= s->n_chrom || s->n_ind < 2u)
    return -1;
  M = s->chrom[chromosome].length;
  start = s->chrom[chromosome].start;
  H = s->n_hap;
  R = 2u * (s->n_ind - 1u);
  if (H != 2u * s->n_ind || 0u == M || M > SIZE_MAX / R)
    return -1;
  ref = malloc(M * R);
  di = malloc(R * sizeof(*di));
  if (!ref || !di)
  {
    free(ref);
    free(di);
    return -1;
  }
  for (size_t j = 0u; j < s->n_ind; j++)
    if (j != target)
      for (size_t side = 0u; side < 2u; side++)
      {
        di[r] = j;
        for (size_t m = 0u; m < M; m++)
          ref[m * R + r] = s->hap[(start + m) * H + 2u * j + side];
        r++;
      }
  if (r != R)
  {
    free(ref);
    free(di);
    return -1;
  }
  *ref_out = ref;
  *donor_ind_out = di;
  *R_out = R;
  return 0;
}

static int st_phase_two_stage_state_posterior(
    const st_phase_data *d, const st_phase_cohort_snapshot *s,
    size_t target, size_t chromosome, size_t side, const uint8_t *ref,
    size_t R, double e, double br, double **gamma_out)
{
  double *f = NULL, *b = NULL, *g = NULL, *z = NULL, *q = NULL;
  size_t M, start;
  if (!d || !s || !ref || !gamma_out || target >= d->n_ind ||
      chromosome >= d->n_chrom || side > 1u || 0u == R ||
      !st_phase_hmm_error_valid(e) || !isfinite(br) || br < 0.0)
    return -1;
  M = d->chrom_length[chromosome];
  start = d->chrom_start[chromosome];
  if (0u == M || M > SIZE_MAX / R || M * R > SIZE_MAX / sizeof(double))
    return -1;
  f = malloc(M * R * sizeof(*f));
  b = malloc(M * R * sizeof(*b));
  g = malloc(M * R * sizeof(*g));
  z = malloc(R * sizeof(*z));
  if (M > 1u)
    q = malloc((M - 1u) * sizeof(*q));
  if (!f || !b || !g || !z || (M > 1u && !q))
    goto fail;
  for (size_t m = 1u; m < M; m++)
    if (0 != st_phase_hmm_break_probability(d->pos_cm[start + m] -
                                                d->pos_cm[start + m - 1u],
                                            br, &q[m - 1u]))
      goto fail;

  {
    double sum = 0.0;
    uint8_t obs = s->hap[start * s->n_hap + 2u * target + side],
            A = d->n_alleles[start];
    for (size_t r = 0u; r < R; r++)
    {
      double ep;
      if (0 != st_phase_haploid_emission_probability_posterior(obs, ref[r],
                                                               A, e, &ep))
        goto fail;
      f[r] = ep / (double)R;
      sum += f[r];
    }
    if (!(sum > 0.0))
      goto fail;
    for (size_t r = 0u; r < R; r++)
      f[r] /= sum;
  }
  for (size_t m = 1u; m < M; m++)
  {
    double sumprev = 0.0, sum = 0.0;
    uint8_t obs = s->hap[(start + m) * s->n_hap + 2u * target + side],
            A = d->n_alleles[start + m];
    for (size_t r = 0u; r < R; r++)
      sumprev += f[(m - 1u) * R + r];
    for (size_t r = 0u; r < R; r++)
    {
      double ep, pr = (1.0 - q[m - 1u]) * f[(m - 1u) * R + r] +
                      (q[m - 1u] / (double)R) * sumprev;
      if (0 != st_phase_haploid_emission_probability_posterior(
                   obs, ref[m * R + r], A, e, &ep))
        goto fail;
      f[m * R + r] = pr * ep;
      sum += f[m * R + r];
    }
    if (!(sum > 0.0))
      goto fail;
    for (size_t r = 0u; r < R; r++)
      f[m * R + r] /= sum;
  }

  for (size_t r = 0u; r < R; r++)
    b[(M - 1u) * R + r] = 1.0 / (double)R;
  for (size_t m = M - 1u; m > 0u; m--)
  {
    double sumz = 0.0, sumb = 0.0;
    uint8_t obs = s->hap[(start + m) * s->n_hap + 2u * target + side],
            A = d->n_alleles[start + m];
    for (size_t r = 0u; r < R; r++)
    {
      double ep;
      if (0 != st_phase_haploid_emission_probability_posterior(
                   obs, ref[m * R + r], A, e, &ep))
        goto fail;
      z[r] = ep * b[m * R + r];
      sumz += z[r];
    }
    if (!(sumz > 0.0))
      goto fail;
    for (size_t r = 0u; r < R; r++)
    {
      b[(m - 1u) * R + r] =
          (1.0 - q[m - 1u]) * z[r] + (q[m - 1u] / (double)R) * sumz;
      sumb += b[(m - 1u) * R + r];
    }
    if (!(sumb > 0.0))
      goto fail;
    for (size_t r = 0u; r < R; r++)
      b[(m - 1u) * R + r] /= sumb;
  }
  for (size_t m = 0u; m < M; m++)
  {
    double sum = 0.0;
    for (size_t r = 0u; r < R; r++)
    {
      g[m * R + r] = f[m * R + r] * b[m * R + r];
      sum += g[m * R + r];
    }
    if (!(sum > 0.0) || !isfinite(sum))
      goto fail;
    for (size_t r = 0u; r < R; r++)
      g[m * R + r] /= sum;
  }
  *gamma_out = g;
  g = NULL;
  free(f);
  free(b);
  free(z);
  free(q);
  return 0;
fail:
  free(f);
  free(b);
  free(g);
  free(z);
  free(q);
  return -1;
}

static int
st_phase_two_stage_find_scaffold_chrom(const st_phase_data *common,
                                       int label, size_t *c_out)
{
  if (!common || !c_out)
    return -1;
  for (size_t c = 0u; c < common->n_chrom; c++)
    if (common->chrom_label[common->chrom_start[c]] == label)
    {
      *c_out = c;
      return 0;
    }
  return -1;
}

static double
st_phase_two_stage_interp_weight(const st_phase_data *common, size_t cc,
                                 const size_t *common_to_full,
                                 const st_phase_data *full, size_t fm,
                                 size_t *left_out, size_t *right_out)
{
  size_t cs = common->chrom_start[cc], M = common->chrom_length[cc],
         q = 0u;
  while (q < M && common_to_full[cs + q] < fm)
    q++;
  if (0u == q)
  {
    *left_out = *right_out = cs;
    return 0.0;
  }
  if (q == M)
  {
    *left_out = *right_out = cs + M - 1u;
    return 0.0;
  }
  *left_out = cs + q - 1u;
  *right_out = cs + q;
  {
    double pl = common->pos_cm[*left_out], pr = common->pos_cm[*right_out],
           p = full->pos_cm[fm];
    if (!(pr > pl))
      return 0.0;
    if (p <= pl)
      return 0.0;
    if (p >= pr)
      return 1.0;
    return (p - pl) / (pr - pl);
  }
}

static int st_phase_two_stage_tie_coin(size_t individual, size_t marker,
                                       int chrom_label)
{
  uint64_t x = (uint64_t)individual + UINT64_C(0x9e3779b97f4a7c15);
  x ^= (uint64_t)marker + UINT64_C(0xbf58476d1ce4e5b9) + (x << 6) +
       (x >> 2);
  x ^= (uint64_t)(uint32_t)chrom_label + UINT64_C(0x94d049bb133111eb) +
       (x << 6) + (x >> 2);
  x ^= x >> 30;
  x *= UINT64_C(0xbf58476d1ce4e5b9);
  x ^= x >> 27;
  x *= UINT64_C(0x94d049bb133111eb);
  x ^= x >> 31;
  return (int)(x & UINT64_C(1));
}

static int st_phase_two_stage_lowfreq_insert(
    const st_phase_data *full, const uint8_t *high_mask,
    const st_phase_data *common, const size_t *common_to_full,
    const st_phase_cohort_snapshot *scaffold, double e, double br,
    double confidence_threshold, st_phase_cohort_snapshot *out,
    uint8_t *orientation_conf)
{
  st_phase_cohort_snapshot tmp;
  int failed = 0;
  st_phase_cohort_snapshot_init(&tmp);
  if (!full || !high_mask || !common || !common_to_full || !scaffold ||
      !out || !orientation_conf || !st_phase_hmm_error_valid(e) ||
      !isfinite(br) || br < 0.0 || !isfinite(confidence_threshold) ||
      confidence_threshold < 0.5 || confidence_threshold > 1.0 ||
      !st_phase_cohort_snapshot_valid(common, scaffold))
    return -1;
  memset(orientation_conf, 0,
         full->n_marker * full->n_ind * sizeof(*orientation_conf));
  if (0 != st_phase_cohort_snapshot_build(full, &tmp))
    return -1;

  for (size_t cm = 0u; cm < common->n_marker; cm++)
  {
    size_t fm = common_to_full[cm];
    for (size_t h = 0u; h < tmp.n_hap; h++)
      tmp.hap[fm * tmp.n_hap + h] =
          scaffold->hap[cm * scaffold->n_hap + h];
  }

#ifdef MT
#pragma omp parallel for schedule(dynamic) reduction(| : failed)
#endif
  for (size_t ii = 0u; ii < full->n_ind; ii++)
  {
    size_t i = (size_t)ii;
    for (size_t fc = 0u; fc < full->n_chrom; fc++)
    {
      size_t cc;
      uint8_t *ref = NULL;
      size_t *donor_ind = NULL;
      size_t R = 0u;
      double *g0 = NULL, *g1 = NULL;
      int label = full->chrom_label[full->chrom_start[fc]];
      int needs_phase = 0;

      for (size_t fm = full->chrom_start[fc];
           fm < full->chrom_start[fc] + full->chrom_length[fc]; fm++)
      {
        size_t z = fm * full->n_ind + i;
        if (!high_mask[fm] && full->n_alleles[fm] == 2u &&
            full->geno1[z] != ST_PHASE_STATE_MISSING &&
            full->geno1[z] != full->geno2[z])
        {
          needs_phase = 1;
          break;
        }
      }
      if (!needs_phase)
        continue;

      if (0 != st_phase_two_stage_find_scaffold_chrom(common, label, &cc))
        continue;
      if (0 != st_phase_two_stage_scaffold_reference(scaffold, i, cc, &ref,
                                                     &donor_ind, &R) ||
          0 != st_phase_two_stage_state_posterior(
                   common, scaffold, i, cc, 0u, ref, R, e, br, &g0) ||
          0 != st_phase_two_stage_state_posterior(common, scaffold, i, cc,
                                                  1u, ref, R, e, br, &g1))
      {
        free(ref);
        free(donor_ind);
        free(g0);
        free(g1);
        failed = 1;
        break;
      }
      size_t fs = full->chrom_start[fc], fe = fs + full->chrom_length[fc],
             cs = common->chrom_start[cc];
      for (size_t fm = fs; fm < fe; fm++)
        if (!high_mask[fm] && full->n_alleles[fm] == 2u)
        {
          size_t z = fm * full->n_ind + i;
          uint8_t a = full->geno1[z], b = full->geno2[z];
          if (ST_PHASE_STATE_MISSING == a || a == b)
            continue;
          size_t c0 = 0u, c1 = 0u;
          for (size_t j = 0u; j < full->n_ind; j++)
          {
            size_t q = fm * full->n_ind + j;
            uint8_t x = full->geno1[q], y = full->geno2[q];
            if (ST_PHASE_STATE_MISSING == x)
              continue;
            c0 += (0u == x) + (0u == y);
            c1 += (1u == x) + (1u == y);
          }
          uint8_t minor = (c1 < c0) ? 1u : 0u;
          size_t Lc, Rc;
          double t = st_phase_two_stage_interp_weight(
              common, cc, common_to_full, full, fm, &Lc, &Rc);
          size_t l = Lc - cs, r = Rc - cs;
          double pm[2] = {0.0, 0.0};
          for (size_t state = 0u; state < R; state++)
          {
            double w0 = (Lc == Rc) ? g0[l * R + state]
                                   : (1.0 - t) * g0[l * R + state] +
                                         t * g0[r * R + state];
            double w1 = (Lc == Rc) ? g1[l * R + state]
                                   : (1.0 - t) * g1[l * R + state] +
                                         t * g1[r * R + state];
            size_t j = donor_ind[state], q = fm * full->n_ind + j;
            uint8_t x = full->geno1[q], y = full->geno2[q];
            double carries;
            if (ST_PHASE_STATE_MISSING == x)
              carries = 0.5;
            else if (x == y)
              carries = (x == minor) ? 1.0 : 0.0;
            else
              carries = 1.0;
            pm[0] += w0 * carries;
            pm[1] += w1 * carries;
          }
          {
            double pa0 = (a == minor) ? pm[0] : (1.0 - pm[0]),
                   pb0 = (b == minor) ? pm[0] : (1.0 - pm[0]);
            double pa1 = (a == minor) ? pm[1] : (1.0 - pm[1]),
                   pb1 = (b == minor) ? pm[1] : (1.0 - pm[1]);
            double lab = pa0 * pb1, lba = pb0 * pa1, sum = lab + lba;
            size_t h = fm * tmp.n_hap + 2u * i;
            if (lab > lba)
            {
              tmp.hap[h] = a;
              tmp.hap[h + 1u] = b;
            }
            else if (lba > lab)
            {
              tmp.hap[h] = b;
              tmp.hap[h + 1u] = a;
            }
            else if (st_phase_two_stage_tie_coin(i, fm, label))
            {
              tmp.hap[h] = a;
              tmp.hap[h + 1u] = b;
            }
            else
            {
              tmp.hap[h] = b;
              tmp.hap[h + 1u] = a;
            }
            if (sum > 0.0 && isfinite(sum) && lab != lba)
            {
              double pc = ((lab > lba) ? lab : lba) / sum;
              if (pc > 0.5 && pc >= confidence_threshold)
                orientation_conf[z] = 1u;
            }
          }
        }
      free(ref);
      free(donor_ind);
      free(g0);
      free(g1);
    }
  }
  if (failed || !st_phase_cohort_snapshot_valid(full, &tmp))
  {
    st_phase_cohort_snapshot_free(&tmp);
    return -1;
  }
  st_phase_cohort_snapshot_free(out);
  *out = tmp;
  st_phase_cohort_snapshot_init(&tmp);
  return 0;
}

static size_t st_phase_two_stage_relative_links(const st_phase_data *d)
{
  size_t links = 0u;
  if (!st_phase_data_adapter_valid(d))
    return 0u;
  for (size_t i = 0u; i < d->n_ind; i++)
    for (size_t c = 0u; c < d->n_chrom; c++)
    {
      size_t nhet = 0u, start = d->chrom_start[c],
             end = start + d->chrom_length[c];
      for (size_t m = start; m < end; m++)
        if (st_phase_cohort_observed_heterozygous(d, m, i))
          nhet++;
      if (nhet > 1u)
        links += nhet - 1u;
    }
  return links;
}

static int st_phase_two_stage_stage_final(
    const gs_varset_type *gv, const st_phase_data *d,
    const st_phase_cohort_snapshot *phase, const uint8_t *high_mask,
    const st_phase_data *common, const size_t *common_to_full,
    const uint8_t *common_link_conf,
    const uint8_t *lowfreq_orientation_conf, st_phase_dataset_stage *stage,
    st_phase_two_stage_pipeline_diagnostics *diag)
{
  size_t *full_to_common = NULL;
  if ((NULL == gv) || (NULL == d) || (NULL == phase) ||
      (NULL == high_mask) || (NULL == common) ||
      (NULL == common_to_full) || (NULL == common_link_conf) ||
      (NULL == lowfreq_orientation_conf) || (NULL == stage) ||
      (NULL == diag) || !st_phase_cohort_snapshot_valid(d, phase))
    return -1;
  if (ST_PHASE_STAGE_OK != st_phase_dataset_stage_build(gv, d, stage))
    return -1;

  full_to_common = malloc(d->n_marker * sizeof(*full_to_common));
  if (!full_to_common)
    return -1;
  for (size_t m = 0u; m < d->n_marker; m++)
    full_to_common[m] = SIZE_MAX;
  for (size_t cm = 0u; cm < common->n_marker; cm++)
  {
    size_t fm = common_to_full[cm];
    if (fm >= d->n_marker || full_to_common[fm] != SIZE_MAX)
    {
      free(full_to_common);
      return -1;
    }
    full_to_common[fm] = cm;
  }

  for (size_t m = 0u; m < d->n_marker; m++)
    for (size_t i = 0u; i < d->n_ind; i++)
    {
      size_t raw = (size_t)d->source_marker[m] * d->n_ind + i;
      stage->covered[raw] = 1u;
    }
  stage->n_covered = stage->n_genotype;

  for (size_t i = 0u; i < d->n_ind; i++)
    for (size_t c = 0u; c < d->n_chrom; c++)
    {
      size_t start = d->chrom_start[c], end = start + d->chrom_length[c];
      size_t nhet = 0u, nlink = 0u, nfixed = 0u;
      size_t previous_common_het = SIZE_MAX;
      int ready = 1;
      diag->n_target_chromosome++;

      for (size_t m = start; m < end; m++)
        if (st_phase_cohort_observed_heterozygous(d, m, i))
        {
          nhet++;
          if (high_mask[m])
          {
            size_t cm = full_to_common[m];
            if (SIZE_MAX == cm)
            {
              free(full_to_common);
              return -1;
            }
            if (SIZE_MAX != previous_common_het)
            {
              if (common_link_conf[cm * d->n_ind + i])
                nfixed++;
              else
                ready = 0;
            }
            previous_common_het = cm;
          }
          else
          {
            if (lowfreq_orientation_conf[m * d->n_ind + i])
              nfixed++;
            else
              ready = 0;
          }
        }

      nlink = (nhet > 1u) ? nhet - 1u : 0u;
      diag->n_observed_heterozygous += nhet;
      diag->n_valid_link += nlink;
      if (0u == nhet)
      {
        diag->n_no_heterozygote_chromosome++;
        diag->n_phase_ready_chromosome++;
        continue;
      }

      {
        size_t left = SIZE_MAX;
        diag->n_phase_ready_chromosome++;
        diag->n_confident_heterozygous += nhet;
        diag->n_confident_link += nlink;
        if (!ready && nfixed > 0u)
          diag->n_partial_fixed_chromosome++;
        for (size_t m = start; m < end; m++)
          if (st_phase_cohort_observed_heterozygous(d, m, i))
          {
            size_t z = m * d->n_ind + i;
            size_t h = m * phase->n_hap + 2u * i;
            uint8_t a = d->geno1[z], b = d->geno2[z];
            uint8_t x = phase->hap[h], y = phase->hap[h + 1u];
            size_t raw = (size_t)d->source_marker[m] * d->n_ind + i;
            size_t off = d->allele_offset[m];
            if (!(((x == a) && (y == b)) || ((x == b) && (y == a))) ||
                x >= d->n_alleles[m] || y >= d->n_alleles[m])
            {
              free(full_to_common);
              return -1;
            }
            if (SIZE_MAX != left)
            {
              int same = 0;
              if (0 != st_phase_cohort_relative_relation_from_snapshot(
                           d, phase, i, left, m, &same))
              {
                free(full_to_common);
                return -1;
              }
              if (same)
                diag->n_fixed_same++;
              else
                diag->n_fixed_switch++;
            }
            left = m;
            stage->staged1[raw] = d->allele_value[off + (size_t)x];
            stage->staged2[raw] = d->allele_value[off + (size_t)y];
            stage->claimed[raw] = 1u;
            stage->n_written++;
          }
      }
    }

  diag->n_unresolved_link = diag->n_valid_link - diag->n_confident_link;
  free(full_to_common);
  return 0;
}

static int st_phase_two_stage_impute_stage(
    const st_phase_data *d, const st_phase_cohort_snapshot *phase,
    size_t max_states, double error_rate, double break_rate1_per_morgan,
    double break_rate2_per_morgan, double confidence_threshold,
    double imputation_threshold, st_phase_dataset_stage *stage,
    st_phase_two_stage_pipeline_diagnostics *diag)
{
  typedef struct
  {
    st_phase_step0_imputation_result result;
    size_t missing;
    int status, has_reference;
  } st_phase_impute_job;
  st_phase_impute_job *jobs = NULL;
  size_t n_job, batch_size = 1u, n_worker, longest = 0u;
  int status = -1;

  if ((NULL == d) || (NULL == phase) || (NULL == stage) ||
      (NULL == diag) || !st_phase_cohort_snapshot_valid(d, phase) ||
      !st_phase_dataset_stage_structurally_valid(stage) ||
      0u == max_states ||
      !st_phase_step0_imputation_threshold_valid(imputation_threshold) ||
      0 != st_phase_size_mul(d->n_ind, d->n_chrom, &n_job))
    return -1;

#ifdef MT
  {
    int requested = omp_get_max_threads();
    if (requested > 1)
      batch_size = (size_t)requested;
  }
#endif
  for (size_t c = 0u; c < d->n_chrom; c++)
    if (d->chrom_length[c] > longest)
      longest = d->chrom_length[c];
  {

    double bytes = (double)longest *
        (128.0 * (double)max_states * (double)max_states +
         8.0 * (double)d->n_ind + 128.0 * (double)max_states + 1024.0);
    double allowed = bytes > 0.0 ? floor(2147483648.0 / bytes) : 1.0;
    if (allowed < 1.0)
      allowed = 1.0;
    if ((double)batch_size > allowed)
      batch_size = (size_t)allowed;
  }
  if (batch_size > n_job)
    batch_size = n_job;
  n_worker = batch_size;

  if (n_worker > 1u)
    batch_size = n_worker > n_job / 4u ? n_job : 4u * n_worker;
  if (0u == batch_size || batch_size > SIZE_MAX / sizeof(*jobs))
    return -1;
  jobs = calloc(batch_size, sizeof(*jobs));
  if (NULL == jobs)
    return -1;

  for (size_t first = 0u; first < n_job; )
  {
    size_t count = n_job - first;
    if (count > batch_size)
      count = batch_size;
#ifdef MT
#pragma omp parallel for schedule(dynamic, 1) num_threads((int)n_worker) if(n_worker > 1u)
#endif
    for (size_t j = 0u; j < count; j++)
    {
      size_t index = first + j, i = index / d->n_chrom,
             c = index % d->n_chrom, start = d->chrom_start[c],
             len = d->chrom_length[c], R = 0u;
      uint8_t *reference = NULL;
      st_phase_step0_hmm_work w;
      st_phase_step0_hmm_work_init(&w);
      st_phase_step0_imputation_result_init(&jobs[j].result);
      jobs[j].missing = 0u;
      jobs[j].status = 0;
      jobs[j].has_reference = 0;
      for (size_t m = start; m < start + len; m++)
        if (ST_PHASE_STATE_MISSING == d->geno1[m * d->n_ind + i])
          jobs[j].missing++;
      if (0u == jobs[j].missing)
        continue;
      if (0 != st_phase_cohort_build_loo_reference(phase, i, c, &reference, &R))
        jobs[j].status = -1;
      else if (R > 0u)
      {
        jobs[j].has_reference = 1;
        if (0 != st_phase_step0_hmm_build(
                     d, i, c, reference, R, max_states, NULL, 0u, error_rate,
                     break_rate1_per_morgan, break_rate2_per_morgan, &w) ||
            0 != st_phase_step0_decode_missing(
                     d, reference, R, &w, confidence_threshold,
                     imputation_threshold, &jobs[j].result))
          jobs[j].status = -1;
      }
      free(reference);
      st_phase_step0_hmm_work_free(&w);
    }

    for (size_t j = 0u; j < count; j++)
    {
      size_t index = first + j, i = index / d->n_chrom,
             c = index % d->n_chrom;
      st_phase_step0_imputation_result *imp = &jobs[j].result;
      diag->n_missing_requested += jobs[j].missing;
      if (jobs[j].status != 0)
        goto cleanup;
      if (!jobs[j].has_reference)
        diag->n_imputation_unresolved += jobs[j].missing;
      else
      {
        if (0 != st_phase_step0_stage_apply_imputation(d, i, c, imp, stage))
          goto cleanup;
        diag->n_imputed += imp->n_imputed;
        diag->n_imputation_unresolved += imp->n_unresolved;
        diag->n_imputed_hom0 += imp->n_imputed_hom0;
        diag->n_imputed_hom1 += imp->n_imputed_hom1;
        diag->n_imputed_heterozygous += imp->n_imputed_heterozygous;
        diag->n_imputed_heterozygous_phased += imp->n_imputed_heterozygous_phased;
        diag->n_imputed_heterozygous_unphased += imp->n_imputed_heterozygous_unphased;
      }
      st_phase_step0_imputation_result_free(imp);
    }
    first += count;
  }
  status = 0;
cleanup:
  for (size_t j = 0u; j < batch_size; j++)
    st_phase_step0_imputation_result_free(&jobs[j].result);
  free(jobs);
  return status;
}

static int st_phase_dataset_pipeline_cohort_two_stage(
    gs_varset_type *gv, size_t burnin_max, size_t burnin_min,
    double burnin_convergence, size_t progressive_max,
    size_t progressive_min, double progressive_convergence,
    size_t burnin_states, size_t local_states, size_t max_states,
    size_t stage1_minor_allele_count, double provisional_threshold,
    double error_rate, double break_rate1_per_morgan,
    double break_rate2_per_morgan, double confidence_threshold,
    int impute_requested, double imputation_threshold,
    int commit_requested,
    st_phase_two_stage_pipeline_diagnostics *diagnostics)
{
  st_phase_two_stage_pipeline_diagnostics diag;
  st_phase_data full, common;
  st_phase_cohort_snapshot initial, current, next, final;
  st_phase_dataset_stage stage;
  st_phase_step0_stage_audit audit;
  uint8_t *mask = NULL, *conf_a = NULL, *conf_b = NULL;
  uint8_t *link_conf = NULL, *lowfreq_conf = NULL;
  size_t *common_to_full = NULL;
  size_t stage1 = 0u, stage2 = 0u, H = 0u, nconf = 0u, nlinkconf = 0u,
         nlowconf = 0u;
  int status = ST_PHASE_COHORT_PIPELINE_INVALID;

  st_phase_two_stage_pipeline_diagnostics_init(&diag);
  st_phase_data_init(&full);
  st_phase_data_init(&common);
  st_phase_cohort_snapshot_init(&initial);
  st_phase_cohort_snapshot_init(&current);
  st_phase_cohort_snapshot_init(&next);
  st_phase_cohort_snapshot_init(&final);
  st_phase_dataset_stage_init(&stage);
  st_phase_step0_stage_audit_init(&audit);
  if (NULL == diagnostics)
    goto cleanup;
  diag.commit_requested = commit_requested;
  diag.impute_requested = impute_requested;
  diag.imputation_threshold = imputation_threshold;
  diag.burnin_states = burnin_states;
  diag.local_states = local_states;
  diag.max_states = max_states;
  diag.stage1_minor_allele_count = stage1_minor_allele_count;
  diag.provisional_threshold = provisional_threshold;

  if ((NULL == gv) || (0u == burnin_max) || (0u == burnin_min) ||
      (burnin_min > burnin_max) || (0u == progressive_max) ||
      (0u == progressive_min) || (progressive_min > progressive_max) ||
      (0u == burnin_states) || (0u == local_states) ||
      (0u == max_states) || (0u == stage1_minor_allele_count) ||
      !isfinite(burnin_convergence) || burnin_convergence < 0.0 ||
      burnin_convergence > 1.0 || !isfinite(progressive_convergence) ||
      progressive_convergence < 0.0 || progressive_convergence > 1.0 ||
      !isfinite(provisional_threshold) || provisional_threshold < 0.5 ||
      provisional_threshold > 1.0 ||
      !st_phase_hmm_error_valid(error_rate) ||
      !isfinite(break_rate1_per_morgan) || break_rate1_per_morgan < 0.0 ||
      !isfinite(break_rate2_per_morgan) || break_rate2_per_morgan < 0.0 ||
      confidence_threshold < 0.5 || confidence_threshold > 1.0 ||
      ((0 != impute_requested) && (1 != impute_requested)) ||
      !st_phase_step0_imputation_threshold_valid(imputation_threshold) ||
      ((0 != commit_requested) && (1 != commit_requested)))
    goto cleanup;

  if (0 != st_phase_build_data(gv, &full))
  {
    status = ST_PHASE_COHORT_PIPELINE_DATA_FAILED;
    goto cleanup;
  }
  diag.n_ind = full.n_ind;
  diag.n_marker = full.n_marker;
  diag.n_chrom = full.n_chrom;
  if (0 != st_phase_size_mul(full.n_ind, full.n_marker, &diag.n_genotype))
    goto cleanup;
  if (0 != st_phase_two_stage_highfreq_mask(&full,
                                            stage1_minor_allele_count,
                                            &mask, &stage1, &stage2) ||
      0 != st_phase_data_subset_marker_mask(&full, mask, &common,
                                            &common_to_full) ||
      0 != st_phase_cohort_snapshot_build(&common, &initial))
  {
    status = ST_PHASE_COHORT_PIPELINE_PHASING_FAILED;
    goto cleanup;
  }
  diag.stage1_markers = stage1;
  diag.stage2_markers = stage2;

  {
    const st_phase_cohort_snapshot *frozen = &initial;
    size_t links = st_phase_two_stage_relative_links(&common);
    diag.burnin_relative_link = links;
    for (size_t it = 0u; it < burnin_max; it++)
    {
      size_t changed = 0u;
      st_phase_cohort_snapshot_init(&next);
      if (0 != st_phase_cohort_rephase_once_viterbi_neutral(
                   &common, frozen, burnin_states, error_rate,
                   break_rate1_per_morgan, &next, &changed))
      {
        status = ST_PHASE_COHORT_PIPELINE_PHASING_FAILED;
        goto cleanup;
      }
      if (frozen == &current)
        st_phase_cohort_snapshot_free(&current);
      current = next;
      st_phase_cohort_snapshot_init(&next);
      frozen = &current;
      diag.burnin_iterations_completed = it + 1u;
      diag.burnin_last_changed_link = changed;
      diag.burnin_last_change_fraction =
          links ? (double)changed / (double)links : 0.0;
      if (diag.burnin_iterations_completed >= burnin_min &&
          (0u == links ||
           diag.burnin_last_change_fraction <= burnin_convergence))
        break;
    }
  }
  if (0u == diag.burnin_iterations_completed)
  {
    status = ST_PHASE_COHORT_PIPELINE_PHASING_FAILED;
    goto cleanup;
  }

  H = current.n_hap;
  if ((0u == H) || (common.n_marker > SIZE_MAX / H))
  {
    status = ST_PHASE_COHORT_PIPELINE_PHASING_FAILED;
    goto cleanup;
  }
  nconf = common.n_marker * H;
  conf_a = malloc(nconf * sizeof(*conf_a));
  conf_b = malloc(nconf * sizeof(*conf_b));
  if (0 != st_phase_size_mul(common.n_marker, common.n_ind, &nlinkconf) ||
      0 != st_phase_size_mul(full.n_marker, full.n_ind, &nlowconf))
  {
    status = ST_PHASE_COHORT_PIPELINE_PHASING_FAILED;
    goto cleanup;
  }
  link_conf = malloc(nlinkconf * sizeof(*link_conf));
  lowfreq_conf = malloc(nlowconf * sizeof(*lowfreq_conf));
  if ((NULL == conf_a) || (NULL == conf_b) || (NULL == link_conf) ||
      (NULL == lowfreq_conf))
  {
    status = ST_PHASE_COHORT_PIPELINE_PHASING_FAILED;
    goto cleanup;
  }
  memset(link_conf, 0, nlinkconf * sizeof(*link_conf));
  memset(lowfreq_conf, 0, nlowconf * sizeof(*lowfreq_conf));
  st_phase_cohort_confidence_initialize(&initial, conf_a);
  diag.progressive_relative_link =
      st_phase_two_stage_relative_links(&common);

  for (size_t it = 0u; it < progressive_max; it++)
  {
    size_t changed = 0u, conf_missing = 0u;
    uint8_t *swap;
    st_phase_cohort_snapshot_init(&next);
    if (0 != st_phase_cohort_rephase_once_haploid_subset_posterior(
                 &common, &current, conf_a, local_states, max_states,
                 error_rate, break_rate1_per_morgan, provisional_threshold,
                 confidence_threshold, &next, conf_b, link_conf, &changed,
                 &conf_missing))
    {
      status = ST_PHASE_COHORT_PIPELINE_PHASING_FAILED;
      goto cleanup;
    }
    st_phase_cohort_snapshot_free(&current);
    current = next;
    st_phase_cohort_snapshot_init(&next);
    swap = conf_a;
    conf_a = conf_b;
    conf_b = swap;
    diag.progressive_iterations_completed = it + 1u;
    diag.progressive_last_changed_link = changed;
    diag.n_confident_missing_haplotype = conf_missing;
    diag.progressive_last_change_fraction =
        diag.progressive_relative_link
            ? (double)changed / (double)diag.progressive_relative_link
            : 0.0;
    if (diag.progressive_iterations_completed >= progressive_min &&
        (0u == diag.progressive_relative_link ||
         diag.progressive_last_change_fraction <= progressive_convergence))
      break;
  }

  if (0 != st_phase_two_stage_lowfreq_insert(
               &full, mask, &common, common_to_full, &current, error_rate,
               break_rate2_per_morgan, confidence_threshold, &final,
               lowfreq_conf))
  {
    status = ST_PHASE_COHORT_PIPELINE_PHASING_FAILED;
    goto cleanup;
  }
  diag.phasing_succeeded = 1;

  if (0 != st_phase_two_stage_stage_final(gv, &full, &final, mask, &common,
                                          common_to_full, link_conf,
                                          lowfreq_conf, &stage, &diag))
  {
    status = ST_PHASE_COHORT_PIPELINE_STAGING_FAILED;
    goto cleanup;
  }
  if (impute_requested &&
      0 != st_phase_two_stage_impute_stage(
               &full, &final, max_states, error_rate,
               break_rate1_per_morgan, break_rate2_per_morgan,
               confidence_threshold, imputation_threshold, &stage, &diag))
  {
    status = ST_PHASE_COHORT_PIPELINE_STAGING_FAILED;
    goto cleanup;
  }

  diag.audit_status =
      st_phase_step0_stage_audit_readonly(gv, &full, &stage, &audit);
  if (ST_PHASE_AUDIT_OK != diag.audit_status)
  {
    status = ST_PHASE_COHORT_PIPELINE_STAGING_FAILED;
    goto cleanup;
  }
  diag.staging_succeeded = 1;
  diag.n_covered = audit.audit.n_covered;
  diag.n_uncovered = audit.audit.n_uncovered;
  diag.n_written = audit.audit.n_written;
  diag.n_covered_unwritten = audit.audit.n_covered_unwritten;
  diag.n_changed_genotypes =
      st_phase_dataset_stage_changed_genotypes(&stage);
  diag.n_missing = audit.audit.n_missing;
  diag.n_homozygous = audit.audit.n_homozygous;
  diag.n_heterozygous = audit.audit.n_heterozygous;
  diag.n_covered_heterozygous_unresolved =
      audit.audit.n_covered_heterozygous_unresolved;
  diag.n_written_same_source_order =
      audit.audit.n_written_same_source_order;
  diag.n_written_swapped_source_order =
      audit.audit.n_written_swapped_source_order;
  diag.audit_complete = audit.audit.complete;
  diag.audit_source_bytes_current = audit.audit.source_bytes_current;
  diag.audit_workspace_current = audit.audit.workspace_current;

  if (commit_requested)
  {
    diag.commit_attempted = 1;
    diag.commit_status =
        st_phase_step0_dataset_stage_commit(gv, &full, &stage);
    if (ST_PHASE_COMMIT_OK != diag.commit_status)
    {
      status = ST_PHASE_COHORT_PIPELINE_COMMIT_FAILED;
      goto cleanup;
    }
    diag.committed = 1;
  }
  status = ST_PHASE_COHORT_PIPELINE_OK;

cleanup:
  diag.pipeline_status = status;
  if (NULL != diagnostics)
    *diagnostics = diag;
  free(mask);
  free(common_to_full);
  free(conf_a);
  free(conf_b);
  free(link_conf);
  free(lowfreq_conf);
  st_phase_dataset_stage_free(&stage);
  st_phase_cohort_snapshot_free(&final);
  st_phase_cohort_snapshot_free(&next);
  st_phase_cohort_snapshot_free(&current);
  st_phase_cohort_snapshot_free(&initial);
  st_phase_data_free(&common);
  st_phase_data_free(&full);
  return status;
}

/* st.phase */
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
  const char *stages[] = {"phasing", "parameter estimation",
                           "common-marker phasing", "haplotype coding",
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
        fprintf(
            trace, " %d",
            (int)((w->comp[(size_t)c * in->words + (m >> 6)] >> (m & 63)) & 1));
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
    float joint = 0.0f, state_sum = 0.0f, fwd_sum = 0.0f, mismatch_sum = 0.0f;
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
    if (rows[i].count > 0 && rows[i].mismatch > 0 && isfinite(rows[i].mismatch))
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
      ci.donors[2 * s + h] = ctx->in->donors[(size_t)s * ci.H + 2 * sample + h];
  int K = cs_run(&ci, cw, sample, ctx->in->seed, NULL);
  for (int h = 0; h < 2; ++h)
    for (int m = 0; m < ci.M; ++m)
    {
      int allele = (ci.haps[(size_t)(2 * sample + h) * ci.words + (m >> 6)] >>
                    (m & 63)) &
                   1;
      for (int k = 0; k < K; ++k)
      {
        int ref = (cw->comp[(size_t)k * ci.words + (m >> 6)] >> (m & 63)) & 1;
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
    fprintf(trace, "FINAL %d %u %u\n", out.rounds, pe_float_bits(out.mismatch),
            pe_float_bits(out.intensity));
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
    if (prev.sample == next.sample && pos[next.start] - pos[prev.end] <= 4.0)
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
            ibs_cluster child = {ibs_calloc(parent.n, sizeof(int)), 0,
                                 parent.hom && code != 1};
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
            ibs_append(out.sample + s, (ibs_segment){t, start, end - 1});
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
      while (v->start > 0 && ibs_compatible(g, N, v->start - 1, s, v->sample))
        --v->start;
      while (v->end + 1 < M && ibs_compatible(g, N, v->end + 1, s, v->sample))
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
      while (v - u < candidates && (backward ? (step <= left || step <= right)
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
          if (index != i && !ibs_query(ibs, hap / 2, a[index] / 2, start, end))
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
            (int)((refs[(size_t)k * words + (start >> 6)] >> (start & 63)) & 1);
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
        mismatch[2 * plane + z] = same ? x : !hi_equal(b, ref, start, end);
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
      ba_backward(w->bwd + t * in->K, in->recomb[c], em, w->rows[t * in->C + c],
                  in->K);
    else
      w->sums[t] = ba_forward(w->fwd + t * in->K, w->sums[t], in->recomb[c], em,
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
  float *b1 = w->bhet1 + (size_t)u * in->K, *b2 = w->bhet2 + (size_t)u * in->K;
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
    fprintf(trace, "H %d %" PRIu32 " %" PRIu32 " %d %d %d\n", c, ba_bits(num),
            ba_bits(den), w->swap, w->swaps, w->types[c]);
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
        for (int m = c ? in->ends[c - 1] : 0; m < in->ends[end - 1]; ++m)
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
    for (int s = backward ? first - 1 : first; backward ? s >= start : s < end;
         s += backward ? -1 : 1)
    {
      pbwt_update(&w, out.codes + (size_t)s * H, out.sizes[s], s, backward, a,
                  d);
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
    st_donor_task t = {in, ibs, out, a, d, steps, tid, count, backward, team};
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
    for (int s = backward ? first - 1 : first; backward ? s >= start : s < end;
         s += backward ? -1 : 1)
    {
      pbwt_update(&w, out.codes + (size_t)s * H, out.sizes[s], s, backward, a,
                  d);
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
          memcpy(snap_a + (size_t)count * H, a, (size_t)H * sizeof(int));
          memcpy(snap_d + (size_t)count * (H + 1), d,
                 (size_t)(H + 1) * sizeof(int));
          steps[count++] = s;
          if (count == BLOCK)
          {
            st_donors_block(in, ibs, &out, snap_a, snap_d, steps, count,
                            backward, workers);
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
          ci.donors[2 * s + j] = donor.donors[(size_t)s * H + 2 * sample + j];
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
      next[(size_t)2 * sample * words + (m >> 6)] |= (uint64_t)a[m] << (m & 63);
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

  memcpy(next + (size_t)2 * in->T * words, in->haps + (size_t)2 * in->T * words,
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
    sq_task task = {in, donor, ci, next, trace, tid, in->T, team, counts + tid};
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
      if ((int)((in->haps[(size_t)h * words + (m >> 6)] >> (m & 63)) & 1) !=
          in->gt[(size_t)m * H + h])
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
      fprintf(trace, "ITER %d %u\n", in->it, pe_float_bits(in->threshold));
    it_donors donor = it_find_donors(in, ibs, trace);
    if (seq->em && (in->it == 0 || in->it < in->burnin))
    {
      int min_steps = (int)ceil((double)(1.0f / in->step));
      if (min_steps < 200)
        min_steps = 200;
      int *starts = ba_alloc(donor.S, sizeof(int));
      for (int s = 1; s < donor.S; ++s)
        starts[s] = donor.ends[s - 1];
      er_input est = {{in->M, 2 * in->N, donor.S, in->max_states, min_steps,
                       (in->M + 63) / 64, starts, NULL, in->haps},
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
        in->recomb[m] = (float)-expm1(-(double)e.intensity * seq->distance[m]);
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
        int count = vote(h, a[h], a[h + 1], order, rank, a, unph, 2 * N) +
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
            a[h + side] = fill(h + side, order, rank, a, unph, 2 * N);
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
  for (int list = in->list_start[step]; list < in->list_start[step + 1]; ++list)
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
    if (v + 1 == ni && (in->backward ? d[v + 1] >= limit : d[v + 1] <= limit))
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
        stop = in->backward ? start - 1 : end, delta = in->backward ? -1 : 1;
    for (int s = first; s != stop; s += delta)
    {
      pbwt_update(&w, in->codes + (size_t)s * H, in->sizes[s], s, in->backward,
                  a, d);
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
  ls_input li = {C, H, rp.S, Kmax, min_steps, 0, in->seed, starts, pair, phase};
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
        memcpy(mm + (size_t)m * K[bit], lw.mismatch + (size_t)m * Kmax, K[bit]);
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
                raw->gt[(size_t)m * H + h], raw->gt[(size_t)m * H + (h ^ 1)],
                h >= 2 * in->T, prob[bit][(size_t)left * K[bit] + j],
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
    memcpy(out + (size_t)m * H + 2 * in->T, raw->gt + (size_t)m * H + 2 * in->T,
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
        it_ints(trace, b ? "BWD" : "FWD", s, donors[b] + (size_t)s * rp.T,
                rp.T);
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
              size_t q = (size_t)prev * H + h, z = (size_t)m * H + h;
              int sameRev = (rev[q] < rev[q + 1]) == (rev[z] < rev[z + 1]);
              int sameFwd = (fwd[q] < fwd[q + 1]) == (a[h] < a[h + 1]);
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
                      cdf[1] ? (next_int((uint32_t)cdf[1]) >= (uint32_t)cdf[0])
                             : 0;
                else
                {
                  int aligned = side;
                  if (prev >= 0)
                  {
                    size_t q = (size_t)prev * H + h;
                    if ((rev[q] < rev[q + 1]) != (fwd[q] < fwd[q + 1]))
                      aligned ^= 1;
                  }
                  a[h + side] = rev[(size_t)m * H + h + aligned];
                }
              }
          }
        }
      memcpy((pass ? fwd : rev) + (size_t)m * H, a, (size_t)H * sizeof(int));
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
  double first = n ? gm_position(bp, cm, n, positions[0]) : positions[0] * 1e-6;
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
      it_ints(trace, "INITENDS", s, in->clusters[s].end, in->clusters[s].n);
      it_ints(trace, "INITTYPES", s, in->clusters[s].type, in->clusters[s].n);
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
    in->recomb[m] = (float)-expm1(-(double)seq->intensity * seq->distance[m]);
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
  float window = (float)options.window_cm, overlap = (float)options.overlap_cm;
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
      double cm = map_n ? gm_position(map_bp, map_cm, map_n, positions[next])
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
      double cm = map_n ? gm_position(map_bp, map_cm, map_n, positions[end - 1])
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
      while (keep > 0 && positions[start + keep - 1] == positions[start + keep])
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
      !isfinite(c->rare) || c->rare < 0 || c->rare > 1 || !isfinite(c->scale) ||
      c->scale <= 0 || !isfinite(c->ne) || c->ne <= 0 || !isfinite(c->buffer) ||
      c->buffer <= 0 || !isfinite(c->lr) || c->lr <= 0)
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

typedef struct
{
  SEXP chrom, cm, bp, a, b, na, options, result;
  int seed, windows, targets;
} st_context;
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
static SEXP st_compute(void *data)
{
  st_context *c = data;
  int M = Rf_length(c->chrom), N = INTEGER(Rf_getAttrib(c->a, R_DimSymbol))[1];
  atomic_store(&st_max_team, 1);
  j_random master;
  jr_seed(&master, (uint64_t)(int64_t)c->seed);
  int column = 0;
  for (int start = 0; start < M; ++column)
  {
    R_CheckUserInterrupt();
    const double *o = REAL(c->options) +
                      (XLENGTH(c->options) == 14 ? 0 : (R_xlen_t)14 * column);
    er_set_threads((int)o[13]);
    st_phase_threads = (int)o[13];
    int chrom = INTEGER(c->chrom)[start], end = start + 1;
    while (end < M && INTEGER(c->chrom)[end] == chrom)
      ++end;
    int m = end - start;
    np_input raw = {0};
    raw.M = m;
    raw.N = N;
    raw.T = c->targets;
    raw.seed = (uint64_t)(int64_t)c->seed;
    raw.burn = (int)o[2];
    raw.iterations = (int)o[3];
    raw.K = (int)o[4];
    raw.rare = (float)o[5];
    raw.scale = (float)o[6];
    raw.ne = (float)o[7];
    raw.buffer = (float)o[8];
    raw.lr = (float)o[9];
    raw.em = (int)o[10];
    raw.bp = INTEGER(c->bp) + start;
    raw.gt = st_zero((size_t)m * 2 * N, sizeof(int));
    int *output = st_zero((size_t)m * 2 * N, sizeof(int));
    for (int j = 0; j < m; ++j)
      for (int s = 0; s < N; ++s)
      {
        size_t q = start + j + (size_t)M * s, z = (size_t)j * 2 * N + 2 * s;
        raw.gt[z] = RAW(c->a)[q] == 255 ? -1 : RAW(c->a)[q];
        raw.gt[z + 1] = RAW(c->b)[q] == 255 ? -1 : RAW(c->b)[q];
      }
    int before = c->windows;
    int ok = st_phase_windowed(
        &raw, INTEGER(c->na) + start, raw.bp, REAL(c->cm) + start, m,
        (st_window_options){o[0], o[1], (int)o[12]}, o[11], &master, chrom,
        &c->windows, output, NULL);
    if (!ok)
      st_phase_error("Invalid data or options for chromosome %d", chrom);
    INTEGER(VECTOR_ELT(c->result, 4))[column] = c->windows - before;
    for (int j = 0; j < m; ++j)
      for (int s = 0; s < N; ++s)
      {
        size_t q = start + j + (size_t)M * s, z = (size_t)j * 2 * N + 2 * s;
        RAW(VECTOR_ELT(c->result, 0))[q] = (Rbyte)output[z];
        RAW(VECTOR_ELT(c->result, 1))[q] = (Rbyte)output[z + 1];
      }
    st_release(raw.gt);
    st_release(output);
    start = end;
  }
  INTEGER(VECTOR_ELT(c->result, 2))[0] = c->windows;
  INTEGER(VECTOR_ELT(c->result, 3))[0] = atomic_load(&st_max_team);
  return c->result;
}
static SEXP st_hd_phase_impl(SEXP chrom, SEXP cm, SEXP bp, SEXP a, SEXP b,
                             SEXP na, SEXP options, SEXP seed, int targets)
{
  if (TYPEOF(chrom) != INTSXP || TYPEOF(cm) != REALSXP ||
      TYPEOF(bp) != INTSXP || TYPEOF(a) != RAWSXP || TYPEOF(b) != RAWSXP ||
      TYPEOF(na) != INTSXP || TYPEOF(options) != REALSXP ||
      XLENGTH(options) < 14 || XLENGTH(options) % 14 ||
      TYPEOF(seed) != INTSXP || XLENGTH(seed) != 1 ||
      INTEGER(seed)[0] == NA_INTEGER)
    st_phase_error("Invalid phasing argument types");
  int M = Rf_length(chrom);
  SEXP dim = Rf_getAttrib(a, R_DimSymbol), dim2 = Rf_getAttrib(b, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 || TYPEOF(dim2) != INTSXP ||
      XLENGTH(dim2) != 2 || M < 2 || INTEGER(dim)[0] != M ||
      INTEGER(dim2)[0] != M || INTEGER(dim)[1] < 2 || INTEGER(dim)[1] > 5000 ||
      INTEGER(dim2)[1] != INTEGER(dim)[1] ||
      XLENGTH(a) != (R_xlen_t)M * INTEGER(dim)[1] || XLENGTH(b) != XLENGTH(a) ||
      XLENGTH(cm) != M || XLENGTH(bp) != M || XLENGTH(na) != M)
    st_phase_error("Invalid phasing dimensions");
  if (targets < 0)
    targets = INTEGER(dim)[1];
  if (targets < 1 || targets > INTEGER(dim)[1])
    st_phase_error("Invalid target/reference count");
  int chromosomes = 0;
  for (int m = 0; m < M; ++m)
  {
    if (INTEGER(chrom)[m] < 1 ||
        (m && INTEGER(chrom)[m] < INTEGER(chrom)[m - 1]) ||
        !R_FINITE(REAL(cm)[m]) || INTEGER(bp)[m] < 1 || INTEGER(na)[m] < 1 ||
        INTEGER(na)[m] > 100)
      st_phase_error("Invalid marker map or allele count");
    if (!m || INTEGER(chrom)[m] != INTEGER(chrom)[m - 1])
      ++chromosomes;
  }
  if (XLENGTH(options) != 14 && XLENGTH(options) != (R_xlen_t)14 * chromosomes)
    st_phase_error("Options need one column or one column per chromosome");
  for (R_xlen_t column = 0; column < XLENGTH(options) / 14; ++column)
  {
    const double *o = REAL(options) + 14 * column;
    for (int i = 0; i < 14; ++i)
      if (i != 0 && i != 11 && !R_FINITE(o[i]))
        st_phase_error("Nonfinite phasing option");
    int indices[] = {2, 3, 4, 10, 12, 13};
    for (int j = 0; j < 6; ++j)
    {
      double v = o[indices[j]];
      if (v < 0 || v > INT_MAX || v != floor(v))
        st_phase_error("Invalid integer phasing option");
    }
    if (o[2] < 1 || o[3] < 1 || o[4] < 2 || o[4] > 1000000 || o[5] <= 0 ||
        o[5] > .5 || o[6] <= 0 || o[7] <= 0 || o[8] <= 0 || o[9] < 1 ||
        o[10] > 1 || o[12] < 4 || o[13] < 1)
      st_phase_error("Invalid phasing option value");
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, Rf_allocMatrix(RAWSXP, M, INTEGER(dim)[1]));
  SET_VECTOR_ELT(out, 1, Rf_allocMatrix(RAWSXP, M, INTEGER(dim)[1]));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(0));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(1));
  SET_VECTOR_ELT(out, 4, Rf_allocVector(INTSXP, chromosomes));
  st_context context = {
      chrom, cm, bp, a, b, na, options, out, INTEGER(seed)[0], 0, targets};
  SEXP result = st_compute(&context);
  UNPROTECT(1);
  return result;
}

extern gs_varset_type *gs_fdta(char *set_name);

static SEXP phase_field(SEXP x, const char *name)
{
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (TYPEOF(x) != VECSXP || TYPEOF(names) != STRSXP)
    st_phase_error("Invalid phase snapshot");
  for (R_xlen_t j = 0; j < XLENGTH(x); ++j)
    if (!strcmp(CHAR(STRING_ELT(names, j)), name))
      return VECTOR_ELT(x, j);
  st_phase_error("Missing phase snapshot field: %s", name);
  return R_NilValue;
}
static SEXP st_phase_apply(SEXP before, SEXP after, SEXP dataset)
{
  if (TYPEOF(dataset) != STRSXP || XLENGTH(dataset) != 1 ||
      STRING_ELT(dataset, 0) == NA_STRING)
    st_phase_error("Invalid phase destination");
  gs_varset_type *gv = gs_fdta((char *)CHAR(STRING_ELT(dataset, 0)));
  SEXP id = phase_field(before, "ID"), ind = phase_field(before, "individual");
  SEXP a = phase_field(before, "allele1"), b = phase_field(before, "allele2");
  SEXP x = phase_field(after, "allele1"), y = phase_field(after, "allele2");
  size_t M = gv->NoMar, N = gv->NoInd;
  if (!M || !N || M > SIZE_MAX / N || TYPEOF(id) != STRSXP ||
      TYPEOF(ind) != STRSXP || XLENGTH(id) != (R_xlen_t)M ||
      XLENGTH(ind) != (R_xlen_t)N || !gv->All1 || !gv->All2)
    st_phase_error("Phase destination dimensions changed");
  size_t n = M * N;
  if (TYPEOF(a) != RAWSXP || TYPEOF(b) != RAWSXP || TYPEOF(x) != RAWSXP ||
      TYPEOF(y) != RAWSXP || XLENGTH(a) != (R_xlen_t)n ||
      XLENGTH(b) != (R_xlen_t)n || XLENGTH(x) != (R_xlen_t)n ||
      XLENGTH(y) != (R_xlen_t)n)
    st_phase_error("Invalid phase result dimensions");
  for (size_t i = 0; i < N; ++i)
    if (strcmp(gv->IndName[i], CHAR(STRING_ELT(ind, i))))
      st_phase_error("Phase destination individuals changed");
  GS_ALLELE_TYPE *s1 = (GS_ALLELE_TYPE *)R_alloc(n, sizeof(*s1));
  GS_ALLELE_TYPE *s2 = (GS_ALLELE_TYPE *)R_alloc(n, sizeof(*s2));
  for (size_t m = 0; m < M; ++m)
  {
    R_CheckUserInterrupt();
    if (strcmp(gv->MarName[m], CHAR(STRING_ELT(id, m))))
      st_phase_error("Phase destination markers changed");
    int coding[100], used[255] = {0};
    for (int j = 0; j < 100; ++j)
      coding[j] = -1;
    int maximum = -1;
    for (size_t i = 0; i < N; ++i)
    {
      size_t r = m + M * i, c = m * N + i;
      int old[2] = {RAW(a)[r], RAW(b)[r]};
      int actual[2] = {gv->All1[c], gv->All2[c]};
      for (int h = 0; h < 2; ++h)
      {
        if (old[h] == 255)
        {
          if (actual[h] != 255)
            st_phase_error("Phase destination calls changed");
          continue;
        }
        if (old[h] >= 100 || actual[h] == 255)
          st_phase_error("Invalid phase source coding");
        if (coding[old[h]] >= 0 && coding[old[h]] != actual[h])
          st_phase_error("Phase destination calls changed");
        coding[old[h]] = actual[h];
        used[actual[h]] = 1;
        if (actual[h] > maximum)
          maximum = actual[h];
      }
      int p = RAW(x)[r], q = RAW(y)[r];
      if ((p != 255 && p >= 100) || (q != 255 && q >= 100))
        st_phase_error("Invalid phase result allele");
      if (old[0] != 255 && old[1] != 255 &&
          !((p == old[0] && q == old[1]) || (p == old[1] && q == old[0])))
        st_phase_error("Phasing attempted to change an observed genotype");
    }

    for (size_t i = 0; i < N; ++i)
    {
      size_t r = m + M * i, c = m * N + i;
      int out[2] = {RAW(x)[r], RAW(y)[r]};
      for (int h = 0; h < 2; ++h)
        if (out[h] != 255 && coding[out[h]] < 0)
        {
          int code = maximum < 1 ? 1 : maximum + 1;
          if (code >= 255)
          {
            code = 1;
            while (code < 255 && used[code])
              ++code;
          }
          if (code >= 255)
            st_phase_error("No available internal allele code");
          coding[out[h]] = code;
          used[code] = 1;
          if (code > maximum)
            maximum = code;
        }
      s1[c] =
          out[0] == 255 ? GS_ALLELE_MISSING : (GS_ALLELE_TYPE)coding[out[0]];
      s2[c] =
          out[1] == 255 ? GS_ALLELE_MISSING : (GS_ALLELE_TYPE)coding[out[1]];
    }
  }

  memcpy(gv->All1, s1, n * sizeof(*s1));
  memcpy(gv->All2, s2, n * sizeof(*s2));
  gv->phased = 1;
  return R_NilValue;
}

extern SEXP st_markerdata_to_STvcf(SEXP dataset);

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

static void stp_number(SEXP options, int index, double lower, double upper,
                       int whole, int vector, int infinity)
{
  SEXP x = VECTOR_ELT(options, index);
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

static SEXP stp_options(SEXP input)
{
  SEXP names = Rf_getAttrib(input, R_NamesSymbol);
  if (TYPEOF(input) != VECSXP ||
      (XLENGTH(input) &&
       (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(input))))
    st_phase_error("phasing.options must be a uniquely named list");

  SEXP out = PROTECT(stp_list(STP_OPTIONS, stp_option_names));
  const double defaults[] = {0, 2, 4000000, 3, 12, 280, 100000, 0,
                             1, 100000, 3, .002, 1, 0};
  for (int i = 0; i < STP_OPTIONS; i++)
  {
    if (i == STP_ERR || i == STP_SEED)
      continue;
    if (i == STP_WINDOW)
      SET_VECTOR_ELT(out, i, Rf_mkString("auto"));
    else if (i == STP_EM)
      SET_VECTOR_ELT(out, i, Rf_ScalarLogical(1));
    else if (i >= STP_MARKERS && i <= STP_STATES)
      SET_VECTOR_ELT(out, i, Rf_ScalarInteger((int)defaults[i]));
    else
      SET_VECTOR_ELT(out, i, Rf_ScalarReal(defaults[i]));
  }
  int supplied[STP_OPTIONS] = {0};
  if (XLENGTH(input) && Rf_any_duplicated(names, FALSE))
    st_phase_error("phasing.options must be a uniquely named list");
  for (R_xlen_t i = 0; i < XLENGTH(input); i++)
  {
    if (STRING_ELT(names, i) == NA_STRING || !CHAR(STRING_ELT(names, i))[0])
      st_phase_error("phasing.options must be a uniquely named list");
    const char *name = CHAR(STRING_ELT(names, i));
    size_t size = strlen(name);
    char *normalized = (char *)R_alloc(size + 1, 1);
    for (size_t j = 0; j <= size; j++)
      normalized[j] = name[j] == '-' ? '.' : name[j];
    if (!strcmp(normalized, "nthreads"))
      st_phase_error("Set threads with st.set.num.threads(); nthreads is not a "
                     "phasing.options entry");
    if (!strcmp(normalized, "ref"))
      st_phase_error("Use the top-level ref argument");
    if (!strcmp(normalized, "impute"))
      st_phase_error("Use the top-level impute argument for missing calls at typed "
                     "markers");
    int index = -1;
    for (int j = 0; j < STP_OPTIONS; j++)
      if (!strcmp(normalized, stp_option_names[j]))
        index = j;
    if (index < 0)
      st_phase_error("Unsupported phasing option: %s. See help('st.phase').",
                     normalized);
    if (supplied[index])
      st_phase_error("Duplicate phasing option after name normalization");
    supplied[index] = 1;
    SET_VECTOR_ELT(out, index, VECTOR_ELT(input, i));
  }

  int automatic = stp_automatic(VECTOR_ELT(out, STP_WINDOW));
  if (automatic && !supplied[STP_OVERLAP])
    SET_VECTOR_ELT(out, STP_OVERLAP, R_NilValue);
  if (!automatic)
    stp_number(out, STP_WINDOW, DBL_MIN, R_PosInf, 0, 1, 1);
  if (!automatic || VECTOR_ELT(out, STP_OVERLAP) != R_NilValue)
    stp_number(out, STP_OVERLAP, 0, R_PosInf, 0, 1, 0);
  if (!automatic)
  {
    SEXP window = VECTOR_ELT(out, STP_WINDOW),
         overlap = VECTOR_ELT(out, STP_OVERLAP);
    R_xlen_t nw = XLENGTH(window), no = XLENGTH(overlap);
    if (nw > 1 && no > 1 && nw != no)
      st_phase_error("Window/overlap lengths do not match");
    for (R_xlen_t i = 0; i < (nw > no ? nw : no); i++)
    {
      double w = stp_number_at(window, i % nw),
             o = stp_number_at(overlap, i % no);
      if (R_FINITE(w) && w < 1.1 * o)
        st_phase_error("window must be at least 1.1 times overlap");
    }
  }
  stp_number(out, STP_BURNIN, 1, INT_MAX, 1, 0, 0);
  stp_number(out, STP_ITERATIONS, 1, INT_MAX, 1, 0, 0);
  stp_number(out, STP_STATES, 2, 1000000, 1, 0, 0);
  stp_number(out, STP_MARKERS, 4, INT_MAX, 1, 0, 0);
  stp_number(out, STP_NE, DBL_MIN, 3.4e38, 0, 0, 0);
  stp_number(out, STP_SCALE, DBL_MIN, 3.4e38, 0, 0, 0);
  stp_number(out, STP_BUFFER, DBL_MIN, 3.4e38, 0, 0, 0);
  stp_number(out, STP_LR, 1, 3.4e38, 0, 0, 0);
  stp_number(out, STP_RARE, DBL_MIN, .5, 0, 0, 0);
  if (VECTOR_ELT(out, STP_ERR) != R_NilValue)
    stp_number(out, STP_ERR, 0, 1, 0, 0, 0);
  stp_logical(VECTOR_ELT(out, STP_EM), "em");
  if (VECTOR_ELT(out, STP_SEED) != R_NilValue)
    stp_number(out, STP_SEED, -INT_MAX, INT_MAX, 1, 0, 0);
  UNPROTECT(1);
  return out;
}

static void stp_names(SEXP x, const char *role, const char *field)
{
  if (TYPEOF(x) != STRSXP || Rf_any_duplicated(x, FALSE))
    st_phase_error("Invalid or duplicated %s %s", role, field);
  for (R_xlen_t i = 0; i < XLENGTH(x); i++)
    if (STRING_ELT(x, i) == NA_STRING || !CHAR(STRING_ELT(x, i))[0])
      st_phase_error("Invalid or duplicated %s %s", role, field);
}

static SEXP stp_allele_labels(SEXP ref, SEXP alt, const char *role)
{
  const char *r = CHAR(ref), *a = CHAR(alt);
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
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, n));
  SET_STRING_ELT(labels, 0, ref);
  const char *start = a;
  for (int j = 1; j < n; j++)
  {
    const char *end = strchr(start, ',');
    size_t length = end ? (size_t)(end - start) : strlen(start);
    if (!length || length > INT_MAX || (length == 1 && *start == '.'))
      st_phase_error("Invalid %s allele label", role);
    SET_STRING_ELT(labels, j,
                   Rf_mkCharLenCE(start, (int)length, Rf_getCharCE(alt)));
    start += length + (end != NULL);
  }
  if (Rf_any_duplicated(labels, FALSE))
    st_phase_error("Invalid or duplicated %s alleles", role);
  UNPROTECT(1);
  return labels;
}

static SEXP stp_validate(SEXP x, const char *role)
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
      (TYPEOF(pos) != INTSXP && TYPEOF(pos) != REALSXP) || XLENGTH(chr) != m ||
      XLENGTH(pos) != m || TYPEOF(a) != RAWSXP || TYPEOF(b) != RAWSXP ||
      TYPEOF(da) != INTSXP || TYPEOF(db) != INTSXP || XLENGTH(da) != 2 ||
      XLENGTH(db) != 2 || INTEGER(da)[0] != m || INTEGER(da)[1] != n ||
      INTEGER(db)[0] != m || INTEGER(db)[1] != n ||
      XLENGTH(a) != (R_xlen_t)m * n || XLENGTH(b) != XLENGTH(a))
    st_phase_error("Invalid %s dimensions, map or allele matrices", role);
  stp_names(id, role, "ID");
  stp_names(ind, role, "individual");
  stp_names(levels, role, "CHROM.levels");
  for (int j = 0; j < m; j++)
    if (INTEGER(chr)[j] < 1 || INTEGER(chr)[j] > XLENGTH(levels) ||
        !R_FINITE(stp_number_at(pos, j)) || stp_number_at(pos, j) < 0)
      st_phase_error("Invalid %s dimensions, map or allele matrices", role);
  SEXP ref = stp_field(x, "REF"), alt = stp_field(x, "ALT");
  if (TYPEOF(ref) != STRSXP || TYPEOF(alt) != STRSXP || XLENGTH(ref) != m ||
      XLENGTH(alt) != m)
    st_phase_error("Invalid %s allele definitions", role);
  SEXP definitions = PROTECT(Rf_allocVector(VECSXP, m));
  for (int j = 0; j < m; j++)
  {
    if (STRING_ELT(ref, j) == NA_STRING || STRING_ELT(alt, j) == NA_STRING)
      st_phase_error("Invalid %s allele definitions", role);
    SET_VECTOR_ELT(
        definitions, j,
        stp_allele_labels(STRING_ELT(ref, j), STRING_ELT(alt, j), role));
    int count = Rf_length(VECTOR_ELT(definitions, j));
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
  UNPROTECT(1);
  return definitions;
}

static SEXP stp_align(SEXP target, SEXP definitions, SEXP reference,
                      int *ignored)
{
  SEXP ref = reference;
  int nprotect = 0;
  if (TYPEOF(reference) == STRSXP && XLENGTH(reference) == 1 &&
      STRING_ELT(reference, 0) != NA_STRING)
  {
    if (!st_phase_step1_find_dataset(CHAR(STRING_ELT(reference, 0))))
      st_phase_error("Unknown reference data set: %s",
                     CHAR(STRING_ELT(reference, 0)));
    ref = PROTECT(st_markerdata_to_STvcf(reference));
    nprotect++;
  }
  SEXP refdefs = PROTECT(stp_validate(ref, "reference"));
  SEXP ix = PROTECT(Rf_match(stp_field(ref, "ID"), stp_field(target, "ID"), 0));
  SEXP a = stp_field(target, "allele1"), b = stp_field(target, "allele2");
  SEXP ra = stp_field(ref, "allele1"), rb = stp_field(ref, "allele2");
  int m = Rf_length(stp_field(target, "ID")),
      rm = Rf_length(stp_field(ref, "ID"));
  int n = Rf_length(stp_field(target, "individual")),
      rn = Rf_length(stp_field(ref, "individual"));
  if (n > 5000 || rn > 5000 - n)
    st_phase_error("Phasing supports at most 5000 target and reference "
                   "samples combined");
  *ignored = rm - m;
  SEXP input = PROTECT(Rf_duplicate(target));
  SEXP ca = PROTECT(Rf_allocMatrix(RAWSXP, m, n + rn));
  SEXP cb = PROTECT(Rf_allocMatrix(RAWSXP, m, n + rn));
  memcpy(RAW(ca), RAW(a), (size_t)m * n);
  memcpy(RAW(cb), RAW(b), (size_t)m * n);
  SEXP rphase = Rf_getAttrib(ref, Rf_install("phased"));
  int phased = TYPEOF(rphase) == LGLSXP && XLENGTH(rphase) == 1 &&
               LOGICAL(rphase)[0] == TRUE;
  SEXP chr = stp_field(target, "CHROM"), rchr = stp_field(ref, "CHROM");
  SEXP levels = stp_field(target, "CHROM.levels"),
       rlevels = stp_field(ref, "CHROM.levels");
  SEXP pos = stp_field(target, "POS"), rpos = stp_field(ref, "POS");
  SEXP alt = stp_field(target, "ALT");
  for (int j = 0; j < m; j++)
  {
    int r = INTEGER(ix)[j] - 1;
    const char *id = CHAR(STRING_ELT(stp_field(target, "ID"), j));
    if (r < 0)
      st_phase_error("Reference panel lacks target marker: %s", id);
    if (strcmp(CHAR(STRING_ELT(levels, INTEGER(chr)[j] - 1)),
               CHAR(STRING_ELT(rlevels, INTEGER(rchr)[r] - 1))))
      st_phase_error("Reference chromosome labels disagree with target markers");
    if (fabs(stp_number_at(pos, j) - stp_number_at(rpos, r)) > 1e-6)
      st_phase_error("Reference genetic positions disagree with target markers");
    SEXP td = VECTOR_ELT(definitions, j), rd = VECTOR_ELT(refdefs, r);
    SEXP labels[100];
    int codes[100], count = Rf_length(td), shared = 0;
    for (int k = 0; k < count; k++)
      labels[k] = STRING_ELT(td, k);
    for (int k = 0; k < Rf_length(rd); k++)
    {
      int code = -1;
      for (int l = 0; l < count; l++)
        if (!strcmp(CHAR(labels[l]), CHAR(STRING_ELT(rd, k))))
          code = l;
      if (code >= 0)
        shared = 1;
      else
      {
        if (count == 100)
          st_phase_error("Too many combined alleles at %s", id);
        code = count++;
        labels[code] = STRING_ELT(rd, k);
      }
      codes[k] = code;
    }
    if (!shared)
      st_phase_error("No shared allele definition at reference marker: %s", id);
    size_t length = 1;
    for (int k = 1; k < count; k++)
    {
      size_t add = strlen(CHAR(labels[k])) + 1;
      if (add > SIZE_MAX - length)
        st_phase_error("Reference allele label size overflow");
      length += add;
    }
    char *joined = (char *)R_alloc(length + 1, 1);
    char *end = joined;
    for (int k = 1; k < count; k++)
    {
      if (k > 1)
        *end++ = ',';
      size_t size = strlen(CHAR(labels[k]));
      memcpy(end, CHAR(labels[k]), size);
      end += size;
    }
    *end = '\0';
    SET_STRING_ELT(alt, j, Rf_mkChar(count == 1 ? "." : joined));
    for (int i = 0; i < rn; i++)
    {
      R_xlen_t from = r + (R_xlen_t)rm * i, to = j + (R_xlen_t)m * (n + i);
      int x = RAW(ra)[from], y = RAW(rb)[from];
      if (x == 255 || y == 255)
        st_phase_error("Reference calls must be complete at every target marker");
      if (x != y && !phased)
        st_phase_error("Heterozygous reference genotypes must be declared phased");
      RAW(ca)[to] = (Rbyte)codes[x];
      RAW(cb)[to] = (Rbyte)codes[y];
    }
    if (!(j % 1024))
      R_CheckUserInterrupt();
  }
  SET_VECTOR_ELT(input, stp_field_index(input, "ALT"), alt);
  SET_VECTOR_ELT(input, stp_field_index(input, "allele1"), ca);
  SET_VECTOR_ELT(input, stp_field_index(input, "allele2"), cb);
  SEXP individuals = PROTECT(Rf_allocVector(STRSXP, n + rn));
  for (int i = 0; i < n; i++)
    SET_STRING_ELT(individuals, i,
                   STRING_ELT(stp_field(target, "individual"), i));
  for (int i = 0; i < rn; i++)
  {
    char name[32];
    snprintf(name, sizeof(name), "Reference%d", i + 1);
    SET_STRING_ELT(individuals, n + i, Rf_mkChar(name));
  }
  SET_VECTOR_ELT(input, stp_field_index(input, "individual"), individuals);
  UNPROTECT(nprotect + 6);
  return input;
}

static int stp_double_compare(const void *left, const void *right)
{
  double x = *(const double *)left, y = *(const double *)right;
  return (x > y) - (x < y);
}

static double stp_expanded(SEXP value, int index, int chromosomes)
{
  if (XLENGTH(value) != 1 && XLENGTH(value) != chromosomes)
    st_phase_error("Window/overlap must be scalar or have one value per chromosome");
  return stp_number_at(value, XLENGTH(value) == 1 ? 0 : index);
}

static SEXP stp_windows(SEXP input, SEXP options, int chromosomes)
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
  const char *explicit_names[] = {"chromosome", "window", "overlap", "reason"};
  int automatic = stp_automatic(VECTOR_ELT(options, STP_WINDOW));
  int fields = automatic ? 9 : 4, wi = automatic ? 6 : 1, oi = wi + 1,
      ri = wi + 2;
  SEXP out = PROTECT(stp_list(fields, automatic ? auto_names : explicit_names));
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
  SEXP chr = stp_field(input, "CHROM"), pos = stp_field(input, "POS");
  SEXP a = stp_field(input, "allele1"), b = stp_field(input, "allele2");
  int m = Rf_length(chr), n = INTEGER(Rf_getAttrib(a, R_DimSymbol))[1];
  double *common = automatic ? (double *)R_alloc(m, sizeof(double)) : NULL;
  double *gaps = automatic ? (double *)R_alloc(m, sizeof(double)) : NULL;
  float carrier_product =
      (float)n * (float)stp_number_at(VECTOR_ELT(options, STP_RARE), 0);
  int limit = (int)floor(carrier_product);
  if (limit < 3)
    limit = 3;
  int start = 0;
  for (int c = 0; c < chromosomes; c++)
  {
    int end = start + 1;
    while (end < m && INTEGER(chr)[end] == INTEGER(chr)[start])
      end++;
    INTEGER(VECTOR_ELT(out, 0))[c] = INTEGER(chr)[start];
    double window, overlap;
    const char *reason = "explicit";
    char reason_buffer[128];
    if (automatic)
    {
      int count = 0;
      double span = stp_number_at(pos, end - 1) - stp_number_at(pos, start);
      if (end - start < 2 || span <= 0)
        st_phase_error("Window selection requires an ordered, nonzero genetic map on "
                       "chromosome %d",
                       INTEGER(chr)[start]);
      for (int j = start; j < end; j++)
      {
        if (j > start && stp_number_at(pos, j) < stp_number_at(pos, j - 1))
          st_phase_error("Window selection requires an ordered genetic map");
        int carriers[100] = {0}, common_alleles = 0;
        for (int i = 0; i < n; i++)
        {
          R_xlen_t q = j + (R_xlen_t)m * i;
          int x = RAW(a)[q], y = RAW(b)[q];
          if (x != 255)
            carriers[x]++;
          if (y != 255 && y != x)
            carriers[y]++;
        }
        for (int k = 0; k < 100; k++)
          common_alleles += carriers[k] > limit;
        if (common_alleles >= 2)
          common[count++] = stp_number_at(pos, j);
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
          overlap_span = fmax(overlap_span, common[j] - common[j - 40]);
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
        median = ngaps % 2 ? gaps[ngaps / 2]
                           : (gaps[ngaps / 2 - 1] + gaps[ngaps / 2]) / 2;
      }
      INTEGER(VECTOR_ELT(out, 1))[c] = end - start;
      INTEGER(VECTOR_ELT(out, 2))[c] = count;
      REAL(VECTOR_ELT(out, 3))[c] = span;
      REAL(VECTOR_ELT(out, 4))[c] = median;
      INTEGER(VECTOR_ELT(out, 5))[c] = 400;
      if (VECTOR_ELT(options, STP_OVERLAP) != R_NilValue)
      {
        overlap =
            stp_expanded(VECTOR_ELT(options, STP_OVERLAP), c, chromosomes);
        window = fmax(window, 1.1 * overlap);
        snprintf(reason_buffer, sizeof(reason_buffer), "%s;explicit_overlap",
                 reason);
        reason = reason_buffer;
      }
      if (!R_FINITE(window) &&
          end - start > stp_number_at(VECTOR_ELT(options, STP_MARKERS), 0))
        st_phase_error("Automatic whole-chromosome choice exceeds window.markers; "
                       "supply explicit windows");
    }
    else
    {
      window = stp_expanded(VECTOR_ELT(options, STP_WINDOW), c, chromosomes);
      overlap = stp_expanded(VECTOR_ELT(options, STP_OVERLAP), c, chromosomes);
    }
    REAL(VECTOR_ELT(out, wi))[c] = window;
    REAL(VECTOR_ELT(out, oi))[c] = overlap;
    SET_STRING_ELT(VECTOR_ELT(out, ri), c, Rf_mkChar(reason));
    start = end;
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
  SEXP dataset = c->dataset, reference = c->reference, requested = c->requested;
  SEXP impute_s = c->impute, commit_s = c->commit;
  const char *name = NULL;
  st_busy = 1;
  if (st_phase_step1_dataset_name(dataset, &name))
    st_phase_dataset = name;
  int impute = stp_logical(impute_s, "impute");
  int commit = stp_logical(commit_s, "commit");
  if (!st_phase_step1_dataset_name(dataset, &name) ||
      !st_phase_step1_find_dataset(name))
    st_phase_error("Unknown source data set");
  SEXP original = PROTECT(st_markerdata_to_STvcf(dataset));
  SEXP definitions = PROTECT(stp_validate(original, "target"));
  SEXP options = PROTECT(stp_options(requested));
  int with_reference = reference != R_NilValue;
  if (commit && TYPEOF(reference) == STRSXP && XLENGTH(reference) == 1 &&
      STRING_ELT(reference, 0) != NA_STRING &&
      !strcmp(name, CHAR(STRING_ELT(reference, 0))))
    st_phase_error(
        "Reference and target must be different data sets when commit=TRUE");

  SEXP output = PROTECT(Rf_duplicate(original));
  int ignored = 0;
  SEXP input = PROTECT(with_reference
                           ? stp_align(output, definitions, reference, &ignored)
                           : output);
  SEXP chr = stp_field(input, "CHROM"), a = stp_field(input, "allele1"),
       b = stp_field(input, "allele2");
  int m = Rf_length(chr), n = INTEGER(Rf_getAttrib(a, R_DimSymbol))[1];
  int targets = Rf_length(stp_field(original, "individual"));
  if (m < 2 || n < 2 || n > 5000)
    st_phase_error("Phasing requires at least 2 markers and 2 to 5000 target "
                   "and reference samples combined");
  int chromosomes = 0;
  for (int j = 0; j < m; j++)
  {
    if (j && INTEGER(chr)[j] < INTEGER(chr)[j - 1])
      st_phase_error("Phasing requires markers ordered by chromosome");
    if (!j || INTEGER(chr)[j] != INTEGER(chr)[j - 1])
      chromosomes++;
  }
  int threads = 1;
#ifdef MT
  threads = omp_get_max_threads();
#endif
  if (threads < 1)
    st_phase_error("Invalid SelectionTools thread setting; use st.set.num.threads()");
  SEXP windows = PROTECT(stp_windows(input, options, chromosomes));
  SEXP seed = VECTOR_ELT(options, STP_SEED);
  if (seed == R_NilValue)
  {
    GetRNGstate();
    int drawn = 1 + (int)R_unif_index((double)INT_MAX);
    PutRNGstate();
    seed = PROTECT(Rf_ScalarInteger(drawn));
  }
  else
    PROTECT(seed);
  SEXP seed_int = PROTECT(Rf_coerceVector(seed, INTSXP));
  SEXP pos = PROTECT(Rf_coerceVector(stp_field(input, "POS"), REALSXP));
  SEXP bp = PROTECT(Rf_allocVector(INTSXP, m));
  SEXP na = PROTECT(Rf_allocVector(INTSXP, m));
  for (int j = 0; j < m; j++)
  {
    double coordinate = floor(REAL(pos)[j] * 1000000 + .5) + 1;
    if (j && INTEGER(chr)[j] == INTEGER(chr)[j - 1])
      coordinate = fmax(coordinate, (double)INTEGER(bp)[j - 1] + 1);
    if (!R_FINITE(coordinate) || coordinate < 1 || coordinate > INT_MAX)
      st_phase_error("Map cannot be represented by integer positions");
    INTEGER(bp)[j] = (int)coordinate;
    const char *alt = CHAR(STRING_ELT(stp_field(input, "ALT"), j));
    int count = strcmp(alt, ".") ? 2 : 1;
    for (const char *p = alt; *p; p++)
      count += *p == ',';
    INTEGER(na)[j] = count;
  }
  SEXP numeric = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)14 * chromosomes));
  int automatic = stp_automatic(VECTOR_ELT(options, STP_WINDOW));
  const int indices[] = {STP_BURNIN, STP_ITERATIONS, STP_STATES, STP_RARE,
                         STP_SCALE,  STP_NE,         STP_BUFFER, STP_LR};
  for (int c = 0; c < chromosomes; c++)
  {
    double *o = REAL(numeric) + (R_xlen_t)14 * c;
    o[0] = REAL(VECTOR_ELT(windows, automatic ? 6 : 1))[c];
    o[1] = REAL(VECTOR_ELT(windows, automatic ? 7 : 2))[c];
    for (int j = 0; j < 8; j++)
      o[j + 2] = stp_number_at(VECTOR_ELT(options, indices[j]), 0);
    o[10] = LOGICAL(VECTOR_ELT(options, STP_EM))[0];
    o[11] = VECTOR_ELT(options, STP_ERR) == R_NilValue
                ? R_NaN
                : stp_number_at(VECTOR_ELT(options, STP_ERR), 0);
    o[12] = stp_number_at(VECTOR_ELT(options, STP_MARKERS), 0);
    o[13] = threads;
  }
  SEXP phase = PROTECT(
      st_hd_phase_impl(chr, pos, bp, a, b, na, numeric, seed_int, targets));
  SEXP pa = VECTOR_ELT(phase, 0), pb = VECTOR_ELT(phase, 1);
  if (with_reference)
  {
    size_t offset = (size_t)m * targets, size = (size_t)m * (n - targets);
    if (memcmp(RAW(pa) + offset, RAW(a) + offset, size) ||
        memcmp(RAW(pb) + offset, RAW(b) + offset, size))
      st_phase_error("Phasing changed a reference haplotype");
  }
  SEXP oa = PROTECT(Rf_allocMatrix(RAWSXP, m, targets));
  SEXP ob = PROTECT(Rf_allocMatrix(RAWSXP, m, targets));
  R_xlen_t length = (R_xlen_t)m * targets, filled = 0;
  for (R_xlen_t q = 0; q < length; q++)
  {
    int x = RAW(a)[q], y = RAW(b)[q], p = RAW(pa)[q], r = RAW(pb)[q];
    int missing = x == 255 || y == 255;
    if (!missing && !((p == x && r == y) || (p == y && r == x)))
      st_phase_error("Phasing attempted to change an observed genotype");
    if (missing && !impute)
    {
      p = x;
      r = y;
    }
    filled += missing && p != 255 && r != 255;
    RAW(oa)[q] = (Rbyte)p;
    RAW(ob)[q] = (Rbyte)r;
  }
  Rf_setAttrib(oa, R_DimNamesSymbol,
               Rf_getAttrib(stp_field(output, "allele1"), R_DimNamesSymbol));
  Rf_setAttrib(ob, R_DimNamesSymbol,
               Rf_getAttrib(stp_field(output, "allele2"), R_DimNamesSymbol));
  SET_VECTOR_ELT(output, stp_field_index(output, "allele1"), oa);
  SET_VECTOR_ELT(output, stp_field_index(output, "allele2"), ob);
  SEXP phased = PROTECT(Rf_ScalarLogical(TRUE));
  Rf_setAttrib(output, Rf_install("phased"), phased);
  const char *setting_names[] = {
      "method",      "impute",          "phasing.options", "nthreads",
      "actual.seed", "window.settings", "automatic.window", "ref"};
  SEXP settings = PROTECT(stp_list(8, setting_names));
  SET_VECTOR_ELT(settings, 0,
                 Rf_mkString(with_reference ? "reference" : "cohort"));
  SET_VECTOR_ELT(settings, 1, impute_s);
  SET_VECTOR_ELT(settings, 2, options);
  SET_VECTOR_ELT(settings, 3, Rf_ScalarInteger(threads));
  SET_VECTOR_ELT(settings, 4, seed);
  SET_VECTOR_ELT(settings, 5, windows);
  SET_VECTOR_ELT(settings, 6, Rf_ScalarLogical(automatic));
  SET_VECTOR_ELT(settings, 7, reference);
  const char *diag_names[] = {
      "windows",          "maximum.worker.team", "chromosome.windows",
      "filled.genotypes", "reference",           "committed"};
  const char *cohort_names[] = {"windows", "maximum.worker.team",
                                "chromosome.windows", "filled.genotypes",
                                "committed"};
  SEXP diagnostics = PROTECT(stp_list(
      with_reference ? 6 : 5, with_reference ? diag_names : cohort_names));
  SET_VECTOR_ELT(diagnostics, 0, VECTOR_ELT(phase, 2));
  SET_VECTOR_ELT(diagnostics, 1, VECTOR_ELT(phase, 3));
  SET_VECTOR_ELT(diagnostics, 2, VECTOR_ELT(phase, 4));
  SEXP chromosome_names =
      PROTECT(Rf_coerceVector(VECTOR_ELT(windows, 0), STRSXP));
  Rf_setAttrib(VECTOR_ELT(diagnostics, 2), R_NamesSymbol, chromosome_names);
  SET_VECTOR_ELT(diagnostics, 3,
                 filled <= INT_MAX ? Rf_ScalarInteger((int)filled)
                                   : Rf_ScalarReal((double)filled));
  if (with_reference)
  {
    const char *ref_names[] = {"samples", "haplotypes", "unchanged",
                               "ignored.markers"};
    SEXP refdiag = PROTECT(stp_list(4, ref_names));
    SET_VECTOR_ELT(refdiag, 0, Rf_ScalarInteger(n - targets));
    SET_VECTOR_ELT(refdiag, 1, Rf_ScalarInteger(2 * (n - targets)));
    SET_VECTOR_ELT(refdiag, 2, Rf_ScalarLogical(TRUE));
    SET_VECTOR_ELT(refdiag, 3, Rf_ScalarInteger(ignored));
    SET_VECTOR_ELT(diagnostics, 4, refdiag);
    UNPROTECT(1);
  }
  SET_VECTOR_ELT(diagnostics, with_reference ? 5 : 4, commit_s);
  const char *result_names[] = {"data", "settings", "diagnostics"};
  SEXP result = PROTECT(stp_list(3, result_names));
  SET_VECTOR_ELT(result, 0, output);
  SET_VECTOR_ELT(result, 1, settings);
  SET_VECTOR_ELT(result, 2, diagnostics);
  if (commit)
    st_phase_apply(original, output, dataset);
  UNPROTECT(20);
  return result;
}

SEXP st_phase(SEXP dataset, SEXP reference, SEXP requested,
              SEXP impute, SEXP commit)
{
  if (st_busy)
    Rf_error("st.phase is already running");
  st_phase_call c = {dataset, reference, requested, impute, commit};
  return R_ExecWithCleanup(st_phase_body, &c, st_cleanup, NULL);
}
