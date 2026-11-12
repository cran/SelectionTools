#ifndef SELECTIONTOOLS_ALLELE_H
#define SELECTIONTOOLS_ALLELE_H
#include <stddef.h>
#include <stdint.h>
#include <limits.h>
#include <stdatomic.h>
#include <string.h>
#define ST_ALLELE_AUTO_PACK 1
#ifndef GS_ALLELE_TYPE
#define GS_ALLELE_TYPE uint8_t
#define GS_ALLELE_MISSING ((GS_ALLELE_TYPE) - 1)
#define GS_ALLELE_MAX (GS_ALLELE_MISSING - 1)
#endif

#ifndef gs_MAXSTR
#define gs_MAXSTR 4096
#endif
#ifndef gs_MAXNME
#define gs_MAXNME 256
#endif

typedef struct st_io_stream st_io_stream;
void st_io_begin(void);
void st_io_end(void);
const char *st_io_memory(const char *name);
const char *st_io_text(const char *name);
st_io_stream *st_io_open(const char *, const char *);
int st_io_close(st_io_stream *);
int st_io_printf(st_io_stream *, const char *, ...);
size_t st_io_write(const void *, size_t, size_t, st_io_stream *);
int st_io_getc(st_io_stream *);
char *st_io_gets(char *, int, st_io_stream *);
int st_io_eof(st_io_stream *);
int st_io_error(st_io_stream *);
int st_io_scan(st_io_stream *, const char *, ...);
int st_io_remove(const char *);
void st_legacy_C_call(const char *, void **, int);

void *st_work_alloc(size_t count, size_t size);
void *st_work_mark(void);
void st_work_reset(void *mark);
void st_work_adopt(void *pointer, void (*destroy)(void *));
void st_work_detach(void *pointer);
void st_work_native_interface(void (*body)(void *), void *payload);
struct SEXPREC *st_work_interface(struct SEXPREC *(*body)(void *),
                                  void *payload);

int st_identifier_valid(const char *name);

typedef enum
{
    ST_DATA_NULL,
    ST_DATA_REAL,
    ST_DATA_INTEGER,
    ST_DATA_LOGICAL,
    ST_DATA_RAW,
    ST_DATA_STRING,
    ST_DATA_LIST
} st_data_kind;
typedef struct st_data_value st_data_value;
typedef struct st_data_attribute
{
    char *name;
    st_data_value *value;
    struct st_data_attribute *next;
} st_data_attribute;
struct st_data_value
{
    st_data_kind kind;
    size_t length, references;
    union
    {
        double *real;
        int *integer;
        unsigned char *raw;
        char **string;
        st_data_value **list;
    } data;
    st_data_attribute *attributes;
};
st_data_value *st_data_from_R(struct SEXPREC *input);
st_data_value *st_data_from_native(const st_data_value *input);
struct SEXPREC *st_data_to_R(const st_data_value *input);
void st_data_retain(st_data_value *value);
void st_data_release(st_data_value *value);
size_t st_data_bytes(const st_data_value *value);
const st_data_value *st_data_get(const st_data_value *value, const char *name);
const st_data_value *st_data_attribute_get(const st_data_value *value,
                                           const char *name);
const char *st_data_string(const st_data_value *value, size_t index);
st_data_value *st_data_with_real_column(const st_data_value *, const char *,
                                        const double *, size_t);

typedef unsigned long int gs_UNSIGNED;
#define gs_UNSIGNED_MAX ULONG_MAX
typedef double gs_FLOAT;

typedef struct gs_mappoint_type
{
    int chrom;
    double pos;
    const char *name;
    const char *class;
    gs_UNSIGNED nr;
} gs_mappoint_type;

typedef struct gs_cross_type
{
    gs_UNSIGNED p1;
    gs_UNSIGNED p2;
    gs_FLOAT gd;
    gs_FLOAT mu;
    gs_FLOAT mi;
    gs_FLOAT ma;
    gs_FLOAT va;
    gs_FLOAT es;
} gs_cross_type;

typedef struct gs_chrom_type
{
    gs_UNSIGNED NLocC;
    gs_mappoint_type **first_loc;
    gs_UNSIGNED first_loc_i;
    gs_FLOAT *rf;
    gs_FLOAT *q;
    gs_FLOAT *ld;
} gs_chrom_type;

typedef struct gs_varset_type
{
    char name[gs_MAXSTR];
    char msg[3 * gs_MAXSTR];
    int info_level;
    int MaxIL;
    int MaxML;
    gs_UNSIGNED NoInd, NoMar;
    int phased;
    char **IndName;
    char **Pedigree;
    char (*MarName)[gs_MAXNME];
    GS_ALLELE_TYPE *allele_raw1;
    GS_ALLELE_TYPE *allele_raw2;
    unsigned char *allele_packed;
    unsigned char *allele_dictionary;
    size_t allele_packed_bytes;
    unsigned int allele_bits;
    gs_FLOAT *genv;
    gs_FLOAT *genv1;
    gs_FLOAT *genv2;
    gs_UNSIGNED *NoAll;
    GS_ALLELE_TYPE *All_S;
    GS_ALLELE_TYPE **All;
    gs_UNSIGNED *AllCnt_S;
    gs_UNSIGNED **AllCnt;
    double *Pic;
    double *Mis;
    double *IndMis;
    double *Z;
    double *V;
    double *y;
    gs_UNSIGNED NoEff;
    char (*EffNme)[gs_MAXNME];
    int *EffAll;
    int *CmpAll;
    double aa, aA, AA;
    double an, nn, nA;
    gs_UNSIGNED dimu;
    gs_UNSIGNED dimA;
    double *A;

    struct st_marker_state *prediction_state;
    int prediction_role;
    st_data_value *physical_map;
    st_data_value *genomic_relationship_A;
    uint64_t genomic_relationship_signature;
    char genomic_relationship_coding[64];
    double *u;
    double *b;
    double *As;
    double *bs;
    double *GI;
    double *a;
    double *p;
    double *ladi;
    double *est;
    double *SMsqm;
    double *SMsqe;
    double *SMF;
    double *SMp;
    double *SMvar;
    double *sigma;
    double sigma_e;
    gs_UNSIGNED NLoc;
    gs_mappoint_type *map_S;
    char *map_str_S;
    size_t map_str_N;
    gs_mappoint_type **map;
    gs_UNSIGNED *mdl;
    GS_ALLELE_TYPE *Ref1;
    GS_ALLELE_TYPE *Ref2;
    gs_UNSIGNED NHap;
    gs_mappoint_type *Hap_S;
    char *Hap_str_S;
    size_t Hap_str_N;
    gs_UNSIGNED *h_u;
    gs_UNSIGNED *h_o;
    gs_UNSIGNED *h_n;
    gs_UNSIGNED NCrs;
    unsigned int cross_ready;
    uint64_t cross_markers_signature;
    uint64_t cross_effects_signature;
    uint64_t cross_map_signature;
    gs_cross_type *crs_S;
    gs_cross_type **crs;
    gs_UNSIGNED NChrom;
    gs_chrom_type *chrom;
    double *t0, *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *t9;
    gs_UNSIGNED *u0;
} gs_varset_type;

double *st_make_A_values(gs_varset_type *x, const char *coding);
#include <stdlib.h>

static inline char *st_text_dup(const char *text)
{
    if (!text)
        return NULL;
    size_t n = strlen(text);
    if (n >= INT_MAX || n == SIZE_MAX)
        return NULL;
    char *copy = malloc(n + 1);
    if (copy)
        memcpy(copy, text, n + 1);
    return copy;
}

static inline int st_text_set(char **slot, const char *text)
{
    char *copy = text ? st_text_dup(text) : NULL;
    if (text && !copy)
        return -1;
    free(*slot);
    *slot = copy;
    return 0;
}

static inline void st_text_array_free(char **strings, size_t n)
{
    if (strings)
    {
        for (size_t i = 0; i < n; i++)
            free(strings[i]);
        free(strings);
    }
}

static inline int st_identity_set(gs_varset_type *x, size_t i, const char *name,
                                  const char *pedigree)
{
    if (!x || i >= x->NoInd || !name)
        return -1;
    char *n = st_text_dup(name), *p = pedigree ? st_text_dup(pedigree) : NULL;
    if (!n || (pedigree && !p))
    {
        free(n);
        free(p);
        return -1;
    }
    if (!x->IndName)
        x->IndName = calloc((size_t)x->NoInd, sizeof(*x->IndName));
    if (pedigree && !x->Pedigree)
        x->Pedigree = calloc((size_t)x->NoInd, sizeof(*x->Pedigree));
    if (!x->IndName || (pedigree && !x->Pedigree))
    {
        free(n);
        free(p);
        return -1;
    }
    free(x->IndName[i]);
    x->IndName[i] = n;
    if (x->Pedigree)
    {
        free(x->Pedigree[i]);
        x->Pedigree[i] = p;
    }
    return 0;
}

static inline int st_identity_copy(gs_varset_type *dst, size_t di,
                                   const gs_varset_type *src, size_t si)
{
    if (!src || si >= src->NoInd || !src->IndName)
        return -1;
    return st_identity_set(dst, di, src->IndName[si],
                           src->Pedigree ? src->Pedigree[si] : NULL);
}

void st_marker_prediction_clear(gs_varset_type *);
void st_marker_prediction_copy(gs_varset_type *, const gs_varset_type *);
int st_marker_prediction_import(gs_varset_type *, struct SEXPREC *);
size_t st_marker_prediction_bytes(const gs_varset_type *);
const char *st_marker_prediction_id(const gs_varset_type *);
struct SEXPREC *st_marker_prediction_export(const gs_varset_type *);
gs_varset_type *st_prediction_find_source(const char *, int);

void st_marker_cache_clear(gs_varset_type *);
int st_marker_cache_current(gs_varset_type *);
gs_varset_type *st_marker_find(const char *);
void st_hybrid_require_unreferenced(const gs_varset_type *);
int st_hybrid_name_exists(const char *);

static inline unsigned int st_allele_bits(const gs_varset_type *x)
{
    return x->allele_bits ? x->allele_bits : 16;
}
static inline int st_allele_present(const gs_varset_type *x)
{
    return st_allele_bits(x) == 16 ? x->allele_raw1 && x->allele_raw2
                                   : x->allele_packed && x->allele_dictionary;
}

typedef struct
{
    const GS_ALLELE_TYPE *raw1, *raw2;
    const unsigned char *packed, *dictionary;
    unsigned int bits, offset;
} st_allele_marker_view;
typedef struct
{
    GS_ALLELE_TYPE a, b;
} st_allele_pair;

static inline st_allele_marker_view st_allele_marker(const gs_varset_type *x,
                                                     size_t marker)
{
    st_allele_marker_view v = {0};
    size_t base = marker * (size_t)x->NoInd;
    v.bits = st_allele_bits(x);
    if (v.bits == 16)
    {
        v.raw1 = x->allele_raw1 + base;
        v.raw2 = x->allele_raw2 + base;
    }
    else
    {
        v.packed = x->allele_packed + (v.bits == 4 ? base >> 1 : base);
        v.offset = v.bits == 4 ? (unsigned int)(base & 1) : 0;
        v.dictionary = x->allele_dictionary + marker * (v.bits == 4 ? 3 : 15);
    }
    return v;
}
static inline st_allele_pair st_allele_read_pair(const st_allele_marker_view *v,
                                                 size_t individual)
{
    st_allele_pair pair;
    if (v->bits == 16)
    {
        pair.a = v->raw1[individual];
        pair.b = v->raw2[individual];
    }
    else
    {
        size_t at = individual + v->offset;
        unsigned int value = atomic_load_explicit(
            (const _Atomic unsigned char *)(v->packed +
                                            (v->bits == 4 ? at >> 1 : at)),
            memory_order_relaxed);
        unsigned int mask = v->bits == 4 ? 3 : 15;
        if (v->bits == 4)
            value >>= 4 * (at & 1);
        unsigned int a = value & mask, b = (value >> (v->bits / 2)) & mask;
        pair.a = a == mask ? GS_ALLELE_MISSING : v->dictionary[a];
        pair.b = b == mask ? GS_ALLELE_MISSING : v->dictionary[b];
    }
    return pair;
}

static inline GS_ALLELE_TYPE
st_allele_at(const gs_varset_type *x, size_t marker, size_t individual, int hom)
{
    size_t index = marker * (size_t)x->NoInd + individual;
    unsigned int bits = st_allele_bits(x);
    if (bits == 16)
        return hom ? x->allele_raw2[index] : x->allele_raw1[index];
    size_t at = bits == 4 ? index >> 1 : index;
    unsigned int value = atomic_load_explicit(
        (const _Atomic unsigned char *)(x->allele_packed + at),
        memory_order_relaxed);
    unsigned int mask = bits == 4 ? 3 : 15;
    if (bits == 4)
        value >>= 4 * (index & 1);
    unsigned int code = (value >> (hom * (bits / 2))) & mask;
    return code == mask ? GS_ALLELE_MISSING
                        : x->allele_dictionary[marker * mask + code];
}
static inline GS_ALLELE_TYPE st_allele_get(const gs_varset_type *x,
                                           size_t index, int hom)
{
    if (st_allele_bits(x) == 16)
        return hom ? x->allele_raw2[index] : x->allele_raw1[index];
    return st_allele_at(x, index / x->NoInd, index % x->NoInd, hom);
}
static inline int st_allele_write_known_pair(gs_varset_type *x, size_t index,
                                             GS_ALLELE_TYPE a, GS_ALLELE_TYPE b)
{
    unsigned int bits = st_allele_bits(x);
    if (bits == 16)
    {
        x->allele_raw1[index] = a;
        x->allele_raw2[index] = b;
        return 0;
    }
    unsigned int mask = bits == 4 ? 3 : 15, code[2] = {mask, mask};
    const unsigned char *dict =
        x->allele_dictionary + (index / x->NoInd) * mask;
    GS_ALLELE_TYPE value[2] = {a, b};
    for (int h = 0; h < 2; h++)
        if (value[h] != GS_ALLELE_MISSING)
        {
            unsigned int k = 0;
            while (k < mask && dict[k] != value[h])
                k++;
            if (k == mask)
                return -1;
            code[h] = k;
        }
    unsigned char pair = (unsigned char)(code[0] | code[1] << (bits / 2));
    _Atomic unsigned char *cell =
        (_Atomic unsigned char *)(x->allele_packed +
                                  (bits == 4 ? index >> 1 : index));
    if (bits == 8)
        atomic_store_explicit(cell, pair, memory_order_relaxed);
    else
    {
        unsigned int shift = 4 * (index & 1);
        unsigned char old = atomic_load_explicit(cell, memory_order_relaxed),
                      next;
        do
        {
            next = (unsigned char)((old & ~(15u << shift)) |
                                   ((unsigned int)pair << shift));
        } while (!atomic_compare_exchange_weak_explicit(
            cell, &old, next, memory_order_relaxed, memory_order_relaxed));
    }
    return 0;
}

static inline void st_allele_copy(const gs_varset_type *x, int hom,
                                  size_t start, size_t count,
                                  GS_ALLELE_TYPE *out)
{
    if (st_allele_bits(x) == 16)
        memcpy(out, (hom ? x->allele_raw2 : x->allele_raw1) + start,
               count * sizeof(*out));
    else
        for (size_t i = 0; i < count; i++)
            out[i] = st_allele_get(x, start + i, hom);
}
static inline int st_allele_array_equal(const gs_varset_type *x, int hom,
                                        const GS_ALLELE_TYPE *values,
                                        size_t count)
{
    if (!st_allele_present(x))
        return 0;
    if (st_allele_bits(x) == 16)
        return !memcmp(values, hom ? x->allele_raw2 : x->allele_raw1,
                       count * sizeof(*values));
    for (size_t i = 0; i < count; i++)
        if (st_allele_get(x, i, hom) != values[i])
            return 0;
    return 1;
}
void st_allele_replace(gs_varset_type *, const GS_ALLELE_TYPE *,
                       const GS_ALLELE_TYPE *);
int st_allele_statistics_prepare(gs_varset_type *);
void st_allele_statistics_clear(gs_varset_type *);
void st_allele_statistics_move(gs_varset_type *, gs_varset_type *);

void st_allele_storage_release(gs_varset_type *);
void st_allele_pack_native(gs_varset_type *, int);
int st_allele_pack_try(gs_varset_type *, int);
void st_allele_replace_owned(gs_varset_type *, GS_ALLELE_TYPE *,
                             GS_ALLELE_TYPE *);
void st_allele_make_mutable(gs_varset_type *);
struct SEXPREC *st_allele_repack(struct SEXPREC *, struct SEXPREC *);

#endif
