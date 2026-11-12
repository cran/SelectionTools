

#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <ctype.h>
#include <limits.h>
#include <stddef.h>
#include <errno.h>

#include <stdint.h>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined(__linux__) || defined(__APPLE__) || defined(__FreeBSD__)
#define ST_THREAD_DLOPEN
#include <dlfcn.h>
#endif

#ifdef _OPENMP
#ifndef ST
#define MT
#include <omp.h>
#endif
#endif

#include <R.h>
#include <Rversion.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "SelectionToolsAllele.h"

typedef struct st_work_block
{
    void *data;
    void (*destroy)(void *);
    struct st_work_block *previous;
} st_work_block;

static st_work_block *st_work_top = NULL;
static unsigned int st_work_depth = 0;

void *st_work_mark(void)
{
    return st_work_top;
}

void st_work_reset(void *mark)
{
    while (st_work_top != (st_work_block *)mark)
    {
        st_work_block *block = st_work_top;
        if (!block)
            error("Invalid native workspace mark");
        st_work_top = block->previous;
        if (block->destroy)
            block->destroy(block->data);
        else
            free(block->data);
        free(block);
    }
}

void *st_work_alloc(size_t count, size_t size)
{
    if (!st_work_depth)
        error("Native workspace allocation outside an interface boundary");
    if (size && count > SIZE_MAX / size)
        error("Native workspace allocation overflow");
    size_t bytes = count * size;
    st_work_block *block = malloc(sizeof(*block));
    if (!block)
        error("Cannot allocate native workspace record");
    block->data = malloc(bytes ? bytes : 1);
    if (!block->data)
    {
        free(block);
        error("Cannot allocate native workspace");
    }
    block->destroy = NULL;
    block->previous = st_work_top;
    st_work_top = block;
    return block->data;
}

void st_work_adopt(void *pointer, void (*destroy)(void *))
{
    st_work_block *block = malloc(sizeof(*block));
    if (!block)
    {
        destroy(pointer);
        error("Cannot register native allocation cleanup");
    }
    block->data = pointer;
    block->destroy = destroy;
    block->previous = st_work_top;
    st_work_top = block;
}

void st_work_detach(void *pointer)
{
    st_work_block **link = &st_work_top;
    while (*link)
    {
        st_work_block *block = *link;
        if (block->data == pointer)
        {
            *link = block->previous;
            free(block);
            return;
        }
        link = &block->previous;
    }
    error("Native allocation is not owned by the current workspace");
}

typedef struct
{
    SEXP (*body)(void *);
    void *payload, *mark;
    unsigned int depth;
} st_work_scope;

static SEXP st_work_body(void *payload)
{
    st_work_scope *scope = payload;
    st_work_depth = scope->depth + 1;
    return scope->body(scope->payload);
}

static void st_work_cleanup(void *payload, Rboolean jump)
{
    (void)jump;
    st_work_scope *scope = payload;
    st_work_reset(scope->mark);
    st_work_depth = scope->depth;
}

SEXP st_work_interface(SEXP (*body)(void *), void *payload)
{
    st_work_scope scope = {body, payload, st_work_mark(), st_work_depth};
    SEXP token = PROTECT(R_MakeUnwindCont());
    SEXP result =
        R_UnwindProtect(st_work_body, &scope, st_work_cleanup, &scope, token);
    UNPROTECT(1);
    return result;
}

typedef struct st_io_buffer
{
    char *text;
    size_t used, capacity;
    const char *name;
    struct st_io_buffer *next;
} st_io_buffer;
typedef struct st_io_context
{
    st_io_buffer *buffers;
    struct st_io_context *previous;
} st_io_context;
struct st_io_stream
{
    FILE *file;
    st_io_buffer *buffer;
    size_t position;
    int ended;
};
static st_io_context *st_io_active;
static void st_io_context_destroy(void *p)
{
    st_io_context *context = p;
    st_io_active = context->previous;
    free(context);
}
void st_io_begin(void)
{
    st_io_context *context = calloc(1, sizeof(*context));
    if (!context)
        error("Cannot allocate memory stream context");
    context->previous = st_io_active;
    st_work_adopt(context, st_io_context_destroy);
    st_io_active = context;
}
void st_io_end(void)
{
    if (st_io_active)
        st_io_active = st_io_active->previous;
}
static void st_io_buffer_destroy(void *p)
{
    st_io_buffer *buffer = p;
    free(buffer->text);
    free(buffer);
}
const char *st_io_memory(const char *name)
{
    if (!st_io_active)
        error("Missing memory stream context");
    st_io_buffer *buffer = calloc(1, sizeof(*buffer));
    if (!buffer)
        error("Cannot allocate memory stream");
    st_work_adopt(buffer, st_io_buffer_destroy);
    size_t n = strlen(name) + 18;
    char *key = st_work_alloc(n, 1);
    snprintf(key, n, "@SelectionTools:%s", name);
    buffer->name = key;
    buffer->next = st_io_active->buffers;
    st_io_active->buffers = buffer;
    return key;
}
static st_io_buffer *st_io_find(const char *name)
{
    for (st_io_context *c = st_io_active; c; c = c->previous)
        for (st_io_buffer *b = c->buffers; b; b = b->next)
            if (!strcmp(name, b->name))
                return b;
    return NULL;
}
const char *st_io_text(const char *name)
{
    st_io_buffer *buffer = st_io_find(name);
    if (!buffer)
        error("Unknown memory stream");
    return buffer->text ? buffer->text : "";
}
static void st_io_stream_destroy(void *p)
{
    st_io_stream *stream = p;
    if (stream->file)
        fclose(stream->file);
    free(stream);
}
st_io_stream *st_io_open(const char *name, const char *mode)
{
    st_io_buffer *buffer = st_io_find(name);
    FILE *file = NULL;
    if (!buffer)
    {
        if (!strncmp(name, "@SelectionTools:", strlen("@SelectionTools:")))
            error("Unbound memory stream");
        file = fopen(name, mode);
        if (!file)
            return NULL;
    }
    st_io_stream *stream = calloc(1, sizeof(*stream));
    if (!stream)
    {
        if (file)
            fclose(file);
        error("Cannot allocate stream handle");
    }
    stream->file = file;
    stream->buffer = buffer;
    st_work_adopt(stream, st_io_stream_destroy);
    if (buffer && mode[0] == 'w')
    {
        buffer->used = 0;
        if (buffer->text)
            buffer->text[0] = '\0';
    }
    if (buffer && mode[0] == 'a')
        stream->position = buffer->used;
    return stream;
}
int st_io_close(st_io_stream *stream)
{
    int status = 0;
    if (stream && stream->file)
    {
        status = fclose(stream->file);
        stream->file = NULL;
    }
    return status;
}
size_t st_io_write(const void *source, size_t size, size_t count,
                   st_io_stream *stream)
{
    if (stream->file)
        return fwrite(source, size, count, stream->file);
    if (size && count > (SIZE_MAX - stream->position - 1) / size)
        error("Memory stream size overflow");
    size_t bytes = size * count, needed = stream->position + bytes + 1;
    st_io_buffer *buffer = stream->buffer;
    if (needed > buffer->capacity)
    {
        size_t capacity = needed > SIZE_MAX / 2 ? needed : needed * 2;
        char *text = realloc(buffer->text, capacity);
        if (!text)
            error("Cannot grow memory stream");
        buffer->text = text;
        buffer->capacity = capacity;
    }
    memcpy(buffer->text + stream->position, source, bytes);
    stream->position += bytes;
    if (buffer->used < stream->position)
        buffer->used = stream->position;
    buffer->text[buffer->used] = '\0';
    return count;
}
int st_io_printf(st_io_stream *stream, const char *format, ...)
{
    va_list args, copy;
    va_start(args, format);
    if (stream->file)
    {
        int result = vfprintf(stream->file, format, args);
        va_end(args);
        return result;
    }
    va_copy(copy, args);
    int n = vsnprintf(NULL, 0, format, copy);
    va_end(copy);
    if (n < 0)
    {
        va_end(args);
        error("Cannot format memory stream");
    }
    size_t bytes = (size_t)n;
    if (bytes > SIZE_MAX - stream->position - 1)
    {
        va_end(args);
        error("Memory stream size overflow");
    }
    size_t needed = stream->position + bytes + 1;
    st_io_buffer *buffer = stream->buffer;
    if (needed > buffer->capacity)
    {
        size_t capacity = needed > SIZE_MAX / 2 ? needed : needed * 2;
        char *text = realloc(buffer->text,capacity);
        if (!text)
        {
            va_end(args);
            error("Cannot grow memory stream");
        }
        buffer->text = text; buffer->capacity = capacity;
    }
    vsnprintf(buffer->text + stream->position,bytes + 1,format,args);
    va_end(args);
    stream->position += bytes;
    if (buffer->used < stream->position) buffer->used = stream->position;
    buffer->text[buffer->used] = '\0';
    return n;
}
int st_io_getc(st_io_stream *stream)
{
    if (stream->file)
        return fgetc(stream->file);
    if (stream->position >= stream->buffer->used)
    {
        stream->ended = 1;
        return EOF;
    }
    return (unsigned char)stream->buffer->text[stream->position++];
}
char *st_io_gets(char *text, int size, st_io_stream *stream)
{
    if (stream->file)
        return fgets(text, size, stream->file);
    int n = 0, c = EOF;
    while (n < size - 1 && (c = st_io_getc(stream)) != EOF)
    {
        text[n++] = c;
        if (c == '\n')
            break;
    }
    if (!n && c == EOF)
        return NULL;
    text[n] = '\0';
    return text;
}
int st_io_eof(st_io_stream *stream)
{
    return stream->file ? feof(stream->file) : stream->ended;
}
int st_io_error(st_io_stream *stream)
{
    return stream->file ? ferror(stream->file) : 0;
}
int st_io_scan(st_io_stream *stream, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    if (stream->file)
    {
        int result = vfscanf(stream->file, format, args);
        va_end(args);
        return result;
    }
    if (strcmp(format, "%4095s"))
    {
        va_end(args);
        error("Unsupported memory stream scan");
    }
    char *out = va_arg(args, char *);
    va_end(args);
    int c;
    do { c = st_io_getc(stream); } while (c != EOF && isspace((unsigned char)c));
    if (c == EOF)
        return EOF;
    size_t n = 0;
    do
    {
        out[n++] = c;
        if (n == 4095)
            break;
        c = st_io_getc(stream);
    } while (c != EOF && !isspace((unsigned char)c));
    out[n] = '\0';
    return 1;
}
int st_io_remove(const char *name)
{
    st_io_buffer *buffer = st_io_find(name);
    if (!buffer)
        return remove(name);
    buffer->used = 0;
    if (buffer->text)
        buffer->text[0] = '\0';
    return 0;
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

#define MAXTHREADS 256

#include "MathLib12.h"

#if MAXTHREADS != ST_RNG_MAXTHREADS
#error "MAXTHREADS and ST_RNG_MAXTHREADS differ"
#endif

#define GS_FREE(ptr) (free((ptr)), (ptr) = NULL)

SEXP st_has_openmp(void)
{
#ifdef MT
    return ScalarLogical(1);
#else
    return ScalarLogical(0);
#endif
}

#ifndef MT
#define omp_get_thread_num() 0
#define omp_set_num_threads(thr)                                               \
    {                                                                          \
    }
#define omp_get_num_threads() 1
#define omp_set_dynamic(thr)                                                   \
    {                                                                          \
    }
#endif

#ifdef MultPrec
#include <gmp.h>
#endif

#define gs_MAXSTR 4096

#define gs_MAXNME 256

#define gs_MAXALL (GS_ALLELE_MAX + 1)

#define gs_MAXLNE 65536

#define ERR_F                                                                  \
    {                                                                          \
        gs_info(gv, -2, "File open error");                                    \
        *retval = -2;                                                          \
        return;                                                                \
    }
#define ERR_M                                                                  \
    {                                                                          \
        gs_info(gv, -2, "Memory allocation error");                            \
        *retval = -2;                                                          \
        return;                                                                \
    }

typedef struct
{
    gs_UNSIGNED a;
    gs_UNSIGNED n;
} gs_usp;

gs_varset_type GV = {.name = "default",
                     .msg = "",
                     .info_level = 0,
                     .MaxIL = 0,
                     .MaxML = 0,
                     .NoInd = 0,
                     .NoMar = 0,
                     .phased = 0,
                     .IndName = NULL,
                     .MarName = NULL,
                     .allele_raw1 = NULL,
                     .allele_raw2 = NULL,
                     .genv = NULL,
                     .genv1 = NULL,
                     .genv2 = NULL,
                     .NoAll = NULL,
                     .All_S = NULL,
                     .All = NULL,
                     .AllCnt_S = NULL,
                     .AllCnt = NULL,
                     .Pic = NULL,
                     .Mis = NULL,
                     .IndMis = NULL,
                     .Z = NULL,
                     .V = NULL,
                     .y = NULL,
                     .NoEff = 0,
                     .EffNme = NULL,
                     .EffAll = NULL,
                     .CmpAll = NULL,
                     .aa = 0,
                     .aA = 1,
                     .AA = 2,
                     .an = 0.5,
                     .nn = 1,
                     .nA = 1.5,
                     .dimu = 0,
                     .dimA = 0,
                     .A = NULL,
                     .u = NULL,
                     .b = NULL,
                     .As = NULL,
                     .bs = NULL,
                     .GI = NULL,
                     .a = NULL,
                     .p = NULL,
                     .ladi = NULL,
                     .est = NULL,
                     .SMsqm = NULL,
                     .SMsqe = NULL,
                     .SMF = NULL,
                     .SMp = NULL,
                     .SMvar = NULL,
                     .sigma = NULL,
                     .sigma_e = 0,
                     .NLoc = 0,
                     .map_S = NULL,
                     .map_str_S = NULL,
                     .map_str_N = 0,
                     .map = NULL,
                     .mdl = NULL,
                     .Ref1 = NULL,
                     .Ref2 = NULL,
                     .NHap = 0,
                     .Hap_S = NULL,
                     .Hap_str_S = NULL,
                     .Hap_str_N = 0,
                     .h_u = NULL,
                     .h_o = NULL,
                     .h_n = NULL,
                     .NCrs = 0,
                     .cross_ready = 0,
                     .cross_markers_signature = 0,
                     .cross_effects_signature = 0,
                     .cross_map_signature = 0,
                     .crs_S = NULL,
                     .crs = NULL,
                     .NChrom = 0,
                     .chrom = NULL,
                     .t0 = NULL,
                     .t1 = NULL,
                     .t2 = NULL,
                     .t3 = NULL,
                     .t4 = NULL,
                     .t5 = NULL,
                     .t6 = NULL,
                     .t7 = NULL,
                     .t8 = NULL,
                     .t9 = NULL,
                     .u0 = NULL};

gs_varset_type *gv = &GV;
gs_varset_type *sh = &GV;

gs_varset_type ES;
gs_varset_type *es = &ES;
int esi = 0;
gs_varset_type VS;
gs_varset_type *vs = &VS;
int vsi = 0;

gs_varset_type **gs_data = NULL;
gs_UNSIGNED gs_ndta = 0;

const char *gs_tok_del = "; \r\n\t\"";
char *nix = "";
char **nixp = &nix;
int nein = 0;
int *neinp = &nein;
double null = 0;
double *nullp = &null;

int Srv = 0;
int *Sretval = &Srv;
char *empty = "";
char **Sout_filename = &empty;
int auf = 0;
int *Sauxfiles = &auf;
int mz = -2;
int *imz = &mz;

typedef struct
{
    void (*body)(void *);
    void *payload;
} st_work_C_scope;

static SEXP st_work_C_body(void *payload)
{
    st_work_C_scope *scope = payload;
    scope->body(scope->payload);
    return R_NilValue;
}

void st_work_native_interface(void (*body)(void *), void *payload)
{
    st_work_C_scope scope = {body, payload};
    st_work_interface(st_work_C_body, &scope);
}

static void *st_data_allocate(size_t count, size_t size)
{
    if (size && count > SIZE_MAX / size)
        error("Native data allocation overflow");
    void *out = calloc(count ? count : 1, size ? size : 1);
    if (!out)
        error("Cannot allocate native data");
    return out;
}

static char *st_data_copy_string(const char *text)
{
    if (!text)
        return NULL;
    size_t n = strlen(text);
    char *out = st_data_allocate(n + 1, 1);
    memcpy(out, text, n + 1);
    return out;
}

void st_data_retain(st_data_value *value)
{
    if (value)
    {
        if (value->references == SIZE_MAX)
            error("Native data reference overflow");
        value->references++;
    }
}

void st_data_release(st_data_value *value)
{
    if (!value || --value->references)
        return;
    st_data_attribute *attribute = value->attributes;
    while (attribute)
    {
        st_data_attribute *next = attribute->next;
        free(attribute->name);
        st_data_release(attribute->value);
        free(attribute);
        attribute = next;
    }
    switch (value->kind)
    {
    case ST_DATA_STRING:
        if (value->data.string)
            for (size_t i = 0; i < value->length; i++)
                free(value->data.string[i]);
        free(value->data.string);
        break;
    case ST_DATA_LIST:
        if (value->data.list)
            for (size_t i = 0; i < value->length; i++)
                st_data_release(value->data.list[i]);
        free(value->data.list);
        break;
    case ST_DATA_REAL:
        free(value->data.real);
        break;
    case ST_DATA_INTEGER:
    case ST_DATA_LOGICAL:
        free(value->data.integer);
        break;
    case ST_DATA_RAW:
        free(value->data.raw);
        break;
    default:
        break;
    }
    free(value);
}

static st_data_value *st_data_new(void)
{
    st_data_value *value = st_data_allocate(1, sizeof(*value));
    value->references = 1;
    return value;
}

static SEXP st_data_attributes(SEXP input)
{
    SEXP quoted = PROTECT(lang2(install("quote"), input)),
         call =
             PROTECT(lang2(findFun(install("attributes"), R_BaseEnv), quoted));
    SEXP result = eval(call, R_BaseEnv);
    UNPROTECT(2);
    return result;
}

static void st_data_read_value(st_data_value *out, SEXP input, int depth)
{
    if (depth > 100)
        error("Native interface data are too deeply nested or cyclic");
    int type = TYPEOF(input);
    if (type == NILSXP)
    {
        out->kind = ST_DATA_NULL;
        return;
    }
    if (type == REALSXP)
        out->kind = ST_DATA_REAL;
    else if (type == INTSXP)
        out->kind = ST_DATA_INTEGER;
    else if (type == LGLSXP)
        out->kind = ST_DATA_LOGICAL;
    else if (type == RAWSXP)
        out->kind = ST_DATA_RAW;
    else if (type == STRSXP)
        out->kind = ST_DATA_STRING;
    else if (type == VECSXP)
        out->kind = ST_DATA_LIST;
    else
        error("Persistent native data cannot contain R objects of type '%s'",
              type2char(type));
    out->length = (size_t)XLENGTH(input);
    switch (out->kind)
    {
    case ST_DATA_REAL:
        out->data.real = st_data_allocate(out->length, sizeof(double));
        memcpy(out->data.real, REAL(input), out->length * sizeof(double));
        break;
    case ST_DATA_INTEGER:
    case ST_DATA_LOGICAL:
        out->data.integer = st_data_allocate(out->length, sizeof(int));
        memcpy(out->data.integer,
               type == INTSXP ? INTEGER(input) : LOGICAL(input),
               out->length * sizeof(int));
        break;
    case ST_DATA_RAW:
        out->data.raw = st_data_allocate(out->length, 1);
        memcpy(out->data.raw, RAW(input), out->length);
        break;
    case ST_DATA_STRING:
        out->data.string = st_data_allocate(out->length, sizeof(char *));
        for (size_t i = 0; i < out->length; i++)
            if (STRING_ELT(input, (R_xlen_t)i) != NA_STRING)
                out->data.string[i] = st_data_copy_string(
                    translateCharUTF8(STRING_ELT(input, (R_xlen_t)i)));
        break;
    case ST_DATA_LIST:
        out->data.list = st_data_allocate(out->length, sizeof(st_data_value *));
        for (size_t i = 0; i < out->length; i++)
        {
            out->data.list[i] = st_data_new();
            st_data_read_value(out->data.list[i],
                               VECTOR_ELT(input, (R_xlen_t)i), depth + 1);
        }
        break;
    default:
        break;
    }
    SEXP attributes = PROTECT(st_data_attributes(input)),
         names = getAttrib(attributes, R_NamesSymbol);
    st_data_attribute **tail = &out->attributes;
    for (R_xlen_t i = 0; i < (isNull(attributes) ? 0 : XLENGTH(attributes));
         i++)
    {
        st_data_attribute *a = st_data_allocate(1, sizeof(*a));
        *tail = a;
        tail = &a->next;
        a->name = st_data_copy_string(CHAR(STRING_ELT(names, i)));
        a->value = st_data_new();
        st_data_read_value(a->value, VECTOR_ELT(attributes, i), depth + 1);
    }
    UNPROTECT(1);
}

typedef struct
{
    SEXP input;
    st_data_value *result;
    int complete;
} st_data_read_context;

static SEXP st_data_read_body(void *payload)
{
    st_data_read_context *c = payload;
    c->result = st_data_new();
    st_data_read_value(c->result, c->input, 0);
    c->complete = 1;
    return R_NilValue;
}

static void st_data_read_cleanup(void *payload, Rboolean jump)
{
    st_data_read_context *c = payload;
    if (jump || !c->complete)
    {
        st_data_release(c->result);
        c->result = NULL;
    }
}

static void st_data_destroy(void *value)
{
    st_data_release((st_data_value *)value);
}

st_data_value *st_data_from_R(SEXP input)
{
    st_data_read_context c = {input, NULL, 0};
    SEXP token = PROTECT(R_MakeUnwindCont());
    R_UnwindProtect(st_data_read_body, &c, st_data_read_cleanup, &c, token);
    UNPROTECT(1);
    st_work_adopt(c.result, st_data_destroy);
    return c.result;
}

static SEXP st_data_write_value(const st_data_value *value, int depth)
{
    if (!value || value->kind == ST_DATA_NULL)
        return R_NilValue;
    if (depth > 100 || value->length > R_XLEN_T_MAX)
        error("Invalid native interface data dimensions");
    int type = value->kind == ST_DATA_REAL      ? REALSXP
               : value->kind == ST_DATA_INTEGER ? INTSXP
               : value->kind == ST_DATA_LOGICAL ? LGLSXP
               : value->kind == ST_DATA_RAW     ? RAWSXP
               : value->kind == ST_DATA_STRING  ? STRSXP
                                                : VECSXP;
    SEXP out = PROTECT(allocVector(type, (R_xlen_t)value->length));
    switch (value->kind)
    {
    case ST_DATA_REAL:
        memcpy(REAL(out), value->data.real, value->length * sizeof(double));
        break;
    case ST_DATA_INTEGER:
        memcpy(INTEGER(out), value->data.integer, value->length * sizeof(int));
        break;
    case ST_DATA_LOGICAL:
        memcpy(LOGICAL(out), value->data.integer, value->length * sizeof(int));
        break;
    case ST_DATA_RAW:
        memcpy(RAW(out), value->data.raw, value->length);
        break;
    case ST_DATA_STRING:
        for (size_t i = 0; i < value->length; i++)
            SET_STRING_ELT(out, (R_xlen_t)i,
                           value->data.string[i]
                               ? mkCharCE(value->data.string[i], CE_UTF8)
                               : NA_STRING);
        break;
    case ST_DATA_LIST:
        for (size_t i = 0; i < value->length; i++)
        {
            SEXP child =
                PROTECT(st_data_write_value(value->data.list[i], depth + 1));
            SET_VECTOR_ELT(out, (R_xlen_t)i, child);
            UNPROTECT(1);
        }
        break;
    default:
        break;
    }
    for (st_data_attribute *a = value->attributes; a; a = a->next)
    {
        SEXP data = PROTECT(st_data_write_value(a->value, depth + 1));
        setAttrib(out, install(a->name), data);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return out;
}

SEXP st_data_to_R(const st_data_value *value)
{
    return st_data_write_value(value, 0);
}

size_t st_data_bytes(const st_data_value *value)
{
    if (!value)
        return 0;
    size_t bytes = sizeof(*value);
    if (value->kind == ST_DATA_LIST)
    {
        bytes += value->length * sizeof(st_data_value *);
        for (size_t i = 0; i < value->length; i++)
            bytes += st_data_bytes(value->data.list[i]);
    }
    else if (value->kind == ST_DATA_STRING)
    {
        bytes += value->length * sizeof(char *);
        for (size_t i = 0; i < value->length; i++)
            if (value->data.string[i])
                bytes += strlen(value->data.string[i]) + 1;
    }
    else
        bytes += value->length * (value->kind == ST_DATA_REAL  ? sizeof(double)
                                  : value->kind == ST_DATA_RAW ? 1
                                                               : sizeof(int));
    for (st_data_attribute *a = value->attributes; a; a = a->next)
        bytes += sizeof(*a) + strlen(a->name) + 1 + st_data_bytes(a->value);
    return bytes;
}

const st_data_value *st_data_attribute_get(const st_data_value *value,
                                           const char *name)
{
    if (value)
    {
        for (st_data_attribute *a = value->attributes; a; a = a->next)
            if (!strcmp(a->name, name))
                return a->value;
    }
    return NULL;
}

const char *st_data_string(const st_data_value *value, size_t index)
{
    return value && value->kind == ST_DATA_STRING && index < value->length
               ? value->data.string[index]
               : NULL;
}

const st_data_value *st_data_get(const st_data_value *value, const char *name)
{
    if (!value || value->kind != ST_DATA_LIST)
        return NULL;
    const st_data_value *names = st_data_attribute_get(value, "names");
    for (size_t i = 0; i < value->length; i++)
    {
        const char *key = st_data_string(names, i);
        if (key && !strcmp(key, name))
            return value->data.list[i];
    }
    return NULL;
}

static void st_data_copy_value(st_data_value *out, const st_data_value *input,
                               int depth)
{
    if (!input)
        return;
    if (depth > 100)
        error("Native interface data are too deeply nested or cyclic");
    out->kind = input->kind;
    out->length = input->length;
    switch (input->kind)
    {
    case ST_DATA_REAL:
        out->data.real = st_data_allocate(out->length, sizeof(double));
        memcpy(out->data.real, input->data.real, out->length * sizeof(double));
        break;
    case ST_DATA_INTEGER:
    case ST_DATA_LOGICAL:
        out->data.integer = st_data_allocate(out->length, sizeof(int));
        memcpy(out->data.integer, input->data.integer,
               out->length * sizeof(int));
        break;
    case ST_DATA_RAW:
        out->data.raw = st_data_allocate(out->length, 1);
        memcpy(out->data.raw, input->data.raw, out->length);
        break;
    case ST_DATA_STRING:
        out->data.string = st_data_allocate(out->length, sizeof(char *));
        for (size_t i = 0; i < out->length; i++)
            if (input->data.string[i])
                out->data.string[i] =
                    st_data_copy_string(input->data.string[i]);
        break;
    case ST_DATA_LIST:
        out->data.list = st_data_allocate(out->length, sizeof(st_data_value *));
        for (size_t i = 0; i < out->length; i++)
        {
            out->data.list[i] = st_data_new();
            st_data_copy_value(out->data.list[i], input->data.list[i],
                               depth + 1);
        }
        break;
    default:
        break;
    }
    st_data_attribute **tail = &out->attributes;
    for (const st_data_attribute *a = input->attributes; a; a = a->next)
    {
        st_data_attribute *copy = st_data_allocate(1, sizeof(*copy));
        *tail = copy;
        tail = &copy->next;
        copy->name = st_data_copy_string(a->name);
        copy->value = st_data_new();
        st_data_copy_value(copy->value, a->value, depth + 1);
    }
}

st_data_value *st_data_from_native(const st_data_value *input)
{
    st_data_value *out = st_data_new();
    st_work_adopt(out, st_data_destroy);
    st_data_copy_value(out, input, 0);
    return out;
}

st_data_value *st_data_with_real_column(const st_data_value *table,
                                        const char *name, const double *values,
                                        size_t n)
{
    if (!table || table->kind != ST_DATA_LIST)
        error("Native performance table is missing");
    const st_data_value *names = st_data_attribute_get(table, "names");
    if (!names || names->kind != ST_DATA_STRING ||
        names->length != table->length)
        error("Native performance table has invalid column names");
    size_t at = table->length;
    for (size_t i = 0; i < table->length; i++)
    {
        if (table->data.list[i]->length != n)
            error("Native performance row counts differ");
        if (!strcmp(names->data.string[i], name))
            at = i;
    }
    st_data_value *out = st_data_new();
    st_work_adopt(out, st_data_destroy);
    out->kind = ST_DATA_LIST;
    out->length = table->length + (at == table->length);
    out->data.list = st_data_allocate(out->length, sizeof(st_data_value *));
    for (size_t i = 0; i < table->length; i++)
        if (i != at)
        {
            out->data.list[i] = table->data.list[i];
            st_data_retain(out->data.list[i]);
        }
    st_data_value *column = st_data_new();
    out->data.list[at] = column;
    column->kind = ST_DATA_REAL;
    column->length = n;
    column->data.real = st_data_allocate(n, sizeof(double));
    memcpy(column->data.real, values, n * sizeof(double));
    st_data_attribute **tail = &out->attributes;
    for (st_data_attribute *a = table->attributes; a; a = a->next)
    {
        st_data_attribute *copy = st_data_allocate(1, sizeof(*copy));
        *tail = copy;
        tail = &copy->next;
        copy->name = st_data_copy_string(a->name);
        if (strcmp(a->name, "names"))
        {
            copy->value = a->value;
            st_data_retain(copy->value);
        }
        else
        {
            copy->value = st_data_new();
            st_data_value *v = copy->value;
            v->kind = ST_DATA_STRING;
            v->length = out->length;
            v->data.string = st_data_allocate(v->length, sizeof(char *));
            for (size_t i = 0; i < table->length; i++)
                v->data.string[i] = st_data_copy_string(names->data.string[i]);
            if (at == table->length)
                v->data.string[at] = st_data_copy_string(name);
        }
    }
    return out;
}

double **st_dmatv(double *a, gs_UNSIGNED nr, gs_UNSIGNED nc)
{
    double **p = malloc(nr * sizeof(double *));

    p[0] = &a[0];
    for (gs_UNSIGNED i = 1; i < nr; i++)
        p[i] = p[i - 1] + nc;

    return (p);
}

void gs_info(gs_varset_type *, int, char *);

double st_fprob(double f, gs_UNSIGNED df1, gs_UNSIGNED df2)
{

    return pf(f, (double)df1, (double)df2, 0, 0);
}

double st_tprob(double t, gs_UNSIGNED df)
{
    return 2.0 * pt(fabs(t), (double)df, 0, 0);
}

void st_fprob_GV(double *f, int *df1, int *df2, double *prob)
{
    *prob = st_fprob(*f, (gs_UNSIGNED)(*df1), (gs_UNSIGNED)(*df2));
}

void st_tprob_GV(double *t, int *df, double *prob)
{
    *prob = st_tprob(*t, (gs_UNSIGNED)(*df));
}

#include <R.h>
#include <Rversion.h>
#include <R_ext/RS.h>
#if defined(_WIN32)

#ifndef BLAS_extern
#define BLAS_extern extern __declspec(dllimport)
#endif
#ifndef La_extern
#define La_extern extern __declspec(dllimport)
#endif
#endif
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "MathLib11.h"

void gs_mm_add(gs_UNSIGNED nr, gs_UNSIGNED nc, double *A, double *B, double *C)
{
    gs_UNSIGNED i, j;

    for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
            C[i * nc + j] = A[i * nc + j] + B[i * nc + j];
}

void gs_sm_add(gs_UNSIGNED nr, gs_UNSIGNED nc, double *A, double *B, double *C)
{
    gs_UNSIGNED i, j;

    for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
            C[i * nc + j] = *A + B[i * nc + j];
}

static int st_threads_initialized = 0;
static int st_threads_requested = 1;
static int st_threads_effective = 1;
static int st_threads_explicit = 0;
static const char *st_blas_backend = "unavailable";
static char st_blas_library[4096] = "";
static int (*st_blas_get_threads)(void) = NULL;
static void (*st_blas_set_threads)(int) = NULL;
static int (*st_blas_get_parallel)(void) = NULL;

#if defined(_WIN32)
typedef HMODULE st_thread_library;
#else
typedef void *st_thread_library;
#endif

static void st_thread_symbol(st_thread_library library, const void *base,
                             const char *name, void *target, size_t size)
{
#if defined(_WIN32)
    FARPROC symbol = GetProcAddress(library, name);
    HMODULE owner = NULL;
    const void *address = NULL;

    (void)base;
    if (sizeof(address) != sizeof(symbol) || sizeof(symbol) != size)
        return;
    memcpy(&address, &symbol, sizeof(address));
    if (symbol != NULL &&
        GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                               GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                           (LPCSTR)address, &owner) &&
        owner == library)
        memcpy(target, &symbol, size);
#elif defined(ST_THREAD_DLOPEN)
    void *symbol = dlsym(library, name);
    Dl_info owner;

    if (symbol != NULL && sizeof(symbol) == size &&
        dladdr(symbol, &owner) != 0 && owner.dli_fbase == base)
        memcpy(target, &symbol, size);
#else
    (void)library;
    (void)base;
    (void)name;
    (void)target;
    (void)size;
#endif
}

static void st_thread_blas_initialize(void)
{
    st_thread_library library = NULL;
    const void *base = NULL;
    const void *address = NULL;
    void (*function)(void) = (void (*)(void))F77_CALL(dgemm);

    if (sizeof(address) != sizeof(function))
        return;
    memcpy(&address, &function, sizeof(address));
#if defined(_WIN32)
    if (!GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                                GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                            (LPCSTR)address, &library))
        return;
    {
        WCHAR path[4096];
        DWORD length = GetModuleFileNameW(library, path, 4096);

        if (length > 0 && length < 4096)
        {
            if (WideCharToMultiByte(CP_UTF8, 0, path, -1, st_blas_library,
                                    sizeof(st_blas_library), NULL, NULL) == 0)
                st_blas_library[0] = '\0';
        }
    }
#elif defined(ST_THREAD_DLOPEN) && defined(RTLD_NOLOAD)
    {
        Dl_info owner;

        if (dladdr(address, &owner) == 0 || owner.dli_fname == NULL)
            return;
        base = owner.dli_fbase;
        snprintf(st_blas_library, sizeof(st_blas_library), "%s",
                 owner.dli_fname);
        library = dlopen(owner.dli_fname, RTLD_LAZY | RTLD_NOLOAD);
        if (library == NULL)
            return;
    }
#else
    return;
#endif

    st_thread_symbol(library, base, "flexiblas_get_num_threads",
                     &st_blas_get_threads, sizeof(st_blas_get_threads));
    st_thread_symbol(library, base, "flexiblas_set_num_threads",
                     &st_blas_set_threads, sizeof(st_blas_set_threads));
    if (st_blas_get_threads != NULL && st_blas_set_threads != NULL)
        st_blas_backend = "FlexiBLAS";
    else
    {
        st_blas_get_threads = NULL;
        st_blas_set_threads = NULL;
        st_thread_symbol(library, base, "openblas_get_num_threads",
                         &st_blas_get_threads, sizeof(st_blas_get_threads));
        st_thread_symbol(library, base, "openblas_set_num_threads",
                         &st_blas_set_threads, sizeof(st_blas_set_threads));
        if (st_blas_get_threads != NULL && st_blas_set_threads != NULL)
        {
            st_blas_backend = "OpenBLAS";
            st_thread_symbol(library, base, "openblas_get_parallel",
                             &st_blas_get_parallel,
                             sizeof(st_blas_get_parallel));
        }
        else
        {
            st_blas_get_threads = NULL;
            st_blas_set_threads = NULL;
            st_thread_symbol(library, base, "MKL_Get_Max_Threads",
                             &st_blas_get_threads, sizeof(st_blas_get_threads));
            st_thread_symbol(library, base, "MKL_Set_Num_Threads",
                             &st_blas_set_threads, sizeof(st_blas_set_threads));
            if (st_blas_get_threads != NULL && st_blas_set_threads != NULL)
                st_blas_backend = "MKL";
            else
            {
                st_blas_get_threads = NULL;
                st_blas_set_threads = NULL;
            }
        }
    }
#if defined(ST_THREAD_DLOPEN) && defined(RTLD_NOLOAD)

    dlclose(library);
#endif
}

static void st_thread_initialize(void)
{
    if (st_threads_initialized)
        return;
#ifdef MT
    st_threads_requested = omp_get_max_threads();
    if (st_threads_requested > MAXTHREADS)
    {
        st_threads_requested = MAXTHREADS;

        omp_set_num_threads(MAXTHREADS);
    }
    if (st_threads_requested < 1)
        st_threads_requested = 1;
    st_threads_effective = st_threads_requested;
    int limit = omp_get_thread_limit();
    if (limit > 0 && st_threads_effective > limit)
        st_threads_effective = limit;
#endif
    st_thread_blas_initialize();
    st_threads_initialized = 1;
}

int st_num_threads(void)
{
    return st_threads_effective;
}

void gs_set_num_threads(int *activate)
{
    st_thread_initialize();
    if (*activate > MAXTHREADS)
        *activate = MAXTHREADS;
    if (*activate < 1)
        *activate = 1;
    st_threads_requested = *activate;
    st_threads_explicit = 1;
    int blas_threads = *activate;
#ifdef MT
    st_threads_effective = *activate;
    int limit = omp_get_thread_limit();
    if (limit > 0 && st_threads_effective > limit)
        st_threads_effective = limit;
    blas_threads = st_threads_effective;

    omp_set_dynamic(0);
#else

    if ((st_blas_get_parallel && st_blas_get_parallel() == 2) ||
        !strcmp(st_blas_backend, "FlexiBLAS"))
        blas_threads = 1;
#endif
    if (st_blas_set_threads != NULL)
        st_blas_set_threads(blas_threads);

    omp_set_num_threads(st_threads_effective);
}

void gs_get_num_threads(int *activate)
{
    st_thread_initialize();
    *activate = st_num_threads();
    snprintf(gv->msg, sizeof(gv->msg), "Number of threads %i", *activate);
    gs_info(gv, 0, gv->msg);
}

SEXP st_thread_control_native(int complete)
{
    const char *names[] = {
        "requested",           "package.threads",    "explicit",
        "openmp.available",    "openmp.max.threads", "openmp.dynamic",
        "openmp.thread.limit", "blas.backend",       "blas.library",
        "blas.threads",        "blas.control",       "blas.parallel",
        "legacy.outer.threads"};
    int available = 0, maximum = 1, dynamic = 0, limit = 1;
    int blas_threads = NA_INTEGER, blas_parallel = NA_INTEGER;
    SEXP answer, labels;

    st_thread_initialize();
    if (!complete)
        return ScalarInteger(st_num_threads());
#ifdef MT
    available = 1;
    maximum = omp_get_max_threads();
    dynamic = omp_get_dynamic();
    limit = omp_get_thread_limit();
#endif
    if (st_blas_get_threads != NULL)
    {
        blas_threads = st_blas_get_threads();
        if (blas_threads < 1)
            blas_threads = NA_INTEGER;
    }
    if (st_blas_get_parallel != NULL)
        blas_parallel = st_blas_get_parallel();
    answer = PROTECT(allocVector(VECSXP, 13));
    labels = PROTECT(allocVector(STRSXP, 13));
    for (int i = 0; i < 13; i++)
        SET_STRING_ELT(labels, i, mkChar(names[i]));
    setAttrib(answer, R_NamesSymbol, labels);
    SET_VECTOR_ELT(answer, 0, ScalarInteger(st_threads_requested));
    SET_VECTOR_ELT(answer, 1, ScalarInteger(st_num_threads()));
    SET_VECTOR_ELT(answer, 2, ScalarLogical(st_threads_explicit));
    SET_VECTOR_ELT(answer, 3, ScalarLogical(available));
    SET_VECTOR_ELT(answer, 4, ScalarInteger(maximum));
    SET_VECTOR_ELT(answer, 5, ScalarLogical(dynamic));
    SET_VECTOR_ELT(answer, 6, ScalarInteger(limit));
    SET_VECTOR_ELT(answer, 7, mkString(st_blas_backend));
    SET_VECTOR_ELT(answer, 8, mkString(st_blas_library));
    SET_VECTOR_ELT(answer, 9, ScalarInteger(blas_threads));
    SET_VECTOR_ELT(answer, 10, ScalarLogical(st_blas_set_threads != NULL));
    SET_VECTOR_ELT(answer, 11, ScalarInteger(blas_parallel));
    SET_VECTOR_ELT(answer, 12, ScalarInteger(1));
    UNPROTECT(2);
    return answer;
}

SEXP st_thread_control(SEXP complete)
{
    if (TYPEOF(complete) != LGLSXP || XLENGTH(complete) != 1 ||
        LOGICAL(complete)[0] == NA_LOGICAL)
        error("'complete' must be TRUE or FALSE");
    return st_thread_control_native(LOGICAL(complete)[0]);
}

void gs_reset(gs_varset_type *gv)
{
    gs_UNSIGNED i;
    st_marker_prediction_clear(gv);
    st_marker_cache_clear(gv);
    st_data_release(gv->physical_map);
    gv->physical_map = NULL;

    gv->msg[0] = '\0';
    gv->MaxIL = 0;
    gv->MaxML = 0;
    st_text_array_free(gv->IndName, gv->NoInd);
    st_text_array_free(gv->Pedigree, gv->NoInd);
    gv->IndName = gv->Pedigree = NULL;
    gv->NoInd = 0;
    gv->NoMar = 0;
    gv->phased = 0;

    GS_FREE(gv->IndName);
    GS_FREE(gv->MarName);
    st_allele_storage_release(gv);
    GS_FREE(gv->genv);
    GS_FREE(gv->genv1);
    GS_FREE(gv->genv2);
    GS_FREE(gv->NoAll);
    GS_FREE(gv->All_S);
    GS_FREE(gv->All);
    GS_FREE(gv->AllCnt_S);
    GS_FREE(gv->AllCnt);
    GS_FREE(gv->Pic);
    GS_FREE(gv->Mis);
    GS_FREE(gv->IndMis);
    GS_FREE(gv->Z);
    GS_FREE(gv->V);
    GS_FREE(gv->y);

    gv->NoEff = 0;

    GS_FREE(gv->EffNme);
    GS_FREE(gv->EffAll);
    GS_FREE(gv->CmpAll);

    gv->aa = 0;
    gv->aA = 1;
    gv->AA = 2;
    gv->an = 0.5;
    gv->nn = 1;
    gv->nA = 1.5;

    gv->dimu = 0;
    gv->dimA = 0;

    GS_FREE(gv->A);
    GS_FREE(gv->u);
    GS_FREE(gv->b);
    GS_FREE(gv->As);
    GS_FREE(gv->bs);
    GS_FREE(gv->GI);
    GS_FREE(gv->a);
    GS_FREE(gv->p);
    GS_FREE(gv->ladi);
    GS_FREE(gv->est);
    GS_FREE(gv->SMsqm);
    GS_FREE(gv->SMsqe);
    GS_FREE(gv->SMF);
    GS_FREE(gv->SMp);
    GS_FREE(gv->SMvar);
    GS_FREE(gv->sigma);

    gv->sigma_e = 0;
    gv->NLoc = 0;

    GS_FREE(gv->map_S);
    GS_FREE(gv->map_str_S);
    gv->map_str_N = 0;
    GS_FREE(gv->map);
    GS_FREE(gv->mdl);
    GS_FREE(gv->Ref1);
    GS_FREE(gv->Ref2);

    gv->NHap = 0;

    GS_FREE(gv->Hap_S);
    GS_FREE(gv->Hap_str_S);
    gv->Hap_str_N = 0;
    GS_FREE(gv->h_u);
    GS_FREE(gv->h_o);
    GS_FREE(gv->h_n);

    gv->NCrs = 0;
    gv->cross_ready = 0;
    gv->cross_markers_signature = 0;
    gv->cross_effects_signature = 0;
    gv->cross_map_signature = 0;

    GS_FREE(gv->crs_S);
    GS_FREE(gv->crs);

    for (i = 0; i < gv->NChrom; i++)
    {
        GS_FREE(gv->chrom[i].rf);
        GS_FREE(gv->chrom[i].q);
        GS_FREE(gv->chrom[i].ld);
    }

    gv->NChrom = 0;
    GS_FREE(gv->chrom);

    GS_FREE(gv->t0);
    GS_FREE(gv->t1);
    GS_FREE(gv->t2);
    GS_FREE(gv->t3);
    GS_FREE(gv->t4);
    GS_FREE(gv->t5);
    GS_FREE(gv->t6);
    GS_FREE(gv->t7);
    GS_FREE(gv->t8);
    GS_FREE(gv->t9);
    GS_FREE(gv->u0);
}

void gs_init(gs_varset_type *gv)
{
    gv->prediction_state = NULL;
    gv->prediction_role = 0;
    gv->allele_packed = NULL;
    gv->allele_dictionary = NULL;
    gv->allele_packed_bytes = 0;
    gv->allele_bits = 16;
    gv->physical_map = NULL;
    gv->genomic_relationship_A = NULL;
    gv->genomic_relationship_signature = 0;
    gv->genomic_relationship_coding[0] = '\0';
    gv->name[0] = '\0';
    gv->msg[0] = '\0';
    gv->info_level = 0;
    gv->MaxIL = 0;
    gv->MaxML = 0;
    gv->NoInd = 0;
    gv->NoMar = 0;
    gv->phased = 0;
    gv->IndName = NULL;
    gv->Pedigree = NULL;
    gv->MarName = NULL;
    gv->allele_raw1 = NULL;
    gv->allele_raw2 = NULL;
    gv->genv = NULL;
    gv->genv1 = NULL;
    gv->genv2 = NULL;
    gv->NoAll = NULL;
    gv->All_S = NULL;
    gv->All = NULL;
    gv->AllCnt_S = NULL;
    gv->AllCnt = NULL;
    gv->Pic = NULL;
    gv->Mis = NULL;
    gv->IndMis = NULL;
    gv->Z = NULL;
    gv->V = NULL;
    gv->y = NULL;
    gv->NoEff = 0;
    gv->EffNme = NULL;
    gv->EffAll = NULL;
    gv->CmpAll = NULL;
    gv->aa = 0;
    gv->aA = 1;
    gv->AA = 2;
    gv->an = 0.5;
    gv->nn = 1;
    gv->nA = 1.5;
    gv->dimu = 0;
    gv->dimA = 0;
    gv->A = NULL;
    gv->u = NULL;
    gv->b = NULL;
    gv->As = NULL;
    gv->bs = NULL;
    gv->GI = NULL;
    gv->a = NULL;
    gv->p = NULL;
    gv->ladi = NULL;
    gv->est = NULL;
    gv->SMsqm = NULL;
    gv->SMsqe = NULL;
    gv->SMF = NULL;
    gv->SMp = NULL;
    gv->SMvar = NULL;
    gv->sigma = NULL;
    gv->sigma_e = 0;
    gv->NLoc = 0;
    gv->map_S = NULL;
    gv->map_str_S = NULL;
    gv->map_str_N = 0;
    gv->map = NULL;
    gv->mdl = NULL;
    gv->Ref1 = NULL;
    gv->Ref2 = NULL;
    gv->NHap = 0;
    gv->Hap_S = NULL;
    gv->Hap_str_S = NULL;
    gv->Hap_str_N = 0;
    gv->h_u = NULL;
    gv->h_o = NULL;
    gv->h_n = NULL;
    gv->NCrs = 0;
    gv->cross_ready = 0;
    gv->cross_markers_signature = 0;
    gv->cross_effects_signature = 0;
    gv->cross_map_signature = 0;
    gv->crs_S = NULL;
    gv->crs = NULL;
    gv->NChrom = 0;
    gv->chrom = NULL;
    gv->t0 = NULL;
    gv->t1 = NULL;
    gv->t2 = NULL;
    gv->t3 = NULL;
    gv->t4 = NULL;
    gv->t5 = NULL;
    gv->t6 = NULL;
    gv->t7 = NULL;
    gv->t8 = NULL;
    gv->t9 = NULL;
    gv->u0 = NULL;
}

static int st_console_output_allowed(void)
{
#ifdef MT
    if (omp_in_parallel())
        return 0;
#endif

    return 1;
}

void gs_info(gs_varset_type *gv, int level, char *message)
{
    if (gv->info_level >= level)
    {
        if (!st_console_output_allowed())
            return;

        if (1 <= level)
            Rprintf("I (data set '%s'): %s\n", gv->name, message);

        else if (0 == level)
            Rprintf("M (data set '%s'): %s\n", gv->name, message);

        else if (-1 == level)
            Rprintf("W (data set '%s'): %s\n", gv->name, message);

        else if (-2 == level)
            Rprintf("E (data set '%s'): %s\n", gv->name, message);
    }
}

void gs_info_r(gs_varset_type *gv, int level, char *message)
{
    if (gv->info_level >= level)
    {
        if (!st_console_output_allowed())
            return;

        if (1 <= level)
            Rprintf("I (data set '%s'): %s\r", gv->name, message);

        else if (0 == level)
            Rprintf("M (data set '%s'): %s\r", gv->name, message);

        else if (-1 == level)
            Rprintf("W (data set '%s'): %s\r", gv->name, message);

        else if (-2 == level)
            Rprintf("E (data set '%s'): %s\r", gv->name, message);
    }
}

void gs_r_info(int *level, char **message)
{
    gs_info(gv, *level, *message);
}

void gs_set_info_level(gs_varset_type *gv, int *level)
{
    if (*level < -2)
        *level = -2;
    if (*level > 3)
        *level = 3;
    gv->info_level = *level;
    snprintf(gv->msg, sizeof(gv->msg), "Info level %i", gv->info_level);
    gs_info(gv, 3, gv->msg);
}

void gs_get_info_level(gs_varset_type *gv, int *level)
{
    *level = gv->info_level;
}

gs_varset_type *gs_fdta(char *set_name)
{
    gs_UNSIGNED i;
    gs_varset_type *sp = NULL;

    if (0 == strcmp(set_name, "default"))
        sp = gv;

    if (NULL == sp)
    {

        for (i = 0; i < gs_ndta; i++)
            if (0 == strcmp(set_name, gs_data[i]->name))
                sp = gs_data[i];
    }

    if (NULL == sp)
    {

        if (st_hybrid_name_exists(set_name))
            Rf_error("Data set '%s' is already a hybrid data set", set_name);
        gs_ndta++;
        gs_data = realloc(gs_data, gs_ndta * sizeof(gs_varset_type *));
        gs_data[gs_ndta - 1] = malloc(sizeof(gs_varset_type));
        sp = gs_data[gs_ndta - 1];

        gs_init(sp);
        strcpy(sp->name, set_name);
        sp->info_level = gv->info_level;

        snprintf(sp->msg, sizeof(sp->msg), "New data set generated '%s'",
                 sp->name);

        gs_info(sp, 0, sp->msg);
    }

    snprintf(sp->msg, sizeof(sp->msg), "Working on data set '%s'", sp->name);

    gs_info(sp, 2, sp->msg);

    return sp;
}

void gs_get_info_level_GV(int *level, char **set_name)
{
    gs_get_info_level(gs_fdta(*set_name), level);
}

void gs_set_info_level_GV(int *level, char **set_name)
{
    gs_set_info_level(gs_fdta(*set_name), level);
}

void gs_set_all_info_levels(int *level)
{
    gs_set_info_level(gv, level);
    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        gs_set_info_level(gs_data[i], level);
}

int gs_string_cmp(const void *a, const void *b)
{
    const char **ia = (const char **)a;
    const char **ib = (const char **)b;
    return strcmp(*ia, *ib);
}

typedef struct
{
    const char *name;
    gs_UNSIGNED idx;
} gs_name_index_type;

static int gs_name_index_cmp(const void *a, const void *b)
{
    const gs_name_index_type *ia = (const gs_name_index_type *)a;
    const gs_name_index_type *ib = (const gs_name_index_type *)b;
    int cmp = strcmp(ia->name, ib->name);

    if (0 != cmp)
        return cmp;

    return (ia->idx > ib->idx) - (ia->idx < ib->idx);
}

static int gs_size_mul(size_t a, size_t b, size_t *out)
{
    if ((0 != a) && (b > ((size_t)-1) / a))
        return -1;

    *out = a * b;
    return 0;
}

static void *gs_malloc_array(size_t count, size_t element_size)
{
    size_t bytes;

    if (0 != gs_size_mul(count, element_size, &bytes))
        return NULL;

    return malloc(bytes);
}

static void *gs_realloc_array(void *ptr, size_t count, size_t element_size)
{
    size_t bytes;

    if (0 != gs_size_mul(count, element_size, &bytes))
        return NULL;

    return realloc(ptr, bytes);
}

static int gs_size_add(size_t a, size_t b, size_t *result)
{
    if (a > SIZE_MAX - b)
        return -1;

    *result = a + b;
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

static int gs_pool_add_string_size(size_t *total, const char *text)
{
    size_t len;
    size_t next;

    if (NULL == text)
        return -1;

    len = strlen(text);
    if (len >= gs_MAXSTR)
        return -1;

    if (0 != gs_size_add(len, 1, &next))
        return -1;
    if (0 != gs_size_add(*total, next, total))
        return -1;

    return 0;
}

static int gs_pool_put_string(char *pool, size_t pool_n, size_t *offset,
                              const char *text, const char **stored)
{
    size_t len;
    size_t needed;

    if ((NULL == pool) || (NULL == offset) || (NULL == text) ||
        (NULL == stored))
        return -1;

    len = strlen(text);
    if (len >= gs_MAXSTR)
        return -1;

    if (0 != gs_size_add(len, 1, &needed))
        return -1;
    if ((*offset > pool_n) || (needed > pool_n - *offset))
        return -1;

    memcpy(pool + *offset, text, needed);
    *stored = pool + *offset;
    *offset += needed;

    return 0;
}

static int gs_copy_map_point_to_pool(gs_mappoint_type *dst,
                                     const gs_mappoint_type *src,
                                     gs_UNSIGNED nr, char *pool, size_t pool_n,
                                     size_t *offset)
{
    dst->chrom = src->chrom;
    dst->pos = src->pos;
    dst->nr = nr;

    if ((0 !=
         gs_pool_put_string(pool, pool_n, offset, src->name, &(dst->name))) ||
        (0 !=
         gs_pool_put_string(pool, pool_n, offset, src->class, &(dst->class))))
        return -1;

    return 0;
}

typedef struct
{
    int chrom;
    double pos;
    size_t name_offset;
    size_t class_offset;
    gs_UNSIGNED nr;
} gs_mappoint_build_type;

static int gs_string_builder_append(char **pool, size_t *used, size_t *capacity,
                                    const char *text, size_t *stored_offset)
{
    size_t len;
    size_t needed;
    size_t required;

    if ((NULL == pool) || (NULL == used) || (NULL == capacity) ||
        (NULL == text) || (NULL == stored_offset))
        return -1;

    len = strlen(text);
    if (len >= gs_MAXSTR)
        return -1;

    if ((0 != gs_size_add(len, 1, &needed)) ||
        (0 != gs_size_add(*used, needed, &required)))
        return -1;

    if (required > *capacity)
    {
        size_t new_capacity = (0 == *capacity) ? 256 : *capacity;
        char *new_pool;

        while (new_capacity < required)
        {
            if (new_capacity > SIZE_MAX / 2)
            {
                new_capacity = required;
                break;
            }
            new_capacity *= 2;
        }

        new_pool = gs_realloc_array(*pool, new_capacity, sizeof(char));
        if (NULL == new_pool)
            return -1;

        *pool = new_pool;
        *capacity = new_capacity;
    }

    *stored_offset = *used;
    memcpy(*pool + *used, text, needed);
    *used = required;

    return 0;
}

static int gs_make_name_index(char (*names)[gs_MAXNME], gs_UNSIGNED n,
                              gs_name_index_type **index)
{
    size_t bytes;
    gs_name_index_type *p;

    *index = NULL;

    if (0 == n)
        return 0;

    if (0 != gs_size_mul((size_t)n, sizeof(*p), &bytes))
        return -1;

    p = malloc(bytes);
    if (NULL == p)
        return -1;

    for (gs_UNSIGNED i = 0; i < n; i++)
    {
        p[i].name = names[i];
        p[i].idx = i;
    }

    qsort(p, n, sizeof(*p), gs_name_index_cmp);
    *index = p;
    return 0;
}

static int gs_make_ind_name_index(char *const *names, gs_UNSIGNED n,
                                  gs_name_index_type **index)
{
    *index = NULL;
    if (!n)
        return 0;
    gs_name_index_type *p = gs_malloc_array((size_t)n, sizeof(*p));
    if (!p)
        return -1;
    for (gs_UNSIGNED i = 0; i < n; i++)
    {
        p[i].name = names[i];
        p[i].idx = i;
    }
    qsort(p, n, sizeof(*p), gs_name_index_cmp);
    *index = p;
    return 0;
}

static gs_UNSIGNED gs_find_sorted_ind_name(char *const *names, gs_UNSIGNED n,
                                           const char *name)
{
    gs_UNSIGNED lo = 0, hi = n;
    while (lo < hi)
    {
        gs_UNSIGNED mid = lo + (hi - lo) / 2;
        int cmp = strcmp(name, names[mid]);
        if (!cmp)
            return mid;
        if (cmp < 0)
            hi = mid;
        else
            lo = mid + 1;
    }
    return gs_UNSIGNED_MAX;
}

static int gs_name_index_has_duplicates(const gs_name_index_type *index,
                                        gs_UNSIGNED n)
{
    for (gs_UNSIGNED i = 1; i < n; i++)
        if (0 == strcmp(index[i - 1].name, index[i].name))
            return 1;

    return 0;
}

static gs_UNSIGNED gs_name_index_find(const gs_name_index_type *index,
                                      gs_UNSIGNED n, const char *name)
{
    gs_UNSIGNED left = 0;
    gs_UNSIGNED right = n;

    while (left < right)
    {
        gs_UNSIGNED mid = left + (right - left) / 2;
        int cmp = strcmp(name, index[mid].name);

        if (0 == cmp)
            return index[mid].idx;

        if (cmp < 0)
            right = mid;
        else
            left = mid + 1;
    }

    return gs_UNSIGNED_MAX;
}

static gs_UNSIGNED gs_find_sorted_fixed_name(char (*names)[gs_MAXNME],
                                             gs_UNSIGNED n, const char *name)
{
    gs_UNSIGNED left = 0;
    gs_UNSIGNED right = n;

    while (left < right)
    {
        gs_UNSIGNED mid = left + (right - left) / 2;
        int cmp = strcmp(name, names[mid]);

        if (0 == cmp)
            return mid;

        if (cmp < 0)
            right = mid;
        else
            left = mid + 1;
    }

    return gs_UNSIGNED_MAX;
}

static void gs_tmp_init_from(gs_varset_type *tmp, const gs_varset_type *src)
{
    gs_init(tmp);
    snprintf(tmp->name, sizeof(tmp->name), "%s", src->name);
    tmp->info_level = src->info_level;
}

static int gs_prepare_allele_commit(gs_varset_type *tmp)
{
    if (ST_ALLELE_AUTO_PACK && tmp->NoInd && tmp->NoMar &&
        st_allele_bits(tmp) == 16)
        return st_allele_pack_try(tmp, 0);
    return 0;
}

static void gs_tmp_install(gs_varset_type *dst, gs_varset_type *tmp)
{
    gs_reset(dst);
    *dst = *tmp;
    gs_init(tmp);
}

extern void ps_identity_observe(const char *);

static int gs_tmp_commit(gs_varset_type *dst, gs_varset_type *tmp)
{
    if (gs_prepare_allele_commit(tmp))
        return -1;
    gs_tmp_install(dst, tmp);
    for (gs_UNSIGNED i = 0; i < dst->NoInd; i++)
        ps_identity_observe(dst->IndName[i]);
    return 0;
}

static int gs_alloc_genotype_matrix(gs_varset_type *dst, gs_UNSIGNED no_mar,
                                    gs_UNSIGNED no_ind)
{
    size_t n, bytes;

    if ((0 == no_mar) || (0 == no_ind))
        return -1;

    if ((0 != gs_size_mul((size_t)no_mar, (size_t)no_ind, &n)) ||
        (0 != gs_size_mul(n, sizeof(GS_ALLELE_TYPE), &bytes)))
        return -1;

    dst->allele_raw1 = malloc(bytes);
    dst->allele_raw2 = malloc(bytes);

    if (!st_allele_present(dst))
        return -1;

    return 0;
}

static void gs_set_genotypes_missing(gs_varset_type *dst)
{
    size_t n = (size_t)dst->NoMar * (size_t)dst->NoInd;

    for (size_t i = 0; i < n; i++)
    {
        dst->allele_raw1[i] = GS_ALLELE_MISSING;
        dst->allele_raw2[i] = GS_ALLELE_MISSING;
    }
}

static int gs_copy_map_primary(gs_varset_type *dst, const gs_varset_type *src)
{
    size_t bytes;
    size_t offset = 0;
    gs_mappoint_type *map_S = NULL;
    char *map_str_S = NULL;
    gs_mappoint_type **map = NULL;
    gs_UNSIGNED *mdl = NULL;
    gs_UNSIGNED *storage = NULL;
    int result = -2;

    if (dst->physical_map != src->physical_map)
    {
        st_data_release(dst->physical_map);
        dst->physical_map = src->physical_map;
        st_data_retain(dst->physical_map);
    }
    if ((0 == src->NLoc) && (NULL == src->map_S) && (NULL == src->map_str_S) &&
        (0 == src->map_str_N) && (NULL == src->map) && (NULL == src->mdl))
        return 0;

    if ((0 == src->NLoc) || (0 == src->NoMar) || (src->NLoc > src->NoMar) ||
        (NULL == src->map_S) || (NULL == src->map_str_S) ||
        (0 == src->map_str_N) || (NULL == src->map) || (NULL == src->mdl) ||
        !gs_map_pool_valid(src->map_S, src->NLoc, src->map_str_S,
                           src->map_str_N))
        return -2;

    storage = gs_malloc_array((size_t)src->NoMar, sizeof(*storage));
    if (NULL == storage)
        return -1;

    for (gs_UNSIGNED m = 0; m < src->NoMar; m++)
        storage[m] = gs_UNSIGNED_MAX;

    for (gs_UNSIGNED k = 0; k < src->NLoc; k++)
    {
        gs_UNSIGNED m = src->map_S[k].nr;

        if ((m >= src->NoMar) || (gs_UNSIGNED_MAX != storage[m]))
            goto cleanup;

        storage[m] = k;
    }

    for (gs_UNSIGNED lne = 0; lne < src->NLoc; lne++)
    {
        gs_UNSIGNED m = src->mdl[lne];

        if ((NULL == src->map[lne]) || (m >= src->NoMar) ||
            (gs_UNSIGNED_MAX == storage[m]) ||
            (src->map[lne] != &(src->map_S[storage[m]])) ||
            (src->map[lne]->nr != m))
            goto cleanup;
    }

    if (0 != gs_size_mul((size_t)src->NLoc, sizeof(gs_mappoint_type), &bytes))
    {
        result = -1;
        goto cleanup;
    }
    map_S = malloc(bytes);

    map_str_S = malloc(src->map_str_N);

    if (0 != gs_size_mul((size_t)src->NLoc, sizeof(gs_mappoint_type *), &bytes))
    {
        result = -1;
        goto cleanup;
    }
    map = malloc(bytes);

    if (0 != gs_size_mul((size_t)src->NLoc, sizeof(gs_UNSIGNED), &bytes))
    {
        result = -1;
        goto cleanup;
    }
    mdl = malloc(bytes);

    if ((NULL == map_S) || (NULL == map_str_S) || (NULL == map) ||
        (NULL == mdl))
    {
        result = -1;
        goto cleanup;
    }

    memcpy(map_str_S, src->map_str_S, src->map_str_N);
    memcpy(mdl, src->mdl, (size_t)src->NLoc * sizeof(gs_UNSIGNED));

    for (gs_UNSIGNED k = 0; k < src->NLoc; k++)
        if (0 != gs_copy_map_point_to_pool(&(map_S[k]), &(src->map_S[k]),
                                           src->map_S[k].nr, map_str_S,
                                           src->map_str_N, &offset))
            goto cleanup;

    if (offset != src->map_str_N)
        goto cleanup;

    for (gs_UNSIGNED lne = 0; lne < src->NLoc; lne++)
    {
        gs_UNSIGNED m = src->mdl[lne];
        map[lne] = &(map_S[storage[m]]);
    }

    GS_FREE(dst->map_S);
    GS_FREE(dst->map_str_S);
    GS_FREE(dst->map);
    GS_FREE(dst->mdl);

    dst->NLoc = src->NLoc;
    dst->map_S = map_S;
    dst->map_str_S = map_str_S;
    dst->map_str_N = src->map_str_N;
    dst->map = map;
    dst->mdl = mdl;

    map_S = NULL;
    map_str_S = NULL;
    map = NULL;
    mdl = NULL;
    result = 0;

cleanup:
    free(map_S);
    free(map_str_S);
    free(map);
    free(mdl);
    free(storage);
    return result;
}

static int gs_parse_long_strict(const char *text, long *value)
{
    char *end;
    long v;

    if ((NULL == text) || ('\0' == *text))
        return -1;

    errno = 0;
    v = strtol(text, &end, 10);
    if ((ERANGE == errno) || (end == text) || ('\0' != *end))
        return -1;

    *value = v;
    return 0;
}

static int gs_parse_double_strict(char *text, double *value)
{
    char *end;
    double v;

    if ((NULL == text) || ('\0' == *text))
        return -1;

    {
        char *comma = strchr(text, ',');
        if (NULL != comma)
            *comma = '.';
    }
    v = strtod(text, &end);
    if ((end == text) || ('\0' != *end) || !isfinite(v))
        return -1;

    *value = v;
    return 0;
}

static int gs_allele_missing_span(const char *text, size_t len)
{
    while ((0 < len) && isspace((unsigned char)*text))
    {
        text++;
        len--;
    }

    while ((0 < len) && isspace((unsigned char)text[len - 1]))
        len--;

    return ((1 == len && ('-' == text[0] || '.' == text[0])) ||
            (2 == len && '-' == text[0] && '-' == text[1]) ||
            (2 == len && '-' == text[0] && '1' == text[1]) ||
            (2 == len && 'N' == text[0] && 'A' == text[1]));
}

static int gs_parse_allele_span(const char *text, size_t len,
                                GS_ALLELE_TYPE *allele)
{
    char buf[64];
    char *end;
    long v;

    if ((NULL == text) || (NULL == allele))
        return -1;

    while ((0 < len) && isspace((unsigned char)*text))
    {
        text++;
        len--;
    }

    while ((0 < len) && isspace((unsigned char)text[len - 1]))
        len--;

    if (0 == len)
        return -1;

    if (gs_allele_missing_span(text, len))
    {
        *allele = GS_ALLELE_MISSING;
        return 0;
    }

    if (1 == len)
    {
        switch (text[0])
        {
        case 'A':
        case 'a':
            *allele = (GS_ALLELE_TYPE)1;
            return 0;
        case 'C':
        case 'c':
            *allele = (GS_ALLELE_TYPE)2;
            return 0;
        case 'G':
        case 'g':
            *allele = (GS_ALLELE_TYPE)3;
            return 0;
        case 'T':
        case 't':
            *allele = (GS_ALLELE_TYPE)4;
            return 0;
        case 'D':
        case 'd':
            *allele = (GS_ALLELE_TYPE)5;
            return 0;
        default:
            break;
        }
    }

    if (len >= sizeof(buf))
        return -1;

    memcpy(buf, text, len);
    buf[len] = '\0';

    errno = 0;
    v = strtol(buf, &end, 10);
    if ((ERANGE == errno) || (end == buf) || ('\0' != *end) || (v < 1) ||
        (v > GS_ALLELE_MAX))
        return -1;

    *allele = (GS_ALLELE_TYPE)v;
    return 0;
}

static int gs_parse_single_allele(const char *text, GS_ALLELE_TYPE *allele)
{
    if (NULL == text)
        return -1;

    return gs_parse_allele_span(text, strlen(text), allele);
}

static int gs_allele_cmp(const void *a, const void *b)
{
    GS_ALLELE_TYPE aa = *(const GS_ALLELE_TYPE *)a;
    GS_ALLELE_TYPE bb = *(const GS_ALLELE_TYPE *)b;

    return (aa > bb) - (aa < bb);
}

static int gs_aggregate_genotype(GS_ALLELE_TYPE *a1, GS_ALLELE_TYPE *a2,
                                 unsigned char *state, GS_ALLELE_TYPE a,
                                 GS_ALLELE_TYPE b, int pair)
{
    if ((NULL == a1) || (NULL == a2) || (NULL == state))
        return -1;

    if ((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b))
    {
        *state = 2;
        *a1 = GS_ALLELE_MISSING;
        *a2 = GS_ALLELE_MISSING;
        return 0;
    }

    if (2 == *state)
        return 0;

    if (pair)
    {
        if (0 == *state)
        {
            *a1 = a;
            *a2 = b;
            *state = 1;
            return 0;
        }

        if (((*a1 == a) || (*a1 == b)) && ((*a2 == a) || (*a2 == b)))
        {
            *a1 = a;
            *a2 = b;
            return 0;
        }

        return -1;
    }

    if (0 == *state)
    {
        *a1 = a;
        *a2 = a;
        *state = 1;
        return 0;
    }

    if ((*a1 == a) || (*a2 == a))
        return 0;

    if (*a1 == *a2)
    {
        *a2 = a;
        return 0;
    }

    return -1;
}

static int gs_marker_genotypes_valid(const gs_varset_type *gv)
{
    size_t n;

    if ((NULL == gv) || !st_allele_present(gv) || (0 == gv->NoInd) ||
        (0 == gv->NoMar) ||
        (0 != gs_size_mul((size_t)gv->NoInd, (size_t)gv->NoMar, &n)))
        return 0;

    for (size_t m = 0; m < gv->NoMar; m++)
    {
        st_allele_marker_view view = st_allele_marker(gv, m);
        for (size_t i = 0; i < gv->NoInd; i++)
        {
            st_allele_pair pair = st_allele_read_pair(&view, i);
            GS_ALLELE_TYPE a = pair.a, b = pair.b;
            int ma = (GS_ALLELE_MISSING == a);
            int mb = (GS_ALLELE_MISSING == b);

            if (ma != mb)
                return 0;

            if (!ma && ((0 == a) || (0 == b)))
                return 0;
        }
    }
    return 1;
}

static size_t gs_count_sep_fields(const char *line, char sep)
{
    size_t n = 1;

    for (const char *p = line; '\0' != *p; p++)
        if (sep == *p)
            n++;

    return n;
}

static size_t gs_split_sep_fields(char *line, char sep, char **fields,
                                  size_t capacity)
{
    size_t n = 0;
    char *p = line;

    for (;;)
    {
        char *q;

        if (n >= capacity)
            return capacity + 1;

        fields[n++] = p;
        q = strchr(p, sep);
        if (NULL == q)
            break;

        *q = '\0';
        p = q + 1;
    }

    return n;
}

static char *gs_trim_field(char *text)
{
    char *end;

    while (isspace((unsigned char)*text))
        text++;

    end = text + strlen(text);
    while ((end > text) && isspace((unsigned char)end[-1]))
        *--end = '\0';

    if ((end - text >= 2) && ('"' == text[0]) && ('"' == end[-1]))
    {
        text++;
        end[-1] = '\0';
    }

    return text;
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
    {
        for (m = 0; m < gv->NoMar; m++)
        {
            if (!seen[m])
            {
                gs_info(gv, -2,
                        "Marker data are not completely covered by map");
                goto cleanup;
            }
        }
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
        {
            if (0 == strcmp(names[m - 1], names[m]))
            {
                snprintf(gv->msg, sizeof(gv->msg),
                         "Duplicate marker name %s in marker data", names[m]);
                gs_info(gv, -2, gv->msg);
                goto cleanup;
            }
        }
    }

    ok = 1;

cleanup:

    free(seen);
    free(names);

    if (!ok)
    {
        free(storage);
        storage = NULL;
    }

    return storage;
}

int gs_strcmp(const void *a, const void *b)
{
    const char *ia = (const char *)a;
    const char *ib = (const char *)b;
    return strcmp(ia, ib);
}

int gs_intcmp(const void *a, const void *b)
{
    const int *ia = (const int *)a;
    const int *ib = (const int *)b;
    return *ia - *ib;
}

int gs_shortcmp(const void *a, const void *b)
{
    const short *ia = (const short *)a;
    const short *ib = (const short *)b;
    return *ia - *ib;
}

int gs_dbl_cmp(const void *a, const void *b)
{
    const double *ia = (const double *)a;
    const double *ib = (const double *)b;
    return (*ia > *ib) - (*ia < *ib);
}

int gs_dbl_cmp_i(const void *a, const void *b)
{
    const double *ia = (const double *)a;
    const double *ib = (const double *)b;
    return (*ia < *ib) - (*ia > *ib);
}

void gs_marker_stats_01(gs_varset_type *, char **, char **, char **, int *,
                        int *, int *, int *, int *);

long long int st_getline(char **, size_t *, FILE *);

void gs_remove_newline(char *);

int gs_determine_alleles(GS_ALLELE_TYPE *, GS_ALLELE_TYPE *, const char *);

static gs_varset_type *st_find_data_set_no_create(const char *set_name)
{
    if ((NULL == set_name) || ('\0' == set_name[0]))
        return NULL;

    if (0 == strcmp(set_name, "default"))
        return gv;

    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        if ((NULL != gs_data[i]) && (0 == strcmp(set_name, gs_data[i]->name)))
            return gs_data[i];

    return NULL;
}

extern int ps_population_exists(const char *);
extern SEXP ps_population_pedigree(const char *);

static SEXP st_marker_pedigree_object(const gs_varset_type *x)
{
    if (x->NoInd > INT_MAX)
        Rf_error("Dataset exceeds R data-frame limits");
    int n = (int)x->NoInd;
    SEXP out = PROTECT(allocVector(VECSXP, 2)),
         ids = PROTECT(allocVector(STRSXP, n)),
         ped = PROTECT(allocVector(STRSXP, n)),
         nm = PROTECT(allocVector(STRSXP, 2)),
         rn = PROTECT(allocVector(INTSXP, n ? 2 : 0));
    for (int i = 0; i < n; i++)
    {
        SET_STRING_ELT(ids, i, mkCharCE(x->IndName[i], CE_UTF8));
        SET_STRING_ELT(ped, i,
                       mkCharCE(x->Pedigree && x->Pedigree[i] ? x->Pedigree[i]
                                                              : x->IndName[i],
                                CE_UTF8));
    }
    SET_VECTOR_ELT(out, 0, ids);
    SET_VECTOR_ELT(out, 1, ped);
    SET_STRING_ELT(nm, 0, mkChar("name"));
    SET_STRING_ELT(nm, 1, mkChar("pedigree"));
    setAttrib(out, R_NamesSymbol, nm);
    if (n)
    {
        INTEGER(rn)[0] = NA_INTEGER;
        INTEGER(rn)[1] = -n;
    }
    setAttrib(out, R_RowNamesSymbol, rn);
    setAttrib(out, R_ClassSymbol, mkString("data.frame"));
    SEXP version = PROTECT(ScalarInteger(1));
    setAttrib(out, install("version"), version);
    UNPROTECT(6);
    return out;
}

static int st_marker_pedigree_import(gs_varset_type *x, SEXP metadata)
{
    if (isNull(metadata))
        return 0;
    SEXP nm = getAttrib(metadata, R_NamesSymbol),
         version = getAttrib(metadata, install("version"));
    if (TYPEOF(metadata) != VECSXP || LENGTH(metadata) != 2 ||
        TYPEOF(nm) != STRSXP || strcmp(CHAR(STRING_ELT(nm, 0)), "name") ||
        strcmp(CHAR(STRING_ELT(nm, 1)), "pedigree") ||
        TYPEOF(version) != INTSXP || LENGTH(version) != 1 ||
        INTEGER(version)[0] != 1)
        return -1;
    SEXP ids = VECTOR_ELT(metadata, 0), ped = VECTOR_ELT(metadata, 1);
    if (TYPEOF(ids) != STRSXP || TYPEOF(ped) != STRSXP ||
        (uintmax_t)XLENGTH(ids) != (uintmax_t)x->NoInd || (uintmax_t)XLENGTH(ped) != (uintmax_t)x->NoInd)
        return -1;
    gs_name_index_type *index = NULL;
    char **strings =
        calloc((size_t)x->NoInd ? (size_t)x->NoInd : 1, sizeof(*strings));
    if (!strings || gs_make_ind_name_index(x->IndName, x->NoInd, &index))
    {
        free(strings);
        return -1;
    }
    for (R_xlen_t i = 0; i < XLENGTH(ids); i++)
    {
        if (STRING_ELT(ids, i) == NA_STRING || STRING_ELT(ped, i) == NA_STRING)
            goto fail;
        const char *id = CHAR(STRING_ELT(ids, i)),
                   *value = CHAR(STRING_ELT(ped, i));
        gs_UNSIGNED at = gs_name_index_find(index, x->NoInd, id);
        if (at == gs_UNSIGNED_MAX || strings[at] || !*value)
            goto fail;
        strings[at] = st_text_dup(value);
        if (!strings[at])
            goto fail;
    }
    free(index);
    st_text_array_free(x->Pedigree, x->NoInd);
    x->Pedigree = strings;
    return 0;
fail:
    free(index);
    st_text_array_free(strings, x->NoInd);
    return -1;
}

SEXP st_population_pedigree(SEXP population, SEXP source)
{
    if (TYPEOF(population) != STRSXP || XLENGTH(population) != 1 ||
        STRING_ELT(population, 0) == NA_STRING || TYPEOF(source) != STRSXP ||
        XLENGTH(source) != 1 || STRING_ELT(source, 0) == NA_STRING)
        Rf_error("name and source must each be one nonmissing string");
    const char *name = CHAR(STRING_ELT(population, 0)),
               *from = CHAR(STRING_ELT(source, 0));
    if (!*name)
        Rf_error("name must not be empty");
    int simulation = ps_population_exists(name);
    gs_varset_type *marker = st_find_data_set_no_create(name);
    if (!strcmp(from, "auto"))
    {
        if (simulation && marker)
            Rf_error("Name exists in both registries; specify source");
        if (simulation)
            return ps_population_pedigree(name);
        if (marker)
            return st_marker_pedigree_object(marker);
    }
    else if (!strcmp(from, "simulation"))
        return ps_population_pedigree(name);
    else if (!strcmp(from, "marker"))
    {
        if (marker)
            return st_marker_pedigree_object(marker);
    }
    else
        Rf_error("source must be auto, simulation or marker");
    Rf_error("Population or marker dataset '%s' does not exist", name);
    return R_NilValue;
}

gs_varset_type *st_marker_find(const char *name)
{
    return st_find_data_set_no_create(name);
}

gs_varset_type *st_prediction_find_source(const char *modelid, int role)
{
    for (gs_UNSIGNED i = 0; i <= gs_ndta; i++)
    {
        gs_varset_type *x = i == 0 ? gv : gs_data[i - 1];
        if (!x || !x->prediction_state || x->prediction_role != role)
            continue;
        const char *id = st_marker_prediction_id(x);
        if (id && !strcmp(modelid, id))
            return x;
    }
    return NULL;
}

static int st_ram_add(size_t *total, size_t count, size_t element_size)
{
    size_t bytes;
    if ((NULL == total) || (0 != gs_size_mul(count, element_size, &bytes)) ||
        (0 != gs_size_add(*total, bytes, total)))
        return -1;
    return 0;
}

static int st_ram_add_square(size_t *total, size_t n, size_t element_size)
{
    size_t count;
    if (0 != gs_size_mul(n, n, &count))
        return -1;
    return st_ram_add(total, count, element_size);
}

static double st_data_set_ram_core_mb(const gs_varset_type *x)
{
    size_t bytes = 0, n = 0;
    if (NULL == x)
        return NA_REAL;

#define ST_RAM_ADD_IF(PTR, COUNT, TYPE)                                        \
    do                                                                         \
    {                                                                          \
        if ((NULL != (PTR)) &&                                                 \
            (0 != st_ram_add(&bytes, (size_t)(COUNT), sizeof(TYPE))))          \
            return NA_REAL;                                                    \
    } while (0)

    if (st_ram_add(&bytes, st_marker_prediction_bytes(x), 1))
        return NA_REAL;
    ST_RAM_ADD_IF(x->y, x->NoInd, *x->y);
    ST_RAM_ADD_IF(x->IndName, x->NoInd, *x->IndName);
    ST_RAM_ADD_IF(x->Pedigree, x->NoInd, *x->Pedigree);
    for (gs_UNSIGNED i = 0; i < x->NoInd; i++)
    {
        if (x->IndName && x->IndName[i] &&
            st_ram_add(&bytes, strlen(x->IndName[i]) + 1, 1))
            return NA_REAL;
        if (x->Pedigree && x->Pedigree[i] &&
            st_ram_add(&bytes, strlen(x->Pedigree[i]) + 1, 1))
            return NA_REAL;
    }
    ST_RAM_ADD_IF(x->MarName, x->NoMar, *x->MarName);
    if (0 != gs_size_mul((size_t)x->NoMar, (size_t)x->NoInd, &n))
        return NA_REAL;
    ST_RAM_ADD_IF(x->allele_raw1, n, *x->allele_raw1);
    ST_RAM_ADD_IF(x->allele_raw2, n, *x->allele_raw2);
    ST_RAM_ADD_IF(x->allele_packed, x->allele_packed_bytes, unsigned char);
    ST_RAM_ADD_IF(x->allele_dictionary,
                  x->NoMar * (st_allele_bits(x) == 4 ? 3 : 15), unsigned char);

    ST_RAM_ADD_IF(x->map_S, x->NLoc, *x->map_S);
    ST_RAM_ADD_IF(x->map_str_S, x->map_str_N, *x->map_str_S);
    ST_RAM_ADD_IF(x->map, x->NLoc, *x->map);
    ST_RAM_ADD_IF(x->mdl, x->NLoc, *x->mdl);

#undef ST_RAM_ADD_IF
    return (double)bytes / (1024.0 * 1024.0);
}

static double st_data_set_ram_work_mb(const gs_varset_type *x)
{
    size_t bytes = 0, n = 0, total_all = 0;
    size_t dim_max, sigma_n;
    if (NULL == x)
        return NA_REAL;

#define ST_RAM_ADD_IF(PTR, COUNT, TYPE)                                        \
    do                                                                         \
    {                                                                          \
        if ((NULL != (PTR)) &&                                                 \
            (0 != st_ram_add(&bytes, (size_t)(COUNT), sizeof(TYPE))))          \
            return NA_REAL;                                                    \
    } while (0)
#define ST_RAM_ADD_SQ_IF(PTR, N, TYPE)                                         \
    do                                                                         \
    {                                                                          \
        if ((NULL != (PTR)) &&                                                 \
            (0 != st_ram_add_square(&bytes, (size_t)(N), sizeof(TYPE))))       \
            return NA_REAL;                                                    \
    } while (0)

    if (0 != gs_size_mul((size_t)x->NoMar, (size_t)x->NoInd, &n))
        return NA_REAL;
    if (x->physical_map &&
        st_ram_add(&bytes, st_data_bytes(x->physical_map), 1))
        return NA_REAL;
    if (x->genomic_relationship_A &&
        st_ram_add(&bytes, st_data_bytes(x->genomic_relationship_A), 1))
        return NA_REAL;
    ST_RAM_ADD_IF(x->genv, n, *x->genv);
    ST_RAM_ADD_IF(x->genv1, n, *x->genv1);
    ST_RAM_ADD_IF(x->genv2, n, *x->genv2);

    ST_RAM_ADD_IF(x->NoAll, x->NoMar, *x->NoAll);
    ST_RAM_ADD_IF(x->All, x->NoMar, *x->All);
    ST_RAM_ADD_IF(x->AllCnt, x->NoMar, *x->AllCnt);
    if (NULL != x->NoAll)
        for (gs_UNSIGNED m = 0; m < x->NoMar; m++)
            if (0 != gs_size_add(total_all, (size_t)x->NoAll[m], &total_all))
                return NA_REAL;
    ST_RAM_ADD_IF(x->All_S, total_all, *x->All_S);
    ST_RAM_ADD_IF(x->AllCnt_S, total_all, *x->AllCnt_S);
    ST_RAM_ADD_IF(x->Pic, x->NoMar, *x->Pic);
    ST_RAM_ADD_IF(x->Mis, x->NoMar, *x->Mis);
    ST_RAM_ADD_IF(x->IndMis, x->NoInd, *x->IndMis);

    ST_RAM_ADD_IF(x->Z, n, *x->Z);
    ST_RAM_ADD_SQ_IF(x->V, x->NoInd, *x->V);

    ST_RAM_ADD_IF(x->EffNme, x->NoEff, *x->EffNme);
    ST_RAM_ADD_IF(x->EffAll, x->NoEff, *x->EffAll);
    ST_RAM_ADD_IF(x->CmpAll, x->NoEff, *x->CmpAll);

    ST_RAM_ADD_SQ_IF(x->A, x->dimA, *x->A);
    ST_RAM_ADD_IF(x->u, x->dimu, *x->u);
    ST_RAM_ADD_IF(x->b, x->dimA, *x->b);
    ST_RAM_ADD_SQ_IF(x->As, x->dimA, *x->As);
    ST_RAM_ADD_IF(x->bs, x->dimA, *x->bs);
    ST_RAM_ADD_SQ_IF(x->GI, x->NoInd, *x->GI);
    ST_RAM_ADD_IF(x->a, x->dimA, *x->a);
    ST_RAM_ADD_IF(x->p, x->dimu, *x->p);
    ST_RAM_ADD_IF(x->ladi, x->NoMar, *x->ladi);
    ST_RAM_ADD_IF(x->est, x->NoInd, *x->est);
    ST_RAM_ADD_IF(x->SMsqm, x->NoMar, *x->SMsqm);
    ST_RAM_ADD_IF(x->SMsqe, x->NoMar, *x->SMsqe);
    ST_RAM_ADD_IF(x->SMF, x->NoMar, *x->SMF);
    ST_RAM_ADD_IF(x->SMp, x->NoMar, *x->SMp);
    ST_RAM_ADD_IF(x->SMvar, x->NoMar, *x->SMvar);
    sigma_n = (size_t)x->dimA;
    if (sigma_n < (size_t)x->NoMar + 1u)
        sigma_n = (size_t)x->NoMar + 1u;
    ST_RAM_ADD_IF(x->sigma, sigma_n, *x->sigma);

    ST_RAM_ADD_IF(x->Ref1, n, *x->Ref1);
    ST_RAM_ADD_IF(x->Ref2, n, *x->Ref2);
    ST_RAM_ADD_IF(x->Hap_S, x->NHap, *x->Hap_S);
    ST_RAM_ADD_IF(x->Hap_str_S, x->Hap_str_N, *x->Hap_str_S);
    ST_RAM_ADD_IF(x->h_u, x->NHap, *x->h_u);
    ST_RAM_ADD_IF(x->h_o, x->NHap, *x->h_o);
    ST_RAM_ADD_IF(x->h_n, x->NHap, *x->h_n);
    ST_RAM_ADD_IF(x->crs_S, x->NCrs, *x->crs_S);
    ST_RAM_ADD_IF(x->crs, x->NCrs, *x->crs);

    ST_RAM_ADD_IF(x->chrom, x->NChrom, *x->chrom);
    if (NULL != x->chrom)
        for (gs_UNSIGNED c = 0; c < x->NChrom; c++)
        {
            ST_RAM_ADD_SQ_IF(x->chrom[c].rf, x->chrom[c].NLocC,
                             *x->chrom[c].rf);
            ST_RAM_ADD_SQ_IF(x->chrom[c].q, x->chrom[c].NLocC, *x->chrom[c].q);
            ST_RAM_ADD_SQ_IF(x->chrom[c].ld, x->chrom[c].NLocC,
                             *x->chrom[c].ld);
        }

    dim_max = (x->NoInd > x->NoMar) ? (size_t)x->NoInd : (size_t)x->NoMar;
    ST_RAM_ADD_IF(x->t0, x->dimA, *x->t0);
    ST_RAM_ADD_IF(x->t1, 1, *x->t1);
    ST_RAM_ADD_IF(x->t2, dim_max, *x->t2);
    ST_RAM_ADD_IF(x->t3, dim_max, *x->t3);
    ST_RAM_ADD_SQ_IF(x->t4, x->dimA, *x->t4);
    ST_RAM_ADD_SQ_IF(x->t5, x->NoInd, *x->t5);
    ST_RAM_ADD_SQ_IF(x->t6, x->NoInd, *x->t6);

    ST_RAM_ADD_SQ_IF(x->t9, x->NoInd, *x->t9);
    {
        size_t u0_n = (x->dimA > x->NoInd) ? (size_t)x->dimA : (size_t)x->NoInd;
        ST_RAM_ADD_IF(x->u0, u0_n, *x->u0);
    }

#undef ST_RAM_ADD_SQ_IF
#undef ST_RAM_ADD_IF
    return (double)bytes / (1024.0 * 1024.0);
}

static void st_compact_data_set_native(gs_varset_type *x)
{
    st_marker_cache_clear(x);
    if (NULL == x)
        return;

    GS_FREE(x->genv);
    GS_FREE(x->genv1);
    GS_FREE(x->genv2);
    GS_FREE(x->NoAll);
    GS_FREE(x->All_S);
    GS_FREE(x->All);
    GS_FREE(x->AllCnt_S);
    GS_FREE(x->AllCnt);
    GS_FREE(x->Pic);
    GS_FREE(x->Mis);
    GS_FREE(x->IndMis);
    GS_FREE(x->Z);
    GS_FREE(x->V);
    GS_FREE(x->y);

    x->NoEff = 0;
    GS_FREE(x->EffNme);
    GS_FREE(x->EffAll);
    GS_FREE(x->CmpAll);
    x->aa = 0;
    x->aA = 1;
    x->AA = 2;
    x->an = 0.5;
    x->nn = 1;
    x->nA = 1.5;

    x->dimu = 0;
    x->dimA = 0;
    GS_FREE(x->A);
    GS_FREE(x->u);
    GS_FREE(x->b);
    GS_FREE(x->As);
    GS_FREE(x->bs);
    GS_FREE(x->GI);
    GS_FREE(x->a);
    GS_FREE(x->p);
    GS_FREE(x->ladi);
    GS_FREE(x->est);
    GS_FREE(x->SMsqm);
    GS_FREE(x->SMsqe);
    GS_FREE(x->SMF);
    GS_FREE(x->SMp);
    GS_FREE(x->SMvar);
    GS_FREE(x->sigma);
    x->sigma_e = 0;

    GS_FREE(x->Ref1);
    GS_FREE(x->Ref2);
    x->NHap = 0;
    GS_FREE(x->Hap_S);
    GS_FREE(x->Hap_str_S);
    x->Hap_str_N = 0;
    GS_FREE(x->h_u);
    GS_FREE(x->h_o);
    GS_FREE(x->h_n);
    x->NCrs = 0;
    x->cross_ready = 0;
    x->cross_markers_signature = 0;
    x->cross_effects_signature = 0;
    x->cross_map_signature = 0;
    GS_FREE(x->crs_S);
    GS_FREE(x->crs);

    if (NULL != x->chrom)
        for (gs_UNSIGNED c = 0; c < x->NChrom; c++)
        {
            GS_FREE(x->chrom[c].rf);
            GS_FREE(x->chrom[c].q);
            GS_FREE(x->chrom[c].ld);
        }
    x->NChrom = 0;
    GS_FREE(x->chrom);

    GS_FREE(x->t0);
    GS_FREE(x->t1);
    GS_FREE(x->t2);
    GS_FREE(x->t3);
    GS_FREE(x->t4);
    GS_FREE(x->t5);
    GS_FREE(x->t6);
    GS_FREE(x->t7);
    GS_FREE(x->t8);
    GS_FREE(x->t9);
    GS_FREE(x->u0);
}

static void st_data_set_summary_row(const gs_varset_type *x, R_xlen_t row,
                                    SEXP name_s, SEXP no_mar_s, SEXP no_ind_s,
                                    SEXP hzy_s, SEXP mm_s, SEXP phased_s,
                                    SEXP map_s, SEXP no_chr_s, SEXP eff_s,
                                    SEXP pv_s, SEXP gv_s, SEXP ram_core_s,
                                    SEXP ram_work_s)
{
    double complete = 0.0, heterozygous = 0.0, missing = 0.0, total = 0.0;
    gs_UNSIGNED no_chr = 0;
    st_marker_cache_current((gs_varset_type *)x);

    SET_STRING_ELT(name_s, row, mkChar(x->name));
    REAL(no_mar_s)[row] = (double)x->NoMar;
    REAL(no_ind_s)[row] = (double)x->NoInd;

    if ((0 < x->NoMar) && (0 < x->NoInd) && st_allele_present(x))
    {
        size_t n;
        if (0 == gs_size_mul((size_t)x->NoMar, (size_t)x->NoInd, &n))
        {
            total = (double)n;
            unsigned int bits = st_allele_bits(x);
            if (bits == 16)
            {
                for (size_t i = 0; i < n; i++)
                {
                    GS_ALLELE_TYPE a = st_allele_get(x, i, 0),
                                   b = st_allele_get(x, i, 1);
                    if (a == GS_ALLELE_MISSING || b == GS_ALLELE_MISSING)
                        missing += 1;
                    else
                    {
                        complete += 1;
                        if (a != b)
                            heterozygous += 1;
                    }
                }
            }
            else if (bits == 4)
            {

                for (size_t byte = 0; byte < x->allele_packed_bytes; byte++)
                {
                    unsigned int value = atomic_load_explicit(
                        (const _Atomic unsigned char *)(x->allele_packed +
                                                        byte),
                        memory_order_relaxed);
                    for (unsigned int h = 0; h < 2 && byte * 2 + h < n; h++)
                    {
                        unsigned int pair = (value >> (4 * h)) & 15,
                                     a = pair & 3, b = pair >> 2;
                        if (a == 3 || b == 3)
                            missing += 1;
                        else
                        {
                            complete += 1;
                            if (a != b)
                                heterozygous += 1;
                        }
                    }
                }
            }
            else
            {
                for (size_t i = 0; i < n; i++)
                {
                    unsigned int value = atomic_load_explicit(
                        (const _Atomic unsigned char *)(x->allele_packed + i),
                        memory_order_relaxed);
                    unsigned int a = value & 15, b = value >> 4;
                    if (a == 15 || b == 15)
                        missing += 1;
                    else
                    {
                        complete += 1;
                        if (a != b)
                            heterozygous += 1;
                    }
                }
            }
        }
    }

    REAL(hzy_s)[row] = (complete > 0.0) ? heterozygous / complete : NA_REAL;
    REAL(mm_s)[row] = (total > 0.0) ? missing / total : NA_REAL;
    LOGICAL(phased_s)[row] = x->phased ? TRUE : FALSE;
    LOGICAL(map_s)[row] = ((0 < x->NLoc) && (NULL != x->map)) ? TRUE : FALSE;

    if ((0 < x->NLoc) && (NULL != x->map))
    {
        int have_previous = 0;
        int previous = 0;
        for (gs_UNSIGNED l = 0; l < x->NLoc; l++)
            if (NULL != x->map[l])
            {
                int chrom = x->map[l]->chrom;
                if (!have_previous || (chrom != previous))
                {
                    no_chr++;
                    previous = chrom;
                    have_previous = 1;
                }
            }
    }
    REAL(no_chr_s)[row] = (double)no_chr;

    LOGICAL(eff_s)[row] = ((0 < x->NoEff) && (NULL != x->u)) ? TRUE : FALSE;
    LOGICAL(pv_s)[row] = (NULL != x->y) ? TRUE : FALSE;
    LOGICAL(gv_s)
    [row] =
        ((NULL != x->genv) || (NULL != x->genv1) || (NULL != x->genv2)) ? TRUE
                                                                        : FALSE;
    REAL(ram_core_s)[row] = st_data_set_ram_core_mb(x);
    REAL(ram_work_s)[row] = st_data_set_ram_work_mb(x);
}

void st_allele_statistics_clear(gs_varset_type *x)
{
    GS_FREE(x->NoAll);
    GS_FREE(x->All_S);
    GS_FREE(x->All);
    GS_FREE(x->AllCnt_S);
    GS_FREE(x->AllCnt);
    GS_FREE(x->Pic);
    GS_FREE(x->Mis);
    GS_FREE(x->IndMis);
}

int st_allele_statistics_prepare(gs_varset_type *x)
{
    char *empty = "";
    int zero = 0, status = -2;
    gs_marker_stats_01(x, &empty, &empty, &empty, &zero, &zero, &zero, &zero,
                       &status);
    return status < 0 ? -1 : 0;
}

void st_allele_statistics_move(gs_varset_type *dst, gs_varset_type *src)
{
    st_allele_statistics_clear(dst);
    dst->NoAll = src->NoAll;
    src->NoAll = NULL;
    dst->All_S = src->All_S;
    src->All_S = NULL;
    dst->All = src->All;
    src->All = NULL;
    dst->AllCnt_S = src->AllCnt_S;
    src->AllCnt_S = NULL;
    dst->AllCnt = src->AllCnt;
    src->AllCnt = NULL;
    dst->Pic = src->Pic;
    src->Pic = NULL;
    dst->Mis = src->Mis;
    src->Mis = NULL;
    dst->IndMis = src->IndMis;
    src->IndMis = NULL;

    GS_FREE(dst->Z);
    GS_FREE(dst->V);
    GS_FREE(dst->A);
}

SEXP st_data_set_list(void)
{
    static const char *names[] = {"data.set",
                                  "no.markers",
                                  "no.individuals",
                                  "heterozygosity",
                                  "missing.markers",
                                  "phased",
                                  "map",
                                  "no.chromosomes",
                                  "effects",
                                  "phenotypic.value",
                                  "genotypic.value",
                                  "ram.core.MB",
                                  "ram.work.MB",
                                  "A_matrix",
                                  "BitsPerMarker",
                                  "Z_matrix"};

    enum
    {
        NCOL = 16
    };

    R_xlen_t n = 0, row = 0;

    if ((0 < gv->NoMar) || (0 < gv->NoInd))
        n++;
    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        if ((0 < gs_data[i]->NoMar) || (0 < gs_data[i]->NoInd))
            n++;
    SEXP ans = PROTECT(allocVector(VECSXP, NCOL));
    SEXP nms = PROTECT(allocVector(STRSXP, NCOL));
    SEXP row_names = PROTECT(allocVector(INTSXP, 2));
    SEXP name_s = PROTECT(allocVector(STRSXP, n));
    SEXP no_mar_s = PROTECT(allocVector(REALSXP, n));
    SEXP no_ind_s = PROTECT(allocVector(REALSXP, n));
    SEXP hzy_s = PROTECT(allocVector(REALSXP, n));
    SEXP mm_s = PROTECT(allocVector(REALSXP, n));
    SEXP phased_s = PROTECT(allocVector(LGLSXP, n));
    SEXP map_s = PROTECT(allocVector(LGLSXP, n));
    SEXP no_chr_s = PROTECT(allocVector(REALSXP, n));
    SEXP eff_s = PROTECT(allocVector(LGLSXP, n));
    SEXP pv_s = PROTECT(allocVector(LGLSXP, n));
    SEXP gv_s = PROTECT(allocVector(LGLSXP, n));
    SEXP ram_core_s = PROTECT(allocVector(REALSXP, n));
    SEXP ram_work_s = PROTECT(allocVector(REALSXP, n));
    SEXP genomic_s = PROTECT(allocVector(LGLSXP, n));
    SEXP z_s = PROTECT(allocVector(LGLSXP, n));

    if ((0 < gv->NoMar) || (0 < gv->NoInd))
        st_data_set_summary_row(gv, row++, name_s, no_mar_s, no_ind_s, hzy_s,
                                mm_s, phased_s, map_s, no_chr_s, eff_s, pv_s,
                                gv_s, ram_core_s, ram_work_s);
    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        if ((0 < gs_data[i]->NoMar) || (0 < gs_data[i]->NoInd))
            st_data_set_summary_row(gs_data[i], row++, name_s, no_mar_s,
                                    no_ind_s, hzy_s, mm_s, phased_s, map_s,
                                    no_chr_s, eff_s, pv_s, gv_s, ram_core_s,
                                    ram_work_s);

    SET_VECTOR_ELT(ans, 0, name_s);
    SET_VECTOR_ELT(ans, 1, no_mar_s);
    SET_VECTOR_ELT(ans, 2, no_ind_s);
    SET_VECTOR_ELT(ans, 3, hzy_s);
    SET_VECTOR_ELT(ans, 4, mm_s);
    SET_VECTOR_ELT(ans, 5, phased_s);
    SET_VECTOR_ELT(ans, 6, map_s);
    SET_VECTOR_ELT(ans, 7, no_chr_s);
    SET_VECTOR_ELT(ans, 8, eff_s);
    SET_VECTOR_ELT(ans, 9, pv_s);
    SET_VECTOR_ELT(ans, 10, gv_s);
    SET_VECTOR_ELT(ans, 11, ram_core_s);
    SET_VECTOR_ELT(ans, 12, ram_work_s);
    row = 0;
    if ((0 < gv->NoMar) || (0 < gv->NoInd))
    {
        LOGICAL(genomic_s)[row] = gv->genomic_relationship_A != NULL;
        LOGICAL(z_s)[row++] = gv->Z != NULL;
    }
    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        if ((0 < gs_data[i]->NoMar) || (0 < gs_data[i]->NoInd))
        {
            LOGICAL(genomic_s)
            [row] = gs_data[i]->genomic_relationship_A != NULL;
            LOGICAL(z_s)[row++] = gs_data[i]->Z != NULL;
        }
    SET_VECTOR_ELT(ans, 13, genomic_s);
    SET_VECTOR_ELT(ans, 15, z_s);
    SET_VECTOR_ELT(ans, 14, allocVector(INTSXP, n));
    row = 0;
    if ((0 < gv->NoMar) || (0 < gv->NoInd))
        INTEGER(VECTOR_ELT(ans, 14))[row++] = (int)st_allele_bits(gv);
    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        if ((0 < gs_data[i]->NoMar) || (0 < gs_data[i]->NoInd))
            INTEGER(VECTOR_ELT(ans, 14))
            [row++] = (int)st_allele_bits(gs_data[i]);
    for (int j = 0; j < NCOL; j++)
        SET_STRING_ELT(nms, j, mkChar(names[j]));
    setAttrib(ans, R_NamesSymbol, nms);
    if (0 < n)
    {
        INTEGER(row_names)[0] = NA_INTEGER;
        INTEGER(row_names)[1] = -(int)n;
        setAttrib(ans, R_RowNamesSymbol, row_names);
    }
    else
    {
        SEXP empty_row_names = PROTECT(allocVector(INTSXP, 0));
        setAttrib(ans, R_RowNamesSymbol, empty_row_names);
        UNPROTECT(1);
    }
    {
        SEXP cls = PROTECT(mkString("data.frame"));
        setAttrib(ans, R_ClassSymbol, cls);
        UNPROTECT(1);
    }

    UNPROTECT(18);
    return ans;
}

SEXP st_compact_data_set(SEXP dataset)
{
    const char *name;
    gs_varset_type *x;
    double before, after;

    if ((STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
        (NA_STRING == STRING_ELT(dataset, 0)) ||
        ('\0' == CHAR(STRING_ELT(dataset, 0))[0]))
        Rf_error("data.set must be one non-empty character string");

    name = CHAR(STRING_ELT(dataset, 0));
    x = st_find_data_set_no_create(name);
    if (NULL == x)
        Rf_error("Data set '%s' does not exist", name);

    before = st_data_set_ram_work_mb(x);
    st_compact_data_set_native(x);
    after = st_data_set_ram_work_mb(x);
    if (ISNA(before) || ISNA(after))
        return ScalarReal(NA_REAL);
    return ScalarReal(before - after);
}

SEXP st_data_set_remove(SEXP dataset)
{
    const char *name;

    if ((STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
        (NA_STRING == STRING_ELT(dataset, 0)) ||
        ('\0' == CHAR(STRING_ELT(dataset, 0))[0]))
        Rf_error("data.set must be one non-empty character string");

    name = CHAR(STRING_ELT(dataset, 0));

    if (0 == strcmp(name, "default"))
    {
        st_hybrid_require_unreferenced(gv);
        gs_reset(gv);
        return ScalarLogical(1);
    }

    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        if ((NULL != gs_data[i]) && (0 == strcmp(name, gs_data[i]->name)))
        {
            gs_varset_type **smaller;

            st_hybrid_require_unreferenced(gs_data[i]);
            gs_reset(gs_data[i]);
            free(gs_data[i]);

            for (gs_UNSIGNED j = i + 1; j < gs_ndta; j++)
                gs_data[j - 1] = gs_data[j];

            gs_ndta--;
            if (0 == gs_ndta)
            {
                free(gs_data);
                gs_data = NULL;
            }
            else
            {
                smaller = gs_realloc_array(gs_data, (size_t)gs_ndta,
                                           sizeof(*gs_data));
                if (NULL != smaller)
                    gs_data = smaller;
            }

            return ScalarLogical(1);
        }

    return ScalarLogical(0);
}

SEXP st_data_set_remove_all(void)
{
    gs_UNSIGNED removed = gs_ndta;

    st_hybrid_require_unreferenced(gv);
    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        st_hybrid_require_unreferenced(gs_data[i]);
    gs_reset(gv);
    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        if (NULL != gs_data[i])
        {
            gs_reset(gs_data[i]);
            free(gs_data[i]);
        }

    free(gs_data);
    gs_data = NULL;
    gs_ndta = 0;

    return ScalarReal((double)removed);
}

static int st_build_marker_subset(const gs_varset_type *src,
                                  const gs_UNSIGNED *keep, gs_UNSIGNED nkeep,
                                  gs_varset_type *tmp)
{
    gs_name_index_type *map_index = NULL;
    gs_name_index_type *selected_index = NULL;
    size_t map_str_offset = 0;
    gs_UNSIGNED nmap_selected = 0;
    int srv = 0;
    char *empty = "";
    char **empty_p = &empty;
    int zero = 0;
    int *zero_p = &zero;

    if ((NULL == src) || (NULL == keep) || (0 == nkeep) || (0 == src->NoInd) ||
        (0 == src->NoMar) || (NULL == src->IndName) || (NULL == src->MarName) ||
        !st_allele_present(src))
        return -1;

    gs_tmp_init_from(tmp, src);
    tmp->NoInd = src->NoInd;
    tmp->NoMar = nkeep;
    tmp->phased = src->phased;
    tmp->MaxIL = src->MaxIL;
    tmp->MaxML = src->MaxML;

    tmp->IndName = calloc((size_t)tmp->NoInd, sizeof(*tmp->IndName));
    tmp->MarName = gs_malloc_array((size_t)tmp->NoMar, sizeof(*tmp->MarName));
    if ((NULL == tmp->IndName) || (NULL == tmp->MarName) ||
        (0 != gs_alloc_genotype_matrix(tmp, tmp->NoMar, tmp->NoInd)))
        goto fail;

    for (gs_UNSIGNED i = 0; i < tmp->NoInd; i++)
        if (st_identity_copy(tmp, i, src, i))
            goto fail;

    for (gs_UNSIGNED m = 0; m < nkeep; m++)
    {
        gs_UNSIGNED old = keep[m];
        if (old >= src->NoMar)
            goto fail;

        strcpy(tmp->MarName[m], src->MarName[old]);
        st_allele_copy(src, 0, (size_t)old * src->NoInd, (size_t)tmp->NoInd,
                       &(tmp->allele_raw1[(size_t)m * tmp->NoInd]));
        st_allele_copy(src, 1, (size_t)old * src->NoInd, (size_t)tmp->NoInd,
                       &(tmp->allele_raw2[(size_t)m * tmp->NoInd]));
    }

    if (NULL != src->y)
    {
        tmp->y = gs_malloc_array((size_t)tmp->NoInd, sizeof(*tmp->y));
        if (NULL == tmp->y)
            goto fail;
        memcpy(tmp->y, src->y, (size_t)tmp->NoInd * sizeof(*tmp->y));
    }

    if (0 < src->NLoc)
    {
        if ((NULL == src->map) || (NULL == src->map_S))
            goto fail;

        map_index = gs_malloc_array((size_t)src->NLoc, sizeof(*map_index));
        if (NULL == map_index)
            goto fail;

        for (gs_UNSIGNED l = 0; l < src->NLoc; l++)
        {
            if ((NULL == src->map[l]) || (NULL == src->map[l]->name) ||
                (NULL == src->map[l]->class))
                goto fail;
            map_index[l].name = src->map[l]->name;
            map_index[l].idx = l;
        }
        qsort(map_index, src->NLoc, sizeof(*map_index), gs_name_index_cmp);
        if (gs_name_index_has_duplicates(map_index, src->NLoc))
            goto fail;

        if ((0 !=
             gs_make_name_index(tmp->MarName, tmp->NoMar, &selected_index)) ||
            gs_name_index_has_duplicates(selected_index, tmp->NoMar))
            goto fail;

        for (gs_UNSIGNED m = 0; m < nkeep; m++)
            if (gs_UNSIGNED_MAX ==
                gs_name_index_find(map_index, src->NLoc, tmp->MarName[m]))
                goto fail;

        for (gs_UNSIGNED l = 0; l < src->NLoc; l++)
        {
            gs_UNSIGNED m = gs_name_index_find(selected_index, tmp->NoMar,
                                               src->map[l]->name);
            if (gs_UNSIGNED_MAX == m)
                continue;

            if ((0 != gs_pool_add_string_size(&(tmp->map_str_N),
                                              src->map[l]->name)) ||
                (0 != gs_pool_add_string_size(&(tmp->map_str_N),
                                              src->map[l]->class)))
                goto fail;
            nmap_selected++;
        }

        if (nmap_selected != nkeep)
            goto fail;

        tmp->NLoc = nkeep;
        tmp->map_S = gs_malloc_array((size_t)tmp->NLoc, sizeof(*tmp->map_S));
        tmp->map_str_S =
            gs_malloc_array(tmp->map_str_N, sizeof(*tmp->map_str_S));
        tmp->map = gs_malloc_array((size_t)tmp->NLoc, sizeof(*tmp->map));
        tmp->mdl = gs_malloc_array((size_t)tmp->NLoc, sizeof(*tmp->mdl));
        if ((NULL == tmp->map_S) || (NULL == tmp->map_str_S) ||
            (NULL == tmp->map) || (NULL == tmp->mdl))
            goto fail;

        nmap_selected = 0;
        for (gs_UNSIGNED l = 0; l < src->NLoc; l++)
        {
            gs_UNSIGNED m = gs_name_index_find(selected_index, tmp->NoMar,
                                               src->map[l]->name);
            if (gs_UNSIGNED_MAX == m)
                continue;

            if (0 != gs_copy_map_point_to_pool(&(tmp->map_S[nmap_selected]),
                                               src->map[l], nmap_selected,
                                               tmp->map_str_S, tmp->map_str_N,
                                               &map_str_offset))
                goto fail;
            tmp->map[nmap_selected] = &(tmp->map_S[nmap_selected]);
            tmp->mdl[nmap_selected] = m;
            nmap_selected++;
        }

        if ((nmap_selected != tmp->NLoc) || (map_str_offset != tmp->map_str_N))
            goto fail;
    }

    gs_marker_stats_01(tmp, empty_p, empty_p, empty_p, zero_p, zero_p, zero_p,
                       zero_p, &srv);
    if (-2 == srv)
        goto fail;

    tmp->physical_map = src->physical_map;
    st_data_retain(tmp->physical_map);
    free(map_index);
    free(selected_index);
    return 0;

fail:
    free(map_index);
    free(selected_index);
    gs_reset(tmp);
    return -1;
}

SEXP st_keep_intersection_markers(SEXP dataset1, SEXP dataset2)
{
    const char *name1;
    const char *name2;
    gs_varset_type *set1;
    gs_varset_type *set2;
    gs_name_index_type *index1 = NULL;
    gs_name_index_type *index2 = NULL;
    gs_UNSIGNED *match2 = NULL;
    gs_UNSIGNED *keep1 = NULL;
    gs_UNSIGNED *keep2 = NULL;
    gs_UNSIGNED ncommon = 0;
    int need1 = 0;
    int need2 = 0;
    int have_tmp1 = 0;
    int have_tmp2 = 0;
    gs_varset_type tmp1;
    gs_varset_type tmp2;
    const char *errmsg = NULL;

    if ((STRSXP != TYPEOF(dataset1)) || (1 != XLENGTH(dataset1)) ||
        (NA_STRING == STRING_ELT(dataset1, 0)) ||
        ('\0' == CHAR(STRING_ELT(dataset1, 0))[0]) ||
        (STRSXP != TYPEOF(dataset2)) || (1 != XLENGTH(dataset2)) ||
        (NA_STRING == STRING_ELT(dataset2, 0)) ||
        ('\0' == CHAR(STRING_ELT(dataset2, 0))[0]))
        Rf_error("data.set1 and data.set2 must be non-empty character strings");

    name1 = CHAR(STRING_ELT(dataset1, 0));
    name2 = CHAR(STRING_ELT(dataset2, 0));
    set1 = st_find_data_set_no_create(name1);
    set2 = st_find_data_set_no_create(name2);

    if (NULL == set1)
        Rf_error("Unknown data set '%s'", name1);
    if (NULL == set2)
        Rf_error("Unknown data set '%s'", name2);
    if ((0 == set1->NoMar) || (NULL == set1->MarName) ||
        !st_allele_present(set1))
        Rf_error("No marker data available in data set '%s'", name1);
    if ((0 == set2->NoMar) || (NULL == set2->MarName) ||
        !st_allele_present(set2))
        Rf_error("No marker data available in data set '%s'", name2);

    if (set1 == set2)
        return ScalarReal((double)set1->NoMar);

    if ((0 != gs_make_name_index(set1->MarName, set1->NoMar, &index1)) ||
        (0 != gs_make_name_index(set2->MarName, set2->NoMar, &index2)))
    {
        errmsg = "Memory allocation error";
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index1, set1->NoMar) ||
        gs_name_index_has_duplicates(index2, set2->NoMar))
    {
        errmsg = "Duplicate marker names are not supported";
        goto cleanup;
    }

    match2 = gs_malloc_array((size_t)set1->NoMar, sizeof(*match2));
    if (NULL == match2)
    {
        errmsg = "Memory allocation error";
        goto cleanup;
    }
    for (gs_UNSIGNED i = 0; i < set1->NoMar; i++)
        match2[i] = gs_UNSIGNED_MAX;

    {
        gs_UNSIGNED a = 0;
        gs_UNSIGNED b = 0;
        while ((a < set1->NoMar) && (b < set2->NoMar))
        {
            int cmp = strcmp(index1[a].name, index2[b].name);
            if (0 == cmp)
            {
                match2[index1[a].idx] = index2[b].idx;
                ncommon++;
                a++;
                b++;
            }
            else if (cmp < 0)
                a++;
            else
                b++;
        }
    }

    if (0 == ncommon)
    {
        errmsg = "The two data sets have no markers in common";
        goto cleanup;
    }

    keep1 = gs_malloc_array((size_t)ncommon, sizeof(*keep1));
    keep2 = gs_malloc_array((size_t)ncommon, sizeof(*keep2));
    if ((NULL == keep1) || (NULL == keep2))
    {
        errmsg = "Memory allocation error";
        goto cleanup;
    }

    {
        gs_UNSIGNED k = 0;
        for (gs_UNSIGNED i = 0; i < set1->NoMar; i++)
            if (gs_UNSIGNED_MAX != match2[i])
            {
                keep1[k] = i;
                keep2[k] = match2[i];
                k++;
            }
    }

    need1 = (ncommon != set1->NoMar);
    need2 = (ncommon != set2->NoMar);
    if (!need2)
        for (gs_UNSIGNED k = 0; k < ncommon; k++)
            if (keep2[k] != k)
            {
                need2 = 1;
                break;
            }

    if (need1)
    {
        if (0 != st_build_marker_subset(set1, keep1, ncommon, &tmp1))
        {
            errmsg = "Failed to restrict first data set to common markers";
            goto cleanup;
        }
        have_tmp1 = 1;
    }

    if (need2)
    {
        if (0 != st_build_marker_subset(set2, keep2, ncommon, &tmp2))
        {
            errmsg = "Failed to restrict second data set to common markers";
            goto cleanup;
        }
        have_tmp2 = 1;
    }

    if ((have_tmp1 && gs_prepare_allele_commit(&tmp1)) ||
        (have_tmp2 && gs_prepare_allele_commit(&tmp2)))
    {
        errmsg = "Cannot prepare packed intersection data";
        goto cleanup;
    }
    if (have_tmp1)
    {
        gs_tmp_install(set1, &tmp1);
        have_tmp1 = 0;
    }
    if (have_tmp2)
    {
        gs_tmp_install(set2, &tmp2);
        have_tmp2 = 0;
    }

cleanup:
    if (have_tmp1)
        gs_reset(&tmp1);
    if (have_tmp2)
        gs_reset(&tmp2);
    free(index1);
    free(index2);
    free(match2);
    free(keep1);
    free(keep2);

    if (NULL != errmsg)
        Rf_error("%s", errmsg);

    return ScalarReal((double)ncommon);
}

static int gs_finish_marker_read(gs_varset_type *gv, gs_varset_type *tmp,
                                 char **ind_filename, char **mar_filename,
                                 char **gen_filename, int *auxfiles,
                                 int *retval)
{
    int srv = 0;
    int *sretval = &srv;

    gs_marker_stats_01(tmp, ind_filename, mar_filename, gen_filename, auxfiles,
                       auxfiles, auxfiles, auxfiles, sretval);

    if (-2 == *sretval)
    {
        *retval = -2;
        return -1;
    }

    if (gs_tmp_commit(gv, tmp))
    {
        *retval = -2;
        return -1;
    }
    *retval = 0;
    return 0;
}

void gs_read_marker_data_01(gs_varset_type *gv, char **dta_filename,
                            char **ind_filename, char **mar_filename,
                            char **gen_filename, int *auxfiles, int *retval)
{
    gs_varset_type tmp;
    FILE *fp = NULL;
    char *line = NULL;
    size_t n = 0;
    long long int rc;
    gs_UNSIGNED NoLines = 0;
    gs_UNSIGNED IndSize = 0;
    gs_UNSIGNED MarSize = 0;
    char **ind = NULL;
    char(*mar)[gs_MAXNME] = NULL;
    GS_ALLELE_TYPE *rec_a = NULL;
    GS_ALLELE_TYPE *rec_b = NULL;
    unsigned char *rec_pair = NULL;
    unsigned char *geno_state = NULL;
    char **names = NULL;
    int failed = 0;

    gs_tmp_init_from(&tmp, gv);

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        *retval = -2;
        return;
    }

    rc = st_getline(&line, &n, fp);
    if (rc < 0)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (0 == rc)
    {
        gs_info(gv, -2, "Marker data file contains no header");
        failed = 1;
        goto cleanup;
    }

    while (0 < (rc = st_getline(&line, &n, fp)))
    {
        char *p;
        size_t len;

        gs_remove_newline(line);
        p = strtok(line, gs_tok_del);
        if (NULL == p)
        {
            gs_info(gv, -2, "Empty line in marker data file");
            failed = 1;
            goto cleanup;
        }
        len = strlen(p);
        if (len >= INT_MAX)
        {
            gs_info(gv, -2, "Individual name too long");
            failed = 1;
            goto cleanup;
        }
        if (len > IndSize)
            IndSize = (gs_UNSIGNED)len;

        p = strtok(NULL, gs_tok_del);
        if (NULL == p)
        {
            gs_info(gv, -2, "Marker name missing in marker data file");
            failed = 1;
            goto cleanup;
        }
        len = strlen(p);
        if (len >= gs_MAXNME)
        {
            gs_info(gv, -2, "Marker name too long");
            failed = 1;
            goto cleanup;
        }
        if (len > MarSize)
            MarSize = (gs_UNSIGNED)len;

        if (gs_UNSIGNED_MAX == NoLines)
        {
            gs_info(gv, -2, "Marker data file too large");
            failed = 1;
            goto cleanup;
        }
        NoLines++;
    }

    if (rc < 0)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    fclose(fp);
    fp = NULL;

    if (0 == NoLines)
    {
        gs_info(gv, -2, "Marker data file contains no data lines");
        failed = 1;
        goto cleanup;
    }

    ind = calloc((size_t)NoLines, sizeof(*ind));
    mar = gs_malloc_array((size_t)NoLines, sizeof(*mar));
    rec_a = gs_malloc_array((size_t)NoLines, sizeof(*rec_a));
    rec_b = gs_malloc_array((size_t)NoLines, sizeof(*rec_b));
    rec_pair = gs_malloc_array((size_t)NoLines, sizeof(*rec_pair));
    names = gs_malloc_array((size_t)NoLines, sizeof(*names));

    if ((NULL == ind) || (NULL == mar) || (NULL == rec_a) || (NULL == rec_b) ||
        (NULL == rec_pair) || (NULL == names))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        failed = 1;
        goto cleanup;
    }

    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2, "Error reading marker data file");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED j = 0; j < NoLines; j++)
    {
        char *p;

        rc = st_getline(&line, &n, fp);
        if (rc <= 0)
        {
            gs_info(gv, -2, "Unexpected end of marker data file");
            failed = 1;
            goto cleanup;
        }

        gs_remove_newline(line);
        p = strtok(line, gs_tok_del);
        if (NULL == p)
        {
            gs_info(gv, -2, "Error in marker data file");
            failed = 1;
            goto cleanup;
        }
        if (st_text_set(&ind[j], p))
        {
            failed = 1;
            goto cleanup;
        }

        p = strtok(NULL, gs_tok_del);
        if (NULL == p)
        {
            gs_info(gv, -2, "Error in marker data file");
            failed = 1;
            goto cleanup;
        }
        strcpy(mar[j], p);

        p = strtok(NULL, gs_tok_del);
        if (NULL == p)
        {
            gs_info(gv, -2, "Missing allele field in marker data file");
            failed = 1;
            goto cleanup;
        }

        rec_pair[j] = (NULL != strchr(p, '/'));
        if (rec_pair[j])
        {
            if (0 != gs_determine_alleles(&(rec_a[j]), &(rec_b[j]), p))
            {
                gs_info(gv, -2, "Invalid allele in marker data file");
                failed = 1;
                goto cleanup;
            }
        }
        else
        {
            if (0 != gs_parse_single_allele(p, &(rec_a[j])))
            {
                gs_info(gv, -2, "Invalid allele in marker data file");
                failed = 1;
                goto cleanup;
            }
            rec_b[j] = rec_a[j];
        }
    }

    fclose(fp);
    fp = NULL;

    for (gs_UNSIGNED j = 0; j < NoLines; j++)
        names[j] = ind[j];
    qsort(names, NoLines, sizeof(*names), gs_string_cmp);
    tmp.NoInd = 1;
    for (gs_UNSIGNED j = 1; j < NoLines; j++)
        if (0 != strcmp(names[j - 1], names[j]))
            tmp.NoInd++;

    if ((0 == tmp.NoInd) ||
        ((size_t)tmp.NoInd > ((size_t)-1) / sizeof(*tmp.IndName)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    tmp.IndName = calloc((size_t)tmp.NoInd, sizeof(*tmp.IndName));
    if (NULL == tmp.IndName)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    {
        gs_UNSIGNED k = 0;
        if (st_identity_set(&tmp, k++, names[0], NULL))
        {
            failed = 1;
            goto cleanup;
        }
        for (gs_UNSIGNED j = 1; j < NoLines; j++)
            if (0 != strcmp(names[j - 1], names[j]))
                if (st_identity_set(&tmp, k++, names[j], NULL))
                {
                    failed = 1;
                    goto cleanup;
                }
    }

    for (gs_UNSIGNED j = 0; j < NoLines; j++)
        names[j] = mar[j];
    qsort(names, NoLines, sizeof(*names), gs_string_cmp);
    tmp.NoMar = 1;
    for (gs_UNSIGNED j = 1; j < NoLines; j++)
        if (0 != strcmp(names[j - 1], names[j]))
            tmp.NoMar++;

    if ((0 == tmp.NoMar) ||
        ((size_t)tmp.NoMar > ((size_t)-1) / sizeof(*tmp.MarName)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    tmp.MarName = malloc((size_t)tmp.NoMar * sizeof(*tmp.MarName));
    if (NULL == tmp.MarName)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    {
        gs_UNSIGNED k = 0;
        strcpy(tmp.MarName[k++], names[0]);
        for (gs_UNSIGNED j = 1; j < NoLines; j++)
            if (0 != strcmp(names[j - 1], names[j]))
                strcpy(tmp.MarName[k++], names[j]);
    }

    tmp.MaxIL = (int)IndSize;
    tmp.MaxML = (int)MarSize;

    if (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    gs_set_genotypes_missing(&tmp);

    {
        size_t n_genotypes;
        if (0 !=
            gs_size_mul((size_t)tmp.NoMar, (size_t)tmp.NoInd, &n_genotypes))
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }
        geno_state = calloc(n_genotypes, sizeof(*geno_state));
        if (NULL == geno_state)
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }
    }

    for (gs_UNSIGNED k = 0; k < NoLines; k++)
    {
        gs_UNSIGNED m =
            gs_find_sorted_fixed_name(tmp.MarName, tmp.NoMar, mar[k]);
        gs_UNSIGNED i = gs_find_sorted_ind_name(tmp.IndName, tmp.NoInd, ind[k]);
        size_t idx;
        GS_ALLELE_TYPE *a1;
        GS_ALLELE_TYPE *a2;

        if ((gs_UNSIGNED_MAX == m) || (gs_UNSIGNED_MAX == i))
        {
            gs_info(gv, -2, "Internal marker-data name lookup error");
            failed = 1;
            goto cleanup;
        }

        idx = (size_t)m * (size_t)tmp.NoInd + (size_t)i;
        a1 = &(tmp.allele_raw1[idx]);
        a2 = &(tmp.allele_raw2[idx]);

        if (0 != gs_aggregate_genotype(a1, a2, &(geno_state[idx]), rec_a[k],
                                       rec_b[k], rec_pair[k]))
        {
            gs_info(gv, -2,
                    rec_pair[k]
                        ? "Conflicting genotype records in list marker data"
                        : "More than two alleles for one list genotype");
            failed = 1;
            goto cleanup;
        }
    }

    if (0 != gs_finish_marker_read(gv, &tmp, ind_filename, mar_filename,
                                   gen_filename, auxfiles, retval))
        failed = 1;

cleanup:
    if (NULL != fp)
        fclose(fp);
    free(line);
    st_text_array_free(ind, NoLines);
    free(mar);
    free(rec_a);
    free(rec_b);
    free(rec_pair);
    free(geno_state);
    free(names);

    if (failed)
    {
        gs_reset(&tmp);
        *retval = -2;
        return;
    }

    *retval = 0;
}

void gs_read_marker_data_01_GV(char **dta_filename, char **ind_filename,
                               char **mar_filename, char **gen_filename,
                               int *auxfiles, int *retval, char **set_name)
{
    gs_read_marker_data_01(gs_fdta(*set_name), dta_filename, ind_filename,
                           mar_filename, gen_filename, auxfiles, retval);
}

static void gs_build_Z_impl(gs_varset_type *gv, char **out_filename,
                            int *auxfiles, int *retval, int allow_monomorphic)
{

    gs_UNSIGNED i, m, idx_Z;
    FILE *outfp;

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2, "No marker data available");
        *retval = -2;
        return;
    }

    gv->Z = realloc(gv->Z, gv->NoInd * gv->NoMar * sizeof(double));

    gv->NoEff = 1 + gv->NoMar;
    gv->EffNme = realloc(gv->EffNme, gv->NoEff * sizeof(char[gs_MAXNME]));
    gv->EffAll = realloc(gv->EffAll, gv->NoEff * sizeof(int));
    gv->CmpAll = realloc(gv->CmpAll, gv->NoEff * sizeof(int));

    if ((NULL == gv->Z) || (NULL == gv->EffNme) || (NULL == gv->EffAll) ||
        (NULL == gv->CmpAll))
        ERR_M;

    strcpy(gv->EffNme[0], "mu");
    gv->EffAll[0] = 0;
    gv->CmpAll[0] = 0;
    for (m = 1; m < gv->NoEff; m++)
        strcpy(gv->EffNme[m], gv->MarName[m - 1]);

    int mono_warning = 0, threads = st_num_threads();
    gs_UNSIGNED bad_marker = gv->NoMar;
    (void)threads;

#ifdef MT
#pragma omp parallel for schedule(static) num_threads(                         \
        threads) if ((size_t)gv->NoInd * gv->NoMar >= 65536 && threads > 1)    \
    reduction(| : mono_warning) reduction(min : bad_marker)
#endif
    for (gs_UNSIGNED first = 0; first < gv->NoMar; first += 32)
    {
        gs_UNSIGNED last = gv->NoMar - first < 32 ? gv->NoMar : first + 32;
        for (gs_UNSIGNED marker = first; marker < last; marker++)
        {
            st_allele_marker_view view = st_allele_marker(gv, marker);
            int lar_all = 0, sml_all = INT_MAX, invalid = 0;
            if (NULL == gv->u)
            {
                for (gs_UNSIGNED ind = 0; ind < gv->NoInd; ind++)
                {
                    st_allele_pair pair = st_allele_read_pair(&view, ind);
                    if (pair.a != GS_ALLELE_MISSING)
                    {
                        if (pair.a > lar_all)
                            lar_all = pair.a;
                        if (pair.a > 0 && pair.a < sml_all)
                            sml_all = pair.a;
                    }
                    if (pair.b != GS_ALLELE_MISSING)
                    {
                        if (pair.b > lar_all)
                            lar_all = pair.b;
                        if (pair.b > 0 && pair.b < sml_all)
                            sml_all = pair.b;
                    }
                }
                if (sml_all == INT_MAX)
                    sml_all = lar_all;
                mono_warning |= lar_all == sml_all || lar_all < 2;
            }
            else
            {
                lar_all = gv->EffAll[marker + 1];
                sml_all = gv->CmpAll[marker + 1];
            }
            for (gs_UNSIGNED ind = 0; ind < gv->NoInd; ind++)
            {
                st_allele_pair pair = st_allele_read_pair(&view, ind);
                if (NULL == gv->u && ((pair.a != lar_all && pair.a != sml_all &&
                                       pair.a != GS_ALLELE_MISSING) ||
                                      (pair.b != lar_all && pair.b != sml_all &&
                                       pair.b != GS_ALLELE_MISSING)))
                    invalid = 1;
                gv->Z[(size_t)ind * gv->NoMar + marker] =
                    (allow_monomorphic && lar_all == sml_all) ||
                            (pair.a == lar_all && pair.b == lar_all)
                        ? gv->AA
                    : pair.a == sml_all && pair.b == sml_all ? gv->aa
                                                             : gv->aA;
            }
            if (invalid)
            {
                if (marker < bad_marker)
                    bad_marker = marker;
            }
            else
            {
                gv->EffAll[marker + 1] = lar_all;
                gv->CmpAll[marker + 1] = sml_all;
            }
        }
    }
    if (bad_marker < gv->NoMar)
    {
        snprintf(gv->msg, sizeof(gv->msg), "More than two alleles: %s",
                 gv->MarName[bad_marker]);
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    if (1 == *auxfiles)
    {

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s.%i ", gv->EffNme[m + 1], gv->EffAll[m + 1]);
        fprintf(outfp, "\n");

        for (i = 0; i < gv->NoInd; i++)
        {
            fprintf(outfp, "%s ", gv->IndName[i]);
            for (m = 0; m < gv->NoMar; m++)
            {
                idx_Z = i * gv->NoMar + m;
                fprintf(outfp, "%f ", gv->Z[idx_Z]);
            }
            fprintf(outfp, "\n");
        }

        fclose(outfp);
    }

    if (1 == mono_warning && !allow_monomorphic)
    {

        snprintf(gv->msg, sizeof(gv->msg),
                 "Monomorphic markers in the design matrix");

        gs_info(gv, -1, gv->msg);
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_build_Z_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                   int *retval)
{
    gs_build_Z_impl(gv, out_filename, auxfiles, retval, 0);
}

void st_build_Z_for_A(gs_varset_type *gv, int *retval)
{
    char *empty = "";
    int aux = 0;
    gs_build_Z_impl(gv, &empty, &aux, retval, 1);
}

void gs_build_Z_01_GV(char **out_filename, int *auxfiles, int *retval,
                      char **set_name)
{
    gs_build_Z_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_neg_effall_01(gs_varset_type *gv, int *retval)
{
    gs_UNSIGNED m;
    int c;

    if (NULL == gv->u)
    {
        gs_info(gv, -2, "No effects estimated");
        *retval = -2;
        return;
    }

    for (m = 1; m < gv->NoEff; m++)
    {
        if (gv->u[m] > 0)
        {
            c = gv->EffAll[m];
            gv->EffAll[m] = gv->CmpAll[m];
            gv->CmpAll[m] = c;
        }
        gv->u[m] = 0;
    }

    *retval = 0;
    return;
}

void gs_neg_effall_01_GV(int *retval, char **set_name)
{
    gs_neg_effall_01(gs_fdta(*set_name), retval);
}

void gs_pos_effall_01(gs_varset_type *gv, int *retval)
{
    gs_UNSIGNED m;
    int c;

    if (NULL == gv->u)
    {
        gs_info(gv, -2, "No effects estimated");
        *retval = -2;
        return;
    }

    for (m = 1; m < gv->NoEff; m++)
    {
        if (gv->u[m] < 0)
        {
            c = gv->EffAll[m];
            gv->EffAll[m] = gv->CmpAll[m];
            gv->CmpAll[m] = c;
        }
        gv->u[m] = 0;
    }

    *retval = 0;
    return;
}

void gs_pos_effall_01_GV(int *retval, char **set_name)
{
    gs_pos_effall_01(gs_fdta(*set_name), retval);
}

long long int st_getline(char **, size_t *, FILE *);

void gs_remove_newline(char *);

void gs_read_performance_data_01(gs_varset_type *gv, char **dta_filename,
                                 char **out_filename, int *auxfiles,
                                 int *retval)
{
    FILE *fp = NULL;
    char *line = NULL;
    size_t n = 0;
    long long int rc;
    gs_name_index_type *index = NULL;
    double *values = NULL;
    unsigned char *seen = NULL;
    gs_varset_type tmp;
    gs_UNSIGNED no_seen = 0;
    int failed = 0;

    if ((0 == gv->NoInd) || (0 == gv->NoMar) || (NULL == gv->IndName) ||
        (NULL == gv->MarName) || !st_allele_present(gv))
    {
        gs_info(gv, -2,
                "No marker data available, run 'read.marker.data' first");
        *retval = -2;
        return;
    }

    gs_tmp_init_from(&tmp, gv);

    if (0 != gs_make_ind_name_index(gv->IndName, gv->NoInd, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        *retval = -2;
        return;
    }

    if (gs_name_index_has_duplicates(index, gv->NoInd))
    {
        gs_info(gv, -2, "Duplicate individual name in marker data");
        free(index);
        *retval = -2;
        return;
    }

    values = gs_malloc_array((size_t)gv->NoInd, sizeof(*values));
    seen = calloc((size_t)gv->NoInd, sizeof(*seen));
    if ((NULL == values) || (NULL == seen))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        failed = 1;
        goto cleanup;
    }

    rc = st_getline(&line, &n, fp);
    if (rc < 0)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (0 == rc)
    {
        gs_info(gv, -2, "Performance data file contains no header");
        failed = 1;
        goto cleanup;
    }

    while (0 < (rc = st_getline(&line, &n, fp)))
    {
        char *name;
        char *value_text;
        gs_UNSIGNED i;

        gs_remove_newline(line);
        name = strtok(line, gs_tok_del);
        if (NULL == name)
        {
            gs_info(gv, -2, "Empty line in performance data file");
            failed = 1;
            goto cleanup;
        }

        i = gs_name_index_find(index, gv->NoInd, name);
        if (gs_UNSIGNED_MAX == i)
            continue;

        value_text = strtok(NULL, gs_tok_del);
        if (NULL == value_text)
        {
            gs_info(gv, -2, "Performance value missing");
            failed = 1;
            goto cleanup;
        }

        {
            char *end;
            double value;

            if (strcmp(value_text, "NA") == 0)
                value = NA_REAL;
            else
            {
                value = strtod(value_text, &end);
                if ((end == value_text) || ('\0' != *end) || !isfinite(value))
                {
                    gs_info(gv, -2, "Invalid performance value");
                    failed = 1;
                    goto cleanup;
                }
            }

            values[i] = value;
            if (!seen[i])
            {
                seen[i] = 1;
                no_seen++;
            }
        }
    }

    if (rc < 0)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    fclose(fp);
    fp = NULL;

    if (0 == no_seen)
    {
        gs_info(gv, -2, "No performance data match marker-data individuals");
        failed = 1;
        goto cleanup;
    }

    tmp.NoInd = no_seen;
    tmp.NoMar = gv->NoMar;
    tmp.phased = gv->phased;
    tmp.MaxIL = gv->MaxIL;
    tmp.MaxML = gv->MaxML;

    tmp.IndName = calloc((size_t)tmp.NoInd, sizeof(*tmp.IndName));
    tmp.MarName = gs_malloc_array((size_t)tmp.NoMar, sizeof(*tmp.MarName));
    tmp.y = gs_malloc_array((size_t)tmp.NoInd, sizeof(*tmp.y));

    if ((NULL == tmp.IndName) || (NULL == tmp.MarName) || (NULL == tmp.y) ||
        (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED m = 0; m < tmp.NoMar; m++)
        strcpy(tmp.MarName[m], gv->MarName[m]);

    {
        gs_UNSIGNED j = 0;

        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            if (seen[i])
            {
                if (st_identity_copy(&tmp, j, gv, i))
                {
                    failed = 1;
                    goto cleanup;
                }
                tmp.y[j] = values[i];

                for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
                {
                    size_t old_idx = (size_t)m * (size_t)gv->NoInd + i;
                    size_t new_idx = (size_t)m * (size_t)tmp.NoInd + j;
                    tmp.allele_raw1[new_idx] = st_allele_get(gv, old_idx, 0);
                    tmp.allele_raw2[new_idx] = st_allele_get(gv, old_idx, 1);
                }

                j++;
            }
    }

    {
        int map_rc = gs_copy_map_primary(&tmp, gv);
        if (0 != map_rc)
        {
            gs_info(gv, -2,
                    (-1 == map_rc) ? "Memory allocation error"
                                   : "Invalid linkage map state");
            failed = 1;
            goto cleanup;
        }
    }

    {
        int srv = 0;
        char *empty = "";
        char **empty_p = &empty;
        int zero = 0;
        int *zero_p = &zero;

        gs_marker_stats_01(&tmp, empty_p, empty_p, empty_p, zero_p, zero_p,
                           zero_p, zero_p, &srv);
        if (-2 == srv)
        {
            failed = 1;
            goto cleanup;
        }
    }

    snprintf(gv->msg, sizeof(gv->msg), "No. of performance data: %lu",
             tmp.NoInd);
    gs_info(gv, 0, gv->msg);

    if (1 == *auxfiles)
    {
        FILE *outfp = fopen(*out_filename, "w");
        if (NULL == outfp)
        {
            gs_info(gv, -2, "File open error");
            failed = 1;
            goto cleanup;
        }

        fprintf(outfp, "y\n");
        for (gs_UNSIGNED i = 0; i < tmp.NoInd; i++)
        {
            if (ISNA(tmp.y[i]))
                fprintf(outfp, "%s NA\n", tmp.IndName[i]);
            else
                fprintf(outfp, "%s %f\n", tmp.IndName[i], tmp.y[i]);
        }
        fclose(outfp);
    }

    if (gs_tmp_commit(gv, &tmp))
    {
        failed = 1;
        goto cleanup;
    }
    *retval = 0;

cleanup:
    if (NULL != fp)
        fclose(fp);
    free(line);
    free(index);
    free(values);
    free(seen);

    if (failed)
    {
        gs_reset(&tmp);
        *retval = -2;
    }
}

void gs_read_performance_data_01_GV(char **dta_filename, char **out_filename,
                                    int *auxfiles, int *retval, char **set_name)
{
    gs_read_performance_data_01(gs_fdta(*set_name), dta_filename, out_filename,
                                auxfiles, retval);
}

void gs_return_performance_data_01(gs_varset_type *gv, char **out_filename,
                                   int *auxfiles, int *retval)
{
    gs_UNSIGNED i;

    if (NULL == gv->y)
    {
        gs_info(
            gv, -2,
            "No performance data available, run 'read.performance.data' first");
        *retval = -2;
        return;
    }

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2,
                "No marker data available, run 'read.marker.data' first");
        *retval = -2;
        return;
    }

    if (1 == *auxfiles)
    {

        FILE *outfp;

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "y\n");
        for (i = 0; i < gv->NoInd; i++)
            fprintf(outfp, "%s %f\n", gv->IndName[i], gv->y[i]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_return_performance_data_01_GV(char **out_filename, int *auxfiles,
                                      int *retval, char **set_name)
{
    gs_return_performance_data_01(gs_fdta(*set_name), out_filename, auxfiles,
                                  retval);
}

void gs_lambda_const_01(gs_varset_type *gv, double *lambda, char **out_filename,
                        int *auxfiles, int *retval)
{
    gs_UNSIGNED i, m;

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available");
        *retval = -2;
        return;
    }

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;

    for (i = 0; i < gv->NoMar; i++)
        gv->ladi[i] = *lambda;

    if (1 == *auxfiles)
    {
        FILE *outfp;

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_lambda_const_01_GV(double *lambda, char **out_filename, int *auxfiles,
                           int *retval, char **set_name)

{
    gs_lambda_const_01(gs_fdta(*set_name), lambda, out_filename, auxfiles,
                       retval);
}

void gs_lambda_var_01(gs_varset_type *gv, double *lambda, int *n,
                      char **out_filename, int *auxfiles, int *retval)
{
    gs_UNSIGNED m;

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available");
        *retval = -2;
        return;
    }

    if ((gs_UNSIGNED)*n != gv->NoMar)
    {
        gs_info(gv, -2, "Incorrect number of lambda values");
        *retval = -2;
        return;
    }

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;

    memcpy(gv->ladi, lambda, gv->NoMar * sizeof(double));

    if (1 == *auxfiles)
    {
        FILE *outfp;

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_lambda_var_01_GV(double *lambda, int *n, char **out_filename,
                         int *auxfiles, int *retval, char **set_name)

{
    gs_lambda_var_01(gs_fdta(*set_name), lambda, n, out_filename, auxfiles,
                     retval);
}

void gs_build_V_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                   int *retval)
{

    gs_UNSIGNED i, j, idx_A;
    FILE *outfp;

    if (NULL == gv->Z)
    {
        int Srv;
        int *Sretval = &Srv;
        gs_build_Z_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }

    gv->V = realloc(gv->V, gv->NoInd * gv->NoInd * sizeof(double));
    if (NULL == gv->V)
        ERR_M;

    st_dgemm("N", "T", gv->NoInd, gv->NoMar, gv->NoInd, gv->NoMar, gv->Z, gv->Z,
             gv->V);

    if (1 == *auxfiles)
    {
        if (NULL == (outfp = fopen(*out_filename, "w")))
            ERR_F;

        for (i = 0; i < gv->NoInd; i++)
        {
            for (j = 0; j < gv->NoInd; j++)
            {
                idx_A = i * gv->NoInd + j;
                fprintf(outfp, "%f ", gv->V[idx_A]);
            }
            fprintf(outfp, "\n");
        }

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_build_V_01_GV(char **out_filename, int *auxfiles, int *retval,
                      char **set_name)
{
    gs_build_V_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

static int gs_restrict_read_name_file(gs_varset_type *gv, const char *filename,
                                      const gs_name_index_type *index,
                                      gs_UNSIGNED n, unsigned char *keep)
{
    FILE *fp = NULL;
    char *line = NULL;
    size_t line_size = 0;
    long long int rc;

    if ((NULL == filename) || ('\0' == filename[0]))
        return 0;

    memset(keep, 0, (size_t)n * sizeof(*keep));

    fp = fopen(filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        return -1;
    }

    while (0 < (rc = st_getline(&line, &line_size, fp)))
    {
        char *p;

        gs_remove_newline(line);
        p = strtok(line, gs_tok_del);

        while (NULL != p)
        {
            gs_UNSIGNED idx = gs_name_index_find(index, n, p);

            if (gs_UNSIGNED_MAX != idx)
                keep[idx] = 1;

            p = strtok(NULL, gs_tok_del);
        }
    }

    fclose(fp);
    free(line);

    if (rc < 0)
    {
        gs_info(gv, -2, "Memory allocation error");
        return -1;
    }

    return 0;
}

static int gs_restrict_marker_effects_valid(const gs_varset_type *gv)
{
    gs_UNSIGNED m;

    if ((gs_UNSIGNED_MAX == gv->NoMar) || (gv->NoEff != gv->NoMar + 1) ||
        (gv->dimu != gv->NoEff) || (NULL == gv->EffNme) ||
        (NULL == gv->EffAll) || (NULL == gv->CmpAll) || (NULL == gv->u))
        return 0;

    if ((NULL == memchr(gv->EffNme[0], '\0', gs_MAXNME)) ||
        (0 != strcmp(gv->EffNme[0], "mu")))
        return 0;

    for (m = 0; m < gv->NoMar; m++)
    {
        if ((NULL == memchr(gv->EffNme[m + 1], '\0', gs_MAXNME)) ||
            (0 != strcmp(gv->EffNme[m + 1], gv->MarName[m])))
            return 0;
    }

    return 1;
}

void gs_restrict_marker_data_03(gs_varset_type *gv, char **ind_fname,
                                char **mar_fname, int *NoAllMAX,
                                double *MaMisMAX, double *ExHetMIN,
                                double *InMisMAX, int *retval)
{
    unsigned char *keep_i = NULL;
    unsigned char *keep_m = NULL;
    gs_UNSIGNED *old_to_new_m = NULL;
    gs_name_index_type *ind_index = NULL;
    gs_name_index_type *mar_index = NULL;
    gs_UNSIGNED *map_storage = NULL;
    gs_varset_type tmp;
    gs_UNSIGNED CountInd = 0;
    gs_UNSIGNED CountMar = 0;
    int has_map;
    int no_marker_reduction;
    int effects_valid = 0;
    int failed = 0;

    *retval = -2;

    if ((0 == gv->NoInd) || (0 == gv->NoMar) || (NULL == gv->IndName) ||
        (NULL == gv->MarName) || !st_allele_present(gv))
    {
        gs_info(gv, -2,
                "No marker data available, run 'read.marker.data' first");
        return;
    }

    if ((NULL == gv->NoAll) || (NULL == gv->All) || (NULL == gv->Mis) ||
        (NULL == gv->Pic) || (NULL == gv->IndMis))
    {
        gs_info(gv, -2, "Marker statistics are not available");
        return;
    }

    if ((*NoAllMAX < 0) || !isfinite(*MaMisMAX) || !isfinite(*ExHetMIN) ||
        !isfinite(*InMisMAX))
    {
        gs_info(gv, -2, "Invalid marker restriction criterion");
        return;
    }

    keep_i = calloc((size_t)gv->NoInd, sizeof(*keep_i));
    keep_m = calloc((size_t)gv->NoMar, sizeof(*keep_m));

    if ((NULL == keep_i) || (NULL == keep_m))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    memset(keep_i, 1, (size_t)gv->NoInd * sizeof(*keep_i));
    memset(keep_m, 1, (size_t)gv->NoMar * sizeof(*keep_m));

    if ('\0' != (*ind_fname)[0])
    {
        if (0 != gs_make_ind_name_index(gv->IndName, gv->NoInd, &ind_index))
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }

        if (0 != gs_restrict_read_name_file(gv, *ind_fname, ind_index,
                                            gv->NoInd, keep_i))
        {
            failed = 1;
            goto cleanup;
        }
    }

    if ('\0' != (*mar_fname)[0])
    {
        if (0 != gs_make_name_index(gv->MarName, gv->NoMar, &mar_index))
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }

        if (0 != gs_restrict_read_name_file(gv, *mar_fname, mar_index,
                                            gv->NoMar, keep_m))
        {
            failed = 1;
            goto cleanup;
        }
    }

    for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
    {
        if ((*InMisMAX < gv->IndMis[i]) || !keep_i[i])
            keep_i[i] = 0;
        else
            CountInd++;
    }

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
        gs_UNSIGNED no_all = gv->NoAll[m];

        if ((0 < no_all) && (GS_ALLELE_MISSING == gv->All[m][no_all - 1]))
            no_all--;

        if (((gs_UNSIGNED)(*NoAllMAX) < no_all) || (*MaMisMAX < gv->Mis[m]) ||
            (*ExHetMIN > gv->Pic[m]) || !keep_m[m])
            keep_m[m] = 0;
        else
            CountMar++;
    }

    if (0 == CountInd)
    {
        gs_info(gv, -2, "No individuals remain after marker-data restriction");
        failed = 1;
        goto cleanup;
    }

    if (0 == CountMar)
    {
        gs_info(gv, -2, "No markers remain after marker-data restriction");
        failed = 1;
        goto cleanup;
    }

    if ((CountInd == gv->NoInd) && (CountMar == gv->NoMar))
    {
        *retval = 0;
        goto cleanup;
    }

    no_marker_reduction = (int)(CountMar == gv->NoMar);
    effects_valid = no_marker_reduction && gs_restrict_marker_effects_valid(gv);

    has_map =
        ((0 != gv->NLoc) || (NULL != gv->map_S) || (NULL != gv->map_str_S) ||
         (0 != gv->map_str_N) || (NULL != gv->map) || (NULL != gv->mdl));

    if (has_map)
    {
        map_storage = gs_map_storage_index(gv, 0, 0);
        if (NULL == map_storage)
        {
            failed = 1;
            goto cleanup;
        }
    }

    gs_tmp_init_from(&tmp, gv);
    tmp.NoInd = CountInd;
    tmp.NoMar = CountMar;
    tmp.phased = gv->phased;

    tmp.IndName = calloc((size_t)tmp.NoInd, sizeof(*tmp.IndName));
    tmp.MarName = gs_malloc_array((size_t)tmp.NoMar, sizeof(*tmp.MarName));

    if ((NULL == tmp.IndName) || (NULL == tmp.MarName) ||
        (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup_tmp;
    }

    old_to_new_m = gs_malloc_array((size_t)gv->NoMar, sizeof(*old_to_new_m));
    if (NULL == old_to_new_m)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup_tmp;
    }

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
        old_to_new_m[m] = gs_UNSIGNED_MAX;

    {
        gs_UNSIGNED i_new = 0;

        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            if (keep_i[i])
            {
                size_t len = strlen(gv->IndName[i]);
                if (st_identity_copy(&tmp, i_new, gv, i))
                {
                    failed = 1;
                    goto cleanup_tmp;
                }
                if ((int)len > tmp.MaxIL)
                    tmp.MaxIL = (int)len;
                i_new++;
            }
    }

    {
        gs_UNSIGNED m_new = 0;

        for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            if (keep_m[m])
            {
                gs_UNSIGNED i_new = 0;
                size_t len = strlen(gv->MarName[m]);

                strcpy(tmp.MarName[m_new], gv->MarName[m]);
                if ((int)len > tmp.MaxML)
                    tmp.MaxML = (int)len;
                old_to_new_m[m] = m_new;

                for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
                    if (keep_i[i])
                    {
                        size_t old_idx = (size_t)m * gv->NoInd + i;
                        size_t new_idx = (size_t)m_new * tmp.NoInd + i_new;

                        tmp.allele_raw1[new_idx] =
                            st_allele_get(gv, old_idx, 0);
                        tmp.allele_raw2[new_idx] =
                            st_allele_get(gv, old_idx, 1);
                        i_new++;
                    }

                m_new++;
            }
    }

    if (NULL != gv->y)
    {
        gs_UNSIGNED i_new = 0;

        tmp.y = gs_malloc_array((size_t)tmp.NoInd, sizeof(*tmp.y));
        if (NULL == tmp.y)
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup_tmp;
        }

        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            if (keep_i[i])
                tmp.y[i_new++] = gv->y[i];
    }

    if (has_map)
    {
        gs_UNSIGNED nloc = 0;

        for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
            if (keep_m[gv->mdl[lne]])
                nloc++;

        if (0 < nloc)
        {
            size_t map_str_offset = 0;

            for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
                if (keep_m[gv->mdl[lne]])
                    if ((0 != gs_pool_add_string_size(&(tmp.map_str_N),
                                                      gv->map[lne]->name)) ||
                        (0 != gs_pool_add_string_size(&(tmp.map_str_N),
                                                      gv->map[lne]->class)))
                    {
                        gs_info(gv, -2, "Linkage map string storage too large");
                        failed = 1;
                        goto cleanup_tmp;
                    }

            tmp.map_S = gs_malloc_array((size_t)nloc, sizeof(*tmp.map_S));
            tmp.map_str_S =
                gs_malloc_array(tmp.map_str_N, sizeof(*tmp.map_str_S));
            tmp.map = gs_malloc_array((size_t)nloc, sizeof(*tmp.map));
            tmp.mdl = gs_malloc_array((size_t)nloc, sizeof(*tmp.mdl));

            if ((NULL == tmp.map_S) || (NULL == tmp.map_str_S) ||
                (NULL == tmp.map) || (NULL == tmp.mdl))
            {
                gs_info(gv, -2, "Memory allocation error");
                failed = 1;
                goto cleanup_tmp;
            }

            tmp.NLoc = nloc;
            nloc = 0;

            for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
            {
                gs_UNSIGNED old_m = gv->mdl[lne];

                if (keep_m[old_m])
                {
                    gs_UNSIGNED new_m = old_to_new_m[old_m];

                    if (0 != gs_copy_map_point_to_pool(
                                 &(tmp.map_S[nloc]), gv->map[lne], new_m,
                                 tmp.map_str_S, tmp.map_str_N, &map_str_offset))
                    {
                        gs_info(gv, -2, "Invalid linkage map string storage");
                        failed = 1;
                        goto cleanup_tmp;
                    }
                    tmp.map[nloc] = &(tmp.map_S[nloc]);
                    tmp.mdl[nloc] = new_m;
                    nloc++;
                }
            }

            if (map_str_offset != tmp.map_str_N)
            {
                gs_info(gv, -2, "Invalid linkage map string storage");
                failed = 1;
                goto cleanup_tmp;
            }
        }
    }

    {
        int srv = 0;
        char *empty = "";
        char **empty_p = &empty;
        int zero = 0;
        int *zero_p = &zero;

        gs_marker_stats_01(&tmp, empty_p, empty_p, empty_p, zero_p, zero_p,
                           zero_p, zero_p, &srv);
        if (-2 == srv)
        {
            failed = 1;
            goto cleanup_tmp;
        }
    }

    if (effects_valid)
    {
        tmp.EffNme = gs_malloc_array((size_t)gv->NoEff, sizeof(*tmp.EffNme));
        tmp.EffAll = gs_malloc_array((size_t)gv->NoEff, sizeof(*tmp.EffAll));
        tmp.CmpAll = gs_malloc_array((size_t)gv->NoEff, sizeof(*tmp.CmpAll));
        tmp.u = gs_malloc_array((size_t)gv->dimu, sizeof(*tmp.u));

        if ((NULL == tmp.EffNme) || (NULL == tmp.EffAll) ||
            (NULL == tmp.CmpAll) || (NULL == tmp.u))
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup_tmp;
        }

        tmp.NoEff = gv->NoEff;
        tmp.dimu = gv->dimu;
        memcpy(tmp.EffNme, gv->EffNme, (size_t)gv->NoEff * sizeof(*tmp.EffNme));
        memcpy(tmp.EffAll, gv->EffAll, (size_t)gv->NoEff * sizeof(*tmp.EffAll));
        memcpy(tmp.CmpAll, gv->CmpAll, (size_t)gv->NoEff * sizeof(*tmp.CmpAll));
        memcpy(tmp.u, gv->u, (size_t)gv->dimu * sizeof(*tmp.u));

        tmp.aa = gv->aa;
        tmp.aA = gv->aA;
        tmp.AA = gv->AA;
        tmp.an = gv->an;
        tmp.nn = gv->nn;
        tmp.nA = gv->nA;

        gs_info(gv, 1, "Effects retained");
    }
    else if (!no_marker_reduction &&
             ((NULL != gv->u) || (NULL != gv->EffNme) || (NULL != gv->EffAll) ||
              (NULL != gv->CmpAll)))
    {
        gs_info(gv, 1, "Effects removed");
    }

    tmp.physical_map = gv->physical_map;
    st_data_retain(tmp.physical_map);
    if (gs_tmp_commit(gv, &tmp))
    {
        failed = 1;
        goto cleanup_tmp;
    }
    *retval = 0;

cleanup_tmp:
    if (failed)
        gs_reset(&tmp);

cleanup:
    free(keep_i);
    free(keep_m);
    free(old_to_new_m);
    free(ind_index);
    free(mar_index);
    free(map_storage);
}

void gs_restrict_marker_data_03_GV(char **ind_fname, char **mar_fname,
                                   int *NoAllMAX, double *MaMisMAX,
                                   double *ExHetMIN, double *InMisMAX,
                                   int *retval, char **set_name)
{
    gs_restrict_marker_data_03(gs_fdta(*set_name), ind_fname, mar_fname,
                               NoAllMAX, MaMisMAX, ExHetMIN, InMisMAX, retval);
}

static void gs_copy_marker_data_01_raw(gs_varset_type *new, gs_varset_type *gv,
                                       int *retval)
{
    gs_UNSIGNED i, m;

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2,
                "No marker data available, run 'read.marker.data' first");
        *retval = -2;
        return;
    }

    if (new == gv)
    {
        *retval = 0;
        return;
    }

    gs_reset(new);

    new->NoInd = gv->NoInd;
    new->NoMar = gv->NoMar;
    new->phased = gv->phased;

    new->IndName = calloc(new->NoInd, sizeof(*new->IndName));
    new->MarName = realloc(new->MarName, new->NoMar * sizeof(char[gs_MAXNME]));
    new->allele_raw1 = realloc(new->allele_raw1, new->NoMar *new->NoInd *
                                                     sizeof(*new->allele_raw1));
    new->allele_raw2 = realloc(new->allele_raw2, new->NoMar *new->NoInd *
                                                     sizeof(*new->allele_raw2));

    if ((NULL == new->allele_raw1) || (NULL == new->allele_raw2) ||
        (NULL == new->IndName) || (NULL == new->MarName))
        ERR_M;

    for (i = 0; i < gv->NoInd; i++)
        if (st_identity_copy(new, i, gv, i))
        {
            ERR_M;
        }

    for (m = 0; m < gv->NoMar; m++)
        strcpy(new->MarName[m], gv->MarName[m]);

    for (m = 0; m < gv->NoMar; m++)
    {
        st_allele_copy(gv, 0, m * (size_t)gv->NoInd, new->NoInd,
                       new->allele_raw1 + m *(size_t)new->NoInd);
        st_allele_copy(gv, 1, m * (size_t)gv->NoInd, new->NoInd,
                       new->allele_raw2 + m *(size_t)new->NoInd);
    }

    int Srv = 0;
    int *Sretval = &Srv;
    char *empty = "";
    char **ind_filename = &empty;
    char **mar_filename = &empty;
    char **gen_filename = &empty;
    int auf = 0;
    int *Sauxfiles = &auf;
    gs_marker_stats_01(new, ind_filename, mar_filename, gen_filename, Sauxfiles,
                       Sauxfiles, Sauxfiles, Sauxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    {
        int map_rc = gs_copy_map_primary(new, gv);

        if (-1 == map_rc)
            ERR_M;
        if (-2 == map_rc)
        {
            *retval = -2;
            return;
        }
    }

    if (NULL != gv->y)
    {
        new->y = realloc(new->y, new->NoInd * sizeof(double));
        if (NULL == new->y)
            ERR_M;
        memcpy(new->y, gv->y, (size_t)new->NoInd * sizeof(double));
    }

    st_marker_prediction_copy(new, gv);
    *retval = 0;
    return;
}

void gs_copy_marker_data_01(gs_varset_type *dst, gs_varset_type *src,
                            int *retval)
{
    if (dst == src)
    {
        *retval = 0;
        return;
    }
    gs_varset_type tmp;
    gs_tmp_init_from(&tmp, dst);
    gs_copy_marker_data_01_raw(&tmp, src, retval);
    if (*retval >= 0 && gs_tmp_commit(dst, &tmp))
        *retval = -2;
    gs_reset(&tmp);
}

void gs_copy_marker_data_01_GV(int *retval, char **set_name_t,
                               char **set_name_s)
{
    gs_copy_marker_data_01(gs_fdta(*set_name_t), gs_fdta(*set_name_s), retval);
}

static void st_combine_destroy(void *pointer)
{
    gs_varset_type *tmp = pointer;
    gs_reset(tmp);
    free(tmp);
}

static const char *st_combine_source_name(SEXP sources, R_xlen_t i)
{
    SEXP value;
    if (TYPEOF(sources) == VECSXP)
    {
        SEXP item = VECTOR_ELT(sources, i);
        if (TYPEOF(item) != STRSXP || XLENGTH(item) != 1)
            error("Each source.data.sets entry must be one dataset name");
        value = STRING_ELT(item, 0);
    }
    else
        value = STRING_ELT(sources, i);
    if (value == NA_STRING)
        error("source.data.sets must not contain missing dataset names");
    const char *name = translateCharUTF8(value);
    if (!*name || strlen(name) >= gs_MAXNME)
        error("Source dataset names must be nonempty and shorter than %d "
              "bytes",
              gs_MAXNME);
    return name;
}

SEXP st_data_combine(SEXP target, SEXP sources)
{
    if (TYPEOF(target) != STRSXP || XLENGTH(target) != 1 ||
        STRING_ELT(target, 0) == NA_STRING)
        error("target.data.set must be one nonempty dataset name");
    const char *name = translateCharUTF8(STRING_ELT(target, 0));
    if (!*name || strlen(name) >= gs_MAXNME)
        error("target.data.set must be nonempty and shorter than %d bytes",
              gs_MAXNME);
    if ((TYPEOF(sources) != STRSXP && TYPEOF(sources) != VECSXP) ||
        XLENGTH(sources) < 1)
        error("source.data.sets must be a nonempty list or character vector "
              "of dataset names");
    if (st_hybrid_name_exists(name))
        error("Target '%s' is a hybrid data set", name);
    gs_varset_type *destination = st_find_data_set_no_create(name);
    if (destination)
        st_hybrid_require_unreferenced(destination);

    size_t nsource = (size_t)XLENGTH(sources);
    gs_varset_type **input = st_work_alloc(nsource, sizeof(*input));
    gs_UNSIGNED *offset = st_work_alloc(nsource, sizeof(*offset));
    gs_UNSIGNED individuals = 0, markers = 0;
    int phenotype = 0, phased = 1;
    for (size_t s = 0; s < nsource; s++)
    {
        const char *source = st_combine_source_name(sources, (R_xlen_t)s);
        input[s] = st_find_data_set_no_create(source);
        if (!input[s])
            error("Unknown source marker data set '%s'", source);
        if (input[s] == destination)
            error("The target must differ from every source data set");
        for (size_t previous = 0; previous < s; previous++)
            if (input[s] == input[previous])
                error("Source data set '%s' is listed more than once", source);
        gs_varset_type *x = input[s];
        if (!x->NoInd || !x->NoMar || !x->IndName || !x->MarName ||
            !gs_marker_genotypes_valid(x))
            error("Source '%s' has no valid marker data", source);
        if (!s)
            markers = x->NoMar;
        else if (x->NoMar != markers)
            error("Source '%s' does not have exactly the same marker set "
                  "as '%s'",
                  source, input[0]->name);
        if (x->NoInd > gs_UNSIGNED_MAX - individuals)
            error("Too many combined individuals");
        offset[s] = individuals;
        individuals += x->NoInd;
        phenotype |= x->y != NULL;
        phased &= x->phased != 0;
    }

    size_t mappings, cells;
    if (gs_size_mul(nsource, (size_t)markers, &mappings) ||
        gs_size_mul((size_t)individuals, (size_t)markers, &cells))
        error("Combined marker dimensions exceed native address limits");
    gs_UNSIGNED *order = st_work_alloc(mappings, sizeof(*order));
    gs_name_index_type *index = st_work_alloc(markers, sizeof(*index));
    for (size_t s = 0; s < nsource; s++)
    {
        for (gs_UNSIGNED m = 0; m < markers; m++)
        {
            if (!input[s]->MarName[m][0])
                error("Source '%s' has an empty marker name", input[s]->name);
            index[m].name = input[s]->MarName[m];
            index[m].idx = m;
        }
        qsort(index, markers, sizeof(*index), gs_name_index_cmp);
        if (gs_name_index_has_duplicates(index, markers))
            error("Source '%s' has duplicate marker names", input[s]->name);
        for (gs_UNSIGNED m = 0; m < markers; m++)
        {
            gs_UNSIGNED at =
                gs_name_index_find(index, markers, input[0]->MarName[m]);
            if (at == gs_UNSIGNED_MAX)
                error("Source '%s' does not have exactly the same marker "
                      "set: marker '%s' is missing",
                      input[s]->name, input[0]->MarName[m]);
            order[s * (size_t)markers + m] = at;
        }
    }

    gs_name_index_type *identities =
        st_work_alloc(individuals, sizeof(*identities));
    for (size_t s = 0; s < nsource; s++)
        for (gs_UNSIGNED i = 0; i < input[s]->NoInd; i++)
        {
            const char *id = input[s]->IndName[i];
            if (!id || !*id)
                error("Source '%s' has an invalid individual name",
                      input[s]->name);
            identities[offset[s] + i].name = id;
            identities[offset[s] + i].idx = offset[s] + i;
        }
    qsort(identities, individuals, sizeof(*identities), gs_name_index_cmp);
    for (gs_UNSIGNED i = 1; i < individuals; i++)
        if (!strcmp(identities[i - 1].name, identities[i].name))
            error("Duplicate individual name '%s' in source data sets",
                  identities[i].name);

    gs_varset_type *tmp = malloc(sizeof(*tmp));
    if (!tmp)
        error("Cannot allocate combined marker data set");
    gs_init(tmp);
    st_work_adopt(tmp, st_combine_destroy);
    snprintf(tmp->name, sizeof(tmp->name), "%s", name);
    tmp->info_level = destination ? destination->info_level : gv->info_level;
    tmp->NoInd = individuals;
    tmp->NoMar = markers;
    tmp->phased = phased;
    tmp->MaxML = input[0]->MaxML;
    tmp->IndName = calloc((size_t)individuals, sizeof(*tmp->IndName));
    tmp->MarName = gs_malloc_array((size_t)markers, sizeof(*tmp->MarName));
    if (!tmp->IndName || !tmp->MarName ||
        gs_alloc_genotype_matrix(tmp, markers, individuals))
        error("Cannot allocate combined marker matrix");
    memcpy(tmp->MarName, input[0]->MarName,
           (size_t)markers * sizeof(*tmp->MarName));
    if (phenotype)
    {
        tmp->y = gs_malloc_array((size_t)individuals, sizeof(*tmp->y));
        if (!tmp->y)
            error("Cannot allocate combined phenotypes");
    }
    for (size_t s = 0; s < nsource; s++)
        for (gs_UNSIGNED i = 0; i < input[s]->NoInd; i++)
        {
            gs_UNSIGNED to = offset[s] + i;
            if (st_identity_copy(tmp, to, input[s], i))
                error("Cannot copy combined individual identities");
            size_t length = strlen(tmp->IndName[to]);
            if ((int)length > tmp->MaxIL)
                tmp->MaxIL = (int)length;
            if (phenotype)
                tmp->y[to] = input[s]->y ? input[s]->y[i] : NA_REAL;
        }
    if (gs_copy_map_primary(tmp, input[0]))
        error("Cannot copy the first source's linkage map");

    R_CheckUserInterrupt();
#ifdef MT
    int threads = st_num_threads();
#pragma omp parallel for num_threads(threads)                                  \
    schedule(static) if (cells >= 1000000)
#endif
    for (gs_UNSIGNED m = 0; m < markers; m++)
        for (size_t s = 0; s < nsource; s++)
        {
            size_t from =
                (size_t)order[s * (size_t)markers + m] * input[s]->NoInd;
            size_t to = (size_t)m * individuals + offset[s];
            st_allele_copy(input[s], 0, from, input[s]->NoInd,
                           tmp->allele_raw1 + to);
            st_allele_copy(input[s], 1, from, input[s]->NoInd,
                           tmp->allele_raw2 + to);
        }
    R_CheckUserInterrupt();
    int status = 0, zero = 0;
    char *empty = "";
    gs_marker_stats_01(tmp, &empty, &empty, &empty, &zero, &zero, &zero, &zero,
                       &status);
    if (status == -2 || gs_prepare_allele_commit(tmp))
        error("Cannot prepare combined marker statistics or packed storage");

    SEXP result = PROTECT(ScalarInteger(0));
    if (!destination)
    {
        if (gs_ndta == gs_UNSIGNED_MAX)
            error("Too many named marker data sets");
        destination = malloc(sizeof(*destination));
        if (!destination)
            error("Cannot allocate the target marker data set");
        gs_varset_type **registry =
            gs_realloc_array(gs_data, (size_t)gs_ndta + 1, sizeof(*registry));
        if (!registry)
        {
            free(destination);
            error("Cannot register the target marker data set");
        }
        gs_init(destination);
        gs_data = registry;
        gs_data[gs_ndta++] = destination;
    }
    gs_tmp_install(destination, tmp);
    for (gs_UNSIGNED i = 0; i < destination->NoInd; i++)
        ps_identity_observe(destination->IndName[i]);
    UNPROTECT(1);
    return result;
}

void gs_estimate_bv_01(gs_varset_type *effect_set, gs_varset_type *marker_set,
                       char **out_filename, int *auxfiles, int *retval)
{
    if (NULL == effect_set->u)
    {
        gs_info(effect_set, -2, "No estimated effects available");
        *retval = -2;
        return;
    }

    if ((0 == marker_set->NoInd) || (0 == marker_set->NoMar))
    {
        gs_info(marker_set, -2, "No marker data available");
        *retval = -2;
        return;
    }
    gs_UNSIGNED i, m, idx_all, idx_Z;

    marker_set->est =
        realloc(marker_set->est, marker_set->NoInd * sizeof(double));
    if (NULL == marker_set->est)
        ERR_M;

    int by_matrix_multiplication = 1;

    if (marker_set->NoMar != effect_set->NoMar)
    {
        by_matrix_multiplication = 0;
    }
    else
    {
        for (m = 0; m < marker_set->NoMar; m++)
        {
            if (0 != strcmp(marker_set->MarName[m], effect_set->MarName[m]))
            {
                by_matrix_multiplication = 0;
                break;
            }
        }
    }

    if (1 == by_matrix_multiplication)
    {
        gs_info(effect_set, 1, "Marker sets are identical");

        if (NULL == marker_set->Z)
        {

            gs_info(effect_set, 1, "Constructing design matrix");

            marker_set->Z =
                realloc(marker_set->Z,
                        marker_set->NoInd * marker_set->NoMar * sizeof(double));
            int a, b, lar_all, sml_all;
            for (m = 0; m < marker_set->NoMar; m++)
            {

                lar_all = effect_set->EffAll[m + 1];
                sml_all = effect_set->CmpAll[m + 1];

                for (i = 0; i < marker_set->NoInd; i++)
                {
                    idx_all = m * marker_set->NoInd + i;

                    a = st_allele_get(marker_set, idx_all, 0);
                    b = st_allele_get(marker_set, idx_all, 1);

                    idx_Z = i * marker_set->NoMar + m;

                    if ((a == lar_all) && (b == lar_all))
                    {
                        marker_set->Z[idx_Z] = marker_set->AA;
                    }
                    else
                    {
                        if ((a == sml_all) && (b == sml_all))
                        {
                            marker_set->Z[idx_Z] = marker_set->aa;
                        }
                        else
                        {
                            if (((a == lar_all) && (b == sml_all)) ||
                                ((a == sml_all) && (b == lar_all)))
                            {
                                marker_set->Z[idx_Z] = marker_set->aA;
                            }
                            else
                            {
                                if ((a == lar_all) || (b == lar_all))
                                {
                                    marker_set->Z[idx_Z] = marker_set->nA;
                                }
                                else
                                {
                                    if ((a == sml_all) || (b == sml_all))
                                    {
                                        marker_set->Z[idx_Z] = marker_set->an;
                                    }
                                    else
                                    {
                                        marker_set->Z[idx_Z] = marker_set->nn;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        st_dgemm("N", "N", marker_set->NoInd, marker_set->NoMar,
                 marker_set->NoMar, 1, marker_set->Z, &(effect_set->u[1]),
                 marker_set->est);

        for (i = 0; i < marker_set->NoInd; i++)
            marker_set->est[i] += effect_set->u[0];
    }
    else
    {
        gs_info(effect_set, 1, "Different marker sets");

#ifdef MT
#pragma omp parallel
#endif
        {
            gs_UNSIGNED i, j, m, idx_all;
            int a, b, lar_all, sml_all;
            double count;

#ifdef MT
#pragma omp for
#endif
            for (i = 0; i < marker_set->NoInd; i++)
            {
                marker_set->est[i] = effect_set->u[0];
                for (j = 1; j < effect_set->dimu; j++)
                {
                    int found = 0;
                    for (m = 0; m < marker_set->NoMar; m++)
                    {
                        if (0 == strcmp(marker_set->MarName[m],
                                        effect_set->EffNme[j]))
                        {
                            lar_all = effect_set->EffAll[j];
                            sml_all = effect_set->CmpAll[j];

                            idx_all = m * marker_set->NoInd + i;
                            a = st_allele_get(marker_set, idx_all, 0);
                            b = st_allele_get(marker_set, idx_all, 1);

                            if ((a == lar_all) && (b == lar_all))
                            {
                                count = effect_set->AA;
                            }
                            else
                            {
                                if ((a == sml_all) && (b == sml_all))
                                {
                                    count = effect_set->aa;
                                }
                                else
                                {
                                    if (((a == lar_all) && (b == sml_all)) ||
                                        ((a == sml_all) && (b == lar_all)))
                                    {
                                        count = effect_set->aA;
                                    }
                                    else
                                    {
                                        if ((a == lar_all) || (b == lar_all))
                                        {
                                            count = effect_set->nA;
                                        }
                                        else
                                        {
                                            if ((a == sml_all) ||
                                                (b == sml_all))
                                            {
                                                count = effect_set->an;
                                            }
                                            else
                                            {
                                                count = effect_set->nn;
                                            }
                                        }
                                    }
                                }
                            }

                            marker_set->est[i] += count * effect_set->u[j];
                            found = 1;
                            break;
                        }
                    }
                    if (0 == found)
                    {
                        count = effect_set->nn;
                        marker_set->est[i] += count * effect_set->u[j];
                    }
                }
            }
        }
    }

    if (1 == *auxfiles)
    {
        gs_UNSIGNED i;
        FILE *fp;
        if ((fp = fopen(*out_filename, "w")) == NULL)
            ERR_F;
        if (NULL == marker_set->y)
        {
            fprintf(fp, "yhat \n");
            for (i = 0; i < marker_set->NoInd; i++)
                fprintf(fp, "%s %f\n", marker_set->IndName[i],
                        marker_set->est[i]);
        }
        else
        {
            fprintf(fp, "y yhat \n");
            for (i = 0; i < marker_set->NoInd; i++)
                fprintf(fp, "%s %f %f\n", marker_set->IndName[i],
                        marker_set->y[i], marker_set->est[i]);
        }

        fclose(fp);
    }

    *retval = 0;
    return;
}

void gs_estimate_bv_01_GV(char **effect_set, char **marker_set,
                          char **out_filename, int *auxfiles, int *retval)
{
    gs_estimate_bv_01(gs_fdta(*effect_set), gs_fdta(*marker_set), out_filename,
                      auxfiles, retval);
}

void gs_check_model_fit_01_GV(char **out_filename, int *auxfiles, int *retval)
{
    gs_estimate_bv_01(gv, gv, out_filename, auxfiles, retval);
}

void gs_estimate_breeding_values_01_GV(char **in_marker, char **out_filename,
                                       int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_varset_type VS;
    gs_varset_type *valset = &VS;
    gs_init(valset);

    gs_read_marker_data_01(valset, in_marker, nixp, nixp, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        gs_reset(valset);
        *retval = -2;
        return;
    }

    gs_estimate_bv_01(gv, valset, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        gs_reset(valset);
        *retval = -2;
        return;
    }

    gs_reset(valset);

    *retval = 0;
    return;
}

void gs_validation_01_GV(char **in_marker, char **in_performance,
                         char **out_filename, int *auxfiles, int *retval)
{

    int Srv;
    int *Sretval = &Srv;

    gs_varset_type VS;
    gs_varset_type *valset = &VS;
    gs_init(valset);

    gs_read_marker_data_01(valset, in_marker, nixp, nixp, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        gs_reset(valset);
        *retval = -2;
        return;
    }

    gs_read_performance_data_01(valset, in_performance, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        gs_reset(valset);
        *retval = -2;
        return;
    }

    gs_estimate_bv_01(gv, valset, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        gs_reset(valset);
        *retval = -2;
        return;
    }

    if (1 == *auxfiles)
    {
        FILE *fp;
        if (NULL == (fp = fopen(*out_filename, "w")))
        {
            gs_reset(valset);
            ERR_F;
        }
        fprintf(fp, "y yhat \n");
        for (gs_UNSIGNED i = 0; i < valset->NoInd; i++)
            fprintf(fp, "%s %f %f\n", valset->IndName[i], valset->y[i],
                    valset->est[i]);
        fclose(fp);
    }

    gs_reset(valset);

    *retval = 0;
    return;
}

static void gs_build_esvs_01_raw(gs_varset_type *gv, gs_varset_type *es,
                                 gs_varset_type *vs, int *n_estimation,
                                 char **out_es_m, char **out_es_p,
                                 char **out_vs_m, char **out_vs_p,
                                 int *auxfiles, int *retval)
{

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2,
                "No marker data available, run 'read.marker.data' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->y)
    {
        gs_info(
            gv, -2,
            "No performance data available, run 'read.performance.data' first");
        *retval = -2;
        return;
    }

    gs_reset(es);
    gs_reset(vs);

    gs_UNSIGNED CVSiEst = (gs_UNSIGNED)*n_estimation;

    gs_UNSIGNED i, j, k, m, idx_g, idx_e, idx_v, idx;

    gs_UNSIGNED *CVes = malloc(gv->NoInd * sizeof(gs_UNSIGNED));
    if (NULL == CVes)
        ERR_M;

    int *p = malloc(gv->NoInd * sizeof(int));
    if (NULL == p)
    {
        free(CVes);
        ERR_M;
    }

    p[0] = 0;
    for (int i = 1; i < (int)gv->NoInd; ++i)
    {
        int j = st_rand_index(i + 1);
        p[i] = p[j];
        p[j] = i;
    }

    for (i = 0; i < gv->NoInd; i++)
        CVes[i] = 0;
    for (i = 0; i < CVSiEst; i++)
        CVes[p[i]] = 1;
    free(p);

    es->NoInd = CVSiEst;
    vs->NoInd = gv->NoInd - es->NoInd;

    es->IndName = calloc(es->NoInd, sizeof(*es->IndName));
    vs->IndName = calloc(vs->NoInd, sizeof(*vs->IndName));
    if ((NULL == es->IndName) || (NULL == vs->IndName))
    {
        free(CVes);
        ERR_M;
    }

    for (i = 0, j = 0, k = 0; i < gv->NoInd; i++)
    {
        if (1 == CVes[i])
        {
            if (st_identity_copy(es, j, gv, i))
            {
                free(CVes);
                ERR_M;
            }
            j++;
        }
        else
        {
            if (st_identity_copy(vs, k, gv, i))
            {
                free(CVes);
                ERR_M;
            }
            k++;
        }
    }

    es->NoMar = gv->NoMar;
    vs->NoMar = gv->NoMar;

    es->MarName = realloc(es->MarName, es->NoMar * sizeof(char[gs_MAXNME]));
    vs->MarName = realloc(vs->MarName, vs->NoMar * sizeof(char[gs_MAXNME]));
    if ((NULL == es->MarName) || (NULL == vs->MarName))
    {
        free(CVes);
        ERR_M;
    }

    for (m = 0; m < gv->NoMar; m++)
    {
        strcpy(es->MarName[m], gv->MarName[m]);
        strcpy(vs->MarName[m], gv->MarName[m]);
    }

    es->allele_raw1 = realloc(es->allele_raw1,
                              es->NoMar * es->NoInd * sizeof(*es->allele_raw1));
    es->allele_raw2 = realloc(es->allele_raw2,
                              es->NoMar * es->NoInd * sizeof(*es->allele_raw1));
    if ((NULL == es->allele_raw1) || (NULL == es->allele_raw2))
    {
        free(CVes);
        ERR_M;
    }

    vs->allele_raw1 = realloc(vs->allele_raw1,
                              vs->NoMar * vs->NoInd * sizeof(*vs->allele_raw1));
    vs->allele_raw2 = realloc(vs->allele_raw2,
                              vs->NoMar * vs->NoInd * sizeof(*vs->allele_raw1));
    if ((NULL == vs->allele_raw1) || (NULL == vs->allele_raw2))
    {
        free(CVes);
        ERR_M;
    }

    for (m = 0; m < gv->NoMar; m++)
    {
        for (i = 0, j = 0, k = 0; i < gv->NoInd; i++)
        {
            idx_g = m * gv->NoInd + i;
            idx_e = m * es->NoInd + j;
            idx_v = m * vs->NoInd + k;

            if (1 == CVes[i])
            {
                es->allele_raw1[idx_e] = st_allele_get(gv, idx_g, 0);
                es->allele_raw2[idx_e] = st_allele_get(gv, idx_g, 1);
                j++;
            }
            else
            {
                vs->allele_raw1[idx_v] = st_allele_get(gv, idx_g, 0);
                vs->allele_raw2[idx_v] = st_allele_get(gv, idx_g, 1);
                k++;
            }
        }
    }

    es->y = realloc(es->y, es->NoInd * sizeof(double));
    vs->y = realloc(vs->y, vs->NoInd * sizeof(double));
    if ((NULL == es->y) || (NULL == vs->y))
    {
        free(CVes);
        ERR_M;
    }

    for (i = 0, j = 0, k = 0; i < gv->NoInd; i++)
    {
        if (1 == CVes[i])
        {
            es->y[j] = gv->y[i];
            j++;
        }
        else
        {
            vs->y[k] = gv->y[i];
            k++;
        }
    }

    free(CVes);

    if (1 == *auxfiles)
    {

        FILE *fp1 = fopen(*out_es_m, "w");
        FILE *fp2 = fopen(*out_es_p, "w");
        FILE *fp3 = fopen(*out_vs_m, "w");
        FILE *fp4 = fopen(*out_vs_p, "w");
        if ((NULL == fp1) || (NULL == fp2) || (NULL == fp3) || (NULL == fp4))
            ERR_F;

        fprintf(fp1, "I;M;A\n");
        for (i = 0; i < es->NoInd; i++)
            for (m = 0; m < es->NoMar; m++)
            {
                idx = m * es->NoInd + i;
                if (GS_ALLELE_MISSING != st_allele_get(es, idx, 0))
                {
                    fprintf(fp1, "%s;%s;%i\n", es->IndName[i], es->MarName[m],
                            st_allele_get(es, idx, 0));
                    if (st_allele_get(es, idx, 0) != st_allele_get(es, idx, 1))
                        fprintf(fp1, "%s;%s;%i\n", es->IndName[i],
                                es->MarName[m], st_allele_get(es, idx, 1));
                }
            }

        fprintf(fp2, "I;P\n");
        for (i = 0; i < es->NoInd; i++)
            fprintf(fp2, "%s;%f\n", es->IndName[i], es->y[i]);

        fprintf(fp3, "I;M;A\n");
        for (i = 0; i < vs->NoInd; i++)
            for (m = 0; m < vs->NoMar; m++)
            {
                idx = m * vs->NoInd + i;
                if (GS_ALLELE_MISSING != st_allele_get(vs, idx, 0))
                {
                    fprintf(fp3, "%s;%s;%i\n", vs->IndName[i], vs->MarName[m],
                            st_allele_get(vs, idx, 0));
                    if (st_allele_get(vs, idx, 0) != st_allele_get(vs, idx, 1))
                        fprintf(fp3, "%s;%s;%i\n", vs->IndName[i],
                                vs->MarName[m], st_allele_get(vs, idx, 1));
                }
            }

        fprintf(fp4, "I;P\n");
        for (i = 0; i < vs->NoInd; i++)
            fprintf(fp4, "%s;%f\n", vs->IndName[i], vs->y[i]);

        fclose(fp1);
        fclose(fp2);
        fclose(fp3);
        fclose(fp4);
    }

    *retval = -0;
    return;
}

void gs_build_esvs_01(gs_varset_type *gv, gs_varset_type *es,
                      gs_varset_type *vs, int *n_estimation, char **out_es_m,
                      char **out_es_p, char **out_vs_m, char **out_vs_p,
                      int *auxfiles, int *retval)
{
    gs_varset_type tmp_es, tmp_vs;
    *retval = -2;
    if (gv == es || gv == vs || es == vs || *n_estimation < 1 ||
        (gs_UNSIGNED)*n_estimation >= gv->NoInd)
    {
        gs_info(
            gv, -2,
            "Estimation and validation sets must be distinct and non-empty");
        return;
    }
    gs_tmp_init_from(&tmp_es, es);
    gs_tmp_init_from(&tmp_vs, vs);
    gs_build_esvs_01_raw(gv, &tmp_es, &tmp_vs, n_estimation, out_es_m, out_es_p,
                         out_vs_m, out_vs_p, auxfiles, retval);
    if (*retval >= 0)
    {
        tmp_es.phased = tmp_vs.phased = gv->phased;
        if (gs_copy_map_primary(&tmp_es, gv) ||
            gs_copy_map_primary(&tmp_vs, gv) ||
            st_allele_statistics_prepare(&tmp_es) ||
            st_allele_statistics_prepare(&tmp_vs))
            *retval = -2;
    }
    if (*retval >= 0)
    {
        if (gs_prepare_allele_commit(&tmp_es) ||
            gs_prepare_allele_commit(&tmp_vs))
            *retval = -2;
        else
        {
            gs_tmp_install(es, &tmp_es);
            gs_tmp_install(vs, &tmp_vs);
        }
    }
    gs_reset(&tmp_es);
    gs_reset(&tmp_vs);
}

void gs_build_esvs_01_GV(int *n_estimation, char **base_set,
                         char **estimation_set, char **validation_set,
                         char **out_es_m, char **out_es_p, char **out_vs_m,
                         char **out_vs_p, int *auxfiles, int *retval)
{

    if (0 == esi)
    {
        gs_init(es);
        esi = 1;
    }
    if (0 == vsi)
    {
        gs_init(vs);
        vsi = 1;
    }

    gs_varset_type *ees = es;
    gs_varset_type *vvs = vs;

    if (0 != strcmp("none", *estimation_set))
        ees = gs_fdta(*estimation_set);

    if (0 != strcmp("none", *validation_set))
        vvs = gs_fdta(*validation_set);

    gs_build_esvs_01(gs_fdta(*base_set), ees, vvs, n_estimation, out_es_m,
                     out_es_p, out_vs_m, out_vs_p, auxfiles, retval);
}

void gs_return_effects_01(gs_varset_type *gv, char **out_filename, int *retval)
{
    gs_UNSIGNED i;
    FILE *outfp;

    if (NULL == gv->u)
    {
        gs_info(gv, -2, "No effects estimated");
        *retval = -2;
        return;
    }

    if ((outfp = fopen(*out_filename, "w")) == NULL)
        ERR_F;

    fprintf(outfp, "marker allele cmpallele effect\n");
    for (i = 0; i < gv->dimu; i++)
        fprintf(outfp, "%s %i %i %20.16f \n", gv->EffNme[i], gv->EffAll[i],
                gv->CmpAll[i], gv->u[i]);

    fclose(outfp);

    *retval = 0;
    return;
}

void gs_return_effects_01_GV(char **out_filename, int *retval, char **set_name)
{
    gs_return_effects_01(gs_fdta(*set_name), out_filename, retval);
}

void gs_set_effects_01(gs_varset_type *gv, char **filename, int *retval)
{

    gs_UNSIGNED i, j, k;

    if (0 == gv->NoMar)
    {
        gs_info(gv, -2, "No marker data available");
        *retval = -2;
        return;
    }

    GS_FREE(gv->Z);
    GS_FREE(gv->est);

    gv->NoEff = 1 + gv->NoMar;
    gv->EffNme = realloc(gv->EffNme, gv->NoEff * sizeof(char[gs_MAXNME]));
    gv->EffAll = realloc(gv->EffAll, gv->NoEff * sizeof(int));
    gv->CmpAll = realloc(gv->CmpAll, gv->NoEff * sizeof(int));
    if ((NULL == gv->EffNme) || (NULL == gv->EffAll) || (NULL == gv->CmpAll))
        ERR_M;

    strcpy(gv->EffNme[0], "mu");
    for (i = 1; i < gv->NoEff; i++)
        strcpy(gv->EffNme[i], gv->MarName[i - 1]);

    gv->dimu = 1 + gv->NoMar;
    gv->u = realloc(gv->u, gv->dimu * sizeof(double));
    if (NULL == gv->u)
        ERR_M;

    for (i = 0; i < gv->dimu; i++)
    {
        gv->EffAll[i] = 0;
        gv->u[i] = 0;
    }

    {
        gs_UNSIGNED NoLines = 0;

        FILE *fp;

        char *line = NULL;
        char **linep = &line;
        size_t n = 0;
        size_t *np = &n;

        char *p;

        if (NULL == (fp = fopen(*filename, "r")))
            ERR_F;

        st_getline(linep, np, fp);

        while (!feof(fp))
        {
            if (0 < st_getline(linep, np, fp))
            {

                p = strtok(line, gs_tok_del);
                if (NULL == p)
                {
                    gs_info(gv, -2, "Error in map file");
                    *retval = -2;
                    return;
                }
                if (strlen(p) > gs_MAXNME)
                {
                    gs_info(gv, -2, "Marker name to long");
                    *retval = -2;
                    return;
                }

                p = strtok(NULL, gs_tok_del);
                if (NULL == p)
                {
                    gs_info(gv, -2, "Error in map file");
                    *retval = -2;
                    return;
                }

                p = strtok(NULL, gs_tok_del);
                if (NULL == p)
                {
                    gs_info(gv, -2, "Error in map file");
                    *retval = -2;
                    return;
                }

                p = strtok(NULL, gs_tok_del);
                if (NULL == p)
                {
                    gs_info(gv, -2, "Error in map file");
                    *retval = -2;
                    return;
                }

                strcpy(line, "\0");
                NoLines++;
            }
        }
        fclose(fp);

        char(*EffNme)[1 + gs_MAXNME] = NULL;
        int *EffAll = NULL;
        int *CmpAll = NULL;
        double *Eff = NULL;

        EffNme = realloc(EffNme, NoLines * sizeof(char[1 + gs_MAXNME]));
        EffAll = realloc(EffAll, NoLines * sizeof(int));
        CmpAll = realloc(CmpAll, NoLines * sizeof(int));
        Eff = realloc(Eff, NoLines * sizeof(double));
        if ((NULL == EffNme) || (NULL == EffAll) || (NULL == CmpAll) ||
            (NULL == Eff))
            ERR_M;

        if (NULL == (fp = fopen(*filename, "r")))
            ERR_F;
        st_getline(linep, np, fp);
        for (j = 0; j < NoLines; j++)
        {
            st_getline(linep, np, fp);
            gs_remove_newline(line);

            p = strtok(line, gs_tok_del);
            strcpy(EffNme[j], p);

            p = strtok(NULL, gs_tok_del);
            {
                long allele_value;
                if ((NULL == p) ||
                    (0 != gs_parse_long_strict(p, &allele_value)) ||
                    (allele_value < 0) || (allele_value > GS_ALLELE_MAX))
                {
                    fclose(fp);
                    free(line);
                    free(EffNme);
                    free(EffAll);
                    free(CmpAll);
                    free(Eff);
                    gs_info(gv, -2, "Invalid effect allele code");
                    *retval = -2;
                    return;
                }
                EffAll[j] = (int)allele_value;
            }

            p = strtok(NULL, gs_tok_del);
            {
                long allele_value;
                if ((NULL == p) ||
                    (0 != gs_parse_long_strict(p, &allele_value)) ||
                    (allele_value < 0) || (allele_value > GS_ALLELE_MAX))
                {
                    fclose(fp);
                    free(line);
                    free(EffNme);
                    free(EffAll);
                    free(CmpAll);
                    free(Eff);
                    gs_info(gv, -2, "Invalid complementary allele code");
                    *retval = -2;
                    return;
                }
                CmpAll[j] = (int)allele_value;
            }

            p = strtok(NULL, gs_tok_del);
            Eff[j] = (double)atof(p);
        }
        fclose(fp);
        free(line);

        for (j = 0; j < NoLines; j++)
        {
            if (0 == strcmp(EffNme[j], gv->EffNme[j]))
            {
                gv->EffAll[j] = EffAll[j];
                gv->CmpAll[j] = CmpAll[j];
                gv->u[j] = Eff[j];
            }
            else
            {
                for (k = 0; k < (1 + gv->NoMar); k++)
                {
                    if (0 == strcmp(EffNme[j], gv->EffNme[k]))
                    {
                        gv->EffAll[k] = EffAll[j];
                        gv->CmpAll[k] = CmpAll[j];
                        gv->u[k] = Eff[j];
                        break;
                    }
                }
            }
        }

        free(EffNme);
        free(EffAll);
        free(CmpAll);
        free(Eff);
    }

    *retval = 0;
    return;
}

void gs_set_effects_01_GV(char **filename, int *retval, char **set_name)
{
    gs_set_effects_01(gs_fdta(*set_name), filename, retval);
}

void gs_return_pvals_01(gs_varset_type *gv, char **out_filename, int *retval)
{
    gs_UNSIGNED i;
    FILE *outfp;

    if (NULL == gv->p)
    {
        gs_info(gv, -2, "No p-values estimated");
        *retval = -2;
        return;
    }

    if ((outfp = fopen(*out_filename, "w")) == NULL)
        ERR_F;

    fprintf(outfp, " effect pvalue \n");
    for (i = 0; i < gv->dimu; i++)
        fprintf(outfp, "%s.%i %f %f \n", gv->EffNme[i], gv->EffAll[i], gv->u[i],
                gv->p[i]);

    fclose(outfp);

    *retval = 0;
    return;
}

void gs_return_pvals_01_GV(char **out_filename, int *retval, char **set_name)
{
    gs_return_pvals_01(gs_fdta(*set_name), out_filename, retval);
}

void gs_check_pvals_01(gs_varset_type *gv, int *retval)
{
    if (NULL == gv->p)
        *retval = 0;
    else
        *retval = 1;
    return;
}

void gs_check_pvals_01_GV(int *retval, char **set_name)
{
    gs_check_pvals_01(gs_fdta(*set_name), retval);
}

void gs_set_allele_codes_01(gs_varset_type *gv, double *aa, double *aA,
                            double *AA, double *an, double *nn, double *nA)
{
    gv->aa = *aa;
    gv->aA = *aA;
    gv->AA = *AA;
    gv->an = *an;
    gv->nn = *nn;
    gv->nA = *nA;

    GS_FREE(gv->Z);
    GS_FREE(gv->est);
}

void gs_set_allele_codes_01_GV(double *aa, double *aA, double *AA, double *an,
                               double *nn, double *nA, char **set_name)
{

    gs_set_allele_codes_01(gs_fdta(*set_name), aa, aA, AA, an, nn, nA);
}

int gs_mappoint_cmp(const void *a, const void *b)
{
    const gs_mappoint_type **ia = (const gs_mappoint_type **)a;
    const gs_mappoint_type **ib = (const gs_mappoint_type **)b;
    int cmp;

    if ((*ia)->chrom > (*ib)->chrom)
        return 1;

    if ((*ia)->chrom < (*ib)->chrom)
        return -1;

    if ((*ia)->pos > (*ib)->pos)
        return 1;

    if ((*ia)->pos < (*ib)->pos)
        return -1;

    cmp = strcmp((*ia)->name, (*ib)->name);

    if (0 != cmp)
        return cmp;

    return strcmp((*ia)->class, (*ib)->class);
}

void gs_reset_all_but_marker_matrix(gs_varset_type *gv)
{
    gs_UNSIGNED NoInd = gv->NoInd;
    gs_UNSIGNED NoMar = gv->NoMar;
    int phased = gv->phased;
    struct st_marker_state *prediction = gv->prediction_state;
    int prediction_role = gv->prediction_role;
    gv->prediction_state = NULL;
    gv->prediction_role = 0;
    double *performance = gv->y;
    gv->y = NULL;

    char **IndName = gv->IndName;
    char **Pedigree = gv->Pedigree;
    gv->Pedigree = NULL;
    gv->IndName = NULL;
    char(*MarName)[gs_MAXNME] = gv->MarName;
    gv->MarName = NULL;

    GS_ALLELE_TYPE *All1 = gv->allele_raw1;
    gv->allele_raw1 = NULL;
    GS_ALLELE_TYPE *All2 = gv->allele_raw2;
    gv->allele_raw2 = NULL;
    unsigned char *packed = gv->allele_packed,
                  *dictionary = gv->allele_dictionary;
    size_t packed_bytes = gv->allele_packed_bytes;
    unsigned int bits = gv->allele_bits;
    gv->allele_packed = gv->allele_dictionary = NULL;
    gv->allele_packed_bytes = 0;
    gv->allele_bits = 16;

    gs_reset(gv);

    gv->NoInd = NoInd;
    gv->NoMar = NoMar;
    gv->phased = phased;
    gv->prediction_state = prediction;
    gv->prediction_role = prediction_role;
    gv->y = performance;

    gv->IndName = IndName;
    gv->Pedigree = Pedigree;
    gv->MarName = MarName;

    gv->allele_raw1 = All1;
    gv->allele_raw2 = All2;
    gv->allele_packed = packed;
    gv->allele_dictionary = dictionary;
    gv->allele_packed_bytes = packed_bytes;
    gv->allele_bits = bits;
}

void gs_reset_chrom_stats(gs_varset_type *gv)
{
    gs_UNSIGNED c;

    if (gv->chrom != NULL)
    {
        for (c = 0; c < gv->NChrom; c++)
        {
            GS_FREE(gv->chrom[c].rf);
            GS_FREE(gv->chrom[c].q);
            GS_FREE(gv->chrom[c].ld);
        }

        GS_FREE(gv->chrom);
        gv->NChrom = 0;
    }
}

void gs_chrom_stats_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                       int *retval)
{
    gs_UNSIGNED c, i;

    if (0 == gv->NLoc)
    {

        snprintf(gv->msg, sizeof(gv->msg), "No linkage map has been loaded");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    gs_reset_chrom_stats(gv);

    gv->NChrom = 1;
    for (i = 1; i < gv->NLoc; i++)
        if (gv->map[i - 1]->chrom != gv->map[i]->chrom)
            gv->NChrom++;

    gv->chrom = realloc(gv->chrom, gv->NChrom * sizeof(gs_chrom_type));

    for (c = 0; c < gv->NChrom; c++)
    {
        gv->chrom[c].rf = NULL;
        gv->chrom[c].q = NULL;
        gv->chrom[c].ld = NULL;
    }

    snprintf(gv->msg, sizeof(gv->msg), "No. of chromosomes: %lu", gv->NChrom);

    gs_info(gv, 1, gv->msg);

    gv->chrom[0].NLocC = 1;
    gv->chrom[0].first_loc = &(gv->map[0]);
    gv->chrom[0].first_loc_i = 0;
    for (c = 0, i = 1; i < gv->NLoc; i++)
    {
        if (gv->map[i - 1]->chrom == gv->map[i]->chrom)
        {
            gv->chrom[c].NLocC++;
        }
        else
        {
            c++;
            gv->chrom[c].first_loc = &(gv->map[i]);
            gv->chrom[c].first_loc_i = i;
            gv->chrom[c].NLocC = 1;
        }
    }

    if (1 == *auxfiles)
    {
        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        gs_mappoint_type **p;

        for (c = 0; c < gv->NChrom; c++)
        {
            p = gv->chrom[c].first_loc;
            snprintf(gv->msg, sizeof(gv->msg), "%lu %lu %f \n", 1 + c,
                     gv->chrom[c].NLocC, (p[gv->chrom[c].NLocC - 1])->pos);
            fprintf(outfp, "%s", gv->msg);
        }

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_chrom_stats_01_GV(char **out_filename, int *auxfiles, int *retval,
                          char **set_name)
{
    gs_chrom_stats_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

gs_FLOAT gs_rf_haldane(gs_FLOAT d)
{
    return ((1 - exp(-2 * d)) / 2);
}

gs_FLOAT gs_rf_kosambi(gs_FLOAT d)
{
    return (exp(4 * d) - 1) / (exp(4 * d) + 1) / 2;
}

void gs_calc_rf_01(gs_varset_type *gv, char **map_function, char **out_filename,
                   int *auxfiles, int *retval)
{
    gs_UNSIGNED c, i, j, idx;
    gs_FLOAT d;
    gs_mappoint_type **p;

    if (0 == gv->NChrom)
    {
        snprintf(gv->msg, sizeof(gv->msg),
                 "No chromosome statistics available");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    gs_FLOAT (*gs_rf)(gs_FLOAT d);

    if (0 == strcmp(*map_function, "Haldane"))
    {
        gs_rf = gs_rf_haldane;
    }

    else if (0 == strcmp(*map_function, "Kosambi"))
    {
        gs_rf = gs_rf_kosambi;
    }

    else
    {
        gs_info(gv, -2, "Unknown map function");
        *retval = -2;
        return;
    }

    for (c = 0; c < gv->NChrom; c++)
    {
        gv->chrom[c].rf =
            realloc(gv->chrom[c].rf,
                    gv->chrom[c].NLocC * gv->chrom[c].NLocC * sizeof(gs_FLOAT));
        gv->chrom[c].q =
            realloc(gv->chrom[c].q,
                    gv->chrom[c].NLocC * gv->chrom[c].NLocC * sizeof(gs_FLOAT));
    }

    for (c = 0; c < gv->NChrom; c++)
    {
        p = gv->chrom[c].first_loc;
        for (i = 0; i < gv->chrom[c].NLocC; i++)
            for (j = 0; j < gv->chrom[c].NLocC; j++)
            {
                d = fabs((p[i])->pos - (p[j])->pos);

                idx = i * gv->chrom[c].NLocC + j;
                (gv->chrom[c].rf)[idx] = gs_rf(d / 100);
            }
    }

    if (1 == *auxfiles)
    {
        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        for (c = 0; c < gv->NChrom; c++)
        {
            for (i = 0; i < (gv->chrom[c].NLocC - 1); i++)
                for (j = i + 1; j < gv->chrom[c].NLocC; j++)
                {
                    idx = i * gv->chrom[c].NLocC + j;

                    snprintf(gv->msg, sizeof(gv->msg), "%lu %lu %lu %f \n",
                             1 + c, 1 + i, 1 + j,
                             (double)(gv->chrom[c].rf)[idx]);
                    fprintf(outfp, "%s", gv->msg);
                }
        }

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_calc_rf_01_GV(char **map_function, char **out_filename, int *auxfiles,
                      int *retval, char **set_name)
{
    gs_calc_rf_01(gs_fdta(*set_name), map_function, out_filename, auxfiles,
                  retval);
}

gs_FLOAT gs_q_dh(int *t, gs_FLOAT r)
{
    return 0.5 + (1 - 2 * r) / 2 * pow((double)(1 - r), (double)*t);
}

gs_FLOAT gs_q_ssd(int *t, gs_FLOAT r)
{

    return 0.5 + (1.0 - 2.0 * r) / (2.0 + 4.0 * r) *
                     pow((double)(1.0 - r), (double)(*t));
}

gs_FLOAT gs_q_crs(int *t, gs_FLOAT r)
{
    int tt = *t;
    tt = 1;
    return tt * (1.0 - r);
}

void gs_calc_q_01(gs_varset_type *gv, char **pop_type, int *t,
                  char **out_filename, int *auxfiles, int *retval)
{
    gs_UNSIGNED c, i, j, idx;

    if ((0 == gv->NChrom) || (NULL == gv->chrom[0].rf))
    {
        snprintf(gv->msg, sizeof(gv->msg),
                 "No recombination frequencies available");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    gs_FLOAT (*gs_q)(int *t, gs_FLOAT r);

    if (0 == strcmp(*pop_type, "DH"))
    {
        gs_q = gs_q_dh;
    }

    else if (0 == strcmp(*pop_type, "SSD"))
    {
        gs_q = gs_q_ssd;
    }

    else if (0 == strcmp(*pop_type, "CRS"))
    {
        gs_q = gs_q_crs;
    }

    else
    {
        gs_info(gv, -2, "Unknown population type");
        *retval = -2;
        return;
    }

    for (c = 0; c < gv->NChrom; c++)
    {
        for (i = 0; i < gv->chrom[c].NLocC; i++)
            for (j = 0; j < gv->chrom[c].NLocC; j++)
            {
                idx = i * gv->chrom[c].NLocC + j;
                (gv->chrom[c].q)[idx] = gs_q(t, (gv->chrom[c].rf)[idx]);
            }
    }

    if (1 == *auxfiles)
    {
        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        for (c = 0; c < gv->NChrom; c++)
        {
            for (i = 0; i < (gv->chrom[c].NLocC - 1); i++)
                for (j = i + 1; j < gv->chrom[c].NLocC; j++)
                {
                    idx = i * gv->chrom[c].NLocC + j;

                    snprintf(gv->msg, sizeof(gv->msg), "%lu %lu %lu %f \n",
                             1 + c, 1 + i, 1 + j,
                             (double)(gv->chrom[c].q)[idx]);
                    fprintf(outfp, "%s", gv->msg);
                }
        }

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_calc_q_01_GV(char **pop_type, int *t, char **out_filename,
                     int *auxfiles, int *retval, char **set_name)
{
    gs_calc_q_01(gs_fdta(*set_name), pop_type, t, out_filename, auxfiles,
                 retval);
}

void gs_replace_comma(char *p)
{
    char *q = strchr(p, ',');
    if (NULL != q)
        *q = '.';
}

void gs_read_map_01(gs_varset_type *gv, char **map_filename, int *skip,
                    char **fmt, double *map_stretch, char **out_filename,
                    int *auxfiles, int *retval)
{
    FILE *fp = NULL;
    char *line = NULL;
    size_t n = 0;
    long long int rc;
    gs_mappoint_type *map_s = NULL;
    gs_mappoint_build_type *map_build = NULL;
    char *map_str_s = NULL;
    size_t map_str_n = 0;
    size_t map_str_capacity = 0;
    gs_mappoint_type **map = NULL;
    gs_name_index_type *map_names = NULL;
    gs_name_index_type *marker_names = NULL;
    gs_UNSIGNED *match = NULL;
    gs_UNSIGNED nloc = 0;
    gs_UNSIGNED capacity = 0;
    gs_UNSIGNED no_matching = 0;
    gs_UNSIGNED mar_size = 0;
    gs_UNSIGNED cls_size = 0;
    gs_varset_type tmp;
    int failed = 0;

    if (!((0 == strcmp(*fmt, "cpms")) || (0 == strcmp(*fmt, "mcp"))))
    {
        gs_info(gv, -2, "Implemented map formats: cpms, mcp");
        *retval = -2;
        return;
    }

    if (*skip < 0)
    {
        gs_info(gv, -2, "Invalid number of map lines to skip");
        *retval = -2;
        return;
    }

    if (!isfinite(*map_stretch))
    {
        gs_info(gv, -2, "Invalid map stretch");
        *retval = -2;
        return;
    }

    if ((0 == gv->NoInd) || (0 == gv->NoMar) || (NULL == gv->IndName) ||
        (NULL == gv->MarName) || !st_allele_present(gv))
    {
        gs_info(gv, -2, "No marker data available");
        *retval = -2;
        return;
    }

    gs_tmp_init_from(&tmp, gv);

    fp = fopen(*map_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        *retval = -2;
        return;
    }

    for (int s = 0; s < *skip; s++)
    {
        rc = st_getline(&line, &n, fp);
        if (rc < 0)
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }
        if (0 == rc)
        {
            gs_info(gv, -2, "Map file shorter than skip value");
            failed = 1;
            goto cleanup;
        }
    }

    while (0 < (rc = st_getline(&line, &n, fp)))
    {
        char *p;
        char *name;
        char *class_name;
        long chrom_l;
        double pos;
        size_t name_len;
        size_t class_len;

        gs_remove_newline(line);

        if (0 == strcmp(*fmt, "cpms"))
        {
            p = strtok(line, gs_tok_del);
            if ((NULL == p) || (0 != gs_parse_long_strict(p, &chrom_l)) ||
                (0 == chrom_l) || (chrom_l < INT_MIN) || (chrom_l > INT_MAX))
            {
                gs_info(gv, -2, "Error in map file chromosome");
                failed = 1;
                goto cleanup;
            }

            p = strtok(NULL, gs_tok_del);
            if ((NULL == p) || (0 != gs_parse_double_strict(p, &pos)))
            {
                gs_info(gv, -2, "Error in map file position");
                failed = 1;
                goto cleanup;
            }

            name = strtok(NULL, gs_tok_del);
            class_name = strtok(NULL, gs_tok_del);
            if ((NULL == name) || (NULL == class_name))
            {
                gs_info(gv, -2, "Error in map file");
                failed = 1;
                goto cleanup;
            }
        }
        else
        {
            name = strtok(line, gs_tok_del);
            p = strtok(NULL, gs_tok_del);
            if ((NULL == name) || (NULL == p) ||
                (0 != gs_parse_long_strict(p, &chrom_l)) || (0 == chrom_l) ||
                (chrom_l < INT_MIN) || (chrom_l > INT_MAX))
            {
                gs_info(gv, -2, "Error in map file chromosome");
                failed = 1;
                goto cleanup;
            }

            p = strtok(NULL, gs_tok_del);
            if ((NULL == p) || (0 != gs_parse_double_strict(p, &pos)))
            {
                gs_info(gv, -2, "Error in map file position");
                failed = 1;
                goto cleanup;
            }

            class_name = "m";
        }

        if (NULL != strtok(NULL, gs_tok_del))
        {
            gs_info(gv, -2, "Too many fields in map file line");
            failed = 1;
            goto cleanup;
        }

        pos *= *map_stretch;
        if (!isfinite(pos))
        {
            gs_info(gv, -2, "Invalid stretched map position");
            failed = 1;
            goto cleanup;
        }

        name_len = strlen(name);
        class_len = strlen(class_name);
        if ((name_len >= gs_MAXSTR) || (class_len >= gs_MAXSTR))
        {
            gs_info(gv, -2, "Description too long");
            failed = 1;
            goto cleanup;
        }
        if (name_len > mar_size)
            mar_size = (gs_UNSIGNED)name_len;
        if (class_len > cls_size)
            cls_size = (gs_UNSIGNED)class_len;

        if (nloc == capacity)
        {
            gs_UNSIGNED new_capacity = (0 == capacity) ? 128 : capacity * 2;
            gs_mappoint_build_type *new_map_build;

            if (new_capacity < capacity)
            {
                gs_info(gv, -2, "Map file too large");
                failed = 1;
                goto cleanup;
            }

            new_map_build = gs_realloc_array(map_build, (size_t)new_capacity,
                                             sizeof(*new_map_build));
            if (NULL == new_map_build)
            {
                gs_info(gv, -2, "Memory allocation error");
                failed = 1;
                goto cleanup;
            }

            map_build = new_map_build;
            capacity = new_capacity;
        }

        map_build[nloc].chrom = (int)chrom_l;
        map_build[nloc].pos = pos;
        map_build[nloc].nr = gs_UNSIGNED_MAX;

        if ((0 != gs_string_builder_append(&map_str_s, &map_str_n,
                                           &map_str_capacity, name,
                                           &(map_build[nloc].name_offset))) ||
            (0 != gs_string_builder_append(&map_str_s, &map_str_n,
                                           &map_str_capacity, class_name,
                                           &(map_build[nloc].class_offset))))
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }

        nloc++;
    }

    if (rc < 0)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    fclose(fp);
    fp = NULL;

    if (0 == nloc)
    {
        gs_info(gv, -2, "Linkage map contains no loci");
        failed = 1;
        goto cleanup;
    }

    map_s = gs_malloc_array((size_t)nloc, sizeof(*map_s));
    if (NULL == map_s)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
    {
        map_s[lne].chrom = map_build[lne].chrom;
        map_s[lne].pos = map_build[lne].pos;
        map_s[lne].name = map_str_s + map_build[lne].name_offset;
        map_s[lne].class = map_str_s + map_build[lne].class_offset;
        map_s[lne].nr = map_build[lne].nr;
    }

    free(map_build);
    map_build = NULL;

    snprintf(gv->msg, sizeof(gv->msg), "Description length: Mar %lu, Class %lu",
             mar_size, cls_size);
    gs_info(gv, 1, gv->msg);
    snprintf(gv->msg, sizeof(gv->msg), "No. of markers in map file: %lu", nloc);
    gs_info(gv, 1, gv->msg);

    map = gs_malloc_array((size_t)nloc, sizeof(*map));
    map_names = gs_malloc_array((size_t)nloc, sizeof(*map_names));
    match = gs_malloc_array((size_t)nloc, sizeof(*match));
    if ((NULL == map) || (NULL == map_names) || (NULL == match))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
        map[lne] = &(map_s[lne]);

    qsort(map, nloc, sizeof(*map), gs_mappoint_cmp);

    for (gs_UNSIGNED first = 0; first < nloc;)
    {
        gs_UNSIGNED last = first + 1;
        double base = map[first]->pos;
        int chrom = map[first]->chrom;

        while ((last < nloc) && (map[last]->chrom == chrom) &&
               (map[last]->pos == base))
            last++;

        if (last - first > 1)
        {
            gs_UNSIGNED k = last - first;
            double step = 0.00123456;
            int has_next = ((last < nloc) && (map[last]->chrom == chrom));

            if (has_next)
            {
                double gap = map[last]->pos - base;
                double preferred = step * (double)(k - 1);

                if (!(gap > 0.0))
                {
                    gs_info(gv, -2, "Invalid linkage-map position ordering");
                    failed = 1;
                    goto cleanup;
                }

                if (preferred >= gap)
                    step = gap / (double)k;
            }

            for (gs_UNSIGNED j = 1; j < k; j++)
            {
                double pos = base + step * (double)j;
                double prev = map[first + j - 1]->pos;

                if (!(pos > prev))
                    pos = nextafter(prev, INFINITY);

                if (!isfinite(pos))
                {
                    gs_info(gv, -2,
                            "Cannot separate identical linkage-map positions");
                    failed = 1;
                    goto cleanup;
                }

                if (has_next && !(pos < map[last]->pos))
                {
                    gs_info(gv, -2,
                            "Cannot separate identical linkage-map positions");
                    failed = 1;
                    goto cleanup;
                }

                map[first + j]->pos = pos;
            }
        }

        first = last;
    }

    for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
    {
        map_names[lne].name = map[lne]->name;
        map_names[lne].idx = lne;
        match[lne] = gs_UNSIGNED_MAX;
    }
    qsort(map_names, nloc, sizeof(*map_names), gs_name_index_cmp);

    if (gs_name_index_has_duplicates(map_names, nloc))
    {
        gs_info(gv, -2, "Duplicate marker name in linkage map");
        failed = 1;
        goto cleanup;
    }

    if (0 != gs_make_name_index(gv->MarName, gv->NoMar, &marker_names))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    if (gs_name_index_has_duplicates(marker_names, gv->NoMar))
    {
        gs_info(gv, -2, "Duplicate marker name in marker data");
        failed = 1;
        goto cleanup;
    }

    {
        gs_UNSIGNED a = 0;
        gs_UNSIGNED b = 0;

        while ((a < nloc) && (b < gv->NoMar))
        {
            int cmp = strcmp(map_names[a].name, marker_names[b].name);

            if (0 == cmp)
            {
                match[map_names[a].idx] = marker_names[b].idx;
                no_matching++;
                a++;
                b++;
            }
            else if (cmp < 0)
                a++;
            else
                b++;
        }
    }

    if (0 == no_matching)
    {
        gs_info(gv, -2, "No linkage-map markers match marker data");
        failed = 1;
        goto cleanup;
    }

    tmp.NoInd = gv->NoInd;
    tmp.NoMar = no_matching;
    tmp.phased = gv->phased;
    tmp.NLoc = no_matching;
    tmp.MaxIL = gv->MaxIL;
    tmp.MaxML = gv->MaxML;

    for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
        if (gs_UNSIGNED_MAX != match[lne])
            if ((0 !=
                 gs_pool_add_string_size(&(tmp.map_str_N), map[lne]->name)) ||
                (0 !=
                 gs_pool_add_string_size(&(tmp.map_str_N), map[lne]->class)))
            {
                gs_info(gv, -2, "Linkage map string storage too large");
                failed = 1;
                goto cleanup;
            }

    tmp.IndName = calloc((size_t)tmp.NoInd, sizeof(*tmp.IndName));
    tmp.MarName = gs_malloc_array((size_t)tmp.NoMar, sizeof(*tmp.MarName));
    tmp.map_S = gs_malloc_array((size_t)tmp.NLoc, sizeof(*tmp.map_S));
    tmp.map_str_S = gs_malloc_array(tmp.map_str_N, sizeof(*tmp.map_str_S));
    tmp.map = gs_malloc_array((size_t)tmp.NLoc, sizeof(*tmp.map));
    tmp.mdl = gs_malloc_array((size_t)tmp.NLoc, sizeof(*tmp.mdl));

    if ((NULL == tmp.IndName) || (NULL == tmp.MarName) || (NULL == tmp.map_S) ||
        (NULL == tmp.map_str_S) || (NULL == tmp.map) || (NULL == tmp.mdl) ||
        (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED i = 0; i < tmp.NoInd; i++)
        if (st_identity_copy(&tmp, i, gv, i))
        {
            failed = 1;
            goto cleanup;
        }

    {
        gs_UNSIGNED m_new = 0;
        size_t map_str_offset = 0;

        for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
            if (gs_UNSIGNED_MAX != match[lne])
            {
                gs_UNSIGNED m_old = match[lne];

                strcpy(tmp.MarName[m_new], gv->MarName[m_old]);
                st_allele_copy(gv, 0, (size_t)m_old * gv->NoInd,
                               (size_t)tmp.NoInd,
                               &(tmp.allele_raw1[(size_t)m_new * tmp.NoInd]));
                st_allele_copy(gv, 1, (size_t)m_old * gv->NoInd,
                               (size_t)tmp.NoInd,
                               &(tmp.allele_raw2[(size_t)m_new * tmp.NoInd]));

                if (0 != gs_copy_map_point_to_pool(
                             &(tmp.map_S[m_new]), map[lne], m_new,
                             tmp.map_str_S, tmp.map_str_N, &map_str_offset))
                {
                    gs_info(gv, -2, "Invalid linkage map string storage");
                    failed = 1;
                    goto cleanup;
                }
                tmp.map[m_new] = &(tmp.map_S[m_new]);
                tmp.mdl[m_new] = m_new;
                m_new++;
            }

        if (map_str_offset != tmp.map_str_N)
        {
            gs_info(gv, -2, "Invalid linkage map string storage");
            failed = 1;
            goto cleanup;
        }
    }

    if (NULL != gv->y)
    {
        tmp.y = gs_malloc_array((size_t)tmp.NoInd, sizeof(*tmp.y));
        if (NULL == tmp.y)
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }
        memcpy(tmp.y, gv->y, (size_t)tmp.NoInd * sizeof(*tmp.y));
    }

    {
        int srv = 0;
        char *empty = "";
        char **empty_p = &empty;
        int zero = 0;
        int *zero_p = &zero;

        gs_marker_stats_01(&tmp, empty_p, empty_p, empty_p, zero_p, zero_p,
                           zero_p, zero_p, &srv);
        if (-2 == srv)
        {
            failed = 1;
            goto cleanup;
        }
    }

    if (1 == *auxfiles)
    {
        FILE *outfp = fopen(*out_filename, "w");
        if (NULL == outfp)
        {
            gs_info(gv, -2, "File open error");
            failed = 1;
            goto cleanup;
        }

        for (gs_UNSIGNED lne = 0; lne < tmp.NLoc; lne++)
            fprintf(outfp, "%i %15.12f %s %s\n", tmp.map[lne]->chrom,
                    tmp.map[lne]->pos, tmp.map[lne]->name, tmp.map[lne]->class);
        fclose(outfp);
    }

    snprintf(gv->msg, sizeof(gv->msg), "No. of matching markers: %lu",
             tmp.NLoc);
    gs_info(gv, 1, gv->msg);

    if (gs_tmp_commit(gv, &tmp))
    {
        failed = 1;
        goto cleanup;
    }
    *retval = 0;

cleanup:
    if (NULL != fp)
        fclose(fp);
    free(line);
    free(map_build);
    free(map_s);
    free(map_str_s);
    free(map);
    free(map_names);
    free(marker_names);
    free(match);

    if (failed)
    {
        gs_reset(&tmp);
        *retval = -2;
    }
}

void gs_read_map_01_GV(char **map_filename, int *skip, char **fmt,
                       double *map_stretch, char **out_filename, int *auxfiles,
                       int *retval, char **set_name)
{

    gs_read_map_01(gs_fdta(*set_name), map_filename, skip, fmt, map_stretch,
                   out_filename, auxfiles, retval);
}

static SEXP st_stvcf_element(SEXP list, const char *name)
{
    SEXP names;

    if (VECSXP != TYPEOF(list))
        return R_NilValue;

    names = getAttrib(list, R_NamesSymbol);
    if (STRSXP != TYPEOF(names))
        return R_NilValue;

    for (R_xlen_t i = 0; i < XLENGTH(list); i++)
        if ((NA_STRING != STRING_ELT(names, i)) &&
            (0 == strcmp(CHAR(STRING_ELT(names, i)), name)))
            return VECTOR_ELT(list, i);

    return R_NilValue;
}

static int st_stvcf_components(SEXP list)
{
    static const char *required[9] = {"CHROM",   "POS",     "ID",
                                      "REF",     "ALT",     "individual",
                                      "allele1", "allele2", "CHROM.levels"};
    int seen[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    SEXP names;

    if ((VECSXP != TYPEOF(list)) || (9 != XLENGTH(list)))
        return 0;

    names = getAttrib(list, R_NamesSymbol);
    if ((STRSXP != TYPEOF(names)) || (9 != XLENGTH(names)))
        return 0;

    for (R_xlen_t i = 0; i < 9; i++)
    {
        SEXP s = STRING_ELT(names, i);
        int found = -1;

        if (NA_STRING == s)
            return 0;

        for (int j = 0; j < 9; j++)
            if (0 == strcmp(CHAR(s), required[j]))
            {
                found = j;
                break;
            }

        if ((found < 0) || seen[found])
            return 0;
        seen[found] = 1;
    }

    return 1;
}

static int st_stvcf_phased_attribute(SEXP stvcf, int *phased)
{
    SEXP x;

    if (NULL == phased)
        return -1;

    *phased = 0;
    x = getAttrib(stvcf, install("phased"));
    if (R_NilValue == x)
        return 0;

    if ((LGLSXP != TYPEOF(x)) || (1 != XLENGTH(x)) ||
        (NA_LOGICAL == LOGICAL(x)[0]))
        return -1;

    *phased = LOGICAL(x)[0] ? 1 : 0;
    return 0;
}

static int st_stvcf_matrix_dim(SEXP x, R_xlen_t nr, R_xlen_t nc)
{
    SEXP dim = getAttrib(x, R_DimSymbol);

    if ((INTSXP != TYPEOF(dim)) || (2 != XLENGTH(dim)))
        return 0;

    return ((INTEGER(dim)[0] == nr) && (INTEGER(dim)[1] == nc));
}

static int st_stvcf_valid_name(const char *name)
{
    const unsigned char *p = (const unsigned char *)name;

    if ((NULL == name) || ('\0' == name[0]))
        return 0;

    for (; '\0' != *p; p++)
        if ((*p <= 32u) || (127u == *p))
            return 0;

    return 1;
}

static int st_stvcf_cstr_cmp(const void *a, const void *b)
{
    const char *const *aa = (const char *const *)a;
    const char *const *bb = (const char *const *)b;

    return strcmp(*aa, *bb);
}

static int st_stvcf_alleles(const char *ref, const char *alt,
                            gs_UNSIGNED *n_alt)
{
    char *copy = NULL;
    char **tokens = NULL;
    size_t len;
    gs_UNSIGNED n = 1;
    gs_UNSIGNED k = 0;
    int rc = -1;

    if ((NULL == ref) || ('\0' == ref[0]) || (0 == strcmp(ref, ".")) ||
        (NULL != strchr(ref, ',')) || (NULL == alt) || ('\0' == alt[0]))
        return -1;

    if (0 == strcmp(alt, "."))
    {
        *n_alt = 0;
        return 0;
    }

    for (const char *p = alt; '\0' != *p; p++)
        if (',' == *p)
        {
            if (n >= (gs_UNSIGNED)(GS_ALLELE_MAX - 1))
                return -1;
            n++;
        }

    if (1 == n)
    {
        if ((0 == strcmp(alt, ref)) || ('\0' == alt[0]))
            return -1;
        *n_alt = 1;
        return 0;
    }

    len = strlen(alt);
    copy = malloc(len + 1);
    tokens = gs_malloc_array((size_t)n, sizeof(*tokens));
    if ((NULL == copy) || (NULL == tokens))
        goto cleanup;

    memcpy(copy, alt, len + 1);
    tokens[k++] = copy;
    for (size_t i = 0; i < len; i++)
        if (',' == copy[i])
        {
            copy[i] = '\0';
            tokens[k++] = copy + i + 1;
        }

    if (k != n)
        goto cleanup;

    for (gs_UNSIGNED i = 0; i < n; i++)
        if (('\0' == tokens[i][0]) || (0 == strcmp(tokens[i], ".")) ||
            (0 == strcmp(tokens[i], ref)))
            goto cleanup;

    qsort(tokens, n, sizeof(*tokens), st_stvcf_cstr_cmp);
    for (gs_UNSIGNED i = 1; i < n; i++)
        if (0 == strcmp(tokens[i - 1], tokens[i]))
            goto cleanup;

    *n_alt = n;
    rc = 0;

cleanup:
    free(copy);
    free(tokens);
    return rc;
}

static int st_stvcf_get_allele(SEXP x, R_xlen_t idx, gs_UNSIGNED n_alt,
                               int *missing, GS_ALLELE_TYPE *allele)
{
    long value;

    if ((NULL == missing) || (NULL == allele) ||
        (n_alt > (gs_UNSIGNED)(GS_ALLELE_MAX - 1)))
        return -1;

    if (RAWSXP == TYPEOF(x))
    {
        unsigned char raw = RAW(x)[idx];

        if (255 == raw)
        {
            *missing = 1;
            *allele = GS_ALLELE_MISSING;
            return 0;
        }
        value = raw;
    }
    else if (INTSXP == TYPEOF(x))
    {
        int integer = INTEGER(x)[idx];

        if (NA_INTEGER == integer)
        {
            *missing = 1;
            *allele = GS_ALLELE_MISSING;
            return 0;
        }
        value = integer;
    }
    else
        return -1;

    if ((value < 0) || ((gs_UNSIGNED)value > n_alt) || (value >= GS_ALLELE_MAX))
        return -1;

    *missing = 0;
    *allele = (GS_ALLELE_TYPE)(value + 1);
    return 0;
}

static int st_stvcf_map_cmp(const void *a, const void *b)
{
    const gs_mappoint_type *const *ia = (const gs_mappoint_type *const *)a;
    const gs_mappoint_type *const *ib = (const gs_mappoint_type *const *)b;

    if ((*ia)->chrom < (*ib)->chrom)
        return -1;
    if ((*ia)->chrom > (*ib)->chrom)
        return 1;
    if ((*ia)->pos < (*ib)->pos)
        return -1;
    if ((*ia)->pos > (*ib)->pos)
        return 1;
    if ((*ia)->nr < (*ib)->nr)
        return -1;
    if ((*ia)->nr > (*ib)->nr)
        return 1;
    return 0;
}

static int st_stvcf_separate_map_ties(gs_varset_type *gv,
                                      gs_mappoint_type **map, gs_UNSIGNED nloc)
{
    for (gs_UNSIGNED first = 0; first < nloc;)
    {
        gs_UNSIGNED last = first + 1;
        double base = map[first]->pos;
        int chrom = map[first]->chrom;

        while ((last < nloc) && (map[last]->chrom == chrom) &&
               (map[last]->pos == base))
            last++;

        if (last - first > 1)
        {
            gs_UNSIGNED k = last - first;
            double step = 0.00123456;
            int has_next = ((last < nloc) && (map[last]->chrom == chrom));

            if (has_next)
            {
                double gap = map[last]->pos - base;
                double preferred = step * (double)(k - 1);

                if (!(gap > 0.0))
                {
                    gs_info(gv, -2, "Invalid linkage-map position ordering");
                    return -1;
                }

                if (preferred >= gap)
                    step = gap / (double)k;
            }

            for (gs_UNSIGNED j = 1; j < k; j++)
            {
                double pos = base + step * (double)j;
                double prev = map[first + j - 1]->pos;

                if (!(pos > prev))
                    pos = nextafter(prev, INFINITY);

                if (!isfinite(pos) || (has_next && !(pos < map[last]->pos)))
                {
                    gs_info(gv, -2,
                            "Cannot separate identical linkage-map positions");
                    return -1;
                }

                map[first + j]->pos = pos;
            }
        }

        first = last;
    }

    return 0;
}

static int st_markerdata_observed_alleles(const gs_varset_type *gv,
                                          gs_UNSIGNED marker,
                                          GS_ALLELE_TYPE *alleles,
                                          gs_UNSIGNED *n_alleles,
                                          int *has_missing)
{
    gs_UNSIGNED n = 0;
    int missing = 0;
    st_allele_marker_view view = st_allele_marker(gv, marker);

    for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
    {
        st_allele_pair pair = st_allele_read_pair(&view, i);
        GS_ALLELE_TYPE a1 = pair.a, a2 = pair.b;
        GS_ALLELE_TYPE values[2] = {a1, a2};

        if ((GS_ALLELE_MISSING == a1) || (GS_ALLELE_MISSING == a2))
        {
            missing = 1;
            continue;
        }

        for (int h = 0; h < 2; h++)
        {
            int found = 0;

            for (gs_UNSIGNED k = 0; k < n; k++)
                if (alleles[k] == values[h])
                {
                    found = 1;
                    break;
                }

            if (!found)
            {
                if (n >= gs_MAXALL)
                    return -1;
                alleles[n++] = values[h];
            }
        }
    }

    if ((gs_UNSIGNED)missing + n > gs_MAXALL)
        return -1;

    qsort(alleles, n, sizeof(*alleles), gs_allele_cmp);
    *n_alleles = n;
    *has_missing = missing;
    return 0;
}

static int st_markerdata_allele_index(const GS_ALLELE_TYPE *alleles,
                                      gs_UNSIGNED n, int allele)
{
    gs_UNSIGNED left = 0;
    gs_UNSIGNED right = n;

    while (left < right)
    {
        gs_UNSIGNED mid = left + (right - left) / 2;

        if (alleles[mid] < allele)
            left = mid + 1;
        else
            right = mid;
    }

    if ((left < n) && (alleles[left] == allele))
        return (int)left;

    return -1;
}

#define ST_VCF_POSITION_SCALE 1000000.0

int st_identifier_valid(const char *name)
{
    const unsigned char *p = (const unsigned char *)name;

    if ((NULL == name) || ('\0' == name[0]))
        return 0;

    for (; '\0' != *p; p++)
        if ((*p <= 32u) || (127u == *p))
            return 0;

    return 1;
}

static char *st_vcf_strdup(const char *s)
{
    size_t n;
    char *out;

    if (NULL == s)
        return NULL;

    n = strlen(s) + 1u;
    out = malloc(n);
    if (NULL != out)
        memcpy(out, s, n);
    return out;
}

static int st_vcf_getline(FILE *fp, char **line, size_t *capacity)
{
    size_t used = 0u;

    if ((NULL == fp) || (NULL == line) || (NULL == capacity))
        return -1;

    if ((*capacity < 8192u) || (NULL == *line))
    {
        char *tmp = realloc(*line, 8192u);
        if (NULL == tmp)
            return -1;
        *line = tmp;
        *capacity = 8192u;
    }

    (*line)[0] = '\0';
    for (;;)
    {
        size_t avail_size;
        int avail;
        char *z;
        size_t added;

        if (used + 2u >= *capacity)
        {
            size_t next = (*capacity > SIZE_MAX / 2u) ? 0u : 2u * (*capacity);
            char *tmp;

            if ((0u == next) || (next > (size_t)INT_MAX))
                next = (size_t)INT_MAX;
            if (next <= *capacity)
                return -1;
            tmp = realloc(*line, next);
            if (NULL == tmp)
                return -1;
            *line = tmp;
            *capacity = next;
        }

        avail_size = *capacity - used;
        avail = (avail_size > (size_t)INT_MAX) ? INT_MAX : (int)avail_size;
        z = fgets(*line + used, avail, fp);
        if (NULL == z)
        {
            if (used > 0u)
                break;
            if (feof(fp))
                return 0;
            return -1;
        }

        added = strlen(*line + used);
        used += added;
        if ((used > 0u) && ('\n' == (*line)[used - 1u]))
            break;
        if (feof(fp))
            break;
    }

    while ((used > 0u) &&
           (('\n' == (*line)[used - 1u]) || ('\r' == (*line)[used - 1u])))
        (*line)[--used] = '\0';

    return 1;
}

static int st_vcf_split_tabs(char *line, char ***field, size_t *capacity,
                             size_t *nfield)
{
    size_t n = 1u;
    size_t len;

    if ((NULL == line) || (NULL == field) || (NULL == capacity) ||
        (NULL == nfield))
        return -1;

    len = strlen(line);
    for (size_t i = 0u; i < len; i++)
        if ('\t' == line[i])
            n++;

    if (n > *capacity)
    {
        char **tmp;
        if (n > SIZE_MAX / sizeof(**field))
            return -1;
        tmp = realloc(*field, n * sizeof(**field));
        if (NULL == tmp)
            return -1;
        *field = tmp;
        *capacity = n;
    }

    n = 0u;
    (*field)[n++] = line;
    for (size_t i = 0u; i < len; i++)
        if ('\t' == line[i])
        {
            line[i] = '\0';
            (*field)[n++] = line + i + 1u;
        }

    *nfield = n;
    return 0;
}

static int st_vcf_format_gt_index(const char *format, size_t *gt_index)
{
    const char *start;
    const char *p;
    size_t index = 0u;
    int found = 0;

    if ((NULL == format) || ('\0' == format[0]) || (0 == strcmp(format, ".")) ||
        (NULL == gt_index))
        return -1;

    start = format;
    for (p = format;; p++)
        if ((':' == *p) || ('\0' == *p))
        {
            size_t len = (size_t)(p - start);
            if (0u == len)
                return -1;
            if ((2u == len) && ('G' == start[0]) && ('T' == start[1]))
            {
                if (found)
                    return -1;
                *gt_index = index;
                found = 1;
            }
            if ('\0' == *p)
                break;
            index++;
            start = p + 1;
        }

    return found ? 0 : -1;
}

static int st_vcf_subfield(const char *sample, size_t wanted,
                           const char **start_out, size_t *len_out)
{
    const char *start;
    const char *p;
    size_t index = 0u;

    if ((NULL == sample) || (NULL == start_out) || (NULL == len_out))
        return -1;

    start = sample;
    for (p = sample;; p++)
        if ((':' == *p) || ('\0' == *p))
        {
            if (index == wanted)
            {
                *start_out = start;
                *len_out = (size_t)(p - start);
                return 0;
            }
            if ('\0' == *p)
                break;
            index++;
            start = p + 1;
        }

    return -1;
}

static int st_vcf_parse_allele_index(const char *s, size_t n, gs_UNSIGNED n_alt,
                                     int *missing, gs_UNSIGNED *value)
{
    unsigned long v = 0ul;

    if ((NULL == s) || (0u == n) || (NULL == missing) || (NULL == value))
        return -1;

    if ((1u == n) && ('.' == s[0]))
    {
        *missing = 1;
        *value = 0u;
        return 0;
    }

    for (size_t i = 0u; i < n; i++)
    {
        unsigned int d;
        if ((s[i] < '0') || (s[i] > '9'))
            return -1;
        d = (unsigned int)(s[i] - '0');
        if (v > (ULONG_MAX - d) / 10ul)
            return -1;
        v = 10ul * v + d;
    }

    if ((gs_UNSIGNED)v > n_alt)
        return -1;

    *missing = 0;
    *value = (gs_UNSIGNED)v;
    return 0;
}

static int st_vcf_parse_gt(const char *gt, size_t len, gs_UNSIGNED n_alt,
                           GS_ALLELE_TYPE *a1, GS_ALLELE_TYPE *a2)
{
    size_t sep = SIZE_MAX;
    int miss1, miss2;
    gs_UNSIGNED v1, v2;

    if ((NULL == gt) || (NULL == a1) || (NULL == a2) || (0u == len))
        return -1;

    if ((1u == len) && ('.' == gt[0]))
    {
        *a1 = GS_ALLELE_MISSING;
        *a2 = GS_ALLELE_MISSING;
        return 0;
    }

    for (size_t i = 0u; i < len; i++)
        if (('/' == gt[i]) || ('|' == gt[i]))
        {
            if (SIZE_MAX != sep)
                return -1;
            sep = i;
        }

    if ((SIZE_MAX == sep) || (0u == sep) || (sep + 1u >= len))
        return -1;

    if ((0 != st_vcf_parse_allele_index(gt, sep, n_alt, &miss1, &v1)) ||
        (0 != st_vcf_parse_allele_index(gt + sep + 1u, len - sep - 1u, n_alt,
                                        &miss2, &v2)))
        return -1;

    if (miss1 || miss2)
    {
        *a1 = GS_ALLELE_MISSING;
        *a2 = GS_ALLELE_MISSING;
        return 0;
    }

    if ((v1 >= (gs_UNSIGNED)GS_ALLELE_MAX) ||
        (v2 >= (gs_UNSIGNED)GS_ALLELE_MAX))
        return -1;

    *a1 = (GS_ALLELE_TYPE)(v1 + 1u);
    *a2 = (GS_ALLELE_TYPE)(v2 + 1u);
    return 0;
}

static int st_vcf_parse_pos(const char *text, unsigned long long *pos)
{
    char *end = NULL;
    unsigned long long value;

    if ((NULL == text) || ('\0' == text[0]) || (NULL == pos) ||
        ('-' == text[0]) || ('+' == text[0]))
        return -1;

    errno = 0;
    value = strtoull(text, &end, 10);
    if ((0 != errno) || (end == text) || ('\0' != *end) || (0ull == value))
        return -1;

    *pos = value;
    return 0;
}

static int st_vcf_allele_label(gs_UNSIGNED index, char *buffer, size_t capacity)
{
    static const char base[4] = {'A', 'C', 'G', 'T'};
    char rev[32];
    size_t n = 0u;
    gs_UNSIGNED value;

    if ((NULL == buffer) || (capacity < 2u))
        return -1;

    if (index < 4u)
    {
        buffer[0] = base[index];
        buffer[1] = '\0';
        return 0;
    }

    value = index - 4u;
    do
    {
        if (n >= sizeof(rev))
            return -1;
        rev[n++] = base[value & 3u];
        value >>= 2u;
    } while (value > 0u);

    if (n + 2u > capacity)
        return -1;
    buffer[0] = 'A';
    for (size_t i = 0u; i < n; i++)
        buffer[1u + i] = rev[n - i - 1u];
    buffer[n + 1u] = '\0';
    return 0;
}

static int st_vcf_generated_marker_name(char *buffer, size_t capacity,
                                        gs_UNSIGNED index)
{
    int n;

    if ((NULL == buffer) || (0u == capacity))
        return -1;

    n = snprintf(buffer, capacity, "M%lu", (unsigned long)(index + 1u));
    if ((n < 1) || ((size_t)n >= capacity))
        return -1;
    return 0;
}

static int st_vcf_map_pool_add(size_t *total, const char *s)
{
    size_t n;

    if ((NULL == total) || (NULL == s))
        return -1;
    n = strlen(s) + 1u;
    if (*total > SIZE_MAX - n)
        return -1;
    *total += n;
    return 0;
}

typedef struct st_vcf_pedigree_record
{
    char *id, *pedigree;
    struct st_vcf_pedigree_record *next;
} st_vcf_pedigree_record;

static char *st_vcf_quoted_read(const char **cursor)
{
    const char *from = *cursor;
    if ((unsigned char)*from != 34)
        return NULL;
    from++;
    char *out = st_work_alloc(strlen(from) + 1, 1), *to = out;
    while (*from && (unsigned char)*from != 34)
    {
        if (*from == '\\')
        {
            from++;
            if (*from != '\\' && (unsigned char)*from != 34)
                return NULL;
        }
        if ((unsigned char)*from < 32)
            return NULL;
        *to++ = *from++;
    }
    if ((unsigned char)*from != 34)
        return NULL;
    *to = 0;
    *cursor = from + 1;
    return out;
}

static int st_vcf_pedigree_parse(const char *line,
                                 st_vcf_pedigree_record **records)
{
    const char *from = line + strlen("##SelectionToolsPedigree=<ID=");
    st_vcf_pedigree_record *record = st_work_alloc(1, sizeof(*record));
    record->id = st_vcf_quoted_read(&from);
    if (!record->id || !*record->id || strncmp(from, ",Pedigree=", 10))
        return -1;
    from += 10;
    record->pedigree = st_vcf_quoted_read(&from);
    if (!record->pedigree || !*record->pedigree || strcmp(from, ">"))
        return -1;
    record->next = *records;
    *records = record;
    return 0;
}

static int st_vcf_pedigree_apply(gs_varset_type *x,
                                 st_vcf_pedigree_record *records)
{
    if (!records)
        return 0;
    gs_name_index_type *index = NULL;
    if (gs_make_ind_name_index(x->IndName, x->NoInd, &index))
        return -1;
    char **strings = calloc(x->NoInd ? x->NoInd : 1, sizeof(*strings));
    if (!strings)
    {
        free(index);
        return -1;
    }
    gs_UNSIGNED count = 0;
    for (st_vcf_pedigree_record *record = records; record;
         record = record->next)
    {
        gs_UNSIGNED at = gs_name_index_find(index, x->NoInd, record->id);
        if (at == gs_UNSIGNED_MAX || strings[at])
            goto fail;
        strings[at] = st_text_dup(record->pedigree);
        if (!strings[at])
            goto fail;
        count++;
    }
    if (count != x->NoInd)
        goto fail;
    free(index);
    st_text_array_free(x->Pedigree, x->NoInd);
    x->Pedigree = strings;
    return 0;
fail:
    free(index);
    st_text_array_free(strings, x->NoInd);
    return -1;
}

static GS_ALLELE_TYPE st_vcf_actg_code(const char *label, size_t length)
{
    if (length == 3)
    {
        static const char *const abbreviation[4] = {"ADE", "CYT", "GUA", "THY"};
        for (int code = 0; code < 4; code++)
            if (toupper((unsigned char)label[0]) == abbreviation[code][0] &&
                toupper((unsigned char)label[1]) == abbreviation[code][1] &&
                toupper((unsigned char)label[2]) == abbreviation[code][2])
                return (GS_ALLELE_TYPE)(code + 1);
    }
    if (length != 1)
        return GS_ALLELE_MISSING;
    switch (label[0])
    {
    case 'A':
    case 'a':
        return 1;
    case 'C':
    case 'c':
        return 2;
    case 'G':
    case 'g':
        return 3;
    case 'T':
    case 't':
        return 4;
    default:
        return GS_ALLELE_MISSING;
    }
}

static int st_vcf_actg_codes(const char *ref, const char *alt,
                             GS_ALLELE_TYPE codes[256])
{
    memset(codes, GS_ALLELE_MISSING, 256 * sizeof(*codes));
    codes[1] = st_vcf_actg_code(ref, strlen(ref));
    int count = codes[1] != GS_ALLELE_MISSING;
    if (!strcmp(alt, "."))
        return count;
    unsigned int index = 2;
    for (const char *start = alt; *start; index++)
    {
        const char *end = strchr(start, ',');
        size_t length = end ? (size_t)(end - start) : strlen(start);
        codes[index] = st_vcf_actg_code(start, length);
        count += codes[index] != GS_ALLELE_MISSING;
        if (!end)
            break;
        start = end + 1;
    }
    return count;
}

static void st_vcf_actg_pair(GS_ALLELE_TYPE *a, GS_ALLELE_TYPE *b,
                             const GS_ALLELE_TYPE codes[256])
{
    *a = codes[*a];
    *b = codes[*b];
    if (*a == GS_ALLELE_MISSING || *b == GS_ALLELE_MISSING)
        *a = *b = GS_ALLELE_MISSING;
}

static void *st_vcf_shrink(void *data, size_t count, size_t size)
{
    void *smaller = realloc(data, count * size);
    return smaller ? smaller : data;
}

static int st_vcf_actg_compact(gs_varset_type *tmp, const unsigned char *keep)
{
    gs_UNSIGNED retained = 0;
    size_t pool_size = 0, offset = 0;
    int has_map = tmp->NLoc != 0;
    char *pool = NULL;
    for (gs_UNSIGNED m = 0; m < tmp->NoMar; m++)
        if (keep[m])
        {
            retained++;
            if (has_map &&
                (st_vcf_map_pool_add(&pool_size, tmp->map_S[m].name) ||
                 st_vcf_map_pool_add(&pool_size, tmp->map_S[m].class)))
                return -1;
        }
    if (!retained)
    {
        gs_info(tmp, -2, "No markers with A/C/G/T allele labels remain");
        return -1;
    }
    if (retained == tmp->NoMar)
        return 0;
    if (has_map)
    {
        pool = gs_malloc_array(pool_size, sizeof(*pool));
        if (!pool)
            return -1;
    }

    gs_UNSIGNED output = 0;
    int previous_chrom = 0, chrom = 0;
    tmp->MaxML = 0;
    for (gs_UNSIGNED m = 0; m < tmp->NoMar; m++)
        if (keep[m])
        {
            if (output != m)
            {
                size_t from = (size_t)m * tmp->NoInd;
                size_t to = (size_t)output * tmp->NoInd;
                size_t bytes = (size_t)tmp->NoInd * sizeof(GS_ALLELE_TYPE);
                memmove(tmp->allele_raw1 + to, tmp->allele_raw1 + from, bytes);
                memmove(tmp->allele_raw2 + to, tmp->allele_raw2 + from, bytes);
                memmove(tmp->MarName[output], tmp->MarName[m],
                        sizeof(*tmp->MarName));
            }
            size_t length = strlen(tmp->MarName[output]);
            if ((int)length > tmp->MaxML)
                tmp->MaxML = (int)length;
            if (has_map)
            {
                gs_mappoint_type point = tmp->map_S[m];
                if (point.chrom != previous_chrom)
                {
                    previous_chrom = point.chrom;
                    chrom++;
                }
                point.chrom = chrom;
                if (gs_copy_map_point_to_pool(&tmp->map_S[output], &point,
                                              output, pool, pool_size, &offset))
                {
                    free(pool);
                    return -1;
                }
            }
            output++;
        }
    if (has_map)
    {
        free(tmp->map_str_S);
        tmp->map_str_S = pool;
        tmp->map_str_N = pool_size;
        tmp->NLoc = retained;
        tmp->map_S = st_vcf_shrink(tmp->map_S, retained, sizeof(*tmp->map_S));
        tmp->map = st_vcf_shrink(tmp->map, retained, sizeof(*tmp->map));
        tmp->mdl = st_vcf_shrink(tmp->mdl, retained, sizeof(*tmp->mdl));
        for (gs_UNSIGNED m = 0; m < retained; m++)
        {
            tmp->map[m] = &tmp->map_S[m];
            tmp->mdl[m] = m;
        }
    }
    size_t elements = (size_t)retained * tmp->NoInd;
    tmp->allele_raw1 =
        st_vcf_shrink(tmp->allele_raw1, elements, sizeof(*tmp->allele_raw1));
    tmp->allele_raw2 =
        st_vcf_shrink(tmp->allele_raw2, elements, sizeof(*tmp->allele_raw2));
    tmp->MarName = st_vcf_shrink(tmp->MarName, retained, sizeof(*tmp->MarName));
    tmp->NoMar = retained;
    return 0;
}

SEXP st_read_vcf(SEXP filename, SEXP dataset, SEXP strictACTG)
{
    gs_varset_type *gv;
    gs_varset_type tmp;
    FILE *fp = NULL;
    char *line = NULL;
    size_t line_capacity = 0u;
    char **field = NULL;
    size_t field_capacity = 0u;
    char **sample_name = NULL;
    gs_UNSIGNED no_ind = 0u, no_mar = 0u;
    size_t map_string_size = 0u;
    double position_scale = ST_VCF_POSITION_SCALE;
    int have_header = 0;
    int failed = 0;
    int status = -2;
    unsigned char *actg_keep = NULL;
    int strict;
    gs_name_index_type *index = NULL;
    char **chrom_name = NULL;
    size_t chrom_capacity = 0u;
    gs_UNSIGNED no_chrom = 0u;
    size_t map_string_offset = 0u;
    unsigned long long previous_raw_pos = 0ull;
    int previous_chrom = 0;
    int saw_phased_gt = 0;
    int saw_unphased_gt = 0;
    st_vcf_pedigree_record *pedigrees = NULL;
    int import_pedigree =
        asLogical(getAttrib(filename, install("pedigree.enabled"))) != FALSE;

    if ((STRSXP != TYPEOF(filename)) || (1 != XLENGTH(filename)) ||
        (NA_STRING == STRING_ELT(filename, 0)) ||
        ('\0' == CHAR(STRING_ELT(filename, 0))[0]) ||
        (STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
        (NA_STRING == STRING_ELT(dataset, 0)) ||
        ('\0' == CHAR(STRING_ELT(dataset, 0))[0]))
        return ScalarInteger(-2);

    if (TYPEOF(strictACTG) != LGLSXP || XLENGTH(strictACTG) != 1 ||
        LOGICAL(strictACTG)[0] == NA_LOGICAL)
        return ScalarInteger(-2);
    strict = LOGICAL(strictACTG)[0];

    gv = gs_fdta((char *)CHAR(STRING_ELT(dataset, 0)));
    gs_tmp_init_from(&tmp, gv);

    fp = fopen(CHAR(STRING_ELT(filename, 0)), "rb");
    if (NULL == fp)
    {
        gs_info(gv, -2, "Cannot open VCF file");
        failed = 1;
        goto cleanup;
    }

    for (;;)
    {
        int lr = st_vcf_getline(fp, &line, &line_capacity);
        size_t nfield = 0u;

        if (0 == lr)
            break;
        if (lr < 0)
        {
            gs_info(gv, -2, "Error reading VCF file");
            failed = 1;
            goto cleanup;
        }
        if ('\0' == line[0])
            continue;

        if (0 == strncmp(line, "##SelectionToolsPositionScale=", 30))
        {
            char *end = NULL;
            double scale = strtod(line + 30, &end);
            if ((end != line + 30) && ('\0' == *end) && isfinite(scale) &&
                (scale > 0.0))
                position_scale = scale;
            continue;
        }

        if (import_pedigree &&
            !strncmp(line, "##SelectionToolsPedigree=<ID=",
                     strlen("##SelectionToolsPedigree=<ID=")))
        {
            if (st_vcf_pedigree_parse(line, &pedigrees))
            {
                gs_info(gv, -2, "Malformed VCF pedigree metadata");
                failed = 1;
                goto cleanup;
            }
            continue;
        }
        if (0 == strncmp(line, "##", 2))
            continue;

        if ('#' == line[0])
        {
            if (have_header)
            {
                gs_info(gv, -2, "VCF contains more than one column header");
                failed = 1;
                goto cleanup;
            }
            if (0 != st_vcf_split_tabs(line, &field, &field_capacity, &nfield))
            {
                failed = 1;
                goto cleanup;
            }
            if ((nfield < 10u) || (0 != strcmp(field[0], "#CHROM")) ||
                (0 != strcmp(field[1], "POS")) ||
                (0 != strcmp(field[2], "ID")) ||
                (0 != strcmp(field[3], "REF")) ||
                (0 != strcmp(field[4], "ALT")) ||
                (0 != strcmp(field[5], "QUAL")) ||
                (0 != strcmp(field[6], "FILTER")) ||
                (0 != strcmp(field[7], "INFO")) ||
                (0 != strcmp(field[8], "FORMAT")))
            {
                gs_info(gv, -2,
                        "VCF header must contain standard columns through "
                        "FORMAT and at least one sample");
                failed = 1;
                goto cleanup;
            }
            if (nfield - 9u > (size_t)gs_UNSIGNED_MAX)
            {
                gs_info(gv, -2, "Too many VCF samples");
                failed = 1;
                goto cleanup;
            }
            no_ind = (gs_UNSIGNED)(nfield - 9u);
            sample_name = calloc((size_t)no_ind, sizeof(*sample_name));
            if (NULL == sample_name)
            {
                gs_info(gv, -2, "Memory allocation error");
                failed = 1;
                goto cleanup;
            }
            for (gs_UNSIGNED i = 0u; i < no_ind; i++)
            {
                const char *name = field[9u + i];
                if (!st_identifier_valid(name) || (strlen(name) >= INT_MAX))
                {
                    gs_info(gv, -2, "Invalid or overlong VCF sample name");
                    failed = 1;
                    goto cleanup;
                }
                sample_name[i] = st_vcf_strdup(name);
                if (NULL == sample_name[i])
                {
                    gs_info(gv, -2, "Memory allocation error");
                    failed = 1;
                    goto cleanup;
                }
            }
            have_header = 1;
            continue;
        }

        if (!have_header)
        {
            gs_info(gv, -2, "VCF data encountered before #CHROM header");
            failed = 1;
            goto cleanup;
        }

        if (0 != st_vcf_split_tabs(line, &field, &field_capacity, &nfield))
        {
            failed = 1;
            goto cleanup;
        }
        if (nfield != (size_t)no_ind + 9u)
        {
            gs_info(gv, -2, "VCF record has the wrong number of columns");
            failed = 1;
            goto cleanup;
        }
        {
            unsigned long long pos;
            gs_UNSIGNED n_alt;
            size_t gt_index;
            char generated[64];
            const char *id = field[2];

            if (('\0' == field[0][0]) ||
                (0 != st_vcf_parse_pos(field[1], &pos)) ||
                ('\0' == field[3][0]) || (0 == strcmp(field[3], ".")) ||
                (0 != st_stvcf_alleles(field[3], field[4], &n_alt)) ||
                (0 != st_vcf_format_gt_index(field[8], &gt_index)))
            {
                gs_info(gv, -2,
                        "Invalid CHROM/POS/REF/ALT/FORMAT field in VCF");
                failed = 1;
                goto cleanup;
            }

            if (0 == strcmp(id, "."))
            {
                if (0 != st_vcf_generated_marker_name(
                             generated, sizeof(generated), no_mar))
                {
                    failed = 1;
                    goto cleanup;
                }
                id = generated;
            }
            else if (!st_identifier_valid(id) || (strlen(id) >= gs_MAXNME))
            {
                gs_info(gv, -2, "Invalid or overlong VCF marker ID");
                failed = 1;
                goto cleanup;
            }

            if ((0 != st_vcf_map_pool_add(&map_string_size, id)) ||
                (0 != st_vcf_map_pool_add(&map_string_size, "m")))
            {
                gs_info(gv, -2, "VCF map string storage too large");
                failed = 1;
                goto cleanup;
            }
        }

        if (no_mar == gs_UNSIGNED_MAX)
        {
            gs_info(gv, -2, "Too many VCF markers");
            failed = 1;
            goto cleanup;
        }
        no_mar++;
    }

    if (!have_header || (0u == no_ind) || (0u == no_mar))
    {
        gs_info(gv, -2, "VCF contains no samples or no markers");
        failed = 1;
        goto cleanup;
    }

    if (0 != fclose(fp))
    {
        fp = NULL;
        gs_info(gv, -2, "Error closing VCF file");
        failed = 1;
        goto cleanup;
    }
    fp = NULL;

    tmp.NoInd = no_ind;
    tmp.NoMar = no_mar;
    tmp.NLoc = no_mar;
    tmp.map_str_N = map_string_size;
    tmp.IndName = calloc((size_t)no_ind, sizeof(*tmp.IndName));
    tmp.MarName = gs_malloc_array((size_t)no_mar, sizeof(*tmp.MarName));
    tmp.map_S = gs_malloc_array((size_t)no_mar, sizeof(*tmp.map_S));
    tmp.map_str_S = gs_malloc_array(tmp.map_str_N, sizeof(*tmp.map_str_S));
    tmp.map = gs_malloc_array((size_t)no_mar, sizeof(*tmp.map));
    tmp.mdl = gs_malloc_array((size_t)no_mar, sizeof(*tmp.mdl));
    if (strict)
        actg_keep = gs_malloc_array((size_t)no_mar, sizeof(*actg_keep));
    if ((strict && !actg_keep) || (NULL == tmp.IndName) ||
        (NULL == tmp.MarName) || (NULL == tmp.map_S) ||
        (NULL == tmp.map_str_S) || (NULL == tmp.map) || (NULL == tmp.mdl) ||
        (0 != gs_alloc_genotype_matrix(&tmp, no_mar, no_ind)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED i = 0u; i < no_ind; i++)
    {
        size_t len = strlen(sample_name[i]);
        if (st_identity_set(&tmp, i, sample_name[i], NULL))
        {
            failed = 1;
            goto cleanup;
        }
        if ((int)len > tmp.MaxIL)
            tmp.MaxIL = (int)len;
    }
    if ((0 != gs_make_ind_name_index(tmp.IndName, tmp.NoInd, &index)) ||
        gs_name_index_has_duplicates(index, tmp.NoInd))
    {
        gs_info(gv, -2, "Duplicate or invalid VCF sample names");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    fp = fopen(CHAR(STRING_ELT(filename, 0)), "rb");
    if (NULL == fp)
    {
        gs_info(gv, -2, "Cannot reopen VCF file");
        failed = 1;
        goto cleanup;
    }

    {
        gs_UNSIGNED marker = 0u;
        int saw_header = 0;

        for (;;)
        {
            int lr = st_vcf_getline(fp, &line, &line_capacity);
            size_t nfield = 0u;

            if (0 == lr)
                break;
            if (lr < 0)
            {
                gs_info(gv, -2, "Error reading VCF file");
                failed = 1;
                goto cleanup;
            }
            if ('\0' == line[0])
                continue;
            if (0 == strncmp(line, "##", 2))
                continue;
            if ('#' == line[0])
            {
                saw_header = 1;
                continue;
            }
            if (!saw_header || (marker >= no_mar) ||
                (0 !=
                 st_vcf_split_tabs(line, &field, &field_capacity, &nfield)) ||
                (nfield != (size_t)no_ind + 9u))
            {
                gs_info(gv, -2, "Invalid VCF structure on second pass");
                failed = 1;
                goto cleanup;
            }

            {
                unsigned long long raw_pos;
                gs_UNSIGNED n_alt;
                size_t gt_index;
                char generated[64];
                const char *id = field[2];
                int chrom_code;
                double genetic_pos;
                gs_mappoint_type source_point;
                GS_ALLELE_TYPE actg[256];

                if ((0 != st_vcf_parse_pos(field[1], &raw_pos)) ||
                    (0 != st_stvcf_alleles(field[3], field[4], &n_alt)) ||
                    (0 != st_vcf_format_gt_index(field[8], &gt_index)))
                {
                    gs_info(gv, -2, "Invalid VCF record on second pass");
                    failed = 1;
                    goto cleanup;
                }

                if (0 == strcmp(id, "."))
                {
                    if (0 != st_vcf_generated_marker_name(
                                 generated, sizeof(generated), marker))
                    {
                        failed = 1;
                        goto cleanup;
                    }
                    id = generated;
                }

                if ((0u == no_chrom) ||
                    (0 != strcmp(chrom_name[no_chrom - 1u], field[0])))
                {
                    for (gs_UNSIGNED c = 0u; c < no_chrom; c++)
                        if (0 == strcmp(chrom_name[c], field[0]))
                        {
                            gs_info(
                                gv, -2,
                                "VCF chromosome records must be contiguous");
                            failed = 1;
                            goto cleanup;
                        }

                    if ((size_t)no_chrom >= chrom_capacity)
                    {
                        size_t next =
                            (0u == chrom_capacity) ? 8u : 2u * chrom_capacity;
                        char **z;
                        if ((next <= chrom_capacity) ||
                            (next > SIZE_MAX / sizeof(*chrom_name)))
                        {
                            failed = 1;
                            goto cleanup;
                        }
                        z = realloc(chrom_name, next * sizeof(*chrom_name));
                        if (NULL == z)
                        {
                            failed = 1;
                            goto cleanup;
                        }
                        chrom_name = z;
                        chrom_capacity = next;
                    }
                    chrom_name[no_chrom] = st_vcf_strdup(field[0]);
                    if (NULL == chrom_name[no_chrom])
                    {
                        failed = 1;
                        goto cleanup;
                    }
                    no_chrom++;
                    previous_raw_pos = 0ull;
                    previous_chrom = (int)no_chrom;
                }

                chrom_code = (int)no_chrom;
                if ((raw_pos < previous_raw_pos) &&
                    (chrom_code == previous_chrom))
                {
                    gs_info(gv, -2,
                            "VCF positions must be nondecreasing within "
                            "chromosome");
                    failed = 1;
                    goto cleanup;
                }
                previous_raw_pos = raw_pos;
                previous_chrom = chrom_code;

                genetic_pos = (position_scale > 0.0)
                                  ? (((double)raw_pos - 1.0) / position_scale)
                                  : (double)raw_pos;
                if (!isfinite(genetic_pos) || (genetic_pos < 0.0))
                {
                    gs_info(gv, -2,
                            "VCF position cannot be represented as a genetic "
                            "map position");
                    failed = 1;
                    goto cleanup;
                }

                strcpy(tmp.MarName[marker], id);
                if ((int)strlen(id) > tmp.MaxML)
                    tmp.MaxML = (int)strlen(id);
                source_point.chrom = chrom_code;
                source_point.pos = genetic_pos;
                source_point.name = id;
                source_point.class = "m";
                source_point.nr = marker;
                if (0 != gs_copy_map_point_to_pool(
                             &(tmp.map_S[marker]), &source_point, marker,
                             tmp.map_str_S, tmp.map_str_N, &map_string_offset))
                {
                    gs_info(gv, -2, "Invalid VCF map string storage");
                    failed = 1;
                    goto cleanup;
                }
                tmp.map[marker] = &(tmp.map_S[marker]);
                tmp.mdl[marker] = marker;

                if (strict)
                    actg_keep[marker] =
                        st_vcf_actg_codes(field[3], field[4], actg) > 0;

                for (gs_UNSIGNED i = 0u; i < no_ind; i++)
                {
                    const char *gt = NULL;
                    size_t gt_len = 0u;
                    GS_ALLELE_TYPE a1, a2;
                    size_t idx = (size_t)marker * (size_t)no_ind + i;

                    if ((0 != st_vcf_subfield(field[9u + i], gt_index, &gt,
                                              &gt_len)) ||
                        (0u == gt_len) ||
                        (0 != st_vcf_parse_gt(gt, gt_len, n_alt, &a1, &a2)))
                    {
                        gs_info(gv, -2, "Invalid diploid GT field in VCF");
                        failed = 1;
                        goto cleanup;
                    }
                    if ((!strict || actg_keep[marker]) &&
                        NULL != memchr(gt, '|', gt_len))
                        saw_phased_gt = 1;
                    if ((!strict || actg_keep[marker]) &&
                        NULL != memchr(gt, '/', gt_len))
                        saw_unphased_gt = 1;
                    if (strict)
                        st_vcf_actg_pair(&a1, &a2, actg);
                    tmp.allele_raw1[idx] = a1;
                    tmp.allele_raw2[idx] = a2;
                }
            }
            marker++;
        }

        if (marker != no_mar)
        {
            gs_info(gv, -2, "VCF marker count changed between passes");
            failed = 1;
            goto cleanup;
        }
    }

    tmp.phased = (saw_phased_gt && !saw_unphased_gt) ? 1 : 0;

    if (map_string_offset != tmp.map_str_N)
    {
        gs_info(gv, -2, "Invalid VCF map string storage");
        failed = 1;
        goto cleanup;
    }

    if (0 != st_stvcf_separate_map_ties(gv, tmp.map, tmp.NLoc))
    {
        failed = 1;
        goto cleanup;
    }

    if ((0 != gs_make_name_index(tmp.MarName, tmp.NoMar, &index)) ||
        gs_name_index_has_duplicates(index, tmp.NoMar))
    {
        gs_info(gv, -2, "Duplicate VCF marker IDs");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    {
        if (strict && st_vcf_actg_compact(&tmp, actg_keep))
        {
            failed = 1;
            goto cleanup;
        }
        int srv = 0;
        char *empty = "";
        char **empty_p = &empty;
        int zero = 0;
        int *zero_p = &zero;

        gs_marker_stats_01(&tmp, empty_p, empty_p, empty_p, zero_p, zero_p,
                           zero_p, zero_p, &srv);
        if (-2 == srv)
        {
            failed = 1;
            goto cleanup;
        }
    }

    if (st_vcf_pedigree_apply(&tmp, pedigrees) ||
        st_marker_pedigree_import(&tmp,
                                  getAttrib(filename, install("pedigree"))))
    {
        gs_info(gv, -2, "Invalid or inconsistent pedigree metadata");
        failed = 1;
        goto cleanup;
    }
    if (gs_tmp_commit(gv, &tmp))
    {
        failed = 1;
        goto cleanup;
    }
    status = 0;

cleanup:
    free(actg_keep);
    if (NULL != fp)
        fclose(fp);
    free(line);
    free(field);
    free(index);
    if (NULL != sample_name)
    {
        for (gs_UNSIGNED i = 0u; i < no_ind; i++)
            free(sample_name[i]);
        free(sample_name);
    }
    if (NULL != chrom_name)
    {
        for (gs_UNSIGNED c = 0u; c < no_chrom; c++)
            free(chrom_name[c]);
        free(chrom_name);
    }
    if (failed)
        gs_reset(&tmp);
    return ScalarInteger(status);
}

typedef struct
{
    FILE *plain;
} st_vcf_output;

static int st_vcf_output_open(st_vcf_output *out, const char *filename)
{
    memset(out, 0, sizeof(*out));
    out->plain = fopen(filename, "wb");
    return (NULL == out->plain) ? -1 : 0;
}

static int st_vcf_output_write(st_vcf_output *out, const char *text, size_t n)
{
    return (fwrite(text, 1u, n, out->plain) == n) ? 0 : -1;
}

static int st_vcf_output_printf(st_vcf_output *out, const char *format, ...)
{
    char stack[4096];
    char *buffer = stack;
    size_t capacity = sizeof(stack);
    va_list ap;
    int n;
    int rc;

    va_start(ap, format);
    n = vsnprintf(buffer, capacity, format, ap);
    va_end(ap);
    if (n < 0)
        return -1;

    if ((size_t)n >= capacity)
    {
        capacity = (size_t)n + 1u;
        buffer = malloc(capacity);
        if (NULL == buffer)
            return -1;
        va_start(ap, format);
        n = vsnprintf(buffer, capacity, format, ap);
        va_end(ap);
        if ((n < 0) || ((size_t)n >= capacity))
        {
            free(buffer);
            return -1;
        }
    }

    rc = st_vcf_output_write(out, buffer, (size_t)n);
    if (buffer != stack)
        free(buffer);
    return rc;
}

static int st_vcf_output_close(st_vcf_output *out)
{
    int rc = (NULL == out->plain) ? 0 : fclose(out->plain);
    out->plain = NULL;
    return (0 == rc) ? 0 : -1;
}

static int st_vcf_quoted_write(st_vcf_output *out, const char *text)
{
    const char quote = 34, slash = 92;
    if (st_vcf_output_write(out, &quote, 1))
        return -1;
    for (const char *p = text; *p; p++)
    {
        if ((unsigned char)*p < 32)
            return -1;
        if (((unsigned char)*p == 34 || (unsigned char)*p == 92) &&
            st_vcf_output_write(out, &slash, 1))
            return -1;
        if (st_vcf_output_write(out, p, 1))
            return -1;
    }
    return st_vcf_output_write(out, &quote, 1);
}

SEXP st_write_vcf(SEXP filename, SEXP dataset)
{
    gs_varset_type *gv;
    st_vcf_output out;
    gs_UNSIGNED *order = NULL;
    gs_mappoint_type **point_by_marker = NULL;
    unsigned char *seen = NULL;
    GS_ALLELE_TYPE *alleles = NULL;
    gs_UNSIGNED norder = 0u;
    long long previous_bp = 0ll;
    int previous_chrom = -1;
    int opened = 0;
    int failed = 0;
    int status = -2;

    memset(&out, 0, sizeof(out));

    if ((STRSXP != TYPEOF(filename)) || (1 != XLENGTH(filename)) ||
        (NA_STRING == STRING_ELT(filename, 0)) ||
        ('\0' == CHAR(STRING_ELT(filename, 0))[0]) ||
        (STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
        (NA_STRING == STRING_ELT(dataset, 0)) ||
        ('\0' == CHAR(STRING_ELT(dataset, 0))[0]))
        return ScalarInteger(-2);

    gv = gs_fdta((char *)CHAR(STRING_ELT(dataset, 0)));
    if ((0u == gv->NoInd) || (0u == gv->NoMar) || (NULL == gv->IndName) ||
        (NULL == gv->MarName) || !st_allele_present(gv) || (0u == gv->NLoc) ||
        (NULL == gv->map) || (NULL == gv->mdl))
    {
        gs_info(gv, -2,
                "No complete marker data/linkage map available for VCF output");
        return ScalarInteger(-2);
    }

    order = gs_malloc_array((size_t)gv->NoMar, sizeof(*order));
    point_by_marker =
        gs_malloc_array((size_t)gv->NoMar, sizeof(*point_by_marker));
    seen = calloc((size_t)gv->NoMar, sizeof(*seen));
    alleles = gs_malloc_array((size_t)gs_MAXALL, sizeof(*alleles));
    if ((NULL == order) || (NULL == point_by_marker) || (NULL == seen) ||
        (NULL == alleles))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED l = 0u; l < gv->NLoc; l++)
    {
        gs_UNSIGNED m = gv->mdl[l];
        if (gs_UNSIGNED_MAX == m)
            continue;
        if ((m >= gv->NoMar) || (NULL == gv->map[l]) || seen[m] ||
            (0 != strcmp(gv->map[l]->name, gv->MarName[m])) ||
            (gv->map[l]->chrom < 1) || !isfinite(gv->map[l]->pos) ||
            (gv->map[l]->pos < 0.0))
        {
            gs_info(
                gv, -2,
                "Marker data do not have one valid map position per marker");
            failed = 1;
            goto cleanup;
        }
        seen[m] = 1u;
        point_by_marker[m] = gv->map[l];
        order[norder++] = m;
    }
    if (norder != gv->NoMar)
    {
        gs_info(gv, -2, "Marker data do not have one map position per marker");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED i = 0u; i < gv->NoInd; i++)
        if (!st_identifier_valid(gv->IndName[i]))
        {
            gs_info(gv, -2, "Invalid sample name for VCF output");
            failed = 1;
            goto cleanup;
        }
    for (gs_UNSIGNED m = 0u; m < gv->NoMar; m++)
        if (!st_identifier_valid(gv->MarName[m]) ||
            (strlen(gv->MarName[m]) >= gs_MAXNME))
        {
            gs_info(gv, -2, "Invalid marker name for VCF output");
            failed = 1;
            goto cleanup;
        }

    if (0 != st_vcf_output_open(&out, CHAR(STRING_ELT(filename, 0))))
    {
        gs_info(gv, -2, "Cannot open VCF output file");
        failed = 1;
        goto cleanup;
    }
    opened = 1;

    if ((0 != st_vcf_output_printf(&out, "##fileformat=VCFv4.3\n")) ||
        (0 != st_vcf_output_printf(&out, "##source=SelectionTools\n")) ||
        (0 != st_vcf_output_printf(
                  &out, "##SelectionToolsPositionScale=1000000\n")) ||
        (0 != st_vcf_output_printf(&out, "##FORMAT=<ID=GT,Number=1,Type=String,"
                                         "Description=\"Genotype\">\n")))
    {
        failed = 1;
        goto cleanup;
    }
    if (asLogical(getAttrib(filename, install("pedigree.enabled"))) != FALSE)
        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        {
            const char *pedigree = gv->Pedigree && gv->Pedigree[i]
                                       ? gv->Pedigree[i]
                                       : gv->IndName[i];
            if (st_vcf_output_printf(&out, "##SelectionToolsPedigree=<ID=") ||
                st_vcf_quoted_write(&out, gv->IndName[i]) ||
                st_vcf_output_printf(&out, ",Pedigree=") ||
                st_vcf_quoted_write(&out, pedigree) ||
                st_vcf_output_printf(&out, ">\n"))
            {
                failed = 1;
                goto cleanup;
            }
        }
    if (st_vcf_output_printf(
            &out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"))
    {
        failed = 1;
        goto cleanup;
    }
    for (gs_UNSIGNED i = 0u; i < gv->NoInd; i++)
        if (0 != st_vcf_output_printf(&out, "\t%s", gv->IndName[i]))
        {
            failed = 1;
            goto cleanup;
        }
    if (0 != st_vcf_output_printf(&out, "\n"))
    {
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED k = 0u; k < norder; k++)
    {
        gs_UNSIGNED m = order[k];
        gs_mappoint_type *point = point_by_marker[m];
        gs_UNSIGNED n_alleles = 0u;
        int has_missing = 0;
        long long bp;

        if (NULL == point)
        {
            failed = 1;
            goto cleanup;
        }

        if (0 != st_markerdata_observed_alleles(gv, m, alleles, &n_alleles,
                                                &has_missing))
        {
            gs_info(gv, -2, "Invalid marker alleles for VCF output");
            failed = 1;
            goto cleanup;
        }
        (void)has_missing;

        if (point->pos > ((double)LLONG_MAX - 2.0) / ST_VCF_POSITION_SCALE)
        {
            gs_info(gv, -2,
                    "Genetic position is too large for Beagle VCF scaling");
            failed = 1;
            goto cleanup;
        }
        bp = llround(point->pos * ST_VCF_POSITION_SCALE) + 1ll;
        if ((point->chrom == previous_chrom) && (bp <= previous_bp))
            bp = previous_bp + 1ll;
        if (bp < 1ll)
        {
            failed = 1;
            goto cleanup;
        }
        previous_chrom = point->chrom;
        previous_bp = bp;

        {
            char label[40];
            if ((0 != st_vcf_allele_label(0u, label, sizeof(label))) ||
                (0 != st_vcf_output_printf(&out, "%d\t%lld\t%s\t%s\t",
                                           point->chrom, bp, gv->MarName[m],
                                           label)))
            {
                failed = 1;
                goto cleanup;
            }
        }
        if (n_alleles <= 1u)
        {
            if (0 != st_vcf_output_printf(&out, "."))
            {
                failed = 1;
                goto cleanup;
            }
        }
        else
            for (gs_UNSIGNED a = 1u; a < n_alleles; a++)
            {
                char label[40];
                if ((0 != st_vcf_allele_label(a, label, sizeof(label))) ||
                    (0 != st_vcf_output_printf(&out, (1u == a) ? "%s" : ",%s",
                                               label)))
                {
                    failed = 1;
                    goto cleanup;
                }
            }

        if (0 != st_vcf_output_printf(&out, "\t.\tPASS\t.\tGT"))
        {
            failed = 1;
            goto cleanup;
        }

        for (gs_UNSIGNED i = 0u; i < gv->NoInd; i++)
        {
            size_t idx = (size_t)m * (size_t)gv->NoInd + i;
            GS_ALLELE_TYPE a1 = st_allele_get(gv, idx, 0);
            GS_ALLELE_TYPE a2 = st_allele_get(gv, idx, 1);

            if ((GS_ALLELE_MISSING == a1) || (GS_ALLELE_MISSING == a2))
            {
                if (0 !=
                    st_vcf_output_printf(&out, gv->phased ? "\t.|." : "\t./."))
                {
                    failed = 1;
                    goto cleanup;
                }
            }
            else
            {
                int i1 = st_markerdata_allele_index(alleles, n_alleles, a1);
                int i2 = st_markerdata_allele_index(alleles, n_alleles, a2);
                if ((i1 < 0) || (i2 < 0) ||
                    (0 !=
                     st_vcf_output_printf(
                         &out, gv->phased ? "\t%d|%d" : "\t%d/%d", i1, i2)))
                {
                    gs_info(gv, -2,
                            "Allele normalization failed during VCF output");
                    failed = 1;
                    goto cleanup;
                }
            }
        }
        if (0 != st_vcf_output_printf(&out, "\n"))
        {
            failed = 1;
            goto cleanup;
        }
    }

    if (0 != st_vcf_output_close(&out))
    {
        opened = 0;
        gs_info(gv, -2, "Error closing VCF output file");
        failed = 1;
        goto cleanup;
    }
    opened = 0;
    status = 0;

cleanup:
    if (opened)
        st_vcf_output_close(&out);
    free(order);
    free(point_by_marker);
    free(seen);
    free(alleles);
    if (failed && (0 == status))
        status = -2;
    return ScalarInteger(status);
}

static int st_physical_map_valid(const st_data_value *map, SEXP ids)
{
    const st_data_value *names = st_data_get(map, "marker"),
                        *chrom = st_data_get(map, "chromosome"),
                        *pos = st_data_get(map, "position"),
                        *unit = st_data_attribute_get(map, "position.unit");
    if (!names || names->kind != ST_DATA_STRING || !chrom || !pos ||
        chrom->length != names->length || pos->length != names->length ||
        (pos->kind != ST_DATA_REAL && pos->kind != ST_DATA_INTEGER) ||
        (chrom->kind != ST_DATA_STRING && chrom->kind != ST_DATA_INTEGER &&
         chrom->kind != ST_DATA_REAL) || !unit || unit->length != 1 ||
        !st_data_string(unit, 0) || strcmp(st_data_string(unit, 0), "bp"))
        return 0;
    size_t n = names->length;
    if (n > gs_UNSIGNED_MAX)
        return 0;
    gs_name_index_type *index = st_work_alloc(n, sizeof(*index));
    for (size_t i = 0; i < n; i++)
    {
        const char *name = st_data_string(names, i);
        double value = pos->kind == ST_DATA_REAL ? pos->data.real[i]
                                                : pos->data.integer[i];
        if (!name || !name[0] || !R_FINITE(value) || value < 0)
            return 0;
        if (chrom->kind == ST_DATA_STRING)
        {
            const char *label = st_data_string(chrom, i);
            if (!label || !label[0])
                return 0;
        }
        else
        {
            double label = chrom->kind == ST_DATA_REAL ? chrom->data.real[i]
                                                       : chrom->data.integer[i];
            if (!R_FINITE(label) || label <= 0)
                return 0;
        }
        index[i].name = name;
        index[i].idx = (gs_UNSIGNED)i;
    }
    qsort(index, n, sizeof(*index), gs_name_index_cmp);
    if (gs_name_index_has_duplicates(index, (gs_UNSIGNED)n))
        return 0;
    for (R_xlen_t i = 0; i < XLENGTH(ids); i++)
        if (STRING_ELT(ids, i) == NA_STRING ||
            gs_name_index_find(index, (gs_UNSIGNED)n,
                               CHAR(STRING_ELT(ids, i))) == gs_UNSIGNED_MAX)
            return 0;
    return 1;
}

SEXP st_load_vcf_data(SEXP stvcf, SEXP dataset, SEXP strictACTG)
{
    gs_varset_type *gv;
    gs_varset_type tmp;
    SEXP chrom_s, pos_s, id_s, ref_s, alt_s, ind_s;
    SEXP allele1_s, allele2_s, chrom_levels_s, allele_codes;
    gs_mappoint_type *map_s = NULL;
    gs_mappoint_type **map = NULL;
    gs_UNSIGNED *n_alt = NULL;
    gs_name_index_type *index = NULL;
    unsigned char *chrom_seen = NULL;
    unsigned char *actg_keep = NULL;
    int strict;
    R_xlen_t no_mar_x, no_ind_x;
    gs_UNSIGNED no_mar, no_ind, no_chrom;
    int input_phased = 0;
    int has_map = 0;
    int failed = 0;
    int status = -2;

    if ((STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
        (NA_STRING == STRING_ELT(dataset, 0)))
        return ScalarInteger(-2);

    if (TYPEOF(strictACTG) != LGLSXP || XLENGTH(strictACTG) != 1 ||
        LOGICAL(strictACTG)[0] == NA_LOGICAL)
        return ScalarInteger(-2);
    strict = LOGICAL(strictACTG)[0];

    gv = gs_fdta((char *)translateCharUTF8(STRING_ELT(dataset, 0)));
    gs_tmp_init_from(&tmp, gv);

    if (!st_stvcf_components(stvcf))
    {
        gs_info(gv, -2,
                "STvcf must contain exactly the nine required components");
        failed = 1;
        goto cleanup;
    }
    if (0 != st_stvcf_phased_attribute(stvcf, &input_phased))
    {
        gs_info(gv, -2, "Invalid STvcf phased attribute");
        failed = 1;
        goto cleanup;
    }

    chrom_s = st_stvcf_element(stvcf, "CHROM");
    pos_s = st_stvcf_element(stvcf, "POS");
    id_s = st_stvcf_element(stvcf, "ID");
    ref_s = st_stvcf_element(stvcf, "REF");
    alt_s = st_stvcf_element(stvcf, "ALT");
    ind_s = st_stvcf_element(stvcf, "individual");
    allele1_s = st_stvcf_element(stvcf, "allele1");
    allele2_s = st_stvcf_element(stvcf, "allele2");
    chrom_levels_s = st_stvcf_element(stvcf, "CHROM.levels");
    allele_codes = getAttrib(stvcf, install("allele.codes"));

    if ((INTSXP != TYPEOF(chrom_s)) || (REALSXP != TYPEOF(pos_s)) ||
        (STRSXP != TYPEOF(id_s)) || (STRSXP != TYPEOF(ref_s)) ||
        (STRSXP != TYPEOF(alt_s)) || (STRSXP != TYPEOF(ind_s)) ||
        (STRSXP != TYPEOF(chrom_levels_s)) ||
        !(((RAWSXP == TYPEOF(allele1_s)) && (RAWSXP == TYPEOF(allele2_s))) ||
          ((INTSXP == TYPEOF(allele1_s)) && (INTSXP == TYPEOF(allele2_s)))))
    {
        gs_info(gv, -2, "Invalid STvcf component types");
        failed = 1;
        goto cleanup;
    }

    no_mar_x = XLENGTH(id_s);
    no_ind_x = XLENGTH(ind_s);
    if ((0 == no_mar_x) || (0 == no_ind_x) ||
        ((R_xlen_t)(gs_UNSIGNED)no_mar_x != no_mar_x) ||
        ((R_xlen_t)(gs_UNSIGNED)no_ind_x != no_ind_x))
    {
        gs_info(gv, -2, "Invalid STvcf dimensions");
        failed = 1;
        goto cleanup;
    }

    no_mar = (gs_UNSIGNED)no_mar_x;
    no_ind = (gs_UNSIGNED)no_ind_x;
    no_chrom = (gs_UNSIGNED)XLENGTH(chrom_levels_s);
    has_map = no_chrom > 0;
    tmp.phased = input_phased;

    if ((XLENGTH(chrom_s) != no_mar_x) || (XLENGTH(pos_s) != no_mar_x) ||
        (XLENGTH(ref_s) != no_mar_x) || (XLENGTH(alt_s) != no_mar_x) ||
        !st_stvcf_matrix_dim(allele1_s, no_mar_x, no_ind_x) ||
        !st_stvcf_matrix_dim(allele2_s, no_mar_x, no_ind_x))
    {
        gs_info(gv, -2, "Inconsistent STvcf component lengths or dimensions");
        failed = 1;
        goto cleanup;
    }

    SEXP physical = getAttrib(stvcf, install("physical.map"));
    if (!isNull(physical))
    {
        st_data_value *value = st_data_from_R(physical);
        if (!st_physical_map_valid(value, id_s))
        {
            gs_info(gv, -2, "Invalid physical map or marker alignment");
            failed = 1;
            goto cleanup;
        }
        tmp.physical_map = value;
        st_work_detach(value);
    }

    map_s = gs_malloc_array((size_t)no_mar, sizeof(*map_s));
    map = gs_malloc_array((size_t)no_mar, sizeof(*map));
    n_alt = gs_malloc_array((size_t)no_mar, sizeof(*n_alt));
    chrom_seen = calloc(has_map ? (size_t)no_chrom : 1u, sizeof(*chrom_seen));
    if (strict)
        actg_keep = gs_malloc_array((size_t)no_mar, sizeof(*actg_keep));
    if ((strict && !actg_keep) || (NULL == map_s) || (NULL == map) ||
        (NULL == n_alt) || (NULL == chrom_seen))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED c = 0; c < no_chrom; c++)
    {
        SEXP s = STRING_ELT(chrom_levels_s, c);

        if ((NA_STRING == s) || ('\0' == translateCharUTF8(s)[0]))
        {
            gs_info(gv, -2, "Invalid CHROM.levels in STvcf");
            failed = 1;
            goto cleanup;
        }

        for (gs_UNSIGNED d = 0; d < c; d++)
            if (0 == strcmp(translateCharUTF8(s),
                            translateCharUTF8(STRING_ELT(chrom_levels_s, d))))
            {
                gs_info(gv, -2, "Duplicate chromosome label in STvcf");
                failed = 1;
                goto cleanup;
            }
    }

    for (gs_UNSIGNED m = 0; m < no_mar; m++)
    {
        int chrom = INTEGER(chrom_s)[m];
        double pos = REAL(pos_s)[m];
        SEXP id = STRING_ELT(id_s, m);
        SEXP ref = STRING_ELT(ref_s, m);
        SEXP alt = STRING_ELT(alt_s, m);
        size_t id_len;

        if ((has_map && ((NA_INTEGER == chrom) || (chrom < 1) ||
                         ((gs_UNSIGNED)chrom > no_chrom) || !isfinite(pos) ||
                         (pos < 0.0))) ||
            (!has_map && (chrom != NA_INTEGER || !ISNA(pos))))
        {
            gs_info(gv, -2, "Invalid chromosome or position in STvcf");
            failed = 1;
            goto cleanup;
        }
        if (has_map)
            chrom_seen[chrom - 1] = 1;

        if ((NA_STRING == id) || !st_stvcf_valid_name(translateCharUTF8(id)))
        {
            gs_info(gv, -2, "Invalid marker ID in STvcf");
            failed = 1;
            goto cleanup;
        }
        id_len = strlen(translateCharUTF8(id));
        if (id_len >= gs_MAXNME)
        {
            gs_info(gv, -2, "Marker ID too long in STvcf");
            failed = 1;
            goto cleanup;
        }

        if ((NA_STRING == ref) || ('\0' == translateCharUTF8(ref)[0]) ||
            (0 == strcmp(translateCharUTF8(ref), ".")) || (NA_STRING == alt) ||
            (0 != st_stvcf_alleles(translateCharUTF8(ref),
                                   translateCharUTF8(alt), &(n_alt[m]))))
        {
            gs_info(gv, -2, "Invalid REF or ALT field in STvcf");
            failed = 1;
            goto cleanup;
        }

        map_s[m].chrom = chrom;
        map_s[m].pos = pos;
        map_s[m].name = translateCharUTF8(id);
        map_s[m].class = "m";
        map_s[m].nr = m;
        map[m] = &(map_s[m]);

        if (has_map &&
            ((0 != gs_pool_add_string_size(&(tmp.map_str_N), map_s[m].name)) ||
             (0 != gs_pool_add_string_size(&(tmp.map_str_N), map_s[m].class))))
        {
            gs_info(gv, -2, "STvcf map string storage too large");
            failed = 1;
            goto cleanup;
        }
    }

    for (gs_UNSIGNED c = 0; c < no_chrom; c++)
        if (!chrom_seen[c])
        {
            gs_info(gv, -2, "STvcf chromosome codes must be consecutive");
            failed = 1;
            goto cleanup;
        }

    if (has_map)
        qsort(map, no_mar, sizeof(*map), st_stvcf_map_cmp);
    if (has_map && 0 != st_stvcf_separate_map_ties(gv, map, no_mar))
    {
        failed = 1;
        goto cleanup;
    }

    tmp.NoInd = no_ind;
    tmp.NoMar = no_mar;
    tmp.NLoc = has_map ? no_mar : 0;
    tmp.IndName = calloc((size_t)no_ind, sizeof(*tmp.IndName));
    tmp.MarName = gs_malloc_array((size_t)no_mar, sizeof(*tmp.MarName));
    if (has_map)
    {
        tmp.map_S = gs_malloc_array((size_t)no_mar, sizeof(*tmp.map_S));
        tmp.map_str_S = gs_malloc_array(tmp.map_str_N, sizeof(*tmp.map_str_S));
        tmp.map = gs_malloc_array((size_t)no_mar, sizeof(*tmp.map));
        tmp.mdl = gs_malloc_array((size_t)no_mar, sizeof(*tmp.mdl));
    }

    if ((NULL == tmp.IndName) || (NULL == tmp.MarName) ||
        (has_map && ((NULL == tmp.map_S) || (NULL == tmp.map_str_S) ||
                     (NULL == tmp.map) || (NULL == tmp.mdl))) ||
        (0 != gs_alloc_genotype_matrix(&tmp, no_mar, no_ind)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED i = 0; i < no_ind; i++)
    {
        SEXP s = STRING_ELT(ind_s, i);
        size_t len;

        if ((NA_STRING == s) || !st_stvcf_valid_name(translateCharUTF8(s)))
        {
            gs_info(gv, -2, "Invalid individual name in STvcf");
            failed = 1;
            goto cleanup;
        }
        len = strlen(translateCharUTF8(s));
        if (len >= INT_MAX)
        {
            gs_info(gv, -2, "Individual name exceeds R string limit in STvcf");
            failed = 1;
            goto cleanup;
        }
        if (st_identity_set(&tmp, i, translateCharUTF8(s), NULL))
        {
            failed = 1;
            goto cleanup;
        }
        if ((int)len > tmp.MaxIL)
            tmp.MaxIL = (int)len;
    }

    if (0 != gs_make_ind_name_index(tmp.IndName, tmp.NoInd, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, tmp.NoInd))
    {
        gs_info(gv, -2, "Duplicate individual name in STvcf");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    {
        size_t map_str_offset = 0;

        for (gs_UNSIGNED m_new = 0; m_new < no_mar; m_new++)
        {
            gs_UNSIGNED m_old = map[m_new]->nr;
            size_t name_len = strlen(map[m_new]->name);

            strcpy(tmp.MarName[m_new], map[m_new]->name);
            if ((int)name_len > tmp.MaxML)
                tmp.MaxML = (int)name_len;

            if (has_map &&
                0 != gs_copy_map_point_to_pool(&(tmp.map_S[m_new]), map[m_new],
                                               m_new, tmp.map_str_S,
                                               tmp.map_str_N, &map_str_offset))
            {
                gs_info(gv, -2, "Invalid STvcf map string storage");
                failed = 1;
                goto cleanup;
            }
            if (has_map)
            {
                tmp.map[m_new] = &(tmp.map_S[m_new]);
                tmp.mdl[m_new] = m_new;
            }

            GS_ALLELE_TYPE actg[256];
            if (strict)
                actg_keep[m_new] =
                    st_vcf_actg_codes(
                        translateCharUTF8(STRING_ELT(ref_s, m_old)),
                        translateCharUTF8(STRING_ELT(alt_s, m_old)), actg) > 0;

            for (gs_UNSIGNED i = 0; i < no_ind; i++)
            {
                R_xlen_t ridx =
                    (R_xlen_t)m_old + (R_xlen_t)no_mar * (R_xlen_t)i;
                size_t idx = (size_t)m_new * (size_t)no_ind + i;
                int miss1, miss2;
                GS_ALLELE_TYPE a1, a2;

                if ((0 != st_stvcf_get_allele(allele1_s, ridx, n_alt[m_old],
                                              &miss1, &a1)) ||
                    (0 != st_stvcf_get_allele(allele2_s, ridx, n_alt[m_old],
                                              &miss2, &a2)))
                {
                    gs_info(gv, -2, "Invalid allele index in STvcf");
                    failed = 1;
                    goto cleanup;
                }

                if (strict)
                    st_vcf_actg_pair(&a1, &a2, actg);
                if (miss1 || miss2 || a1 == GS_ALLELE_MISSING ||
                    a2 == GS_ALLELE_MISSING)
                {
                    tmp.allele_raw1[idx] = GS_ALLELE_MISSING;
                    tmp.allele_raw2[idx] = GS_ALLELE_MISSING;
                }
                else
                {
                    tmp.allele_raw1[idx] = a1;
                    tmp.allele_raw2[idx] = a2;
                }
            }
        }

        if (map_str_offset != tmp.map_str_N)
        {
            gs_info(gv, -2, "Invalid STvcf map string storage");
            failed = 1;
            goto cleanup;
        }
    }

    if (0 != gs_make_name_index(tmp.MarName, tmp.NoMar, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, tmp.NoMar))
    {
        gs_info(gv, -2, "Duplicate marker ID in STvcf");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    {
        if (strict && st_vcf_actg_compact(&tmp, actg_keep))
        {
            failed = 1;
            goto cleanup;
        }
        int srv = 0;
        char *empty = "";
        char **empty_p = &empty;
        int zero = 0;
        int *zero_p = &zero;

        gs_marker_stats_01(&tmp, empty_p, empty_p, empty_p, zero_p, zero_p,
                           zero_p, zero_p, &srv);
        if (-2 == srv)
        {
            failed = 1;
            goto cleanup;
        }
    }

    if (!strict && !isNull(allele_codes))
    {
        SEXP labels = getAttrib(allele_codes, R_NamesSymbol);
        if (TYPEOF(allele_codes) != VECSXP || (uintmax_t)XLENGTH(allele_codes) != (uintmax_t)no_mar ||
            !R_compute_identical(labels, id_s, 0))
        {
            failed = 1;
            goto cleanup;
        }
        for (gs_UNSIGNED m_new = 0; m_new < no_mar; m_new++)
        {
            gs_UNSIGNED m_old = map[m_new]->nr;
            SEXP codes = VECTOR_ELT(allele_codes, m_old);
            if (TYPEOF(codes) != INTSXP || (uintmax_t)XLENGTH(codes) != (uintmax_t)n_alt[m_old] + 1)
            {
                failed = 1;
                goto cleanup;
            }
            unsigned char seen[256] = {0};
            for (int j = 0; j < LENGTH(codes); j++)
            {
                int code = INTEGER(codes)[j];
                if (code < 0 || code > GS_ALLELE_MAX || seen[code])
                {
                    failed = 1;
                    goto cleanup;
                }
                seen[code] = 1;
            }
            for (gs_UNSIGNED i = 0; i < no_ind; i++)
            {
                size_t k = (size_t)m_new * no_ind + i;
                if (tmp.allele_raw1[k] != GS_ALLELE_MISSING)
                {
                    tmp.allele_raw1[k] =
                        (GS_ALLELE_TYPE)INTEGER(codes)[tmp.allele_raw1[k] - 1];
                    tmp.allele_raw2[k] =
                        (GS_ALLELE_TYPE)INTEGER(codes)[tmp.allele_raw2[k] - 1];
                }
            }
        }

        int srv = 0, zero = 0;
        char *empty = "";
        gs_marker_stats_01(&tmp, &empty, &empty, &empty, &zero, &zero, &zero,
                           &zero, &srv);
        if (srv == -2)
        {
            failed = 1;
            goto cleanup;
        }
    }
    if (st_marker_pedigree_import(&tmp, getAttrib(stvcf, install("pedigree"))))
    {
        gs_info(gv, -2, "Invalid or inconsistent pedigree metadata");
        failed = 1;
        goto cleanup;
    }
    if (!strict && st_marker_prediction_import(
                       &tmp, getAttrib(stvcf, install("prediction.state"))))
    {
        gs_info(gv, -2, "Invalid stored marker prediction state");
        failed = 1;
        goto cleanup;
    }
    if (gs_tmp_commit(gv, &tmp))
    {
        failed = 1;
        goto cleanup;
    }
    status = 0;

cleanup:
    free(actg_keep);
    free(map_s);
    free(map);
    free(n_alt);
    free(chrom_seen);
    free(index);

    if (failed)
        gs_reset(&tmp);

    return ScalarInteger(status);
}

SEXP st_markerdata_to_STvcf(SEXP dataset)
{
    gs_varset_type *gv;
    gs_name_index_type *index = NULL;
    int *map_chrom = NULL;
    double *map_pos = NULL;
    unsigned char *map_seen = NULL;
    int *chrom_levels = NULL;
    GS_ALLELE_TYPE *alleles = NULL;
    gs_UNSIGNED no_chrom = 0;
    int failed = 0;
    int nprot = 0;
    SEXP out = R_NilValue;
    SEXP chrom_s, pos_s, id_s, ref_s, alt_s, ind_s;
    SEXP allele1_s, allele2_s, chrom_levels_s;
    SEXP names_s, class_s, allele_codes;

    if ((STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
        (NA_STRING == STRING_ELT(dataset, 0)))
        return R_NilValue;

    gv = gs_fdta((char *)CHAR(STRING_ELT(dataset, 0)));

    if ((0 == gv->NoInd) || (0 == gv->NoMar) || (NULL == gv->IndName) ||
        (NULL == gv->MarName) || !st_allele_present(gv))
    {
        gs_info(gv, -2, "No marker data available");
        return R_NilValue;
    }

    if ((gv->NoInd > INT_MAX) || (gv->NoMar > INT_MAX) ||
        ((size_t)gv->NoMar > SIZE_MAX / (size_t)gv->NoInd))
    {
        gs_info(gv, -2, "Marker data are too large for STvcf conversion");
        return R_NilValue;
    }

    if ((!gv->NLoc && !gv->physical_map) ||
        (gv->NLoc && (!gv->map_S || !gv->map_str_S || !gv->map_str_N ||
                     !gv->map || !gv->mdl)))
    {
        gs_info(gv, -2, "No linkage map available");
        return R_NilValue;
    }

    if (gv->NLoc &&
        !gs_map_pool_valid(gv->map_S, gv->NLoc, gv->map_str_S, gv->map_str_N))
    {
        gs_info(gv, -2, "Invalid linkage map string storage");
        return R_NilValue;
    }

    map_chrom = gs_malloc_array((size_t)gv->NoMar, sizeof(*map_chrom));
    map_pos = gs_malloc_array((size_t)gv->NoMar, sizeof(*map_pos));
    map_seen = calloc((size_t)gv->NoMar, sizeof(*map_seen));
    chrom_levels = gs_malloc_array((size_t)gv->NoMar, sizeof(*chrom_levels));
    alleles = gs_malloc_array((size_t)gs_MAXALL, sizeof(*alleles));
    if ((NULL == map_chrom) || (NULL == map_pos) || (NULL == map_seen) ||
        (NULL == chrom_levels) || (NULL == alleles))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
    {
        const char *name = gv->IndName[i];

        if (!st_stvcf_valid_name(name) || (strlen(name) >= INT_MAX))
        {
            gs_info(gv, -2,
                    "Individual names are not valid for STvcf conversion");
            failed = 1;
            goto cleanup;
        }
    }

    if (0 != gs_make_ind_name_index(gv->IndName, gv->NoInd, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, gv->NoInd))
    {
        gs_info(gv, -2, "Duplicate individual name in marker data");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
        const char *name = gv->MarName[m];

        if (!st_stvcf_valid_name(name) || (strlen(name) >= gs_MAXNME))
        {
            gs_info(gv, -2, "Marker names are not valid for STvcf conversion");
            failed = 1;
            goto cleanup;
        }
    }

    if (0 != gs_make_name_index(gv->MarName, gv->NoMar, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, gv->NoMar))
    {
        gs_info(gv, -2, "Duplicate marker name in marker data");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    for (gs_UNSIGNED l = 0; l < gv->NLoc; l++)
    {
        gs_UNSIGNED m = gv->mdl[l];
        gs_mappoint_type *point;

        if (gs_UNSIGNED_MAX == m)
            continue;

        if ((m >= gv->NoMar) || (NULL == gv->map[l]))
        {
            gs_info(gv, -2, "Invalid linkage map state");
            failed = 1;
            goto cleanup;
        }

        point = gv->map[l];
        if (map_seen[m] || (point->chrom < 1) || !isfinite(point->pos) ||
            (point->pos < 0.0) || (0 != strcmp(point->name, gv->MarName[m])))
        {
            gs_info(
                gv, -2,
                "Marker data do not have one valid map position per marker");
            failed = 1;
            goto cleanup;
        }

        map_seen[m] = 1;
        map_chrom[m] = point->chrom;
        map_pos[m] = point->pos;
    }

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
        if (gv->NLoc && !map_seen[m])
        {
            gs_info(
                gv, -2,
                "Marker data do not have one valid map position per marker");
            failed = 1;
            goto cleanup;
        }

    for (gs_UNSIGNED m = 0; gv->NLoc && m < gv->NoMar; m++)
    {
        int found = -1;

        for (gs_UNSIGNED c = 0; c < no_chrom; c++)
            if (chrom_levels[c] == map_chrom[m])
            {
                found = (int)c;
                break;
            }

        if (found < 0)
            chrom_levels[no_chrom++] = map_chrom[m];
    }

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
        gs_UNSIGNED n_alleles = 0;
        int has_missing = 0;

        if (0 != st_markerdata_observed_alleles(gv, m, alleles, &n_alleles,
                                                &has_missing))
        {
            gs_info(gv, -2,
                    "Marker alleles are not valid for STvcf conversion");
            failed = 1;
            goto cleanup;
        }
    }

    allele_codes = PROTECT(allocVector(VECSXP, gv->NoMar));
    nprot++;
    out = PROTECT(allocVector(VECSXP, 9));
    nprot++;
    chrom_s = PROTECT(allocVector(INTSXP, gv->NoMar));
    nprot++;
    pos_s = PROTECT(allocVector(REALSXP, gv->NoMar));
    nprot++;
    id_s = PROTECT(allocVector(STRSXP, gv->NoMar));
    nprot++;
    ref_s = PROTECT(allocVector(STRSXP, gv->NoMar));
    nprot++;
    alt_s = PROTECT(allocVector(STRSXP, gv->NoMar));
    nprot++;
    ind_s = PROTECT(allocVector(STRSXP, gv->NoInd));
    nprot++;

    allele1_s = PROTECT(allocMatrix(RAWSXP, gv->NoMar, gv->NoInd));
    nprot++;
    allele2_s = PROTECT(allocMatrix(RAWSXP, gv->NoMar, gv->NoInd));
    nprot++;
    memset(RAW(allele1_s), 255, (size_t)gv->NoMar * gv->NoInd);
    memset(RAW(allele2_s), 255, (size_t)gv->NoMar * gv->NoInd);

    chrom_levels_s = PROTECT(allocVector(STRSXP, no_chrom));
    nprot++;
    names_s = PROTECT(allocVector(STRSXP, 9));
    nprot++;
    class_s = PROTECT(mkString("STvcf"));
    nprot++;

    for (gs_UNSIGNED c = 0; c < no_chrom; c++)
    {
        char label[64];
        snprintf(label, sizeof(label), "%d", chrom_levels[c]);
        SET_STRING_ELT(chrom_levels_s, c, mkChar(label));
    }

    for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        SET_STRING_ELT(ind_s, i, mkChar(gv->IndName[i]));

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
        gs_UNSIGNED n_alleles = 0;
        int has_missing = 0;
        int chrom_code = -1;

        if (0 != st_markerdata_observed_alleles(gv, m, alleles, &n_alleles,
                                                &has_missing))
        {
            gs_info(gv, -2,
                    "Marker alleles are not valid for STvcf conversion");
            failed = 1;
            goto cleanup;
        }

        for (gs_UNSIGNED c = 0; c < no_chrom; c++)
            if (chrom_levels[c] == map_chrom[m])
            {
                chrom_code = (int)c + 1;
                break;
            }

        if (gv->NLoc && chrom_code < 1)
        {
            gs_info(gv, -2, "Invalid chromosome mapping");
            failed = 1;
            goto cleanup;
        }

        INTEGER(chrom_s)[m] = gv->NLoc ? chrom_code : NA_INTEGER;
        REAL(pos_s)[m] = gv->NLoc ? map_pos[m] : NA_REAL;
        SET_STRING_ELT(id_s, m, mkChar(gv->MarName[m]));
        SET_STRING_ELT(ref_s, m, mkChar("1"));
        SEXP codes = PROTECT(allocVector(INTSXP, n_alleles ? n_alleles : 1));
        if (!n_alleles)
            INTEGER(codes)[0] = 1;
        for (gs_UNSIGNED k = 0; k < n_alleles; k++)
            INTEGER(codes)[k] = (int)alleles[k];
        SET_VECTOR_ELT(allele_codes, m, codes);
        UNPROTECT(1);

        if (n_alleles <= 1)
            SET_STRING_ELT(alt_s, m, mkChar("."));
        else
        {
            size_t capacity = (size_t)n_alleles * 12 + 1;
            char *alt = malloc(capacity);
            size_t used = 0;

            if (NULL == alt)
            {
                gs_info(gv, -2, "Memory allocation error");
                failed = 1;
                goto cleanup;
            }

            alt[0] = '\0';
            for (gs_UNSIGNED a = 2; a <= n_alleles; a++)
            {
                int n = snprintf(alt + used, capacity - used,
                                 (2 == a) ? "%lu" : ",%lu", a);

                if ((n < 0) || ((size_t)n >= capacity - used))
                {
                    free(alt);
                    gs_info(gv, -2, "ALT field construction failed");
                    failed = 1;
                    goto cleanup;
                }
                used += (size_t)n;
            }

            SET_STRING_ELT(alt_s, m, mkChar(alt));
            free(alt);
        }

        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        {
            size_t old_idx = (size_t)m * (size_t)gv->NoInd + i;
            R_xlen_t ridx = (R_xlen_t)m + (R_xlen_t)gv->NoMar * (R_xlen_t)i;
            GS_ALLELE_TYPE a1 = st_allele_get(gv, old_idx, 0);
            GS_ALLELE_TYPE a2 = st_allele_get(gv, old_idx, 1);

            if ((GS_ALLELE_MISSING == a1) || (GS_ALLELE_MISSING == a2))
                continue;

            {
                int i1 = st_markerdata_allele_index(alleles, n_alleles, a1);
                int i2 = st_markerdata_allele_index(alleles, n_alleles, a2);

                if ((i1 < 0) || (i2 < 0))
                {
                    gs_info(gv, -2, "Allele normalization failed");
                    failed = 1;
                    goto cleanup;
                }

                RAW(allele1_s)[ridx] = (unsigned char)i1;
                RAW(allele2_s)[ridx] = (unsigned char)i2;
            }
        }
    }

    SET_VECTOR_ELT(out, 0, chrom_s);
    SET_VECTOR_ELT(out, 1, pos_s);
    SET_VECTOR_ELT(out, 2, id_s);
    SET_VECTOR_ELT(out, 3, ref_s);
    SET_VECTOR_ELT(out, 4, alt_s);
    SET_VECTOR_ELT(out, 5, ind_s);
    SET_VECTOR_ELT(out, 6, allele1_s);
    SET_VECTOR_ELT(out, 7, allele2_s);
    SET_VECTOR_ELT(out, 8, chrom_levels_s);

    {
        const char *names[9] = {"CHROM",   "POS",     "ID",
                                "REF",     "ALT",     "individual",
                                "allele1", "allele2", "CHROM.levels"};

        for (int i = 0; i < 9; i++)
            SET_STRING_ELT(names_s, i, mkChar(names[i]));
    }

    setAttrib(allele_codes, R_NamesSymbol, id_s);
    setAttrib(out, install("allele.codes"), allele_codes);
    setAttrib(out, R_NamesSymbol, names_s);
    setAttrib(out, R_ClassSymbol, class_s);
    if (gv->physical_map)
    {
        SEXP physical = PROTECT(st_data_to_R(gv->physical_map));
        setAttrib(out, install("physical.map"), physical);
        UNPROTECT(1);
    }
    SEXP pedigree_s = PROTECT(st_marker_pedigree_object(gv));
    nprot++;
    setAttrib(out, install("pedigree"), pedigree_s);
    if (gv->prediction_state)
    {
        SEXP state = PROTECT(st_marker_prediction_export(gv));
        setAttrib(out, install("prediction.state"), state);
        UNPROTECT(1);
    }
    {
        SEXP phased_s = PROTECT(ScalarLogical(gv->phased ? TRUE : FALSE));
        nprot++;
        setAttrib(out, install("phased"), phased_s);
    }

cleanup:
    free(index);
    free(map_chrom);
    free(map_pos);
    free(map_seen);
    free(chrom_levels);
    free(alleles);

    if (failed)
        out = R_NilValue;

    UNPROTECT(nprot);
    return out;
}

SEXP st_load_performance_data(SEXP id_s, SEXP y_s, SEXP dataset)
{
    gs_varset_type *gv;
    gs_varset_type tmp;
    gs_name_index_type *index = NULL;
    double *values = NULL;
    unsigned char *seen = NULL;
    gs_UNSIGNED no_seen = 0;
    int failed = 0;
    int status = -2;

    if ((STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
        (NA_STRING == STRING_ELT(dataset, 0)))
        return ScalarInteger(-2);

    gv = gs_fdta((char *)translateCharUTF8(STRING_ELT(dataset, 0)));
    gs_tmp_init_from(&tmp, gv);

    if ((0 == gv->NoInd) || (0 == gv->NoMar) || (NULL == gv->IndName) ||
        (NULL == gv->MarName) || !st_allele_present(gv))
    {
        gs_info(gv, -2,
                "No marker data available, run marker-data input first");
        failed = 1;
        goto cleanup;
    }

    if ((STRSXP != TYPEOF(id_s)) || (REALSXP != TYPEOF(y_s)) ||
        (0 == XLENGTH(id_s)) || (XLENGTH(id_s) != XLENGTH(y_s)))
    {
        gs_info(gv, -2, "Invalid performance data");
        failed = 1;
        goto cleanup;
    }

    if (0 != gs_make_ind_name_index(gv->IndName, gv->NoInd, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, gv->NoInd))
    {
        gs_info(gv, -2, "Duplicate individual name in marker data");
        failed = 1;
        goto cleanup;
    }

    values = gs_malloc_array((size_t)gv->NoInd, sizeof(*values));
    seen = calloc((size_t)gv->NoInd, sizeof(*seen));
    if ((NULL == values) || (NULL == seen))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (R_xlen_t j = 0; j < XLENGTH(id_s); j++)
    {
        SEXP id = STRING_ELT(id_s, j);
        gs_UNSIGNED i;
        double value;

        if ((NA_STRING == id) || !st_stvcf_valid_name(translateCharUTF8(id)) ||
            (strlen(translateCharUTF8(id)) >= INT_MAX))
        {
            gs_info(gv, -2, "Invalid individual ID in performance data");
            failed = 1;
            goto cleanup;
        }

        i = gs_name_index_find(index, gv->NoInd, translateCharUTF8(id));
        if (gs_UNSIGNED_MAX == i)
            continue;

        value = REAL(y_s)[j];
        if (!isfinite(value) && !R_IsNA(value))
        {
            gs_info(gv, -2, "Invalid performance value");
            failed = 1;
            goto cleanup;
        }

        values[i] = value;
        if (!seen[i])
        {
            seen[i] = 1;
            no_seen++;
        }
    }

    if (0 == no_seen)
    {
        gs_info(gv, -2, "No performance data match marker-data individuals");
        failed = 1;
        goto cleanup;
    }

    tmp.NoInd = no_seen;
    tmp.NoMar = gv->NoMar;
    tmp.phased = gv->phased;
    tmp.MaxIL = gv->MaxIL;
    tmp.MaxML = gv->MaxML;
    tmp.IndName = calloc((size_t)tmp.NoInd, sizeof(*tmp.IndName));
    tmp.MarName = gs_malloc_array((size_t)tmp.NoMar, sizeof(*tmp.MarName));
    tmp.y = gs_malloc_array((size_t)tmp.NoInd, sizeof(*tmp.y));

    if ((NULL == tmp.IndName) || (NULL == tmp.MarName) || (NULL == tmp.y) ||
        (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED m = 0; m < tmp.NoMar; m++)
        strcpy(tmp.MarName[m], gv->MarName[m]);

    {
        gs_UNSIGNED k = 0;

        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            if (seen[i])
            {
                if (st_identity_copy(&tmp, k, gv, i))
                {
                    failed = 1;
                    goto cleanup;
                }
                tmp.y[k] = values[i];

                for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
                {
                    size_t old_idx = (size_t)m * (size_t)gv->NoInd + i;
                    size_t new_idx = (size_t)m * (size_t)tmp.NoInd + k;
                    tmp.allele_raw1[new_idx] = st_allele_get(gv, old_idx, 0);
                    tmp.allele_raw2[new_idx] = st_allele_get(gv, old_idx, 1);
                }
                k++;
            }
    }

    {
        int map_rc = gs_copy_map_primary(&tmp, gv);
        if (0 != map_rc)
        {
            gs_info(gv, -2,
                    (-1 == map_rc) ? "Memory allocation error"
                                   : "Invalid linkage map state");
            failed = 1;
            goto cleanup;
        }
    }

    {
        int srv = 0;
        char *empty = "";
        char **empty_p = &empty;
        int zero = 0;
        int *zero_p = &zero;

        gs_marker_stats_01(&tmp, empty_p, empty_p, empty_p, zero_p, zero_p,
                           zero_p, zero_p, &srv);
        if (-2 == srv)
        {
            failed = 1;
            goto cleanup;
        }
    }

    snprintf(gv->msg, sizeof(gv->msg), "No. of performance data: %lu",
             tmp.NoInd);
    gs_info(gv, 0, gv->msg);

    if (gs_tmp_commit(gv, &tmp))
    {
        failed = 1;
        goto cleanup;
    }
    status = 0;

cleanup:
    free(index);
    free(values);
    free(seen);

    if (failed)
        gs_reset(&tmp);

    return ScalarInteger(status);
}

void gs_get_map_01(gs_varset_type *gv, char **out_filename, int *retval)
{
    gs_UNSIGNED lne = 0;
    FILE *outfp;

    if ((outfp = fopen(*out_filename, "w")) == NULL)
        ERR_F;

    if (!gv->NLoc && gv->physical_map)
    {
        const st_data_value *names = st_data_get(gv->physical_map, "marker"),
                            *chrom = st_data_get(gv->physical_map, "chromosome"),
                            *pos = st_data_get(gv->physical_map, "position");
        gs_name_index_type *index = NULL;
        if (gs_make_name_index(gv->MarName, gv->NoMar, &index))
        {
            fclose(outfp);
            *retval = -2;
            return;
        }
        for (size_t i = 0; i < names->length; i++)
        {
            const char *name = st_data_string(names, i);
            if (gs_name_index_find(index, gv->NoMar, name) == gs_UNSIGNED_MAX)
                continue;
            char label[64];
            const char *group = chrom->kind == ST_DATA_STRING
                                    ? st_data_string(chrom, i) : label;
            if (chrom->kind != ST_DATA_STRING)
                snprintf(label, sizeof(label), "%.15g",
                         chrom->kind == ST_DATA_REAL ? chrom->data.real[i]
                                                    : chrom->data.integer[i]);
            double value = pos->kind == ST_DATA_REAL ? pos->data.real[i]
                                                    : pos->data.integer[i];
            fprintf(outfp, "%s %.15g %s bp\n", group, value, name);
        }
        free(index);
    }

    for (lne = 0; lne < gv->NLoc; lne++)
        if (gv->mdl[lne] < gs_UNSIGNED_MAX)
        {
            snprintf(gv->msg, sizeof(gv->msg), "%i %15.12f %s %s\n",
                     (gv->map[lne])->chrom, (gv->map[lne])->pos,
                     (gv->map[lne])->name, (gv->map[lne])->class);
            fprintf(outfp, "%s", gv->msg);
        }

    fclose(outfp);
    *retval = 0;
    return;
}

void gs_get_map_01_GV(char **out_filename, int *retval, char **set_name)
{

    gs_get_map_01(gs_fdta(*set_name), out_filename, retval);
}

int minimum(int a, int b)
{
    return (a < b) ? a : b;
}

void gs_str_rem_us(char *strg)
{
    gs_UNSIGNED i, ii = strlen(strg);
    for (i = 0; i < ii; i++)
    {
        if (('.' == strg[i]) || ('-' == strg[i]) || ('/' == strg[i]) ||
            ('_' == strg[i]))
            strg[i] = 'x';
    }
}

void gs_recode_ref_01(gs_varset_type *gv, int *reference, int *descending,
                      int *retval)

{
    gs_UNSIGNED i, j, idx, idx_R;
    int a, b, c, d;

    if ((*reference < 1) || ((gs_UNSIGNED)*reference > gv->NoInd))
    {
        gs_info(gv, -2, "Invalid reference individual");
        *retval = -2;
        return;
    }

    GS_ALLELE_TYPE *Ref1 = malloc(gv->NoMar * gv->NoInd * sizeof(*Ref1));
    GS_ALLELE_TYPE *Ref2 = malloc(gv->NoMar * gv->NoInd * sizeof(*Ref2));
    if ((NULL == Ref1) || (NULL == Ref2))
        ERR_M;

    gs_UNSIGNED R = (gs_UNSIGNED)*reference - 1;

    for (i = 0; i < gv->NoMar; i++)
    {
        idx_R = i * gv->NoInd + R;
        a = st_allele_get(gv, idx_R, 0);
        b = st_allele_get(gv, idx_R, 1);

        for (j = 0; j < gv->NoInd; j++)
        {
            idx = i * gv->NoInd + j;
            c = st_allele_get(gv, idx, 0);
            d = st_allele_get(gv, idx, 1);

            if ((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b) ||
                (GS_ALLELE_MISSING == c) || (GS_ALLELE_MISSING == d))
            {
                Ref1[idx] = GS_ALLELE_MISSING;
                Ref2[idx] = GS_ALLELE_MISSING;
            }
            else
            {
                if (0 == *descending)
                {
                    if (((c == a) && (d == b)) || ((c == b) && (d == a)))
                    {
                        Ref1[idx] = 1;
                        Ref2[idx] = 1;
                    }
                    else if (((c == a) || (c == b)) || ((d == a) || (d == b)))
                    {
                        Ref1[idx] = 1;
                        Ref2[idx] = 2;
                    }
                    else
                    {
                        Ref1[idx] = 2;
                        Ref2[idx] = 2;
                    }
                }
                else
                {
                    if (((c == a) || (c == b)) && ((d == a) || (d == b)))
                    {
                        Ref1[idx] = 1;
                        Ref2[idx] = 1;
                    }
                    else if (((c == a) || (c == b)) || ((d == a) || (d == b)))
                    {
                        Ref1[idx] = 1;
                        Ref2[idx] = 2;
                    }
                    else
                    {
                        Ref1[idx] = 2;
                        Ref2[idx] = 2;
                    }
                }
            }
        }
    }

    st_allele_replace_owned(gv, Ref1, Ref2);

    int Srv = 0;
    int *Sretval = &Srv;
    char *empty = "";
    char **ind_filename = &empty;
    char **mar_filename = &empty;
    char **gen_filename = &empty;
    int auf = 0;
    int *Sauxfiles = &auf;
    gs_marker_stats_01(gv, ind_filename, mar_filename, gen_filename, Sauxfiles,
                       Sauxfiles, Sauxfiles, Sauxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_recode_ref_01_GV(int *reference, int *descending, int *retval,
                         char **set_name)
{
    gs_recode_ref_01(gs_fdta(*set_name), reference, descending, retval);
}

void gs_recode_ref_name_01_GV(char **reference_name, int *descending,
                              int *retval, char **set_name)
{
    gs_varset_type *gv = gs_fdta(*set_name);
    gs_UNSIGNED i;
    int reference = 0;

    if ((NULL == reference_name) || (NULL == *reference_name) ||
        ('\0' == (*reference_name)[0]))
    {
        gs_info(gv, -2, "Invalid reference individual");
        *retval = -2;
        return;
    }

    for (i = 0; i < gv->NoInd; i++)
        if (0 == strcmp(gv->IndName[i], *reference_name))
        {
            reference = (int)i + 1;
            break;
        }

    if (0 == reference)
    {
        snprintf(gv->msg, sizeof(gv->msg), "Reference individual not found: %s",
                 *reference_name);
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    gs_recode_ref_01(gv, &reference, descending, retval);
}

gs_FLOAT gs_avg_ld(gs_chrom_type *chrom, gs_UNSIGNED mm_lft, gs_UNSIGNED mm_rgt)
{
    gs_FLOAT ld = 0;
    gs_FLOAT n;
    gs_UNSIGNED idx, i, j;

    for (i = mm_lft; i < mm_rgt; i++)
        for (j = i + 1; j <= mm_rgt; j++)
        {
            idx = i * chrom->NLocC + j;
            ld += chrom->ld[idx];
        }
    n = ((gs_FLOAT)(mm_rgt - mm_lft) + 1) * ((gs_FLOAT)mm_rgt - mm_lft) / 2;

    return (ld / n);
}

gs_FLOAT gs_check_extend(gs_chrom_type *chrom, gs_UNSIGNED mm_lft,
                         gs_UNSIGNED mm_rgt, gs_UNSIGNED tol_lft,
                         gs_UNSIGNED tol_rgt, char **ld_criterion)
{

    gs_UNSIGNED t, idx;
    gs_FLOAT tld, ld = 0;

    if (0 == strcmp("flanking_fast", *ld_criterion))
    {
        idx = mm_lft * chrom->NLocC + mm_rgt;
        ld = chrom->ld[idx];
    }

    if (0 == strcmp("flanking", *ld_criterion))
    {

        idx = mm_lft * chrom->NLocC + mm_rgt;
        ld = chrom->ld[idx];

        if (tol_lft > 0)
        {
            if (tol_lft > mm_lft)
                tol_lft = mm_lft;
            for (t = 1; t <= tol_lft; t++)
            {
                idx = (mm_lft - t) * chrom->NLocC + mm_rgt;
                if (ld < chrom->ld[idx])
                    ld = chrom->ld[idx];
            }
        }

        if (tol_rgt > 0)
        {
            if (tol_rgt > ((chrom->NLocC) - 1 - mm_rgt))
                tol_rgt = ((chrom->NLocC) - 1 - mm_rgt);
            for (t = 1; t <= tol_rgt; t++)
            {
                idx = mm_lft * chrom->NLocC + (mm_rgt + t);
                if (ld < chrom->ld[idx])
                    ld = chrom->ld[idx];
            }
        }
    }

    if (0 == strcmp("average", *ld_criterion))
    {

        ld = gs_avg_ld(chrom, mm_lft, mm_rgt);

        if (tol_lft > 0)
        {
            if (tol_lft > mm_lft)
                tol_lft = mm_lft;
            for (t = 1; t <= tol_lft; t++)
            {
                tld = gs_avg_ld(chrom, mm_lft - t, mm_rgt);
                if (ld < tld)
                    ld = tld;
            }
        }

        if (tol_rgt > 0)
        {
            if (tol_rgt > ((chrom->NLocC) - 1 - mm_rgt))
                tol_rgt = ((chrom->NLocC) - 1 - mm_rgt);
            for (t = 1; t <= tol_rgt; t++)
            {
                tld = gs_avg_ld(chrom, mm_lft, mm_rgt + t);
                if (ld < tld)
                    ld = tld;
            }
        }
    }

    return (ld);
}

static int gs_build_haplotype_map(gs_varset_type *gv, char **hap_symbol)
{
    gs_mappoint_type *Hap_S = NULL;
    char *Hap_str_S = NULL;
    size_t Hap_str_N = 0;
    size_t offset = 0;

    if ((0 == gv->NHap) || (NULL == hap_symbol) || (NULL == *hap_symbol))
        return -2;

    for (gs_UNSIGNED l = 0; l < gv->NHap; l++)
    {
        char name[gs_MAXNME];
        int n = snprintf(name, sizeof(name), "%1s%06lu", *hap_symbol,
                         (unsigned long)l);
        size_t needed;

        if ((n < 0) || ((size_t)n >= sizeof(name)) || ((size_t)n >= gs_MAXSTR))
            return -2;

        needed = (size_t)n + 1;
        if ((0 != gs_size_add(Hap_str_N, needed, &Hap_str_N)) ||
            (0 != gs_size_add(Hap_str_N, 2, &Hap_str_N)))
            return -1;
    }

    Hap_S = gs_malloc_array((size_t)gv->NHap, sizeof(*Hap_S));
    Hap_str_S = gs_malloc_array(Hap_str_N, sizeof(*Hap_str_S));
    if ((NULL == Hap_S) || (NULL == Hap_str_S))
    {
        free(Hap_S);
        free(Hap_str_S);
        return -1;
    }

    for (gs_UNSIGNED l = 0; l < gv->NHap; l++)
    {
        char name[gs_MAXNME];
        int n = snprintf(name, sizeof(name), "%1s%06lu", *hap_symbol,
                         (unsigned long)l);

        if ((n < 0) || ((size_t)n >= sizeof(name)))
        {
            free(Hap_S);
            free(Hap_str_S);
            return -2;
        }

        Hap_S[l].chrom = gv->map[gv->h_u[l]]->chrom;
        Hap_S[l].pos =
            (gv->map[gv->h_o[l]]->pos + gv->map[gv->h_u[l]]->pos) / 2;
        Hap_S[l].nr = l;

        if ((0 != gs_pool_put_string(Hap_str_S, Hap_str_N, &offset, name,
                                     &(Hap_S[l].name))) ||
            (0 != gs_pool_put_string(Hap_str_S, Hap_str_N, &offset, "b",
                                     &(Hap_S[l].class))))
        {
            free(Hap_S);
            free(Hap_str_S);
            return -1;
        }
    }

    if (offset != Hap_str_N)
    {
        free(Hap_S);
        free(Hap_str_S);
        return -1;
    }

    GS_FREE(gv->Hap_S);
    GS_FREE(gv->Hap_str_S);
    gv->Hap_S = Hap_S;
    gv->Hap_str_S = Hap_str_S;
    gv->Hap_str_N = Hap_str_N;

    return 0;
}

void gs_def_hblocks_01(gs_varset_type *gv, int *hap, int *hap_unit,
                       double *ld_threshold, char **ld_criterion,
                       int *tolerance, char **hap_symbol, char **out_filename,
                       int *auxfiles, int *retval)

{

    if (0 == gv->NLoc)
    {
        snprintf(gv->msg, sizeof(gv->msg), "No linkage map has been loaded");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    gs_UNSIGNED i, l, o, n, fib, lib;

    gv->h_u = realloc(gv->h_u, gv->NLoc * sizeof(gs_UNSIGNED));
    gv->h_o = realloc(gv->h_o, gv->NLoc * sizeof(gs_UNSIGNED));
    gv->h_n = realloc(gv->h_n, gv->NLoc * sizeof(gs_UNSIGNED));
    if ((NULL == gv->h_u) || (NULL == gv->h_o) || (gv->h_n == 0))
        ERR_M;

    gv->NHap = 0;
    o = 0;
    n = 0;

    fib = gs_UNSIGNED_MAX;
    lib = gs_UNSIGNED_MAX;

    if ((0 == strcmp("average", *ld_criterion)) ||
        (0 == strcmp("flanking", *ld_criterion)))
    {

        if (0 == strcmp("average", *ld_criterion))
            snprintf(gv->msg, sizeof(gv->msg),
                     "Average LD in haplotype block > %4.2f", *ld_threshold);
        if (0 == strcmp("flanking", *ld_criterion))
            snprintf(gv->msg, sizeof(gv->msg),
                     "LD of markers flanking haplotype blocks > %4.2f",
                     *ld_threshold);
        gs_info(gv, 0, gv->msg);

        if ((NULL == gv->chrom) || (NULL == gv->chrom[0].ld))
        {
            snprintf(gv->msg, sizeof(gv->msg), "No LD values available");
            gs_info(gv, -2, gv->msg);
            *retval = -2;
            return;
        }

        gs_UNSIGNED tol;
        if (*tolerance < 0)
        {
            snprintf(gv->msg, sizeof(gv->msg), "Incorrect tolerance value");
            gs_info(gv, -2, gv->msg);
            *retval = -2;
            return;
        }

        tol = (gs_UNSIGNED)*tolerance;

        gs_UNSIGNED *blk = malloc(gv->NLoc * sizeof(gs_UNSIGNED));
        if (NULL == blk)
            ERR_M;
        for (i = 0; i < gv->NLoc; i++)
            blk[i] = 0;

        gs_UNSIGNED c, mm, idx, lfdnr;
        gs_UNSIGNED locus1, locus2;
        gs_FLOAT ld;
        int unasigned_loci, block_extendable;
        gs_UNSIGNED greatest_ld_locus;
        gs_FLOAT greatest_ld_value, extend_lft, extend_rgt;
        gs_UNSIGNED mm_rgt, mm_lft, loc_rgt, loc_lft;

        lfdnr = 0;

        for (c = 0; c < gv->NChrom; c++)
        {

            do
            {

                greatest_ld_locus = 0;
                greatest_ld_value = 0;

                for (mm = 1; mm < gv->chrom[c].NLocC; mm++)
                {

                    locus2 = gv->chrom[c].first_loc_i + mm;
                    locus1 = locus2 - 1;

                    idx = (mm - 1) * gv->chrom[c].NLocC + mm;
                    ld = (gv->chrom[c].ld)[idx];

                    if ((0 == blk[locus1]) && (0 == blk[locus2]) &&
                        ld > *ld_threshold)
                    {
                        if (ld > greatest_ld_value)
                        {
                            greatest_ld_locus = mm;
                            greatest_ld_value = ld;
                        }
                    }
                }

                if (0 == greatest_ld_locus)

                {
                    for (mm = 0; mm < gv->chrom[c].NLocC; mm++)
                        if (0 == blk[gv->chrom[c].first_loc_i + mm])
                            blk[gv->chrom[c].first_loc_i + mm] = ++lfdnr;
                }
                else
                {

                    mm_rgt = greatest_ld_locus;
                    mm_lft = greatest_ld_locus - 1;

                    loc_rgt = gv->chrom[c].first_loc_i + mm_rgt;
                    loc_lft = gv->chrom[c].first_loc_i + mm_lft;

                    blk[loc_rgt] = ++lfdnr;
                    blk[loc_lft] = lfdnr;

                    do
                    {
                        block_extendable = 0;

                        if (mm_lft > 0)
                            extend_lft =
                                gs_check_extend(&(gv->chrom[c]), mm_lft - 1,
                                                mm_rgt, tol, 0, ld_criterion);
                        else
                            extend_lft = 0;

                        if (mm_rgt < ((gv->chrom[c].NLocC) - 1))
                            extend_rgt = gs_check_extend(&(gv->chrom[c]),
                                                         mm_lft, mm_rgt + 1, 0,
                                                         tol, ld_criterion);
                        else
                            extend_rgt = 0;

                        if ((extend_lft >= extend_rgt) &&
                            (extend_lft >= *ld_threshold))

                        {
                            mm_lft = mm_lft - 1;
                            loc_lft = gv->chrom[c].first_loc_i + mm_lft;
                            blk[loc_lft] = lfdnr;
                            block_extendable = 1;
                        }
                        else if ((extend_rgt > extend_lft) &&
                                 (extend_rgt >= *ld_threshold))

                        {
                            mm_rgt = mm_rgt + 1;
                            loc_rgt = gv->chrom[c].first_loc_i + mm_rgt;
                            blk[loc_rgt] = lfdnr;
                            block_extendable = 1;
                        }

                    } while (1 == block_extendable);
                }

                unasigned_loci = 0;
                for (mm = 0; mm < gv->chrom[c].NLocC; mm++)
                    if (0 == blk[gv->chrom[c].first_loc_i + mm])
                        unasigned_loci = 1;

            } while (1 == unasigned_loci);
        }

        gv->NHap = 0;
        fib = 0;

        for (i = 1; i < gv->NLoc; i++)
        {
            if (blk[i - 1] != blk[i])
            {
                gv->h_u[gv->NHap] = fib;
                gv->h_o[gv->NHap] = i - 1;
                gv->h_n[gv->NHap] = i - fib;
                gv->NHap++;
                fib = i;
            }
        }
        gv->h_u[gv->NHap] = fib;
        gv->h_o[gv->NHap] = i - 1;
        gv->h_n[gv->NHap] = i - fib;
        gv->NHap++;

        free(blk);
    }

    else if (0 == *hap_unit)
    {

        snprintf(gv->msg, sizeof(gv->msg),
                 "Haplotype borders form disjunct LD segments ");
        gs_info(gv, 0, gv->msg);

        gs_chrom_stats_01(gv, Sout_filename, Sauxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_UNSIGNED i, m, cc, l1, l2, idx1, idx2;
        for (cc = 0; cc < gv->NChrom; cc++)
        {
            fib = gv->chrom[cc].first_loc_i;

            for (m = 1; m < gv->chrom[cc].NLocC; m++)
            {
                n++;
                for (i = 0; i < gv->NoInd; i++)
                {
                    l2 = gv->chrom[cc].first_loc_i + m;
                    l1 = l2 - 1;

                    idx1 = l1 * gv->NoInd + i;
                    idx2 = l2 * gv->NoInd + i;

                    if ((st_allele_get(gv, idx1, 0) !=
                         st_allele_get(gv, idx2, 0)) ||
                        (st_allele_get(gv, idx1, 1) !=
                         st_allele_get(gv, idx2, 1)))
                    {
                        lib = fib + n - 1;
                        gv->h_u[gv->NHap] = fib;
                        gv->h_o[gv->NHap] = lib;
                        gv->h_n[gv->NHap] = n;
                        gv->NHap++;
                        n = 0;
                        fib = lib + 1;
                        break;
                    }
                }
            }

            lib = fib + n;
            gv->h_u[gv->NHap] = fib;
            gv->h_o[gv->NHap] = lib;
            gv->h_n[gv->NHap] = n;
            gv->NHap++;
            n = 0;
        }
    }

    else if (1 == *hap_unit)
    {
        snprintf(gv->msg, sizeof(gv->msg), "Haplotyp borders every %i markers",
                 *hap);
        gs_info(gv, 0, gv->msg);

        do
        {

            do
            {
                if (gv->mdl[o] < gs_UNSIGNED_MAX)
                {
                    if (gs_UNSIGNED_MAX == fib)
                        fib = o;
                    lib = o;
                    n++;
                }
                o++;
            } while ((o < gv->NLoc) && !((int)n > (*hap - 1)) &&
                     !(gv->map[o - 1]->chrom != gv->map[o]->chrom));

            if (n > 0)
            {
                gv->h_u[gv->NHap] = fib;
                gv->h_o[gv->NHap] = lib;
                gv->h_n[gv->NHap] = n;
                gv->NHap++;

                n = 0;
                fib = gs_UNSIGNED_MAX;
                lib = gs_UNSIGNED_MAX;
            }

        } while (o < (gv->NLoc));
    }

    else if (2 == *hap_unit)
    {

        snprintf(gv->msg, sizeof(gv->msg),
                 "Haplotyp borders every %i centiMorgans", *hap);

        gs_info(gv, 0, gv->msg);

        do
        {

            do
            {
                if (gv->mdl[o] < gs_UNSIGNED_MAX)
                {
                    if (gs_UNSIGNED_MAX == fib)
                        fib = o;
                    lib = o;
                    n++;
                }
                o++;
            } while (
                (o < gv->NLoc) &&
                !((double)*hap < (gv->map[lib]->pos - gv->map[fib]->pos)) &&
                !(gv->map[o - 1]->chrom != gv->map[o]->chrom));

            if (n > 0)
            {
                gv->h_u[gv->NHap] = fib;
                gv->h_o[gv->NHap] = lib;
                gv->h_n[gv->NHap] = n;
                gv->NHap++;

                n = 0;
                fib = gs_UNSIGNED_MAX;
                lib = gs_UNSIGNED_MAX;
            }
        } while (o < (gv->NLoc));
    }

    else

    {

        snprintf(gv->msg, sizeof(gv->msg),
                 "Provide either hap.unit or ld.criterion");

        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    {
        int map_rc = gs_build_haplotype_map(gv, hap_symbol);

        if (0 != map_rc)
        {
            gv->NHap = 0;
            GS_FREE(gv->Hap_S);
            GS_FREE(gv->Hap_str_S);
            gv->Hap_str_N = 0;

            if (-1 == map_rc)
                ERR_M;

            gs_info(gv, -2, "Haplotype marker name too long");
            *retval = -2;
            return;
        }
    }

    if (1 == *auxfiles)
    {

        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        for (l = 0; l < gv->NHap; l++)
        {
            char *lne = (char *)malloc(gs_MAXLNE * sizeof(char));

            snprintf(lne, gs_MAXLNE, "%i %f %s %s ", gv->Hap_S[l].chrom,
                     gv->Hap_S[l].pos, gv->Hap_S[l].name, gv->Hap_S[l].class);

            for (i = gv->h_u[l]; i <= gv->h_o[l]; i++)
            {
                if (gv->mdl[i] < gs_UNSIGNED_MAX)
                {
                    size_t a = strlen(lne);
                    size_t b = strlen(gv->MarName[gv->mdl[i]]);
                    lne = realloc(lne, 2 + a + b);
                    strcat(lne, gv->MarName[gv->mdl[i]]);
                    strcat(lne, ";");
                }
            }

            fprintf(outfp, "%s\n", lne);
            free(lne);
        }

        fclose(outfp);
    }
}

void gs_def_hblocks_01_GV(int *hap, int *hap_unit, double *ld_threshold,
                          char **ld_criterion, int *tolerance,
                          char **hap_symbol, char **out_filename, int *auxfiles,
                          int *retval, char **set_name)
{
    gs_def_hblocks_01(gs_fdta(*set_name), hap, hap_unit, ld_threshold,
                      ld_criterion, tolerance, hap_symbol, out_filename,
                      auxfiles, retval);
}

int gs_swap_block_marker(gs_varset_type *gv)
{

    gs_UNSIGNED i;

    gs_UNSIGNED NoInd = gv->NoInd;
    gs_UNSIGNED NoMar = gv->NHap;
    int phased = gv->phased;

    char(*MarName)[gs_MAXNME] =
        gs_malloc_array((size_t)NoMar, sizeof(*MarName));
    gs_mappoint_type **map = gs_malloc_array((size_t)NoMar, sizeof(*map));
    gs_UNSIGNED *mdl = gs_malloc_array((size_t)NoMar, sizeof(*mdl));

    if ((NULL == MarName) || (NULL == map) || (NULL == mdl))
    {
        free(MarName);
        free(map);
        free(mdl);
        return -1;
    }

    for (i = 0; i < gv->NHap; i++)
    {
        strcpy(MarName[i], gv->Hap_S[i].name);
        map[i] = &(gv->Hap_S[i]);
        mdl[i] = i;
    }

    gs_varset_type storage;
    gs_tmp_init_from(&storage, gv);
    storage.NoInd = NoInd;
    storage.NoMar = NoMar;
    if (gs_alloc_genotype_matrix(&storage, NoMar, NoInd))
    {
        gs_reset(&storage);
        free(MarName);
        free(map);
        free(mdl);
        return -1;
    }
    memcpy(storage.allele_raw1, gv->Ref1, (size_t)NoMar * NoInd);
    memcpy(storage.allele_raw2, gv->Ref2, (size_t)NoMar * NoInd);
    if (gs_prepare_allele_commit(&storage))
    {
        gs_reset(&storage);
        free(MarName);
        free(map);
        free(mdl);
        return -1;
    }

    char **IndName = gv->IndName;
    char **Pedigree = gv->Pedigree;
    gv->Pedigree = NULL;
    gv->IndName = NULL;

    gs_UNSIGNED NLoc = gv->NHap;
    gs_mappoint_type *map_S = gv->Hap_S;
    gv->Hap_S = NULL;
    char *map_str_S = gv->Hap_str_S;
    gv->Hap_str_S = NULL;
    size_t map_str_N = gv->Hap_str_N;
    gv->Hap_str_N = 0;

    double *y = gv->y;
    gv->y = NULL;

    gs_reset(gv);

    gv->NoInd = NoInd;
    gv->NoMar = NoMar;
    gv->phased = phased;

    gv->IndName = IndName;
    gv->Pedigree = Pedigree;
    gv->MarName = MarName;

    gv->allele_raw1 = storage.allele_raw1;
    gv->allele_raw2 = storage.allele_raw2;
    gv->allele_packed = storage.allele_packed;
    gv->allele_dictionary = storage.allele_dictionary;
    gv->allele_packed_bytes = storage.allele_packed_bytes;
    gv->allele_bits = storage.allele_bits;

    gv->NLoc = NLoc;
    gv->map_S = map_S;
    gv->map_str_S = map_str_S;
    gv->map_str_N = map_str_N;
    gv->map = map;
    gv->mdl = mdl;

    gv->y = y;

    return 0;
}

void gs_recode_hbc_01(gs_varset_type *gv, int *reference, int *retval)

{

    if (0 == gv->NHap)
    {
        snprintf(gv->msg, sizeof(gv->msg), "No haplotype blocks were defined");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    gs_UNSIGNED l, i, j, idx, idx_R;

    if ((*reference < 1) || ((gs_UNSIGNED)*reference > gv->NoInd))
    {
        snprintf(gv->msg, sizeof(gv->msg), "Invalid reference individual");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    gv->Ref1 = realloc(gv->Ref1, gv->NoMar * gv->NoInd * sizeof(*gv->Ref1));
    gv->Ref2 = realloc(gv->Ref2, gv->NoMar * gv->NoInd * sizeof(*gv->Ref2));
    if ((NULL == gv->Ref1) || (NULL == gv->Ref2))
        ERR_M;

    gs_UNSIGNED R = (gs_UNSIGNED)*reference - 1;
    gs_UNSIGNED idx_H, n0, n1, n2;
    int a, b, c, d;

    (void)n2;

    for (l = 0; l < gv->NHap; l++)
    {

        for (j = 0; j < gv->NoInd; j++)
        {
            n0 = 0;
            n1 = 0;
            n2 = 0;
            idx_H = l * gv->NoInd + j;

            for (i = gv->h_u[l]; i <= gv->h_o[l]; i++)
                if (gv->mdl[i] < gs_UNSIGNED_MAX)
                {

                    idx_R = i * gv->NoInd + R;
                    a = st_allele_get(gv, idx_R, 0);
                    b = st_allele_get(gv, idx_R, 1);

                    idx = i * gv->NoInd + j;
                    c = st_allele_get(gv, idx, 0);
                    d = st_allele_get(gv, idx, 1);

                    if (!((GS_ALLELE_MISSING == a) ||
                          (GS_ALLELE_MISSING == b) ||
                          (GS_ALLELE_MISSING == c) || (GS_ALLELE_MISSING == d)))
                    {
                        if (((c == a) || (c == b)) && ((d == a) || (d == b)))
                        {
                            n2++;
                        }
                        else if (((c == a) || (c == b)) ||
                                 ((d == a) || (d == b)))
                        {
                            n1++;
                        }
                        else
                        {
                            n0++;
                        }
                    }
                }

            gv->Ref1[idx_H] = 1;
            gv->Ref2[idx_H] = 1;

            if (n0 > 0)
            {
                gv->Ref1[idx_H] = 2;
                gv->Ref2[idx_H] = 2;
            }
            else if (n1 > 0)
            {
                gv->Ref1[idx_H] = 1;
                gv->Ref2[idx_H] = 2;
            }
        }
    }

    if (0 != gs_swap_block_marker(gv))
        ERR_M;

    int Srv = 0;
    int *Sretval = &Srv;
    char *empty = "";
    char **ind_filename = &empty;
    char **mar_filename = &empty;
    char **gen_filename = &empty;
    int auf = 0;
    int *Sauxfiles = &auf;
    gs_marker_stats_01(gv, ind_filename, mar_filename, gen_filename, Sauxfiles,
                       Sauxfiles, Sauxfiles, Sauxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_recode_hbc_01_GV(int *reference, int *retval, char **set_name)
{
    gs_recode_hbc_01(gs_fdta(*set_name), reference, retval);
}

int gs_hall_cmp(const void *a, const void *b)
{

    const gs_usp *ia = (const gs_usp *)a;
    const gs_usp *ib = (const gs_usp *)b;

    if ((ia)->n < (ib)->n)
        return (1);
    if ((ia)->n > (ib)->n)
        return (-1);
    return (0);
}

void gs_recode_hil_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                      int *retval)

{
    GS_ALLELE_TYPE *AT1 = NULL;
    GS_ALLELE_TYPE *AT2 = NULL;
    GS_ALLELE_TYPE **Hall = NULL;
    gs_UNSIGNED *OcHall = NULL;
    gs_UNSIGNED *Hidx1 = NULL;
    gs_UNSIGNED *Hidx2 = NULL;
    gs_usp *ar = NULL;
    char **hal = NULL;
    gs_UNSIGNED n_hal = 0;
    size_t genotype_count;
    size_t max_hall;
    int failed = 0;

    if (0 == gv->NHap)
    {
        snprintf(gv->msg, sizeof(gv->msg), "No haplotype blocks were defined");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    if ((0 !=
         gs_size_mul((size_t)gv->NoMar, (size_t)gv->NoInd, &genotype_count)) ||
        ((size_t)gv->NoInd > SIZE_MAX / 2))
    {
        gs_info(gv, -2, "Memory allocation error");
        *retval = -2;
        return;
    }
    max_hall = 2 * (size_t)gv->NoInd;

    gv->Ref1 = realloc(gv->Ref1, genotype_count * sizeof(*gv->Ref1));
    gv->Ref2 = realloc(gv->Ref2, genotype_count * sizeof(*gv->Ref2));
    AT1 = gs_malloc_array(genotype_count, sizeof(*AT1));
    AT2 = gs_malloc_array(genotype_count, sizeof(*AT2));
    Hall = gs_malloc_array(max_hall, sizeof(*Hall));
    OcHall = gs_malloc_array(max_hall, sizeof(*OcHall));
    Hidx1 = gs_malloc_array((size_t)gv->NoInd, sizeof(*Hidx1));
    Hidx2 = gs_malloc_array((size_t)gv->NoInd, sizeof(*Hidx2));

    if ((NULL == gv->Ref1) || (NULL == gv->Ref2) || (NULL == AT1) ||
        (NULL == AT2) || (NULL == Hall) || (NULL == OcHall) ||
        (NULL == Hidx1) || (NULL == Hidx2))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED i = 0; i < gv->NoMar; i++)
        for (gs_UNSIGNED j = 0; j < gv->NoInd; j++)
        {
            size_t idx = (size_t)i * (size_t)gv->NoInd + j;
            size_t idx_T = (size_t)j * (size_t)gv->NoMar + i;
            AT1[idx_T] = st_allele_get(gv, idx, 0);
            AT2[idx_T] = st_allele_get(gv, idx, 1);
        }

    for (gs_UNSIGNED l = 0; l < gv->NHap; l++)
    {
        gs_UNSIGNED i = gv->h_u[l];
        gs_UNSIGNED NoHall = 0;
        size_t n = sizeof(GS_ALLELE_TYPE) * (size_t)gv->h_n[l];

        memset(OcHall, 0, max_hall * sizeof(*OcHall));

        Hall[NoHall] = &(AT1[(size_t)i]);
        OcHall[NoHall]++;
        NoHall++;

        for (gs_UNSIGNED j = 0; j < gv->NoInd; j++)
        {
            size_t idx_T = (size_t)j * (size_t)gv->NoMar + i;
            GS_ALLELE_TYPE *P1 = &(AT1[idx_T]);
            int found = 0;

            for (gs_UNSIGNED k = 0; k < NoHall; k++)
                if (0 == memcmp(Hall[k], P1, n))
                {
                    Hidx1[j] = k;
                    OcHall[k]++;
                    found = 1;
                    break;
                }

            if (!found)
            {
                if (NoHall >= GS_ALLELE_MAX)
                {
                    gs_info(gv, -2,
                            "Haplotyping would result in too many alleles");
                    failed = 1;
                    goto cleanup;
                }
                Hall[NoHall] = P1;
                Hidx1[j] = NoHall;
                OcHall[NoHall]++;
                NoHall++;
            }

            P1 = &(AT2[idx_T]);
            found = 0;
            for (gs_UNSIGNED k = 0; k < NoHall; k++)
                if (0 == memcmp(Hall[k], P1, n))
                {
                    Hidx2[j] = k;
                    OcHall[k]++;
                    found = 1;
                    break;
                }

            if (!found)
            {
                if (NoHall >= GS_ALLELE_MAX)
                {
                    gs_info(gv, -2,
                            "Haplotyping would result in too many alleles");
                    failed = 1;
                    goto cleanup;
                }
                Hall[NoHall] = P1;
                Hidx2[j] = NoHall;
                OcHall[NoHall]++;
                NoHall++;
            }
        }

        ar = gs_realloc_array(ar, (size_t)NoHall, sizeof(*ar));
        if (NULL == ar)
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }

        for (gs_UNSIGNED k = 0; k < NoHall; k++)
        {
            ar[k].a = k;
            ar[k].n = OcHall[k];
        }
        qsort(ar, NoHall, sizeof(*ar), gs_hall_cmp);

        for (gs_UNSIGNED j = 0; j < gv->NoInd; j++)
        {
            size_t idx_H = (size_t)l * (size_t)gv->NoInd + j;

            for (gs_UNSIGNED k = 0; k < NoHall; k++)
                if (Hidx1[j] == ar[k].a)
                {
                    gv->Ref1[idx_H] = (GS_ALLELE_TYPE)(1 + k);
                    break;
                }
            for (gs_UNSIGNED k = 0; k < NoHall; k++)
                if (Hidx2[j] == ar[k].a)
                {
                    gv->Ref2[idx_H] = (GS_ALLELE_TYPE)(1 + k);
                    break;
                }
        }

        for (gs_UNSIGNED j = 0; j < gv->NoInd; j++)
        {
            size_t idx_H = (size_t)l * (size_t)gv->NoInd + j;
            gs_UNSIGNED rank = (gs_UNSIGNED)gv->Ref1[idx_H] - 1;

            for (gs_UNSIGNED kk = 0; kk < gv->h_n[l]; kk++)
                if (GS_ALLELE_MISSING == Hall[ar[rank].a][kk])
                {
                    gv->Ref1[idx_H] = GS_ALLELE_MISSING;
                    break;
                }

            rank = (gs_UNSIGNED)gv->Ref2[idx_H] - 1;
            for (gs_UNSIGNED kk = 0; kk < gv->h_n[l]; kk++)
                if (GS_ALLELE_MISSING == Hall[ar[rank].a][kk])
                {
                    gv->Ref2[idx_H] = GS_ALLELE_MISSING;
                    break;
                }
        }

        {
            char **new_hal = gs_realloc_array(hal, (size_t)(n_hal + NoHall),
                                              sizeof(*new_hal));
            if (NULL == new_hal)
            {
                gs_info(gv, -2, "Memory allocation error");
                failed = 1;
                goto cleanup;
            }
            hal = new_hal;
        }

        for (gs_UNSIGNED k = 0; k < NoHall; k++)
        {
            size_t cap = 6 * (size_t)gv->h_n[l] + gs_MAXNME + 10;
            char *tmp_str = malloc(cap);
            char tmp_str2[gs_MAXSTR];

            if (NULL == tmp_str)
            {
                gs_info(gv, -2, "Memory allocation error");
                failed = 1;
                goto cleanup;
            }

            snprintf(tmp_str, cap, "%s %lu ", gv->Hap_S[l].name, k + 1);
            for (gs_UNSIGNED kk = 0; kk < gv->h_n[l]; kk++)
            {
                int value = (GS_ALLELE_MISSING == Hall[ar[k].a][kk])
                                ? -1
                                : (int)Hall[ar[k].a][kk];
                snprintf(tmp_str2, sizeof(tmp_str2),
                         (kk + 1 < gv->h_n[l]) ? "%i;" : "%i", value);
                strcat(tmp_str, tmp_str2);
            }
            hal[n_hal + k] = tmp_str;
        }

        n_hal += NoHall;
    }

    free(AT1);
    AT1 = NULL;
    free(AT2);
    AT2 = NULL;
    free(Hall);
    Hall = NULL;
    free(OcHall);
    OcHall = NULL;
    free(Hidx1);
    Hidx1 = NULL;
    free(Hidx2);
    Hidx2 = NULL;
    free(ar);
    ar = NULL;

    if (0 != gs_swap_block_marker(gv))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    {
        int srv = 0;
        char *empty = "";
        char **empty_p = &empty;
        int zero = 0;

        gs_marker_stats_01(gv, empty_p, empty_p, empty_p, &zero, &zero, &zero,
                           &zero, &srv);
        if (-2 == srv)
        {
            failed = 1;
            goto cleanup;
        }
    }

    if (1 == *auxfiles)
    {
        FILE *outfp = fopen(*out_filename, "w");
        if (NULL == outfp)
        {
            gs_info(gv, -2, "File open error");
            failed = 1;
            goto cleanup;
        }
        for (gs_UNSIGNED k = 0; k < n_hal; k++)
            fprintf(outfp, "%s\n", hal[k]);
        fclose(outfp);
    }

cleanup:
    free(AT1);
    free(AT2);
    free(Hall);
    free(OcHall);
    free(Hidx1);
    free(Hidx2);
    free(ar);
    for (gs_UNSIGNED k = 0; k < n_hal; k++)
        free(hal[k]);
    free(hal);

    if (failed)
    {
        GS_FREE(gv->Ref1);
        GS_FREE(gv->Ref2);
        *retval = -2;
        return;
    }

    *retval = 0;
}

void gs_recode_hil_01_GV(char **out_filename, int *auxfiles, int *retval,
                         char **set_name)
{
    gs_recode_hil_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_write_marker_data_02(gs_varset_type *gv, char **type,
                             char **lpo_filename, char **mmp_filename,
                             char **mpo_filename, char **ngp_filename,
                             char **nmp_filename, char **npo_filename,
                             char **ind_filename, int *f_ind, int *l_ind,
                             char **add_inum, int *auxfiles, int *retval)
{
    gs_UNSIGNED *pli_idx = NULL;
    gs_UNSIGNED npli = 0;
    gs_UNSIGNED capacity;
    char **n_names = NULL;
    FILE *outfp = NULL;
    char line[gs_MAXSTR];
    char t;
    int failed = 0;

    (void)auxfiles;

    if (!st_allele_present(gv) || (0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2, "No marker data available");
        *retval = -2;
        return;
    }

    t = (*type)[0];
    if (!(('l' == t) || ('a' == t) || ('m' == t) || ('n' == t)) ||
        ('\0' != (*type)[1]))
    {
        gs_info(gv, -2, "Implemented output formats: l, a, m, n");
        *retval = -2;
        return;
    }

    capacity = (0 == gv->NoInd) ? 1 : gv->NoInd;
    pli_idx = gs_malloc_array((size_t)capacity, sizeof(*pli_idx));
    if (NULL == pli_idx)
    {
        gs_info(gv, -2, "Memory allocation error");
        *retval = -2;
        return;
    }

    npli = gv->NoInd;
    for (gs_UNSIGNED j = 0; j < npli; j++)
        pli_idx[j] = j;

    if ((0 != *l_ind) || (0 != *f_ind))
    {
        long f = *f_ind;
        long l = *l_ind;

        if ((long)gv->NoInd < f)
            f = (long)gv->NoInd;
        if ((long)gv->NoInd < l)
            l = (long)gv->NoInd;
        if (0 < f)
            f--;
        else if (0 > f)
            f = 0;
        if ((1 > l) || ((long)gv->NoInd < l))
            l = (long)gv->NoInd;

        npli = 0;
        for (long j = f; j < l; j++)
            pli_idx[npli++] = (gs_UNSIGNED)j;
    }

    if (0 != strcmp(*ind_filename, ""))
    {
        FILE *fp = fopen(*ind_filename, "r");
        if (NULL == fp)
        {
            gs_info(gv, -2, "File open error");
            failed = 1;
            goto cleanup;
        }

        npli = 0;
        while (1 == fscanf(fp, "%4095s", line))
        {
            for (gs_UNSIGNED k = 0; k < gv->NoInd; k++)
                if (0 == strcmp(gv->IndName[k], line))
                {
                    if (npli == capacity)
                    {
                        gs_UNSIGNED new_capacity = capacity * 2;
                        gs_UNSIGNED *new_idx;

                        if (new_capacity <= capacity)
                        {
                            fclose(fp);
                            gs_info(gv, -2, "Too many output individuals");
                            failed = 1;
                            goto cleanup;
                        }

                        new_idx = gs_realloc_array(
                            pli_idx, (size_t)new_capacity, sizeof(*new_idx));
                        if (NULL == new_idx)
                        {
                            fclose(fp);
                            gs_info(gv, -2, "Memory allocation error");
                            failed = 1;
                            goto cleanup;
                        }
                        pli_idx = new_idx;
                        capacity = new_capacity;
                    }

                    pli_idx[npli++] = k;
                    break;
                }
        }

        fclose(fp);
    }

    if ((('m' == t) || ('n' == t)) && (NULL != gv->map))
    {
        gs_UNSIGNED *storage = gs_map_storage_index(gv, 0, 0);

        if (NULL == storage)
        {
            failed = 1;
            goto cleanup;
        }
        free(storage);
    }

    if (('l' == t) || ('a' == t))
    {
        if ('a' == t)
            outfp = fopen(*lpo_filename, "a");
        else
            outfp = fopen(*lpo_filename, "w");

        if (NULL == outfp)
        {
            gs_info(gv, -2, "File open error");
            failed = 1;
            goto cleanup;
        }

        if ('l' == t)
            fprintf(outfp, "I;M;A\n");

        for (gs_UNSIGNED ii = 0; ii < npli; ii++)
        {
            gs_UNSIGNED i = pli_idx[ii];

            for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            {
                size_t idx = (size_t)m * (size_t)gv->NoInd + i;

                if (GS_ALLELE_MISSING != st_allele_get(gv, idx, 0))
                {
                    fprintf(outfp, "%s;%s;%u\n", gv->IndName[i], gv->MarName[m],
                            (unsigned)st_allele_get(gv, idx, 0));
                    if (st_allele_get(gv, idx, 0) != st_allele_get(gv, idx, 1))
                        fprintf(outfp, "%s;%s;%u\n", gv->IndName[i],
                                gv->MarName[m],
                                (unsigned)st_allele_get(gv, idx, 1));
                }
            }
        }

        fclose(outfp);
        outfp = NULL;
    }

    if ('m' == t)
    {
        if (NULL != gv->map)
        {
            outfp = fopen(*mmp_filename, "w");
            if (NULL == outfp)
            {
                gs_info(gv, -2, "File open error");
                failed = 1;
                goto cleanup;
            }

            for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
                if (gv->mdl[lne] < gs_UNSIGNED_MAX)
                    fprintf(outfp, "%i %15.12f %s %s\n", gv->map[lne]->chrom,
                            gv->map[lne]->pos, gv->map[lne]->name,
                            gv->map[lne]->class);

            fclose(outfp);
            outfp = NULL;
        }

        outfp = fopen(*mpo_filename, "w");
        if (NULL == outfp)
        {
            gs_info(gv, -2, "File open error");
            failed = 1;
            goto cleanup;
        }

        for (gs_UNSIGNED jj = 0; jj < npli; jj++)
            fprintf(outfp, "%s ", gv->IndName[pli_idx[jj]]);
        fprintf(outfp, "\n");

        if (NULL != gv->map)
        {
            for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
                if (gv->mdl[lne] < gs_UNSIGNED_MAX)
                {
                    gs_UNSIGNED m = gv->mdl[lne];
                    fprintf(outfp, "%s ", gv->MarName[m]);
                    for (gs_UNSIGNED jj = 0; jj < npli; jj++)
                    {
                        gs_UNSIGNED i = pli_idx[jj];
                        size_t idx = (size_t)m * (size_t)gv->NoInd + i;
                        if (GS_ALLELE_MISSING == st_allele_get(gv, idx, 0))
                            fprintf(outfp, "-- ");
                        else
                            fprintf(outfp, "%u/%u ",
                                    (unsigned)st_allele_get(gv, idx, 0),
                                    (unsigned)st_allele_get(gv, idx, 1));
                    }
                    fprintf(outfp, "\n");
                }
        }
        else
        {
            for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            {
                fprintf(outfp, "%s ", gv->MarName[m]);
                for (gs_UNSIGNED jj = 0; jj < npli; jj++)
                {
                    gs_UNSIGNED i = pli_idx[jj];
                    size_t idx = (size_t)m * (size_t)gv->NoInd + i;
                    if (GS_ALLELE_MISSING == st_allele_get(gv, idx, 0))
                        fprintf(outfp, "-- ");
                    else
                        fprintf(outfp, "%u/%u ",
                                (unsigned)st_allele_get(gv, idx, 0),
                                (unsigned)st_allele_get(gv, idx, 1));
                }
                fprintf(outfp, "\n");
            }
        }

        fclose(outfp);
        outfp = NULL;
    }

    if ('n' == t)
    {
        gs_name_index_type *name_index = NULL;

        n_names = calloc((size_t)((0 == npli) ? 1 : npli), sizeof(*n_names));
        if (NULL == n_names)
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }

        for (gs_UNSIGNED jj = 0; jj < npli; jj++)
        {
            if (st_text_set(&n_names[jj], gv->IndName[pli_idx[jj]]))
            {
                failed = 1;
                goto cleanup;
            }
            gs_str_rem_us(n_names[jj]);
        }

        if (0 < npli)
        {
            name_index = gs_malloc_array((size_t)npli, sizeof(*name_index));
            if (NULL == name_index)
            {
                gs_info(gv, -2, "Memory allocation error");
                failed = 1;
                goto cleanup;
            }

            for (gs_UNSIGNED jj = 0; jj < npli; jj++)
            {
                name_index[jj].name = n_names[jj];
                name_index[jj].idx = pli_idx[jj];
            }
            qsort(name_index, npli, sizeof(*name_index), gs_name_index_cmp);

            for (gs_UNSIGNED jj = 1; jj < npli; jj++)
                if ((0 ==
                     strcmp(name_index[jj - 1].name, name_index[jj].name)) &&
                    (name_index[jj - 1].idx != name_index[jj].idx))
                {
                    free(name_index);
                    gs_info(gv, -2, "Individual names collide in NTSys output");
                    failed = 1;
                    goto cleanup;
                }
            free(name_index);
        }

        if ((NULL == gv->All) || (NULL == gv->NoAll))
        {
            gs_info(gv, -2, "Marker statistics unavailable");
            failed = 1;
            goto cleanup;
        }

        if (NULL != gv->map)
        {
            for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
                if ((gv->mdl[lne] < gs_UNSIGNED_MAX) &&
                    (NULL != strchr(gv->MarName[gv->mdl[lne]], '.')))
                {
                    gs_info(gv, -2, "Marker name contains '.' in NTSys output");
                    failed = 1;
                    goto cleanup;
                }
        }
        else
        {
            for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
                if (NULL != strchr(gv->MarName[m], '.'))
                {
                    gs_info(gv, -2, "Marker name contains '.' in NTSys output");
                    failed = 1;
                    goto cleanup;
                }
        }

        if (NULL != gv->map)
        {
            outfp = fopen(*nmp_filename, "w");
            if (NULL == outfp)
            {
                gs_info(gv, -2, "File open error");
                failed = 1;
                goto cleanup;
            }

            for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
                if (gv->mdl[lne] < gs_UNSIGNED_MAX)
                    fprintf(outfp, "%i %15.12f %s %s\n", gv->map[lne]->chrom,
                            gv->map[lne]->pos / 100, gv->map[lne]->name,
                            gv->map[lne]->class);
            fclose(outfp);
            outfp = NULL;

            outfp = fopen(*ngp_filename, "w");
            if (NULL == outfp)
            {
                gs_info(gv, -2, "File open error");
                failed = 1;
                goto cleanup;
            }

            {
                gs_UNSIGNED last_lne = gs_UNSIGNED_MAX;

                for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
                    if (gv->mdl[lne] < gs_UNSIGNED_MAX)
                    {
                        if ((gs_UNSIGNED_MAX != last_lne) &&
                            (gv->map[lne]->chrom > gv->map[last_lne]->chrom))
                            fprintf(outfp, "%9.4f\n",
                                    0.01 + gv->map[last_lne]->pos / 100);
                        last_lne = lne;
                    }

                if (gs_UNSIGNED_MAX != last_lne)
                    fprintf(outfp, "%9.4f\n",
                            0.01 + gv->map[last_lne]->pos / 100);
            }

            fclose(outfp);
            outfp = NULL;
        }

        outfp = fopen(*npo_filename, "w");
        if (NULL == outfp)
        {
            gs_info(gv, -2, "File open error");
            failed = 1;
            goto cleanup;
        }

        for (gs_UNSIGNED jj = 0; jj < npli; jj++)
            fprintf(outfp, "%s%s ", n_names[jj], *add_inum);
        fprintf(outfp, "\n");

        if (NULL != gv->map)
        {
            for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
                if (gv->mdl[lne] < gs_UNSIGNED_MAX)
                {
                    gs_UNSIGNED m = gv->mdl[lne];

                    if (NULL != strchr(gv->MarName[m], '.'))
                    {
                        fclose(outfp);
                        outfp = NULL;
                        gs_info(gv, -2,
                                "Marker name contains '.' in NTSys output");
                        failed = 1;
                        goto cleanup;
                    }

                    for (gs_UNSIGNED a = 0; a < gv->NoAll[m]; a++)
                        if (GS_ALLELE_MISSING != gv->All[m][a])
                        {
                            fprintf(outfp, "%s.%u ", gv->MarName[m],
                                    (unsigned)gv->All[m][a]);
                            for (gs_UNSIGNED jj = 0; jj < npli; jj++)
                            {
                                gs_UNSIGNED i = pli_idx[jj];
                                size_t idx = (size_t)m * (size_t)gv->NoInd + i;
                                int present = ((gv->All[m][a] ==
                                                st_allele_get(gv, idx, 0)) ||
                                               (gv->All[m][a] ==
                                                st_allele_get(gv, idx, 1)))
                                                  ? 1
                                                  : 0;
                                fprintf(outfp, "%i ", present);
                            }
                            fprintf(outfp, "\n");
                        }
                }
        }
        else
        {
            for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            {
                for (gs_UNSIGNED a = 0; a < gv->NoAll[m]; a++)
                    if (GS_ALLELE_MISSING != gv->All[m][a])
                    {
                        fprintf(outfp, "%s.%u ", gv->MarName[m],
                                (unsigned)gv->All[m][a]);
                        for (gs_UNSIGNED jj = 0; jj < npli; jj++)
                        {
                            gs_UNSIGNED i = pli_idx[jj];
                            size_t idx = (size_t)m * (size_t)gv->NoInd + i;
                            int present =
                                ((gv->All[m][a] == st_allele_get(gv, idx, 0)) ||
                                 (gv->All[m][a] == st_allele_get(gv, idx, 1)))
                                    ? 1
                                    : 0;
                            fprintf(outfp, "%i ", present);
                        }
                        fprintf(outfp, "\n");
                    }
            }
        }

        fclose(outfp);
        outfp = NULL;
    }

cleanup:
    if (NULL != outfp)
        fclose(outfp);
    free(pli_idx);
    st_text_array_free(n_names, npli);

    if (failed)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
}

void gs_write_marker_data_02_GV(char **type, char **lpo_filename,
                                char **mmp_filename, char **mpo_filename,
                                char **ngp_filename, char **nmp_filename,
                                char **npo_filename, char **ind_filename,
                                int *f_ind, int *l_ind, char **add_inum,
                                int *auxfiles, int *retval, char **set_name)
{

    gs_write_marker_data_02(gs_fdta(*set_name), type, lpo_filename,
                            mmp_filename, mpo_filename, ngp_filename,
                            nmp_filename, npo_filename, ind_filename, f_ind,
                            l_ind, add_inum, auxfiles, retval);
}

void gs_write_map_01(gs_varset_type *gv, char **mmp_filename, int *retval)
{

    if (NULL == gv->map)
    {
        snprintf(gv->msg, sizeof(gv->msg), "No linkage map was loaded");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    gs_UNSIGNED lne;

    FILE *outfp;

    if ((outfp = fopen(*mmp_filename, "w")) == NULL)
        ERR_F;

    for (lne = 0; lne < gv->NLoc; lne++)
    {
        if (gv->mdl[lne] < gs_UNSIGNED_MAX)
        {
            snprintf(gv->msg, sizeof(gv->msg), "%i %15.12f %s %s\n",
                     (gv->map[lne])->chrom, (gv->map[lne])->pos,
                     (gv->map[lne])->name, (gv->map[lne])->class);

            fprintf(outfp, "%s", gv->msg);
        }
    }
    fclose(outfp);

    *retval = 0;
    return;
}

void gs_write_map_01_GV(char **mmp_filename, int *retval, char **set_name)
{

    gs_write_map_01(gs_fdta(*set_name), mmp_filename, retval);
}

int gs_determine_alleles(GS_ALLELE_TYPE *a, GS_ALLELE_TYPE *b, const char *p)
{
    const char *sep;
    size_t len;

    if ((NULL == a) || (NULL == b) || (NULL == p) || ('\0' == *p))
        return -1;

    len = strlen(p);

    if (gs_allele_missing_span(p, len))
    {
        *a = GS_ALLELE_MISSING;
        *b = GS_ALLELE_MISSING;
        return 0;
    }

    sep = strchr(p, '/');
    if (NULL != sep)
    {
        if (NULL != strchr(sep + 1, '/'))
            return -1;

        if ((0 != gs_parse_allele_span(p, (size_t)(sep - p), a)) ||
            (0 !=
             gs_parse_allele_span(sep + 1, len - (size_t)(sep - p) - 1, b)))
            return -1;

        if ((GS_ALLELE_MISSING == *a) || (GS_ALLELE_MISSING == *b))
        {
            *a = GS_ALLELE_MISSING;
            *b = GS_ALLELE_MISSING;
        }
        return 0;
    }

    if (2 == len)
    {
        if ((0 != gs_parse_allele_span(p, 1, a)) ||
            (0 != gs_parse_allele_span(p + 1, 1, b)))
            return -1;

        if ((GS_ALLELE_MISSING == *a) || (GS_ALLELE_MISSING == *b))
        {
            *a = GS_ALLELE_MISSING;
            *b = GS_ALLELE_MISSING;
        }
        return 0;
    }

    if (0 != gs_parse_single_allele(p, a))
        return -1;

    *b = *a;
    return 0;
}

void gs_remove_newline(char *p)
{
    char *pp = strchr(p, '\n');
    if (NULL != pp)
        *pp = '\0';
    pp = strchr(p, '\r');
    if (NULL != pp)
        *pp = '\0';
}

long long int st_getline(char **lineptr, size_t *n, FILE *stream)
{
    const size_t chunk = 32000;
    size_t len = 0;
    int c;

    if ((NULL == lineptr) || (NULL == n) || (NULL == stream))
        return -1;

    if ((NULL == *lineptr) || (0 == *n))
    {
        char *tmp = realloc(*lineptr, chunk);
        if (NULL == tmp)
            return -1;
        *lineptr = tmp;
        *n = chunk;
    }

    while (EOF != (c = fgetc(stream)))
    {
        if (len + 2 > *n)
        {
            size_t new_n;
            char *tmp;

            if (*n > ((size_t)-1) - chunk)
                return -1;

            new_n = *n + chunk;
            tmp = realloc(*lineptr, new_n);
            if (NULL == tmp)
                return -1;

            *lineptr = tmp;
            *n = new_n;
        }

        (*lineptr)[len++] = (char)c;
        if ('\n' == c)
            break;
    }

    if (0 == len)
    {
        (*lineptr)[0] = '\0';
        return 0;
    }

    if ('\r' == (*lineptr)[len - 1])
        (*lineptr)[len - 1] = '\n';

    (*lineptr)[len] = '\0';
    return (long long int)len;
}

void gs_read_marker_data_02(gs_varset_type *gv, char **dta_filename,
                            char **ind_filename, char **mar_filename,
                            char **gen_filename, int *auxfiles, int *retval)
{
    FILE *fp = NULL;
    char *line = NULL;
    size_t n = 0;
    long long int rc;
    gs_UNSIGNED header_fields = 0;
    gs_UNSIGNED row_fields = 0;
    gs_UNSIGNED no_ind = 0;
    gs_UNSIGNED no_mar = 0;
    int corner = 0;
    gs_varset_type tmp;
    gs_name_index_type *index = NULL;
    int failed = 0;

    gs_tmp_init_from(&tmp, gv);

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        *retval = -2;
        return;
    }

    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2,
                (rc < 0) ? "Memory allocation error"
                         : "Marker data file contains no header");
        failed = 1;
        goto cleanup;
    }
    gs_remove_newline(line);
    {
        char *p = strtok(line, gs_tok_del);
        while (NULL != p)
        {
            header_fields++;
            p = strtok(NULL, gs_tok_del);
        }
    }

    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2, "Marker data file contains no data rows");
        failed = 1;
        goto cleanup;
    }
    gs_remove_newline(line);
    {
        char *p = strtok(line, gs_tok_del);
        while (NULL != p)
        {
            row_fields++;
            p = strtok(NULL, gs_tok_del);
        }
    }

    if ((0 == header_fields) || (row_fields < 2))
    {
        gs_info(gv, -2, "Invalid marker matrix dimensions");
        failed = 1;
        goto cleanup;
    }

    if (row_fields == header_fields + 1)
    {
        corner = 0;
        no_ind = header_fields;
    }
    else if (row_fields == header_fields)
    {
        if (header_fields < 2)
        {
            gs_info(gv, -2, "Invalid marker matrix header");
            failed = 1;
            goto cleanup;
        }
        corner = 1;
        no_ind = header_fields - 1;
    }
    else
    {
        gs_info(gv, -2,
                "Marker matrix header and data row have different widths");
        failed = 1;
        goto cleanup;
    }

    fclose(fp);
    fp = NULL;

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        failed = 1;
        goto cleanup;
    }
    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2, "Error reading marker matrix header");
        failed = 1;
        goto cleanup;
    }

    while (0 < (rc = st_getline(&line, &n, fp)))
    {
        gs_UNSIGNED fields = 0;
        char *p;

        gs_remove_newline(line);
        p = strtok(line, gs_tok_del);
        while (NULL != p)
        {
            fields++;
            p = strtok(NULL, gs_tok_del);
        }

        if (fields != no_ind + 1)
        {
            gs_info(gv, -2, "Marker matrix row has wrong number of fields");
            failed = 1;
            goto cleanup;
        }
        no_mar++;
    }

    if (rc < 0)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    fclose(fp);
    fp = NULL;

    if (0 == no_mar)
    {
        gs_info(gv, -2, "Marker data file contains no marker rows");
        failed = 1;
        goto cleanup;
    }

    tmp.NoInd = no_ind;
    tmp.NoMar = no_mar;
    tmp.IndName = calloc((size_t)tmp.NoInd, sizeof(*tmp.IndName));
    tmp.MarName = gs_malloc_array((size_t)tmp.NoMar, sizeof(*tmp.MarName));
    if ((NULL == tmp.IndName) || (NULL == tmp.MarName) ||
        (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        failed = 1;
        goto cleanup;
    }

    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2, "Error reading marker matrix header");
        failed = 1;
        goto cleanup;
    }
    gs_remove_newline(line);
    {
        char *p = strtok(line, gs_tok_del);
        gs_UNSIGNED i = 0;

        if (corner)
            p = strtok(NULL, gs_tok_del);
        while (NULL != p)
        {
            size_t len = strlen(p);
            if ((i >= tmp.NoInd) || (len >= INT_MAX))
            {
                gs_info(gv, -2, "Individual name too long or invalid header");
                failed = 1;
                goto cleanup;
            }
            if (st_identity_set(&tmp, i, p, NULL))
            {
                failed = 1;
                goto cleanup;
            }
            if ((int)len > tmp.MaxIL)
                tmp.MaxIL = (int)len;
            i++;
            p = strtok(NULL, gs_tok_del);
        }

        if (i != tmp.NoInd)
        {
            gs_info(gv, -2, "Invalid marker matrix header");
            failed = 1;
            goto cleanup;
        }
    }

    if (0 != gs_make_ind_name_index(tmp.IndName, tmp.NoInd, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, tmp.NoInd))
    {
        gs_info(gv, -2, "Duplicate individual name in marker matrix");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    for (gs_UNSIGNED m = 0; m < tmp.NoMar; m++)
    {
        char *p;
        size_t len;

        rc = st_getline(&line, &n, fp);
        if (rc <= 0)
        {
            gs_info(gv, -2, "Unexpected end of marker matrix");
            failed = 1;
            goto cleanup;
        }
        gs_remove_newline(line);
        p = strtok(line, gs_tok_del);
        if (NULL == p)
        {
            gs_info(gv, -2, "Missing marker name");
            failed = 1;
            goto cleanup;
        }
        len = strlen(p);
        if (len >= gs_MAXNME)
        {
            gs_info(gv, -2, "Marker name too long");
            failed = 1;
            goto cleanup;
        }
        strcpy(tmp.MarName[m], p);
        if ((int)len > tmp.MaxML)
            tmp.MaxML = (int)len;

        for (gs_UNSIGNED i = 0; i < tmp.NoInd; i++)
        {
            GS_ALLELE_TYPE a, b;
            p = strtok(NULL, gs_tok_del);
            if (NULL == p)
            {
                gs_info(gv, -2, "Marker matrix row is too short");
                failed = 1;
                goto cleanup;
            }
            if (0 != gs_determine_alleles(&a, &b, p))
            {
                gs_info(gv, -2, "Invalid allele in marker matrix");
                failed = 1;
                goto cleanup;
            }
            tmp.allele_raw1[(size_t)m * tmp.NoInd + i] = a;
            tmp.allele_raw2[(size_t)m * tmp.NoInd + i] = b;
        }

        if (NULL != strtok(NULL, gs_tok_del))
        {
            gs_info(gv, -2, "Marker matrix row is too long");
            failed = 1;
            goto cleanup;
        }
    }

    fclose(fp);
    fp = NULL;

    if (0 != gs_make_name_index(tmp.MarName, tmp.NoMar, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, tmp.NoMar))
    {
        gs_info(gv, -2, "Duplicate marker name in marker matrix");
        failed = 1;
        goto cleanup;
    }

    if (0 != gs_finish_marker_read(gv, &tmp, ind_filename, mar_filename,
                                   gen_filename, auxfiles, retval))
        failed = 1;

cleanup:
    if (NULL != fp)
        fclose(fp);
    free(line);
    free(index);
    if (failed)
    {
        gs_reset(&tmp);
        *retval = -2;
    }
}

void gs_read_marker_data_02_GV(char **dta_filename, char **ind_filename,
                               char **mar_filename, char **gen_filename,
                               int *auxfiles, int *retval, char **set_name)
{
    gs_read_marker_data_02(gs_fdta(*set_name), dta_filename, ind_filename,
                           mar_filename, gen_filename, auxfiles, retval);
}

void gs_read_marker_data_04(gs_varset_type *gv, char **dta_filename,
                            char **ind_filename, char **mar_filename,
                            char **gen_filename, int *auxfiles, int *retval)
{
    const char sep = ';';
    FILE *fp = NULL;
    char *line = NULL;
    size_t n = 0;
    long long int rc;
    size_t header_fields;
    size_t data_fields;
    gs_UNSIGNED no_mar;
    gs_UNSIGNED no_ind = 0;
    int corner;
    char **fields = NULL;
    gs_varset_type tmp;
    gs_name_index_type *index = NULL;
    int failed = 0;

    gs_tmp_init_from(&tmp, gv);

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        *retval = -2;
        return;
    }

    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2, "Marker data file contains no header");
        failed = 1;
        goto cleanup;
    }
    gs_remove_newline(line);
    header_fields = gs_count_sep_fields(line, sep);

    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2, "Marker data file contains no data rows");
        failed = 1;
        goto cleanup;
    }
    gs_remove_newline(line);
    data_fields = gs_count_sep_fields(line, sep);

    if (data_fields == header_fields + 1)
    {
        corner = 0;
        no_mar = (gs_UNSIGNED)header_fields;
    }
    else if (data_fields == header_fields)
    {
        if (header_fields < 2)
        {
            gs_info(gv, -2, "Invalid transposed marker matrix header");
            failed = 1;
            goto cleanup;
        }
        corner = 1;
        no_mar = (gs_UNSIGNED)(header_fields - 1);
    }
    else
    {
        gs_info(gv, -2, "Transposed marker matrix has inconsistent width");
        failed = 1;
        goto cleanup;
    }

    if (0 == no_mar)
    {
        gs_info(gv, -2, "Marker matrix contains no markers");
        failed = 1;
        goto cleanup;
    }

    no_ind = 1;
    while (0 < (rc = st_getline(&line, &n, fp)))
    {
        gs_remove_newline(line);
        if (gs_count_sep_fields(line, sep) != (size_t)no_mar + 1)
        {
            gs_info(gv, -2, "Transposed marker matrix row has wrong width");
            failed = 1;
            goto cleanup;
        }
        no_ind++;
    }
    if (rc < 0)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    fclose(fp);
    fp = NULL;

    tmp.NoMar = no_mar;
    tmp.NoInd = no_ind;
    tmp.MarName = gs_malloc_array((size_t)tmp.NoMar, sizeof(*tmp.MarName));
    tmp.IndName = calloc((size_t)tmp.NoInd, sizeof(*tmp.IndName));
    fields = gs_malloc_array((size_t)tmp.NoMar + 1, sizeof(*fields));

    if ((NULL == tmp.MarName) || (NULL == tmp.IndName) || (NULL == fields) ||
        (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        failed = 1;
        goto cleanup;
    }

    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2, "Error reading transposed marker matrix header");
        failed = 1;
        goto cleanup;
    }
    gs_remove_newline(line);
    {
        size_t cap = (size_t)tmp.NoMar + (corner ? 1 : 0);
        char **header = gs_malloc_array(cap, sizeof(*header));
        size_t nf;
        size_t start = corner ? 1 : 0;

        if (NULL == header)
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }

        nf = gs_split_sep_fields(line, sep, header, cap);
        if (nf != cap)
        {
            free(header);
            gs_info(gv, -2, "Invalid transposed marker matrix header");
            failed = 1;
            goto cleanup;
        }

        for (gs_UNSIGNED m = 0; m < tmp.NoMar; m++)
        {
            char *name = gs_trim_field(header[start + m]);
            size_t len = strlen(name);

            if ((0 == len) || (len >= gs_MAXNME))
            {
                free(header);
                gs_info(gv, -2, "Invalid marker name in transposed matrix");
                failed = 1;
                goto cleanup;
            }
            strcpy(tmp.MarName[m], name);
            if ((int)len > tmp.MaxML)
                tmp.MaxML = (int)len;
        }
        free(header);
    }

    if (0 != gs_make_name_index(tmp.MarName, tmp.NoMar, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, tmp.NoMar))
    {
        gs_info(gv, -2, "Duplicate marker name in transposed matrix");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    for (gs_UNSIGNED i = 0; i < tmp.NoInd; i++)
    {
        size_t nf;
        char *name;
        size_t len;

        rc = st_getline(&line, &n, fp);
        if (rc <= 0)
        {
            gs_info(gv, -2, "Unexpected end of transposed marker matrix");
            failed = 1;
            goto cleanup;
        }
        gs_remove_newline(line);
        nf = gs_split_sep_fields(line, sep, fields, (size_t)tmp.NoMar + 1);
        if (nf != (size_t)tmp.NoMar + 1)
        {
            gs_info(gv, -2, "Transposed marker matrix row has wrong width");
            failed = 1;
            goto cleanup;
        }

        name = gs_trim_field(fields[0]);
        len = strlen(name);
        if ((0 == len) || (len >= INT_MAX))
        {
            gs_info(gv, -2, "Invalid individual name in transposed matrix");
            failed = 1;
            goto cleanup;
        }
        if (st_identity_set(&tmp, i, name, NULL))
        {
            failed = 1;
            goto cleanup;
        }
        if ((int)len > tmp.MaxIL)
            tmp.MaxIL = (int)len;

        for (gs_UNSIGNED m = 0; m < tmp.NoMar; m++)
        {
            char *geno = gs_trim_field(fields[m + 1]);
            GS_ALLELE_TYPE a, b;
            if (0 != gs_determine_alleles(&a, &b, geno))
            {
                gs_info(gv, -2, "Invalid allele in transposed marker matrix");
                failed = 1;
                goto cleanup;
            }
            tmp.allele_raw1[(size_t)m * tmp.NoInd + i] = a;
            tmp.allele_raw2[(size_t)m * tmp.NoInd + i] = b;
        }
    }

    fclose(fp);
    fp = NULL;

    if (0 != gs_make_ind_name_index(tmp.IndName, tmp.NoInd, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, tmp.NoInd))
    {
        gs_info(gv, -2, "Duplicate individual name in transposed matrix");
        failed = 1;
        goto cleanup;
    }

    if (0 != gs_finish_marker_read(gv, &tmp, ind_filename, mar_filename,
                                   gen_filename, auxfiles, retval))
        failed = 1;

cleanup:
    if (NULL != fp)
        fclose(fp);
    free(line);
    free(fields);
    free(index);
    if (failed)
    {
        gs_reset(&tmp);
        *retval = -2;
    }
}

void gs_read_marker_data_04_GV(char **dta_filename, char **ind_filename,
                               char **mar_filename, char **gen_filename,
                               int *auxfiles, int *retval, char **set_name)
{
    gs_read_marker_data_04(gs_fdta(*set_name), dta_filename, ind_filename,
                           mar_filename, gen_filename, auxfiles, retval);
}

void gs_marker_stats_01(gs_varset_type *gv, char **ind_filename,
                        char **mar_filename, char **gen_filename, int *ind,
                        int *mar, int *gen, int *auxfiles, int *retval)
{
    gs_UNSIGNED *NoAll = NULL;
    GS_ALLELE_TYPE *All_S = NULL;
    GS_ALLELE_TYPE **All = NULL;
    gs_UNSIGNED *AllCnt_S = NULL;
    gs_UNSIGNED **AllCnt = NULL;
    double *Pic = NULL;
    double *Mis = NULL;
    double *IndMis = NULL;
    GS_ALLELE_TYPE *all_work = NULL;
    FILE *outfp;
    size_t total_all = 0;
    size_t offset = 0;

    if ((0 == gv->NoInd) || (0 == gv->NoMar) || !st_allele_present(gv) ||
        (NULL == gv->IndName) || (NULL == gv->MarName))
    {
        gs_info(gv, -2, "No marker data available");
        *retval = -2;
        return;
    }

    if (!gs_marker_genotypes_valid(gv))
    {
        gs_info(gv, -2, "Invalid marker genotype representation");
        *retval = -2;
        return;
    }

    snprintf(gv->msg, sizeof(gv->msg),
             "No. of individuals: %lu, no. of markers: %lu", gv->NoInd,
             gv->NoMar);
    gs_info(gv, 0, gv->msg);

    NoAll = gs_malloc_array((size_t)gv->NoMar, sizeof(*NoAll));
    Pic = gs_malloc_array((size_t)gv->NoMar, sizeof(*Pic));
    Mis = gs_malloc_array((size_t)gv->NoMar, sizeof(*Mis));
    IndMis = gs_malloc_array((size_t)gv->NoInd, sizeof(*IndMis));
    all_work = gs_malloc_array((size_t)gs_MAXALL, sizeof(*all_work));

    if ((NULL == NoAll) || (NULL == Pic) || (NULL == Mis) || (NULL == IndMis) ||
        (NULL == all_work))
    {
        free(NoAll);
        free(Pic);
        free(Mis);
        free(IndMis);
        free(all_work);
        gs_info(gv, -2, "Memory allocation error");
        *retval = -2;
        return;
    }

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
        st_allele_marker_view view = st_allele_marker(gv, m);
        NoAll[m] = 0;

        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        {
            st_allele_pair pair = st_allele_read_pair(&view, i);
            GS_ALLELE_TYPE alleles[2] = {pair.a, pair.b};

            for (int h = 0; h < 2; h++)
            {
                int found = 0;

                for (gs_UNSIGNED k = 0; k < NoAll[m]; k++)
                    if (alleles[h] == all_work[k])
                    {
                        found = 1;
                        break;
                    }

                if (!found)
                {
                    if (NoAll[m] >= gs_MAXALL)
                    {
                        free(NoAll);
                        free(Pic);
                        free(Mis);
                        free(IndMis);
                        free(all_work);
                        gs_info(gv, -2, "Too many alleles per marker");
                        *retval = -2;
                        return;
                    }
                    all_work[NoAll[m]++] = alleles[h];
                }
            }
        }

        if ((size_t)NoAll[m] > ((size_t)-1) - total_all)
        {
            free(NoAll);
            free(Pic);
            free(Mis);
            free(IndMis);
            free(all_work);
            gs_info(gv, -2, "Memory allocation error");
            *retval = -2;
            return;
        }
        total_all += (size_t)NoAll[m];
    }

    All = gs_malloc_array((size_t)gv->NoMar, sizeof(*All));
    All_S = gs_malloc_array(total_all, sizeof(*All_S));
    AllCnt = gs_malloc_array((size_t)gv->NoMar, sizeof(*AllCnt));
    AllCnt_S = gs_malloc_array(total_all, sizeof(*AllCnt_S));

    if ((NULL == All) || (NULL == All_S) || (NULL == AllCnt) ||
        (NULL == AllCnt_S))
    {
        free(NoAll);
        free(All_S);
        free(All);
        free(AllCnt_S);
        free(AllCnt);
        free(Pic);
        free(Mis);
        free(IndMis);
        free(all_work);
        gs_info(gv, -2, "Memory allocation error");
        *retval = -2;
        return;
    }

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
        All[m] = All_S + offset;
        AllCnt[m] = AllCnt_S + offset;
        offset += (size_t)NoAll[m];
    }

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
        st_allele_marker_view view = st_allele_marker(gv, m);
        gs_UNSIGNED no_all = 0;

        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        {
            st_allele_pair pair = st_allele_read_pair(&view, i);
            GS_ALLELE_TYPE alleles[2] = {pair.a, pair.b};

            for (int h = 0; h < 2; h++)
            {
                int found = 0;

                for (gs_UNSIGNED k = 0; k < no_all; k++)
                    if (alleles[h] == all_work[k])
                    {
                        found = 1;
                        break;
                    }

                if (!found)
                {
                    if (no_all >= gs_MAXALL)
                    {
                        free(NoAll);
                        free(All_S);
                        free(All);
                        free(AllCnt_S);
                        free(AllCnt);
                        free(Pic);
                        free(Mis);
                        free(IndMis);
                        free(all_work);
                        gs_info(gv, -2, "Too many alleles per marker");
                        *retval = -2;
                        return;
                    }
                    all_work[no_all++] = alleles[h];
                }
            }
        }

        if (no_all != NoAll[m])
        {
            free(NoAll);
            free(All_S);
            free(All);
            free(AllCnt_S);
            free(AllCnt);
            free(Pic);
            free(Mis);
            free(IndMis);
            free(all_work);
            gs_info(gv, -2, "Marker statistics error");
            *retval = -2;
            return;
        }

        memcpy(All[m], all_work, (size_t)NoAll[m] * sizeof(*All[m]));
        qsort(All[m], NoAll[m], sizeof(*All[m]), gs_allele_cmp);

        for (gs_UNSIGNED k = 0; k < NoAll[m]; k++)
        {
            AllCnt[m][k] = 0;
            for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            {
                st_allele_pair pair = st_allele_read_pair(&view, i);
                if (pair.a == All[m][k])
                    AllCnt[m][k]++;
                if (pair.b == All[m][k])
                    AllCnt[m][k]++;
            }
        }

        {
            double observed = 2.0 * (double)gv->NoInd;
            gs_UNSIGNED real_no_all = NoAll[m];

            if ((0 < NoAll[m]) && (GS_ALLELE_MISSING == All[m][NoAll[m] - 1]))
            {
                gs_UNSIGNED km = NoAll[m] - 1;
                Mis[m] = (double)AllCnt[m][km] / 2.0 / (double)gv->NoInd;
                observed -= AllCnt[m][km];
                real_no_all--;
            }
            else
                Mis[m] = 0;

            if (observed > 0.0)
            {
                Pic[m] = 1.0;
                for (gs_UNSIGNED k = 0; k < real_no_all; k++)
                {
                    double f = (double)AllCnt[m][k] / observed;
                    Pic[m] -= f * f;
                }
            }
            else
                Pic[m] = 0.0;
        }
    }

    free(all_work);

    for (size_t i = 0; i < gv->NoInd; i++)
        IndMis[i] = 0;
    for (size_t m = 0; m < gv->NoMar; m++)
    {
        st_allele_marker_view view = st_allele_marker(gv, m);
        for (size_t i = 0; i < gv->NoInd; i++)
            if (st_allele_read_pair(&view, i).a == GS_ALLELE_MISSING)
                IndMis[i] += 1;
    }
    for (size_t i = 0; i < gv->NoInd; i++)
        IndMis[i] /= (double)gv->NoMar;

    free(gv->NoAll);
    gv->NoAll = NoAll;
    free(gv->All_S);
    gv->All_S = All_S;
    free(gv->All);
    gv->All = All;
    free(gv->AllCnt_S);
    gv->AllCnt_S = AllCnt_S;
    free(gv->AllCnt);
    gv->AllCnt = AllCnt;
    free(gv->Pic);
    gv->Pic = Pic;
    free(gv->Mis);
    gv->Mis = Mis;
    free(gv->IndMis);
    gv->IndMis = IndMis;

    if ((1 == *auxfiles) && (1 == *ind))
    {
        outfp = fopen(*ind_filename, "w");
        if (NULL == outfp)
        {
            gs_info(gv, -2, "File open error");
            *retval = -2;
            return;
        }

        fprintf(outfp, "Name InMis\n");
        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            fprintf(outfp, "%s %f\n", gv->IndName[i], gv->IndMis[i]);
        fclose(outfp);
    }

    if ((1 == *auxfiles) && (1 == *mar))
    {
        int *obs_all = NULL;
        gs_UNSIGNED no_obs_all = 0;
        gs_UNSIGNED capacity = 0;
        int have_missing = 0;

        for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            for (gs_UNSIGNED k = 0; k < gv->NoAll[m]; k++)
            {
                int allele = (int)gv->All[m][k];
                int found = 0;

                if (GS_ALLELE_MISSING == gv->All[m][k])
                {
                    have_missing = 1;
                    continue;
                }

                for (gs_UNSIGNED l = 0; l < no_obs_all; l++)
                    if (allele == obs_all[l])
                    {
                        found = 1;
                        break;
                    }

                if (!found)
                {
                    if (no_obs_all == capacity)
                    {
                        gs_UNSIGNED new_capacity =
                            (0 == capacity) ? 16 : capacity * 2;
                        int *new_obs;

                        if (new_capacity <= capacity)
                        {
                            free(obs_all);
                            gs_info(gv, -2, "Too many observed allele codes");
                            *retval = -2;
                            return;
                        }

                        new_obs = gs_realloc_array(
                            obs_all, (size_t)new_capacity, sizeof(*new_obs));
                        if (NULL == new_obs)
                        {
                            free(obs_all);
                            gs_info(gv, -2, "Memory allocation error");
                            *retval = -2;
                            return;
                        }
                        obs_all = new_obs;
                        capacity = new_capacity;
                    }

                    obs_all[no_obs_all++] = allele;
                }
            }

        if (no_obs_all > 1)
            qsort(obs_all, no_obs_all, sizeof(*obs_all), gs_intcmp);

        outfp = fopen(*mar_filename, "w");
        if (NULL == outfp)
        {
            free(obs_all);
            gs_info(gv, -2, "File open error");
            *retval = -2;
            return;
        }

        fprintf(outfp, "Name NoAll MaMis ExHet ");
        if (have_missing)
            fprintf(outfp, "AM ");
        for (gs_UNSIGNED l = 0; l < no_obs_all; l++)
            fprintf(outfp, " A%i", obs_all[l]);
        fprintf(outfp, "\n");

        for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
        {
            gs_UNSIGNED real_no_all = gv->NoAll[m];
            gs_UNSIGNED miss_count = 0;

            if ((0 < gv->NoAll[m]) &&
                (GS_ALLELE_MISSING == gv->All[m][gv->NoAll[m] - 1]))
            {
                real_no_all--;
                miss_count = gv->AllCnt[m][gv->NoAll[m] - 1];
            }

            fprintf(outfp, "%s %lu %5.3f %5.3f ", gv->MarName[m], real_no_all,
                    gv->Mis[m], gv->Pic[m]);

            if (have_missing)
                fprintf(outfp, " %lu", miss_count);

            for (gs_UNSIGNED l = 0; l < no_obs_all; l++)
            {
                gs_UNSIGNED ac = 0;
                for (gs_UNSIGNED k = 0; k < real_no_all; k++)
                    if (obs_all[l] == (int)gv->All[m][k])
                        ac = gv->AllCnt[m][k];
                fprintf(outfp, " %lu", ac);
            }
            fprintf(outfp, "\n");
        }

        fclose(outfp);
        free(obs_all);
    }

    if ((1 == *auxfiles) && (1 == *gen))
    {
        outfp = fopen(*gen_filename, "w");
        if (NULL == outfp)
        {
            gs_info(gv, -2, "File open error");
            *retval = -2;
            return;
        }

        fprintf(outfp, "Mar/Ind ");
        for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            fprintf(outfp, "%s ", gv->IndName[i]);
        fprintf(outfp, "\n");

        for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
        {
            fprintf(outfp, "%s ", gv->MarName[m]);
            for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            {
                size_t idx = (size_t)m * (size_t)gv->NoInd + i;
                if (GS_ALLELE_MISSING == st_allele_get(gv, idx, 0))
                    fprintf(outfp, "-- ");
                else
                    fprintf(outfp, "%u/%u ",
                            (unsigned)st_allele_get(gv, idx, 0),
                            (unsigned)st_allele_get(gv, idx, 1));
            }
            fprintf(outfp, "\n");
        }
        fclose(outfp);
    }

    *retval = 0;
}

void gs_marker_stats_01_GV(char **ind_filename, char **mar_filename,
                           char **gen_filename, int *ind, int *mar, int *gen,
                           int *auxfiles, int *retval, char **set_name)
{
    gs_marker_stats_01(gs_fdta(*set_name), ind_filename, mar_filename,
                       gen_filename, ind, mar, gen, auxfiles, retval);
}

void gs_read_marker_data_03(gs_varset_type *gv, char **dta_filename,
                            char **ind_filename, char **mar_filename,
                            char **gen_filename, int *auxfiles, int *retval)
{
    FILE *fp = NULL;
    char *line = NULL;
    size_t n = 0;
    long long int rc;
    gs_varset_type tmp;
    gs_UNSIGNED marker_capacity = 0;
    gs_name_index_type *index = NULL;
    int failed = 0;

    gs_tmp_init_from(&tmp, gv);

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        *retval = -2;
        return;
    }

    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2, "NTSys marker file contains no header");
        failed = 1;
        goto cleanup;
    }
    gs_remove_newline(line);

    {
        char *p = strtok(line, gs_tok_del);
        while (NULL != p)
        {
            size_t len = strlen(p);
            char **new_names;

            if (len >= INT_MAX)
            {
                gs_info(gv, -2, "Individual name too long");
                failed = 1;
                goto cleanup;
            }

            if (tmp.NoInd == marker_capacity)
            {
                gs_UNSIGNED new_capacity =
                    (0 == marker_capacity) ? 64 : marker_capacity * 2;
                if (new_capacity <= marker_capacity)
                {
                    gs_info(gv, -2, "Too many individuals");
                    failed = 1;
                    goto cleanup;
                }
                new_names = gs_realloc_array(tmp.IndName, (size_t)new_capacity,
                                             sizeof(*new_names));
                if (NULL == new_names)
                {
                    gs_info(gv, -2, "Memory allocation error");
                    failed = 1;
                    goto cleanup;
                }
                memset(new_names + marker_capacity, 0,
                       (size_t)(new_capacity - marker_capacity) *
                           sizeof(*new_names));
                tmp.IndName = new_names;
                marker_capacity = new_capacity;
            }

            if (st_identity_set(&tmp, tmp.NoInd++, p, NULL))
            {
                failed = 1;
                goto cleanup;
            }
            if ((int)len > tmp.MaxIL)
                tmp.MaxIL = (int)len;
            p = strtok(NULL, gs_tok_del);
        }
    }

    if (0 == tmp.NoInd)
    {
        gs_info(gv, -2, "NTSys marker file contains no individuals");
        failed = 1;
        goto cleanup;
    }

    if (0 != gs_make_ind_name_index(tmp.IndName, tmp.NoInd, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, tmp.NoInd))
    {
        gs_info(gv, -2, "Duplicate individual name in NTSys marker file");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    marker_capacity = 0;
    while (0 < (rc = st_getline(&line, &n, fp)))
    {
        char *p;
        char *dot;
        char marker[gs_MAXNME];
        size_t len;
        gs_UNSIGNED fields = 0;

        gs_remove_newline(line);
        p = strtok(line, gs_tok_del);
        if (NULL == p)
        {
            gs_info(gv, -2, "Empty line in NTSys marker file");
            failed = 1;
            goto cleanup;
        }

        dot = strchr(p, '.');
        if ((NULL == dot) || (dot == p) || ('\0' == dot[1]))
        {
            gs_info(gv, -2, "Invalid marker.allele field in NTSys file");
            failed = 1;
            goto cleanup;
        }

        len = (size_t)(dot - p);
        if (len >= gs_MAXNME)
        {
            gs_info(gv, -2, "Marker name too long");
            failed = 1;
            goto cleanup;
        }
        memcpy(marker, p, len);
        marker[len] = '\0';

        {
            long allele;
            if ((0 != gs_parse_long_strict(dot + 1, &allele)) || (allele < 1) ||
                (allele > GS_ALLELE_MAX))
            {
                gs_info(gv, -2, "Invalid allele number in NTSys file");
                failed = 1;
                goto cleanup;
            }
        }

        while (NULL != (p = strtok(NULL, gs_tok_del)))
        {
            long value;
            if ((0 != gs_parse_long_strict(p, &value)) ||
                !((0 == value) || (1 == value)))
            {
                gs_info(gv, -2, "NTSys incidence value must be 0 or 1");
                failed = 1;
                goto cleanup;
            }
            fields++;
        }

        if (fields != tmp.NoInd)
        {
            gs_info(gv, -2, "NTSys marker row has wrong number of fields");
            failed = 1;
            goto cleanup;
        }

        if ((0 == tmp.NoMar) ||
            (0 != strcmp(tmp.MarName[tmp.NoMar - 1], marker)))
        {
            char(*new_names)[gs_MAXNME];

            if (tmp.NoMar == marker_capacity)
            {
                gs_UNSIGNED new_capacity =
                    (0 == marker_capacity) ? 64 : marker_capacity * 2;
                if (new_capacity <= marker_capacity)
                {
                    gs_info(gv, -2, "Too many markers");
                    failed = 1;
                    goto cleanup;
                }
                new_names = gs_realloc_array(tmp.MarName, (size_t)new_capacity,
                                             sizeof(*new_names));
                if (NULL == new_names)
                {
                    gs_info(gv, -2, "Memory allocation error");
                    failed = 1;
                    goto cleanup;
                }
                tmp.MarName = new_names;
                marker_capacity = new_capacity;
            }

            strcpy(tmp.MarName[tmp.NoMar++], marker);
            if ((int)len > tmp.MaxML)
                tmp.MaxML = (int)len;
        }
    }

    if (rc < 0)
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    fclose(fp);
    fp = NULL;

    if (0 == tmp.NoMar)
    {
        gs_info(gv, -2, "NTSys marker file contains no marker rows");
        failed = 1;
        goto cleanup;
    }

    if (0 != gs_make_name_index(tmp.MarName, tmp.NoMar, &index))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    if (gs_name_index_has_duplicates(index, tmp.NoMar))
    {
        gs_info(gv, -2, "Marker block reappears in NTSys marker file");
        failed = 1;
        goto cleanup;
    }
    free(index);
    index = NULL;

    if (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    gs_set_genotypes_missing(&tmp);

    fp = fopen(*dta_filename, "r");
    if (NULL == fp)
    {
        gs_info(gv, -2, "File open error");
        failed = 1;
        goto cleanup;
    }
    rc = st_getline(&line, &n, fp);
    if (rc <= 0)
    {
        gs_info(gv, -2, "Error reading NTSys header");
        failed = 1;
        goto cleanup;
    }

    {
        gs_UNSIGNED current_m = gs_UNSIGNED_MAX;
        char current_name[gs_MAXNME] = "";

        while (0 < (rc = st_getline(&line, &n, fp)))
        {
            char *p;
            char *dot;
            char marker[gs_MAXNME];
            size_t len;
            long allele_l;
            GS_ALLELE_TYPE allele;

            gs_remove_newline(line);
            p = strtok(line, gs_tok_del);
            if (NULL == p)
            {
                gs_info(gv, -2, "Error in NTSys marker file");
                failed = 1;
                goto cleanup;
            }
            dot = strchr(p, '.');
            if ((NULL == dot) || (dot == p) || ('\0' == dot[1]))
            {
                gs_info(gv, -2, "Invalid marker.allele field in NTSys file");
                failed = 1;
                goto cleanup;
            }
            len = (size_t)(dot - p);
            if (len >= gs_MAXNME)
            {
                gs_info(gv, -2, "Marker name too long");
                failed = 1;
                goto cleanup;
            }
            memcpy(marker, p, len);
            marker[len] = '\0';
            if ((0 != gs_parse_long_strict(dot + 1, &allele_l)) ||
                (allele_l < 1) || (allele_l > GS_ALLELE_MAX))
            {
                gs_info(gv, -2, "Invalid allele number in NTSys file");
                failed = 1;
                goto cleanup;
            }
            allele = (GS_ALLELE_TYPE)allele_l;

            if ((gs_UNSIGNED_MAX == current_m) ||
                (0 != strcmp(current_name, marker)))
            {
                current_m++;
                if ((current_m >= tmp.NoMar) ||
                    (0 != strcmp(tmp.MarName[current_m], marker)))
                {
                    gs_info(gv, -2, "Invalid marker block order in NTSys file");
                    failed = 1;
                    goto cleanup;
                }
                strcpy(current_name, marker);
            }

            for (gs_UNSIGNED i = 0; i < tmp.NoInd; i++)
            {
                long present;
                size_t idx = (size_t)current_m * tmp.NoInd + i;

                p = strtok(NULL, gs_tok_del);
                if ((NULL == p) || (0 != gs_parse_long_strict(p, &present)) ||
                    !((0 == present) || (1 == present)))
                {
                    gs_info(gv, -2, "Invalid NTSys incidence row");
                    failed = 1;
                    goto cleanup;
                }

                if (1 == present)
                {
                    if ((GS_ALLELE_MISSING == st_allele_get(&(tmp), idx, 0)) &&
                        (GS_ALLELE_MISSING == st_allele_get(&(tmp), idx, 1)))
                    {
                        tmp.allele_raw1[idx] = allele;
                        tmp.allele_raw2[idx] = allele;
                    }
                    else if ((st_allele_get(&(tmp), idx, 0) == allele) ||
                             (st_allele_get(&(tmp), idx, 1) == allele))
                    {
                        continue;
                    }
                    else if (st_allele_get(&(tmp), idx, 0) ==
                             st_allele_get(&(tmp), idx, 1))
                    {
                        tmp.allele_raw2[idx] = allele;
                    }
                    else
                    {
                        gs_info(
                            gv, -2,
                            "More than two present alleles in NTSys genotype");
                        failed = 1;
                        goto cleanup;
                    }
                }
            }

            if (NULL != strtok(NULL, gs_tok_del))
            {
                gs_info(gv, -2, "NTSys marker row has wrong number of fields");
                failed = 1;
                goto cleanup;
            }
        }

        if (rc < 0)
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }
    }

    fclose(fp);
    fp = NULL;

    if (0 != gs_finish_marker_read(gv, &tmp, ind_filename, mar_filename,
                                   gen_filename, auxfiles, retval))
        failed = 1;

cleanup:
    if (NULL != fp)
        fclose(fp);
    free(line);
    free(index);
    if (failed)
    {
        gs_reset(&tmp);
        *retval = -2;
    }
}

void gs_read_marker_data_03_GV(char **dta_filename, char **ind_filename,
                               char **mar_filename, char **gen_filename,
                               int *auxfiles, int *retval, char **set_name)
{
    gs_read_marker_data_03(gs_fdta(*set_name), dta_filename, ind_filename,
                           mar_filename, gen_filename, auxfiles, retval);
}

double st_time_store[5];

void st_start_timer(int *i, double *time)
{
    st_time_store[*i] = *time;
}

void st_stop_timer(int *i, double *time)
{
    *time = st_time_store[*i];
}

int st_info_level = 0;

extern int INFOLEVEL;

void st_set_info_level(int *i)
{
    st_info_level = *i;
    INFOLEVEL = *i;
    gs_set_all_info_levels(i);
}

void st_get_info_level(int *i)
{
    *i = st_info_level;
}

void st_info(int level, char *message)
{
    if (st_info_level >= level)
    {
        if (!st_console_output_allowed())
            return;

        if (1 <= level)
            Rprintf("I: %s\n", message);

        else if (0 == level)
            Rprintf("M: %s\n", message);

        else if (-1 == level)
            Rprintf("W: %s\n", message);

        else if (-2 == level)
            Rprintf("E: %s\n", message);
    }
}

void st_r_info(int *level, char **message)
{
    st_info(*level, *message);
}

enum
{
    ST_CROSS_GD = 1, ST_CROSS_MU = 2, ST_CROSS_MI = 4,
    ST_CROSS_MA = 8, ST_CROSS_VA = 16, ST_CROSS_ES = 32
};

static uint64_t st_cross_hash(uint64_t hash, const void *data, size_t n)
{
    const unsigned char *bytes = data;
    for (size_t i = 0; i < n; i++)
        hash = (hash ^ bytes[i]) * UINT64_C(1099511628211);
    return hash;
}

static void st_cross_current(gs_varset_type *gv)
{
    uint64_t markers = UINT64_C(14695981039346656037),
             effects = markers, map = markers;
    markers = st_cross_hash(markers, &gv->NoInd, sizeof(gv->NoInd));
    markers = st_cross_hash(markers, &gv->NoMar, sizeof(gv->NoMar));
    for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        markers = st_cross_hash(markers, gv->IndName[i], strlen(gv->IndName[i]) + 1);
    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
        markers = st_cross_hash(markers, gv->MarName[m], strlen(gv->MarName[m]) + 1);
    if (st_allele_present(gv))
        for (size_t i = 0; i < (size_t)gv->NoInd * gv->NoMar; i++)
        {
            GS_ALLELE_TYPE a = st_allele_get(gv, i, 0), b = st_allele_get(gv, i, 1);
            markers = st_cross_hash(markers, &a, sizeof(a));
            markers = st_cross_hash(markers, &b, sizeof(b));
        }
    effects = st_cross_hash(effects, &gv->NoEff, sizeof(gv->NoEff));
    effects = st_cross_hash(effects, &gv->dimu, sizeof(gv->dimu));
    if (gv->u) effects = st_cross_hash(effects, gv->u, gv->dimu * sizeof(*gv->u));
    if (gv->EffAll) effects = st_cross_hash(effects, gv->EffAll, gv->NoEff * sizeof(*gv->EffAll));
    map = st_cross_hash(map, &gv->NLoc, sizeof(gv->NLoc));
    if (gv->map)
        for (gs_UNSIGNED i = 0; i < gv->NLoc; i++)
        {
            map = st_cross_hash(map, &gv->map[i]->chrom, sizeof(gv->map[i]->chrom));
            map = st_cross_hash(map, &gv->map[i]->pos, sizeof(gv->map[i]->pos));
            map = st_cross_hash(map, gv->map[i]->name, strlen(gv->map[i]->name) + 1);
            if (gv->mdl) map = st_cross_hash(map, &gv->mdl[i], sizeof(gv->mdl[i]));
        }
    if (markers != gv->cross_markers_signature)
    {
        GS_FREE(gv->crs_S);
        GS_FREE(gv->crs);
        gv->NCrs = 0;
        gv->cross_ready = 0;
    }
    if (effects != gv->cross_effects_signature)
        gv->cross_ready &= ST_CROSS_GD;
    if (map != gv->cross_map_signature)
        gv->cross_ready &= ~(ST_CROSS_VA | ST_CROSS_ES);
    gv->cross_markers_signature = markers;
    gv->cross_effects_signature = effects;
    gv->cross_map_signature = map;
}

void gs_cross_init_01(gs_varset_type *gv, int *retval)
{
    st_cross_current(gv);
    if (gv->NoInd < 2 || !gv->NoMar || !st_allele_present(gv))
    {
        gs_info(gv, -2, "At least two parents with marker data are required");
        *retval = -2;
        return;
    }
    if (!gv->crs || !gv->crs_S)
    {
        size_t pairs, bytes;
        if (gs_size_mul(gv->NoInd, gv->NoInd - 1, &pairs) ||
            pairs / 2 > gs_UNSIGNED_MAX ||
            gs_size_mul(pairs / 2, sizeof(*gv->crs_S), &bytes))
        {
            gs_info(gv, -2, "Too many parental pairs");
            *retval = -2;
            return;
        }
        pairs /= 2;
        gs_cross_type *storage = calloc(pairs, sizeof(*storage));
        gs_cross_type **order = gs_malloc_array(pairs, sizeof(*order));
        if (!storage || !order)
        {
            free(storage);
            free(order);
            gs_info(gv, -2, "Memory allocation error");
            *retval = -2;
            return;
        }
        gs_UNSIGNED k = 0;
        for (gs_UNSIGNED i = 0; i < gv->NoInd - 1; i++)
            for (gs_UNSIGNED j = i + 1; j < gv->NoInd; j++)
            {
                storage[k].p1 = i;
                storage[k].p2 = j;
                order[k] = &storage[k];
                k++;
            }
        GS_FREE(gv->crs_S);
        GS_FREE(gv->crs);
        gv->crs_S = storage;
        gv->crs = order;
        gv->NCrs = pairs;
        gv->cross_ready = 0;
    }
    *retval = 0;
}

void gs_calc_rd(gs_varset_type *gv, gs_UNSIGNED k)
{
    gs_UNSIGNED m, idx_i, idx_j;
    int a, b, c, d;

    gs_UNSIGNED i = (gv->crs_S[k]).p1;
    gs_UNSIGNED j = (gv->crs_S[k]).p2;

    gs_UNSIGNED n = 0;
    gs_FLOAT di = 0;

    for (m = 0; m < gv->NoMar; m++)
    {
        idx_i = m * gv->NoInd + i;
        a = st_allele_get(gv, idx_i, 0);
        b = st_allele_get(gv, idx_i, 1);

        idx_j = m * gv->NoInd + j;
        c = st_allele_get(gv, idx_j, 0);
        d = st_allele_get(gv, idx_j, 1);

        if (!((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b) ||
              (GS_ALLELE_MISSING == c) || (GS_ALLELE_MISSING == d)))
        {
            n++;

            if (!(((c == a) && (d == b)) || ((c == b) && (d == a))))
            {
                if (((c == a) || (c == b)) || ((d == a) || (d == b)))
                {
                    di += 0.5;
                }
                else
                {
                    di += 1;
                }
            }
        }
    }

    (gv->crs_S[k]).gd = (gs_FLOAT)di / (gs_FLOAT)n;
}

void gs_calc_mrd(gs_varset_type *gv, gs_UNSIGNED k)
{
    gs_UNSIGNED m, idx_i, idx_j;
    int a, b, c, d;

    gs_UNSIGNED i = (gv->crs_S[k]).p1;
    gs_UNSIGNED j = (gv->crs_S[k]).p2;

    gs_UNSIGNED n = 0;
    gs_FLOAT di = 0;

    for (m = 0; m < gv->NoMar; m++)
    {
        idx_i = m * gv->NoInd + i;
        a = st_allele_get(gv, idx_i, 0);
        b = st_allele_get(gv, idx_i, 1);

        idx_j = m * gv->NoInd + j;
        c = st_allele_get(gv, idx_j, 0);
        d = st_allele_get(gv, idx_j, 1);

        if (!((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b) ||
              (GS_ALLELE_MISSING == c) || (GS_ALLELE_MISSING == d)))
        {
            n++;

            if (!(((c == a) && (d == b)) || ((c == b) && (d == a))))
            {
                if (((c == a) || (c == b)) || ((d == a) || (d == b)))
                {
                    di += 0.5;
                }
                else
                {
                    di += 2;
                }
            }
        }
    }

    (gv->crs_S[k]).gd = (gs_FLOAT)sqrt(di / (gs_FLOAT)n / 2);
}

void gs_calc_euc(gs_varset_type *gv, gs_UNSIGNED k)
{
    gs_UNSIGNED m, idx_i, idx_j;
    int a, b, c, d;

    gs_UNSIGNED i = (gv->crs_S[k]).p1;
    gs_UNSIGNED j = (gv->crs_S[k]).p2;

    gs_UNSIGNED n = 0;
    gs_FLOAT di = 0;

    (void)n;

    for (m = 0; m < gv->NoMar; m++)
    {
        idx_i = m * gv->NoInd + i;
        a = st_allele_get(gv, idx_i, 0);
        b = st_allele_get(gv, idx_i, 1);

        idx_j = m * gv->NoInd + j;
        c = st_allele_get(gv, idx_j, 0);
        d = st_allele_get(gv, idx_j, 1);

        if (!((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b) ||
              (GS_ALLELE_MISSING == c) || (GS_ALLELE_MISSING == d)))
        {
            n++;

            if (!(((c == a) && (d == b)) || ((c == b) && (d == a))))
            {
                if (((c == a) || (c == b)) || ((d == a) || (d == b)))
                {
                    di += 0.5;
                }
                else
                {
                    di += 2;
                }
            }
        }
    }

    (gv->crs_S[k]).gd = (gs_FLOAT)sqrt(di);
}

void gs_calc_mu(gs_varset_type *gv, gs_UNSIGNED k)
{
    gs_UNSIGNED m, idx_i, idx_j;
    int a, b, c, d, e, count;

    gs_UNSIGNED i = (gv->crs_S[k]).p1;
    gs_UNSIGNED j = (gv->crs_S[k]).p2;

    gs_FLOAT *mu = &((gv->crs_S[k]).mu);

    *mu = gv->u[0];

    for (m = 0; m < gv->NoMar; m++)
    {
        idx_i = m * gv->NoInd + i;
        a = st_allele_get(gv, idx_i, 0);
        b = st_allele_get(gv, idx_i, 1);

        idx_j = m * gv->NoInd + j;
        c = st_allele_get(gv, idx_j, 0);
        d = st_allele_get(gv, idx_j, 1);

        e = gv->EffAll[1 + m];

        count = (int)(e == a) + (int)(e == b) + (int)(e == c) + (int)(e == d);
        *mu += (gs_FLOAT)count / (gs_FLOAT)2 * (gs_FLOAT)(gv->u[1 + m]);
    }
}

void gs_cross_eval_mu_01(gs_varset_type *gv, int *retval)
{

    if (NULL == gv->u)
    {
        gs_info(gv, -2, "No estimated effects available");
        *retval = -2;
        return;
    }

    int Srv = 0;
    int *Sretval = &Srv;
    gs_cross_init_01(gv, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_UNSIGNED k;

#ifdef MT
#pragma omp parallel for
#endif

    for (k = 0; k < gv->NCrs; k++)
    {
        gs_calc_mu(gv, k);
    }

    gv->cross_ready &= ~ST_CROSS_ES;
    gv->cross_ready |= ST_CROSS_MU;
    *retval = 0;
    return;
}

void gs_cross_eval_mu_01_GV(int *retval, char **set_name)
{
    gs_cross_eval_mu_01(gs_fdta(*set_name), retval);
}

void gs_calc_ma(gs_varset_type *gv, gs_UNSIGNED k)
{
    gs_UNSIGNED m, idx_i, idx_j;
    int a, b, c, d, e, count;

    gs_UNSIGNED i = (gv->crs_S[k]).p1;
    gs_UNSIGNED j = (gv->crs_S[k]).p2;

    gs_FLOAT *ma = &((gv->crs_S[k]).ma);

    *ma = gv->u[0];

    for (m = 0; m < gv->NoMar; m++)
    {
        idx_i = m * gv->NoInd + i;
        a = st_allele_get(gv, idx_i, 0);
        b = st_allele_get(gv, idx_i, 1);

        idx_j = m * gv->NoInd + j;
        c = st_allele_get(gv, idx_j, 0);
        d = st_allele_get(gv, idx_j, 1);

        e = gv->EffAll[1 + m];

        count = (int)(e == a) + (int)(e == b) + (int)(e == c) + (int)(e == d);

        if ((0 < count) && (0 < gv->u[1 + m]))
            *ma += (gs_FLOAT)2 * (gs_FLOAT)(gv->u[1 + m]);
    }
}

void gs_cross_eval_ma_01(gs_varset_type *gv, int *retval)
{

    if (NULL == gv->u)
    {
        gs_info(gv, -2, "No estimated effects available");
        *retval = -2;
        return;
    }

    int Srv = 0;
    int *Sretval = &Srv;
    gs_cross_init_01(gv, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_UNSIGNED k;

#ifdef MT
#pragma omp parallel for
#endif

    for (k = 0; k < gv->NCrs; k++)
    {
        gs_calc_ma(gv, k);
    }

    gv->cross_ready |= ST_CROSS_MA;
    *retval = 0;
    return;
}

void gs_cross_eval_ma_01_GV(int *retval, char **set_name)
{
    gs_cross_eval_ma_01(gs_fdta(*set_name), retval);
}

void gs_calc_mi(gs_varset_type *gv, gs_UNSIGNED k)
{
    gs_UNSIGNED m, idx_i, idx_j;
    int a, b, c, d, e, count;

    gs_UNSIGNED i = (gv->crs_S[k]).p1;
    gs_UNSIGNED j = (gv->crs_S[k]).p2;

    gs_FLOAT *mi = &((gv->crs_S[k]).mi);

    *mi = gv->u[0];

    for (m = 0; m < gv->NoMar; m++)
    {
        idx_i = m * gv->NoInd + i;
        a = st_allele_get(gv, idx_i, 0);
        b = st_allele_get(gv, idx_i, 1);

        idx_j = m * gv->NoInd + j;
        c = st_allele_get(gv, idx_j, 0);
        d = st_allele_get(gv, idx_j, 1);

        e = gv->EffAll[1 + m];

        count = (int)(e == a) + (int)(e == b) + (int)(e == c) + (int)(e == d);

        if ((0 < count) && (0 > gv->u[1 + m]))
            *mi += (gs_FLOAT)2 * (gs_FLOAT)(gv->u[1 + m]);
    }
}

void gs_cross_eval_mi_01(gs_varset_type *gv, int *retval)
{

    if (NULL == gv->u)
    {
        gs_info(gv, -2, "No estimated effects available");
        *retval = -2;
        return;
    }

    int Srv = 0;
    int *Sretval = &Srv;
    gs_cross_init_01(gv, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_UNSIGNED k;

#ifdef MT
#pragma omp parallel for
#endif

    for (k = 0; k < gv->NCrs; k++)
    {
        gs_calc_mi(gv, k);
    }

    gv->cross_ready |= ST_CROSS_MI;
    *retval = 0;
    return;
}

void gs_cross_eval_mi_01_GV(int *retval, char **set_name)
{
    gs_cross_eval_mi_01(gs_fdta(*set_name), retval);
}

void gs_cross_eval_gd_01(gs_varset_type *gv, char **dist, int *retval)
{
    if (strcmp(*dist, "rd") && strcmp(*dist, "mrd") && strcmp(*dist, "euc"))
    {
        gs_info(gv, -2, "Unknown genetic distance");
        *retval = -2;
        return;
    }


    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2,
                "No marker data available, run 'read.marker.data' first");
        *retval = -2;
        return;
    }

    int Srv = 0;
    int *Sretval = &Srv;
    gs_cross_init_01(gv, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_UNSIGNED k;

    if (0 == strcmp("mrd", *dist))
    {
#ifdef MT
#pragma omp parallel for
#endif

        for (k = 0; k < gv->NCrs; k++)
            gs_calc_mrd(gv, k);
    }
    else if (0 == strcmp("rd", *dist))
    {
#ifdef MT
#pragma omp parallel for
#endif

        for (k = 0; k < gv->NCrs; k++)
            gs_calc_rd(gv, k);
    }
    else
    {
#ifdef MT
#pragma omp parallel for
#endif

        for (k = 0; k < gv->NCrs; k++)
            gs_calc_euc(gv, k);
    }

    gv->cross_ready |= ST_CROSS_GD;
    *retval = 0;
    return;
}

void gs_cross_eval_gd_01_GV(char **dist, int *retval, char **set_name)
{
    gs_cross_eval_gd_01(gs_fdta(*set_name), dist, retval);
}

void gs_calc_genv_01(gs_varset_type *gv)

{
    gs_UNSIGNED ii, jj, idx;
    int a, b, e, count;

    gv->genv = realloc(gv->genv, gv->NoMar * gv->NoInd * sizeof(gs_FLOAT));

    for (ii = 0; ii < gv->NoMar; ii++)
    {
        e = gv->EffAll[1 + ii];

        for (jj = 0; jj < gv->NoInd; jj++)
        {
            idx = ii * gv->NoInd + jj;
            a = st_allele_get(gv, idx, 0);
            b = st_allele_get(gv, idx, 1);
            count = (int)(e == a) + (int)(e == b);
            gv->genv[idx] = (gs_FLOAT)count * (gs_FLOAT)(gv->u[1 + ii]) / 2;
        }
    }
}

void gs_calc_va_linked(gs_varset_type *gv, gs_UNSIGNED k)
{
    gs_UNSIGNED c, ii, jj, m, idx, idx_i, idx_j;

    gs_UNSIGNED i = (gv->crs_S[k]).p1;
    gs_UNSIGNED j = (gv->crs_S[k]).p2;

    gs_FLOAT q, ggu, ggv, hhu, hhv;

    gs_FLOAT *va = &((gv->crs_S[k]).va);

    *va = 0;

    for (c = 0; c < gv->NChrom; c++)
    {
        for (ii = 0; ii < gv->chrom[c].NLocC; ii++)
            for (jj = 0; jj < gv->chrom[c].NLocC; jj++)
            {
                idx = ii * gv->chrom[c].NLocC + jj;
                q = (gv->chrom[c].q)[idx];

                m = ii + gv->chrom[c].first_loc_i;
                idx_i = m * gv->NoInd + i;
                idx_j = m * gv->NoInd + j;
                ggu = gv->genv[idx_i];
                hhu = gv->genv[idx_j];

                m = jj + gv->chrom[c].first_loc_i;
                idx_i = m * gv->NoInd + i;
                idx_j = m * gv->NoInd + j;
                ggv = gv->genv[idx_i];
                hhv = gv->genv[idx_j];

                *va += (q / 2 - 0.25) *
                       (ggu * ggv + hhu * hhv - ggu * hhv - hhu * ggv);
            }
    }

    *va = *va * 4;
}

void gs_calc_genv_02(gs_varset_type *gv)
{
    gs_UNSIGNED ii, jj, idx;
    int a, b, e, count1, count2;

    gv->genv1 = realloc(gv->genv1, gv->NoMar * gv->NoInd * sizeof(gs_FLOAT));
    gv->genv2 = realloc(gv->genv2, gv->NoMar * gv->NoInd * sizeof(gs_FLOAT));

    for (ii = 0; ii < gv->NoMar; ii++)
    {
        e = gv->EffAll[1 + ii];

        for (jj = 0; jj < gv->NoInd; jj++)
        {
            idx = ii * gv->NoInd + jj;
            a = st_allele_get(gv, idx, 0);
            b = st_allele_get(gv, idx, 1);
            count1 = (int)(e == a);
            count2 = (int)(e == b);
            gv->genv1[idx] = (gs_FLOAT)count1 * (gs_FLOAT)(gv->u[1 + ii]);
            gv->genv2[idx] = (gs_FLOAT)count2 * (gs_FLOAT)(gv->u[1 + ii]);
        }
    }
}

void gs_calc_va_crs(gs_varset_type *gv, gs_UNSIGNED k)
{
    gs_UNSIGNED c, ii, jj, m, idx, idx_i, idx_j;

    gs_UNSIGNED i = (gv->crs_S[k]).p1;
    gs_UNSIGNED j = (gv->crs_S[k]).p2;

    gs_FLOAT q, ggu, ggv, hhu, hhv;

    gs_FLOAT *va = &((gv->crs_S[k]).va);

    *va = 0;

    for (c = 0; c < gv->NChrom; c++)
    {
        for (ii = 0; ii < gv->chrom[c].NLocC; ii++)
            for (jj = 0; jj < gv->chrom[c].NLocC; jj++)
            {
                idx = ii * gv->chrom[c].NLocC + jj;
                q = (gv->chrom[c].q)[idx];

                m = ii + gv->chrom[c].first_loc_i;
                idx_i = m * gv->NoInd + i;
                idx_j = m * gv->NoInd + j;
                ggu = gv->genv1[idx_i];
                hhu = gv->genv2[idx_i];

                m = jj + gv->chrom[c].first_loc_i;
                idx_i = m * gv->NoInd + i;
                idx_j = m * gv->NoInd + j;
                ggv = gv->genv1[idx_i];
                hhv = gv->genv2[idx_i];

                *va += (q / 2 - 0.25) *
                       (ggu * ggv + hhu * hhv - ggu * hhv - hhu * ggv);

                m = ii + gv->chrom[c].first_loc_i;
                idx_i = m * gv->NoInd + i;
                idx_j = m * gv->NoInd + j;
                ggu = gv->genv1[idx_j];
                hhu = gv->genv2[idx_j];

                m = jj + gv->chrom[c].first_loc_i;
                idx_i = m * gv->NoInd + i;
                idx_j = m * gv->NoInd + j;
                ggv = gv->genv1[idx_j];
                hhv = gv->genv2[idx_j];

                *va += (q / 2 - 0.25) *
                       (ggu * ggv + hhu * hhv - ggu * hhv - hhu * ggv);
            }
    }
}

void gs_calc_va(gs_varset_type *gv, gs_UNSIGNED k)
{
    gs_UNSIGNED m, idx_i, idx_j;
    int a, b, c, d, e, count;

    gs_UNSIGNED i = (gv->crs_S[k]).p1;
    gs_UNSIGNED j = (gv->crs_S[k]).p2;

    gs_FLOAT EZ2 = 0;
    gs_FLOAT EZ = 0;

    gs_FLOAT *va = &((gv->crs_S[k]).va);

    *va = 0;

    for (m = 0; m < gv->NoMar; m++)
    {
        idx_i = m * gv->NoInd + i;
        a = st_allele_get(gv, idx_i, 0);
        b = st_allele_get(gv, idx_i, 1);

        idx_j = m * gv->NoInd + j;
        c = st_allele_get(gv, idx_j, 0);
        d = st_allele_get(gv, idx_j, 1);

        e = gv->EffAll[1 + m];

        count = (int)(e == a) + (int)(e == b) + (int)(e == c) + (int)(e == d);

        EZ2 = (gs_FLOAT)count * (gs_FLOAT)(gv->u[1 + m]) *
              (gs_FLOAT)(gv->u[1 + m]) / 4;
        EZ = (gs_FLOAT)count * (gs_FLOAT)(gv->u[1 + m]) / 4;

        *va += (EZ2 - EZ * EZ) * 4;
    }
}

void gs_cross_eval_va_01(gs_varset_type *gv, char **pop_type, int *t,
                         char **map_function, int *retval)
{
    if (*t < 0 ||
        (strcmp(*pop_type, "unlinked") && strcmp(*pop_type, "DH") &&
         strcmp(*pop_type, "SSD") && strcmp(*pop_type, "CRS")) ||
        (strcmp(*pop_type, "unlinked") && strcmp(*map_function, "Haldane") &&
         strcmp(*map_function, "Kosambi")))
    {
        gs_info(gv, -2, "Invalid cross population, generation or map function");
        *retval = -2;
        return;
    }


    if (NULL == gv->u)
    {
        gs_info(gv, -2, "No estimated effects available");
        *retval = -2;
        return;
    }

    int Srv = 0;
    int *Sretval = &Srv;
    gs_cross_init_01(gv, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    char *empty = "";
    char **Sout_filename = &empty;
    int auf = 0;
    int *Sauxfiles = &auf;

    gs_UNSIGNED k;

    if (0 == strcmp("unlinked", *pop_type))
    {
#ifdef MT
#pragma omp parallel for
#endif

        for (k = 0; k < gv->NCrs; k++)
        {
            gs_calc_va(gv, k);
        }
    }

    else if ((0 == strcmp(*pop_type, "DH")) || (0 == strcmp(*pop_type, "SSD")))
    {

        gs_chrom_stats_01(gv, Sout_filename, Sauxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_calc_rf_01(gv, map_function, Sout_filename, Sauxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_calc_q_01(gv, pop_type, t, Sout_filename, Sauxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_calc_genv_01(gv);

#ifdef MT
#pragma omp parallel for
#endif

        for (k = 0; k < gv->NCrs; k++)
        {
            gs_calc_va_linked(gv, k);
        }
    }

    else if (0 == strcmp(*pop_type, "CRS"))
    {

        gs_chrom_stats_01(gv, Sout_filename, Sauxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_calc_rf_01(gv, map_function, Sout_filename, Sauxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_calc_q_01(gv, pop_type, t, Sout_filename, Sauxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_calc_genv_02(gv);

#ifdef MT
#pragma omp parallel for
#endif

        for (k = 0; k < gv->NCrs; k++)
        {
            gs_calc_va_crs(gv, k);
        }
    }

    else
    {
        gs_info(gv, -2, "Unknown population type");
        *retval = -2;
        return;
    }

    gv->cross_ready &= ~ST_CROSS_ES;
    gv->cross_ready |= ST_CROSS_VA;
    *retval = 0;
    return;
}

void gs_cross_eval_va_01_GV(char **pop_type, int *t, char **map_function,
                            int *retval, char **set_name)
{
    gs_cross_eval_va_01(gs_fdta(*set_name), pop_type, t, map_function, retval);
}

void gs_calc_es(gs_varset_type *gv, gs_UNSIGNED k, double *i)
{
    gs_FLOAT *mu = &((gv->crs_S[k]).mu);
    gs_FLOAT *va = &((gv->crs_S[k]).va);
    gs_FLOAT *es = &((gv->crs_S[k]).es);

    *es = *mu + *i * sqrt(*va);
}

void gs_cross_eval_es_01(gs_varset_type *gv, double *i, int *retval)
{
    if (!R_FINITE(*i) || *i < 0)
    {
        gs_info(gv, -2, "Selection intensity must be finite and nonnegative");
        *retval = -2;
        return;
    }
    st_cross_current(gv);
    if (!gv->crs || !gv->crs_S ||
        (gv->cross_ready & (ST_CROSS_MU | ST_CROSS_VA)) != (ST_CROSS_MU | ST_CROSS_VA))
    {
        gs_info(gv, -2, "Estimated cross means and variances are required");
        *retval = -2;
        return;
    }
    for (gs_UNSIGNED k = 0; k < gv->NCrs; k++)
        if (!R_FINITE(gv->crs_S[k].mu) || !R_FINITE(gv->crs_S[k].va) ||
            gv->crs_S[k].va < 0 ||
            !R_FINITE(gv->crs_S[k].mu + *i * sqrt(gv->crs_S[k].va)))
        {
            gs_info(gv, -2, "Cross means and nonnegative variances must be finite");
            *retval = -2;
            return;
        }
    snprintf(gv->msg, sizeof(gv->msg),
             "Standardized selection differential %5.2f", *i);
    gs_info(gv, 0, gv->msg);
#ifdef MT
#pragma omp parallel for
#endif
    for (gs_UNSIGNED k = 0; k < gv->NCrs; k++)
        gs_calc_es(gv, k, i);
    gv->cross_ready |= ST_CROSS_ES;
    *retval = 0;
}

void gs_cross_eval_es_01_GV(double *i, int *retval, char **set_name)
{
    gs_cross_eval_es_01(gs_fdta(*set_name), i, retval);
}

static int st_cross_cmp_parents(const gs_cross_type *a, const gs_cross_type *b)
{
    if (a->p1 != b->p1) return a->p1 < b->p1 ? -1 : 1;
    if (a->p2 != b->p2) return a->p2 < b->p2 ? -1 : 1;
    return 0;
}

static int st_cross_cmp_value(const gs_cross_type *a, const gs_cross_type *b,
                              double x, double y, int ascending)
{
    if (ISNAN(x) != ISNAN(y)) return ISNAN(x) ? 1 : -1;
    if (x < y) return ascending ? -1 : 1;
    if (x > y) return ascending ? 1 : -1;
    return st_cross_cmp_parents(a, b);
}

int gs_cross_cmp_gd(const void *a, const void *b)
{
    const gs_cross_type *x = *(const gs_cross_type *const *)a,
                        *y = *(const gs_cross_type *const *)b;
    return st_cross_cmp_value(x, y, x->gd, y->gd, 0);
}

int gs_cross_cmp_mu(const void *a, const void *b)
{
    const gs_cross_type *x = *(const gs_cross_type *const *)a,
                        *y = *(const gs_cross_type *const *)b;
    return st_cross_cmp_value(x, y, x->mu, y->mu, 0);
}

int gs_cross_cmp_ma(const void *a, const void *b)
{
    const gs_cross_type *x = *(const gs_cross_type *const *)a,
                        *y = *(const gs_cross_type *const *)b;
    return st_cross_cmp_value(x, y, x->ma, y->ma, 0);
}

int gs_cross_cmp_mi(const void *a, const void *b)
{
    const gs_cross_type *x = *(const gs_cross_type *const *)a,
                        *y = *(const gs_cross_type *const *)b;
    return st_cross_cmp_value(x, y, x->mi, y->mi, 1);
}

int gs_cross_cmp_va(const void *a, const void *b)
{
    const gs_cross_type *x = *(const gs_cross_type *const *)a,
                        *y = *(const gs_cross_type *const *)b;
    return st_cross_cmp_value(x, y, x->va, y->va, 0);
}

int gs_cross_cmp_es(const void *a, const void *b)
{
    const gs_cross_type *x = *(const gs_cross_type *const *)a,
                        *y = *(const gs_cross_type *const *)b;
    return st_cross_cmp_value(x, y, x->es, y->es, 0);
}

int gs_cross_cmp_index(const void *a, const void *b)
{
    return st_cross_cmp_parents(*(const gs_cross_type *const *)a,
                                *(const gs_cross_type *const *)b);
}

void gs_cross_info_01(gs_varset_type *gv, int *bestn, char **sortby,
                      char **out_filename, int *retval)
{
    gs_UNSIGNED i, m;
    FILE *outfp;

    if (NULL == gv->crs)
    {
        gs_info(gv, -2, "No cross information available");
        *retval = -2;
        return;
    }

    if ((outfp = fopen(*out_filename, "w")) == NULL)
        ERR_F;

    if (0 == strcmp(*sortby, "index"))
        qsort(gv->crs, gv->NCrs, sizeof(gs_cross_type *), gs_cross_cmp_index);

    if (0 == strcmp(*sortby, "gd"))
        qsort(gv->crs, gv->NCrs, sizeof(gs_cross_type *), gs_cross_cmp_gd);

    if (0 == strcmp(*sortby, "mu"))
        qsort(gv->crs, gv->NCrs, sizeof(gs_cross_type *), gs_cross_cmp_mu);

    if (0 == strcmp(*sortby, "mi"))
        qsort(gv->crs, gv->NCrs, sizeof(gs_cross_type *), gs_cross_cmp_mi);

    if (0 == strcmp(*sortby, "ma"))
        qsort(gv->crs, gv->NCrs, sizeof(gs_cross_type *), gs_cross_cmp_ma);

    if (0 == strcmp(*sortby, "va"))
        qsort(gv->crs, gv->NCrs, sizeof(gs_cross_type *), gs_cross_cmp_va);

    if (0 == strcmp(*sortby, "es"))
        qsort(gv->crs, gv->NCrs, sizeof(gs_cross_type *), gs_cross_cmp_es);

    if (0 == *bestn)
        m = gv->NCrs;
    else
        m = ((gs_UNSIGNED)(*bestn) < gv->NCrs) ? (gs_UNSIGNED)(*bestn)
                                               : gv->NCrs;

    fprintf(outfp, "P1No P2No P1Name P2Name gd mu mi ma va es\n");
    for (i = 0; i < m; i++)
        fprintf(outfp, "%lu %lu %s %s %.17g %.17g %.17g %.17g %.17g %.17g\n",
                1 + (gv->crs[i])->p1, 1 + (gv->crs[i])->p2,
                gv->IndName[(gv->crs[i])->p1], gv->IndName[(gv->crs[i])->p2],
                (double)(gv->crs[i])->gd, (double)(gv->crs[i])->mu,
                (double)(gv->crs[i])->mi, (double)(gv->crs[i])->ma,
                (double)(gv->crs[i])->va, (double)(gv->crs[i])->es);

    fclose(outfp);

    *retval = 0;
    return;
}

void gs_cross_info_01_GV(int *bestn, char **sortby, char **out_filename,
                         int *retval, char **set_name)
{
    gs_cross_info_01(gs_fdta(*set_name), bestn, sortby, out_filename, retval);
}

extern void ps_get_simpop(char **pop_name, gs_UNSIGNED *NoInd,
                          gs_UNSIGNED *NoMar, char ***IndName, char ***Pedigree,
                          char (**MarName)[gs_MAXNME], GS_ALLELE_TYPE *(*All1),
                          GS_ALLELE_TYPE *(*All2), gs_UNSIGNED *NLoc,
                          gs_mappoint_type **map_S, char **map_str_S,
                          size_t *map_str_N, gs_mappoint_type ***map,
                          gs_UNSIGNED **mdl, int *bg, int *Sretval);

static void st_get_simpop_raw(gs_varset_type *gv, char **pop_name, int *bg,
                              int *retval)
{
    gs_reset(gv);

    gs_UNSIGNED *NoInd = &(gv->NoInd);
    gs_UNSIGNED *NoMar = &(gv->NoMar);

    char ***IndName = &(gv->IndName);
    char ***Pedigree = &(gv->Pedigree);
    char(**MarName)[gs_MAXNME] = &(gv->MarName);

    GS_ALLELE_TYPE *(*All1) = &(gv->allele_raw1);
    GS_ALLELE_TYPE *(*All2) = &(gv->allele_raw2);

    gs_UNSIGNED *NLoc = &(gv->NLoc);
    gs_mappoint_type **map_S = &(gv->map_S);
    char **map_str_S = &(gv->map_str_S);
    size_t *map_str_N = &(gv->map_str_N);
    gs_mappoint_type ***map = &(gv->map);
    gs_UNSIGNED **mdl = &(gv->mdl);

    int Srv = 0;
    int *Sretval = &Srv;
    ps_get_simpop(pop_name, NoInd, NoMar, IndName, Pedigree, MarName, All1,
                  All2, NLoc, map_S, map_str_S, map_str_N, map, mdl, bg,
                  Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    {
        int Srv = 0;
        int *Sretval = &Srv;
        char *empty = "";
        char **ind_filename = &empty;
        char **mar_filename = &empty;
        char **gen_filename = &empty;
        int auf = 0;
        int *Sauxfiles = &auf;
        gs_marker_stats_01(gv, ind_filename, mar_filename, gen_filename,
                           Sauxfiles, Sauxfiles, Sauxfiles, Sauxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }

    gv->phased = 1;

    *retval = 1;
    return;
}

void st_get_simpop(gs_varset_type *dst, char **pop_name, int *bg, int *retval)
{
    gs_varset_type tmp;
    gs_tmp_init_from(&tmp, dst);
    st_get_simpop_raw(&tmp, pop_name, bg, retval);
    if (*retval >= 0 && gs_tmp_commit(dst, &tmp))
        *retval = -2;
    gs_reset(&tmp);
}

void st_get_simpop_GV(char **pop_name, int *bg, int *retval, char **set_name)
{
    st_get_simpop(gs_fdta(*set_name), pop_name, bg, retval);
}

extern void ps_get_simpop2(char **pop_name, gs_UNSIGNED *NoInd,
                           gs_UNSIGNED *NoMar, char ***IndName,
                           char ***Pedigree, char (**MarName)[gs_MAXNME],
                           GS_ALLELE_TYPE *(*All1), GS_ALLELE_TYPE *(*All2),
                           gs_UNSIGNED *NLoc, gs_mappoint_type **map_S,
                           gs_mappoint_type ***map, gs_UNSIGNED **mdl, int *bg,
                           int *Sretval);

static void st_get_simpop2_raw(gs_varset_type *gv, char **pop_name, int *bg,
                               int *retval)
{
    gs_reset(gv);

    gs_UNSIGNED *NoInd = &(gv->NoInd);
    gs_UNSIGNED *NoMar = &(gv->NoMar);

    char ***IndName = &(gv->IndName);
    char ***Pedigree = &(gv->Pedigree);
    char(**MarName)[gs_MAXNME] = &(gv->MarName);

    GS_ALLELE_TYPE *(*All1) = &(gv->allele_raw1);
    GS_ALLELE_TYPE *(*All2) = &(gv->allele_raw2);

    gs_UNSIGNED *NLoc = &(gv->NLoc);
    gs_mappoint_type **map_S = &(gv->map_S);
    gs_mappoint_type ***map = &(gv->map);
    gs_UNSIGNED **mdl = &(gv->mdl);

    int Srv = 0;
    int *Sretval = &Srv;
    ps_get_simpop2(pop_name, NoInd, NoMar, IndName, Pedigree, MarName, All1,
                   All2, NLoc, map_S, map, mdl, bg, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    {
        int Srv = 0;
        int *Sretval = &Srv;
        char *empty = "";
        char **ind_filename = &empty;
        char **mar_filename = &empty;
        char **gen_filename = &empty;
        int auf = 0;
        int *Sauxfiles = &auf;
        gs_marker_stats_01(gv, ind_filename, mar_filename, gen_filename,
                           Sauxfiles, Sauxfiles, Sauxfiles, Sauxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }

    gv->phased = 1;

    *retval = 1;
    return;
}

void st_get_simpop2(gs_varset_type *dst, char **pop_name, int *bg, int *retval)
{
    gs_varset_type tmp;
    gs_tmp_init_from(&tmp, dst);
    st_get_simpop2_raw(&tmp, pop_name, bg, retval);
    if (*retval >= 0 && gs_tmp_commit(dst, &tmp))
        *retval = -2;
    gs_reset(&tmp);
}

void st_get_simpop2_GV(char **pop_name, int *bg, int *retval, char **set_name)
{
    st_get_simpop2(gs_fdta(*set_name), pop_name, bg, retval);
}

void il_get_population_size(gs_varset_type *gv, int *NoLines)
{
    *NoLines = (int)gv->NoInd;
}

void il_get_population_size_GV(int *NoLines, char **set_name)
{
    il_get_population_size(gs_fdta(*set_name), NoLines);
}

void il_calc_segment_length(gs_varset_type *gv)
{
    gs_UNSIGNED l, c, m;
    gs_UNSIGNED i1, i2, i3;
    gs_FLOAT p1, p2, p3;

    gv->dimu = gv->NLoc;
    gv->u = realloc(gv->u, gv->dimu * sizeof(double));

    for (c = 0, l = 0; c < gv->NChrom; c++)
    {
        m = 0;
        i2 = gv->chrom[c].first_loc_i + m;
        p2 = gv->map[i2]->pos;
        i3 = gv->chrom[c].first_loc_i + m + 1;
        p3 = gv->map[i3]->pos;
        gv->u[l] = (p3 - p2) / 2;
        l++;

        for (m = 1; m < (-1 + gv->chrom[c].NLocC); m++, l++)
        {
            i1 = gv->chrom[c].first_loc_i + m - 1;
            p1 = gv->map[i1]->pos;
            i2 = gv->chrom[c].first_loc_i + m;
            p2 = gv->map[i2]->pos;
            i3 = gv->chrom[c].first_loc_i + m + 1;
            p3 = gv->map[i3]->pos;
            gv->u[l] = (p2 - p1) / 2 + (p3 - p2) / 2;
        }

        i1 = gv->chrom[c].first_loc_i + m - 1;
        p1 = gv->map[i1]->pos;
        i2 = gv->chrom[c].first_loc_i + m;
        p2 = gv->map[i2]->pos;
        gv->u[l] = (p2 - p1) / 2;
        l++;
    }
}

void il_mark_loci(gs_varset_type *gv, int *chrom, double *lower, double *upper,
                  int *exclude)
{

    gs_UNSIGNED l, c, m;
    gs_UNSIGNED i2;
    gs_FLOAT p2;

    int hit = 0;

    if (0 == *chrom)
        return;

    for (c = 0, l = 0; c < gv->NChrom; c++)
    {
        for (m = 0; m < gv->chrom[c].NLocC; m++, l++)
        {
            i2 = gv->chrom[c].first_loc_i + m;
            p2 = gv->map[i2]->pos;

            hit = ((c == (gs_UNSIGNED)(-1 + *chrom)) &&
                   ((0 == *upper) || ((*lower <= p2) && (p2 <= *upper))));

            if (1 == *exclude)
                hit = !hit;

            if (!hit)
                gv->u[l] = 0;
        }
    }
}

double il_calc_tot(gs_varset_type *gv)
{
    gs_UNSIGNED l, m;

    int cc, s_c, e_c;
    s_c = 0;
    e_c = (int)gv->NChrom;

    double tot = 0;

    for (cc = s_c; cc < e_c; cc++)
    {
        for (m = 0; m < gv->chrom[cc].NLocC; m++)
        {
            l = gv->chrom[cc].first_loc_i + m;

            tot += gv->u[l];
        }
    }
    tot = tot * 2;
    return (tot);
}

void il_eval_lines(gs_varset_type *gv, int *allele, int *chrom, double *lower,
                   double *upper, int *exclude, double *abs, double *tot,
                   double *nseg, double *lseg, int *retval)
{
    gs_UNSIGNED i, l, m, idx;
    gs_UNSIGNED l1, l2, idx1, idx2;

    int cc, s_c, e_c;

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2, "No marker data available");
        *retval = -2;
        return;
    }

    if ((*allele < 1) || (*allele > GS_ALLELE_MAX))
    {
        gs_info(gv, -2, "Allele code must be in the GS allele range");
        *retval = -2;
        return;
    }

    gs_chrom_stats_01(gv, Sout_filename, Sauxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    s_c = 0;
    e_c = (int)gv->NChrom;

    if ((*chrom < 0) || (*chrom > (int)gv->NChrom))
    {
        gs_info(gv, -2, "Incorrect chromosome number");
        *retval = -2;
        return;
    }

    for (i = 0; i < gv->NoInd; i++)
    {

        il_calc_segment_length(gv);

        il_mark_loci(gv, &(chrom[i]), &(lower[i]), &(upper[i]), exclude);

        tot[i] = il_calc_tot(gv);

        abs[i] = 0;

        for (cc = s_c; cc < e_c; cc++)
        {
            for (m = 0; m < gv->chrom[cc].NLocC; m++)
            {
                l = gv->chrom[cc].first_loc_i + m;

                idx = l * gv->NoInd + i;
                abs[i] += (st_allele_get(gv, idx, 0) == *allele) * gv->u[l] +
                          (st_allele_get(gv, idx, 1) == *allele) * gv->u[l];
            }
        }

        nseg[i] = 0;
        for (cc = s_c; cc < e_c; cc++)
        {
            m = 0;
            l2 = gv->chrom[cc].first_loc_i + m;
            idx2 = l2 * gv->NoInd + i;
            if (0 != gv->u[l2])
            {
                if (*allele == st_allele_get(gv, idx2, 0))
                    nseg[i]++;
                if (*allele == st_allele_get(gv, idx2, 1))
                    nseg[i]++;
            }

            for (m = 1; m < gv->chrom[cc].NLocC; m++)
            {
                l2 = gv->chrom[cc].first_loc_i + m;
                l1 = l2 - 1;

                idx1 = l1 * gv->NoInd + i;
                idx2 = l2 * gv->NoInd + i;

                if ((0 != gv->u[l1]) && ((0 != gv->u[l2])))
                {
                    if ((st_allele_get(gv, idx1, 0) !=
                         st_allele_get(gv, idx2, 0)) &&
                        (*allele == st_allele_get(gv, idx2, 0)))
                        nseg[i]++;

                    if ((st_allele_get(gv, idx1, 1) !=
                         st_allele_get(gv, idx2, 1)) &&
                        (*allele == st_allele_get(gv, idx2, 1)))
                        nseg[i]++;
                }
                else if ((0 == gv->u[l1]) && ((0 != gv->u[l2])))
                {
                    if (*allele == st_allele_get(gv, idx2, 0))
                        nseg[i]++;
                    if (*allele == st_allele_get(gv, idx2, 1))
                        nseg[i]++;
                }
            }
        }
        if (0 != nseg[i])
            lseg[i] = abs[i] / nseg[i];
    }

    *retval = 0;
}

void il_eval_lines_GV(int *allele, int *chrom, double *lower, double *upper,
                      int *exclude, double *abs, double *tot, double *nseg,
                      double *lseg, int *retval, char **set_name)
{
    il_eval_lines(gs_fdta(*set_name), allele, chrom, lower, upper, exclude, abs,
                  tot, nseg, lseg, retval);
}

void il_eval_library(gs_varset_type *gv, double *cov, double *dep,
                     double *seglib, int *allele, int *retval)
{
    gs_UNSIGNED i, j, m, cc, l1, l2, idx, idx1, idx2, cnt;
    double tot;

    if ((*allele < 1) || (*allele > GS_ALLELE_MAX))
    {
        gs_info(gv, -2, "Allele code must be in the GS allele range");
        *retval = -2;
        return;
    }

    gs_chrom_stats_01(gv, Sout_filename, Sauxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    il_calc_segment_length(gv);

    tot = 0;
    *cov = 0;
    for (i = 0; i < gv->NoMar; i++)
    {
        tot += gv->u[i];
        for (j = 0; j < gv->NoInd; j++)
        {
            idx = i * gv->NoInd + j;
            if ((*allele == st_allele_get(gv, idx, 0)) ||
                (*allele == st_allele_get(gv, idx, 1)))
            {
                *cov += gv->u[i];
                break;
            }
        }
    }
    *cov = *cov / tot;

    gs_UNSIGNED *depth = malloc(gv->NoMar * sizeof(gs_UNSIGNED));
    if (NULL == depth)
        ERR_M;

    for (i = 0; i < gv->NoMar; i++)
    {
        depth[i] = 0;

        for (j = 0; j < gv->NoInd; j++)
        {
            idx = i * gv->NoInd + j;
            if ((*allele == st_allele_get(gv, idx, 0)) ||
                (*allele == st_allele_get(gv, idx, 1)))
            {
                depth[i]++;
            }
        }
    }

    *dep = 0;
    cnt = 0;
    for (i = 0; i < gv->NoMar; i++)
    {
        if (depth[i] > 0)
        {
            *dep += (double)depth[i];
            cnt++;
        }
    }
    *dep = *dep / (double)cnt;

    free(depth);

    *seglib = (double)gv->NChrom;

    for (cc = 0; cc < gv->NChrom; cc++)
        for (m = 1; m < gv->chrom[cc].NLocC; m++)
            for (i = 0; i < gv->NoInd; i++)
            {
                l2 = gv->chrom[cc].first_loc_i + m;
                l1 = l2 - 1;

                idx1 = l1 * gv->NoInd + i;
                idx2 = l2 * gv->NoInd + i;

                if ((st_allele_get(gv, idx1, 0) !=
                     st_allele_get(gv, idx2, 0)) ||
                    (st_allele_get(gv, idx1, 1) != st_allele_get(gv, idx2, 1)))
                {
                    (*seglib)++;
                    break;
                }
            }
}

void il_eval_library_GV(double *cov, double *dep, double *seglib, int *allele,
                        int *retval, char **set_name)
{
    il_eval_library(gs_fdta(*set_name), cov, dep, seglib, allele, retval);
}

void gs_build_GI_01(gs_varset_type *gv, int *retval)
{
    const double small = 1e-5;

    gs_UNSIGNED ii, jj;
    double d;
    double det;

    const gs_UNSIGNED n = gv->NoInd;
    const gs_UNSIGNED m = gv->NoMar;

    gv->t9 = realloc(gv->t9, (size_t)n * (size_t)n * sizeof(double));
    gv->GI = realloc(gv->GI, (size_t)n * (size_t)n * sizeof(double));
    gv->u0 = realloc(gv->u0, (size_t)n * sizeof(gs_UNSIGNED));
    if ((NULL == gv->t9) || (NULL == gv->GI) || (NULL == gv->u0))
        ERR_M;

    st_dsymm_ZZt(n, m, gv->Z, gv->t9);

    memcpy(gv->GI, gv->t9, (size_t)n * (size_t)n * sizeof(double));

    {

        int info = st_dpotrf_rm_sym(gv->t9, (int)n);
        if (info == 0)
        {
            info = st_dpotri_rm_sym(gv->t9, (int)n);
            if (info == 0)
            {

                symmetrize_from_lower_rm(gv->t9, (int)n);

                memcpy(gv->GI, gv->t9, (size_t)n * (size_t)n * sizeof(double));

                *retval = 0;
                return;
            }
        }
    }

    memcpy(gv->t9, gv->GI, (size_t)n * (size_t)n * sizeof(double));
    st_dgetrf(gv->t9, n, gv->u0, &d);

    for (ii = 0, det = 1.0; ii < n; ++ii)
    {
        if (fabs(gv->t9[(size_t)ii * (size_t)n + ii]) < 1e-10)
        {
            det = 0.0;
            break;
        }
    }

    if (det != 0.0)
    {
        st_dgetri(gv->t9, n, gv->u0, gv->GI);
        *retval = 0;
        return;
    }

    gs_info(gv, -1, "ZZ' has not full rank, using approximation)");

    double *w = (double *)malloc((size_t)n * sizeof(double));
    double *z = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    if (!w || !z)
    {
        free(z);
        free(w);
        *retval = -1;
        return;
    }

    memcpy(gv->t9, gv->GI, (size_t)n * (size_t)n * sizeof(double));
    st_dsyevr_v(gv->t9, w, z, n);

    for (jj = 0; jj < n; ++jj)
    {
        double denom = (w[jj] < small) ? small : w[jj];
        double dinv = 1.0 / denom;
        for (ii = 0; ii < n; ++ii)
        {
            gv->t9[(size_t)ii * (size_t)n + jj] =
                z[(size_t)ii * (size_t)n + jj] * dinv;
        }
    }

    st_dgemm("N", "T", n, n, n, n, gv->t9, z, gv->GI);

    free(z);
    free(w);

    *retval = 0;
    return;
}

void gs_build_GI_02(gs_varset_type *gv, int *retval)
{

    double small = 1e-5;

    gv->t9 = realloc(gv->t9, gv->NoInd * gv->NoInd * sizeof(double));
    gv->GI = realloc(gv->GI, gv->NoInd * gv->NoInd * sizeof(double));
    if ((NULL == gv->t9) || (NULL == gv->GI))
        ERR_M;

    st_dgemm("N", "T", gv->NoInd, gv->NoMar, gv->NoInd, gv->NoMar, gv->Z, gv->Z,
             gv->t9);

    double d;
    gs_UNSIGNED ii, jj;
    double det = 0;
    gv->u0 = realloc(gv->u0, gv->NoInd * sizeof(gs_UNSIGNED));

    st_dgetrf(gv->t9, gv->NoInd, gv->u0, &d);

    for (ii = 0, det = 1; ii < gv->NoInd; ii++)
        if (fabs(gv->t9[ii * gv->NoInd + ii]) < 1e-10)
            det = 0.0;

    if (0.0 == det)
    {
        gs_info(gv, -1, "ZZ' has not full rank, using approximation)");

        double *w = malloc(gv->NoInd * sizeof(double));
        double *z = malloc(gv->NoInd * gv->NoInd * sizeof(double));

        double *t6 = malloc(gv->NoInd * gv->NoInd * sizeof(double));
        double *t7 = malloc(gv->NoInd * gv->NoInd * sizeof(double));
        double *t8 = malloc(gv->NoInd * gv->NoInd * sizeof(double));
        double *t9 = malloc(gv->NoInd * gv->NoInd * sizeof(double));

        st_dgemm("N", "T", gv->NoInd, gv->NoMar, gv->NoInd, gv->NoMar, gv->Z,
                 gv->Z, t9);

        while (0.0 == det)

        {

            memcpy(gv->t9, t9, gv->NoInd * gv->NoInd * sizeof(double));
            st_dsyevr_v(gv->t9, w, z, gv->NoInd);

            for (ii = 0; ii < gv->NoInd; ii++)
                for (jj = 0; jj < gv->NoInd; jj++)
                    t6[ii * gv->NoInd + jj] = 0;

            for (ii = 0; ii < gv->NoInd; ii++)
                t6[ii * gv->NoInd + ii] = (w[ii] < small) ? small : w[ii];

            st_dgemm("N", "N", gv->NoInd, gv->NoInd, gv->NoInd, gv->NoInd, z,
                     t6, t7);

            st_dgemm("N", "T", gv->NoInd, gv->NoInd, gv->NoInd, gv->NoInd, t7,
                     z, t8);

            for (ii = 0; ii < gv->NoInd; ii++)
                for (jj = 0; jj < gv->NoInd; jj++)
                    t9[ii * gv->NoInd + jj] =
                        t9[ii * gv->NoInd + jj] -
                        (t9[ii * gv->NoInd + jj] - t8[ii * gv->NoInd + jj]);

            memcpy(gv->t9, t9, gv->NoInd * gv->NoInd * sizeof(double));
            st_dgetrf(gv->t9, gv->NoInd, gv->u0, &d);

            for (ii = 0, det = 1; ii < gv->NoInd; ii++)
                if (fabs(gv->t9[ii * gv->NoInd + ii]) < 1e-10)
                    det = 0.0;
        }

        free(w);
        free(z);
        free(t6);
        free(t7);
        free(t8);
        free(t9);
    }

    st_dgetri(gv->t9, gv->NoInd, gv->u0, gv->GI);

    *retval = 0;
    return;
}

static void check_len_rowmajor(SEXP x, R_xlen_t need, const char *name)
{
    if (!isReal(x))
        error("%s must be a double vector", name);
    if (XLENGTH(x) < need)
        error("%s too short: need %lld, got %lld", name, (long long)need,
              (long long)XLENGTH(x));
}

SEXP call_gs_build_GI_01_R(SEXP Z_rm_, SEXP NoInd_, SEXP NoMar_)
{
    if (!isInteger(NoInd_) || LENGTH(NoInd_) != 1)
        error("NoInd must be integer scalar");
    if (!isInteger(NoMar_) || LENGTH(NoMar_) != 1)
        error("NoMar must be integer scalar");

    int NoInd = asInteger(NoInd_);
    int NoMar = asInteger(NoMar_);
    if (NoInd <= 0 || NoMar <= 0)
        error("NoInd/NoMar must be > 0");

    R_xlen_t lenZ = (R_xlen_t)NoInd * (R_xlen_t)NoMar;
    check_len_rowmajor(Z_rm_, lenZ, "Z_rm");

    gs_varset_type gv;
    memset(&gv, 0, sizeof(gv));

    gv.NoInd = (gs_UNSIGNED)NoInd;
    gv.NoMar = (gs_UNSIGNED)NoMar;

    gv.Z = (double *)malloc((size_t)lenZ * sizeof(double));
    if (!gv.Z)
        error("Out of memory allocating gv.Z");
    memcpy(gv.Z, REAL(Z_rm_), (size_t)lenZ * sizeof(double));

    gv.t9 = NULL;
    gv.GI = NULL;
    gv.u0 = NULL;

    int retval = -999;
    gs_build_GI_01(&gv, &retval);
    if (retval != 0)
    {
        free(gv.Z);
        free(gv.t9);
        free(gv.GI);
        free(gv.u0);
        error("gs_build_GI_01 failed (retval=%d)", retval);
    }

    R_xlen_t lenGI = (R_xlen_t)NoInd * (R_xlen_t)NoInd;
    SEXP out = PROTECT(allocVector(REALSXP, lenGI));
    memcpy(REAL(out), gv.GI, (size_t)lenGI * sizeof(double));

    free(gv.Z);
    free(gv.t9);
    free(gv.GI);
    free(gv.u0);

    UNPROTECT(1);
    return out;
}

SEXP call_gs_build_GI_02_R(SEXP Z_rm_, SEXP NoInd_, SEXP NoMar_)
{
    if (!isInteger(NoInd_) || LENGTH(NoInd_) != 1)
        error("NoInd must be integer scalar");
    if (!isInteger(NoMar_) || LENGTH(NoMar_) != 1)
        error("NoMar must be integer scalar");

    int NoInd = asInteger(NoInd_);
    int NoMar = asInteger(NoMar_);
    if (NoInd <= 0 || NoMar <= 0)
        error("NoInd/NoMar must be > 0");

    R_xlen_t lenZ = (R_xlen_t)NoInd * (R_xlen_t)NoMar;
    check_len_rowmajor(Z_rm_, lenZ, "Z_rm");

    gs_varset_type gv;
    memset(&gv, 0, sizeof(gv));

    gv.NoInd = (gs_UNSIGNED)NoInd;
    gv.NoMar = (gs_UNSIGNED)NoMar;

    gv.Z = (double *)malloc((size_t)lenZ * sizeof(double));
    if (!gv.Z)
        error("Out of memory allocating gv.Z");
    memcpy(gv.Z, REAL(Z_rm_), (size_t)lenZ * sizeof(double));

    gv.t9 = NULL;
    gv.GI = NULL;
    gv.u0 = NULL;

    int retval = -999;
    gs_build_GI_02(&gv, &retval);
    if (retval != 0)
    {
        free(gv.Z);
        free(gv.t9);
        free(gv.GI);
        free(gv.u0);
        error("gs_build_GI_02 failed (retval=%d)", retval);
    }

    R_xlen_t lenGI = (R_xlen_t)NoInd * (R_xlen_t)NoInd;
    SEXP out = PROTECT(allocVector(REALSXP, lenGI));
    memcpy(REAL(out), gv.GI, (size_t)lenGI * sizeof(double));

    free(gv.Z);
    free(gv.t9);
    free(gv.GI);
    free(gv.u0);

    UNPROTECT(1);
    return out;
}

static int gs_rr_boundary_01(gs_varset_type *gv, const double *weight,
                             int *retval)
{
    int all_zero = 1, all_infinite = 1;
    if (!gv->ladi || !gv->NoMar)
        return 0;
    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
        if (gv->ladi[m] != 0)
            all_zero = 0;
        if (R_FINITE(gv->ladi[m]) || gv->ladi[m] <= 0 || ISNAN(gv->ladi[m]))
            all_infinite = 0;
    }
    if (!all_zero && !all_infinite)
        return 0;
    if (!gv->y || !gv->Z || !gv->NoInd || gv->NoInd >= INT_MAX)
    {
        gs_info(gv, -2, "Invalid data for a boundary variance fit");
        *retval = -2;
        return -1;
    }
    size_t n = gv->NoInd, m = gv->NoMar, d = n + 1;
    if (m >= SIZE_MAX / sizeof(double) || d > SIZE_MAX / d / sizeof(double))
    {
        gs_info(gv, -2, "Boundary variance workspace is too large");
        *retval = -2;
        return -1;
    }
    double *u = realloc(gv->u, (m + 1) * sizeof(double));
    if (!u)
    {
        gs_info(gv, -2, "Memory allocation error in boundary variance fitting");
        *retval = -2;
        return -1;
    }
    gv->u = u;
    gv->dimu = m + 1;
    if (all_infinite)
    {
        long double mean = 0;
        for (size_t i = 0; i < n; i++)
            mean += gv->y[i];
        gv->u[0] = (double)(mean / n);
        for (size_t j = 1; j <= m; j++)
            gv->u[j] = 0;
    }
    else
    {
        double *A = calloc(d * d, sizeof(double));
        double *G = malloc(n * n * sizeof(double));
        double *rhs = calloc(d, sizeof(double));
        double *sol = malloc(d * sizeof(double));
        if (!A || !G || !rhs || !sol)
        {
            free(A);
            free(G);
            free(rhs);
            free(sol);
            gs_info(gv, -2,
                    "Memory allocation error in boundary variance fitting");
            *retval = -2;
            return -1;
        }
        double scale = 1;
        if (weight)
        {
            scale = 0;
            for (size_t j = 0; j < m; j++)
            {
                if (!R_FINITE(weight[j]) || weight[j] < 0)
                {
                    free(A);
                    free(G);
                    free(rhs);
                    free(sol);
                    gs_info(
                        gv, -2,
                        "Invalid marker weights for boundary variance fitting");
                    *retval = -2;
                    return -1;
                }
                if (weight[j] > scale)
                    scale = weight[j];
            }
            if (!(scale > 0))
            {
                free(A);
                free(G);
                free(rhs);
                free(sol);
                gs_info(
                    gv, -2,
                    "No positive marker weights for boundary variance fitting");
                *retval = -2;
                return -1;
            }
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j <= i; j++)
                {
                    double value = 0;
                    for (size_t h = 0; h < m; h++)
                        value += gv->Z[i * m + h] * (weight[h] / scale) *
                                 gv->Z[j * m + h];
                    G[i * n + j] = G[j * n + i] = value;
                }
        }
        else
            st_dgemm("N", "T", n, m, n, m, gv->Z, gv->Z, G);
        for (size_t i = 0; i < n; i++)
        {
            A[i + 1] = A[(i + 1) * d] = 1;
            rhs[i + 1] = gv->y[i];
            for (size_t j = 0; j < n; j++)
                A[(i + 1) * d + j + 1] = G[i * n + j];
        }
        int status = 0;
        st_dgesv((int)d, (int)d, A, sol, rhs, &status);
        if (status != -2)
        {
            gv->u[0] = sol[0];
            st_dgemm("T", "N", n, m, n, 1, gv->Z, sol + 1, gv->u + 1);
            if (weight)
                for (size_t j = 0; j < m; j++)
                    gv->u[j + 1] *= weight[j] / scale;
            for (size_t j = 0; j <= m; j++)
                if (!R_FINITE(gv->u[j]))
                    status = -2;
        }
        free(A);
        free(G);
        free(rhs);
        free(sol);
        if (status == -2)
        {
            gs_info(gv, -2,
                    "Singular equations for a zero residual variance fit");
            *retval = -2;
            return -1;
        }
    }
    *retval = 0;
    return 1;
}

void gs_mmet_coeff_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                      int *retval)
{

    gs_UNSIGNED i, j, idx_A;
    FILE *outfp;

    if (NULL == gv->y)
    {
        gs_info(
            gv, -2,
            "No performance data available, run 'read.performance.data' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available, run 'build.Z' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->ladi)
    {
        gs_info(gv, -2, "No shrinkage factor defined");
        *retval = -2;
        return;
    }

    gv->dimA = 1 + gv->NoInd;
    gv->A = realloc(gv->A, gv->dimA * gv->dimA * sizeof(double));
    if (NULL == gv->A)
        ERR_M;

    for (i = 0; i < gv->dimA; i++)
        for (j = 0; j < gv->dimA; j++)
            gv->A[i * gv->dimA + j] = 0;

    gv->A[0] = (double)gv->NoInd;

    for (j = 0; j < gv->NoInd; j++)
    {
        gv->A[0 * gv->dimA + (1 + j)] = 1;
        gv->A[(j + 1) * gv->dimA + 0] = 1;
    }

    if (!R_FINITE(gv->ladi[0]) && gv->ladi[0] > 0)
    {
        for (i = 0; i < gv->dimA; i++)
            for (j = 0; j < gv->dimA; j++)
                gv->A[i * gv->dimA + j] =
                    (i == j) ? (i == 0 ? gv->NoInd : 1) : 0;
    }
    else
    {

        gs_build_GI_01(gv, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        for (i = 0; i < gv->NoInd; i++)
            for (j = 0; j < gv->NoInd; j++)
                gv->A[(1 + i) * gv->dimA + (1 + j)] =
                    gv->GI[i * gv->NoInd + j] * gv->ladi[0];

        for (i = 0; i < gv->NoInd; i++)
            gv->A[(1 + i) * gv->dimA + (1 + i)] += 1;
    }

    if (1 == *auxfiles)
    {
        if (NULL == (outfp = fopen(*out_filename, "w")))
            ERR_F;

        for (i = 0; i < gv->dimA; i++)
        {
            for (j = 0; j < gv->dimA; j++)
            {
                idx_A = i * gv->dimA + j;
                fprintf(outfp, "%f ", gv->A[idx_A]);
            }
            fprintf(outfp, "\n");
        }

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_mmet_coeff_01_GV(char **out_filename, int *auxfiles, int *retval,
                         char **set_name)
{
    gs_mmet_coeff_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_mmet_coeff_02(gs_varset_type *gv, char **out_filename, int *auxfiles,
                      int *retval)
{

    gs_UNSIGNED i, j, idx_A;
    FILE *outfp;

    if (NULL == gv->y)
    {
        gs_info(
            gv, -2,
            "No performance data available, run 'read.performance.data' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available, run 'build.Z' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->ladi)
    {
        gs_info(gv, -2, "No shrinkage factor defined");
        *retval = -2;
        return;
    }

    gv->dimA = 1 + gv->NoInd;
    gv->A = realloc(gv->A, gv->dimA * gv->dimA * sizeof(double));
    if (NULL == gv->A)
        ERR_M;

    for (i = 0; i < gv->dimA; i++)
        for (j = 0; j < gv->dimA; j++)
            gv->A[i * gv->dimA + j] = 0;

    gv->A[0] = (double)gv->NoInd;

    for (j = 0; j < gv->NoInd; j++)
    {
        gv->A[0 * gv->dimA + (1 + j)] = 1;
        gv->A[(j + 1) * gv->dimA + 0] = 1;
    }

    if (!R_FINITE(gv->ladi[0]) && gv->ladi[0] > 0)
    {
        for (i = 0; i < gv->dimA; i++)
            for (j = 0; j < gv->dimA; j++)
                gv->A[i * gv->dimA + j] =
                    (i == j) ? (i == 0 ? gv->NoInd : 1) : 0;
    }
    else
    {

        for (i = 0; i < gv->NoInd; i++)
            for (j = 0; j < gv->NoInd; j++)
                gv->A[(1 + i) * gv->dimA + (1 + j)] =
                    gv->ladi[0] * gv->GI[i * gv->NoInd + j];

        for (i = 0; i < gv->NoInd; i++)
            gv->A[(1 + i) * gv->dimA + (1 + i)] += 1;
    }

    if (1 == *auxfiles)
    {
        if (NULL == (outfp = fopen(*out_filename, "w")))
            ERR_F;

        for (i = 0; i < gv->dimA; i++)
        {
            for (j = 0; j < gv->dimA; j++)
            {
                idx_A = i * gv->dimA + j;
                fprintf(outfp, "%f ", gv->A[idx_A]);
            }
            fprintf(outfp, "\n");
        }

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_mmet_coeff_02_GV(char **out_filename, int *auxfiles, int *retval,
                         char **set_name)
{
    gs_mmet_coeff_02(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_mmet_rhs_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                    int *retval)
{
    gs_UNSIGNED i, j, dim_b;
    FILE *outfp;

    if (NULL == gv->y)
    {
        gs_info(
            gv, -2,
            "No performance data available, run 'read.performance.data' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available, run 'build.Z' first");
        *retval = -2;
        return;
    }

    dim_b = 1 + gv->NoInd;
    gv->b = realloc(gv->b, dim_b * sizeof(double));
    if (NULL == gv->b)
        ERR_M;

    gv->b[0] = 0;
    for (j = 0; j < gv->NoInd; j++)
        gv->b[0] += gv->y[j];

    for (j = 0; j < gv->NoInd; j++)
        gv->b[1 + j] = (gv->ladi && !R_FINITE(gv->ladi[0]) && gv->ladi[0] > 0)
                           ? 0
                           : gv->y[j];

    if (1 == *auxfiles)
    {

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        for (i = 0; i < dim_b; i++)
            fprintf(outfp, "%f \n", gv->b[i]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_mmet_rhs_01_GV(char **out_filename, int *auxfiles, int *retval,
                       char **set_name)
{
    gs_mmet_rhs_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_mmet_solve_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                      int *retval)
{
    gs_UNSIGNED i;
    FILE *outfp;

    int Srv = 0;
    int *Sretval = &Srv;

    if (NULL == gv->A)
    {
        gs_info(gv, -2, "No coefficient matrix available");
        *retval = -2;
        return;
    }

    if (NULL == gv->b)
    {
        gs_info(gv, -2, "No right hand side available");
        *retval = -2;
        return;
    }

    int boundary = gs_rr_boundary_01(gv, NULL, retval);
    if (boundary < 0)
        return;
    if (!boundary)
    {

        gv->a = realloc(gv->a, (1 + gv->NoInd) * sizeof(double));
        if (NULL == gv->a)
            ERR_M;

        st_dgesv(gv->dimA, gv->dimA, gv->A, gv->a, gv->b, Sretval);
        if (*Sretval == -2)
        {
            *retval = -2;
            return;
        }

        gv->dimu = 1 + gv->NoMar;
        gv->u = realloc(gv->u, gv->dimu * sizeof(double));
        if (NULL == gv->u)
            ERR_M;

        double *tmp1;
        tmp1 = malloc(gv->NoInd * gv->NoMar * sizeof(double));
        if (NULL == tmp1)
            ERR_M;

        st_dgemm("T", "N", gv->NoInd, gv->NoMar, gv->NoInd, gv->NoInd, gv->Z,
                 gv->GI, tmp1);

        st_dgemm("N", "N", gv->NoMar, gv->NoInd, gv->NoInd, 1, tmp1, 1 + gv->a,
                 1 + gv->u);

        gv->u[0] = gv->a[0];

        free(tmp1);
    }
    free(gv->A);
    gv->A = NULL;

    if (1 == *auxfiles)
    {

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "effect\n");
        for (i = 0; i < gv->dimu; i++)
            fprintf(outfp, "%s.%i %f \n", gv->EffNme[i], gv->EffAll[i],
                    gv->u[i]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_mmet_solve_01_GV(char **out_filename, int *auxfiles, int *retval,
                         char **set_name)
{
    gs_mmet_solve_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_lambda_rrblup_01(gs_varset_type *gv, int *maxiter, double *precision,
                         double *hsq, char **out_filename, int *auxfiles,
                         int *retval);

void gs_lambda_rmla_02(gs_varset_type *gv, double *alpha, int *maxiter,
                       double *precision, double *hsq, char **out_filename,
                       int *auxfiles, int *retval);

void gs_lambda_rmlc_01(gs_varset_type *gv, int *maxiter, double *precision,
                       double *hsq, char **out_filename, int *auxfiles,
                       int *retval);

void gs_mme_coeff_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                     int *retval);

void gs_lambda_rmla_01(gs_varset_type *gv, double *alpha, int *maxiter,
                       double *precision, double *hsq, char **out_filename,
                       int *auxfiles, int *retval);

void gs_mme_rhs_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                   int *retval);

void gs_mme_solve_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                     int *retval);

void gs_lambda_rml_02(gs_varset_type *gv, double *alpha, int *maxiter,
                      double *precision, double *hsq, char **out_filename,
                      int *auxfiles, int *retval);

void gs_lambda_rmlv_01(gs_varset_type *gv, int *maxiter, double *precision,
                       double *hsq, char **out_filename, int *auxfiles,
                       int *retval);

void gs_vc_rrblup_01(gs_varset_type *gv, int *maxiter, double *precision,
                     double *hsq, char **out_filename, int *auxfiles,
                     int *retval);

void gs_w_aov_01(gs_varset_type *gv, double *alpha, char **out_filename,
                 int *auxfiles, int *retval);

void gs_w_I_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
               int *retval);

void gs_single_marker_aov_01(gs_varset_type *gv, double *alpha,
                             char **out_filename, int *auxfiles, int *retval);

void gs_esteff_rr_01(gs_varset_type *gv, char **methode, int *maxiter,
                     double *precision, double *hsq, double *alpha,
                     char **out_filename, int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;
    int fit_status = 0;

    if (strcmp(*methode, "BLUP") == 0)
    {

        gs_build_Z_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_vc_rrblup_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
        fit_status = *Sretval;

        gs_lambda_rrblup_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mmet_coeff_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mmet_rhs_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mmet_solve_01(gv, out_filename, auxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }

    else if (strcmp(*methode, "RMLA") == 0)
    {
        gs_build_Z_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_vc_rrblup_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
        fit_status = *Sretval;

        gs_lambda_rrblup_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_single_marker_aov_01(gv, alpha, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        if (gv->sigma[0] == 0)
        {
            if (gs_rr_boundary_01(gv, gv->SMvar, Sretval) < 0)
            {
                *retval = -2;
                return;
            }
            if (1 == *auxfiles)
            {
                FILE *outfp = fopen(*out_filename, "w");
                if (!outfp)
                    ERR_F;
                fprintf(outfp, "effect\n");
                for (gs_UNSIGNED i = 0; i < gv->dimu; i++)
                    fprintf(outfp, "%s.%i %f \n", gv->EffNme[i], gv->EffAll[i],
                            gv->u[i]);
                fclose(outfp);
            }
            *retval = fit_status;
            return;
        }

        gs_lambda_rmla_02(gv, alpha, maxiter, precision, hsq, nixp, neinp,
                          Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_coeff_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_rhs_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }

    else if (strcmp(*methode, "rmlc") == 0)
    {
        gs_build_Z_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_lambda_rmlc_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_coeff_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_rhs_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }

    else if (strcmp(*methode, "rmla") == 0)
    {
        gs_build_Z_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_lambda_rmla_01(gv, alpha, maxiter, precision, hsq, nixp, neinp,
                          Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_coeff_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_rhs_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }

    else if (strcmp(*methode, "rmlv") == 0)
    {
        gs_build_Z_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_lambda_rmlv_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_coeff_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_rhs_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }

    else
    {
        gs_info(gv, -2, "Unknown estimation method");
        *retval = -2;
        return;
    }

    *retval = fit_status;
    return;
}

void gs_esteff_rr_01_GV(char **methode, int *maxiter, double *precision,
                        double *hsq, double *alpha, char **out_filename,
                        int *auxfiles, int *retval, char **set_name)
{
    gs_esteff_rr_01(gs_fdta(*set_name), methode, maxiter, precision, hsq, alpha,
                    out_filename, auxfiles, retval);
}

void gs_cross_validation_02(gs_varset_type *gv, char **methode,
                            int *n_estimation, int *n_runs, double *hsq,
                            double *alpha, int *maxiter, double *precision,
                            char **out_filename_r, char **out_filename_u,
                            int *auxfiles, int *retval)
{
    int ntr = 1;
    int prev_info_level;
    int no_info[1] = {-2};
    int CVNoRuns;
    int parallel_error = 0;
    int fit_status = 0;

    double *CVcor = NULL;
    gs_FLOAT *CVest = NULL;

    size_t no_mar;
    size_t cv_cols;
    size_t cv_count;

    *retval = -2;

    if (gv->NoInd > (gs_UNSIGNED)INT_MAX)
    {
        gs_info(gv, -2, "More individuals than INT_MAX");
        return;
    }

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2, "No marker data available");
        return;
    }

    if (NULL == gv->y)
    {
        gs_info(gv, -2, "No performance data available");
        return;
    }

    if ((*n_estimation < 1) || (*n_estimation >= (int)gv->NoInd))
    {
        gs_info(gv, -2, "Invalid size of training set");
        return;
    }

    if (*n_runs < 1)
    {
        gs_info(gv, -2, "Number of cross-validation runs must be positive");
        return;
    }

    if (strcmp(*methode, "BLUP") == 0)
    {
    }
    else if (strcmp(*methode, "RMLA") == 0)
    {
    }
    else if (strcmp(*methode, "rmlc") == 0)
    {
    }
    else if (strcmp(*methode, "rmlv") == 0)
    {
    }
    else if (strcmp(*methode, "rmla") == 0)
    {
    }
    else
    {
        gs_info(gv, -2, "Unknown estimation method");
        return;
    }

    CVNoRuns = *n_runs;

    no_mar = (size_t)gv->NoMar;
    if (((gs_UNSIGNED)no_mar != gv->NoMar) || (no_mar == ((size_t)-1)))
    {
        gs_info(gv, -2, "Cross-validation result size is too large");
        return;
    }

    cv_cols = no_mar + 1;

    if ((size_t)CVNoRuns > ((size_t)-1) / cv_cols)
    {
        gs_info(gv, -2, "Cross-validation result size is too large");
        return;
    }

    cv_count = (size_t)CVNoRuns * cv_cols;

    if (((size_t)CVNoRuns > ((size_t)-1) / sizeof(*CVcor)) ||
        (cv_count > ((size_t)-1) / sizeof(*CVest)))
    {
        gs_info(gv, -2, "Cross-validation result size is too large");
        return;
    }

    CVcor = malloc((size_t)CVNoRuns * sizeof(*CVcor));
    CVest = malloc(cv_count * sizeof(*CVest));

    if ((NULL == CVcor) || (NULL == CVest))
    {
        free(CVcor);
        free(CVest);
        gs_info(gv, -2, "Memory allocation error");
        return;
    }

    for (int c = 0; c < CVNoRuns; c++)
        CVcor[c] = 0;

    ntr = 1;
    (void)ntr;

    st_get_info_level(&prev_info_level);
    st_set_info_level(no_info);

#ifdef MT
#pragma omp parallel num_threads(ntr) if (ntr > 1)                             \
    reduction(| : parallel_error) reduction(max : fit_status)
#endif
    {
        int Srv;
        int *Sretval = &Srv;

        gs_varset_type ES;
        gs_varset_type *es = &ES;

        gs_varset_type VS;
        gs_varset_type *vs = &VS;

        gs_UNSIGNED i;
        double sy, syy, sx, sxx, sxy, k;

        gs_init(es);
        gs_init(vs);

        es->info_level = -2;
        vs->info_level = -2;

#ifdef MT
        omp_set_schedule(omp_sched_static, 0);
#pragma omp for schedule(runtime)
#endif
        for (int r = 0; r < CVNoRuns; r++)
        {
            int success_rep = 0;
            int no_unsucc_rep = 0;

            while ((0 == success_rep) && (no_unsucc_rep < 1000))
            {
                *Sretval = 0;

                gs_build_esvs_01(gv, es, vs, n_estimation, nixp, nixp, nixp,
                                 nixp, neinp, Sretval);

                if (-2 == *Sretval)
                {
                    parallel_error = 1;
                    break;
                }

                *Sretval = 0;

                gs_esteff_rr_01(es, methode, maxiter, precision, hsq, alpha,
                                nixp, neinp, Sretval);

                if (-2 == *Sretval)
                {
                    no_unsucc_rep++;
                    continue;
                }

                if (*Sretval > fit_status)
                    fit_status = *Sretval;
                success_rep = 1;
            }

            if (parallel_error)
                continue;

            if (0 == success_rep)
            {
                parallel_error = 1;
                continue;
            }

            if ((NULL == es->u) || ((size_t)es->dimu != cv_cols))
            {
                parallel_error = 1;
                continue;
            }

            for (i = 0; i < es->dimu; i++)
                CVest[(size_t)r * cv_cols + (size_t)i] = es->u[i];

            *Sretval = 0;

            gs_estimate_bv_01(es, vs, nixp, neinp, Sretval);

            if ((-2 == *Sretval) || (NULL == vs->est))
            {
                parallel_error = 1;
                continue;
            }

            sy = 0;
            syy = 0;
            sx = 0;
            sxx = 0;
            sxy = 0;
            k = (double)vs->NoInd;

            for (i = 0; i < vs->NoInd; i++)
            {
                sx += vs->est[i];
                sxx += vs->est[i] * vs->est[i];
                sy += vs->y[i];
                syy += vs->y[i] * vs->y[i];
                sxy += vs->est[i] * vs->y[i];
            }

            CVcor[r] = (k * sxy - sx * sy) /
                       sqrt((k * sxx - sx * sx) * (k * syy - sy * sy));
        }

        gs_reset(es);
        gs_reset(vs);
    }

    if (parallel_error)
    {
        st_set_info_level(&prev_info_level);
        free(CVcor);
        free(CVest);
        gs_info(gv, -2, "Cross validation failed in parallel calculation");
        return;
    }

    {
        int Srv = 0;
        int *Sretval = &Srv;
        char *meth_ = "BLUP";
        char **meth = &meth_;

        gs_esteff_rr_01(gv, meth, maxiter, precision, hsq, alpha, nixp, neinp,
                        Sretval);

        if (-2 == *Sretval)
        {
            st_set_info_level(&prev_info_level);
            free(CVcor);
            free(CVest);
            *retval = -2;
            return;
        }
        if (*Sretval > fit_status)
            fit_status = *Sretval;
    }

    if ((NULL == gv->u) || ((size_t)gv->dimu != cv_cols))
    {
        st_set_info_level(&prev_info_level);
        free(CVcor);
        free(CVest);
        gs_info(gv, -2, "Inconsistent number of effects in cross validation");
        return;
    }

    for (gs_UNSIGNED i = 0; i < gv->dimu; i++)
    {
        gs_FLOAT sum = 0;

        for (int r = 0; r < CVNoRuns; r++)
            sum += CVest[(size_t)r * cv_cols + (size_t)i];

        gv->u[i] = (double)(sum / (double)CVNoRuns);
    }

    st_set_info_level(&prev_info_level);

    if (1 == *auxfiles)
    {
        FILE *fp;

        fp = fopen(*out_filename_r, "w");
        if (NULL == fp)
        {
            free(CVcor);
            free(CVest);
            gs_info(gv, -2, "File open error");
            return;
        }

        fprintf(fp, "cor \n");
        for (int c = 0; c < CVNoRuns; c++)
            fprintf(fp, "%i %f\n", c, CVcor[c]);
        fclose(fp);

        fp = fopen(*out_filename_u, "w");
        if (NULL == fp)
        {
            free(CVcor);
            free(CVest);
            gs_info(gv, -2, "File open error");
            return;
        }

        fprintf(fp, "effect\n");
        for (gs_UNSIGNED i = 0; i < gv->dimu; i++)
            fprintf(fp, "%s.%i %f \n", gv->EffNme[i], gv->EffAll[i], gv->u[i]);
        fclose(fp);
    }

    free(CVcor);
    free(CVest);

    *retval = fit_status;
}

void gs_cross_validation_02_GV(char **scheme, int *n_estimation, int *n_runs,
                               double *hsq, double *alpha, int *maxiter,
                               double *precision, char **out_filename_r,
                               char **out_filename_u, int *auxfiles,
                               int *retval, char **set_name)
{
    gs_cross_validation_02(gs_fdta(*set_name), scheme, n_estimation, n_runs,
                           hsq, alpha, maxiter, precision, out_filename_r,
                           out_filename_u, auxfiles, retval);
}

gs_FLOAT gs_ld_rsquare(gs_varset_type *gv, gs_UNSIGNED ii, gs_UNSIGNED jj)
{
    gs_UNSIGNED i, k, l;

    gs_UNSIGNED NoAllii = gv->NoAll[ii];
    gs_UNSIGNED NoAlljj = gv->NoAll[jj];
    GS_ALLELE_TYPE *Allii = gv->All[ii];
    GS_ALLELE_TYPE *Alljj = gv->All[jj];
    gs_UNSIGNED *AllCntii = gv->AllCnt[ii];
    gs_UNSIGNED *AllCntjj = gv->AllCnt[jj];

    if ((0 < NoAllii) && (GS_ALLELE_MISSING == Allii[NoAllii - 1]))
        NoAllii--;
    if ((0 < NoAlljj) && (GS_ALLELE_MISSING == Alljj[NoAlljj - 1]))
        NoAlljj--;

    st_allele_marker_view vi = st_allele_marker(gv, ii),
                          vj = st_allele_marker(gv, jj);
    gs_UNSIGNED CNT;
    double N = 2 * (double)gv->NoInd;
    double r2 = 0;

    for (k = 0; k < NoAllii; k++)
        for (l = 0; l < NoAlljj; l++)
        {
            CNT = 0;
            for (i = 0; i < gv->NoInd; i++)
            {
                st_allele_pair pi = st_allele_read_pair(&vi, i),
                               pj = st_allele_read_pair(&vj, i);
                if ((pi.a == Allii[k]) && (pj.a == Alljj[l]))
                    CNT++;
                if ((pi.b == Allii[k]) && (pj.b == Alljj[l]))
                    CNT++;
            }
            double fx = (double)AllCntii[k] / N;
            double fy = (double)AllCntjj[l] / N;
            double Dxy = (double)CNT / N - fx * fy;
            double r2xy = Dxy * Dxy / (fx * (1 - fx) * fy * (1 - fy));
            r2 += r2xy * fx * fy;
        }
    return (r2);
}

gs_FLOAT gs_ld_dprime(gs_varset_type *gv, gs_UNSIGNED ii, gs_UNSIGNED jj)
{
    gs_UNSIGNED i, k, l;

    gs_UNSIGNED NoAllii = gv->NoAll[ii];
    gs_UNSIGNED NoAlljj = gv->NoAll[jj];
    GS_ALLELE_TYPE *Allii = gv->All[ii];
    GS_ALLELE_TYPE *Alljj = gv->All[jj];
    gs_UNSIGNED *AllCntii = gv->AllCnt[ii];
    gs_UNSIGNED *AllCntjj = gv->AllCnt[jj];

    if ((0 < NoAllii) && (GS_ALLELE_MISSING == Allii[NoAllii - 1]))
        NoAllii--;
    if ((0 < NoAlljj) && (GS_ALLELE_MISSING == Alljj[NoAlljj - 1]))
        NoAlljj--;

    st_allele_marker_view vi = st_allele_marker(gv, ii),
                          vj = st_allele_marker(gv, jj);
    gs_UNSIGNED CNT;
    double N = 2 * (double)gv->NoInd;
    double Dp = 0;
    double fx, fy, Dxy, c1, c2, Dmax, Dpxy;

    for (k = 0; k < NoAllii; k++)
        for (l = 0; l < NoAlljj; l++)
        {
            CNT = 0;
            for (i = 0; i < gv->NoInd; i++)
            {
                st_allele_pair pi = st_allele_read_pair(&vi, i),
                               pj = st_allele_read_pair(&vj, i);
                if ((pi.a == Allii[k]) && (pj.a == Alljj[l]))
                    CNT++;
                if ((pi.b == Allii[k]) && (pj.b == Alljj[l]))
                    CNT++;
            }
            fx = (double)(AllCntii[k]) / N;
            fy = (double)(AllCntjj[l]) / N;
            Dxy = (double)CNT / N - fx * fy;
            if (Dxy < 0)
            {
                c1 = fx * fy;
                c2 = (1 - fx) * (1 - fy);
            }
            else
            {
                c1 = fx * (1 - fy);
                c2 = (1 - fx) * fy;
            }
            if (c1 < c2)
                Dmax = c1;
            else
                Dmax = c2;
            Dpxy = fabs(Dxy / Dmax);
            Dp += Dpxy * fx * fy;
        }
    return (Dp);
}

static int gs_ld_unphased_marker_alleles(const gs_varset_type *gv,
                                         gs_UNSIGNED marker,
                                         GS_ALLELE_TYPE allele[2])
{
    gs_UNSIGNED n = 0;

    if ((NULL == gv) || (NULL == gv->NoAll) || (NULL == gv->All) ||
        (marker >= gv->NoMar) || (NULL == gv->All[marker]))
        return 0;

    for (gs_UNSIGNED k = 0; k < gv->NoAll[marker]; k++)
    {
        GS_ALLELE_TYPE a = gv->All[marker][k];

        if (GS_ALLELE_MISSING == a)
            continue;

        if (2 <= n)
            return 0;

        allele[n++] = a;
    }

    return (2 == n);
}

static int gs_ld_unphased_data_are_biallelic(const gs_varset_type *gv)
{
    GS_ALLELE_TYPE allele[2];

    if ((NULL == gv) || (0 == gv->NoMar))
        return 0;

    for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
        if (!gs_ld_unphased_marker_alleles(gv, m, allele))
            return 0;

    return 1;
}

static int gs_ld_unphased_dosage(GS_ALLELE_TYPE a1, GS_ALLELE_TYPE a2,
                                 const GS_ALLELE_TYPE allele[2])
{
    int d = 0;

    if ((GS_ALLELE_MISSING == a1) || (GS_ALLELE_MISSING == a2))
        return -1;

    if ((a1 != allele[0]) && (a1 != allele[1]))
        return -2;
    if ((a2 != allele[0]) && (a2 != allele[1]))
        return -2;

    if (a1 == allele[1])
        d++;
    if (a2 == allele[1])
        d++;

    return d;
}

static int gs_ld_unphased_table(const gs_varset_type *gv, gs_UNSIGNED ii,
                                gs_UNSIGNED jj, gs_UNSIGNED table[9])
{
    GS_ALLELE_TYPE ai[2];
    GS_ALLELE_TYPE aj[2];

    if (!gs_ld_unphased_marker_alleles(gv, ii, ai) ||
        !gs_ld_unphased_marker_alleles(gv, jj, aj))
        return 0;

    for (int k = 0; k < 9; k++)
        table[k] = 0;

    st_allele_marker_view vi = st_allele_marker(gv, ii),
                          vj = st_allele_marker(gv, jj);
    for (gs_UNSIGNED ind = 0; ind < gv->NoInd; ind++)
    {
        st_allele_pair pi = st_allele_read_pair(&vi, ind),
                       pj = st_allele_read_pair(&vj, ind);
        int di = gs_ld_unphased_dosage(pi.a, pi.b, ai);
        int dj = gs_ld_unphased_dosage(pj.a, pj.b, aj);

        if ((-1 == di) || (-1 == dj))
            continue;

        if ((di < 0) || (dj < 0))
            return 0;

        table[3 * di + dj]++;
    }

    return 1;
}

static double gs_ld_unphased_loglik(const gs_UNSIGNED table[9],
                                    const double h[4])
{
    double p[9];
    double ll = 0.0;

    p[0] = h[0] * h[0];
    p[1] = 2.0 * h[0] * h[1];
    p[2] = h[1] * h[1];
    p[3] = 2.0 * h[0] * h[2];
    p[4] = 2.0 * (h[0] * h[3] + h[1] * h[2]);
    p[5] = 2.0 * h[1] * h[3];
    p[6] = h[2] * h[2];
    p[7] = 2.0 * h[2] * h[3];
    p[8] = h[3] * h[3];

    for (int k = 0; k < 9; k++)
        if (0 < table[k])
        {
            if (!(p[k] > 0.0) || !isfinite(p[k]))
                return -HUGE_VAL;
            ll += (double)table[k] * log(p[k]);
        }

    return ll;
}

static int gs_ld_unphased_em(const gs_UNSIGNED table[9],
                             double initial_coupling, double h[4],
                             double *loglik)
{
    const int max_iter = 10000;
    const double tol = 1.0e-13;
    double fixed[4];
    double total_ind = 0.0;
    double double_het = (double)table[4];
    double total_hap;
    double p = initial_coupling;

    for (int k = 0; k < 9; k++)
        total_ind += (double)table[k];

    if (!(total_ind > 0.0))
        return 0;

    total_hap = 2.0 * total_ind;

    fixed[0] = 2.0 * (double)table[0] + (double)table[1] + (double)table[3];
    fixed[1] = 2.0 * (double)table[2] + (double)table[1] + (double)table[5];
    fixed[2] = 2.0 * (double)table[6] + (double)table[3] + (double)table[7];
    fixed[3] = 2.0 * (double)table[8] + (double)table[5] + (double)table[7];

    if (p < 0.0)
        p = 0.0;
    if (p > 1.0)
        p = 1.0;

    h[0] = (fixed[0] + double_het * p) / total_hap;
    h[1] = (fixed[1] + double_het * (1.0 - p)) / total_hap;
    h[2] = (fixed[2] + double_het * (1.0 - p)) / total_hap;
    h[3] = (fixed[3] + double_het * p) / total_hap;

    if (0.0 == double_het)
    {
        *loglik = gs_ld_unphased_loglik(table, h);
        return isfinite(*loglik);
    }

    for (int iter = 0; iter < max_iter; iter++)
    {
        double coupling = h[0] * h[3];
        double repulsion = h[1] * h[2];
        double denom = coupling + repulsion;
        double next[4];
        double q;
        double delta = 0.0;

        if (!(denom > 0.0) || !isfinite(denom))
            return 0;

        q = coupling / denom;

        next[0] = (fixed[0] + double_het * q) / total_hap;
        next[1] = (fixed[1] + double_het * (1.0 - q)) / total_hap;
        next[2] = (fixed[2] + double_het * (1.0 - q)) / total_hap;
        next[3] = (fixed[3] + double_het * q) / total_hap;

        for (int k = 0; k < 4; k++)
        {
            double d = fabs(next[k] - h[k]);
            if (d > delta)
                delta = d;
            h[k] = next[k];
        }

        if (delta <= tol)
        {
            *loglik = gs_ld_unphased_loglik(table, h);
            return isfinite(*loglik);
        }
    }

    return 0;
}

static int gs_ld_unphased_haplotypes(const gs_UNSIGNED table[9], double h[4],
                                     int *ambiguous_sign)
{

    static const double starts[] = {0.0,  0.01, 0.10, 0.25, 0.50,
                                    0.75, 0.90, 0.99, 1.0};
    double best_ll = -HUGE_VAL;
    double best_h[4] = {0.0, 0.0, 0.0, 0.0};
    double candidate_ll[sizeof(starts) / sizeof(starts[0])];
    double candidate_d[sizeof(starts) / sizeof(starts[0])];
    size_t candidate_count = 0;
    int found = 0;

    if (ambiguous_sign)
        *ambiguous_sign = 0;

    for (size_t s = 0; s < sizeof(starts) / sizeof(starts[0]); s++)
    {
        double candidate[4];
        double ll;

        if (!gs_ld_unphased_em(table, starts[s], candidate, &ll))
            continue;

        if (ambiguous_sign)
        {
            candidate_ll[candidate_count] = ll;
            candidate_d[candidate_count] =
                candidate[0] * candidate[3] - candidate[1] * candidate[2];
            candidate_count++;
        }

        if ((!found) || (ll > best_ll))
        {
            best_ll = ll;
            for (int k = 0; k < 4; k++)
                best_h[k] = candidate[k];
            found = 1;
        }
    }

    if (!found)
        return 0;

    for (int k = 0; k < 4; k++)
        h[k] = best_h[k];

    if (ambiguous_sign)
    {
        double best_d = h[0] * h[3] - h[1] * h[2];
        double tolerance = 1.0e-10 * (1.0 + fabs(best_ll));
        for (size_t s = 0; s < candidate_count; s++)
            if (fabs(candidate_ll[s] - best_ll) <= tolerance &&
                fabs(best_d) > 1.0e-12 && fabs(candidate_d[s]) > 1.0e-12 &&
                ((best_d < 0.0) != (candidate_d[s] < 0.0)))
                *ambiguous_sign = 1;
    }

    return 1;
}

static int gs_ld_unphased_measures(const gs_varset_type *gv, gs_UNSIGNED ii,
                                   gs_UNSIGNED jj, double *r2, double *dprime,
                                   double *signed_r)
{
    gs_UNSIGNED table[9];
    double h[4];
    double pa0, pa1, pb0, pb1;
    double D, denom, dmax;
    int ambiguous_sign = 0;

    *r2 = 0.0;
    *dprime = 0.0;
    if (signed_r)
        *signed_r = NA_REAL;

    if (!gs_ld_unphased_table(gv, ii, jj, table))
        return 0;
    if (!gs_ld_unphased_haplotypes(table, h, signed_r ? &ambiguous_sign : NULL))
        return 0;

    pa0 = h[0] + h[1];
    pa1 = h[2] + h[3];
    pb0 = h[0] + h[2];
    pb1 = h[1] + h[3];

    denom = pa0 * pa1 * pb0 * pb1;
    if (!(denom > 0.0))
        return 2;

    D = h[0] - pa0 * pb0;
    *r2 = D * D / denom;
    if (signed_r && !ambiguous_sign)
    {
        *signed_r = D / sqrt(denom);
        if (*signed_r > 1.0 && *signed_r < 1.0 + 1.0e-12)
            *signed_r = 1.0;
        if (*signed_r < -1.0 && *signed_r > -1.0 - 1.0e-12)
            *signed_r = -1.0;
    }

    if (D < 0.0)
        dmax = fmin(pa0 * pb0, pa1 * pb1);
    else
        dmax = fmin(pa0 * pb1, pa1 * pb0);

    if (dmax > 0.0)
        *dprime = fabs(D) / dmax;

    if (*r2 < 0.0)
        *r2 = 0.0;
    if (*r2 > 1.0 && *r2 < 1.0 + 1.0e-12)
        *r2 = 1.0;
    if (*dprime < 0.0)
        *dprime = 0.0;
    if (*dprime > 1.0 && *dprime < 1.0 + 1.0e-12)
        *dprime = 1.0;

    return 1;
}

static gs_FLOAT gs_ld_rsquare_unphased(gs_varset_type *gv, gs_UNSIGNED ii,
                                       gs_UNSIGNED jj)
{
    double r2, dprime;

    if (!gs_ld_unphased_measures(gv, ii, jj, &r2, &dprime, NULL))
        return 0.0;

    return (gs_FLOAT)r2;
}

static gs_FLOAT gs_ld_dprime_unphased(gs_varset_type *gv, gs_UNSIGNED ii,
                                      gs_UNSIGNED jj)
{
    double r2, dprime;

    if (!gs_ld_unphased_measures(gv, ii, jj, &r2, &dprime, NULL))
        return 0.0;

    return (gs_FLOAT)dprime;
}

static gs_FLOAT gs_ld_r_unphased(gs_varset_type *gv, gs_UNSIGNED ii,
                                 gs_UNSIGNED jj)
{
    double r2, dprime, signed_r;

    if (gs_ld_unphased_measures(gv, ii, jj, &r2, &dprime, &signed_r) != 1)
        return NA_REAL;

    return signed_r;
}

void gs_calc_ld_01(gs_varset_type *gv, char **ld_measure, char **out_filename,
                   int *auxfiles, int *retval)

{
    gs_UNSIGNED c, i, j, idx;
    gs_UNSIGNED ii, jj;
    gs_mappoint_type **p;
    int Srv = 0;
    int *Sretval = &Srv;

    gs_chrom_stats_01(gv, Sout_filename, Sauxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_FLOAT (*gs_ld)(gs_varset_type *gv, gs_UNSIGNED i, gs_UNSIGNED j);

    if (0 == strcmp(*ld_measure, "r2"))
    {
        gs_ld = gs_ld_rsquare;
    }

    else if (0 == strcmp(*ld_measure, "r-estimate-phases"))
    {
        if (gs_ld_unphased_data_are_biallelic(gv))
            gs_ld = gs_ld_r_unphased;
        else
        {
            gs_info(gv, -2,
                    "Phase-estimating LD requires exactly two biological "
                    "alleles at every marker");
            *retval = -2;
            return;
        }
    }

    else if (0 == strcmp(*ld_measure, "r2-estimate-phases"))
    {
        if (gs_ld_unphased_data_are_biallelic(gv))
            gs_ld = gs_ld_rsquare_unphased;
        else
        {
            gs_info(gv, -2,
                    "Phase-estimating LD requires exactly two biological "
                    "alleles at every marker");
            *retval = -2;
            return;
        }
    }

    else if (0 == strcmp(*ld_measure, "Dp"))
    {
        gs_ld = gs_ld_dprime;
    }

    else if (0 == strcmp(*ld_measure, "Dp-estimate-phases"))
    {
        if (gs_ld_unphased_data_are_biallelic(gv))
            gs_ld = gs_ld_dprime_unphased;
        else
        {
            gs_info(gv, -2,
                    "Phase-estimating LD requires exactly two biological "
                    "alleles at every marker");
            *retval = -2;
            return;
        }
    }

    else
    {
        gs_info(gv, -2, "LD measure not (yet) implemented");
        *retval = -2;
        return;
    }

    for (c = 0; c < gv->NChrom; c++)
    {
        gv->chrom[c].ld =
            realloc(gv->chrom[c].ld,
                    gv->chrom[c].NLocC * gv->chrom[c].NLocC * sizeof(gs_FLOAT));
    }

    for (c = 0; c < gv->NChrom; c++)
    {
        p = gv->chrom[c].first_loc;
#ifdef MT
#pragma omp parallel for private(i, j, ii, jj, idx)
#endif
        for (i = 0; i < (gv->chrom[c].NLocC - 1); i++)
            for (j = i + 1; j < gv->chrom[c].NLocC; j++)
            {
                ii = (p[i])->nr;
                jj = (p[j])->nr;
                idx = i * gv->chrom[c].NLocC + j;
                (gv->chrom[c].ld)[idx] = gs_ld(gv, ii, jj);
            }
    }

    if (1 == *auxfiles)
    {
        gs_UNSIGNED loc_rgt, loc_lft;

        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        for (c = 0; c < gv->NChrom; c++)
        {
            for (i = 0; i < (gv->chrom[c].NLocC - 1); i++)
                for (j = i + 1; j < gv->chrom[c].NLocC; j++)
                {

                    idx = i * gv->chrom[c].NLocC + j;

                    loc_lft = gv->chrom[c].first_loc_i + i;
                    loc_rgt = gv->chrom[c].first_loc_i + j;

                    snprintf(gv->msg, sizeof gv->msg, "%lu %lu %lu %s %s %f \n",
                             1UL + (unsigned long)c, 1UL + (unsigned long)i,
                             1UL + (unsigned long)j, gv->map[loc_lft]->name,
                             gv->map[loc_rgt]->name,
                             (double)(gv->chrom[c].ld)[idx]);

                    fprintf(outfp, "%s", gv->msg);
                }
        }

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_calc_ld_01_GV(char **ld_measure, char **out_filename, int *auxfiles,
                      int *retval, char **set_name)
{
    gs_calc_ld_01(gs_fdta(*set_name), ld_measure, out_filename, auxfiles,
                  retval);
}

void gs_set_hblocks_01(gs_varset_type *gv, int *n_hap, int *blk,
                       char **hap_symbol, char **out_filename, int *auxfiles,
                       int *retval)

{

    if (0 == gv->NLoc)
    {
        snprintf(gv->msg, sizeof(gv->msg), "No linkage map has been loaded");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    if ((gs_UNSIGNED)(*n_hap) != gv->NLoc)
    {
        snprintf(gv->msg, sizeof(gv->msg),
                 "Incorrect number of markers: %lu - %lu",
                 (gs_UNSIGNED)(*n_hap), gv->NLoc);

        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    gs_UNSIGNED i, l, fib;

    gv->h_u = realloc(gv->h_u, gv->NLoc * sizeof(gs_UNSIGNED));
    gv->h_o = realloc(gv->h_o, gv->NLoc * sizeof(gs_UNSIGNED));
    gv->h_n = realloc(gv->h_n, gv->NLoc * sizeof(gs_UNSIGNED));
    if ((NULL == gv->h_u) || (NULL == gv->h_o) || (gv->h_n == 0))
        ERR_M;

    gv->NHap = 0;

    gv->NHap = 0;
    fib = 0;

    for (i = 1; i < gv->NLoc; i++)
    {
        if (blk[i - 1] != blk[i])
        {
            gv->h_u[gv->NHap] = fib;
            gv->h_o[gv->NHap] = i - 1;
            gv->h_n[gv->NHap] = i - fib;
            gv->NHap++;
            fib = i;
        }
    }
    gv->h_u[gv->NHap] = fib;
    gv->h_o[gv->NHap] = i - 1;
    gv->h_n[gv->NHap] = i - fib;
    gv->NHap++;

    {
        int map_rc = gs_build_haplotype_map(gv, hap_symbol);

        if (0 != map_rc)
        {
            gv->NHap = 0;
            GS_FREE(gv->Hap_S);
            GS_FREE(gv->Hap_str_S);
            gv->Hap_str_N = 0;

            if (-1 == map_rc)
                ERR_M;

            gs_info(gv, -2, "Haplotype marker name too long");
            *retval = -2;
            return;
        }
    }

    if (1 == *auxfiles)
    {

        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        for (l = 0; l < gv->NHap; l++)
        {
            char *lne = (char *)malloc(gs_MAXLNE * sizeof(char));
            snprintf(lne, gs_MAXLNE, "%i %f %s %s ", gv->Hap_S[l].chrom,
                     gv->Hap_S[l].pos, gv->Hap_S[l].name, gv->Hap_S[l].class);

            for (i = gv->h_u[l]; i <= gv->h_o[l]; i++)
            {
                if (gv->mdl[i] < gs_UNSIGNED_MAX)
                {
                    size_t a = strlen(lne);
                    size_t b = strlen(gv->MarName[gv->mdl[i]]);
                    lne = realloc(lne, 2 + a + b);
                    strcat(lne, gv->MarName[gv->mdl[i]]);
                    strcat(lne, ";");
                }
            }

            fprintf(outfp, "%s\n", lne);
            free(lne);
        }

        fclose(outfp);
    }
}

void gs_set_hblocks_01_GV(int *n_hap, int *hap_list, char **hap_symbol,
                          char **out_filename, int *auxfiles, int *retval,
                          char **set_name)
{
    gs_set_hblocks_01(gs_fdta(*set_name), n_hap, hap_list, hap_symbol,
                      out_filename, auxfiles, retval);
}

void gs_cross_eval_gd_fct_01(gs_varset_type *gv, char **dist, int *split,
                             int *retval)
{
    if (strcmp(*dist, "rd") && strcmp(*dist, "mrd") && strcmp(*dist, "euc"))
    {
        gs_info(gv, -2, "Unknown genetic distance");
        *retval = -2;
        return;
    }


    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2,
                "No marker data available, run 'read.marker.data' first");
        *retval = -2;
        return;
    }

    gs_UNSIGNED sp = (gs_UNSIGNED)*split;

    if ((sp < 1) || (sp > (gv->NoInd) - 1))
    {
        gs_info(gv, -2, "Incorrect split value");
        *retval = -2;
        return;
    }

    gs_UNSIGNED i, j, k, l;

    int Srv = 0;
    int *Sretval = &Srv;
    gs_cross_init_01(gv, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    for (k = 0; k < gv->NCrs; k++)
        gv->crs[k]->gd = NA_REAL;

    gs_UNSIGNED *fact = malloc(sp * (gv->NoInd - sp) * sizeof(gs_UNSIGNED));

    k = 0;
    l = 0;
    for (i = 0; i < (gv->NoInd - 1); i++)
        for (j = i + 1; j < gv->NoInd; j++)
        {
            if ((i < sp) && (j >= sp))
            {
                fact[l] = k;
                l++;
            }
            k++;
        }

    if (0 == strcmp("mrd", *dist))
    {
#ifdef MT
#pragma omp parallel for
#endif

        for (gs_UNSIGNED l = 0; l < (sp * (gv->NoInd - sp)); l++)
            gs_calc_mrd(gv, fact[l]);
    }
    else if (0 == strcmp("rd", *dist))
    {
#ifdef MT
#pragma omp parallel for
#endif

        for (gs_UNSIGNED l = 0; l < (sp * (gv->NoInd - sp)); l++)
            gs_calc_rd(gv, fact[l]);
    }
    else
    {
#ifdef MT
#pragma omp parallel for
#endif

        for (gs_UNSIGNED l = 0; l < (sp * (gv->NoInd - sp)); l++)
            gs_calc_euc(gv, fact[l]);
    }

    free(fact);
    gv->cross_ready |= ST_CROSS_GD;
    *retval = 0;
    return;
}

void gs_cross_eval_gd_fct_01_GV(char **dist, int *split, int *retval,
                                char **set_name)
{
    gs_cross_eval_gd_fct_01(gs_fdta(*set_name), dist, split, retval);
}

void gs_recode_ref_02(gs_varset_type *gv, int *r1, int *r2, int *descending,
                      int *retval)

{
    gs_UNSIGNED i, j, idx, idx_R1, idx_R2;
    int a, b, c, d, x, y;

    if ((*r1 < 1) || (*r2 < 1) || ((gs_UNSIGNED)*r1 > gv->NoInd) ||
        ((gs_UNSIGNED)*r2 > gv->NoInd))
    {
        gs_info(gv, -2, "Invalid reference individual");
        *retval = -2;
        return;
    }

    GS_ALLELE_TYPE *Ref1 = malloc(gv->NoMar * gv->NoInd * sizeof(*Ref1));
    GS_ALLELE_TYPE *Ref2 = malloc(gv->NoMar * gv->NoInd * sizeof(*Ref2));
    if ((NULL == Ref1) || (NULL == Ref2))
        ERR_M;

    gs_UNSIGNED R1 = (gs_UNSIGNED)*r1 - 1;
    gs_UNSIGNED R2 = (gs_UNSIGNED)*r2 - 1;

    for (i = 0; i < gv->NoMar; i++)
    {
        idx_R1 = i * gv->NoInd + R1;
        a = st_allele_get(gv, idx_R1, 0);
        b = st_allele_get(gv, idx_R1, 1);

        idx_R2 = i * gv->NoInd + R2;
        x = st_allele_get(gv, idx_R2, 0);
        y = st_allele_get(gv, idx_R2, 1);

        for (j = 0; j < gv->NoInd; j++)
        {
            idx = i * gv->NoInd + j;
            c = st_allele_get(gv, idx, 0);
            d = st_allele_get(gv, idx, 1);

            if ((GS_ALLELE_MISSING == a) || (GS_ALLELE_MISSING == b) ||
                (GS_ALLELE_MISSING == c) || (GS_ALLELE_MISSING == d))
            {
                Ref1[idx] = GS_ALLELE_MISSING;
                Ref2[idx] = GS_ALLELE_MISSING;
            }
            else
            {
                if (0 == *descending)
                {

                    if (((c == a) && (d == b)) || ((c == b) && (d == a)))
                    {
                        Ref1[idx] = 1;
                        Ref2[idx] = 1;
                    }
                    else if (((c == x) && (d == y)) || ((c == y) && (d == x)))
                    {
                        Ref1[idx] = 2;
                        Ref2[idx] = 2;
                    }
                    else if (((c == a) || (c == b)) || ((d == a) || (d == b)))
                    {
                        Ref1[idx] = 1;
                        Ref2[idx] = 3;
                    }
                    else if (((c == x) || (c == y)) || ((d == x) || (d == y)))
                    {
                        Ref1[idx] = 2;
                        Ref2[idx] = 3;
                    }
                    else
                    {
                        Ref1[idx] = 3;
                        Ref2[idx] = 3;
                    }
                }
                else
                {
                    if (((c == a) || (c == b)) && ((d == a) || (d == b)))
                    {
                        Ref1[idx] = 1;
                        Ref2[idx] = 1;
                    }
                    else if (((c == a) || (c == b)) || ((d == a) || (d == b)))
                    {
                        Ref1[idx] = 1;
                        Ref2[idx] = 2;
                    }
                    else
                    {
                        Ref1[idx] = 2;
                        Ref2[idx] = 2;
                    }
                }
            }
        }
    }

    st_allele_replace_owned(gv, Ref1, Ref2);

    int Srv = 0;
    int *Sretval = &Srv;
    char *empty = "";
    char **ind_filename = &empty;
    char **mar_filename = &empty;
    char **gen_filename = &empty;
    int auf = 0;
    int *Sauxfiles = &auf;
    gs_marker_stats_01(gv, ind_filename, mar_filename, gen_filename, Sauxfiles,
                       Sauxfiles, Sauxfiles, Sauxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_recode_ref_02_GV(int *r1, int *r2, int *descending, int *retval,
                         char **set_name)
{
    gs_recode_ref_02(gs_fdta(*set_name), r1, r2, descending, retval);
}

extern int st_mxd_mv_scalar(int n, const double *y, int p, const double *X,
                            int k, double *const *V, double *b,
                            double *variance, int maxit, double tolerance,
                            char *message, size_t message_size);

void st_mxd(int *n, double *y, int *cX, double *X, int *kk, double **Z,
            double **VV, double *b, double *u, double *varcmp, int *maxcyc,
            double *tol, int *retval)
{
    (void)u;
    if (!retval)
        return;
    *retval = -2;
    if (!n || !y || !cX || !X || !kk || !varcmp || !maxcyc || !tol || Z ||
        *n < 2 || *cX < 0 || *kk < 0 || *kk == INT_MAX || (*kk && !VV) ||
        *maxcyc < 0 || *maxcyc == INT_MAX || !R_FINITE(*tol) || *tol <= 0)
    {
        st_info(-1, "Invalid arguments to st_mxd");
        return;
    }
    char message[512] = {0};
    *retval = st_mxd_mv_scalar(*n, y, *cX, X, *kk, VV, b, varcmp, *maxcyc,
                               *tol, message, sizeof(message));
    if (*retval < 0)
        st_info(-1, message[0] ? message : "Mixed-model fitting failed");
    else if (*retval)
        st_info(-1, "No convergence in variance component estimation");
}

SEXP st_get_Z(SEXP dataSet)
{
    char *set_name = (char *)CHAR(STRING_ELT(dataSet, 0));
    gs_varset_type *gv = gs_fdta(set_name);

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2, "No marker data available");
        return (R_NilValue);
    }

    if (NULL == gv->Z)
    {
        int Srv;
        int *Sretval = &Srv;
        gs_build_Z_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
            return (R_NilValue);
    }

    SEXP ZZ;
    R_len_t n, m;

    n = (R_len_t)gv->NoInd;
    m = (R_len_t)gv->NoMar;

    PROTECT(ZZ = allocMatrix(REALSXP, n, m));

    double *Z = REAL(ZZ);

    size_t ii = (size_t)gv->NoInd;
    size_t mm = (size_t)gv->NoMar;

    for (size_t i = 0; i < ii; i++)
    {
        for (size_t m = 0; m < mm; m++)
        {
            Z[i + m * ii] = gv->Z[i * mm + m];
        }
    }

    UNPROTECT(1);
    return (ZZ);
}

SEXP st_get_Z_named(SEXP dataSet, SEXP names)
{
    if (TYPEOF(dataSet) != STRSXP || XLENGTH(dataSet) != 1 ||
        STRING_ELT(dataSet, 0) == NA_STRING ||
        CHAR(STRING_ELT(dataSet, 0))[0] == '\0')
    {
        st_info(-2, "gs.get.Z: data.set must be a single nonempty string");
        return R_NilValue;
    }
    if (TYPEOF(names) != LGLSXP || XLENGTH(names) != 1 ||
        LOGICAL(names)[0] == NA_LOGICAL)
    {
        st_info(-2, "gs.get.Z: names must be TRUE or FALSE");
        return R_NilValue;
    }
    if (!LOGICAL(names)[0])
        return st_get_Z(dataSet);

    SEXP Z;
    PROTECT(Z = st_get_Z(dataSet));
    if (Z == R_NilValue)
    {
        UNPROTECT(1);
        return R_NilValue;
    }
    gs_varset_type *gv = gs_fdta((char *)CHAR(STRING_ELT(dataSet, 0)));
    if (gv->IndName == NULL)
    {
        gs_info(gv, -2, "No individual identifiers available");
        UNPROTECT(1);
        return R_NilValue;
    }

    int n = INTEGER(getAttrib(Z, R_DimSymbol))[0];
    SEXP labels, dimnames;
    PROTECT(labels = allocVector(STRSXP, n));
    for (int i = 0; i < n; i++)
        SET_STRING_ELT(labels, i, mkChar(gv->IndName[i]));
    PROTECT(dimnames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, labels);
    SET_VECTOR_ELT(dimnames, 1, R_NilValue);
    setAttrib(Z, R_DimNamesSymbol, dimnames);
    UNPROTECT(3);
    return Z;
}

SEXP st_get_V(SEXP dataSet)
{
    char *set_name = (char *)CHAR(STRING_ELT(dataSet, 0));
    gs_varset_type *gv = gs_fdta(set_name);

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2, "No marker data available");
        return (R_NilValue);
    }

    if (NULL == gv->Z)
    {
        int Srv;
        int *Sretval = &Srv;
        gs_build_Z_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
            return (R_NilValue);
    }

    if (NULL == gv->V)
    {
        int Srv;
        int *Sretval = &Srv;
        gs_build_V_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
            return (R_NilValue);
    }

    SEXP VV;

    R_len_t n;

    n = (R_len_t)gv->NoInd;

    PROTECT(VV = allocMatrix(REALSXP, n, n));

    double *V = REAL(VV);

    size_t ii = (size_t)gv->NoInd;

    memcpy(V, gv->V, ii * ii * sizeof(double));

    UNPROTECT(1);
    return (VV);
}

SEXP st_get_y(SEXP dataSet)
{
    char *set_name = (char *)CHAR(STRING_ELT(dataSet, 0));
    gs_varset_type *gv = gs_fdta(set_name);

    if (NULL == gv->y)
    {
        gs_info(gv, -2, "No performance data available");
        return (R_NilValue);
    }

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2, "No marker data available");
        return (R_NilValue);
    }

    SEXP yy;
    R_len_t n = (R_len_t)gv->NoInd;

    PROTECT(yy = allocVector(REALSXP, n));

    double *y = REAL(yy);

    size_t ii = (size_t)gv->NoInd;

    memcpy(y, gv->y, ii * sizeof(double));

    UNPROTECT(1);
    return (yy);
}

static SEXP st_get_y_named_values(SEXP dataSet, int names)
{
    if (!names)
        return st_get_y(dataSet);

    SEXP y;
    PROTECT(y = st_get_y(dataSet));
    if (y == R_NilValue)
    {
        UNPROTECT(1);
        return R_NilValue;
    }
    gs_varset_type *gv = gs_fdta((char *)CHAR(STRING_ELT(dataSet, 0)));
    if (gv->IndName == NULL)
    {
        gs_info(gv, -2, "No individual identifiers available");
        UNPROTECT(1);
        return R_NilValue;
    }

    SEXP labels;
    PROTECT(labels = allocVector(STRSXP, XLENGTH(y)));
    for (R_xlen_t i = 0; i < XLENGTH(y); i++)
        SET_STRING_ELT(labels, i, mkChar(gv->IndName[i]));
    setAttrib(y, R_NamesSymbol, labels);
    UNPROTECT(2);
    return y;
}

SEXP st_get_y_named(SEXP dataSet, SEXP names)
{
    if (TYPEOF(dataSet) != STRSXP || XLENGTH(dataSet) != 1 ||
        STRING_ELT(dataSet, 0) == NA_STRING ||
        CHAR(STRING_ELT(dataSet, 0))[0] == '\0')
    {
        st_info(-2, "gs.get.y: data.set must be a single nonempty string");
        return R_NilValue;
    }
    if (TYPEOF(names) != LGLSXP || XLENGTH(names) != 1 ||
        LOGICAL(names)[0] == NA_LOGICAL)
    {
        st_info(-2, "gs.get.y: names must be TRUE or FALSE");
        return R_NilValue;
    }
    return st_get_y_named_values(dataSet, LOGICAL(names)[0]);
}

SEXP st_lmm_marker_data(SEXP dataSet)
{
    if (TYPEOF(dataSet) != STRSXP || XLENGTH(dataSet) != 1 ||
        STRING_ELT(dataSet, 0) == NA_STRING ||
        CHAR(STRING_ELT(dataSet, 0))[0] == '\0')
        error("st.lmm: marker data set must be a single nonempty string");
    const char *name = CHAR(STRING_ELT(dataSet, 0));
    gs_varset_type *x = st_find_data_set_no_create(name);
    if (x == NULL)
        error("st.lmm: marker data set '%s' does not exist", name);
    if (!x->NoInd || !x->NoMar || x->IndName == NULL)
        error("st.lmm: no marker data available in data set '%s'", name);
    if (x->y == NULL)
        error("st.lmm: no performance data available in data set '%s'", name);
    if (x->NoInd > INT_MAX)
        error("st.lmm: marker data set has too many individuals");

    int n = (int)x->NoInd;
    SEXP y = PROTECT(st_get_y_named_values(dataSet, 1));
    SEXP gen = PROTECT(allocVector(INTSXP, n));
    for (int j = 0; j < n; j++)
        INTEGER(gen)[j] = j + 1;
    setAttrib(gen, R_LevelsSymbol, getAttrib(y, R_NamesSymbol));
    SEXP factor_class = PROTECT(mkString("factor"));
    setAttrib(gen, R_ClassSymbol, factor_class);
    SEXP data = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(data, 0, y);
    SET_VECTOR_ELT(data, 1, gen);
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("y"));
    SET_STRING_ELT(names, 1, mkChar("gen"));
    setAttrib(data, R_NamesSymbol, names);
    setAttrib(data, R_RowNamesSymbol, getAttrib(y, R_NamesSymbol));
    SEXP data_class = PROTECT(mkString("data.frame"));
    setAttrib(data, R_ClassSymbol, data_class);
    UNPROTECT(6);
    return data;
}

SEXP st_get_im(SEXP i, SEXP m, SEXP dataSet)
{
    char *set_name = (char *)CHAR(STRING_ELT(dataSet, 0));
    gs_varset_type *gv = gs_fdta(set_name);

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2, "No marker data available");
        return (R_NilValue);
    }

    int *ii = INTEGER(i);
    int *mm = INTEGER(m);

    *ii = (int)gv->NoInd;
    *mm = (int)gv->NoMar;

    return (R_NilValue);
}

void gs_lambda_rrblup_01(gs_varset_type *gv, int *maxiter, double *precision,
                         double *hsq, char **out_filename, int *auxfiles,
                         int *retval)
{

    (void)maxiter;
    (void)precision;
    (void)hsq;

    if (!gv->sigma || !gv->NoMar || !R_FINITE(gv->sigma[0]) ||
        !R_FINITE(gv->sigma[1]) || gv->sigma[0] < 0 || gv->sigma[1] < 0 ||
        (gv->sigma[0] == 0 && gv->sigma[1] == 0))
    {
        gs_info(gv, -2, "No valid variance components available");
        *retval = -2;
        return;
    }

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;

    double lambda = gv->sigma[1] == 0 ? R_PosInf : gv->sigma[0] / gv->sigma[1];

    for (gs_UNSIGNED i = 0; i < gv->NoMar; i++)
        gv->ladi[i] = lambda;

    if (1 == *auxfiles)
    {
        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_lambda_rrblup_01_GV(int *maxiter, double *precision, double *hsq,
                            char **out_filename, int *auxfiles, int *retval,
                            char **set_name)
{
    gs_lambda_rrblup_01(gs_fdta(*set_name), maxiter, precision, hsq,
                        out_filename, auxfiles, retval);
}

void gs_vc_rrblup_01(gs_varset_type *gv, int *maxiter, double *precision,
                     double *hsq, char **out_filename, int *auxfiles,
                     int *retval)
{

    char msg[gs_MAXSTR];
    int fit_status = 0;

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available");
        *retval = -2;
        return;
    }

    if (NULL == gv->y)
    {
        gs_info(gv, -2, "No performance data available");
        *retval = -2;
        return;
    }

    int k = 1;
    int n[1] = {(int)gv->NoInd};
    int cX[1] = {1};

    double *y = malloc(gv->NoInd * sizeof(double));
    double *X = malloc(gv->NoInd * sizeof(double));
    double **V = calloc(k, sizeof(double *));
    if (V)
        V[0] = malloc((size_t)gv->NoInd * gv->NoInd * sizeof(double));
    if (!y || !X || !V || !V[0])
    {
        free(y);
        free(X);
        if (V)
        {
            free(V[0]);
            free(V);
        }
        gs_info(gv, -2,
                "Memory allocation error in variance-component fitting");
        *retval = -2;
        return;
    }

    st_dgemm("N", "T", gv->NoInd, gv->NoMar, gv->NoInd, gv->NoMar, gv->Z, gv->Z,
             V[0]);

    for (int i = 0; i < *n; i++)
        X[i] = 1.0;

    st_dcopy(*n, y, gv->y);

    double varcmp[2] = {0, 0};

    double mean_y = 0;
    double var_y = 0;
    for (int i = 0; i < *n; i++)
        mean_y += y[i], var_y += (y[i]) * (y[i]);
    mean_y = mean_y / *n;
    var_y = (var_y - *n * mean_y * mean_y) / (*n - 1);

    double vg = *hsq * var_y;
    double ve = var_y - vg;

    varcmp[0] = vg / (double)gv->NoMar;
    varcmp[1] = ve;

    if (0 < *maxiter)
    {
        snprintf(msg, sizeof(msg),
                 "Start values from hsq = %4.2f: var(m) = %f, var(e) = %5f",
                 *hsq, varcmp[0], varcmp[1]);
        gs_info(gv, 2, msg);

        st_mxd(n, y, cX, X, &k, NULL, V, NULL, NULL, varcmp, maxiter, precision,
               retval);
        fit_status = *retval;
        if (*retval < 0)
        {
            free(y);
            free(X);
            free(V[0]);
            free(V);
            *retval = -2;
            return;
        }
    }

    snprintf(msg, sizeof(msg), "Variance components: var(m) = %f, var(e) = %5f",
             varcmp[0], varcmp[1]);

    gs_info(gv, 2, msg);

    double *sigma =
        realloc(gv->sigma, ((size_t)1 + gv->NoMar) * sizeof(double));
    if (!sigma)
    {
        free(y);
        free(X);
        free(V[0]);
        free(V);
        gs_info(gv, -2,
                "Memory allocation error in variance-component fitting");
        *retval = -2;
        return;
    }
    gv->sigma = sigma;

    gv->sigma[0] = varcmp[1];

    for (gs_UNSIGNED i = 1; i < (1 + gv->NoMar); i++)
        gv->sigma[i] = varcmp[0];

    free(y);
    free(X);
    free(V[0]);
    free(V);

    if (1 == *auxfiles)
    {
        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "vc \n");
        for (gs_UNSIGNED m = 0; m < 2; m++)
            fprintf(outfp, "%f\n", gv->sigma[m]);

        fclose(outfp);
    }
    *retval = fit_status;
    return;
}

void gs_vc_rrblup_01_GV(int *maxiter, double *precision, double *hsq,
                        char **out_filename, int *auxfiles, int *retval,
                        char **set_name)
{
    gs_vc_rrblup_01(gs_fdta(*set_name), maxiter, precision, hsq, out_filename,
                    auxfiles, retval);
}

void gs_single_marker_aov_01(gs_varset_type *gv, double *alpha,
                             char **out_filename, int *auxfiles, int *retval)
{

    (void)alpha;

    gs_UNSIGNED i, m, idx_Z;
    FILE *outfp;

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available, run 'build.Z' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->y)
    {
        gs_info(
            gv, -2,
            "No performance data available, run 'read.performance.data' first");
        *retval = -2;
        return;
    }

    int ext = 0;
    double zz;
    for (i = 0; i < gv->NoInd; i++)
        for (m = 0; m < gv->NoMar; m++)
        {
            idx_Z = i * gv->NoMar + m;
            zz = gv->Z[idx_Z];
            if ((0 != zz) && (1 != zz) && (2 != zz))
            {
                ext = 1;
                break;
            };
        }

    if (1 == ext)
    {
        gs_info(gv, -2, "Incorrect design matrix for single marker anova");
        *retval = -2;
        return;
    }

    gv->SMp = realloc(gv->SMp, gv->NoMar * sizeof(double));
    gv->SMF = realloc(gv->SMF, gv->NoMar * sizeof(double));
    gv->SMsqm = realloc(gv->SMsqm, gv->NoMar * sizeof(double));
    gv->SMsqe = realloc(gv->SMsqe, gv->NoMar * sizeof(double));
    gv->SMvar = realloc(gv->SMvar, gv->NoMar * sizeof(double));
    if ((NULL == gv->SMp) || (NULL == gv->SMF) || (NULL == gv->SMsqm) ||
        (NULL == gv->SMsqe) || (NULL == gv->SMvar))
        ERR_M;

    for (m = 0; m < gv->NoMar; m++)
        gv->SMp[m] = 0;

    double k = (double)(gv->NoInd);
    double sy = 0, syy = 0, SQT = 0, SQM = 0, SQE = 0, mm, kk;

    double sy0, sy1, sy2, ny0, ny1, ny2;

    for (i = 0; i < gv->NoInd; i++)
    {
        sy += gv->y[i];
        syy += gv->y[i] * gv->y[i];
    }
    SQT = syy - sy * sy / k;

    for (m = 0; m < gv->NoMar; m++)
    {
        sy0 = 0;
        sy1 = 0;
        sy2 = 0;
        ny0 = 0;
        ny1 = 0;
        ny2 = 0;

        for (i = 0; i < gv->NoInd; i++)
        {
            idx_Z = i * gv->NoMar + m;
            if (0 == gv->Z[idx_Z])
            {
                ny0++;
                sy0 += gv->y[i];
            }
            else if (2 == gv->Z[idx_Z])
            {
                ny2++;
                sy2 += gv->y[i];
            }
            else if (1 == gv->Z[idx_Z])
            {
                ny1++;
                sy1 += gv->y[i];
            }
        }

        mm = (0 == ny0 ? 0 : 1) + (0 == ny1 ? 0 : 1) + (0 == ny2 ? 0 : 1);

        SQE = syy - (0 == ny0 ? 0 : sy0 * sy0 / ny0) -
              (0 == ny1 ? 0 : sy1 * sy1 / ny1) -
              (0 == ny2 ? 0 : sy2 * sy2 / ny2);

        SQM = SQT - SQE;

        gv->SMsqm[m] = SQM;
        gv->SMsqe[m] = SQE;

        if ((2 > mm) || ((k - mm) < 1))
        {
            gv->SMvar[m] = FLT_MIN;
            gv->SMF[m] = 0;
            gv->SMp[m] = 1;
        }
        else
        {
            kk = (2 / (mm - 1)) * (k - (ny0 * ny0 + ny1 * ny1 + ny2 * ny2) / k);
            gv->SMvar[m] = (SQM / (mm - 1) - SQE / (k - mm)) / kk;

            gv->SMF[m] = (SQM / (mm - 1)) / (SQT / (k - mm));
            gv->SMp[m] = pf(gv->SMF[m], (int)(mm - 1), (int)(k - mm), 0, 0);
        }
    }

    double vcmin;
    vcmin = fabs(gv->SMvar[0]);
    for (m = 1; m < gv->NoMar; m++)
        if (fabs(gv->SMvar[m]) < vcmin)
            vcmin = fabs(gv->SMvar[m]);
    vcmin = vcmin / 1000.0;

    for (m = 0; m < gv->NoMar; m++)
        if (gv->SMvar[m] < 0)
            gv->SMvar[m] = vcmin;

    if (1 == *auxfiles)
    {

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "SQM SQE var F p\n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f %f %f %f %f\n", gv->MarName[m], gv->SMsqm[m],
                    gv->SMsqe[m], gv->SMvar[m], gv->SMF[m], gv->SMp[m]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_single_marker_aov_01_GV(double *alpha, char **out_filename,
                                int *auxfiles, int *retval, char **set_name)
{
    gs_single_marker_aov_01(gs_fdta(*set_name), alpha, out_filename, auxfiles,
                            retval);
}

void gs_lambda_rmla_02(gs_varset_type *gv, double *alpha, int *maxiter,
                       double *precision, double *hsq, char **out_filename,
                       int *auxfiles, int *retval)
{
    gs_UNSIGNED m;
    FILE *outfp;

    (void)alpha;
    (void)maxiter;
    (void)precision;
    (void)hsq;

    if (NULL == gv->ladi)
    {
        gs_info(gv, -2,
                "No shrinkage factors available that could be adjusted");
        *retval = -2;
        return;
    }

    long double Svar = 0.0;
    for (m = 0; m < gv->NoMar; m++)
        Svar += gv->SMvar[m];

    Svar = Svar / (long double)gv->NoMar;

    for (m = 0; m < gv->NoMar; m++)
        if (R_FINITE(gv->ladi[m]) && gv->ladi[m] != 0)
            gv->ladi[m] = (Svar / gv->SMvar[m]) * gv->ladi[m];

    if (1 == *auxfiles)
    {
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_lambda_rmla_02_GV(double *alpha, int *maxiter, double *precision,
                          double *hsq, char **out_filename, int *auxfiles,
                          int *retval, char **set_name)
{
    gs_lambda_rmla_02(gs_fdta(*set_name), alpha, maxiter, precision, hsq,
                      out_filename, auxfiles, retval);
}

void gs_w_aov_01(gs_varset_type *gv, double *alpha, char **out_filename,
                 int *auxfiles, int *retval)

{
    gs_UNSIGNED m;
    FILE *outfp;
    int Srv;
    int *Sretval = &Srv;

    gs_single_marker_aov_01(gv, alpha, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    long double Svar = 0.0;
    for (m = 0; m < gv->NoMar; m++)
        Svar += gv->SMvar[m];

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;

    for (m = 0; m < gv->NoMar; m++)
        gv->ladi[m] = (gv->SMvar[m] / Svar) * (double)gv->NoMar;

    if (1 == *auxfiles)
    {
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "W \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_w_aov_01_GV(double *alpha, char **out_filename, int *auxfiles,
                    int *retval, char **set_name)
{
    gs_w_aov_01(gs_fdta(*set_name), alpha, out_filename, auxfiles, retval);
}

void gs_w_I_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
               int *retval)

{
    gs_UNSIGNED m;
    FILE *outfp;

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;

    for (m = 0; m < gv->NoMar; m++)
        gv->ladi[m] = 1;

    if (1 == *auxfiles)
    {
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "W \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_w_I_01_GV(char **out_filename, int *auxfiles, int *retval,
                  char **set_name)
{
    gs_w_I_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_lambda_const_02(gs_varset_type *gv, double *hsq, char **out_filename,
                        int *auxfiles, int *retval)
{
    gs_UNSIGNED i, m;

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available, run 'build.Z' first");
        *retval = -2;
        return;
    }

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;

    double lambda = (1 / *hsq - 1) * (double)gv->NoMar;

    for (i = 0; i < gv->NoMar; i++)
        gv->ladi[i] = lambda;

    if (1 == *auxfiles)
    {
        FILE *outfp;

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_lambda_const_02_GV(double *hsq, char **out_filename, int *auxfiles,
                           int *retval, char **set_name)
{
    gs_lambda_const_02(gs_fdta(*set_name), hsq, out_filename, auxfiles, retval);
}

void gs_lambda_aov_01(gs_varset_type *gv, double *hsq, double *alpha,
                      char **out_filename, int *auxfiles, int *retval)
{
    gs_UNSIGNED m;
    FILE *outfp;
    int Srv;
    int *Sretval = &Srv;

    gs_single_marker_aov_01(gv, alpha, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;

    double Svar = 0;
    for (m = 0; m < gv->NoMar; m++)
        Svar += gv->SMvar[m];

    for (m = 0; m < gv->NoMar; m++)
    {
        gv->ladi[m] = ((1 / *hsq) - 1) * Svar / gv->SMvar[m];
    }

    if (1 == *auxfiles)
    {
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_lambda_aov_01_GV(double *hsq, double *alpha, char **out_filename,
                         int *auxfiles, int *retval, char **set_name)
{
    gs_lambda_aov_01(gs_fdta(*set_name), hsq, alpha, out_filename, auxfiles,
                     retval);
}

void gs_single_marker_reg_01(gs_varset_type *gv, char **out_filename,
                             int *auxfiles, int *retval)
{
    gs_UNSIGNED i, m, idx_Z;
    FILE *outfp;

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available, run 'build.Z' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->y)
    {
        gs_info(
            gv, -2,
            "No performance data available, run 'read.performance.data' first");
        *retval = -2;
        return;
    }

    gv->SMp = realloc(gv->SMp, gv->NoMar * sizeof(double));
    gv->SMF = realloc(gv->SMF, gv->NoMar * sizeof(double));
    gv->SMsqm = realloc(gv->SMsqm, gv->NoMar * sizeof(double));
    gv->SMsqe = realloc(gv->SMsqe, gv->NoMar * sizeof(double));
    if ((NULL == gv->SMp) || (NULL == gv->SMF) || (NULL == gv->SMsqm) ||
        (NULL == gv->SMsqe))
        ERR_M;

    for (m = 0; m < gv->NoMar; m++)
        gv->SMp[m] = 0;

    double k = (double)(gv->NoInd);
    double sy = 0, syy = 0, sx = 0, sxx = 0, sxy = 0, b, SQT = 0, SQM = 0,
           SQE = 0, F, p;

    for (i = 0; i < gv->NoInd; i++)
    {
        sy += gv->y[i];
        syy += gv->y[i] * gv->y[i];
    }
    SQT = syy - sy * sy / k;

    for (m = 0; m < gv->NoMar; m++)
    {
        sx = 0;
        sxx = 0;
        sxy = 0;
        for (i = 0; i < gv->NoInd; i++)
        {
            idx_Z = i * gv->NoMar + m;
            sx += gv->Z[idx_Z];
            sxx += gv->Z[idx_Z] * gv->Z[idx_Z];
            sxy += gv->Z[idx_Z] * gv->y[i];
        }

        b = (sxy - sx * sy / k) / (sxx - sx * sx / k);

        SQM = b * (sxy - sx * sy / k);
        SQE = SQT - SQM;
        F = SQM / (SQE / (k - 2));
        p = pf(F, 1, (int)gv->NoInd, 0, 0);

        gv->SMsqm[m] = SQM;
        gv->SMsqe[m] = SQE;
        gv->SMF[m] = F;
        gv->SMp[m] = p;
    }

    if (1 == *auxfiles)
    {

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "SQM SQE F p\n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f %f %f %f\n", gv->MarName[m], gv->SMsqm[m],
                    gv->SMsqe[m], gv->SMF[m], gv->SMp[m]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_single_marker_reg_01_GV(char **out_filename, int *auxfiles, int *retval,
                                char **set_name)
{
    gs_single_marker_reg_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_lambda_reg_01(gs_varset_type *gv, double *hsq, char **out_filename,
                      int *auxfiles, int *retval)
{
    gs_UNSIGNED m;
    FILE *outfp;
    int Srv;
    int *Sretval = &Srv;

    gs_single_marker_reg_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;

    double SSQM = 0;
    for (m = 0; m < gv->NoMar; m++)
        SSQM += gv->SMsqm[m];

    for (m = 0; m < gv->NoMar; m++)
        gv->ladi[m] = ((1 / *hsq) - 1) * SSQM / gv->SMsqm[m];

    if (1 == *auxfiles)
    {
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_lambda_reg_01_GV(double *hsq, char **out_filename, int *auxfiles,
                         int *retval, char **set_name)
{
    gs_lambda_reg_01(gs_fdta(*set_name), hsq, out_filename, auxfiles, retval);
}

void gs_mme_coeff_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                     int *retval)
{

    gs_UNSIGNED i, j, idx_A;
    FILE *outfp;

    if (NULL == gv->y)
    {
        gs_info(
            gv, -2,
            "No performance data available, run 'read.performance.data' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available, run 'build.Z' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->ladi)
    {
        gs_info(gv, -2, "No shrinkage factor defined");
        *retval = -2;
        return;
    }

    gv->dimA = 1 + gv->NoMar;
    gv->A = realloc(gv->A, gv->dimA * gv->dimA * sizeof(double));
    if (NULL == gv->A)
        ERR_M;

    for (i = 0; i < gv->dimA; i++)
        for (j = 0; j < gv->dimA; j++)
            gv->A[i * gv->dimA + j] = 0;

    double *one_vec;

    one_vec = malloc(gv->NoInd * sizeof(double));
    if (NULL == one_vec)
        ERR_M;
    for (i = 0; i < gv->NoInd; i++)
        one_vec[i] = 1;

    gv->A[0] = (double)gv->NoInd;

    double *tmp1;
    tmp1 = malloc(gv->NoMar * sizeof(double));
    if (NULL == tmp1)
        ERR_M;

    st_dgemm("T", "N", gv->NoInd, 1, gv->NoInd, gv->NoMar, one_vec, gv->Z,
             tmp1);

    for (j = 0; j < gv->NoMar; j++)
    {
        gv->A[0 * gv->dimA + (1 + j)] = tmp1[j];
        gv->A[(j + 1) * gv->dimA + 0] = tmp1[j];
    }

    double *tmp2;
    tmp2 = malloc(gv->NoMar * gv->NoMar * sizeof(double));
    if (NULL == tmp2)
        ERR_M;

    st_dgemm("T", "N", gv->NoInd, gv->NoMar, gv->NoInd, gv->NoMar, gv->Z, gv->Z,
             tmp2);

    for (i = 0; i < gv->NoMar; i++)
        for (j = 0; j < gv->NoMar; j++)
        {
            gv->A[(1 + i) * gv->dimA + (1 + j)] = tmp2[i * gv->NoMar + j];
        }

    if (1)
    {
        gv->As = realloc(gv->As, gv->dimA * gv->dimA * sizeof(double));
        if (NULL == gv->As)
            ERR_M;
        for (i = 0; i < gv->dimA; i++)
            for (j = 0; j < gv->dimA; j++)
            {
                gv->As[i * gv->dimA + j] = gv->A[i * gv->dimA + j];
            }
    }

    for (i = 1; i < gv->dimA; i++)
    {
        if (!R_FINITE(gv->ladi[i - 1]) && gv->ladi[i - 1] > 0)
        {
            for (j = 0; j < gv->dimA; j++)
                gv->A[i * gv->dimA + j] = gv->A[j * gv->dimA + i] = 0;
            gv->A[i * gv->dimA + i] = 1;
        }
        else
            gv->A[i * gv->dimA + i] += gv->ladi[i - 1];
    }

    free(one_vec);
    free(tmp1);
    free(tmp2);

    if (1 == *auxfiles)
    {
        if (NULL == (outfp = fopen(*out_filename, "w")))
            ERR_F;

        for (i = 0; i < gv->dimA; i++)
        {
            for (j = 0; j < gv->dimA; j++)
            {
                idx_A = i * gv->dimA + j;
                fprintf(outfp, "%f ", gv->A[idx_A]);
            }
            fprintf(outfp, "\n");
        }

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_mme_coeff_01_GV(char **out_filename, int *auxfiles, int *retval,
                        char **set_name)
{
    gs_mme_coeff_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_mme_restcoeff_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                         int *retval)
{

    gs_UNSIGNED i, j, idx_A;
    FILE *outfp;

    int Srv;
    int *Sretval = &Srv;

    if (NULL == gv->ladi)
    {
        gs_info(gv, -2, "No shrinkage factor defined");
        *retval = -2;
        return;
    }

    if (NULL == gv->As)
    {

        gs_mme_coeff_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }
    else
    {

        gv->dimA = 1 + gv->NoMar;
        gv->A = realloc(gv->A, gv->dimA * gv->dimA * sizeof(double));
        if (NULL == gv->A)
            ERR_M;

        memcpy(gv->A, gv->As,
               (size_t)gv->dimA * (size_t)gv->dimA * sizeof(double));

        for (i = 1; i < gv->dimA; i++)
            gv->A[i * gv->dimA + i] += gv->ladi[i - 1];
    }

    if (1 == *auxfiles)
    {
        if (NULL == (outfp = fopen(*out_filename, "w")))
            ERR_F;

        for (i = 0; i < gv->dimA; i++)
        {
            for (j = 0; j < gv->dimA; j++)
            {
                idx_A = i * gv->dimA + j;
                fprintf(outfp, "%f ", gv->A[idx_A]);
            }
            fprintf(outfp, "\n");
        }

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_mme_restcoeff_01_GV(char **out_filename, int *auxfiles, int *retval)
{
    gs_mme_restcoeff_01(gv, out_filename, auxfiles, retval);
}

void gs_mme_rhs_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                   int *retval)
{
    gs_UNSIGNED i, j, dim_b;
    FILE *outfp;

    if (NULL == gv->y)
    {
        gs_info(
            gv, -2,
            "No performance data available, run 'read.performance.data' first");
        *retval = -2;
        return;
    }

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available, run 'build.Z' first");
        *retval = -2;
        return;
    }

    dim_b = 1 + gv->NoMar;
    gv->b = realloc(gv->b, dim_b * sizeof(double));
    if (NULL == gv->b)
        ERR_M;
    for (i = 0; i < dim_b; i++)
        gv->b[i] = 0;

    double *one_vec;
    one_vec = malloc(gv->NoInd * sizeof(double));
    if (NULL == one_vec)
        ERR_M;
    for (i = 0; i < gv->NoInd; i++)
        one_vec[i] = 1;

    double *tmp;
    tmp = malloc(1 * sizeof(double));
    if (NULL == tmp)
        ERR_M;

    st_dgemm("T", "N", gv->NoInd, 1, gv->NoInd, 1, one_vec, gv->y, tmp);

    gv->b[0] = tmp[0];

    double *tmp2;
    tmp2 = malloc(gv->NoMar * sizeof(double));
    if (NULL == tmp2)
        ERR_M;

    st_dgemm("T", "N", gv->NoInd, gv->NoMar, gv->NoInd, 1, gv->Z, gv->y, tmp2);

    for (j = 0; j < gv->NoMar; j++)
        gv->b[1 + j] = (gv->ladi && !R_FINITE(gv->ladi[j]) && gv->ladi[j] > 0)
                           ? 0
                           : tmp2[j];

    if (1)
    {
        gv->bs = realloc(gv->bs, dim_b * sizeof(double));
        if (NULL == gv->bs)
            ERR_M;
        for (j = 0; j < dim_b; j++)
            gv->bs[j] = gv->b[j];
    }

    free(one_vec);
    free(tmp);
    free(tmp2);

    if (1 == *auxfiles)
    {

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        for (i = 0; i < dim_b; i++)
            fprintf(outfp, "%f \n", gv->b[i]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_mme_rhs_01_GV(char **out_filename, int *auxfiles, int *retval,
                      char **set_name)
{
    gs_mme_rhs_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_mme_restrhs_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                       int *retval)
{
    gs_UNSIGNED i, dim_b;
    FILE *outfp;

    int Srv;
    int *Sretval = &Srv;

    dim_b = 1 + gv->NoMar;

    if (NULL == gv->bs)
    {
        gs_mme_rhs_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }
    else
    {
        gv->b = realloc(gv->b, dim_b * sizeof(double));
        if (NULL == gv->b)
            ERR_M;

        memcpy(gv->b, gv->bs, (size_t)dim_b);
    }

    if (1 == *auxfiles)
    {

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        for (i = 0; i < dim_b; i++)
            fprintf(outfp, "%f \n", gv->b[i]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_mme_restrhs_01_GV(char **out_filename, int *auxfiles, int *retval)
{
    gs_mme_restrhs_01(gv, out_filename, auxfiles, retval);
}

void gs_mme_solve_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                     int *retval)
{
    gs_UNSIGNED i;
    FILE *outfp;

    int Srv = 0;
    int *Sretval = &Srv;

    if (NULL == gv->A)
    {
        gs_info(gv, -2, "No coefficient matrix available");
        *retval = -2;
        return;
    }

    if (NULL == gv->b)
    {
        gs_info(gv, -2, "No right hand side available");
        *retval = -2;
        return;
    }

    gv->dimu = 1 + gv->NoMar;
    gv->u = realloc(gv->u, gv->dimu * sizeof(double));
    if (NULL == gv->u)
        ERR_M;

    int boundary = gs_rr_boundary_01(gv, NULL, retval);
    if (boundary < 0)
        return;
    if (!boundary)
        st_dgesv(gv->dimA, gv->dimA, gv->A, gv->u, gv->b, Sretval);

    free(gv->A);
    gv->A = NULL;

    if (*Sretval == -2)
    {
        *retval = -2;
        return;
    }

    if (1 == *auxfiles)
    {

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "effect\n");
        for (i = 0; i < gv->dimu; i++)
            fprintf(outfp, "%s.%i %f \n", gv->EffNme[i], gv->EffAll[i],
                    gv->u[i]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_mme_solve_01_GV(char **out_filename, int *auxfiles, int *retval,
                        char **set_name)
{
    gs_mme_solve_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_esteff_rrwr_01(gs_varset_type *gv, double *hsq, char **out_filename,
                       int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_build_Z_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_reg_01(gv, hsq, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_esteff_rrwr_01_GV(double *hsq, char **out_filename, int *auxfiles,
                          int *retval, char **set_name)
{
    gs_esteff_rrwr_01(gs_fdta(*set_name), hsq, out_filename, auxfiles, retval);
}

void gs_esteff_rrwa_01(gs_varset_type *gv, double *hsq, double *alpha,
                       char **out_filename, int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_build_Z_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_aov_01(gv, hsq, alpha, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_esteff_rrwa_01_GV(double *hsq, double *alpha, char **out_filename,
                          int *auxfiles, int *retval, char **set_name)
{
    gs_esteff_rrwa_01(gs_fdta(*set_name), hsq, alpha, out_filename, auxfiles,
                      retval);
}

void gs_esteff_rrwe_01(gs_varset_type *gv, double *lambda, char **out_filename,
                       int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_build_Z_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_const_01(gv, lambda, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_esteff_rrwe_01_GV(double *lambda, char **out_filename, int *auxfiles,
                          int *retval, char **set_name)
{
    gs_esteff_rrwe_01(gs_fdta(*set_name), lambda, out_filename, auxfiles,
                      retval);
}

void gs_esteff_rrwe_02(gs_varset_type *gv, double *hsq, char **out_filename,
                       int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_build_Z_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_const_02(gv, hsq, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_esteff_rrwe_02_GV(double *hsq, char **out_filename, int *auxfiles,
                          int *retval, char **set_name)
{
    gs_esteff_rrwe_02(gs_fdta(*set_name), hsq, out_filename, auxfiles, retval);
}

void gs_lambda_emstep_01(gs_varset_type *gv, int *constvar, char **out_filename,
                         int *auxfiles, int *retval)
{
    gs_UNSIGNED i, j, m;

    if (NULL == gv->ladi)
    {
        gs_info(gv, -2, "No lambda defined");
        *retval = -2;
        return;
    }

    int Srv;
    int *Sretval = &Srv;

    gs_mme_restcoeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_restrhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gv->dimu = 1 + gv->NoMar;
    gv->u = realloc(gv->u, gv->dimu * sizeof(double));
    if (NULL == gv->u)
        ERR_M;

    for (i = 0; i < gv->dimu; i++)
        gv->u[i] = gv->b[i];

    gv->t0 = realloc(gv->t0, gv->dimA * sizeof(double));

    {
        double d;

        gv->u0 = realloc(gv->u0, gv->dimA * sizeof(gs_UNSIGNED));

        st_dgetrf(gv->A, gv->dimA, gv->u0, &d);

        st_dgetrs(gv->A, gv->dimA, gv->u0, gv->u);

        gv->t4 = realloc(gv->t4, gv->dimA * gv->dimA * sizeof(double));

        st_dgetri(gv->A, gv->dimA, gv->u0, gv->t4);

        for (j = 0; j < gv->dimA; j++)
            gv->t0[j] = gv->t4[gv->dimA * j + j];
    }

    gv->sigma = realloc(gv->sigma, gv->dimA * sizeof(double));
    if (NULL == gv->sigma)
        ERR_M;

    gv->t1 = realloc(gv->t1, sizeof(double));
    gv->t2 = realloc(gv->t2, ((gv->NoInd > gv->NoMar) ? gv->NoInd : gv->NoMar) *
                                 sizeof(double));
    gv->t3 = realloc(gv->t3, ((gv->NoInd > gv->NoMar) ? gv->NoInd : gv->NoMar) *
                                 sizeof(double));

    for (i = 0; i < gv->NoInd; i++)
        gv->t2[i] = 1;

    st_dgemm("T", "N", gv->NoInd, 1, gv->NoInd, 1, gv->y, gv->y, gv->sigma);

    st_dgemm("T", "N", gv->NoInd, 1, gv->NoInd, 1, gv->t2, gv->y, gv->t1);

    st_dgemm("T", "N", gv->NoInd, gv->NoMar, gv->NoInd, 1, gv->Z, gv->y,
             gv->t2);

    st_dgemm("T", "N", gv->NoMar, 1, gv->NoMar, 1, 1 + (gv->u), gv->t2, gv->t3);

    if (0 == *constvar)
    {

        gv->sigma[0] = gv->sigma_e;
        for (i = 1; i < gv->dimA; i++)
        {
            gv->sigma[i] = gv->u[i] * gv->u[i] + gv->sigma[0] * gv->t0[i];
        }
    }
    else
    {
        gv->sigma[0] -= gv->u[0] * gv->t1[0];
        gv->sigma[0] -= gv->t3[0];
        gv->sigma[0] = gv->sigma[0] / (double)gv->NoInd;
        gv->sigma_e = gv->sigma[0];

        st_dgemm("T", "N", gv->NoMar, 1, gv->NoMar, 1, 1 + (gv->u), 1 + (gv->u),
                 gv->t2);

        double tr = 0;
        for (i = 1; i < gv->dimA; i++)
            tr += gv->t0[i];

        for (i = 1; i < gv->dimA; i++)
            gv->sigma[i] = (gv->t2[0] + gv->sigma[0] * tr) / (double)gv->NoMar;
    }

    for (i = 1; i < gv->dimA; i++)
        gv->ladi[i - 1] = gv->sigma[0] / gv->sigma[i];

    if (1 == *auxfiles)
    {
        FILE *outfp;

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda var \n");
        fprintf(outfp, "err %f %f \n", 0.0, gv->sigma[0]);
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f %f\n", gv->MarName[m], gv->ladi[m],
                    gv->sigma[1 + m]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_lambda_emstep_01_GV(int *constvar, char **out_filename, int *auxfiles,
                            int *retval, char **set_name)

{
    gs_lambda_emstep_01(gs_fdta(*set_name), constvar, out_filename, auxfiles,
                        retval);
}

void gs_lambda_rml_01(gs_varset_type *gv, int *maxiter, double *precision,
                      int *constvar, int *init_lambda, double *hsq,
                      char **out_filename, int *auxfiles, int *retval)
{
    gs_UNSIGNED m;
    int i;
    double qu, maxqu = 0;
    int Srv;
    int *Sretval = &Srv;

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available");
        *retval = -2;
        return;
    }

    double *vc_save = malloc((1 + gv->NoMar) * sizeof(double));
    if (NULL == vc_save)
        ERR_M;

    for (m = 0; m < 1 + gv->NoMar; m++)
        vc_save[m] = 0;

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;

    if (1 == *init_lambda)
    {

        double lambda = (1 / *hsq - 1) * (double)gv->NoMar;
        for (m = 0; m < gv->NoMar; m++)
            gv->ladi[m] = lambda;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    for (i = 0; i < *maxiter; i++)
    {
        if (i > 0)
            for (m = 0; m < 1 + gv->NoMar; m++)
                vc_save[m] = gv->sigma[m];

        gs_lambda_emstep_01(gv, constvar, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        maxqu = 0;
        for (m = 0; m < 1 + gv->NoMar; m++)
        {
            qu = fabs((vc_save[m] - gv->sigma[m]) / gv->sigma[m]);
            if (qu > maxqu)
                maxqu = qu;
        }
        if (maxqu < *precision)
            break;
    }

    snprintf(gv->msg, sizeof(gv->msg),
             "Number of iterations: %i, precision: %f", i, maxqu);

    gs_info(gv, 1, gv->msg);
    if (1 == *constvar)
    {
        snprintf(gv->msg, sizeof(gv->msg),
                 "Variance components: var(m) = %f, var(e) = %5f", gv->sigma[1],
                 gv->sigma[0]);

        gs_info(gv, 1, gv->msg);
    }

    free(vc_save);

    if (1 == *auxfiles)
    {
        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_lambda_rmlc_01(gs_varset_type *gv, int *maxiter, double *precision,
                       double *hsq, char **out_filename, int *auxfiles,
                       int *retval)
{
    int _cc = 1;
    int *constvar = &_cc;
    int _il = 1;
    int *init_lambda = &_il;

    gs_lambda_rml_01(gv, maxiter, precision, constvar, init_lambda, hsq,
                     out_filename, auxfiles, retval);
}

void gs_lambda_rmlv_01(gs_varset_type *gv, int *maxiter, double *precision,
                       double *hsq, char **out_filename, int *auxfiles,
                       int *retval)
{

    int Srv = 0;
    int *Sretval = &Srv;

    int _cc = 1;
    int *constvar = &_cc;
    int _il = 1;
    int *init_lambda = &_il;

    _cc = 1;
    _il = 1;
    gs_lambda_rml_01(gv, maxiter, precision, constvar, init_lambda, hsq, nixp,
                     neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    _cc = 0;
    _il = 1;
    gs_lambda_rml_01(gv, maxiter, precision, constvar, init_lambda, hsq,
                     out_filename, auxfiles, retval);
}

void gs_lambda_rmlv_01_GV(int *maxiter, double *precision, double *hsq,
                          char **out_filename, int *auxfiles, int *retval,
                          char **set_name)
{
    gs_lambda_rmlv_01(gs_fdta(*set_name), maxiter, precision, hsq, out_filename,
                      auxfiles, retval);
}

void gs_lambda_rmlc_01_GV(int *maxiter, double *precision, double *hsq,
                          char **out_filename, int *auxfiles, int *retval,
                          char **set_name)
{
    gs_lambda_rmlc_01(gs_fdta(*set_name), maxiter, precision, hsq, out_filename,
                      auxfiles, retval);
}

void gs_lambda_rmlr_01(gs_varset_type *gv, int *maxiter, double *precision,
                       double *hsq, char **out_filename, int *auxfiles,
                       int *retval)
{
    gs_UNSIGNED m;
    FILE *outfp;
    int Srv;
    int *Sretval = &Srv;

    gs_single_marker_reg_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_rmlc_01(gv, maxiter, precision, hsq, out_filename, auxfiles,
                      Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    long double SSQM = 0.0;
    for (m = 0; m < gv->NoMar; m++)
        SSQM += gv->SMsqm[m];

    SSQM = SSQM / (long double)gv->NoMar;

    for (m = 0; m < gv->NoMar; m++)
        gv->ladi[m] = (SSQM / gv->SMsqm[m]) * gv->ladi[m];

    if (1 == *auxfiles)
    {
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_lambda_rmlr_01_GV(int *maxiter, double *precision, double *hsq,
                          char **out_filename, int *auxfiles, int *retval,
                          char **set_name)
{
    gs_lambda_rmlr_01(gs_fdta(*set_name), maxiter, precision, hsq, out_filename,
                      auxfiles, retval);
}

void gs_lambda_rmla_01(gs_varset_type *gv, double *alpha, int *maxiter,
                       double *precision, double *hsq, char **out_filename,
                       int *auxfiles, int *retval)
{
    gs_UNSIGNED m;
    FILE *outfp;
    int Srv;
    int *Sretval = &Srv;

    gs_single_marker_aov_01(gv, alpha, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_rmlc_01(gv, maxiter, precision, hsq, out_filename, auxfiles,
                      Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    long double Svar = 0.0;
    for (m = 0; m < gv->NoMar; m++)
        Svar += gv->SMvar[m];

    Svar = Svar / (long double)gv->NoMar;

    for (m = 0; m < gv->NoMar; m++)
        gv->ladi[m] = (Svar / gv->SMvar[m]) * gv->ladi[m];

    if (1 == *auxfiles)
    {
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_lambda_rmla_01_GV(double *alpha, int *maxiter, double *precision,
                          double *hsq, char **out_filename, int *auxfiles,
                          int *retval, char **set_name)
{
    gs_lambda_rmla_01(gs_fdta(*set_name), alpha, maxiter, precision, hsq,
                      out_filename, auxfiles, retval);
}

void gs_esteff_rmlc_01(gs_varset_type *, int *, double *, double *, char **,
                       int *, int *);

void gs_esteff_lsq_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                      int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    if ((1 + gv->NoMar) > gv->NoInd)
    {
        gs_info(gv, -2, "p > n");
        *retval = -2;
        return;
    }

    gs_build_Z_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    double lambda = 0;
    double *lambdaptr = &lambda;

    gs_lambda_const_01(gv, lambdaptr, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

#include <R_ext/Applic.h>

static int gs_lsq_rank_ready(gs_varset_type *x)
{
    int status = 0;
    gs_build_Z_01(x, nixp, neinp, &status);
    if (status == -2)
        return 0;
    if (x->NoInd > INT_MAX || x->NoMar > INT_MAX)
        Rf_error("Design matrix exceeds the QR integer range");
    int n = (int)x->NoInd, p = (int)x->NoMar, rank = 0;
    size_t cells;
    if (gs_size_mul((size_t)n, (size_t)p, &cells))
        Rf_error("Design matrix dimensions overflow");
    void *mark = st_work_mark();
    double *matrix = st_work_alloc(cells, sizeof(double));
    double *auxiliary = st_work_alloc((size_t)p, sizeof(double));
    double *work = st_work_alloc((size_t)p * 2, sizeof(double));
    int *pivot = st_work_alloc((size_t)p, sizeof(int));
    double tolerance = 1e-7;
    char formatted[DBL_MAX_10_EXP + 32];
    double values[3] = {x->AA, x->aa, x->aA};
    for (int j = 0; j < 3; j++)
    {
        snprintf(formatted, sizeof(formatted), "%.6f", values[j]);
        values[j] = strtod(formatted, NULL);
    }
    for (int j = 0; j < p; j++)
    {
        pivot[j] = j + 1;
        for (int i = 0; i < n; i++)
        {
            double value = x->Z[(size_t)i * p + j];
            matrix[i + (size_t)n * j] = value == x->AA   ? values[0]
                                        : value == x->aa ? values[1]
                                                         : values[2];
        }
    }
    F77_CALL(dqrdc2)(matrix, &n, &n, &p, &tolerance, &rank, auxiliary, pivot,
                     work);
    st_work_reset(mark);
    if (rank < p - 1)
    {
        snprintf(x->msg, sizeof(x->msg),
                 "Design matrix is not of full rank (%i < %i)", rank, p - 1);
        gs_info(x, -2, x->msg);
        return 0;
    }
    return 1;
}

void gs_esteff_lsq_01_GV(char **out_filename, int *auxfiles, int *retval,
                         char **set_name)
{
    gs_varset_type *x = gs_fdta(*set_name);
    if (!gs_lsq_rank_ready(x))
    {
        *retval = -2;
        return;
    }
    gs_esteff_lsq_01(x, out_filename, auxfiles, retval);
}

void gs_lm_ftest_01(gs_varset_type *gv, char **out_filename, int *auxfiles,
                    int *retval)

{
    int Srv;
    int *Sretval = &Srv;
    gs_UNSIGNED i, j;
    FILE *outfp;

    double *Q_val = malloc(gv->dimu * sizeof(double));
    double *F_val = malloc(gv->dimu * sizeof(double));
    gv->p = realloc(gv->p, gv->dimu * sizeof(double));
    if (NULL == Q_val || NULL == F_val || NULL == gv->p)
        ERR_M;

    double *k = malloc(gv->dimu * sizeof(double));
    if (NULL == k)
        ERR_M;

    gs_mme_restcoeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    double *ixpx = malloc(gv->dimu * gv->dimu * sizeof(double));
    if (NULL == ixpx)
        ERR_M;

    st_dgetri_inv(gv->A, ixpx, gv->dimA);

    free(gv->A);
    gv->A = NULL;

    gs_mme_restrhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_UNSIGNED dfe;
    double kb, kxxk, yy, bxy, sse;
    double *kixtx = malloc(gv->dimA * sizeof(double));

    st_dgemm("T", "N", gv->NoInd, 1, gv->NoInd, 1, gv->y, gv->y, &yy);
    st_dgemm("T", "N", gv->dimu, 1, gv->dimu, 1, gv->u, gv->b, &bxy);
    sse = yy - bxy;
    dfe = gv->NoInd - (gv->NoMar + 1);

    snprintf(gv->msg, sizeof(gv->msg), "SSE = %9f, DFE = %lu", sse, dfe);
    gs_info(gv, 1, gv->msg);

    for (i = 0; i < gv->dimu; i++)
    {
        for (j = 0; j < gv->dimu; j++)
            k[j] = 0;
        k[i] = 1;

        st_dgemm("T", "N", gv->dimA, 1, gv->dimA, 1, k, gv->u, &kb);
        st_dgemm("T", "N", gv->dimA, 1, gv->dimA, gv->dimA, k, ixpx, kixtx);
        st_dgemm("T", "N", gv->dimA, 1, gv->dimA, 1, kixtx, k, &kxxk);
        Q_val[i] = kb * (1 / kxxk) * kb;
        F_val[i] = Q_val[i] / (sse / (double)dfe);
        gv->p[i] = st_fprob(F_val[i], 1, dfe);
    }

    if (1 == *auxfiles)
    {

        char sig[5];

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "effect ssh f pvalue sig\n");
        for (i = 0; i < gv->dimu; i++)
        {
            if (gv->p[i] < 0.001)
                strcpy(sig, "***");
            else if (gv->p[i] < 0.01)
                strcpy(sig, "** ");
            else if (gv->p[i] < 0.05)
                strcpy(sig, "*  ");
            else if (gv->p[i] < 0.10)
                strcpy(sig, ".  ");
            else
                strcpy(sig, "");

            fprintf(outfp, "%s.%i %f %f %f %f \"%s\"\n", gv->EffNme[i],
                    gv->EffAll[i], gv->u[i], Q_val[i], F_val[i], gv->p[i], sig);
        }

        fclose(outfp);
    }

    free(Q_val);
    free(F_val);
    free(k);
    free(ixpx);
    free(kixtx);

    *retval = 0;
    return;
}

void gs_esteff_rmlc_01(gs_varset_type *gv, int *maxiter, double *precision,
                       double *hsq, char **out_filename, int *auxfiles,
                       int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_build_Z_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_rmlc_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_esteff_rmlc_01_GV(int *maxiter, double *precision, double *hsq,
                          char **out_filename, int *auxfiles, int *retval,
                          char **set_name)
{
    gs_esteff_rmlc_01(gs_fdta(*set_name), maxiter, precision, hsq, out_filename,
                      auxfiles, retval);
}

void gs_esteff_rmlv_01(gs_varset_type *gv, int *maxiter, double *precision,
                       double *hsq, char **out_filename, int *auxfiles,
                       int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_build_Z_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_rmlv_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_esteff_rmlv_01_GV(int *maxiter, double *precision, double *hsq,
                          char **out_filename, int *auxfiles, int *retval,
                          char **set_name)
{
    gs_esteff_rmlv_01(gs_fdta(*set_name), maxiter, precision, hsq, out_filename,
                      auxfiles, retval);
}

void gs_permtest_01(gs_varset_type *gv, char *scheme, int *nperm, int *maxiter,
                    double *precision, double *hsq, char **out_filename,
                    int *auxfiles, int *retval)
{
    int ntr = 1;
    int parallel_error = 0;
    int sh_info_level;
    size_t dimu;
    size_t eff_count;

    double *eff_st = NULL;
    double *p_seq = NULL;
    double *p_leq = NULL;
    gs_varset_type **thread_es = NULL;

    *retval = -2;

    if (gv->NoInd > (gs_UNSIGNED)INT_MAX)
    {
        gs_info(gv, -2, "More individuals than INT_MAX");
        return;
    }

    if ((0 == gv->NoInd) || (0 == gv->NoMar) || !st_allele_present(gv))
    {
        gs_info(gv, -2, "No marker data available");
        return;
    }

    if (*nperm < 1)
    {
        gs_info(gv, -2, "Number of permutations must be positive");
        return;
    }

    if ((0 != strcmp(scheme, "lsq")) && (0 != strcmp(scheme, "rmlc")) &&
        (0 != strcmp(scheme, "rmlv")))
    {
        gs_info(gv, -2, "Unknown permutation-test scheme");
        return;
    }

    if ((0 == gv->dimu) || (NULL == gv->u))
    {
        gs_info(gv, -2, "No estimated effects available");
        return;
    }

    dimu = (size_t)gv->dimu;
    if ((gs_UNSIGNED)dimu != gv->dimu)
    {
        gs_info(gv, -2, "Permutation result size is too large");
        return;
    }

    if ((size_t)*nperm > ((size_t)-1) / dimu)
    {
        gs_info(gv, -2, "Permutation result size is too large");
        return;
    }

    eff_count = (size_t)*nperm * dimu;

    if (eff_count > ((size_t)-1) / sizeof(*eff_st))
    {
        gs_info(gv, -2, "Permutation result size is too large");
        return;
    }

    eff_st = calloc(eff_count, sizeof(*eff_st));
    if (NULL == eff_st)
    {
        gs_info(gv, -2, "Memory allocation error");
        return;
    }

    ntr = 1;

    if ((size_t)ntr > ((size_t)-1) / sizeof(*thread_es))
    {
        free(eff_st);
        gs_info(gv, -2, "Memory allocation error");
        return;
    }

    thread_es = malloc((size_t)ntr * sizeof(*thread_es));
    if (NULL == thread_es)
    {
        free(eff_st);
        gs_info(gv, -2, "Memory allocation error");
        return;
    }

    gs_get_info_level(sh, &sh_info_level);
    gs_set_info_level(sh, imz);

    for (int t = 0; t < ntr; t++)
    {
        int Srv = 0;
        int *Sretval = &Srv;
        char esname[24];

        snprintf(esname, sizeof(esname), "dta_t%03i", t);
        thread_es[t] = gs_fdta(esname);

        if (NULL == thread_es[t])
        {
            gs_set_info_level(sh, &sh_info_level);
            free(thread_es);
            free(eff_st);
            gs_info(gv, -2, "Memory allocation error");
            return;
        }

        thread_es[t]->info_level = -2;

        gs_copy_marker_data_01(thread_es[t], gv, Sretval);

        if (-2 == *Sretval || (st_allele_bits(thread_es[t]) != 16 &&
                               st_allele_pack_try(thread_es[t], 16)))
        {
            gs_set_info_level(sh, &sh_info_level);
            free(thread_es);
            free(eff_st);
            return;
        }
    }

    gs_set_info_level(sh, &sh_info_level);

#ifdef MT
#pragma omp parallel num_threads(ntr) if (ntr > 1) reduction(| : parallel_error)
#endif
    {
        gs_UNSIGNED i, j, idx1, idx2;
        int Srv;
        int *Sretval = &Srv;
        int t = omp_get_thread_num();
        gs_varset_type *es = NULL;
        int *p = NULL;
        int local_error = 0;

        if (t >= ntr)
        {
            local_error = 1;
            parallel_error = 1;
        }
        else
        {
            es = thread_es[t];
        }

        if (gv->NoInd > (gs_UNSIGNED)(((size_t)-1) / sizeof(*p)))
        {
            local_error = 1;
            parallel_error = 1;
        }
        else if (!local_error)
        {
            p = malloc((size_t)gv->NoInd * sizeof(*p));
            if (NULL == p)
            {
                local_error = 1;
                parallel_error = 1;
            }
        }

#ifdef MT

        omp_set_schedule(st_rng_is_seeded() ? omp_sched_static
                                            : omp_sched_dynamic,
                         st_rng_is_seeded() ? 0 : 1);
#pragma omp for schedule(runtime)
#endif
        for (int r = 0; r < *nperm; r++)
        {
            gs_UNSIGNED no_unsucc_rep = 0;
            int success_rep = 0;

            if (local_error)
                continue;

            while ((0 == success_rep) && (no_unsucc_rep < 1000))
            {
                success_rep = 1;

                for (i = 1; i < gv->dimu; i++)
                {

                    p[0] = 0;
                    for (int ii = 1; ii < (int)gv->NoInd; ii++)
                    {
                        int jj = st_rand_index(ii + 1);
                        p[ii] = p[jj];
                        p[jj] = ii;
                    }

                    for (j = 0; j < gv->NoInd; j++)
                    {
                        idx1 = (i - 1) * gv->NoInd + j;
                        idx2 = (i - 1) * gv->NoInd + p[j];
                        es->allele_raw1[idx1] = st_allele_get(gv, idx2, 0);
                        es->allele_raw2[idx1] = st_allele_get(gv, idx2, 1);
                    }

                    *Sretval = 0;

                    if (0 == strcmp(scheme, "lsq"))
                        gs_esteff_lsq_01(es, nixp, neinp, Sretval);
                    else if (0 == strcmp(scheme, "rmlc"))
                        gs_esteff_rmlc_01(es, maxiter, precision, hsq, nixp,
                                          neinp, Sretval);
                    else
                        gs_esteff_rmlv_01(es, maxiter, precision, hsq, nixp,
                                          neinp, Sretval);

                    if (-2 == *Sretval)
                    {
                        success_rep = 0;
                        break;
                    }

                    eff_st[(size_t)i * (size_t)*nperm + (size_t)r] = es->u[i];
                }

                if (0 == success_rep)
                    no_unsucc_rep++;
            }

            if (0 == success_rep)
                parallel_error = 1;
        }

        free(p);
    }

    for (int t = 0; t < ntr; t++)
        if (ST_ALLELE_AUTO_PACK && st_allele_pack_try(thread_es[t], 0))
            parallel_error = 1;
    free(thread_es);
    thread_es = NULL;

    if (parallel_error)
    {
        free(eff_st);
        gs_info(gv, -2, "Permutation test failed in parallel calculation");
        return;
    }

#ifdef MT
#pragma omp parallel for
#endif
    for (gs_UNSIGNED i = 0; i < gv->dimu; i++)
        qsort(&(eff_st[(size_t)i * (size_t)*nperm]), (size_t)*nperm,
              sizeof(double), gs_dbl_cmp);

    if (dimu > ((size_t)-1) / sizeof(double))
    {
        free(eff_st);
        gs_info(gv, -2, "Permutation result size is too large");
        return;
    }

    p_seq = malloc(dimu * sizeof(*p_seq));
    p_leq = malloc(dimu * sizeof(*p_leq));

    if ((NULL == p_seq) || (NULL == p_leq))
    {
        free(p_seq);
        free(p_leq);
        free(eff_st);
        gs_info(gv, -2, "Memory allocation error");
        return;
    }

    {
        double *tmp_p = realloc(gv->p, dimu * sizeof(*gv->p));

        if (NULL == tmp_p)
        {
            free(p_seq);
            free(p_leq);
            free(eff_st);
            gs_info(gv, -2, "Memory allocation error");
            return;
        }

        gv->p = tmp_p;
    }

#ifdef MT
#pragma omp parallel for
#endif
    for (gs_UNSIGNED i = 0; i < gv->dimu; i++)
    {
        gs_UNSIGNED r;

        for (r = 0; r < (gs_UNSIGNED)*nperm; r++)
            if (gv->u[i] < eff_st[(size_t)i * (size_t)*nperm + (size_t)r])
                break;
        p_seq[i] = (double)r;

        for (r = 0; r < (gs_UNSIGNED)*nperm; r++)
            if (gv->u[i] > eff_st[(size_t)i * (size_t)*nperm +
                                  (size_t)(*nperm - 1) - (size_t)r])
                break;
        p_leq[i] = (double)r;

        if (p_seq[i] < p_leq[i])
            gv->p[i] = p_seq[i] / (double)(*nperm);
        else
            gv->p[i] = p_leq[i] / (double)(*nperm);
    }

    if (1 == *auxfiles)
    {
        FILE *outfp;
        char sig[5];

        outfp = fopen(*out_filename, "w");
        if (NULL == outfp)
        {
            free(p_seq);
            free(p_leq);
            free(eff_st);
            gs_info(gv, -2, "File open error");
            return;
        }

        fprintf(outfp, "effect  p.seq p.leq pvalue sig\n");
        for (gs_UNSIGNED i = 0; i < gv->dimu; i++)
        {
            if (gv->p[i] < 0.001)
                strcpy(sig, "***");
            else if (gv->p[i] < 0.01)
                strcpy(sig, "** ");
            else if (gv->p[i] < 0.05)
                strcpy(sig, "*  ");
            else if (gv->p[i] < 0.10)
                strcpy(sig, ".  ");
            else
                strcpy(sig, "");

            fprintf(outfp, "%s.%i %f %f %f %f \"%s\"\n", gv->EffNme[i],
                    gv->EffAll[i], gv->u[i], p_seq[i], p_leq[i], gv->p[i], sig);
        }

        fclose(outfp);
    }

    free(p_seq);
    free(p_leq);
    free(eff_st);

    *retval = 0;
}

void gs_testeff_lsq_01(gs_varset_type *gv, int *nperm, char **out_filename,
                       int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_esteff_lsq_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    if (0 == *nperm)
    {
        gs_lm_ftest_01(gv, out_filename, auxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }
    else
    {
        char scheme[] = "lsq";
        int mi = 0;
        int *maxiter = &mi;
        double pr = 0;
        double *precision = &pr;
        double _hsq = 0.9;
        double *hsq = &_hsq;
        gs_permtest_01(gv, scheme, nperm, maxiter, precision, hsq, out_filename,
                       auxfiles, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }
    }

    *retval = 0;
    return;
}

void gs_testeff_lsq_01_GV(int *nperm, char **out_filename, int *auxfiles,
                          int *retval, char **set_name)
{
    gs_varset_type *x = gs_fdta(*set_name);
    if (!gs_lsq_rank_ready(x))
    {
        *retval = -2;
        return;
    }
    gs_testeff_lsq_01(x, nperm, out_filename, auxfiles, retval);
}

void gs_testeff_rmlc_01(gs_varset_type *gv, int *nperm, int *maxiter,
                        double *precision, double *hsq, char **out_filename,
                        int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_esteff_rmlc_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    char scheme[] = "rmlc";
    gs_permtest_01(gv, scheme, nperm, maxiter, precision, hsq, out_filename,
                   auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_testeff_rmlc_01_GV(int *nperm, int *maxiter, double *precision,
                           double *hsq, char **out_filename, int *auxfiles,
                           int *retval, char **set_name)
{
    gs_testeff_rmlc_01(gs_fdta(*set_name), nperm, maxiter, precision, hsq,
                       out_filename, auxfiles, retval);
}

void gs_testeff_rmlv_01(gs_varset_type *gv, int *nperm, int *maxiter,
                        double *precision, double *hsq, char **out_filename,
                        int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_esteff_rmlv_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    char scheme[] = "rmlv";
    gs_permtest_01(gv, scheme, nperm, maxiter, precision, hsq, out_filename,
                   auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_testeff_rmlv_01_GV(int *nperm, int *maxiter, double *precision,
                           double *hsq, char **out_filename, int *auxfiles,
                           int *retval, char **set_name)
{
    gs_testeff_rmlv_01(gs_fdta(*set_name), nperm, maxiter, precision, hsq,
                       out_filename, auxfiles, retval);
}

void gs_esteff_rmlr_01(gs_varset_type *gv, int *maxiter, double *precision,
                       double *hsq, char **out_filename, int *auxfiles,
                       int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_build_Z_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_rmlr_01(gv, maxiter, precision, hsq, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_esteff_rmlr_01_GV(int *maxiter, double *precision, double *hsq,
                          char **out_filename, int *auxfiles, int *retval,
                          char **set_name)
{
    gs_esteff_rmlr_01(gs_fdta(*set_name), maxiter, precision, hsq, out_filename,
                      auxfiles, retval);
}

void gs_esteff_rmla_01(gs_varset_type *gv, double *alpha, int *maxiter,
                       double *precision, double *hsq, char **out_filename,
                       int *auxfiles, int *retval)
{
    int Srv;
    int *Sretval = &Srv;

    gs_build_Z_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_rmla_01(gv, alpha, maxiter, precision, hsq, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_coeff_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mme_solve_01(gv, out_filename, auxfiles, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    *retval = 0;
    return;
}

void gs_esteff_rmla_01_GV(double *alpha, int *maxiter, double *precision,
                          double *hsq, char **out_filename, int *auxfiles,
                          int *retval, char **set_name)
{
    gs_esteff_rmla_01(gs_fdta(*set_name), alpha, maxiter, precision, hsq,
                      out_filename, auxfiles, retval);
}

void gs_lambda_emstep_tc_01(gs_varset_type *gv, char **out_filename,
                            int *auxfiles, int *retval)
{

    gs_UNSIGNED i, j, m;

    if (NULL == gv->ladi)
    {
        gs_info(gv, -2, "No lambda defined");
        *retval = -2;
        return;
    }

    int Srv;
    int *Sretval = &Srv;

    gs_mmet_coeff_02(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_mmet_rhs_01(gv, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gv->a = realloc(gv->a, gv->dimA * sizeof(double));
    if (NULL == gv->a)
        ERR_M;

    for (i = 0; i < gv->dimA; i++)
        gv->a[i] = gv->b[i];

    {
        double d;

        gv->u0 = realloc(gv->u0, gv->dimA * sizeof(gs_UNSIGNED));

        st_dgetrf(gv->A, gv->dimA, gv->u0, &d);

        st_dgetrs(gv->A, gv->dimA, gv->u0, gv->a);

        gv->t4 = realloc(gv->t4, gv->dimA * gv->dimA * sizeof(double));

        st_dgetri(gv->A, gv->dimA, gv->u0, gv->t4);
    }

    gv->sigma = realloc(gv->sigma, gv->NoMar * sizeof(double));
    if (NULL == gv->sigma)
        ERR_M;

    gv->t1 = realloc(gv->t1, sizeof(double));

    gv->t2 = realloc(gv->t2, ((gv->NoInd > gv->NoMar) ? gv->NoInd : gv->NoMar) *
                                 sizeof(double));

    gv->t3 = realloc(gv->t3, ((gv->NoInd > gv->NoMar) ? gv->NoInd : gv->NoMar) *
                                 sizeof(double));

    for (i = 0; i < gv->NoInd; i++)
        gv->t2[i] = 1;

    st_dgemm("T", "N", gv->NoInd, 1, gv->NoInd, 1, gv->y, gv->y, gv->sigma);

    st_dgemm("T", "N", gv->NoInd, 1, gv->NoInd, 1, gv->t2, gv->y, gv->t1);

    st_dgemm("T", "N", gv->NoInd, 1, gv->NoInd, 1, 1 + (gv->a), gv->y, gv->t3);

    gv->sigma[0] -= gv->a[0] * gv->t1[0];
    gv->sigma[0] -= gv->t3[0];
    gv->sigma[0] = gv->sigma[0] / (double)gv->NoInd;
    gv->sigma_e = gv->sigma[0];

    st_dgemm("T", "N", gv->NoInd, 1, gv->NoInd, gv->NoInd, 1 + (gv->a), gv->GI,
             gv->t3);

    st_dgemm("N", "T", 1, gv->NoInd, 1, gv->NoInd, gv->t3, 1 + (gv->a), gv->t2);

    gv->t5 = realloc(gv->t5, gv->NoInd * gv->NoInd * sizeof(double));
    for (j = 1; j < gv->dimA; j++)
        memcpy(&gv->t5[(j - 1) * gv->NoInd], &gv->t4[j * gv->dimA + 1],
               (gv->NoInd) * sizeof(double));

    gv->t6 = realloc(gv->t6, gv->NoInd * gv->NoInd * sizeof(double));
    st_dgemm("N", "N", gv->NoInd, gv->NoInd, gv->NoInd, gv->NoInd, gv->GI,
             gv->t5, gv->t6);

    double tr = 0;
    for (i = 0; i < gv->NoInd; i++)
        tr += gv->t6[gv->NoInd * i + i];

    gv->sigma[1] = (gv->t2[0] + gv->sigma[0] * tr) / (double)gv->NoInd;

    for (i = 0; i < gv->NoMar; i++)
        gv->ladi[i] = gv->sigma[0] / gv->sigma[1];

    if (1 == *auxfiles)
    {
        FILE *outfp;

        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda var \n");
        fprintf(outfp, "err %f %f \n", 0.0, gv->sigma[0]);
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f %f\n", gv->MarName[m], gv->ladi[m],
                    gv->sigma[1 + m]);

        fclose(outfp);
    }

    *retval = 0;
    return;
}

void gs_lambda_emstep_tc_01_GV(char **out_filename, int *auxfiles, int *retval,
                               char **set_name)

{
    gs_lambda_emstep_tc_01(gs_fdta(*set_name), out_filename, auxfiles, retval);
}

void gs_lambda_rmlct_01(gs_varset_type *gv, int *maxiter, double *precision,
                        double *hsq, char **out_filename, int *auxfiles,
                        int *retval)
{

    gs_UNSIGNED m;
    int i;
    double qu, maxqu = 0;
    int Srv;
    int *Sretval = &Srv;

    if (NULL == gv->Z)
    {
        gs_info(gv, -2, "No Z matrix available");
        *retval = -2;
        return;
    }

    double vc_save[2] = {0, 0};

    gv->ladi = realloc(gv->ladi, gv->NoMar * sizeof(double));
    if (NULL == gv->ladi)
        ERR_M;
    for (m = 0; m < gv->NoMar; m++)
        gv->ladi[m] = (1 / *hsq - 1) * (double)gv->NoInd;

    gs_build_GI_01(gv, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    for (i = 0; i < *maxiter; i++)
    {

        if (i > 0)
            memcpy(vc_save, gv->sigma, 2 * sizeof(double));

        gs_lambda_emstep_tc_01(gv, nixp, neinp, Sretval);
        if (-2 == *Sretval)
        {
            *retval = -2;
            return;
        }

        maxqu = 0;
        for (m = 0; m < 2; m++)
        {
            qu = fabs((vc_save[m] - gv->sigma[m]) / gv->sigma[m]);
            if (qu > maxqu)
                maxqu = qu;
        }
        if (maxqu < *precision)
            break;
    }

    snprintf(gv->msg, sizeof(gv->msg),
             "Number of EM iterations: %i, precision: %f", i, maxqu);
    gs_info(gv, 1, gv->msg);
    snprintf(gv->msg, sizeof(gv->msg),
             "Variance components: var(m) = %f, var(e) = %5f", gv->sigma[1],
             gv->sigma[0]);
    gs_info(gv, 1, gv->msg);

    if (1 == *auxfiles)
    {
        FILE *outfp;
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_lambda_rmlct_01_GV(int *maxiter, double *precision, double *hsq,
                           char **out_filename, int *auxfiles, int *retval,
                           char **set_name)
{
    gs_lambda_rmlct_01(gs_fdta(*set_name), maxiter, precision, hsq,
                       out_filename, auxfiles, retval);
}

void gs_lambda_rmlat_01(gs_varset_type *gv, double *alpha, int *maxiter,
                        double *precision, double *hsq, char **out_filename,
                        int *auxfiles, int *retval)
{
    gs_UNSIGNED m;
    FILE *outfp;
    int Srv;
    int *Sretval = &Srv;

    gs_single_marker_aov_01(gv, alpha, nixp, neinp, Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    gs_lambda_rmlct_01(gv, maxiter, precision, hsq, out_filename, auxfiles,
                       Sretval);
    if (-2 == *Sretval)
    {
        *retval = -2;
        return;
    }

    long double Svar = 0.0;
    for (m = 0; m < gv->NoMar; m++)
        Svar += gv->SMvar[m];

    Svar = Svar / (long double)gv->NoMar;

    for (m = 0; m < gv->NoMar; m++)
        gv->ladi[m] = (Svar / gv->SMvar[m]) * gv->ladi[m];

    if (1 == *auxfiles)
    {
        if ((outfp = fopen(*out_filename, "w")) == NULL)
            ERR_F;

        fprintf(outfp, "lambda \n");
        for (m = 0; m < gv->NoMar; m++)
            fprintf(outfp, "%s %f\n", gv->MarName[m], gv->ladi[m]);

        fclose(outfp);
    }
    *retval = 0;
    return;
}

void gs_lambda_rmlat_01_GV(double *alpha, int *maxiter, double *precision,
                           double *hsq, char **out_filename, int *auxfiles,
                           int *retval, char **set_name)
{
    gs_lambda_rmlat_01(gs_fdta(*set_name), alpha, maxiter, precision, hsq,
                       out_filename, auxfiles, retval);
}

static int st_sim_marker_data_valid(gs_varset_type *gv)
{
    if ((0 == gv->NoInd) || (0 == gv->NoMar) || !st_allele_present(gv) ||
        (NULL == gv->MarName))
    {
        gs_info(gv, -2, "No marker data available");
        return 0;
    }

    if (!gs_marker_genotypes_valid(gv))
    {
        gs_info(gv, -2, "Invalid marker genotype representation");
        return 0;
    }

    return 1;
}

static int st_sim_data_valid(gs_varset_type *gv)
{
    gs_UNSIGNED *storage;

    if (!st_sim_marker_data_valid(gv))
        return 0;

    storage = gs_map_storage_index(gv, 1, 1);

    if (NULL == storage)
        return 0;

    free(storage);

    return 1;
}

void st_set_sim_gp_provide_data(int *n, double **chrlen, int *retval,
                                char **set_name)
{
    gs_varset_type *gv = gs_fdta(*set_name);
    gs_UNSIGNED lne;
    double *maxpos = NULL;
    unsigned char *seen = NULL;
    int maxchrom = 0;
    int chrom;

    *retval = -2;
    *n = 0;
    *chrlen = NULL;

    if (!st_sim_data_valid(gv))
        return;

    for (lne = 0; lne < gv->NLoc; lne++)
    {
        chrom = gv->map[lne]->chrom;

        if ((chrom < 1) || !isfinite(gv->map[lne]->pos) ||
            (gv->map[lne]->pos < 0))
        {
            gs_info(gv, -2, "Invalid chromosome or position in linkage map");
            return;
        }

        if (chrom > maxchrom)
            maxchrom = chrom;
    }

    if (maxchrom < 1)
    {
        gs_info(gv, -2, "No valid chromosomes in linkage map");
        return;
    }

    if ((size_t)maxchrom > ((size_t)-1) / sizeof(double))
    {
        gs_info(gv, -2, "Memory allocation error");
        return;
    }

    *chrlen = malloc((size_t)maxchrom * sizeof(double));
    maxpos = calloc((size_t)maxchrom, sizeof(double));
    seen = calloc((size_t)maxchrom, sizeof(unsigned char));

    if ((NULL == *chrlen) || (NULL == maxpos) || (NULL == seen))
    {
        gs_info(gv, -2, "Memory allocation error");
        goto fail;
    }

    for (lne = 0; lne < gv->NLoc; lne++)
    {
        chrom = gv->map[lne]->chrom - 1;
        seen[chrom] = 1;

        if (gv->map[lne]->pos > maxpos[chrom])
            maxpos[chrom] = gv->map[lne]->pos;
    }

    for (chrom = 0; chrom < maxchrom; chrom++)
    {
        if (!seen[chrom])
        {
            gs_info(gv, -2,
                    "Chromosome numbering in linkage map is not contiguous");
            goto fail;
        }

        (*chrlen)[chrom] = 0.01 + maxpos[chrom] / 100.0;

        if (!isfinite((*chrlen)[chrom]) || ((*chrlen)[chrom] <= 0))
        {
            gs_info(gv, -2,
                    "Invalid chromosome length derived from linkage map");
            goto fail;
        }
    }

    *n = maxchrom;

    free(maxpos);
    free(seen);

    *retval = 0;
    return;

fail:

    free(*chrlen);
    free(maxpos);
    free(seen);

    *chrlen = NULL;
    *n = 0;
}

void st_set_sim_mp_provide_data(gs_UNSIGNED *NLoc, gs_mappoint_type ***map,
                                int *retval, char **set_name)
{

    gs_varset_type *gv = gs_fdta(*set_name);

    if (NULL == gv->map)
    {
        snprintf(gv->msg, sizeof(gv->msg), "No map data available");
        gs_info(gv, -2, gv->msg);
        *retval = -2;
        return;
    }

    {
        gs_UNSIGNED *storage = gs_map_storage_index(gv, 0, 0);

        if (NULL == storage)
        {
            *retval = -2;
            return;
        }
        free(storage);
    }

    *NLoc = gv->NLoc;
    *map = gv->map;

    *retval = 0;
    return;
}

static int st_sim_effect_data_state(gs_varset_type *gv)
{
    if ((0 == gv->NoEff) || (0 == gv->dimu) || (NULL == gv->u))
        return 0;

    if ((0 == gv->NoMar) || (NULL == gv->MarName) ||
        !gs_restrict_marker_effects_valid(gv))
        return -1;

    return 1;
}

void st_set_sim_ef_provide_data(gs_UNSIGNED *NoEff, char (**EffNme)[gs_MAXNME],
                                int **EffAll, int **CmpAll, double **u,
                                int *retval, char **set_name)
{
    gs_varset_type *gv = gs_fdta(*set_name);
    int state;

    *retval = -2;
    *NoEff = 0;
    *EffNme = NULL;
    *EffAll = NULL;
    *CmpAll = NULL;
    *u = NULL;

    state = st_sim_effect_data_state(gv);

    if (state < 0)
    {
        gs_info(gv, -2, "Invalid marker effect data");
        return;
    }

    if (0 == state)
    {
        *retval = 0;
        return;
    }

    *NoEff = gv->NoEff;
    *EffNme = gv->EffNme;
    *EffAll = gv->EffAll;
    *CmpAll = gv->CmpAll;
    *u = gv->u;

    *retval = 0;
    return;
}

void st_set_sim_pp_provide_data(gs_UNSIGNED *NoInd, gs_UNSIGNED *NoMar,
                                char ***IndName, char ***Pedigree,
                                char (**MarName)[gs_MAXNME],
                                GS_ALLELE_TYPE **All1, GS_ALLELE_TYPE **All2,
                                int *retval, char **set_name)
{
    gs_varset_type *gv = gs_fdta(*set_name);

    *retval = -2;
    *NoInd = 0;
    *NoMar = 0;
    *IndName = *Pedigree = NULL;
    *MarName = NULL;
    *All1 = NULL;
    *All2 = NULL;

    if (!st_sim_marker_data_valid(gv))
        return;

    *IndName = gv->IndName;
    *Pedigree = gv->Pedigree;
    *NoInd = gv->NoInd;
    *NoMar = gv->NoMar;
    *MarName = gv->MarName;
    if (st_allele_bits(gv) == 16)
    {
        *All1 = gv->allele_raw1;
        *All2 = gv->allele_raw2;
    }
    else
    {
        size_t n;
        if (gs_size_mul((size_t)gv->NoMar, (size_t)gv->NoInd, &n))
            return;
        *All1 = (GS_ALLELE_TYPE *)st_work_alloc(n, sizeof(GS_ALLELE_TYPE));
        *All2 = (GS_ALLELE_TYPE *)st_work_alloc(n, sizeof(GS_ALLELE_TYPE));
        st_allele_copy(gv, 0, 0, n, *All1);
        st_allele_copy(gv, 1, 0, n, *All2);
    }

    *retval = 0;
    return;
}

void gs_set_performance_data_01(gs_varset_type *gv, int *n, double *y,
                                int *retval)
{

    if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
        gs_info(gv, -2,
                "No marker data available, run 'read.marker.data' first");
        *retval = -2;
        return;
    }

    if (*n != (int)gv->NoInd)
    {
        gs_info(gv, -2, "Incorrect numer of observations");
        *retval = -2;
        return;
    }

    gv->y = realloc(gv->y, gv->NoInd * sizeof(double));
    if (NULL == gv->y)
        ERR_M;

    st_dcopy(gv->NoInd, gv->y, y);

    *retval = 0;
    return;
}

void gs_set_performance_data_01_GV(int *n, double *y, int *retval,
                                   char **set_name)
{
    gs_set_performance_data_01(gs_fdta(*set_name), n, y, retval);
}

typedef enum
{
    K_STD = 0,
    K_VANRADEN = 1,
    K_SOMMER = 2
} KMethod;

typedef enum
{
    CODING_012 = 0,
    CODING_M101 = 1,
    CODING_001 = 2,
    CODING_OTHER = 3
} ZCoding;

typedef enum
{
    TRANS_NONE = 0,
    TRANS_ADD1 = 1,
    TRANS_MUL2 = 2
} ZTransform;

static int is_in_set_012(double v)
{
    return (v == 0.0 || v == 1.0 || v == 2.0);
}

static int is_in_set_m101(double v)
{
    return (v == -1.0 || v == 0.0 || v == 1.0);
}

static int is_in_set_001(double v)
{
    return (v == 0.0 || v == 0.5 || v == 1.0);
}

static int Z_is_compatible_with(const double *Zrm, int nG, int m,
                                ZCoding coding)
{
    const size_t nm = (size_t)nG * (size_t)m;
    for (size_t t = 0; t < nm; ++t)
    {
        double v = Zrm[t];
        if (!(v == v))
            continue;
        int ok = 0;
        if (coding == CODING_012)
            ok = is_in_set_012(v);
        if (coding == CODING_M101)
            ok = is_in_set_m101(v);
        if (coding == CODING_001)
            ok = is_in_set_001(v);
        if (!ok)
            return 0;
    }
    return 1;
}

static int detect_coding_strict_or_ambiguous(const double *Zrm, int nG, int m,
                                             ZCoding *out)
{
    int ok012 = Z_is_compatible_with(Zrm, nG, m, CODING_012);
    int okm101 = Z_is_compatible_with(Zrm, nG, m, CODING_M101);
    int ok001 = Z_is_compatible_with(Zrm, nG, m, CODING_001);
    int n_ok = ok012 + okm101 + ok001;

    if (n_ok == 0)
        return -2001;
    if (n_ok > 1)
        return -2004;

    if (ok012)
        *out = CODING_012;
    else if (okm101)
        *out = CODING_M101;
    else
        *out = CODING_001;
    return 0;
}

static int parse_coding_text(const char *s, int *is_detect, ZCoding *out)
{
    if (!s)
        return -2005;
    if (strcmp(s, "detect") == 0)
    {
        *is_detect = 1;
        return 0;
    }
    *is_detect = 0;
    if (strcmp(s, "0 1 2") == 0)
    {
        *out = CODING_012;
        return 0;
    }
    if (strcmp(s, "-1 0 1") == 0)
    {
        *out = CODING_M101;
        return 0;
    }
    if (strcmp(s, "0 0.5 1") == 0)
    {
        *out = CODING_001;
        return 0;
    }
    return -2005;
}

static int build_K_from_Z_text(const double *Zrm, int nG, int m, KMethod method,
                               const char *coding, double *Krm_out,
                               double *out_z_mean, double *out_z_sd,
                               ZCoding *out_coding_used,
                               ZTransform *out_transform)
{
    const size_t nm = (size_t)nG * (size_t)m;

    int is_detect = 0;
    ZCoding coding_forced = CODING_OTHER;
    int prc = parse_coding_text(coding, &is_detect, &coding_forced);
    if (prc != 0)
        return -2005;

    ZCoding coding_used = CODING_OTHER;
    if (is_detect)
    {
        int drc = detect_coding_strict_or_ambiguous(Zrm, nG, m, &coding_used);
        if (drc != 0)
            return drc;
    }
    else
    {

        if (!Z_is_compatible_with(Zrm, nG, m, coding_forced))
            return -2006;
        coding_used = coding_forced;
    }
    if (out_coding_used)
        *out_coding_used = coding_used;

    double *Zwork = (double *)malloc(nm * sizeof(double));
    if (!Zwork)
        return -1010;
    memcpy(Zwork, Zrm, nm * sizeof(double));

    if (method == K_STD)
    {
    }

    ZTransform tr = TRANS_NONE;

    if (coding_used == CODING_M101)
    {

        tr = TRANS_ADD1;
        for (size_t t = 0; t < nm; ++t)
            if (Zwork[t] == Zwork[t])
                Zwork[t] += 1.0;
    }
    else if (coding_used == CODING_001)
    {

        tr = TRANS_MUL2;
        for (size_t t = 0; t < nm; ++t)
            if (Zwork[t] == Zwork[t])
                Zwork[t] *= 2.0;
    }
    else
    {

        tr = TRANS_NONE;
    }
    if (out_transform)
        *out_transform = tr;

    if (out_z_mean || out_z_sd)
    {
        double sum = 0.0, sumsq = 0.0;
        size_t cnt = 0;
        for (size_t t = 0; t < nm; ++t)
        {
            double v = Zwork[t];
            if (!(v == v))
                continue;
            sum += v;
            sumsq += v * v;
            cnt++;
        }
        double mean = (cnt > 0) ? (sum / (double)cnt) : NAN;
        double var = (cnt > 1) ? fmax(0.0, (sumsq - (double)cnt * mean * mean) /
                                               (double)(cnt - 1))
                               : NAN;
        double sd = (cnt > 1) ? sqrt(var) : NAN;
        if (out_z_mean)
            *out_z_mean = mean;
        if (out_z_sd)
            *out_z_sd = sd;
    }

    for (size_t t = 0; t < nm; ++t)
    {
        double v = Zwork[t];
        if (!(v == v))
            continue;
        if (!is_in_set_012(v))
        {
            free(Zwork);
            return -2002;
        }
    }

    double *p = (double *)malloc((size_t)m * sizeof(double));
    if (!p)
    {
        free(Zwork);
        return -1010;
    }

    int threads = st_num_threads();
    (void)threads;
#ifdef MT
#pragma omp parallel for schedule(static)                                      \
    num_threads(threads) if (nm >= 65536 && threads > 1)
#endif
    for (int j = 0; j < m; ++j)
    {
        double s = 0.0;
        int c = 0;
        for (int i = 0; i < nG; ++i)
        {
            double v = Zwork[(size_t)i * (size_t)m + (size_t)j];
            if (!(v == v))
                continue;
            s += v;
            c++;
        }

        p[j] = (c > 0) ? ((s / (double)c) / 2.0) : 0.5;
    }

    double denom = 0.0;
    for (int j = 0; j < m; ++j)
        denom += p[j] * (1.0 - p[j]);
    denom = 2.0 * denom;

    if (!(denom > 0.0) || !(denom == denom))
    {
        free(p);
        free(Zwork);
        return -2003;
    }

#ifdef MT
#pragma omp parallel for schedule(static)                                      \
    num_threads(threads) if (nm >= 65536 && threads > 1)
#endif
    for (int i = 0; i < nG; ++i)
        for (int j = 0; j < m; ++j)
        {
            double *cell = &Zwork[(size_t)i * (size_t)m + (size_t)j];
            if (*cell == *cell)
                *cell -= 2.0 * p[j];
        }

    st_dsymm_ZZt((gs_UNSIGNED)nG, (gs_UNSIGNED)m, Zwork, Krm_out);

    const double invd = 1.0 / denom;
    for (size_t t = 0; t < (size_t)nG * (size_t)nG; ++t)
        Krm_out[t] *= invd;

    free(p);
    free(Zwork);
    return 0;
}

static int build_K_from_Z(const double *Zrm, int nG, int m, KMethod method,
                          SEXP coding, double *Krm_out, double *out_z_mean,
                          double *out_z_sd, ZCoding *out_coding_used,
                          ZTransform *out_transform)
{
    if (!isString(coding) || XLENGTH(coding) != 1)
        return -2005;
    return build_K_from_Z_text(Zrm, nG, m, method, CHAR(STRING_ELT(coding, 0)),
                               Krm_out, out_z_mean, out_z_sd, out_coding_used,
                               out_transform);
}

double *st_make_A_values(gs_varset_type *x, const char *coding)
{
    if (!x || !x->Z || !x->NoInd || !x->NoMar || x->NoInd > INT_MAX ||
        x->NoMar > INT_MAX ||
        (size_t)x->NoInd > SIZE_MAX / sizeof(double) / (size_t)x->NoMar ||
        (size_t)x->NoInd > SIZE_MAX / sizeof(double) / (size_t)x->NoInd ||
        (R_xlen_t)x->NoInd > R_XLEN_T_MAX / (R_xlen_t)x->NoInd)
        error("Cannot construct oversized native relationship matrix");
    int n = (int)x->NoInd, m = (int)x->NoMar;
    size_t count = (size_t)n * m;
    for (size_t k = 0; k < count; k++)
        if (!R_FINITE(x->Z[k]))
            error("gs.make.A: Z must not contain missing or infinite values");
    double *native = st_work_alloc((size_t)n * n, sizeof(double));
    int rc = build_K_from_Z_text(x->Z, n, m, K_VANRADEN, coding, native, NULL,
                                 NULL, NULL, NULL);
    if (rc)
        error("Cannot construct native relationship matrix (code %d)", rc);
    return native;
}

SEXP st_make_A_native(gs_varset_type *x, SEXP coding)
{
    const char *code = isString(coding) && XLENGTH(coding) == 1
                           ? CHAR(STRING_ELT(coding, 0))
                           : NULL;
    double *native = st_make_A_values(x, code);
    int n = (int)x->NoInd;
    SEXP A = PROTECT(allocMatrix(REALSXP, n, n)),
         labels = PROTECT(allocVector(STRSXP, n)),
         dn = PROTECT(allocVector(VECSXP, 2));
    memcpy(REAL(A), native, (size_t)n * n * sizeof(double));
    for (int i = 0; i < n; i++)
        SET_STRING_ELT(labels, i, mkChar(x->IndName[i]));
    SET_VECTOR_ELT(dn, 0, labels);
    SET_VECTOR_ELT(dn, 1, labels);
    setAttrib(A, R_DimNamesSymbol, dn);
    UNPROTECT(3);
    return A;
}

SEXP gs_make_A(SEXP Z, SEXP coding)
{
    int nprot = 0;
    char *diagnostic = NULL;
    SEXP A = R_NilValue;

    if (!isMatrix(Z) || (TYPEOF(Z) != REALSXP && TYPEOF(Z) != INTSXP))
    {
        diagnostic = "gs.make.A: Z must be a numeric matrix";
        goto cleanup;
    }
    if (TYPEOF(coding) != STRSXP || XLENGTH(coding) != 1 ||
        STRING_ELT(coding, 0) == NA_STRING)
    {
        diagnostic = "gs.make.A: coding must be a single nonmissing string";
        goto cleanup;
    }

    SEXP dim = getAttrib(Z, R_DimSymbol);
    int n = INTEGER(dim)[0], m = INTEGER(dim)[1];
    if (n <= 0 || m <= 0)
    {
        diagnostic = "gs.make.A: Z must have at least one row and one column";
        goto cleanup;
    }
    if ((size_t)n > SIZE_MAX / sizeof(double) / (size_t)m ||
        (size_t)n > SIZE_MAX / sizeof(double) / (size_t)n ||
        (R_xlen_t)n > R_XLEN_T_MAX / (R_xlen_t)n)
    {
        diagnostic = "gs.make.A: matrix dimensions exceed supported storage";
        goto cleanup;
    }

    size_t count = (size_t)n * (size_t)m;
    for (size_t k = 0; k < count; k++)
    {
        double value = TYPEOF(Z) == REALSXP          ? REAL(Z)[k]
                       : INTEGER(Z)[k] == NA_INTEGER ? NA_REAL
                                                     : INTEGER(Z)[k];
        if (!R_FINITE(value))
        {
            diagnostic =
                "gs.make.A: Z must not contain missing or infinite values";
            goto cleanup;
        }
    }

    double *Zrm = (double *)st_work_alloc(count, sizeof(double));
    for (int j = 0; j < m; j++)
    {
        if (j % 1024 == 0)
            R_CheckUserInterrupt();
        for (int i = 0; i < n; i++)
        {
            size_t k = (size_t)i + (size_t)n * (size_t)j;
            Zrm[(size_t)i * (size_t)m + (size_t)j] =
                TYPEOF(Z) == REALSXP ? REAL(Z)[k] : INTEGER(Z)[k];
        }
    }

    double *native_A = st_work_alloc((size_t)n * n, sizeof(double));
    int rc = build_K_from_Z(Zrm, n, m, K_VANRADEN, coding, native_A, NULL, NULL,
                            NULL, NULL);
    if (rc != 0)
    {
        switch (rc)
        {
        case -1010:
            diagnostic = "gs.make.A: memory allocation error";
            break;
        case -2001:
            diagnostic = "gs.make.A: marker coding was not recognized";
            break;
        case -2002:
            diagnostic = "gs.make.A: invalid diploid marker dosage";
            break;
        case -2003:
            diagnostic = "gs.make.A: allele-frequency scale must be positive";
            break;
        case -2004:
            diagnostic =
                "gs.make.A: marker coding is ambiguous; specify coding";
            break;
        case -2005:
            diagnostic = "gs.make.A: unsupported marker coding";
            break;
        case -2006:
            diagnostic = "gs.make.A: marker values do not match coding";
            break;
        default:
            diagnostic = "gs.make.A: relationship-matrix construction failed";
        }
        goto cleanup;
    }

    PROTECT(A = allocMatrix(REALSXP, n, n));
    nprot++;
    memcpy(REAL(A), native_A, (size_t)n * n * sizeof(double));
    SEXP names = getAttrib(Z, R_DimNamesSymbol);
    if (TYPEOF(names) == VECSXP && XLENGTH(names) == 2 &&
        VECTOR_ELT(names, 0) != R_NilValue)
    {
        SEXP labels;
        PROTECT(labels = allocVector(VECSXP, 2));
        nprot++;
        SET_VECTOR_ELT(labels, 0, VECTOR_ELT(names, 0));
        SET_VECTOR_ELT(labels, 1, VECTOR_ELT(names, 0));
        setAttrib(A, R_DimNamesSymbol, labels);
    }

cleanup:
    UNPROTECT(nprot);
    if (diagnostic != NULL)
    {
        st_info(-2, diagnostic);
        return R_NilValue;
    }
    return A;
}

static SEXP getListElement(SEXP list, const char *str)
{
    SEXP names = getAttrib(list, R_NamesSymbol);
    if (names == R_NilValue)
        return R_NilValue;

    for (R_xlen_t i = 0; i < XLENGTH(list); i++)
    {
        if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
        {
            return VECTOR_ELT(list, i);
        }
    }
    return R_NilValue;
}

#include <stdint.h>

typedef struct
{
    const char *name;
    size_t index;
} ste_name;

typedef struct
{
    const char *name;
    size_t length;
    size_t *base, *donor;
} ste_block;

typedef struct ste_variant
{
    size_t block;
    int hidden;
    GS_ALLELE_TYPE *sequence;
    uint64_t *carriers;
    struct ste_variant *next;
} ste_variant;

typedef struct
{
    SEXP base_arg, donor_arg, options, blocks_arg, exclude_arg, n_arg,
        method_arg;
    gs_varset_type tmp, *base, *donor;
    size_t nm, nb, nd, nblocks, nevents, *base_order, *donor_order;
    ste_name *marker_names, *donor_names;
    ste_block *blocks;
    int *main_count, exclude, select_n, forward;
    ste_variant *first, *last;
} ste_context;

static size_t ste_product(size_t a, size_t b)
{
    if (b && a > SIZE_MAX / b)
        Rf_error("Enrichment allocation size overflow");
    return a * b;
}

static void *ste_alloc(size_t n, size_t size)
{
    ste_product(n, size);
    return st_work_alloc(n ? n : 1, size);
}

static int ste_name_cmp(const void *a, const void *b)
{
    return strcmp(((const ste_name *)a)->name, ((const ste_name *)b)->name);
}

static size_t ste_lookup(const ste_name *a, size_t n, const char *name)
{
    size_t lo = 0, hi = n;
    while (lo < hi)
    {
        size_t m = lo + (hi - lo) / 2;
        int c = strcmp(name, a[m].name);
        if (!c)
            return a[m].index;
        if (c < 0)
            hi = m;
        else
            lo = m + 1;
    }
    return SIZE_MAX;
}

static const char *ste_string(SEXP x, const char *what)
{
    if (TYPEOF(x) != STRSXP || XLENGTH(x) != 1 ||
        STRING_ELT(x, 0) == NA_STRING || !LENGTH(STRING_ELT(x, 0)))
        Rf_error("%s must be one nonempty string", what);
    return CHAR(STRING_ELT(x, 0));
}

static int ste_integer(SEXP x, int minimum, const char *what)
{
    if ((TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP) || XLENGTH(x) != 1)
        Rf_error("%s must be an integer >= %d", what, minimum);
    double v = Rf_asReal(x);
    if (!R_FINITE(v) || v < minimum || v > INT_MAX || v != floor(v))
        Rf_error("%s must be an integer >= %d", what, minimum);
    return (int)v;
}

static SEXP ste_column(SEXP x, const char *name)
{
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(x) != VECSXP || TYPEOF(names) != STRSXP ||
        XLENGTH(names) != XLENGTH(x))
        return R_NilValue;
    for (R_xlen_t i = 0; i < XLENGTH(x); i++)
        if (STRING_ELT(names, i) != NA_STRING &&
            !strcmp(CHAR(STRING_ELT(names, i)), name))
            return VECTOR_ELT(x, i);
    return R_NilValue;
}

static SEXP ste_characters(SEXP x)
{
    return Rf_isFactor(x) ? Rf_asCharacterFactor(x)
                          : Rf_coerceVector(x, STRSXP);
}

static SEXP ste_frame(int nr, int nc, const char **names, const int *types)
{
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, nc)),
         nms = PROTECT(Rf_allocVector(STRSXP, nc));
    SEXP rn = PROTECT(Rf_allocVector(INTSXP, 2)),
         cls = PROTECT(Rf_mkString("data.frame"));
    for (int i = 0; i < nc; i++)
    {
        SET_STRING_ELT(nms, i, Rf_mkChar(names[i]));
        SET_VECTOR_ELT(ans, i, Rf_allocVector(types[i], nr));
    }
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -nr;
    Rf_setAttrib(ans, R_NamesSymbol, nms);
    Rf_setAttrib(ans, R_RowNamesSymbol, rn);
    Rf_setAttrib(ans, R_ClassSymbol, cls);
    UNPROTECT(4);
    return ans;
}

static gs_varset_type *ste_pool(SEXP arg)
{
    const char *name = ste_string(arg, "Pool name");
    gs_varset_type *x = NULL;
    if (!strcmp(name, "default"))
        x = gv;
    else
        for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
            if (!strcmp(name, gs_data[i]->name))
            {
                x = gs_data[i];
                break;
            }
    if (!x || !x->NoMar || !x->NoInd || !st_allele_present(x) || !x->MarName ||
        !x->IndName || !x->NLoc || !x->map || !x->mdl)
        Rf_error("Pool '%s' has no usable marker data/map", name);
    if (x->NoInd > INT_MAX / 2 || x->NoMar > INT_MAX)
        Rf_error("Pool dimensions exceed supported enrichment limits");
    return x;
}

static size_t *ste_map_order(gs_varset_type *x)
{
    size_t *order = ste_alloc(x->NoMar, sizeof(size_t));
    unsigned char *seen = ste_alloc(x->NoMar, 1);
    memset(seen, 0, x->NoMar);
    size_t n = 0;
    for (gs_UNSIGNED i = 0; i < x->NLoc; i++)
    {
        if (x->mdl[i] == gs_UNSIGNED_MAX)
            continue;
        if (x->mdl[i] >= x->NoMar || seen[x->mdl[i]] || !x->map[i] ||
            x->map[i]->chrom < 1 || !R_FINITE(x->map[i]->pos) ||
            x->map[i]->pos < 0 ||
            strcmp(x->map[i]->name, x->MarName[x->mdl[i]]))
            Rf_error("Invalid or ambiguous marker map");
        seen[x->mdl[i]] = 1;
        order[n++] = i;
    }
    if (n != x->NoMar)
        Rf_error("Every marker must have exactly one map location");
    return order;
}

static void ste_validate(ste_context *c)
{
    c->base = ste_pool(c->base_arg);
    c->donor = ste_pool(c->donor_arg);
    c->nb = c->base->NoInd;
    c->nd = c->donor->NoInd;
    c->nm = c->base->NoMar;
    if (c->nm != c->donor->NoMar)
        Rf_error("Pools must have the same marker set");
    if (TYPEOF(c->exclude_arg) != LGLSXP || XLENGTH(c->exclude_arg) != 1 ||
        LOGICAL(c->exclude_arg)[0] == NA_LOGICAL)
        Rf_error("exclude.compatible must be TRUE or FALSE");
    c->exclude = LOGICAL(c->exclude_arg)[0];
    if (c->n_arg != R_NilValue)
    {
        c->select_n = ste_integer(c->n_arg, 1, "n");
        if ((size_t)c->select_n > c->nd)
            Rf_error("n exceeds the enrichment pool size");
    }
    const char *method = ste_string(c->method_arg, "method");
    if (strcmp(method, "best.subset") && strcmp(method, "forward"))
        Rf_error("method must be 'best.subset' or 'forward'");
    c->forward = !strcmp(method, "forward");
    if ((c->options == R_NilValue) == (c->blocks_arg == R_NilValue))
        Rf_error("Supply exactly one of haplotype.options or blocks");
    c->base_order = ste_map_order(c->base);
    size_t *dm = ste_map_order(c->donor);
    c->donor_order = ste_alloc(c->nm, sizeof(size_t));
    c->marker_names = ste_alloc(c->nm, sizeof(ste_name));
    ste_name *dn = ste_alloc(c->nm, sizeof(ste_name));
    for (size_t i = 0; i < c->nm; i++)
    {
        c->marker_names[i] =
            (ste_name){c->base->map[c->base_order[i]]->name, i};
        dn[i] = (ste_name){c->donor->map[dm[i]]->name, dm[i]};
    }
    qsort(c->marker_names, c->nm, sizeof(ste_name), ste_name_cmp);
    qsort(dn, c->nm, sizeof(ste_name), ste_name_cmp);
    for (size_t i = 0; i < c->nm; i++)
    {
        if (i &&
            (!strcmp(c->marker_names[i - 1].name, c->marker_names[i].name) ||
             !strcmp(dn[i - 1].name, dn[i].name)))
            Rf_error("Marker names must be unique");
        if (strcmp(c->marker_names[i].name, dn[i].name))
            Rf_error("Pools must have the same marker set");
        size_t b = c->marker_names[i].index, d = dn[i].index;
        gs_mappoint_type *bp = c->base->map[c->base_order[b]],
                         *dp = c->donor->map[d];
        if (bp->chrom != dp->chrom || fabs(bp->pos - dp->pos) > 1e-8)
            Rf_error("Pools must have matching chromosome assignments and "
                     "positions");
        c->donor_order[b] = d;
    }
    c->donor_names = ste_alloc(c->nd, sizeof(ste_name));
    for (size_t i = 0; i < c->nd; i++)
        c->donor_names[i] = (ste_name){c->donor->IndName[i], i};
    qsort(c->donor_names, c->nd, sizeof(ste_name), ste_name_cmp);
    for (size_t i = 1; i < c->nd; i++)
        if (!strcmp(c->donor_names[i - 1].name, c->donor_names[i].name))
            Rf_error("Donor genotype names must be unique");
    ste_name *bn = ste_alloc(c->nb, sizeof(ste_name));
    for (size_t i = 0; i < c->nb; i++)
        bn[i] = (ste_name){c->base->IndName[i], i};
    qsort(bn, c->nb, sizeof(ste_name), ste_name_cmp);
    for (size_t i = 1; i < c->nb; i++)
        if (!strcmp(bn[i - 1].name, bn[i].name))
            Rf_error("Base genotype names must be unique");
}

static void ste_define_blocks(ste_context *c)
{
    SEXP op = c->options, names = Rf_getAttrib(op, R_NamesSymbol);
    if (TYPEOF(op) != VECSXP ||
        (XLENGTH(op) &&
         (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(op))))
        Rf_error("haplotype.options must be a named list");
    int hap = 5, unit = 1, tolerance = 0;
    double threshold = .8;
    const char *criterion = "none", *symbol = "b", *measure = "r2";
    const char *allowed[] = {"hap",          "hap.unit",  "ld.threshold",
                             "ld.criterion", "tolerance", "hap.symbol",
                             "ld.measure"};
    unsigned int seen = 0;
    for (R_xlen_t i = 0; i < XLENGTH(op); i++)
    {
        if (STRING_ELT(names, i) == NA_STRING)
            Rf_error("Invalid haplotyping option name");
        const char *key = CHAR(STRING_ELT(names, i));
        int k;
        for (k = 0; k < 7; k++)
            if (!strcmp(key, allowed[k]))
                break;
        if (k == 7 || (seen & (1U << k)))
            Rf_error("Unknown or repeated haplotyping option '%s'", key);
        seen |= 1U << k;
        SEXP v = VECTOR_ELT(op, i);
        switch (k)
        {
        case 0:
            hap = ste_integer(v, 1, key);
            break;
        case 1:
            unit = ste_integer(v, 0, key);
            if (unit > 2)
                Rf_error("hap.unit must be 0, 1 or 2");
            break;
        case 2:
            if ((TYPEOF(v) != INTSXP && TYPEOF(v) != REALSXP) ||
                XLENGTH(v) != 1)
                Rf_error("ld.threshold must lie between zero and one");
            threshold = Rf_asReal(v);
            if (!R_FINITE(threshold) || threshold < 0 || threshold > 1)
                Rf_error("ld.threshold must lie between zero and one");
            break;
        case 3:
            criterion = ste_string(v, key);
            break;
        case 4:
            tolerance = ste_integer(v, 0, key);
            break;
        case 5:
            symbol = ste_string(v, key);
            break;
        case 6:
            measure = ste_string(v, key);
            break;
        }
    }
    if (strcmp(criterion, "none") && strcmp(criterion, "flanking") &&
        strcmp(criterion, "average"))
        Rf_error("Invalid LD criterion");
    if (strcmp(measure, "r2") && strcmp(measure, "Dp") &&
        strcmp(measure, "r2-estimate-phases") &&
        strcmp(measure, "Dp-estimate-phases"))
        Rf_error("Invalid LD measure");
    if (strlen(symbol) > 20)
        Rf_error("hap.symbol must be 1 to 20 letters or digits");
    for (const unsigned char *p = (const unsigned char *)symbol; *p; p++)
        if (!((*p >= 'a' && *p <= 'z') || (*p >= 'A' && *p <= 'Z') ||
              (*p >= '0' && *p <= '9')))
            Rf_error("hap.symbol must be 1 to 20 letters or digits");
    int rv = 0, aux = 0;
    char *empty = "", *crit = (char *)criterion, *sym = (char *)symbol,
         *met = (char *)measure;
    c->tmp.info_level = c->base->info_level;
    snprintf(c->tmp.name, sizeof(c->tmp.name), "enrichment work");
    gs_copy_marker_data_01(&c->tmp, c->base, &rv);
    if (rv == -2)
        Rf_error("Cannot create temporary base data");
    if (strcmp(criterion, "none"))
    {
        gs_calc_ld_01(&c->tmp, &met, &empty, &aux, &rv);
        if (rv == -2)
            Rf_error("LD calculation failed");
    }
    gs_def_hblocks_01(&c->tmp, &hap, &unit, &threshold, &crit, &tolerance, &sym,
                      &empty, &aux, &rv);
    if (rv == -2 || !c->tmp.NHap || c->tmp.NHap > INT_MAX)
        Rf_error("No usable haplotype blocks were defined");
    c->nblocks = c->tmp.NHap;
    c->blocks = ste_alloc(c->nblocks, sizeof(ste_block));
    c->main_count = ste_alloc(c->nblocks, sizeof(int));
    for (size_t b = 0; b < c->nblocks; b++)
    {
        ste_block *block = &c->blocks[b];
        block->name = c->tmp.Hap_S[b].name;
        block->length = 0;
        for (size_t j = c->tmp.h_u[b]; j <= c->tmp.h_o[b]; j++)
            if (c->tmp.mdl[j] != gs_UNSIGNED_MAX)
                block->length++;
        if (!block->length)
            Rf_error("Empty haplotype block");
        block->base = ste_alloc(block->length, sizeof(size_t));
        block->donor = ste_alloc(block->length, sizeof(size_t));
        size_t k = 0;
        for (size_t j = c->tmp.h_u[b]; j <= c->tmp.h_o[b]; j++)
        {
            if (c->tmp.mdl[j] == gs_UNSIGNED_MAX)
                continue;
            size_t idx = ste_lookup(c->marker_names, c->nm,
                                    c->tmp.MarName[c->tmp.mdl[j]]);
            if (idx == SIZE_MAX)
                Rf_error("Unknown marker in generated haplotype block");
            block->base[k] = c->base->mdl[c->base_order[idx]];
            block->donor[k] = c->donor->mdl[c->donor_order[idx]];
            k++;
        }
    }
}

static SEXP ste_export_blocks(ste_context *c)
{
    const char *cols[] = {"Chrom", "Pos", "Name", "Class", "Markers"};
    int types[] = {INTSXP, REALSXP, STRSXP, STRSXP, STRSXP};
    SEXP ans = PROTECT(ste_frame((int)c->tmp.NHap, 5, cols, types));
    size_t capacity = 1;
    for (size_t i = 0; i < c->nm; i++)
    {
        size_t len = strlen(c->base->map[c->base_order[i]]->name) + 1;
        if (capacity > SIZE_MAX - len)
            Rf_error("Block text size overflow");
        capacity += len;
    }
    char *text = ste_alloc(capacity, 1);
    for (size_t b = 0; b < c->tmp.NHap; b++)
    {
        gs_mappoint_type *p = &c->tmp.Hap_S[b];
        INTEGER(VECTOR_ELT(ans, 0))[b] = p->chrom;

        char rounded[DBL_MAX_10_EXP + 32];
        snprintf(rounded, sizeof(rounded), "%.6f", p->pos);
        REAL(VECTOR_ELT(ans, 1))[b] = strtod(rounded, NULL);
        SET_STRING_ELT(VECTOR_ELT(ans, 2), b, Rf_mkChar(p->name));
        SET_STRING_ELT(VECTOR_ELT(ans, 3), b, Rf_mkChar(p->class));
        size_t used = 0;
        for (size_t j = c->tmp.h_u[b]; j <= c->tmp.h_o[b]; j++)
        {
            if (c->tmp.mdl[j] == gs_UNSIGNED_MAX)
                continue;
            const char *m = c->tmp.MarName[c->tmp.mdl[j]];
            size_t len = strlen(m);
            memcpy(text + used, m, len);
            used += len;
            text[used++] = ';';
        }
        text[used] = '\0';
        SET_STRING_ELT(VECTOR_ELT(ans, 4), b, Rf_mkChar(text));
    }
    UNPROTECT(1);
    return ans;
}

static void ste_parse_blocks(ste_context *c, SEXP table)
{
    SEXP n = ste_column(table, "Name"), m = ste_column(table, "Markers"),
         ch = ste_column(table, "Chrom");
    if (!Rf_inherits(table, "data.frame") || n == R_NilValue ||
        m == R_NilValue || ch == R_NilValue || !XLENGTH(n) ||
        XLENGTH(n) > INT_MAX || XLENGTH(m) != XLENGTH(n) ||
        XLENGTH(ch) != XLENGTH(n))
        Rf_error("blocks must contain Name, Chrom and Markers columns");
    SEXP names = PROTECT(ste_characters(n)),
         markers = PROTECT(ste_characters(m));
    SEXP chrom = PROTECT(Rf_isFactor(ch) ? Rf_asCharacterFactor(ch) : ch);
    c->nblocks = XLENGTH(n);
    c->blocks = ste_alloc(c->nblocks, sizeof(ste_block));
    c->main_count = ste_alloc(c->nblocks, sizeof(int));
    ste_name *bn = ste_alloc(c->nblocks, sizeof(ste_name));
    unsigned char *seen = ste_alloc(c->nm, 1);
    memset(seen, 0, c->nm);
    for (size_t i = 0; i < c->nblocks; i++)
    {
        if (STRING_ELT(names, i) == NA_STRING ||
            !LENGTH(STRING_ELT(names, i)) ||
            STRING_ELT(markers, i) == NA_STRING)
            Rf_error("Block names and marker lists must be nonmissing");
        const char *sname = CHAR(STRING_ELT(names, i));
        char *name = ste_alloc(strlen(sname) + 1, 1);
        strcpy(name, sname);
        c->blocks[i].name = name;
        bn[i] = (ste_name){name, i};
        double chromosome;
        if (TYPEOF(chrom) == STRSXP)
        {
            if (STRING_ELT(chrom, i) == NA_STRING)
                Rf_error("Invalid block chromosome");
            const char *s = CHAR(STRING_ELT(chrom, i));
            char *end;
            chromosome = strtod(s, &end);
            if (end == s || *end)
                Rf_error("Invalid block chromosome");
        }
        else if (TYPEOF(chrom) == INTSXP)
            chromosome = INTEGER(chrom)[i];
        else if (TYPEOF(chrom) == REALSXP)
            chromosome = REAL(chrom)[i];
        else
            Rf_error("Invalid block chromosome");
        if (!R_FINITE(chromosome) || chromosome < 1 || chromosome > INT_MAX ||
            chromosome != floor(chromosome))
            Rf_error("Invalid block chromosome");
        const char *s = CHAR(STRING_ELT(markers, i));
        size_t len = strlen(s);
        if (len && s[len - 1] == ';')
            len--;
        if (!len)
            Rf_error("Empty haplotype block");
        size_t count = 1;
        for (size_t j = 0; j < len; j++)
            if (s[j] == ';')
                count++;
        c->blocks[i].length = count;
        c->blocks[i].base = ste_alloc(count, sizeof(size_t));
        c->blocks[i].donor = ste_alloc(count, sizeof(size_t));
        char *copy = ste_alloc(len + 1, 1);
        memcpy(copy, s, len);
        copy[len] = '\0';
        size_t start = 0, k = 0, previous = SIZE_MAX;
        for (size_t j = 0; j <= len; j++)
            if (j == len || copy[j] == ';')
            {
                copy[j] = '\0';
                size_t idx = ste_lookup(c->marker_names, c->nm, copy + start);
                if (idx == SIZE_MAX || seen[idx])
                    Rf_error("Block markers must be known and nonoverlapping");
                if ((k && idx != previous + 1) ||
                    c->base->map[c->base_order[idx]]->chrom != (int)chromosome)
                    Rf_error("Each block must contain consecutive markers on "
                             "one chromosome");
                seen[idx] = 1;
                previous = idx;
                c->blocks[i].base[k] = c->base->mdl[c->base_order[idx]];
                c->blocks[i].donor[k] = c->donor->mdl[c->donor_order[idx]];
                k++;
                start = j + 1;
            }
    }
    qsort(bn, c->nblocks, sizeof(ste_name), ste_name_cmp);
    for (size_t i = 1; i < c->nblocks; i++)
        if (!strcmp(bn[i - 1].name, bn[i].name))
            Rf_error("Block names must be unique");
    UNPROTECT(3);
}

static uint64_t ste_hash(const GS_ALLELE_TYPE *s, size_t n)
{
    uint64_t h = UINT64_C(1469598103934665603);
    for (size_t i = 0; i < n; i++)
    {
        h ^= (uint64_t)s[i];
        h *= UINT64_C(1099511628211);
    }
    return h;
}

static void ste_catalogue(ste_context *c)
{
    size_t copies = 2 * (c->nb + c->nd), maxlen = 0;
    for (size_t i = 0; i < c->nblocks; i++)
        if (c->blocks[i].length > maxlen)
            maxlen = c->blocks[i].length;
    GS_ALLELE_TYPE *buffer =
        ste_alloc(ste_product(copies, maxlen), sizeof(GS_ALLELE_TYPE));
    unsigned char *complete = ste_alloc(copies, 1);
    size_t cap = 1;
    while (cap < ste_product(copies, 2))
    {
        if (cap > SIZE_MAX / 2)
            Rf_error("Haplotype hash size overflow");
        cap *= 2;
    }
    size_t *slots = ste_alloc(cap, sizeof(size_t));
    ste_variant **event = ste_alloc(cap, sizeof(ste_variant *));
    size_t donor_words = (c->nd + 63) / 64;
    for (size_t b = 0; b < c->nblocks; b++)
    {
        R_CheckUserInterrupt();
        ste_block *block = &c->blocks[b];
        size_t len = block->length;
        memset(slots, 0, cap * sizeof(size_t));
        memset(event, 0, cap * sizeof(ste_variant *));
        memset(complete, 1, copies);
        c->main_count[b] = 0;
        for (size_t p = 0; p < copies; p++)
        {
            int is_donor = p >= 2 * c->nb;
            gs_varset_type *pool = is_donor ? c->donor : c->base;
            size_t copy = is_donor ? p - 2 * c->nb : p,
                   ind = copy % pool->NoInd;
            int homologue = copy < pool->NoInd ? 0 : 1;
            size_t *rows = is_donor ? block->donor : block->base;
            GS_ALLELE_TYPE *seq = buffer + p * len;
            for (size_t j = 0; j < len; j++)
            {
                seq[j] = st_allele_at(pool, rows[j], ind, homologue);
                if (seq[j] == GS_ALLELE_MISSING)
                    complete[p] = 0;
            }
            if (!complete[p])
                continue;
            size_t slot = (size_t)ste_hash(seq, len) & (cap - 1);
            while (slots[slot] && memcmp(buffer + (slots[slot] - 1) * len, seq,
                                         len * sizeof(*seq)))
                slot = (slot + 1) & (cap - 1);
            if (!slots[slot])
            {
                slots[slot] = p + 1;
                if (!is_donor)
                {
                    c->main_count[b]++;
                    continue;
                }
                if (c->nevents >= INT_MAX)
                    Rf_error("Too many additional haplotypes");
                ste_variant *v = ste_alloc(1, sizeof(*v));
                v->block = b;
                v->hidden = 0;
                v->next = NULL;
                v->sequence = ste_alloc(len, sizeof(*seq));
                memcpy(v->sequence, seq, len * sizeof(*seq));
                v->carriers = ste_alloc(donor_words, sizeof(uint64_t));
                memset(v->carriers, 0, donor_words * sizeof(uint64_t));
                for (size_t j = 0; j < 2 * c->nb; j++)
                    if (!complete[j])
                    {
                        const GS_ALLELE_TYPE *base = buffer + j * len;
                        int match = 1;
                        for (size_t k = 0; k < len; k++)
                            if (base[k] != GS_ALLELE_MISSING &&
                                base[k] != seq[k])
                            {
                                match = 0;
                                break;
                            }
                        if (match)
                        {
                            v->hidden = 1;
                            break;
                        }
                    }
                if (c->last)
                    c->last->next = v;
                else
                    c->first = v;
                c->last = v;
                c->nevents++;
                event[slot] = v;
            }
            if (is_donor && event[slot])
                event[slot]->carriers[ind / 64] |= UINT64_C(1) << (ind % 64);
        }
    }
}

static SEXP ste_calculate(ste_context *c, SEXP blocks)
{
    const char *sn[] = {"Base.lines", "Enrichment.lines", "Main",
                        "Additional", "Hidden",           "No.match",
                        "Score",      "Enriched",         "Gaining.blocks"};
    int st[] = {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
                INTSXP, INTSXP, INTSXP, INTSXP};
    const char *bn[] = {"Block",      "Main",  "Hidden",  "No.match",
                        "Additional", "Score", "Enriched"};
    int bt[] = {STRSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};
    const char *an[] = {"Block", "Sequence", "Hidden"};
    int at[] = {STRSXP, STRSXP, LGLSXP};
    SEXP summary = PROTECT(ste_frame(1, 9, sn, st)),
         by = PROTECT(ste_frame((int)c->nblocks, 7, bn, bt));
    SEXP alleles = PROTECT(ste_frame((int)c->nevents, 3, an, at));
    for (size_t b = 0; b < c->nblocks; b++)
    {
        SET_STRING_ELT(VECTOR_ELT(by, 0), b, Rf_mkChar(c->blocks[b].name));
        INTEGER(VECTOR_ELT(by, 1))[b] = c->main_count[b];
        INTEGER(VECTOR_ELT(by, 2))[b] = INTEGER(VECTOR_ELT(by, 3))[b] = 0;
    }
    size_t maxlen = 0;
    for (size_t b = 0; b < c->nblocks; b++)
        if (c->blocks[b].length > maxlen)
            maxlen = c->blocks[b].length;
    size_t text_capacity = ste_product(maxlen, 4) + 1;
    char *text = ste_alloc(text_capacity, 1);
    size_t i = 0;
    for (ste_variant *v = c->first; v; v = v->next, i++)
    {
        INTEGER(VECTOR_ELT(by, v->hidden ? 2 : 3))[v->block]++;
        SET_STRING_ELT(VECTOR_ELT(alleles, 0), i,
                       Rf_mkChar(c->blocks[v->block].name));
        size_t used = 0, len = c->blocks[v->block].length;
        for (size_t j = 0; j < len; j++)
            used +=
                (size_t)snprintf(text + used, text_capacity - used,
                                 j ? ";%u" : "%u", (unsigned)v->sequence[j]);
        SET_STRING_ELT(VECTOR_ELT(alleles, 1), i, Rf_mkChar(text));
        LOGICAL(VECTOR_ELT(alleles, 2))[i] = v->hidden;
    }
    int64_t totals[9] = {(int)c->nb, (int)c->nd, 0, 0, 0, 0, 0, 0, 0};
    for (size_t b = 0; b < c->nblocks; b++)
    {
        int main = c->main_count[b], hidden = INTEGER(VECTOR_ELT(by, 2))[b],
            no = INTEGER(VECTOR_ELT(by, 3))[b];
        int additional = hidden + no, score = c->exclude ? no : additional;
        if ((int64_t)main + additional > INT_MAX)
            Rf_error("Allele count exceeds integer limit");
        INTEGER(VECTOR_ELT(by, 4))[b] = additional;
        INTEGER(VECTOR_ELT(by, 5))[b] = score;
        INTEGER(VECTOR_ELT(by, 6))[b] = main + additional;
        totals[2] += main;
        totals[3] += additional;
        totals[4] += hidden;
        totals[5] += no;
        totals[6] += score;
        totals[7] += main + additional;
        totals[8] += (additional > 0);
    }
    for (int j = 0; j < 9; j++)
    {
        if (totals[j] > INT_MAX)
            Rf_error("Allele count exceeds integer limit");
        INTEGER(VECTOR_ELT(summary, j))[0] = (int)totals[j];
    }
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 5)),
         names = PROTECT(Rf_allocVector(STRSXP, 5));
    const char *rn[] = {"summary", "by.block", "alleles", "blocks",
                        "exclude.compatible"};
    for (int j = 0; j < 5; j++)
        SET_STRING_ELT(names, j, Rf_mkChar(rn[j]));
    SET_VECTOR_ELT(ans, 0, summary);
    SET_VECTOR_ELT(ans, 1, by);
    SET_VECTOR_ELT(ans, 2, alleles);
    SET_VECTOR_ELT(ans, 3, blocks);
    SET_VECTOR_ELT(ans, 4, Rf_ScalarLogical(c->exclude));
    Rf_setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(5);
    return ans;
}

static int ste_popcount(uint64_t x)
{
    int n = 0;
    while (x)
    {
        x &= x - 1;
        n++;
    }
    return n;
}

static int ste_count(const uint64_t *x, size_t words)
{
    int n = 0;
    for (size_t w = 0; w < words; w++)
        n += ste_popcount(x[w]);
    return n;
}

static int ste_lex_less(const int *a, const int *b, int n)
{
    for (int i = 0; i < n; i++)
        if (a[i] != b[i])
            return a[i] < b[i];
    return 0;
}

typedef struct
{
    uint64_t *covered;
    int *chosen, *next, *task_best, *best, *gains, best_score;
} ste_worker;

static int ste_gain_cmp(const void *a, const void *b)
{
    int x = *(const int *)a, y = *(const int *)b;
    return (x < y) - (x > y);
}

static void ste_search_task(int first, int n, int nd, size_t words,
                            const uint64_t *inc, const uint64_t *suffix,
                            ste_worker *s)
{
    int level = 1, task_score = -1;
    s->chosen[0] = first;
    s->next[1] = first + 1;
    memcpy(s->covered + words, inc + (size_t)first * words,
           words * sizeof(uint64_t));
    while (level >= 1)
    {
        uint64_t *covered = s->covered + (size_t)level * words;
        if (level == n)
        {
            int score = ste_count(covered, words);
            if (score > task_score)
            {
                task_score = score;
                memcpy(s->task_best, s->chosen, n * sizeof(int));
            }
            level--;
            continue;
        }
        int start = s->next[level], left = n - level;
        if (start > nd - left)
        {
            level--;
            continue;
        }
        int upper = 0;
        for (size_t w = 0; w < words; w++)
            upper +=
                ste_popcount(covered[w] | suffix[(size_t)start * words + w]);
        if (upper <= task_score)
        {
            level--;
            continue;
        }
        int current = ste_count(covered, words), ng = nd - start;
        for (int j = start; j < nd; j++)
        {
            int gain = 0;
            for (size_t w = 0; w < words; w++)
                gain += ste_popcount(inc[(size_t)j * words + w] & ~covered[w]);
            s->gains[j - start] = gain;
        }
        qsort(s->gains, (size_t)ng, sizeof(int), ste_gain_cmp);
        int64_t bound = current;
        for (int j = 0; j < left; j++)
            bound += s->gains[j];
        if (bound <= task_score)
        {
            level--;
            continue;
        }
        int chosen = start;
        s->next[level] = start + 1;
        s->chosen[level] = chosen;
        uint64_t *next = s->covered + (size_t)(level + 1) * words;
        for (size_t w = 0; w < words; w++)
            next[w] = covered[w] | inc[(size_t)chosen * words + w];
        level++;
        s->next[level] = chosen + 1;
    }
    if (task_score > s->best_score ||
        (task_score == s->best_score && ste_lex_less(s->task_best, s->best, n)))
    {
        s->best_score = task_score;
        memcpy(s->best, s->task_best, n * sizeof(int));
    }
}

static int ste_threads(int nd, int tasks)
{
    (void)nd;
    (void)tasks;
    int threads = 1;
#ifdef MT
    if (nd >= 16 && tasks > 1)
    {
        threads = omp_get_max_threads();
        if (threads > MAXTHREADS)
            threads = MAXTHREADS;
        if (threads > tasks)
            threads = tasks;
    }
#endif
    return threads;
}

static void ste_forward(ste_context *c, const uint64_t *inc, size_t words,
                        SEXP ans)
{
    int nd = (int)c->nd, n = c->select_n, threads = ste_threads(nd, nd);
    (void)threads;
    uint64_t *covered = ste_alloc(words, sizeof(uint64_t));
    memset(covered, 0, words * sizeof(uint64_t));
    unsigned char *selected = ste_alloc(nd, 1);
    memset(selected, 0, nd);
    int *gains = ste_alloc(nd, sizeof(int));

    for (int step = 0; step < n; step++)
    {
        R_CheckUserInterrupt();
#ifdef MT
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
#endif
        for (int j = 0; j < nd; j++)
        {
            int gain = -1;
            if (!selected[j])
            {
                gain = 0;
                for (size_t w = 0; w < words; w++)
                    gain +=
                        ste_popcount(inc[(size_t)j * words + w] & ~covered[w]);
            }
            gains[j] = gain;
        }

        int best = -1, best_gain = -1;
        for (int j = 0; j < nd; j++)
            if (gains[j] > best_gain)
            {
                best = j;
                best_gain = gains[j];
            }
        selected[best] = 1;
        for (size_t w = 0; w < words; w++)
            covered[w] |= inc[(size_t)best * words + w];
        SET_STRING_ELT(ans, step, Rf_mkChar(c->donor_names[best].name));
    }
}

static SEXP ste_select(ste_context *c)
{
    int n = c->select_n, nd = (int)c->nd;
    size_t qualifying = 0;
    for (ste_variant *v = c->first; v; v = v->next)
        if (!c->exclude || !v->hidden)
            qualifying++;
    SEXP ans = PROTECT(Rf_allocVector(STRSXP, n));
    if (!qualifying || (!c->forward && n == nd))
    {
        for (int i = 0; i < n; i++)
            SET_STRING_ELT(ans, i, Rf_mkChar(c->donor_names[i].name));
        UNPROTECT(1);
        return ans;
    }
    size_t words = (qualifying + 63) / 64;
    uint64_t *inc = ste_alloc(ste_product(c->nd, words), sizeof(uint64_t));
    memset(inc, 0, ste_product(c->nd, words) * sizeof(uint64_t));
    size_t k = 0;
    for (ste_variant *v = c->first; v; v = v->next)
        if (!c->exclude || !v->hidden)
        {
            for (size_t j = 0; j < c->nd; j++)
            {
                size_t original = c->donor_names[j].index;
                if (v->carriers[original / 64] &
                    (UINT64_C(1) << (original % 64)))
                    inc[j * words + k / 64] |= UINT64_C(1) << (k % 64);
            }
            k++;
        }
    if (c->forward)
    {
        ste_forward(c, inc, words, ans);
        UNPROTECT(1);
        return ans;
    }
    if (n == 1)
    {
        int best = -1, index = 0;
        for (int j = 0; j < nd; j++)
        {
            int score = ste_count(inc + (size_t)j * words, words);
            if (score > best)
            {
                best = score;
                index = j;
            }
        }
        SET_STRING_ELT(ans, 0, Rf_mkChar(c->donor_names[index].name));
        UNPROTECT(1);
        return ans;
    }
    uint64_t *suffix =
        ste_alloc(ste_product(c->nd + 1, words), sizeof(uint64_t));
    memset(suffix + c->nd * words, 0, words * sizeof(uint64_t));
    for (int j = nd - 1; j >= 0; j--)
        for (size_t w = 0; w < words; w++)
            suffix[(size_t)j * words + w] =
                suffix[(size_t)(j + 1) * words + w] |
                inc[(size_t)j * words + w];
    int tasks = nd - n + 1, threads = ste_threads(nd, tasks);
    ste_worker *workers = ste_alloc(threads, sizeof(ste_worker));
    for (int t = 0; t < threads; t++)
    {
        ste_worker *s = &workers[t];
        s->best_score = -1;
        s->covered =
            ste_alloc(ste_product((size_t)n + 1, words), sizeof(uint64_t));
        s->chosen = ste_alloc(n, sizeof(int));
        s->next = ste_alloc((size_t)n + 1, sizeof(int));
        s->task_best = ste_alloc(n, sizeof(int));
        s->best = ste_alloc(n, sizeof(int));
        s->gains = ste_alloc(nd, sizeof(int));
    }
    for (int begin = 0; begin < tasks; begin += 32)
    {
        int end = tasks - begin < 32 ? tasks : begin + 32;
        R_CheckUserInterrupt();
#ifdef MT
#pragma omp parallel for schedule(dynamic, 1)                                  \
    num_threads(threads) if (threads > 1)
#endif
        for (int task = begin; task < end; task++)
        {
            int t = 0;
#ifdef MT
            t = omp_get_thread_num();
#endif
            ste_search_task(task, n, nd, words, inc, suffix, &workers[t]);
        }
    }
    int winner = 0;
    for (int t = 1; t < threads; t++)
        if (workers[t].best_score > workers[winner].best_score ||
            (workers[t].best_score == workers[winner].best_score &&
             workers[t].best_score >= 0 &&
             ste_lex_less(workers[t].best, workers[winner].best, n)))
            winner = t;
    for (int i = 0; i < n; i++)
        SET_STRING_ELT(ans, i,
                       Rf_mkChar(c->donor_names[workers[winner].best[i]].name));
    UNPROTECT(1);
    return ans;
}

static SEXP ste_body(void *data)
{
    ste_context *c = (ste_context *)data;
    ste_validate(c);
    if (c->blocks_arg != R_NilValue)
        ste_parse_blocks(c, c->blocks_arg);
    else
        ste_define_blocks(c);
    ste_catalogue(c);
    if (c->select_n)
        return ste_select(c);
    SEXP blocks = PROTECT(c->blocks_arg != R_NilValue ? c->blocks_arg
                                                      : ste_export_blocks(c));
    SEXP ans = PROTECT(ste_calculate(c, blocks));
    UNPROTECT(2);
    return ans;
}

static void ste_cleanup(void *data)
{
    ste_context *c = (ste_context *)data;
    gs_reset(&c->tmp);
}

SEXP st_enrichment(SEXP base, SEXP donor, SEXP options, SEXP blocks,
                   SEXP exclude, SEXP n, SEXP method)
{
    ste_context c;
    memset(&c, 0, sizeof(c));
    gs_init(&c.tmp);
    c.base_arg = base;
    c.donor_arg = donor;
    c.options = options;
    c.blocks_arg = blocks;
    c.exclude_arg = exclude;
    c.n_arg = n;
    c.method_arg = method;
    return R_ExecWithCleanup(ste_body, &c, ste_cleanup, &c);
}

static void ste_core_catalogue(ste_context *c)
{
    size_t copies = 2 * c->nd, maxlen = 0;
    for (size_t i = 0; i < c->nblocks; i++)
        if (c->blocks[i].length > maxlen)
            maxlen = c->blocks[i].length;
    GS_ALLELE_TYPE *buffer =
        ste_alloc(ste_product(copies, maxlen), sizeof(GS_ALLELE_TYPE));
    unsigned char *complete = ste_alloc(copies, 1);
    size_t cap = 1;
    while (cap < ste_product(copies, 2))
    {
        if (cap > SIZE_MAX / 2)
            Rf_error("Haplotype hash size overflow");
        cap *= 2;
    }
    size_t *slots = ste_alloc(cap, sizeof(size_t));
    ste_variant **event = ste_alloc(cap, sizeof(ste_variant *));
    size_t donor_words = (c->nd + 63) / 64;
    for (size_t b = 0; b < c->nblocks; b++)
    {
        R_CheckUserInterrupt();
        ste_block *block = &c->blocks[b];
        size_t len = block->length;
        memset(slots, 0, cap * sizeof(size_t));
        memset(event, 0, cap * sizeof(ste_variant *));
        memset(complete, 1, copies);
        c->main_count[b] = 0;
        for (size_t p = 0; p < copies; p++)
        {
            size_t ind = p % c->nd;
            int homologue = p < c->nd ? 0 : 1;
            GS_ALLELE_TYPE *seq = buffer + p * len;
            for (size_t j = 0; j < len; j++)
            {
                seq[j] =
                    st_allele_at(c->donor, block->donor[j], ind, homologue);
                if (seq[j] == GS_ALLELE_MISSING)
                    complete[p] = 0;
            }
            if (!complete[p])
                continue;
            size_t slot = (size_t)ste_hash(seq, len) & (cap - 1);
            while (slots[slot] && memcmp(buffer + (slots[slot] - 1) * len, seq,
                                         len * sizeof(*seq)))
                slot = (slot + 1) & (cap - 1);
            if (!slots[slot])
            {
                slots[slot] = p + 1;
                if (c->nevents >= INT_MAX)
                    Rf_error("Too many core haplotypes");
                ste_variant *v = ste_alloc(1, sizeof(*v));
                v->block = b;
                v->hidden = 0;
                v->next = NULL;
                v->sequence = ste_alloc(len, sizeof(*seq));
                memcpy(v->sequence, seq, len * sizeof(*seq));
                v->carriers = ste_alloc(donor_words, sizeof(uint64_t));
                memset(v->carriers, 0, donor_words * sizeof(uint64_t));
                if (c->last)
                    c->last->next = v;
                else
                    c->first = v;
                c->last = v;
                c->nevents++;
                c->main_count[b]++;
                event[slot] = v;
            }
            event[slot]->carriers[ind / 64] |= UINT64_C(1) << (ind % 64);
        }
    }
}

static void ste_core_add_compatible(ste_context *c, size_t original,
                                    uint64_t *hidden, size_t words)
{
    (void)words;
    size_t k = 0;
    for (ste_variant *v = c->first; v; v = v->next, k++)
    {
        ste_block *block = &c->blocks[v->block];
        for (int h = 0; h < 2; h++)
        {

            int incomplete = 0, match = 1;
            for (size_t j = 0; j < block->length; j++)
            {
                GS_ALLELE_TYPE a =
                    st_allele_at(c->donor, block->donor[j], original, h);
                if (a == GS_ALLELE_MISSING)
                    incomplete = 1;
                else if (a != v->sequence[j])
                {
                    match = 0;
                    break;
                }
            }
            if (incomplete && match)
            {
                hidden[k / 64] |= UINT64_C(1) << (k % 64);
                break;
            }
        }
    }
}

typedef struct
{
    uint64_t *covered, *hidden;
    int *chosen, *next, *task_best, *best, *gains;
    int best_score, best_segments;
} ste_core_worker;

static int ste_core_score(const uint64_t *covered, const uint64_t *hidden,
                          size_t words, int exclude)
{
    int n = 0;
    for (size_t w = 0; w < words; w++)
        n += ste_popcount(exclude ? covered[w] & ~hidden[w] : covered[w]);
    return n;
}

static void ste_core_search_task(int first, int n, int nd, size_t words,
                                 const uint64_t *inc,
                                 const uint64_t *compatible,
                                 const uint64_t *suffix, int exclude,
                                 ste_core_worker *s)
{
    int level = 1, task_score = -1, task_segments = 0;
    s->chosen[0] = first;
    s->next[1] = first + 1;
    memcpy(s->covered + words, inc + (size_t)first * words,
           words * sizeof(uint64_t));
    if (exclude)
        memcpy(s->hidden + words, compatible + (size_t)first * words,
               words * sizeof(uint64_t));
    else
        memset(s->hidden + words, 0, words * sizeof(uint64_t));

    while (level >= 1)
    {
        uint64_t *covered = s->covered + (size_t)level * words;
        uint64_t *hidden = s->hidden + (size_t)level * words;
        if (level == n)
        {
            int score = ste_core_score(covered, hidden, words, exclude);
            if (score > task_score)
            {
                task_score = score;
                task_segments = ste_count(covered, words);
                memcpy(s->task_best, s->chosen, n * sizeof(int));
            }
            level--;
            continue;
        }

        int start = s->next[level], left = n - level;
        if (start > nd - left)
        {
            level--;
            continue;
        }

        int upper = 0;
        for (size_t w = 0; w < words; w++)
        {
            uint64_t possible = covered[w] | suffix[(size_t)start * words + w];
            if (exclude)
                possible &= ~hidden[w];
            upper += ste_popcount(possible);
        }
        if (upper <= task_score)
        {
            level--;
            continue;
        }

        int current = ste_core_score(covered, hidden, words, exclude);
        int ng = nd - start;
        for (int j = start; j < nd; j++)
        {
            int gain = 0;
            for (size_t w = 0; w < words; w++)
            {
                uint64_t add = inc[(size_t)j * words + w] & ~covered[w];
                if (exclude)
                    add &= ~hidden[w];
                gain += ste_popcount(add);
            }
            s->gains[j - start] = gain;
        }
        qsort(s->gains, (size_t)ng, sizeof(int), ste_gain_cmp);
        int64_t bound = current;
        for (int j = 0; j < left; j++)
            bound += s->gains[j];
        if (bound <= task_score)
        {
            level--;
            continue;
        }

        int chosen = start;
        s->next[level] = start + 1;
        s->chosen[level] = chosen;
        uint64_t *next_covered = s->covered + (size_t)(level + 1) * words;
        uint64_t *next_hidden = s->hidden + (size_t)(level + 1) * words;
        for (size_t w = 0; w < words; w++)
        {
            next_covered[w] = covered[w] | inc[(size_t)chosen * words + w];
            next_hidden[w] =
                exclude ? hidden[w] | compatible[(size_t)chosen * words + w]
                        : 0;
        }
        level++;
        s->next[level] = chosen + 1;
    }

    if (task_score > s->best_score ||
        (task_score == s->best_score && task_score >= 0 &&
         ste_lex_less(s->task_best, s->best, n)))
    {
        s->best_score = task_score;
        s->best_segments = task_segments;
        memcpy(s->best, s->task_best, n * sizeof(int));
    }
}

static SEXP ste_core_forward(ste_context *c)
{
    int nd = (int)c->nd, nmax = c->select_n, threads = ste_threads(nd, nd);
    (void)threads;
    size_t words = (c->nevents + 63) / 64;
    int *chosen = ste_alloc(nmax, sizeof(int));
    int *segments = ste_alloc(nmax, sizeof(int));
    int nsel = 0;

    if (c->nevents)
    {
        uint64_t *inc = ste_alloc(ste_product(c->nd, words), sizeof(uint64_t));
        memset(inc, 0, ste_product(c->nd, words) * sizeof(uint64_t));
        size_t k = 0;
        for (ste_variant *v = c->first; v; v = v->next, k++)
            for (size_t j = 0; j < c->nd; j++)
            {
                size_t original = c->donor_names[j].index;
                if (v->carriers[original / 64] &
                    (UINT64_C(1) << (original % 64)))
                    inc[j * words + k / 64] |= UINT64_C(1) << (k % 64);
            }

        uint64_t *covered = ste_alloc(words, sizeof(uint64_t));
        uint64_t *hidden = ste_alloc(words, sizeof(uint64_t));
        memset(covered, 0, words * sizeof(uint64_t));
        memset(hidden, 0, words * sizeof(uint64_t));
        unsigned char *selected = ste_alloc(nd, 1);
        memset(selected, 0, nd);
        int *gains = ste_alloc(nd, sizeof(int));

        for (int step = 0; step < nmax; step++)
        {
            R_CheckUserInterrupt();
#ifdef MT
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
#endif
            for (int j = 0; j < nd; j++)
            {
                int gain = -1;
                if (!selected[j])
                {
                    gain = 0;
                    for (size_t w = 0; w < words; w++)
                    {
                        uint64_t add = inc[(size_t)j * words + w] & ~covered[w];
                        if (c->exclude)
                            add &= ~hidden[w];
                        gain += ste_popcount(add);
                    }
                }
                gains[j] = gain;
            }
            int best = -1, best_gain = 0;
            for (int j = 0; j < nd; j++)
                if (gains[j] > best_gain)
                {
                    best = j;
                    best_gain = gains[j];
                }
            if (best < 0)
                break;

            selected[best] = 1;
            chosen[nsel] = best;
            for (size_t w = 0; w < words; w++)
                covered[w] |= inc[(size_t)best * words + w];
            segments[nsel] = ste_count(covered, words);
            if (c->exclude)
                ste_core_add_compatible(c, c->donor_names[best].index, hidden,
                                        words);
            nsel++;
        }
    }

    SEXP sets = PROTECT(Rf_allocVector(VECSXP, nsel));
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, nsel));
    for (int i = 0; i < nsel; i++)
    {
        SEXP set = PROTECT(Rf_allocVector(STRSXP, i + 1));
        for (int j = 0; j <= i; j++)
            SET_STRING_ELT(set, j, Rf_mkChar(c->donor_names[chosen[j]].name));
        SET_VECTOR_ELT(sets, i, set);
        INTEGER(counts)[i] = segments[i];
        UNPROTECT(1);
    }
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("genotypes"));
    SET_STRING_ELT(names, 1, Rf_mkChar("segments"));
    SET_VECTOR_ELT(ans, 0, sets);
    SET_VECTOR_ELT(ans, 1, counts);
    Rf_setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(4);
    return ans;
}

static SEXP ste_core_best_subset(ste_context *c)
{
    int nd = (int)c->nd, nmax = c->select_n;
    size_t words = (c->nevents + 63) / 64;
    if (!c->nevents)
    {
        SEXP ans = PROTECT(Rf_allocVector(VECSXP, 2));
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(names, 0, Rf_mkChar("genotypes"));
        SET_STRING_ELT(names, 1, Rf_mkChar("segments"));
        SET_VECTOR_ELT(ans, 0, Rf_allocVector(VECSXP, 0));
        SET_VECTOR_ELT(ans, 1, Rf_allocVector(INTSXP, 0));
        Rf_setAttrib(ans, R_NamesSymbol, names);
        UNPROTECT(2);
        return ans;
    }

    uint64_t *inc = ste_alloc(ste_product(c->nd, words), sizeof(uint64_t));
    memset(inc, 0, ste_product(c->nd, words) * sizeof(uint64_t));
    size_t event = 0;
    for (ste_variant *v = c->first; v; v = v->next, event++)
        for (size_t j = 0; j < c->nd; j++)
        {
            size_t original = c->donor_names[j].index;
            if (v->carriers[original / 64] & (UINT64_C(1) << (original % 64)))
                inc[j * words + event / 64] |= UINT64_C(1) << (event % 64);
        }

    uint64_t *compatible = NULL;
    int maxthreads = ste_threads(nd, nd);
    if (c->exclude)
    {
        compatible = ste_alloc(ste_product(c->nd, words), sizeof(uint64_t));
        memset(compatible, 0, ste_product(c->nd, words) * sizeof(uint64_t));
#ifdef MT
#pragma omp parallel for schedule(static)                                      \
    num_threads(maxthreads) if (maxthreads > 1)
#endif
        for (int j = 0; j < nd; j++)
            ste_core_add_compatible(c, c->donor_names[j].index,
                                    compatible + (size_t)j * words, words);
    }

    uint64_t *suffix =
        ste_alloc(ste_product(c->nd + 1, words), sizeof(uint64_t));
    memset(suffix + c->nd * words, 0, words * sizeof(uint64_t));
    for (int j = nd - 1; j >= 0; j--)
        for (size_t w = 0; w < words; w++)
            suffix[(size_t)j * words + w] =
                suffix[(size_t)(j + 1) * words + w] |
                inc[(size_t)j * words + w];

    ste_core_worker *workers = ste_alloc(maxthreads, sizeof(ste_core_worker));
    for (int t = 0; t < maxthreads; t++)
    {
        ste_core_worker *s = &workers[t];
        s->covered =
            ste_alloc(ste_product((size_t)nmax + 1, words), sizeof(uint64_t));
        s->hidden =
            ste_alloc(ste_product((size_t)nmax + 1, words), sizeof(uint64_t));
        s->chosen = ste_alloc(nmax, sizeof(int));
        s->next = ste_alloc((size_t)nmax + 1, sizeof(int));
        s->task_best = ste_alloc(nmax, sizeof(int));
        s->best = ste_alloc(nmax, sizeof(int));
        s->gains = ste_alloc(nd, sizeof(int));
    }

    int **allsets = ste_alloc(nmax, sizeof(*allsets));
    int *allcounts = ste_alloc(nmax, sizeof(*allcounts));
    int nout = 0;

    for (int n = 1; n <= nmax; n++)
    {
        int tasks = nd - n + 1, threads = ste_threads(nd, tasks);
        for (int t = 0; t < threads; t++)
        {
            workers[t].best_score = -1;
            workers[t].best_segments = 0;
        }

        for (int begin = 0; begin < tasks; begin += 32)
        {
            int end = tasks - begin < 32 ? tasks : begin + 32;
            R_CheckUserInterrupt();
#ifdef MT
#pragma omp parallel for schedule(dynamic, 1)                                  \
    num_threads(threads) if (threads > 1)
#endif
            for (int task = begin; task < end; task++)
            {
                int t = 0;
#ifdef MT
                t = omp_get_thread_num();
#endif
                ste_core_search_task(task, n, nd, words, inc, compatible,
                                     suffix, c->exclude, &workers[t]);
            }
        }

        int winner = 0;
        for (int t = 1; t < threads; t++)
            if (workers[t].best_score > workers[winner].best_score ||
                (workers[t].best_score == workers[winner].best_score &&
                 workers[t].best_score >= 0 &&
                 ste_lex_less(workers[t].best, workers[winner].best, n)))
                winner = t;

        if (workers[winner].best_score <= 0)
            break;

        allsets[nout] = ste_alloc(n, sizeof(int));
        memcpy(allsets[nout], workers[winner].best, (size_t)n * sizeof(int));
        allcounts[nout] = workers[winner].best_segments;
        nout++;

        if (workers[winner].best_segments == (int)c->nevents)
            break;
    }

    SEXP sets = PROTECT(Rf_allocVector(VECSXP, nout));
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, nout));
    for (int i = 0; i < nout; i++)
    {
        SEXP set = PROTECT(Rf_allocVector(STRSXP, i + 1));
        for (int j = 0; j <= i; j++)
            SET_STRING_ELT(set, j,
                           Rf_mkChar(c->donor_names[allsets[i][j]].name));
        SET_VECTOR_ELT(sets, i, set);
        INTEGER(counts)[i] = allcounts[i];
        UNPROTECT(1);
    }
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("genotypes"));
    SET_STRING_ELT(names, 1, Rf_mkChar("segments"));
    SET_VECTOR_ELT(ans, 0, sets);
    SET_VECTOR_ELT(ans, 1, counts);
    Rf_setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(4);
    return ans;
}

static SEXP ste_core_body(void *data)
{
    ste_context *c = (ste_context *)data;
    ste_validate(c);
    if (c->blocks_arg != R_NilValue)
        ste_parse_blocks(c, c->blocks_arg);
    else
        ste_define_blocks(c);
    ste_core_catalogue(c);
    return c->forward ? ste_core_forward(c) : ste_core_best_subset(c);
}

SEXP st_detect_core(SEXP pool, SEXP nmax, SEXP options, SEXP blocks,
                    SEXP exclude, SEXP method)
{
    ste_context c;
    memset(&c, 0, sizeof(c));
    gs_init(&c.tmp);
    c.base_arg = pool;
    c.donor_arg = pool;
    c.options = options;
    c.blocks_arg = blocks;
    c.exclude_arg = exclude;
    c.n_arg = nmax;
    c.method_arg = method;
    return R_ExecWithCleanup(ste_core_body, &c, ste_cleanup, &c);
}

#include "SelectionTools01a.h"

extern void st_adjust_pvalues(double *, int, const char *);

static const char *st_legacy_marker_string(SEXP value, const char *name)
{
    if (TYPEOF(value) != STRSXP || XLENGTH(value) != 1 ||
        STRING_ELT(value, 0) == NA_STRING)
        error("%s must be one character string", name);
    return CHAR(STRING_ELT(value, 0));
}

static const char *st_legacy_adjust_method(SEXP value)
{
    const char *name = st_legacy_marker_string(value, "p.adjust.method");
    const char *methods[] = {"holm", "hochberg", "hommel", "bonferroni",
                             "BH",   "BY",       "fdr",    "none"};
    int selected = -1, matches = 0;
    for (int i = 0; i < 8; i++)
        if (!strcmp(name, methods[i]))
            return methods[i];
        else if (!strncmp(name, methods[i], strlen(name)))
            selected = i, matches++;
    if (matches != 1)
        error("Invalid or ambiguous p.adjust.method");
    return methods[selected];
}

static double st_legacy_export_precision(double value, int digits)
{
    char text[DBL_MAX_10_EXP + 64];
    snprintf(text, sizeof(text), "%.*f", digits, value);
    return strtod(text, NULL);
}

static const char *st_legacy_input_string(SEXP x, R_xlen_t i)
{
    if (TYPEOF(x) == STRSXP)
        return STRING_ELT(x, i) == NA_STRING ? NULL : CHAR(STRING_ELT(x, i));
    if (Rf_isFactor(x))
    {
        int level = INTEGER(x)[i];
        SEXP levels = getAttrib(x, R_LevelsSymbol);
        if (level == NA_INTEGER || level < 1 || level > LENGTH(levels))
            return NULL;
        return CHAR(STRING_ELT(levels, level - 1));
    }
    error("Marker names must be character strings or factors");
    return NULL;
}

static void st_legacy_mark_name(const gs_name_index_type *index, gs_UNSIGNED n,
                                unsigned char *selected, const char *name,
                                size_t length)
{
    if (length >= gs_MAXNME)
        return;
    char label[gs_MAXNME];
    memcpy(label, name, length);
    label[length] = '\0';
    gs_UNSIGNED marker = gs_name_index_find(index, n, label);
    if (marker != gs_UNSIGNED_MAX)
        selected[marker] = 1;
}

static SEXP st_legacy_mark_alleles(SEXP args)
{
    const char *source_name =
        st_legacy_marker_string(getListElement(args, "source"), "data.set");
    const char *target_name =
        st_legacy_marker_string(getListElement(args, "target"), "target.set");
    SEXP replacement = getListElement(args, "replace");
    const char *replacement_text =
        st_legacy_marker_string(replacement, "replace");
    for (const char *p = replacement_text; *p; p++)
        if (*p < '0' || *p > '9')
            error("replace must be an allele code from 1 to 254");
    char *replacement_end;
    long replace = strtol(replacement_text, &replacement_end, 10);
    if (!*replacement_text || *replacement_end || replace < 1 ||
        replace > GS_ALLELE_MAX)
        error("replace must be an allele code from 1 to 254");
    gs_varset_type *source = st_find_data_set_no_create(source_name);
    if (!source || !source->NoMar || !source->NoInd ||
        !st_allele_present(source))
        error("Source marker data are not available");
    unsigned char *selected = st_work_alloc(source->NoMar, 1);
    memset(selected, 0, source->NoMar);
    gs_name_index_type *marker_index =
        st_work_alloc(source->NoMar, sizeof(*marker_index));
    for (gs_UNSIGNED i = 0; i < source->NoMar; i++)
    {
        marker_index[i].name = source->MarName[i];
        marker_index[i].idx = i;
    }
    qsort(marker_index, source->NoMar, sizeof(*marker_index),
          gs_name_index_cmp);
    SEXP markers = getListElement(args, "markers");
    if (markers != R_NilValue)
    {
        for (R_xlen_t i = 0; i < XLENGTH(markers); i++)
        {
            const char *name = st_legacy_input_string(markers, i);
            if (name)
                st_legacy_mark_name(marker_index, source->NoMar, selected, name,
                                    strlen(name));
        }
    }
    else
    {
        SEXP effect_arg = getListElement(args, "effect");
        const char *effect_name =
            isNull(effect_arg)
                ? source_name
                : st_legacy_marker_string(effect_arg, "effect.set");
        gs_varset_type *effect = st_find_data_set_no_create(effect_name);
        if (!effect || !effect->p || !effect->u || !effect->EffNme ||
            !effect->EffAll)
        {
            gs_info(source, -2, "No p-values estimated");
            return R_NilValue;
        }
        if (effect->dimu > INT_MAX)
            error("Too many p-values for multiplicity adjustment");
        const char *method =
            st_legacy_adjust_method(getListElement(args, "adjust"));
        double alpha = asReal(getListElement(args, "alpha"));
        double *probability = st_work_alloc(effect->dimu, sizeof(double));
        gs_UNSIGNED *index = st_work_alloc(effect->dimu, sizeof(gs_UNSIGNED));
        unsigned char *missing = st_work_alloc(effect->dimu, 1);
        int n = 0;
        for (gs_UNSIGNED i = 1; i < effect->dimu; i++)
            if (fabs(st_legacy_export_precision(effect->u[i], 6)) > 1e-8)
            {
                index[n] = i;
                probability[n] = st_legacy_export_precision(effect->p[i], 6);
                missing[n] = ISNAN(probability[n]);
                if (missing[n])
                    probability[n] = 1;
                n++;
            }
        st_adjust_pvalues(probability, n, method);
        SEXP haplotypes = getListElement(args, "haplotypes");
        SEXP hapnames = isNull(haplotypes) ? R_NilValue
                                           : getListElement(haplotypes, "Name");
        SEXP membership = isNull(haplotypes)
                              ? R_NilValue
                              : getListElement(haplotypes, "Markers");
        if (!isNull(haplotypes) && (isNull(hapnames) || isNull(membership) ||
                                    XLENGTH(hapnames) != XLENGTH(membership)))
            error("hap.lst must contain matching Name and Markers columns");
        R_xlen_t nhap = isNull(hapnames) ? 0 : XLENGTH(hapnames), valid_hap = 0;
        gs_name_index_type *hap_index =
            st_work_alloc(nhap ? nhap : 1, sizeof(*hap_index));
        for (R_xlen_t i = 0; i < nhap; i++)
        {
            const char *name = st_legacy_input_string(hapnames, i);
            if (name)
            {
                hap_index[valid_hap].name = name;
                hap_index[valid_hap++].idx = i;
            }
        }
        qsort(hap_index, valid_hap, sizeof(*hap_index), gs_name_index_cmp);
        for (int i = 0; i < n; i++)
        {
            if (missing[i] || !(probability[i] < alpha))
                continue;
            char label[gs_MAXNME + 32];
            snprintf(label, sizeof(label), "%s.%i", effect->EffNme[index[i]],
                     effect->EffAll[index[i]]);
            size_t bytes = 0;
            for (int chars = 0; label[bytes] && chars < 7; chars++)
            {
                bytes++;
                while (((unsigned char)label[bytes] & 0xc0) == 0x80)
                    bytes++;
            }
            label[bytes] = '\0';
            if (isNull(haplotypes))
                st_legacy_mark_name(marker_index, source->NoMar, selected,
                                    label, bytes);
            else
            {
                R_xlen_t low = 0, high = valid_hap;
                while (low < high)
                {
                    R_xlen_t mid = low + (high - low) / 2;
                    if (strcmp(hap_index[mid].name, label) < 0)
                        low = mid + 1;
                    else
                        high = mid;
                }
                for (R_xlen_t h = low;
                     h < valid_hap && !strcmp(hap_index[h].name, label); h++)
                {
                    R_xlen_t j = hap_index[h].idx;
                    const char *list = st_legacy_input_string(membership, j);
                    if (!list)
                        continue;
                    const char *start = list;
                    for (const char *p = list;; p++)
                        if (!*p || *p == ';')
                        {
                            if (p > start)
                                st_legacy_mark_name(marker_index, source->NoMar,
                                                    selected, start, p - start);
                            start = p + 1;
                            if (!*p)
                                break;
                        }
                }
            }
        }
    }
    gs_varset_type *target = gs_fdta((char *)target_name);
    gs_varset_type *tmp = malloc(sizeof(*tmp));
    if (!tmp)
        error("Cannot allocate allele-replacement staging");
    gs_tmp_init_from(tmp, target);
    st_work_adopt(tmp, st_combine_destroy);
    int status = 0;
    gs_copy_marker_data_01_raw(tmp, source, &status);
    if (status == -2)
        error("Cannot copy marker data for allele replacement");
    free(tmp->y);
    tmp->y = NULL;
    st_marker_prediction_clear(tmp);
    for (gs_UNSIGNED m = 0; m < tmp->NoMar; m++)
        for (gs_UNSIGNED i = 0; i < tmp->NoInd; i++)
        {
            size_t at = (size_t)m * tmp->NoInd + i;
            GS_ALLELE_TYPE a = tmp->allele_raw1[at], b = tmp->allele_raw2[at];
            if (a == GS_ALLELE_MISSING)
                tmp->allele_raw2[at] = GS_ALLELE_MISSING;
            else if (selected[m] && b == 2 && (a == 1 || a == 2))
            {
                tmp->allele_raw2[at] = replace;
                if (a == 2)
                    tmp->allele_raw1[at] = replace;
            }
        }
    tmp->phased = 0;
    for (gs_UNSIGNED i = 0; i < tmp->NLoc; i++)
        tmp->map[i]->pos = st_legacy_export_precision(tmp->map[i]->pos, 12);
    char *empty = "";
    int no = 0;
    gs_marker_stats_01(tmp, &empty, &empty, &empty, &no, &no, &no, &no,
                       &status);
    if (status == -2 || gs_tmp_commit(target, tmp))
        error("Cannot commit marker allele replacement");
    return ScalarInteger(0);
}

static double st_legacy_marker_mean(const double *x, size_t n)
{
    if (!n)
        return R_NaN;
    long double sum = 0;
    for (size_t i = 0; i < n; i++)
        sum += x[i];
    double mean = (double)(sum / n);
    if (R_FINITE(mean))
    {
        sum = 0;
        for (size_t i = 0; i < n; i++)
            sum += x[i] - mean;
        mean += (double)(sum / n);
    }
    return mean;
}

static SEXP st_legacy_il_summary(SEXP args)
{
    const char *name =
        st_legacy_marker_string(getListElement(args, "source"), "data.set");
    gs_varset_type *gv = gs_fdta((char *)name);
    int donor = asInteger(getListElement(args, "donor"));
    int recipient = asInteger(getListElement(args, "recipient"));
    SEXP chromosome = getListElement(args, "chromosome");
    SEXP lower_arg = getListElement(args, "lower");
    SEXP upper_arg = getListElement(args, "upper");
    if (donor < 1 || donor > GS_ALLELE_MAX || recipient < 1 ||
        recipient > GS_ALLELE_MAX)
        error("rp.allele and dp.allele must be integers from 1 to 254");
    if (TYPEOF(chromosome) != INTSXP || !XLENGTH(chromosome) ||
        TYPEOF(lower_arg) != REALSXP || TYPEOF(upper_arg) != REALSXP)
        error("Invalid target region borders");
    int regions = INTEGER(chromosome)[0] != 0, status = 0;
    double values[10] = {0};
    il_eval_library(gv, values, values + 1, values + 2, &donor, &status);
    if (status == -2)
        return R_NilValue;
    size_t n = gv->NoInd;
    if (regions && ((uintmax_t)XLENGTH(chromosome) != (uintmax_t)n || (uintmax_t)XLENGTH(lower_arg) != (uintmax_t)n ||
                    (uintmax_t)XLENGTH(upper_arg) != (uintmax_t)n))
    {
        gs_info(gv, -2, "Incorrect length of target region borders");
        return R_NilValue;
    }
    double *absolute = st_work_alloc(n ? n : 1, sizeof(double));
    double *total = st_work_alloc(n ? n : 1, sizeof(double));
    double *segments = st_work_alloc(n ? n : 1, sizeof(double));
    double *length = st_work_alloc(n ? n : 1, sizeof(double));
    int *chrom = st_work_alloc(n ? n : 1, sizeof(int));
    double *lower = st_work_alloc(n ? n : 1, sizeof(double));
    double *upper = st_work_alloc(n ? n : 1, sizeof(double));
    memset(chrom, 0, n * sizeof(int));
    memset(lower, 0, n * sizeof(double));
    memset(upper, 0, n * sizeof(double));
    for (int pass = 0; pass < (regions ? 4 : 1); pass++)
    {
        memset(absolute, 0, n * sizeof(double));
        memset(total, 0, n * sizeof(double));
        memset(segments, 0, n * sizeof(double));
        memset(length, 0, n * sizeof(double));
        if (pass == 1)
        {
            memcpy(chrom, INTEGER(chromosome), n * sizeof(int));
            memcpy(lower, REAL(lower_arg), n * sizeof(double));
            memcpy(upper, REAL(upper_arg), n * sizeof(double));
        }
        int allele = pass == 2 ? recipient : donor, exclude = pass > 1;
        il_eval_lines(gv, &allele, chrom, lower, upper, &exclude, absolute,
                      total, segments, length, &status);
        if (status == -2)
            return R_NilValue;
        for (size_t i = 0; i < n; i++)
            absolute[i] /= total[i];
        if (!pass)
        {
            values[3] = st_legacy_marker_mean(absolute, n);
            values[4] = st_legacy_marker_mean(segments, n);
            values[5] = st_legacy_marker_mean(length, n);
        }
        else if (pass < 3)
            values[pass + 5] = st_legacy_marker_mean(absolute, n);
        else
        {
            values[8] = st_legacy_marker_mean(segments, n);
            values[9] = st_legacy_marker_mean(length, n);
        }
    }
    int count = regions ? 10 : 6;
    const char *labels[] = {"cov",  "dep",    "seglib", "dpg",     "ndps",
                            "ldps", "tr.dpg", "nt.rpg", "nt.ndps", "nt.ldps"};
    SEXP out = PROTECT(allocVector(REALSXP, count));
    SEXP names = PROTECT(allocVector(STRSXP, count));
    for (int i = 0; i < count; i++)
    {
        REAL(out)[i] = values[i];
        SET_STRING_ELT(names, i, mkChar(labels[i]));
    }
    setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

static SEXP st_legacy_marker_frame(int columns, const char **labels, int rows,
                                   int *np)
{
    SEXP out = PROTECT(allocVector(VECSXP, columns));
    SEXP names = PROTECT(allocVector(STRSXP, columns));
    SEXP rn = PROTECT(allocVector(INTSXP, rows ? 2 : 0));
    *np += 3;
    for (int i = 0; i < columns; i++)
        SET_STRING_ELT(names, i, mkChar(labels[i]));
    if (rows)
    {
        INTEGER(rn)[0] = NA_INTEGER;
        INTEGER(rn)[1] = -rows;
    }
    setAttrib(out, R_NamesSymbol, names);
    setAttrib(out, R_RowNamesSymbol, rn);
    setAttrib(out, R_ClassSymbol, mkString("data.frame"));
    return out;
}

static SEXP st_legacy_il_lines(SEXP args)
{
    const char *name =
        st_legacy_marker_string(getListElement(args, "source"), "data.set");
    gs_varset_type *gv = gs_fdta((char *)name);
    int allele = asInteger(getListElement(args, "allele"));
    int exclude = asInteger(getListElement(args, "exclude"));
    SEXP chr_arg = getListElement(args, "chromosome");
    SEXP lower_arg = getListElement(args, "lower");
    SEXP upper_arg = getListElement(args, "upper");
    if (allele < 1 || allele > GS_ALLELE_MAX)
        error("allele must be an integer from 1 to 254");
    R_xlen_t provided = XLENGTH(chr_arg);
    if (!provided || provided != XLENGTH(lower_arg) ||
        provided != XLENGTH(upper_arg))
    {
        gs_info(gv, -2, "Incorrect target region borders");
        return R_NilValue;
    }
    if (gv->NoInd > INT_MAX)
        error("Too many lines to export");
    int n = (int)gv->NoInd, status = 0;
    int *chromosome = st_work_alloc(n ? n : 1, sizeof(int));
    double *lower = st_work_alloc(n ? n : 1, sizeof(double));
    double *upper = st_work_alloc(n ? n : 1, sizeof(double));
    for (int i = 0; i < n; i++)
    {
        R_xlen_t j = provided < n ? 0 : i;
        chromosome[i] = INTEGER(chr_arg)[j];
        lower[i] = REAL(lower_arg)[j];
        upper[i] = REAL(upper_arg)[j];
    }
    const char *labels[] = {"name", "nseg", "lseg", "abs", "tot", "rel"};
    int np = 0;
    SEXP out = st_legacy_marker_frame(6, labels, n, &np);
    for (int j = 0; j < 6; j++)
    {
        SEXP col = PROTECT(allocVector(j ? REALSXP : STRSXP, n));
        np++;
        SET_VECTOR_ELT(out, j, col);
        if (j)
            memset(REAL(col), 0, n * sizeof(double));
    }
    il_eval_lines(gv, &allele, chromosome, lower, upper, &exclude,
                  REAL(VECTOR_ELT(out, 3)), REAL(VECTOR_ELT(out, 4)),
                  REAL(VECTOR_ELT(out, 1)), REAL(VECTOR_ELT(out, 2)), &status);
    if (status == -2)
    {
        UNPROTECT(np);
        return R_NilValue;
    }
    for (int i = 0; i < n; i++)
    {
        SET_STRING_ELT(VECTOR_ELT(out, 0), i, mkChar(gv->IndName[i]));
        REAL(VECTOR_ELT(out, 5))
        [i] = REAL(VECTOR_ELT(out, 3))[i] / REAL(VECTOR_ELT(out, 4))[i];
    }
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_il_create(SEXP args)
{
    double nc = asReal(getListElement(args, "chromosomes"));
    double lc = asReal(getListElement(args, "length"));
    double sl = asReal(getListElement(args, "segment"));
    int overlap = asLogical(getListElement(args, "overlap"));
    if (!R_FINITE(nc) || !R_FINITE(lc) || !R_FINITE(sl) || nc < 1 || lc < 1 ||
        sl < 1 || nc != floor(nc) || lc != floor(lc) || sl != floor(sl) ||
        nc * lc > INT_MAX || fmod(nc * lc, sl) != 0 ||
        nc * lc / sl + overlap * nc > INT_MAX)
        error(
            "Chromosome and segment sizes must be positive integers with complete segments");
    int chromosomes = (int)nc, length = (int)lc, segment = (int)sl;
    int markers = chromosomes * length;
    int individuals = markers / segment + overlap * chromosomes;
    if (individuals < 1)
        error("Library must contain at least one line");
    const char *name =
        st_legacy_marker_string(getListElement(args, "source"), "data.set");
    gs_varset_type *destination = gs_fdta((char *)name),
                   *tmp = malloc(sizeof(*tmp));
    if (!tmp)
        error("Cannot allocate introgression-library staging");
    gs_tmp_init_from(tmp, destination);
    st_work_adopt(tmp, st_combine_destroy);
    tmp->NoInd = individuals;
    tmp->NoMar = tmp->NLoc = markers;
    tmp->IndName = calloc(individuals, sizeof(*tmp->IndName));
    tmp->MarName = gs_malloc_array(markers, sizeof(*tmp->MarName));
    tmp->map_S = gs_malloc_array(markers, sizeof(*tmp->map_S));
    tmp->map = gs_malloc_array(markers, sizeof(*tmp->map));
    tmp->mdl = gs_malloc_array(markers, sizeof(*tmp->mdl));
    if (!tmp->IndName || !tmp->MarName || !tmp->map_S || !tmp->map ||
        !tmp->mdl || gs_alloc_genotype_matrix(tmp, markers, individuals))
        error("Cannot allocate introgression library");
    for (int i = 0; i < individuals; i++)
    {
        char id[64];
        snprintf(id, sizeof(id), "IL%d", i + 1);
        if (st_identity_set(tmp, i, id, NULL))
            error("Cannot initialize introgression-line identity");
        if ((int)strlen(id) > tmp->MaxIL)
            tmp->MaxIL = (int)strlen(id);
    }
    for (int c = 0, at = 0; c < chromosomes; c++)
        for (int p = 0; p < length; p++, at++)
        {
            char class_name[64];
            snprintf(tmp->MarName[at], sizeof(*tmp->MarName), "C%dMA%d", c + 1,
                     p + 1);
            snprintf(class_name, sizeof(class_name), "C%d", c + 1);
            if (gs_pool_add_string_size(&tmp->map_str_N, tmp->MarName[at]) ||
                gs_pool_add_string_size(&tmp->map_str_N, class_name))
                error("Introgression map names are too large");
            if ((int)strlen(tmp->MarName[at]) > tmp->MaxML)
                tmp->MaxML = (int)strlen(tmp->MarName[at]);
        }
    tmp->map_str_S = malloc(tmp->map_str_N);
    if (!tmp->map_str_S)
        error("Cannot allocate introgression map labels");
    size_t offset = 0;
    for (int c = 0, at = 0; c < chromosomes; c++)
        for (int p = 0; p < length; p++, at++)
        {
            char class_name[64];
            snprintf(class_name, sizeof(class_name), "C%d", c + 1);
            gs_mappoint_type point = {0};
            point.name = tmp->MarName[at];
            point.class = class_name;
            point.chrom = c + 1;
            point.pos = p + 1;
            if (gs_copy_map_point_to_pool(&tmp->map_S[at], &point, at,
                                          tmp->map_str_S, tmp->map_str_N,
                                          &offset))
                error("Cannot construct introgression map");
            tmp->map[at] = &tmp->map_S[at];
            tmp->mdl[at] = at;
            int first = at / segment + (overlap ? at / length : 0);
            for (int i = 0; i < individuals; i++)
            {
                size_t index = (size_t)at * individuals + i;
                int allele = i == first || (overlap && i == first + 1) ? 2 : 1;
                tmp->allele_raw1[index] = tmp->allele_raw2[index] = allele;
            }
        }
    int status = 0, zero = 0;
    char *empty = "";
    gs_marker_stats_01(tmp, &empty, &empty, &empty, &zero, &zero, &zero, &zero,
                       &status);
    if (status == -2 || gs_tmp_commit(destination, tmp))
        error("Cannot commit introgression library");
    return ScalarInteger(0);
}

typedef struct
{
    double key, value;
    int index;
} st_legacy_phenotype;

static int st_legacy_phenotype_compare(const void *a, const void *b)
{
    const st_legacy_phenotype *x = a, *y = b;
    int mx = ISNAN(x->key), my = ISNAN(y->key);
    if (mx != my)
        return mx - my;
    if (!mx && x->key != y->key)
        return x->key < y->key ? -1 : 1;
    return (x->index > y->index) - (x->index < y->index);
}

static double st_legacy_marker_number(SEXP input, R_xlen_t i)
{
    if (TYPEOF(input) == REALSXP)
        return REAL(input)[i];
    if (TYPEOF(input) == INTSXP || TYPEOF(input) == LGLSXP)
        return INTEGER(input)[i] == NA_INTEGER ? NA_REAL : INTEGER(input)[i];
    error("Numeric values required");
    return NA_REAL;
}

static void st_legacy_subset_cell(SEXP out, R_xlen_t to, SEXP input,
                                  R_xlen_t from)
{
    switch (TYPEOF(out))
    {
    case LGLSXP:
    case INTSXP:
        INTEGER(out)[to] = from < 0 ? NA_INTEGER : INTEGER(input)[from];
        break;
    case REALSXP:
        REAL(out)
        [to] = from < 0 ? NA_REAL : st_legacy_marker_number(input, from);
        break;
    case STRSXP:
        SET_STRING_ELT(out, to, from < 0 ? NA_STRING : STRING_ELT(input, from));
        break;
    case VECSXP:
        SET_VECTOR_ELT(out, to,
                       from < 0 ? R_NilValue : VECTOR_ELT(input, from));
        break;
    case RAWSXP:
        RAW(out)[to] = from < 0 ? 0 : RAW(input)[from];
        break;
    case CPLXSXP:
        if (from < 0)
            COMPLEX(out)[to].r = COMPLEX(out)[to].i = NA_REAL;
        else
            COMPLEX(out)[to] = COMPLEX(input)[from];
        break;
    default:
        error("Unsupported phenotype table column");
    }
}

static SEXP st_legacy_select_phenotype(SEXP args)
{
    SEXP input = getListElement(args, "data"), y = R_NilValue;
    gs_varset_type *source = NULL;
    if (isNull(input))
    {
        const char *name =
            st_legacy_marker_string(getListElement(args, "source"), "data.set");
        source = st_find_data_set_no_create(name);
        if (!source || !source->y)
            error("No performance data available for data set '%s'", name);
    }
    else
    {
        if (!inherits(input, "data.frame"))
            error("p must be a data frame");
        y = getListElement(input, "y");
        if (isNull(y))
            y = getListElement(input, "yy");
        if (isNull(getListElement(input, "i")) &&
            isNull(getListElement(input, "ii")))
            error("p must contain an individual column named 'i' or 'ii'");
        if (!isNumeric(y))
            error("performance values must be numeric");
    }
    R_xlen_t rows = source ? (R_xlen_t)source->NoInd : XLENGTH(y);
    if (rows > INT_MAX)
        error("Phenotype table is too large");
    int n = (int)rows,
        decreasing = asLogical(getListElement(args, "decreasing"));
    double threshold = asReal(getListElement(args, "threshold"));
    double keep = asReal(getListElement(args, "n"));
    if (ISNAN(threshold) || ISNAN(keep) || decreasing == NA_LOGICAL)
        error("Selection options must not be missing");
    st_legacy_phenotype *order = st_work_alloc(n ? n : 1, sizeof(*order));
    for (int i = 0; i < n; i++)
    {
        double value = source ? st_legacy_export_precision(source->y[i], 6)
                              : st_legacy_marker_number(y, i);
        order[i] = (st_legacy_phenotype){decreasing ? -value : value, value, i};
    }
    qsort(order, n, sizeof(*order), st_legacy_phenotype_compare);
    int selected = 0, missing_row = 0;
    for (int i = 0; i < n; i++)
    {
        if (threshold > 0)
        {
            if (ISNAN(order[i].value))
                order[i].index = -1, missing_row = 1;
            else if (!(decreasing ? order[i].value > threshold
                                  : order[i].value < threshold))
                continue;
        }
        order[selected++] = order[i];
    }
    if (keep > 0 && keep < selected)
        selected = (int)keep;
    int np = 0, columns = source ? 2 : LENGTH(input);
    const char *labels[] = {"i", "y"};
    SEXP out;
    if (source)
        out = st_legacy_marker_frame(2, labels, selected, &np);
    else
    {
        out = PROTECT(allocVector(VECSXP, columns));
        np++;
        copyMostAttrib(input, out);
        setAttrib(out, R_NamesSymbol, getAttrib(input, R_NamesSymbol));
    }
    for (int j = 0; j < columns; j++)
    {
        SEXP col = source ? R_NilValue : VECTOR_ELT(input, j);
        int type = source ? (j ? REALSXP : STRSXP) : TYPEOF(col);
        int matrix = !source && isMatrix(col);
        SEXP dimensions = matrix ? getAttrib(col, R_DimSymbol) : R_NilValue;
        int width = matrix ? INTEGER(dimensions)[1] : 1;
        SEXP result = PROTECT(allocVector(type, (R_xlen_t)selected * width));
        np++;
        if (!source)
        {
            if (matrix ? INTEGER(dimensions)[0] != n : XLENGTH(col) != n)
                error(
                    "Phenotype table columns must be vectors of equal length");
            copyMostAttrib(col, result);
        }
        for (int block = 0; block < width; block++)
            for (int i = 0; i < selected; i++)
            {
                int from = order[i].index;
                if (!source)
                    st_legacy_subset_cell(
                        result, i + (R_xlen_t)selected * block, col,
                        from < 0 ? -1 : from + (R_xlen_t)n * block);
                else if (j)
                    REAL(result)[i] = from < 0 ? NA_REAL : order[i].value;
                else
                    SET_STRING_ELT(result, i,
                                   from < 0 ? NA_STRING
                                            : mkChar(source->IndName[from]));
            }
        SEXP original_names =
            source ? R_NilValue : getAttrib(col, R_NamesSymbol);
        if (matrix)
        {
            SEXP dim = PROTECT(duplicate(dimensions));
            np++;
            INTEGER(dim)[0] = selected;
            setAttrib(result, R_DimSymbol, dim);
            SEXP original_dn = getAttrib(col, R_DimNamesSymbol);
            if (!isNull(original_dn))
            {
                SEXP dn = PROTECT(allocVector(VECSXP, 2));
                np++;
                SET_VECTOR_ELT(dn, 1, VECTOR_ELT(original_dn, 1));
                original_names = VECTOR_ELT(original_dn, 0);
                if (!isNull(original_names))
                {
                    SEXP labels = PROTECT(allocVector(STRSXP, selected));
                    np++;
                    for (int i = 0; i < selected; i++)
                        SET_STRING_ELT(
                            labels, i,
                            order[i].index < 0
                                ? NA_STRING
                                : STRING_ELT(original_names, order[i].index));
                    SET_VECTOR_ELT(dn, 0, labels);
                }
                setAttrib(result, R_DimNamesSymbol, dn);
            }
        }
        else if (!isNull(original_names))
        {
            SEXP labels = PROTECT(allocVector(STRSXP, selected));
            np++;
            for (int i = 0; i < selected; i++)
                SET_STRING_ELT(labels, i,
                               order[i].index < 0 ? NA_STRING
                                                  : STRING_ELT(original_names,
                                                               order[i].index));
            setAttrib(result, R_NamesSymbol, labels);
        }
        SET_VECTOR_ELT(out, j, result);
    }
    SEXP original_rows =
        source ? R_NilValue : getAttrib(input, R_RowNamesSymbol);
    int integer_rows =
        !missing_row && (source || TYPEOF(original_rows) == INTSXP);
    SEXP rownames =
        PROTECT(allocVector(integer_rows ? INTSXP : STRSXP, selected));
    np++;
    for (int i = 0; i < selected; i++)
    {
        int from = order[i].index;
        if (integer_rows)
            INTEGER(rownames)
            [i] = source ? from + 1 : INTEGER(original_rows)[from];
        else
        {
            char number[64];
            const char *name;
            if (from < 0)
                name = "NA";
            else if (!source && TYPEOF(original_rows) == STRSXP)
                name = CHAR(STRING_ELT(original_rows, from));
            else
            {
                snprintf(number, sizeof(number), "%d",
                         source ? from + 1 : INTEGER(original_rows)[from]);
                name = number;
            }
            size_t size = strlen(name) + 32;
            char *unique = st_work_alloc(size, 1);
            snprintf(unique, size, "%s", name);
            for (int suffix = 0, duplicate = 1; duplicate;)
            {
                duplicate = 0;
                for (int j = 0; j < i; j++)
                    if (!strcmp(unique, CHAR(STRING_ELT(rownames, j))))
                    {
                        snprintf(unique, size, "%s.%d", name, ++suffix);
                        duplicate = 1;
                        break;
                    }
            }
            SET_STRING_ELT(rownames, i, mkChar(unique));
        }
    }
    setAttrib(out, R_RowNamesSymbol, rownames);
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_pseff_text(SEXP column, int row, int labels)
{
    if (TYPEOF(column) == STRSXP)
        return STRING_ELT(column, row);
    if (labels && inherits(column, "factor"))
    {
        int index = INTEGER(column)[row];
        SEXP levels = getAttrib(column, R_LevelsSymbol);
        return index == NA_INTEGER ? NA_STRING : STRING_ELT(levels, index - 1);
    }
    double value = st_legacy_marker_number(column, row);
    if (ISNA(value))
        return NA_STRING;
    char text[64];
    snprintf(text, sizeof(text), "%.15g", value);
    return mkChar(text);
}

static SEXP st_legacy_pseff(SEXP args)
{
    SEXP e = getListElement(args, "data");
    if (TYPEOF(e) != VECSXP || LENGTH(e) < 4)
        error("e must be an effect data frame with four columns");
    SEXP effect = VECTOR_ELT(e, 3);
    R_xlen_t count = XLENGTH(effect);
    if (count < 2 || count > INT_MAX / 2)
        error("e must contain an intercept and marker effects");
    int n = (int)count - 1, np = 0;
    double intercept = asReal(getListElement(args, "intercept"));
    if (ISNAN(intercept))
        intercept = st_legacy_marker_number(effect, 0);
    long double sum = 0;
    for (int i = 1; i <= n; i++)
        sum += st_legacy_marker_number(effect, i);
    const char *keys[] = {"marker", "allele", "effect"};
    SEXP effects = st_legacy_marker_frame(3, keys, 2 * n, &np);
    int second_character = TYPEOF(VECTOR_ELT(e, 0)) == STRSXP ||
                           TYPEOF(VECTOR_ELT(e, 2)) == STRSXP;
    for (int j = 0; j < 3; j++)
    {
        SEXP first = VECTOR_ELT(e, j == 2 ? 3 : j);
        SEXP second = VECTOR_ELT(e, j == 1 ? 2 : j == 2 ? 3 : 0);
        int factor = j < 2 && inherits(first, "factor");
        int type = j == 2 ? REALSXP
                   : factor || TYPEOF(first) == STRSXP || second_character
                       ? STRSXP
                       : REALSXP;
        if (XLENGTH(first) != count || XLENGTH(second) != count)
            error("Effect table columns must have equal lengths");
        SEXP col = PROTECT(allocVector(type, 2 * n));
        np++;
        int invalid_factor = 0;
        for (int i = 0; i < n; i++)
        {
            if (j == 2)
            {
                REAL(col)[i] = st_legacy_marker_number(effect, i + 1) / 2;
                REAL(col)[i + n] = -REAL(col)[i];
            }
            else if (type == STRSXP)
            {
                SET_STRING_ELT(col, i, st_legacy_pseff_text(first, i + 1, 1));
                SEXP text = st_legacy_pseff_text(second, i + 1, 0);
                if (factor && !second_character && text != NA_STRING)
                {
                    SEXP levels = getAttrib(first, R_LevelsSymbol);
                    int found = 0;
                    for (R_xlen_t k = 0; k < XLENGTH(levels); k++)
                        if (!strcmp(CHAR(text), CHAR(STRING_ELT(levels, k))))
                        {
                            found = 1;
                            break;
                        }
                    if (!found)
                    {
                        text = NA_STRING;
                        invalid_factor = 1;
                    }
                }
                SET_STRING_ELT(col, i + n, text);
            }
            else
            {
                REAL(col)[i] = st_legacy_marker_number(first, i + 1);
                REAL(col)[i + n] = st_legacy_marker_number(second, i + 1);
            }
        }
        if (invalid_factor)
            warning("invalid factor level, NA generated");
        SET_VECTOR_ELT(effects, j, col);
    }
    SEXP out = PROTECT(allocVector(VECSXP, 2));
    np++;
    SET_VECTOR_ELT(out, 0, ScalarReal(intercept + (double)sum));
    SET_VECTOR_ELT(out, 1, effects);
    UNPROTECT(np);
    return out;
}

static void st_legacy_presence_error(gs_varset_type *source,
                                     const char *message)
{
    gs_info(source, -2, (char *)message);
    error("st.genetic.distances.02: %s", message);
}

void st_legacy_marker_presence(const char *dataset, int *rows, int *columns,
                               int *markers, int **counts,
                               const char ***marker_names,
                               const char ***individual_names, double **values)
{
    gs_varset_type *source = gs_fdta((char *)dataset);
    if (!st_allele_present(source) || !source->NoInd || !source->NoMar)
        st_legacy_presence_error(source, "No marker data available");
    if (source->NoInd > INT_MAX || source->NoMar > INT_MAX ||
        source->NLoc > INT_MAX)
        st_legacy_presence_error(source,
                                 "Marker data exceed supported dimensions");
    if (source->map)
    {
        gs_UNSIGNED *storage = gs_map_storage_index(source, 0, 0);
        if (!storage)
            error("st.genetic.distances.02: invalid marker map");
        free(storage);
    }
    int nc = (int)source->NoInd;
    const char **names = st_work_alloc(nc, sizeof(*names));
    gs_name_index_type *name_index = st_work_alloc(nc, sizeof(*name_index));
    for (int i = 0; i < nc; i++)
    {
        size_t n = strlen(source->IndName[i]);
        char *name = st_work_alloc(n + 1, 1);
        memcpy(name, source->IndName[i], n + 1);
        gs_str_rem_us(name);
        if (strpbrk(name, "\"'#"))
            st_legacy_presence_error(
                source,
                "Invalid quoted or commented individual name in NTSys output");
        names[i] = name;
        name_index[i].name = name;
        name_index[i].idx = i;
    }
    qsort(name_index, nc, sizeof(*name_index), gs_name_index_cmp);
    for (int i = 1; i < nc; i++)
        if (!strcmp(name_index[i - 1].name, name_index[i].name))
            st_legacy_presence_error(
                source, "Individual names collide in NTSys output");
    if (!source->All || !source->NoAll)
        st_legacy_presence_error(source, "Marker statistics unavailable");
    int available = (int)(source->map ? source->NLoc : source->NoMar);
    int *index = st_work_alloc(available, sizeof(*index));
    int *nall = st_work_alloc(available, sizeof(*nall));
    const char **marker = st_work_alloc(available, sizeof(*marker));
    int *seen = st_work_alloc(source->NoMar, sizeof(*seen));
    memset(seen, 0, source->NoMar * sizeof(*seen));
    int nr = 0, nm = 0;
    for (int j = 0; j < available; j++)
    {
        gs_UNSIGNED m = source->map ? source->mdl[j] : (gs_UNSIGNED)j;
        if (m == gs_UNSIGNED_MAX)
            continue;
        if (strchr(source->MarName[m], '.'))
            st_legacy_presence_error(
                source, "Marker name contains '.' in NTSys output");
        if (strpbrk(source->MarName[m], "\"'#"))
            st_legacy_presence_error(
                source,
                "Invalid quoted or commented marker name in NTSys output");
        int count = 0;
        for (gs_UNSIGNED a = 0; a < source->NoAll[m]; a++)
            count += source->All[m][a] != GS_ALLELE_MISSING;
        if (!count)
            continue;
        if (seen[m]++)
            error("duplicate 'row.names' are not allowed");
        if (count > INT_MAX - nr)
            st_legacy_presence_error(source,
                                     "Marker data exceed supported dimensions");
        index[nm] = (int)m;
        nall[nm] = count;
        marker[nm++] = source->MarName[m];
        nr += count;
    }
    double *data = st_work_alloc((size_t)nr * nc, sizeof(*data));
    for (int j = 0, row = 0; j < nm; j++)
    {
        int m = index[j];
        for (gs_UNSIGNED a = 0; a < source->NoAll[m]; a++)
            if (source->All[m][a] != GS_ALLELE_MISSING)
            {
                for (int i = 0; i < nc; i++)
                {
                    size_t at = (size_t)m * nc + i;
                    data[row + (size_t)nr * i] =
                        source->All[m][a] == st_allele_get(source, at, 0) ||
                        source->All[m][a] == st_allele_get(source, at, 1);
                }
                row++;
            }
    }
    *rows = nr;
    *columns = nc;
    *markers = nm;
    *counts = nall;
    *marker_names = marker;
    *individual_names = names;
    *values = data;
}

static SEXP st_buffer_get(SEXP x, const char *key)
{
    SEXP names = getAttrib(x, R_NamesSymbol);
    if (TYPEOF(x) != VECSXP || TYPEOF(names) != STRSXP)
        return R_NilValue;
    for (R_xlen_t i = 0; i < XLENGTH(x); i++)
        if (!strcmp(CHAR(STRING_ELT(names, i)), key))
            return VECTOR_ELT(x, i);
    return R_NilValue;
}
static int st_buffer_flag(SEXP x, const char *key, int fallback)
{
    SEXP value = st_buffer_get(x, key);
    return value == R_NilValue ? fallback : asLogical(value);
}
typedef struct
{
    char **value;
    size_t *offset, *width;
    size_t rows;
} st_buffer_table;
static st_buffer_table st_buffer_tokenize(const char *source, int skip)
{
    size_t n = strlen(source);
    size_t fields = 1, lines = 1;
    int space = 1;
    for (size_t i = 0; i < n; i++)
    {
        int next = isspace((unsigned char)source[i]) != 0;
        if (!next && space) fields++;
        if (source[i] == '\n') lines++;
        space = next;
    }
    char *text = st_work_alloc(n + 1, 1);
    char **value = st_work_alloc(fields, sizeof(*value));
    size_t *offset = st_work_alloc(lines, sizeof(*offset));
    size_t *width = st_work_alloc(lines, sizeof(*width));
    size_t used = 0, count = 0, rows = 0, start = 0;
    const char *p = source;
    while (skip-- > 0 && *p)
    {
        while (*p && *p != '\n') p++;
        if (*p) p++;
    }
    for (;;)
    {
        while (*p && *p != '\n' && isspace((unsigned char)*p)) p++;
        if (*p == '#')
            while (*p && *p != '\n') p++;
        if (!*p || *p == '\n')
        {
            if (count > start)
            {
                offset[rows] = start;
                width[rows++] = count - start;
                start = count;
            }
            if (!*p) break;
            p++;
            continue;
        }
        value[count++] = text + used;
        int quote = 0;
        while (*p)
        {
            unsigned char c = *p;
            if (!quote && (isspace(c) || c == '#')) break;
            if (c == '\'' || c == '"')
            {
                if (!quote) { quote = c; p++; continue; }
                if (quote == c) { quote = 0; p++; continue; }
            }
            text[used++] = *p++;
        }
        if (quote) error("Unterminated quoted output field");
        text[used++] = '\0';
    }
    st_buffer_table result = {value, offset, width, rows};
    return result;
}
static int st_buffer_missing(const char *s)
{
    return !strcmp(s, "NA");
}
static int st_buffer_logical(const char *s)
{
    if (!strcmp(s, "TRUE") || !strcmp(s, "T")) return 1;
    if (!strcmp(s, "FALSE") || !strcmp(s, "F")) return 0;
    return NA_LOGICAL;
}
static double st_buffer_number(const char *s)
{
    if (st_buffer_missing(s)) return NA_REAL;
    char *end;
    double value = R_strtod(s, &end);
    if (end == s || *end) error("Invalid numerical output field '%s'", s);
    return value;
}
static char *st_buffer_name(const char *source, int check)
{
    size_t n = strlen(source);
    char *name = st_work_alloc(2 * n + 40, 1);
    if (!check)
    {
        memcpy(name, source, n + 1);
        return name;
    }
    size_t used = 0;
    if (!n || !(isalpha((unsigned char)source[0]) ||
        (source[0] == '.' && (n < 2 || !isdigit((unsigned char)source[1]))) ||
        (unsigned char)source[0] >= 128)) name[used++] = 'X';
    for (size_t i = 0; i < n; i++)
    {
        unsigned char c = source[i];
        name[used++] = isalnum(c) || c == '.' || c == '_' || c >= 128 ? c : '.';
    }
    name[used] = '\0';
    const char *reserved[] = {"if","else","repeat","while","function","for","in","next","break","TRUE","FALSE","NULL","Inf","NaN","NA","NA_integer_","NA_real_","NA_complex_","NA_character_"};
    for (size_t i = 0; i < sizeof(reserved) / sizeof(*reserved); i++)
        if (!strcmp(name, reserved[i])) { name[used++] = '.'; name[used] = '\0'; break; }
    return name;
}
static SEXP st_buffer_result(const char *source, SEXP config)
{
    SEXP skip_arg = st_buffer_get(config, "skip");
    st_buffer_table table = st_buffer_tokenize(source, skip_arg == R_NilValue ? 0 : asInteger(skip_arg));
    if (!table.rows) error("No lines available in result");
    int header = st_buffer_flag(config, "header", 0);
    size_t width = table.width[0];
    int implicit_rows = table.rows > 1 && table.width[1] == width + 1;
    if (header == NA_LOGICAL) header = implicit_rows;
    size_t first = header ? 1 : 0;
    size_t nr = table.rows - first;
    size_t nc = width + (implicit_rows && !header ? 1 : 0);
    int row_names = implicit_rows && header;
    size_t columns = nc;
    if (row_names) columns = width;
    for (size_t r = first; r < table.rows; r++)
        if (table.width[r] != columns + row_names) error("Inconsistent result column count");
    if (nr > INT_MAX || columns > INT_MAX) error("Result dimensions exceed integer range");
    int *kind = st_work_alloc(columns, sizeof(*kind));
    char **labels = st_work_alloc(columns, sizeof(*labels));
    SEXP classes = st_buffer_get(config, "colClasses"), names_arg = st_buffer_get(config, "names");
    int check = st_buffer_flag(config, "check.names", 1);
    for (size_t c = 0; c < columns; c++)
    {
        char default_name[40]; snprintf(default_name, sizeof(default_name), "V%zu", c + 1);
        const char *label = header ? table.value[table.offset[0] + c] : default_name;
        labels[c] = st_buffer_name(label, check);
        if (check)
        {
            int suffix = 0, collision;
            char *base = labels[c];
            do
            {
                collision = 0;
                for (size_t k = 0; k < c; k++)
                    if (!strcmp(labels[k], labels[c])) { collision = 1; break; }
                if (collision)
                {
                    size_t n = strlen(base) + 40;
                    labels[c] = st_work_alloc(n, 1);
                    snprintf(labels[c], n, "%s.%d", base, ++suffix);
                }
            } while (collision);
        }
        kind[c] = LGLSXP;
        for (size_t r = 0; r < nr; r++)
        {
            const char *s = table.value[table.offset[r + first] + c + row_names];
            if (st_buffer_missing(s)) continue;
            if (st_buffer_logical(s) != NA_LOGICAL) continue;
            char *end;
            double value = R_strtod(s, &end);
            if (end == s || *end) { kind[c] = STRSXP; break; }
            if (kind[c] == LGLSXP) kind[c] = INTSXP;
            if (!R_FINITE(value) || value <= INT_MIN || value > INT_MAX || value != floor(value)) kind[c] = REALSXP;
        }
        if (TYPEOF(classes) == STRSXP && XLENGTH(classes))
        {
            const char *type = CHAR(STRING_ELT(classes, c % XLENGTH(classes)));
            if (!strcmp(type,"character")) kind[c] = STRSXP;
            else if (!strcmp(type,"numeric") || !strcmp(type,"double")) kind[c] = REALSXP;
            else if (!strcmp(type,"integer")) kind[c] = INTSXP;
            else if (!strcmp(type,"logical")) kind[c] = LGLSXP;
        }
        if (TYPEOF(names_arg) == STRSXP && c < (size_t)XLENGTH(names_arg))
            labels[c] = (char *)CHAR(STRING_ELT(names_arg, c));
    }
    SEXP matrix_arg = st_buffer_get(config, "matrix");
    int matrix = matrix_arg == R_NilValue ? 0 : asInteger(matrix_arg);
    int output_kind = LGLSXP;
    if (matrix)
        for (size_t c = 0; c < columns; c++)
        {
            int k = matrix == 1 && kind[c] == STRSXP ? INTSXP : kind[c];
            if (matrix == 1 && k == LGLSXP) k = INTSXP;
            if (k > output_kind) output_kind = k;
        }
    int np = 0;
    SEXP out;
    PROTECT(out = matrix ? allocMatrix(output_kind, nr, columns) : allocVector(VECSXP, columns)); np++;
    for (size_t c = 0; c < columns; c++)
    {
        SEXP column = out;
        if (!matrix)
        {
            PROTECT(column = allocVector(kind[c], nr)); np++;
            SET_VECTOR_ELT(out, c, column);
        }
        const char **levels = NULL;
        size_t nl = 0;
        if (matrix == 1 && kind[c] == STRSXP)
        {
            levels = st_work_alloc(nr ? nr : 1, sizeof(*levels));
            for (size_t r = 0; r < nr; r++)
            {
                const char *s = table.value[table.offset[r + first] + c + row_names];
                if (st_buffer_missing(s)) continue;
                size_t k = 0; while (k < nl && strcmp(levels[k], s)) k++;
                if (k == nl) levels[nl++] = s;
            }
            for (size_t i = 1; i < nl; i++)
            {
                const char *s = levels[i]; size_t k = i;
                while (k && strcoll(levels[k - 1], s) > 0) { levels[k] = levels[k - 1]; k--; }
                levels[k] = s;
            }
        }
        for (size_t r = 0; r < nr; r++)
        {
            const char *s = table.value[table.offset[r + first] + c + row_names];
            size_t at = matrix ? r + nr * c : r;
            int k = matrix ? output_kind : kind[c];
            if (k == STRSXP)
                SET_STRING_ELT(column, at, st_buffer_missing(s) ? NA_STRING : mkChar(s));
            else
            {
                double value;
                if (st_buffer_missing(s)) value = NA_REAL;
                else if (levels)
                {
                    size_t j = 0; while (j < nl && strcmp(levels[j], s)) j++;
                    value = j + 1;
                }
                else if (kind[c] == LGLSXP) value = st_buffer_logical(s);
                else value = st_buffer_number(s);
                if (k == REALSXP) REAL(column)[at] = value;
                else if (k == INTSXP) INTEGER(column)[at] = ISNAN(value) ? NA_INTEGER : (int)value;
                else LOGICAL(column)[at] = ISNAN(value) ? NA_LOGICAL : (int)value;
            }
        }
    }
    SEXP labels_out = R_NilValue, rows_out = R_NilValue;
    int keep_names = st_buffer_flag(config, "keep.names", 1);
    if (!matrix || keep_names)
    {
        PROTECT(labels_out = allocVector(STRSXP, columns)); np++;
        for (size_t c = 0; c < columns; c++) SET_STRING_ELT(labels_out, c, mkChar(labels[c]));
    }
    if (row_names)
    {
        PROTECT(rows_out = allocVector(STRSXP, nr)); np++;
        for (size_t r = 0; r < nr; r++)
        {
            const char *s = table.value[table.offset[r + first]];
            if (st_buffer_missing(s)) error("Missing result row name");
            for (size_t k = 0; k < r; k++)
                if (!strcmp(s, CHAR(STRING_ELT(rows_out, k)))) error("Duplicate result row name");
            SET_STRING_ELT(rows_out, r, mkChar(s));
        }
    }
    if (matrix)
    {
        if (matrix)
        {
            SEXP dn; PROTECT(dn = allocVector(VECSXP, 2)); np++;
            SET_VECTOR_ELT(dn, 0, rows_out); SET_VECTOR_ELT(dn, 1, labels_out);
            setAttrib(out, R_DimNamesSymbol, dn);
        }
    }
    else
    {
        setAttrib(out, R_NamesSymbol, labels_out);
        if (!row_names)
        {
            PROTECT(rows_out = allocVector(INTSXP, 2)); np++;
            INTEGER(rows_out)[0] = NA_INTEGER; INTEGER(rows_out)[1] = -(int)nr;
        }
        setAttrib(out, R_RowNamesSymbol, rows_out);
        classgets(out, mkString("data.frame"));
    }
    UNPROTECT(np);
    return out;
}
static void st_buffer_input(st_io_stream *stream, SEXP input, int header)
{
    int frame = TYPEOF(input) == VECSXP;
    R_xlen_t nc = frame ? XLENGTH(input) : 1;
    R_xlen_t nr = frame ? (nc ? XLENGTH(VECTOR_ELT(input,0)) : 0) : XLENGTH(input);
    if (header && frame)
    {
        SEXP names = getAttrib(input, R_NamesSymbol);
        for (R_xlen_t c = 0; c < nc; c++)
            st_io_printf(stream, "%s%s", c ? " " : "", CHAR(STRING_ELT(names,c)));
        st_io_printf(stream, "\n");
    }
    for (R_xlen_t r = 0; r < nr; r++)
    {
        for (R_xlen_t c = 0; c < nc; c++)
        {
            SEXP column = frame ? VECTOR_ELT(input,c) : input;
            if (r >= XLENGTH(column)) error("Inconsistent input column length");
            if (c) st_io_printf(stream, " ");
            if (Rf_isFactor(column))
            {
                int at = INTEGER(column)[r];
                SEXP levels = getAttrib(column,R_LevelsSymbol);
                st_io_printf(stream,"%s",at == NA_INTEGER ? "NA" : CHAR(STRING_ELT(levels,at-1)));
            }
            else if (TYPEOF(column) == STRSXP)
                st_io_printf(stream, "%s", STRING_ELT(column,r) == NA_STRING ? "NA" : CHAR(STRING_ELT(column,r)));
            else if (TYPEOF(column) == REALSXP)
            {
                if (ISNAN(REAL(column)[r])) st_io_printf(stream,"NA");
                else st_io_printf(stream, "%.15g", REAL(column)[r]);
            }
            else if (TYPEOF(column) == INTSXP || TYPEOF(column) == LGLSXP)
            {
                int v = TYPEOF(column) == INTSXP ? INTEGER(column)[r] : LOGICAL(column)[r];
                if (v == NA_INTEGER) st_io_printf(stream,"NA");
                else st_io_printf(stream,"%d",v);
            }
            else error("Unsupported input column type");
        }
        st_io_printf(stream, "\n");
    }
}
static void st_legacy_foreign_validate(SEXP value, int argument, int naok)
{
    R_xlen_t n = XLENGTH(value);
    switch (TYPEOF(value))
    {
    case INTSXP:
    case LGLSXP:
        if (!naok)
            for (R_xlen_t i = 0; i < n; i++)
                if ((TYPEOF(value) == INTSXP ? INTEGER(value)[i] : LOGICAL(value)[i]) == NA_INTEGER)
                    error("NAs in foreign function call (arg %d)", argument);
        break;
    case REALSXP:
        if (!naok)
            for (R_xlen_t i = 0; i < n; i++)
                if (!R_FINITE(REAL(value)[i]))
                    error("NA/NaN/Inf in foreign function call (arg %d)", argument);
        break;
    case CPLXSXP:
        if (!naok)
            for (R_xlen_t i = 0; i < n; i++)
                if (!R_FINITE(COMPLEX(value)[i].r) || !R_FINITE(COMPLEX(value)[i].i))
                    error("NA/NaN/Inf in foreign function call (arg %d)", argument);
        break;
    case STRSXP:
    case RAWSXP:
    case NILSXP:
        break;
    default:
        error("Unsupported legacy interface argument type");
    }
}

static SEXP st_legacy_result_call(SEXP args)
{
    SEXP input = st_buffer_get(args, "arguments"),
         routine = st_buffer_get(args, "routine"),
         selected = st_buffer_get(args, "result");
    if (TYPEOF(input) != VECSXP || XLENGTH(input) > 17 ||
        TYPEOF(routine) != STRSXP || XLENGTH(routine) != 1 ||
        STRING_ELT(routine, 0) == NA_STRING)
        error("Invalid legacy interface request");
    int count = (int)XLENGTH(input), full = isNull(selected),
        slot = full ? 0 : asInteger(selected),
        naok = st_buffer_flag(args, "naok", 0),
        logical = st_buffer_flag(args, "logical", 0), np = 0;
    if (!full && (slot < 0 || slot > count || slot == NA_INTEGER))
        error("Invalid legacy result position");
    if (logical && (full || !slot))
        error("A logical legacy result needs one selected position");
    for (int i = 0; i < count; i++)
        st_legacy_foreign_validate(VECTOR_ELT(input, i), i + 1, naok);
    SEXP result = R_NilValue;
    if (full)
    {
        PROTECT(result = allocVector(VECSXP, count)); np++;
        setAttrib(result, R_NamesSymbol, getAttrib(input, R_NamesSymbol));
    }
    void **value = st_work_alloc(count ? count : 1, sizeof(void *));
    for (int i = 0; i < count; i++)
    {
        SEXP source = VECTOR_ELT(input, i), output = R_NilValue;
        int retain = full || (i + 1 == slot && !logical);
        size_t n = (size_t)XLENGTH(source), bytes = 0;
        const void *data = NULL;
        if (retain)
        {
            PROTECT(output = duplicate(source)); np++;
            if (full) SET_VECTOR_ELT(result, i, output);
            else result = output;
        }
        if (!n)
        {
            value[i] = NULL;
            continue;
        }
        switch (TYPEOF(source))
        {
        case INTSXP:
        case LGLSXP:
            bytes = sizeof(int);
            data = TYPEOF(source) == INTSXP ? INTEGER(source) : LOGICAL(source);
            if (retain)
                value[i] = TYPEOF(output) == INTSXP ? INTEGER(output) : LOGICAL(output);
            break;
        case REALSXP:
            bytes = sizeof(double); data = REAL(source);
            if (retain) value[i] = REAL(output);
            break;
        case CPLXSXP:
            bytes = sizeof(Rcomplex); data = COMPLEX(source);
            if (retain) value[i] = COMPLEX(output);
            break;
        case RAWSXP:
            bytes = sizeof(Rbyte); data = RAW(source);
            if (retain) value[i] = RAW(output);
            break;
        case STRSXP:
        {
            char **strings = st_work_alloc(n ? n : 1, sizeof(char *));
            for (size_t j = 0; j < n; j++)
            {
                const char *text = translateChar(STRING_ELT(source, j));
                strings[j] = st_work_alloc(strlen(text) + 1, 1);
                strcpy(strings[j], text);
            }
            value[i] = strings;
            break;
        }
        case NILSXP:
            value[i] = NULL;
            break;
        }
        if (bytes && !retain)
        {
            value[i] = st_work_alloc(n ? n : 1, bytes);
            if (n) memcpy(value[i], data, n * bytes);
        }
    }
    st_legacy_C_call(CHAR(STRING_ELT(routine, 0)), value, count);
    if (logical)
    {
        SEXP source = VECTOR_ELT(input, slot - 1);
        R_xlen_t n = XLENGTH(source);
        PROTECT(result = allocVector(LGLSXP, n)); np++;
        for (R_xlen_t i = 0; i < n; i++)
        {
            double v;
            if (TYPEOF(source) == INTSXP || TYPEOF(source) == LGLSXP)
            {
                int integer = ((const int *)value[slot - 1])[i];
                v = integer == NA_INTEGER ? NA_REAL : integer;
            }
            else if (TYPEOF(source) == REALSXP)
                v = ((const double *)value[slot - 1])[i];
            else
                error("A logical legacy result needs a numeric argument");
            LOGICAL(result)[i] = ISNAN(v) ? NA_LOGICAL : v != 0;
        }
    }
    else
        for (int i = 0; i < count; i++)
            if ((full || i + 1 == slot) && TYPEOF(VECTOR_ELT(input, i)) == STRSXP)
            {
                SEXP output = full ? VECTOR_ELT(result, i) : result;
                const char *const *strings = value[i];
                for (R_xlen_t j = 0; j < XLENGTH(output); j++)
                    SET_STRING_ELT(output, j, mkChar(strings[j]));
            }
    UNPROTECT(np);
    return result;
}

static SEXP st_legacy_map_write(SEXP args)
{
    SEXP files = st_buffer_get(args, "file"),
         datasets = st_buffer_get(args, "dataset"),
         environment = st_buffer_get(args, "environment");
    if (TYPEOF(files) != STRSXP || !XLENGTH(files) ||
        TYPEOF(datasets) != STRSXP || !XLENGTH(datasets) ||
        TYPEOF(environment) != ENVSXP)
        error("Invalid map-output interface request");
    char *path = NULL, *dataset = NULL;
    SEXP inputs[] = {files, datasets};
    char **destinations[] = {&path, &dataset};
    for (int k = 0; k < 2; k++)
        for (R_xlen_t j = 0; j < XLENGTH(inputs[k]); j++)
        {
            const char *text = translateChar(STRING_ELT(inputs[k], j));
            if (!j)
            {
                *destinations[k] = st_work_alloc(strlen(text) + 1, 1);
                strcpy(*destinations[k], text);
            }
        }
    int status = 0;
    gs_varset_type *source = gs_fdta(dataset);
    gs_write_map_01(source, &path, &status);
    int np = 0;
    SEXP available = PROTECT(ScalarLogical(status != -2)); np++;
    SEXP condition = PROTECT(lang3(install("&&"), install("auxfiles"), available)); np++;
    SEXP yes = PROTECT(ScalarLogical(TRUE)); np++;
    SEXP no = PROTECT(ScalarLogical(FALSE)); np++;
    SEXP choice = PROTECT(lang4(install("if"), condition, yes, no)); np++;
    SEXP enabled = PROTECT(eval(choice, environment)); np++;
    if (asLogical(enabled) != TRUE)
    {
        UNPROTECT(np);
        return R_NilValue;
    }
    if (XLENGTH(files) != 1)
        error("invalid 'description' argument");
    st_io_begin();
    char *memory = (char *)st_io_memory("map-output");
    gs_write_map_01(source, &memory, &status);
    const char *text = st_io_text(memory);
    if (!*text) error("no lines available in input");
    SEXP result = st_buffer_result(text, st_buffer_get(args, "output"));
    st_io_end();
    UNPROTECT(np);
    return result;
}

static SEXP st_legacy_buffer_call(SEXP args)
{
    SEXP routine = st_buffer_get(args,"routine"), input = st_buffer_get(args,"arguments");
    if (TYPEOF(routine) != STRSXP || XLENGTH(routine) != 1 || STRING_ELT(routine,0) == NA_STRING || TYPEOF(input) != VECSXP || XLENGTH(input) > 17)
        error("Invalid memory interface request");
    int count = XLENGTH(input);
    void **value = st_work_alloc(17,sizeof(*value));
    for (int i = 0; i < count; i++)
    {
        SEXP v = VECTOR_ELT(input,i);
        size_t n = v == R_NilValue ? 0 : XLENGTH(v);
        if (!n)
        {
            if (v != R_NilValue && TYPEOF(v) != INTSXP && TYPEOF(v) != LGLSXP &&
                TYPEOF(v) != REALSXP && TYPEOF(v) != STRSXP)
                error("Unsupported memory interface argument");
            value[i] = NULL;
            continue;
        }
        if (TYPEOF(v) == INTSXP || TYPEOF(v) == LGLSXP)
        {
            int *p = st_work_alloc(n ? n : 1,sizeof(int));
            if (n) memcpy(p, TYPEOF(v) == INTSXP ? INTEGER(v) : LOGICAL(v), n*sizeof(int));
            value[i] = p;
        }
        else if (TYPEOF(v) == REALSXP)
        {
            double *p = st_work_alloc(n ? n : 1,sizeof(double));
            if (n) memcpy(p,REAL(v),n*sizeof(double));
            value[i] = p;
        }
        else if (TYPEOF(v) == STRSXP)
        {
            char **p = st_work_alloc(n ? n : 1,sizeof(*p));
            for (size_t j = 0; j < n; j++)
            {
                const char *s = CHAR(STRING_ELT(v,j));
                p[j] = st_work_alloc(strlen(s)+1,1); strcpy(p[j],s);
            }
            value[i] = p;
        }
        else if (v == R_NilValue) value[i] = NULL;
        else error("Unsupported memory interface argument");
    }
    SEXP retval = st_buffer_get(args,"retval");
    int status_at = -1;
    if (retval != R_NilValue)
    {
        int at = asInteger(retval);
        if (at == NA_INTEGER || at < 1 || at > count)
            error("Invalid memory status position");
        SEXP status = VECTOR_ELT(input,at-1);
        if ((TYPEOF(status) != INTSXP && TYPEOF(status) != LGLSXP) ||
            XLENGTH(status) < 1)
            error("Memory status argument must contain an integer");
        status_at = at-1;
    }
    st_io_begin();
    SEXP outputs = st_buffer_get(args,"outputs"), inputs = st_buffer_get(args,"inputs");
    if ((outputs != R_NilValue && TYPEOF(outputs) != VECSXP) ||
        (inputs != R_NilValue && TYPEOF(inputs) != VECSXP))
        error("Invalid memory bindings");
    size_t no = outputs == R_NilValue ? 0 : XLENGTH(outputs);
    const char **paths = st_work_alloc(no ? no : 1,sizeof(*paths));
    for (size_t i = 0; i < no; i++)
    {
        SEXP config = VECTOR_ELT(outputs,i);
        int at = asInteger(st_buffer_get(config,"argument"));
        if (at == NA_INTEGER || at < 1 || at > count || at-1 == status_at)
            error("Invalid output binding");
        at--;
        char name[40]; snprintf(name,sizeof(name),"output-%zu",i);
        paths[i] = st_io_memory(name);
        char **slot = st_work_alloc(1,sizeof(*slot)); *slot = (char *)paths[i]; value[at] = slot;
    }
    for (R_xlen_t i = 0; i < (inputs == R_NilValue ? 0 : XLENGTH(inputs)); i++)
    {
        SEXP config = VECTOR_ELT(inputs,i);
        int at = asInteger(st_buffer_get(config,"argument"));
        if (at == NA_INTEGER || at < 1 || at > count || at-1 == status_at)
            error("Invalid input binding");
        at--;
        char name[40]; snprintf(name,sizeof(name),"input-%lld",(long long)i);
        const char *path = st_io_memory(name);
        char **slot = st_work_alloc(1,sizeof(*slot)); *slot = (char *)path; value[at] = slot;
        st_io_stream *stream = st_io_open(path,"w");
        st_buffer_input(stream,st_buffer_get(config,"data"),st_buffer_flag(config,"header",0));
        st_io_close(stream);
    }
    st_legacy_C_call(CHAR(STRING_ELT(routine,0)),value,count);
    int status = status_at < 0 ? 0 : *(int *)value[status_at];
    if (status > 0 && st_buffer_flag(args,"warn.nonconvergence",0))
        warning("Variance-component estimation did not converge");
    if (status == -2 || !st_buffer_flag(args,"enabled",1))
    {
        st_io_end();
        return R_NilValue;
    }
    size_t active = 0;
    for (size_t i = 0; i < no; i++)
        active += st_buffer_flag(VECTOR_ELT(outputs,i),"enabled",1) != 0;
    if (!active)
    {
        st_io_end();
        return R_NilValue;
    }
    if (no == 1 && !st_buffer_flag(args,"list",0))
    {
        SEXP result = st_buffer_result(st_io_text(paths[0]), VECTOR_ELT(outputs,0));
        st_io_end();
        return result;
    }
    int np = 0;
    SEXP out; PROTECT(out = allocVector(VECSXP,active)); np++;
    SEXP names = st_buffer_get(args,"output.names"), final_names = R_NilValue;
    if (names != R_NilValue)
    {
        PROTECT(final_names = allocVector(STRSXP,active)); np++;
        setAttrib(out,R_NamesSymbol,final_names);
    }
    size_t at = 0;
    for (size_t i = 0; i < no; i++)
        if (st_buffer_flag(VECTOR_ELT(outputs,i),"enabled",1))
        {
            SET_VECTOR_ELT(out,at,st_buffer_result(st_io_text(paths[i]),VECTOR_ELT(outputs,i)));
            if (names != R_NilValue) SET_STRING_ELT(final_names,at,STRING_ELT(names,i));
            at++;
        }
    st_io_end();
    UNPROTECT(np);
    return out;
}

static SEXP st_legacy_hblocks(SEXP args)
{
    SEXP markers = getListElement(args,"markers");
    if (TYPEOF(markers) != STRSXP) error("Markers must be character strings");
    size_t total = 0;
    for (R_xlen_t i = 0; i < XLENGTH(markers); i++)
    {
        const char *s = CHAR(STRING_ELT(markers,i));
        if (*s) total++;
        for (const char *p = s; *p; p++) if (*p == ';' && p[1]) total++;
    }
    if (total > INT_MAX) error("Too many haplotype markers");
    int *index = st_work_alloc(total ? total : 1,sizeof(int)); size_t at = 0;
    for (R_xlen_t i = 0; i < XLENGTH(markers); i++)
    {
        const char *s = CHAR(STRING_ELT(markers,i));
        if (*s) index[at++] = i + 1;
        for (const char *p = s; *p; p++) if (*p == ';' && p[1]) index[at++] = i + 1;
    }
    int n = total, status = 0, aux = asLogical(getListElement(args,"auxfiles"));
    char *symbol = (char *)st_legacy_marker_string(getListElement(args,"symbol"),"hap.symbol");
    char *dataset = (char *)st_legacy_marker_string(getListElement(args,"dataset"),"data.set");
    st_io_begin(); char *path = (char *)st_io_memory("blocks");
    gs_set_hblocks_01_GV(&n,index,&symbol,&path,&aux,&status,&dataset);
    SEXP result = status == -2 || !aux ? R_NilValue : st_buffer_result(st_io_text(path),getListElement(args,"output"));
    st_io_end(); return result;
}
static double st_cross_scalar(SEXP value, const char *name)
{
    if (XLENGTH(value) != 1 || (TYPEOF(value) != REALSXP && TYPEOF(value) != INTSXP))
        error("%s must be one finite number", name);
    double out = asReal(value);
    if (!R_FINITE(out)) error("%s must be one finite number", name);
    return out;
}

static gs_varset_type *st_cross_dataset(SEXP value)
{
    const char *name = st_legacy_marker_string(value, "data.set");
    if (!strcmp(name, "default")) return gv;
    for (gs_UNSIGNED i = 0; i < gs_ndta; i++)
        if (!strcmp(name, gs_data[i]->name)) return gs_data[i];
    error("Marker data set '%s' does not exist", name);
    return NULL;
}

double st_cross_propsd_validate(SEXP matrix, int nparents, const char **ids,
                                const int *slot_counts, double maximum);

static double st_cross_intensity(double alpha, double n, double g)
{
    double intensity = 0;
    if (alpha < 1)
    {
        intensity = exp(dnorm4(qnorm5(alpha, 0, 1, 0, 0), 0, 1, 1) - log(alpha));
        if (n > 0) intensity -= (1 - n / g) / (2 * n * (1 + 1 / g) * intensity);
    }
    if (!R_FINITE(intensity) || intensity < 0)
        error("The selection request has no finite nonnegative intensity");
    return intensity;
}

static SEXP st_allocate_keep(SEXP value, int *np)
{
    PROTECT(value); (*np)++; return value;
}

static SEXP st_allocate_list(int n, const char **names, int *np)
{
    SEXP out = st_allocate_keep(allocVector(VECSXP, n), np),
         labels = st_allocate_keep(allocVector(STRSXP, n), np);
    for (int i = 0; i < n; i++) SET_STRING_ELT(labels, i, mkChar(names[i]));
    setAttrib(out, R_NamesSymbol, labels);
    return out;
}

static int st_allocate_numeric(SEXP value)
{
    return (TYPEOF(value) == REALSXP || TYPEOF(value) == INTSXP) &&
        !inherits(value, "factor") && isNull(getAttrib(value, R_DimSymbol));
}

static int st_allocate_count(SEXP value, int row, int n, const char *name)
{
    if (!st_allocate_numeric(value) || (XLENGTH(value) != 1 && XLENGTH(value) != n))
        error("st.cross.allocate: %s must be numeric with length one or the number of input rows", name);
    double x = st_legacy_marker_number(value, XLENGTH(value) == 1 ? 0 : row);
    if (!R_FINITE(x) || x < 1 || x > INT_MAX || x != floor(x))
        error("st.cross.allocate: %s must contain positive supported whole counts", name);
    return (int)x;
}

static double st_allocate_accuracy(SEXP value, int row, int n)
{
    if (!st_allocate_numeric(value) || (XLENGTH(value) != 1 && XLENGTH(value) != n))
        error("st.cross.allocate: accuracy must be numeric with length one or the number of input rows");
    double x = st_legacy_marker_number(value, XLENGTH(value) == 1 ? 0 : row);
    if (!R_FINITE(x) || x < 0 || x > 1)
        error("st.cross.allocate: accuracy must be between zero and one");
    return x;
}

static const char *st_allocate_id(SEXP value, int i)
{
    if (STRING_ELT(value, i) == NA_STRING || !CHAR(STRING_ELT(value, i))[0])
        error("st.cross.allocate: parental IDs must be complete nonempty strings");
    return translateCharUTF8(STRING_ELT(value, i));
}

static int st_allocate_find(const char *name, const char **ids, int n)
{
    for (int i = 0; i < n; i++) if (!strcmp(name, ids[i])) return i;
    return -1;
}

static double st_allocate_response(int retained, int produced, double scale)
{
    double out = scale * st_cross_intensity((double)retained / produced, retained, produced);
    if (!R_FINITE(out)) error("st.cross.allocate: a family response is not finite");
    return out;
}

static SEXP st_legacy_cross_allocate(SEXP args)
{
    SEXP plan = getListElement(args, "plan"), scores = getListElement(args, "scores"),
         control = getListElement(args, "control");
    int inherited = inherits(plan, "st.cross.plan");
    SEXP table = inherited ? getListElement(plan, "crosses") : plan;
    if (!inherits(table, "data.frame"))
        error("st.cross.allocate: plan must be an st.cross.plan result or a data frame");
    SEXP first = getListElement(table, "P1Name"), second = getListElement(table, "P2Name");
    if (TYPEOF(first) != STRSXP || TYPEOF(second) != STRSXP ||
        XLENGTH(first) < 1 || XLENGTH(first) > INT_MAX || XLENGTH(second) != XLENGTH(first))
        error("st.cross.allocate: plan needs nonempty P1Name and P2Name character columns");
    int nr = (int)XLENGTH(first);
    if (!st_allocate_numeric(getListElement(args, "total.progeny")))
        error("st.cross.allocate: total.progeny must be a numeric scalar");
    double budget_value = st_cross_scalar(getListElement(args, "total.progeny"), "total.progeny");
    if (budget_value < 1 || budget_value > INT_MAX || budget_value != floor(budget_value))
        error("st.cross.allocate: total.progeny must be a positive supported whole count");
    int budget = (int)budget_value;
    double max_cells = 1e7, max_transitions = 1e8, tolerance = 1e-10;
    if (TYPEOF(control) != VECSXP || Rf_isObject(control))
        error("st.cross.allocate: control must be a named list");
    SEXP cnames = getAttrib(control, R_NamesSymbol);
    if (XLENGTH(control) && (TYPEOF(cnames) != STRSXP || XLENGTH(cnames) != XLENGTH(control)))
        error("st.cross.allocate: control must be a named list");
    unsigned int seen = 0;
    for (R_xlen_t i = 0; i < XLENGTH(control); i++)
    {
        if (STRING_ELT(cnames, i) == NA_STRING) error("st.cross.allocate: missing control name");
        const char *name = CHAR(STRING_ELT(cnames, i));
        int key = !strcmp(name, "max.cells") ? 0 : !strcmp(name, "max.transitions") ? 1 :
            !strcmp(name, "tolerance") ? 2 : -1;
        if (key < 0 || (seen & (1u << key)))
            error("st.cross.allocate: unknown or repeated control '%s'", name);
        seen |= 1u << key;
        if (!st_allocate_numeric(VECTOR_ELT(control, i)))
            error("st.cross.allocate: control values must be numeric scalars");
        double x = st_cross_scalar(VECTOR_ELT(control, i), name);
        if (key < 2 && (x < 1 || x != floor(x) || x > 9007199254740991.0))
            error("st.cross.allocate: resource limits must be positive whole counts at most 2^53-1");
        if (key == 2 && x < 0) error("st.cross.allocate: tolerance must be nonnegative");
        if (key == 0) max_cells = x;
        if (key == 1) max_transitions = x;
        if (key == 2) tolerance = x;
    }
    SEXP retained_arg = getListElement(args, "n.selected"), lower_arg = getListElement(args, "min.progeny"),
         upper_arg = getListElement(args, "max.progeny"), accuracy_arg = getListElement(args, "accuracy");
    int *p1 = st_work_alloc(nr, sizeof(int)), *p2 = st_work_alloc(nr, sizeof(int)),
        *retained = st_work_alloc(nr, sizeof(int)), *lower = st_work_alloc(nr, sizeof(int)),
        *upper = st_work_alloc(nr, sizeof(int)), *produced = st_work_alloc(nr, sizeof(int)),
        *input_family = st_work_alloc(nr, sizeof(int));
    double *accuracy = st_work_alloc(nr, sizeof(double)), *mu = st_work_alloc(nr, sizeof(double)),
        *va = st_work_alloc(nr, sizeof(double));
    const char **ids;
    int nparents = 0;
    SEXP inputs = inherited ? getListElement(plan, "inputs") : R_NilValue,
         matrix = inherited ? getListElement(inputs, "K") : R_NilValue,
         constraints_in = inherited ? getListElement(plan, "constraints") : R_NilValue;
    if (inherited)
    {
        SEXP parent_names = getListElement(getListElement(plan, "parents"), "parent");
        if (TYPEOF(parent_names) != STRSXP || XLENGTH(parent_names) < 2 || XLENGTH(parent_names) > INT_MAX)
            error("st.cross.allocate: the inherited parent table is invalid");
        nparents = (int)XLENGTH(parent_names);
        ids = st_work_alloc(nparents, sizeof(*ids));
        for (int i = 0; i < nparents; i++)
        {
            ids[i] = st_allocate_id(parent_names, i);
            if (st_allocate_find(ids[i], ids, i) >= 0)
                error("st.cross.allocate: inherited parental IDs must be unique");
        }
    }
    else ids = st_work_alloc((size_t)nr * 2, sizeof(*ids));
    int nf = 0;
    for (int i = 0; i < nr; i++)
    {
        const char *a = st_allocate_id(first, i), *b = st_allocate_id(second, i);
        if (!strcmp(a, b)) error("st.cross.allocate: self crosses are not supported");
        int ia = st_allocate_find(a, ids, nparents), ib;
        if (ia < 0)
        {
            if (inherited) error("st.cross.allocate: an input pair has an unknown parent");
            if (nparents == INT_MAX) error("st.cross.allocate: too many distinct parents");
            ia = nparents; ids[nparents++] = a;
        }
        ib = st_allocate_find(b, ids, nparents);
        if (ib < 0)
        {
            if (inherited) error("st.cross.allocate: an input pair has an unknown parent");
            if (nparents == INT_MAX) error("st.cross.allocate: too many distinct parents");
            ib = nparents; ids[nparents++] = b;
        }
        int r = st_allocate_count(retained_arg, i, nr, "n.selected"),
            lo = st_allocate_count(lower_arg, i, nr, "min.progeny"),
            hi = st_allocate_count(upper_arg, i, nr, "max.progeny");
        double a0 = st_allocate_accuracy(accuracy_arg, i, nr);
        if (lo < r || hi < lo)
            error("st.cross.allocate: bounds must satisfy n.selected <= min.progeny <= max.progeny in every input row");
        int e = 0;
        while (e < nf && !((p1[e] == ia && p2[e] == ib) || (p1[e] == ib && p2[e] == ia))) e++;
        if (e == nf)
        {
            p1[e] = ia; p2[e] = ib; retained[e] = lower[e] = upper[e] = 0;
            accuracy[e] = a0; nf++;
        }
        else if (accuracy[e] != a0)
            error("st.cross.allocate: repeated unordered pairs require the same accuracy");
        if (r > INT_MAX - retained[e] || lo > INT_MAX - lower[e] || hi > INT_MAX - upper[e])
            error("st.cross.allocate: pooled family counts exceed the supported integer range");
        retained[e] += r; lower[e] += lo; upper[e] += hi; input_family[i] = e;
        if (!(i % 1024)) R_CheckUserInterrupt();
    }
    int64_t sum_lower = 0, sum_upper = 0, sum_retained = 0;
    for (int e = 0; e < nf; e++)
    {
        sum_lower += lower[e]; sum_upper += upper[e]; sum_retained += retained[e];
    }
    if (sum_lower > budget || sum_upper < budget)
        error("st.cross.allocate: total.progeny is infeasible for the pooled family bounds");
    if (!inherits(scores, "data.frame")) error("st.cross.allocate: scores must be a data frame");
    SEXP sf = getListElement(scores, "P1Name"), ss = getListElement(scores, "P2Name"),
         sm = getListElement(scores, "mu"), sv = getListElement(scores, "va");
    if (TYPEOF(sf) != STRSXP || TYPEOF(ss) != STRSXP || XLENGTH(sf) > INT_MAX ||
        XLENGTH(ss) != XLENGTH(sf) || !st_allocate_numeric(sm) || !st_allocate_numeric(sv) ||
        XLENGTH(sm) != XLENGTH(sf) || XLENGTH(sv) != XLENGTH(sf))
        error("st.cross.allocate: scores need matching character P1Name/P2Name and numeric mu/va columns");
    int nscores = (int)XLENGTH(sf), *found = st_work_alloc(nf, sizeof(int));
    memset(found, 0, (size_t)nf * sizeof(int));
    for (int i = 0; i < nscores; i++)
    {
        const char *a = st_allocate_id(sf, i), *b = st_allocate_id(ss, i);
        double m = st_legacy_marker_number(sm, i), v = st_legacy_marker_number(sv, i);
        if (!R_FINITE(m) || !R_FINITE(v) || v < 0)
            error("st.cross.allocate: scores must have finite mu and nonnegative finite va");
        for (int e = 0; e < nf; e++)
            if ((!strcmp(a, ids[p1[e]]) && !strcmp(b, ids[p2[e]])) ||
                (!strcmp(a, ids[p2[e]]) && !strcmp(b, ids[p1[e]])))
            {
                if (found[e]) error("st.cross.allocate: scores contain a repeated planned unordered pair");
                found[e] = 1; mu[e] = m; va[e] = v;
            }
        if (!(i % 1024)) R_CheckUserInterrupt();
    }
    for (int e = 0; e < nf; e++) if (!found[e])
        error("st.cross.allocate: scores must cover every pooled planned pair");
    double *cret = st_work_alloc(nparents, sizeof(double)), *cprod = st_work_alloc(nparents, sizeof(double));
    memset(cret, 0, (size_t)nparents * sizeof(double));
    memset(cprod, 0, (size_t)nparents * sizeof(double));
    for (int e = 0; e < nf; e++)
    {
        double c = retained[e] / (2.0 * sum_retained);
        cret[p1[e]] += c; cret[p2[e]] += c;
    }
    double qret = NA_REAL, qprod = NA_REAL;
    if (inherited)
    {
        SEXP dim = getAttrib(matrix, R_DimSymbol);
        if ((TYPEOF(matrix) != REALSXP && TYPEOF(matrix) != INTSXP) || TYPEOF(dim) != INTSXP ||
            XLENGTH(dim) != 2 || INTEGER(dim)[0] != nparents || INTEGER(dim)[1] != nparents)
            error("st.cross.allocate: inherited K has invalid dimensions");
        SEXP dn = getAttrib(matrix, R_DimNamesSymbol);
        if (!isNull(dn))
        {
            if (TYPEOF(dn) != VECSXP || XLENGTH(dn) != 2)
                error("st.cross.allocate: inherited K names are invalid");
            for (int axis = 0; axis < 2; axis++)
            {
                SEXP names = VECTOR_ELT(dn, axis);
                if (TYPEOF(names) != STRSXP || XLENGTH(names) != nparents)
                    error("st.cross.allocate: inherited K names are invalid");
                for (int j = 0; j < nparents; j++)
                    if (strcmp(st_allocate_id(names, j), ids[j]))
                        error("st.cross.allocate: inherited K names must agree with the parental order");
            }
        }
        if ((size_t)nparents > SIZE_MAX / (size_t)nparents / sizeof(double))
            error("st.cross.allocate: inherited K exceeds the supported workspace size");
        double *copy = st_work_alloc((size_t)nparents * nparents, sizeof(double)), scale = 0;
        for (size_t i = 0; i < (size_t)nparents * nparents; i++)
        {
            double x = st_legacy_marker_number(matrix, i);
            if (!R_FINITE(x)) error("st.cross.allocate: inherited K must be finite");
            scale = fmax(scale, fabs(x));
        }
        for (int i = 0; i < nparents; i++) for (int j = 0; j < nparents; j++)
        {
            double x = scale ? st_legacy_marker_number(matrix, i + (size_t)nparents * j) / scale : 0,
                   y = scale ? st_legacy_marker_number(matrix, j + (size_t)nparents * i) / scale : 0;
            if (fabs(x - y) > 256 * DBL_EPSILON * nparents)
                error("st.cross.allocate: inherited K must be symmetric");
            copy[i + (size_t)nparents * j] = x * .5 + y * .5;
        }
        double *eigen = st_work_alloc(nparents, sizeof(double)), query;
        int lwork = -1, info = 0;
        F77_CALL(dsyev)("N", "U", &nparents, copy, &nparents, eigen, &query, &lwork, &info FCONE FCONE);
        if (info || !R_FINITE(query) || query > INT_MAX)
            error("st.cross.allocate: inherited K eigendecomposition failed");
        lwork = (int)query;
        double *work = st_work_alloc(lwork, sizeof(double));
        F77_CALL(dsyev)("N", "U", &nparents, copy, &nparents, eigen, work, &lwork, &info FCONE FCONE);
        if (info || eigen[0] < -256 * DBL_EPSILON * nparents * fmax(1, eigen[nparents - 1]))
            error("st.cross.allocate: inherited K must be positive semidefinite");
        SEXP lo = getListElement(constraints_in, "lower"), hi = getListElement(constraints_in, "upper"),
             fixed = getListElement(constraints_in, "requested.contributions");
        if (!st_allocate_numeric(lo) || !st_allocate_numeric(hi) || XLENGTH(lo) != nparents || XLENGTH(hi) != nparents)
            error("st.cross.allocate: inherited contribution bounds are invalid");
        for (int i = 0; i < nparents; i++)
        {
            double l = st_legacy_marker_number(lo, i), u = st_legacy_marker_number(hi, i);
            if (!R_FINITE(l) || !R_FINITE(u) || l < 0 || u > 1 || l > u)
                error("st.cross.allocate: inherited contribution bounds are invalid");
            if (cret[i] < l - 1e-12 || cret[i] > u + 1e-12)
                error("st.cross.allocate: nominal retained contributions violate inherited bounds");
            if (!isNull(fixed) && (!st_allocate_numeric(fixed) || XLENGTH(fixed) != nparents ||
                !R_FINITE(st_legacy_marker_number(fixed, i)) || fabs(cret[i] - st_legacy_marker_number(fixed, i)) > 1e-12))
                error("st.cross.allocate: nominal retained contributions violate inherited fixed contributions");
        }
        double theta = st_cross_scalar(getListElement(constraints_in, "theta"), "inherited theta"),
            inherited_tolerance = st_cross_scalar(getListElement(getListElement(plan, "control"), "tolerance"), "inherited tolerance");
        if (inherited_tolerance <= 0 || inherited_tolerance >= 1)
            error("st.cross.allocate: inherited tolerance must lie between zero and one");
        int64_t *uses = st_work_alloc(nparents, sizeof(int64_t));
        memset(uses, 0, (size_t)nparents * sizeof(int64_t));
        for (int e = 0; e < nf; e++) { uses[p1[e]] += retained[e]; uses[p2[e]] += retained[e]; }
        long double denominator = 2 * (long double)sum_retained, square = denominator * denominator,
            numerator = 0, absolute = 0, row_residual = 0, active_maximum = 0;
        for (int i = 0; i < nparents; i++) if (uses[i])
            for (int j = 0; j < nparents; j++) if (uses[j])
                active_maximum = fmaxl(active_maximum,
                    fabsl(st_legacy_marker_number(matrix, i + (size_t)nparents * j)));
        int integer_product = active_maximum <= LDBL_MAX / square;
        for (int i = 0; i < nparents; i++)
        {
            if (theta == 0)
            {
                long double row_scale = 0;
                for (int j = 0; j < nparents; j++) row_scale = fmaxl(row_scale,
                    fabsl(st_legacy_marker_number(matrix, i + (size_t)nparents * j)));
                if (row_scale)
                {
                    long double row_value = 0, norm = 0;
                    for (int j = 0; j < nparents; j++)
                    {
                        long double x = st_legacy_marker_number(matrix, i + (size_t)nparents * j) / row_scale;
                        row_value += x * uses[j]; norm += x * x;
                    }
                    row_residual = fmaxl(row_residual, fabsl(row_value) / (denominator * sqrtl(norm)));
                }
            }
            if (!uses[i]) continue;
            for (int j = 0; j < nparents; j++)
            {
                if (!uses[j]) continue;
                long double x = st_legacy_marker_number(matrix, i + (size_t)nparents * j);
                if (integer_product) x *= (long double)uses[i] * uses[j];
                else x *= ((long double)uses[i] / denominator) * ((long double)uses[j] / denominator);
                numerator += x; absolute += fabsl(x);
            }
        }
        if (integer_product) { numerator /= square; absolute /= square; }
        long double operations = 3 * (long double)nparents * nparents + 8,
            gamma = operations * LDBL_EPSILON / (1 - operations * LDBL_EPSILON),
            arithmetic_error = gamma * absolute + operations * nextafterl(0.0L, 1.0L),
            limit = (long double)theta + inherited_tolerance * fmaxl(fabsl(theta), absolute);
        qret = (double)numerator;
        if (!isfinite(numerator) || !isfinite(absolute) || !isfinite(arithmetic_error) ||
            !R_FINITE(qret) || operations * LDBL_EPSILON >= .5)
            error("st.cross.allocate: inherited retained relatedness cannot be certified");
        if ((theta == 0 && row_residual > 128 * DBL_EPSILON * nparents) ||
            (active_maximum == 0 ? theta < 0 : nextafterl(numerator + arithmetic_error, INFINITY) > limit))
            error("st.cross.allocate: nominal retained relatedness violates or cannot certify inherited theta");
    }
    SEXP diversity = inherited ? getListElement(inputs, "diversity") : R_NilValue;
    double diversity_loss = NA_REAL;
    if (inherited && isNull(diversity) &&
        (!isNull(getListElement(plan, "diversity")) || !isNull(getListElement(constraints_in, "diversity"))))
        error("st.cross.allocate: inherited diversity inputs are missing");
    if (!isNull(diversity))
    {
        SEXP labels = getAttrib(diversity, R_NamesSymbol);
        if (TYPEOF(diversity) != VECSXP || XLENGTH(diversity) != 3 ||
            TYPEOF(labels) != STRSXP || XLENGTH(labels) != 3)
            error("st.cross.allocate: inherited diversity needs G, max and sample components");
        unsigned int names_seen = 0;
        for (int i = 0; i < 3; i++)
        {
            if (STRING_ELT(labels, i) == NA_STRING) error("st.cross.allocate: invalid inherited diversity names");
            const char *name = CHAR(STRING_ELT(labels, i));
            int key = !strcmp(name, "G") ? 0 : !strcmp(name, "max") ? 1 : !strcmp(name, "sample") ? 2 : -1;
            if (key < 0 || (names_seen & (1u << key))) error("st.cross.allocate: invalid inherited diversity names");
            names_seen |= 1u << key;
        }
        if (strcmp(st_legacy_marker_string(getListElement(diversity, "sample"), "inherited diversity sample"), "family.slots"))
            error("st.cross.allocate: inherited diversity sample must be family.slots");
        SEXP maximum_arg = getListElement(diversity, "max");
        if (!st_allocate_numeric(maximum_arg)) error("st.cross.allocate: inherited diversity max must be numeric");
        double maximum = st_cross_scalar(maximum_arg, "inherited diversity max");
        int *slot_counts = st_work_alloc(nparents, sizeof(int));
        memset(slot_counts, 0, (size_t)nparents * sizeof(int));
        for (int e = 0; e < nf; e++) { slot_counts[p1[e]]++; slot_counts[p2[e]]++; }
        diversity_loss = st_cross_propsd_validate(getListElement(diversity, "G"), nparents, ids, slot_counts, maximum);
    }
    int extra = budget - (int)sum_lower;
    size_t utilities = 0;
    double estimated_transitions = 0;
    for (int e = 0; e < nf; e++)
    {
        int width = upper[e] - lower[e];
        if (width > extra) width = extra;
        if ((size_t)width + 1 > SIZE_MAX - utilities)
            error("st.cross.allocate: utility-table dimensions exceed the supported workspace size");
        utilities += (size_t)width + 1;
        estimated_transitions += ((double)extra + 1) * ((double)width + 1);
    }
    double cells = (double)nf * ((double)extra + 1) + utilities + 4 * ((double)extra + 1);
    int complete = cells <= max_cells && estimated_transitions <= max_transitions;
    const char *reason = complete ? "complete.fixed.pair.dp" : cells > max_cells ? "max.cells" : "max.transitions";
    int remain = extra;
    for (int e = 0; e < nf; e++)
    {
        int add = upper[e] - lower[e]; if (add > remain) add = remain;
        produced[e] = lower[e] + add; remain -= add;
    }
    double transitions = 0, response_upper = R_PosInf;
    if (complete)
    {
        size_t stride = (size_t)extra + 1;
        if ((size_t)nf > SIZE_MAX / stride / sizeof(int))
            error("st.cross.allocate: traceback dimensions exceed the supported workspace size");
        int *trace = st_work_alloc((size_t)nf * stride, sizeof(int)),
            *widths = st_work_alloc(nf, sizeof(int));
        size_t *offset = st_work_alloc(nf, sizeof(size_t));
        double *utility = st_work_alloc(utilities, sizeof(double)), *previous = st_work_alloc(stride, sizeof(double)),
            *current = st_work_alloc(stride, sizeof(double)), *up_previous = st_work_alloc(stride, sizeof(double)),
            *up_current = st_work_alloc(stride, sizeof(double));
        size_t at = 0;
        for (int e = 0; e < nf; e++)
        {
            widths[e] = upper[e] - lower[e]; if (widths[e] > extra) widths[e] = extra;
            offset[e] = at;
            double s = (retained[e] / (double)sum_retained) * accuracy[e] * sqrt(va[e]);
            if (!R_FINITE(s)) error("st.cross.allocate: a response scale is not finite");
            for (int x = 0; x <= widths[e]; x++)
            {
                utility[at++] = st_allocate_response(retained[e], lower[e] + x, s);
                if (!(x % 4096)) R_CheckUserInterrupt();
            }
        }
        for (int b = 0; b <= extra; b++) previous[b] = up_previous[b] = R_NegInf;
        previous[0] = up_previous[0] = 0;
        int reached = 0;
        for (int e = 0; e < nf; e++)
        {
            int next_reached = widths[e] > extra - reached ? extra : reached + widths[e];
            for (int b = 0; b <= next_reached; b++)
            {
                double best = R_NegInf, bound = R_NegInf;
                int bestx = -1, lo = b > reached ? b - reached : 0,
                    hi = b < widths[e] ? b : widths[e];
                for (int x = lo; x <= hi; x++)
                {
                    double u = utility[offset[e] + x], v = previous[b - x] + u,
                        z = nextafter(up_previous[b - x] + u, R_PosInf);
                    if (!R_FINITE(v) || !R_FINITE(z))
                        error("st.cross.allocate: accumulated response is not finite");
                    if (v > best) { best = v; bestx = x; }
                    if (z > bound) bound = z;
                    transitions++;
                }
                current[b] = best; up_current[b] = bound;
                trace[(size_t)e * stride + b] = bestx;
                if (!(b % 128)) R_CheckUserInterrupt();
            }
            double *swap = previous; previous = current; current = swap;
            swap = up_previous; up_previous = up_current; up_current = swap;
            reached = next_reached;
        }
        response_upper = up_previous[extra];
        int b = extra;
        for (int e = nf - 1; e >= 0; e--)
        {
            int x = trace[(size_t)e * stride + b];
            if (x < 0) error("st.cross.allocate: invalid allocation traceback");
            produced[e] = lower[e] + x; b -= x;
        }
        if (b) error("st.cross.allocate: incomplete allocation traceback");
    }
    long double mean_sum = 0, mean_upper = 0, response_sum = 0;
    double response_lower = 0;
    for (int e = 0; e < nf; e++)
    {
        double s = (retained[e] / (double)sum_retained) * accuracy[e] * sqrt(va[e]),
            response = st_allocate_response(retained[e], produced[e], s),
            c = produced[e] / (2.0 * budget);
        long double weight = (long double)retained[e] / sum_retained,
            weight_bound = nextafterl(weight, mu[e] < 0 ? -INFINITY : INFINITY),
            term_upper = nextafterl(weight_bound * mu[e], INFINITY);
        mean_sum += weight * mu[e];
        mean_upper = nextafterl(mean_upper + term_upper, INFINITY);
        response_sum += response;
        response_lower = nextafter(response_lower + response, R_NegInf);
        cprod[p1[e]] += c; cprod[p2[e]] += c;
    }
    double base = (double)mean_sum, response = (double)response_sum, value = (double)(mean_sum + response_sum);
    if (!R_FINITE(base) || !R_FINITE(response) || !R_FINITE(value))
        error("st.cross.allocate: the weighted objective is not finite");
    double base_upper = (double)mean_upper;
    if ((long double)base_upper < mean_upper) base_upper = nextafter(base_upper, R_PosInf);
    double gap = complete ? fmax(0, nextafter(response_upper - response_lower, R_PosInf)) : R_PosInf,
        bound = complete ? fmax(value, nextafter(base_upper + response_upper, R_PosInf)) : R_PosInf;
    const char *status = !complete ? "resource.limit" : gap <= tolerance * fmax(1, fabs(response)) ?
        "optimal" : "numerically.unresolved";
    if (inherited)
    {
        long double q = 0;
        for (int i = 0; i < nparents; i++) for (int j = 0; j < nparents; j++)
            q += (long double)cprod[i] * st_legacy_marker_number(matrix, i + (size_t)nparents * j) * cprod[j];
        qprod = (double)q;
        if (!R_FINITE(qprod)) error("st.cross.allocate: produced relatedness is not finite");
    }
    int np = 0;
    const char *out_names[] = {"crosses", "contributions.retained", "contributions.produced", "relatedness.retained",
        "relatedness.produced", "status", "objective", "constraints", "diagnostics", "inputs", "control"};
    SEXP out = st_allocate_list(11, out_names, &np);
    const char *columns[] = {"family", "P1Name", "P2Name", "n.progeny", "n.selected", "mu", "va", "accuracy",
        "intensity", "selected.mean", "min.progeny", "max.progeny"};
    SEXP crosses = st_legacy_marker_frame(12, columns, nf, &np);
    for (int j = 0; j < 12; j++) SET_VECTOR_ELT(crosses, j, st_allocate_keep(allocVector(
        j == 1 || j == 2 ? STRSXP : j == 0 || j == 3 || j == 4 || j > 9 ? INTSXP : REALSXP, nf), &np));
    for (int e = 0; e < nf; e++)
    {
        INTEGER(VECTOR_ELT(crosses, 0))[e] = e + 1;
        SET_STRING_ELT(VECTOR_ELT(crosses, 1), e, mkCharCE(ids[p1[e]], CE_UTF8));
        SET_STRING_ELT(VECTOR_ELT(crosses, 2), e, mkCharCE(ids[p2[e]], CE_UTF8));
        INTEGER(VECTOR_ELT(crosses, 3))[e] = produced[e];
        INTEGER(VECTOR_ELT(crosses, 4))[e] = retained[e];
        REAL(VECTOR_ELT(crosses, 5))[e] = mu[e]; REAL(VECTOR_ELT(crosses, 6))[e] = va[e];
        REAL(VECTOR_ELT(crosses, 7))[e] = accuracy[e];
        double intensity = st_cross_intensity((double)retained[e] / produced[e], retained[e], produced[e]),
            selected = mu[e] + accuracy[e] * sqrt(va[e]) * intensity;
        if (!R_FINITE(selected)) error("st.cross.allocate: a selected family mean is not finite");
        REAL(VECTOR_ELT(crosses, 8))[e] = intensity;
        REAL(VECTOR_ELT(crosses, 9))[e] = selected;
        INTEGER(VECTOR_ELT(crosses, 10))[e] = lower[e]; INTEGER(VECTOR_ELT(crosses, 11))[e] = upper[e];
    }
    SET_VECTOR_ELT(out, 0, crosses);
    SEXP parent_names = st_allocate_keep(allocVector(STRSXP, nparents), &np);
    for (int i = 0; i < nparents; i++) SET_STRING_ELT(parent_names, i, mkCharCE(ids[i], CE_UTF8));
    for (int j = 1; j <= 2; j++)
    {
        SEXP x = st_allocate_keep(allocVector(REALSXP, nparents), &np);
        memcpy(REAL(x), j == 1 ? cret : cprod, (size_t)nparents * sizeof(double));
        setAttrib(x, R_NamesSymbol, parent_names); SET_VECTOR_ELT(out, j, x);
    }
    SET_VECTOR_ELT(out, 3, st_allocate_keep(ScalarReal(qret), &np));
    SET_VECTOR_ELT(out, 4, st_allocate_keep(ScalarReal(qprod), &np));
    SET_VECTOR_ELT(out, 5, st_allocate_keep(mkString(status), &np));
    const char *objective_names[] = {"value", "response", "unselected.mean", "bound", "gap", "bound.type"};
    SEXP objective = st_allocate_list(6, objective_names, &np);
    SET_VECTOR_ELT(objective, 0, st_allocate_keep(ScalarReal(value), &np));
    SET_VECTOR_ELT(objective, 1, st_allocate_keep(ScalarReal(response), &np));
    SET_VECTOR_ELT(objective, 2, st_allocate_keep(ScalarReal(base), &np));
    SET_VECTOR_ELT(objective, 3, st_allocate_keep(ScalarReal(bound), &np));
    SET_VECTOR_ELT(objective, 4, st_allocate_keep(ScalarReal(gap), &np));
    SET_VECTOR_ELT(objective, 5, st_allocate_keep(mkString(complete ?
        "outward.addition.enclosure.of.evaluated.response.table" : "unbounded.resource.limit"), &np));
    SET_VECTOR_ELT(out, 6, objective);
    const char *constraint_names[] = {"total.progeny", "total.selected", "n.families", "contribution.sample",
        "inherited", "family.slots.fixed", "diversity.loss"};
    SEXP constraints = st_allocate_list(7, constraint_names, &np);
    SET_VECTOR_ELT(constraints, 0, st_allocate_keep(ScalarInteger(budget), &np));
    SET_VECTOR_ELT(constraints, 1, st_allocate_keep(ScalarReal((double)sum_retained), &np));
    SET_VECTOR_ELT(constraints, 2, st_allocate_keep(ScalarInteger(nf), &np));
    SET_VECTOR_ELT(constraints, 3, st_allocate_keep(mkString("nominal.retained.family.shares"), &np));
    SET_VECTOR_ELT(constraints, 4, constraints_in);
    SET_VECTOR_ELT(constraints, 5, st_allocate_keep(ScalarLogical(inherited), &np));
    SET_VECTOR_ELT(constraints, 6, st_allocate_keep(ScalarReal(diversity_loss), &np));
    SET_VECTOR_ELT(out, 7, constraints);
    const char *diagnostic_names[] = {"reason", "complete.search", "input.rows", "pooled.rows", "workspace.cells",
        "transition.upper.estimate", "transitions", "numerical.qualification", "ancestry.qualification"};
    SEXP diagnostics = st_allocate_list(9, diagnostic_names, &np);
    SET_VECTOR_ELT(diagnostics, 0, st_allocate_keep(mkString(reason), &np));
    SET_VECTOR_ELT(diagnostics, 1, st_allocate_keep(ScalarLogical(complete), &np));
    SET_VECTOR_ELT(diagnostics, 2, st_allocate_keep(ScalarInteger(nr), &np));
    SET_VECTOR_ELT(diagnostics, 3, st_allocate_keep(ScalarInteger(nr - nf), &np));
    SET_VECTOR_ELT(diagnostics, 4, st_allocate_keep(ScalarReal(cells), &np));
    SET_VECTOR_ELT(diagnostics, 5, st_allocate_keep(ScalarReal(estimated_transitions), &np));
    SET_VECTOR_ELT(diagnostics, 6, st_allocate_keep(ScalarReal(transitions), &np));
    SET_VECTOR_ELT(diagnostics, 7, st_allocate_keep(mkString(
        "Gap encloses addition rounding for the evaluated binary64 response table; it does not bound the normal or finite-family approximation error."), &np));
    SET_VECTOR_ELT(diagnostics, 8, st_allocate_keep(mkString(
        "Produced and retained contributions are nominal parental shares, not realized genomic ancestry after within-family selection."), &np));
    SET_VECTOR_ELT(out, 8, diagnostics);
    const char *input_names[] = {"plan", "scores", "input.family"};
    SEXP input_out = st_allocate_list(3, input_names, &np);
    SET_VECTOR_ELT(input_out, 0, plan); SET_VECTOR_ELT(input_out, 1, scores);
    SEXP mapping = st_allocate_keep(allocVector(INTSXP, nr), &np);
    for (int i = 0; i < nr; i++) INTEGER(mapping)[i] = input_family[i] + 1;
    SET_VECTOR_ELT(input_out, 2, mapping); SET_VECTOR_ELT(out, 9, input_out);
    const char *control_names[] = {"max.cells", "max.transitions", "tolerance"};
    SEXP control_out = st_allocate_list(3, control_names, &np);
    SET_VECTOR_ELT(control_out, 0, st_allocate_keep(ScalarReal(max_cells), &np));
    SET_VECTOR_ELT(control_out, 1, st_allocate_keep(ScalarReal(max_transitions), &np));
    SET_VECTOR_ELT(control_out, 2, st_allocate_keep(ScalarReal(tolerance), &np));
    SET_VECTOR_ELT(out, 10, control_out);
    setAttrib(out, R_ClassSymbol, st_allocate_keep(mkString("st.cross.allocation"), &np));
    UNPROTECT(np);
    return out;
}

static int st_allocate_print_setting(SEXP args, SEXP options, const char *name,
                                      const char *fallback, int min, int max, int infinity)
{
    SEXP value = getListElement(args, name);
    if (isNull(value)) value = getListElement(options, fallback);
    if (!st_allocate_numeric(value) || XLENGTH(value) != 1)
        error("print: %s must be one number", name);
    double x = asReal(value);
    if (infinity && x == R_PosInf) return INT_MAX;
    if (!R_FINITE(x) || x != floor(x) || x < min || x > max)
        error("print: %s is outside its supported whole-number range", name);
    return (int)x;
}

static void st_allocate_print_line(int width, const char *format, ...)
{
    va_list ap, copy; va_start(ap, format); va_copy(copy, ap);
    int n = vsnprintf(NULL, 0, format, copy); va_end(copy);
    if (n < 0) error("print: unable to format allocation");
    char *text = st_work_alloc((size_t)n + 1, sizeof(char));
    vsnprintf(text, (size_t)n + 1, format, ap); va_end(ap);
    for (int i = 0; i < n; i++) if ((unsigned char)text[i] < 32 || text[i] == 127) text[i] = ' ';
    char *start = text;
    while ((int)strlen(start) > width)
    {
        int split = width;
        while (split > 0 && start[split] != ' ') split--;
        if (!split) { split = width; while (split > 0 && ((unsigned char)start[split] & 0xc0) == 0x80) split--; }
        if (!split) split = width;
        Rprintf("%.*s\n", split, start); start += split;
        while (*start == ' ') start++;
    }
    Rprintf("%s\n", start);
}

static SEXP st_legacy_cross_allocate_print(SEXP args)
{
    SEXP x = getListElement(args, "object"), options = getListElement(args, "options");
    int digits = st_allocate_print_setting(args, options, "digits", "print.digits", 1, 16, 0),
        top = st_allocate_print_setting(args, options, "top", "print.top", 0, INT_MAX, 1),
        bottom = st_allocate_print_setting(args, options, "bottom", "print.bottom", 0, INT_MAX, 1),
        width = st_allocate_print_setting(args, options, "char.per.line", "print.char.perline", 40, 10000, 0);
    SEXP crosses = getListElement(x, "crosses"), objective = getListElement(x, "objective"),
         constraints = getListElement(x, "constraints");
    if (!inherits(x, "st.cross.allocation") || !inherits(crosses, "data.frame"))
        error("print: invalid st.cross.allocation object");
    SEXP p1 = getListElement(crosses, "P1Name"), p2 = getListElement(crosses, "P2Name"),
         produced = getListElement(crosses, "n.progeny"), selected = getListElement(crosses, "n.selected"),
         means = getListElement(crosses, "selected.mean");
    if (TYPEOF(p1) != STRSXP || TYPEOF(p2) != STRSXP || XLENGTH(p1) > INT_MAX ||
        TYPEOF(produced) != INTSXP || TYPEOF(selected) != INTSXP || TYPEOF(means) != REALSXP ||
        XLENGTH(p2) != XLENGTH(p1) || XLENGTH(produced) != XLENGTH(p1) ||
        XLENGTH(selected) != XLENGTH(p1) || XLENGTH(means) != XLENGTH(p1))
        error("print: invalid st.cross.allocation family columns");
    int n = (int)XLENGTH(p1);
    for (int i = 0; i < n; i++)
    {
        st_allocate_id(p1, i); st_allocate_id(p2, i);
        if (INTEGER(produced)[i] < 1 || INTEGER(selected)[i] < 1 ||
            INTEGER(selected)[i] > INTEGER(produced)[i] || !R_FINITE(REAL(means)[i]))
            error("print: invalid st.cross.allocation family values");
    }
    st_allocate_print_line(width, "Fixed-pair progeny allocation");
    st_allocate_print_line(width, "Status: %s", st_legacy_marker_string(getListElement(x, "status"), "status"));
    st_allocate_print_line(width, "Families: %d; produced: %.0f; retained: %.0f", n,
        asReal(getListElement(constraints, "total.progeny")), asReal(getListElement(constraints, "total.selected")));
    st_allocate_print_line(width, "Predicted retained mean: %.*g; response: %.*g", digits,
        asReal(getListElement(objective, "value")), digits, asReal(getListElement(objective, "response")));
    st_allocate_print_line(width, "Objective upper bound: %.*g; response gap: %.*g", digits,
        asReal(getListElement(objective, "bound")), digits, asReal(getListElement(objective, "gap")));
    st_allocate_print_line(width, "Group relatedness, nominal retained: %.*g; produced: %.*g", digits,
        asReal(getListElement(x, "relatedness.retained")), digits, asReal(getListElement(x, "relatedness.produced")));
    st_allocate_print_line(width, "Family: parents; produced / retained; predicted retained mean");
    int ntop = top < n ? top : n, nbot = bottom < n - ntop ? bottom : n - ntop;
    for (int i = 0; i < n; i++)
    {
        if (i == ntop && n - ntop - nbot > 0)
            st_allocate_print_line(width, "... %d families omitted ...", n - ntop - nbot);
        if (i >= ntop && i < n - nbot) continue;
        st_allocate_print_line(width, "%d: %s x %s; %d / %d; %.*g", i + 1,
            st_allocate_id(p1, i), st_allocate_id(p2, i), INTEGER(produced)[i], INTEGER(selected)[i],
            digits, REAL(means)[i]);
    }
    return x;
}

static SEXP st_legacy_cross_es(SEXP args)
{
    double alpha = st_cross_scalar(getListElement(args, "alpha"), "alpha"),
           n = st_cross_scalar(getListElement(args, "N"), "N"),
           g = st_cross_scalar(getListElement(args, "G"), "G");
    if (alpha <= 0 || alpha > 1) error("alpha must be greater than zero and at most one");
    if (n < 0 || g < 0 || n != floor(n) || g != floor(g) ||
        n > 9007199254740991.0 || g > 9007199254740991.0 ||
        ((n == 0) != (g == 0)) || n > g)
        error("N and G must both be zero or whole counts satisfying 1 <= N <= G");
    if (n > 0) alpha = n / g;
    double intensity = st_cross_intensity(alpha, n, g);
    gs_varset_type *dataset = st_cross_dataset(getListElement(args, "dataset"));
    int status = 0;
    gs_cross_eval_es_01(dataset, &intensity, &status);
    if (status == -2) error("Cross means and nonnegative variances must be evaluated before selection");
    return R_NilValue;
}

static SEXP st_legacy_cross_va(SEXP args)
{
    const char *population = st_legacy_marker_string(getListElement(args, "pop.type"), "pop.type"),
               *map_function = st_legacy_marker_string(getListElement(args, "map.function"), "map.function");
    double generation = st_cross_scalar(getListElement(args, "t"), "t");
    if (generation < 0 || generation > INT_MAX || generation != floor(generation))
        error("t must be a nonnegative whole number at most the maximum integer");
    if (strcmp(population, "unlinked") && strcmp(population, "DH") &&
        strcmp(population, "SSD") && strcmp(population, "CRS"))
        error("pop.type must be unlinked, DH, SSD or CRS");
    if (strcmp(map_function, "Haldane") && strcmp(map_function, "Kosambi"))
        error("map.function must be Haldane or Kosambi");
    gs_varset_type *dataset = st_cross_dataset(getListElement(args, "dataset"));
    if (dataset->NoInd < 2 || !dataset->NoMar || !st_allele_present(dataset))
        error("At least two parents with marker data are required");
    if (!dataset->u || !dataset->EffAll || dataset->dimu != dataset->NoMar + 1 ||
        dataset->NoEff != dataset->NoMar + 1)
        error("Marker effects are required for every marker");
    if (strcmp(population, "unlinked") && (!dataset->NLoc || !dataset->map))
        error("A linkage map is required for linked cross variances");
    for (gs_UNSIGNED i = 0; i < dataset->dimu; i++)
        if (!R_FINITE(dataset->u[i])) error("Marker effects must be finite");
    int t = (int)generation, status = 0;
    char *pop = (char *)population, *map = (char *)map_function;
    gs_cross_eval_va_01(dataset, &pop, &t, &map, &status);
    if (status == -2) error("Cross variances could not be evaluated");
    return R_NilValue;
}

static SEXP st_legacy_cross_info(SEXP args)
{
    double bestn = st_cross_scalar(getListElement(args, "bestn"), "bestn");
    if (bestn < 0 || bestn != floor(bestn) || bestn > INT_MAX)
        error("bestn must be a whole number between zero and the maximum integer");
    const char *sortby = st_legacy_marker_string(getListElement(args, "sortby"), "sortby");
    int (*compare)(const void *, const void *) = NULL;
    unsigned int criterion = 0;
    if (!strcmp(sortby, "index")) compare = gs_cross_cmp_index;
    else if (!strcmp(sortby, "gd")) { compare = gs_cross_cmp_gd; criterion = ST_CROSS_GD; }
    else if (!strcmp(sortby, "mu")) { compare = gs_cross_cmp_mu; criterion = ST_CROSS_MU; }
    else if (!strcmp(sortby, "mi")) { compare = gs_cross_cmp_mi; criterion = ST_CROSS_MI; }
    else if (!strcmp(sortby, "ma")) { compare = gs_cross_cmp_ma; criterion = ST_CROSS_MA; }
    else if (!strcmp(sortby, "va")) { compare = gs_cross_cmp_va; criterion = ST_CROSS_VA; }
    else if (!strcmp(sortby, "es")) { compare = gs_cross_cmp_es; criterion = ST_CROSS_ES; }
    else error("sortby must be index, gd, mu, mi, ma, va or es");
    gs_varset_type *dataset = st_cross_dataset(getListElement(args, "dataset"));
    st_cross_current(dataset);
    if (!dataset->crs || !dataset->crs_S) error("No evaluated cross information is available");
    if (criterion && !(dataset->cross_ready & criterion))
        error("The requested sorting criterion has not been evaluated");
    if (dataset->NoInd > INT_MAX || dataset->NCrs > INT_MAX)
        error("Cross information exceeds the supported data-frame dimensions");
    gs_cross_type **order = st_work_alloc(dataset->NCrs, sizeof(*order));
    memcpy(order, dataset->crs, dataset->NCrs * sizeof(*order));
    qsort(order, dataset->NCrs, sizeof(*order), compare);
    int n = bestn > 0 && bestn < dataset->NCrs ? (int)bestn : (int)dataset->NCrs;
    const char *names[] = {"P1No", "P2No", "P1Name", "P2Name", "gd", "mu", "mi", "ma", "va", "es"};
    unsigned int flags[] = {ST_CROSS_GD, ST_CROSS_MU, ST_CROSS_MI, ST_CROSS_MA, ST_CROSS_VA, ST_CROSS_ES};
    int np = 0;
    SEXP out = PROTECT(allocVector(VECSXP, 10)); np++;
    SEXP labels = PROTECT(allocVector(STRSXP, 10)); np++;
    for (int j = 0; j < 10; j++)
    {
        SET_STRING_ELT(labels, j, mkChar(names[j]));
        SET_VECTOR_ELT(out, j, allocVector(j < 2 ? INTSXP : (j < 4 ? STRSXP : REALSXP), n));
    }
    for (int i = 0; i < n; i++)
    {
        const gs_cross_type *cross = order[i];
        INTEGER(VECTOR_ELT(out, 0))[i] = (int)cross->p1 + 1;
        INTEGER(VECTOR_ELT(out, 1))[i] = (int)cross->p2 + 1;
        SET_STRING_ELT(VECTOR_ELT(out, 2), i, mkChar(dataset->IndName[cross->p1]));
        SET_STRING_ELT(VECTOR_ELT(out, 3), i, mkChar(dataset->IndName[cross->p2]));
        double values[] = {cross->gd, cross->mu, cross->mi, cross->ma, cross->va, cross->es};
        for (int j = 0; j < 6; j++)
            REAL(VECTOR_ELT(out, j + 4))[i] = dataset->cross_ready & flags[j] ? values[j] : NA_REAL;
    }
    SEXP rows = PROTECT(allocVector(INTSXP, 2)); np++;
    INTEGER(rows)[0] = NA_INTEGER; INTEGER(rows)[1] = -n;
    setAttrib(out, R_NamesSymbol, labels);
    setAttrib(out, R_RowNamesSymbol, rows);
    SEXP classes = PROTECT(mkString("data.frame")); np++;
    setAttrib(out, R_ClassSymbol, classes);
    UNPROTECT(np);
    return out;
}

static SEXP st_vcf_keep(SEXP value, int *np)
{
    PROTECT(value);
    ++*np;
    return value;
}

static SEXP st_vcf_coerce(const char *name, SEXP value, int quiet, int *np)
{
    SEXP quoted = st_vcf_keep(lang2(install("quote"), value), np);
    SEXP call = st_vcf_keep(lang2(findFun(install(name), R_BaseEnv), quoted), np);
    if (quiet)
        call = st_vcf_keep(lang2(findFun(install("suppressWarnings"), R_BaseEnv), call), np);
    return st_vcf_keep(eval(call, R_GlobalEnv), np);
}

static const char *st_vcf_text(SEXP value)
{
    return getCharCE(value) == CE_BYTES ? CHAR(value) : translateCharUTF8(value);
}

static int st_vcf_bytes_equal(SEXP a, SEXP b)
{
    if (a == NA_STRING || b == NA_STRING) return a == b;
    return (getCharCE(a) == CE_BYTES) == (getCharCE(b) == CE_BYTES) && !strcmp(CHAR(a), CHAR(b));
}

static int st_vcf_identifier(SEXP value, double maximum)
{
    if (value == NA_STRING || !*CHAR(value) || strlen(CHAR(value)) > maximum)
        return 0;
    const unsigned char *text = (const unsigned char *)st_vcf_text(value);
    for (; *text; text++)
        if (*text <= 32 || *text == 127) return 0;
    return 1;
}

static SEXP st_vcf_identifiers(SEXP args)
{
    SEXP value = getListElement(args, "data"), maximum = getListElement(args, "max.bytes");
    R_xlen_t n = isNull(value) ? 0 : XLENGTH(value);
    int np = 0;
    SEXP text = value;
    if (n && TYPEOF(text) != STRSXP) text = st_vcf_coerce("as.character", value, 0, &np);
    SEXP out = st_vcf_keep(allocVector(LGLSXP, n), &np);
    for (R_xlen_t i = 0; i < n; i++)
    {
        SEXP x = STRING_ELT(text, i);
        if (x == NA_STRING || !*CHAR(x)) { LOGICAL(out)[i] = FALSE; continue; }
        int length = isNull(maximum) ? 0 : LENGTH(maximum);
        if (length > 1) error("'length = %d' in coercion to 'logical(1)'", length);
        double limit = length ? asReal(maximum) : NA_REAL;
        if (ISNAN(limit)) error("missing value where TRUE/FALSE needed");
        if (strlen(CHAR(x)) > limit) { LOGICAL(out)[i] = FALSE; continue; }
        if (TYPEOF(value) != STRSXP) error("argument is not a character vector");
        LOGICAL(out)[i] = st_vcf_identifier(x, limit);
    }
    SEXP names = getAttrib(value, R_NamesSymbol);
    if (isNull(names) && TYPEOF(value) == STRSXP) names = value;
    if (!isNull(names)) setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(np);
    return out;
}

static int st_vcf_same(SEXP a, SEXP b)
{
    if (a == NA_STRING || b == NA_STRING) return a == b;
    if (getCharCE(a) == CE_BYTES || getCharCE(b) == CE_BYTES) return st_vcf_bytes_equal(a, b);
    return !strcmp(translateCharUTF8(a), translateCharUTF8(b));
}

static int st_vcf_unique(SEXP value)
{
    R_xlen_t n = XLENGTH(value);
    int bytes = 0;
    for (R_xlen_t i = 0; i < n; i++) bytes |= getCharCE(STRING_ELT(value, i)) == CE_BYTES;
    size_t capacity = 8;
    while (capacity < (size_t)n * 2)
    {
        if (capacity > SIZE_MAX / 2) error("Too many VCF identifiers");
        capacity *= 2;
    }
    R_xlen_t *index = st_work_alloc(capacity, sizeof(*index));
    memset(index, 0, capacity * sizeof(*index));
    for (R_xlen_t i = 0; i < n; i++)
    {
        SEXP element = STRING_ELT(value, i);
        const char *text = element == NA_STRING ? "NA" : bytes ? CHAR(element) : translateCharUTF8(element);
        uint64_t hash = UINT64_C(14695981039346656037);
        for (const unsigned char *p = (const unsigned char *)text; *p; p++)
            hash = (hash ^ *p) * UINT64_C(1099511628211);
        size_t at = (size_t)hash & (capacity - 1);
        while (index[at] && !(bytes ? st_vcf_bytes_equal(element, STRING_ELT(value, index[at] - 1))
                                  : st_vcf_same(element, STRING_ELT(value, index[at] - 1))))
            at = (at + 1) & (capacity - 1);
        if (index[at]) return 0;
        index[at] = i + 1;
    }
    return 1;
}

static int st_vcf_fields(const char *text, char separator, int *gt)
{
    int n = 0, matches = 0;
    const char *p = text;
    if (!*p || !strcmp(p, ".")) return 0;
    for (;;)
    {
        const char *end = strchr(p, separator);
        size_t size = end ? (size_t)(end - p) : strlen(p);
        if (!size) return 0;
        n++;
        if (size == 2 && p[0] == 'G' && p[1] == 'T')
            *gt = n, matches++;
        if (!end) break;
        p = end + 1;
    }
    return matches == 1 ? n : 0;
}

static int st_vcf_alt(SEXP alt, SEXP ref)
{
    int bytes = getCharCE(alt) == CE_BYTES || getCharCE(ref) == CE_BYTES;
    int comparable = (getCharCE(alt) == CE_BYTES) == (getCharCE(ref) == CE_BYTES);
    const char *text = bytes ? CHAR(alt) : translateCharUTF8(alt),
               *reference = bytes ? CHAR(ref) : translateCharUTF8(ref);
    if (!strcmp(text, ".")) return 0;
    const char *start = text;
    int count = 0;
    for (;;)
    {
        const char *end = strchr(start, ',');
        size_t n = end ? (size_t)(end - start) : strlen(start);
        if (!n || (n == 1 && *start == '.') ||
            (comparable && strlen(reference) == n && !memcmp(reference, start, n))) return -1;
        const char *earlier = text;
        while (earlier < start)
        {
            const char *last = strchr(earlier, ',');
            if (last && (size_t)(last - earlier) == n && !memcmp(earlier, start, n))
                return -1;
            earlier = last ? last + 1 : start;
        }
        count++;
        if (!end) break;
        start = end + 1;
    }
    return count;
}

static SEXP st_vcf_numeric(SEXP value, int *np)
{
    if (inherits(value, "factor")) value = st_vcf_coerce("as.character", value, 0, np);
    return st_vcf_coerce("as.numeric", value, 1, np);
}

static SEXP st_vcf_cell(SEXP column, int row, int index, int *np)
{
    if (TYPEOF(column) == STRSXP && !isObject(column)) return STRING_ELT(column, row);
    SEXP quoted = st_vcf_keep(lang2(install("quote"), column), np),
         call = st_vcf_keep(lang3(findFun(install("["), R_BaseEnv), quoted,
                                  st_vcf_keep(ScalarInteger(row + 1), np)), np),
         sliced = st_vcf_keep(eval(call, R_GlobalEnv), np);
    quoted = st_vcf_keep(lang2(install("quote"), sliced), np);
    call = st_vcf_keep(lang3(findFun(install("[["), R_BaseEnv), quoted,
                             st_vcf_keep(ScalarInteger(1), np)), np);
    SEXP scalar = st_vcf_keep(eval(call, R_GlobalEnv), np),
         text = st_vcf_coerce("as.character", scalar, 0, np);
    if (TYPEOF(text) != STRSXP || LENGTH(text) != 1)
        error("values must be length 1,\n but FUN(X[[%d]]) result is length %d", index + 1, LENGTH(text));
    return STRING_ELT(text, 0);
}

static int st_vcf_allele(const char *text, size_t n, int maximum)
{
    if ((n == 1 && *text == '.') || (n == 2 && !memcmp(text, "-1", 2))) return -1;
    int value = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (text[i] < '0' || text[i] > '9') return -2;
        if (value <= maximum) value = value * 10 + text[i] - '0';
    }
    return value > maximum ? -3 : value;
}

static void st_vcf_metadata(SEXP from, SEXP to)
{
    const char *keys[] = {"pedigree", "prediction.state", "allele.codes"};
    for (int i = 0; i < 3; i++)
        setAttrib(to, install(keys[i]), getAttrib(from, install(keys[i])));
}

static SEXP st_vcf_from_frame(SEXP args)
{
    SEXP data = getListElement(args, "data");
    int np = 0;
    if (!inherits(data, "data.frame")) error("data must be a data frame");
    SEXP names = getAttrib(data, R_NamesSymbol);
    for (int j = 0; j < LENGTH(names); j++)
        if (STRING_ELT(names, j) == NA_STRING || !*CHAR(STRING_ELT(names, j)))
            error("data frame column names must be non-missing and non-empty");
    int nc = LENGTH(data);
    if (nc < 6) error("data must contain fixed fields and at least one sample column");
    if (TYPEOF(names) != STRSXP || LENGTH(names) < 5) error("first five columns must be CHROM/#CHROM, POS, ID, REF, ALT");
    const char *fixed[] = {"POS", "ID", "REF", "ALT"};
    int valid = !strcmp(CHAR(STRING_ELT(names, 0)), "CHROM") ||
                !strcmp(CHAR(STRING_ELT(names, 0)), "#CHROM");
    for (int j = 1; j < 5; j++) valid &= !strcmp(CHAR(STRING_ELT(names, j)), fixed[j - 1]);
    if (!valid) error("first five columns must be CHROM/#CHROM, POS, ID, REF, ALT");
    int nr = LENGTH(VECTOR_ELT(data, 0));
    if (nr <= 0) error("data contains no markers");
    SEXP chromosomes = st_vcf_coerce("as.character", VECTOR_ELT(data, 0), 0, &np);
    int *chrom = st_work_alloc(nr, sizeof(int)), nchrom = 0, chromosome_bytes = 0;
    for (int i = 0; i < nr; i++) chromosome_bytes |= getCharCE(STRING_ELT(chromosomes, i)) == CE_BYTES;
    SEXP *levels = st_work_alloc(nr, sizeof(*levels));
    size_t message_size = 32;
    for (int i = 0; i < nr; i++)
    {
        SEXP label = STRING_ELT(chromosomes, i);
        if (label == NA_STRING || !*CHAR(label)) error("CHROM contains missing or empty values");
        int at = 0;
        while (at < nchrom && !(chromosome_bytes ? st_vcf_bytes_equal(label, levels[at]) : st_vcf_same(label, levels[at]))) at++;
        if (at == nchrom)
        {
            levels[nchrom++] = label;

        }
        chrom[i] = at + 1;
    }
    for (int i = 0; i < nchrom; i++) message_size += strlen(translateCharUTF8(levels[i])) + 40;
    char *message = st_work_alloc(message_size, 1);
    strcpy(message, "Chromosome coding: ");
    size_t used = strlen(message);
    for (int i = 0; i < nchrom; i++)
        used += (size_t)snprintf(message + used, message_size - used, "%s%d=%s", i ? ", " : "", i + 1, translateCharUTF8(levels[i]));
    st_info(0, message);
    int format = -1, cm = -1, formats = 0, cms = 0;
    SEXP format_values = R_NilValue;
    for (int j = 5; j < nc; j++)
        if (!strcmp(CHAR(STRING_ELT(names, j)), "FORMAT"))
        {
            SEXP value = st_vcf_coerce("as.character", VECTOR_ELT(data, j), 0, &np);
            int good = 1;
            for (int i = 0; i < nr; i++)
            {
                int gt = 0;
                SEXP x = STRING_ELT(value, i);
                if (x == NA_STRING || !st_vcf_fields(CHAR(x), ':', &gt)) good = 0;
            }
            if (good)
            {
                formats++; format = j; format_values = value;
            }
        }
    if (formats > 1) error("more than one column can be interpreted as FORMAT");
    if (format >= 0)
    {
        const char *metadata[] = {"QUAL", "FILTER", "INFO", "CM"};
        int seen[4] = {0};
        for (int j = 5; j < format; j++)
        {
            int at = 0;
            while (at < 4 && strcmp(CHAR(STRING_ELT(names, j)), metadata[at])) at++;
            if (at == 4) error("unrecognized column before FORMAT");
            seen[at]++;
            if (at == 3) cm = j;
        }
        for (int j = 0; j < 4; j++)
            if (seen[j] > 1) error("VCF metadata columns before FORMAT must be unique");
        if (format == nc - 1) error("FORMAT is present but no sample columns follow it");
    }
    else
        for (int j = 5; j < nc; j++)
            if (!strcmp(CHAR(STRING_ELT(names, j)), "CM"))
            {
                SEXP value = st_vcf_numeric(VECTOR_ELT(data, j), &np);
                int good = 1;
                for (int i = 0; i < LENGTH(value); i++)
                    if (!R_FINITE(REAL(value)[i]) || REAL(value)[i] < 0) good = 0;
                if (good)
                {
                    cms++; cm = j;
                }
            }
    if (cms > 1) error("more than one column can be interpreted as CM");
    int ns = 0, *samples = st_work_alloc(nc, sizeof(int));
    for (int j = format < 0 ? 5 : format + 1; j < nc; j++)
        if (j != cm) samples[ns++] = j;
    if (ns <= 0) error("data contains no sample genotype columns");
    SEXP individuals = st_vcf_keep(allocVector(STRSXP, ns), &np);
    for (int j = 0; j < ns; j++) SET_STRING_ELT(individuals, j, STRING_ELT(names, samples[j]));
    if (!st_vcf_unique(individuals)) error("sample names must be unique");
    for (int j = 0; j < ns; j++)
        if (!st_vcf_identifier(STRING_ELT(individuals, j), R_PosInf))
            error("sample names are invalid for SelectionTools");
    SEXP positions = st_vcf_numeric(VECTOR_ELT(data, cm < 0 ? 1 : cm), &np);
    const char *position_name = cm < 0 ? "POS" : "CM";
    for (int i = 0; i < LENGTH(positions); i++)
        if (!R_FINITE(REAL(positions)[i])) error("%s contains invalid positions", position_name);
    for (int i = 0; i < LENGTH(positions); i++)
        if (REAL(positions)[i] < 0) error("%s contains negative genetic positions", position_name);
    SEXP ids = st_vcf_coerce("as.character", VECTOR_ELT(data, 2), 0, &np);
    valid = st_vcf_unique(ids);
    for (int i = 0; i < LENGTH(ids); i++)
        if (STRING_ELT(ids, i) == NA_STRING || !*CHAR(STRING_ELT(ids, i))) valid = 0;
    if (!valid) error("marker IDs are required and must be unique");
    for (int i = 0; i < LENGTH(ids); i++)
        if (!st_vcf_identifier(STRING_ELT(ids, i), 255)) error("marker IDs are invalid for SelectionTools");
    SEXP refs = st_vcf_coerce("as.character", VECTOR_ELT(data, 3), 0, &np);
    for (int i = 0; i < LENGTH(refs); i++)
        if (STRING_ELT(refs, i) == NA_STRING || !*CHAR(STRING_ELT(refs, i)) ||
            !strcmp(CHAR(STRING_ELT(refs, i)), ".") || strchr(CHAR(STRING_ELT(refs, i)), ','))
            error("REF must contain one non-missing allele label per marker");
    SEXP alts = st_vcf_coerce("as.character", VECTOR_ELT(data, 4), 0, &np);
    for (int i = 0; i < LENGTH(alts); i++)
        if (STRING_ELT(alts, i) == NA_STRING || !*CHAR(STRING_ELT(alts, i)))
            error("ALT may be '.' but may not be missing or empty");
    int *nalt = st_work_alloc(nr, sizeof(int)), too_many = 0;
    for (int i = 0; i < nr; i++)
    {
        nalt[i] = st_vcf_alt(STRING_ELT(alts, i), STRING_ELT(refs, i));
        if (nalt[i] < 0) error("invalid or duplicated ALT field at marker %s", CHAR(STRING_ELT(ids, i)));
        too_many |= nalt[i] > 253;
    }
    if (too_many) error("STvcf supports at most 253 ALT alleles per marker");
    if ((size_t)nr > (size_t)PTRDIFF_MAX / (size_t)ns)
        error("VCF genotype storage exceeds the maximum object size");
    size_t cells = (size_t)nr * (size_t)ns;
    unsigned char *first = st_work_alloc(cells, sizeof(*first));
    unsigned char *second = st_work_alloc(cells, sizeof(*second));
    memset(first, GS_ALLELE_MISSING, cells);
    memset(second, GS_ALLELE_MISSING, cells);
    int phased = 0, unphased = 0;
    for (int i = 0; i < nr; i++)
    {
        int gt = 1;
        if (format >= 0) st_vcf_fields(CHAR(STRING_ELT(format_values, i)), ':', &gt);
        int before = np;
        void *row_mark = st_work_mark();
        SEXP *values = st_work_alloc(ns, sizeof(*values));
        const char **texts = st_work_alloc(ns, sizeof(*texts));
        for (int j = 0; j < ns; j++) values[j] = st_vcf_cell(VECTOR_ELT(data, samples[j]), i, j, &np);
        unsigned char seen[254] = {0};
        int missing = 0;
        for (int j = 0; j < ns; j++)
        {
            const char *x = values[j] == NA_STRING ? NULL : CHAR(values[j]);
            if (format >= 0 && x && *x && strcmp(x, "."))
            {
                for (int k = 1; k < gt; k++)
                {
                    x = strchr(x, ':');
                    if (!x) break;
                    x++;
                }
                if (!x || !*x || *x == ':')
                    error("sample FORMAT field does not contain GT at marker %s sample %s", CHAR(STRING_ELT(ids, i)), CHAR(STRING_ELT(individuals, j)));
                const char *end = strchr(x, ':');
                if (end)
                {
                    size_t n = (size_t)(end - x); char *copy = st_work_alloc(n + 1, 1);
                    memcpy(copy, x, n); copy[n] = '\0'; x = copy;
                }
            }
            texts[j] = x;
        }
        for (int j = 0; j < ns; j++)
        {
            const char *x = texts[j];
            if (x) phased |= strchr(x, '|') != NULL, unphased |= strchr(x, '/') != NULL;
            if (!x || !*x || !strcmp(x, ".") || !strcmp(x, "-1") ||
                !strcmp(x, "./.") || !strcmp(x, ".|.") || !strcmp(x, "-1/-1") ||
                !strcmp(x, "-1|-1") || !strcmp(x, "NA"))
            { missing = 1; continue; }
            const char *separator = strpbrk(x, "/|");
            if (!separator || strpbrk(separator + 1, "/|"))
                error("non-diploid or invalid GT at marker %s sample %s", CHAR(STRING_ELT(ids, i)), CHAR(STRING_ELT(individuals, j)));
            size_t n1 = (size_t)(separator - x), n2 = strlen(separator + 1);
            if (!n1 || !n2) error("invalid GT at marker %s sample %s", CHAR(STRING_ELT(ids, i)), CHAR(STRING_ELT(individuals, j)));
            int a = st_vcf_allele(x, n1, nalt[i]), b = st_vcf_allele(separator + 1, n2, nalt[i]);
            if (a == -2 || b == -2) error("invalid allele index at marker %s sample %s", CHAR(STRING_ELT(ids, i)), CHAR(STRING_ELT(individuals, j)));
            if (a == -3 || b == -3) error("allele index exceeds REF/ALT definition at marker %s sample %s", CHAR(STRING_ELT(ids, i)), CHAR(STRING_ELT(individuals, j)));
            if (a < 0 || b < 0) { missing = 1; continue; }
            first[i + (size_t)nr * j] = (unsigned char)a;
            second[i + (size_t)nr * j] = (unsigned char)b;
            seen[a] = seen[b] = 1;
        }
        int states = missing;
        for (int j = 0; j < 254; j++) states += seen[j];
        if (states > (int)GS_ALLELE_MAX + 1) error("too many observed allele states at marker %s", CHAR(STRING_ELT(ids, i)));
        UNPROTECT(np - before); np = before;
        st_work_reset(row_mark);
    }
    const char *keys[] = {"CHROM", "POS", "ID", "REF", "ALT", "individual", "allele1", "allele2", "CHROM.levels"};
    SEXP out = st_vcf_keep(allocVector(VECSXP, 9), &np),
         out_names = st_vcf_keep(allocVector(STRSXP, 9), &np),
         codes = st_vcf_keep(allocVector(INTSXP, nr), &np),
         chromosome_levels = st_vcf_keep(allocVector(STRSXP, nchrom), &np),
         a1 = st_vcf_keep(allocMatrix(RAWSXP, nr, ns), &np),
         a2 = st_vcf_keep(allocMatrix(RAWSXP, nr, ns), &np);
    memcpy(INTEGER(codes), chrom, (size_t)nr * sizeof(int));
    memcpy(RAW(a1), first, (size_t)nr * ns); memcpy(RAW(a2), second, (size_t)nr * ns);
    for (int i = 0; i < nchrom; i++) SET_STRING_ELT(chromosome_levels, i, levels[i]);
    SEXP fields[] = {codes, positions, ids, refs, alts, individuals, a1, a2, chromosome_levels};
    for (int j = 0; j < 9; j++)
    {
        SET_VECTOR_ELT(out, j, fields[j]); SET_STRING_ELT(out_names, j, mkChar(keys[j]));
    }
    setAttrib(out, R_NamesSymbol, out_names);
    setAttrib(out, R_ClassSymbol, st_vcf_keep(mkString("STvcf"), &np));
    setAttrib(out, install("phased"), st_vcf_keep(ScalarLogical(phased && !unphased), &np));
    st_vcf_metadata(data, out);
    UNPROTECT(np);
    return out;
}

static int st_vcf_attributes(SEXP value)
{
#if R_VERSION >= R_Version(4, 6, 0)
    return R_getAttribCount(value) != 0;
#else
    return ATTRIB(value) != R_NilValue;
#endif
}

static SEXP st_vcf_to_frame(SEXP args)
{
    SEXP value = getListElement(args, "data");
    int np = 0;
    if (TYPEOF(value) == LISTSXP)
    {
        SEXP original = value;
        value = st_vcf_keep(coerceVector(value, VECSXP), &np);
        copyMostAttrib(original, value);
    }
    if (!st_stvcf_components(value)) error("STvcf must contain exactly the nine required components");
    SEXP phase = getAttrib(value, install("phased"));
    if (!isNull(phase) && (TYPEOF(phase) != LGLSXP || LENGTH(phase) != 1 || LOGICAL(phase)[0] == NA_LOGICAL))
        error("invalid STvcf phased attribute");
    int phased = !isNull(phase) && LOGICAL(phase)[0] == TRUE;
    SEXP chrom = st_stvcf_element(value, "CHROM"), pos = st_stvcf_element(value, "POS"),
         ids = st_stvcf_element(value, "ID"), refs = st_stvcf_element(value, "REF"),
         alts = st_stvcf_element(value, "ALT"), individuals = st_stvcf_element(value, "individual"),
         levels = st_stvcf_element(value, "CHROM.levels"), a1 = st_stvcf_element(value, "allele1"),
         a2 = st_stvcf_element(value, "allele2");
    int nr = LENGTH(ids), ns = LENGTH(individuals);
    if (!nr || !ns) error("STvcf contains no markers or no individuals");
    if (TYPEOF(chrom) != INTSXP || LENGTH(chrom) != nr || TYPEOF(pos) != REALSXP || LENGTH(pos) != nr ||
        TYPEOF(ids) != STRSXP || TYPEOF(refs) != STRSXP || LENGTH(refs) != nr ||
        TYPEOF(alts) != STRSXP || LENGTH(alts) != nr || TYPEOF(individuals) != STRSXP || TYPEOF(levels) != STRSXP)
        error("invalid STvcf component types or lengths");
    for (int i = 0; i < nr; i++)
        if (!R_FINITE(REAL(pos)[i]) || REAL(pos)[i] < 0) error("invalid genetic position in STvcf");
    int valid = st_vcf_unique(ids);
    for (int i = 0; i < nr; i++) valid &= st_vcf_identifier(STRING_ELT(ids, i), 255);
    if (!valid) error("invalid or duplicated marker ID in STvcf");
    valid = st_vcf_unique(individuals);
    for (int i = 0; i < ns; i++) valid &= st_vcf_identifier(STRING_ELT(individuals, i), R_PosInf);
    if (!valid) error("invalid or duplicated individual ID in STvcf");
    for (int i = 0; i < nr; i++)
        if (STRING_ELT(refs, i) == NA_STRING || !*CHAR(STRING_ELT(refs, i)) ||
            !strcmp(CHAR(STRING_ELT(refs, i)), ".") || strchr(CHAR(STRING_ELT(refs, i)), ',') ||
            STRING_ELT(alts, i) == NA_STRING || !*CHAR(STRING_ELT(alts, i))) error("invalid REF or ALT field in STvcf");
    SEXP d1 = getAttrib(a1, R_DimSymbol), d2 = getAttrib(a2, R_DimSymbol);
    if (st_vcf_attributes(d1) || st_vcf_attributes(d2) || TYPEOF(d1) != INTSXP || LENGTH(d1) != 2 || INTEGER(d1)[0] != nr || INTEGER(d1)[1] != ns ||
        TYPEOF(d2) != INTSXP || LENGTH(d2) != 2 || INTEGER(d2)[0] != nr || INTEGER(d2)[1] != ns)
        error("STvcf allele matrices have invalid dimensions");
    if (!((TYPEOF(a1) == RAWSXP && TYPEOF(a2) == RAWSXP) || (TYPEOF(a1) == INTSXP && TYPEOF(a2) == INTSXP)))
        error("STvcf allele matrices must both use raw or both use integer storage");
    int nlevels = LENGTH(levels);
    valid = nlevels > 0 && st_vcf_unique(levels);
    for (int i = 0; i < nlevels; i++)
        if (STRING_ELT(levels, i) == NA_STRING || !*CHAR(STRING_ELT(levels, i))) valid = 0;
    for (int i = 0; i < nr; i++)
        if (INTEGER(chrom)[i] == NA_INTEGER || INTEGER(chrom)[i] < 1 || INTEGER(chrom)[i] > nlevels) valid = 0;
    if (!valid) error("invalid chromosome coding in STvcf");
    unsigned char *used = st_work_alloc(nlevels, 1); memset(used, 0, (size_t)nlevels);
    for (int i = 0; i < nr; i++) used[INTEGER(chrom)[i] - 1] = 1;
    for (int i = 0; i < nlevels; i++)
        if (!used[i]) error("STvcf chromosome codes must be consecutive and used");
    int *nalt = st_work_alloc(nr, sizeof(int)), too_many = 0;
    for (int i = 0; i < nr; i++)
    {
        nalt[i] = st_vcf_alt(STRING_ELT(alts, i), STRING_ELT(refs, i));
        if (nalt[i] < 0) error("invalid or duplicated ALT field in STvcf marker %s", CHAR(STRING_ELT(ids, i)));
        too_many |= nalt[i] > 253;
    }
    if (too_many) error("STvcf supports at most 253 ALT alleles per marker");
    int raw = TYPEOF(a1) == RAWSXP;
    for (int i = 0; i < nr; i++)
    {
        unsigned char seen[254] = {0}; int missing = 0;
        for (int j = 0; j < ns; j++)
        {
            size_t at = i + (size_t)nr * j;
            int first = raw ? RAW(a1)[at] : INTEGER(a1)[at], second = raw ? RAW(a2)[at] : INTEGER(a2)[at],
                m1 = raw ? first == 255 : first == NA_INTEGER, m2 = raw ? second == 255 : second == NA_INTEGER;
            if ((!m1 && (first < 0 || first > nalt[i])) || (!m2 && (second < 0 || second > nalt[i])))
                error("invalid allele index in STvcf marker %s", CHAR(STRING_ELT(ids, i)));
            if (m1 || m2) missing = 1;
            else seen[first] = seen[second] = 1;
        }
        int states = missing;
        for (int j = 0; j < 254; j++) states += seen[j];
        if (states > (int)GS_ALLELE_MAX + 1) error("too many observed allele states at marker %s", CHAR(STRING_ELT(ids, i)));
    }
    SEXP chromosome = st_vcf_keep(allocVector(STRSXP, nr), &np),
         format = st_vcf_keep(allocVector(STRSXP, nr), &np),
         level_names = getAttrib(levels, R_NamesSymbol), chromosome_names = R_NilValue;
    if (!isNull(level_names)) chromosome_names = st_vcf_keep(allocVector(STRSXP, nr), &np);
    for (int i = 0; i < nr; i++)
    {
        SET_STRING_ELT(chromosome, i, STRING_ELT(levels, INTEGER(chrom)[i] - 1));
        SET_STRING_ELT(format, i, mkChar("GT"));
        if (!isNull(chromosome_names)) SET_STRING_ELT(chromosome_names, i, STRING_ELT(level_names, INTEGER(chrom)[i] - 1));
    }
    if (!isNull(chromosome_names)) setAttrib(chromosome, R_NamesSymbol, chromosome_names);
    SEXP tail = st_vcf_keep(allocList(8 + ns), &np),
         call = st_vcf_keep(lcons(findFun(install("data.frame"), R_BaseEnv), tail), &np);
    const char *keys[] = {"CHROM", "POS", "ID", "REF", "ALT", "FORMAT"};
    SEXP fields[] = {chromosome, pos, ids, refs, alts, format};
    for (int j = 0; j < 6; j++, tail = CDR(tail))
    {
        SETCAR(tail, fields[j]); SET_TAG(tail, install(keys[j]));
    }
    for (int j = 0; j < ns; j++, tail = CDR(tail))
    {
        SEXP column = allocVector(STRSXP, nr);
        SETCAR(tail, column); SET_TAG(tail, install("st.vcf.sample"));
        for (int i = 0; i < nr; i++)
        {
            size_t at = i + (size_t)nr * j;
            int first = raw ? RAW(a1)[at] : INTEGER(a1)[at], second = raw ? RAW(a2)[at] : INTEGER(a2)[at];
            char text[32];
            if (raw ? first == 255 || second == 255 : first == NA_INTEGER || second == NA_INTEGER)
                snprintf(text, sizeof(text), ".%c.", phased ? '|' : '/');
            else snprintf(text, sizeof(text), "%d%c%d", first, phased ? '|' : '/', second);
            SET_STRING_ELT(column, i, mkChar(text));
        }
    }
    SETCAR(tail, st_vcf_keep(ScalarLogical(FALSE), &np)); SET_TAG(tail, install("stringsAsFactors"));
    tail = CDR(tail);
    SETCAR(tail, st_vcf_keep(ScalarLogical(FALSE), &np)); SET_TAG(tail, install("check.names"));
    SEXP out = st_vcf_keep(eval(call, R_BaseEnv), &np),
         out_names = st_vcf_keep(duplicate(getAttrib(out, R_NamesSymbol)), &np);
    for (int j = 0; j < ns; j++) SET_STRING_ELT(out_names, LENGTH(out_names) - ns + j, STRING_ELT(individuals, j));
    setAttrib(out, R_NamesSymbol, out_names);
    st_vcf_metadata(value, out);
    UNPROTECT(np);
    return out;
}

SEXP st_legacy_marker(SEXP args)
{
    const char *action =
        st_legacy_marker_string(getListElement(args, "action"), "action");
    if (!strcmp(action, "vcf.identifiers")) return st_vcf_identifiers(args);
    if (!strcmp(action, "vcf.from.frame")) return st_vcf_from_frame(args);
    if (!strcmp(action, "vcf.to.frame")) return st_vcf_to_frame(args);
    if (!strcmp(action, "map.write"))
        return st_legacy_map_write(args);
    if (!strcmp(action, "legacy.result"))
        return st_legacy_result_call(args);
    if (!strcmp(action, "hblocks.set"))
        return st_legacy_hblocks(args);
    if (!strcmp(action, "cross.va"))
        return st_legacy_cross_va(args);
    if (!strcmp(action, "cross.info"))
        return st_legacy_cross_info(args);
    if (!strcmp(action, "cross.es"))
        return st_legacy_cross_es(args);
    if (!strcmp(action, "cross.allocate"))
        return st_legacy_cross_allocate(args);
    if (!strcmp(action, "cross.allocate.print"))
        return st_legacy_cross_allocate_print(args);
    if (!strcmp(action, "call"))
        return st_legacy_buffer_call(args);
    if (!strcmp(action, "il.lines"))
        return st_legacy_il_lines(args);
    if (!strcmp(action, "il.create"))
        return st_legacy_il_create(args);
    if (!strcmp(action, "select.phen"))
        return st_legacy_select_phenotype(args);
    if (!strcmp(action, "pseff"))
        return st_legacy_pseff(args);
    if (!strcmp(action, "mark.alleles"))
        return st_legacy_mark_alleles(args);
    if (!strcmp(action, "il.summary"))
        return st_legacy_il_summary(args);
    error("Unknown legacy marker operation");
    return R_NilValue;
}
