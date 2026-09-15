/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   Simulation routines by Hans Peter Mauerer and Matthias Frisch
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#define AUTHOR "$Author: frisch-m $"
#define DATE "$Date: 2019/06/10 00:00:00 $"
#define REVISION "$Revision: 19.2 $"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include <stdint.h>
#include <limits.h>
#include <float.h>

#include "SelectionToolsAllele.h"

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Random.h>

#define MAXTHREADS 256

#include "MathLib12.h"

#if MAXTHREADS != ST_RNG_MAXTHREADS
#error "MAXTHREADS and ST_RNG_MAXTHREADS differ"
#endif


#ifdef _OPENMP
  #ifndef ST
     #define MT
  #endif
#endif

#ifdef MT
   #include <omp.h> 
#else
  #define omp_get_thread_num() 0
  #define omp_set_num_threads(thr) {}
  #define omp_get_num_threads() 0
  #define omp_set_dynamic(thr) {}
#endif

/* not used

#define EXPYEAR1 "2017"
#define EXPYEAR2 "2015"
#define EXPYEAR3 "2016"

*/

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* FUNCTION DECLARATIONS */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void info( int level, char*info, ... );
void remove_genotype_population( char **PopNames );
void remove_evaluate_population( char **PopNames );
void flush_population( char**PopName, int mode );
void effmap_remove_all( void );

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* RANDOM NUMBERS */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int rng_initialized = 0;

#ifdef use_gsl

   #include <gsl/gsl_rng.h> 
   #include <gsl/gsl_randist.h>

   gsl_rng *rng[MAXTHREADS];

   const gsl_rng_type **rng_default; // pointer to array element
   const gsl_rng_type **all_types;   // array

   char rng_def_nme [] ="taus";
   char * pp_rng_def_nme = rng_def_nme;
   char ** p_rng_def_nme = &pp_rng_def_nme;

   void rng_init (void);

   void rng_info(void)
   {
     info(0,"Random number generator type: '%s'\n",gsl_rng_name(rng[0]));
   }

   void rng_choose ( char ** name)
   {
     for ( rng_default=all_types; *rng_default != 0; rng_default++)
       {
	 if( 0== strcmp ( (*rng_default)->name, *name) ) break;
       }
     if (*rng_default == 0) rng_choose( p_rng_def_nme);
     rng_init();
     rng_info();
   }


   void rng_list(void)
   {
     const gsl_rng_type**t;
     
     info(0,"Available generators:\n");
     for( t = all_types; *t != 0; t++)
       {
       if (0==strcmp(gsl_rng_name(rng[0]),(*t)->name))
	 {
	   info(0,"--> %s <--\n",(*t)->name);
	 }
       else
	 {
	   info(0,"    %s\n",(*t)->name);
	 }
     }
   }

#endif

#ifndef use_gsl

/* -------------------------------------------------------------------------
   Clean-room, license-safe drop-in replacement for the old NR-derived code.

   Compatibility goals (as requested):
   - same macros (IA..RNMX) and same #undef list
   - same typedef gsl_rng { long seed; long iyy; long ivv[NTAB]; }
   - same global: gsl_rng* rng[MAXTHREADS];
   - same function names + signatures, all present
   - keep "seed negative forces (re)initialization" semantics
   - thread-safe: no shared mutable static state

   Notes:
   - RNG core uses Park–Miller "minimal standard" with Schrage method (widely
     published; implemented here independently).
   - gammln() is implemented via standard C lgamma() (no NR code).
   - Poisson sampler:
       * Knuth method for small mean
       * PTRS (Hoermann) for large mean
     No static caching, so thread-safe.

   Required headers:
     #include <math.h>
     #include <stdlib.h>
   ------------------------------------------------------------------------- */

   #define IA 16807
   #define IM 2147483647
   #define AM (1.0/IM)
   #define IQ 127773
   #define IR 2836
   #define NTAB 32
   #define NDIV (1+(IM-1)/NTAB)
   #define EPS 1.2e-7
   #define RNMX (1.0-EPS)

   typedef struct {
     long seed       ;
     long iyy        ;
     long ivv[NTAB]  ;
   } gsl_rng;

   gsl_rng* rng[MAXTHREADS];

   /* one Park–Miller step using Schrage method; returns new state in (0, IM) */
   static inline long pm_step(long x)
   {
     /* x is assumed positive */
     const long k  = x / IQ;
     long xnew = IA * (x - k * IQ) - IR * k;
     if (xnew < 0) xnew += IM;
     return xnew;
   }

   double gsl_rng_uniform (gsl_rng* rng)
   {
     /* This preserves the original state layout/behavior:
        - rng->seed holds the evolving Park–Miller state
        - rng->ivv[] is a shuffle table
        - rng->iyy holds the last shuffle output
        - negative seed or iyy==0 triggers initialization
     */
     int j;
     long *idum = &(rng->seed);
     long *iy   = &(rng->iyy);
     long *iv   = &(rng->ivv[0]);

     if ((*idum) <= 0 || (*iy) == 0)
     {
       long x = *idum;

       /* seed normalization: treat 0 as 1; flip negative to positive */
       if (x == 0) x = 1;
       if (x < 0)  x = -x;

       /* warm up and fill shuffle table */
       for (j = NTAB + 7; j >= 0; --j)
       {
         x = pm_step(x);
         if (j < NTAB) iv[j] = x;
       }

       *iy   = iv[0];
       *idum = x;
     }

     /* advance */
     *idum = pm_step(*idum);

     /* shuffle */
     j = (int)((*iy) / NDIV);
     if (j < 0) j = 0;
     if (j >= NTAB) j = NTAB - 1;

     *iy   = iv[j];
     iv[j] = *idum;

     /* scale to (0,1); clamp like original */
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
    /* Thread-safe, standards-based replacement */
    return lgamma(xx);
  }

  #define rPI 3.141592654

  double gsl_rng_uniform_pos (gsl_rng* rng);

  unsigned int gsl_ran_poisson (gsl_rng* rng, double xm)
  {
    /* Thread-safe: no static caches.

       Small mean: Knuth multiplicative method.
       Large mean: PTRS (Hörmann, transformed rejection).
    */
    if (!(xm >= 0.0) || !isfinite(xm)) return 0U;
    if (xm == 0.0) return 0U;

    /* Small mean: Knuth */
    if (xm < 12.0)
    {
      const double L = exp(-xm);
      unsigned int k = 0U;
      double p = 1.0;

      do {
        ++k;
        p *= gsl_rng_uniform_pos(rng);
      } while (p > L);

      return k - 1U;
    }

    /* Large mean: PTRS */
    {
      const double mu = xm;
      const double s  = sqrt(mu);
      const double b  = 0.931 + 2.53 * s;
      const double a  = -0.059 + 0.02483 * b;
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

        /* acceptance test on log scale */
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

  double gsl_rng_uniform_pos (gsl_rng* rng)
  {
    /* ensure strict (0,1) */
    double r;
    do {
      r = gsl_rng_uniform(rng);
    } while (r <= 0.0 || r >= 1.0);
    return r;
  }

unsigned long int gsl_rng_uniform_int (gsl_rng* rng, unsigned long int n)
{
  const unsigned long int range = 2147483646UL;
  unsigned long int x, limit, q, r, h;

  if (n <= 1UL) return 0UL;

  if (n <= range)
  {
    q = range / n;
    limit = q * n;

    do
    {
      (void) gsl_rng_uniform(rng);
      x = (unsigned long int)(rng->iyy - 1L);
    }
    while (x >= limit);

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
  }
  while ((h == q) && (x >= r));

  return h * range + x;
}


  gsl_rng* gsl_rng_alloc(void)
  {
    /* allocate and initialize to a valid default state */
    gsl_rng *r = (gsl_rng*) malloc(sizeof(gsl_rng));
    if (!r) return NULL;
    r->seed = -1; /* negative => init on first use */
    r->iyy  = 0;
    for (int i = 0; i < (int)(sizeof(r->ivv)/sizeof(r->ivv[0])); ++i) r->ivv[i] = 0;
    return r;
  }

  void gsl_rng_set (gsl_rng * rng, unsigned long int s)
  {
    /* preserve original semantics: store negative seed to force init */
    if (s == 0UL) s = 1UL;
    rng->seed = -(long) s;
    rng->iyy  = 0;
  }

  void gsl_rng_free (gsl_rng * r)
  {
    free(r);
  }

  void rng_choose(char**name_rng)
  {
    (void)name_rng;
    info(0, "Random number generator: Park-Miller (clean-room) + shuffle table\n" );
  }

  void rng_info(void)
  {
    info(0, "Random number generator type: Park-Miller (clean-room) + shuffle table\n");
  }

  void rng_list(void)
  {
    info(0, "Available generators: Park-Miller (clean-room) + shuffle table\n" );
  }

  #define p_rng_default

#endif

void rng_init (void) {
 
  int i, j;
  unsigned long int seed;
  gsl_rng *tmp_rng[MAXTHREADS] = { NULL };

  #ifdef use_gsl
    if ( 0 == rng_initialized )
      {
	all_types = gsl_rng_types_setup();
	
	for ( rng_default=all_types; *rng_default != 0; rng_default++)
	  {
	    if( 0== strcmp ( (*rng_default)->name, rng_def_nme) ) break;
	  }
      }
    const gsl_rng_type * p_rng_default = *rng_default;
  #endif

  for (i =0; i < MAXTHREADS; i++ )
    {
      seed = (unsigned long)st_rand() + 1UL; // Using MathLib12
      tmp_rng[i] = gsl_rng_alloc (p_rng_default);

      if ( NULL == tmp_rng[i] )
	goto fail;

      gsl_rng_set (tmp_rng[i],seed);
    }

  if ( 1 == rng_initialized)
    {
      for (i =0; i < MAXTHREADS; i++ )
	gsl_rng_free (rng[i]);
      info(1, "Re-initializing random number generator \n" );
    }
  else
    {
      info(1, "Initializing random number generator \n" );
    }

  for (i =0; i < MAXTHREADS; i++ )
    rng[i] = tmp_rng[i];

  rng_initialized = 1;
  return;

fail:

  for (j =0; j < i; j++ )
    gsl_rng_free (tmp_rng[j]);

  info(-2, "Memory allocation failed in rng_init\n" );
  return;

}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* PLABSIM MAIN  */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#define MAXNAME 256
#define MAXLINE 100
#define CPOS double
#define POPULATION long int
#define MAP_DISPERSE_FACTOR_DEFAULT 100.0

#define CCTIME __TIME__
#define CCDATE __DATE__

/* #define xCalloc(...)Calloc(__VA_ARGS__) */
/* #define xRealloc(...)Realloc(__VA_ARGS__) */
/* #define xFree(...)Free(__VA_ARGS__) */

/* Alternative memory allocation without R wrappers */

#define yCalloc(n,t)(t*)malloc((size_t)(n)*sizeof(t))

static void *yRealloc_impl(void *ptr, size_t count, size_t size)
{
   void *new_ptr;

   if ((count == 0) || (size == 0))  goto fail;
   if (count > SIZE_MAX / size)      goto fail;

   new_ptr = realloc(ptr, count * size);

   if (new_ptr == NULL)      goto fail;

   return new_ptr;

fail:
   free(ptr);
   return NULL;
}

#define yRealloc(p,n,t) \
   (t*)yRealloc_impl((void*)(p),(size_t)(n),sizeof(t))

#define yFree(p)(free((void*)(p)),(p)= NULL)


/* Type definitions */

typedef struct
{
    CPOS Pos;
    int  chrom;
    char *name;
    char *class;
    int  counter;
} 
TMappoint;

typedef struct
{
    double    ChromLen;
    int       NrMappoints;
    TMappoint *Mappoint;
} 
TChromDescr;

typedef struct
{
  CPOS           Pos;
  unsigned short All;
}
TLoc;

typedef struct
{
  int  NoLociAlloc;
  int  NoLociUsed;
  TLoc *Loc;
}
THom;

typedef struct
{
  THom *Hom;
}
TChrom;

typedef struct
{
  TChrom *Chrom;
  float  *GValue;
  float  *PValue;
  int    *genotype;
  char   *info;
}
TInd;

typedef struct
{
  char       Name[ MAXNAME ];
  POPULATION NoIndsAlloc;
  POPULATION NoInds;
  TInd *Ind;
}
TPop;

typedef struct
{
    int  sizeRows, sizeCols;
    int  nPop, nMarker;
    int  *nAll, *nInd;
    int  *nAll_off, *nInd_off;
    int  *data;
}
TPopMatrix;

struct TIntLinkList
{
    int number;
    struct TIntLinkList *child, *parent;
    struct TIntLinkList *next,  *prev;
};

typedef struct TIntLinkList IntNode;

typedef struct SPopListEl
{
    TPop* Pop;
    struct SPopListEl *Prev;
    struct SPopListEl *Next;
}
TPopListEl;

typedef struct
{
    TMappoint *MP;
    int  NrAll;
    int  *All;
}
TGen;

typedef struct
{
    float Eff;
    int   NrLoci;
    TGen  *Loc;
}
TEff;

typedef struct
{
    float GValue_avg;
    float weight;
    int   NrEff;
    int   Pos;
    TEff  *Eff;
    char  Name[ 3*MAXNAME ];
    int   MPs;
    TMappoint**MPset;
}
TEffects;


typedef struct
{
    int  betrag;
    int  *member;
}
TSet;

/* Variable declarations  */

int     NoChroms = 0;
int     NoHoms = 0;
int     NoLociInit = 10;
IntNode *startNode = NULL;
double  COFreq = 1;

TPopListEl *PopListEl = NULL;
int        bg, missing, doubleHetero;
int        dimLA, dimLB, locusA, locusB, nrGeno;
int        **table, *genotype;


TPopListEl  *PopListStart = NULL;
TChromDescr *PopDescr = NULL;
int         cMP = 0;
int         cMPhom = 0;
int         INFOLEVEL = 0;

TEffects*Effects = NULL;
int NrEffects = 0;
TSet*AllSet;
unsigned long int MDP = 0;

/* Generic functions */

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

void info( int level, char*info, ... )
{
    if ( INFOLEVEL >= level )
    {
        /*
         * Never call the R console from an OpenMP region.
         *
         * Ordinary information and warnings are simply suppressed.
         * Errors are remembered so that the calling routine can report
         * one error after leaving the parallel region.
         */
        if ( !st_console_output_allowed() )
        {
            if ( -2 >= level )
                st_parallel_info_error_set();

            return;
        }

        va_list ap;
        char Info[ 200 ];

        va_start( ap, info );
        vsnprintf( Info, sizeof( Info ), info, ap );
        va_end( ap );

        if ( 0 <= level )
            Rprintf( "M: %s", Info );
        else if ( -1 == level )
            Rprintf( "W: %s", Info );
        else if ( -2 == level )
            Rprintf( "E: %s", Info );
    }
}


void r_info(int*level, char**message)
{
  info( *level, "%s\n", *message );
}

void plabsim_info_cctime_get( char**message )
{
  snprintf(*message,
	   strlen(CCTIME) + strlen(CCDATE) + 4,
	   "%s - %s", CCTIME, CCDATE);
}

 
void print_version(void)
{
    info( 0, "%s (%s): Rev. %s \n", AUTHOR, DATE, REVISION );
    info( 0, "compiled: %s - %s\n", CCTIME, CCDATE );
}

void print_version_2( char **message )
{
  char * co_date_complete = DATE;
  char co_date[20];
  strncpy(co_date,co_date_complete+7,10);
  co_date[10] = '\0';
  info( 0, "plabsim c %s\n", co_date );
  info( 0, "plabsim r %s\n", *message );
}

FILE *fileopen ( char* path,
		 char* mode)
{
    FILE*fp;

    if ( ( fp = fopen( path, mode ) ) == NULL )
        info( -2, "Can't open the file %s\n", path );

    return fp;
}

void fileclose( FILE*fp )
{
    if ( fclose( fp ) != 0 )
        info( -2, "Can't close the file stream\n");
}


void set_info_level( int*level )
{
  INFOLEVEL = *level;
  if ( INFOLEVEL < -2 )  INFOLEVEL = -2;
  if ( INFOLEVEL > 10 )  INFOLEVEL = 10;
  info(2,"Info level set to %d\n",INFOLEVEL);
}

void get_info_level (int* i)
{
  *i = INFOLEVEL;
}

/* Simulations functions */

void find_list_element ( char        *Name, 
			 TPopListEl  **PopListEl )
{
    *PopListEl = PopListStart;

    while ( ( *PopListEl != NULL ) && 
	    ( strcmp( ( *PopListEl ) ->Pop->Name, Name ) != 0 ) )
      {
	( *PopListEl ) = ( *PopListEl ) ->Next;
      }
}


void remove_individuum ( TInd *DelInd )
{
    TChrom  VDelChrom;
    TChrom *DelChrom = &VDelChrom;
    int chrom, hom;

    for ( chrom = 0; chrom < NoChroms; chrom++ ) 
    {
      DelChrom = &(DelInd->Chrom[chrom]);
      for ( hom = 0; hom < NoHoms; hom++ ) yFree (DelChrom->Hom[hom].Loc);
      yFree (DelChrom->Hom);
    }

    yFree (DelInd->Chrom);

    if ( NULL != DelInd->GValue )   yFree (DelInd->GValue);
    if ( NULL != DelInd->PValue )   yFree (DelInd->PValue);
    if ( NULL != DelInd->genotype ) yFree (DelInd->genotype);
    if ( NULL != DelInd->info )     yFree (DelInd->info);
}

int new_individuum ( TInd *NewInd )
{
    int    chrom, hom;
    THom   VNewHom;
    THom   *NewHom = &VNewHom;
    TChrom VNewChrom;
    TChrom *NewChrom = &VNewChrom;

    NewInd->Chrom    = NULL;
    NewInd->GValue   = NULL;
    NewInd->PValue   = NULL;
    NewInd->genotype = NULL;
    NewInd->info     = NULL;

    NewInd->Chrom = (TChrom*) yCalloc( NoChroms, TChrom );

    if ( NULL == NewInd->Chrom )
        goto fail;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
        NewInd->Chrom[ chrom ].Hom = NULL;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        NewChrom = &( NewInd->Chrom[ chrom ] );
        NewChrom->Hom = ( THom* ) yCalloc( NoHoms, THom );

        if ( NULL == NewChrom->Hom )
            goto fail;

        for ( hom = 0;hom < NoHoms;hom++ )
            NewChrom->Hom[ hom ].Loc = NULL;

        for ( hom = 0;hom < NoHoms;hom++ )
        {
            NewHom = &( NewChrom->Hom[ hom ] );
            NewHom->NoLociAlloc = NoLociInit;
            NewHom->NoLociUsed = 0;
            NewHom->Loc = ( TLoc* ) yCalloc( NoLociInit, TLoc );

            if ( NULL == NewHom->Loc )
                goto fail;
        }
    }

    return 0;

fail:

    if ( NULL != NewInd->Chrom )
    {
        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            NewChrom = &( NewInd->Chrom[ chrom ] );

            if ( NULL != NewChrom->Hom )
            {
                for ( hom = 0;hom < NoHoms;hom++ )
                    yFree( NewChrom->Hom[ hom ].Loc );
                yFree( NewChrom->Hom );
            }
        }
        yFree( NewInd->Chrom );
    }

    return -1;
}


TInd* get_individual ( TPop* Pop, POPULATION nrInd )
{
    return &( Pop->Ind[ nrInd ] );
}


static void st_clear_population_genotype(TPop *Pop)
{
  POPULATION ind;

  if (NULL == Pop)
    return;

  for (ind = 0; ind < Pop->NoIndsAlloc; ind++)
    {
      free(Pop->Ind[ind].genotype);
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
      free(Pop->Ind[ind].PValue);
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
      free(Pop->Ind[ind].GValue);
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
      free(*target);
      *target = NULL;
      return 1;
    }

  n = (size_t)NrEffects + 1;

  if (n > SIZE_MAX / sizeof(*tmp))
    return 0;

  tmp = (float *)malloc(n * sizeof(*tmp));

  if (NULL == tmp)
    return 0;

  memcpy(tmp,source,n * sizeof(*tmp));
  free(*target);
  *target = tmp;

  return 1;
}


static int st_effect_position(char **EffName, int *eff);
static int st_selection_criterion(char **criterion);
static float *st_individual_value_array(TInd *Ind, int criterion);
static int st_population_selection_ready(
    TPop *Pop,
    int criterion,
    int eff);



void set_nolociinit ( int* NewLociInit)
{
  if ( 0 < *NewLociInit )  NoLociInit = *NewLociInit;
}


void remove_map( void )
{
  TPopListEl *PopListEl;
  POPULATION ind;
  int chrom, i;
  int mapRemoved = 0;

  PopListEl = PopListStart;

  while ( NULL != PopListEl )
    {
      for ( ind = 0;
	    ind < PopListEl->Pop->NoIndsAlloc;
	    ind++ )
	{
	  free( PopListEl->Pop->Ind[ ind ].GValue );
	  PopListEl->Pop->Ind[ ind ].GValue = NULL;

	  free( PopListEl->Pop->Ind[ ind ].PValue );
	  PopListEl->Pop->Ind[ ind ].PValue = NULL;

	  free( PopListEl->Pop->Ind[ ind ].genotype );
	  PopListEl->Pop->Ind[ ind ].genotype = NULL;
	}

      PopListEl = PopListEl->Next;
    }

  effmap_remove_all();

  if ( NULL != PopDescr )
    {
      for ( chrom = 0;chrom < NoChroms;chrom++ )
	{
	  if ( 0 < PopDescr[ chrom ].NrMappoints )
	    {
	      mapRemoved = 1;

	      for ( i = 0;
		    i < PopDescr[ chrom ].NrMappoints;
		    i++ )
		{
		  free(
		    PopDescr[ chrom ].Mappoint[ i ].name );

		  free(
		    PopDescr[ chrom ].Mappoint[ i ].class );
		}

	      free( PopDescr[ chrom ].Mappoint );

	      PopDescr[ chrom ].Mappoint = NULL;
	      PopDescr[ chrom ].NrMappoints = 0;
	    }
	}
    }

  cMP = 0;
  cMPhom = 0;

  if ( mapRemoved )
    info( 1, "Linkage map removed\n" );
}

static int set_genome_par_impl ( int*    _NoChroms,
		      int*    _NoHoms,
		      double* _ChromLen )
{
  time_t ltime;
  char tme[ 64 ];
  TChromDescr *tmpPopDescr;
  int newNoChroms = *_NoChroms;
  int newNoHoms = *_NoHoms;
  int chrom;

  time( &ltime );
  strcpy( tme, ctime( &ltime ) );

  /* Not used
  if (
      (strstr(tme,EXPYEAR1)==NULL)&&
      (strstr(tme,EXPYEAR2)==NULL)&&
      (strstr(tme,EXPYEAR3)==NULL)
      )
    {
      info(-2,"Compilation date: %s - %s\n", CCTIME, CCDATE );
      info(-2,"Software update required: %s\n",tme);
      return ;
    }
  */

  if ( newNoChroms <= 0 )
    {
      info( -2,
	    "Number of chromosomes must be positive\n" );
      return 0;
    }

  if ( newNoHoms <= 0 )
    {
      info( -2,
	    "Number of homologues must be positive\n" );
      return 0;
    }

  for ( chrom = 0;chrom < newNoChroms;chrom++ )
    {
      if ( !isfinite( _ChromLen[ chrom ] ) ||
	   ( _ChromLen[ chrom ] <= 0 ) )
	{
	  info(
	    -2,
	    "Chromosome length must be finite and positive\n" );
	  return 0;
	}
    }

  if ( NULL != PopListStart )
    {
      info( -2, "First remove all populations\n" );
      return 0;
    }

  if ( NrEffects > 0 )
    {
      info( -2, "First remove all effects\n" );
      return 0;
    }

  if ( ( PopDescr != NULL ) && ( cMP > 0 ) )
    {
      info(
	-2,
	"First remove previously loaded linkage map\n" );
      return 0;
    }

  if ( ( size_t ) newNoChroms >
       SIZE_MAX / sizeof( TChromDescr ) )
    {
      info(
	-2,
	"Memory allocation failed in set_genome_par\n" );
      return 0;
    }

  if ( NULL == PopDescr )
    tmpPopDescr = ( TChromDescr* ) malloc(
      ( size_t ) newNoChroms *
      sizeof( TChromDescr ) );
  else
    tmpPopDescr = ( TChromDescr* ) realloc(
      PopDescr,
      ( size_t ) newNoChroms *
      sizeof( TChromDescr ) );

  if ( NULL == tmpPopDescr )
    {
      info(
	-2,
	"Memory allocation failed in set_genome_par\n" );
      return 0;
    }

  PopDescr = tmpPopDescr;
  NoChroms = newNoChroms;
  NoHoms = newNoHoms;

  for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
      PopDescr[ chrom ].ChromLen =
	_ChromLen[ chrom ];

      PopDescr[ chrom ].NrMappoints = 0;
      PopDescr[ chrom ].Mappoint = NULL;
    }

  info( 0,
	"%d chromosomes defined\n",
	NoChroms );

  return 1;
}


void set_genome_par ( int*    _NoChroms,
		      int*    _NoHoms,
		      double* _ChromLen )
{
  ( void ) set_genome_par_impl(
    _NoChroms,
    _NoHoms,
    _ChromLen );
}


void optimize_population ( char** PopNames )
{
  int         chrom, hom;
  POPULATION  ind;
  TPopListEl* PopListEl;
  TPop*       Pop = NULL;
  char*       PopName;

  PopName = strtok( *PopNames, "[' ','\n']" );
  
  while ( NULL != PopName )
    {
      find_list_element( PopName, &PopListEl );
      
      if ( PopListEl == NULL ) {
	info( -2, "Population not defined: %s\n", PopName );
	return;
      }
      
      Pop = PopListEl->Pop;
      
      if ( NULL != Pop->Ind )
	{
	  for ( ind = 0; ind < Pop->NoInds; ind++ )
	    for ( chrom = 0; chrom < NoChroms; chrom++ )
	      for ( hom = 0; hom < NoHoms; hom++ )
		{
		  THom*Hom = &( ( *Pop ).Ind[ind].Chrom[chrom].Hom[hom] );
		  
		  if ( ( Hom->NoLociUsed + NoLociInit ) 
		       <= Hom->NoLociAlloc                )
		    {
		      int newNoLociAlloc;
		      TLoc*tmpLoc;

		      newNoLociAlloc = ( Hom->NoLociUsed / NoLociInit ) * 
			                   NoLociInit;
		      
		      if ( 0 == newNoLociAlloc )
			newNoLociAlloc += NoLociInit;
		      
		      if ( newNoLociAlloc < Hom->NoLociUsed )
			newNoLociAlloc += NoLociInit;
		      
		      tmpLoc = ( TLoc* ) realloc( Hom->Loc,
			( size_t ) newNoLociAlloc * sizeof( TLoc ) );

		      if ( NULL != tmpLoc )
			{
			  Hom->Loc = tmpLoc;
			  Hom->NoLociAlloc = newNoLociAlloc;
			}
		    }
		}
	  
	  if ( ( Pop->NoInds ) != ( Pop->NoIndsAlloc ) )
	    {
	      info( -1, "NoIndsAlloc is greater than NoInds\n" );
	      info( -1, "use resize to remove unused individuals\n" );
	      
	      for ( ind = Pop->NoInds;ind < Pop->NoIndsAlloc;ind++ )
		{
		  for ( chrom = 0;chrom < NoChroms;chrom++ )
		    {
		      for ( hom = 0;hom < NoHoms;hom++ )
			{
			  THom*Hom = &( ( *Pop ).Ind[ ind ].Chrom[ chrom ].Hom[ hom ] );
			  TLoc*tmpLoc;

			  tmpLoc = ( TLoc* ) realloc( Hom->Loc,
			    ( size_t ) NoLociInit * sizeof( TLoc ) );

			  if ( NULL != tmpLoc )
			    {
			      Hom->Loc = tmpLoc;
			      Hom->NoLociAlloc = NoLociInit;
			    }
			}
		    }
		}
	    }
	}
      
      PopName = strtok( NULL, "[' ','\n']" );
    }
}

void int_resize_population( char**PopName, POPULATION newSize )
{
  POPULATION ind;
  TPopListEl*PopListEl;
  TPop*Pop = NULL;
  TInd*tmpInd;
  
  find_list_element( *PopName, &PopListEl );

  if ( PopListEl == NULL ) {
    info( -2, "Population not defined: %s\n", *PopName );
    return;
  }

  if ( newSize < 0 ){
    info( -2, "Population size must be a positive number \n" );
    return;
  }
  
  Pop = PopListEl->Pop;

  if ( newSize > Pop->NoIndsAlloc )
    {
      tmpInd = ( TInd* ) realloc( Pop->Ind,
				    ( size_t ) newSize * sizeof( TInd ) );

      if ( NULL == tmpInd )
	{
	  info( -2, "Memory allocation failed in int_resize_population\n" );
	  return;
	}

      Pop->Ind = tmpInd;
      
      for ( ind = Pop->NoIndsAlloc; 
	    ind < newSize; 
	    ind++ )
	if ( 0 != new_individuum( &( Pop->Ind[ ind ] ) ) )
	  goto fail;
      
      Pop->NoIndsAlloc = newSize;
    }
  
  if ( newSize < Pop->NoIndsAlloc )
    {
      if ( Pop->NoInds > newSize )
	Pop->NoInds = newSize;

      for ( ind = newSize; 
	    ind < Pop->NoIndsAlloc;
	    ind++ )
	remove_individuum( &( Pop->Ind[ ind ] ) );
      
      if ( 0 == newSize )
	{
	  yFree( Pop->Ind );
	}
      else
	{
	  tmpInd = ( TInd* ) realloc( Pop->Ind,
				      ( size_t ) newSize * sizeof( TInd ) );

	  if ( NULL != tmpInd )
	    Pop->Ind = tmpInd;
	}

      Pop->NoIndsAlloc = newSize;
    }

  info( 1, "Population %s resized to %ld individuals\n",	*PopName, newSize );
  return;

fail:

  while ( ind > Pop->NoIndsAlloc )
    {
      ind--;
      remove_individuum( &( Pop->Ind[ ind ] ) );
    }

  if ( 0 == Pop->NoIndsAlloc )
    {
      yFree( Pop->Ind );
    }
  else
    {
      tmpInd = ( TInd* ) realloc( Pop->Ind,
				  ( size_t ) Pop->NoIndsAlloc *
				  sizeof( TInd ) );

      if ( NULL != tmpInd )
	Pop->Ind = tmpInd;
    }

  info( -2, "Memory allocation failed in int_resize_population\n" );
  return;
}

void resize_population ( char** PopName, char** cnewSize )
{
  POPULATION newSize = atol( *cnewSize );
  POPULATION oldSize;
  TPopListEl *PopListEl;

  find_list_element( *PopName, &PopListEl );

  if ( NULL == PopListEl )
    {
      int_resize_population( PopName, newSize );
      return;
    }

  oldSize = PopListEl->Pop->NoIndsAlloc;

  int_resize_population( PopName, newSize );

  if ( ( oldSize != newSize ) &&
       ( PopListEl->Pop->NoIndsAlloc == newSize ) )
    {
      st_clear_population_genotype( PopListEl->Pop );
      st_clear_population_gvalue( PopListEl->Pop );
    }
}

TPop* new_population ( char** Name, POPULATION* NoInds )
{
  POPULATION  ind = 0;
  TPop*       NewPop = NULL;
  TPopListEl* NewPopListEl;
  
  if ( 2 > strlen( *Name ) ) {
    info( -2, "Population name is too short\n" );
    return NULL;
  }

  if ( strlen( *Name ) >= MAXNAME ) {
    info( -2, "Population name is too long\n" );
    return NULL;
  }
  
  if ( 0 > *NoInds ) {
    info( -2, "Incorrect number of individuals\n" );
    return NULL;
  }
  
  find_list_element( *Name, &NewPopListEl );
  
  if ( NewPopListEl != NULL )  
    NewPop = NewPopListEl->Pop;
  
  if ( NULL != NewPop ) {
    info( -2, "Population already exists\n" );
    return NULL;
  }

  NewPop = (TPop*) yCalloc( 1, TPop );

  if ( NULL == NewPop )
    goto fail_alloc;
 
  NewPopListEl = (TPopListEl*) yCalloc( 1, TPopListEl );

if ( NULL == NewPopListEl )
    goto fail_alloc;
  
  strcpy( NewPop->Name, *Name );
  
  NewPop->NoIndsAlloc = *NoInds;
  NewPop->NoInds      = 0;
  NewPopListEl->Pop   = NewPop;
  NewPopListEl->Prev  = NULL;
  NewPopListEl->Next  = PopListStart;

  if ( NewPopListEl->Next != NULL ) 
    NewPopListEl->Next->Prev = NewPopListEl;

  PopListStart = NewPopListEl;

  if ( NewPop->NoIndsAlloc > 0 )
    {
      NewPop->Ind = (TInd*) yCalloc( NewPop->NoIndsAlloc, TInd );

      if ( NULL == NewPop->Ind )
        goto fail;

      for ( ind = 0; ind < NewPop->NoIndsAlloc; ind++ )
        {
          if ( 0 != new_individuum( &( NewPop->Ind[ ind ] ) ) )
            goto fail;
        }
    }
  else
    {
      NewPop->Ind = NULL;
    }

    info( 1, "Population %s with %ld individuals generated\n", *Name, *NoInds );
    return ( NewPop );

fail_alloc:

    yFree( NewPopListEl );
    yFree( NewPop );

    info( -2, "Memory allocation failed in new_population\n" );
    return NULL;
    
fail:

    while ( ind > 0 )
      {
        ind--;
        remove_individuum( &( NewPop->Ind[ ind ] ) );
      }

    yFree( NewPop->Ind );

    PopListStart = NewPopListEl->Next;

    if ( NULL != PopListStart )
      PopListStart->Prev = NULL;

    yFree( NewPop );
    yFree( NewPopListEl );

    info( -2, "Memory allocation failed in new_population\n" );
    return NULL;
}

void remove_population( char** PopNames )
{
  char*       token;
  TPopListEl* DelPopListEl;
  char*       _PopNames;
  
  _PopNames = ( char* ) yCalloc( strlen( *PopNames ) + 2, char );

  if ( NULL == _PopNames )
    {
      info( -2, "Memory allocation failed in remove_population\n" );
      return;
    }
  
  strcpy( _PopNames, *PopNames );
  token = strtok( _PopNames, "[' ','\n']" );
  
  while ( NULL != token )
    {
      find_list_element( token, &DelPopListEl );
      if ( DelPopListEl != NULL )
	{
	  POPULATION ind;
	  for ( ind = 0; 
		ind < DelPopListEl->Pop->NoIndsAlloc; 
		ind++ )
	    {
	      remove_individuum( &( DelPopListEl->Pop->Ind[ ind ] ) );
	    }
	  
	  if ( ( DelPopListEl->Prev == NULL ) && 
	       ( DelPopListEl->Next == NULL )     )
	    {
	      PopListStart = NULL;
	    }
	  else
            {
	      if ( DelPopListEl->Prev == NULL )
                {
		  DelPopListEl->Next->Prev = NULL;
		  PopListStart = DelPopListEl->Next;
                }
	      else
                {
		  if ( DelPopListEl->Next == NULL )
                    {
		      DelPopListEl->Prev->Next = NULL;
                    }
		  else
                    {
		      DelPopListEl->Next->Prev = DelPopListEl->Prev;
		      DelPopListEl->Prev->Next = DelPopListEl->Next;
                    }
                }
            }
	  
	  if ( NULL != DelPopListEl->Pop->Ind )
	    yFree( DelPopListEl->Pop->Ind );

	  yFree( DelPopListEl->Pop );
	  yFree( DelPopListEl );

	  info( 1, "Population %s deleted\n", token );
        }
      
      token = strtok( NULL, "[' ','\n']" );
    }
  yFree( _PopNames );
}


TPop *exist_population ( char*      PopName,
			 POPULATION NoInds,
			 int        empty )
{
  TPopListEl *PopListEl;

  find_list_element( PopName, &PopListEl );

  if ( PopListEl == NULL )
    {
      return new_population( &PopName, &NoInds );
    }

  if ( PopListEl->Pop->NoIndsAlloc != NoInds )
    {
      int_resize_population( &PopName, NoInds );

      if ( PopListEl->Pop->NoIndsAlloc != NoInds )
        return NULL;
    }

  if ( 1 == empty )
    {
      remove_genotype_population( &PopName );
      remove_evaluate_population( &PopName );
    }
  
  return PopListEl->Pop;
}

int insertLocus( unsigned short All, 
		 CPOS  Pos, 
		 THom* Hom )
{
  if ( ( 0 == Hom->NoLociUsed ) ||
       ( ( Hom->Loc[ Hom->NoLociUsed - 1 ].Pos < Pos ) &&
	 ( Hom->Loc[ Hom->NoLociUsed - 1 ].All != All ) ) )
    {
      TLoc * Loc;
      
      if ( Hom->NoLociUsed == Hom->NoLociAlloc )
        {
	  int newNoLociAlloc = Hom->NoLociAlloc + NoLociInit;
	  TLoc *tmpLoc;

	  tmpLoc = ( TLoc* ) realloc(
	    Hom->Loc,
	    ( size_t ) newNoLociAlloc * sizeof( TLoc ) );

	  if ( NULL == tmpLoc )
	    {
	      info( -2, "Memory allocation failed in insertLocus\n" );
	      return -1;
	    }

	  Hom->Loc = tmpLoc;
	  Hom->NoLociAlloc = newNoLociAlloc;
        }
      
      Loc = &( Hom->Loc[ Hom->NoLociUsed ] );
      Loc->Pos = Pos;
      Loc->All = All;
      Hom->NoLociUsed++;
      return 0;
    }

  return 1;
}

void init_population( char** Name,
		      char** cCount,
		      int*_Ind,
		      int*_Chrom,
		      int*_Hom,
		      double*_Pos,
		      int*_All        )
{
  int min_chrom = 0, 
    max_chrom = 0,
    min_hom   = 0, 
    max_hom    = 0;
  POPULATION row, 
    Count = atol( *cCount );
  POPULATION max_Ind = 0, 
    min_Ind = 0;
  TPopListEl* PopListEl;
  TPop* Pop;
  THom* Hom;

  if ( Count < 1 )
    {
      info( -2, "Population data contain no loci\n" );
      return;
    }
  
  find_list_element( *Name, &PopListEl );
  
  if ( PopListEl != NULL ) {
    info( -2, "First remove population %s\n", *Name );
    return;
  }

  for ( row = 0; row < Count; row++ )
    {
      max_chrom = ( *( _Chrom + row ) - 1 > max_chrom ) ? *( _Chrom + row ) - 1 : max_chrom;
      min_chrom = ( *( _Chrom + row ) - 1 < min_chrom ) ? *( _Chrom + row ) - 1 : min_chrom;
      max_hom = ( *( _Hom + row ) - 1 > max_hom ) ? *( _Hom + row ) - 1 : max_hom;
      min_hom = ( *( _Hom + row ) - 1 < min_hom ) ? *( _Hom + row ) - 1 : min_hom;
      max_Ind = ( *( _Ind + row ) - 1 > max_Ind ) ? *( _Ind + row ) - 1 : max_Ind;
      min_Ind = ( *( _Ind + row ) - 1 < min_Ind ) ? *( _Ind + row ) - 1 : min_Ind;
    }
  
  if ( ( min_chrom < 0 ) || ( max_chrom >= NoChroms ) ) {
    info( -2, "Error in the chromosome column\n" );
    return;
  }
 
  if ( ( min_hom < 0 ) || ( max_hom >= NoHoms ) ) {
    info( -2, "Error in the homologue column\n" );
    return;
  }
  
  if ( min_Ind < 0 ) {
    info( -2, "Error in the individual column\n" );
    return;
  }
  
  Pop = exist_population( *Name, max_Ind + 1, 1 );
  if ( NULL == Pop )  return;
  
  Pop->NoInds = max_Ind + 1;
  
  for ( row = 0; row < Count; row++ )
    {
      TInd*Ind = get_individual( Pop, *( _Ind + row ) - 1 );
      Hom = &( Ind->Chrom[*(_Chrom+row)-1].Hom[*(_Hom+row)-1] );

      if ( -1 == insertLocus( _All[ row ], _Pos[ row ], Hom ) )
        {
          remove_population( Name );
          return;
        }
      
      if ( 0 == _All[ row ] )
	info( -1, "Don't use 0 to code alleles\n" );
    }
}


void get_allele_number( char** Name, 
			int*Count  )
{
    int         chrom, hom;
    POPULATION  ind;
    TPopListEl* PopListEl;
    TPop*       Pop;
    TInd*       Ind;
    TChrom*     Chrom;
    THom*       Hom;

    find_list_element( *Name, &PopListEl );

    if ( PopListEl != NULL )
      {
        *Count = 0;
        Pop = PopListEl->Pop;
        for ( ind = 0; ind < Pop->NoInds; ind++ )
	  {
            Ind = get_individual( PopListEl->Pop, ind );
            for ( chrom = 0;chrom < NoChroms;chrom++ )
	      {
                Chrom = &( Ind->Chrom[ chrom ] );
                for ( hom = 0; hom < NoHoms; hom++ )
		  {
                    Hom = &( Chrom->Hom[ hom ] );
                    *Count += Hom->NoLociUsed;
		  }
	      }
	  }
      }
}


void get_population( 
		    char**  PopName,
		    int*    _Ind,
		    int*    _Chrom,
		    int*    _Hom,
		    double* _Pos,
		    int*    _All      )
{
  int chrom, hom, loc;
  POPULATION ind;
  POPULATION row = 0;
  TPopListEl*PopListEl;
  TPop*Pop;
  TInd*Ind;
  TChrom*Chrom;
  THom*Hom;
  TLoc*Loc;
  
  find_list_element( *PopName, &PopListEl );

  if ( NULL == PopListEl ) {
    info( -2, "Population not defined: %s\n", *PopName );
    return;
  }
  
  Pop = PopListEl->Pop;
  
  for ( ind = 0;ind < Pop->NoInds;ind++ )
    {
      Ind = get_individual( PopListEl->Pop, ind );
      for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
	  Chrom = &( Ind->Chrom[chrom] );
	  for ( hom = 0;hom < NoHoms;hom++ )
            {
	      Hom = &( Chrom->Hom[hom] );
	      for ( loc = 0; loc < Hom->NoLociUsed; loc++ )
                {
		  Loc = &( Hom->Loc[ loc ] );
		  _Ind[ row ] = ind + 1;
		  _Chrom[ row ] = chrom + 1;
		  _Hom[ row ] = hom + 1;
		  _Pos[ row ] = ( double ) Loc->Pos;
		  _All[ row++ ] = Loc->All;
                }
            }
        }
    }
}


void list_populations_parameter ( int* lPopName, int* NoPops )
{
  TPopListEl * PopListEl;
  PopListEl = PopListStart;

  while ( PopListEl != NULL )
    {
      (*NoPops) ++ ;
      (*lPopName) += strlen( PopListEl->Pop->Name ) + 1;
      PopListEl = PopListEl->Next;
    }
  ( *lPopName ) += 10;
}


void list_populations ( char **all_pops, int *inds )
{
  TPopListEl * PopListEl;
  PopListEl = PopListStart;
  
  strcpy( *all_pops, "" );
  
  while ( PopListEl != NULL )
    {
      *inds = PopListEl->Pop->NoInds;
      inds = inds + 1;
      strcat( *all_pops, PopListEl->Pop->Name );
      strcat( *all_pops, " " );
      PopListEl = PopListEl->Next;
    }
}


void set_co_freq( double*cofreq )
{
  if ( 0 <= *cofreq )
    {
      COFreq = *cofreq;
      info( 1, "Crossoverfrequency set to %f\n", *cofreq );
    }
  else
    info( -1, "Recombination Frequency must be a positive number\n" );
}


int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}


static int meiosis_impl( THom*   PgHom,
	      THom*   PHom, 
	      double  Length, 
	      double* significant,
	      int     r)
{
  unsigned int NuCO = gsl_ran_poisson(rng[r], COFreq * Length );
  CPOS posCO [NuCO + 2] ;
  unsigned int i, aHom;
  int aPos[] = {0, 0};
  CPOS pos;

  if ( -1 < *significant )
    {
      for ( i = 1;i <= NuCO;i++ )
        {
	  posCO[i] = fround( gsl_rng_uniform_pos(rng[r]) * Length, *significant );
	  posCO[i] = ( *(posCO+i) > Length ) ? Length : *(posCO+i);
        }
    }
  else
    {
      for ( i = 1;i < NuCO + 1;i++ )
        {
	  posCO[i] = gsl_rng_uniform_pos(rng[r]) * Length;
        }
    }

  posCO[ 0 ] = 0;
  posCO[ NuCO + 1 ] = Length;

  switch ( NuCO )
    {

    case 0:

    case 1:
      break;

    case 2:
      if ( ( posCO[ 1 ] > posCO[ 2 ] ) )
        {
	  pos = posCO[ 1 ];
	  posCO[ 1 ] = posCO[ 2 ];
	  posCO[ 2 ] = pos;
        }
      break;

    case 3:
      if ( ( posCO[ 1 ] < posCO[ 2 ] ) && ( posCO[ 2 ] < posCO[ 3 ] ) )
	break;

    /* fall through */
    default:
      // R_rsort( posCO + 1, NuCO );
      qsort ( posCO+1 , NuCO, sizeof (double), compare_doubles);

    }

  aHom = ( int ) gsl_rng_uniform_int( rng[r], 2 );

    for ( i = 0; i <= NuCO; i++ )
      {
        while ( ( aPos[ aHom ] + 1 < PHom[ aHom ].NoLociUsed ) &&
                ( PHom[ aHom ].Loc[ aPos[ aHom ] ].Pos < posCO[ i ] ) &&
                ( PHom[ aHom ].Loc[ aPos[ aHom ] + 1 ].Pos <= posCO[ i ] ) )
	  {
            aPos[ aHom ] ++ ;
	  }
        while ( ( aPos[ aHom ] < PHom[ aHom ].NoLociUsed ) &&
                ( PHom[ aHom ].Loc[ aPos[ aHom ] ].Pos < posCO[ i + 1 ] ) )
	  {
            pos = PHom[ aHom ].Loc[ aPos[ aHom ] ].Pos;
            pos = ( pos > posCO[ i ] ) ? pos : posCO[ i ];
            if ( -1 == insertLocus(
                   PHom[ aHom ].Loc[ aPos[ aHom ] ].All,
                   pos,
                   PgHom ) )
              return -1;

            aPos[ aHom ] ++;
	  }
        if ( aPos[ aHom ] > 0 )
	  aPos[ aHom ] -- ;
        aHom = ( aHom == 1 ) ? 0 : 1;
      }
  return 0;
}


void meiosis( THom*   PgHom,
              THom*   PHom,
              double  Length,
              double* significant,
              int     r)
{
  (void)meiosis_impl( PgHom, PHom, Length, significant, r );
}


/* qsort int comparison function */ 
int intcmp(const void *a, const void *b) 
{ 
    const int *ia = (const int *)a; 
    const int *ib = (const int *)b;
    return *ia  - *ib; 
} 

/* qsort unsigned short comparison function */ 
int unsignedshortcmp(const void *a, const void *b) 
{ 
    const unsigned short *ia = (const unsigned short *)a; 
    const unsigned short *ib = (const unsigned short *)b;
    return *ia  - *ib; 
} 


void cross ( 
	   char**  NamePg,  
	   char**  NameP1,
	   char**  NameP2,
	   double* Self,
	   char**  cNoPg,
	   double* significant,
	   int*    mode,
	   int*    crossing_scheme )
{
  POPULATION ind, 
             NoPg = atol( *cNoPg );

  TPopListEl *PopListEl;

  TPop *PgPop = NULL, 
       *P1Pop = NULL, 
       *P2Pop = NULL;
  int parallel_error = 0;


  if ( 2 != NoHoms ) {
    info( -2, "Implemented for diploids only\n" );
    return;
  }

  find_list_element( *NamePg, &PopListEl );
  if ( PopListEl != NULL )  PgPop = PopListEl->Pop;

  find_list_element( *NameP1, &PopListEl );
  if ( PopListEl != NULL )  P1Pop = PopListEl->Pop;

  find_list_element( *NameP2, &PopListEl );
  if ( PopListEl != NULL )  P2Pop = PopListEl->Pop;

  if ( -1 == NoPg )
    {
      if ( ( NULL == PgPop ) || ( PgPop->NoIndsAlloc < 1 ) )
        {
          info( -2,
                "Progeny population must exist when NoPg is -1\n" );
          return;
        }

      NoPg = PgPop->NoIndsAlloc;
    }
  else if ( 1 > NoPg )
    {
      info( -2, "Number of progenies must be positive\n" );
      return;
    }


  if ( ( P1Pop == NULL ) ||  ( P2Pop == NULL ) ) {
    info( -2, "Parent population not defined\n" );
    return;
  }
  
  if ( ( P1Pop == PgPop ) || ( P2Pop == PgPop ) ) {
    info( -2, "Progeny population incorrect\n" );
    return;
  }

  if ( ( 0 == P1Pop->NoInds ) || ( 0 == P2Pop->NoInds ) ){
    info( -2, "Parental population of size 0\n" );
    return;
  }

  PgPop = exist_population( *NamePg, NoPg, 1 );
  if ( NULL == PgPop )  return;

  PgPop->NoInds = NoPg;

  /*
   * Clear an error possibly remembered from an earlier parallel
   * operation.  info() suppresses R console output inside OpenMP
   * and records level -2 errors in this flag.
   */
  st_parallel_info_error_clear();

#ifdef MT
  #pragma omp parallel reduction(|:parallel_error)
#endif
  
  {
    int chrom;
    double r;
    POPULATION ind1, ind2; 
    TInd *P1Ind, *P2Ind, *PgInd;
    int local_error = 0;

    int i = omp_get_thread_num();

#ifdef MT
  #pragma omp for 
#endif
 
  for ( ind = 0; ind < PgPop->NoInds; ind++ )
  {

    if ( ( 0 == *crossing_scheme ) || ( 1 == *crossing_scheme ) ) {
      P1Ind = get_individual( P1Pop, 
			      ind1 = (POPULATION)gsl_rng_uniform_int ( rng[i], 
                                     (unsigned long int)P1Pop->NoInds ) );
    }
    else {
      P1Ind = get_individual( P1Pop, 
			      ind1 = ind % P1Pop->NoInds );
    }

    if ( ( 0 == *crossing_scheme ) || ( 2 == *crossing_scheme ) ) {
      P2Ind = get_individual( P2Pop,
			      ind2 = (POPULATION) gsl_rng_uniform_int ( rng[i], 
                                     (unsigned long int)P2Pop->NoInds ) );
    }
    else {
      P2Ind = get_individual( P2Pop, 
			      ind2 = ind % P2Pop->NoInds );
    }

    if ( ( 0 == *mode ) && 
	 ( P1Ind == P2Ind ) && 
	 ( 1 < P1Pop->NoInds ) )
    {
      while ( P1Ind == P2Ind )
	{
	  P2Ind = get_individual( P2Pop,
				  ind2 = (POPULATION)gsl_rng_uniform_int( rng[i],
	                          (unsigned long int)P2Pop->NoInds) );
	}
    }

    r = gsl_rng_uniform_pos(rng[i]);
    PgInd = get_individual( PgPop, ind );

    if ( *Self / 2 > r )
      {
	P2Ind = P1Ind;
	info( 2, "  %s (%ld) selfed -> %s (%ld)\n", *NameP1, ind1, *NamePg, ind );
      }
    else
      {
	if ( *Self > r )
	  {
	    P1Ind = P2Ind;
	    info( 2, "  %s (%ld) selfed -> %s (%ld)\n", *NameP2, ind2, *NamePg, ind );
	  }
	else
	  {
	    info( 2, "  %s (%ld) x %s (%ld) -> %s (%ld)\n",
		  *NameP1,
		  ind1,
		  *NameP2,
		  ind2,
		  *NamePg,
		  ind );
	  }
      }
    
    if ( 2 == NoHoms )
      {
	for ( chrom = 0;chrom < NoChroms;chrom++ )
	  {
	    PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed = 0;
	    PgInd->Chrom[ chrom ].Hom[ 1 ].NoLociUsed = 0;
	    
	    if ( -1 == meiosis_impl(
                   PgInd->Chrom[ chrom ].Hom,
                   ( *P1Ind ).Chrom[ chrom ].Hom,
                   PopDescr[ chrom ].ChromLen,
                   significant,
                   i ) )
              local_error = 1;

            if ( !local_error &&
                 ( -1 == meiosis_impl(
                   &( PgInd->Chrom[ chrom ].Hom[ 1 ] ),
                   ( *P2Ind ).Chrom[ chrom ].Hom,
                   PopDescr[ chrom ].ChromLen,
                   significant,
                   i ) ) )
              local_error = 1;

            if ( local_error )
              {
                parallel_error = 1;
                break;
              }
	  }
      }
  }
  } //parallel

  /*
   * Back on the serial R thread: report one generic error if a
   * worker encountered a level -2 error.  Routine worker output
   * is intentionally not retained or replayed.
   */
  if ( parallel_error || st_parallel_info_error_get() )
    {
      st_parallel_info_error_clear();
      remove_population( NamePg );
      info( -2, "Error during parallel crossing\n" );
      return;
    }

  st_parallel_info_error_clear();
}


void ps_set_num_threads (int * activate)
{
  if ( *activate > MAXTHREADS)  *activate = MAXTHREADS ;
  if ( *activate < 1 )          *activate = 1 ;
  omp_set_num_threads (*activate);
}

void ssd_mating( char**  NamePg,
		 char**  NameP1,
		 char**  cNoPg,
		 int*    maxcycles,
		 double* significant )
{
  POPULATION ind, NoPg = atol( *cNoPg ), nPg = atol( *cNoPg ), prog;
  int chrom;

  TPopListEl *PopListEl;
  TPop       *PgPop = NULL, 
             *P1Pop = NULL;
  TInd       *PgInd, 
             *P1Ind;
  THom       *tmpHomA = NULL, 
             *tmpHomB;
  int cycle, 
      nocycles = *maxcycles, 
      hom, 
      loc;
  int destination_prepared = 0;

  if ( 2 != NoHoms ){
    info( -2, "Implemented only for diploids\n" );
    return;
  }
  
  if ( 1 > NoPg ){
    info( -2, "Number of progenies must be a positive number\n" );
    return;
  }

  if ( 1 > nocycles ){
    info( -2, "Number of mating cycles must be positive\n" );
    return;
  }

  find_list_element( *NamePg, &PopListEl );

  if ( PopListEl != NULL )
    PgPop = PopListEl->Pop;

  find_list_element( *NameP1, &PopListEl );

  if ( PopListEl != NULL )
    P1Pop = PopListEl->Pop;

  if ( P1Pop == NULL )
    {
      info( -2, "Population not defined: %s\n", *NameP1 );
      return;
    }

  if ( P1Pop == PgPop )
    {
      info(
	   -2,
	   "Progeny population name (%s) is identical to parent population name (%s)\n",
	   *NamePg, *NameP1
	   );
      return;
    }

  if ( 0 == P1Pop->NoInds )
    {
      info( -2, "Parental population of size 0\n" );
      return;
    }

  if ( P1Pop->NoInds > LONG_MAX / NoPg )
    {
      info( -2, "Number of progenies is too large\n" );
      return;
    }

  NoPg = P1Pop->NoInds * NoPg;

  tmpHomA = ( THom* ) malloc(
    ( size_t ) NoHoms * sizeof( THom ) );

  if ( NULL == tmpHomA )
    {
      info( -2, "Memory allocation failed in ssd_mating\n" );
      return;
    }

  for ( hom = 0;hom < NoHoms;hom++ )
    tmpHomA[ hom ].Loc = NULL;

  for ( hom = 0;hom < NoHoms;hom++ )
    {
      tmpHomA[ hom ].NoLociAlloc = NoLociInit;
      tmpHomA[ hom ].NoLociUsed = 0;
      tmpHomA[ hom ].Loc = ( TLoc* ) malloc(
	( size_t ) NoLociInit * sizeof( TLoc ) );

      if ( NULL == tmpHomA[ hom ].Loc )
	{
	  info( -2, "Memory allocation failed in ssd_mating\n" );
	  goto cleanup;
	}
    }

  PgPop = exist_population( *NamePg, NoPg, 1 );
  if ( NULL == PgPop ) goto cleanup;

  destination_prepared = 1;
  PgPop->NoInds = NoPg;

  for ( ind = 0; ind < P1Pop->NoInds; ind++ )
    {
      P1Ind = get_individual( P1Pop, ind );
      for ( prog = 0;prog < nPg;prog++ )
	{
	  PgInd = get_individual( PgPop, ind * nPg + prog );
	  for ( chrom = 0;chrom < NoChroms;chrom++ )
	    {
	      int diff = 0;
	      cycle = nocycles;
	      tmpHomB = P1Ind->Chrom[ chrom ].Hom;

	      do
		{
		  diff = 0;

		  for ( hom = 0;hom < NoHoms;hom++ )
		    {
		      if ( tmpHomA[ hom ].NoLociAlloc <
			   tmpHomB[ hom ].NoLociUsed )
			{
			  int newNoLociAlloc =
			    tmpHomB[ hom ].NoLociUsed;
			  TLoc *tmpLoc;

			  tmpLoc = ( TLoc* ) realloc(
			    tmpHomA[ hom ].Loc,
			    ( size_t ) newNoLociAlloc *
			    sizeof( TLoc ) );

			  if ( NULL == tmpLoc )
			    {
			      info(
				   -2,
				   "Memory allocation failed in ssd_mating\n"
				   );
			      goto fail_population;
			    }

			  tmpHomA[ hom ].Loc = tmpLoc;
			  tmpHomA[ hom ].NoLociAlloc =
			    newNoLociAlloc;
			}

		      tmpHomA[ hom ].NoLociUsed =
			tmpHomB[ hom ].NoLociUsed;

		      for ( loc = 0;
			    loc < tmpHomA[ hom ].NoLociUsed;
			    loc++ )
			{
			  tmpHomA[ hom ].Loc[ loc ].All =
			    tmpHomB[ hom ].Loc[ loc ].All;
			  tmpHomA[ hom ].Loc[ loc ].Pos =
			    tmpHomB[ hom ].Loc[ loc ].Pos;
			}
		    }

		  PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed = 0;
		  PgInd->Chrom[ chrom ].Hom[ 1 ].NoLociUsed = 0;

		  if ( -1 == meiosis_impl(
                         PgInd->Chrom[ chrom ].Hom,
                         tmpHomA,
                         PopDescr[ chrom ].ChromLen,
                         significant,
                         0 ) )
                    goto fail_population;

		  if ( -1 == meiosis_impl(
                         &( PgInd->Chrom[ chrom ].Hom[ 1 ] ),
                         tmpHomA,
                         PopDescr[ chrom ].ChromLen,
                         significant,
                         0 ) )
                    goto fail_population;

		  if ( PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed !=
		       PgInd->Chrom[ chrom ].Hom[ 1 ].NoLociUsed )
		    {
		      diff++;
		    }
		  else
		    {
		      for ( loc = 0;
			    loc <
			      PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed;
			    loc++ )
			{
			  if (
			      ( PgInd->Chrom[ chrom ].Hom[ 0 ].
				Loc[ loc ].All !=
				PgInd->Chrom[ chrom ].Hom[ 1 ].
				Loc[ loc ].All ) ||
			      ( PgInd->Chrom[ chrom ].Hom[ 0 ].
				Loc[ loc ].Pos !=
				PgInd->Chrom[ chrom ].Hom[ 1 ].
				Loc[ loc ].Pos )
			      )
			    {
			      diff++;
			      break;
			    }
			}
		    }

		  tmpHomB = PgInd->Chrom[ chrom ].Hom;
		}
	      while ( ( 0 != --cycle ) && ( diff ) );
	    }
	}
    }

  destination_prepared = 0;
  goto cleanup;

 fail_population:

  if ( destination_prepared )
    {
      remove_population( NamePg );
      PgPop = NULL;
      destination_prepared = 0;
    }

 cleanup:

  if ( NULL != tmpHomA )
    {
      for ( hom = 0;hom < NoHoms;hom++ )
	free( tmpHomA[ hom ].Loc );

      free( tmpHomA );
    }
}


void dh( char**  NamePg,
	 char**  NameP1,
	 char**  cNoPg,
	 double* significant )
{
    POPULATION ind, 
               NoPg = atol( *cNoPg );

    int chrom;

    TPopListEl*PopListEl;

    TPop *PgPop = NULL, 
         *P1Pop = NULL;

    TInd *PgInd, 
         *P1Ind;
    int destination_prepared = 0;

    if ( 2 != NoHoms ) {
      info( -2, "Implemented for diploids only\n" );
      return;
    }

    find_list_element( *NamePg, &PopListEl );
    if ( PopListEl != NULL ) PgPop = PopListEl->Pop;

    find_list_element( *NameP1, &PopListEl );
    if ( PopListEl != NULL ) P1Pop = PopListEl->Pop;

    if ( -1 == NoPg )
      {
        if ( ( NULL == PgPop ) || ( PgPop->NoIndsAlloc < 1 ) )
          {
            info( -2,
                  "Progeny population must exist when NoPg is -1\n" );
            return;
          }

        NoPg = PgPop->NoIndsAlloc;
      }
    else if ( 1 > NoPg )
      {
        info( -2, "Number of progenies must be positive\n" );
        return;
      }


    if ( P1Pop == NULL )  {
      info( -2, "Parent population not defined\n" );
      return;
    }

    if ( P1Pop == PgPop ) {
      info( -2, "Progeny population incorrect\n" );
      return;
    }

    if  ( 0 == P1Pop->NoInds ) {
      info( -2, "Parental population of size 0\n" );
      return;
    }

    PgPop = exist_population( *NamePg, NoPg, 1 );
    if ( NULL == PgPop )  return;

    destination_prepared = 1;
    PgPop->NoInds = NoPg;

    for ( ind = 0;ind < PgPop->NoInds;ind++ )
    {
        P1Ind = get_individual( P1Pop,
           (POPULATION)gsl_rng_uniform_int(rng[0], (unsigned long int)P1Pop->NoInds ) );
        PgInd = get_individual( PgPop, ind );

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            int loc, hom;

            PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed = 0;
            PgInd->Chrom[ chrom ].Hom[ 1 ].NoLociUsed = 0;

            if ( -1 == meiosis_impl(
                   PgInd->Chrom[ chrom ].Hom,
                   P1Ind->Chrom[ chrom ].Hom,
                   PopDescr[ chrom ].ChromLen,
                   significant,
                   0 ) )
              goto fail_population;

            for ( hom = 1;hom < NoHoms;hom++ )
            {
                THom*P1Hom = &( PgInd->Chrom[ chrom ] ).Hom[ hom ];
                THom*P2Hom = &( PgInd->Chrom[ chrom ] ).Hom[ 0 ];
                PgInd->Chrom[ chrom ].Hom[ hom ].NoLociUsed = 0;

                for ( loc = 0;loc < P2Hom->NoLociUsed;loc++ )
                {
                    if ( -1 == insertLocus(
                           P2Hom->Loc[ loc ].All,
                           P2Hom->Loc[ loc ].Pos,
                           P1Hom ) )
                      goto fail_population;
                }
            }
        }
    }

    return;

fail_population:

    if ( destination_prepared )
      remove_population( NamePg );

    return;
}




static int parse_int_token( const char *token,
			    int        *value )
{
  char *end;
  long x;

  if ( ( NULL == token ) || ( '\0' == *token ) )
    return 0;

  errno = 0;
  x = strtol( token, &end, 10 );

  if ( ( ERANGE == errno ) ||
       ( end == token ) ||
       ( '\0' != *end ) ||
       ( x < INT_MIN ) ||
       ( x > INT_MAX ) )
    return 0;

  *value = ( int ) x;

  return 1;
}

static int parse_double_token( const char *token,
			       CPOS       *value )
{
  char *end;
  double x;

  if ( ( NULL == token ) || ( '\0' == *token ) )
    return 0;

  errno = 0;
  x = strtod( token, &end );

  if ( ( ERANGE == errno ) ||
       ( end == token ) ||
       ( '\0' != *end ) ||
       ( !isfinite( x ) ) )
    return 0;

  *value = ( CPOS ) x;

  return 1;
}

static int parse_chrom_pos_token( const char *token,
				  int        *chrom,
				  CPOS       *pos )
{
  char *endChrom;
  char *endPos;
  long c;
  double p;

  if ( ( NULL == token ) || ( '\0' == *token ) )
    return 0;

  errno = 0;
  c = strtol( token, &endChrom, 10 );

  if ( ( ERANGE == errno ) ||
       ( endChrom == token ) ||
       ( '/' != *endChrom ) ||
       ( c < INT_MIN ) ||
       ( c > INT_MAX ) )
    return 0;

  errno = 0;
  p = strtod( endChrom + 1, &endPos );

  if ( ( ERANGE == errno ) ||
       ( endPos == endChrom + 1 ) ||
       ( '\0' != *endPos ) ||
       ( !isfinite( p ) ) )
    return 0;

  *chrom = ( int ) c;
  *pos = ( CPOS ) p;

  return 1;
}

static void free_mappoints( TMappoint *points,
			    int         n )
{
  int i;

  if ( NULL != points )
    {
      for ( i = 0;i < n;i++ )
	{
	  free( points[ i ].name );
	  free( points[ i ].class );
	}

      free( points );
    }
}

static int append_mappoint( TMappoint **points,
			    int        *n,
			    int         chrom,
			    CPOS        pos,
			    const char *name,
			    const char *class )
{
  TMappoint *tmp;
  TMappoint *point;

  if ( *n == INT_MAX )
    return 0;

  if ( ( size_t )( *n + 1 ) >
       SIZE_MAX / sizeof( TMappoint ) )
    return 0;

  tmp = ( TMappoint* ) realloc(
    *points,
    ( size_t )( *n + 1 ) * sizeof( TMappoint ) );

  if ( NULL == tmp )
    return 0;

  *points = tmp;

  point = &( tmp[ *n ] );

  point->chrom = chrom;
  point->Pos = pos;
  point->counter = *n;
  point->name = NULL;
  point->class = NULL;

  point->name = ( char* ) malloc( strlen( name ) + 1 );

  if ( NULL == point->name )
    return 0;

  point->class = ( char* ) malloc( strlen( class ) + 1 );

  if ( NULL == point->class )
    {
      free( point->name );
      point->name = NULL;
      return 0;
    }

  strcpy( point->name, name );
  strcpy( point->class, class );

  ( *n )++;

  return 1;
}

static int mnamecmp( const void *i,
		     const void *j )
{
  const TMappoint *m0 = ( const TMappoint* ) i;
  const TMappoint *m1 = ( const TMappoint* ) j;

  return strcmp( m0->name, m1->name );
}


int mccmp( const void *i,
	   const void *j )
{
  const TMappoint *m0 = ( const TMappoint* ) i;
  const TMappoint *m1 = ( const TMappoint* ) j;

  if ( m0->chrom > m1->chrom ) return 1;
  if ( m0->chrom < m1->chrom ) return -1;

  if ( m0->Pos > m1->Pos ) return 1;
  if ( m0->Pos < m1->Pos ) return -1;

  if ( m0->counter > m1->counter ) return 1;
  if ( m0->counter < m1->counter ) return -1;

  return 0;
}



static int validate_map_genome( void )
{
  int chrom;

  if ( ( NULL == PopDescr ) ||
       ( NoChroms < 1 ) ||
       ( NoHoms < 1 ) )
    {
      info( -2, "No valid genome defined\n" );
      return 0;
    }

  for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
      if ( !isfinite( PopDescr[ chrom ].ChromLen ) ||
           ( PopDescr[ chrom ].ChromLen <= 0 ) )
        {
          info( -2, "Invalid chromosome length\n" );
          return 0;
        }
    }

  return 1;
}


static int disperse_mappoints( TMappoint *points,
                               int         n,
                               double      factor,
                               int        *modified )
{
  int start = 0;
  int i;

  *modified = 0;

  if ( !isfinite( factor ) ||
       ( factor < 1 ) )
    {
      info(
        -2,
        "disperse.factor must be finite and >= 1\n" );

      return 0;
    }

  while ( start < n )
    {
      int end = start;
      int size;
      int chrom;
      int j;
      int backward = 0;
      CPOS pos;
      CPOS distance;
      CPOS delta;

      while (
        ( end + 1 < n ) &&
        ( points[ end + 1 ].chrom ==
          points[ start ].chrom ) &&
        ( points[ end + 1 ].Pos ==
          points[ start ].Pos ) )
        end++;

      size = end - start + 1;

      if ( size > 1 )
        {
          chrom = points[ start ].chrom;
          pos = points[ start ].Pos;

          if ( ( end + 1 < n ) &&
               ( points[ end + 1 ].chrom == chrom ) )
            {
              distance =
                points[ end + 1 ].Pos - pos;
            }

          else if (
            pos == PopDescr[ chrom - 1 ].ChromLen )
            {
              backward = 1;

              if ( ( start > 0 ) &&
                   ( points[ start - 1 ].chrom == chrom ) )
                {
                  distance =
                    pos - points[ start - 1 ].Pos;
                }

              else
                {
                  distance = pos;
                }
            }

          else
            {
              distance =
                PopDescr[ chrom - 1 ].ChromLen - pos;
            }

          delta =
            distance /
            ( factor * ( double ) size );

          if ( !( delta > 0 ) ||
               !isfinite( delta ) )
            {
              info(
                -2,
                "Unable to disperse loci at chromosome "
                "%d position %f\n",
                chrom,
                pos );

              return 0;
            }

          if ( backward )
            {
              for ( j = 0;j < size;j++ )
                {
                  points[ start + j ].Pos -=
                    ( size - 1 - j ) * delta;
                }
            }

          else
            {
              for ( j = 0;j < size;j++ )
                {
                  points[ start + j ].Pos +=
                    j * delta;
                }
            }

          *modified += size - 1;
        }

      start = end + 1;
    }

  /*
   * A dispersed map must be strictly increasing within
   * every chromosome. This also catches cases where the
   * requested delta is lost through floating-point rounding.
   */
  for ( i = 0;i < n;i++ )
    {
      int chrom = points[ i ].chrom;

      if ( !isfinite( points[ i ].Pos ) ||
           ( points[ i ].Pos < 0 ) ||
           ( points[ i ].Pos >
             PopDescr[ chrom - 1 ].ChromLen ) )
        {
          info(
            -2,
            "Invalid dispersed locus position on chromosome %d\n",
            chrom );

          return 0;
        }

      if ( ( i > 0 ) &&
           ( points[ i - 1 ].chrom == chrom ) &&
           !( points[ i ].Pos > points[ i - 1 ].Pos ) )
        {
          info(
            -2,
            "Unable to obtain distinct locus positions on "
            "chromosome %d\n",
            chrom );

          return 0;
        }
    }

  return 1;
}


static int prepare_mappoints( TMappoint *points,
                              int         n,
                              int         disperse,
                              double      factor,
                              int        *modified )
{
  int i;

  *modified = 0;

  if ( !validate_map_genome() )
    return 0;

  if ( ( NULL == points ) ||
       ( n < 1 ) )
    {
      info( -2, "Linkage map contains no loci\n" );
      return 0;
    }

  if ( n > INT_MAX / NoHoms )
    {
      info( -2, "Linkage map contains too many loci\n" );
      return 0;
    }

  if ( ( 0 != disperse ) &&
       ( 1 != disperse ) )
    {
      info( -2, "Invalid disperse argument\n" );
      return 0;
    }

  for ( i = 0;i < n;i++ )
    {
      int chrom = points[ i ].chrom;

      if ( ( chrom < 1 ) ||
           ( chrom > NoChroms ) )
        {
          info(
            -2,
            "Wrong chromosome (%d)\n",
            chrom );

          return 0;
        }

      if ( !isfinite( points[ i ].Pos ) ||
           ( points[ i ].Pos < 0 ) ||
           ( points[ i ].Pos >
             PopDescr[ chrom - 1 ].ChromLen ) )
        {
          info(
            -2,
            "Chromosome %d: Locusposition outside of the chromosome "
            "(0 <= %f <= %f)\n",
            chrom,
            points[ i ].Pos,
            PopDescr[ chrom - 1 ].ChromLen );

          return 0;
        }

      if ( ( NULL == points[ i ].name ) ||
           ( NULL == points[ i ].class ) ||
           ( '\0' == points[ i ].name[ 0 ] ) ||
           ( '\0' == points[ i ].class[ 0 ] ) )
        {
          info( -2, "Empty locus name or class\n" );
          return 0;
        }

      if ( ( strlen( points[ i ].name ) >= MAXNAME ) ||
           ( strlen( points[ i ].class ) >= MAXNAME ) )
        {
          info( -2, "Locus name or class is too long\n" );
          return 0;
        }

      if ( ( NULL !=
             strpbrk( points[ i ].name, " \t\r\n" ) ) ||
           ( NULL !=
             strpbrk( points[ i ].class, " \t\r\n" ) ) )
        {
          info(
            -2,
            "Whitespace not allowed in locus names or classes\n" );

          return 0;
        }
    }

  qsort(
    points,
    n,
    sizeof( TMappoint ),
    mnamecmp );

  for ( i = 1;i < n;i++ )
    {
      if ( 0 == strcmp(
             points[ i - 1 ].name,
             points[ i ].name ) )
        {
          info(
            -2,
            "Duplicate locus name: %s\n",
            points[ i ].name );

          return 0;
        }
    }

  qsort(
    points,
    n,
    sizeof( TMappoint ),
    mccmp );

  if ( disperse )
    {
      if ( !disperse_mappoints(
             points,
             n,
             factor,
             modified ) )
        return 0;
    }

  return 1;
}


static int read_mappoint_file( char       *fName,
                               TMappoint **points,
                               int        *nMP )
{
  FILE *fp = NULL;
  TMappoint *newMPs = NULL;
  char line[ 4 * MAXNAME + 128 ];
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

  if ( ( NULL == fName ) ||
       ( '\0' == *fName ) )
    {
      info( -2, "Invalid linkage map file\n" );
      return 0;
    }

  fp = fileopen( fName, "r" );

  if ( NULL == fp )
    goto fail;

  while ( NULL != fgets( line, sizeof( line ), fp ) )
    {
      int ch;

      lineNo++;

      if ( NULL == strchr( line, '\n' ) )
        {
          ch = fgetc( fp );

          if ( EOF != ch )
            {
              while ( ( '\n' != ch ) &&
                      ( EOF != ch ) )
                ch = fgetc( fp );

              info(
                -2,
                "Map file line %d is too long\n",
                lineNo );

              goto fail;
            }
        }

      cChrom = strtok( line, " \t\r\n" );
      cPos = strtok( NULL, " \t\r\n" );
      lociname = strtok( NULL, " \t\r\n" );
      lociclass = strtok( NULL, " \t\r\n" );
      extra = strtok( NULL, " \t\r\n" );

      if ( ( NULL == cChrom ) ||
           ( NULL == cPos ) ||
           ( NULL == lociname ) ||
           ( NULL == lociclass ) ||
           ( NULL != extra ) )
        {
          info(
            -2,
            "Error in map file line %d\n",
            lineNo );

          goto fail;
        }

      if ( !parse_int_token(
             cChrom,
             &chrom ) )
        {
          info(
            -2,
            "Invalid chromosome in map file line %d\n",
            lineNo );

          goto fail;
        }

      if ( !parse_double_token(
             cPos,
             &pos ) )
        {
          info(
            -2,
            "Invalid locus position in map file line %d\n",
            lineNo );

          goto fail;
        }

      if ( !append_mappoint(
             &newMPs,
             nMP,
             chrom,
             pos,
             lociname,
             lociclass ) )
        goto fail_alloc;
    }

  if ( ferror( fp ) )
    {
      info(
        -2,
        "Error reading map file %s\n",
        fName );

      goto fail;
    }

  if ( 0 != fclose( fp ) )
    {
      fp = NULL;

      info(
        -2,
        "Can't close the file stream\n" );

      goto fail;
    }

  fp = NULL;

  if ( 0 == *nMP )
    {
      info( -2, "Linkage map contains no loci\n" );
      goto fail;
    }

  *points = newMPs;

  return 1;

fail_alloc:

  info(
    -2,
    "Memory allocation failed while reading linkage map\n" );

fail:

  if ( NULL != fp )
    fileclose( fp );

  free_mappoints( newMPs, *nMP );

  *nMP = 0;

  return 0;
}


static int write_mappoint_file( char            *fName,
                                const TMappoint *points,
                                int              n )
{
  FILE *fp;
  int i;

  fp = fileopen( fName, "w" );

  if ( NULL == fp )
    return 0;

  for ( i = 0;i < n;i++ )
    {
      if ( 0 > fprintf(
             fp,
             "%d\t%.17g\t%s\t%s\n",
             points[ i ].chrom,
             points[ i ].Pos,
             points[ i ].name,
             points[ i ].class ) )
        {
          info(
            -2,
            "Error writing map file %s\n",
            fName );

          goto fail;
        }
    }

  if ( 0 != fclose( fp ) )
    {
      info(
        -2,
        "Can't close the file stream\n" );

      remove( fName );

      return 0;
    }

  return 1;

fail:

  fclose( fp );
  remove( fName );

  return 0;
}


static int install_mappoints( TMappoint *points,
			      int         n )
{
  TMappoint **newMap = NULL;
  int *count = NULL;
  int chrom, i, j;
  int counter = 0;

  if ( ( NULL == PopDescr ) ||
       ( NoChroms < 1 ) ||
       ( NoHoms < 1 ) ||
       ( n < 1 ) )
    return 0;

  if ( n > INT_MAX / NoHoms )
    return 0;

  if ( ( size_t ) NoChroms >
       SIZE_MAX / sizeof( TMappoint* ) )
    return 0;

  newMap = ( TMappoint** ) malloc(
    ( size_t ) NoChroms * sizeof( TMappoint* ) );

  if ( NULL == newMap )
    return 0;

  for ( chrom = 0;chrom < NoChroms;chrom++ )
    newMap[ chrom ] = NULL;

  if ( ( size_t ) NoChroms >
       SIZE_MAX / sizeof( int ) )
    goto fail;

  count = ( int* ) malloc(
    ( size_t ) NoChroms * sizeof( int ) );

  if ( NULL == count )
    goto fail;

  for ( chrom = 0;chrom < NoChroms;chrom++ )
    count[ chrom ] = 0;

  for ( i = 0;i < n;i++ )
    {
      if ( ( points[ i ].chrom < 1 ) ||
	   ( points[ i ].chrom > NoChroms ) )
	goto fail;

      count[ points[ i ].chrom - 1 ]++;
    }

  for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
      if ( 0 < count[ chrom ] )
	{
	  if ( ( size_t ) count[ chrom ] >
	       SIZE_MAX / sizeof( TMappoint ) )
	    goto fail;

	  newMap[ chrom ] = ( TMappoint* ) malloc(
	    ( size_t ) count[ chrom ] *
	    sizeof( TMappoint ) );

	  if ( NULL == newMap[ chrom ] )
	    goto fail;
	}
    }

  i = 0;

  for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
      for ( j = 0;j < count[ chrom ];j++ )
	{
	  if ( ( i >= n ) ||
	       ( points[ i ].chrom != chrom + 1 ) )
	    goto fail;

	  newMap[ chrom ][ j ] = points[ i ];
	  newMap[ chrom ][ j ].counter = counter++;

	  i++;
	}
    }

  if ( i != n )
    goto fail;

  /*
   * Ownership of name/class now passes from points
   * to newMap. There can be no failing operation below.
   */
  for ( i = 0;i < n;i++ )
    {
      points[ i ].name = NULL;
      points[ i ].class = NULL;
    }

  remove_map();

  for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
      PopDescr[ chrom ].Mappoint =
	newMap[ chrom ];

      PopDescr[ chrom ].NrMappoints =
	count[ chrom ];
    }

  cMP = n;
  cMPhom = NoHoms * n;

  free( newMap );
  free( count );

  return 1;

fail:

  if ( NULL != newMap )
    {
      for ( chrom = 0;chrom < NoChroms;chrom++ )
	free( newMap[ chrom ] );
    }

  free( newMap );
  free( count );

  return 0;
}

void linkage_map_load( char   **fName,
                       int     *disperse,
                       char   **fNameDisperse,
                       double  *factor,
                       int     *retval )
{
  TMappoint *points = NULL;
  int nMP = 0;
  int modified = 0;

  if ( NULL == retval )
    {
      info( -2, "Invalid linkage map arguments\n" );
      return;
    }

  *retval = -2;

  if ( ( NULL == fName ) ||
       ( NULL == *fName ) ||
       ( NULL == disperse ) ||
       ( NULL == factor ) )
    {
      info( -2, "Invalid linkage map arguments\n" );
      return;
    }

  if ( !validate_map_genome() )
    return;

  if ( !read_mappoint_file(
         *fName,
         &points,
         &nMP ) )
    return;

  if ( !prepare_mappoints(
         points,
         nMP,
         *disperse,
         *factor,
         &modified ) )
    goto fail;

  /*
   * file.disperse is retained only for compatibility.
   * No additional file is written for the normal empty name.
   */
  if ( *disperse &&
       ( NULL != fNameDisperse ) &&
       ( NULL != *fNameDisperse ) &&
       ( '\0' != ( *fNameDisperse )[ 0 ] ) )
    {
      if ( !write_mappoint_file(
             *fNameDisperse,
             points,
             nMP ) )
        goto fail;
    }

  if ( !install_mappoints(
         points,
         nMP ) )
    {
      info(
        -2,
        "Memory allocation failed in linkage_map_load\n" );

      goto fail;
    }

  free_mappoints( points, nMP );

  if ( *disperse )
    {
      if (0<modified) info(
        0,
        "%d loci positions modified\n",
        modified );
    }

  if ( *disperse &&
       ( NULL != fNameDisperse ) &&
       ( NULL != *fNameDisperse ) &&
       ( '\0' != ( *fNameDisperse )[ 0 ] ) )
    {
      info(
        0,
        "Dispersed linkage map saved to file: %s\n",
        *fNameDisperse );
    }

  info(
    0,
    "Map file %s loaded\n",
    *fName );

  info(
    0,
    "%d loci defined\n",
    nMP );

  *retval = 0;

  return;

fail:

  free_mappoints( points, nMP );

  return;
}


void define_map( char **fName,
                 int  *retval )
{
  int disperse = 1;
  double factor = MAP_DISPERSE_FACTOR_DEFAULT;

  linkage_map_load(
    fName,
    &disperse,
    NULL,
    &factor,
    retval );
}


void define_linkage_map(
    int    *dim,
    int    *chrom,
    double *pos,
    char  **names,
    char  **class
)
{
  TMappoint *newMPs = NULL;
  int nMP = 0;
  int modified = 0;
  int i;

  if ( !validate_map_genome() )
    return;

  if ( ( NULL == dim ) ||
       ( *dim < 1 ) )
    {
      info( -2, "Linkage map contains no loci\n" );
      return;
    }

  if ( ( NULL == chrom ) ||
       ( NULL == pos ) ||
       ( NULL == names ) ||
       ( NULL == class ) )
    {
      info( -2, "Invalid linkage map data\n" );
      return;
    }

  if ( *dim > INT_MAX / NoHoms )
    {
      info( -2, "Linkage map contains too many loci\n" );
      return;
    }

  for ( i = 0;i < *dim;i++ )
    {
      if ( ( NULL == names[ i ] ) ||
           ( NULL == class[ i ] ) )
        {
          info( -2, "Invalid linkage map data\n" );
          goto fail;
        }

      if ( !append_mappoint(
             &newMPs,
             &nMP,
             chrom[ i ],
             pos[ i ],
             names[ i ],
             class[ i ] ) )
        goto fail_alloc;
    }

  if ( !prepare_mappoints(
         newMPs,
         nMP,
         1,
         MAP_DISPERSE_FACTOR_DEFAULT,
         &modified ) )
    goto fail;

  if ( !install_mappoints(
         newMPs,
         nMP ) )
    goto fail_alloc;

  free_mappoints( newMPs, nMP );

  if (0<modified)info(
    0,
    "%d loci positions modified\n",
    modified );

  info(
    0,
    "%d loci defined\n",
    nMP );

  return;

fail_alloc:

  info(
    -2,
    "Memory allocation failed in define_linkage_map\n" );

fail:

  free_mappoints( newMPs, nMP );

  return;
}


void generate_map_file( char **fName,
			char **Description,
			int  *retval )
{
  TMappoint *newMPs = NULL;
  char *description = NULL;
  char *token;
  char name[ MAXNAME ];
  char classname[ MAXNAME ];
  int status = 0;
  int number = 1;
  int cMPs = 0;
  int maxMP;
  int i;
  FILE *fp = NULL;
  int outputOpened = 0;

  *retval = -2;

  if ( ( NULL == PopDescr ) ||
       ( NoChroms < 1 ) ||
       ( NoHoms < 1 ) )
    {
      info( -2, "No valid genome defined\n" );
      return;
    }

  for ( i = 0;i < NoChroms;i++ )
    {
      if ( !isfinite( PopDescr[ i ].ChromLen ) ||
	   ( PopDescr[ i ].ChromLen <= 0 ) )
	{
	  info( -2, "Invalid chromosome length\n" );
	  return;
	}
    }

  maxMP = INT_MAX / NoHoms;

  description = ( char* ) malloc(
    strlen( *Description ) + 1 );

  if ( NULL == description )
    goto fail_alloc;

  strcpy( description, *Description );

  token = strtok(
    description,
    "[]', \t\r\n" );

  while ( NULL != token )
    {
      switch ( status )
	{
	case 0:

	  if ( '*' == token[ 0 ] )
	    {
	      if ( !parse_int_token(
		     token + 1,
		     &number ) ||
		   ( number < 1 ) )
		{
		  info(
		    -2,
		    "Invalid number of loci: %s\n",
		    token );

		  goto fail;
		}

	      status = 1;
	    }

	  else
	    {
	      if ( strlen( token ) >= MAXNAME )
		{
		  info(
		    -2,
		    "Locus name too long: %s\n",
		    token );

		  goto fail;
		}

	      strcpy( name, token );
	      number = 1;
	      status = 2;
	    }

	  break;

	case 1:

	  if ( strlen( token ) >= MAXNAME )
	    {
	      info(
		-2,
		"Locus name too long: %s\n",
		token );

	      goto fail;
	    }

	  strcpy( name, token );
	  status = 2;
	  break;

	case 2:

	  if ( strlen( token ) >= MAXNAME )
	    {
	      info(
		-2,
		"Locus class too long: %s\n",
		token );

	      goto fail;
	    }

	  strcpy( classname, token );
	  status = 3;
	  break;

	case 3:

	  if ( number > maxMP - cMPs )
	    {
	      info(
		-2,
		"Linkage map contains too many loci\n" );

	      goto fail;
	    }

	  if ( 0 == strcmp( token, "random" ) )
	    {
	      double genomelength = 0;
	      int lNR = 0;

	      for ( i = 0;i < NoChroms;i++ )
		genomelength +=
		  PopDescr[ i ].ChromLen;

	      if ( !isfinite( genomelength ) ||
		   ( genomelength <= 0 ) )
		{
		  info( -2, "Invalid genome length\n" );
		  goto fail;
		}

	      for ( i = 0;i < number;i++ )
		{
		  char tmpstring[ MAXNAME ];
		  int nchars;
		  int chrom = -1;
		  CPOS pos =
		    - gsl_rng_uniform_pos( rng[ 0 ] ) *
		    genomelength;

		  while ( ( pos < 0 ) &&
			  ( chrom + 1 < NoChroms ) )
		    pos +=
		      PopDescr[ ++chrom ].ChromLen;

		  if ( ( chrom < 0 ) ||
		       ( chrom >= NoChroms ) ||
		       ( pos < 0 ) ||
		       ( pos >
			 PopDescr[ chrom ].ChromLen ) )
		    {
		      info(
			-2,
			"Error generating random map position\n" );

		      goto fail;
		    }

		  nchars = snprintf(
		    tmpstring,
		    sizeof( tmpstring ),
		    "%s%d",
		    name,
		    lNR++ );

		  if ( ( nchars < 0 ) ||
		       ( nchars >= MAXNAME ) )
		    {
		      info(
			-2,
			"Generated locus name is too long\n" );

		      goto fail;
		    }

		  if ( !append_mappoint(
			 &newMPs,
			 &cMPs,
			 chrom + 1,
			 pos,
			 tmpstring,
			 classname ) )
		    goto fail_alloc;
		}

	      number = 1;
	      status = 0;
	    }

	  else if ( 0 ==
		    strcmp( token, "terminal" ) )
	    {
	      int chrom, j;
	      int lNR = 0;

	      if ( 0 != number % NoChroms )
		{
		  info(
		    -2,
		    "Number of terminal loci must be a multiple "
		    "of the number of chromosomes\n" );

		  goto fail;
		}

	      for ( chrom = 0;
		    chrom < NoChroms;
		    chrom++ )
		{
		  int nPerChrom =
		    number / NoChroms;

		  CPOS delta =
		    ( nPerChrom > 1 ) ?
		    PopDescr[ chrom ].ChromLen /
		    ( nPerChrom - 1 ) : 0;

		  for ( j = 0;j < nPerChrom;j++ )
		    {
		      char tmpstring[ MAXNAME ];
		      int nchars;
		      CPOS pos = delta * j;

		      nchars = snprintf(
			tmpstring,
			sizeof( tmpstring ),
			"%s%d",
			name,
			lNR++ );

		      if ( ( nchars < 0 ) ||
			   ( nchars >= MAXNAME ) )
			{
			  info(
			    -2,
			    "Generated locus name is too long\n" );

			  goto fail;
			}

		      if ( !append_mappoint(
			     &newMPs,
			     &cMPs,
			     chrom + 1,
			     pos,
			     tmpstring,
			     classname ) )
			goto fail_alloc;
		    }
		}

	      number = 1;
	      status = 0;
	    }

	  else if ( 0 ==
		    strcmp( token, "nonterminal" ) )
	    {
	      int chrom, j;
	      int lNR = 0;

	      if ( 0 != number % NoChroms )
		{
		  info(
		    -2,
		    "Number of nonterminal loci must be a multiple "
		    "of the number of chromosomes\n" );

		  goto fail;
		}

	      for ( chrom = 0;
		    chrom < NoChroms;
		    chrom++ )
		{
		  int nPerChrom =
		    number / NoChroms;

		  CPOS delta =
		    PopDescr[ chrom ].ChromLen /
		    ( nPerChrom + 1 );

		  for ( j = 1;j <= nPerChrom;j++ )
		    {
		      char tmpstring[ MAXNAME ];
		      int nchars;
		      CPOS pos = delta * j;

		      nchars = snprintf(
			tmpstring,
			sizeof( tmpstring ),
			"%s%d",
			name,
			lNR++ );

		      if ( ( nchars < 0 ) ||
			   ( nchars >= MAXNAME ) )
			{
			  info(
			    -2,
			    "Generated locus name is too long\n" );

			  goto fail;
			}

		      if ( !append_mappoint(
			     &newMPs,
			     &cMPs,
			     chrom + 1,
			     pos,
			     tmpstring,
			     classname ) )
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

	      if ( ( 1 != number ) ||
		   !parse_chrom_pos_token(
		     token,
		     &chrom,
		     &pos ) )
		{
		  info(
		    -2,
		    "Invalid map specification: %s\n",
		    token );

		  goto fail;
		}

	      if ( ( chrom < 1 ) ||
		   ( chrom > NoChroms ) )
		{
		  info(
		    -2,
		    "Wrong chromosome (%d)\n",
		    chrom );

		  goto fail;
		}

	      if ( ( pos < 0 ) ||
		   ( pos >
		     PopDescr[ chrom - 1 ].ChromLen ) )
		{
		  info(
		    -2,
		    "Chromosome %d: Locusposition outside of "
		    "the chromosome (0 <= %f <= %f)\n",
		    chrom,
		    pos,
		    PopDescr[ chrom - 1 ].ChromLen );

		  goto fail;
		}

	      if ( !append_mappoint(
		     &newMPs,
		     &cMPs,
		     chrom,
		     pos,
		     name,
		     classname ) )
		goto fail_alloc;

	      number = 1;
	      status = 0;
	    }

	  break;
	}

      token = strtok(
	NULL,
	"[]', \t\r\n" );
    }

  if ( 0 != status )
    {
      info( -2, "Incomplete map description\n" );
      goto fail;
    }

  if ( 0 == cMPs )
    {
      info( -2, "Linkage map contains no loci\n" );
      goto fail;
    }

  qsort(
    newMPs,
    cMPs,
    sizeof( TMappoint ),
    mnamecmp );

  for ( i = 1;i < cMPs;i++ )
    {
      if ( 0 == strcmp(
	     newMPs[ i - 1 ].name,
	     newMPs[ i ].name ) )
	{
	  info(
	    -2,
	    "Duplicate locus name: %s\n",
	    newMPs[ i ].name );

	  goto fail;
	}
    }

  qsort(
    newMPs,
    cMPs,
    sizeof( TMappoint ),
    mccmp );

  /*
   * Only now open the output file. All description and
   * allocation errors above therefore leave an existing
   * file untouched.
   */
  fp = fileopen( *fName, "w" );

  if ( NULL == fp )
    goto fail;

  outputOpened = 1;

  for ( i = 0;i < cMPs;i++ )
    {
      if ( 0 > fprintf(
	     fp,
	     "%d\t%.17g\t%s\t%s\n",
	     newMPs[ i ].chrom,
	     newMPs[ i ].Pos,
	     newMPs[ i ].name,
	     newMPs[ i ].class ) )
	{
	  info(
	    -2,
	    "Error writing map file %s\n",
	    *fName );

	  goto fail;
	}
    }

  if ( 0 != fclose( fp ) )
    {
      fp = NULL;

      info(
	-2,
	"Can't close the file stream\n" );

      goto fail;
    }

  fp = NULL;

  free_mappoints( newMPs, cMPs );
  free( description );

  *retval = 0;
  return;

fail_alloc:

  info(
    -2,
    "Memory allocation failed in generate_map_file\n" );

fail:

  if ( NULL != fp )
    {
      fclose( fp );
      fp = NULL;
    }

  if ( outputOpened )
    remove( *fName );

  free_mappoints( newMPs, cMPs );
  free( description );

  return;
}


void generate_effect_file(char**fName, char**Description )
{

    int status = 0;
    char*token, classname[ MAXNAME * 3 ], all0[ MAXNAME * 3 ], all1[ MAXNAME * 3 ];

    FILE*fp = fileopen( *fName, "w" );

    if ( NULL == fp )
        return;

    fprintf( fp, "# Automatisch erzeugt\n" );
    fprintf( fp, "# %s\n", *Description );

    GetRNGstate();

    token = strtok( *Description, "[' ',',']" );

    while ( NULL != token )
    {
        switch ( status )
        {

        case 0:
            {
                int tmpi;

                if ( 1 == sscanf( token, "%d", &tmpi ) )
                {
                    fprintf( fp, "%d\n", tmpi );

                }

                else
                {
                    if ( strlen( token ) >= sizeof( classname ) )
                    {
                        info( -2, "Class name too long: %s\n", token );
                        goto fail;
                    }

                    strcpy( classname, token );
                    strcpy( all0, "" );
                    strcpy( all1, "" );
                    status = 2;
                }

                break;
            }

        case 2:
            if ( strlen( token ) >= sizeof( all0 ) )
            {
                info( -2, "Allele name too long: %s\n", token );
                goto fail;
            }

            strcpy( all0, token );
            status = 4;
            break;

        case 3:
            {
                float mean, var;
                int chrom, i;

                if ( ( 2 == sscanf( token, "%f/%f", &mean, &var ) ) &&
                     ( 0.0 <= var ) )
                {
                    for ( chrom = 0;chrom < NoChroms;chrom++ )
                    {
                        for ( i = 0;i < PopDescr[ chrom ].NrMappoints;i++ )
                        {
                            if ( NULL != strstr( PopDescr[ chrom ].Mappoint[ i ].class, classname ) )
                            {
                                if ( 0 == strcmp( all1, "" ) )
                                    fprintf( fp, "%s\t%s\t%f\n",
                                             PopDescr[ chrom ].Mappoint[ i ].name,
                                             all0,
                                             norm_rand() * sqrt( var ) + mean );
                                else
                                    fprintf( fp, "%s\t%s\t%s\t%s\t%f\n",
                                             PopDescr[ chrom ].Mappoint[ i ].name,
                                             all0,
                                             PopDescr[ chrom ].Mappoint[ i ].name,
                                             all1,
                                             norm_rand() * sqrt( var ) + mean );
                            }
                        }
                    }

                    status = 0;
                }

                else 
                {
                    info( -2, "Error\n" );
                    goto fail;
                }
                break;
            }

        case 4:

            if ( 0 == strcmp( token, "normal" ) )
                status = 3;
            else if ( 0 == strcmp( token, "uniform" ) )
                status = 6;
            else
            {
                if ( strlen( token ) >= sizeof( all1 ) )
                {
                    info( -2, "Allele name too long: %s\n", token );
                    goto fail;
                }

                strcpy( all1, token );
                status = 5;
            }

            break;

        case 5:

            if ( 0 == strcmp( token, "normal" ) )
                status = 3;
            else if ( 0 == strcmp( token, "uniform" ) )
                status = 6;
            else {
                info( -2, "Fehler Verteilungsart erwartet: %s\n", token );
                goto fail;
            }

            break;

        case 6:
            {
                float eff;
                int chrom, i;

                if ( 1 == sscanf( token, "%f", &eff ) )
                {
                    for ( chrom = 0;chrom < NoChroms;chrom++ )
                    {
                        for ( i = 0;i < PopDescr[ chrom ].NrMappoints;i++ )
                        {
                            if ( NULL != strstr( PopDescr[ chrom ].Mappoint[ i ].class, classname ) )
                            {
                                if ( 0 == strcmp( all1, "" ) )
                                    fprintf( fp, "%s\t%s\t%f\n",
                                             PopDescr[ chrom ].Mappoint[ i ].name,
                                             all0,
                                             eff );
                                else
                                    fprintf( fp, "%s\t%s\t%s\t%s\t%f\n",
                                             PopDescr[ chrom ].Mappoint[ i ].name,
                                             all0,
                                             PopDescr[ chrom ].Mappoint[ i ].name,
                                             all1,
                                             eff );
                            }
                        }
                    }

                    status = 0;
                }

                else {
                    info( -2, "Error\n" );
                    goto fail;
                }
                break;
            }
        }

        token = strtok( NULL, "[' ',',']" );
    }

fail:

    PutRNGstate();
    fileclose( fp );
}


unsigned short get_allele( THom*Hom, CPOS Pos )
{
    int start, end, pos;

    if ( ( NULL == Hom ) ||
         ( NULL == Hom->Loc ) ||
         ( 1 > Hom->NoLociUsed ) )
        return 0;

    start = 0;
    end = Hom->NoLociUsed - 1;

    while ( start <= end )
      {
        pos = start + ( end - start ) / 2;

        if ( Pos < Hom->Loc[ pos ].Pos )
	  {
            end = pos - 1;
	  }
        else
	  {
            if ( ( pos + 1 < Hom->NoLociUsed ) &&
                 ( Hom->Loc[ pos + 1 ].Pos <= Pos ) )
	      {
                start = pos + 1;
	      }
            else
	      {
                return Hom->Loc[ pos ].All;
	      }
	  }
      }

    return Hom->Loc[ 0 ].All;
}


IntNode* insertIntLinkList( 
			   IntNode* Genotype,
			   int*     Phaenotyp,
			   int      counter,
			   int*     nrGt )
{
  IntNode*nodepos = Genotype;
  
  if ( NULL == Genotype )
    {
      Genotype = ( IntNode* ) yCalloc( 1, IntNode );
      Genotype->child = NULL;
      Genotype->parent = NULL;
      Genotype->prev = NULL;
      Genotype->next = NULL;

      if ( 0 == counter )
        {
	  Genotype->number = 1;
	  ( *nrGt ) ++;
        }
      else
        {
	  Genotype->number = Phaenotyp[ 0 ];
	  Genotype->child = insertIntLinkList( Genotype->child, &Phaenotyp[ 1 ], counter - 1, nrGt );
	  Genotype->child->parent = Genotype;
        }

        return Genotype;
    }

  if ( 0 == counter )
    {
      ( Genotype->number ) ++;
      return Genotype;
    }
  do
    {
      if ( nodepos->number != Phaenotyp[ 0 ] )
        {
	  if ( nodepos->number > Phaenotyp[ 0 ] )
            {
	      
	      IntNode * tmp;
	      tmp = ( IntNode* ) yCalloc( 1, IntNode );
	      tmp->child = NULL;
	      tmp->parent = nodepos->parent;
	      tmp->prev = NULL;
	      tmp->next = nodepos;
	      nodepos->prev = tmp;
	      
	      if ( nodepos->parent != NULL )
		nodepos->parent->child = tmp;
	      
	      tmp->number = Phaenotyp[ 0 ];
	      
	      Genotype = nodepos = tmp;
            }
	  else
            {
	      if ( ( NULL != nodepos->next ) && ( nodepos->next->number <= Phaenotyp[ 0 ] ) )
                {
		  nodepos = nodepos->next;
                }
	      else
                {
		  IntNode*tmp;
		  tmp = ( IntNode* ) yCalloc( 1, IntNode );
		  tmp->child = NULL;
		  tmp->parent = nodepos->parent;
		  tmp->prev = nodepos;
		  tmp->next = nodepos->next;
		  tmp->number = Phaenotyp[ 0 ];
		  nodepos->next = tmp;
		  
		  if ( tmp->next != NULL )
                    {
		      tmp->next->prev = tmp;
                    }
		  
		  nodepos = tmp;
                }
            }
        }
    }

  while ( nodepos->number < Phaenotyp[ 0 ] );

  nodepos->child = insertIntLinkList( nodepos->child,
				      &Phaenotyp[ 1 ],
				      counter - 1,
				      nrGt );
  
  nodepos->child->parent = nodepos;
  
  return Genotype;

}

IntNode* go_down ( IntNode *Node )
{
  while ( Node->child != NULL ) Node = Node->child;
  return Node;
}

void print_path( IntNode*Node,
		 int*_nrLoci,
		 int*_nrGt,
		 int*m )
{
    int i = *_nrLoci;

    do {
      m[ ( i-- ) * ( *_nrGt ) ] = Node->number;
    }
    while ( NULL != ( Node = Node->parent ) );
}


IntNode* find_next_leaf( IntNode* Node)
{
    IntNode*tmp;

    while ( NULL != Node->parent )
      {
        tmp = Node;
        Node = Node->parent;
        yFree( tmp );
	
        if ( Node->next != NULL )
	  {
            tmp = Node->next;
            yFree( Node );
            return ( go_down( tmp ) );
	  }
      }

    yFree( Node );
    return ( NULL );
}


SEXP genotype_evaluate( SEXP _PopName,
			SEXP _Loci,
			SEXP _mode )
{
    TPopListEl*PopListEl = NULL;
    int mode = 0;
    char*PopName = NULL;

    if ( 0 == cMP ){
      info( -2, "No linkage map defined\n" );
      return R_NilValue;
    }

    if ( !isInteger( _mode ) || ( length( _mode ) != 1 ) ) 
      {
        info( -2, "Value of mode is not an integer\n" );
	return R_NilValue;
      }
    else
      {
        mode = INTEGER( _mode ) [ 0 ];
      }

    if ( !isString( _PopName ) || ( length( _PopName ) != 1 ) )
      {
        info( -2, "Name of population is not a single string\n" );
	return R_NilValue;
      }

    PopName = ( char* ) yCalloc(
        strlen( CHAR( STRING_ELT( _PopName, 0 ) ) ) + 2,
        char );

    strcpy( PopName, CHAR( STRING_ELT( _PopName, 0 ) ) );

    if ( !isString( _Loci ) )
      {
        info( -2, "Name of Loci is not a string\n" );
	goto fail;
      }

    find_list_element( PopName, &PopListEl );

    if ( NULL == PopListEl )
      {
	info( -2, "Population in not defined: %s\n", PopName );
	goto fail;
      }

    else
    {

        typedef struct
        {
            int nr_marker;
            TMappoint*Mappoint;
            int*order;
        }

        TEvChrom;

        TEvChrom*EvChrom = ( TEvChrom* ) yCalloc( NoChroms, TEvChrom );
        int chrom, i, j, l;

        int nrLoci = 0;
        int*pnrLoci = &nrLoci;
        int nrGt = 0;
        int*pnrGt = &nrGt;

        SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;
        int*locus_idx = ( int* ) yCalloc( length( _Loci ), int );

        POPULATION ind;
        CPOS Pos;

        startNode = NULL;


        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            EvChrom[ chrom ].nr_marker = 0;
            EvChrom[ chrom ].Mappoint = NULL;
            EvChrom[ chrom ].order = NULL;
        }

        for ( l = 0;l < length( _Loci );l++ )
        {
            locus_idx[ l ] = 0;

            for ( chrom = 0;( ( chrom < NoChroms ) &&
                              ( 0 == locus_idx[ l ] ) );chrom++ )
            {
                for ( j = 0;( ( j < PopDescr[ chrom ].NrMappoints ) &&
                              ( 0 == locus_idx[ l ] ) );j++ )
                {
                    if ( 0 == strcmp(
                                PopDescr[ chrom ].Mappoint[ j ].name,
                                CHAR( STRING_ELT( _Loci, l ) ) ) )
                    {
                        locus_idx[ l ] = 1;
                        EvChrom[ chrom ].nr_marker++;

                        if ( EvChrom[ chrom ].nr_marker == 1 )
                        {
                            EvChrom[ chrom ].Mappoint =
                                ( TMappoint* ) yCalloc( 1, TMappoint );

                            EvChrom[ chrom ].order =
                                ( int* ) yCalloc( 1, int );
                        }

                        else
                        {
                            EvChrom[ chrom ].Mappoint =
                                ( TMappoint* ) yRealloc(
                                    EvChrom[ chrom ].Mappoint,
                                    EvChrom[ chrom ].nr_marker,
                                    TMappoint );

                            EvChrom[ chrom ].order =
                                ( int* ) yRealloc(
                                    EvChrom[ chrom ].order,
                                    EvChrom[ chrom ].nr_marker,
                                    int );
                        }

                        EvChrom[ chrom ].Mappoint[
                            EvChrom[ chrom ].nr_marker - 1 ] =
                            PopDescr[ chrom ].Mappoint[ j ];

                        EvChrom[ chrom ].order[
                            EvChrom[ chrom ].nr_marker - 1 ] = nrLoci;

                        nrLoci += NoHoms;
                    }
                }
            }
        }

        if ( 0 == nrLoci )
	  {
            info( -2, "Loci not found\n" );
	    goto fail_data;
	  }

        if ( 0 == PopListEl->Pop->NoInds )
	  {
            info( -2, "Empty population: %s\n", PopName );
	    goto fail_data;
	  }

        {
            TChrom*Chrom;
            TInd*Ind;
            IntNode*Node;
            int*alleles = ( int* ) yCalloc( NoHoms, int );
            int*IndGT = ( int* ) yCalloc( nrLoci + 1, int );
            int hom, tmp_counter = 0;

            for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
            {
                Ind = get_individual( PopListEl->Pop, ind );

                for ( chrom = 0;chrom < NoChroms;chrom++ )
                {
                    for ( l = 0;l < EvChrom[ chrom ].nr_marker;l++ )
                    {
                        int offset = EvChrom[ chrom ].order[ l ];

                        if ( NULL == Ind->genotype )
                        {
                            Chrom = &( Ind->Chrom[ chrom ] );
                            Pos = EvChrom[ chrom ].Mappoint[ l ].Pos;

                            for ( hom = 0;hom < NoHoms;hom++ )
                            {
                                alleles[ hom ] =
                                    get_allele( &( Chrom->Hom[ hom ] ),
                                                Pos );
                            }
                        }

                        else
                        {
                            for ( hom = 0;hom < NoHoms;hom++ )
                            {
                                alleles[ hom ] =
                                    Ind->genotype[
                                        EvChrom[ chrom ].
                                        Mappoint[ l ].counter *
                                        NoHoms + hom ];
                            }
                        }

                        switch ( mode )
                        {

                        case 0:
                            break;

                        case 1:

                            if ( EvChrom[ chrom ].nr_marker == 1 )
                                R_isort( &alleles[ 0 ], NoHoms );

                            break;

                        default:
                            R_isort( &alleles[ 0 ], NoHoms );

                            break;
                        }

                        for ( hom = 0;hom < NoHoms;hom++ )
                        {
                            IndGT[ offset + hom ] = alleles[ hom ];
                        }
                    }
                }

                startNode =
                    insertIntLinkList( startNode,
                                       IndGT,
                                       nrLoci,
                                       pnrGt );
            }

            {
                int*_matrix =
                    ( int* ) yCalloc( ( nrLoci + 1 ) * nrGt, int );
                int offset;
                char tmp[ 1024 ];
                int N = 0;

                Node = startNode;

                if ( NULL != ( Node = go_down( Node ) ) )
                {
                    do
                    {
                        print_path( Node,
                                    pnrLoci,
                                    pnrGt,
                                    &_matrix[ tmp_counter++ ] );
                    }

                    while ( NULL !=
                            ( Node = find_next_leaf( Node ) ) );

                    startNode = NULL;
                }

                PROTECT( matrix =
                         allocMatrix( INTSXP, nrGt, nrLoci + 2 ) );
                PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
                INTEGER( matrixDim ) [ 0 ] = nrGt;
                INTEGER( matrixDim ) [ 1 ] = nrLoci + 2;
                setAttrib( matrix, R_DimSymbol, matrixDim );

                PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
                PROTECT( xAxis = allocVector( STRSXP, nrLoci + 2 ) );
                PROTECT( yAxis = allocVector( STRSXP, nrGt ) );

                for ( i = offset = 0;i < nrLoci + 1;i++ )
                {
                    for ( j = 0;j < nrGt;j++, offset++ )
                    {
                        INTEGER( matrix ) [ offset ] =
                            _matrix[ offset ];

                        if ( i == nrLoci )
                        {
                            N += _matrix[ offset ];
                        }
                    }
                }

                for ( j = 0;j < nrGt;j++, offset++ )
                {
                    INTEGER( matrix ) [ offset ] = N;
                }

                for ( i = 0, offset = 0;
                      i < length( _Loci );
                      i++ )
                {
                    if ( locus_idx[ i ] )
                    {
                        for ( hom = 1;hom <= NoHoms;hom++ )
                        {
			  snprintf( tmp,
                                    sizeof( tmp ),
				    "%s-%d ",
				    CHAR( STRING_ELT( _Loci, i ) ),
                                    hom );

                            SET_STRING_ELT( xAxis,
                                            offset++,
                                            mkChar( tmp ) );
                        }
                    }
                }

                SET_STRING_ELT( xAxis,
                                offset++,
                                mkChar( "count" ) );

                SET_STRING_ELT( xAxis,
                                offset++,
                                mkChar( "frequency" ) );

                for ( i = 0;i < nrGt;i++ )
                {
		  snprintf( tmp,
                            sizeof( tmp ),
			    "%d ",
                            i + 1 );

		  SET_STRING_ELT( yAxis,
                                  i,
                                  mkChar( tmp ) );
                }

                SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
                SET_VECTOR_ELT( matrixDimNames, 0, yAxis );
                setAttrib( matrix,
                           R_DimNamesSymbol,
                           matrixDimNames );

                UNPROTECT( 5 );

		yFree( PopName );

		for ( chrom = 0;chrom < NoChroms;chrom++ )
                {
		    yFree( EvChrom[ chrom ].Mappoint );
		    yFree( EvChrom[ chrom ].order );
                }

		yFree( EvChrom );
		yFree( locus_idx );
		yFree( alleles );
		yFree( IndGT );
		yFree( _matrix );

                return ( matrix );
            }
        }

fail_data:

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            yFree( EvChrom[ chrom ].Mappoint );
            yFree( EvChrom[ chrom ].order );
        }

        yFree( EvChrom );
        yFree( locus_idx );
    }

fail:

    yFree( PopName );
    return(R_NilValue);
}

SEXP allele_evaluate(
		     SEXP _PopName,
		     SEXP _Loci,
		     SEXP _missing,
		     SEXP _bg)

{
    char*PopName = NULL;
    int bg = 0;
    int missing = 0;
    TPopListEl*PopListEl = NULL;

    if ( 0 == cMP ){
        info( -2, "No linkage map defined\n" );
	return R_NilValue;
    }

    if ( !isString( _PopName ) || ( length( _PopName ) != 1 ) )
    {
        info( -2, "Name of population is not a single string\n" );
	return R_NilValue;
    }

    PopName = ( char* ) yCalloc(
        strlen( CHAR( STRING_ELT( _PopName, 0 ) ) ) + 2,
        char );

    strcpy( PopName, CHAR( STRING_ELT( _PopName, 0 ) ) );

    if ( !isString( _Loci ) ){
        info( -2, "Name of Loci is not a string\n" );
	goto fail;
    }

    if ( !isInteger( _bg ) || ( length( _bg ) != 1 ) ) {
        info( -2, "Value of background allele is not an integer\n" );
	goto fail;
    }
    else
        bg = INTEGER( _bg ) [ 0 ];

    if ( !isInteger( _missing ) || ( length( _missing ) != 1 ) ) {
        info( -2, "Value of missing is not an integer\n" );
	goto fail;
    }
    else
        missing = INTEGER( _missing ) [ 0 ];

    find_list_element( PopName, &PopListEl );

    if ( NULL == PopListEl )
    {
        info( -2, "Population not defined: %s\n", PopName );
	goto fail;
    }

    else
    {

        int*locus_idx = ( int* ) yCalloc( length( _Loci ), int );

        int nrGt = 0;
        int*pnrGt = &nrGt;

        typedef struct
        {
            int nr_marker;
            TMappoint*Mappoint;
            int*order;
        }

        TEvChrom;

        TEvChrom*EvChrom = ( TEvChrom* ) yCalloc( NoChroms, TEvChrom );
        int j, l, chrom;
        int nrLoci = 0;
        int*pnrLoci = &nrLoci;

        startNode = NULL;


        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            EvChrom[ chrom ].nr_marker = 0;
            EvChrom[ chrom ].Mappoint = NULL;
            EvChrom[ chrom ].order = NULL;
        }


        for ( l = 0;l < length( _Loci );l++ )
        {
            locus_idx[ l ] = 0;

            for ( chrom = 0;( ( chrom < NoChroms ) &&
                              ( 0 == locus_idx[ l ] ) );chrom++ )
            {
                for ( j = 0;( ( j < PopDescr[ chrom ].NrMappoints ) &&
                              ( 0 == locus_idx[ l ] ) );j++ )
                {
                    if ( 0 == strcmp(
                                PopDescr[ chrom ].Mappoint[ j ].name,
                                CHAR( STRING_ELT( _Loci, l ) ) ) )
                    {
                        locus_idx[ l ] = 1;
                        EvChrom[ chrom ].nr_marker++;

                        if ( 1 == EvChrom[ chrom ].nr_marker )
                        {
                            EvChrom[ chrom ].Mappoint =
                                ( TMappoint* ) yCalloc( 1, TMappoint );
                            EvChrom[ chrom ].order =
                                ( int* ) yCalloc( 1, int );
                        }

                        else
                        {
                            EvChrom[ chrom ].Mappoint =
                                ( TMappoint* ) yRealloc(
                                    EvChrom[ chrom ].Mappoint,
                                    EvChrom[ chrom ].nr_marker,
                                    TMappoint );

                            EvChrom[ chrom ].order =
                                ( int* ) yRealloc(
                                    EvChrom[ chrom ].order,
                                    EvChrom[ chrom ].nr_marker,
                                    int );
                        }

                        EvChrom[ chrom ].Mappoint[
                            EvChrom[ chrom ].nr_marker - 1 ] =
                            PopDescr[ chrom ].Mappoint[ j ];

                        EvChrom[ chrom ].order[
                            EvChrom[ chrom ].nr_marker - 1 ] = nrLoci;

                        nrLoci++;
                    }
                }
            }
        }


        if ( 0 == nrLoci )
        {
            info( -2, "Loci not found\n" );
	    goto fail_data;
        }

        else if ( 0 == PopListEl->Pop->NoInds )
        {
            info( -2, "Empty population %s\n", PopName );
	    goto fail_data;
        }

        else
        {
            CPOS Pos;
            int chrom, hom, l, temp, i;
            int*all = ( int* ) yCalloc( nrLoci, int );
            POPULATION ind;
            TInd*Ind;
            TChrom*Chrom;

            for ( i = 0;i < nrLoci;i++ )
            {
                all[ i ] = 0;
            }

            for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
            {
                Ind = get_individual( PopListEl->Pop, ind );

                for ( chrom = 0;chrom < NoChroms;chrom++ )
                {
                    for ( l = 0;l < EvChrom[ chrom ].nr_marker;l++ )
                    {
                        temp = EvChrom[ chrom ].order[ l ];

                        if ( NULL == Ind->genotype )
                        {
                            Chrom = &( Ind->Chrom[ chrom ] );
                            Pos = EvChrom[ chrom ].Mappoint[ l ].Pos;

                            for ( hom = 0;hom < NoHoms;hom++ )
                            {
                                if ( 0 == ( all[ temp ] =
                                          get_allele(
                                              &( Chrom->Hom[ hom ] ),
                                              Pos ) ) )
                                    info( -2, "Allele 0 not allowed\n" );

                                if ( ( 1 == missing ) &&
                                     ( all[ temp ] == bg ) )
                                {
                                    all[ temp ] = 0;
                                }

                                else
                                {
                                    startNode = insertIntLinkList(
                                        startNode,
                                        all,
                                        nrLoci,
                                        pnrGt );
                                }
                            }
                        }

                        else
                        {
                            for ( hom = 0;hom < NoHoms;hom++ )
                            {
                                if ( 0 == ( all[ temp ] =
                                          Ind->genotype[
                                              EvChrom[ chrom ].
                                              Mappoint[ l ].counter *
                                              NoHoms + hom ] ) )
                                    info( -2, "Allele 0 not allowed\n" );

                                if ( ( 1 == missing ) &&
                                     ( all[ temp ] == bg ) )
                                {
                                    all[ temp ] = 0;
                                }

                                else
                                {
                                    startNode = insertIntLinkList(
                                        startNode,
                                        all,
                                        nrLoci,
                                        pnrGt );
                                }
                            }
                        }

                        all[ temp ] = 0;
                    }
                }
            }


            for ( chrom = 0;chrom < NoChroms;chrom++ )
            {
                yFree( EvChrom[ chrom ].Mappoint );
                yFree( EvChrom[ chrom ].order );
            }

            {
                int*_matrix =
                    ( int* ) yCalloc( ( nrLoci + 1 ) * nrGt, int );
                int offset;
                char tmp[ 1024 ];
                int N = 0;
                IntNode*Node;
                int tmp_counter = 0;
                SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;


                Node = startNode;

                if ( ( NULL != Node ) &&
                     ( NULL != ( Node = go_down( Node ) ) ) )
                {
                    do
                    {
                        print_path( Node,
                                    pnrLoci,
                                    pnrGt,
                                    &_matrix[ tmp_counter++ ] );
                    }

                    while ( NULL !=
                            ( Node = find_next_leaf( Node ) ) );

                    startNode = NULL;
                }

                PROTECT( matrix =
                         allocMatrix( INTSXP, nrGt, nrLoci + 2 ) );
                PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
                INTEGER( matrixDim ) [ 0 ] = nrGt;
                INTEGER( matrixDim ) [ 1 ] = nrLoci + 2;
                setAttrib( matrix, R_DimSymbol, matrixDim );
                PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
                PROTECT( xAxis = allocVector( STRSXP, nrLoci + 2 ) );
                PROTECT( yAxis = allocVector( STRSXP, nrGt ) );

                for ( i = offset = 0;i < nrLoci + 1;i++ )
                {
                    for ( j = 0;j < nrGt;j++, offset++ )
                    {
                        INTEGER( matrix ) [ offset ] =
                            _matrix[ offset ];

                        if ( i == nrLoci )
                        {
                            N += _matrix[ offset ];
                        }
                    }
                }

                for ( j = 0;j < nrGt;j++, offset++ )
                {
                    INTEGER( matrix ) [ offset ] = N;
                }

                for ( i = 0, offset = 0;
                      i < length( _Loci );
                      i++ )
                {
                    if ( locus_idx[ i ] )
                    {
		      snprintf( tmp,
                                sizeof( tmp ),
			        "%s",
			        CHAR( STRING_ELT( _Loci, i ) ) );

		      SET_STRING_ELT( xAxis,
                                      offset++,
                                      mkChar( tmp ) );
                    }
                }

                SET_STRING_ELT( xAxis,
                                offset++,
                                mkChar( "count" ) );
                SET_STRING_ELT( xAxis,
                                offset++,
                                mkChar( "frequency" ) );

                for ( i = 0;i < nrGt;i++ )
                {
		  snprintf( tmp,
                            sizeof( tmp ),
			    "%d ",
                            i + 1 );

		  SET_STRING_ELT( yAxis,
                                  i,
                                  mkChar( tmp ) );
                }

                SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
                SET_VECTOR_ELT( matrixDimNames, 0, yAxis );
                setAttrib( matrix,
                           R_DimNamesSymbol,
                           matrixDimNames );

                UNPROTECT( 5 );

                yFree( PopName );
                yFree( locus_idx );
                yFree( EvChrom );
                yFree( all );
                yFree( _matrix );

                return ( matrix );
            }
        }

fail_data:

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            yFree( EvChrom[ chrom ].Mappoint );
            yFree( EvChrom[ chrom ].order );
        }

        yFree( locus_idx );
        yFree( EvChrom );
    }

fail:

    yFree( PopName );
    return R_NilValue;
}

void evaluate_locus( 
		    char**PopName,
		    char**Loci,
		    int*_matrix,
		    int*noLoci,
		    int*bg )

{
    TPopListEl*PopListEl;
    char*pch;
    int chrom, i = 0, j;
    POPULATION noHomo, noHetero, noMissing, ind;

    if ( 2 != NoHoms ){
        info( -2, "This function is only defined for diploid organisms!\n" );
	return;
    }


    if ( 0 == cMP ){
        info( -2, "No linkage map defined\n" );
 	return;
   }


    find_list_element( *PopName, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Population not defined: %s\n", *PopName );
	return;
    }

    pch = strtok( *Loci, " " );

    while ( pch != NULL )
    {
        noHomo = 0;
        noHetero = 0;
        noMissing = 0;

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            if ( PopDescr[ chrom ].NrMappoints > 0 )
            {
                for ( j = 0;j < PopDescr[ chrom ].NrMappoints;j++ )
                {
                    if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ j ].name, pch ) )
                    {
                        for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
                        {
                            TInd*Ind = get_individual( PopListEl->Pop, ind );

                            if ( NULL == Ind->genotype )
                            {
                                info( -2, "first genotype_population (chrom=%d, j=%d, ind=%ld, pch=%s)",
                                      chrom, j, ind, pch );
				return;
                            }

                            if ( ( Ind->genotype[ 2 * PopDescr[ chrom ].Mappoint[ j ].counter ] == *bg ) ||
                                    ( Ind->genotype[ 2 * PopDescr[ chrom ].Mappoint[ j ].counter + 1 ] == *bg ) )
                            {
                                noMissing++;
                            }

                            else if ( Ind->genotype[ 2 * PopDescr[ chrom ].Mappoint[ j ].counter ] !=
                                    Ind->genotype[ 2 * PopDescr[ chrom ].Mappoint[ j ].counter + 1 ] )
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

        _matrix[ i ] = noHomo;
        _matrix[ i + ( *noLoci ) ] = noHetero;
        _matrix[ i + 2 * ( *noLoci ) ] = noMissing;
        i++;
        pch = strtok( NULL, " " );
    }
}


SEXP evaluate_haplotype( 
			SEXP _PopName,
			SEXP _Loci,
			SEXP _missing,
			SEXP _bg)
{
    int bg = 0, missing = 0;
    SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;


    if ( !isString( _PopName ) || ( length( _PopName ) != 1 ) )
      {
        info( -2, "Name of population is not a single string\n" );
	return  R_NilValue;
      }


    if ( !isString( _Loci ) || ( length( _Loci ) != 1 ) )
      {
        info( -2, "Name of Loci is not a single string\n" );
	return R_NilValue;
      }

    if ( !isInteger( _bg ) || ( length( _bg ) != 1 ) )
      {
        info( -2, "Value of background allele is not an integer\n" );
	return  R_NilValue;
      }
    else
      {
        bg = INTEGER( _bg ) [ 0 ];
      }

    if ( !isInteger( _missing ) || ( length( _missing ) != 1 ) )
      {
        info( -2, "Value of missing variable is not an integer\n" );
	return R_NilValue;
      }
    else
      {
        missing = INTEGER( _missing ) [ 0 ];
      }

    {
        TPopListEl*PopListEl;
        POPULATION ind;
        TMappoint*Mappoint = NULL;

        int nLoci = 0;
        int nGt = 0;
        int*pGt = &nGt;
        int*hap = NULL;
        int free_startNode = 0;

        char*PopName = ( char* ) yCalloc( strlen( CHAR( STRING_ELT( _PopName, 0 ) ) ) + 2, char );
        char*Loci = ( char* ) yCalloc( strlen( CHAR( STRING_ELT( _Loci, 0 ) ) ) + 2, char );

        strcpy( PopName, CHAR( STRING_ELT( _PopName, 0 ) ) );
        strcpy( Loci, CHAR( STRING_ELT( _Loci, 0 ) ) );

        find_list_element( PopName, &PopListEl );

        if ( !PopListEl ){
            info( -2, "Population not defined: %s\n", PopName );
	    goto fail;
	}

        if ( 1 > PopListEl->Pop->NoInds ){
            info( -2, "Empty population %s\n", PopName );
	    goto fail;
	}

        {
            int lastLoci;
            int chrom, loc;
            char*locus;

            locus = strtok( Loci, " " );

            while ( locus != NULL )
            {
                lastLoci = nLoci;

                for ( chrom = 0;chrom < NoChroms;chrom++ )
                {
                    if ( PopDescr[ chrom ].NrMappoints > 0 )
                    {
                        for ( loc = 0;loc < PopDescr[ chrom ].NrMappoints;loc++ )
                        {
                            if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ loc ].name, locus ) )
                            {
                                if ( 0 == nLoci )
                                {
                                    Mappoint = ( TMappoint* ) yCalloc( 1, TMappoint );
                                }

                                else
                                {
                                    Mappoint = ( TMappoint* ) yRealloc( Mappoint, nLoci + 1, TMappoint );
                                }

                                Mappoint[ nLoci ] = ( PopDescr[ chrom ].Mappoint[ loc ] );
                                nLoci++;
                                break;
                            }
                        }
                    }
                }

                if ( lastLoci == nLoci )
                {
                    info( -2, "Cant find locus %s\n ", locus );
		    goto fail;
                }

                locus = strtok( NULL, " " );
            }
        }


        {
            TInd*Ind;
            int loc, hom, nrmissing;

            hap = ( int* ) yCalloc( nLoci, int );


            startNode = NULL;

            for ( loc = 0;loc < nLoci;loc++ )
            {
                hap[ loc ] = 0;
            }

            startNode = insertIntLinkList( startNode, hap, nLoci, pGt );
            free_startNode = 1;

            for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
            {
                Ind = get_individual( PopListEl->Pop, ind );

                if ( NULL == Ind->genotype ) {
                    info( -2, "First genotype_population\n" );
		    goto fail;
		}

                for ( hom = 0;hom < NoHoms;hom++ )
                {
                    nrmissing = 0;

                    for ( loc = 0;loc < nLoci;loc++ )
                    {
                        hap[ loc ] = Ind->genotype[
                            Mappoint[ loc ].counter * NoHoms + hom ];

                        if ( ( 1 == missing ) && ( bg == hap[ loc ] ) )
                        {
                            nrmissing++;
                        }
                    }

                    if ( 0 == nrmissing )
                    {
                        startNode = insertIntLinkList( startNode, hap, nLoci, pGt );
                    }
                }
            }

            yFree( hap );
        }

        {
            int y_dim = nGt;
            int x_dim = nLoci + 1;
            IntNode*Node = startNode;
            int tmpcounter = 0;
            char*locus;

            PROTECT( matrix = allocMatrix( INTSXP, y_dim, x_dim ) );
            PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
            INTEGER( matrixDim ) [ 0 ] = y_dim;
            INTEGER( matrixDim ) [ 1 ] = x_dim;
            setAttrib( matrix, R_DimSymbol, matrixDim );
            PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
            PROTECT( xAxis = allocVector( STRSXP, x_dim ) );
            PROTECT( yAxis = allocVector( STRSXP, y_dim ) );

            if ( ( NULL != Node ) && ( NULL != ( Node = go_down( Node ) ) ) )
            {
                int * path = ( int* ) yCalloc( nGt * x_dim, int );
                int ttt = 0;
                int iii;

                do
                {
                    print_path( Node, &nLoci, &nGt, &path[ ttt++ ] );
                }

                while ( NULL != ( Node = find_next_leaf( Node ) ) );

                startNode = NULL;

                for ( iii = 0;iii < x_dim;iii++ )
                {
                    for ( ttt = 0;ttt < y_dim;ttt++, tmpcounter++ )
                    {
                        INTEGER( matrix ) [ tmpcounter ] = path[ tmpcounter ];
                    }
                }

                yFree( path );
            }

            strcpy( Loci, CHAR( STRING_ELT( _Loci, 0 ) ) );
            locus = strtok( Loci, " " );
            {
                int col = 0;

                while ( locus != NULL )
                {
                    SET_STRING_ELT( xAxis, col++, mkChar( locus ) );
                    locus = strtok( NULL, " " );
                }

                SET_STRING_ELT( xAxis, col, mkChar( "count" ) );
            }

            SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
            setAttrib( matrix, R_DimNamesSymbol, matrixDimNames );
            UNPROTECT( 5 );
        }

        if ( 2 <= nLoci )
        {
            int x;
            int**table = NULL;
            int dimLA = 0, dimLB = 0;
            int lx, ly, xPos, yPos;

            for ( x = 0;x < nGt;x++ )
            {
                xPos = 0;
                yPos = 0;

                if ( ( 0 == dimLA ) && ( 0 == dimLB ) )
                {
                    dimLA = 1;
                    dimLB = 1;
                    table =  yCalloc( dimLA + 1, int* );
                    table[ 0 ] = ( int* ) yCalloc( dimLB + 1, int );
                    table[ 1 ] = ( int* ) yCalloc( dimLB + 1, int );
                    table[ 0 ][ 0 ] = INTEGER( matrix ) [ x + 2 * nGt ];
                    table[ 1 ][ 0 ] = INTEGER( matrix ) [ x ];
                    table[ 0 ][ 1 ] = INTEGER( matrix ) [ x + nGt ];
                    table[ 1 ][ 1 ] = INTEGER( matrix ) [ x + 2 * nGt ];
                }

                else
                {
                    for ( lx = 1;lx <= dimLB;lx++ )
                    {
                        if ( table[ 0 ][ lx ] == INTEGER( matrix ) [ x + nGt ] )
                        {
                            xPos = lx;
                        }
                    }

                    if ( 0 == xPos )
                    {
                        dimLB++;

                        for ( ly = 0;ly <= dimLA;ly++ )
                        {
                            table[ ly ] = ( int* ) yRealloc( table[ ly ], dimLB + 1, int );
                            table[ ly ][ dimLB ] = 0;
                        }

                        lx = 1;

                        while ( table[ 0 ][ lx ] < INTEGER( matrix ) [ x + nGt ] )
                        {
                            lx++;

                            if ( lx >= dimLB )
                                break;
                        }

                        xPos = lx;

                        if ( lx != dimLB )
                        {
                            for ( ly = 0;ly <= dimLA;ly++ )
                            {
                                for ( lx = dimLB;lx > xPos;lx-- )
                                {
                                    table[ ly ][ lx ] = table[ ly ][ lx - 1 ];
                                }
                            }

                            for ( ly = 0;ly <= dimLA;ly++ )
                            {
                                table[ ly ][ xPos ] = 0;
                            }
                        }

                        table[ 0 ][ xPos ] = INTEGER( matrix ) [ x + nGt ];
                    }

                    for ( ly = 1;ly <= dimLA;ly++ )
                    {
                        if ( table[ ly ][ 0 ] == INTEGER( matrix ) [ x ] )
                        {
                            yPos = ly;
                        }
                    }

                    if ( 0 == yPos )
                    {
                        dimLA++;
                        table = ( int** ) yRealloc( table, dimLA + 1, int* );
                        table[ dimLA ] = ( int* ) yCalloc( dimLB + 1, int );

                        table[ dimLA ][ 0 ] = INTEGER( matrix ) [ x ];

                        for ( lx = 1;lx <= dimLB;lx++ )
                        {
                            table[ dimLA ][ lx ] = 0;
                        }

                        yPos = dimLA;
                    }

                    if ( ( xPos != 0 ) && ( yPos != 0 ) )
                    {
                        table[ 0 ][ 0 ] += INTEGER( matrix ) [ x + 2 * nGt ];
                        table[ yPos ][ xPos ] = INTEGER( matrix ) [ x + 2 * nGt ];
                    }
                }
            }

            for ( ly = 0;ly <= dimLA;ly++ )
            {
                yFree( table[ ly ] );
            }

            yFree( table );
        }

	yFree( Mappoint );
        yFree( PopName );
        yFree( Loci );

        return ( matrix );

fail:

        if ( free_startNode )
        {
            IntNode*Node = go_down( startNode );

            while ( NULL != ( Node = find_next_leaf( Node ) ) )
            {
            }

            startNode = NULL;
        }

        yFree( hap );
	yFree( Mappoint );
        yFree( PopName );
        yFree( Loci );

        return R_NilValue;
    }
}

void haplotype_init( 
		    char**PopName,
		    int*_missing,
		    int*_doubleHetero,
		    int*_bg)

{
    PopListEl = NULL;
    find_list_element( *PopName, &PopListEl );

    if ( !PopListEl ) {
      info( -2, "Population not defined: %s\n", *PopName );
      return ;
    }

    if ( 1 > PopListEl->Pop->NoInds ){
      info( -2, "Population %s is empty\n", *PopName );
      return ;
    }

    bg = *_bg;

    missing = *_missing;

    doubleHetero = *_doubleHetero;
}



void haplotype_next(
		    int*_locusA,
		    int*_locusB,
		    int*_noTI,
		    int*_noI,
		    int*_nodH,
		    int*_dLocA,
		    int*_dLocB)

{
    int nGt = 0;
    int*pGt = &nGt;
    TInd*Ind;
    POPULATION ind;
    int nLoci = 2, hap[ 2 ], dH = 0;

    *_noTI = 0;
    *_noI = 0;
    *_nodH = 0;
    *_dLocA = 0;
    *_dLocB = 0;

    if ( ( NULL == PopListEl ) || ( NULL == PopListEl->Pop ) )
    {
        info( -2, "First initialize haplotype evaluation\n" );
	return;
    }

    if ( 2 != NoHoms )
    {
        info( -2, "Haplotype evaluation is implemented for diploids only\n" );
	return;
    }

    if ( ( 0 > *_locusA ) || ( cMP <= *_locusA ) ||
         ( 0 > *_locusB ) || ( cMP <= *_locusB ) )
    {
        info( -2, "Locus index out of range\n" );
	return;
    }

    if ( NULL != table )
    {
        int ly;

        for ( ly = 0;ly <= dimLA;ly++ )
        {
            yFree( table[ ly ] );
        }

        yFree( table );
    }

    yFree( genotype );
    nrGeno = 0;
    dimLA = 0;
    dimLB = 0;

    locusA = *_locusA;
    locusA *= 2;
    locusB = *_locusB;
    locusB *= 2;

    startNode = NULL;
    *_noTI = PopListEl->Pop->NoInds;

    for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
    {
        Ind = get_individual( PopListEl->Pop, ind );

        if ( NULL == Ind->genotype ){
            IntNode*Node;

            info( -2, "first genotype_population\n" );

            if ( NULL != startNode )
            {
                Node = go_down( startNode );

                while ( NULL != ( Node = find_next_leaf( Node ) ) )
                {
                }

                startNode = NULL;
            }

            yFree( genotype );
            nrGeno = 0;
            *_noI = 0;

	    return ;
	}

        if ( doubleHetero &&
                ( Ind->genotype[ locusA ] != Ind->genotype[ locusA + 1 ] ) &&
                ( Ind->genotype[ locusB ] != Ind->genotype[ locusB + 1 ] ) )
        {
            dH++;
        }

        else
        {
            if ( missing ||
                    ( ( Ind->genotype[ locusA ] != bg ) &&
                      ( Ind->genotype[ locusA + 1 ] != bg ) &&
                      ( Ind->genotype[ locusB ] != bg ) &&
                      ( Ind->genotype[ locusB + 1 ] != bg ) ) )
            {
                hap[ 0 ] = Ind->genotype[ locusA ];
                hap[ 1 ] = Ind->genotype[ locusB ];
                startNode = insertIntLinkList( startNode, hap, 2, pGt );
                hap[ 0 ] = Ind->genotype[ locusA + 1 ];
                hap[ 1 ] = Ind->genotype[ locusB + 1 ];
                ( *_noI ) ++;
                nrGeno += 4;

                if ( 4 == nrGeno )
		  genotype = ( int* ) yCalloc( nrGeno, int );
                else
		  genotype = ( int* ) yRealloc( genotype, nrGeno, int );

                genotype[ nrGeno - 4 ] = Ind->genotype[ locusA ];
                genotype[ nrGeno - 3 ] = Ind->genotype[ locusB ];
                genotype[ nrGeno - 2 ] = Ind->genotype[ locusA + 1 ];
                genotype[ nrGeno - 1 ] = Ind->genotype[ locusB + 1 ];

                startNode = insertIntLinkList( startNode, hap, 2, pGt );
            }
        }
    }


    {
      IntNode*Node = startNode;

      if ( ( NULL != Node ) && ( NULL != ( Node = go_down( Node ) ) ) )
        {
	  int ttt = 0;
	  int x, xPos = 0, yPos = 0;
	  int*path = ( int* ) yCalloc( nGt * 3, int );
	  int lx, ly;
	  
	  for ( x = 0;x < nGt*3;x++ )
            {
	      path[ x ] = 0;
            }
	  
	  do
	    {
	      print_path( Node, &nLoci, &nGt, &path[ ttt++ ] );
	    }
	  while ( NULL != ( Node = find_next_leaf( Node ) ) );

          startNode = NULL;

	  for ( x = 0;x < nGt;x++ )
            {
	      if ( ( 0 == dimLA ) && ( 0 == dimLB ) )
                {
		  dimLA = 1;
		  dimLB = 1;
		  table =  yCalloc( dimLA + 1, int*);
		  table[ 0 ] =  yCalloc( dimLB + 1, int );
		  table[ 1 ] =  yCalloc( dimLB + 1, int );
		  table[ 0 ][ 0 ] = path[ x + 2 * nGt ];
		  table[ 1 ][ 0 ] = path[ x ];
		  table[ 0 ][ 1 ] = path[ x + nGt ];
		  table[ 1 ][ 1 ] = path[ x + 2 * nGt ];
                }
	      else
                {
		  for ( lx = 1, xPos = 0;lx <= dimLB;lx++ )
                    {
		      if ( table[ 0 ][ lx ] == path[ x + nGt ] )
                        {
			  xPos = lx;
                        }
                    }
		  if ( 0 == xPos )
                    {
		      dimLB++;
		      
		      for ( ly = 0;ly <= dimLA;ly++ )
                        {
			  table[ ly ] = ( int* ) yRealloc( table[ ly ],
							   dimLB + 1,
							   int );
			  table[ ly ][ dimLB ] = 0;
                        }
		      
		      lx = 1;
		      
		      while ( table[ 0 ][ lx ] < path[ x + nGt ] )
                        {
			  lx++;
			  
			  if ( lx >= dimLB )
			    break;
                        }
		      
		      xPos = lx;
		      
		      if ( lx != dimLB )
                        {
			  for ( ly = 0;ly <= dimLA;ly++ )
                            {
			      for ( lx = dimLB;lx > xPos;lx-- )
                                {
				  table[ ly ][ lx ] = table[ ly ][ lx - 1 ];
                                }
                            }
			  
			  for ( ly = 0;ly <= dimLA;ly++ )
                            {
			      table[ ly ][ xPos ] = 0;
                            }
                        }
		      
		      table[ 0 ][ xPos ] = path[ x + nGt ];
                    }
		  
		  for ( ly = 1, yPos = 0;ly <= dimLA;ly++ )
                    {
		      if ( table[ ly ][ 0 ] == path[ x ] )
                        {
			  yPos = ly;
                        }
                    }
		  
		  if ( 0 == yPos )
                    {
		      dimLA++;
		      table = ( int** ) yRealloc( table,
						  dimLA + 1,
						  int* );
		      table[ dimLA ] = ( int* ) yCalloc( dimLB + 1,
							 int );
		      
		      table[ dimLA ][ 0 ] = path[ x ];
		      
		      for ( lx = 1;lx <= dimLB;lx++ )
                        {
			  table[ dimLA ][ lx ] = 0;
                        }
		      
		      yPos = dimLA;
                    }
		  
		  if ( ( xPos != 0 ) && ( yPos != 0 ) )
                    {
		      table[ 0 ][ 0 ] += path[ x + 2 * nGt ];
		      table[ yPos ][ xPos ] = path[ x + 2 * nGt ];
                    }
                }
            }
	  
	  yFree( path );
        }
    }
    
    *_nodH = dH;
    *_dLocA = dimLA;
    *_dLocB = dimLB;

    {
        int maxAllelA = 0;
        int maxAllelB = 0;

        int*transAllelA;
        int*transAllelB;

        int i;

        for ( i = 1;i <= dimLA;i++ )
        {
            if ( maxAllelA < table[ i ][ 0 ] )
            {
                maxAllelA = table[ i ][ 0 ];
            }
        }

        for ( i = 1;i <= dimLB;i++ )
        {
            if ( maxAllelB < table[ 0 ][ i ] )
            {
                maxAllelB = table[ 0 ][ i ];
            }
        }

        transAllelA = ( int* ) yCalloc( maxAllelA + 1, int );
        transAllelB = ( int* ) yCalloc( maxAllelB + 1, int );

        for ( i = 0;i <= maxAllelA;i++ )
        {
            transAllelA[ i ] = 0;
        }

        for ( i = 0;i <= maxAllelB;i++ )
        {
            transAllelB[ i ] = 0;
        }

        for ( i = 1;i <= dimLA;i++ )
        {
            transAllelA[ table[ i ][ 0 ] ] = i - 1;
        }

        for ( i = 1;i <= dimLB;i++ )
        {
            transAllelB[ table[ 0 ][ i ] ] = ( i - 1 ) * dimLA;
        }


        for ( i = 0;i < nrGeno;i += 2 )
        {
            genotype[ i ] = transAllelA[ genotype[ i ] ];
            genotype[ i + 1 ] = transAllelB[ genotype[ i + 1 ] ];
        }

        yFree( transAllelA );
        yFree( transAllelB );
    }
}


SEXP haplotype_return_table(void)

{
  SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;

  PROTECT( matrix = allocMatrix( INTSXP, dimLA, dimLB ) );
  PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
  INTEGER( matrixDim ) [ 0 ] = dimLA;
  INTEGER( matrixDim ) [ 1 ] = dimLB;
  setAttrib( matrix, R_DimSymbol, matrixDim );
  PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
  PROTECT( xAxis = allocVector( STRSXP, dimLB ) );
  PROTECT( yAxis = allocVector( STRSXP, dimLA ) );
  
  {
    int tmpcounter, x, y;
    char tmp[ 200 ];
    
    for ( x = 1, tmpcounter = 0;x <= dimLB;x++ )
      {
	snprintf(tmp, sizeof(tmp),
		 "locusB.%d ", table[0][x]);
	SET_STRING_ELT( xAxis, x - 1, mkChar( tmp ) );
	
	for ( y = 1;y <= dimLA;y++, tmpcounter++ )
	  {
	    if ( x == 1 )
	      {
		snprintf(tmp, sizeof(tmp),
			 "locusA.%d ", table[y][0]);

		SET_STRING_ELT( yAxis, y - 1, mkChar( tmp ) );
	      }
	    
	    INTEGER( matrix ) [ tmpcounter ] = table[ y ][ x ];
	  }
      }
  }
  
  SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
  SET_VECTOR_ELT( matrixDimNames, 0, yAxis );
  setAttrib( matrix, R_DimNamesSymbol, matrixDimNames );
  
  UNPROTECT( 5 );
  return ( matrix );
}


SEXP haplotype_return_genotype(void)
  
{
  SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;
  
  PROTECT( matrix = allocMatrix( INTSXP, nrGeno / 2, 2 ) );
  PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
  INTEGER( matrixDim ) [ 0 ] = nrGeno / 2;
  INTEGER( matrixDim ) [ 1 ] = 2;
  setAttrib( matrix, R_DimSymbol, matrixDim );
  PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
  PROTECT( xAxis = allocVector( STRSXP, 2 ) );
  PROTECT( yAxis = allocVector( STRSXP, nrGeno / 2 ) );
  
  {
    int tmpcounter, x, y;
    char tmp[ 200 ];
    
    SET_STRING_ELT( xAxis, 0, mkChar( "LocusA" ) );
    SET_STRING_ELT( xAxis, 1, mkChar( "LocusB" ) );
    
    for ( x = 1, tmpcounter = 0;x <= ( nrGeno / 4 );x++ )
      {
	for ( y = 1;y <= 2;y++, tmpcounter++ )
	  {
	    snprintf(tmp, sizeof(tmp),
		     "Ind%d.%d", x, y);

	    SET_STRING_ELT( yAxis, 2 * ( x - 1 ) + ( y - 1 ), mkChar( tmp ) );
	    INTEGER( matrix ) [ tmpcounter ] = genotype[ 4 * ( x - 1 ) + 2 * ( y - 1 ) ];
	    INTEGER( matrix ) [ tmpcounter + nrGeno / 2 ] = genotype[ 4 * ( x - 1 ) + 2 * ( y - 1 ) + 1 ];
	  }
      }
  }
  
  SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
  SET_VECTOR_ELT( matrixDimNames, 0, yAxis );
  setAttrib( matrix, R_DimNamesSymbol, matrixDimNames );
  
  UNPROTECT( 5 );
  return ( matrix );
}


void haplotype_calc_d ( double *_dd,
		        double *_rr,
		        double *_dprime )

{
  int x, y, xx, yy, sumxx, sumyy, sumTotal;
  double D, DD = 0, RA = 1, RB = 1, dPrime = 0;
  double p, q;

  RA = 1;
  RB = 1;

  for ( x = 1, sumTotal = 0; x <= dimLB; x++ ) {
    for ( y = 1; y <= dimLA; y++ ) {
      sumTotal += table[ y ][ x ];
    }
  }

  for ( x = 1;x <= dimLB;x++ ) {

    for ( yy = 1, sumyy = 0;yy <= dimLA;yy++ ) {
      sumyy += table[ yy ][ x ];
    }
      
    for ( y = 1; y <= dimLA; y++ ) {

      for ( xx = 1, sumxx = 0;xx <= dimLB;xx++ ) {
	sumxx += table[ y ][ xx ];
      }

      D = ( ( (double) table[y][x] - ( (double) sumxx * (double) sumyy ) /
	      ( (double) sumTotal ) ) / ( (double) sumTotal ) );

      DD += D * D;
      
      if ( 1 == x ) {
	RA -= (double) ( ( (double) sumxx / ( (double) sumTotal ) ) *
			 ( (double) sumxx / ( (double) sumTotal ) ) );
      }
      
      if ( 1 == y ) {
	RB -= (double) ( ( (double) sumyy / ( (double) sumTotal ) ) *
			 ( (double) sumyy / ( (double) sumTotal ) ) );
      }
      
      p = (double) sumxx / ( (double) sumTotal );
      q = (double) sumyy / ( ( double ) sumTotal );

      if ( D < 0 ) {

	  if ( (p*q) < ( (1-p)*(1-q) ) ) {
	    dPrime += ( p * q * D * ( -1 ) / ( p * q ) );
	  }
	  else {
	    dPrime += ( p * q * D * ( -1 ) / ( ( 1 - p ) * ( 1 - q ) ) );
	  }

      }
      else {

          if ( ( (1-p)*q ) < (p*(1-q) ) ) {
	    dPrime += ( p * q * D / ( ( 1 - p ) * q ) );
	  }
	  else {
	    dPrime += ( p * q * D / ( p * ( 1 - q ) ) );
	  }

      }
    }
  }

  *_dd = DD;
  *_rr = DD / ( RA * RB );
  *_dprime = dPrime;
}


void haplotype_free(void)

{
    int ly;

    if ( NULL != table )
    {
        for ( ly = 0;ly <= dimLA;ly++ )
        {
            yFree( table[ ly ] );
        }

        yFree( table );
    }

    yFree( genotype );
    nrGeno = 0;
    dimLA = 0;
    dimLB = 0;
}


void concat_population( 
		       char**NameP1,
		       char**NameP2)
{
    TInd tmpInd;
    TInd *tmpP1Ind;
    POPULATION ind = 0, j, oldP1Alloc, newP1Alloc;
    TPopListEl*PopListEl;
    TPop*P1Pop = NULL;
    TPop*P2Pop = NULL;

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl ){
      info( -2, "Population not defined: %s!", *NameP2 );
      return ;
    }

    P2Pop = PopListEl->Pop;

    if ( 0 == strcmp( *NameP1, *NameP2 ) )
    {
        info( -2, "Source and destination population must be different\n" );
        return;
    }

    find_list_element( *NameP1, &PopListEl );

    if ( NULL == PopListEl )
    {
        P1Pop = new_population( NameP1, &ind );
	if ( NULL == P1Pop ) return ;
    }

    else
    {
        P1Pop = PopListEl->Pop;
    }

    ind = P1Pop->NoInds, j = 0, oldP1Alloc = P1Pop->NoIndsAlloc;
    newP1Alloc = oldP1Alloc + P2Pop->NoIndsAlloc;

    if ( NULL == P1Pop->Ind )
    {
        if ( 0 < newP1Alloc )
        {
            tmpP1Ind = ( TInd* ) malloc(
                ( size_t ) newP1Alloc * sizeof( *tmpP1Ind ) );

            if ( NULL == tmpP1Ind )
            {
                info( -2, "Memory allocation failed in concat_population\n" );
                return;
            }

            P1Pop->Ind = tmpP1Ind;
        }
    }

    else if ( oldP1Alloc != newP1Alloc )
    {
        tmpP1Ind = ( TInd* ) realloc(
            P1Pop->Ind,
            ( size_t ) newP1Alloc * sizeof( *tmpP1Ind ) );

        if ( NULL == tmpP1Ind )
        {
            info( -2, "Memory allocation failed in concat_population\n" );
            return;
        }

        P1Pop->Ind = tmpP1Ind;
    }

    remove_evaluate_population( NameP2 );
    remove_genotype_population( NameP2 );

    if ( NULL != PopListEl )
    {
        remove_evaluate_population( NameP1 );
        remove_genotype_population( NameP1 );
    }

    P1Pop->NoInds += P2Pop->NoInds;
    P1Pop->NoIndsAlloc = newP1Alloc;

    for ( ;j < P2Pop->NoIndsAlloc;ind++, j++ )
    {
        tmpInd = P1Pop->Ind[ ind ];
        P1Pop->Ind[ ind ] = P2Pop->Ind[ j ];

        if ( ind < oldP1Alloc )
            P1Pop->Ind[ P1Pop->NoIndsAlloc - j - 1 ] = tmpInd;
    }

    P2Pop->NoIndsAlloc = 0;
    remove_population( NameP2 );

    info( 1, "Population %s concatenated to population %s\n",
          *NameP2, *NameP1 );
}


void sample_population( 
                       char**NameP1,
                       char**NameP2,
                       char**cNoI,
                       int*rep)
{
    int chrom, hom, loc, r;
    int copy_genotype, copy_gvalue, copy_pvalue;
    int genotype_alloc_error = 0, gvalue_alloc_error = 0;
    int pvalue_alloc_error = 0, info_alloc_error = 0;
    POPULATION ind, NoI = atol( *cNoI ), *samples = NULL;
    POPULATION *tmp_samples;
    TInd*P1Ind, *P2Ind;
    TPopListEl*PopListEl;

    TPop*P1Pop = NULL, *P2Pop = NULL;

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl )
      {
        info( -2, "Population not defined: %s \n", *NameP2 );
        return ;
      }

    P2Pop = PopListEl->Pop;

    if ( 0 == P2Pop->NoInds ){
        info( -2, "Population >%s< is empty\n", *NameP2 );
        return;
      }

    if ( 0 == strcmp( *NameP1, *NameP2 ) )
      {
        info( -2, "Source and destination population must be different\n" );
        return;
      }

    if ( -1 == NoI )
        NoI = P2Pop->NoInds;

    if ( 0 > NoI )
      {
        info( -2, "Incorrect number of individuals\n" );
        return;
      }

    if ( ( 0 == *rep ) && ( P2Pop->NoInds < NoI ) )
      {
        info( -2, "sampling only with repetition possible\n" );
        return;
      }

    copy_genotype = st_population_all_genotyped( P2Pop );
    copy_gvalue = st_population_all_gvalue( P2Pop );
    copy_pvalue = st_population_all_pvalue( P2Pop );

    P1Pop = exist_population( *NameP1, NoI, 1 );
    if ( NULL == P1Pop )  return;

    if ( NoI != P1Pop->NoIndsAlloc ) return;

    P1Pop->NoInds = NoI;

    for ( ind = 0;ind < P1Pop->NoInds;ind++ )
    {

      r = (int) gsl_rng_uniform_int( rng[0],
                                     (unsigned long int)P2Pop->NoInds );


        if ( 0 == *rep )
        {
            int j = 0;
            POPULATION ind_i;

            do
            {
                j = 0;

                for ( ind_i = 0;ind_i < ind;ind_i++ )
                {
                    if ( r == *( samples + ind_i ) )
                    {
                        j = 1;
                        break;
                    }
                }

                if ( 0 == j )
                {
                    if ( 0 == ind )
                      tmp_samples = ( POPULATION* ) malloc(
                        sizeof( *tmp_samples ) );
                    else
                      tmp_samples = ( POPULATION* ) realloc(
                        samples,
                        ( size_t )( ind + 1 ) * sizeof( *tmp_samples ) );

                    if ( NULL == tmp_samples )
                      {
                        info( -2,
                              "Memory allocation failed in sample_population\n" );
                        goto fail_population;
                      }

                    samples = tmp_samples;
                    *( samples + ind ) = r;
                }

                else
                  r = (int) gsl_rng_uniform_int(
                                                    rng[0],
                                                    (unsigned long int)
                                                    P2Pop->NoInds );
            }

            while ( 1 == j );
        }

        P1Ind = get_individual( P1Pop, ind );
        P2Ind = get_individual( P2Pop, r );

        if ( copy_genotype )
        {
            int i;

            P1Ind->genotype = ( int* ) yCalloc( cMPhom, int );

            if ( NULL == P1Ind->genotype )
              {
                genotype_alloc_error = 1;
                copy_genotype = 0;
              }
            else
              {
                for ( i = 0;i < cMPhom;i++ )
                  P1Ind->genotype[ i ] = P2Ind->genotype[ i ];
              }
        }

        else
        {
            P1Ind->genotype = NULL;
        }

        if ( !st_copy_value_array(
               &(P1Ind->GValue),
               copy_gvalue ? P2Ind->GValue : NULL ) )
          {
            gvalue_alloc_error = 1;
            copy_gvalue = 0;
          }

        if ( !st_copy_value_array(
               &(P1Ind->PValue),
               copy_pvalue ? P2Ind->PValue : NULL ) )
          {
            pvalue_alloc_error = 1;
            copy_pvalue = 0;
          }

        if ( NULL != P1Ind->info )
            yFree( P1Ind->info );

        if ( NULL != P2Ind->info )
        {
            P1Ind->info = ( char* ) malloc( MAXNAME * sizeof( char ) );

            if ( NULL == P1Ind->info )
              info_alloc_error = 1;
            else
              strcpy( P1Ind->info, P2Ind->info );
        }

        else
        {
            P1Ind->info = NULL;
        }

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            for ( hom = 0;hom < NoHoms;hom++ )
            {
                THom*P1Hom = &( P1Ind->Chrom[ chrom ] ).Hom[ hom ];
                THom*P2Hom = &( P2Ind->Chrom[ chrom ] ).Hom[ hom ];
                P1Hom->NoLociUsed = 0;

                for ( loc = 0;loc < P2Hom->NoLociUsed;loc++ )
                {
                    if ( -1 == insertLocus(
                           P2Hom->Loc[ loc ].All,
                           P2Hom->Loc[ loc ].Pos,
                           P1Hom ) )
                      goto fail_population;
                }
            }
        }

    }

    if ( genotype_alloc_error )
      st_clear_population_genotype( P1Pop );

    if ( gvalue_alloc_error )
      st_clear_population_gvalue( P1Pop );

    if ( pvalue_alloc_error )
      st_clear_population_pvalue( P1Pop );

    free( samples );

    if ( genotype_alloc_error ||
         gvalue_alloc_error ||
         pvalue_alloc_error ||
         info_alloc_error )
      info( -2, "Memory allocation failed in sample_population\n" );

    return;

fail_population:

    free( samples );

    if ( NULL != P1Pop )
      remove_population( NameP1 );

    return;
}


void copy_population( 
                     char**NameP1,
                     char**NameP2,
                     char**cstart,
                     char**cn)

{
    int chrom, hom, loc;
    int copy_genotype, copy_gvalue, copy_pvalue;
    int genotype_alloc_error = 0, gvalue_alloc_error = 0;
    int pvalue_alloc_error = 0, info_alloc_error = 0;
    POPULATION ind, start = atol( *cstart ), n = atol( *cn );
    TInd*P1Ind, *P2Ind;
    TPop*P1Pop = NULL, *P2Pop = NULL;
    TPopListEl*PopListEl;

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl )
      {
        info( -2, "Population not defined: %s !\n", *NameP2 );
        return;
      }

    P2Pop = PopListEl->Pop;

    if ( 0 == P2Pop->NoInds )
      {
        info( -2, "Population >%s< is empty\n", *NameP2 );
        return;
      }

    if ( 0 == strcmp( *NameP1, *NameP2 ) )
      {
        info( -2, "Source and destination population must be different\n" );
        return;
      }

    start--;

    if ( 0 > start )
        start = 0;

    if ( P2Pop->NoInds <= start )
      {
        info( -2, "Start index outside population %s\n", *NameP2 );
        return;
      }

    if ( 0 > n )
        n = P2Pop->NoInds;

    if ( P2Pop->NoInds < start + n )
        n = P2Pop->NoInds - start;

    copy_genotype = st_population_all_genotyped( P2Pop );
    copy_gvalue = st_population_all_gvalue( P2Pop );
    copy_pvalue = st_population_all_pvalue( P2Pop );

    P1Pop = exist_population( *NameP1, n, 1 );
    if ( NULL == P1Pop ) return;

    if ( n != P1Pop->NoIndsAlloc ) return;

    P1Pop->NoInds = n;

    for ( ind = 0;ind < n;ind++ )
    {
        P1Ind = get_individual( P1Pop, ind );
        P2Ind = get_individual( P2Pop, start + ind );


        if ( copy_genotype )
        {
            int i;

            P1Ind->genotype = ( int* ) yCalloc( cMPhom, int );

            if ( NULL == P1Ind->genotype )
              {
                genotype_alloc_error = 1;
                copy_genotype = 0;
              }
            else
              {
                for ( i = 0;i < cMPhom;i++ )
                  P1Ind->genotype[ i ] = P2Ind->genotype[ i ];
              }
        }

        else
        {
            P1Ind->genotype = NULL;
        }

        if ( !st_copy_value_array(
               &(P1Ind->GValue),
               copy_gvalue ? P2Ind->GValue : NULL ) )
          {
            gvalue_alloc_error = 1;
            copy_gvalue = 0;
          }

        if ( !st_copy_value_array(
               &(P1Ind->PValue),
               copy_pvalue ? P2Ind->PValue : NULL ) )
          {
            pvalue_alloc_error = 1;
            copy_pvalue = 0;
          }

        if ( NULL != P1Ind->info )
            yFree( P1Ind->info );

        if ( NULL != P2Ind->info )
        {
            P1Ind->info = ( char* ) malloc( MAXNAME * sizeof( char ) );

            if ( NULL == P1Ind->info )
              info_alloc_error = 1;
            else
              strcpy( P1Ind->info, P2Ind->info );
        }

        else
        {
            P1Ind->info = NULL;
        }

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            for ( hom = 0;hom < NoHoms;hom++ )
            {
                THom*P1Hom = &( P1Ind->Chrom[ chrom ] ).Hom[ hom ];
                THom*P2Hom = &( P2Ind->Chrom[ chrom ] ).Hom[ hom ];
                P1Hom->NoLociUsed = 0;

                for ( loc = 0;loc < P2Hom->NoLociUsed;loc++ )
                {
                    if ( -1 == insertLocus(
                           P2Hom->Loc[ loc ].All,
                           P2Hom->Loc[ loc ].Pos,
                           P1Hom ) )
                      goto fail_population;
                }
            }
        }

    }

    if ( genotype_alloc_error )
      st_clear_population_genotype( P1Pop );

    if ( gvalue_alloc_error )
      st_clear_population_gvalue( P1Pop );

    if ( pvalue_alloc_error )
      st_clear_population_pvalue( P1Pop );

    if ( genotype_alloc_error ||
         gvalue_alloc_error ||
         pvalue_alloc_error ||
         info_alloc_error )
      {
        info( -2, "Memory allocation failed in copy_population\n" );
        return;
      }

    info( 1, "Copied %ld individuals (Index %ld to %ld) from population %s to population %s.\n",
          n, start + 1, start + n, *NameP2, *NameP1 );
    return;

fail_population:

    if ( NULL != P1Pop )
      remove_population( NameP1 );

    return;
}


void append_population(
                       char**NameP1,
                       char**NameP2)
{
    int chrom, hom, loc;
    int copy_genotype = 0, copy_gvalue = 0, copy_pvalue = 0;
    int genotype_alloc_error = 0, gvalue_alloc_error = 0;
    int pvalue_alloc_error = 0, info_alloc_error = 0;
    int source_genotyped, source_gvalue, source_pvalue;
    int destination_existed = 0;
    POPULATION P1Inds = 0, P1Alloc = 0, NoIndsNew, ind;
    TInd*P1Ind, *P2Ind;
    TPopListEl*PopListEl;
    TPop*P1Pop = NULL, *P2Pop = NULL;

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Population not defined: %s \n", *NameP2 );
        return;
    }

    P2Pop = PopListEl->Pop;

    if ( 0 == P2Pop->NoInds ) {
        info( -2, "Population >%s< is empty\n", *NameP2 );
        return;
    }

    if ( 0 == strcmp( *NameP1, *NameP2 ) ) {
        info( -2, "Source and destination population must be different\n" );
        return;
    }

    source_genotyped = st_population_all_genotyped( P2Pop );
    source_gvalue = st_population_all_gvalue( P2Pop );
    source_pvalue = st_population_all_pvalue( P2Pop );

    find_list_element( *NameP1, &PopListEl );

    if ( NULL != PopListEl )
      {
        destination_existed = 1;
        P1Inds = PopListEl->Pop->NoInds;
        P1Alloc = PopListEl->Pop->NoIndsAlloc;
      }

    if ( P2Pop->NoInds > LONG_MAX - P1Inds )
      {
        info( -2, "Resulting population is too large\n" );
        return;
      }

    NoIndsNew = P1Inds + P2Pop->NoInds;

    P1Pop = exist_population( *NameP1, NoIndsNew, 0 );
    if ( NULL == P1Pop ) return;

    if ( NoIndsNew != P1Pop->NoIndsAlloc ) {
        info( -2, "Population %s could not be resized\n", *NameP1 );
        return;
    }

    if ( 0 == P1Inds )
      {
        copy_gvalue = source_gvalue;
        copy_pvalue = source_pvalue;
        copy_genotype = source_genotyped;
      }
    else
      {
        copy_gvalue =
          st_population_all_gvalue( P1Pop ) && source_gvalue;

        copy_pvalue =
          st_population_all_pvalue( P1Pop ) && source_pvalue;

        copy_genotype =
          st_population_all_genotyped( P1Pop ) && source_genotyped;
      }

    for ( ind = 0;ind < P2Pop->NoInds;ind++ )
    {
        P1Ind = get_individual( P1Pop, ind + P1Inds );
        P2Ind = get_individual( P2Pop, ind );

        if ( NULL != P1Ind->genotype )
            yFree( P1Ind->genotype );

        if ( NULL != P1Ind->GValue )
            yFree( P1Ind->GValue );

        if ( NULL != P1Ind->PValue )
            yFree( P1Ind->PValue );

        if ( NULL != P1Ind->info )
            yFree( P1Ind->info );

        if ( copy_genotype )
        {
            int i;

            P1Ind->genotype = ( int* ) yCalloc( cMPhom, int );

            if ( NULL == P1Ind->genotype )
              {
                genotype_alloc_error = 1;
                copy_genotype = 0;
              }
            else
              {
                for ( i = 0;i < cMPhom;i++ )
                  P1Ind->genotype[ i ] = P2Ind->genotype[ i ];
              }
        }

        if ( !st_copy_value_array(
               &(P1Ind->GValue),
               copy_gvalue ? P2Ind->GValue : NULL ) )
          {
            gvalue_alloc_error = 1;
            copy_gvalue = 0;
          }

        if ( !st_copy_value_array(
               &(P1Ind->PValue),
               copy_pvalue ? P2Ind->PValue : NULL ) )
          {
            pvalue_alloc_error = 1;
            copy_pvalue = 0;
          }

        if ( NULL != P2Ind->info )
        {
            P1Ind->info = ( char* ) malloc( MAXNAME * sizeof( char ) );

            if ( NULL == P1Ind->info )
              info_alloc_error = 1;
            else
              strcpy( P1Ind->info, P2Ind->info );
        }

        else
        {
            P1Ind->info = NULL;
        }

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            for ( hom = 0;hom < NoHoms;hom++ )
            {
                THom*P1Hom = &( P1Ind->Chrom[ chrom ] ).Hom[ hom ];
                THom*P2Hom = &( P2Ind->Chrom[ chrom ] ).Hom[ hom ];
                P1Hom->NoLociUsed = 0;

                for ( loc = 0;loc < P2Hom->NoLociUsed;loc++ )
                {
                    if ( -1 == insertLocus(
                           P2Hom->Loc[ loc ].All,
                           P2Hom->Loc[ loc ].Pos,
                           P1Hom ) )
                      goto fail_population;
                }
            }
        }

    }

    P1Pop->NoInds = NoIndsNew;

    if ( !copy_genotype || genotype_alloc_error )
      st_clear_population_genotype( P1Pop );

    if ( !copy_gvalue || gvalue_alloc_error )
      st_clear_population_gvalue( P1Pop );

    if ( !copy_pvalue || pvalue_alloc_error )
      st_clear_population_pvalue( P1Pop );

    if ( genotype_alloc_error ||
         gvalue_alloc_error ||
         pvalue_alloc_error ||
         info_alloc_error )
      {
        info( -2, "Memory allocation failed in append_population\n" );
        return;
      }

    info( 1, "Population %s (%ld individuals) appended to population %s.\n",
          *NameP2, P2Pop->NoInds, *NameP1 );
    return;

fail_population:

    for ( ind = P1Inds;
          ( ind < NoIndsNew ) && ( ind < P1Pop->NoIndsAlloc );
          ind++ )
      {
        P1Ind = get_individual( P1Pop, ind );

        free( P1Ind->genotype );
        P1Ind->genotype = NULL;
        free( P1Ind->GValue );
        P1Ind->GValue = NULL;
        free( P1Ind->PValue );
        P1Ind->PValue = NULL;
        free( P1Ind->info );
        P1Ind->info = NULL;

        for ( chrom = 0;chrom < NoChroms;chrom++ )
          for ( hom = 0;hom < NoHoms;hom++ )
            P1Ind->Chrom[ chrom ].Hom[ hom ].NoLociUsed = 0;
      }

    P1Pop->NoInds = P1Inds;

    if ( !destination_existed )
      remove_population( NameP1 );
    else if ( P1Pop->NoIndsAlloc != P1Alloc )
      int_resize_population( NameP1, P1Alloc );

    return;
}


void select_all_n_best(
    char **NamePop,
    char **indexbest,
    int *n,
    double *minscore,
    double *maxscore,
    char **criterion_name,
    char **EffName,
    int *retval)
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

  if (!st_effect_position(EffName,&eff))
    return;

  find_list_element(*NamePop,&PopListEl);

  if (NULL == PopListEl)
    {
      info(-2,"Population not defined: %s\n",*NamePop);
      return;
    }

  Pop = PopListEl->Pop;

  if (0 == Pop->NoInds)
    {
      info(-2,"Population %s is empty\n",*NamePop);
      return;
    }

  if (1 > _n)
    {
      info(-2,
           "Number of individuals to select must be at least 1\n");
      return;
    }

  if (!st_population_selection_ready(Pop,criterion,eff))
    {
      info(-2,
           (0 == criterion) ?
           "Population not evaluated\n" :
           "Population not phenotyped\n");
      return;
    }

  Ind0 = get_individual(Pop,ind);
  value0 = st_individual_value_array(Ind0,criterion);
  *maxscore = value0[eff];

  for (ind = 1; ind < Pop->NoInds; ind++)
    {
      Ind1 = get_individual(Pop,ind);
      value1 = st_individual_value_array(Ind1,criterion);

      if (0.0f != (value0[eff] - value1[eff]))
        {
          if (0 == --_n)
            break;
        }

      Ind0 = Ind1;
      value0 = value1;
    }

  *minscore = value0[eff];

  snprintf(*indexbest,strlen(*indexbest),"%ld",ind);

  info(1,
       "%ld individuals selected, (score %f .. %f)\n",
       ind,*minscore,*maxscore);

  *retval = 0;
}


void divide_population ( char **NameP1,
			 char **NameP2,
			 char **cNoI   )
{
  POPULATION NoI = atol( *cNoI );
  TPopListEl *PopListEl;
  TPop *P1Pop, *P2Pop = NULL;

  find_list_element( *NameP2, &PopListEl );

  if ( PopListEl != NULL ) 
    {
      P2Pop = PopListEl->Pop;
      P1Pop = exist_population( *NameP1, NoI, 1 );
      if ( NULL == P1Pop )  return;

      if ( NoI > P2Pop->NoInds )
	{
	  info( -1, 
		"population.divide - population size adjusted %ld -> %ld \n", 
		NoI, 
		P2Pop->NoInds );

	  NoI = P2Pop->NoInds;
	}
      
      {
        TInd tmpInd;
        POPULATION ind, ind_j;

        P1Pop->NoInds = NoI;

        for ( ind = 0; 
	      ind < P1Pop->NoInds; 
	      ind++ )
        {
            tmpInd = P1Pop->Ind[ ind ];
            P1Pop->Ind[ ind ] = P2Pop->Ind[ ind ];
            P2Pop->Ind[ ind ] = tmpInd;
        }

        for ( ind = NoI, ind_j = 0;
	      ind < P2Pop->NoInds;
	      ind++, ind_j++ )
        {
            tmpInd = P2Pop->Ind[ ind ];
            P2Pop->Ind[ ind ] = P2Pop->Ind[ ind_j ];
            P2Pop->Ind[ ind_j ] = tmpInd;
        }

        P2Pop->NoInds -= NoI;
      }
    }   else
    {
      info( -2, "Population %s is not defined.", *NameP2 );
      return;
    }

}

void rename_population(
		       char**OldName,
		       char**NewName)

{
    TPopListEl*PopListEl1, *PopListEl2;
    TPop*Pop = NULL;

    find_list_element( *OldName, &PopListEl1 );

    if ( PopListEl1 == NULL ){
        info( -2, "Population not defined: %s\n", *OldName );
	return;
    }

    find_list_element( *NewName, &PopListEl2 );

    if ( PopListEl2 != NULL ){
        info( -2, "Populationname >%s< is already used\n", *NewName );
	return;
    }

    if ( 2 > strlen( *NewName ) ){
        info( -2, "New populationname is too short\n" );
	return;
    }

    if ( strlen( *NewName ) >= MAXNAME ){
        info( -2, "New populationname is too long\n" );
	return;
    }

    Pop = PopListEl1->Pop;

    strcpy( Pop->Name, *NewName );
}


void swap_population_name( 
			  char**NameP1,
			  char**NameP2)

{
    TPopListEl*PopListEl;
    find_list_element( *NameP1, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Population not defined: %s\n", *NameP1 );
	return;
    }

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Population not defined:  %s\n", *NameP2 );
	return;
    }

    if ( 0 == strcmp( *NameP1, *NameP2 ) ){
        info( -2, "You have to specify two different populations\n" );
	return;
    }

    {

        {
            TPopListEl*PopListEl;
            char tmpName[ MAXNAME ];
            TPop*P1Pop = NULL, *P2Pop = NULL;

            find_list_element( *NameP1, &PopListEl );

            if ( PopListEl != NULL )
                P1Pop = PopListEl->Pop;

            if ( NULL == P1Pop ){
                info( -2, "Population not defined: %s\n", *NameP1 );
		return;
	    }

            find_list_element( *NameP2, &PopListEl );

            if ( PopListEl != NULL )
                P2Pop = PopListEl->Pop;

            if ( NULL == P2Pop ){
                info( -2, "Population not defined:  %s\n", *NameP2 );
		return;
	    }

            strcpy( tmpName, P1Pop->Name );

            strcpy( P1Pop->Name, P2Pop->Name );

            strcpy( P2Pop->Name, tmpName );
        }

    }
}

void get_population_size( 
			 char**PopNames,
			 char**cNoInds,
			 char**cNoIndsAlloc)

{
    TPopListEl*PopListEl;
    POPULATION NoInds = 0, NoIndsAlloc = 0;
    char*PopName;

    PopName = strtok( *PopNames, " " );

    while ( NULL != PopName )
    {
        find_list_element( PopName, &PopListEl );

        if ( NULL == PopListEl ){
	  info( -2, "Population not defined: %s\n", PopName );
	  return;
	}

        NoInds += PopListEl->Pop->NoInds;

        NoIndsAlloc += PopListEl->Pop->NoIndsAlloc;

        PopName = strtok( NULL, " " );
    }

    snprintf(*cNoInds, strlen(*cNoInds),
	     "%ld", NoInds);
    snprintf(*cNoIndsAlloc, strlen(*cNoIndsAlloc),
	     "%ld", NoIndsAlloc);

}


void get_genome_par( int*_NoChroms, int*_NoHoms )

{
    *_NoChroms = NoChroms;
    *_NoHoms = NoHoms;
}

void get_genome_par_b( double*_ChromLen )

{
    int chrom;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        _ChromLen[ chrom ] = PopDescr[ chrom ].ChromLen;
    }

}


void get_map_a(
    int*_NoMappoints,
    int*_maxname,
    char**_classlgth,
    char**_namelgth
)
{
    int chrom;
    int point;
    unsigned long int class_lgth = 0;
    unsigned long int name_lgth = 0;

    *_NoMappoints = 0;
    *_maxname = MAXNAME;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        *_NoMappoints += PopDescr[ chrom ].NrMappoints;

        for ( point = 0;point < PopDescr[ chrom ].NrMappoints;point++ )
        {
            name_lgth += ( 1 + strlen( PopDescr[ chrom ].Mappoint[ point ].name ) );
            class_lgth += ( 1 + strlen( PopDescr[ chrom ].Mappoint[ point ].class ) );
        }
    }

    snprintf(*_classlgth, strlen(*_classlgth),
         "%lu", class_lgth);

    snprintf(*_namelgth, strlen(*_namelgth),
         "%lu", name_lgth);


}

void get_map_b(
    double*_chrom,
    double*_pos,
    char**_name,
    char**_class
)
{
    int chrom;
    int point;
    long counter = 0;
    char*p = *_name;
    char*q = *_class;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        for ( point = 0;point < PopDescr[ chrom ].NrMappoints;point++ )
        {
            size_t len;

            _chrom[ counter ] = chrom + 1;
            _pos[ counter ] = PopDescr[ chrom ].Mappoint[ point ].Pos;

            len = strlen( PopDescr[ chrom ].Mappoint[ point ].name );
            memcpy( p, PopDescr[ chrom ].Mappoint[ point ].name, len );
            p += len;
            *p++ = ' ';

            len = strlen( PopDescr[ chrom ].Mappoint[ point ].class );
            memcpy( q, PopDescr[ chrom ].Mappoint[ point ].class, len );
            q += len;
            *q++ = ' ';

            counter++;
        }
    }
}


TMappoint*find_mappoint( char*pch )
{
    int chrom, i;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        if ( PopDescr[ chrom ].NrMappoints > 0 )
        {
            for ( i = 0;i < PopDescr[ chrom ].NrMappoints;i++ )
            {
                if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ i ].name, pch ) )
                    return & ( PopDescr[ chrom ].Mappoint[ i ] );
            }
        }
    }

    return NULL;
}


int insertmember( 
		  TSet*Set,
		  int newmember)
{
    int i;
    int*tmp;
    size_t count;

    for ( i = 0;i < Set->betrag;i++ )
    {
        if ( Set->member[ i ] == newmember )
            return 1;
    }

    if ( ( Set->betrag < 0 ) || ( Set->betrag == INT_MAX ) )
        return 0;

    count = ( size_t ) Set->betrag + 1;

    if ( count > SIZE_MAX / sizeof( *Set->member ) )
        return 0;

    tmp = ( int* ) realloc(
		      Set->member,
		      count * sizeof( *Set->member )
		      );

    if ( NULL == tmp )
        return 0;

    Set->member = tmp;
    Set->member[ Set->betrag ] = newmember;
    Set->betrag++;

    return 1;
}

TEffects* find_Effect( char*Name )
{
    int i;

    for ( i = 0; i < NrEffects; i++ )
    {
        if ( 0 == strcmp( Name, Effects[ i ].Name ) )
            return & ( Effects[ i ] );
    }

    return NULL;
}


void list_eff_parameter( int* lEffName, int* nrEffects )
{
    int i;
    *nrEffects = NrEffects;

    for ( i = 0; i < NrEffects; i++ )
    {
        *lEffName += strlen( Effects[ i ].Name ) + 1;
    }

    ( *lEffName ) += 10;
}



void list_eff( char**all_eff, double*weigth_eff )
{
    int i;
    strcpy( *all_eff, "" );

    for ( i = 0;i < NrEffects;i++ )
    {
        strcat( *all_eff, Effects[ i ].Name );
        strcat( *all_eff, " " );
        weigth_eff[ i ] = Effects[ i ].weight;
    }
}



static int st_effect_position(char **EffName, int *eff)
{
  TEffects *Effect;

  if (NULL == eff)
    return 0;

  *eff = 0;

  if ((NULL == EffName) ||
      (NULL == *EffName) ||
      ('\0' == **EffName))
    return 1;

  Effect = find_Effect(*EffName);

  if (NULL == Effect)
    {
      info(-2,"Effect %s not defined\n",*EffName);
      return 0;
    }

  *eff = Effect->Pos;

  return 1;
}


static int st_selection_criterion(char **criterion)
{
  if ((NULL == criterion) ||
      (NULL == *criterion) ||
      (0 == strcmp(*criterion,"G")))
    return 0;

  if (0 == strcmp(*criterion,"P"))
    return 1;

  info(-2,"Selection criterion must be G or P\n");
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


static int st_population_selection_ready(
    TPop *Pop,
    int criterion,
    int eff)
{
  POPULATION ind;

  if ((NULL == Pop) ||
      (0 > eff) ||
      (eff > NrEffects))
    return 0;

  for (ind = 0; ind < Pop->NoInds; ind++)
    {
      float *value =
        st_individual_value_array(&(Pop->Ind[ind]),criterion);

      if ((NULL == value) ||
          !isfinite((double)value[eff]))
        return 0;
    }

  return 1;
}


void get_population_value(
    char **PopName,
    double *value,
    char **EffName,
    char **criterion_name,
    int *retval)
{
  int criterion;
  int eff;
  POPULATION ind;
  TPopListEl *PopListEl;

  *retval = -2;

  criterion = st_selection_criterion(criterion_name);

  if (0 > criterion)
    return;

  if (!st_effect_position(EffName,&eff))
    return;

  find_list_element(*PopName,&PopListEl);

  if (NULL == PopListEl)
    {
      info(-2,"Population %s not defined\n",*PopName);
      return;
    }

  if (!st_population_selection_ready(
         PopListEl->Pop,criterion,eff))
    {
      info(-2,
           (0 == criterion) ?
           "Population not evaluated\n" :
           "Population not phenotyped\n");
      return;
    }

  for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
      TInd *Ind = get_individual(PopListEl->Pop,ind);
      float *v = st_individual_value_array(Ind,criterion);

      value[ind] = (double)v[eff];
    }

  *retval = 0;
}


void set_population_gvalue(
    char **PopName,
    double *gvalue,
    int *ngvalue,
    int *retval)
{
  POPULATION ind;
  TPopListEl *PopListEl;
  float **new_value = NULL;
  int need_allocate = 0;

  *retval = -2;

  find_list_element(*PopName,&PopListEl);

  if (NULL == PopListEl)
    {
      info(-2,"Population %s not defined\n",*PopName);
      return;
    }

  if (0 == PopListEl->Pop->NoInds)
    {
      info(-2,"Population %s is empty\n",*PopName);
      return;
    }

  if ((0 > *ngvalue) ||
      ((POPULATION)*ngvalue != PopListEl->Pop->NoInds))
    {
      info(-2,
           "Number of genotypic values does not match population size\n");
      return;
    }

  for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
      if (!isfinite(gvalue[ind]) ||
          (fabs(gvalue[ind]) > FLT_MAX))
        {
          info(-2,"Genotypic values must be finite float values\n");
          return;
        }
    }

  if (!st_population_all_gvalue(PopListEl->Pop))
    {
      if (0 != NrEffects)
        {
          info(-2,"First evaluate population %s\n",*PopName);
          return;
        }

      need_allocate = 1;
      new_value = (float **)calloc(
        (size_t)PopListEl->Pop->NoInds,
        sizeof(*new_value));

      if (NULL == new_value)
        goto fail_alloc;

      for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
          new_value[ind] = (float *)malloc(sizeof(**new_value));

          if (NULL == new_value[ind])
            goto fail_alloc;

          new_value[ind][0] = (float)gvalue[ind];
        }
    }

  if (need_allocate)
    {
      for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
          TInd *Ind = get_individual(PopListEl->Pop,ind);

          free(Ind->GValue);
          Ind->GValue = new_value[ind];
          new_value[ind] = NULL;
        }
    }
  else
    {
      for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
          TInd *Ind = get_individual(PopListEl->Pop,ind);

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
        free(new_value[ind]);

      free(new_value);
    }

  info(-2,"Memory allocation failed in set_population_gvalue\n");
}


void phenotype_population(
    char **PopName,
    double *hsq,
    int *nhsq,
    int *retval)
{
  POPULATION ind;
  TPopListEl *PopListEl;
  float **new_value = NULL;
  double *mean = NULL;
  double *variance = NULL;
  double *sd = NULL;
  int i;

  *retval = -2;

  find_list_element(*PopName,&PopListEl);

  if (NULL == PopListEl)
    {
      info(-2,"Population %s not defined\n",*PopName);
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
      info(-2,"No effects defined\n");
      return;
    }

  if (*nhsq != NrEffects)
    {
      info(-2,
           "One heritability must be supplied for each effect\n");
      return;
    }

  if (!st_population_all_gvalue(PopListEl->Pop))
    {
      info(-2,"First evaluate population %s\n",*PopName);
      return;
    }

  for (i = 0; i < NrEffects; i++)
    {
      if (!isfinite(hsq[i]) ||
          (0.0 >= hsq[i]) ||
          (1.0 < hsq[i]))
        {
          info(-2,"Heritabilities must be in (0,1]\n");
          return;
        }
    }

  mean = (double *)calloc((size_t)NrEffects,sizeof(*mean));
  variance = (double *)calloc(
    (size_t)NrEffects,sizeof(*variance));
  sd = (double *)calloc((size_t)NrEffects,sizeof(*sd));
  new_value = (float **)calloc(
    (size_t)PopListEl->Pop->NoInds,
    sizeof(*new_value));

  if ((NULL == mean) ||
      (NULL == variance) ||
      (NULL == sd) ||
      (NULL == new_value))
    goto fail_alloc;

  for (i = 1; i < NrEffects + 1; i++)
    {
      double ss = 0.0;

      for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
          TInd *Ind = get_individual(PopListEl->Pop,ind);

          if (!isfinite((double)Ind->GValue[i]))
            {
              info(-2,"Genotypic values must be finite\n");
              goto fail;
            }

          mean[i-1] += (double)Ind->GValue[i];
        }

      mean[i-1] /= (double)PopListEl->Pop->NoInds;

      for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        {
          TInd *Ind = get_individual(PopListEl->Pop,ind);
          double d = (double)Ind->GValue[i] - mean[i-1];

          ss += d * d;
        }

      variance[i-1] =
        ss / (double)(PopListEl->Pop->NoInds - 1);

      if (!isfinite(variance[i-1]) ||
          (0.0 > variance[i-1]))
        {
          info(-2,"Invalid genotypic variance\n");
          goto fail;
        }

      sd[i-1] = sqrt(
        variance[i-1] / hsq[i-1] - variance[i-1]);

      if (!isfinite(sd[i-1]))
        {
          info(-2,"Invalid residual standard deviation\n");
          goto fail;
        }
    }

  for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
      if ((size_t)(NrEffects + 1) >
          SIZE_MAX / sizeof(**new_value))
        goto fail_alloc;

      new_value[ind] = (float *)malloc(
        (size_t)(NrEffects + 1) * sizeof(**new_value));

      if (NULL == new_value[ind])
        goto fail_alloc;
    }

  GetRNGstate();

  for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
      TInd *Ind = get_individual(PopListEl->Pop,ind);
      double index = 0.0;

      for (i = 1; i < NrEffects + 1; i++)
        {
          double p =
            (double)Ind->GValue[i] +
            norm_rand() * sd[i-1];

          if (!isfinite(p) || (fabs(p) > FLT_MAX))
            {
              PutRNGstate();
              info(-2,"Generated phenotypic value is not finite\n");
              goto fail;
            }

          new_value[ind][i] = (float)p;
          index += p * (double)Effects[i-1].weight;
        }

      if (!isfinite(index) || (fabs(index) > FLT_MAX))
        {
          PutRNGstate();
          info(-2,"Phenotypic selection index is not finite\n");
          goto fail;
        }

      new_value[ind][0] = (float)index;
    }

  PutRNGstate();

  for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
    {
      TInd *Ind = get_individual(PopListEl->Pop,ind);

      free(Ind->PValue);
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

  info(-2,"Memory allocation failed in phenotype_population\n");

fail:

  if (NULL != new_value)
    {
      for (ind = 0; ind < PopListEl->Pop->NoInds; ind++)
        free(new_value[ind]);
    }

  free(new_value);
  free(mean);
  free(variance);
  free(sd);
}


void set_population_info ( char **PopName,
			   char **_ind,
			   char **_info )
{
  POPULATION  ind = atol( *_ind ) - 1;
  TPopListEl *PopListEl;
  TInd       *Ind;
  
  find_list_element( *PopName, &PopListEl );

  if ( NULL == PopListEl ) 
      {
	info( -2, "Population %s not defined\n", *PopName );
	return;
      }
  
  Ind = ( TInd* ) get_individual( PopListEl->Pop, ind );
  
  if ( NULL == Ind->info )
    Ind->info = ( char* ) yCalloc( MAXNAME, char );
  
  if ( strlen( *_info ) >= MAXNAME )
    {
      info( -2, "Info string is too long\n" );
      return;
    }
  
  strcpy( Ind->info, *_info );
}


void get_population_info ( char **PopName,
			   char **_ind,
			   char **_info )
{
    POPULATION  ind = atol( *_ind ) - 1;
    TPopListEl  *PopListEl;
    TInd        *Ind;

    find_list_element( *PopName, &PopListEl );

    if ( NULL == PopListEl ) 
      {
	info( -2, "Population %s not defined\n", *PopName );
	return;
      }

    Ind = ( TInd* ) get_individual( PopListEl->Pop, ind );
    
    if ( NULL == Ind->info )
      strcpy( *_info, "none" );
    else
      strcpy( *_info, Ind->info );
}



void set_eff_weight( char **Name, double *weight )
{
    TEffects*Effect = find_Effect( *Name );

    if ( NULL == Effect ) {
        info( -2, "Can't find effect %s\n", *Name );
	return;
    }

    Effect->weight = ( float ) * weight;
    st_clear_all_population_gvalue();
}

void effect_weight_set_all ( double *weight )
{
    int i;

    for ( i = 0;i < NrEffects;i++ )
    {
        Effects[ i ].weight = *weight;
    }

    st_clear_all_population_gvalue();
}


void Insert_ALL_Loc( TEff*Eff, TMappoint*MP, int All )
{
    int i;

    for ( i = 0;i < Eff->NrLoci;i++ )
    {
        if ( MP == Eff->Loc[ i ].MP )
        {
            if ( Eff->Loc[ i ].NrAll++ < NoHoms )
            {
                Eff->Loc[ i ].All[ Eff->Loc[ i ].NrAll - 1 ] = All;
            }

            else
            {
                Eff->Loc[ i ].NrAll--;
                info( -2, "Error in the effect-file - too many alleles\n" );
            }

            return ;
        }
    }

    if ( 0 == Eff->NrLoci++ )
        Eff->Loc = ( TGen* ) yCalloc( 1, TGen );
    else
        Eff->Loc = ( TGen* ) yRealloc( Eff->Loc, Eff->NrLoci, TGen );

    Eff->Loc[ Eff->NrLoci - 1 ].MP = MP;
    Eff->Loc[ Eff->NrLoci - 1 ].NrAll = 1;
    Eff->Loc[ Eff->NrLoci - 1 ].All = ( int* ) yCalloc( NoHoms, int );
    Eff->Loc[ Eff->NrLoci - 1 ].All[ 0 ] = All;

    for ( i = 1;i < NoHoms;i++ )
        Eff->Loc[ Eff->NrLoci - 1 ].All[ i ] = -1;
}

void load_effmap_mod ( char**effectName,
		       char**effectSpec,
		       double*effectWeight )

{
    char line[ 10 * MAXLINE ], *token;
    int i, j, k, tmpi, status, maxLoc = 0;
    float tmpf;
    TMappoint*tmpMP = NULL;
    TEffects*Effect;

    if ( NULL != find_Effect( *effectName ) ){
        info( -2, "This effect is already loaded (same name)\n" );
	return;
    }

    if ( strlen( *effectSpec ) >= sizeof( line ) ){
        info( -2, "Effect specification too long\n" );
	return;
    }

    if ( 0 == NrEffects++ )
        Effects = ( TEffects* ) yCalloc( 1, TEffects );
    else
        Effects = ( TEffects* ) yRealloc( Effects, NrEffects, TEffects );

    Effect             = &Effects[ NrEffects - 1 ];
    Effect->GValue_avg = 0;
    Effect->weight     = 1;
    Effect->NrEff      = 0;
    Effect->Pos        = NrEffects;
    Effect->Eff        = NULL;
    Effect->MPset      = NULL;
    strcpy( Effect->Name, *effectName );
    Effect->weight = ( float ) * effectWeight;
    Effect->GValue_avg = 0;

    for ( i = 0; i < 1; i++ )
    {
        strcpy( line, *effectSpec );
        status = 0;
        token = strtok( line, "[ \t]" );

        while ( ( NULL != token ) && ( status != -1 ) && ( status != 1 ) )
        {
            switch ( status )
            {

            case 0:

                if ( 0 == strcmp( line, "\n" ) )
                    status = -1;
                else
                {
                    status = 2;

                    if ( 0 == ( Effect->NrEff ) ++ )
                        Effect->Eff = ( TEff* ) yCalloc( 1, TEff );
                    else
                        Effect->Eff = ( TEff* ) yRealloc( Effect->Eff,
                                                         Effect->NrEff,
                                                         TEff );

                    Effect->Eff[ Effect->NrEff - 1 ].Eff = 0;
                    Effect->Eff[ Effect->NrEff - 1 ].NrLoci = 0;
                    Effect->Eff[ Effect->NrEff - 1 ].Loc = NULL;
                }

                break;

            case 2:

                if ( NULL != ( tmpMP = find_mappoint( token ) ) )
                {
                    status = 3;
                }

                else
                {
                    status = -1;
                    info( -1, "Marker %s not found.\nline %d\n", token, i );
                }

                break;

            case 3:

                if ( 1 == sscanf( token, "%d", &tmpi ) )
                {
                    status = 4;
                    Insert_ALL_Loc( &Effect->Eff[ Effect->NrEff - 1 ],
				    tmpMP, tmpi );
                }

                else
                {
                    status = -1;
                    info( -1, "Integer expected\nline %d\n", i );
                }

                break;

            case 4:

                if ( 1 == sscanf( token, "%f", &tmpf ) )
                {
                    status = 1;
                    Effect->Eff[ Effect->NrEff - 1 ].Eff = tmpf;

                    if ( maxLoc < Effect->Eff[ Effect->NrEff - 1 ].NrLoci )
                    {
                        maxLoc = Effect->Eff[ Effect->NrEff - 1 ].NrLoci;
                    }
                }

                else
                {
                    status = 2;
                }

                break;
            }

            if ( 2 != status )
                token = strtok( NULL, "[ \t]" );
        }
    }

    if ( 1 != status ){
        if ( -1 != status )
            info( -1, "Incomplete effect specification\n" );
	goto fail;
    }

    Effect->MPs = 0;

    for ( i = 0;i < Effect->NrEff;i++ )
    {
        for ( j = 0;j < Effect->Eff[ i ].NrLoci;j++ )
        {
            for ( k = 0;k < Effect->MPs;k++ )
            {
                if ( *( Effect->MPset + k ) ==
                     ( Effect->Eff[ i ].Loc[ j ].MP ) )
                    break;
            }

            if ( k == Effect->MPs )
            {
                Effect->MPs++;

                if ( 1 == Effect->MPs )
                    Effect->MPset =
                        ( TMappoint** ) yCalloc( 1, TMappoint* );
                else
                    Effect->MPset =
                        ( TMappoint** ) yRealloc( Effect->MPset,
                                                 Effect->MPs,
                                                 TMappoint* );

                Effect->MPset[ Effect->MPs - 1 ] =
                    Effect->Eff[ i ].Loc[ j ].MP;
            }
        }
    }

    st_clear_all_population_gvalue();
    return;

 fail:

    for ( i = 0;i < Effect->NrEff;i++ )
    {
        for ( j = 0;j < Effect->Eff[ i ].NrLoci;j++ )
            yFree( Effect->Eff[ i ].Loc[ j ].All );

        yFree( Effect->Eff[ i ].Loc );
    }

    yFree( Effect->Eff );
    yFree( Effect->MPset );

    NrEffects--;

    if ( 0 == NrEffects )
        yFree( Effects );
    else
    {
        TEffects *tmp = realloc( Effects, NrEffects * sizeof( *tmp ) );
        if ( NULL != tmp ) Effects = tmp;
    }
}

void load_effmap(
		 char**file,
		 char**name)
{
    FILE*fp;
    char line[ 10 * MAXLINE ], *token;
    int i, j, k, tmpi, status, maxLoc = 0;
    float tmpf;
    TMappoint*tmpMP = NULL;
    TEffects*Effect;

    if ( NULL != find_Effect( *name ) ){
        info( -2, "This effect is already loaded (same name)\n" );
	return;
    }

    fp = fileopen( *file, "r" );
    
    if ( NULL == fp ){
        info( -2, "Effec file can not be opened\n" );
        return;
    }
    
    if ( 0 == NrEffects++ )
        Effects = ( TEffects* ) yCalloc( 1, TEffects );
    else
        Effects = ( TEffects* ) yRealloc( Effects, NrEffects, TEffects );

    Effect = &Effects[ NrEffects - 1 ];

    Effect->GValue_avg = 0;

    Effect->weight = 1;

    Effect->NrEff = 0;

    Effect->Pos = NrEffects;

    Effect->Eff = NULL;

    Effect->MPset = NULL;

    strcpy( Effect->Name, *name );

    for ( i = 0;;i++ )
    {
        if ( fgets( line, 10 * MAXLINE, fp ) == NULL )
            break;

        if ( (unsigned)( j = strcspn( line, "#" ) ) != strlen( line ) )
        {
            if ( j == 0 )
                strcpy( line, "" );
            else
            {
                line[ j++ ] = '\n';
                line[ j ] = '\0';
            }
        }

        status = 0;
        token = strtok( line, "[ \t]" );

        while ( ( NULL != token ) && ( status != -1 ) && ( status != 1 ) )
        {
            switch ( status )
            {

            case 0:

                if ( 0 == strcmp( token, "Weigth:" ) )
                {
                    status = 5;
                }

                else if ( 1 == sscanf( token, "%f", &Effect->GValue_avg ) )
                {
                    status = 1;
                }

                else
                {
                    if ( 0 == strcmp( line, "\n" ) )
                        status = -1;
                    else
                    {
                        status = 2;

                        if ( 0 == ( Effect->NrEff ) ++ )
                            Effect->Eff = ( TEff* ) yCalloc( 1, TEff );
                        else
                            Effect->Eff = ( TEff* ) yRealloc( Effect->Eff, Effect->NrEff, TEff );

                        Effect->Eff[ ( Effect->NrEff ) - 1 ].Eff = 0;

                        Effect->Eff[ ( Effect->NrEff ) - 1 ].NrLoci = 0;
                        Effect->Eff[ ( Effect->NrEff ) - 1 ].Loc = NULL;
                    }
                }

                break;

            case 2:

                if ( NULL != ( tmpMP = find_mappoint( token ) ) )
                {
                    status = 3;
                }

                else
                {
                    status = -1;
                    info( -1, "Marker %s not found.\nline %d\n", token, i );
                }

                break;

            case 3:

                if ( 1 == sscanf( token, "%d", &tmpi ) )
                {
                    status = 4;
                    Insert_ALL_Loc( &Effect->Eff[ Effect->NrEff - 1 ], tmpMP, tmpi );
                }

                else
                {
                    status = -1;
                    info( -1, "Integer expected\nline %d\n", i );
                }

                break;

            case 4:

                if ( 1 == sscanf( token, "%f", &tmpf ) )
                {
                    status = 1;
                    Effect->Eff[ Effect->NrEff - 1 ].Eff = tmpf;

                    if ( maxLoc < Effect->Eff[ Effect->NrEff - 1 ].NrLoci )
                    {
                        maxLoc = Effect->Eff[ Effect->NrEff - 1 ].NrLoci;
                    }
                }

                else
                {
                    status = 2;
                }

                break;

            case 5:

                if ( 1 == sscanf( token, "%f", &Effect->weight ) )
                {
                    status = 1;
                }

                else
                {
                    status = -1;
                    info( -1, "Error: weight" );
                }

                break;
            }

            if ( 2 != status )
                token = strtok( NULL, "[ \t]" );
        }
    }

    fileclose( fp );

    Effect->MPs = 0;

    for ( i = 0;i < Effect->NrEff;i++ )
    {
        for ( j = 0;j < Effect->Eff[ i ].NrLoci;j++ )
        {
            for ( k = 0;k < Effect->MPs;k++ )
            {
                if ( *( Effect->MPset + k ) == ( Effect->Eff[ i ].Loc[ j ].MP ) )
                    break;
            }

            if ( k == Effect->MPs )
            {
                Effect->MPs++;

                if ( 1 == Effect->MPs )
                    Effect->MPset = ( TMappoint** ) yCalloc( 1, TMappoint* );
                else
                    Effect->MPset = ( TMappoint** ) yRealloc( Effect->MPset, Effect->MPs, TMappoint* );

                Effect->MPset[ Effect->MPs - 1 ] = Effect->Eff[ i ].Loc[ j ].MP;
            }
        }
    }


    {
        int l, m;

        for ( i = 0;i < maxLoc;i++ )
        {
            for ( k = 0;k < NoHoms;k++ )
            {
                int oldEff = Effect->NrEff;
                int n;

                for ( n = k + 1;n < NoHoms;n++ )
                {
                    for ( j = 0;j < oldEff;j++ )
                    {
                        if ( ( i < Effect->Eff[ j ].NrLoci ) && ( k < Effect->Eff[ j ].Loc[ i ].NrAll ) &&
                                ( Effect->Eff[ j ].Loc[ i ].All[ k ] != Effect->Eff[ j ].Loc[ i ].All[ n ] ) )
                        {
                            Effect->NrEff++;
                            Effect->Eff = ( TEff* ) yRealloc( Effect->Eff, Effect->NrEff, TEff );
                            Effect->Eff[ Effect->NrEff - 1 ].Eff = Effect->Eff[ j ].Eff;
                            Effect->Eff[ Effect->NrEff - 1 ].NrLoci = Effect->Eff[ j ].NrLoci;
                            Effect->Eff[ Effect->NrEff - 1 ].Loc =
                                ( TGen* ) yCalloc( Effect->Eff[ j ].NrLoci, TGen );

                            for ( l = 0;l < Effect->Eff[ j ].NrLoci;l++ )
                            {
                                Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].MP =
                                    Effect->Eff[ j ].Loc[ l ].MP;
                                Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].NrAll =
                                    Effect->Eff[ j ].Loc[ l ].NrAll;
                                Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].All =
                                    ( int* ) yCalloc( NoHoms, int );

                                for ( m = 0;m < NoHoms;m++ )
                                {
                                    if ( ( l == i ) && ( m == k ) )
                                    {
                                        Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].All[ m ] = Effect->Eff[ j ].Loc[ l ].All[ n ];
                                    }

                                    else if ( ( l == i ) && ( m == n ) )
                                    {
                                        Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].All[ m ] = Effect->Eff[ j ].Loc[ l ].All[ k ];
                                    }

                                    else
                                    {
                                        Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].All[ m ] = Effect->Eff[ j ].Loc[ l ].All[ m ];
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
{}


void effmap_remove( char**effectName)
{
    if ( NULL == find_Effect( *effectName ) ) {
        info( -2, "Couldn't find effect: %s\n", *effectName );
	return;
    }

    int i, j, k;

    for ( i = NrEffects - 1;i >= 0;i-- ) {
        if ( 0 == strcmp( *effectName, Effects[ i ].Name ) ) {

            for ( j = 0;j < Effects[ i ].NrEff;j++ ) {
                for ( k = 0;k < Effects[ i ].Eff[ j ].NrLoci;k++ )
                    if ( NULL != Effects[ i ].Eff[ j ].Loc[ k ].All )
                        yFree( Effects[ i ].Eff[ j ].Loc[ k ].All );

                yFree( Effects[ i ].Eff[ j ].Loc );
            }

            if ( NULL != Effects[ i ].Eff )
                yFree( Effects[ i ].Eff );

            if ( NULL != Effects[ i ].MPset )
                yFree( Effects[ i ].MPset );

            for ( j = i + 1;j < NrEffects;j++ ) {
                Effects[ j - 1 ] = Effects[ j ];
                Effects[ j - 1 ].Pos = j;
            }

            NrEffects--;

            if ( 0 == NrEffects )
                yFree( Effects );
            else
                Effects = ( TEffects* ) yRealloc( Effects,
                                                  NrEffects,
                                                  TEffects );

            st_clear_all_population_gvalue();
            info( 1, "Effect %s removed\n", *effectName );
            break;
        }
    }
}

void effmap_remove_all(void)
{
    int i, j, k;

    for ( i = 0;i < NrEffects;i++ )
    {
        for ( j = 0; j < Effects[i].NrEff; j++ )
        {
            for ( k = 0;k < Effects[ i ].Eff[ j ].NrLoci;k++ )
            {
                if ( NULL != Effects[i].Eff[j].Loc[k].All )
                {
                    yFree( Effects[i].Eff[j].Loc[k].All );
                }
            }

            yFree( Effects[ i ].Eff[ j ].Loc );
        }

        if ( NULL != Effects[i].Eff )   yFree( Effects[i].Eff );
        if ( NULL != Effects[i].MPset ) yFree( Effects[i].MPset );

    }

    if ( NULL != Effects ) yFree( Effects );

    if ( 0 < NrEffects )  
      info( 1, "%d effects removed\n", NrEffects );

    NrEffects = 0;
    st_clear_all_population_gvalue();
}


void remove_genotype_population ( char **PopNames )
{
  char *PopName;

  PopName = strtok( *PopNames, "[' ','\n']" );

  while ( NULL != PopName )
    {
      TPopListEl *PopListEl;

      find_list_element( PopName, &PopListEl );

      if ( NULL == PopListEl )
	info( -1,
	      "Population not defined: %s\n",
	      PopName );
      else
	{
	  POPULATION ind;

	  for ( ind = 0;
		ind < PopListEl->Pop->NoIndsAlloc;
		ind++ )
	    {
	      if ( NULL !=
		   PopListEl->Pop->Ind[ ind ].genotype )
		{
		  free(
		    PopListEl->Pop->Ind[ ind ].genotype );

		  PopListEl->Pop->Ind[ ind ].genotype = NULL;
		}
	    }
          st_clear_population_pvalue(PopListEl->Pop);


	  info( 1,
		"Genotype removed for population %s\n",
		PopName );
	}

      PopName = strtok( NULL, "[' ','\n']" );
    }
}


void genotype_population( char **PopNames)
{
    char *PopName;

    TPopListEl *PopListEl;
    int alloc_error;

    if ( 0 == cMP )
      {
        info( -2, "No linkage map defined\n" );
	return;
      }

    PopName = strtok( *PopNames, "[' ','\n']" );

    while ( NULL != PopName )
    {
      
      find_list_element( PopName, &PopListEl );

      if ( NULL == PopListEl )
	info( -1, "Population not defined: %s\n", PopName );
      else if ( 0 == PopListEl->Pop->NoInds )
	info( -1, "Population %s is empty\n", PopName );
      else
        {
	  if ( st_population_all_genotyped( PopListEl->Pop ) ) {
	    info( 1, "Population %s already genotyped - genotyping skipped\n", PopName );
	  }  else {

	    alloc_error = 0;

#ifdef MT
  #pragma omp parallel 
#endif

	    {

	      POPULATION ind;
	      int  chrom, hom;

	      TInd   *Ind;
	      TChrom *Chrom;
	      int    i, m ;


#ifdef MT
  #pragma omp for 
#endif
  
	    for ( ind = 0; 
		  ind < PopListEl->Pop->NoInds; 
		  ind++ )
	    {
	      int *tmp;

	      Ind = get_individual( PopListEl->Pop, ind );

	      tmp = ( int* ) realloc(
				    Ind->genotype,
				    cMPhom * sizeof( *tmp )
				    );

	      if ( NULL == tmp )
		{
#ifdef MT
  #pragma omp atomic write
#endif
		  alloc_error = 1;
		  continue;
		}

	      Ind->genotype = tmp;

	      m = 0;

	      for ( chrom = 0; chrom < NoChroms; chrom++ )
		{

		  for ( i = 0; 
			i < PopDescr[chrom].NrMappoints;
			i++, m++ )
		    {
		      
		      Chrom = &( Ind->Chrom[ chrom ] );
		      
		      for ( hom = 0; hom < NoHoms; hom++ )
			{
			  Ind->genotype[m*NoHoms+hom] = 
			    get_allele( &( Chrom->Hom[ hom ] ),
					PopDescr[chrom].Mappoint[i].Pos );
			}

		    }
		}
	    }
	    } //parallel

	    if ( alloc_error )
	      {
		st_clear_population_genotype( PopListEl->Pop );
		info( -2, "Memory allocation failed in genotype_population\n" );
		return;
	      }
	    
	    info( 1, "Population %s genotyped\n", PopName );
	  }
        }

      PopName = strtok( NULL, "[' ','\n']" );
    }
}

void evaluate_population ( char**PopNames,
			   char**EffName   )
{
  char       *PopName;
  int        eff = -1;
  int        alloc_error = 0;

  TPopListEl *PopListEl;
  
  if ( 0 == NrEffects ) 
    {
      info( -2, "No effects defined\n" );
      return;
    }

  PopName = strtok( *PopNames, "[' ','\n']" );

  while ( NULL != PopName )
    {
      find_list_element( PopName, &PopListEl );

      if ( NULL == PopListEl )
	info( -1, "Population not defined: %s\n", PopName );
      else
        {
	  if ( 0 == PopListEl->Pop->NoInds )
	    {
	      info( -2, "Population %s is empty\n", PopName );
	      return;
	    }

	  if ( !st_population_all_genotyped( PopListEl->Pop ) )
	    {
	      info( -2, "First genotype population %s\n", PopName );
	      return;
	    }
	  
	  if ( NULL != EffName )
            {
	      TEffects * Effect;
	      
	      if ( NULL == ( Effect = find_Effect( *EffName ) ) ){
		info( -2, "Can't find the effect %s\n", *EffName );
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
	    TInd       *Ind;
	    int        i, j, hom, loc, status;


#ifdef MT
  #pragma omp for 
#endif
 
	  for ( ind = 0;
		ind < PopListEl->Pop->NoInds;
		ind++ )
          {
	    float *tmp;

	    Ind = ( TInd* ) get_individual( PopListEl->Pop, ind );

	    tmp = ( float* ) realloc(
	      Ind->GValue,
	      ( NrEffects + 1 ) * sizeof( *tmp )
	    );

	    if ( NULL == tmp )
	      {
#ifdef MT
  #pragma omp atomic write
#endif
		alloc_error = 1;
		continue;
	      }

	    Ind->GValue = tmp;

	    Ind->GValue[0] = 0;
	      
	    for ( i = 1; 
		  i < NrEffects + 1; 
		  i++ )
            {
	      Ind->GValue[i] = Effects[i-1].GValue_avg;

	      for ( j = 0; 
		    j < Effects[i-1].NrEff;
		    j++ )
		{
		  status = 0;
		  loc = 0;

		  while ( ( 0 == status ) 
			  && 
			  ( loc < Effects[i-1].Eff[j].NrLoci ) )
		    {
		      hom = 0;
		      
		    while ( ( 0 == status ) && ( hom < NoHoms ) )
		    {
		      if ( ( -1 < Effects[i-1].Eff[j].Loc[loc].All[hom] ) &&
			   ( Effects[i-1].Eff[j].Loc[loc].All[hom] !=
			     Ind->genotype[Effects[i-1].Eff[j].Loc[loc].MP->counter * NoHoms + hom ] ) )
	              {
			status = 1;
		      }
		      
		      hom++;
		    }
		      
		    loc++;
		  }
		  
		  if ( 0 == status )
		    Ind->GValue[i] += Effects[i-1].Eff[j].Eff;
		}
		
		Ind->GValue[0] += Ind->GValue[i] * Effects[i-1].weight;
	      }
	      
	      if ( eff > -1 )
		Ind->GValue[0] = Ind->GValue[eff];
            }
	  } // parallel

	  if ( alloc_error )
	    {
	      st_clear_population_gvalue( PopListEl->Pop );
	      info( -2, "Memory allocation failed in evaluate_population\n" );
	      return;
	    }
	  
	  info( 1, "Population %s evaluated; selection index: %s\n",
		PopName, ( -1 == eff ) ? "weigthed sum" : *EffName );
        }

        PopName = strtok( NULL, "[' ','\n']" );
    }
}


void remove_evaluate_population ( char **PopNames )
{
  char *PopName;

  PopName = strtok( *PopNames, "[' ','\n']" );

  while ( NULL != PopName )
    {
      TPopListEl *PopListEl;

      find_list_element( PopName, &PopListEl );

      if ( NULL == PopListEl )
        info( -1,
              "Population not defined: %s\n",
              PopName );
      else
        {
          st_clear_population_gvalue(PopListEl->Pop);

          info( 1,
                "Evaluation removed for population %s\n",
                PopName );
        }

      PopName = strtok( NULL, "[' ','\n']" );
    }
}


static int st_sort_criterion = 0;
static int st_sort_eff = 0;
static int st_sort_decreasing = 1;


static int st_indcmp(const void *a, const void *b)
{
  const TInd *ia = (const TInd *)a;
  const TInd *ib = (const TInd *)b;
  const float *va =
    (0 == st_sort_criterion) ? ia->GValue : ia->PValue;
  const float *vb =
    (0 == st_sort_criterion) ? ib->GValue : ib->PValue;

  if (va[st_sort_eff] > vb[st_sort_eff])
    return st_sort_decreasing ? -1 : 1;

  if (va[st_sort_eff] < vb[st_sort_eff])
    return st_sort_decreasing ? 1 : -1;

  return 0;
}


void sort_population(
    char **PopNames,
    int *decreasing,
    char **criterion_name,
    char **EffName,
    int *retval)
{
  TPopListEl *PopListEl;
  char *PopName;
  int criterion;
  int eff;

  *retval = -2;

  criterion = st_selection_criterion(criterion_name);

  if (0 > criterion)
    return;

  if (!st_effect_position(EffName,&eff))
    return;

  PopName = strtok(*PopNames," ");

  while (NULL != PopName)
    {
      find_list_element(PopName,&PopListEl);

      if (NULL == PopListEl)
        {
          info(-2,"Population not defined: %s\n",PopName);
          return;
        }

      if (0 == PopListEl->Pop->NoInds)
        {
          info(-2,"Population %s is empty\n",PopName);
          return;
        }

      if (!st_population_selection_ready(
             PopListEl->Pop,criterion,eff))
        {
          info(-2,
               (0 == criterion) ?
               "First evaluate population %s\n" :
               "First phenotype population %s\n",
               PopName);
          return;
        }

      st_sort_criterion = criterion;
      st_sort_eff = eff;
      st_sort_decreasing = *decreasing ? 1 : 0;

      qsort(PopListEl->Pop->Ind,
            PopListEl->Pop->NoInds,
            sizeof(TInd),
            st_indcmp);

      PopName = strtok(NULL," ");
    }

  *retval = 0;
}


int getPos( int all, TSet Set )
{
    int i;

    for ( i = 0;i < Set.betrag;i++ )
    {
        if ( Set.member[ i ] == all )
        {
            return i;
        }
    }

    info( -2, "Fehler in getPos()\n" );
    return 0;
}

SEXP return_population(SEXP _PopNames, SEXP _bg)
{
   int i_bg = 0;
   int nprotect = 0;
   int chrom, hom, marker, cmarker, allele;
   int NoAll = 0, NoInds = 0;
   R_xlen_t col = 0, row = 0, offset;
   R_xlen_t matrix_length;
   POPULATION ind;
   TPopListEl *PopListEl;
   char *PopNames = NULL;
   char *PopName = NULL;
   char tmp[2 * MAXNAME];
   int *bg = &i_bg;
   size_t pop_names_length;
   SEXP result = R_NilValue;
   SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;

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

   i_bg = INTEGER(_bg)[0];

   pop_names_length = strlen(CHAR(STRING_ELT(_PopNames, 0)));

   if (pop_names_length > SIZE_MAX - 2)
   {
      info(-2, "Invalid population name in return_population\n");
      goto cleanup;
   }

   PopNames = (char *)malloc(
      (pop_names_length + 2) * sizeof(*PopNames));

   if (PopNames == NULL)
   {
      info(-2, "Memory allocation failed in return_population\n");
      goto cleanup;
   }

   strcpy(PopNames, CHAR(STRING_ELT(_PopNames, 0)));

   if ((cMP < 0) ||
       ((size_t)cMP > SIZE_MAX / sizeof(*AllSet)))
   {
      info(-2, "Invalid number of markers in return_population\n");
      goto cleanup;
   }

   if (cMP > 0)
   {
      AllSet = (TSet *)malloc((size_t)cMP * sizeof(*AllSet));

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
                  if (!insertmember(
                         AllSet + marker,
                         Ind->genotype[NoHoms * marker + hom]))
                  {
                     info(
                        -2,
                        "Memory allocation failed in return_population\n");
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
         qsort(AllSet[marker].member,
               AllSet[marker].betrag,
               sizeof(int),
               intcmp);
      }

      if (AllSet[marker].betrag > INT_MAX - NoAll)
      {
         info(-2, "Too many alleles in return_population\n");
         goto cleanup;
      }

      NoAll += AllSet[marker].betrag;
   }

   PROTECT(matrix = allocMatrix(INTSXP, NoAll, NoInds));
   nprotect++;

   PROTECT(matrixDim = allocVector(INTSXP, 2));
   nprotect++;

   INTEGER(matrixDim)[0] = NoAll;
   INTEGER(matrixDim)[1] = NoInds;
   setAttrib(matrix, R_DimSymbol, matrixDim);

   PROTECT(matrixDimNames = allocVector(VECSXP, 2));
   nprotect++;

   PROTECT(xAxis = allocVector(STRSXP, NoInds));
   nprotect++;

   PROTECT(yAxis = allocVector(STRSXP, NoAll));
   nprotect++;

   matrix_length = XLENGTH(matrix);

   for (offset = 0; offset < matrix_length; offset++)
      INTEGER(matrix)[offset] = 0;

   strcpy(PopNames, CHAR(STRING_ELT(_PopNames, 0)));

   PopName = strtok(PopNames, " ");
   offset = 0;

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

         snprintf(tmp, sizeof(tmp),
                  "%s.%ld ",
                  PopListEl->Pop->Name, (long)(ind + 1));

         SET_STRING_ELT(xAxis, col++, mkChar(tmp));

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
                  INTEGER(matrix)[
                     offset +
                     getPos(Ind->genotype[marker * NoHoms + hom],
                            AllSet[marker])] = 1;
               }
            }

            offset += AllSet[marker].betrag;
         }
      }

      PopName = strtok(NULL, " ");
   }

   for (chrom = 0, marker = 0, row = 0;
        chrom < NoChroms;
        chrom++)
   {
      for (cmarker = 0;
           cmarker < PopDescr[chrom].NrMappoints;
           cmarker++, marker++)
      {
         for (allele = 0;
              allele < AllSet[marker].betrag;
              allele++)
         {
            snprintf(tmp, sizeof(tmp),
                     "%s.%d ",
                     PopDescr[chrom].Mappoint[cmarker].name,
                     AllSet[marker].member[allele]);

            SET_STRING_ELT(yAxis, row++, mkChar(tmp));
         }
      }
   }

   SET_VECTOR_ELT(matrixDimNames, 1, xAxis);
   SET_VECTOR_ELT(matrixDimNames, 0, yAxis);
   setAttrib(matrix, R_DimNamesSymbol, matrixDimNames);

   result = matrix;

cleanup:

   if (AllSet != NULL)
   {
      for (marker = 0; marker < cMP; marker++)
      {
         free(AllSet[marker].member);
      }

      free(AllSet);
      AllSet = NULL;
   }

   free(PopNames);

   if (nprotect > 0)
      UNPROTECT(nprotect);

   return result;
}

void evaluate_genome ( char**  PopName,
		       int*    all,
		       double* ratio,
		       double* length,
		       double* total,
		       int*    blocks   )
{
  TPopListEl*PopListEl;
  POPULATION ind;
  int chrom, hom, loc;
  double gLength = 0;
  TInd*Ind;

  find_list_element( *PopName, &PopListEl );

  if ( NULL == PopListEl ){
    info( -2, "Can't find Population %s\n", *PopName );
    return;
  }

  for ( chrom = 0;chrom < NoChroms;chrom++ )
    gLength += PopDescr[ chrom ].ChromLen * NoHoms;

  for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ ) {
    Ind = get_individual( PopListEl->Pop, ind);
    for ( chrom = 0; chrom < NoChroms; chrom++ ) {
      for ( hom = 0; hom < NoHoms; hom++ ) {
	for ( loc = 0; loc < Ind->Chrom[chrom].Hom[hom].NoLociUsed; loc++ ) {
	  if ( *all == Ind->Chrom[ chrom ].Hom[ hom ].Loc[ loc ].All )
	    {
	      if ( loc + 1 < Ind->Chrom[ chrom ].Hom[ hom ].NoLociUsed )
		{
		  length[ind] += Ind->Chrom[chrom].Hom[hom].Loc[loc+1].Pos -
		                 Ind->Chrom[chrom].Hom[hom].Loc[loc].Pos;
		}
	      else
		{
		  length[ind] += PopDescr[chrom].ChromLen - 
                                 Ind->Chrom[chrom].Hom[hom].Loc[loc].Pos;
		}
	      blocks[ ind ] ++;
	    }
	}
      }
    }
    ratio[ ind ] = length[ ind ] / gLength;
    total[ ind ] = gLength;
  }
}

void evaluate_genome_region ( char**  PopName,
			      int*    all,
			      int*    chrom,
			      double* b,
			      double* e,
			      double* ratio,
			      double* length,
			      double* total,
			      int*    blocks )
{
  TPopListEl*PopListEl;
  POPULATION ind;
  int hom, loc;
  double gLength = 0;
  double begin, end;
  TInd*Ind;
  
  find_list_element( *PopName, &PopListEl );
  
  if ( NULL == PopListEl ){
    info( -2, "Can't find Population %s\n", *PopName );
    return;
  }

  if ( ( *chrom < 0 ) || ( *chrom >= NoChroms ) ) {
    info( -2, "Chromosome out of range: %d\n", *chrom + 1 );
    return;
  }
  
  gLength += ( *e - *b ) * NoHoms;
  
  for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
    {
      Ind = get_individual( PopListEl->Pop, ind );
      
      for ( hom = 0;hom < NoHoms;hom++ )
	{
	  for ( loc = 0; loc < Ind->Chrom[*chrom].Hom[hom].NoLociUsed; loc++ )
	    {
	      if ( *all == Ind->Chrom[*chrom].Hom[hom].Loc[loc].All )
		{
		  begin = Ind->Chrom[*chrom].Hom[hom].Loc[loc].Pos;
		  
		  if ( loc + 1 < Ind->Chrom[*chrom].Hom[hom].NoLociUsed )
		    {
		      end = Ind->Chrom[*chrom].Hom[hom].Loc[loc+1].Pos;
		    }
		  else
		    {
		      end = PopDescr[ *chrom ].ChromLen;
		    }
		  
		  if ( ( end >= *b ) && ( begin <= *e ) )
		    {
		      if ( begin < *b ) begin = *b;
		      if ( end   > *e ) end   = *e;
		      length[ ind ] += end - begin;
		      blocks[ ind ] ++;
		    }
		}
	    }
	}
      ratio[ ind ] = length[ ind ] / gLength;
      total[ ind ] = gLength;
    }
}


void reset_mdp(void)
{
  MDP = 0;
  info(2, "MDP counter reset\n" );
}

void set_mdp( char** cMDP )
{
  MDP = atol( *cMDP );
}

void evaluate_mdp ( char** progeny,
		    char** parent1,
		    char** parent2,
		    char** effectfile )
{
  TEffects   *Effect = NULL;
  unsigned long int oldMDP = MDP;
  int i;

  if ( 0 == NrEffects ) {
    info( -2, "No Effectfiles loaded\n" );
    return;
  }

  if ( NULL == effectfile ) {
    info( -2, "No effect specified\n" );
    return;
  }

  if ( ( NULL != effectfile ) 
       && ( NULL == ( Effect = find_Effect(*effectfile) ) ) ) {
    info( -2, "Can't find the effect %s\n", *effectfile );
    return;
  }

  {
    TPopListEl* PopListEl;
    POPULATION  ind;
    int* counter = NULL;
    int* all = NULL;
    int  hom;
    TInd*Ind;
    int error = 0;

    if ( Effect->MPs > 0 ) {
      counter = ( int* ) malloc(
                  ( size_t ) Effect->MPs * sizeof( *counter ) );
      all = ( int* ) malloc(
              ( size_t ) Effect->MPs * sizeof( *all ) );

      if ( ( NULL == counter ) || ( NULL == all ) ) {
	info( -2, "Memory allocation failed in evaluate_mdp\n" );
	error = 1;
	goto cleanup;
      }
    }

    for ( i = 0; i < Effect->MPs; i++ ) {
      counter[ i ] = 0;
      all[ i ] = -1;
    }

    find_list_element( *parent1, &PopListEl );

    if ( NULL == PopListEl ) {
      info( -2, "Population not defined: %s\n", *parent1 );
      error = 1;
      goto cleanup;
    }

    if ( !st_population_all_genotyped( PopListEl->Pop ) ) {
      info( -2, "First genotype population %s\n", *parent1 );
      error = 1;
      goto cleanup;
    }

    for ( ind = 0; 
	  ind < PopListEl->Pop->NoInds; 
	  ind++ ) {
      Ind = get_individual( PopListEl->Pop, ind );
      for ( i = 0; i < Effect->MPs; i++ ) {
	if ( counter[ i ] < 2 ) {
	  for ( hom = 0;hom < NoHoms;hom++ ) {
	    if ( all[ i ] 
		 != Ind->genotype[ NoHoms * Effect->MPset[i] ->counter + hom ] )
	      {
		all[ i ] =
		  Ind->genotype[
		    NoHoms * Effect->MPset[i] ->counter + hom
		  ];
		counter[ i ] ++;
	      }
	  }
	}
      }
    }

    for ( i = 0; i < Effect->MPs; i++ ) {
      if ( 1 == counter[ i ] ) {
	counter[ i ] = 0;
	all[ i ] = -1;
      }
    }

    find_list_element( *parent2, &PopListEl );

    if ( NULL == PopListEl ){
      info( -2, "Population not defined: %s\n", *parent2 );
      error = 1;
      goto cleanup;
    }

    if ( !st_population_all_genotyped( PopListEl->Pop ) ) {
      info( -2, "First genotype population %s\n", *parent2 );
      error = 1;
      goto cleanup;
    }

    for ( ind = 0; ind < PopListEl->Pop->NoInds; ind++ ) {
      Ind = get_individual( PopListEl->Pop, ind );
      for ( i = 0; i < Effect->MPs; i++ ) {
	if ( counter[ i ] < 2 ) {
	  for ( hom = 0; hom < NoHoms; hom++ ) {
	    if ( all[ i ] 
		 != Ind->genotype[
		      NoHoms * Effect->MPset[ i ] ->counter + hom
		    ] )
	      {
		all[ i ] =
		  Ind->genotype[
		    NoHoms * Effect->MPset[ i ] ->counter + hom
		  ];
		counter[ i ] ++;
	      }
	  }
	}
      }
    }
    
    find_list_element( *progeny, &PopListEl );

    if ( NULL == PopListEl ){
      info( -2, "Population not defined: %s\n", *progeny );
      error = 1;
      goto cleanup;
    }

    for ( i = 0; i < Effect->MPs; i++ ) {
      if ( counter[ i ] == 2 ) {
	MDP += PopListEl->Pop->NoInds;
      }
    }

  cleanup:

    free( counter );
    free( all );

    if ( error ) return;
  }

  info( 1, "%lu marker data points required in this step (total %lu)\n",
	MDP - oldMDP, MDP );
}


void get_mdp ( char** cmdp )
{
  snprintf(*cmdp, strlen(*cmdp),
	   "%lu", MDP);
}

typedef struct {
    int oldPos;
    TMappoint*Mappoint;
} TGenPop;

TGenPop* PGenPop;

int genpopcmp(  const void* i, const void*j  )
{
    int chrom0, chrom1;
    CPOS pos0, pos1;
    chrom0 = ( ( TGenPop* ) i ) ->Mappoint->chrom;
    chrom1 = ( ( TGenPop* ) j ) ->Mappoint->chrom;
    pos0 = ( ( TGenPop* ) i ) ->Mappoint->Pos;
    pos1 = ( ( TGenPop* ) j ) ->Mappoint->Pos;

    if ( chrom0 > chrom1 )
        return 1;

    if ( chrom0 < chrom1 )
        return -1;

    if ( pos0 > pos1 )
        return 1;

    if ( pos0 < pos1 )
        return -1;

    return 0;
}

void generate_population(
    int*  no_pop,
    int*  no_ind,
    int*  no_mar,
    int*  no_all,
    int*  no_pop_u,
    int*  pop_u,
    int*  no_ind_u,
    int*  ind_u,
    int*  no_mar_u,
    int*  mar_u,
    double* data,
    double* freq ,
    char**  pop_list,
    char**  mar_list,
    int*    al_list,
    int*    bg,
    int*    bc )
{
  int i, j, k, l;
  double*ptr_data, *ptr_dataC, *dataC = NULL;
  double count;
  int rowpos_d;
  size_t colpos, dataC_count;
  int nrow_d, ncol_d;
  int*coffset = NULL;
  int*roff_d = NULL;
  int*roff_f = NULL;
  int*_ind = NULL, *_chrom = NULL, *_hom = NULL, *_all = NULL;
  int*candidates = NULL;
  int*tmp_int;
  double*_pos = NULL, *tmp_double;

  (void)*freq;

  PGenPop = NULL;

  if ( ( *no_pop < 1 ) || ( *no_mar < 1 ) ||
       ( *no_pop_u < 1 ) || ( *no_mar_u < 1 ) ||
       ( NoHoms < 1 ) )
    goto invalid_or_allocation_fail;

  coffset = ( int* ) malloc( ( size_t ) *no_pop * sizeof( *coffset ) );
  roff_d = ( int* ) malloc( ( size_t ) *no_mar * sizeof( *roff_d ) );
  roff_f = ( int* ) malloc( ( size_t ) *no_mar_u * sizeof( *roff_f ) );

  if ( ( NULL == coffset ) || ( NULL == roff_d ) ||
       ( NULL == roff_f ) )
    goto invalid_or_allocation_fail;
  
  coffset[ 0 ] = 0;
  for ( k = 1; k < *no_pop; k++ )
    {
      coffset[ k ] = coffset[ k - 1 ] + no_ind[ k - 1 ];
    }
  
  roff_d[ 0 ] = 0;
  for ( i = 1; i < *no_mar; i++ )
    {
      roff_d[ i ] = roff_d[ i - 1 ] + no_all[ i - 1 ];
    }

  roff_f[ 0 ] = 0;
  for ( i = 1;i < *no_mar_u;i++ )
    {
      roff_f[ i ] = roff_f[ i - 1 ] + no_all[ mar_u[ i - 1 ] - 1 ];
    }
  
    nrow_d = roff_d[ *no_mar - 1 ] + no_all[ *no_mar - 1 ];
    ncol_d = coffset[ *no_pop - 1 ] + no_ind[ *no_pop - 1 ];

    /* Commentend out due to compiler warning 2026-01-17
    int  nrow_f
    nrow_f = roff_f[ *no_mar_u - 1 ] + no_all[ mar_u[ *no_mar_u - 1 ] - 1 ];
    */

    {
      int sum_ind_u = 0;

      if ( ( nrow_d < 1 ) || ( ncol_d < 1 ) ||
           ( ( size_t ) nrow_d > SIZE_MAX / ( size_t ) ncol_d ) )
        goto invalid_or_allocation_fail;

      dataC_count = ( size_t ) nrow_d * ( size_t ) ncol_d;
      if ( dataC_count > SIZE_MAX / sizeof( *dataC ) )
        goto invalid_or_allocation_fail;

      dataC = ( double* ) malloc( dataC_count * sizeof( *dataC ) );
      if ( NULL == dataC )
        goto invalid_or_allocation_fail;

      for ( k = 0;k < *no_pop_u;k++ )
        {
	  for ( l = 0;l < no_ind_u[ k ];l++ )
            {
	      colpos = ( ( size_t ) coffset[ pop_u[ k ] - 1 ] +
                         ( size_t ) ( ind_u[ l + sum_ind_u ] - 1 ) ) *
                       ( size_t ) nrow_d;
	      for ( i = 0;i < *no_mar_u;i++ )
                {
		  count = 0;
		  for ( j = 0;j < no_all[ mar_u[ i ] - 1 ];j++ )
                    {
		      rowpos_d = ( roff_d[ mar_u[ i ] - 1 ] + j );
		      ptr_data = &data[ colpos + rowpos_d ];
		      if ( !ISNA( *ptr_data ) && ( *ptr_data > 0 ) )
                        {
			  *ptr_data = 1;
			  count++;
                        }
                    }
		  for ( j = 0;j < no_all[ mar_u[ i ] - 1 ];j++ )
                    {
		      rowpos_d = ( roff_d[ mar_u[ i ] - 1 ] + j );
		      ptr_data = &data[ colpos + rowpos_d ];
		      ptr_dataC = &dataC[ colpos + rowpos_d ];
		      if ( count > 0 )
			( *ptr_dataC ) = ( *ptr_data ) /
                                         ( double ) count;
		      else
			( *ptr_dataC ) = ( *ptr_data );
                    }
                }
            }
	  sum_ind_u += no_ind_u[ k ];
        }
      
      {
	int chrom, hom, tokencounter = 0;
	char*token;
	
	if ( 0 == cMP ){
	  info( -2, "No linkage map defined\n" );
	  goto fail;
	}

	PGenPop = ( TGenPop* ) malloc(
                     ( size_t ) *no_mar * sizeof( *PGenPop ) );
	if ( NULL == PGenPop )
	  goto invalid_or_allocation_fail;
	token = strtok( *mar_list, " " );
	
	while ( NULL != token )
	  {
	    PGenPop[ tokencounter ].Mappoint = NULL;
	    for ( chrom = 0;chrom < NoChroms;chrom++ )
	      {
		for ( i = 0;i < PopDescr[ chrom ].NrMappoints;i++ )
		  {
		    if ( 0 == strcmp(
                           PopDescr[ chrom ].Mappoint[ i ].name,
                           token ) )
		      {
			PGenPop[ tokencounter ].oldPos = tokencounter;
			PGenPop[ tokencounter ].Mappoint =
                          &( PopDescr[ chrom ].Mappoint[ i ] );
		      }
		  }
	      }

	    if ( NULL == PGenPop[ tokencounter ].Mappoint ){
	      info( -2, "Can't find marker %s!", token );
	      goto fail;
	    }

	    tokencounter++;
	    
	    token = strtok( NULL, " " );
	  }
	
	qsort( PGenPop, *no_mar, sizeof( TGenPop ), genpopcmp );

	{
	  size_t offset = 0;
	  token = strtok( *pop_list, " " );
	  
	  for ( k = 0;k < *no_pop_u;k++ ){
	    int nrow = 0;
	    int ind;

	    _ind = ( int* ) malloc( sizeof( *_ind ) );
	    _chrom = ( int* ) malloc( sizeof( *_chrom ) );
	    _hom = ( int* ) malloc( sizeof( *_hom ) );
	    _pos = ( double* ) malloc( sizeof( *_pos ) );
	    _all = ( int* ) malloc( sizeof( *_all ) );

	    if ( ( NULL == _ind ) || ( NULL == _chrom ) ||
	         ( NULL == _hom ) || ( NULL == _pos ) ||
	         ( NULL == _all ) )
	      goto invalid_or_allocation_fail;
	      
	    for ( ind = 0;ind < no_ind_u[ k ];ind++ ) {
	      for ( chrom = 0;chrom < NoChroms;chrom++ ) {
		int tmpcounter = 0;
		CPOS lastpos = 0;
		for ( i = 0;i < *no_mar;i++ ) {
		  int nr_cand = NoHoms, c, act_cand = 0;

		  candidates = ( int* ) malloc(
                                 ( size_t ) NoHoms *
                                 sizeof( *candidates ) );
		  if ( NULL == candidates )
		    goto invalid_or_allocation_fail;

		  if ( PGenPop[ i ].Mappoint->chrom - 1 == chrom ) {
		    for ( c = 0;c < NoHoms;c++ ){
		      candidates[ c ] = *bg;
		    }
		    for ( j = 0;
                          j < no_all[ PGenPop[ i ].oldPos ];
                          j++ ) {
		      while ( dataC[
                                offset +
                                ( size_t ) ind * ( size_t ) nrow_d +
                                roff_d[ PGenPop[ i ].oldPos ] + j
                              ] > 0 ) {
			if ( nr_cand == act_cand ){
			  ++nr_cand;
			  tmp_int = ( int* ) realloc(
                                     candidates,
                                     ( size_t ) nr_cand *
                                     sizeof( *candidates ) );
			  if ( NULL == tmp_int )
			    goto invalid_or_allocation_fail;
			  candidates = tmp_int;
			}
			candidates [ act_cand++ ] 
			  = ( int ) al_list[
                              roff_d[ PGenPop[ i ].oldPos ] + j ];
			dataC[
                          offset +
                          ( size_t ) ind * ( size_t ) nrow_d +
                          roff_d[ PGenPop[ i ].oldPos ] + j
                        ] -= ( float ) 1 / NoHoms;
		      }
		    }  
			
		    qsort( candidates, nr_cand, sizeof(int), intcmp );
			
		    for ( hom = 0;hom < NoHoms;hom++ )
		      {
			nrow++;
			    
			tmp_int = ( int* ) realloc(
                           _ind, ( size_t ) nrow * sizeof( *_ind ) );
			if ( NULL == tmp_int )
			  goto invalid_or_allocation_fail;
			_ind = tmp_int;

			tmp_int = ( int* ) realloc(
                           _chrom, ( size_t ) nrow * sizeof( *_chrom ) );
			if ( NULL == tmp_int )
			  goto invalid_or_allocation_fail;
			_chrom = tmp_int;

			tmp_int = ( int* ) realloc(
                           _hom, ( size_t ) nrow * sizeof( *_hom ) );
			if ( NULL == tmp_int )
			  goto invalid_or_allocation_fail;
			_hom = tmp_int;

			tmp_double = ( double* ) realloc(
                              _pos, ( size_t ) nrow * sizeof( *_pos ) );
			if ( NULL == tmp_double )
			  goto invalid_or_allocation_fail;
			_pos = tmp_double;

			tmp_int = ( int* ) realloc(
                           _all, ( size_t ) nrow * sizeof( *_all ) );
			if ( NULL == tmp_int )
			  goto invalid_or_allocation_fail;
			_all = tmp_int;
			
			*( _ind + nrow - 1 ) = ind + 1;
			*( _chrom + nrow - 1 ) = chrom + 1;
			*( _hom + nrow - 1 ) = hom + 1;
			  
			
			if ( 0 == tmpcounter ) {
			  * ( _pos + nrow - 1 ) = 0;
			} else {
			  *( _pos + nrow - 1 ) 
			    = ( PGenPop[ i ].Mappoint->Pos +
                                lastpos ) / 2;
			}

			{
			  int sel_all;

			  if ( 1 == *bc ) {
			    sel_all = hom;
			  } else {
			    sel_all = (int) gsl_rng_uniform_int(
                               rng[0],
                               (unsigned long int)nr_cand );
			  }
			  nr_cand--;
			  *( _all + nrow - 1 ) =
                            candidates[ sel_all ];
			  candidates[ sel_all ] = candidates[ nr_cand ];
			}
		      }

		    tmpcounter++;
		    lastpos = PGenPop[ i ].Mappoint->Pos;
			
		    if ( nr_cand > 0 ){
		      info( -1,
			    "More than two alleles. Marker %s Pop: %s Ind: %d\n",
			    PGenPop[ i ].Mappoint->name, token,
			    ind + 1 );
		    }
		  }
		    
		  free( candidates );
		  candidates = NULL;
		}
		  
		if ( 0 == tmpcounter ){
		  if ( 0 == ind ){
		    info( -1,
                          "No marker defined for chromosome %d.\n",
                          chrom + 1 );
		  }

		  for ( hom = 0;hom < NoHoms;hom++ ){
		    nrow++;

		    tmp_int = ( int* ) realloc(
                         _ind, ( size_t ) nrow * sizeof( *_ind ) );
		    if ( NULL == tmp_int )
		      goto invalid_or_allocation_fail;
		    _ind = tmp_int;

		    tmp_int = ( int* ) realloc(
                         _chrom, ( size_t ) nrow * sizeof( *_chrom ) );
		    if ( NULL == tmp_int )
		      goto invalid_or_allocation_fail;
		    _chrom = tmp_int;

		    tmp_int = ( int* ) realloc(
                         _hom, ( size_t ) nrow * sizeof( *_hom ) );
		    if ( NULL == tmp_int )
		      goto invalid_or_allocation_fail;
		    _hom = tmp_int;

		    tmp_double = ( double* ) realloc(
                            _pos, ( size_t ) nrow * sizeof( *_pos ) );
		    if ( NULL == tmp_double )
		      goto invalid_or_allocation_fail;
		    _pos = tmp_double;

		    tmp_int = ( int* ) realloc(
                         _all, ( size_t ) nrow * sizeof( *_all ) );
		    if ( NULL == tmp_int )
		      goto invalid_or_allocation_fail;
		    _all = tmp_int;

		    *( _ind + nrow - 1 ) = ind + 1;
		    *( _chrom + nrow - 1 ) = chrom + 1;
		    *( _hom + nrow - 1 ) = hom + 1;
		    *( _pos + nrow - 1 ) = 0;
		    *( _all + nrow - 1 ) = *bg;
		  }
		}
		
	      }
	    }
	    
	    offset += ( size_t ) ind * ( size_t ) nrow_d;

	    {
	      char cnrow[ 50 ];
	      char*prow = cnrow;

	      snprintf(cnrow, sizeof(cnrow),
		       "%d", nrow);

              {
                TPopListEl *generated = NULL;

                find_list_element( token, &generated );
                if ( NULL != generated )
                  {
                    info( -2, "First remove population %s\n", token );
                    goto fail;
                  }

	        init_population( &token,
			         &prow,
			         &( _ind[ 0 ] ),
			         &( _chrom[ 0 ] ),
			         &( _hom[ 0 ] ),
			         &( _pos[ 0 ] ),
			         &( _all[ 0 ] ) );

                find_list_element( token, &generated );
                if ( NULL == generated )
                  goto fail;
              }
	    }
	    
	    token = strtok( NULL, " " );
	    nrow = 0;
	    free( _ind );
	    free( _chrom );
	    free( _hom );
	    free( _pos );
	    free( _all );
	    _ind = NULL;
	    _chrom = NULL;
	    _hom = NULL;
	    _pos = NULL;
	    _all = NULL;
	  }
	}
	
	free( PGenPop );
	PGenPop = NULL;
      }
      
      free( dataC );
    }

    free( coffset );
    free( roff_d );
    free( roff_f );
    return;

invalid_or_allocation_fail:

    info( -2,
          "Invalid dimensions or memory allocation failed in "
          "generate_population\n" );

fail:

    free( candidates );
    free( _ind );
    free( _chrom );
    free( _hom );
    free( _pos );
    free( _all );
    free( PGenPop );
    PGenPop = NULL;
    free( dataC );
    free( coffset );
    free( roff_d );
    free( roff_f );
    return;
}

SEXP splitdt( SEXP in_string )
{
  const char *ma;
  const char *a;
  int i;
  int n = length( in_string );
  SEXP element;
  SEXP out_string;
  
  PROTECT( out_string = allocMatrix( STRSXP, n, 2 ) );
  for ( i = 0;i < n;i++ )
    {
      element = STRING_ELT( in_string, i );

      if ( NA_STRING == element )
	{
	  SET_STRING_ELT( out_string, i, NA_STRING );
	  SET_STRING_ELT( out_string, n + i, NA_STRING );
	  continue;
	}

      ma = CHAR( element );

      a = strchr( ma, '.' );
      if ( NULL == a )
	a = strchr( ma, '_' );

      if ( NULL == a )
	{
	  SET_STRING_ELT( out_string, i, element );
	  SET_STRING_ELT( out_string, n + i, mkChar( "" ) );
	}
      else
	{
	  SET_STRING_ELT(
			 out_string,
			 i,
			 mkCharLen( ma, (int)( a - ma ) )
			 );
	  SET_STRING_ELT( out_string, n + i, mkChar( a + 1 ) );
	}
    }
  UNPROTECT( 1 );
  return ( out_string );
}

void pq_genotype(
    int*nomar,
    int*noind,
    int*mgt,
    int*k
)
{
    int i, j;
    int*p_m = mgt;
    int*p_k = k;

    for ( j = 0;j < *nomar;j++ )
        for ( i = 0;i < *noind;i++ )
        {
            if ( ( *p_m == 0 ) && ( *( 1 + p_m ) == 1 ) )
                * p_k = 0;
            else if ( ( *p_m == 1 ) && ( *( 1 + p_m ) == 1 ) )
                * p_k = 1;
            else if ( ( *p_m == 1 ) && ( *( 1 + p_m ) == 0 ) )
                * p_k = 2;

            p_k++;

            p_m++;

            p_m++;
        }
}

/**************************************************************************
New code for MABC 2010 BEGIN
**************************************************************************/

void evaluate_ld ( 
		  char** PopName,
		  int* all,
		  int* chr,
		  double* pos,
		  double* ld
		  )

/* length of the linkage drag */
{
    TPopListEl*PopListEl;
    POPULATION ind;
    int chrom, hom, loc;
    double gLength = 0;
    TInd*Ind;
    double  upper, lower;

    (void)gLength;

    if ( ( *chr < 1 ) || ( *chr > NoChroms ) )
      {
	info( -2, "Chromosome out of range: %d\n", *chr );
	return;
      }
    
    find_list_element( *PopName, &PopListEl );

    if ( NULL == PopListEl )
      {
        info( -2, "Can't find Population %s\n", *PopName );
	return;
      }

    for ( chrom = 0;chrom < NoChroms;chrom++ )
        gLength += PopDescr[ chrom ].ChromLen * NoHoms;

    for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
    {
        Ind = get_individual( PopListEl->Pop, ind );

        chrom = *chr - 1;
        {
            for ( hom = 0;hom < NoHoms;hom++ )
            {
                for ( loc = 0;loc < Ind->Chrom[ chrom ].Hom[ hom ].NoLociUsed;loc++ )
                {
		  if ( *all == Ind->Chrom[ chrom ].Hom[ hom ].Loc[ loc ].All ) 
                    {
                        if ( loc + 1 < Ind->Chrom[ chrom ].Hom[ hom ].NoLociUsed )
                        {
			  upper = Ind->Chrom[ chrom ].Hom[ hom ].Loc[ loc + 1 ].Pos;
                        }
                        else
                        {
			  upper = PopDescr[ chrom ].ChromLen;
                        }
			lower = Ind->Chrom[ chrom ].Hom[ hom ].Loc[ loc ].Pos;
			if ( (lower <= *pos) && (*pos < upper) )
			  {
			    if (ld[ ind ] == 0) {         /* first homologue  */
			      ld[ ind ] = upper - lower;
			    } else {                      /* second homologue */
			      ld[ ind ] = ld[ ind ] / 2 + (upper - lower) / 2;
 			    }
			  }
                    }
                }
            }
        }
	/* ld[ ind ] = length[ ind ] ; */
    }
}


/**************************************************************************
New code for MABC 2010 END
**************************************************************************/



/**************************************************************************
Genetic distances BEGIN
**************************************************************************/

#include <R_ext/Random.h> 


static int AllowZeroFrequencies= 0;


void set_allow_zeros(int*allow){
  AllowZeroFrequencies= *allow;
}


void get_allow_zeros(int*allow){
  *allow= AllowZeroFrequencies;
}


SEXP splitdot(SEXP in_string){

  const char*ma,*a;
  int i;
  int n= length(in_string);
  SEXP element;
  SEXP out_string;

  PROTECT(out_string= allocMatrix(STRSXP,n,2));
  for(i= 0;i<n;i++){
    element= STRING_ELT(in_string,i);

    if(element==NA_STRING){
      SET_STRING_ELT(out_string,i,NA_STRING);
      SET_STRING_ELT(out_string,n+i,NA_STRING);
      continue;
    }

    ma= CHAR(element);

    a= strchr(ma,'.');
    if(NULL==a)
      a= strchr(ma,'_');

    if(NULL==a){
      SET_STRING_ELT(out_string,i,element);
      SET_STRING_ELT(out_string,n+i,mkChar(""));
    }
    else{
      SET_STRING_ELT(out_string,i,
		     mkCharLen(ma,(int)(a-ma)));
      SET_STRING_ELT(out_string,n+i,mkChar(a+1));
    }
  }
  UNPROTECT(1);
  return(out_string);
}


void allele_freq(
		 int*no_pop,
		 int*no_ind,
		 int*no_mar,
		 int*no_all,
		 int*no_ind_u,
		 int*ind_u,
		 double*data,
		 double*freq
		 )
{
  size_t*coff = NULL;
  size_t*roff = NULL;
  int i,j,k,l;
  size_t rowpos,colpos,nrow,ncol_d,data_count;
  double*ptr_data,*ptr_freq,*ptr_dataC;
  double*dataC = NULL;
  double count;

  if ( ( *no_pop < 1 ) || ( *no_mar < 1 ) ) goto fail;

  if ( ( size_t ) *no_pop > SIZE_MAX / sizeof( *coff ) ) goto fail;
  if ( ( size_t ) *no_mar > SIZE_MAX / sizeof( *roff ) ) goto fail;

  coff= (size_t*)malloc((size_t)(*no_pop)*sizeof(*coff));
  if ( NULL == coff ) goto fail;

  roff= (size_t*)malloc((size_t)(*no_mar)*sizeof(*roff));
  if ( NULL == roff ) goto fail;

  roff[0]= coff[0]= 0;

  for(k= 1;k<*no_pop;k++){
    if ( ( no_ind[k-1] < 0 ) ||
	 ( coff[k-1] > SIZE_MAX-(size_t)no_ind[k-1] ) )
      goto fail;
    coff[k]= coff[k-1]+(size_t)no_ind[k-1];
  }

  for(i= 1;i<*no_mar;i++){
    if ( ( no_all[i-1] < 0 ) ||
	 ( roff[i-1] > SIZE_MAX-(size_t)no_all[i-1] ) )
      goto fail;
    roff[i]= roff[i-1]+(size_t)no_all[i-1];
  }

  if ( ( no_all[*no_mar-1] < 0 ) ||
       ( roff[*no_mar-1] >
	 SIZE_MAX-(size_t)no_all[*no_mar-1] ) )
    goto fail;

  nrow= roff[*no_mar-1]+(size_t)no_all[*no_mar-1];

  if ( ( no_ind[*no_pop-1] < 0 ) ||
       ( coff[*no_pop-1] >
	 SIZE_MAX-(size_t)no_ind[*no_pop-1] ) )
    goto fail;

  ncol_d= coff[*no_pop-1]+(size_t)no_ind[*no_pop-1];

  if ( ( nrow != 0 ) && ( ncol_d > SIZE_MAX/nrow ) ) goto fail;
  data_count= nrow*ncol_d;

  if ( data_count > SIZE_MAX/sizeof(*dataC) ) goto fail;

  if ( ( nrow != 0 ) &&
       ( (size_t)*no_pop > SIZE_MAX/nrow ) )
    goto fail;

  if ( data_count > 0 ){
    dataC= (double*)malloc(data_count*sizeof(*dataC));
    if ( NULL == dataC ) goto fail;
  }

  for(k= 0;k<*no_pop;k++){
    for(l= 0;l<no_ind_u[k];l++){
      colpos=
	(coff[k]+
	 (size_t)(ind_u[coff[k]+(size_t)l]-1))*nrow;
      for(i= 0;i<*no_mar;i++){
	count= 0;
	for(j= 0;j<no_all[i];j++){
	  rowpos= roff[i]+(size_t)j;
	  ptr_data= &data[colpos+rowpos];
	  if(*ptr_data> 0){
	    *ptr_data= 1;
	    count+= 1;
	  }
	}
	for(j= 0;j<no_all[i];j++){
	  rowpos= roff[i]+(size_t)j;
	  ptr_data= &data[colpos+rowpos];
	  ptr_dataC= &dataC[colpos+rowpos];

	  if(AllowZeroFrequencies){*ptr_dataC= *ptr_data;}
	  else{
	    if(count> 0){
	      if(*ptr_data> 0)(*ptr_dataC)= (*ptr_data)/(double)count;
	      else(*ptr_dataC)= 0;
	    }
	    else{*ptr_dataC= -1;}
	  }
	}
      }
    }
  }

  for(i= 0;i<*no_mar;i++)
    for(j= 0;j<no_all[i];j++){
      rowpos= roff[i]+(size_t)j;
      for(k= 0;k<*no_pop;k++){
	ptr_freq= &freq[(size_t)k*nrow+rowpos];
	*ptr_freq= -1;
	count= 0;
	for(l= 0;l<no_ind_u[k];l++){
	  colpos=
	    (coff[k]+
	     (size_t)(ind_u[coff[k]+(size_t)l]-1))*nrow;
	  ptr_dataC= &dataC[colpos+rowpos];
	  if(*ptr_dataC!=-1){
	    count+= 1;
	    if(*ptr_freq==-1)(*ptr_freq)= 0;
	    (*ptr_freq)+= (*ptr_dataC);
	  }
	}
	if(count> 0){
	  if(*ptr_freq> 0)*ptr_freq= *ptr_freq/(double)count;
	  else*ptr_freq= 0;
	}
	else*ptr_freq= -1;
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

void correct_missing(
		     int*no_pop,
		     int*no_ind,
		     int*no_mar,
		     int*no_all,
		     int*no_ind_u,
		     int*ind_u,
		     double*data,
		     double*freq
		     )
{

  (void)freq;

  size_t*coff = NULL;
  size_t*roff = NULL;
  int i,j,k,l;
  size_t rowpos,colpos,nrow,ncol_d,data_count;
  double*ptr_data,*ptr_dataC;
  double*dataC = NULL;
  
  double count;

  if ( ( *no_pop < 1 ) || ( *no_mar < 1 ) ) goto fail;

  if ( ( size_t ) *no_pop > SIZE_MAX / sizeof( *coff ) ) goto fail;
  if ( ( size_t ) *no_mar > SIZE_MAX / sizeof( *roff ) ) goto fail;

  coff= (size_t*)malloc((size_t)(*no_pop)*sizeof(*coff));
  if ( NULL == coff ) goto fail;

  roff= (size_t*)malloc((size_t)(*no_mar)*sizeof(*roff));
  if ( NULL == roff ) goto fail;

  roff[0]= coff[0]= 0;

  for(k= 1;k<*no_pop;k++){
    if ( ( no_ind[k-1] < 0 ) ||
	 ( coff[k-1] > SIZE_MAX-(size_t)no_ind[k-1] ) )
      goto fail;
    coff[k]= coff[k-1]+(size_t)no_ind[k-1];
  }

  for(i= 1;i<*no_mar;i++){
    if ( ( no_all[i-1] < 0 ) ||
	 ( roff[i-1] > SIZE_MAX-(size_t)no_all[i-1] ) )
      goto fail;
    roff[i]= roff[i-1]+(size_t)no_all[i-1];
  }

  if ( ( no_all[*no_mar-1] < 0 ) ||
       ( roff[*no_mar-1] >
	 SIZE_MAX-(size_t)no_all[*no_mar-1] ) )
    goto fail;

  nrow= roff[*no_mar-1]+(size_t)no_all[*no_mar-1];

  if ( ( no_ind[*no_pop-1] < 0 ) ||
       ( coff[*no_pop-1] >
	 SIZE_MAX-(size_t)no_ind[*no_pop-1] ) )
    goto fail;

  ncol_d= coff[*no_pop-1]+(size_t)no_ind[*no_pop-1];

  if ( ( nrow != 0 ) && ( ncol_d > SIZE_MAX/nrow ) ) goto fail;
  data_count= nrow*ncol_d;

  if ( data_count > SIZE_MAX/sizeof(*dataC) ) goto fail;

  if ( data_count > 0 ){
    dataC= (double*)malloc(data_count*sizeof(*dataC));
    if ( NULL == dataC ) goto fail;
  }

  for(k= 0;k<*no_pop;k++){
    for(l= 0;l<no_ind_u[k];l++){
      colpos= (coff[k]+(size_t)(ind_u[coff[k]+l]-1))*nrow;
      for(i= 0;i<*no_mar;i++){
	count= 0;
	for(j= 0;j<no_all[i];j++){
	  rowpos= roff[i]+(size_t)j;
	  ptr_data= &data[colpos+rowpos];
	  if(*ptr_data> 0){
	    *ptr_data= 1;
	    count+= 1;
	  }
	}
	for(j= 0;j<no_all[i];j++){
	  rowpos= roff[i]+(size_t)j;
	  ptr_data= &data[colpos+rowpos];
	  ptr_dataC= &dataC[colpos+rowpos];
	  
	  if(AllowZeroFrequencies){*ptr_dataC= *ptr_data;}
	  else{
	    if(count> 0){
	      if(*ptr_data> 0)(*ptr_dataC)= (*ptr_data)/(double)count;
	      else(*ptr_dataC)= 0;
	    }
	    else{*ptr_dataC= -1;}
	  }
	}
      }
    }
  }

  for(k= 0;k<*no_pop;k++){
    for(l= 0;l<no_ind_u[k];l++){
      colpos= (coff[k]+(size_t)(ind_u[coff[k]+l]-1))*nrow;
      for(i= 0;i<*no_mar;i++){
	for(j= 0;j<no_all[i];j++){
	  rowpos= roff[i]+(size_t)j;
	  ptr_data= &data[colpos+rowpos];
	  ptr_dataC= &dataC[colpos+rowpos];
	  (*ptr_data)= (*ptr_dataC);
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

  info(-2, "Invalid dimensions or memory allocation failed in correct_missing\n");
  return;
}




void rd(
	int*no_pop,
	int*no_mar,
	int*no_all,
	int*no_mar_u,
	int*mar_u,
	double*freq,
	double*dist
	)
{
  
  long i,j,k1,k2,m,observed,nrow,rowpos;
  double t1,t2;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int*roff;

  if ( *no_mar < 1 )
    {
      info( -2, "Invalid number of markers in rd\n" );
      return;
    }

  if ( ( size_t ) *no_mar > SIZE_MAX / sizeof( *roff ) )
    {
      info( -2, "Memory allocation failed in rd\n" );
      return;
    }

  roff = ( int* ) malloc( ( size_t ) *no_mar * sizeof( *roff ) );

  if ( NULL == roff )
    {
      info( -2, "Memory allocation failed in rd\n" );
      return;
    }
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;

  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      t2= 0;m= 0;
      for(i= 0;i<*no_mar_u;i++){
	observed= 1;
	t1= 0;
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1))
	    t1+= pow((*ptr_freq1-*ptr_freq2),2);
	  else observed= 0;
	}
	if(observed){
	  t2+= sqrt(t1/2);
	  m+= 1;
	}
      }
      (*ptr_dist)= t2/(double)m;
      ptr_dist++;
    }

  free( roff );
}


void mrd(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,m,observed,nrow,rowpos;
  double t1;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int*roff;

  if ( *no_mar < 1 )
    {
      info( -2, "Invalid number of markers in mrd\n" );
      return;
    }

  if ( ( size_t ) *no_mar > SIZE_MAX / sizeof( *roff ) )
    {
      info( -2, "Memory allocation failed in mrd\n" );
      return;
    }

  roff = ( int* ) malloc( ( size_t ) *no_mar * sizeof( *roff ) );

  if ( NULL == roff )
    {
      info( -2, "Memory allocation failed in mrd\n" );
      return;
    }
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;
  
  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      t1= 0;m= 0;
      for(i= 0;i<*no_mar_u;i++){
	observed= 1;
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1))
	    t1+= pow((*ptr_freq1-*ptr_freq2),2);
	  else observed= 0;
	}
	if(observed)m+= 1;
      }
      (*ptr_dist)= sqrt(t1/2/(double)m);
      ptr_dist++;
    }

  free( roff );
}


void euc(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,m,observed,nrow,rowpos;
  double t1;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int*roff;

  if ( *no_mar < 1 )
    {
      info( -2, "Invalid number of markers in euc\n" );
      return;
    }

  if ( ( size_t ) *no_mar > SIZE_MAX / sizeof( *roff ) )
    {
      info( -2, "Memory allocation failed in euc\n" );
      return;
    }

  roff = ( int* ) malloc( ( size_t ) *no_mar * sizeof( *roff ) );

  if ( NULL == roff )
    {
      info( -2, "Memory allocation failed in euc\n" );
      return;
    }

  (void)m;
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;

  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      t1= 0;m= 0;
      for(i= 0;i<*no_mar_u;i++){
	observed= 1;
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1))
	    t1+= pow((*ptr_freq1-*ptr_freq2),2);
	  else observed= 0;
	}
	if(observed)m+= 1;
      }
      (*ptr_dist)= sqrt(t1);
      ptr_dist++;
    }

  free( roff );
}


void jac(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,nrow,rowpos;
  double zaehler,nenner,fa,fb;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int*roff;

  if ( *no_mar < 1 )
    {
      info( -2, "Invalid number of markers in jac\n" );
      return;
    }

  if ( ( size_t ) *no_mar > SIZE_MAX / sizeof( *roff ) )
    {
      info( -2, "Memory allocation failed in jac\n" );
      return;
    }

  roff = ( int* ) malloc( ( size_t ) *no_mar * sizeof( *roff ) );

  if ( NULL == roff )
    {
      info( -2, "Memory allocation failed in jac\n" );
      return;
    }
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;
  
  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      zaehler= nenner= 0;
      for(i= 0;i<*no_mar_u;i++){
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1)){
	    fa= (double)(*ptr_freq1> 0?1:0);
	    fb= (double)(*ptr_freq2> 0?1:0);
	    zaehler+= fa*fb;
	    nenner+= fa+fb-fa*fb;
	  }
	}
      }
      if(nenner> 0)
	(*ptr_dist)= zaehler/nenner;
      else
	(*ptr_dist)= 0;
      ptr_dist++;
    }

  free( roff );
}


void dic(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,nrow,rowpos;
  double zaehler,nenner,fa,fb;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int*roff;

  if ( *no_mar < 1 )
    {
      info( -2, "Invalid number of markers in dic\n" );
      return;
    }

  if ( ( size_t ) *no_mar > SIZE_MAX / sizeof( *roff ) )
    {
      info( -2, "Memory allocation failed in dic\n" );
      return;
    }

  roff = ( int* ) malloc( ( size_t ) *no_mar * sizeof( *roff ) );

  if ( NULL == roff )
    {
      info( -2, "Memory allocation failed in dic\n" );
      return;
    }
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;
  
  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      zaehler= nenner= 0;
      for(i= 0;i<*no_mar_u;i++){
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1)){
	    fa= (double)(*ptr_freq1> 0?1:0);
	    fb= (double)(*ptr_freq2> 0?1:0);
	    zaehler+= fa*fb;
	    nenner+= fa+fb;
	  }
	}
      }
      if(nenner> 0)
	(*ptr_dist)= 2*zaehler/nenner;
      else
	(*ptr_dist)= 0;
      ptr_dist++;
    }

  free( roff );
}


void sma(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,nrow,rowpos;
  double zaehler,nenner,fa,fb;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int*roff;

  if ( *no_mar < 1 )
    {
      info( -2, "Invalid number of markers in sma\n" );
      return;
    }

  if ( ( size_t ) *no_mar > SIZE_MAX / sizeof( *roff ) )
    {
      info( -2, "Memory allocation failed in sma\n" );
      return;
    }

  roff = ( int* ) malloc( ( size_t ) *no_mar * sizeof( *roff ) );

  if ( NULL == roff )
    {
      info( -2, "Memory allocation failed in sma\n" );
      return;
    }
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;
  
  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      zaehler= nenner= 0;
      for(i= 0;i<*no_mar_u;i++){
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1)){
	    fa= (double)(*ptr_freq1> 0?1:0);
	    fb= (double)(*ptr_freq2> 0?1:0);
	    zaehler+= fa*fb+(1-fa)*(1-fb);
	    nenner+= 1;
	  }
	}
      }
      if(nenner> 0)
	(*ptr_dist)= zaehler/nenner;
      else
	(*ptr_dist)= 0;
      ptr_dist++;
    }

  free( roff );
}

void gen_dist(
	      int*measure,
	      int*no_pop,
	      int*no_ind,
	      int*no_mar,
	      int*no_all,
	      int*no_pop_u,
	      int*pop_u,
	      int*no_ind_u,
	      int*ind_u,
	      int*no_mar_u,
	      int*mar_u,
	      double*data,
	      double*freq,
	      double*dist,
	      int*s,
	      double*sdev,
	      char**par
	      )
{

  void (*distance)(int*, int*, int*, int*, int*, double*, double*) = NULL;
  double *dist_j = NULL, *dist_b = NULL;
  double *sum_dist_sq = NULL, *sum_dist = NULL;
  int *mar_u_j = NULL, *ind_u_boot = NULL, *mar_u_boot = NULL;

  (void)no_pop_u;
  (void)pop_u;

  if ( ( *no_pop < 1 ) || ( *no_mar < 1 ) ) goto fail;

  switch(*measure){
    case 1:  distance = rd;  break;
    case 2:  distance = mrd; break;
    case 3:  distance = euc; break;
    case 11: distance = dic; break;
    case 12: distance = jac; break;
    case 13: distance = sma; break;
    default: distance = euc;
    }

  allele_freq(no_pop,no_ind,no_mar,no_all,no_ind_u,ind_u,data,freq);
  (*distance)(no_pop,no_mar,no_all,no_mar_u,mar_u,freq,dist);

  if(*s==1){

    int i,j,k,n_dist;
    size_t n_dist_size,n_dist_a,n_dist_b;
    double f= (double)(*no_mar-1)/(double)(*no_mar);
    int no_mar_u[1];
    *no_mar_u= *no_mar-1;

    if ( ( *no_pop < 2 ) || ( *no_mar_u < 1 ) ) goto fail;

    n_dist_a = ( size_t ) *no_pop;
    n_dist_b = ( size_t ) ( *no_pop - 1 );
    if ( 0 == n_dist_a % 2 ) n_dist_a /= 2;
    else n_dist_b /= 2;
    if ( n_dist_a > SIZE_MAX / n_dist_b ) goto fail;
    n_dist_size = n_dist_a * n_dist_b;
    if ( n_dist_size > INT_MAX ) goto fail;
    n_dist = ( int ) n_dist_size;

    {
      dist_j = ( double* ) yRealloc_impl( NULL, n_dist_size, sizeof( double ) );
      if ( NULL == dist_j ) goto fail;
      sum_dist_sq = ( double* ) yRealloc_impl( NULL, n_dist_size, sizeof( double ) );
      if ( NULL == sum_dist_sq ) goto fail;
      sum_dist = ( double* ) yRealloc_impl( NULL, n_dist_size, sizeof( double ) );
      if ( NULL == sum_dist ) goto fail;
      mar_u_j = ( int* ) yRealloc_impl( NULL, ( size_t ) *no_mar_u, sizeof( int ) );
      if ( NULL == mar_u_j ) goto fail;
      
      for(j= 0;j<n_dist;j++)
	sum_dist_sq[j]= sum_dist[j]= 0;
      
      for(i= 0;i<*no_mar;i++){
	for(j= k= 0;j<*no_mar;j++)
	  if(i!=j){
	    mar_u_j[k]= j+1;
	    k++;
	  }
	for(j= 0;j<n_dist;j++)dist_j[j]= 0;

	(*distance)(no_pop,no_mar,no_all,no_mar_u,mar_u_j,freq,dist_j);

	for(j= 0;j<n_dist;j++){
	  sum_dist_sq[j]+= pow(dist_j[j],2);
	  sum_dist[j]+= dist_j[j];
	}

      }

      for(j= 0;j<n_dist;j++)
	sdev[j] = sqrt(f*(sum_dist_sq[j]-pow(sum_dist[j],2) /
			  (double)(*no_mar)));
      
      yFree( dist_j );
      yFree( sum_dist_sq );
      yFree( sum_dist );
      yFree( mar_u_j );
    }
  }

  if(*s> 1){

    int i,j,k,l,lfdnr,n_dist,sumall= 0;
    size_t n_dist_size,n_dist_a,n_dist_b,sumind= 0;
    double f= (double)1/(double)(*s);

    if ( ( *no_pop < 2 ) || ( *no_mar_u < 1 ) ) goto fail;

    for (i = 0; i<*no_pop; i++)
      {
	if ( ( no_ind[i] < 0 ) ||
	     ( ( size_t ) no_ind[i] > SIZE_MAX - sumind ) )
	  goto fail;
	sumind += ( size_t ) no_ind[i];
      }

    if ( ( 0 == sumind ) || ( sumind > INT_MAX ) ) goto fail;

    n_dist_a = ( size_t ) *no_pop;
    n_dist_b = ( size_t ) ( *no_pop - 1 );
    if ( 0 == n_dist_a % 2 ) n_dist_a /= 2;
    else n_dist_b /= 2;
    if ( n_dist_a > SIZE_MAX / n_dist_b ) goto fail;
    n_dist_size = n_dist_a * n_dist_b;
    if ( n_dist_size > INT_MAX ) goto fail;
    n_dist = ( int ) n_dist_size;

    (void)sumall;
    
    {

      ind_u_boot = ( int* ) yRealloc_impl( NULL, sumind, sizeof( int ) );
      if ( NULL == ind_u_boot ) goto fail;
      mar_u_boot = ( int* ) yRealloc_impl( NULL, ( size_t ) *no_mar_u, sizeof( int ) );
      if ( NULL == mar_u_boot ) goto fail;
      dist_b = ( double* ) yRealloc_impl( NULL, n_dist_size, sizeof( double ) );
      if ( NULL == dist_b ) goto fail;
      sum_dist_sq = ( double* ) yRealloc_impl( NULL, n_dist_size, sizeof( double ) );
      if ( NULL == sum_dist_sq ) goto fail;
      sum_dist = ( double* ) yRealloc_impl( NULL, n_dist_size, sizeof( double ) );
      if ( NULL == sum_dist ) goto fail;
      
      for(j= 0;j<n_dist;j++)
	sum_dist_sq[j]= sum_dist[j]= 0;
      
      for(i= 0;i<*s;i++){
	if(strstr(*par,"m")){
	  sumall= 0;
	  for(k= 0;k<*no_mar_u;k++){
	    mar_u_boot[k] = 1 + st_rand_index(*no_mar_u); // Using MathLib12
	    sumall+= no_all[mar_u_boot[k]-1];
	  }
	}
	else{
	  sumall= 0;
	  for(k= 0;k<*no_mar_u;k++){
	    mar_u_boot[k]= mar_u[k];
	    sumall+= no_all[mar_u_boot[k]-1];
	  }
	}
	
	if(strstr(*par,"i")){
	  lfdnr= 0;
	  for(k= 0;k<*no_pop;k++){
	    for(l= 0;l<no_ind[k];l++){
	      ind_u_boot[lfdnr] = 1 + st_rand_index(no_ind[k]);// Using MathLib12
	      lfdnr++;
	    }
	  }
	}
	else{
	  lfdnr= 0;
	  for(k= 0;k<*no_pop;k++){
	    for(l= 0;l<no_ind[k];l++){
	      ind_u_boot[lfdnr]= ind_u[lfdnr];
	      lfdnr++;
	    }
	  }
	}

	for(j= 0;j<n_dist;j++)dist_b[j]= 0;

	allele_freq(no_pop,no_ind,no_mar,no_all,no_ind_u,ind_u_boot,data,freq);
	(*distance)(no_pop,no_mar,no_all,no_mar_u,mar_u_boot,freq,dist_b);

	for(j= 0;j<n_dist;j++){
	  sum_dist_sq[j]+= pow(dist_b[j],2);
	  sum_dist[j]+= dist_b[j];
	}

      }

      for(j= 0;j<n_dist;j++)
	sdev[j]= sqrt(f*(sum_dist_sq[j]-pow(sum_dist[j],2)/(double)(*s)));
      
      yFree( ind_u_boot );
      yFree( mar_u_boot );
      yFree( dist_b );
      yFree( sum_dist_sq );
      yFree( sum_dist );
    }
  }

  return;

fail:

  yFree( dist_j );
  yFree( dist_b );
  yFree( sum_dist_sq );
  yFree( sum_dist );
  yFree( mar_u_j );
  yFree( ind_u_boot );
  yFree( mar_u_boot );

  info( -2, "Invalid dimensions or memory allocation failed in gen_dist\n" );
  return;
}

/**************************************************************************
Genetic distances END
**************************************************************************/

/**************************************************************************
Interface to Simulation Routines / Part 1: get_simpop
**************************************************************************/

#define gs_MAXSTR 4096
#define gs_MAXNME 256

typedef unsigned long int  gs_UNSIGNED ;
typedef long double        gs_FLOAT ;

typedef struct
{
    int          chrom;
    double       pos;
    const char  *name;
    const char  *class;
    gs_UNSIGNED  nr;
} gs_mappoint_type;

static int ps_get_simpop_name_valid(const char *name)
{
  const unsigned char *p = (const unsigned char *)name;

  if ((NULL == name) || ('\0' == name[0]))
    return 0;

  for (; '\0' != *p; p++)
    if (!(('A' <= *p && *p <= 'Z') ||
          ('a' <= *p && *p <= 'z') ||
          ('0' <= *p && *p <= '9')))
      return 0;

  return 1;
}


static int ps_get_simpop_population_name_valid(const char *name)
{
  size_t n;

  if (!ps_get_simpop_name_valid(name))
    return 0;

  n = strlen(name);

  return ((2 <= n) && (n < MAXNAME));
}


static int ps_get_simpop_allele_pair(
    unsigned short a1,
    unsigned short a2,
    unsigned short bg,
    GS_ALLELE_TYPE *out1,
    GS_ALLELE_TYPE *out2)
{
  if ((bg == a1) || (bg == a2) ||
      (USHRT_MAX == a1) || (USHRT_MAX == a2))
    {
      *out1 = GS_ALLELE_MISSING;
      *out2 = GS_ALLELE_MISSING;
      return 1;
    }

  if ((0 == a1) || (0 == a2) ||
      (a1 > GS_ALLELE_MAX) || (a2 > GS_ALLELE_MAX))
    return 0;

  *out1 = (GS_ALLELE_TYPE)a1;
  *out2 = (GS_ALLELE_TYPE)a2;

  return 1;
}


void ps_get_simpop ( char ** pop_name   ,
		     gs_UNSIGNED * NoInd,
		     gs_UNSIGNED * NoMar,
		     char (**IndName) [gs_MAXNME]  ,
		     char (**MarName) [gs_MAXNME]  ,
		     GS_ALLELE_TYPE *(* All1)  ,
		     GS_ALLELE_TYPE *(* All2)  ,
		     gs_UNSIGNED       *  NLoc ,
		     gs_mappoint_type  ** map_S,
                     char              ** map_str_S,
                     size_t             * map_str_N,
		     gs_mappoint_type *** map  ,
		     gs_UNSIGNED       ** mdl  ,
		     int* bg,
		     int* Sretval         )
{

  gs_UNSIGNED i,m;
  size_t data_count;
  size_t map_pool_size = 0;
  size_t map_offset = 0;
  int c,l;
  int allele_error = 0;
  unsigned short bg_allele;
  char tmp_str[2*gs_MAXNME];

  if ((*bg < 0) || (*bg > USHRT_MAX))
    {
      info(-2,"SIM missing allele indicator is invalid\n");
      *Sretval = -2;
      return;
    }
  bg_allele = (unsigned short)*bg;

  /* Find population */

  TPopListEl *PopListEl;
  TPop *Pop = NULL; 

  find_list_element( *pop_name, &PopListEl );

  if ( PopListEl == NULL ) 
    {
      info( -2, "SIM Population is not defined\n" );
      *Sretval = -2;  return;
    }

  Pop = PopListEl->Pop;

  if ( !ps_get_simpop_population_name_valid( Pop->Name ) )
    {
      info( -2,
            "SIM population name must contain 2 to 255 characters "
            "from A-Z, a-z, and 0-9\n" );
      *Sretval = -2;  return;
    }
  
  if ( 1 > Pop->NoInds )
    {
      info( -2, "SIM Population is empty\n " );
      *Sretval = -2;  return;
    }

  if ( 2 != NoHoms )
    {
      info( -2,
            "SIM transfer to GS data set implemented for diploids only\n" );
      *Sretval = -2;  return;
    }

  /* Indidivuals */

  *NoInd  = Pop->NoInds;

  *IndName = ( char (*)[gs_MAXNME] ) yRealloc_impl
    ( *IndName, *NoInd, sizeof(char[gs_MAXNME]) );
  if ( NULL == *IndName ) goto fail;

  for (i = 0 ; i< *NoInd; i++) 
    {
      snprintf(tmp_str, sizeof(tmp_str),
         "%sI%lu", Pop->Name, 1+i);
      if ( strlen( tmp_str ) >= gs_MAXNME )
	{
	  info( -2, "SIM Individual name too long\n" );
	  goto fail;
	}
      strcpy( (*IndName)[i] ,tmp_str  );
    }

  /* Markers */ 

  *NoMar = cMP;
  *MarName = ( char (*)[gs_MAXNME] ) yRealloc_impl
    ( *MarName, *NoMar, sizeof(char[gs_MAXNME]) );
  if ( NULL == *MarName ) goto fail;

  if ( NULL == PopDescr ) 
    {
      info(-2,"SIM No linkage map loaded\n");
      goto fail;
    }

  m = 0;
  for ( c=0; c < NoChroms; c++ ) 
    {
      for ( l = 0; l < PopDescr[c].NrMappoints; l++ )
	{
	  if ( !ps_get_simpop_name_valid(PopDescr[c].Mappoint[l].name) ||
	       strlen(PopDescr[c].Mappoint[l].name) >= gs_MAXNME ) {
	    info(-2,"SIM marker name invalid for SelectionTools transfer\n");
	    goto fail;
	  }
	  strcpy ( (*MarName)[m] , PopDescr[c].Mappoint[l].name );
	  m++;  
	}
    }

  /* Allocate memory for the genotypic data */

  if ( ( size_t ) *NoMar > SIZE_MAX / ( size_t ) *NoInd )
    goto fail;

  data_count = ( size_t ) *NoMar * ( size_t ) *NoInd;

  *All1 = (GS_ALLELE_TYPE*)yRealloc(*All1,data_count,GS_ALLELE_TYPE);
  if ( NULL == *All1 ) goto fail;
  *All2 = (GS_ALLELE_TYPE*)yRealloc(*All2,data_count,GS_ALLELE_TYPE);
  if ( NULL == *All2 ) goto fail;

  /* "Genotype population" to a ST data set */

#ifdef MT
  #pragma omp parallel 
#endif

  {

    POPULATION ind;
    int  chrom;

    TInd   *Ind;
    TChrom *Chrom;
    int    i, m ;

    size_t idx;

#ifdef MT
  #pragma omp for 
#endif
  
    for ( ind = 0; ind < Pop->NoInds; ind++ )
      {
	Ind = get_individual( PopListEl->Pop, ind );

	for ( chrom=0,m=0 ; chrom < NoChroms; chrom++ )
	  {
	    Chrom = &( Ind->Chrom[chrom] );

	    for ( i = 0; i < PopDescr[chrom].NrMappoints; i++,m++ )
	      {

		idx = ( size_t ) m * ( size_t ) *NoInd +
		      ( size_t ) ind;

                {
                  unsigned short a1 = get_allele(&(Chrom->Hom[0]),
                                                  PopDescr[chrom].Mappoint[i].Pos);
                  unsigned short a2 = get_allele(&(Chrom->Hom[1]),
                                                  PopDescr[chrom].Mappoint[i].Pos);

                  if (!ps_get_simpop_allele_pair(a1,a2,bg_allele,
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
  } //parallel

  if (allele_error)
    {
      info(-2,"SIM allele code cannot be represented in GS marker data\n");
      goto fail;
    }

  /* Linkage map */

  *NLoc = cMP;

  for ( c=0; c < NoChroms; c++ )
    {
      for ( l=0; l < PopDescr[c].NrMappoints; l++ )
        {
          const char *name = PopDescr[c].Mappoint[l].name;
          const char *class_name = PopDescr[c].Mappoint[l].class;
          size_t name_len;
          size_t class_len;

          if ((NULL == name) || (NULL == class_name))
            {
              info( -2, "SIM Marker name or class invalid\n" );
              goto fail;
            }

          name_len = strlen(name);
          class_len = strlen(class_name);

          if ( name_len >= gs_MAXSTR )
            {
              info( -2, "SIM Marker name too long\n" );
              goto fail;
            }

          if ( class_len >= gs_MAXSTR )
            {
              info( -2, "SIM Marker class too long\n" );
              goto fail;
            }

          if ((name_len + 1 > SIZE_MAX - map_pool_size) ||
              (class_len + 1 >
               SIZE_MAX - (map_pool_size + name_len + 1)))
            {
              info( -2, "SIM linkage map string storage too large\n" );
              goto fail;
            }

          map_pool_size += name_len + 1;
          map_pool_size += class_len + 1;
        }
    }

  *map_S = ( gs_mappoint_type* ) yRealloc( *map_S, *NLoc,
                                           gs_mappoint_type );
  if ( NULL == *map_S ) goto fail;
  *map_str_S = ( char* ) yRealloc( *map_str_S, map_pool_size, char );
  if ( NULL == *map_str_S ) goto fail;
  *map_str_N = map_pool_size;
  *map = ( gs_mappoint_type** ) yRealloc( *map, *NLoc,
                                          gs_mappoint_type* );
  if ( NULL == *map ) goto fail;
  *mdl = ( gs_UNSIGNED* ) yRealloc( *mdl, *NLoc, gs_UNSIGNED );
  if ( NULL == *mdl ) goto fail;

  for ( c=0,m=0 ; c < NoChroms; c++ )
    {
      for ( l=0; l < PopDescr[c].NrMappoints; l++,m++ )
        {
          const char *name = PopDescr[c].Mappoint[l].name;
          const char *class_name = PopDescr[c].Mappoint[l].class;
          size_t name_len = strlen(name) + 1;
          size_t class_len = strlen(class_name) + 1;

          (*map_S)[m].chrom = (int)PopDescr[c].Mappoint[l].chrom;
          (*map_S)[m].pos = (double)PopDescr[c].Mappoint[l].Pos * 100;

          memcpy(*map_str_S + map_offset,name,name_len);
          (*map_S)[m].name = *map_str_S + map_offset;
          map_offset += name_len;

          memcpy(*map_str_S + map_offset,class_name,class_len);
          (*map_S)[m].class = *map_str_S + map_offset;
          map_offset += class_len;
        }
    }

  if (map_offset != *map_str_N)
    goto fail;

  /* Array of pointers to mappoints and index in the marker matrix*/
  for ( m=0 ; m < *NLoc ; m++ )
    {
      (*map)[m] = &((*map_S)[m]);
      (*mdl)[m] = m;
      (*map_S)[m].nr = m;
    }

  *Sretval = 0;  return;

fail:

  yFree( *IndName );
  yFree( *MarName );
  yFree( *All1 );
  yFree( *All2 );
  yFree( *map_S );
  yFree( *map_str_S );
  *map_str_N = 0;
  yFree( *map );
  yFree( *mdl );

  *NoInd = 0;
  *NoMar = 0;
  *NLoc = 0;
  *Sretval = -2;  return;

}


//+DEBUG++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ps_get_simpop2 ( char ** pop_name   ,
		     gs_UNSIGNED * NoInd,
		     gs_UNSIGNED * NoMar,
		     char (**IndName) [gs_MAXNME]  ,
		     char (**MarName) [gs_MAXNME]  ,
		     GS_ALLELE_TYPE *(* All1)  ,
		     GS_ALLELE_TYPE *(* All2)  ,
		     gs_UNSIGNED       *  NLoc ,
		     gs_mappoint_type  ** map_S,
		     gs_mappoint_type *** map  ,
		     gs_UNSIGNED       ** mdl  ,
		     int* bg,
		     int* Sretval         )
// No linkage map
{
  (void)NLoc;
  (void)map_S;
  (void)map;
  (void)mdl;
  
  gs_UNSIGNED i,m;
  size_t data_count;
  int c,l;
  int allele_error = 0;
  unsigned short bg_allele;
  char tmp_str[2*gs_MAXNME];

  if ((*bg < 0) || (*bg > USHRT_MAX))
    {
      info(-2,"SIM missing allele indicator is invalid\n");
      *Sretval = -2;
      return;
    }
  bg_allele = (unsigned short)*bg;

  /* Find population */

  TPopListEl *PopListEl;
  TPop *Pop = NULL; 

  find_list_element( *pop_name, &PopListEl );

  if ( PopListEl == NULL ) 
    {
      info( -2, "SIM Population is not defined\n" );
      *Sretval = -2;  return;
    }

  Pop = PopListEl->Pop;

  if ( !ps_get_simpop_population_name_valid( Pop->Name ) )
    {
      info( -2,
            "SIM population name must contain 2 to 255 characters "
            "from A-Z, a-z, and 0-9\n" );
      *Sretval = -2;  return;
    }
  
  if ( 1 > Pop->NoInds )
    {
      info( -2, "SIM Population is empty\n " );
      *Sretval = -2;  return;
    }

  if ( 2 != NoHoms )
    {
      info( -2,
            "SIM transfer to GS data set implemented for diploids only\n" );
      *Sretval = -2;  return;
    }

  /* Indidivuals */

  *NoInd  = Pop->NoInds;

  *IndName = ( char (*)[gs_MAXNME] ) yRealloc_impl
    ( *IndName, *NoInd, sizeof(char[gs_MAXNME]) );
  if ( NULL == *IndName ) goto fail;

  for (i = 0 ; i< *NoInd; i++) 
    {
      snprintf(tmp_str, sizeof(tmp_str), "%sI%lu", Pop->Name,
	       (unsigned long)(1+i));
      if ( strlen( tmp_str ) >= gs_MAXNME )
        {
          info( -2, "SIM Individual name too long\n" );
          goto fail;
        }
      strcpy( (*IndName)[i] ,tmp_str  );
    }

  /* Markers */ 

  *NoMar = cMP;
  *MarName = ( char (*)[gs_MAXNME] ) yRealloc_impl
    ( *MarName, *NoMar, sizeof(char[gs_MAXNME]) );
  if ( NULL == *MarName ) goto fail;

  if ( NULL == PopDescr ) 
    {
      info(-2,"SIM No linkage map loaded\n");
      goto fail;
    }

  m = 0;
  for ( c=0; c < NoChroms; c++ ) 
    {
      for ( l = 0; l < PopDescr[c].NrMappoints; l++ )
	{
	  if ( !ps_get_simpop_name_valid(PopDescr[c].Mappoint[l].name) ||
	       strlen(PopDescr[c].Mappoint[l].name) >= gs_MAXNME ) {
	    info(-2,"SIM marker name invalid for SelectionTools transfer\n");
	    goto fail;
	  }
	  strcpy ( (*MarName)[m] , PopDescr[c].Mappoint[l].name );
	  m++;  
	}
    }

  /* Allocate memory for the genotypic data */

  if ( ( size_t ) *NoMar > SIZE_MAX / ( size_t ) *NoInd )
    goto fail;

  data_count = ( size_t ) *NoMar * ( size_t ) *NoInd;

  *All1 = (GS_ALLELE_TYPE*)yRealloc(*All1,data_count,GS_ALLELE_TYPE);
  if ( NULL == *All1 ) goto fail;
  *All2 = (GS_ALLELE_TYPE*)yRealloc(*All2,data_count,GS_ALLELE_TYPE);
  if ( NULL == *All2 ) goto fail;

  /* "Genotype population" to a ST data set */

#ifdef MT
  #pragma omp parallel 
#endif

  {

    POPULATION ind;
    int  chrom;

    TInd   *Ind;
    TChrom *Chrom;
    int    i, m ;

    gs_UNSIGNED idx;

#ifdef MT
  #pragma omp for 
#endif
  
    for ( ind = 0; ind < Pop->NoInds; ind++ )
      {
	Ind = get_individual( PopListEl->Pop, ind );

	for ( chrom=0,m=0 ; chrom < NoChroms; chrom++ )
	  {
	    Chrom = &( Ind->Chrom[chrom] );

	    for ( i = 0; i < PopDescr[chrom].NrMappoints; i++,m++ )
	      {

		idx = m * *NoInd + ind;

                {
                  unsigned short a1 = get_allele(&(Chrom->Hom[0]),
                                                  PopDescr[chrom].Mappoint[i].Pos);
                  unsigned short a2 = get_allele(&(Chrom->Hom[1]),
                                                  PopDescr[chrom].Mappoint[i].Pos);

                  if (!ps_get_simpop_allele_pair(a1,a2,bg_allele,
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
  } //parallel

  if (allele_error)
    {
      info(-2,"SIM allele code cannot be represented in GS marker data\n");
      goto fail;
    }

  *Sretval = 0;  return;

fail:

  yFree( *IndName );
  yFree( *MarName );
  yFree( *All1 );
  yFree( *All2 );

  *NoInd = 0;
  *NoMar = 0;
  *Sretval = -2;  return;

}

//+DEBUG++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/**************************************************************************
Interface to Simulation Routines / Part 1: get_simpop
**************************************************************************/


/**************************************************************************
Begin: Interface to simulation routines / Part2 set_sim            Jun 2019
**************************************************************************/


static int st_copy_sim_mappoints(
    gs_UNSIGNED        NLoc,
    gs_mappoint_type **map,
    TMappoint         **points,
    int                *nMP )
{
  gs_UNSIGNED lne;

  *points = NULL;
  *nMP = 0;

  if ((0 == NLoc) ||
      (NLoc > (gs_UNSIGNED)INT_MAX) ||
      (NULL == map))
    {
      info(-2,"Invalid GS linkage map\n");
      return 0;
    }

  for (lne=0; lne < NLoc; lne++)
    {
      if (NULL == map[lne])
        {
          info(-2,"Invalid GS linkage map pointer\n");
          goto fail;
        }

      if ((NULL == map[lne]->name) ||
          (NULL == map[lne]->class))
        {
          info(-2,"Invalid GS linkage map name or class\n");
          goto fail;
        }

      if (!append_mappoint(points,
                           nMP,
                           map[lne]->chrom,
                           map[lne]->pos / 100.0,
                           map[lne]->name,
                           map[lne]->class))
        {
          info(-2,
               "Memory allocation failed in st_copy_sim_mappoints\n");
          goto fail;
        }
    }

  return 1;

fail:

  free_mappoints(*points,*nMP);
  *points = NULL;
  *nMP = 0;

  return 0;
}


static int st_sim_allele(
    GS_ALLELE_TYPE  allele,
    unsigned short *sim_allele )
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


void st_set_sim_gp_provide_data ( int*      n,           // Number of chromosomes
				  double ** chrlen,      // Chromosome lengths
			          int*      retval,      // Return value
			          char**    set_name     // Name of the data set
				  ) ;

void st_set_sim_gp              ( int*     retval,      // Return value
			          char**   set_name     // Name of the data set
			         )
{
  int n = 0;
  int Sretval = -2;
  int hom = 2;
  double *chrlen = NULL;

  *retval = -2;

  st_set_sim_gp_provide_data(
    &n,
    &chrlen,
    &Sretval,
    set_name );

  if ((0 != Sretval) ||
      (n < 1) ||
      (NULL == chrlen))
    {
      free(chrlen);
      return;
    }

  if (!set_genome_par_impl(
      &n,
      &hom,
      chrlen))
    {
      free(chrlen);
      return;
    }

  free(chrlen);

  *retval = 0;
  return;
}


void st_set_sim_mp_provide_data (  gs_UNSIGNED       *  NLoc, // Number of loci
                                   gs_mappoint_type *** map , // Map
                                   int*     retval,      // Return value
 			           char**   set_name     // Name of the data set
				   );


void st_set_sim_mp(int *retval,
                   char **set_name)
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
      info(-2,"First remove previously loaded linkage map\n");
      return;
    }

  if (NULL == PopDescr)
    {
      info(-2,"No genome parameters defined\n");
      return;
    }

  st_set_sim_mp_provide_data(
    &NLoc,
    &map,
    &Sretval,
    set_name );

  if (0 != Sretval)
    return;

  if (!st_copy_sim_mappoints(
      NLoc,
      map,
      &points,
      &nMP))
    return;

  if (!prepare_mappoints(
      points,
      nMP,
      1,
      MAP_DISPERSE_FACTOR_DEFAULT,
      &modified))
    goto fail;

  if (!install_mappoints(points,nMP))
    {
      info(-2,
           "Memory allocation failed in st_set_sim_mp\n");
      goto fail;
    }

  free_mappoints(points,nMP);

  if (0<modified) info(0,
       "%d loci positions modified\n",
       modified);

  info(0,
       "%d loci defined\n",
       nMP);

  *retval = 0;
  return;

fail:

  free_mappoints(points,nMP);
}


void st_set_sim_ef_provide_data (  gs_UNSIGNED       *  NoEff,
                                   char             (** EffNme)[gs_MAXNME],
                                   int              ** EffAll,
                                   int              ** CmpAll,
                                   double           ** u,
                                   int               * retval,
                                   char             ** set_name );


static void st_sim_ef_free(TEffects *Effect)
{
  int i,j;

  if (NULL == Effect)
    return;

  if (NULL != Effect->Eff)
    {
      for (i=0; i < Effect->NrEff; i++)
        {
          if (NULL != Effect->Eff[i].Loc)
            {
              for (j=0; j < Effect->Eff[i].NrLoci; j++)
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
  int i,j,k,l,m,n;
  int maxLoc = 0;

  Effect->MPs = 0;
  free(Effect->MPset);
  Effect->MPset = NULL;

  for (i=0; i < Effect->NrEff; i++)
    {
      if (maxLoc < Effect->Eff[i].NrLoci)
        maxLoc = Effect->Eff[i].NrLoci;

      for (j=0; j < Effect->Eff[i].NrLoci; j++)
        {
          TMappoint **tmp;

          for (k=0; k < Effect->MPs; k++)
            if (Effect->MPset[k] == Effect->Eff[i].Loc[j].MP)
              break;

          if (k != Effect->MPs)
            continue;

          if (INT_MAX == Effect->MPs)
            return 0;

          tmp = realloc(Effect->MPset,
                        (size_t)(Effect->MPs + 1) *
                        sizeof(*Effect->MPset));

          if (NULL == tmp)
            return 0;

          Effect->MPset = tmp;
          Effect->MPset[Effect->MPs] = Effect->Eff[i].Loc[j].MP;
          Effect->MPs++;
        }
    }

  for (i=0; i < maxLoc; i++)
    {
      for (k=0; k < NoHoms; k++)
        {
          int oldEff = Effect->NrEff;

          for (n=k+1; n < NoHoms; n++)
            {
              for (j=0; j < oldEff; j++)
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

                  tmp = realloc(Effect->Eff,
                                (size_t)(Effect->NrEff + 1) *
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

                      newEff->Loc =
                        malloc((size_t)newEff->NrLoci *
                               sizeof(*newEff->Loc));

                      if (NULL == newEff->Loc)
                        return 0;

                      for (l=0; l < newEff->NrLoci; l++)
                        {
                          newEff->Loc[l].MP = NULL;
                          newEff->Loc[l].NrAll = 0;
                          newEff->Loc[l].All = NULL;
                        }
                    }

                  for (l=0; l < newEff->NrLoci; l++)
                    {
                      newEff->Loc[l].MP = Effect->Eff[j].Loc[l].MP;
                      newEff->Loc[l].NrAll =
                        Effect->Eff[j].Loc[l].NrAll;

                      if ((NoHoms < 1) ||
                          ((size_t)NoHoms >
                           SIZE_MAX / sizeof(*newEff->Loc[l].All)))
                        return 0;

                      newEff->Loc[l].All =
                        malloc((size_t)NoHoms *
                               sizeof(*newEff->Loc[l].All));

                      if (NULL == newEff->Loc[l].All)
                        return 0;

                      for (m=0; m < NoHoms; m++)
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


void st_set_sim_ef(int *retval,
                   char **effect_name,
                   char **set_name)
{
  gs_UNSIGNED NoEff = 0;
  char (*EffNme)[gs_MAXNME] = NULL;
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

  st_set_sim_ef_provide_data(
    &NoEff,
    &EffNme,
    &EffAll,
    &CmpAll,
    &u,
    &Sretval,
    set_name );

  if (0 != Sretval)
    return;

  /*
   * Marker effects are optional.  A data set without estimated effects
   * leaves the simulation unchanged and is not an error.
   */
  if (0 == NoEff)
    {
      *retval = 0;
      return;
    }

  if ((NULL == effect_name) ||
      (NULL == *effect_name) ||
      ('\0' == (*effect_name)[0]) ||
      (strlen(*effect_name) >= sizeof(Effect.Name)))
    {
      info(-2,"Invalid simulation effect name\n");
      return;
    }

  if (NULL != find_Effect(*effect_name))
    {
      info(-2,"This effect is already loaded (same name)\n");
      return;
    }

  if (2 != NoHoms)
    {
      info(-2,
           "GS effect transfer is implemented for diploids only\n");
      return;
    }

  if ((NULL == PopDescr) || (cMP < 1))
    {
      info(-2,"No simulation linkage map defined\n");
      return;
    }

  if ((NULL == EffNme) ||
      (NULL == EffAll) ||
      (NULL == CmpAll) ||
      (NULL == u))
    {
      info(-2,"Invalid marker effect data\n");
      return;
    }

  memset(&Effect,0,sizeof(Effect));
  gvalue_avg = u[0];
  Effect.GValue_avg = 0;
  Effect.weight = 1;
  strcpy(Effect.Name,*effect_name);

  nmarker = NoEff - 1;

  if (nmarker > (gs_UNSIGNED)(INT_MAX / 2))
    {
      info(-2,"Too many marker effects\n");
      return;
    }

  if (0 < nmarker)
    {
      int baseEff = 2 * (int)nmarker;

      if ((size_t)baseEff >
          SIZE_MAX / sizeof(*Effect.Eff))
        {
          info(-2,"Memory allocation failed in st_set_sim_ef\n");
          return;
        }

      Effect.Eff =
        malloc((size_t)baseEff * sizeof(*Effect.Eff));

      if (NULL == Effect.Eff)
        {
          info(-2,"Memory allocation failed in st_set_sim_ef\n");
          return;
        }

      Effect.NrEff = 0;
    }

  for (i=1; i < NoEff; i++)
    {
      TMappoint *tmpMP;
      int e;

      tmpMP = find_mappoint(EffNme[i]);

      if (NULL == tmpMP)
        {
          info(-1,
               "Marker %s not found; effect ignored\n",
               EffNme[i]);
          continue;
        }

      e = Effect.NrEff;

      Effect.Eff[e].Eff = 0;
      Effect.Eff[e].NrLoci = 0;
      Effect.Eff[e].Loc = NULL;

      Effect.Eff[e+1].Eff = 0;
      Effect.Eff[e+1].NrLoci = 0;
      Effect.Eff[e+1].Loc = NULL;

      Insert_ALL_Loc(&(Effect.Eff[e]),tmpMP,EffAll[i]);
      Effect.Eff[e].Eff = (float)(u[i] / 2.0);

      Insert_ALL_Loc(&(Effect.Eff[e+1]),tmpMP,CmpAll[i]);
      Effect.Eff[e+1].Eff = (float)(-u[i] / 2.0);

      /*
       * GS uses the default 0/1/2 dosage coding, whereas the simulation
       * representation above is equivalent to -1/0/1.  The coding shift
       * therefore adds u[i] to the general mean for every transferred
       * marker.
       */
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
      info(-2,"Memory allocation failed in st_set_sim_ef\n");
      goto fail;
    }

  if (INT_MAX == NrEffects)
    {
      info(-2,"Too many simulation effects\n");
      goto fail;
    }

  if ((size_t)(NrEffects + 1) >
      SIZE_MAX / sizeof(*Effects))
    {
      info(-2,"Memory allocation failed in st_set_sim_ef\n");
      goto fail;
    }

  tmpEffects =
    realloc(Effects,
            (size_t)(NrEffects + 1) * sizeof(*Effects));

  if (NULL == tmpEffects)
    {
      info(-2,"Memory allocation failed in st_set_sim_ef\n");
      goto fail;
    }

  Effects = tmpEffects;
  Effect.Pos = NrEffects + 1;
  Effects[NrEffects] = Effect;
  NrEffects++;

  /*
   * Existing GValue arrays were allocated for the previous number of
   * effects and are no longer valid after a successful effect import.
   */
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
  const st_sim_marker_index_type *ia =
    (const st_sim_marker_index_type *)a;
  const st_sim_marker_index_type *ib =
    (const st_sim_marker_index_type *)b;

  return strcmp(ia->name,ib->name);
}


static gs_UNSIGNED st_sim_marker_index_find(
    const st_sim_marker_index_type *index,
    gs_UNSIGNED n,
    const char *name)
{
  st_sim_marker_index_type key;
  const st_sim_marker_index_type *hit;

  key.name = name;
  key.row = 0;

  hit = bsearch(&key,
                index,
                (size_t)n,
                sizeof(*index),
                st_sim_marker_index_cmp);

  if (NULL == hit)
    return ULONG_MAX;

  return hit->row;
}


void st_set_sim_pp_provide_data  (  gs_UNSIGNED       *  NoInd,
                                    gs_UNSIGNED       *  NoMar,
                                    char             (** MarName)[gs_MAXNME],
                                    GS_ALLELE_TYPE   ** All1,
                                    GS_ALLELE_TYPE   ** All2,
                                    int               * retval,
                                    char             ** set_name );


void st_set_sim_pp              ( int*     retval,      // Return value
			          char**   pop_name,    // Name of the population
			          char**   set_name     // Name of the data set
			         )
{
  gs_UNSIGNED NoInd = 0;
  gs_UNSIGNED NoMar = 0;
  gs_UNSIGNED i,j;
  gs_UNSIGNED *marker_row = NULL;
  char (*MarName)[gs_MAXNME] = NULL;
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
      info(-2,
           "Population name must contain 2 to 255 characters "
           "from A-Z, a-z, and 0-9\n");
      return;
    }

  find_list_element(*pop_name,&checkPopListEl);

  if (NULL != checkPopListEl)
    {
      info(-2,
           "First remove population %s\n",
           *pop_name);
      return;
    }

  st_set_sim_pp_provide_data(
    &NoInd,
    &NoMar,
    &MarName,
    &All1,
    &All2,
    &Sretval,
    set_name );

  if (0 != Sretval)
    return;

  if ((2 != NoHoms) ||
      (NULL == PopDescr) ||
      (cMP < 1) ||
      (0 == NoInd) ||
      (0 == NoMar) ||
      (NULL == MarName) ||
      (NULL == All1) ||
      (NULL == All2))
    {
      info(-2,
           "Simulation map or GS marker data are not available\n");
      return;
    }

  if (NoInd > (gs_UNSIGNED)INT_MAX)
    {
      info(-2,"Too many individuals in GS data set\n");
      return;
    }

  if ((size_t)NoMar >
      SIZE_MAX / sizeof(*marker_index))
    {
      info(-2,"GS marker data are too large\n");
      return;
    }

  marker_index =
    malloc((size_t)NoMar * sizeof(*marker_index));

  if (NULL == marker_index)
    goto fail_alloc;

  for (i=0; i < NoMar; i++)
    {
      if (NULL == memchr(MarName[i],'\0',gs_MAXNME))
        {
          info(-2,"Invalid marker name in GS data set\n");
          goto fail;
        }

      marker_index[i].name = MarName[i];
      marker_index[i].row = i;
    }

  qsort(marker_index,
        (size_t)NoMar,
        sizeof(*marker_index),
        st_sim_marker_index_cmp);

  for (i=1; i < NoMar; i++)
    {
      if (0 == strcmp(marker_index[i-1].name,
                      marker_index[i].name))
        {
          info(-2,
               "Duplicate marker name in GS data set: %s\n",
               marker_index[i].name);
          goto fail;
        }
    }

  if ((size_t)cMP >
      SIZE_MAX / sizeof(*points))
    {
      info(-2,"Simulation linkage map is too large\n");
      goto fail;
    }

  points = calloc((size_t)cMP,sizeof(*points));
  marker_row = malloc((size_t)cMP * sizeof(*marker_row));
  PB = malloc((size_t)cMP * sizeof(*PB));

  if ((NULL == points) ||
      (NULL == marker_row) ||
      (NULL == PB))
    goto fail_alloc;

  for (chrom=0; chrom < NoChroms; chrom++)
    {
      for (loc=0; loc < PopDescr[chrom].NrMappoints; loc++)
        {
          TMappoint *point =
            &(PopDescr[chrom].Mappoint[loc]);

          if ((point->counter < 0) ||
              (point->counter >= cMP) ||
              (NULL != points[point->counter]))
            {
              info(-2,"Invalid installed simulation linkage map\n");
              goto fail;
            }

          points[point->counter] = point;
        }
    }

  for (i=0; i < (gs_UNSIGNED)cMP; i++)
    {
      if (NULL == points[i])
        {
          info(-2,"Invalid installed simulation linkage map\n");
          goto fail;
        }

      marker_row[i] =
        st_sim_marker_index_find(marker_index,
                                 NoMar,
                                 points[i]->name);

      if (ULONG_MAX == marker_row[i])
        {
          info(-2,
               "Marker %s not found in GS data set\n",
               points[i]->name);
          goto fail;
        }

      if ((0 == i) ||
          (points[i]->chrom != points[i-1]->chrom))
        {
          PB[i] = 0.0;
        }
      else
        {
          PB[i] =
            points[i-1]->Pos +
            (points[i]->Pos - points[i-1]->Pos) / 2.0;
        }
    }

  if ((size_t)cMP >
      SIZE_MAX / (size_t)NoInd)
    {
      info(-2,"GS population data are too large\n");
      goto fail;
    }

  loci_count =
    (size_t)cMP * (size_t)NoInd;

  if (loci_count > SIZE_MAX / 2)
    {
      info(-2,"GS population data are too large\n");
      goto fail;
    }

  Cnt = 2 * loci_count;

  if ((Cnt > (size_t)LONG_MAX) ||
      (Cnt > SIZE_MAX / sizeof(double)))
    {
      info(-2,"GS population data are too large\n");
      goto fail;
    }

  _Ind = malloc(Cnt * sizeof(int));
  _Chrom = malloc(Cnt * sizeof(int));
  _Hom = malloc(Cnt * sizeof(int));
  _Pos = malloc(Cnt * sizeof(double));
  _All = malloc(Cnt * sizeof(int));

  if ((NULL == _Ind) ||
      (NULL == _Chrom) ||
      (NULL == _Hom) ||
      (NULL == _Pos) ||
      (NULL == _All))
    goto fail_alloc;

  row = 0;

  for (j=0; j < NoInd; j++)
    {
      for (i=0; i < (gs_UNSIGNED)cMP; i++)
        {
          unsigned short a1;
          unsigned short a2;

          idx =
            (size_t)marker_row[i] * (size_t)NoInd +
            (size_t)j;

          if (!st_sim_allele(All1[idx],&a1) ||
              !st_sim_allele(All2[idx],&a2))
            {
              info(-2,
                   "Allele code cannot be represented in simulation\n");
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

  snprintf(d,
           sizeof(d),
           "%ld",
           (long)Cnt);

  init_population(pop_name,
                  cCount,
                  _Ind,
                  _Chrom,
                  _Hom,
                  _Pos,
                  _All);

  find_list_element(*pop_name,&checkPopListEl);

  if (NULL == checkPopListEl)
    {
      info(-2,"Population import from GS data failed\n");
      goto fail;
    }

  population_created = 1;
  Pop = checkPopListEl->Pop;

  if (Pop->NoInds != (POPULATION)NoInd)
    {
      info(-2,"Population import from GS data failed\n");
      goto fail_import;
    }

  /*
   * Verify the imported population against the source marker matrix at
   * every installed simulation marker.
   */
  for (j=0; j < NoInd; j++)
    {
      TInd *Ind =
        get_individual(Pop,(POPULATION)j);

      for (i=0; i < (gs_UNSIGNED)cMP; i++)
        {
          unsigned short expected1;
          unsigned short expected2;

          chrom = points[i]->chrom - 1;

          idx =
            (size_t)marker_row[i] * (size_t)NoInd +
            (size_t)j;

          if (!st_sim_allele(All1[idx],&expected1) ||
              !st_sim_allele(All2[idx],&expected2))
            goto fail_import;

          if ((get_allele(&(Ind->Chrom[chrom].Hom[0]),
                          points[i]->Pos) != expected1) ||
              (get_allele(&(Ind->Chrom[chrom].Hom[1]),
                          points[i]->Pos) != expected2))
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

  info(-2,
       "Memory allocation failed in st_set_sim_pp\n");
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


/**************************************************************************
End:   Interface to simulation routines / Part2 set_sim            Jun 2019
**************************************************************************/

#include "SelectionTools02a.h"
  

