//////////////////////////////////////////////////////////////////////////
// SelectionTools01.c
///////////////////////////////////////////////////////////////////////////

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

#include "SelectionToolsAllele.h"

/* Include the OpenMP header before R headers.  Rinternals.h defines
   match as Rf_match, while newer omp.h versions use match(...) in
   OpenMP directives. */
#ifdef _OPENMP
  #ifndef ST
    #define MT
    #include <omp.h>
  #endif
#endif

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

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
  #define omp_set_num_threads(thr) {}
  #define omp_get_num_threads() 0
  #define omp_set_dynamic(thr) {}
#endif

#ifdef MultPrec
  #include <gmp.h>
#endif

/* Maximum string size for general use */
#define gs_MAXSTR 4096

/* Maximum string size for names of markers and indviduals */
#define gs_MAXNME 256

/* Maximum number of stored allele states per marker, including missing */
#define gs_MAXALL (GS_ALLELE_MAX + 1)

/* Maximum line size in input files */
#define gs_MAXLNE 65536


/* Common errors */
#define ERR_F {gs_info(gv,-2,"File open error");*retval=-2;return;}
#define ERR_M {gs_info(gv,-2,"Memory allocation error");*retval=-2;return;}

/* Global variables and definitions  */
typedef unsigned long int gs_UNSIGNED ;
# define gs_UNSIGNED_MAX ULONG_MAX

// typedef __float128 gs_FLOAT ;
// typedef __float80 gs_FLOAT ;
typedef double gs_FLOAT ;

/* Type linkage map*/
typedef struct gs_mappoint_type
{
  int          chrom;
  double       pos;
  const char  *name;
  const char  *class;
  gs_UNSIGNED  nr;          // Number of the data line
} gs_mappoint_type;

/* Pair of two unsigned values */
typedef struct
{
  gs_UNSIGNED a;
  gs_UNSIGNED n;
} gs_usp;

/* Parameters of a cross of two individuals*/
typedef struct gs_cross_type
{
  gs_UNSIGNED p1;
  gs_UNSIGNED p2;
  gs_FLOAT    gd;  // Genetic distance of the parents
  gs_FLOAT    mu;  // Expected value of the cross
  gs_FLOAT    mi;  // Minimum value of the cross
  gs_FLOAT    ma;  // Maximum value of the cross
  gs_FLOAT    va;  // Variance of the cross
  gs_FLOAT    es;  // Expectation of the selected fraction
} gs_cross_type;
 
/* Chromosomes and recombination frequencies */
typedef struct gs_chrom_type
{
  gs_UNSIGNED        NLocC;       // Number of loci on the chromosome
  gs_mappoint_type** first_loc;   // Pointer to the first locus 
  gs_UNSIGNED        first_loc_i; // Index of the first locus 
  gs_FLOAT*          rf;          // Matrix of recombination frequencies
  gs_FLOAT*          q;           // Matrix of conditional probabilites
  gs_FLOAT*          ld;          // Matrix of LD measures
} gs_chrom_type;


/* Data set */

typedef struct gs_varset_type {
  
  char name[gs_MAXSTR];
  char msg [3*gs_MAXSTR];
  int info_level;

  /* Variables initialized by gs_read_markerdata */

  int MaxIL;        /* Maximum length of an individual name */
  int MaxML;        /* Maximum length of an marker name */

  gs_UNSIGNED  
    NoInd ,   /* Number of individuals                   */
    NoMar ;   /* Number of markers                       */

  char (*IndName) [gs_MAXNME]; /* Individual names*/ 
  char (*MarName) [gs_MAXNME]; /* Marker names*/ 

  GS_ALLELE_TYPE * All1  ; /* Genotypic data */
  GS_ALLELE_TYPE * All2  ; /* Genotypic data */
  gs_FLOAT  * genv  ; /* Genotypic value */
  gs_FLOAT  * genv1  ; /* Genotypic value first homologue chrom*/
  gs_FLOAT  * genv2  ; /* Genotypic value second homologue chrom*/

  /* Variables initiaized by gs_marker_stats*/

  gs_UNSIGNED (*NoAll)  ;             /* Number of allele per Marker */
  GS_ALLELE_TYPE *All_S;              /* Allele names, storage */
  GS_ALLELE_TYPE **All;               /* Allele names, marker access */
  gs_UNSIGNED *AllCnt_S;              /* Allele count, storage */
  gs_UNSIGNED **AllCnt;               /* Allele count, marker access */

  double (*Pic) ; /* Pic value */
  double (*Mis) ; /* Frequency of Missing values per marker */

  double (*IndMis) ; /* Frequency of missing markers per ind */

  /* Variables for gs_build_Z and gs_build_X */

  double* Z ;   /* Design matrix */
  double* V ;   /* ZZ'           */
  double* y ;   /* Phenotype     */

  gs_UNSIGNED NoEff ;          /* Number of effects  */
  char (*EffNme) [gs_MAXNME] ; /* Marker names */ 
  int *EffAll ;                /* Allele number of the effect */
  int *CmpAll ;                /* Allele number complementary allele */

  double aa; double aA; double AA; /* Coding alleles */
  double an; double nn; double nA; /* Coding alleles */

  /* Variables for gs_build_A */
  /* to solve Ax=b */

  gs_UNSIGNED dimu ;
  gs_UNSIGNED dimA ;
  double*  A ;   
  double*  u ;   
  double*  b ;
  double*  As;      // Save the value of A and b when using 
  double*  bs ;     // in place matrix operations

  /* Variables for the transformed model */

  double*  GI ;  // G-star inverted: inv(ZZ')  
  double*  a ;   // Solution vector for the transformed model   

  /* p-value for test of u */
  double*  p;

  /* Shrinkage factor */
  double* ladi ;

  /* Variables for gs_validation */
  double* est ;

  /* Variables for gs_single_marker_reg_01 and _aov_01*/
  double* SMsqm ;
  double* SMsqe ;
  double* SMF   ;
  double* SMp   ;
  double* SMvar ;

  /* Variables for variance component estimation */
  double* sigma;
  double  sigma_e;

  /* Linkage map */  
  gs_UNSIGNED         NLoc;       // Number of loci in the map
  gs_mappoint_type  * map_S;      // Array for storage
  char              * map_str_S;  // Storage for map names and classes
  size_t              map_str_N;  // Number of bytes in map_str_S
  gs_mappoint_type ** map;        // Pointers for access
  gs_UNSIGNED       * mdl;        // Index of the data line for a mappoint

  /* Storage of recoded of haploblock alleles */
  GS_ALLELE_TYPE * Ref1  ; 
  GS_ALLELE_TYPE * Ref2  ; 

  /* Linkage map of haplotypes */
  gs_UNSIGNED         NHap  ;      // Number of haplotype blocks
  gs_mappoint_type  * Hap_S ;      // Array for storage
  char              * Hap_str_S;   // Storage for haplotype names/classes
  size_t              Hap_str_N;   // Number of bytes in Hap_str_S
  gs_UNSIGNED  * h_u;              // Lower border of block 
  gs_UNSIGNED  * h_o;          // Upper border
  gs_UNSIGNED  * h_n;          // No. of loci in the block

  /* Parameters of a cross */
  gs_UNSIGNED      NCrs;   // Number of crosses
  gs_cross_type  * crs_S;  // Array for storage
  gs_cross_type  **crs;    // Array of pointers for access

  /* Chromosomes and recombination frequencies, ld values */
  gs_UNSIGNED NChrom;      // Number of chromosomes
  gs_chrom_type * chrom ;  // Array of chromosome parameters

  /* Temprary storage for iterations without free */
  double *t0,*t1,*t2,*t3,*t4,*t5,*t6,*t7,*t8,*t9;
  gs_UNSIGNED *u0;   
} gs_varset_type;


/* Standard variable set is accessed with pointer */
gs_varset_type GV =  {
  .name       = "default",    
  .msg        = "",
  .info_level = 0 ,
  .MaxIL      = 0, 
  .MaxML      = 0, 
  .NoInd      = 0,  
  .NoMar      = 0,  
  .IndName    = NULL,
  .MarName    = NULL,
  .All1       = NULL,
  .All2       = NULL,
  .genv       = NULL,
  .genv1      = NULL,
  .genv2      = NULL,
  .NoAll      = NULL,
  .All_S      = NULL,
  .All        = NULL,
  .AllCnt_S   = NULL,
  .AllCnt     = NULL,
  .Pic        = NULL,
  .Mis        = NULL,
  .IndMis     = NULL,
  .Z          = NULL, 
  .V          = NULL, 
  .y          = NULL, 
  .NoEff      = 0,   
  .EffNme     = NULL,
  .EffAll     = NULL,
  .CmpAll     = NULL,
  .aa         = 0, 
  .aA         = 1,
  .AA         = 2,
  .an         = 0.5,
  .nn         = 1  ,
  .nA         = 1.5,
  .dimu       = 0,
  .dimA       = 0,
  .A          = NULL,   
  .u          = NULL,   
  .b          = NULL,
  .As         = NULL,
  .bs         = NULL,
  .GI         = NULL,
  .a          = NULL,
  .p          = NULL,
  .ladi       = NULL,
  .est        = NULL,
  .SMsqm      = NULL,
  .SMsqe      = NULL,
  .SMF        = NULL,
  .SMp        = NULL,
  .SMvar      = NULL,
  .sigma      = NULL,
  .sigma_e    = 0,
  .NLoc       = 0,
  .map_S      = NULL,
  .map_str_S  = NULL,
  .map_str_N  = 0,
  .map        = NULL,
  .mdl        = NULL,
  .Ref1       = NULL,
  .Ref2       = NULL,
  .NHap       = 0,
  .Hap_S      = NULL,
  .Hap_str_S  = NULL,
  .Hap_str_N  = 0,
  .h_u        = NULL,
  .h_o        = NULL,
  .h_n        = NULL,
  .NCrs       = 0,
  .crs_S      = NULL,
  .crs        = NULL,
  .NChrom     = 0,
  .chrom      = NULL,
  .t0=NULL, .t1=NULL, .t2=NULL, .t3=NULL, .t4=NULL, 
  .t5=NULL, .t6=NULL, .t7=NULL, .t8=NULL, .t9=NULL,
  .u0=NULL
} ;

gs_varset_type *gv = &GV; 
gs_varset_type *sh = &GV; 

gs_varset_type ES; gs_varset_type *es = &ES; int esi = 0;
gs_varset_type VS; gs_varset_type *vs = &VS; int vsi = 0;

gs_varset_type **gs_data = NULL;
gs_UNSIGNED      gs_ndta = 0;

/* Global constants */
const char * gs_tok_del = "; \r\n\t\""  ;
char * nix  = ""; char ** nixp  = &nix;
int    nein = 0 ; int  *  neinp = & nein;
double null = 0 ; double *  nullp = & null;

int Srv=0; int* Sretval = &Srv;
char*empty=""; char ** Sout_filename = &empty;
int auf = 0; int*Sauxfiles = &auf;
int mz = -2; int*imz = &mz;

///////////////////////////////////////////////////////////////////////////
// Utility functions
///////////////////////////////////////////////////////////////////////////

double**  st_dmatv (double *a, gs_UNSIGNED nr, gs_UNSIGNED nc)
{
  double **p = malloc ( nr * sizeof(double*) );

  p[0] = &a[0];
  for( gs_UNSIGNED i=1; i<nr; i++)
    p[i] = p[i-1] + nc;

  return (p);
}

///////////////////////////////////////////////////////////////////////////
// Random numbers                                                        //
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// Distribution functions                                                //
///////////////////////////////////////////////////////////////////////////

void gs_info ( gs_varset_type *, int , char* );


double st_fprob(double f, gs_UNSIGNED df1, gs_UNSIGNED df2)
{
    /* pf(x, df1, df2, lower_tail, log_p) */
    return pf( f,
	       (double)df1,
	       (double)df2,
	       /*lower_tail=*/ 0,
	       /*log_p=*/      0);
}


double st_tprob(double t, gs_UNSIGNED df)
{
    return 2.0 * pt( fabs(t),
		     (double)df,
		     /*lower_tail=*/ 0,
		     /*log_p=*/      0);
}


void st_fprob_GV ( double* f, int* df1, int* df2, double* prob)
{
  *prob = st_fprob( *f,	(gs_UNSIGNED)(*df1), (gs_UNSIGNED)(*df2) ); 
}


void st_tprob_GV ( double* t, int* df, double* prob)
{
  *prob = st_tprob( *t,	(gs_UNSIGNED)(*df) ); 
}

///////////////////////////////////////////////////////////////////////////
// Distribution functions                                                //
///////////////////////////////////////////////////////////////////////////

//##############################################################################
//# Begin:  Linear algebra                                               2026-02
//##############################################################################

  #include <R.h>
  #include <R_ext/RS.h>
  #include <R_ext/BLAS.h>
  #include <R_ext/Lapack.h>

  #include "MathLib11.h"

//##############################################################################
//# End:  Linear algebra                                                 2026-02
//##############################################################################




///////////////////////////////////////////////////////////////////////////
// Linear algebra without libraries                                      //
///////////////////////////////////////////////////////////////////////////

/* matrix matrix add */
void gs_mm_add ( gs_UNSIGNED nr, gs_UNSIGNED nc,
		 double *A, double* B, double *C )
{
  gs_UNSIGNED i, j;

  for (i = 0; i < nr; i++) 
    for (j = 0; j < nc; j++) 
	C [ i*nc + j ] = A [ i*nc + j ] + B [ i*nc + j ] ;
}

/* scalar matrix add */
void gs_sm_add ( gs_UNSIGNED nr, gs_UNSIGNED nc,
		 double *A, double* B, double *C )
{
  gs_UNSIGNED i, j;

  for (i = 0; i < nr; i++) 
    for (j = 0; j < nc; j++) 
	C [ i*nc + j ] = *A  + B [ i*nc + j ] ;
}

///////////////////////////////////////////////////////////////////////////
// // Linear algebra without libraries                                      //
///////////////////////////////////////////////////////////////////////////

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void gs_set_num_threads (int * activate)
{

  if ( *activate > MAXTHREADS)  *activate = MAXTHREADS ;
  if ( *activate < 1 )          *activate = 1 ;
  omp_set_num_threads (*activate);


#if defined(use_openblas) && defined(MT)
  openblas_set_num_threads(*activate); 
  goto_set_num_threads(*activate); 
#endif

#if defined(use_mkl) && defined(MT)
  mkl_set_num_threads(*activate); 
#endif

  
}

void gs_get_num_threads (int * activate)
{
#ifdef MT
  #pragma omp parallel
  { 
    *activate = omp_get_num_threads ();
  }
#else
  { 
    *activate = omp_get_num_threads ();
  }
#endif  
  snprintf(gv->msg, sizeof(gv->msg),
         "Number of threads %i", *activate);
  gs_info(gv,0,gv->msg);
}


void gs_reset(gs_varset_type *gv)
{
   gs_UNSIGNED i;

   gv->msg[0]     = '\0';
   gv->MaxIL      = 0;
   gv->MaxML      = 0;
   gv->NoInd      = 0;
   gv->NoMar      = 0;

   GS_FREE(gv->IndName);
   GS_FREE(gv->MarName);
   GS_FREE(gv->All1);
   GS_FREE(gv->All2);
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
   gv->NLoc    = 0;

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

void gs_init( gs_varset_type *gv )
{
  gv->name[0]    = '\0';
  gv->msg[0]     = '\0';
  gv->info_level = 0;
  gv->MaxIL      = 0; 
  gv->MaxML      = 0; 
  gv->NoInd      = 0;  
  gv->NoMar      = 0;  
  gv->IndName    = NULL;
  gv->MarName    = NULL;
  gv->All1       = NULL;
  gv->All2       = NULL;
  gv->genv       = NULL;
  gv->genv1      = NULL;
  gv->genv2      = NULL;
  gv->NoAll      = NULL;
  gv->All_S      = NULL;
  gv->All        = NULL;
  gv->AllCnt_S   = NULL;
  gv->AllCnt     = NULL;
  gv->Pic        = NULL;
  gv->Mis        = NULL;
  gv->IndMis     = NULL;
  gv->Z          = NULL; 
  gv->V          = NULL; 
  gv->y          = NULL; 
  gv->NoEff      = 0;   			                       
  gv->EffNme     = NULL;
  gv->EffAll     = NULL;
  gv->CmpAll     = NULL;
  gv->aa         = 0; 
  gv->aA         = 1;
  gv->AA         = 2;
  gv->an         = 0.5; 
  gv->nn         = 1;
  gv->nA         = 1.5;
  gv->dimu       = 0;                   
  gv->dimA       = 0; 
  gv->A          = NULL;   
  gv->u          = NULL;   
  gv->b          = NULL;
  gv->As         = NULL;   
  gv->bs         = NULL;
  gv->GI         = NULL;   
  gv->a          = NULL;   
  gv->p          = NULL;   
  gv->ladi       = NULL;
  gv->est        = NULL;
  gv->SMsqm      = NULL;
  gv->SMsqe      = NULL;
  gv->SMF        = NULL;
  gv->SMp        = NULL;
  gv->SMvar      = NULL;
  gv->sigma      = NULL;
  gv->sigma_e    = 0;
  gv->NLoc       = 0;
  gv->map_S      = NULL;
  gv->map_str_S  = NULL;
  gv->map_str_N  = 0;
  gv->map        = NULL;
  gv->mdl        = NULL;
  gv->Ref1       = NULL;
  gv->Ref2       = NULL;
  gv->NHap       = 0;
  gv->Hap_S      = NULL;
  gv->Hap_str_S  = NULL;
  gv->Hap_str_N  = 0;
  gv->h_u        = NULL;
  gv->h_o        = NULL;
  gv->h_n        = NULL;
  gv->NCrs       = 0;
  gv->crs_S      = NULL;
  gv->crs        = NULL;
  gv->NChrom     = 0;
  gv->chrom      = NULL;
  gv->t0=NULL;  gv->t1=NULL;  gv->t2=NULL;  gv->t3=NULL;  gv->t4=NULL; 
  gv->t5=NULL;  gv->t6=NULL;  gv->t7=NULL;  gv->t8=NULL;  gv->t9=NULL;
  gv->u0=NULL;
}

/*--------------*/

static int st_console_output_allowed(void)
{
#ifdef MT
  if (omp_in_parallel())
    return 0;
#endif

  return 1;
}

void gs_info ( gs_varset_type *gv,
               int level,
               char* message )
{
    if ( gv->info_level >= level )
    {
        if ( !st_console_output_allowed() )
            return;

        if ( 1 <= level )
            Rprintf(
              "I (data set '%s'): %s\n",
              gv->name,
              message );

        else if ( 0 == level )
            Rprintf(
              "M (data set '%s'): %s\n",
              gv->name,
              message );

        else if ( -1 == level )
            Rprintf(
              "W (data set '%s'): %s\n",
              gv->name,
              message );

        else if ( -2 == level )
            Rprintf(
              "E (data set '%s'): %s\n",
              gv->name,
              message );
    }
}

void gs_info_r ( gs_varset_type *gv,
                 int level,
                 char* message )
{
    if ( gv->info_level >= level )
    {
        if ( !st_console_output_allowed() )
            return;

        if ( 1 <= level )
            Rprintf(
              "I (data set '%s'): %s\r",
              gv->name,
              message );

        else if ( 0 == level )
            Rprintf(
              "M (data set '%s'): %s\r",
              gv->name,
              message );

        else if ( -1 == level )
            Rprintf(
              "W (data set '%s'): %s\r",
              gv->name,
              message );

        else if ( -2 == level )
            Rprintf(
              "E (data set '%s'): %s\r",
              gv->name,
              message );
    }
}

void gs_r_info(int* level, char** message)
{
  gs_info(gv, *level, *message );
}

void gs_set_info_level(gs_varset_type *gv , int *level )
{
  if ( *level < -2 )  *level = -2;
  if ( *level > 3 )   *level =  3;
  gv->info_level = *level;
  snprintf(gv->msg, sizeof(gv->msg),
         "Info level %i",
         gv->info_level);
  gs_info(gv,3,gv->msg);
}

void gs_get_info_level ( gs_varset_type *gv , int *level )
{
  *level = gv->info_level;
}

gs_varset_type* gs_fdta (char* set_name)
{
  gs_UNSIGNED i;
  gs_varset_type* sp = NULL;

  /* Default data set */
  if ( 0 == strcmp ( set_name , "default" ) ) 
    sp = gv;

  if (NULL == sp)
    {
      /* User defined data set is available */
      for ( i = 0; i < gs_ndta; i ++ )
	if ( 0 == strcmp ( set_name , gs_data[i]->name ) )
	  sp = gs_data[i];
    }

  if (NULL == sp)
    {
      /* New dataset */  
      gs_ndta++;
      gs_data = realloc(gs_data,gs_ndta*sizeof(gs_varset_type*));
      gs_data[gs_ndta-1] = malloc(sizeof(gs_varset_type));
      sp = gs_data[gs_ndta-1];
      
      gs_init(sp);
      strcpy (sp->name,set_name);
      sp->info_level = gv->info_level;

      snprintf(sp->msg, sizeof(sp->msg),
         "New data set generated '%s'",
         sp->name);

      gs_info(sp,0,sp->msg);
    }

  snprintf(sp->msg, sizeof(sp->msg),
	   "Working on data set '%s'",
	   sp->name);

  gs_info(sp,2,sp->msg);
  
  return sp;

} 


void gs_get_info_level_GV ( int *level,
			    char** set_name )
{
  gs_get_info_level( gs_fdta(*set_name) , level );
}

void gs_set_info_level_GV ( int *level,
			    char** set_name )
{
  gs_set_info_level( gs_fdta(*set_name), level );
}

void gs_set_all_info_levels ( int *level )
{
  gs_set_info_level( gv, level );
  for ( gs_UNSIGNED i = 0 ; i < gs_ndta; i++)
    gs_set_info_level(gs_data[i], level );
}


/* qsort pointers to string comparison function*/
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

static void *gs_realloc_array(void *ptr,
                              size_t count,
                              size_t element_size)
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

static int gs_map_pool_valid(const gs_mappoint_type *map_S,
                             gs_UNSIGNED nloc,
                             const char *pool,
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

      if ((offset >= pool_n) ||
          (map_S[k].name != pool + offset))
        return 0;

      end = memchr(pool + offset, '\0', pool_n - offset);
      if (NULL == end)
        return 0;
      len = (size_t)(end - (pool + offset));
      if (len >= gs_MAXSTR)
        return 0;
      offset += len + 1;

      if ((offset >= pool_n) ||
          (map_S[k].class != pool + offset))
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

static int gs_pool_put_string(char *pool,
                              size_t pool_n,
                              size_t *offset,
                              const char *text,
                              const char **stored)
{
  size_t len;
  size_t needed;

  if ((NULL == pool) || (NULL == offset) ||
      (NULL == text) || (NULL == stored))
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
                                     gs_UNSIGNED nr,
                                     char *pool,
                                     size_t pool_n,
                                     size_t *offset)
{
  dst->chrom = src->chrom;
  dst->pos = src->pos;
  dst->nr = nr;

  if ((0 != gs_pool_put_string(pool, pool_n, offset,
                               src->name, &(dst->name))) ||
      (0 != gs_pool_put_string(pool, pool_n, offset,
                               src->class, &(dst->class))))
    return -1;

  return 0;
}

typedef struct
{
  int          chrom;
  double       pos;
  size_t       name_offset;
  size_t       class_offset;
  gs_UNSIGNED  nr;
} gs_mappoint_build_type;

static int gs_string_builder_append(char **pool,
                                    size_t *used,
                                    size_t *capacity,
                                    const char *text,
                                    size_t *stored_offset)
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

static int gs_make_name_index(char (*names)[gs_MAXNME],
                              gs_UNSIGNED n,
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

static int gs_name_index_has_duplicates(const gs_name_index_type *index,
                                        gs_UNSIGNED n)
{
  for (gs_UNSIGNED i = 1; i < n; i++)
    if (0 == strcmp(index[i-1].name, index[i].name))
      return 1;

  return 0;
}

static gs_UNSIGNED gs_name_index_find(const gs_name_index_type *index,
                                      gs_UNSIGNED n,
                                      const char *name)
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
                                             gs_UNSIGNED n,
                                             const char *name)
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

static void gs_tmp_commit(gs_varset_type *dst, gs_varset_type *tmp)
{
  gs_reset(dst);
  *dst = *tmp;
  gs_init(tmp);
}

static int gs_alloc_genotype_matrix(gs_varset_type *dst,
                                    gs_UNSIGNED no_mar,
                                    gs_UNSIGNED no_ind)
{
  size_t n, bytes;

  if ((0 == no_mar) || (0 == no_ind))
    return -1;

  if ((0 != gs_size_mul((size_t)no_mar, (size_t)no_ind, &n)) ||
      (0 != gs_size_mul(n, sizeof(GS_ALLELE_TYPE), &bytes)))
    return -1;

  dst->All1 = malloc(bytes);
  dst->All2 = malloc(bytes);

  if ((NULL == dst->All1) || (NULL == dst->All2))
    return -1;

  return 0;
}

static void gs_set_genotypes_missing(gs_varset_type *dst)
{
  size_t n = (size_t)dst->NoMar * (size_t)dst->NoInd;

  for (size_t i = 0; i < n; i++)
    {
      dst->All1[i] = GS_ALLELE_MISSING;
      dst->All2[i] = GS_ALLELE_MISSING;
    }
}

static int gs_copy_map_primary(gs_varset_type *dst,
                               const gs_varset_type *src)
{
  size_t bytes;
  size_t offset = 0;
  gs_mappoint_type *map_S = NULL;
  char *map_str_S = NULL;
  gs_mappoint_type **map = NULL;
  gs_UNSIGNED *mdl = NULL;
  gs_UNSIGNED *storage = NULL;
  int result = -2;

  if ((0 == src->NLoc) &&
      (NULL == src->map_S) &&
      (NULL == src->map_str_S) &&
      (0 == src->map_str_N) &&
      (NULL == src->map) &&
      (NULL == src->mdl))
    return 0;

  if ((0 == src->NLoc) ||
      (0 == src->NoMar) ||
      (src->NLoc > src->NoMar) ||
      (NULL == src->map_S) ||
      (NULL == src->map_str_S) ||
      (0 == src->map_str_N) ||
      (NULL == src->map) ||
      (NULL == src->mdl) ||
      !gs_map_pool_valid(src->map_S,
                         src->NLoc,
                         src->map_str_S,
                         src->map_str_N))
    return -2;

  storage = gs_malloc_array((size_t)src->NoMar,sizeof(*storage));
  if (NULL == storage)
    return -1;

  for (gs_UNSIGNED m = 0; m < src->NoMar; m++)
    storage[m] = gs_UNSIGNED_MAX;

  for (gs_UNSIGNED k = 0; k < src->NLoc; k++)
    {
      gs_UNSIGNED m = src->map_S[k].nr;

      if ((m >= src->NoMar) ||
          (gs_UNSIGNED_MAX != storage[m]))
        goto cleanup;

      storage[m] = k;
    }

  for (gs_UNSIGNED lne = 0; lne < src->NLoc; lne++)
    {
      gs_UNSIGNED m = src->mdl[lne];

      if ((NULL == src->map[lne]) ||
          (m >= src->NoMar) ||
          (gs_UNSIGNED_MAX == storage[m]) ||
          (src->map[lne] != &(src->map_S[storage[m]])) ||
          (src->map[lne]->nr != m))
        goto cleanup;
    }

  if (0 != gs_size_mul((size_t)src->NLoc,
                       sizeof(gs_mappoint_type),
                       &bytes))
    {
      result = -1;
      goto cleanup;
    }
  map_S = malloc(bytes);

  map_str_S = malloc(src->map_str_N);

  if (0 != gs_size_mul((size_t)src->NLoc,
                       sizeof(gs_mappoint_type *),
                       &bytes))
    {
      result = -1;
      goto cleanup;
    }
  map = malloc(bytes);

  if (0 != gs_size_mul((size_t)src->NLoc,
                       sizeof(gs_UNSIGNED),
                       &bytes))
    {
      result = -1;
      goto cleanup;
    }
  mdl = malloc(bytes);

  if ((NULL == map_S) ||
      (NULL == map_str_S) ||
      (NULL == map) ||
      (NULL == mdl))
    {
      result = -1;
      goto cleanup;
    }

  memcpy(map_str_S, src->map_str_S, src->map_str_N);
  memcpy(mdl,
         src->mdl,
         (size_t)src->NLoc * sizeof(gs_UNSIGNED));

  for (gs_UNSIGNED k = 0; k < src->NLoc; k++)
    if (0 != gs_copy_map_point_to_pool(&(map_S[k]),
                                       &(src->map_S[k]),
                                       src->map_S[k].nr,
                                       map_str_S,
                                       src->map_str_N,
                                       &offset))
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
    if (NULL != comma) *comma = '.';
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

  while ((0 < len) && isspace((unsigned char)text[len-1]))
    len--;

  return ((1 == len && ('-' == text[0] || '.' == text[0])) ||
          (2 == len && '-' == text[0] && '-' == text[1]) ||
          (2 == len && '-' == text[0] && '1' == text[1]) ||
          (2 == len && 'N' == text[0] && 'A' == text[1]));
}

static int gs_parse_allele_span(const char *text,
                                size_t len,
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

  while ((0 < len) && isspace((unsigned char)text[len-1]))
    len--;

  if (0 == len)
    return -1;

  if (gs_allele_missing_span(text,len))
    {
      *allele = GS_ALLELE_MISSING;
      return 0;
    }

  if (1 == len)
    {
      switch (text[0])
        {
        case 'A': case 'a': *allele = (GS_ALLELE_TYPE)1; return 0;
        case 'C': case 'c': *allele = (GS_ALLELE_TYPE)2; return 0;
        case 'G': case 'g': *allele = (GS_ALLELE_TYPE)3; return 0;
        case 'T': case 't': *allele = (GS_ALLELE_TYPE)4; return 0;
        case 'D': case 'd': *allele = (GS_ALLELE_TYPE)5; return 0;
        default: break;
        }
    }

  if (len >= sizeof(buf))
    return -1;

  memcpy(buf,text,len);
  buf[len] = '\0';

  errno = 0;
  v = strtol(buf,&end,10);
  if ((ERANGE == errno) || (end == buf) || ('\0' != *end) ||
      (v < 1) || (v > GS_ALLELE_MAX))
    return -1;

  *allele = (GS_ALLELE_TYPE)v;
  return 0;
}

static int gs_parse_single_allele(const char *text,
                                  GS_ALLELE_TYPE *allele)
{
  if (NULL == text)
    return -1;

  return gs_parse_allele_span(text,strlen(text),allele);
}

static int gs_allele_cmp(const void *a, const void *b)
{
  GS_ALLELE_TYPE aa = *(const GS_ALLELE_TYPE *)a;
  GS_ALLELE_TYPE bb = *(const GS_ALLELE_TYPE *)b;

  return (aa > bb) - (aa < bb);
}


static int gs_aggregate_genotype(GS_ALLELE_TYPE *a1,
                                 GS_ALLELE_TYPE *a2,
                                 unsigned char *state,
                                 GS_ALLELE_TYPE a,
                                 GS_ALLELE_TYPE b,
                                 int pair)
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

      if (((*a1 == a) || (*a1 == b)) &&
          ((*a2 == a) || (*a2 == b)))
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

  if ((NULL == gv) || (NULL == gv->All1) || (NULL == gv->All2) ||
      (0 == gv->NoInd) || (0 == gv->NoMar) ||
      (0 != gs_size_mul((size_t)gv->NoInd,(size_t)gv->NoMar,&n)))
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

static size_t gs_count_sep_fields(const char *line, char sep)
{
  size_t n = 1;

  for (const char *p = line; '\0' != *p; p++)
    if (sep == *p)
      n++;

  return n;
}

static size_t gs_split_sep_fields(char *line,
                                  char sep,
                                  char **fields,
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

static gs_UNSIGNED *gs_map_storage_index(
    gs_varset_type *gv,
    int             require_complete,
    int             require_unique_names)
{
  gs_UNSIGNED *storage = NULL;
  unsigned char *seen = NULL;
  char **names = NULL;
  gs_UNSIGNED k,lne,m;
  int ok = 0;

  if ((0 == gv->NoMar) ||
      (NULL == gv->MarName))
    {
      gs_info(gv,-2,"No marker data available");
      return NULL;
    }

  if ((0 == gv->NLoc) ||
      (NULL == gv->map_S) ||
      (NULL == gv->map_str_S) ||
      (0 == gv->map_str_N) ||
      (NULL == gv->map) ||
      (NULL == gv->mdl))
    {
      gs_info(gv,-2,"No linkage map defined");
      return NULL;
    }

  if (gv->NLoc > gv->NoMar)
    {
      gs_info(gv,-2,
              "More linkage-map loci than marker-data rows");
      return NULL;
    }

  if (require_complete &&
      (gv->NLoc != gv->NoMar))
    {
      gs_info(gv,-2,
              "Linkage map is not matching the marker data");
      return NULL;
    }

  if (gv->NoMar >
      (gs_UNSIGNED)(((size_t)-1) / sizeof(*storage)))
    {
      gs_info(gv,-2,"Memory allocation error");
      return NULL;
    }

  storage = malloc((size_t)gv->NoMar * sizeof(*storage));
  seen = calloc((size_t)gv->NoMar, sizeof(*seen));

  if ((NULL == storage) ||
      (NULL == seen))
    {
      gs_info(gv,-2,"Memory allocation error");
      goto cleanup;
    }

  for (m=0; m < gv->NoMar; m++)
    storage[m] = gs_UNSIGNED_MAX;

  if (!gs_map_pool_valid(gv->map_S,
                         gv->NLoc,
                         gv->map_str_S,
                         gv->map_str_N))
    {
      gs_info(gv,-2,
              "Invalid linkage map name or class storage");
      goto cleanup;
    }

  /*
   * map_S is storage. nr identifies the marker-data row represented
   * by each storage element. The storage order itself need not equal
   * either map order or marker-data order.
   */
  for (k=0; k < gv->NLoc; k++)
    {
      m = gv->map_S[k].nr;

      if ((m >= gv->NoMar) ||
          (gs_UNSIGNED_MAX != storage[m]))
        {
          gs_info(gv,-2,
                  "Invalid linkage-map data-row index");
          goto cleanup;
        }

      storage[m] = k;
    }

  /*
   * map[] is access order. mdl[] gives the marker-data row for that
   * map point. nr in the pointed-to storage element must agree.
   */
  for (lne=0; lne < gv->NLoc; lne++)
    {
      if (NULL == gv->map[lne])
        {
          gs_info(gv,-2,"Invalid linkage map pointer");
          goto cleanup;
        }

      m = gv->mdl[lne];

      if ((m >= gv->NoMar) ||
          (gs_UNSIGNED_MAX == storage[m]))
        {
          gs_info(gv,-2,
                  "Invalid map to marker-data index");
          goto cleanup;
        }

      if (gv->map[lne] !=
          &(gv->map_S[storage[m]]))
        {
          gs_info(gv,-2,
                  "Linkage map pointer does not match map storage");
          goto cleanup;
        }

      if (gv->map[lne]->nr != m)
        {
          gs_info(gv,-2,
                  "Inconsistent map to marker-data index");
          goto cleanup;
        }

      if (NULL == memchr(gv->MarName[m],
                         '\0',
                         gs_MAXNME))
        {
          gs_info(gv,-2,"Invalid marker name");
          goto cleanup;
        }

      if (0 != strcmp(gv->map[lne]->name,
                      gv->MarName[m]))
        {
          gs_info(gv,-2,
                  "Linkage-map marker is not matching marker data");
          goto cleanup;
        }

      if (seen[m])
        {
          gs_info(gv,-2,
                  "Marker-data row used more than once by linkage map");
          goto cleanup;
        }

      seen[m] = 1;
    }

  if (require_complete)
    {
      for (m=0; m < gv->NoMar; m++)
        {
          if (!seen[m])
            {
              gs_info(gv,-2,
                      "Marker data are not completely covered by map");
              goto cleanup;
            }
        }
    }

  if (require_unique_names)
    {
      if (gv->NoMar >
          (gs_UNSIGNED)(((size_t)-1) / sizeof(*names)))
        {
          gs_info(gv,-2,"Memory allocation error");
          goto cleanup;
        }

      names = malloc((size_t)gv->NoMar * sizeof(*names));

      if (NULL == names)
        {
          gs_info(gv,-2,"Memory allocation error");
          goto cleanup;
        }

      for (m=0; m < gv->NoMar; m++)
        names[m] = gv->MarName[m];

      qsort(names,
            gv->NoMar,
            sizeof(*names),
            gs_string_cmp);

      for (m=1; m < gv->NoMar; m++)
        {
          if (0 == strcmp(names[m-1], names[m]))
            {
              snprintf(gv->msg,
                       sizeof(gv->msg),
                       "Duplicate marker name %s in marker data",
                       names[m]);
              gs_info(gv,-2,gv->msg);
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


/* qsort strings  comparison function*/
int gs_strcmp(const void *a, const void *b) 
{ 
    const char *ia = (const char *)a;
    const char *ib = (const char *)b;
    return strcmp(ia, ib);
} 

/* qsort int comparison function */ 
int gs_intcmp(const void *a, const void *b) 
{ 
    const int *ia = (const int *)a; 
    const int *ib = (const int *)b;
    return *ia  - *ib; 
} 

/* qsort short comparison function */ 
int gs_shortcmp(const void *a, const void *b) 
{ 
    const short *ia = (const short *)a; 
    const short *ib = (const short *)b;
    return *ia  - *ib; 
} 

/* qsort double comparison function */ 
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


void gs_marker_stats_01 ( gs_varset_type* , char** , char** , char** ,
			  int* , int* , int* , int* , int* );

long long int st_getline (char** , size_t* , FILE* ) ;

void gs_remove_newline(char*);

int gs_determine_alleles (GS_ALLELE_TYPE*, GS_ALLELE_TYPE*, const char*);


static int gs_finish_marker_read(gs_varset_type *gv,
                                 gs_varset_type *tmp,
                                 char **ind_filename,
                                 char **mar_filename,
                                 char **gen_filename,
                                 int *auxfiles,
                                 int *retval)
{
  int srv = 0;
  int *sretval = &srv;

  gs_marker_stats_01(tmp,
                     ind_filename,
                     mar_filename,
                     gen_filename,
                     auxfiles, auxfiles, auxfiles,
                     auxfiles,
                     sretval);

  if (-2 == *sretval)
    {
      *retval = -2;
      return -1;
    }

  gs_tmp_commit(gv, tmp);
  *retval = 0;
  return 0;
}

void gs_read_marker_data_01 ( gs_varset_type *gv ,
                              char** dta_filename,
                              char** ind_filename,
                              char** mar_filename,
                              char** gen_filename,
                              int*   auxfiles,
                              int*   retval        )
{
  gs_varset_type tmp;
  FILE *fp = NULL;
  char *line = NULL;
  size_t n = 0;
  long long int rc;
  gs_UNSIGNED NoLines = 0;
  gs_UNSIGNED IndSize = 0;
  gs_UNSIGNED MarSize = 0;
  char (*ind)[gs_MAXNME] = NULL;
  char (*mar)[gs_MAXNME] = NULL;
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
      gs_info(gv,-2,"File open error");
      *retval = -2;
      return;
    }

  rc = st_getline(&line, &n, fp);
  if (rc < 0)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (0 == rc)
    {
      gs_info(gv,-2,"Marker data file contains no header");
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
          gs_info(gv,-2,"Empty line in marker data file");
          failed = 1;
          goto cleanup;
        }
      len = strlen(p);
      if (len >= gs_MAXNME)
        {
          gs_info(gv,-2,"Individual name too long");
          failed = 1;
          goto cleanup;
        }
      if (len > IndSize) IndSize = (gs_UNSIGNED)len;

      p = strtok(NULL, gs_tok_del);
      if (NULL == p)
        {
          gs_info(gv,-2,"Marker name missing in marker data file");
          failed = 1;
          goto cleanup;
        }
      len = strlen(p);
      if (len >= gs_MAXNME)
        {
          gs_info(gv,-2,"Marker name too long");
          failed = 1;
          goto cleanup;
        }
      if (len > MarSize) MarSize = (gs_UNSIGNED)len;

      if (gs_UNSIGNED_MAX == NoLines)
        {
          gs_info(gv,-2,"Marker data file too large");
          failed = 1;
          goto cleanup;
        }
      NoLines++;
    }

  if (rc < 0)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  fclose(fp);
  fp = NULL;

  if (0 == NoLines)
    {
      gs_info(gv,-2,"Marker data file contains no data lines");
      failed = 1;
      goto cleanup;
    }

  ind = gs_malloc_array((size_t)NoLines,sizeof(*ind));
  mar = gs_malloc_array((size_t)NoLines,sizeof(*mar));
  rec_a = gs_malloc_array((size_t)NoLines,sizeof(*rec_a));
  rec_b = gs_malloc_array((size_t)NoLines,sizeof(*rec_b));
  rec_pair = gs_malloc_array((size_t)NoLines,sizeof(*rec_pair));
  names = gs_malloc_array((size_t)NoLines,sizeof(*names));

  if ((NULL == ind) || (NULL == mar) ||
      (NULL == rec_a) || (NULL == rec_b) ||
      (NULL == rec_pair) || (NULL == names))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  fp = fopen(*dta_filename, "r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      failed = 1;
      goto cleanup;
    }

  rc = st_getline(&line, &n, fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,"Error reading marker data file");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED j = 0; j < NoLines; j++)
    {
      char *p;

      rc = st_getline(&line, &n, fp);
      if (rc <= 0)
        {
          gs_info(gv,-2,"Unexpected end of marker data file");
          failed = 1;
          goto cleanup;
        }

      gs_remove_newline(line);
      p = strtok(line, gs_tok_del);
      if (NULL == p)
        {
          gs_info(gv,-2,"Error in marker data file");
          failed = 1;
          goto cleanup;
        }
      strcpy(ind[j], p);

      p = strtok(NULL, gs_tok_del);
      if (NULL == p)
        {
          gs_info(gv,-2,"Error in marker data file");
          failed = 1;
          goto cleanup;
        }
      strcpy(mar[j], p);

      p = strtok(NULL, gs_tok_del);
      if (NULL == p)
        {
          gs_info(gv,-2,"Missing allele field in marker data file");
          failed = 1;
          goto cleanup;
        }

      rec_pair[j] = (NULL != strchr(p,'/'));
      if (rec_pair[j])
        {
          if (0 != gs_determine_alleles(&(rec_a[j]), &(rec_b[j]), p))
            {
              gs_info(gv,-2,"Invalid allele in marker data file");
              failed = 1;
              goto cleanup;
            }
        }
      else
        {
          if (0 != gs_parse_single_allele(p,&(rec_a[j])))
            {
              gs_info(gv,-2,"Invalid allele in marker data file");
              failed = 1;
              goto cleanup;
            }
          rec_b[j] = rec_a[j];
        }
    }

  fclose(fp);
  fp = NULL;

  for (gs_UNSIGNED j = 0; j < NoLines; j++) names[j] = ind[j];
  qsort(names, NoLines, sizeof(*names), gs_string_cmp);
  tmp.NoInd = 1;
  for (gs_UNSIGNED j = 1; j < NoLines; j++)
    if (0 != strcmp(names[j-1], names[j])) tmp.NoInd++;

  if ((0 == tmp.NoInd) ||
      ((size_t)tmp.NoInd > ((size_t)-1) / sizeof(*tmp.IndName)))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  tmp.IndName = malloc((size_t)tmp.NoInd * sizeof(*tmp.IndName));
  if (NULL == tmp.IndName)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  {
    gs_UNSIGNED k = 0;
    strcpy(tmp.IndName[k++], names[0]);
    for (gs_UNSIGNED j = 1; j < NoLines; j++)
      if (0 != strcmp(names[j-1], names[j]))
        strcpy(tmp.IndName[k++], names[j]);
  }

  for (gs_UNSIGNED j = 0; j < NoLines; j++) names[j] = mar[j];
  qsort(names, NoLines, sizeof(*names), gs_string_cmp);
  tmp.NoMar = 1;
  for (gs_UNSIGNED j = 1; j < NoLines; j++)
    if (0 != strcmp(names[j-1], names[j])) tmp.NoMar++;

  if ((0 == tmp.NoMar) ||
      ((size_t)tmp.NoMar > ((size_t)-1) / sizeof(*tmp.MarName)))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  tmp.MarName = malloc((size_t)tmp.NoMar * sizeof(*tmp.MarName));
  if (NULL == tmp.MarName)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  {
    gs_UNSIGNED k = 0;
    strcpy(tmp.MarName[k++], names[0]);
    for (gs_UNSIGNED j = 1; j < NoLines; j++)
      if (0 != strcmp(names[j-1], names[j]))
        strcpy(tmp.MarName[k++], names[j]);
  }

  tmp.MaxIL = (int)IndSize;
  tmp.MaxML = (int)MarSize;

  if (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  gs_set_genotypes_missing(&tmp);

  {
    size_t n_genotypes;
    if (0 != gs_size_mul((size_t)tmp.NoMar,
                         (size_t)tmp.NoInd,
                         &n_genotypes))
      {
        gs_info(gv,-2,"Memory allocation error");
        failed = 1;
        goto cleanup;
      }
    geno_state = calloc(n_genotypes,sizeof(*geno_state));
    if (NULL == geno_state)
      {
        gs_info(gv,-2,"Memory allocation error");
        failed = 1;
        goto cleanup;
      }
  }

  for (gs_UNSIGNED k = 0; k < NoLines; k++)
    {
      gs_UNSIGNED m = gs_find_sorted_fixed_name(tmp.MarName,
                                                tmp.NoMar,
                                                mar[k]);
      gs_UNSIGNED i = gs_find_sorted_fixed_name(tmp.IndName,
                                                tmp.NoInd,
                                                ind[k]);
      size_t idx;
      GS_ALLELE_TYPE *a1;
      GS_ALLELE_TYPE *a2;

      if ((gs_UNSIGNED_MAX == m) || (gs_UNSIGNED_MAX == i))
        {
          gs_info(gv,-2,"Internal marker-data name lookup error");
          failed = 1;
          goto cleanup;
        }

      idx = (size_t)m * (size_t)tmp.NoInd + (size_t)i;
      a1 = &(tmp.All1[idx]);
      a2 = &(tmp.All2[idx]);

      if (0 != gs_aggregate_genotype(a1,a2,&(geno_state[idx]),
                                     rec_a[k],rec_b[k],rec_pair[k]))
        {
          gs_info(gv,-2,
                  rec_pair[k] ?
                  "Conflicting genotype records in list marker data" :
                  "More than two alleles for one list genotype");
          failed = 1;
          goto cleanup;
        }
    }

  if (0 != gs_finish_marker_read(gv, &tmp,
                                 ind_filename,
                                 mar_filename,
                                 gen_filename,
                                 auxfiles,
                                 retval))
    failed = 1;

cleanup:
  if (NULL != fp) fclose(fp);
  free(line);
  free(ind);
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



void gs_read_marker_data_01_GV ( char** dta_filename,
			      char** ind_filename,
			      char** mar_filename,
			      char** gen_filename,
			      int*   auxfiles,
			      int*   retval,
			      char** set_name)
{
  gs_read_marker_data_01 ( gs_fdta(*set_name),
			dta_filename,
			ind_filename,
			mar_filename,
			gen_filename,
			auxfiles,
			retval         ) ;

}


/* Design matrix of type 1: */
/* biallelic markers */
/* one column per marker */
/* design matrix consists of counts of the greater allele */
/* greater allele is the 'effect allele' */
/* missing equals heterozyous */

/* effect allele:        b */
/* complementray allele: a */
/* everything else:      c */

/* c is missing or alleles that are present in the prediction set but not */
/* the estimation set */

/* Values of the design matrix */

/* a:    0        */
/* c:   1/2       */
/* b:    1        */

/* aa   ab    bb  */
/* 0     1     2  */

/* ac   cc    bc   */
/* 1/2   1    3/2  */

/* Storage of the complementray allele is required to distinguish  */
/* between  */
/* aa ac cc and */
/* ab bc */

void gs_build_Z_01 ( gs_varset_type *gv ,
		     char** out_filename,
		     int*   auxfiles,
		     int*   retval)
{

  gs_UNSIGNED i, m, idx_all, idx_Z;
  FILE *outfp;

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available");
    *retval = -2; return ;
  }

  /* Allocate memory for the design matrix */
  gv->Z =  realloc( gv->Z, gv->NoInd * gv->NoMar * sizeof(double) );

  /* Check wheter EffAlls are already defined */

//  int EffAllDefined = 1;
//  if ( NULL ==  gv->EffAll ) EffAllDefined = 0; 

  /* Allocate memory for the effect descriptions */
  gv->NoEff =  1 + gv->NoMar;
  gv->EffNme = realloc (gv->EffNme, gv->NoEff * sizeof(char[gs_MAXNME]) );
  gv->EffAll = realloc (gv->EffAll, gv->NoEff * sizeof(int) );
  gv->CmpAll = realloc (gv->CmpAll, gv->NoEff * sizeof(int) );

  if ( (NULL==gv->Z) || (NULL==gv->EffNme) || (NULL==gv->EffAll) 
       || (NULL==gv->CmpAll)) ERR_M;

  /* Names of effects */
  strcpy(gv->EffNme[0],"mu");
  gv->EffAll[0] = 0;
  gv->CmpAll[0] = 0;
  for ( m = 1; m < gv->NoEff; m++)
    strcpy(gv->EffNme[m],gv->MarName[m-1]);		       
  
  /* Initialize design matrix */

  for ( i = 0; i < gv->NoInd ; i++) 
    for ( m = 0; m < gv->NoMar; m++) {
      idx_Z   = i*gv->NoMar + m ;
      gv->Z[idx_Z] =  gv->aA;
    }

  int a, b, lar_all, sml_all,
    mono_warning ;

  mono_warning = 0;

  for ( m = 0; m < gv->NoMar; m++)
    { 
      if ( NULL == gv->u ) {
	// Use the larger allele as EffAll
	lar_all = 0;
	for ( i = 0; i < gv->NoInd ; i++){
	  idx_all = m*gv->NoInd + i ;
          if ((GS_ALLELE_MISSING != gv->All1[idx_all]) &&
              (gv->All1[idx_all] > lar_all)) lar_all = gv->All1[idx_all];
          if ((GS_ALLELE_MISSING != gv->All2[idx_all]) &&
              (gv->All2[idx_all] > lar_all)) lar_all = gv->All2[idx_all];
	}
	/* Determine the complementary allele for each marker */
	sml_all = lar_all ;
	for ( i = 0; i < gv->NoInd ; i++){
	  idx_all = m*gv->NoInd + i ;
          if ((GS_ALLELE_MISSING != gv->All1[idx_all]) &&
              (gv->All1[idx_all] < lar_all) && (gv->All1[idx_all] > 0))
            sml_all = gv->All1[idx_all];
          if ((GS_ALLELE_MISSING != gv->All2[idx_all]) &&
              (gv->All2[idx_all] < lar_all) && (gv->All2[idx_all] > 0))
            sml_all = gv->All2[idx_all];
	}
	/* Check whether marker is biallelic */
	if ( (lar_all==sml_all) || (lar_all<2)  ) mono_warning = 1;
        /* Determine genotype */
	for ( i = 0; i < gv->NoInd ; i++){
	  idx_all = m*gv->NoInd + i ;
	  if ( ( (gv->All1[idx_all] != lar_all) && 
		 (gv->All1[idx_all] != sml_all) && 
		 (gv->All1[idx_all] != GS_ALLELE_MISSING) )          ||
	       ( (gv->All2[idx_all] != lar_all) && 
		 (gv->All2[idx_all] != sml_all) && 
		 (gv->All2[idx_all] != GS_ALLELE_MISSING) )              )
	    {
	      snprintf(gv->msg, sizeof(gv->msg),
		       "More than two alleles: %s",
		       gv->MarName[m]);

	      gs_info(gv,-2,gv->msg);
	      *retval=-2;return;
	    }
	}
	/* Save allele for which effects will be determined */
	gv->EffAll[m+1] = lar_all;
	gv->CmpAll[m+1] = sml_all;
      } else {
	// Use EffAll and CmpAll from previous run 
	lar_all = gv->EffAll[m+1] ;
	sml_all = gv->CmpAll[m+1] ;
      }

     /* Construct design matrix */
      for ( i = 0; i < gv->NoInd ; i++)
	{
	  idx_all = m*gv->NoInd + i ;

	  a = gv->All1[idx_all] ;
          b = gv->All2[idx_all] ;

	  idx_Z   = i*gv->NoMar + m ; 

	  if ( (a==lar_all) && (b==lar_all) ) { 
	    gv->Z[idx_Z] = gv->AA; 
	  }
	  else	{
	    if ( (a==sml_all) && (b==sml_all) )  {
	      gv->Z[idx_Z] = gv->aa;
	    }
	    else {
	      gv->Z[idx_Z] = gv->aA ;
	    }
	  }
	}
    }
  
   
  /* Write design matrix to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( m = 0; m < gv->NoMar; m++)
	fprintf(outfp,"%s.%i ",gv->EffNme[m+1],gv->EffAll[m+1]);
      fprintf(outfp,"\n");

      for ( i = 0; i < gv->NoInd ; i++)
	{
	  fprintf(outfp,"%s ",gv->IndName[i]);
	  for ( m = 0; m < gv->NoMar; m++)
	    {
	      idx_Z = i*gv->NoMar + m ;
	      fprintf(outfp,"%f ",gv->Z[idx_Z]);
	    }
	  fprintf(outfp,"\n");
	}

      fclose(outfp);  
  }

  if ( 1 == mono_warning  )
    {

      snprintf(gv->msg, sizeof(gv->msg),
	       "Monomorphic markers in the design matrix");

      gs_info(gv,-1,gv->msg);
      *retval = -2; return;
    }
  
  *retval = 0; return;

}

void gs_build_Z_01_GV ( 
		       char** out_filename        ,
		       int*   auxfiles            ,
		       int*   retval              ,
		       char** set_name                               
		      )
{
  gs_build_Z_01 ( gs_fdta(*set_name) ,
		  out_filename,
		  auxfiles,
		  retval         );
}

void gs_neg_effall_01 (  gs_varset_type *gv ,
			 int*    retval         )
{
  gs_UNSIGNED m;
  int c;

  if ( NULL == gv->u ) {
      gs_info(gv,-2, "No effects estimated");
      *retval = -2; return ;
  }

  for ( m = 1; m < gv->NoEff; m++) {
    if ( gv->u[m] > 0) {
      c = gv->EffAll[m];
      gv->EffAll[m] = gv->CmpAll[m];
      gv->CmpAll[m] = c;
    }
    gv->u[m] = 0;
  }

  *retval = 0; return;
}

void gs_neg_effall_01_GV (  int*    retval,
			    char** set_name         )
{
  gs_neg_effall_01 ( gs_fdta(*set_name),
		     retval         ) ;
}

void gs_pos_effall_01 (  gs_varset_type *gv ,
			 int*    retval         )
{
  gs_UNSIGNED m;
  int c;

  if ( NULL == gv->u ) {
      gs_info(gv,-2, "No effects estimated");
      *retval = -2; return ;
  }

  for ( m = 1; m < gv->NoEff; m++) {
    if ( gv->u[m] < 0) {
      c = gv->EffAll[m];
      gv->EffAll[m] = gv->CmpAll[m];
      gv->CmpAll[m] = c;
    }
    gv->u[m] = 0;
  }

  *retval = 0; return;
}

void gs_pos_effall_01_GV (  int*    retval,
			    char** set_name         )
{
  gs_pos_effall_01 ( gs_fdta(*set_name),
		     retval         ) ;
}


long long int st_getline(char **, size_t *, FILE *);

void gs_remove_newline(char*);

void gs_read_performance_data_01 ( gs_varset_type *gv ,
                                   char** dta_filename,
                                   char** out_filename,
                                   int*   auxfiles,
                                   int*   retval         )
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

  if ((0 == gv->NoInd) || (0 == gv->NoMar) ||
      (NULL == gv->IndName) || (NULL == gv->MarName) ||
      (NULL == gv->All1) || (NULL == gv->All2))
    {
      gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
      *retval = -2;
      return;
    }

  gs_tmp_init_from(&tmp, gv);

  if (0 != gs_make_name_index(gv->IndName, gv->NoInd, &index))
    {
      gs_info(gv,-2,"Memory allocation error");
      *retval = -2;
      return;
    }

  if (gs_name_index_has_duplicates(index, gv->NoInd))
    {
      gs_info(gv,-2,"Duplicate individual name in marker data");
      free(index);
      *retval = -2;
      return;
    }

  values = gs_malloc_array((size_t)gv->NoInd,sizeof(*values));
  seen = calloc((size_t)gv->NoInd, sizeof(*seen));
  if ((NULL == values) || (NULL == seen))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  fp = fopen(*dta_filename, "r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      failed = 1;
      goto cleanup;
    }

  rc = st_getline(&line, &n, fp);
  if (rc < 0)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (0 == rc)
    {
      gs_info(gv,-2,"Performance data file contains no header");
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
          gs_info(gv,-2,"Empty line in performance data file");
          failed = 1;
          goto cleanup;
        }

      i = gs_name_index_find(index, gv->NoInd, name);
      if (gs_UNSIGNED_MAX == i)
        continue;

      value_text = strtok(NULL, gs_tok_del);
      if (NULL == value_text)
        {
          gs_info(gv,-2,"Performance value missing");
          failed = 1;
          goto cleanup;
        }

      {
        char *end;
        double value = strtod(value_text, &end);

        if ((end == value_text) || ('\0' != *end) || !isfinite(value))
          {
            gs_info(gv,-2,"Invalid performance value");
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
    }

  if (rc < 0)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  fclose(fp);
  fp = NULL;

  if (0 == no_seen)
    {
      gs_info(gv,-2,"No performance data match marker-data individuals");
      failed = 1;
      goto cleanup;
    }

  tmp.NoInd = no_seen;
  tmp.NoMar = gv->NoMar;
  tmp.MaxIL = gv->MaxIL;
  tmp.MaxML = gv->MaxML;

  tmp.IndName = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.IndName));
  tmp.MarName = gs_malloc_array((size_t)tmp.NoMar,sizeof(*tmp.MarName));
  tmp.y = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.y));

  if ((NULL == tmp.IndName) ||
      (NULL == tmp.MarName) ||
      (NULL == tmp.y) ||
      (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd)))
    {
      gs_info(gv,-2,"Memory allocation error");
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
          strcpy(tmp.IndName[j], gv->IndName[i]);
          tmp.y[j] = values[i];

          for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            {
              size_t old_idx = (size_t)m * (size_t)gv->NoInd + i;
              size_t new_idx = (size_t)m * (size_t)tmp.NoInd + j;
              tmp.All1[new_idx] = gv->All1[old_idx];
              tmp.All2[new_idx] = gv->All2[old_idx];
            }

          j++;
        }
  }

  {
    int map_rc = gs_copy_map_primary(&tmp, gv);
    if (0 != map_rc)
      {
        gs_info(gv,-2,
                (-1 == map_rc) ?
                "Memory allocation error" :
                "Invalid linkage map state");
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

    gs_marker_stats_01(&tmp,
                       empty_p, empty_p, empty_p,
                       zero_p, zero_p, zero_p,
                       zero_p,
                       &srv);
    if (-2 == srv)
      {
        failed = 1;
        goto cleanup;
      }
  }

  snprintf(gv->msg, sizeof(gv->msg),
           "No. of performance data: %lu", tmp.NoInd);
  gs_info(gv,0,gv->msg);

  if (1 == *auxfiles)
    {
      FILE *outfp = fopen(*out_filename, "w");
      if (NULL == outfp)
        {
          gs_info(gv,-2,"File open error");
          failed = 1;
          goto cleanup;
        }

      fprintf(outfp,"y\n");
      for (gs_UNSIGNED i = 0; i < tmp.NoInd; i++)
        fprintf(outfp,"%s %f\n", tmp.IndName[i], tmp.y[i]);
      fclose(outfp);
    }

  gs_tmp_commit(gv, &tmp);
  *retval = 0;

cleanup:
  if (NULL != fp) fclose(fp);
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


void gs_read_performance_data_01_GV (  char** dta_filename,
				       char** out_filename,
				       int*   auxfiles,
				       int*   retval,
				       char** set_name   )
{
  gs_read_performance_data_01 ( gs_fdta(*set_name),
				dta_filename,
				out_filename,
				auxfiles,
				retval ) ;
}

void gs_return_performance_data_01 ( gs_varset_type *gv ,
				     char** out_filename,
				     int*   auxfiles,
				     int*   retval         )
{
  gs_UNSIGNED i;

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return;
  }

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  /* Write observation vector to output file */
  if (1 == *auxfiles) {

    FILE *outfp;

    if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

    fprintf(outfp,"y\n");
    for ( i = 0; i < gv->NoInd ; i++)
      fprintf(outfp,"%s %f\n",gv->IndName[i], gv->y[i]);
    
    fclose(outfp);  
  }

  *retval = 0; return ;

}

void gs_return_performance_data_01_GV ( char** out_filename,
					int*   auxfiles,
					int*   retval,
					char** set_name   )
{
  gs_return_performance_data_01 ( gs_fdta(*set_name),
				  out_filename,
				  auxfiles,
				  retval ) ;
}

void gs_lambda_const_01 ( gs_varset_type *gv ,
			  double* lambda, 
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  gs_UNSIGNED i,m;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available");
    *retval = -2; return ;
  }

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if ( NULL== gv->ladi )  ERR_M;

  for ( i = 0 ; i < gv->NoMar ; i++)  gv->ladi[i] = * lambda;

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;

      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    } 

  *retval = 0; return ;

}

void gs_lambda_const_01_GV ( double* lambda      , 
			     char**  out_filename,
			     int*    auxfiles    , 
			     int*    retval      ,
			     char**  set_name     )

{
  gs_lambda_const_01 ( gs_fdta(*set_name),
		       lambda, 
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}

void gs_lambda_var_01 ( gs_varset_type *gv ,
			double* lambda,
			int*    n,
			char**  out_filename,
			int*    auxfiles, 
			int*    retval         )
{
  gs_UNSIGNED m;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available");
    *retval = -2; return ;
  }

  if ( (gs_UNSIGNED)*n != gv->NoMar ) {
    gs_info(gv,-2,"Incorrect number of lambda values");
    *retval = -2; return ;
  }
  
  
  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if ( NULL== gv->ladi )  ERR_M;

  memcpy(gv->ladi,lambda,gv->NoMar * sizeof(double));


  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;

      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    } 

  *retval = 0; return ;

}

void gs_lambda_var_01_GV ( double* lambda      ,
			   int*    n           ,
			   char**  out_filename,
			   int*    auxfiles    , 
			   int*    retval      ,
			   char**  set_name     )

{
  gs_lambda_var_01 ( gs_fdta(*set_name),
		     lambda,
		     n,
		     out_filename,
		     auxfiles, 
		     retval         ) ;
}




void gs_build_V_01 ( gs_varset_type* gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  if ( NULL == gv->Z ) {
      int Srv; int* Sretval = &Srv;
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }
  }

  gv->V =  realloc( gv->V, gv->NoInd * gv->NoInd * sizeof(double) );
  if ( NULL == gv->V ) ERR_M;

  st_dgemm ("N","T",
	    gv->NoInd, gv->NoMar,
	    gv->NoInd, gv->NoMar,
	    gv->Z, gv->Z, gv->V);
  
  /* Write V matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->NoInd; i++) {
	  for ( j = 0; j < gv->NoInd; j++) {
	      idx_A   = i*gv->NoInd + j ;
	      fprintf(outfp,"%f ",gv->V [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_build_V_01_GV ( char** out_filename,
			  int*   auxfiles, 
			  int*   retval,
			  char** set_name                )
{
  gs_build_V_01 (gs_fdta(*set_name),
		   out_filename,
		   auxfiles, 
		   retval        ) ;
}


static int gs_restrict_read_name_file(
    gs_varset_type          *gv,
    const char              *filename,
    const gs_name_index_type *index,
    gs_UNSIGNED              n,
    unsigned char           *keep)
{
  FILE *fp = NULL;
  char *line = NULL;
  size_t line_size = 0;
  long long int rc;

  if ((NULL == filename) || ('\0' == filename[0]))
    return 0;

  memset(keep,0,(size_t)n * sizeof(*keep));

  fp = fopen(filename,"r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      return -1;
    }

  while (0 < (rc = st_getline(&line,&line_size,fp)))
    {
      char *p;

      gs_remove_newline(line);
      p = strtok(line,gs_tok_del);

      while (NULL != p)
        {
          gs_UNSIGNED idx = gs_name_index_find(index,n,p);

          if (gs_UNSIGNED_MAX != idx)
            keep[idx] = 1;

          p = strtok(NULL,gs_tok_del);
        }
    }

  fclose(fp);
  free(line);

  if (rc < 0)
    {
      gs_info(gv,-2,"Memory allocation error");
      return -1;
    }

  return 0;
}


static int gs_restrict_marker_effects_valid(const gs_varset_type *gv)
{
  gs_UNSIGNED m;

  if ((gs_UNSIGNED_MAX == gv->NoMar) ||
      (gv->NoEff != gv->NoMar + 1) ||
      (gv->dimu != gv->NoEff) ||
      (NULL == gv->EffNme) ||
      (NULL == gv->EffAll) ||
      (NULL == gv->CmpAll) ||
      (NULL == gv->u))
    return 0;

  if ((NULL == memchr(gv->EffNme[0],'\0',gs_MAXNME)) ||
      (0 != strcmp(gv->EffNme[0],"mu")))
    return 0;

  for (m=0; m < gv->NoMar; m++)
    {
      if ((NULL == memchr(gv->EffNme[m+1],'\0',gs_MAXNME)) ||
          (0 != strcmp(gv->EffNme[m+1],gv->MarName[m])))
        return 0;
    }

  return 1;
}


void gs_restrict_marker_data_03 ( gs_varset_type *gv ,
                                  char**  ind_fname,
                                  char**  mar_fname,
                                  int*    NoAllMAX,
                                  double* MaMisMAX,
                                  double* ExHetMIN,
                                  double* InMisMAX,
                                  int*    retval        )
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

  if ((0 == gv->NoInd) || (0 == gv->NoMar) ||
      (NULL == gv->IndName) || (NULL == gv->MarName) ||
      (NULL == gv->All1) || (NULL == gv->All2))
    {
      gs_info(gv,-2,
              "No marker data available, run 'read.marker.data' first");
      return;
    }

  if ((NULL == gv->NoAll) || (NULL == gv->All) ||
      (NULL == gv->Mis) || (NULL == gv->Pic) ||
      (NULL == gv->IndMis))
    {
      gs_info(gv,-2,"Marker statistics are not available");
      return;
    }

  if ((*NoAllMAX < 0) ||
      !isfinite(*MaMisMAX) ||
      !isfinite(*ExHetMIN) ||
      !isfinite(*InMisMAX))
    {
      gs_info(gv,-2,"Invalid marker restriction criterion");
      return;
    }

  keep_i = calloc((size_t)gv->NoInd,sizeof(*keep_i));
  keep_m = calloc((size_t)gv->NoMar,sizeof(*keep_m));

  if ((NULL == keep_i) || (NULL == keep_m))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  memset(keep_i,1,(size_t)gv->NoInd * sizeof(*keep_i));
  memset(keep_m,1,(size_t)gv->NoMar * sizeof(*keep_m));

  if ('\0' != (*ind_fname)[0])
    {
      if (0 != gs_make_name_index(gv->IndName,gv->NoInd,&ind_index))
        {
          gs_info(gv,-2,"Memory allocation error");
          failed = 1;
          goto cleanup;
        }

      if (0 != gs_restrict_read_name_file(gv,*ind_fname,
                                          ind_index,gv->NoInd,keep_i))
        {
          failed = 1;
          goto cleanup;
        }
    }

  if ('\0' != (*mar_fname)[0])
    {
      if (0 != gs_make_name_index(gv->MarName,gv->NoMar,&mar_index))
        {
          gs_info(gv,-2,"Memory allocation error");
          failed = 1;
          goto cleanup;
        }

      if (0 != gs_restrict_read_name_file(gv,*mar_fname,
                                          mar_index,gv->NoMar,keep_m))
        {
          failed = 1;
          goto cleanup;
        }
    }

  /* Apply all criteria to the statistics of the original data set. */
  for (gs_UNSIGNED i=0; i < gv->NoInd; i++)
    {
      if ((*InMisMAX < gv->IndMis[i]) || !keep_i[i])
        keep_i[i] = 0;
      else
        CountInd++;
    }

  for (gs_UNSIGNED m=0; m < gv->NoMar; m++)
    {
      gs_UNSIGNED no_all = gv->NoAll[m];

      if ((0 < no_all) &&
          (GS_ALLELE_MISSING == gv->All[m][no_all-1]))
        no_all--;

      if (((gs_UNSIGNED)(*NoAllMAX) < no_all) ||
          (*MaMisMAX < gv->Mis[m]) ||
          (*ExHetMIN > gv->Pic[m]) ||
          !keep_m[m])
        keep_m[m] = 0;
      else
        CountMar++;
    }

  if (0 == CountInd)
    {
      gs_info(gv,-2,"No individuals remain after marker-data restriction");
      failed = 1;
      goto cleanup;
    }

  if (0 == CountMar)
    {
      gs_info(gv,-2,"No markers remain after marker-data restriction");
      failed = 1;
      goto cleanup;
    }

  if ((CountInd == gv->NoInd) && (CountMar == gv->NoMar))
    {
      *retval = 0;
      goto cleanup;
    }

  no_marker_reduction = (int)(CountMar == gv->NoMar);
  effects_valid = no_marker_reduction &&
                  gs_restrict_marker_effects_valid(gv);

  has_map = ((0 != gv->NLoc) || (NULL != gv->map_S) ||
             (NULL != gv->map_str_S) || (0 != gv->map_str_N) ||
             (NULL != gv->map) || (NULL != gv->mdl));

  if (has_map)
    {
      map_storage = gs_map_storage_index(gv,0,0);
      if (NULL == map_storage)
        {
          failed = 1;
          goto cleanup;
        }
    }

  gs_tmp_init_from(&tmp,gv);
  tmp.NoInd = CountInd;
  tmp.NoMar = CountMar;

  tmp.IndName = gs_malloc_array((size_t)tmp.NoInd,
                                sizeof(*tmp.IndName));
  tmp.MarName = gs_malloc_array((size_t)tmp.NoMar,
                                sizeof(*tmp.MarName));

  if ((NULL == tmp.IndName) || (NULL == tmp.MarName) ||
      (0 != gs_alloc_genotype_matrix(&tmp,tmp.NoMar,tmp.NoInd)))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup_tmp;
    }

  old_to_new_m = gs_malloc_array((size_t)gv->NoMar,
                                 sizeof(*old_to_new_m));
  if (NULL == old_to_new_m)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup_tmp;
    }

  for (gs_UNSIGNED m=0; m < gv->NoMar; m++)
    old_to_new_m[m] = gs_UNSIGNED_MAX;

  {
    gs_UNSIGNED i_new = 0;

    for (gs_UNSIGNED i=0; i < gv->NoInd; i++)
      if (keep_i[i])
        {
          size_t len = strlen(gv->IndName[i]);
          strcpy(tmp.IndName[i_new],gv->IndName[i]);
          if ((int)len > tmp.MaxIL) tmp.MaxIL = (int)len;
          i_new++;
        }
  }

  {
    gs_UNSIGNED m_new = 0;

    for (gs_UNSIGNED m=0; m < gv->NoMar; m++)
      if (keep_m[m])
        {
          gs_UNSIGNED i_new = 0;
          size_t len = strlen(gv->MarName[m]);

          strcpy(tmp.MarName[m_new],gv->MarName[m]);
          if ((int)len > tmp.MaxML) tmp.MaxML = (int)len;
          old_to_new_m[m] = m_new;

          for (gs_UNSIGNED i=0; i < gv->NoInd; i++)
            if (keep_i[i])
              {
                size_t old_idx = (size_t)m * gv->NoInd + i;
                size_t new_idx = (size_t)m_new * tmp.NoInd + i_new;

                tmp.All1[new_idx] = gv->All1[old_idx];
                tmp.All2[new_idx] = gv->All2[old_idx];
                i_new++;
              }

          m_new++;
        }
  }

  if (NULL != gv->y)
    {
      gs_UNSIGNED i_new = 0;

      tmp.y = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.y));
      if (NULL == tmp.y)
        {
          gs_info(gv,-2,"Memory allocation error");
          failed = 1;
          goto cleanup_tmp;
        }

      for (gs_UNSIGNED i=0; i < gv->NoInd; i++)
        if (keep_i[i])
          tmp.y[i_new++] = gv->y[i];
    }

  if (has_map)
    {
      gs_UNSIGNED nloc = 0;

      for (gs_UNSIGNED lne=0; lne < gv->NLoc; lne++)
        if (keep_m[gv->mdl[lne]])
          nloc++;

      if (0 < nloc)
        {
          size_t map_str_offset = 0;

          for (gs_UNSIGNED lne=0; lne < gv->NLoc; lne++)
            if (keep_m[gv->mdl[lne]])
              if ((0 != gs_pool_add_string_size(&(tmp.map_str_N),
                                                gv->map[lne]->name)) ||
                  (0 != gs_pool_add_string_size(&(tmp.map_str_N),
                                                gv->map[lne]->class)))
                {
                  gs_info(gv,-2,"Linkage map string storage too large");
                  failed = 1;
                  goto cleanup_tmp;
                }

          tmp.map_S = gs_malloc_array((size_t)nloc,sizeof(*tmp.map_S));
          tmp.map_str_S = gs_malloc_array(tmp.map_str_N,
                                          sizeof(*tmp.map_str_S));
          tmp.map = gs_malloc_array((size_t)nloc,sizeof(*tmp.map));
          tmp.mdl = gs_malloc_array((size_t)nloc,sizeof(*tmp.mdl));

          if ((NULL == tmp.map_S) || (NULL == tmp.map_str_S) ||
              (NULL == tmp.map) || (NULL == tmp.mdl))
            {
              gs_info(gv,-2,"Memory allocation error");
              failed = 1;
              goto cleanup_tmp;
            }

          tmp.NLoc = nloc;
          nloc = 0;

          for (gs_UNSIGNED lne=0; lne < gv->NLoc; lne++)
            {
              gs_UNSIGNED old_m = gv->mdl[lne];

              if (keep_m[old_m])
                {
                  gs_UNSIGNED new_m = old_to_new_m[old_m];

                  if (0 != gs_copy_map_point_to_pool(&(tmp.map_S[nloc]),
                                                      gv->map[lne],
                                                      new_m,
                                                      tmp.map_str_S,
                                                      tmp.map_str_N,
                                                      &map_str_offset))
                    {
                      gs_info(gv,-2,"Invalid linkage map string storage");
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
              gs_info(gv,-2,"Invalid linkage map string storage");
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

    gs_marker_stats_01(&tmp,
                       empty_p,empty_p,empty_p,
                       zero_p,zero_p,zero_p,zero_p,
                       &srv);
    if (-2 == srv)
      {
        failed = 1;
        goto cleanup_tmp;
      }
  }

  if (effects_valid)
    {
      tmp.EffNme = gs_malloc_array((size_t)gv->NoEff,
                                   sizeof(*tmp.EffNme));
      tmp.EffAll = gs_malloc_array((size_t)gv->NoEff,
                                   sizeof(*tmp.EffAll));
      tmp.CmpAll = gs_malloc_array((size_t)gv->NoEff,
                                   sizeof(*tmp.CmpAll));
      tmp.u = gs_malloc_array((size_t)gv->dimu,sizeof(*tmp.u));

      if ((NULL == tmp.EffNme) || (NULL == tmp.EffAll) ||
          (NULL == tmp.CmpAll) || (NULL == tmp.u))
        {
          gs_info(gv,-2,"Memory allocation error");
          failed = 1;
          goto cleanup_tmp;
        }

      tmp.NoEff = gv->NoEff;
      tmp.dimu = gv->dimu;
      memcpy(tmp.EffNme,gv->EffNme,
             (size_t)gv->NoEff * sizeof(*tmp.EffNme));
      memcpy(tmp.EffAll,gv->EffAll,
             (size_t)gv->NoEff * sizeof(*tmp.EffAll));
      memcpy(tmp.CmpAll,gv->CmpAll,
             (size_t)gv->NoEff * sizeof(*tmp.CmpAll));
      memcpy(tmp.u,gv->u,(size_t)gv->dimu * sizeof(*tmp.u));

      tmp.aa = gv->aa;
      tmp.aA = gv->aA;
      tmp.AA = gv->AA;
      tmp.an = gv->an;
      tmp.nn = gv->nn;
      tmp.nA = gv->nA;

      gs_info(gv,1,"Effects retained");
    }
  else if (!no_marker_reduction &&
           ((NULL != gv->u) || (NULL != gv->EffNme) ||
            (NULL != gv->EffAll) || (NULL != gv->CmpAll)))
    {
      gs_info(gv,1,"Effects removed");
    }

  gs_tmp_commit(gv,&tmp);
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

void gs_restrict_marker_data_03_GV ( char**  ind_fname,
				     char**  mar_fname,
				     int*    NoAllMAX,
				     double* MaMisMAX,
				     double* ExHetMIN,
				     double* InMisMAX, 
				     int*    retval,
				     char**  set_name       )
{
  gs_restrict_marker_data_03 ( gs_fdta(*set_name) ,
			       ind_fname,
			       mar_fname,
			       NoAllMAX,
			       MaMisMAX,
			       ExHetMIN,
			       InMisMAX, 
			       retval        );
}


void gs_copy_marker_data_01 ( gs_varset_type *new ,
			      gs_varset_type *gv  ,
			      int*   retval         )
{
  gs_UNSIGNED i,m;

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  if (new == gv)
    {
      *retval = 0;
      return;
    }

  gs_reset(new);

  /* Marker matrix */

  new->NoInd = gv->NoInd; 
  new->NoMar = gv->NoMar; 

  new->IndName = realloc( new->IndName, 
			  new->NoInd * sizeof(char[gs_MAXNME]) ); 
  new->MarName = realloc( new->MarName, 
			  new->NoMar * sizeof(char[gs_MAXNME]) ); 
  new->All1   =  realloc( new->All1,
                          new->NoMar * new->NoInd * sizeof(*new->All1) );
  new->All2   =  realloc( new->All2,
                          new->NoMar * new->NoInd * sizeof(*new->All2) );

  if ( (NULL == new->All1)  || (NULL == new->All2) || 
       ( NULL == new->IndName )  || ( NULL == new->MarName ) ) ERR_M;

  for ( i=0; i < gv->NoInd; i++  )
	strcpy ( new->IndName[i] , gv->IndName[i]);

  for ( m=0; m < gv->NoMar; m++ )
	strcpy ( new->MarName[m] , gv->MarName[m] );

  GS_ALLELE_TYPE *p1;
  GS_ALLELE_TYPE *p2;

  for ( m=0; m < gv->NoMar; m++ )
    {
	  p1 = &(gv ->All1[m* gv->NoInd]); 
	  p2 = &(new->All1[m*new->NoInd]); 
	  memcpy(p2,p1,(size_t)new->NoInd * sizeof(*p2));
	  p1 = &(gv ->All2[m* gv->NoInd]);
	  p2 = &(new->All2[m*new->NoInd]);
	  memcpy(p2,p1,(size_t)new->NoInd * sizeof(*p2));
    }
  
  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv;
  char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( new ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }

  {
    int map_rc = gs_copy_map_primary(new,gv);

    if (-1 == map_rc)
      ERR_M;
    if (-2 == map_rc)
      {
        *retval = -2;
        return;
      }
  }

  if (NULL != gv->y)  // Phenotypes
    {
      new->y = realloc( new->y, new->NoInd*sizeof(double) );
      if ( NULL== new->y ) ERR_M;
      memcpy( new->y ,
	      gv->y  , 
	      (size_t)new->NoInd * sizeof(double) );
    }

  *retval = 0; return;
}

void gs_copy_marker_data_01_GV ( int*    retval     ,
				 char**  set_name_t ,
				 char**  set_name_s   )
{
  gs_copy_marker_data_01 ( gs_fdta(*set_name_t) ,
			   gs_fdta(*set_name_s) , 
			   retval        );
}


void gs_estimate_bv_01 ( gs_varset_type *effect_set ,
			 gs_varset_type *marker_set ,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval          )
{
   if ( NULL == effect_set->u ) {
    gs_info(effect_set,-2,"No estimated effects available");
    *retval = -2; return;
   }
 
   if ( ( 0 == marker_set->NoInd) || ( 0 == marker_set->NoMar) ) {
     gs_info(marker_set,-2,"No marker data available");
     *retval = -2; return ;
   }
   gs_UNSIGNED i, m, idx_all, idx_Z;

   // The genetic value
   marker_set->est = realloc ( marker_set->est, 
			       marker_set->NoInd*sizeof(double) );
   if (NULL == marker_set->est) ERR_M;


   int by_matrix_multiplication = 1;

   // Check whether the same marker set is used 
   if ( marker_set->NoMar != effect_set->NoMar)
     {
       by_matrix_multiplication = 0;                // Not the same no. of markers
     }
   else {
     for ( m = 0; m < marker_set->NoMar; m++) {
       if ( 0 != strcmp(marker_set->MarName[m], effect_set->MarName[m]) ) {
	   by_matrix_multiplication = 0;          // Not the same markers
	   break;
	 }
     }
   }

   // if ( marker_set != effect_set) by_matrix_multiplication = 0;

   if ( 1 == by_matrix_multiplication )
     {
       gs_info(effect_set,1,"Marker sets are identical");
       
       // Make design matrix if necessary
       if ( NULL == marker_set->Z) {

	 gs_info(effect_set,1,"Constructing design matrix");

	 marker_set->Z =  realloc( marker_set->Z, marker_set->NoInd * 
				   marker_set->NoMar * sizeof(double) );
	 int a, b, lar_all, sml_all;
	 for ( m = 0; m < marker_set->NoMar; m++)
	   {
	     // Determine the larges allele for the marker 
	     lar_all = effect_set->EffAll[m+1];
	     sml_all = effect_set->CmpAll[m+1] ;
	     
	     // Construct design matrix
	     for ( i = 0; i < marker_set->NoInd ; i++)
	       {
		 idx_all = m*marker_set->NoInd + i ;
		 
		 a = marker_set->All1[idx_all] ;
		 b = marker_set->All2[idx_all] ;
		 
		 idx_Z   = i*marker_set->NoMar + m ;

		 if ( (a==lar_all) && (b==lar_all) ) { 
		   marker_set->Z[idx_Z] = marker_set->AA; 
		 }
		 else	{
		   if ( (a==sml_all) && (b==sml_all) )  {
		     marker_set->Z[idx_Z] = marker_set->aa;
		   }
		   else {
		     if ( ((a==lar_all) && (b==sml_all)) ||
			  ((a==sml_all) && (b==lar_all))   ) {
		       marker_set->Z[idx_Z] = marker_set->aA ;
		     }
		     else {
		       if ( (a==lar_all) || (b==lar_all) ) {
			 marker_set->Z[idx_Z] = marker_set->nA ;
		       }
		       else {
			 if ( (a==sml_all) || (b==sml_all) ) {
			   marker_set->Z[idx_Z] = marker_set->an ;
			 }
			 else {
			   marker_set->Z[idx_Z] = marker_set->nn ;
			 }
		       }
		     }
		   }
		 }
	       }
	   }  // Build design matrix
       } // two different data sets?

       st_dgemm ("N","N",
                  marker_set->NoInd, marker_set->NoMar,
		  marker_set->NoMar,1,
		  marker_set->Z,
		  &(effect_set->u[1]),
		  marker_set->est);		 

       for ( i = 0; i < marker_set->NoInd; i++) 
	 marker_set->est[i] += effect_set->u[0];

     }
   else  // Not by matrix multiplication
     {
       gs_info(effect_set,1,"Different marker sets");

#ifdef MT
  #pragma omp parallel 
#endif
        {
	  gs_UNSIGNED i,j,m,idx_all;
	  int a, b, lar_all, sml_all;
	  double count;
	  
#ifdef MT
          #pragma omp for
#endif
	  for ( i = 0; i < marker_set->NoInd; i++) 
	    {
	      marker_set->est[i] = effect_set->u[0];
	      for (j = 1; j < effect_set->dimu; j++) 
		{
		  int found = 0;
		  for ( m = 0; m < marker_set->NoMar; m++) 
		    {
		      if ( 0 == strcmp (marker_set->MarName[m],
					effect_set->EffNme[j]) )
			{
			  lar_all = effect_set->EffAll[j];
			  sml_all = effect_set->CmpAll[j] ;

			  idx_all = m*marker_set->NoInd + i ;
			  a = marker_set->All1[idx_all] ;
			  b = marker_set->All2[idx_all] ;

			  if ( (a==lar_all) && (b==lar_all) ) { 
			    count = effect_set->AA; 
			  }
			  else	{
			    if ( (a==sml_all) && (b==sml_all) )  {
			      count = effect_set->aa;
			    }
			    else {
			      if ( ((a==lar_all) && (b==sml_all)) ||
				   ((a==sml_all) && (b==lar_all))   ) {
				count = effect_set->aA ;
			      }
			      else {
				if ( (a==lar_all) || (b==lar_all) ) {
				  count = effect_set->nA ;
				}
				else {
				  if ( (a==sml_all) || (b==sml_all) ) {
				    count = effect_set->an ;
				  }
				  else {
				    count = effect_set->nn ;
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
		  if (0 == found) {
		    count = effect_set->nn;
		    marker_set->est[i] += count * effect_set->u[j];
		  }
		}
	    }
	}  // parallel
     } // else by matrix multiplication
     
 /* Write estimated values to output file */
 if (1 == *auxfiles)
   {
     gs_UNSIGNED i;
     FILE *fp;
     if ( (fp = fopen(*out_filename,"w")) == NULL ) ERR_F;
     if (NULL == marker_set->y )
       {
	 fprintf(fp,"yhat \n");
	 for ( i = 0; i < marker_set->NoInd ; i++)
	   fprintf(fp,"%s %f\n",
		   marker_set->IndName[i],
		   marker_set->est[i]);
       }
     else
       {
	 fprintf(fp,"y yhat \n");
	 for ( i = 0; i < marker_set->NoInd ; i++)
	   fprintf(fp,"%s %f %f\n",
		   marker_set->IndName[i],
		   marker_set->y[i],
		   marker_set->est[i]);
       }

     fclose(fp);  
   }

 *retval = 0; return;
}

void gs_estimate_bv_01_GV ( char** effect_set ,
			    char** marker_set ,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval          )
{
  gs_estimate_bv_01 ( gs_fdta(*effect_set),
		      gs_fdta(*marker_set),
		      out_filename,
		      auxfiles, 
		      retval          ) ;

}

/* Estimate breeding values for the loaded marker set */
void gs_check_model_fit_01_GV ( char**  out_filename ,
				int*    auxfiles     , 
				int*    retval       )
{
  gs_estimate_bv_01 ( gv,
		      gv,
		      out_filename,
		      auxfiles, 
		      retval              ) ;

}

/* Estimate breeding values for an new marker set */
void gs_estimate_breeding_values_01_GV ( char**  in_marker,
					 char**  out_filename,
					 int*    auxfiles, 
					 int*    retval    )
{
  int Srv; int* Sretval = &Srv;
 
  gs_varset_type VS; gs_varset_type *valset = &VS;
  gs_init (valset);

  gs_read_marker_data_01 (valset,in_marker,nixp,nixp,nixp,neinp,Sretval);
  if (-2 == *Sretval) {gs_reset(valset); *retval = -2; return; }
                                 
  gs_estimate_bv_01 (gv,valset,out_filename,auxfiles,Sretval);
  if (-2 == *Sretval) {gs_reset(valset); *retval = -2; return; }

  gs_reset(valset);

  *retval = 0; return;

}

void gs_validation_01_GV ( char**  in_marker,
			   char**  in_performance,
			   char**  out_filename,
			   int*    auxfiles, 
			   int*    retval     )
{

  int Srv; int* Sretval = &Srv;
 
  gs_varset_type VS; gs_varset_type *valset = &VS;
  gs_init (valset);

  gs_read_marker_data_01 (valset,in_marker,nixp,nixp,nixp,neinp,Sretval);
  if (-2 == *Sretval) {gs_reset(valset); *retval = -2; return; }

  gs_read_performance_data_01 (valset,in_performance,nixp,neinp,Sretval);
  if (-2 == *Sretval) {gs_reset(valset); *retval = -2; return; }

  gs_estimate_bv_01 (gv,valset,nixp,neinp,Sretval);
  if (-2 == *Sretval) {gs_reset(valset); *retval = -2; return; }


 /* Write estimated values to output file */
  if (1 == *auxfiles)
    {
      FILE *fp;
      if ( NULL == (fp = fopen(*out_filename,"w"))  ) 
	{gs_reset(valset);ERR_F;}
      fprintf(fp,"y yhat \n");
      for (gs_UNSIGNED i = 0; i < valset->NoInd ; i++)
	fprintf(fp,"%s %f %f\n",
		valset->IndName[i],
		valset->y[i],
		valset->est[i]);
      fclose(fp);  
    }

  gs_reset (valset);
  
  *retval = 0; return;
}

void gs_build_esvs_01 ( gs_varset_type *gv ,
			gs_varset_type *es ,
			gs_varset_type *vs ,
			int*    n_estimation ,
			char**  out_es_m,
			char**  out_es_p,
			char**  out_vs_m,
			char**  out_vs_p,
			int*    auxfiles, 
			int*    retval          )
{
  
  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
     gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
     *retval = -2; return ;
   }

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
     *retval = -2; return ;
   }

  gs_reset(es);
  gs_reset(vs);

   gs_UNSIGNED CVSiEst  = (gs_UNSIGNED)*n_estimation;

   gs_UNSIGNED i,j,k,m, idx_g, idx_e, idx_v, idx;

   /* Estimation set */
   gs_UNSIGNED *CVes = malloc ( gv->NoInd * sizeof(gs_UNSIGNED) );
   if (NULL == CVes) ERR_M;
 
   /* Random permutation*/

   int *p = malloc( gv->NoInd * sizeof(int) );
   if (NULL == p) {free(CVes);ERR_M;}

   /* Random permutaion of the indiviudals */
   /* Knuth Stanford Graph Base */
   p[0] = 0;
   for (int i = 1; i < (int)gv->NoInd; ++i) {
     int j = st_rand_index(i + 1); // New implementation using MathLib12
     p[i] = p[j];
     p[j] = i;
   }
   /* Estimation set 1, validation set 0 */
   /* Take the first random individuals as estimation set*/
   for ( i = 0; i < gv->NoInd; i++ ) CVes[i]      = 0;
   for ( i = 0; i < CVSiEst; i++ )   CVes[ p[i] ] = 1;
   free(p);

   /* Copy resampled data */

   /* Individuals */
   es->NoInd = CVSiEst;
   vs->NoInd = gv->NoInd - es->NoInd;

   es->IndName = realloc(es->IndName, es->NoInd * sizeof(char[gs_MAXNME]));
   vs->IndName = realloc(vs->IndName, vs->NoInd * sizeof(char[gs_MAXNME]));
   if ( (NULL==es->IndName) || (NULL==vs->IndName) ) {free(CVes);ERR_M;}

   for ( i=0,j=0,k=0; i < gv->NoInd; i++)
     {
       if ( 1 == CVes[i] )
	 {
	   strcpy( es->IndName[j] , gv->IndName[i] );
	   j++;
	 }
       else
	 {
	   strcpy( vs->IndName[k] , gv->IndName[i] );
	   k++;
	 }
     }

   /* Markers*/
   es->NoMar = gv->NoMar;
   vs->NoMar = gv->NoMar;

   es->MarName = realloc ( es->MarName, es->NoMar * sizeof(char[gs_MAXNME]) );
   vs->MarName = realloc ( vs->MarName, vs->NoMar * sizeof(char[gs_MAXNME]) );
   if ( (NULL==es->MarName) || (NULL==vs->MarName) ) {free(CVes);ERR_M;}
   
   for ( m = 0; m < gv->NoMar; m ++)
     {
	  strcpy( es->MarName[m] , gv->MarName[m] );
	  strcpy( vs->MarName[m] , gv->MarName[m] );
     }

  /* Genotypic data */
  es->All1 =  realloc(es->All1, es->NoMar * es->NoInd * sizeof(*es->All1));
  es->All2 =  realloc(es->All2, es->NoMar * es->NoInd * sizeof(*es->All1));
  if ( (NULL == es->All1)  || (NULL == es->All2) ) {free(CVes);ERR_M;}

  vs->All1 =  realloc(vs->All1, vs->NoMar * vs->NoInd * sizeof(*vs->All1));
  vs->All2 =  realloc(vs->All2, vs->NoMar * vs->NoInd * sizeof(*vs->All1));
  if ( (NULL == vs->All1)  || (NULL == vs->All2) ) {free(CVes);ERR_M;}

  for ( m = 0; m < gv->NoMar; m++ ) 
    {
      for ( i=0,j=0,k=0; i < gv->NoInd; i++)
	{
	   idx_g = m*gv->NoInd + i;
	   idx_e = m*es->NoInd + j;
	   idx_v = m*vs->NoInd + k;

	   if ( 1 == CVes[i] )
	     {
	       es->All1[idx_e] = gv->All1[idx_g];
	       es->All2[idx_e] = gv->All2[idx_g];
	       j++;
	     }
	   else
	     {
	       vs->All1[idx_v] = gv->All1[idx_g];
	       vs->All2[idx_v] = gv->All2[idx_g];
	       k++;
	     }
 	 }
     }

   /* Note: further data as generated by read_marker_data may be required*/

   /* Performance data */
   es->y = realloc( es->y, es->NoInd*sizeof(double) );
   vs->y = realloc( vs->y, vs->NoInd*sizeof(double) );
   if ( ( NULL== es->y ) || ( NULL== vs->y ) ) {free(CVes);ERR_M;}

   for ( i=0,j=0,k=0; i < gv->NoInd; i++)
     {
       if ( 1 == CVes[i] )
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

  /* Write output files */
  if (1 == *auxfiles){

    FILE * fp1 = fopen(*out_es_m,"w") ;
    FILE * fp2 = fopen(*out_es_p,"w") ;
    FILE * fp3 = fopen(*out_vs_m,"w") ;
    FILE * fp4 = fopen(*out_vs_p,"w") ;
    if ( (NULL==fp1) || (NULL==fp2) || (NULL==fp3) || (NULL==fp4)) ERR_F;

      fprintf(fp1, "I;M;A\n");
      for ( i=0; i < es->NoInd; i++  )
	for ( m=0; m < es->NoMar; m++ ) 
	  {
	    idx = m*es->NoInd + i;
	    if  (GS_ALLELE_MISSING != es->All1[idx]) 
	      {
		fprintf(fp1, "%s;%s;%i\n",
			es->IndName[i], es->MarName[m], es->All1[idx] );
		if (es->All1[idx] != es->All2[idx])
		  fprintf(fp1, "%s;%s;%i\n",
			  es->IndName[i], es->MarName[m], es->All2[idx] );
	      }
	  }

      fprintf(fp2, "I;P\n");
      for ( i=0; i < es->NoInd; i++  )
	fprintf(fp2, "%s;%f\n", es->IndName[i], es->y[i] );

      fprintf(fp3, "I;M;A\n");
      for ( i=0; i < vs->NoInd; i++  )
	for ( m=0; m < vs->NoMar; m++ ) 
	  {
	    idx = m*vs->NoInd + i;
	    if  (GS_ALLELE_MISSING != vs->All1[idx]) 
	      {
		fprintf(fp3, "%s;%s;%i\n",
			vs->IndName[i], vs->MarName[m], vs->All1[idx] );
		if (vs->All1[idx] != vs->All2[idx])
		  fprintf(fp3, "%s;%s;%i\n",
			  vs->IndName[i], vs->MarName[m], vs->All2[idx] );
	      }
	  }

      fprintf(fp4, "I;P\n");
      for ( i=0; i < vs->NoInd; i++  )
	fprintf(fp4, "%s;%f\n", vs->IndName[i], vs->y[i] );

      fclose(fp1);  
      fclose(fp2);  
      fclose(fp3);  
      fclose(fp4);  

  }
  
  *retval = -0; return ;

}

void gs_build_esvs_01_GV ( int*    n_estimation ,
			   char**  base_set,
			   char**  estimation_set,
			   char**  validation_set,
			   char**  out_es_m,
			   char**  out_es_p,
			   char**  out_vs_m,
			   char**  out_vs_p,
			   int*    auxfiles, 
			   int*    retval          )
{

  /* Achtung instabil bei gleichen Ausgabefilenamen */

  /* First use of global validation and estimation set */
  if (0==esi) {gs_init(es);esi=1;}
  if (0==vsi) {gs_init(vs);vsi=1;}

  gs_varset_type* ees = es;
  gs_varset_type* vvs = vs;

  if ( 0 != strcmp ( "none" , *estimation_set ) )
    ees = gs_fdta(*estimation_set);

  if ( 0 != strcmp ( "none" , *validation_set ) )
    vvs = gs_fdta(*validation_set);

  gs_build_esvs_01 ( gs_fdta(*base_set) ,
		  ees ,
		  vvs ,
		  n_estimation ,
		  out_es_m,
		  out_es_p,
		  out_vs_m,
		  out_vs_p,
		  auxfiles, 
		  retval          );
}

void gs_return_effects_01 (  gs_varset_type *gv ,
			     char**  out_filename,
			     int*    retval         )
{
  gs_UNSIGNED i;
  FILE *outfp;
   
  if ( NULL == gv->u ) {
      gs_info(gv,-2, "No effects estimated");
      *retval = -2; return ;
  }

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;
  
  fprintf(outfp,"marker allele cmpallele effect\n");
  for ( i = 0; i < gv->dimu ; i++)
    fprintf(outfp,"%s %i %i %20.16f \n",
	    gv->EffNme[i],
	    gv->EffAll[i],
	    gv->CmpAll[i],
	    gv->u [i]);
  
  fclose(outfp);  

  *retval = 0; 
  return ;
}

void gs_return_effects_01_GV ( char**  out_filename,
			       int*    retval,
			       char** set_name         )
{
  gs_return_effects_01 ( gs_fdta(*set_name),
			 out_filename,
			 retval         ) ;
}


void gs_set_effects_01 (  gs_varset_type *gv ,
		          char**  filename,
		          int*    retval         )
{

  gs_UNSIGNED i,j,k;

  /* Check wheter marker data is loaded */
  if ( 0 == gv->NoMar ) {
      gs_info(gv,-2, "No marker data available");
      *retval = -2; return ;
  }

  /* Remove old design matrix and estimated genotypes */
  GS_FREE(gv->Z);
  GS_FREE(gv->est);

  /* Allocate memory for the effect descriptions */
  gv->NoEff =  1 + gv->NoMar;
  gv->EffNme = realloc (gv->EffNme, gv->NoEff * sizeof(char[gs_MAXNME]) );
  gv->EffAll = realloc (gv->EffAll, gv->NoEff * sizeof(int) );
  gv->CmpAll = realloc (gv->CmpAll, gv->NoEff * sizeof(int) );
  if ( (NULL==gv->EffNme) || (NULL==gv->EffAll) || (NULL==gv->CmpAll)) ERR_M;

  /* Names of effects */
  strcpy(gv->EffNme[0],"mu");
  for ( i = 1; i < gv->NoEff; i++)
    strcpy(gv->EffNme[i],gv->MarName[i-1]);		       

  /* Allocate memory for the solution vector */
  gv->dimu = 1 + gv->NoMar;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );
  if ( NULL == gv->u ) ERR_M;

  /* Initialize effects with zeros */
  for ( i = 0; i < gv->dimu ; i++) {
    gv->EffAll[i] = 0;
    gv->u [i]     = 0;

  }

  {
    gs_UNSIGNED NoLines = 0; /* Number of lines in the input file       */

    FILE *fp;

    char*  line = NULL;   char** linep = &line;
    size_t n    = 0;      size_t* np   = &n;
    
    char *p ;

    if ( NULL == (fp = fopen(*filename,"r")) ) ERR_F;

    /* Headline */
    st_getline ( linep , np , fp ) ;

    /* Count lines */
    while (!feof(fp))
      {
	if ( 0 < st_getline ( linep , np , fp ) )
	  { 
            // Read in a data line

	    p = strtok (line,gs_tok_del);
	    if ( NULL == p ) 
	      {gs_info(gv,-2,"Error in map file");*retval=-2;return;}
	    if ( strlen(p) > gs_MAXNME)
	      {gs_info(gv,-2,"Marker name to long");*retval=-2;return;}
	    
	    p = strtok (NULL,gs_tok_del);
	    if  ( NULL == p )  
	      {gs_info(gv,-2,"Error in map file");*retval=-2;return;}
	    
	    p = strtok (NULL,gs_tok_del); 
	    if ( NULL == p ) {gs_info(gv,-2,"Error in map file");*retval=-2;return;}

	    p = strtok (NULL,gs_tok_del); 
	    if ( NULL == p ) {gs_info(gv,-2,"Error in map file");*retval=-2;return;}

	    strcpy(line,"\0");
	    NoLines++; 
	  }
      }
    fclose(fp);

    /* Allocate memory */
    char   (*EffNme) [1+gs_MAXNME] = NULL; /* Effect name */
    int     *EffAll                = NULL;  /* Effect allele */     
    int     *CmpAll                = NULL;  /* Complementary allele */     
    double  *Eff                   = NULL;  /* Effect size */ 
  
    EffNme = realloc ( EffNme, NoLines*sizeof(char[1+gs_MAXNME]) );
    EffAll = realloc ( EffAll, NoLines*sizeof(int) );
    CmpAll = realloc ( CmpAll, NoLines*sizeof(int) );
    Eff    = realloc ( Eff ,   NoLines*sizeof(double) );
    if ( (NULL == EffNme) || (NULL == EffAll) || (NULL == CmpAll)
	 || (NULL == Eff) ) ERR_M;

    /* Read in data lines */
    if ( NULL == (fp = fopen(*filename,"r")) ) ERR_F;
    st_getline ( linep , np , fp ) ;
    for (j = 0; j < NoLines; j++)
      {
	st_getline ( linep , np , fp ) ;
	gs_remove_newline ( line ) ;

	p = strtok (line,gs_tok_del);
        strcpy(EffNme[j],p);

	p = strtok (NULL,gs_tok_del);
        {
          long allele_value;
          if ((NULL == p) ||
              (0 != gs_parse_long_strict(p,&allele_value)) ||
              (allele_value < 0) || (allele_value > GS_ALLELE_MAX))
            {
              fclose(fp);
              free(line);
              free(EffNme);
              free(EffAll);
              free(CmpAll);
              free(Eff);
              gs_info(gv,-2,"Invalid effect allele code");
              *retval = -2;
              return;
            }
          EffAll[j] = (int)allele_value;
        }

	p = strtok (NULL,gs_tok_del);
        {
          long allele_value;
          if ((NULL == p) ||
              (0 != gs_parse_long_strict(p,&allele_value)) ||
              (allele_value < 0) || (allele_value > GS_ALLELE_MAX))
            {
              fclose(fp);
              free(line);
              free(EffNme);
              free(EffAll);
              free(CmpAll);
              free(Eff);
              gs_info(gv,-2,"Invalid complementary allele code");
              *retval = -2;
              return;
            }
          CmpAll[j] = (int)allele_value;
        }

	p = strtok (NULL,gs_tok_del);
	Eff[j] = (double)atof(p);

      }
    fclose(fp);
    free(line);

    for (j = 0; j < NoLines; j++)
      {
	if ( 0 == strcmp(EffNme[j],gv->EffNme[j]) ) 
	  { /* same order of effect names */
	    gv->EffAll[j] = EffAll[j];
	    gv->CmpAll[j] = CmpAll[j];
	    gv->u[j]      = Eff[j] ;
	    
	  }
	else
	  { /* different order of effect names */
	    for (k = 0; k < (1 + gv->NoMar) ; k++)
	      {
		if ( 0 == strcmp(EffNme[j],gv->EffNme[k]) ) {
		  gv->EffAll[k] = EffAll[j];
		  gv->CmpAll[k] = CmpAll[j];
		  gv->u[k]      = Eff[j] ;
		  break;
		}
	      }
	  }
      }

    free(EffNme);
    free(EffAll);
    free(CmpAll);
    free(Eff   );

  }

  *retval = 0; 
  return ;
}

void gs_set_effects_01_GV ( char**  filename,
		            int*    retval,
		            char** set_name         )
{
  gs_set_effects_01 ( gs_fdta(*set_name),
		      filename,
		      retval         ) ;
}




void gs_return_pvals_01 (  gs_varset_type *gv ,
			   char**  out_filename,
			   int*    retval         )
{
  gs_UNSIGNED i;
  FILE *outfp;
   
  if ( NULL == gv->p ) {
      gs_info(gv,-2, "No p-values estimated");
      *retval = -2; return ;
  }

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;
  
  fprintf(outfp," effect pvalue \n");
  for ( i = 0; i < gv->dimu ; i++)
    fprintf(outfp,"%s.%i %f %f \n",
	    gv->EffNme[i],
	    gv->EffAll[i],
	    gv->u [i],
	    gv->p [i]);
  
  fclose(outfp);  

  *retval = 0; 
  return ;
}

void gs_return_pvals_01_GV ( char**  out_filename,
			     int*    retval,
			     char** set_name         )
{
  gs_return_pvals_01 ( gs_fdta(*set_name),
			 out_filename,
			 retval         ) ;
}

void gs_check_pvals_01 (  gs_varset_type *gv ,
			  int*    retval         )
{
  if ( NULL == gv->p ) *retval = 0; 
  else *retval = 1;
  return;
}

void gs_check_pvals_01_GV (  int*    retval,
			     char** set_name         )
{
  gs_check_pvals_01 ( gs_fdta(*set_name),
		      retval         ) ;
}


void gs_set_allele_codes_01 (  gs_varset_type *gv ,
			       double* aa         ,
			       double* aA         ,
			       double* AA         ,
			       double* an         ,
			       double* nn         ,
			       double* nA           )
{
  gv->aa         = * aa ; 
  gv->aA         = * aA ;
  gv->AA         = * AA ;
  gv->an         = * an ; 
  gv->nn         = * nn ;
  gv->nA         = * nA ;

  /* Remove old design matrix and estimated genotypes */
  GS_FREE(gv->Z);
  GS_FREE(gv->est);
  
}

void gs_set_allele_codes_01_GV (  double* aa         ,
				  double* aA         ,
				  double* AA         ,
				  double* an         ,
				  double* nn         ,
				  double* nA         ,
				  char** set_name  )
{

  gs_set_allele_codes_01 ( gs_fdta(*set_name) ,
			   aa                 ,
			   aA                 ,
			   AA                 ,
			   an                 ,
			   nn                 ,
			   nA                  ) ;
}




/* Nov 2011 */
/* Functions for plotting graphical genotypes */

/* qsort pointers to mappoints comparison function*/
int gs_mappoint_cmp(const void *a, const void *b)
{
  const gs_mappoint_type **ia =
    (const gs_mappoint_type **)a;
  const gs_mappoint_type **ib =
    (const gs_mappoint_type **)b;
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


void gs_reset_all_but_marker_matrix (gs_varset_type *gv)
{
  gs_UNSIGNED NoInd = gv->NoInd; 
  gs_UNSIGNED NoMar = gv->NoMar; 

  char (*IndName) [gs_MAXNME] = gv->IndName; gv->IndName = NULL;  
  char (*MarName) [gs_MAXNME] = gv->MarName; gv->MarName = NULL; 

  GS_ALLELE_TYPE * All1  = gv->All1; gv->All1 = NULL;
  GS_ALLELE_TYPE * All2  = gv->All2; gv->All2 = NULL;

  gs_reset(gv);

  gv->NoInd = NoInd; 
  gv->NoMar = NoMar; 

  gv->IndName = IndName;  
  gv->MarName = MarName; 

  gv->All1 = All1;
  gv->All2 = All2;

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

void gs_chrom_stats_01 (gs_varset_type *gv ,
			char** out_filename,
			int* auxfiles,
			int* retval          )
{
  gs_UNSIGNED c, i;

  if ( 0 == gv->NLoc ) {

    snprintf(gv->msg, sizeof(gv->msg),
         "No linkage map has been loaded");
    gs_info(gv,-2,gv->msg);
      *retval = -2;   return;
    }

  gs_reset_chrom_stats(gv);

  /* Count chromosomes and memory allocation*/

  gv->NChrom = 1;
  for ( i = 1 ; i < gv->NLoc; i++ )
    if  (gv->map[i-1]->chrom != gv->map[i]->chrom) gv->NChrom++;

  gv->chrom = realloc( gv->chrom, gv->NChrom * sizeof (gs_chrom_type) );

  for (c = 0; c < gv->NChrom; c ++ ) 
    {
      gv->chrom[c].rf = NULL;
      gv->chrom[c].q  = NULL;
      gv->chrom[c].ld = NULL;
    }

  snprintf(gv->msg, sizeof(gv->msg),
         "No. of chromosomes: %lu",
         gv->NChrom);

  gs_info ( gv, 1, gv->msg );

  /* First locus of each chromosome and number of loci  */

  gv->chrom[0].NLocC = 1;
  gv->chrom[0].first_loc   = &(gv->map[0]);
  gv->chrom[0].first_loc_i = 0;
  for ( c=0, i=1 ; i < gv->NLoc; i++ ) 
    {
      if ( gv->map[i-1]->chrom == gv->map[i]->chrom )
	{
	  gv->chrom[c].NLocC ++ ;
	}
      else
	{ 
	  c ++ ;
	  gv->chrom[c].first_loc   = &(gv->map[i]);
	  gv->chrom[c].first_loc_i = i;
	  gv->chrom[c].NLocC = 1;
	}
    }


  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      gs_mappoint_type ** p;
 
      for ( c=0; c < gv->NChrom; c++ )
	{
	  p =  gv->chrom[c].first_loc ;
	  snprintf(gv->msg, sizeof(gv->msg),
		   "%lu %lu %f \n",
		   1+c,
		   gv->chrom[c].NLocC,
		   (p[gv->chrom[c].NLocC-1])->pos);
 	  fprintf(outfp,"%s",gv->msg);
	}
      
      fclose(outfp);
    }

  *retval=0; return;
}

void gs_chrom_stats_01_GV ( char**  out_filename,
			    int*    auxfiles,
			    int*    retval,
			    char**  set_name)
{
  gs_chrom_stats_01 ( gs_fdta(*set_name),
		      out_filename,
		      auxfiles,
		      retval         ) ;
}

gs_FLOAT gs_rf_haldane( gs_FLOAT d)
{
  return  ( (1 - exp(-2*d) ) /2  );
}

gs_FLOAT gs_rf_kosambi( gs_FLOAT d)
{
  return ( exp(4*d) - 1 ) / ( exp(4*d) + 1 ) / 2;
}

void gs_calc_rf_01 ( gs_varset_type *gv  ,
		     char** map_function ,
		     char** out_filename ,
		     int* auxfiles       ,
		     int* retval           )
{
  gs_UNSIGNED c, i, j, idx;
  gs_FLOAT d;
  gs_mappoint_type ** p;

  if ( 0 == gv->NChrom ) { 
    snprintf(gv->msg, sizeof(gv->msg),
	     "No chromosome statistics available");
    gs_info(gv,-2,gv->msg);
    *retval = -2;   return;
  }

  gs_FLOAT (*gs_rf) (gs_FLOAT d);

  if ( 0 == strcmp(*map_function,"Haldane"))
    {
      gs_rf = gs_rf_haldane;
    }

  else if ( 0 == strcmp(*map_function,"Kosambi"))
    {
      gs_rf = gs_rf_kosambi;
    }

  else 
    {
      gs_info(gv,-2,"Unknown map function");
      *retval = -2; return;
    }

  /* Allocate memory for recombination frequencies*/
  for (c = 0; c < gv->NChrom; c ++ ) 
    {
      gv->chrom[c].rf = realloc (gv->chrom[c].rf, 
	gv->chrom[c].NLocC * gv->chrom[c].NLocC * sizeof(gs_FLOAT) );
      gv->chrom[c].q  = realloc (gv->chrom[c].q , 
	gv->chrom[c].NLocC * gv->chrom[c].NLocC * sizeof(gs_FLOAT) ); 
    }

  for ( c=0; c < gv->NChrom; c++ )
    {
       p =  gv->chrom[c].first_loc ;
       for  ( i = 0 ; i < gv->chrom[c].NLocC ; i ++)
	 for  ( j = 0 ; j < gv->chrom[c].NLocC ; j ++)
	   {
	     d = fabs ( (p[i])->pos - (p[j])->pos );

	     idx = i * gv->chrom[c].NLocC + j;
	     (gv->chrom[c].rf) [idx] = gs_rf (d/100) ;
	   }
    }

  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( c=0; c < gv->NChrom; c++ )
	{
	  for  ( i = 0 ; i < ( gv->chrom[c].NLocC - 1 ) ; i ++)
	    for  ( j = i+1 ; j < gv->chrom[c].NLocC ; j ++)
	      {
		idx = i * gv->chrom[c].NLocC + j;

		snprintf(gv->msg, sizeof(gv->msg),
			 "%lu %lu %lu %f \n",
			 1+c,
			 1+i,
			 1+j,
			 (double)(gv->chrom[c].rf)[idx]);
		fprintf(outfp,"%s",gv->msg);

	      }
	  
	}

      fclose(outfp);
    }

  *retval=0; return;

 }


void gs_calc_rf_01_GV (	char** map_function ,
			char** out_filename ,
			int*   auxfiles     ,
			int*   retval       ,
			char** set_name       )
{
  gs_calc_rf_01 ( gs_fdta(*set_name),
		  map_function      ,
		  out_filename      ,
		  auxfiles          ,
		  retval              ) ;
}

gs_FLOAT gs_q_dh (int* t, gs_FLOAT r)
{
  return 0.5 + ( 1 - 2*r) / 2 * pow ( (double)(1-r), (double)*t);
}

gs_FLOAT gs_q_ssd (int* t, gs_FLOAT r)
{
/*
  old version: check
  return 0.5 + ( 1.0 - 2.0*r) / ( 2.0 - 4.0*r ) *
              pow ( (double)(1.0-r), 1.0 + (double)(*t) );
*/
  return 0.5 + ( 1.0 - 2.0*r) / ( 2.0 - 4.0*r ) *
              pow ( (double)(1.0-r), 4.0 + (double)(*t) );
}

gs_FLOAT gs_q_crs (int* t, gs_FLOAT r)
{
  int tt = *t;   // to aviod compiler warnings
  tt = 1;
  return  tt*(1.0 - r) ;
}

void gs_calc_q_01 ( gs_varset_type *gv   ,
		    char**  pop_type     ,
		    int*    t            ,
		    char** out_filename  ,
		    int* auxfiles        ,
		    int* retval            )
{
  gs_UNSIGNED c, i, j, idx;

  if  ( ( 0 == gv->NChrom ) || ( NULL == gv->chrom[0].rf ) ) { 
    snprintf(gv->msg, sizeof(gv->msg),
         "No recombination frequencies available");
    gs_info(gv,-2,gv->msg);
    *retval = -2;   return;
    }

  gs_FLOAT (*gs_q) (int* t, gs_FLOAT r);

  if ( 0 == strcmp(*pop_type,"DH"))
    {
      gs_q = gs_q_dh;
    }

  else if ( 0 == strcmp(*pop_type,"SSD"))
    {
      gs_q = gs_q_ssd;
    }

  else if ( 0 == strcmp(*pop_type,"CRS"))
    {
      gs_q = gs_q_crs;
    }
  
  else 
    {
      gs_info(gv,-2,"Unknown population type");
      *retval = -2; return;
    }

  
  for ( c=0; c < gv->NChrom; c++ )
    {
      for  ( i = 0 ; i < gv->chrom[c].NLocC ; i ++)
	for  ( j = 0 ; j < gv->chrom[c].NLocC ; j ++)
	  {
	    idx = i * gv->chrom[c].NLocC + j;
	    (gv->chrom[c].q)[idx] =  gs_q (t,(gv->chrom[c].rf)[idx]);
	  }
    }
  
  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( c=0; c < gv->NChrom; c++ )
	{
	  for  ( i = 0 ; i < ( gv->chrom[c].NLocC - 1 ) ; i ++)
	    for  ( j = i+1 ; j < gv->chrom[c].NLocC ; j ++)
	      {
		idx = i * gv->chrom[c].NLocC + j;

		snprintf(gv->msg, sizeof(gv->msg),
			 "%lu %lu %lu %f \n",
			 1+c,
			 1+i,
			 1+j,
			 (double)(gv->chrom[c].q)[idx]);
		fprintf(outfp,"%s",gv->msg);

	      }
	  
	}

      fclose(outfp);
    }

  *retval=0; return;

 }


void gs_calc_q_01_GV ( char**  pop_type     ,
		       int*    t            ,
		       char**  out_filename ,
		       int*    auxfiles     ,
		       int*    retval       ,
		       char**  set_name       )
{
  gs_calc_q_01 ( gs_fdta(*set_name),
		 pop_type,
		 t,
		 out_filename,
		 auxfiles,
		 retval         ) ;
}


void gs_replace_comma (char *p)
{
  char* q = strchr(p,',');
  if (NULL != q) *q = '.';
}

void gs_read_map_01 ( gs_varset_type *gv ,
                      char**  map_filename,
                      int*    skip,
                      char**  fmt,
                      double* map_stretch,
                      char**  out_filename,
                      int*    auxfiles,
                      int*    retval        )
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

  if (!((0 == strcmp(*fmt,"cpms")) ||
        (0 == strcmp(*fmt,"mcp"))))
    {
      gs_info(gv,-2,"Implemented map formats: cpms, mcp");
      *retval = -2;
      return;
    }

  if (*skip < 0)
    {
      gs_info(gv,-2,"Invalid number of map lines to skip");
      *retval = -2;
      return;
    }

  if (!isfinite(*map_stretch))
    {
      gs_info(gv,-2,"Invalid map stretch");
      *retval = -2;
      return;
    }

  if ((0 == gv->NoInd) || (0 == gv->NoMar) ||
      (NULL == gv->IndName) || (NULL == gv->MarName) ||
      (NULL == gv->All1) || (NULL == gv->All2))
    {
      gs_info(gv,-2,"No marker data available");
      *retval = -2;
      return;
    }

  gs_tmp_init_from(&tmp, gv);

  fp = fopen(*map_filename, "r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      *retval = -2;
      return;
    }

  for (int s = 0; s < *skip; s++)
    {
      rc = st_getline(&line, &n, fp);
      if (rc < 0)
        {
          gs_info(gv,-2,"Memory allocation error");
          failed = 1;
          goto cleanup;
        }
      if (0 == rc)
        {
          gs_info(gv,-2,"Map file shorter than skip value");
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

      if (0 == strcmp(*fmt,"cpms"))
        {
          p = strtok(line, gs_tok_del);
          if ((NULL == p) ||
              (0 != gs_parse_long_strict(p, &chrom_l)) ||
              (0 == chrom_l) ||
              (chrom_l < INT_MIN) || (chrom_l > INT_MAX))
            {
              gs_info(gv,-2,"Error in map file chromosome");
              failed = 1;
              goto cleanup;
            }

          p = strtok(NULL, gs_tok_del);
          if ((NULL == p) || (0 != gs_parse_double_strict(p, &pos)))
            {
              gs_info(gv,-2,"Error in map file position");
              failed = 1;
              goto cleanup;
            }

          name = strtok(NULL, gs_tok_del);
          class_name = strtok(NULL, gs_tok_del);
          if ((NULL == name) || (NULL == class_name))
            {
              gs_info(gv,-2,"Error in map file");
              failed = 1;
              goto cleanup;
            }
        }
      else
        {
          name = strtok(line, gs_tok_del);
          p = strtok(NULL, gs_tok_del);
          if ((NULL == name) || (NULL == p) ||
              (0 != gs_parse_long_strict(p, &chrom_l)) ||
              (0 == chrom_l) ||
              (chrom_l < INT_MIN) || (chrom_l > INT_MAX))
            {
              gs_info(gv,-2,"Error in map file chromosome");
              failed = 1;
              goto cleanup;
            }

          p = strtok(NULL, gs_tok_del);
          if ((NULL == p) || (0 != gs_parse_double_strict(p, &pos)))
            {
              gs_info(gv,-2,"Error in map file position");
              failed = 1;
              goto cleanup;
            }

          class_name = "m";
        }

      if (NULL != strtok(NULL, gs_tok_del))
        {
          gs_info(gv,-2,"Too many fields in map file line");
          failed = 1;
          goto cleanup;
        }

      pos *= *map_stretch;
      if (!isfinite(pos))
        {
          gs_info(gv,-2,"Invalid stretched map position");
          failed = 1;
          goto cleanup;
        }

      name_len = strlen(name);
      class_len = strlen(class_name);
      if ((name_len >= gs_MAXSTR) || (class_len >= gs_MAXSTR))
        {
          gs_info(gv,-2,"Description too long");
          failed = 1;
          goto cleanup;
        }
      if (name_len > mar_size) mar_size = (gs_UNSIGNED)name_len;
      if (class_len > cls_size) cls_size = (gs_UNSIGNED)class_len;

      if (nloc == capacity)
        {
          gs_UNSIGNED new_capacity = (0 == capacity) ? 128 : capacity * 2;
          gs_mappoint_build_type *new_map_build;

          if (new_capacity < capacity)
            {
              gs_info(gv,-2,"Map file too large");
              failed = 1;
              goto cleanup;
            }

          new_map_build = gs_realloc_array(map_build,
                                           (size_t)new_capacity,
                                           sizeof(*new_map_build));
          if (NULL == new_map_build)
            {
              gs_info(gv,-2,"Memory allocation error");
              failed = 1;
              goto cleanup;
            }

          map_build = new_map_build;
          capacity = new_capacity;
        }

      map_build[nloc].chrom = (int)chrom_l;
      map_build[nloc].pos = pos;
      map_build[nloc].nr = gs_UNSIGNED_MAX;

      if ((0 != gs_string_builder_append(&map_str_s,
                                         &map_str_n,
                                         &map_str_capacity,
                                         name,
                                         &(map_build[nloc].name_offset))) ||
          (0 != gs_string_builder_append(&map_str_s,
                                         &map_str_n,
                                         &map_str_capacity,
                                         class_name,
                                         &(map_build[nloc].class_offset))))
        {
          gs_info(gv,-2,"Memory allocation error");
          failed = 1;
          goto cleanup;
        }

      nloc++;
    }

  if (rc < 0)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  fclose(fp);
  fp = NULL;

  if (0 == nloc)
    {
      gs_info(gv,-2,"Linkage map contains no loci");
      failed = 1;
      goto cleanup;
    }

  map_s = gs_malloc_array((size_t)nloc,sizeof(*map_s));
  if (NULL == map_s)
    {
      gs_info(gv,-2,"Memory allocation error");
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

  snprintf(gv->msg, sizeof(gv->msg),
           "Description length: Mar %lu, Class %lu",
           mar_size, cls_size);
  gs_info(gv,1,gv->msg);
  snprintf(gv->msg, sizeof(gv->msg),
           "No. of markers in map file: %lu", nloc);
  gs_info(gv,1,gv->msg);

  map = gs_malloc_array((size_t)nloc,sizeof(*map));
  map_names = gs_malloc_array((size_t)nloc,sizeof(*map_names));
  match = gs_malloc_array((size_t)nloc,sizeof(*match));
  if ((NULL == map) || (NULL == map_names) || (NULL == match))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
    map[lne] = &(map_s[lne]);

  qsort(map, nloc, sizeof(*map), gs_mappoint_cmp);

  /* Deterministically separate exact position ties without passing the
     next genuinely distinct map position. */
  for (gs_UNSIGNED first = 0; first < nloc; )
    {
      gs_UNSIGNED last = first + 1;
      double base = map[first]->pos;
      int chrom = map[first]->chrom;

      while ((last < nloc) &&
             (map[last]->chrom == chrom) &&
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
                  gs_info(gv,-2,"Invalid linkage-map position ordering");
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
                  gs_info(gv,-2,
                          "Cannot separate identical linkage-map positions");
                  failed = 1;
                  goto cleanup;
                }

              if (has_next && !(pos < map[last]->pos))
                {
                  gs_info(gv,-2,
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
      gs_info(gv,-2,"Duplicate marker name in linkage map");
      failed = 1;
      goto cleanup;
    }

  if (0 != gs_make_name_index(gv->MarName, gv->NoMar, &marker_names))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  if (gs_name_index_has_duplicates(marker_names, gv->NoMar))
    {
      gs_info(gv,-2,"Duplicate marker name in marker data");
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
      gs_info(gv,-2,"No linkage-map markers match marker data");
      failed = 1;
      goto cleanup;
    }

  tmp.NoInd = gv->NoInd;
  tmp.NoMar = no_matching;
  tmp.NLoc = no_matching;
  tmp.MaxIL = gv->MaxIL;
  tmp.MaxML = gv->MaxML;

  for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
    if (gs_UNSIGNED_MAX != match[lne])
      if ((0 != gs_pool_add_string_size(&(tmp.map_str_N),
                                        map[lne]->name)) ||
          (0 != gs_pool_add_string_size(&(tmp.map_str_N),
                                        map[lne]->class)))
        {
          gs_info(gv,-2,"Linkage map string storage too large");
          failed = 1;
          goto cleanup;
        }

  tmp.IndName = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.IndName));
  tmp.MarName = gs_malloc_array((size_t)tmp.NoMar,sizeof(*tmp.MarName));
  tmp.map_S = gs_malloc_array((size_t)tmp.NLoc,sizeof(*tmp.map_S));
  tmp.map_str_S = gs_malloc_array(tmp.map_str_N,sizeof(*tmp.map_str_S));
  tmp.map = gs_malloc_array((size_t)tmp.NLoc,sizeof(*tmp.map));
  tmp.mdl = gs_malloc_array((size_t)tmp.NLoc,sizeof(*tmp.mdl));

  if ((NULL == tmp.IndName) || (NULL == tmp.MarName) ||
      (NULL == tmp.map_S) || (NULL == tmp.map_str_S) ||
      (NULL == tmp.map) || (NULL == tmp.mdl) ||
      (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd)))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED i = 0; i < tmp.NoInd; i++)
    strcpy(tmp.IndName[i], gv->IndName[i]);

  {
    gs_UNSIGNED m_new = 0;
    size_t map_str_offset = 0;

    for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
      if (gs_UNSIGNED_MAX != match[lne])
        {
          gs_UNSIGNED m_old = match[lne];

          strcpy(tmp.MarName[m_new], gv->MarName[m_old]);
          memcpy(&(tmp.All1[(size_t)m_new * tmp.NoInd]),
                 &(gv->All1[(size_t)m_old * gv->NoInd]),
                 (size_t)tmp.NoInd * sizeof(*tmp.All1));
          memcpy(&(tmp.All2[(size_t)m_new * tmp.NoInd]),
                 &(gv->All2[(size_t)m_old * gv->NoInd]),
                 (size_t)tmp.NoInd * sizeof(*tmp.All2));

          if (0 != gs_copy_map_point_to_pool(&(tmp.map_S[m_new]),
                                              map[lne],
                                              m_new,
                                              tmp.map_str_S,
                                              tmp.map_str_N,
                                              &map_str_offset))
            {
              gs_info(gv,-2,"Invalid linkage map string storage");
              failed = 1;
              goto cleanup;
            }
          tmp.map[m_new] = &(tmp.map_S[m_new]);
          tmp.mdl[m_new] = m_new;
          m_new++;
        }

    if (map_str_offset != tmp.map_str_N)
      {
        gs_info(gv,-2,"Invalid linkage map string storage");
        failed = 1;
        goto cleanup;
      }
  }

  if (NULL != gv->y)
    {
      tmp.y = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.y));
      if (NULL == tmp.y)
        {
          gs_info(gv,-2,"Memory allocation error");
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

    gs_marker_stats_01(&tmp,
                       empty_p, empty_p, empty_p,
                       zero_p, zero_p, zero_p,
                       zero_p,
                       &srv);
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
          gs_info(gv,-2,"File open error");
          failed = 1;
          goto cleanup;
        }

      for (gs_UNSIGNED lne = 0; lne < tmp.NLoc; lne++)
        fprintf(outfp, "%i %15.12f %s %s\n",
                tmp.map[lne]->chrom,
                tmp.map[lne]->pos,
                tmp.map[lne]->name,
                tmp.map[lne]->class);
      fclose(outfp);
    }

  snprintf(gv->msg, sizeof(gv->msg),
           "No. of matching markers: %lu", tmp.NLoc);
  gs_info(gv,1,gv->msg);

  gs_tmp_commit(gv, &tmp);
  *retval = 0;

cleanup:
  if (NULL != fp) fclose(fp);
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


void gs_read_map_01_GV ( char**  map_filename,
			 int*    skip,
			 char**  fmt,
			 double* map_stretch,
			 char**  out_filename,
			 int*    auxfiles,
			 int*    retval,
			 char**  set_name)
{

  gs_read_map_01 ( gs_fdta(*set_name),
		   map_filename,
		   skip,
		   fmt,
		   map_stretch,
		   out_filename,
		   auxfiles,
		   retval         ) ;
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
  static const char *required[9] = {
    "CHROM", "POS", "ID", "REF", "ALT", "individual",
    "allele1", "allele2", "CHROM.levels"
  };
  int seen[9] = { 0,0,0,0,0,0,0,0,0 };
  SEXP names;

  if ((VECSXP != TYPEOF(list)) || (9 != XLENGTH(list)))
    return 0;

  names = getAttrib(list,R_NamesSymbol);
  if ((STRSXP != TYPEOF(names)) || (9 != XLENGTH(names)))
    return 0;

  for (R_xlen_t i = 0; i < 9; i++)
    {
      SEXP s = STRING_ELT(names,i);
      int found = -1;

      if (NA_STRING == s)
        return 0;

      for (int j = 0; j < 9; j++)
        if (0 == strcmp(CHAR(s),required[j]))
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
    if (!(('A' <= *p && *p <= 'Z') ||
          ('a' <= *p && *p <= 'z') ||
          ('0' <= *p && *p <= '9')))
      return 0;

  return 1;
}


static int st_stvcf_cstr_cmp(const void *a, const void *b)
{
  const char * const *aa = (const char * const *)a;
  const char * const *bb = (const char * const *)b;

  return strcmp(*aa,*bb);
}


static int st_stvcf_alleles(const char *ref,
                             const char *alt,
                             gs_UNSIGNED *n_alt)
{
  char *copy = NULL;
  char **tokens = NULL;
  size_t len;
  gs_UNSIGNED n = 1;
  gs_UNSIGNED k = 0;
  int rc = -1;

  if ((NULL == ref) || ('\0' == ref[0]) ||
      (0 == strcmp(ref,".")) || (NULL != strchr(ref,',')) ||
      (NULL == alt) || ('\0' == alt[0]))
    return -1;

  if (0 == strcmp(alt,"."))
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
      if ((0 == strcmp(alt,ref)) || ('\0' == alt[0]))
        return -1;
      *n_alt = 1;
      return 0;
    }

  len = strlen(alt);
  copy = malloc(len + 1);
  tokens = gs_malloc_array((size_t)n,sizeof(*tokens));
  if ((NULL == copy) || (NULL == tokens))
    goto cleanup;

  memcpy(copy,alt,len + 1);
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
    if (('\0' == tokens[i][0]) ||
        (0 == strcmp(tokens[i],".")) ||
        (0 == strcmp(tokens[i],ref)))
      goto cleanup;

  qsort(tokens,n,sizeof(*tokens),st_stvcf_cstr_cmp);
  for (gs_UNSIGNED i = 1; i < n; i++)
    if (0 == strcmp(tokens[i-1],tokens[i]))
      goto cleanup;

  *n_alt = n;
  rc = 0;

 cleanup:
  free(copy);
  free(tokens);
  return rc;
}


static int st_stvcf_get_allele(SEXP x,
                               R_xlen_t idx,
                               gs_UNSIGNED n_alt,
                               int *missing,
                               GS_ALLELE_TYPE *allele)
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

  if ((value < 0) || ((gs_UNSIGNED)value > n_alt) ||
      (value >= GS_ALLELE_MAX))
    return -1;

  *missing = 0;
  *allele = (GS_ALLELE_TYPE)(value + 1);
  return 0;
}


static int st_stvcf_map_cmp(const void *a, const void *b)
{
  const gs_mappoint_type * const *ia =
    (const gs_mappoint_type * const *)a;
  const gs_mappoint_type * const *ib =
    (const gs_mappoint_type * const *)b;

  if ((*ia)->chrom < (*ib)->chrom) return -1;
  if ((*ia)->chrom > (*ib)->chrom) return 1;
  if ((*ia)->pos < (*ib)->pos) return -1;
  if ((*ia)->pos > (*ib)->pos) return 1;
  if ((*ia)->nr < (*ib)->nr) return -1;
  if ((*ia)->nr > (*ib)->nr) return 1;
  return 0;
}


static int st_stvcf_separate_map_ties(gs_varset_type *gv,
                                      gs_mappoint_type **map,
                                      gs_UNSIGNED nloc)
{
  for (gs_UNSIGNED first = 0; first < nloc; )
    {
      gs_UNSIGNED last = first + 1;
      double base = map[first]->pos;
      int chrom = map[first]->chrom;

      while ((last < nloc) &&
             (map[last]->chrom == chrom) &&
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
                  gs_info(gv,-2,"Invalid linkage-map position ordering");
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
                pos = nextafter(prev,INFINITY);

              if (!isfinite(pos) ||
                  (has_next && !(pos < map[last]->pos)))
                {
                  gs_info(gv,-2,
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

  for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
    {
      size_t idx = (size_t)marker * (size_t)gv->NoInd + i;
      GS_ALLELE_TYPE a1 = gv->All1[idx];
      GS_ALLELE_TYPE a2 = gv->All2[idx];
      GS_ALLELE_TYPE values[2] = { a1, a2 };

      if ((GS_ALLELE_MISSING == a1) ||
          (GS_ALLELE_MISSING == a2))
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

  qsort(alleles,n,sizeof(*alleles),gs_allele_cmp);
  *n_alleles = n;
  *has_missing = missing;
  return 0;
}


static int st_markerdata_allele_index(const GS_ALLELE_TYPE *alleles,
                                      gs_UNSIGNED n,
                                      int allele)
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


/* ----------------------------------------------------------------------- */
/* Direct VCF I/O for Beagle interchange                                   */
/* ----------------------------------------------------------------------- */

#define ST_VCF_POSITION_SCALE 1000000.0

static int st_vcf_token_name_valid(const char *name)
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
    memcpy(out,s,n);
  return out;
}


/* Read one complete line from an ordinary VCF file. */
static int st_vcf_getline(FILE *fp, char **line, size_t *capacity)
{
  size_t used = 0u;

  if ((NULL == fp) || (NULL == line) || (NULL == capacity))
    return -1;

  if ((*capacity < 8192u) || (NULL == *line))
    {
      char *tmp = realloc(*line,8192u);
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
          tmp = realloc(*line,next);
          if (NULL == tmp)
            return -1;
          *line = tmp;
          *capacity = next;
        }

      avail_size = *capacity - used;
      avail = (avail_size > (size_t)INT_MAX) ? INT_MAX : (int)avail_size;
      z = fgets(*line + used,avail,fp);
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
      tmp = realloc(*field,n * sizeof(**field));
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

  if ((NULL == format) || ('\0' == format[0]) ||
      (0 == strcmp(format,".")) || (NULL == gt_index))
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


static int st_vcf_parse_allele_index(const char *s, size_t n,
                                     gs_UNSIGNED n_alt,
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

  if ((0 != st_vcf_parse_allele_index(gt,sep,n_alt,&miss1,&v1)) ||
      (0 != st_vcf_parse_allele_index(gt + sep + 1u,len - sep - 1u,
                                      n_alt,&miss2,&v2)))
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
  value = strtoull(text,&end,10);
  if ((0 != errno) || (end == text) || ('\0' != *end) || (0ull == value))
    return -1;

  *pos = value;
  return 0;
}


static int st_vcf_allele_label(gs_UNSIGNED index, char *buffer,
                               size_t capacity)
{
  static const char base[4] = {'A','C','G','T'};
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

  /* Longer alleles begin with REF='A', which is a valid VCF indel form. */
  value = index - 4u;
  do
    {
      if (n >= sizeof(rev))
        return -1;
      rev[n++] = base[value & 3u];
      value >>= 2u;
    }
  while (value > 0u);

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

  n = snprintf(buffer,capacity,"M%lu",(unsigned long)(index + 1u));
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


SEXP st_read_vcf(SEXP filename, SEXP dataset)
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
  gs_name_index_type *index = NULL;
  char **chrom_name = NULL;
  size_t chrom_capacity = 0u;
  gs_UNSIGNED no_chrom = 0u;
  size_t map_string_offset = 0u;
  unsigned long long previous_raw_pos = 0ull;
  int previous_chrom = 0;

  if ((STRSXP != TYPEOF(filename)) || (1 != XLENGTH(filename)) ||
      (NA_STRING == STRING_ELT(filename,0)) ||
      ('\0' == CHAR(STRING_ELT(filename,0))[0]) ||
      (STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
      (NA_STRING == STRING_ELT(dataset,0)) ||
      ('\0' == CHAR(STRING_ELT(dataset,0))[0]))
    return ScalarInteger(-2);

  gv = gs_fdta((char *)CHAR(STRING_ELT(dataset,0)));
  gs_tmp_init_from(&tmp,gv);

  fp = fopen(CHAR(STRING_ELT(filename,0)),"rb");
  if (NULL == fp)
    {
      gs_info(gv,-2,"Cannot open VCF file");
      failed = 1;
      goto cleanup;
    }

  /* First pass: validate the VCF skeleton and size native storage. */
  for (;;)
    {
      int lr = st_vcf_getline(fp,&line,&line_capacity);
      size_t nfield = 0u;

      if (0 == lr)
        break;
      if (lr < 0)
        {
          gs_info(gv,-2,"Error reading VCF file");
          failed = 1;
          goto cleanup;
        }
      if ('\0' == line[0])
        continue;

      if (0 == strncmp(line,"##SelectionToolsPositionScale=",30))
        {
          char *end = NULL;
          double scale = strtod(line + 30,&end);
          if ((end != line + 30) && ('\0' == *end) && isfinite(scale) &&
              (scale > 0.0))
            position_scale = scale;
          continue;
        }

      if (0 == strncmp(line,"##",2))
        continue;

      if ('#' == line[0])
        {
          if (have_header)
            {
              gs_info(gv,-2,"VCF contains more than one column header");
              failed = 1;
              goto cleanup;
            }
          if (0 != st_vcf_split_tabs(line,&field,&field_capacity,&nfield))
            {
              failed = 1;
              goto cleanup;
            }
          if ((nfield < 10u) || (0 != strcmp(field[0],"#CHROM")) ||
              (0 != strcmp(field[1],"POS")) || (0 != strcmp(field[2],"ID")) ||
              (0 != strcmp(field[3],"REF")) || (0 != strcmp(field[4],"ALT")) ||
              (0 != strcmp(field[5],"QUAL")) || (0 != strcmp(field[6],"FILTER")) ||
              (0 != strcmp(field[7],"INFO")) || (0 != strcmp(field[8],"FORMAT")))
            {
              gs_info(gv,-2,"VCF header must contain standard columns through FORMAT and at least one sample");
              failed = 1;
              goto cleanup;
            }
          if (nfield - 9u > (size_t)gs_UNSIGNED_MAX)
            {
              gs_info(gv,-2,"Too many VCF samples");
              failed = 1;
              goto cleanup;
            }
          no_ind = (gs_UNSIGNED)(nfield - 9u);
          sample_name = calloc((size_t)no_ind,sizeof(*sample_name));
          if (NULL == sample_name)
            {
              gs_info(gv,-2,"Memory allocation error");
              failed = 1;
              goto cleanup;
            }
          for (gs_UNSIGNED i = 0u; i < no_ind; i++)
            {
              const char *name = field[9u + i];
              if (!st_vcf_token_name_valid(name) || (strlen(name) >= gs_MAXNME))
                {
                  gs_info(gv,-2,"Invalid or overlong VCF sample name");
                  failed = 1;
                  goto cleanup;
                }
              sample_name[i] = st_vcf_strdup(name);
              if (NULL == sample_name[i])
                {
                  gs_info(gv,-2,"Memory allocation error");
                  failed = 1;
                  goto cleanup;
                }
            }
          have_header = 1;
          continue;
        }

      if (!have_header)
        {
          gs_info(gv,-2,"VCF data encountered before #CHROM header");
          failed = 1;
          goto cleanup;
        }

      if (0 != st_vcf_split_tabs(line,&field,&field_capacity,&nfield))
        {
          failed = 1;
          goto cleanup;
        }
      if (nfield != (size_t)no_ind + 9u)
        {
          gs_info(gv,-2,"VCF record has the wrong number of columns");
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
            (0 != st_vcf_parse_pos(field[1],&pos)) ||
            ('\0' == field[3][0]) || (0 == strcmp(field[3],".")) ||
            (0 != st_stvcf_alleles(field[3],field[4],&n_alt)) ||
            (0 != st_vcf_format_gt_index(field[8],&gt_index)))
          {
            gs_info(gv,-2,"Invalid CHROM/POS/REF/ALT/FORMAT field in VCF");
            failed = 1;
            goto cleanup;
          }

        if (0 == strcmp(id,"."))
          {
            if (0 != st_vcf_generated_marker_name(generated,sizeof(generated),no_mar))
              {
                failed = 1;
                goto cleanup;
              }
            id = generated;
          }
        else if (!st_vcf_token_name_valid(id) || (strlen(id) >= gs_MAXNME))
          {
            gs_info(gv,-2,"Invalid or overlong VCF marker ID");
            failed = 1;
            goto cleanup;
          }

        if ((0 != st_vcf_map_pool_add(&map_string_size,id)) ||
            (0 != st_vcf_map_pool_add(&map_string_size,"m")))
          {
            gs_info(gv,-2,"VCF map string storage too large");
            failed = 1;
            goto cleanup;
          }
      }

      if (no_mar == gs_UNSIGNED_MAX)
        {
          gs_info(gv,-2,"Too many VCF markers");
          failed = 1;
          goto cleanup;
        }
      no_mar++;
    }

  if (!have_header || (0u == no_ind) || (0u == no_mar))
    {
      gs_info(gv,-2,"VCF contains no samples or no markers");
      failed = 1;
      goto cleanup;
    }

  if (0 != fclose(fp))
    {
      fp = NULL;
      gs_info(gv,-2,"Error closing VCF file");
      failed = 1;
      goto cleanup;
    }
  fp = NULL;

  tmp.NoInd = no_ind;
  tmp.NoMar = no_mar;
  tmp.NLoc = no_mar;
  tmp.map_str_N = map_string_size;
  tmp.IndName = gs_malloc_array((size_t)no_ind,sizeof(*tmp.IndName));
  tmp.MarName = gs_malloc_array((size_t)no_mar,sizeof(*tmp.MarName));
  tmp.map_S = gs_malloc_array((size_t)no_mar,sizeof(*tmp.map_S));
  tmp.map_str_S = gs_malloc_array(tmp.map_str_N,sizeof(*tmp.map_str_S));
  tmp.map = gs_malloc_array((size_t)no_mar,sizeof(*tmp.map));
  tmp.mdl = gs_malloc_array((size_t)no_mar,sizeof(*tmp.mdl));
  if ((NULL == tmp.IndName) || (NULL == tmp.MarName) ||
      (NULL == tmp.map_S) || (NULL == tmp.map_str_S) ||
      (NULL == tmp.map) || (NULL == tmp.mdl) ||
      (0 != gs_alloc_genotype_matrix(&tmp,no_mar,no_ind)))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED i = 0u; i < no_ind; i++)
    {
      size_t len = strlen(sample_name[i]);
      strcpy(tmp.IndName[i],sample_name[i]);
      if ((int)len > tmp.MaxIL) tmp.MaxIL = (int)len;
    }
  if ((0 != gs_make_name_index(tmp.IndName,tmp.NoInd,&index)) ||
      gs_name_index_has_duplicates(index,tmp.NoInd))
    {
      gs_info(gv,-2,"Duplicate or invalid VCF sample names");
      failed = 1;
      goto cleanup;
    }
  free(index);
  index = NULL;

  fp = fopen(CHAR(STRING_ELT(filename,0)),"rb");
  if (NULL == fp)
    {
      gs_info(gv,-2,"Cannot reopen VCF file");
      failed = 1;
      goto cleanup;
    }

  {
    gs_UNSIGNED marker = 0u;
    int saw_header = 0;

    for (;;)
      {
        int lr = st_vcf_getline(fp,&line,&line_capacity);
        size_t nfield = 0u;

        if (0 == lr)
          break;
        if (lr < 0)
          {
            gs_info(gv,-2,"Error reading VCF file");
            failed = 1;
            goto cleanup;
          }
        if ('\0' == line[0])
          continue;
        if (0 == strncmp(line,"##",2))
          continue;
        if ('#' == line[0])
          {
            saw_header = 1;
            continue;
          }
        if (!saw_header || (marker >= no_mar) ||
            (0 != st_vcf_split_tabs(line,&field,&field_capacity,&nfield)) ||
            (nfield != (size_t)no_ind + 9u))
          {
            gs_info(gv,-2,"Invalid VCF structure on second pass");
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

          if ((0 != st_vcf_parse_pos(field[1],&raw_pos)) ||
              (0 != st_stvcf_alleles(field[3],field[4],&n_alt)) ||
              (0 != st_vcf_format_gt_index(field[8],&gt_index)))
            {
              gs_info(gv,-2,"Invalid VCF record on second pass");
              failed = 1;
              goto cleanup;
            }

          if (0 == strcmp(id,"."))
            {
              if (0 != st_vcf_generated_marker_name(generated,sizeof(generated),marker))
                {
                  failed = 1;
                  goto cleanup;
                }
              id = generated;
            }

          if ((0u == no_chrom) ||
              (0 != strcmp(chrom_name[no_chrom - 1u],field[0])))
            {
              for (gs_UNSIGNED c = 0u; c < no_chrom; c++)
                if (0 == strcmp(chrom_name[c],field[0]))
                  {
                    gs_info(gv,-2,"VCF chromosome records must be contiguous");
                    failed = 1;
                    goto cleanup;
                  }

              if ((size_t)no_chrom >= chrom_capacity)
                {
                  size_t next = (0u == chrom_capacity) ? 8u : 2u * chrom_capacity;
                  char **z;
                  if ((next <= chrom_capacity) ||
                      (next > SIZE_MAX / sizeof(*chrom_name)))
                    {
                      failed = 1;
                      goto cleanup;
                    }
                  z = realloc(chrom_name,next * sizeof(*chrom_name));
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
          if ((raw_pos < previous_raw_pos) && (chrom_code == previous_chrom))
            {
              gs_info(gv,-2,"VCF positions must be nondecreasing within chromosome");
              failed = 1;
              goto cleanup;
            }
          previous_raw_pos = raw_pos;
          previous_chrom = chrom_code;

          genetic_pos = (position_scale > 0.0) ?
            (((double)raw_pos - 1.0) / position_scale) : (double)raw_pos;
          if (!isfinite(genetic_pos) || (genetic_pos < 0.0))
            {
              gs_info(gv,-2,"VCF position cannot be represented as a genetic map position");
              failed = 1;
              goto cleanup;
            }

          strcpy(tmp.MarName[marker],id);
          if ((int)strlen(id) > tmp.MaxML) tmp.MaxML = (int)strlen(id);
          source_point.chrom = chrom_code;
          source_point.pos = genetic_pos;
          source_point.name = id;
          source_point.class = "m";
          source_point.nr = marker;
          if (0 != gs_copy_map_point_to_pool(&(tmp.map_S[marker]),
                                             &source_point,marker,
                                             tmp.map_str_S,tmp.map_str_N,
                                             &map_string_offset))
            {
              gs_info(gv,-2,"Invalid VCF map string storage");
              failed = 1;
              goto cleanup;
            }
          tmp.map[marker] = &(tmp.map_S[marker]);
          tmp.mdl[marker] = marker;

          for (gs_UNSIGNED i = 0u; i < no_ind; i++)
            {
              const char *gt = NULL;
              size_t gt_len = 0u;
              GS_ALLELE_TYPE a1, a2;
              size_t idx = (size_t)marker * (size_t)no_ind + i;

              if ((0 != st_vcf_subfield(field[9u + i],gt_index,&gt,&gt_len)) ||
                  (0u == gt_len) ||
                  (0 != st_vcf_parse_gt(gt,gt_len,n_alt,&a1,&a2)))
                {
                  gs_info(gv,-2,"Invalid diploid GT field in VCF");
                  failed = 1;
                  goto cleanup;
                }
              tmp.All1[idx] = a1;
              tmp.All2[idx] = a2;
            }
        }
        marker++;
      }

    if (marker != no_mar)
      {
        gs_info(gv,-2,"VCF marker count changed between passes");
        failed = 1;
        goto cleanup;
      }
  }

  if (map_string_offset != tmp.map_str_N)
    {
      gs_info(gv,-2,"Invalid VCF map string storage");
      failed = 1;
      goto cleanup;
    }

  if (0 != st_stvcf_separate_map_ties(gv,tmp.map,tmp.NLoc))
    {
      failed = 1;
      goto cleanup;
    }

  if ((0 != gs_make_name_index(tmp.MarName,tmp.NoMar,&index)) ||
      gs_name_index_has_duplicates(index,tmp.NoMar))
    {
      gs_info(gv,-2,"Duplicate VCF marker IDs");
      failed = 1;
      goto cleanup;
    }
  free(index);
  index = NULL;

  {
    int srv = 0;
    char *empty = "";
    char **empty_p = &empty;
    int zero = 0;
    int *zero_p = &zero;

    gs_marker_stats_01(&tmp,empty_p,empty_p,empty_p,
                       zero_p,zero_p,zero_p,zero_p,&srv);
    if (-2 == srv)
      {
        failed = 1;
        goto cleanup;
      }
  }

  gs_tmp_commit(gv,&tmp);
  status = 0;

 cleanup:
  if (NULL != fp) fclose(fp);
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
  memset(out,0,sizeof(*out));
  out->plain = fopen(filename,"wb");
  return (NULL == out->plain) ? -1 : 0;
}


static int st_vcf_output_write(st_vcf_output *out, const char *text,
                               size_t n)
{
  return (fwrite(text,1u,n,out->plain) == n) ? 0 : -1;
}


static int st_vcf_output_printf(st_vcf_output *out, const char *format, ...)
{
  char stack[4096];
  char *buffer = stack;
  size_t capacity = sizeof(stack);
  va_list ap;
  int n;
  int rc;

  va_start(ap,format);
  n = vsnprintf(buffer,capacity,format,ap);
  va_end(ap);
  if (n < 0)
    return -1;

  if ((size_t)n >= capacity)
    {
      capacity = (size_t)n + 1u;
      buffer = malloc(capacity);
      if (NULL == buffer)
        return -1;
      va_start(ap,format);
      n = vsnprintf(buffer,capacity,format,ap);
      va_end(ap);
      if ((n < 0) || ((size_t)n >= capacity))
        {
          free(buffer);
          return -1;
        }
    }

  rc = st_vcf_output_write(out,buffer,(size_t)n);
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

  memset(&out,0,sizeof(out));

  if ((STRSXP != TYPEOF(filename)) || (1 != XLENGTH(filename)) ||
      (NA_STRING == STRING_ELT(filename,0)) ||
      ('\0' == CHAR(STRING_ELT(filename,0))[0]) ||
      (STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
      (NA_STRING == STRING_ELT(dataset,0)) ||
      ('\0' == CHAR(STRING_ELT(dataset,0))[0]))
    return ScalarInteger(-2);

  gv = gs_fdta((char *)CHAR(STRING_ELT(dataset,0)));
  if ((0u == gv->NoInd) || (0u == gv->NoMar) ||
      (NULL == gv->IndName) || (NULL == gv->MarName) ||
      (NULL == gv->All1) || (NULL == gv->All2) ||
      (0u == gv->NLoc) || (NULL == gv->map) || (NULL == gv->mdl))
    {
      gs_info(gv,-2,"No complete marker data/linkage map available for VCF output");
      return ScalarInteger(-2);
    }

  order = gs_malloc_array((size_t)gv->NoMar,sizeof(*order));
  point_by_marker = gs_malloc_array((size_t)gv->NoMar,
                                    sizeof(*point_by_marker));
  seen = calloc((size_t)gv->NoMar,sizeof(*seen));
  alleles = gs_malloc_array((size_t)gs_MAXALL,sizeof(*alleles));
  if ((NULL == order) || (NULL == point_by_marker) ||
      (NULL == seen) || (NULL == alleles))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED l = 0u; l < gv->NLoc; l++)
    {
      gs_UNSIGNED m = gv->mdl[l];
      if (gs_UNSIGNED_MAX == m)
        continue;
      if ((m >= gv->NoMar) || (NULL == gv->map[l]) || seen[m] ||
          (0 != strcmp(gv->map[l]->name,gv->MarName[m])) ||
          (gv->map[l]->chrom < 1) || !isfinite(gv->map[l]->pos) ||
          (gv->map[l]->pos < 0.0))
        {
          gs_info(gv,-2,"Marker data do not have one valid map position per marker");
          failed = 1;
          goto cleanup;
        }
      seen[m] = 1u;
      point_by_marker[m] = gv->map[l];
      order[norder++] = m;
    }
  if (norder != gv->NoMar)
    {
      gs_info(gv,-2,"Marker data do not have one map position per marker");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED i = 0u; i < gv->NoInd; i++)
    if (!st_vcf_token_name_valid(gv->IndName[i]))
      {
        gs_info(gv,-2,"Invalid sample name for VCF output");
        failed = 1;
        goto cleanup;
      }
  for (gs_UNSIGNED m = 0u; m < gv->NoMar; m++)
    if (!st_vcf_token_name_valid(gv->MarName[m]) ||
        (strlen(gv->MarName[m]) >= gs_MAXNME))
      {
        gs_info(gv,-2,"Invalid marker name for VCF output");
        failed = 1;
        goto cleanup;
      }

  if (0 != st_vcf_output_open(&out,CHAR(STRING_ELT(filename,0))))
    {
      gs_info(gv,-2,"Cannot open VCF output file");
      failed = 1;
      goto cleanup;
    }
  opened = 1;

  if ((0 != st_vcf_output_printf(&out,"##fileformat=VCFv4.3\n")) ||
      (0 != st_vcf_output_printf(&out,"##source=SelectionTools\n")) ||
      (0 != st_vcf_output_printf(&out,
          "##SelectionToolsPositionScale=1000000\n")) ||
      (0 != st_vcf_output_printf(&out,
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")) ||
      (0 != st_vcf_output_printf(&out,
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")))
    {
      failed = 1;
      goto cleanup;
    }
  for (gs_UNSIGNED i = 0u; i < gv->NoInd; i++)
    if (0 != st_vcf_output_printf(&out,"\t%s",gv->IndName[i]))
      {
        failed = 1;
        goto cleanup;
      }
  if (0 != st_vcf_output_printf(&out,"\n"))
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

      if (0 != st_markerdata_observed_alleles(gv,m,alleles,
                                              &n_alleles,&has_missing))
        {
          gs_info(gv,-2,"Invalid marker alleles for VCF output");
          failed = 1;
          goto cleanup;
        }
      (void)has_missing;

      if (point->pos > ((double)LLONG_MAX - 2.0) / ST_VCF_POSITION_SCALE)
        {
          gs_info(gv,-2,"Genetic position is too large for Beagle VCF scaling");
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
        if ((0 != st_vcf_allele_label(0u,label,sizeof(label))) ||
            (0 != st_vcf_output_printf(&out,"%d\t%lld\t%s\t%s\t",
                                        point->chrom,bp,gv->MarName[m],label)))
          {
            failed = 1;
            goto cleanup;
          }
      }
      if (n_alleles <= 1u)
        {
          if (0 != st_vcf_output_printf(&out,"."))
            {
              failed = 1;
              goto cleanup;
            }
        }
      else
        for (gs_UNSIGNED a = 1u; a < n_alleles; a++)
          {
            char label[40];
            if ((0 != st_vcf_allele_label(a,label,sizeof(label))) ||
                (0 != st_vcf_output_printf(&out,(1u == a) ? "%s" : ",%s",
                                            label)))
              {
                failed = 1;
                goto cleanup;
              }
          }

      if (0 != st_vcf_output_printf(&out,"\t.\tPASS\t.\tGT"))
        {
          failed = 1;
          goto cleanup;
        }

      for (gs_UNSIGNED i = 0u; i < gv->NoInd; i++)
        {
          size_t idx = (size_t)m * (size_t)gv->NoInd + i;
          GS_ALLELE_TYPE a1 = gv->All1[idx];
          GS_ALLELE_TYPE a2 = gv->All2[idx];

          if ((GS_ALLELE_MISSING == a1) || (GS_ALLELE_MISSING == a2))
            {
              if (0 != st_vcf_output_printf(&out,"\t./."))
                {
                  failed = 1;
                  goto cleanup;
                }
            }
          else
            {
              int i1 = st_markerdata_allele_index(alleles,n_alleles,a1);
              int i2 = st_markerdata_allele_index(alleles,n_alleles,a2);
              if ((i1 < 0) || (i2 < 0) ||
                  (0 != st_vcf_output_printf(&out,"\t%d/%d",i1,i2)))
                {
                  gs_info(gv,-2,"Allele normalization failed during VCF output");
                  failed = 1;
                  goto cleanup;
                }
            }
        }
      if (0 != st_vcf_output_printf(&out,"\n"))
        {
          failed = 1;
          goto cleanup;
        }
    }

  if (0 != st_vcf_output_close(&out))
    {
      opened = 0;
      gs_info(gv,-2,"Error closing VCF output file");
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
  if (failed && (0 == status)) status = -2;
  return ScalarInteger(status);
}


SEXP st_load_vcf_data(SEXP stvcf, SEXP dataset)
{
  gs_varset_type *gv;
  gs_varset_type tmp;
  SEXP chrom_s, pos_s, id_s, ref_s, alt_s, ind_s;
  SEXP allele1_s, allele2_s, chrom_levels_s;
  gs_mappoint_type *map_s = NULL;
  gs_mappoint_type **map = NULL;
  gs_UNSIGNED *n_alt = NULL;
  gs_name_index_type *index = NULL;
  unsigned char *chrom_seen = NULL;
  R_xlen_t no_mar_x, no_ind_x;
  gs_UNSIGNED no_mar, no_ind, no_chrom;
  int failed = 0;
  int status = -2;

  if ((STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
      (NA_STRING == STRING_ELT(dataset,0)))
    return ScalarInteger(-2);

  gv = gs_fdta((char *)CHAR(STRING_ELT(dataset,0)));
  gs_tmp_init_from(&tmp,gv);

  if (!st_stvcf_components(stvcf))
    {
      gs_info(gv,-2,
              "STvcf must contain exactly the nine required components");
      failed = 1;
      goto cleanup;
    }

  chrom_s = st_stvcf_element(stvcf,"CHROM");
  pos_s = st_stvcf_element(stvcf,"POS");
  id_s = st_stvcf_element(stvcf,"ID");
  ref_s = st_stvcf_element(stvcf,"REF");
  alt_s = st_stvcf_element(stvcf,"ALT");
  ind_s = st_stvcf_element(stvcf,"individual");
  allele1_s = st_stvcf_element(stvcf,"allele1");
  allele2_s = st_stvcf_element(stvcf,"allele2");
  chrom_levels_s = st_stvcf_element(stvcf,"CHROM.levels");

  if ((INTSXP != TYPEOF(chrom_s)) ||
      (REALSXP != TYPEOF(pos_s)) ||
      (STRSXP != TYPEOF(id_s)) ||
      (STRSXP != TYPEOF(ref_s)) ||
      (STRSXP != TYPEOF(alt_s)) ||
      (STRSXP != TYPEOF(ind_s)) ||
      (STRSXP != TYPEOF(chrom_levels_s)) ||
      !(((RAWSXP == TYPEOF(allele1_s)) &&
         (RAWSXP == TYPEOF(allele2_s))) ||
        ((INTSXP == TYPEOF(allele1_s)) &&
         (INTSXP == TYPEOF(allele2_s)))))
    {
      gs_info(gv,-2,"Invalid STvcf component types");
      failed = 1;
      goto cleanup;
    }

  no_mar_x = XLENGTH(id_s);
  no_ind_x = XLENGTH(ind_s);
  if ((0 == no_mar_x) || (0 == no_ind_x) ||
      ((R_xlen_t)(gs_UNSIGNED)no_mar_x != no_mar_x) ||
      ((R_xlen_t)(gs_UNSIGNED)no_ind_x != no_ind_x))
    {
      gs_info(gv,-2,"Invalid STvcf dimensions");
      failed = 1;
      goto cleanup;
    }

  no_mar = (gs_UNSIGNED)no_mar_x;
  no_ind = (gs_UNSIGNED)no_ind_x;
  no_chrom = (gs_UNSIGNED)XLENGTH(chrom_levels_s);

  if ((0 == no_chrom) ||
      (XLENGTH(chrom_s) != no_mar_x) ||
      (XLENGTH(pos_s) != no_mar_x) ||
      (XLENGTH(ref_s) != no_mar_x) ||
      (XLENGTH(alt_s) != no_mar_x) ||
      !st_stvcf_matrix_dim(allele1_s,no_mar_x,no_ind_x) ||
      !st_stvcf_matrix_dim(allele2_s,no_mar_x,no_ind_x))
    {
      gs_info(gv,-2,"Inconsistent STvcf component lengths or dimensions");
      failed = 1;
      goto cleanup;
    }

  map_s = gs_malloc_array((size_t)no_mar,sizeof(*map_s));
  map = gs_malloc_array((size_t)no_mar,sizeof(*map));
  n_alt = gs_malloc_array((size_t)no_mar,sizeof(*n_alt));
  chrom_seen = calloc((size_t)no_chrom,sizeof(*chrom_seen));
  if ((NULL == map_s) || (NULL == map) || (NULL == n_alt) ||
      (NULL == chrom_seen))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED c = 0; c < no_chrom; c++)
    {
      SEXP s = STRING_ELT(chrom_levels_s,c);

      if ((NA_STRING == s) || ('\0' == CHAR(s)[0]))
        {
          gs_info(gv,-2,"Invalid CHROM.levels in STvcf");
          failed = 1;
          goto cleanup;
        }

      for (gs_UNSIGNED d = 0; d < c; d++)
        if (0 == strcmp(CHAR(s),CHAR(STRING_ELT(chrom_levels_s,d))))
          {
            gs_info(gv,-2,"Duplicate chromosome label in STvcf");
            failed = 1;
            goto cleanup;
          }
    }

  for (gs_UNSIGNED m = 0; m < no_mar; m++)
    {
      int chrom = INTEGER(chrom_s)[m];
      double pos = REAL(pos_s)[m];
      SEXP id = STRING_ELT(id_s,m);
      SEXP ref = STRING_ELT(ref_s,m);
      SEXP alt = STRING_ELT(alt_s,m);
      size_t id_len;

      if ((NA_INTEGER == chrom) || (chrom < 1) ||
          ((gs_UNSIGNED)chrom > no_chrom) || !isfinite(pos) ||
          (pos < 0.0))
        {
          gs_info(gv,-2,"Invalid chromosome or position in STvcf");
          failed = 1;
          goto cleanup;
        }
      chrom_seen[chrom - 1] = 1;

      if ((NA_STRING == id) ||
          !st_stvcf_valid_name(CHAR(id)))
        {
          gs_info(gv,-2,"Invalid marker ID in STvcf");
          failed = 1;
          goto cleanup;
        }
      id_len = strlen(CHAR(id));
      if (id_len >= gs_MAXNME)
        {
          gs_info(gv,-2,"Marker ID too long in STvcf");
          failed = 1;
          goto cleanup;
        }

      if ((NA_STRING == ref) || ('\0' == CHAR(ref)[0]) ||
          (0 == strcmp(CHAR(ref),".")) || (NA_STRING == alt) ||
          (0 != st_stvcf_alleles(CHAR(ref),CHAR(alt),&(n_alt[m]))))
        {
          gs_info(gv,-2,"Invalid REF or ALT field in STvcf");
          failed = 1;
          goto cleanup;
        }

      map_s[m].chrom = chrom;
      map_s[m].pos = pos;
      map_s[m].name = CHAR(id);
      map_s[m].class = "m";
      map_s[m].nr = m;
      map[m] = &(map_s[m]);

      if ((0 != gs_pool_add_string_size(&(tmp.map_str_N),
                                        map_s[m].name)) ||
          (0 != gs_pool_add_string_size(&(tmp.map_str_N),
                                        map_s[m].class)))
        {
          gs_info(gv,-2,"STvcf map string storage too large");
          failed = 1;
          goto cleanup;
        }
    }

  for (gs_UNSIGNED c = 0; c < no_chrom; c++)
    if (!chrom_seen[c])
      {
        gs_info(gv,-2,"STvcf chromosome codes must be consecutive");
        failed = 1;
        goto cleanup;
      }

  qsort(map,no_mar,sizeof(*map),st_stvcf_map_cmp);
  if (0 != st_stvcf_separate_map_ties(gv,map,no_mar))
    {
      failed = 1;
      goto cleanup;
    }

  tmp.NoInd = no_ind;
  tmp.NoMar = no_mar;
  tmp.NLoc = no_mar;
  tmp.IndName = gs_malloc_array((size_t)no_ind,sizeof(*tmp.IndName));
  tmp.MarName = gs_malloc_array((size_t)no_mar,sizeof(*tmp.MarName));
  tmp.map_S = gs_malloc_array((size_t)no_mar,sizeof(*tmp.map_S));
  tmp.map_str_S = gs_malloc_array(tmp.map_str_N,sizeof(*tmp.map_str_S));
  tmp.map = gs_malloc_array((size_t)no_mar,sizeof(*tmp.map));
  tmp.mdl = gs_malloc_array((size_t)no_mar,sizeof(*tmp.mdl));

  if ((NULL == tmp.IndName) || (NULL == tmp.MarName) ||
      (NULL == tmp.map_S) || (NULL == tmp.map_str_S) ||
      (NULL == tmp.map) || (NULL == tmp.mdl) ||
      (0 != gs_alloc_genotype_matrix(&tmp,no_mar,no_ind)))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED i = 0; i < no_ind; i++)
    {
      SEXP s = STRING_ELT(ind_s,i);
      size_t len;

      if ((NA_STRING == s) || !st_stvcf_valid_name(CHAR(s)))
        {
          gs_info(gv,-2,"Invalid individual name in STvcf");
          failed = 1;
          goto cleanup;
        }
      len = strlen(CHAR(s));
      if (len >= gs_MAXNME)
        {
          gs_info(gv,-2,"Individual name too long in STvcf");
          failed = 1;
          goto cleanup;
        }
      strcpy(tmp.IndName[i],CHAR(s));
      if ((int)len > tmp.MaxIL) tmp.MaxIL = (int)len;
    }

  if (0 != gs_make_name_index(tmp.IndName,tmp.NoInd,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,tmp.NoInd))
    {
      gs_info(gv,-2,"Duplicate individual name in STvcf");
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

        strcpy(tmp.MarName[m_new],map[m_new]->name);
        if ((int)name_len > tmp.MaxML) tmp.MaxML = (int)name_len;

        if (0 != gs_copy_map_point_to_pool(&(tmp.map_S[m_new]),
                                            map[m_new],
                                            m_new,
                                            tmp.map_str_S,
                                            tmp.map_str_N,
                                            &map_str_offset))
          {
            gs_info(gv,-2,"Invalid STvcf map string storage");
            failed = 1;
            goto cleanup;
          }
        tmp.map[m_new] = &(tmp.map_S[m_new]);
        tmp.mdl[m_new] = m_new;

        for (gs_UNSIGNED i = 0; i < no_ind; i++)
        {
          R_xlen_t ridx = (R_xlen_t)m_old +
                          (R_xlen_t)no_mar * (R_xlen_t)i;
          size_t idx = (size_t)m_new * (size_t)no_ind + i;
          int miss1, miss2;
          GS_ALLELE_TYPE a1, a2;

          if ((0 != st_stvcf_get_allele(allele1_s,ridx,n_alt[m_old],
                                        &miss1,&a1)) ||
              (0 != st_stvcf_get_allele(allele2_s,ridx,n_alt[m_old],
                                        &miss2,&a2)))
            {
              gs_info(gv,-2,"Invalid allele index in STvcf");
              failed = 1;
              goto cleanup;
            }

          if (miss1 || miss2)
            {
              tmp.All1[idx] = GS_ALLELE_MISSING;
              tmp.All2[idx] = GS_ALLELE_MISSING;
            }
          else
            {
              tmp.All1[idx] = a1;
              tmp.All2[idx] = a2;
            }
        }
      }

    if (map_str_offset != tmp.map_str_N)
      {
        gs_info(gv,-2,"Invalid STvcf map string storage");
        failed = 1;
        goto cleanup;
      }
  }

  if (0 != gs_make_name_index(tmp.MarName,tmp.NoMar,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,tmp.NoMar))
    {
      gs_info(gv,-2,"Duplicate marker ID in STvcf");
      failed = 1;
      goto cleanup;
    }
  free(index);
  index = NULL;

  {
    int srv = 0;
    char *empty = "";
    char **empty_p = &empty;
    int zero = 0;
    int *zero_p = &zero;

    gs_marker_stats_01(&tmp,
                       empty_p,empty_p,empty_p,
                       zero_p,zero_p,zero_p,zero_p,
                       &srv);
    if (-2 == srv)
      {
        failed = 1;
        goto cleanup;
      }
  }

  gs_tmp_commit(gv,&tmp);
  status = 0;

 cleanup:
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
  SEXP names_s, class_s;

  if ((STRSXP != TYPEOF(dataset)) || (1 != XLENGTH(dataset)) ||
      (NA_STRING == STRING_ELT(dataset,0)))
    return R_NilValue;

  gv = gs_fdta((char *)CHAR(STRING_ELT(dataset,0)));

  if ((0 == gv->NoInd) || (0 == gv->NoMar) ||
      (NULL == gv->IndName) || (NULL == gv->MarName) ||
      (NULL == gv->All1) || (NULL == gv->All2))
    {
      gs_info(gv,-2,"No marker data available");
      return R_NilValue;
    }

  if ((gv->NoInd > INT_MAX) || (gv->NoMar > INT_MAX) ||
      ((size_t)gv->NoMar > SIZE_MAX / (size_t)gv->NoInd))
    {
      gs_info(gv,-2,"Marker data are too large for STvcf conversion");
      return R_NilValue;
    }

  if ((0 == gv->NLoc) || (NULL == gv->map_S) ||
      (NULL == gv->map_str_S) || (0 == gv->map_str_N) ||
      (NULL == gv->map) || (NULL == gv->mdl))
    {
      gs_info(gv,-2,"No linkage map available");
      return R_NilValue;
    }

  if (!gs_map_pool_valid(gv->map_S,
                         gv->NLoc,
                         gv->map_str_S,
                         gv->map_str_N))
    {
      gs_info(gv,-2,"Invalid linkage map string storage");
      return R_NilValue;
    }

  map_chrom = gs_malloc_array((size_t)gv->NoMar,sizeof(*map_chrom));
  map_pos = gs_malloc_array((size_t)gv->NoMar,sizeof(*map_pos));
  map_seen = calloc((size_t)gv->NoMar,sizeof(*map_seen));
  chrom_levels = gs_malloc_array((size_t)gv->NoMar,
                                 sizeof(*chrom_levels));
  alleles = gs_malloc_array((size_t)gs_MAXALL,sizeof(*alleles));
  if ((NULL == map_chrom) || (NULL == map_pos) || (NULL == map_seen) ||
      (NULL == chrom_levels) || (NULL == alleles))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
    {
      const char *name = gv->IndName[i];

      if (!st_stvcf_valid_name(name) ||
          (strlen(name) >= gs_MAXNME))
        {
          gs_info(gv,-2,
                  "Individual names are not valid for STvcf conversion");
          failed = 1;
          goto cleanup;
        }
    }

  if (0 != gs_make_name_index(gv->IndName,gv->NoInd,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,gv->NoInd))
    {
      gs_info(gv,-2,"Duplicate individual name in marker data");
      failed = 1;
      goto cleanup;
    }
  free(index);
  index = NULL;

  for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
      const char *name = gv->MarName[m];

      if (!st_stvcf_valid_name(name) ||
          (strlen(name) >= gs_MAXNME))
        {
          gs_info(gv,-2,
                  "Marker names are not valid for STvcf conversion");
          failed = 1;
          goto cleanup;
        }
    }

  if (0 != gs_make_name_index(gv->MarName,gv->NoMar,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,gv->NoMar))
    {
      gs_info(gv,-2,"Duplicate marker name in marker data");
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
          gs_info(gv,-2,"Invalid linkage map state");
          failed = 1;
          goto cleanup;
        }

      point = gv->map[l];
      if (map_seen[m] || (point->chrom < 1) ||
          !isfinite(point->pos) || (point->pos < 0.0) ||
          (0 != strcmp(point->name,gv->MarName[m])))
        {
          gs_info(gv,-2,
                  "Marker data do not have one valid map position per marker");
          failed = 1;
          goto cleanup;
        }

      map_seen[m] = 1;
      map_chrom[m] = point->chrom;
      map_pos[m] = point->pos;
    }

  for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    if (!map_seen[m])
      {
        gs_info(gv,-2,
                "Marker data do not have one valid map position per marker");
        failed = 1;
        goto cleanup;
      }

  for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
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

      if (0 != st_markerdata_observed_alleles(gv,m,alleles,
                                              &n_alleles,&has_missing))
        {
          gs_info(gv,-2,
                  "Marker alleles are not valid for STvcf conversion");
          failed = 1;
          goto cleanup;
        }

    }

  out = PROTECT(allocVector(VECSXP,9)); nprot++;
  chrom_s = PROTECT(allocVector(INTSXP,gv->NoMar)); nprot++;
  pos_s = PROTECT(allocVector(REALSXP,gv->NoMar)); nprot++;
  id_s = PROTECT(allocVector(STRSXP,gv->NoMar)); nprot++;
  ref_s = PROTECT(allocVector(STRSXP,gv->NoMar)); nprot++;
  alt_s = PROTECT(allocVector(STRSXP,gv->NoMar)); nprot++;
  ind_s = PROTECT(allocVector(STRSXP,gv->NoInd)); nprot++;

  allele1_s = PROTECT(allocMatrix(RAWSXP,gv->NoMar,gv->NoInd));
  nprot++;
  allele2_s = PROTECT(allocMatrix(RAWSXP,gv->NoMar,gv->NoInd));
  nprot++;
  memset(RAW(allele1_s),255,(size_t)gv->NoMar * gv->NoInd);
  memset(RAW(allele2_s),255,(size_t)gv->NoMar * gv->NoInd);

  chrom_levels_s = PROTECT(allocVector(STRSXP,no_chrom)); nprot++;
  names_s = PROTECT(allocVector(STRSXP,9)); nprot++;
  class_s = PROTECT(mkString("STvcf")); nprot++;

  for (gs_UNSIGNED c = 0; c < no_chrom; c++)
    {
      char label[64];
      snprintf(label,sizeof(label),"%d",chrom_levels[c]);
      SET_STRING_ELT(chrom_levels_s,c,mkChar(label));
    }

  for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
    SET_STRING_ELT(ind_s,i,mkChar(gv->IndName[i]));

  for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
      gs_UNSIGNED n_alleles = 0;
      int has_missing = 0;
      int chrom_code = -1;

      if (0 != st_markerdata_observed_alleles(gv,m,alleles,
                                              &n_alleles,&has_missing))
        {
          gs_info(gv,-2,
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

      if (chrom_code < 1)
        {
          gs_info(gv,-2,"Invalid chromosome mapping");
          failed = 1;
          goto cleanup;
        }

      INTEGER(chrom_s)[m] = chrom_code;
      REAL(pos_s)[m] = map_pos[m];
      SET_STRING_ELT(id_s,m,mkChar(gv->MarName[m]));
      SET_STRING_ELT(ref_s,m,mkChar("1"));

      if (n_alleles <= 1)
        SET_STRING_ELT(alt_s,m,mkChar("."));
      else
        {
          size_t capacity = (size_t)n_alleles * 12 + 1;
          char *alt = malloc(capacity);
          size_t used = 0;

          if (NULL == alt)
            {
              gs_info(gv,-2,"Memory allocation error");
              failed = 1;
              goto cleanup;
            }

          alt[0] = '\0';
          for (gs_UNSIGNED a = 2; a <= n_alleles; a++)
            {
              int n = snprintf(alt + used,capacity - used,
                               (2 == a) ? "%lu" : ",%lu",a);

              if ((n < 0) || ((size_t)n >= capacity - used))
                {
                  free(alt);
                  gs_info(gv,-2,"ALT field construction failed");
                  failed = 1;
                  goto cleanup;
                }
              used += (size_t)n;
            }

          SET_STRING_ELT(alt_s,m,mkChar(alt));
          free(alt);
        }

      for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        {
          size_t old_idx = (size_t)m * (size_t)gv->NoInd + i;
          R_xlen_t ridx = (R_xlen_t)m +
                          (R_xlen_t)gv->NoMar * (R_xlen_t)i;
          GS_ALLELE_TYPE a1 = gv->All1[old_idx];
          GS_ALLELE_TYPE a2 = gv->All2[old_idx];

          if ((GS_ALLELE_MISSING == a1) ||
              (GS_ALLELE_MISSING == a2))
            continue;

          {
            int i1 = st_markerdata_allele_index(alleles,n_alleles,a1);
            int i2 = st_markerdata_allele_index(alleles,n_alleles,a2);

            if ((i1 < 0) || (i2 < 0))
              {
                gs_info(gv,-2,"Allele normalization failed");
                failed = 1;
                goto cleanup;
              }

            RAW(allele1_s)[ridx] = (unsigned char)i1;
            RAW(allele2_s)[ridx] = (unsigned char)i2;
          }
        }
    }

  SET_VECTOR_ELT(out,0,chrom_s);
  SET_VECTOR_ELT(out,1,pos_s);
  SET_VECTOR_ELT(out,2,id_s);
  SET_VECTOR_ELT(out,3,ref_s);
  SET_VECTOR_ELT(out,4,alt_s);
  SET_VECTOR_ELT(out,5,ind_s);
  SET_VECTOR_ELT(out,6,allele1_s);
  SET_VECTOR_ELT(out,7,allele2_s);
  SET_VECTOR_ELT(out,8,chrom_levels_s);

  {
    const char *names[9] = {
      "CHROM", "POS", "ID", "REF", "ALT", "individual",
      "allele1", "allele2", "CHROM.levels"
    };

    for (int i = 0; i < 9; i++)
      SET_STRING_ELT(names_s,i,mkChar(names[i]));
  }

  setAttrib(out,R_NamesSymbol,names_s);
  setAttrib(out,R_ClassSymbol,class_s);

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
      (NA_STRING == STRING_ELT(dataset,0)))
    return ScalarInteger(-2);

  gv = gs_fdta((char *)CHAR(STRING_ELT(dataset,0)));
  gs_tmp_init_from(&tmp,gv);

  if ((0 == gv->NoInd) || (0 == gv->NoMar) ||
      (NULL == gv->IndName) || (NULL == gv->MarName) ||
      (NULL == gv->All1) || (NULL == gv->All2))
    {
      gs_info(gv,-2,"No marker data available, run marker-data input first");
      failed = 1;
      goto cleanup;
    }

  if ((STRSXP != TYPEOF(id_s)) || (REALSXP != TYPEOF(y_s)) ||
      (0 == XLENGTH(id_s)) || (XLENGTH(id_s) != XLENGTH(y_s)))
    {
      gs_info(gv,-2,"Invalid performance data");
      failed = 1;
      goto cleanup;
    }

  if (0 != gs_make_name_index(gv->IndName,gv->NoInd,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,gv->NoInd))
    {
      gs_info(gv,-2,"Duplicate individual name in marker data");
      failed = 1;
      goto cleanup;
    }

  values = gs_malloc_array((size_t)gv->NoInd,sizeof(*values));
  seen = calloc((size_t)gv->NoInd,sizeof(*seen));
  if ((NULL == values) || (NULL == seen))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (R_xlen_t j = 0; j < XLENGTH(id_s); j++)
    {
      SEXP id = STRING_ELT(id_s,j);
      gs_UNSIGNED i;
      double value;

      if ((NA_STRING == id) ||
          !st_stvcf_valid_name(CHAR(id)) ||
          (strlen(CHAR(id)) >= gs_MAXNME))
        {
          gs_info(gv,-2,"Invalid individual ID in performance data");
          failed = 1;
          goto cleanup;
        }

      i = gs_name_index_find(index,gv->NoInd,CHAR(id));
      if (gs_UNSIGNED_MAX == i)
        continue;

      value = REAL(y_s)[j];
      if (!isfinite(value))
        {
          gs_info(gv,-2,"Invalid performance value");
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
      gs_info(gv,-2,"No performance data match marker-data individuals");
      failed = 1;
      goto cleanup;
    }

  tmp.NoInd = no_seen;
  tmp.NoMar = gv->NoMar;
  tmp.MaxIL = gv->MaxIL;
  tmp.MaxML = gv->MaxML;
  tmp.IndName = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.IndName));
  tmp.MarName = gs_malloc_array((size_t)tmp.NoMar,sizeof(*tmp.MarName));
  tmp.y = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.y));

  if ((NULL == tmp.IndName) || (NULL == tmp.MarName) ||
      (NULL == tmp.y) ||
      (0 != gs_alloc_genotype_matrix(&tmp,tmp.NoMar,tmp.NoInd)))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (gs_UNSIGNED m = 0; m < tmp.NoMar; m++)
    strcpy(tmp.MarName[m],gv->MarName[m]);

  {
    gs_UNSIGNED k = 0;

    for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
      if (seen[i])
        {
          strcpy(tmp.IndName[k],gv->IndName[i]);
          tmp.y[k] = values[i];

          for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            {
              size_t old_idx = (size_t)m * (size_t)gv->NoInd + i;
              size_t new_idx = (size_t)m * (size_t)tmp.NoInd + k;
              tmp.All1[new_idx] = gv->All1[old_idx];
              tmp.All2[new_idx] = gv->All2[old_idx];
            }
          k++;
        }
  }

  {
    int map_rc = gs_copy_map_primary(&tmp,gv);
    if (0 != map_rc)
      {
        gs_info(gv,-2,
                (-1 == map_rc) ?
                "Memory allocation error" :
                "Invalid linkage map state");
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

    gs_marker_stats_01(&tmp,
                       empty_p,empty_p,empty_p,
                       zero_p,zero_p,zero_p,zero_p,
                       &srv);
    if (-2 == srv)
      {
        failed = 1;
        goto cleanup;
      }
  }

  snprintf(gv->msg,sizeof(gv->msg),
           "No. of performance data: %lu",tmp.NoInd);
  gs_info(gv,0,gv->msg);

  gs_tmp_commit(gv,&tmp);
  status = 0;

 cleanup:
  free(index);
  free(values);
  free(seen);

  if (failed)
    gs_reset(&tmp);

  return ScalarInteger(status);
}


void gs_get_map_01 ( gs_varset_type *gv ,
		     char**  out_filename,
		     int*    retval        )
{
  gs_UNSIGNED  lne = 0;
  FILE *outfp;

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;
  
  for (lne=0; lne <  gv->NLoc; lne++)
    if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
      {
	snprintf(gv->msg, sizeof(gv->msg),
		 "%i %15.12f %s %s\n",
		 (gv->map[lne])->chrom,
		 (gv->map[lne])->pos,
		 (gv->map[lne])->name,
		 (gv->map[lne])->class);
	fprintf(outfp,"%s",gv->msg);
      }

  fclose(outfp);
  *retval = 0; return ;
}

void gs_get_map_01_GV (	 char**  out_filename,
			 int*    retval,
			 char**  set_name)
{

  gs_get_map_01 ( gs_fdta(*set_name),
		  out_filename,
		  retval         ) ;
}



int minimum (int a, int b) { return (a < b) ? a : b; }

void gs_str_rem_us (char* strg) 
{
  gs_UNSIGNED i, ii = strlen(strg);
  for (i = 0; i< ii; i++) 
    {
      if (
	  ('.' == strg[i]) ||
	  ('-' == strg[i]) ||
	  ('/' == strg[i]) ||
	  ('_' == strg[i]) 
	  )
	strg[i] = 'x';
    }


}

void gs_recode_ref_01 ( gs_varset_type* gv , 
			int* reference,
			int* descending,
			int* retval        )
/* Recoding to a reference indidivudal on individual marker basis */
{
  gs_UNSIGNED  i, j, idx, idx_R;
  int a,b,c,d;

  if ( ( *reference < 1 ) ||
       ( (gs_UNSIGNED) *reference > gv->NoInd ) ) {
      gs_info(gv,-2, "Invalid reference individual");
      *retval = -2; return ;
  }

  GS_ALLELE_TYPE *Ref1 =
    malloc(gv->NoMar * gv->NoInd * sizeof(*Ref1));
  GS_ALLELE_TYPE *Ref2 =
    malloc(gv->NoMar * gv->NoInd * sizeof(*Ref2));
  if ( (NULL == Ref1)  || (NULL == Ref2) ) ERR_M;
  
  gs_UNSIGNED R = (gs_UNSIGNED) *reference - 1;
      
  // Recode alleles in local memory
  for (i=0; i < gv->NoMar; i++)
    {
      idx_R = i*gv->NoInd + R;
      a = gv->All1[idx_R];
      b = gv->All2[idx_R];
      
      for (j=0; j <gv->NoInd; j++)
	{
	  idx = i*gv->NoInd + j;
	  c = gv->All1[idx];
	  d = gv->All2[idx];

          if ((GS_ALLELE_MISSING == a) ||
              (GS_ALLELE_MISSING == b) ||
              (GS_ALLELE_MISSING == c) ||
              (GS_ALLELE_MISSING == d))
            {
              Ref1[idx] = GS_ALLELE_MISSING;
              Ref2[idx] = GS_ALLELE_MISSING;
            } 
	  else
	    {
	      if (0==*descending) {
		if (( (c==a)&&(d==b)) ||
		    ( (c==b)&&(d==a))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 1;
		  } 
		else if (
			 ( (c==a)||(c==b)) ||
			 ( (d==a)||(d==b))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 2;
		  } 
		else
		  {
		    Ref1[idx] = 2;
		    Ref2[idx] = 2;
		  }
	      } else {
		if (( (c==a)||(c==b)) &&
		    ( (d==a)||(d==b))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 1;
		  } 
		else if (
			 ( (c==a)||(c==b)) ||
			 ( (d==a)||(d==b))    ) 
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
  
  // Write back the recoded alleles
  free(gv->All1); gv->All1 = Ref1;
  free(gv->All2); gv->All2 = Ref2;

  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv;
  char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return;

}

void gs_recode_ref_01_GV ( int* reference,
			   int* descending,
			   int* retval,
			   char**  set_name)
{
  gs_recode_ref_01 ( gs_fdta(*set_name),
		     reference,
		     descending,
		     retval            ) ;

}

void gs_recode_ref_name_01_GV ( char** reference_name,
                                int* descending,
                                int* retval,
                                char** set_name )
{
  gs_varset_type *gv = gs_fdta(*set_name);
  gs_UNSIGNED i;
  int reference = 0;

  if ( (NULL == reference_name) ||
       (NULL == *reference_name) ||
       ('\0' == (*reference_name)[0]) )
    {
      gs_info(gv,-2,"Invalid reference individual");
      *retval = -2;
      return;
    }

  for (i=0; i<gv->NoInd; i++)
    if (0 == strcmp(gv->IndName[i],*reference_name))
      {
        reference = (int) i + 1;
        break;
      }

  if (0 == reference)
    {
      snprintf(gv->msg,sizeof(gv->msg),
               "Reference individual not found: %s",
               *reference_name);
      gs_info(gv,-2,gv->msg);
      *retval = -2;
      return;
    }

  gs_recode_ref_01(gv,&reference,descending,retval);
}

gs_FLOAT gs_avg_ld ( gs_chrom_type * chrom,
		     gs_UNSIGNED mm_lft   ,
		     gs_UNSIGNED mm_rgt     )
{
  gs_FLOAT ld=0;
  gs_FLOAT n ;
  gs_UNSIGNED idx,i,j;

  for ( i=mm_lft; i < mm_rgt; i++)
    for ( j=i+1   ; j <= mm_rgt; j++)
    {
      idx = i * chrom->NLocC + j;
      ld += chrom->ld [idx];
    }
  n = ( (gs_FLOAT)(mm_rgt-mm_lft) +1 ) * ( (gs_FLOAT)mm_rgt-mm_lft)    / 2;

  return ( ld /n );

}

gs_FLOAT gs_check_extend ( gs_chrom_type * chrom,
			   gs_UNSIGNED mm_lft,
			   gs_UNSIGNED mm_rgt,
			   gs_UNSIGNED tol_lft,
			   gs_UNSIGNED tol_rgt,
			   char** ld_criterion
			  )
{

  gs_UNSIGNED  t,idx;
  gs_FLOAT tld,ld=0;

  if ( 0 == strcmp ( "flanking_fast" , *ld_criterion) )  
    {
      idx = mm_lft * chrom->NLocC + mm_rgt;
      ld = chrom->ld [idx];
    }

  if ( 0 == strcmp ( "flanking" , *ld_criterion) )  
    {
      // ld between the markers mm_lft and mm_rgt
      idx = mm_lft * chrom->NLocC + mm_rgt;
      ld = chrom->ld [idx];
      // Tolerance on the left is defined
      if ( tol_lft > 0 )
	{
	  if (tol_lft > mm_lft) tol_lft = mm_lft; // Left end of the chromosome
	  for ( t=1; t<=tol_lft; t++)
	  {
	    idx = (mm_lft-t) * chrom->NLocC + mm_rgt;  // LD between flanking is 
	    if (ld < chrom->ld [idx])            // replaced between the largest
	      ld = chrom->ld [idx];              // LD in the tolerance region
	  }
	}
      // Tolerance on the right is defined
      if ( tol_rgt > 0 )
	{
	  if ( tol_rgt > ( (chrom->NLocC)-1-mm_rgt) )  
	    tol_rgt = ((chrom->NLocC)-1-mm_rgt);
	for (t=1;t<=tol_rgt;t++ )
	  {
	    idx = mm_lft * chrom->NLocC + (mm_rgt+t);  // LD between flanking is 
	    if (ld < chrom->ld [idx])            // replaced between the largest
	      ld = chrom->ld [idx];              // LD in the tolerance regionn
	  }
	}
    }

  if  ( 0 == strcmp ( "average" ,  *ld_criterion) )
    {
      // ld between the markers mm_lft and mm_rgt
      ld = gs_avg_ld ( chrom, mm_lft, mm_rgt);

      // Tolerance on the left is defined
      if ( tol_lft > 0 )
	{
	  if (tol_lft > mm_lft) tol_lft = mm_lft; // Left end of the chromosome
	  for ( t=1; t<=tol_lft; t++)
	  {
	    tld = gs_avg_ld ( chrom, mm_lft-t, mm_rgt);
	    if (ld < tld) ld = tld;
	  }
	}
      // Tolerance on the right is defined
      if ( tol_rgt > 0 )
	{
	  if ( tol_rgt > ( (chrom->NLocC)-1-mm_rgt) )  
	    tol_rgt = ((chrom->NLocC)-1-mm_rgt);
	for (t=1;t<=tol_rgt;t++ )
	  {
	    tld = gs_avg_ld ( chrom, mm_lft, mm_rgt+t);
	    if (ld < tld) ld = tld;
	  }
	}
    }

  return(ld);

}

static int gs_build_haplotype_map(gs_varset_type *gv,
                                  char **hap_symbol)
{
  gs_mappoint_type *Hap_S = NULL;
  char *Hap_str_S = NULL;
  size_t Hap_str_N = 0;
  size_t offset = 0;

  if ((0 == gv->NHap) || (NULL == hap_symbol) ||
      (NULL == *hap_symbol))
    return -2;

  for (gs_UNSIGNED l = 0; l < gv->NHap; l++)
    {
      char name[gs_MAXNME];
      int n = snprintf(name,sizeof(name),"%1s%06lu",
                       *hap_symbol,
                       (unsigned long)l);
      size_t needed;

      if ((n < 0) || ((size_t)n >= sizeof(name)) ||
          ((size_t)n >= gs_MAXSTR))
        return -2;

      needed = (size_t)n + 1;
      if ((0 != gs_size_add(Hap_str_N,needed,&Hap_str_N)) ||
          (0 != gs_size_add(Hap_str_N,2,&Hap_str_N)))
        return -1;
    }

  Hap_S = gs_malloc_array((size_t)gv->NHap,sizeof(*Hap_S));
  Hap_str_S = gs_malloc_array(Hap_str_N,sizeof(*Hap_str_S));
  if ((NULL == Hap_S) || (NULL == Hap_str_S))
    {
      free(Hap_S);
      free(Hap_str_S);
      return -1;
    }

  for (gs_UNSIGNED l = 0; l < gv->NHap; l++)
    {
      char name[gs_MAXNME];
      int n = snprintf(name,sizeof(name),"%1s%06lu",
                       *hap_symbol,
                       (unsigned long)l);

      if ((n < 0) || ((size_t)n >= sizeof(name)))
        {
          free(Hap_S);
          free(Hap_str_S);
          return -2;
        }

      Hap_S[l].chrom = gv->map[gv->h_u[l]]->chrom;
      Hap_S[l].pos = (gv->map[gv->h_o[l]]->pos +
                      gv->map[gv->h_u[l]]->pos) / 2;
      Hap_S[l].nr = l;

      if ((0 != gs_pool_put_string(Hap_str_S,Hap_str_N,&offset,
                                   name,&(Hap_S[l].name))) ||
          (0 != gs_pool_put_string(Hap_str_S,Hap_str_N,&offset,
                                   "b",&(Hap_S[l].class))))
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

void gs_def_hblocks_01 ( gs_varset_type* gv , 
			 int* hap,
			 int* hap_unit,
			 double* ld_threshold,
			 char**  ld_criterion,
			 int* tolerance,
			 char** hap_symbol,
			 char** out_filename,
			 int*   auxfiles,
			 int* retval        )
/* Determining haplotype borders */
{

  if ( 0 == gv->NLoc ) 
    {
      snprintf(gv->msg, sizeof(gv->msg),
         "No linkage map has been loaded");
      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }
  
// gs_UNSIGNED i,l,u,o,n, fib, lib;

  gs_UNSIGNED i,l,o,n, fib, lib;

  gv->h_u   = realloc ( gv->h_u , gv->NLoc*sizeof(gs_UNSIGNED) )  ;
  gv->h_o   = realloc ( gv->h_o , gv->NLoc*sizeof(gs_UNSIGNED) )  ;
  gv->h_n   = realloc ( gv->h_n , gv->NLoc*sizeof(gs_UNSIGNED) )  ;
  if ( (NULL == gv->h_u)||(NULL == gv->h_o)||(gv->h_n == 0) ) ERR_M;

  gv->NHap  = 0;                    
  o = 0; // upper block border
  n = 0; // number of (observed) markers
//  u = 0;
  fib = gs_UNSIGNED_MAX; // first in block
  lib = gs_UNSIGNED_MAX; // last in block

if ( ( 0 == strcmp ( "average" ,  *ld_criterion) ) ||
     ( 0 == strcmp ( "flanking" , *ld_criterion) )    )
    {
      // Messages
      
      if ( 0 == strcmp ( "average" ,  *ld_criterion) )
	snprintf(gv->msg, sizeof(gv->msg),
		 "Average LD in haplotype block > %4.2f",
		 *ld_threshold);
      if ( 0 == strcmp ( "flanking" ,  *ld_criterion) )
	snprintf(gv->msg, sizeof(gv->msg),
		 "LD of markers flanking haplotype blocks > %4.2f",
		 *ld_threshold);
      gs_info(gv,0,gv->msg);

   
      if ( (NULL == gv->chrom) || ( NULL == gv->chrom[0].ld ) )
	{
	  snprintf(gv->msg, sizeof(gv->msg),
		   "No LD values available");
	  gs_info(gv,-2,gv->msg);
	  *retval = -2; 
	  return;
	}

      gs_UNSIGNED tol;
      if (*tolerance<0)
	{
	  snprintf(gv->msg, sizeof(gv->msg),
         "Incorrect tolerance value");
	  gs_info(gv,-2,gv->msg);
	  *retval = -2; 
	  return;
	}

      tol = (gs_UNSIGNED) *tolerance; 

      // Assign each locus to block

      gs_UNSIGNED *blk = malloc (gv->NLoc * sizeof(gs_UNSIGNED)); // Block nr
      if (NULL == blk)  ERR_M;                                    // for each 
      for (i=0; i <  gv->NLoc; i++) blk[i] = 0;                   // locus

      gs_UNSIGNED c,mm,idx,lfdnr;
      gs_UNSIGNED locus1,locus2;
      gs_FLOAT ld;
      int unasigned_loci, block_extendable;
      gs_UNSIGNED greatest_ld_locus;
      gs_FLOAT    greatest_ld_value,extend_lft,extend_rgt;
      gs_UNSIGNED mm_rgt, mm_lft, loc_rgt, loc_lft ;

      lfdnr = 0; // Running number for blocks

      for ( c=0 ; c < gv->NChrom ; c ++ ) { // Go through all chromosomes

	do {                                // While there are unasigned loci  

	  // Search for the largest LD value between two adjacent loci 
	  // which were not yet assigned to a block

	  greatest_ld_locus = 0; // right hand locus of pair with geatest ld
	  greatest_ld_value = 0;

	  // Find starting locus pair for the next block
	  for ( mm=1 ; mm < gv->chrom[c].NLocC ; mm++) { 
                                                     
	    locus2 = gv->chrom[c].first_loc_i + mm ;   // Positions map index
	    locus1 = locus2 - 1;
	    
	    idx = (mm-1) * gv->chrom[c].NLocC + mm;     // LD
	    ld = (gv->chrom[c].ld) [idx];

	    if (( 0 == blk[locus1] ) && (0 == blk[locus2]) && // Unassigned and
		ld > *ld_threshold)                           // in LD
	      {
		if ( ld>greatest_ld_value )
		  {
		    greatest_ld_locus = mm ;
		    greatest_ld_value = ld ;
		  }
	      }
	  }

	  if (0 ==greatest_ld_locus)   
	    // No locus pair in LD: each unassingned marker is a block 
	    {                       
	      for ( mm=0 ; mm < gv->chrom[c].NLocC ; mm++) 
		 if ( 0 == blk[ gv->chrom[c].first_loc_i + mm] )
		   blk [gv->chrom[c].first_loc_i + mm] = ++lfdnr;
	    }
	  else // Locus pair in LD is found: build a new block                     
	    {

	      mm_rgt = greatest_ld_locus;        // Index on the chromosome
	      mm_lft = greatest_ld_locus - 1 ;

	      loc_rgt = gv->chrom[c].first_loc_i + mm_rgt ;  // Index on the map
	      loc_lft = gv->chrom[c].first_loc_i + mm_lft ;  

	      blk [loc_rgt] = ++ lfdnr;  // Block consists at least of the
	      blk [loc_lft] =    lfdnr;  // found locus pair in LD

	      // Check whether the block can be extended
	      do {
		block_extendable = 0;

		if ( mm_lft > 0 )
		  extend_lft = gs_check_extend ( &(gv->chrom[c]), 
						 mm_lft-1, mm_rgt,
						 tol, 0,
						 ld_criterion );
		else 
		  extend_lft = 0;

		if ( mm_rgt < (( gv->chrom[c].NLocC) - 1) )
		  extend_rgt = gs_check_extend ( &(gv->chrom[c]), 
						 mm_lft, mm_rgt+1,
						 0, tol,
						 ld_criterion );
		else 
		  extend_rgt = 0;

		if ( ( extend_lft >= extend_rgt ) && ( extend_lft >= *ld_threshold ) )
		  // extension to the left
		  {
		    mm_lft = mm_lft - 1;
		    loc_lft = gv->chrom[c].first_loc_i + mm_lft ;  
		    blk [loc_lft] = lfdnr;
		    block_extendable = 1;
		  }
		else if ( (extend_rgt>extend_lft ) && ( extend_rgt >= *ld_threshold) )
		  // extension to the right
		  {
		    mm_rgt = mm_rgt +1;
		    loc_rgt = gv->chrom[c].first_loc_i + mm_rgt ; 
		    blk [loc_rgt] = lfdnr;
		    block_extendable = 1;
		  }

	      } while (1==block_extendable) ;
	      
	    }

	  // Are there unassigned still loci?
	  unasigned_loci = 0;
	  for ( mm=0 ; mm < gv->chrom[c].NLocC ; mm++) 
	    if ( 0 == blk[ gv->chrom[c].first_loc_i + mm] ) unasigned_loci=1;

	} while (1 == unasigned_loci);

      } // for c = 1 .. NChrom


      
      // Make a list of haplotype blocks from the assignment list of
      // markers to blocks

      gv->NHap = 0;
      fib=0;

      for ( i=1; i<gv->NLoc; i++)
	{
	  if  (blk[i-1] != blk[i])
	    {
	      gv->h_u [gv->NHap] = fib;
	      gv->h_o [gv->NHap] = i-1;
	      gv->h_n [gv->NHap] = i-fib;
	      gv->NHap++;
	      fib=i;
	    }
	}
      gv->h_u [gv->NHap] = fib;
      gv->h_o [gv->NHap] = i-1;
      gv->h_n [gv->NHap] = i-fib;
      gv->NHap++;

      free (blk);

    }

    else if (0 == *hap_unit) // disjunct segments
    {
      
      snprintf(gv->msg, sizeof(gv->msg),
         "Haplotype borders form disjunct LD segments ");
      gs_info(gv,0,gv->msg);

      gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_UNSIGNED i,m,cc,l1,l2,idx1,idx2;
      for ( cc=0 ; cc < gv->NChrom ; cc ++ ) {
	fib = gv->chrom[cc].first_loc_i;

	for ( m=1 ; m < gv->chrom[cc].NLocC ; m++) {
	  n++;
	  for ( i=0; i < gv->NoInd; i++) {
	    l2 = gv->chrom[cc].first_loc_i + m ;
	    l1 = l2-1;
	    
	    idx1 = l1*gv->NoInd + i;
	    idx2 = l2*gv->NoInd + i;
	    
	    if (  (gv->All1[idx1] != gv->All1[idx2]) ||
		  (gv->All2[idx1] != gv->All2[idx2])    )
	      {
		lib = fib + n -1;
		gv->h_u [gv->NHap] = fib;
		gv->h_o [gv->NHap] = lib;
		gv->h_n [gv->NHap] = n;
		gv->NHap++;
		n = 0; 
		fib = lib + 1;
		break;
	      }
	  }
	}

	lib = fib + n;
	gv->h_u [gv->NHap] = fib;
	gv->h_o [gv->NHap] = lib;
	gv->h_n [gv->NHap] = n;
	gv->NHap++;
	n = 0; 
      }
    }

  else if (1 == *hap_unit) // Marker count
    {
      snprintf(gv->msg, sizeof(gv->msg),
         "Haplotyp borders every %i markers",
         *hap);
      gs_info(gv,0,gv->msg);

      do {
      // find blocks
	do {
	  if ( gv->mdl[o] < gs_UNSIGNED_MAX ) {
	    if  ( gs_UNSIGNED_MAX == fib )  fib = o;
	    lib = o;
	    n++;
	  }
	  o++; 
	}  while  (  ( o < gv->NLoc) &&                         
		     !( (int)n > (*hap-1) ) &&
		     !( gv->map[o-1]->chrom !=  gv->map[o]->chrom ) );
	
	if (n>0) {
	  gv->h_u [gv->NHap] = fib;
	  gv->h_o [gv->NHap] = lib;
	  gv->h_n [gv->NHap] = n;
	  gv->NHap++;
	  
	  n = 0; 
	  fib = gs_UNSIGNED_MAX; 
	  lib = gs_UNSIGNED_MAX;
	}
	
      } while  ( o < (gv->NLoc) ) ;
    }

 
  else if (2 == *hap_unit) // Segment length
    {

      snprintf(gv->msg, sizeof(gv->msg),
         "Haplotyp borders every %i centiMorgans",
         *hap);
      
      gs_info(gv,0,gv->msg);

      do {
	// find blocks
	do {
	  if ( gv->mdl[o] < gs_UNSIGNED_MAX ) {
	    if  ( gs_UNSIGNED_MAX == fib )  fib = o;
	    lib = o;
	    n++;
	  }
	  o++; 
	}  while  (  ( o < gv->NLoc) &&                         
		     ! ( (double)*hap < 
			 ( gv->map[lib]->pos -  gv->map[fib]->pos ) ) &&
		     !( gv->map[o-1]->chrom !=  gv->map[o]->chrom ) );
	
	if (n>0) {
	  gv->h_u [gv->NHap] = fib;
	  gv->h_o [gv->NHap] = lib;
	  gv->h_n [gv->NHap] = n;
	  gv->NHap++;

	  n = 0; 
	  fib = gs_UNSIGNED_MAX; 
	  lib = gs_UNSIGNED_MAX;
	}
      } while  ( o < (gv->NLoc) ) ;
    }

  else
    // No valid specification of haplotyping algorithm
    {

      snprintf(gv->msg, sizeof(gv->msg),
         "Provide either hap.unit or ld.criterion");

      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }

  // Determine haploblock map
  {
    int map_rc = gs_build_haplotype_map(gv,hap_symbol);

    if (0 != map_rc)
      {
        gv->NHap = 0;
        GS_FREE(gv->Hap_S);
        GS_FREE(gv->Hap_str_S);
        gv->Hap_str_N = 0;

        if (-1 == map_rc)
          ERR_M;

        gs_info(gv,-2,"Haplotype marker name too long");
        *retval = -2;
        return;
      }
  }
    
  if (1 == *auxfiles){   
    // Write individual names to .hap file
    // New Version 2017-08-31

    FILE *outfp;
    if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

    for ( l=0; l < gv->NHap; l++)
      {
	char* lne = (char *) malloc ( gs_MAXLNE * sizeof(char) );

	snprintf(lne, gs_MAXLNE,
		 "%i %f %s %s ",
		 gv->Hap_S[l].chrom,
		 gv->Hap_S[l].pos,
		 gv->Hap_S[l].name,
		 gv->Hap_S[l].class);

	for ( i = gv->h_u[l] ; i <= gv->h_o[l] ; i++ ) 
	  {
	    if ( gv->mdl[i] < gs_UNSIGNED_MAX )
	      {
		size_t a=strlen(lne);
		size_t b=strlen(gv->MarName[ gv->mdl[i]]);
		lne = realloc (lne, 2+a+b );
		strcat( lne, gv->MarName[ gv->mdl[i] ] );
		strcat( lne, ";" );
	      }
	  }

	fprintf(outfp,"%s\n",lne);
	free (lne);
      }
    
    fclose(outfp);
  }
  
}

void gs_def_hblocks_01_GV ( int*   hap,
			    int*   hap_unit,
			    double* ld_threshold,
			    char**  ld_criterion,
			    int* tolerance,
			    char** hap_symbol,
			    char** out_filename,
			    int*   auxfiles,    
			    int*   retval,
			    char** set_name)
{
  gs_def_hblocks_01 ( gs_fdta(*set_name),
		      hap,
		      hap_unit,
		      ld_threshold,
		      ld_criterion,
		      tolerance,
		      hap_symbol,
		      out_filename,
		      auxfiles,
		      retval            ) ;
}



int gs_swap_block_marker (gs_varset_type* gv )
{

  gs_UNSIGNED i;

  // Marker (=haploblock) data
  gs_UNSIGNED NoInd = gv->NoInd;
  gs_UNSIGNED NoMar = gv->NHap;

  char (*MarName) [gs_MAXNME] =
    gs_malloc_array((size_t)NoMar,sizeof(*MarName));
  gs_mappoint_type **map =
    gs_malloc_array((size_t)NoMar,sizeof(*map));
  gs_UNSIGNED *mdl =
    gs_malloc_array((size_t)NoMar,sizeof(*mdl));

  if ((NULL == MarName) || (NULL == map) || (NULL == mdl))
    {
      free(MarName);
      free(map);
      free(mdl);
      return -1;
    }

  for (i = 0; i < gv->NHap; i++)
    {
      strcpy(MarName[i],gv->Hap_S[i].name);
      map[i] = &(gv->Hap_S[i]);
      mdl[i] = i;
    }

  char (*IndName) [gs_MAXNME] = gv->IndName; gv->IndName = NULL;

  GS_ALLELE_TYPE *All1 = gv->Ref1; gv->Ref1 = NULL;
  GS_ALLELE_TYPE *All2 = gv->Ref2; gv->Ref2 = NULL;

  // Map data
  gs_UNSIGNED NLoc = gv->NHap;
  gs_mappoint_type *map_S = gv->Hap_S; gv->Hap_S = NULL;
  char *map_str_S = gv->Hap_str_S; gv->Hap_str_S = NULL;
  size_t map_str_N = gv->Hap_str_N; gv->Hap_str_N = 0;

  // Phenotypes
  double *y = gv->y;
  gv->y = NULL;       // prevent from being freed

  gs_reset(gv);

  gv->NoInd = NoInd;
  gv->NoMar = NoMar;

  gv->IndName = IndName;
  gv->MarName = MarName;

  gv->All1 = All1;
  gv->All2 = All2;

  gv->NLoc      = NLoc;
  gv->map_S     = map_S;
  gv->map_str_S = map_str_S;
  gv->map_str_N = map_str_N;
  gv->map       = map;
  gv->mdl       = mdl;

  gv->y = y;

  return 0;
}

void gs_recode_hbc_01 ( gs_varset_type* gv , 
			int* reference,
			int* retval        )
/* Recoding to a haplotype blocks for BC individuals */
{

  if ( 0 == gv->NHap ) 
    {
	snprintf(gv->msg, sizeof(gv->msg), 
	       "No haplotype blocks were defined");
      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }

  gs_UNSIGNED  l, i, j, idx, idx_R;

  if ( (*reference < 1) || ((gs_UNSIGNED)*reference > gv->NoInd) )
    {
      snprintf(gv->msg, sizeof(gv->msg),
               "Invalid reference individual");
      gs_info(gv,-2,gv->msg);
      *retval = -2;
      return;
    }

  // Allocate memory

  gv->Ref1 = realloc(gv->Ref1,gv->NoMar * gv->NoInd * sizeof(*gv->Ref1));
  gv->Ref2 = realloc(gv->Ref2,gv->NoMar * gv->NoInd * sizeof(*gv->Ref2));
  if ( (NULL == gv->Ref1)  || (NULL == gv->Ref2) ) ERR_M;

  // Determine hayplotypes
  gs_UNSIGNED R = (gs_UNSIGNED) *reference - 1;
  gs_UNSIGNED idx_H, n0,n1,n2;
  int a,b,c,d;

  (void)n2;
  
  for ( l = 0 ; l < gv->NHap ; l++ ) {

    for (j = 0; j < gv->NoInd ; j ++)
      {
	n0 = 0; n1 = 0; n2 = 0; 
	idx_H = l*gv->NoInd + j;
	
	for ( i = gv->h_u[l] ; i <= gv->h_o[l] ; i++ ) 
	  if ( gv->mdl[i] < gs_UNSIGNED_MAX )
	    {
	      
	      idx_R = i*gv->NoInd + R;
	      a = gv->All1[idx_R];
	      b = gv->All2[idx_R];
	      
	      idx = i*gv->NoInd + j;
	      c = gv->All1[idx];
	      d = gv->All2[idx];
	      
              if (!(                             // not missing
                    (GS_ALLELE_MISSING == a) ||
                    (GS_ALLELE_MISSING == b) ||
                    (GS_ALLELE_MISSING == c) ||
                    (GS_ALLELE_MISSING == d)
                                                    )) 
		{
		  if (                          // two common alleles
		      ( (c==a)||(c==b)) &&
		      ( (d==a)||(d==b))    ) 
		    {
		      n2++ ;
		    } 
		  else if (                     // one common allele
			   ( (c==a)||(c==b)) ||
			   ( (d==a)||(d==b))    ) 
		    {
		      n1++ ;
		    } 
		  else
		    {
		      n0++;
		    }
		} 
 	    }
	
	gv->Ref1[idx_H] =  1; // Can originate from the RP
	gv->Ref2[idx_H] =  1; // 

	if (n0>0)                // At least in one marker in the block
	  {                      // both alleles can not 
	    gv->Ref1[idx_H] = 2; // originate from 
	    gv->Ref2[idx_H] = 2; // the RP
	  }  
	else if  (n1>0)  
	  {
	    gv->Ref1[idx_H] = 1;
	    gv->Ref2[idx_H] = 2;
	  }
 
      }
  }

  /* Use the haplotype data as markers */
  if (0 != gs_swap_block_marker(gv)) ERR_M;

  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv; char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }
  
  *retval = 0; return;
}


void gs_recode_hbc_01_GV ( int*   reference,
			   int*   retval,
			   char** set_name)
{
  gs_recode_hbc_01 ( gs_fdta(*set_name),
		     reference,
		     retval            ) ;

}

/* qsort haplotype alleles compariosn*/
int gs_hall_cmp(const void *a, const void *b) 
{ 

  const gs_usp *ia = (const gs_usp *)a;
  const gs_usp *ib = (const gs_usp *)b;

  if ( (ia)->n < (ib)->n ) return (1);
  if ( (ia)->n > (ib)->n ) return (-1);
  return (0);
} 

void gs_recode_hil_01 ( gs_varset_type* gv , 
                        char** out_filename,
                        int* auxfiles,    
                        int* retval        )
/* Recoding to a haplotype blocks for inbred lines */
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
      snprintf(gv->msg,sizeof(gv->msg),"No haplotype blocks were defined");
      gs_info(gv,-2,gv->msg);
      *retval = -2;
      return;
    }

  if ((0 != gs_size_mul((size_t)gv->NoMar,
                        (size_t)gv->NoInd,
                        &genotype_count)) ||
      ((size_t)gv->NoInd > SIZE_MAX / 2))
    {
      gs_info(gv,-2,"Memory allocation error");
      *retval = -2;
      return;
    }
  max_hall = 2 * (size_t)gv->NoInd;

  gv->Ref1 = realloc(gv->Ref1,genotype_count * sizeof(*gv->Ref1));
  gv->Ref2 = realloc(gv->Ref2,genotype_count * sizeof(*gv->Ref2));
  AT1 = gs_malloc_array(genotype_count,sizeof(*AT1));
  AT2 = gs_malloc_array(genotype_count,sizeof(*AT2));
  Hall = gs_malloc_array(max_hall,sizeof(*Hall));
  OcHall = gs_malloc_array(max_hall,sizeof(*OcHall));
  Hidx1 = gs_malloc_array((size_t)gv->NoInd,sizeof(*Hidx1));
  Hidx2 = gs_malloc_array((size_t)gv->NoInd,sizeof(*Hidx2));

  if ((NULL == gv->Ref1) || (NULL == gv->Ref2) ||
      (NULL == AT1) || (NULL == AT2) ||
      (NULL == Hall) || (NULL == OcHall) ||
      (NULL == Hidx1) || (NULL == Hidx2))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  /* Transpose genotypic data. */
  for (gs_UNSIGNED i=0; i < gv->NoMar; i++)
    for (gs_UNSIGNED j=0; j < gv->NoInd; j++)
      {
        size_t idx = (size_t)i * (size_t)gv->NoInd + j;
        size_t idx_T = (size_t)j * (size_t)gv->NoMar + i;
        AT1[idx_T] = gv->All1[idx];
        AT2[idx_T] = gv->All2[idx];
      }

  for (gs_UNSIGNED l = 0; l < gv->NHap; l++)
    {
      gs_UNSIGNED i = gv->h_u[l];
      gs_UNSIGNED NoHall = 0;
      size_t n = sizeof(GS_ALLELE_TYPE) * (size_t)gv->h_n[l];

      memset(OcHall,0,max_hall * sizeof(*OcHall));

      /* Preserve the original HIL counting/order behavior: initialize the
         first haplotype before iterating over all individuals. */
      Hall[NoHall] = &(AT1[(size_t)i]);
      OcHall[NoHall]++;
      NoHall++;

      for (gs_UNSIGNED j = 0; j < gv->NoInd; j++)
        {
          size_t idx_T = (size_t)j * (size_t)gv->NoMar + i;
          GS_ALLELE_TYPE *P1 = &(AT1[idx_T]);
          int found = 0;

          for (gs_UNSIGNED k = 0; k < NoHall; k++)
            if (0 == memcmp(Hall[k],P1,n))
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
                  gs_info(gv,-2,
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
            if (0 == memcmp(Hall[k],P1,n))
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
                  gs_info(gv,-2,
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

      ar = gs_realloc_array(ar,(size_t)NoHall,sizeof(*ar));
      if (NULL == ar)
        {
          gs_info(gv,-2,"Memory allocation error");
          failed = 1;
          goto cleanup;
        }

      for (gs_UNSIGNED k = 0; k < NoHall; k++)
        {
          ar[k].a = k;
          ar[k].n = OcHall[k];
        }
      qsort(ar,NoHall,sizeof(*ar),gs_hall_cmp);

      /* Recode the original zero-based temporary haplotype indices to the
         same frequency-ranked 1..NoHall codes as in the previous code. */
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

      /* Preserve the existing HIL behavior for haplotypes containing
         missing marker alleles: their final haplotype allele is missing. */
      for (gs_UNSIGNED j = 0; j < gv->NoInd; j++)
        {
          size_t idx_H = (size_t)l * (size_t)gv->NoInd + j;
          gs_UNSIGNED rank = (gs_UNSIGNED)gv->Ref1[idx_H] - 1;

          for (gs_UNSIGNED kk=0; kk < gv->h_n[l]; kk++)
            if (GS_ALLELE_MISSING == Hall[ar[rank].a][kk])
              {
                gv->Ref1[idx_H] = GS_ALLELE_MISSING;
                break;
              }

          rank = (gs_UNSIGNED)gv->Ref2[idx_H] - 1;
          for (gs_UNSIGNED kk=0; kk < gv->h_n[l]; kk++)
            if (GS_ALLELE_MISSING == Hall[ar[rank].a][kk])
              {
                gv->Ref2[idx_H] = GS_ALLELE_MISSING;
                break;
              }
        }

      /* Preserve the auxiliary haplotype list and its previous external
         -1 representation for missing marker alleles. */
      {
        char **new_hal = gs_realloc_array(hal,
                                          (size_t)(n_hal + NoHall),
                                          sizeof(*new_hal));
        if (NULL == new_hal)
          {
            gs_info(gv,-2,"Memory allocation error");
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
              gs_info(gv,-2,"Memory allocation error");
              failed = 1;
              goto cleanup;
            }

          snprintf(tmp_str,cap,"%s %lu ",gv->Hap_S[l].name,k+1);
          for (gs_UNSIGNED kk=0; kk < gv->h_n[l]; kk++)
            {
              int value = (GS_ALLELE_MISSING == Hall[ar[k].a][kk]) ?
                          -1 : (int)Hall[ar[k].a][kk];
              snprintf(tmp_str2,sizeof(tmp_str2),
                       (kk + 1 < gv->h_n[l]) ? "%i;" : "%i",value);
              strcat(tmp_str,tmp_str2);
            }
          hal[n_hal + k] = tmp_str;
        }

      n_hal += NoHall;
    }

  free(AT1); AT1 = NULL;
  free(AT2); AT2 = NULL;
  free(Hall); Hall = NULL;
  free(OcHall); OcHall = NULL;
  free(Hidx1); Hidx1 = NULL;
  free(Hidx2); Hidx2 = NULL;
  free(ar); ar = NULL;

  /* Use the haplotype data as markers. */
  if (0 != gs_swap_block_marker(gv))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  {
    int srv=0;
    char *empty="";
    char **empty_p=&empty;
    int zero=0;

    gs_marker_stats_01(gv,empty_p,empty_p,empty_p,
                       &zero,&zero,&zero,&zero,&srv);
    if (-2 == srv)
      {
        failed = 1;
        goto cleanup;
      }
  }

  if (1 == *auxfiles)
    {
      FILE *outfp = fopen(*out_filename,"w");
      if (NULL == outfp)
        {
          gs_info(gv,-2,"File open error");
          failed = 1;
          goto cleanup;
        }
      for (gs_UNSIGNED k=0; k<n_hal; k++)
        fprintf(outfp,"%s\n",hal[k]);
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
  for (gs_UNSIGNED k=0; k<n_hal; k++) free(hal[k]);
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

void gs_recode_hil_01_GV ( char** out_filename,
			   int*   auxfiles,    
			   int*   retval,
			   char** set_name)
{
  gs_recode_hil_01 ( gs_fdta(*set_name),
		     out_filename,
		     auxfiles,
		     retval            ) ;
}


void gs_write_marker_data_02 ( gs_varset_type *gv ,
                               char** type ,
                               char** lpo_filename,
                               char** mmp_filename,
                               char** mpo_filename,
                               char** ngp_filename,
                               char** nmp_filename,
                               char** npo_filename,
                               char** ind_filename,
                               int*   f_ind,
                               int*   l_ind,
                               char** add_inum,
                               int*   auxfiles,
                               int*   retval )
{
  gs_UNSIGNED *pli_idx = NULL;
  gs_UNSIGNED npli = 0;
  gs_UNSIGNED capacity;
  char (*n_names)[gs_MAXNME] = NULL;
  FILE *outfp = NULL;
  char line[gs_MAXSTR];
  char t;
  int failed = 0;

  (void)auxfiles;

  if ((NULL == gv->All1) || (NULL == gv->All2) ||
      (0 == gv->NoInd) || (0 == gv->NoMar))
    {
      gs_info(gv,-2,"No marker data available");
      *retval = -2;
      return;
    }

  t = (*type)[0];
  if (!(('l' == t) || ('a' == t) || ('m' == t) || ('n' == t)) ||
      ('\0' != (*type)[1]))
    {
      gs_info(gv,-2,"Implemented output formats: l, a, m, n");
      *retval = -2;
      return;
    }

  capacity = (0 == gv->NoInd) ? 1 : gv->NoInd;
  pli_idx = gs_malloc_array((size_t)capacity,sizeof(*pli_idx));
  if (NULL == pli_idx)
    {
      gs_info(gv,-2,"Memory allocation error");
      *retval = -2;
      return;
    }

  npli = gv->NoInd;
  for (gs_UNSIGNED j = 0; j < npli; j++) pli_idx[j] = j;

  if ((0 != *l_ind) || (0 != *f_ind))
    {
      long f = *f_ind;
      long l = *l_ind;

      if ((long)gv->NoInd < f) f = (long)gv->NoInd;
      if ((long)gv->NoInd < l) l = (long)gv->NoInd;
      if (0 < f) f--; else if (0 > f) f = 0;
      if ((1 > l) || ((long)gv->NoInd < l)) l = (long)gv->NoInd;

      npli = 0;
      for (long j = f; j < l; j++)
        pli_idx[npli++] = (gs_UNSIGNED)j;
    }

  if (0 != strcmp(*ind_filename,""))
    {
      FILE *fp = fopen(*ind_filename,"r");
      if (NULL == fp)
        {
          gs_info(gv,-2,"File open error");
          failed = 1;
          goto cleanup;
        }

      npli = 0;
      while (1 == fscanf(fp,"%4095s",line))
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
                        gs_info(gv,-2,"Too many output individuals");
                        failed = 1;
                        goto cleanup;
                      }

                    new_idx = gs_realloc_array(pli_idx,
                                               (size_t)new_capacity,
                                               sizeof(*new_idx));
                    if (NULL == new_idx)
                      {
                        fclose(fp);
                        gs_info(gv,-2,"Memory allocation error");
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
      gs_UNSIGNED *storage = gs_map_storage_index(gv,0,0);

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
        outfp = fopen(*lpo_filename,"a");
      else
        outfp = fopen(*lpo_filename,"w");

      if (NULL == outfp)
        {
          gs_info(gv,-2,"File open error");
          failed = 1;
          goto cleanup;
        }

      if ('l' == t)
        fprintf(outfp,"I;M;A\n");

      for (gs_UNSIGNED ii = 0; ii < npli; ii++)
        {
          gs_UNSIGNED i = pli_idx[ii];

          for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            {
              size_t idx = (size_t)m * (size_t)gv->NoInd + i;

              if (GS_ALLELE_MISSING != gv->All1[idx])
                {
                  fprintf(outfp,"%s;%s;%u\n",
                          gv->IndName[i], gv->MarName[m],
                          (unsigned)gv->All1[idx]);
                  if (gv->All1[idx] != gv->All2[idx])
                    fprintf(outfp,"%s;%s;%u\n",
                            gv->IndName[i], gv->MarName[m],
                            (unsigned)gv->All2[idx]);
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
          outfp = fopen(*mmp_filename,"w");
          if (NULL == outfp)
            {
              gs_info(gv,-2,"File open error");
              failed = 1;
              goto cleanup;
            }

          for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
            if (gv->mdl[lne] < gs_UNSIGNED_MAX)
              fprintf(outfp,"%i %15.12f %s %s\n",
                      gv->map[lne]->chrom,
                      gv->map[lne]->pos,
                      gv->map[lne]->name,
                      gv->map[lne]->class);

          fclose(outfp);
          outfp = NULL;
        }

      outfp = fopen(*mpo_filename,"w");
      if (NULL == outfp)
        {
          gs_info(gv,-2,"File open error");
          failed = 1;
          goto cleanup;
        }

      for (gs_UNSIGNED jj = 0; jj < npli; jj++)
        fprintf(outfp,"%s ",gv->IndName[pli_idx[jj]]);
      fprintf(outfp,"\n");

      if (NULL != gv->map)
        {
          for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
            if (gv->mdl[lne] < gs_UNSIGNED_MAX)
              {
                gs_UNSIGNED m = gv->mdl[lne];
                fprintf(outfp,"%s ",gv->MarName[m]);
                for (gs_UNSIGNED jj = 0; jj < npli; jj++)
                  {
                    gs_UNSIGNED i = pli_idx[jj];
                    size_t idx = (size_t)m * (size_t)gv->NoInd + i;
                    if (GS_ALLELE_MISSING == gv->All1[idx])
                      fprintf(outfp,"-- ");
                    else
                      fprintf(outfp,"%u/%u ",(unsigned)gv->All1[idx],
                              (unsigned)gv->All2[idx]);
                  }
                fprintf(outfp,"\n");
              }
        }
      else
        {
          for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            {
              fprintf(outfp,"%s ",gv->MarName[m]);
              for (gs_UNSIGNED jj = 0; jj < npli; jj++)
                {
                  gs_UNSIGNED i = pli_idx[jj];
                  size_t idx = (size_t)m * (size_t)gv->NoInd + i;
                  if (GS_ALLELE_MISSING == gv->All1[idx])
                      fprintf(outfp,"-- ");
                    else
                      fprintf(outfp,"%u/%u ",(unsigned)gv->All1[idx],
                              (unsigned)gv->All2[idx]);
                }
              fprintf(outfp,"\n");
            }
        }

      fclose(outfp);
      outfp = NULL;
    }

  if ('n' == t)
    {
      gs_name_index_type *name_index = NULL;

      n_names = gs_malloc_array((size_t)((0 == npli) ? 1 : npli),sizeof(*n_names));
      if (NULL == n_names)
        {
          gs_info(gv,-2,"Memory allocation error");
          failed = 1;
          goto cleanup;
        }

      for (gs_UNSIGNED jj = 0; jj < npli; jj++)
        {
          strcpy(n_names[jj], gv->IndName[pli_idx[jj]]);
          gs_str_rem_us(n_names[jj]);
        }

      if (0 < npli)
        {
          name_index = gs_malloc_array((size_t)npli,sizeof(*name_index));
          if (NULL == name_index)
            {
              gs_info(gv,-2,"Memory allocation error");
              failed = 1;
              goto cleanup;
            }

          for (gs_UNSIGNED jj = 0; jj < npli; jj++)
            {
              name_index[jj].name = n_names[jj];
              name_index[jj].idx = pli_idx[jj];
            }
          qsort(name_index,npli,sizeof(*name_index),gs_name_index_cmp);

          for (gs_UNSIGNED jj = 1; jj < npli; jj++)
            if ((0 == strcmp(name_index[jj-1].name,name_index[jj].name)) &&
                (name_index[jj-1].idx != name_index[jj].idx))
              {
                free(name_index);
                gs_info(gv,-2,"Individual names collide in NTSys output");
                failed = 1;
                goto cleanup;
              }
          free(name_index);
        }

      if ((NULL == gv->All) || (NULL == gv->NoAll))
        {
          gs_info(gv,-2,"Marker statistics unavailable");
          failed = 1;
          goto cleanup;
        }

      if (NULL != gv->map)
        {
          for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
            if ((gv->mdl[lne] < gs_UNSIGNED_MAX) &&
                (NULL != strchr(gv->MarName[gv->mdl[lne]], '.')))
              {
                gs_info(gv,-2,"Marker name contains '.' in NTSys output");
                failed = 1;
                goto cleanup;
              }
        }
      else
        {
          for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
            if (NULL != strchr(gv->MarName[m], '.'))
              {
                gs_info(gv,-2,"Marker name contains '.' in NTSys output");
                failed = 1;
                goto cleanup;
              }
        }

      if (NULL != gv->map)
        {
          outfp = fopen(*nmp_filename,"w");
          if (NULL == outfp)
            {
              gs_info(gv,-2,"File open error");
              failed = 1;
              goto cleanup;
            }

          for (gs_UNSIGNED lne = 0; lne < gv->NLoc; lne++)
            if (gv->mdl[lne] < gs_UNSIGNED_MAX)
              fprintf(outfp,"%i %15.12f %s %s\n",
                      gv->map[lne]->chrom,
                      gv->map[lne]->pos / 100,
                      gv->map[lne]->name,
                      gv->map[lne]->class);
          fclose(outfp);
          outfp = NULL;

          outfp = fopen(*ngp_filename,"w");
          if (NULL == outfp)
            {
              gs_info(gv,-2,"File open error");
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
                    fprintf(outfp,"%9.4f\n",
                            0.01 + gv->map[last_lne]->pos / 100);
                  last_lne = lne;
                }

            if (gs_UNSIGNED_MAX != last_lne)
              fprintf(outfp,"%9.4f\n",
                      0.01 + gv->map[last_lne]->pos / 100);
          }

          fclose(outfp);
          outfp = NULL;
        }

      outfp = fopen(*npo_filename,"w");
      if (NULL == outfp)
        {
          gs_info(gv,-2,"File open error");
          failed = 1;
          goto cleanup;
        }

      for (gs_UNSIGNED jj = 0; jj < npli; jj++)
        fprintf(outfp,"%s%s ",n_names[jj],*add_inum);
      fprintf(outfp,"\n");

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
                    gs_info(gv,-2,"Marker name contains '.' in NTSys output");
                    failed = 1;
                    goto cleanup;
                  }

                for (gs_UNSIGNED a = 0; a < gv->NoAll[m]; a++)
                  if (GS_ALLELE_MISSING != gv->All[m][a])
                    {
                      fprintf(outfp,"%s.%u ",gv->MarName[m],(unsigned)gv->All[m][a]);
                      for (gs_UNSIGNED jj = 0; jj < npli; jj++)
                        {
                          gs_UNSIGNED i = pli_idx[jj];
                          size_t idx = (size_t)m * (size_t)gv->NoInd + i;
                          int present = ((gv->All[m][a] == gv->All1[idx]) ||
                                         (gv->All[m][a] == gv->All2[idx])) ? 1 : 0;
                          fprintf(outfp,"%i ",present);
                        }
                      fprintf(outfp,"\n");
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
                    fprintf(outfp,"%s.%u ",gv->MarName[m],(unsigned)gv->All[m][a]);
                    for (gs_UNSIGNED jj = 0; jj < npli; jj++)
                      {
                        gs_UNSIGNED i = pli_idx[jj];
                        size_t idx = (size_t)m * (size_t)gv->NoInd + i;
                        int present = ((gv->All[m][a] == gv->All1[idx]) ||
                                       (gv->All[m][a] == gv->All2[idx])) ? 1 : 0;
                        fprintf(outfp,"%i ",present);
                      }
                    fprintf(outfp,"\n");
                  }
            }
        }

      fclose(outfp);
      outfp = NULL;
    }

cleanup:
  if (NULL != outfp) fclose(outfp);
  free(pli_idx);
  free(n_names);

  if (failed)
    {
      *retval = -2;
      return;
    }

  *retval = 0;
}


void gs_write_marker_data_02_GV ( 
				 char** type ,
				 char** lpo_filename,
				 char** mmp_filename,
				 char** mpo_filename, 
				 char** ngp_filename, 
				 char** nmp_filename, 
				 char** npo_filename, 
				 char** ind_filename,
				 int*   f_ind,
				 int*   l_ind,
				 char** add_inum,
				 int*   auxfiles,
				 int*   retval,         
				 char** set_name)
{

  gs_write_marker_data_02 ( gs_fdta(*set_name),
			    type ,          
			    lpo_filename,   
			    mmp_filename,   
			    mpo_filename,   
			    ngp_filename,   
			    nmp_filename,   
			    npo_filename,   
			    ind_filename,   
			    f_ind,          
			    l_ind,          
			    add_inum,       
			    auxfiles,       
			    retval         ) ;
}

void gs_write_map_01 ( gs_varset_type *gv  ,
		       char** mmp_filename ,
		       int*   retval        )
{

  if (NULL == gv->map)   
    {
      snprintf(gv->msg, sizeof(gv->msg),
         "No linkage map was loaded");
      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }

  gs_UNSIGNED lne;
  
  FILE* outfp;

  if ( (outfp = fopen(*mmp_filename,"w")) == NULL ) ERR_F;

  for (lne=0; lne <  gv->NLoc; lne++)
    {
      if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
	{
	  snprintf(gv->msg, sizeof(gv->msg),
		   "%i %15.12f %s %s\n",
		   (gv->map[lne])->chrom,
		   (gv->map[lne])->pos,
		   (gv->map[lne])->name,
		   (gv->map[lne])->class);

	  fprintf(outfp,"%s",gv->msg);
	}
    }
  fclose(outfp);
  
  *retval = 0;   return ;
}

void gs_write_map_01_GV ( char** mmp_filename,
			  int*   retval,         
			  char** set_name      )
{

  gs_write_map_01 ( gs_fdta(*set_name),
		    mmp_filename,   
		    retval             ) ;
}


int gs_determine_alleles(GS_ALLELE_TYPE *a,
                         GS_ALLELE_TYPE *b,
                         const char *p)
{
  const char *sep;
  size_t len;

  if ((NULL == a) || (NULL == b) || (NULL == p) || ('\0' == *p))
    return -1;

  len = strlen(p);

  if (gs_allele_missing_span(p,len))
    {
      *a = GS_ALLELE_MISSING;
      *b = GS_ALLELE_MISSING;
      return 0;
    }

  sep = strchr(p,'/');
  if (NULL != sep)
    {
      if (NULL != strchr(sep + 1,'/'))
        return -1;

      if ((0 != gs_parse_allele_span(p,(size_t)(sep-p),a)) ||
          (0 != gs_parse_allele_span(sep+1,
                                     len-(size_t)(sep-p)-1,b)))
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
      if ((0 != gs_parse_allele_span(p,1,a)) ||
          (0 != gs_parse_allele_span(p+1,1,b)))
        return -1;

      if ((GS_ALLELE_MISSING == *a) || (GS_ALLELE_MISSING == *b))
        {
          *a = GS_ALLELE_MISSING;
          *b = GS_ALLELE_MISSING;
        }
      return 0;
    }

  if (0 != gs_parse_single_allele(p,a))
    return -1;

  *b = *a;
  return 0;
}


void gs_remove_newline(char*p)
{
  char * pp = strchr(p,'\n');
  if ( NULL != pp ) *pp = '\0'; 
  pp = strchr(p,'\r');
  if ( NULL != pp ) *pp = '\0'; 
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

  if ('\r' == (*lineptr)[len-1])
    (*lineptr)[len-1] = '\n';

  (*lineptr)[len] = '\0';
  return (long long int)len;
}


void gs_read_marker_data_02 ( gs_varset_type *gv ,
                              char** dta_filename,
                              char** ind_filename,
                              char** mar_filename,
                              char** gen_filename,
                              int*   auxfiles,
                              int*   retval        )
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

  fp = fopen(*dta_filename,"r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      *retval = -2;
      return;
    }

  rc = st_getline(&line,&n,fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,(rc < 0) ? "Memory allocation error" :
                             "Marker data file contains no header");
      failed = 1;
      goto cleanup;
    }
  gs_remove_newline(line);
  {
    char *p = strtok(line,gs_tok_del);
    while (NULL != p)
      {
        header_fields++;
        p = strtok(NULL,gs_tok_del);
      }
  }

  rc = st_getline(&line,&n,fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,"Marker data file contains no data rows");
      failed = 1;
      goto cleanup;
    }
  gs_remove_newline(line);
  {
    char *p = strtok(line,gs_tok_del);
    while (NULL != p)
      {
        row_fields++;
        p = strtok(NULL,gs_tok_del);
      }
  }

  if ((0 == header_fields) || (row_fields < 2))
    {
      gs_info(gv,-2,"Invalid marker matrix dimensions");
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
          gs_info(gv,-2,"Invalid marker matrix header");
          failed = 1;
          goto cleanup;
        }
      corner = 1;
      no_ind = header_fields - 1;
    }
  else
    {
      gs_info(gv,-2,"Marker matrix header and data row have different widths");
      failed = 1;
      goto cleanup;
    }

  fclose(fp);
  fp = NULL;

  fp = fopen(*dta_filename,"r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      failed = 1;
      goto cleanup;
    }
  rc = st_getline(&line,&n,fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,"Error reading marker matrix header");
      failed = 1;
      goto cleanup;
    }

  while (0 < (rc = st_getline(&line,&n,fp)))
    {
      gs_UNSIGNED fields = 0;
      char *p;

      gs_remove_newline(line);
      p = strtok(line,gs_tok_del);
      while (NULL != p)
        {
          fields++;
          p = strtok(NULL,gs_tok_del);
        }

      if (fields != no_ind + 1)
        {
          gs_info(gv,-2,"Marker matrix row has wrong number of fields");
          failed = 1;
          goto cleanup;
        }
      no_mar++;
    }

  if (rc < 0)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  fclose(fp);
  fp = NULL;

  if (0 == no_mar)
    {
      gs_info(gv,-2,"Marker data file contains no marker rows");
      failed = 1;
      goto cleanup;
    }

  tmp.NoInd = no_ind;
  tmp.NoMar = no_mar;
  tmp.IndName = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.IndName));
  tmp.MarName = gs_malloc_array((size_t)tmp.NoMar,sizeof(*tmp.MarName));
  if ((NULL == tmp.IndName) || (NULL == tmp.MarName) ||
      (0 != gs_alloc_genotype_matrix(&tmp,tmp.NoMar,tmp.NoInd)))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  fp = fopen(*dta_filename,"r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      failed = 1;
      goto cleanup;
    }

  rc = st_getline(&line,&n,fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,"Error reading marker matrix header");
      failed = 1;
      goto cleanup;
    }
  gs_remove_newline(line);
  {
    char *p = strtok(line,gs_tok_del);
    gs_UNSIGNED i = 0;

    if (corner) p = strtok(NULL,gs_tok_del);
    while (NULL != p)
      {
        size_t len = strlen(p);
        if ((i >= tmp.NoInd) || (len >= gs_MAXNME))
          {
            gs_info(gv,-2,"Individual name too long or invalid header");
            failed = 1;
            goto cleanup;
          }
        strcpy(tmp.IndName[i],p);
        if ((int)len > tmp.MaxIL) tmp.MaxIL = (int)len;
        i++;
        p = strtok(NULL,gs_tok_del);
      }

    if (i != tmp.NoInd)
      {
        gs_info(gv,-2,"Invalid marker matrix header");
        failed = 1;
        goto cleanup;
      }
  }

  if (0 != gs_make_name_index(tmp.IndName,tmp.NoInd,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,tmp.NoInd))
    {
      gs_info(gv,-2,"Duplicate individual name in marker matrix");
      failed = 1;
      goto cleanup;
    }
  free(index);
  index = NULL;

  for (gs_UNSIGNED m = 0; m < tmp.NoMar; m++)
    {
      char *p;
      size_t len;

      rc = st_getline(&line,&n,fp);
      if (rc <= 0)
        {
          gs_info(gv,-2,"Unexpected end of marker matrix");
          failed = 1;
          goto cleanup;
        }
      gs_remove_newline(line);
      p = strtok(line,gs_tok_del);
      if (NULL == p)
        {
          gs_info(gv,-2,"Missing marker name");
          failed = 1;
          goto cleanup;
        }
      len = strlen(p);
      if (len >= gs_MAXNME)
        {
          gs_info(gv,-2,"Marker name too long");
          failed = 1;
          goto cleanup;
        }
      strcpy(tmp.MarName[m],p);
      if ((int)len > tmp.MaxML) tmp.MaxML = (int)len;

      for (gs_UNSIGNED i = 0; i < tmp.NoInd; i++)
        {
          GS_ALLELE_TYPE a,b;
          p = strtok(NULL,gs_tok_del);
          if (NULL == p)
            {
              gs_info(gv,-2,"Marker matrix row is too short");
              failed = 1;
              goto cleanup;
            }
          if (0 != gs_determine_alleles(&a,&b,p))
            {
              gs_info(gv,-2,"Invalid allele in marker matrix");
              failed = 1;
              goto cleanup;
            }
          tmp.All1[(size_t)m*tmp.NoInd+i] = a;
          tmp.All2[(size_t)m*tmp.NoInd+i] = b;
        }

      if (NULL != strtok(NULL,gs_tok_del))
        {
          gs_info(gv,-2,"Marker matrix row is too long");
          failed = 1;
          goto cleanup;
        }
    }

  fclose(fp);
  fp = NULL;

  if (0 != gs_make_name_index(tmp.MarName,tmp.NoMar,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,tmp.NoMar))
    {
      gs_info(gv,-2,"Duplicate marker name in marker matrix");
      failed = 1;
      goto cleanup;
    }

  if (0 != gs_finish_marker_read(gv,&tmp,
                                 ind_filename,mar_filename,gen_filename,
                                 auxfiles,retval))
    failed = 1;

cleanup:
  if (NULL != fp) fclose(fp);
  free(line);
  free(index);
  if (failed)
    {
      gs_reset(&tmp);
      *retval = -2;
    }
}


void gs_read_marker_data_02_GV ( char** dta_filename,
				 char** ind_filename,
				 char** mar_filename,
				 char** gen_filename,
				 int*   auxfiles,
				 int*   retval,
				 char** set_name)
{
  gs_read_marker_data_02 ( gs_fdta(*set_name),
			   dta_filename,
			   ind_filename,
			   mar_filename,
			   gen_filename,
			   auxfiles,
			   retval         ) ;

}

void gs_read_marker_data_04 ( gs_varset_type *gv ,
                              char** dta_filename,
                              char** ind_filename,
                              char** mar_filename,
                              char** gen_filename,
                              int*   auxfiles,
                              int*   retval        )
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

  gs_tmp_init_from(&tmp,gv);

  fp = fopen(*dta_filename,"r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      *retval = -2;
      return;
    }

  rc = st_getline(&line,&n,fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,"Marker data file contains no header");
      failed = 1;
      goto cleanup;
    }
  gs_remove_newline(line);
  header_fields = gs_count_sep_fields(line,sep);

  rc = st_getline(&line,&n,fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,"Marker data file contains no data rows");
      failed = 1;
      goto cleanup;
    }
  gs_remove_newline(line);
  data_fields = gs_count_sep_fields(line,sep);

  if (data_fields == header_fields + 1)
    {
      corner = 0;
      no_mar = (gs_UNSIGNED)header_fields;
    }
  else if (data_fields == header_fields)
    {
      if (header_fields < 2)
        {
          gs_info(gv,-2,"Invalid transposed marker matrix header");
          failed = 1;
          goto cleanup;
        }
      corner = 1;
      no_mar = (gs_UNSIGNED)(header_fields - 1);
    }
  else
    {
      gs_info(gv,-2,"Transposed marker matrix has inconsistent width");
      failed = 1;
      goto cleanup;
    }

  if (0 == no_mar)
    {
      gs_info(gv,-2,"Marker matrix contains no markers");
      failed = 1;
      goto cleanup;
    }

  no_ind = 1;
  while (0 < (rc = st_getline(&line,&n,fp)))
    {
      gs_remove_newline(line);
      if (gs_count_sep_fields(line,sep) != (size_t)no_mar + 1)
        {
          gs_info(gv,-2,"Transposed marker matrix row has wrong width");
          failed = 1;
          goto cleanup;
        }
      no_ind++;
    }
  if (rc < 0)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  fclose(fp);
  fp = NULL;

  tmp.NoMar = no_mar;
  tmp.NoInd = no_ind;
  tmp.MarName = gs_malloc_array((size_t)tmp.NoMar,sizeof(*tmp.MarName));
  tmp.IndName = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.IndName));
  fields = gs_malloc_array((size_t)tmp.NoMar + 1,sizeof(*fields));

  if ((NULL == tmp.MarName) || (NULL == tmp.IndName) ||
      (NULL == fields) ||
      (0 != gs_alloc_genotype_matrix(&tmp,tmp.NoMar,tmp.NoInd)))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  fp = fopen(*dta_filename,"r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      failed = 1;
      goto cleanup;
    }

  rc = st_getline(&line,&n,fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,"Error reading transposed marker matrix header");
      failed = 1;
      goto cleanup;
    }
  gs_remove_newline(line);
  {
    size_t cap = (size_t)tmp.NoMar + (corner ? 1 : 0);
    char **header = gs_malloc_array(cap,sizeof(*header));
    size_t nf;
    size_t start = corner ? 1 : 0;

    if (NULL == header)
      {
        gs_info(gv,-2,"Memory allocation error");
        failed = 1;
        goto cleanup;
      }

    nf = gs_split_sep_fields(line,sep,header,cap);
    if (nf != cap)
      {
        free(header);
        gs_info(gv,-2,"Invalid transposed marker matrix header");
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
            gs_info(gv,-2,"Invalid marker name in transposed matrix");
            failed = 1;
            goto cleanup;
          }
        strcpy(tmp.MarName[m],name);
        if ((int)len > tmp.MaxML) tmp.MaxML = (int)len;
      }
    free(header);
  }

  if (0 != gs_make_name_index(tmp.MarName,tmp.NoMar,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,tmp.NoMar))
    {
      gs_info(gv,-2,"Duplicate marker name in transposed matrix");
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

      rc = st_getline(&line,&n,fp);
      if (rc <= 0)
        {
          gs_info(gv,-2,"Unexpected end of transposed marker matrix");
          failed = 1;
          goto cleanup;
        }
      gs_remove_newline(line);
      nf = gs_split_sep_fields(line,sep,fields,(size_t)tmp.NoMar+1);
      if (nf != (size_t)tmp.NoMar + 1)
        {
          gs_info(gv,-2,"Transposed marker matrix row has wrong width");
          failed = 1;
          goto cleanup;
        }

      name = gs_trim_field(fields[0]);
      len = strlen(name);
      if ((0 == len) || (len >= gs_MAXNME))
        {
          gs_info(gv,-2,"Invalid individual name in transposed matrix");
          failed = 1;
          goto cleanup;
        }
      strcpy(tmp.IndName[i],name);
      if ((int)len > tmp.MaxIL) tmp.MaxIL = (int)len;

      for (gs_UNSIGNED m = 0; m < tmp.NoMar; m++)
        {
          char *geno = gs_trim_field(fields[m+1]);
          GS_ALLELE_TYPE a,b;
          if (0 != gs_determine_alleles(&a,&b,geno))
            {
              gs_info(gv,-2,"Invalid allele in transposed marker matrix");
              failed = 1;
              goto cleanup;
            }
          tmp.All1[(size_t)m*tmp.NoInd+i] = a;
          tmp.All2[(size_t)m*tmp.NoInd+i] = b;
        }
    }

  fclose(fp);
  fp = NULL;

  if (0 != gs_make_name_index(tmp.IndName,tmp.NoInd,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,tmp.NoInd))
    {
      gs_info(gv,-2,"Duplicate individual name in transposed matrix");
      failed = 1;
      goto cleanup;
    }

  if (0 != gs_finish_marker_read(gv,&tmp,
                                 ind_filename,mar_filename,gen_filename,
                                 auxfiles,retval))
    failed = 1;

cleanup:
  if (NULL != fp) fclose(fp);
  free(line);
  free(fields);
  free(index);
  if (failed)
    {
      gs_reset(&tmp);
      *retval = -2;
    }
}


void gs_read_marker_data_04_GV ( char** dta_filename,
				 char** ind_filename,
				 char** mar_filename,
				 char** gen_filename,
				 int*   auxfiles,
				 int*   retval,
				 char** set_name)
{
  gs_read_marker_data_04 ( gs_fdta(*set_name),
			   dta_filename,
			   ind_filename,
			   mar_filename,
			   gen_filename,
			   auxfiles,
			   retval         ) ;

}


void gs_marker_stats_01 ( gs_varset_type *gv ,
                          char** ind_filename,
                          char** mar_filename,
                          char** gen_filename,
                          int*   ind,
                          int*   mar,
                          int*   gen,
                          int*   auxfiles,
                          int*   retval        )
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

  if ((0 == gv->NoInd) || (0 == gv->NoMar) ||
      (NULL == gv->All1) || (NULL == gv->All2) ||
      (NULL == gv->IndName) || (NULL == gv->MarName))
    {
      gs_info(gv,-2,"No marker data available");
      *retval = -2;
      return;
    }

  if (!gs_marker_genotypes_valid(gv))
    {
      gs_info(gv,-2,"Invalid marker genotype representation");
      *retval = -2;
      return;
    }

  snprintf(gv->msg, sizeof(gv->msg),
           "No. of individuals: %lu, no. of markers: %lu",
           gv->NoInd, gv->NoMar);
  gs_info(gv,0,gv->msg);

  NoAll = gs_malloc_array((size_t)gv->NoMar,sizeof(*NoAll));
  Pic = gs_malloc_array((size_t)gv->NoMar,sizeof(*Pic));
  Mis = gs_malloc_array((size_t)gv->NoMar,sizeof(*Mis));
  IndMis = gs_malloc_array((size_t)gv->NoInd,sizeof(*IndMis));
  all_work = gs_malloc_array((size_t)gs_MAXALL,sizeof(*all_work));

  if ((NULL == NoAll) || (NULL == Pic) || (NULL == Mis) ||
      (NULL == IndMis) || (NULL == all_work))
    {
      free(NoAll); free(Pic); free(Mis); free(IndMis); free(all_work);
      gs_info(gv,-2,"Memory allocation error");
      *retval = -2;
      return;
    }

  /* First pass: determine the number of alleles per marker. */
  for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
      NoAll[m] = 0;

      for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        {
          size_t idx = (size_t)m * (size_t)gv->NoInd + i;
          GS_ALLELE_TYPE alleles[2] = { gv->All1[idx], gv->All2[idx] };

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
                      free(NoAll); free(Pic); free(Mis);
                      free(IndMis); free(all_work);
                      gs_info(gv,-2,"Too many alleles per marker");
                      *retval = -2;
                      return;
                    }
                  all_work[NoAll[m]++] = alleles[h];
                }
            }
        }

      if ((size_t)NoAll[m] > ((size_t)-1) - total_all)
        {
          free(NoAll); free(Pic); free(Mis);
          free(IndMis); free(all_work);
          gs_info(gv,-2,"Memory allocation error");
          *retval = -2;
          return;
        }
      total_all += (size_t)NoAll[m];
    }

  All = gs_malloc_array((size_t)gv->NoMar,sizeof(*All));
  All_S = gs_malloc_array(total_all,sizeof(*All_S));
  AllCnt = gs_malloc_array((size_t)gv->NoMar,sizeof(*AllCnt));
  AllCnt_S = gs_malloc_array(total_all,sizeof(*AllCnt_S));

  if ((NULL == All) || (NULL == All_S) ||
      (NULL == AllCnt) || (NULL == AllCnt_S))
    {
      free(NoAll); free(All_S); free(All);
      free(AllCnt_S); free(AllCnt);
      free(Pic); free(Mis); free(IndMis); free(all_work);
      gs_info(gv,-2,"Memory allocation error");
      *retval = -2;
      return;
    }

  for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
      All[m] = All_S + offset;
      AllCnt[m] = AllCnt_S + offset;
      offset += (size_t)NoAll[m];
    }

  /* Second pass: determine alleles and allele counts. */
  for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
    {
      gs_UNSIGNED no_all = 0;

      for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        {
          size_t idx = (size_t)m * (size_t)gv->NoInd + i;
          GS_ALLELE_TYPE alleles[2] = { gv->All1[idx], gv->All2[idx] };

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
                      free(NoAll); free(All_S); free(All);
                      free(AllCnt_S); free(AllCnt);
                      free(Pic); free(Mis); free(IndMis); free(all_work);
                      gs_info(gv,-2,"Too many alleles per marker");
                      *retval = -2;
                      return;
                    }
                  all_work[no_all++] = alleles[h];
                }
            }
        }

      if (no_all != NoAll[m])
        {
          free(NoAll); free(All_S); free(All);
          free(AllCnt_S); free(AllCnt);
          free(Pic); free(Mis); free(IndMis); free(all_work);
          gs_info(gv,-2,"Marker statistics error");
          *retval = -2;
          return;
        }

      memcpy(All[m],all_work,(size_t)NoAll[m] * sizeof(*All[m]));
      qsort(All[m],NoAll[m],sizeof(*All[m]),gs_allele_cmp);

      for (gs_UNSIGNED k = 0; k < NoAll[m]; k++)
        {
          AllCnt[m][k] = 0;
          for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            {
              size_t idx = (size_t)m * (size_t)gv->NoInd + i;
              if (gv->All1[idx] == All[m][k]) AllCnt[m][k]++;
              if (gv->All2[idx] == All[m][k]) AllCnt[m][k]++;
            }
        }

      {
        double observed = 2.0 * (double)gv->NoInd;
        gs_UNSIGNED real_no_all = NoAll[m];

        if ((0 < NoAll[m]) &&
            (GS_ALLELE_MISSING == All[m][NoAll[m]-1]))
          {
            gs_UNSIGNED km = NoAll[m]-1;
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

  for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
    {
      gs_UNSIGNED count_miss = 0;
      for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
        if (GS_ALLELE_MISSING == gv->All1[(size_t)m*gv->NoInd+i]) count_miss++;
      IndMis[i] = (double)count_miss / (double)gv->NoMar;
    }

  free(gv->NoAll);     gv->NoAll = NoAll;
  free(gv->All_S);     gv->All_S = All_S;
  free(gv->All);       gv->All = All;
  free(gv->AllCnt_S);  gv->AllCnt_S = AllCnt_S;
  free(gv->AllCnt);    gv->AllCnt = AllCnt;
  free(gv->Pic);       gv->Pic = Pic;
  free(gv->Mis);       gv->Mis = Mis;
  free(gv->IndMis);    gv->IndMis = IndMis;

  if ((1 == *auxfiles) && (1 == *ind))
    {
      outfp = fopen(*ind_filename,"w");
      if (NULL == outfp)
        {
          gs_info(gv,-2,"File open error");
          *retval = -2;
          return;
        }

      fprintf(outfp,"Name InMis\n");
      for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        fprintf(outfp,"%s %f\n",gv->IndName[i],gv->IndMis[i]);
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
                    gs_UNSIGNED new_capacity = (0 == capacity) ? 16 : capacity * 2;
                    int *new_obs;

                    if (new_capacity <= capacity)
                      {
                        free(obs_all);
                        gs_info(gv,-2,"Too many observed allele codes");
                        *retval = -2;
                        return;
                      }

                    new_obs = gs_realloc_array(obs_all,
                                               (size_t)new_capacity,
                                               sizeof(*new_obs));
                    if (NULL == new_obs)
                      {
                        free(obs_all);
                        gs_info(gv,-2,"Memory allocation error");
                        *retval = -2;
                        return;
                      }
                    obs_all = new_obs;
                    capacity = new_capacity;
                  }

                obs_all[no_obs_all++] = allele;
              }
          }

      qsort(obs_all,no_obs_all,sizeof(*obs_all),gs_intcmp);

      outfp = fopen(*mar_filename,"w");
      if (NULL == outfp)
        {
          free(obs_all);
          gs_info(gv,-2,"File open error");
          *retval = -2;
          return;
        }

      fprintf(outfp,"Name NoAll MaMis ExHet ");
      if (have_missing) fprintf(outfp,"AM ");
      for (gs_UNSIGNED l = 0; l < no_obs_all; l++)
        fprintf(outfp," A%i",obs_all[l]);
      fprintf(outfp,"\n");

      for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
        {
          gs_UNSIGNED real_no_all = gv->NoAll[m];
          gs_UNSIGNED miss_count = 0;

          if ((0 < gv->NoAll[m]) &&
              (GS_ALLELE_MISSING == gv->All[m][gv->NoAll[m]-1]))
            {
              real_no_all--;
              miss_count = gv->AllCnt[m][gv->NoAll[m]-1];
            }

          fprintf(outfp,"%s %lu %5.3f %5.3f ",
                  gv->MarName[m],real_no_all,gv->Mis[m],gv->Pic[m]);

          if (have_missing) fprintf(outfp," %lu",miss_count);

          for (gs_UNSIGNED l = 0; l < no_obs_all; l++)
            {
              gs_UNSIGNED ac = 0;
              for (gs_UNSIGNED k = 0; k < real_no_all; k++)
                if (obs_all[l] == (int)gv->All[m][k])
                  ac = gv->AllCnt[m][k];
              fprintf(outfp," %lu",ac);
            }
          fprintf(outfp,"\n");
        }

      fclose(outfp);
      free(obs_all);
    }

  if ((1 == *auxfiles) && (1 == *gen))
    {
      outfp = fopen(*gen_filename,"w");
      if (NULL == outfp)
        {
          gs_info(gv,-2,"File open error");
          *retval = -2;
          return;
        }

      fprintf(outfp,"Mar/Ind ");
      for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
        fprintf(outfp,"%s ",gv->IndName[i]);
      fprintf(outfp,"\n");

      for (gs_UNSIGNED m = 0; m < gv->NoMar; m++)
        {
          fprintf(outfp,"%s ",gv->MarName[m]);
          for (gs_UNSIGNED i = 0; i < gv->NoInd; i++)
            {
              size_t idx = (size_t)m * (size_t)gv->NoInd + i;
              if (GS_ALLELE_MISSING == gv->All1[idx])
                fprintf(outfp,"-- ");
              else
                fprintf(outfp,"%u/%u ",(unsigned)gv->All1[idx],
                        (unsigned)gv->All2[idx]);
            }
          fprintf(outfp,"\n");
        }
      fclose(outfp);
    }

  *retval = 0;
}


void gs_marker_stats_01_GV ( 
			     char** ind_filename,
			     char** mar_filename,
			     char** gen_filename,
			     int*   ind,
			     int*   mar,
			     int*   gen,
			     int*   auxfiles,
			     int*   retval,
			     char** set_name)
{
  gs_marker_stats_01 (     gs_fdta(*set_name),
			   ind_filename,
			   mar_filename,
			   gen_filename,
			   ind,
			   mar,
			   gen,
			   auxfiles,
			   retval         ) ;

}

void gs_read_marker_data_03 ( gs_varset_type *gv ,
                              char** dta_filename,
                              char** ind_filename,
                              char** mar_filename,
                              char** gen_filename,
                              int*   auxfiles,
                              int*   retval        )
{
  FILE *fp = NULL;
  char *line = NULL;
  size_t n = 0;
  long long int rc;
  gs_varset_type tmp;
  gs_UNSIGNED marker_capacity = 0;
  gs_name_index_type *index = NULL;
  int failed = 0;

  gs_tmp_init_from(&tmp,gv);

  fp = fopen(*dta_filename,"r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      *retval = -2;
      return;
    }

  rc = st_getline(&line,&n,fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,"NTSys marker file contains no header");
      failed = 1;
      goto cleanup;
    }
  gs_remove_newline(line);

  {
    char *p = strtok(line,gs_tok_del);
    while (NULL != p)
      {
        size_t len = strlen(p);
        char (*new_names)[gs_MAXNME];

        if (len >= gs_MAXNME)
          {
            gs_info(gv,-2,"Individual name too long");
            failed = 1;
            goto cleanup;
          }

        if (tmp.NoInd == marker_capacity)
          {
            gs_UNSIGNED new_capacity = (0 == marker_capacity) ? 64 : marker_capacity * 2;
            if (new_capacity <= marker_capacity)
              {
                gs_info(gv,-2,"Too many individuals");
                failed = 1;
                goto cleanup;
              }
            new_names = gs_realloc_array(tmp.IndName,
                                         (size_t)new_capacity,
                                         sizeof(*new_names));
            if (NULL == new_names)
              {
                gs_info(gv,-2,"Memory allocation error");
                failed = 1;
                goto cleanup;
              }
            tmp.IndName = new_names;
            marker_capacity = new_capacity;
          }

        strcpy(tmp.IndName[tmp.NoInd++],p);
        if ((int)len > tmp.MaxIL) tmp.MaxIL = (int)len;
        p = strtok(NULL,gs_tok_del);
      }
  }

  if (0 == tmp.NoInd)
    {
      gs_info(gv,-2,"NTSys marker file contains no individuals");
      failed = 1;
      goto cleanup;
    }

  if (0 != gs_make_name_index(tmp.IndName,tmp.NoInd,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,tmp.NoInd))
    {
      gs_info(gv,-2,"Duplicate individual name in NTSys marker file");
      failed = 1;
      goto cleanup;
    }
  free(index);
  index = NULL;

  marker_capacity = 0;
  while (0 < (rc = st_getline(&line,&n,fp)))
    {
      char *p;
      char *dot;
      char marker[gs_MAXNME];
      size_t len;
      gs_UNSIGNED fields = 0;

      gs_remove_newline(line);
      p = strtok(line,gs_tok_del);
      if (NULL == p)
        {
          gs_info(gv,-2,"Empty line in NTSys marker file");
          failed = 1;
          goto cleanup;
        }

      dot = strchr(p,'.');
      if ((NULL == dot) || (dot == p) || ('\0' == dot[1]))
        {
          gs_info(gv,-2,"Invalid marker.allele field in NTSys file");
          failed = 1;
          goto cleanup;
        }

      len = (size_t)(dot - p);
      if (len >= gs_MAXNME)
        {
          gs_info(gv,-2,"Marker name too long");
          failed = 1;
          goto cleanup;
        }
      memcpy(marker,p,len);
      marker[len] = '\0';

      {
        long allele;
        if ((0 != gs_parse_long_strict(dot+1,&allele)) ||
            (allele < 1) || (allele > GS_ALLELE_MAX))
          {
            gs_info(gv,-2,"Invalid allele number in NTSys file");
            failed = 1;
            goto cleanup;
          }
      }

      while (NULL != (p = strtok(NULL,gs_tok_del)))
        {
          long value;
          if ((0 != gs_parse_long_strict(p,&value)) ||
              !((0 == value) || (1 == value)))
            {
              gs_info(gv,-2,"NTSys incidence value must be 0 or 1");
              failed = 1;
              goto cleanup;
            }
          fields++;
        }

      if (fields != tmp.NoInd)
        {
          gs_info(gv,-2,"NTSys marker row has wrong number of fields");
          failed = 1;
          goto cleanup;
        }

      if ((0 == tmp.NoMar) ||
          (0 != strcmp(tmp.MarName[tmp.NoMar-1],marker)))
        {
          char (*new_names)[gs_MAXNME];

          if (tmp.NoMar == marker_capacity)
            {
              gs_UNSIGNED new_capacity = (0 == marker_capacity) ? 64 : marker_capacity * 2;
              if (new_capacity <= marker_capacity)
                {
                  gs_info(gv,-2,"Too many markers");
                  failed = 1;
                  goto cleanup;
                }
              new_names = gs_realloc_array(tmp.MarName,
                                           (size_t)new_capacity,
                                           sizeof(*new_names));
              if (NULL == new_names)
                {
                  gs_info(gv,-2,"Memory allocation error");
                  failed = 1;
                  goto cleanup;
                }
              tmp.MarName = new_names;
              marker_capacity = new_capacity;
            }

          strcpy(tmp.MarName[tmp.NoMar++],marker);
          if ((int)len > tmp.MaxML) tmp.MaxML = (int)len;
        }
    }

  if (rc < 0)
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  fclose(fp);
  fp = NULL;

  if (0 == tmp.NoMar)
    {
      gs_info(gv,-2,"NTSys marker file contains no marker rows");
      failed = 1;
      goto cleanup;
    }

  if (0 != gs_make_name_index(tmp.MarName,tmp.NoMar,&index))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  if (gs_name_index_has_duplicates(index,tmp.NoMar))
    {
      gs_info(gv,-2,"Marker block reappears in NTSys marker file");
      failed = 1;
      goto cleanup;
    }
  free(index);
  index = NULL;

  if (0 != gs_alloc_genotype_matrix(&tmp,tmp.NoMar,tmp.NoInd))
    {
      gs_info(gv,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  gs_set_genotypes_missing(&tmp);

  fp = fopen(*dta_filename,"r");
  if (NULL == fp)
    {
      gs_info(gv,-2,"File open error");
      failed = 1;
      goto cleanup;
    }
  rc = st_getline(&line,&n,fp);
  if (rc <= 0)
    {
      gs_info(gv,-2,"Error reading NTSys header");
      failed = 1;
      goto cleanup;
    }

  {
    gs_UNSIGNED current_m = gs_UNSIGNED_MAX;
    char current_name[gs_MAXNME] = "";

    while (0 < (rc = st_getline(&line,&n,fp)))
      {
        char *p;
        char *dot;
        char marker[gs_MAXNME];
        size_t len;
        long allele_l;
        GS_ALLELE_TYPE allele;

        gs_remove_newline(line);
        p = strtok(line,gs_tok_del);
        if (NULL == p)
          {
            gs_info(gv,-2,"Error in NTSys marker file");
            failed = 1;
            goto cleanup;
          }
        dot = strchr(p,'.');
        if ((NULL == dot) || (dot == p) || ('\0' == dot[1]))
          {
            gs_info(gv,-2,"Invalid marker.allele field in NTSys file");
            failed = 1;
            goto cleanup;
          }
        len = (size_t)(dot-p);
        if (len >= gs_MAXNME)
          {
            gs_info(gv,-2,"Marker name too long");
            failed = 1;
            goto cleanup;
          }
        memcpy(marker,p,len);
        marker[len] = '\0';
        if ((0 != gs_parse_long_strict(dot+1,&allele_l)) ||
            (allele_l < 1) || (allele_l > GS_ALLELE_MAX))
          {
            gs_info(gv,-2,"Invalid allele number in NTSys file");
            failed = 1;
            goto cleanup;
          }
        allele = (GS_ALLELE_TYPE)allele_l;

        if ((gs_UNSIGNED_MAX == current_m) ||
            (0 != strcmp(current_name,marker)))
          {
            current_m++;
            if ((current_m >= tmp.NoMar) ||
                (0 != strcmp(tmp.MarName[current_m],marker)))
              {
                gs_info(gv,-2,"Invalid marker block order in NTSys file");
                failed = 1;
                goto cleanup;
              }
            strcpy(current_name,marker);
          }

        for (gs_UNSIGNED i = 0; i < tmp.NoInd; i++)
          {
            long present;
            size_t idx = (size_t)current_m * tmp.NoInd + i;

            p = strtok(NULL,gs_tok_del);
            if ((NULL == p) ||
                (0 != gs_parse_long_strict(p,&present)) ||
                !((0 == present) || (1 == present)))
              {
                gs_info(gv,-2,"Invalid NTSys incidence row");
                failed = 1;
                goto cleanup;
              }

            if (1 == present)
              {
                if ((GS_ALLELE_MISSING == tmp.All1[idx]) &&
                    (GS_ALLELE_MISSING == tmp.All2[idx]))
                  {
                    tmp.All1[idx] = allele;
                    tmp.All2[idx] = allele;
                  }
                else if ((tmp.All1[idx] == allele) ||
                         (tmp.All2[idx] == allele))
                  {
                    continue;
                  }
                else if (tmp.All1[idx] == tmp.All2[idx])
                  {
                    tmp.All2[idx] = allele;
                  }
                else
                  {
                    gs_info(gv,-2,
                            "More than two present alleles in NTSys genotype");
                    failed = 1;
                    goto cleanup;
                  }
              }
          }

        if (NULL != strtok(NULL,gs_tok_del))
          {
            gs_info(gv,-2,"NTSys marker row has wrong number of fields");
            failed = 1;
            goto cleanup;
          }
      }

    if (rc < 0)
      {
        gs_info(gv,-2,"Memory allocation error");
        failed = 1;
        goto cleanup;
      }
  }

  fclose(fp);
  fp = NULL;

  if (0 != gs_finish_marker_read(gv,&tmp,
                                 ind_filename,mar_filename,gen_filename,
                                 auxfiles,retval))
    failed = 1;

cleanup:
  if (NULL != fp) fclose(fp);
  free(line);
  free(index);
  if (failed)
    {
      gs_reset(&tmp);
      *retval = -2;
    }
}


void gs_read_marker_data_03_GV ( char** dta_filename,
				 char** ind_filename,
				 char** mar_filename,
				 char** gen_filename,
				 int*   auxfiles,
				 int*   retval,
				 char** set_name)
{
  gs_read_marker_data_03 ( gs_fdta(*set_name),
			   dta_filename,
			   ind_filename,
			   mar_filename,
			   gen_filename,
			   auxfiles,
			   retval         ) ;

}

double st_time_store[5];

void st_start_timer (int* i, double* time)
{
  st_time_store[*i] = *time;
}

void st_stop_timer (int* i, double* time)
{
  *time = st_time_store[*i];
}

int st_info_level = 0;

extern int INFOLEVEL;

void st_set_info_level (int* i)
{
  st_info_level = *i;          // Selection Tools
  INFOLEVEL     = *i;          // Plabsim
  gs_set_all_info_levels (i); // All data sets
}

void st_get_info_level (int* i)
{
  *i = st_info_level;
}

void st_info ( int level,
               char* message )
{
    if ( st_info_level >= level )
    {
        if ( !st_console_output_allowed() )
            return;

        if ( 1 <= level )
            Rprintf( "I: %s\n", message );

        else if ( 0 == level )
            Rprintf( "M: %s\n", message );

        else if ( -1 == level )
            Rprintf( "W: %s\n", message );

        else if ( -2 == level )
            Rprintf( "E: %s\n", message );
    }
}


void st_r_info(int* level, char** message)
{
  st_info( *level, *message );
}

void gs_cross_init_01 ( gs_varset_type *gv  ,
			int*   retval         )
{
  gs_UNSIGNED i,j,k;

  if ( (NULL == gv->crs) || (NULL == gv->crs_S) )
    {
      gv->NCrs  = gv->NoInd * (gv->NoInd-1) / 2;
      gv->crs_S = realloc ( gv->crs_S, gv->NCrs * sizeof( gs_cross_type ) );
      gv->crs   = realloc ( gv->crs,   gv->NCrs * sizeof( gs_cross_type* ) );
      if ( (NULL == gv->crs) || (NULL == gv->crs_S) ) ERR_M;

      k = 0;
      
      for ( i=0; i < (gv->NoInd-1) ; i++ )  
	for ( j = i+1; j < gv->NoInd ; j++ )
	  {
	    (gv->crs_S[k]).p1 = i;
	    (gv->crs_S[k]).p2 = j;
	    gv->crs[k] = &(gv->crs_S[k]);
	    k++;
	  }

      for ( k = 0; k <  gv->NCrs ; k++)
	{
	  gv->crs[k]->gd = 0;
	  gv->crs[k]->mu = 0;
	  gv->crs[k]->mi = 0;
	  gv->crs[k]->ma = 0;
	  gv->crs[k]->va = 0;
	  gv->crs[k]->es = 0;
	}

    }

  *retval = 0; return;
      
}


void gs_calc_rd ( gs_varset_type *gv ,
		  gs_UNSIGNED     k   )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_UNSIGNED  n = 0; // non missing markers
  gs_FLOAT    di = 0; // sum of distances

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];
      
      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
      
      if ( !(  ( GS_ALLELE_MISSING == a ) ||
               ( GS_ALLELE_MISSING == b ) ||
               ( GS_ALLELE_MISSING == c ) ||
               ( GS_ALLELE_MISSING == d )    ) ) 
	{
	  n++;
	  
	  if ( !( ( (c==a)&&(d==b)) ||
		  ( (c==b)&&(d==a))    )  ) 
	    {
	      if ( ( (c==a)||(c==b)) ||
		   ( (d==a)||(d==b))    )  
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

  (gv->crs_S[k]).gd = (gs_FLOAT) di / (gs_FLOAT)n;
  
}

void gs_calc_mrd ( gs_varset_type *gv ,
	           gs_UNSIGNED     k   )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_UNSIGNED  n = 0; // non missing markers
  gs_FLOAT    di = 0; // sum of distances

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];
      
      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
      
      if ( !(  ( GS_ALLELE_MISSING == a ) ||
               ( GS_ALLELE_MISSING == b ) ||
               ( GS_ALLELE_MISSING == c ) ||
               ( GS_ALLELE_MISSING == d )    ) ) 
	{
	  n++;
	  
	  if ( !( ( (c==a)&&(d==b)) ||
		  ( (c==b)&&(d==a))    )  ) 
	    {
	      if ( ( (c==a)||(c==b)) ||
		   ( (d==a)||(d==b))    )  
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

  (gv->crs_S[k]).gd = (gs_FLOAT) sqrt( di / (gs_FLOAT)n / 2 );
  
}

void gs_calc_euc ( gs_varset_type *gv ,
	           gs_UNSIGNED     k   )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_UNSIGNED  n = 0; // non missing markers
  gs_FLOAT    di = 0; // sum of distances

  (void)n;
  
  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];
      
      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
      
      if ( !(  ( GS_ALLELE_MISSING == a ) ||
               ( GS_ALLELE_MISSING == b ) ||
               ( GS_ALLELE_MISSING == c ) ||
               ( GS_ALLELE_MISSING == d )    ) ) 
	{
	  n++;
	  
	  if ( !( ( (c==a)&&(d==b)) ||
		  ( (c==b)&&(d==a))    )  ) 
	    {
	      if ( ( (c==a)||(c==b)) ||
		   ( (d==a)||(d==b))    )  
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

  (gv->crs_S[k]).gd = (gs_FLOAT) sqrt( di );
  
}

void gs_calc_mu ( gs_varset_type *gv,
		  gs_UNSIGNED k )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d, e, count;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT   *mu = &((gv->crs_S[k]).mu); 

  *mu = gv->u[0];
 
  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];

      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
	  
      e = gv->EffAll[ 1+m ];

      count = (int)(e==a) + (int)(e==b) + (int)(e==c) + (int)(e==d);
      *mu  +=  (gs_FLOAT)count / (gs_FLOAT)2 * (gs_FLOAT)(gv->u[1+m]);

    }
}

void gs_cross_eval_mu_01 ( gs_varset_type* gv  ,
			   int* retval          )
{

  if ( NULL == gv->u ) {
    gs_info(gv,-2,"No estimated effects available");
    *retval = -2; return;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  gs_UNSIGNED k;

#ifdef MT
  #pragma omp parallel for 
#endif
  
  for ( k = 0; k <  gv->NCrs ; k++) 
    {
      gs_calc_mu (gv,k);
    }

  *retval = 0; return;
}

void gs_cross_eval_mu_01_GV ( int*    retval  ,
		              char**  set_name   )
{
  gs_cross_eval_mu_01 ( gs_fdta(*set_name) ,
		        retval              );
}

void gs_calc_ma ( gs_varset_type *gv,
		  gs_UNSIGNED k )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d, e, count;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT   *ma = &((gv->crs_S[k]).ma); 

  *ma = gv->u[0];

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];

      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
	  
      e = gv->EffAll[ 1+m ];

      count = (int)(e==a) + (int)(e==b) + (int)(e==c) + (int)(e==d);

      if ( ( 0 < count ) && (0 < gv->u[1+m]) )
	*ma  +=  (gs_FLOAT)2 * (gs_FLOAT)(gv->u[1+m]);
    }
}

void gs_cross_eval_ma_01 ( gs_varset_type* gv  ,
			   int* retval          )
{

  if ( NULL == gv->u ) {
    gs_info(gv,-2,"No estimated effects available");
    *retval = -2; return;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  gs_UNSIGNED k;

#ifdef MT
  #pragma omp parallel for 
#endif
  
  for ( k = 0; k <  gv->NCrs ; k++)
    {
      gs_calc_ma (gv,k);
    }

  *retval = 0; return;
}

void gs_cross_eval_ma_01_GV ( int*    retval  ,
		              char**  set_name   )
{
  gs_cross_eval_ma_01 ( gs_fdta(*set_name) ,
		        retval              );
}
void gs_calc_mi ( gs_varset_type *gv,
		  gs_UNSIGNED k )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d, e, count;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT   *mi = &((gv->crs_S[k]).mi); 

  *mi = gv->u[0];

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];

      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
	  
      e = gv->EffAll[ 1+m ];

      count = (int)(e==a) + (int)(e==b) + (int)(e==c) + (int)(e==d);

      if ( ( 0 < count ) && (0 > gv->u[1+m]) )
	*mi  +=  (gs_FLOAT)2 * (gs_FLOAT)(gv->u[1+m]);
    }
}

void gs_cross_eval_mi_01 ( gs_varset_type* gv  ,
			   int* retval          )
{

  if ( NULL == gv->u ) {
    gs_info(gv,-2,"No estimated effects available");
    *retval = -2; return;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  gs_UNSIGNED k;

#ifdef MT
  #pragma omp parallel for 
#endif
  
  for ( k = 0; k <  gv->NCrs ; k++)
    {
      gs_calc_mi (gv,k);
    }

  *retval = 0; return;
}

void gs_cross_eval_mi_01_GV ( int*    retval  ,
		              char**  set_name   )
{
  gs_cross_eval_mi_01 ( gs_fdta(*set_name) ,
		        retval              );
}

void gs_cross_eval_gd_01 ( gs_varset_type *gv  ,
			   char ** dist ,
			   int*   retval         )
{

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  gs_UNSIGNED k;

  if ( 0 == strcmp("mrd",*dist) ) 
    {
#ifdef MT
  #pragma omp parallel for 
#endif
  
      for ( k = 0; k <  gv->NCrs ; k++)  gs_calc_mrd (gv,k);
    }
  else if ( 0 == strcmp("rd",*dist) ) 
    {
#ifdef MT
  #pragma omp parallel for 
#endif
  
      for ( k = 0; k <  gv->NCrs ; k++)  gs_calc_rd (gv,k);
    }
  else
    {
#ifdef MT
  #pragma omp parallel for 
#endif
  
      for ( k = 0; k <  gv->NCrs ; k++)  gs_calc_euc (gv,k);
    }

  *retval = 0; return;
}

void gs_cross_eval_gd_01_GV ( char**  dist    ,
			      int*    retval  ,
			      char**  set_name   )
{
  gs_cross_eval_gd_01 ( gs_fdta(*set_name) ,
			dist,
			retval              );
}

void gs_calc_genv_01 (gs_varset_type* gv)

{
  gs_UNSIGNED ii, jj, idx;
  int a,b, e, count;

  gv->genv = realloc (gv->genv, gv->NoMar * gv->NoInd * sizeof (gs_FLOAT));

  for (ii=0; ii < gv->NoMar; ii++)
    {
      e = gv->EffAll[ 1+ii ];

      for (jj=0; jj < gv->NoInd; jj++)
	{
	  idx = ii*gv->NoInd + jj;
	  a = gv->All1[idx] ;
	  b = gv->All2[idx] ;
	  count = (int)(e==a) + (int)(e==b);
	  gv->genv[idx] = (gs_FLOAT)count * (gs_FLOAT)(gv->u[1+ii]) /2;
	}
    }
}


void gs_calc_va_linked ( gs_varset_type* gv,
			 gs_UNSIGNED k      )
{
  gs_UNSIGNED c, ii, jj, m, idx, idx_i, idx_j;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT q, ggu, ggv, hhu, hhv;

  gs_FLOAT* va = &((gv->crs_S[k]).va); 

  *va = 0;

  for ( c=0; c < gv->NChrom; c++ )
    {
       for  ( ii = 0 ; ii < gv->chrom[c].NLocC ; ii ++)
	 for  ( jj = 0 ; jj < gv->chrom[c].NLocC ; jj ++)
	   {
	     idx = ii * gv->chrom[c].NLocC + jj;
	     q = (gv->chrom[c].q) [idx]  ;

	     m = ii + gv->chrom[c].first_loc_i;
	     idx_i = m*gv->NoInd + i;
	     idx_j = m*gv->NoInd + j;
	     ggu   = gv->genv[idx_i];
	     hhu   = gv->genv[idx_j];  

	     m = jj + gv->chrom[c].first_loc_i;
	     idx_i = m*gv->NoInd + i;
	     idx_j = m*gv->NoInd + j;
	     ggv   = gv->genv[idx_i];
	     hhv   = gv->genv[idx_j];  
	     
	     *va += (q/2 - 0.25) * (ggu*ggv + hhu*hhv - ggu*hhv - hhu*ggv);
	   }
    }

  *va = *va * 4;

}

void gs_calc_genv_02 (gs_varset_type* gv)
{
  gs_UNSIGNED ii, jj, idx;
  int a,b, e, count1, count2 ;

  gv->genv1 = realloc (gv->genv1, gv->NoMar * gv->NoInd * sizeof (gs_FLOAT));
  gv->genv2 = realloc (gv->genv2, gv->NoMar * gv->NoInd * sizeof (gs_FLOAT));

  for (ii=0; ii < gv->NoMar; ii++)
    {
      e = gv->EffAll[ 1+ii ];

      for (jj=0; jj < gv->NoInd; jj++)
	{
	  idx = ii*gv->NoInd + jj;
	  a = gv->All1[idx] ;
	  b = gv->All2[idx] ;
	  count1 = (int)(e==a) ;
	  count2 = (int)(e==b) ;
	  gv->genv1[idx] = (gs_FLOAT)count1 * (gs_FLOAT)(gv->u[1+ii]);
	  gv->genv2[idx] = (gs_FLOAT)count2 * (gs_FLOAT)(gv->u[1+ii]);
	}
    }
}


void gs_calc_va_crs ( gs_varset_type* gv,
		       gs_UNSIGNED k      )
{
  gs_UNSIGNED c, ii, jj, m, idx, idx_i, idx_j;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT q, ggu, ggv, hhu, hhv;

  gs_FLOAT* va = &((gv->crs_S[k]).va); 

  *va = 0;

  for ( c=0; c < gv->NChrom; c++ )
    {
       for  ( ii = 0 ; ii < gv->chrom[c].NLocC ; ii ++)
	 for  ( jj = 0 ; jj < gv->chrom[c].NLocC ; jj ++)
	   {
	     idx = ii * gv->chrom[c].NLocC + jj;
	     q = (gv->chrom[c].q) [idx]  ;

	     m = ii + gv->chrom[c].first_loc_i;
	     idx_i = m*gv->NoInd + i;
	     idx_j = m*gv->NoInd + j;
	     ggu   = gv->genv1[idx_i];
	     hhu   = gv->genv2[idx_i];  

	     m = jj + gv->chrom[c].first_loc_i;
	     idx_i = m*gv->NoInd + i;
	     idx_j = m*gv->NoInd + j;
	     ggv   = gv->genv1[idx_i];
	     hhv   = gv->genv2[idx_i];  
	     
	     *va += (q/2 - 0.25) * (ggu*ggv + hhu*hhv - ggu*hhv - hhu*ggv);

	     m = ii + gv->chrom[c].first_loc_i;
	     idx_i = m*gv->NoInd + i;
	     idx_j = m*gv->NoInd + j;
	     ggu   = gv->genv1[idx_j];
	     hhu   = gv->genv2[idx_j];  

	     m = jj + gv->chrom[c].first_loc_i;
	     idx_i = m*gv->NoInd + i;
	     idx_j = m*gv->NoInd + j;
	     ggv   = gv->genv1[idx_j];
	     hhv   = gv->genv2[idx_j];  
	     
	     *va += (q/2 - 0.25) * (ggu*ggv + hhu*hhv - ggu*hhv - hhu*ggv);
	     
	   }
    }

   //*va = *va /2;

}



void gs_calc_va ( gs_varset_type* gv,
		  gs_UNSIGNED k )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d, e, count;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT  EZ2 = 0;
  gs_FLOAT  EZ  = 0;

  gs_FLOAT* va = &((gv->crs_S[k]).va); 

  *va = 0;

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];

      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
	  
      e = gv->EffAll[ 1+m ];

      count = (int)(e==a) + (int)(e==b) + (int)(e==c) + (int)(e==d);

      EZ2  = (gs_FLOAT)count * (gs_FLOAT)(gv->u[1+m]) * (gs_FLOAT)(gv->u[1+m]) / 4;
      EZ   = (gs_FLOAT)count * (gs_FLOAT)(gv->u[1+m])  / 4 ;

      *va += ( EZ2 - EZ*EZ ) * 4; // *4 both gametes = one entrire dh line

    }
}

void gs_cross_eval_va_01 ( gs_varset_type* gv  ,
			   char** pop_type     ,
			   int*   t            ,
			   char** map_function ,
			   int*   retval          )
{

  if ( NULL == gv->u ) {
    gs_info(gv,-2,"No estimated effects available");
    *retval = -2; return;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  char*empty=""; char ** Sout_filename = &empty;
  int auf = 0; int*Sauxfiles = &auf;

  gs_UNSIGNED k;

  if ( 0 == strcmp("unlinked",*pop_type) )
    {
#ifdef MT
  #pragma omp parallel for 
#endif
  
       for ( k = 0; k <  gv->NCrs ; k++)
	 {
	   gs_calc_va (gv,k);
	 }
    }

  else if ( ( 0 == strcmp(*pop_type,"DH")) ||
	    ( 0 == strcmp(*pop_type,"SSD"))  )
    {

      gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_calc_rf_01 ( gv, map_function, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_calc_q_01  ( gv, pop_type, t, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_calc_genv_01 (gv);

#ifdef MT
  #pragma omp parallel for 
#endif
  
       for ( k = 0; k <  gv->NCrs ; k++)
	 {
	   gs_calc_va_linked (gv,k);
	 }
    }

  else if ( 0 == strcmp(*pop_type,"CRS") )
    {

      gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_calc_rf_01 ( gv, map_function, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_calc_q_01  ( gv, pop_type, t, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_calc_genv_02 (gv);

#ifdef MT
  #pragma omp parallel for 
#endif
  
       for ( k = 0; k <  gv->NCrs ; k++)
	 {
	   gs_calc_va_crs (gv,k);
	 }
    }
  
  
  else {
    gs_info(gv,-2,"Unknown population type");
    *retval = -2; return;
  }


  *retval = 0; return;
}

void gs_cross_eval_va_01_GV ( char**  pop_type    ,
			      int*    t           ,
			      char**  map_function ,
			      int*    retval      ,
			      char**  set_name      )
{
  gs_cross_eval_va_01 ( gs_fdta(*set_name) ,
			pop_type           ,
			t                  ,
			map_function       ,
		        retval               );
}

void gs_calc_es ( gs_varset_type* gv ,
		  gs_UNSIGNED k      ,  
		  double* i            )
{
  gs_FLOAT* mu = &((gv->crs_S[k]).mu); 
  gs_FLOAT* va = &((gv->crs_S[k]).va); 
  gs_FLOAT* es = &((gv->crs_S[k]).es); 

  *es = *mu + *i * sqrt(*va) ;

}

void gs_cross_eval_es_01 ( gs_varset_type* gv  ,
			   double * i          ,
			   int* retval          )
{
  if ( (NULL == gv->crs) || (NULL == gv->crs_S) 
       || ( 0 == gv->crs[0]->va ) ) {
    gs_info(gv,-2,"No estimated variances available");
    *retval = -2; return;
  }

  snprintf(gv->msg,
         sizeof(gv->msg),
         "Standardized selection differential %5.2f",
         (float)*i);
 
  gs_info(gv,0,gv->msg);

  gs_UNSIGNED k;

#ifdef MT
  #pragma omp parallel for 
#endif
  
  for ( k = 0; k <  gv->NCrs ; k++)
    {
      gs_calc_es (gv,k,i);
    }

  *retval = 0; return;
}

void gs_cross_eval_es_01_GV ( double* i       ,
			      int*    retval  ,
			      char**  set_name   )
{
  gs_cross_eval_es_01 ( gs_fdta(*set_name) ,
			i                  ,
		        retval               );
}

int gs_cross_cmp_gd(const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->gd < (*ib)->gd )  return 1;
  if ( (*ia)->gd > (*ib)->gd )  return -1;
  return 0;

} 

int gs_cross_cmp_mu (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->mu < (*ib)->mu )  return 1;
  if ( (*ia)->mu > (*ib)->mu )  return -1;
  return 0;

} 

int gs_cross_cmp_ma (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->ma < (*ib)->ma )  return 1;
  if ( (*ia)->ma > (*ib)->ma )  return -1;
  return 0;
} 

int gs_cross_cmp_mi (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->ma < (*ib)->ma )  return -1;
  if ( (*ia)->ma > (*ib)->ma )  return 1;
  return 0;
} 

int gs_cross_cmp_va (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->va < (*ib)->va )  return 1;
  if ( (*ia)->va > (*ib)->va )  return -1;
  return 0;

} 

int gs_cross_cmp_es (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->es < (*ib)->es )  return 1;
  if ( (*ia)->es > (*ib)->es )  return -1;
  return 0;

} 

int gs_cross_cmp_index(const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->p1 > (*ib)->p1 )  return 1;
  else {
    if ( (*ia)->p1 < (*ib)->p1 )  return -1;
    else {
      if ( (*ia)->p2 > (*ib)->p2 )  return 1;
      else {
	if ( (*ia)->p1 < (*ib)->p1 )  return -1;
	else return 0;
      }
    }
  }

} 



void gs_cross_info_01 (  gs_varset_type *gv   ,
			 int   * bestn        ,
		         char ** sortby       ,
			 char**  out_filename ,
			 int*    retval         )
{
  gs_UNSIGNED i,m;
  FILE *outfp;
   
  if ( NULL == gv->crs ) {
      gs_info(gv,-2, "No cross information available");
      *retval = -2; return ;
  }

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

  if ( 0 == strcmp(*sortby,"index") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_index           );

  if ( 0 == strcmp(*sortby,"gd") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_gd           );

  if ( 0 == strcmp(*sortby,"mu") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_mu           );

  if ( 0 == strcmp(*sortby,"mi") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_mi           );

  if ( 0 == strcmp(*sortby,"ma") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_ma           );

  if ( 0 == strcmp(*sortby,"va") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_va           );

  if ( 0 == strcmp(*sortby,"es") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_es           );

  if (0 == *bestn) 
    m = gv->NCrs;
  else
    m = ( (gs_UNSIGNED)(*bestn) < gv->NCrs) ? (gs_UNSIGNED)(*bestn) : gv->NCrs;

  fprintf(outfp,"P1No P2No P1Name P2Name gd mu mi ma va es\n");
  for ( i = 0; i < m ; i++)
    fprintf(outfp,"%lu %lu %s %s %f %f %f %f %f %f\n",
	    1+(gv->crs[i])->p1,
	    1+(gv->crs[i])->p2,
	    gv->IndName[ (gv->crs[i])->p1 ],
	    gv->IndName[ (gv->crs[i])->p2 ],
	    (double)(gv->crs[i])->gd,
	    (double)(gv->crs[i])->mu,
	    (double)(gv->crs[i])->mi,
	    (double)(gv->crs[i])->ma,
	    (double)(gv->crs[i])->va,
	    (double)(gv->crs[i])->es
	    );
  
  fclose(outfp);  

  *retval = 0; 
  return ;
}

void gs_cross_info_01_GV (int   * bestn,
			  char ** sortby,
			  char**  out_filename,
			  int*    retval,
			  char**  set_name         )
{
  gs_cross_info_01 ( gs_fdta(*set_name),
		     bestn,
		     sortby,
		     out_filename,
		     retval         ) ;
}


extern void ps_get_simpop ( char ** pop_name   ,
			    gs_UNSIGNED * NoInd,
			    gs_UNSIGNED * NoMar,
			    char (**IndName) [gs_MAXNME] ,
			    char (**MarName) [gs_MAXNME] ,
			    GS_ALLELE_TYPE *(* All1)    ,
			    GS_ALLELE_TYPE *(* All2)    ,
			    gs_UNSIGNED       *  NLoc ,
			    gs_mappoint_type  ** map_S,
                            char              ** map_str_S,
                            size_t             * map_str_N,
			    gs_mappoint_type *** map  ,
			    gs_UNSIGNED       ** mdl  ,
			    int*    bg,
			    int* Sretval                  );

void st_get_simpop (  gs_varset_type *gv   ,
	              char ** pop_name,
		      int*    bg,
		      int*    retval         )
{
  gs_reset (gv);

  gs_UNSIGNED * NoInd = &(gv->NoInd);
  gs_UNSIGNED * NoMar = &(gv->NoMar);

  char   (**IndName) [gs_MAXNME]   = &(gv->IndName); 
  char   (**MarName) [gs_MAXNME]   = &(gv->MarName); 

  GS_ALLELE_TYPE *(* All1)  = &(gv->All1); 
  GS_ALLELE_TYPE *(* All2)  = &(gv->All2); 

  gs_UNSIGNED       *  NLoc      = &(gv->NLoc);
  gs_mappoint_type  ** map_S     = &(gv->map_S);
  char              ** map_str_S = &(gv->map_str_S);
  size_t             * map_str_N = &(gv->map_str_N);
  gs_mappoint_type *** map       = &(gv->map);
  gs_UNSIGNED       ** mdl       = &(gv->mdl);      

 int Srv=0; int* Sretval = &Srv;
 ps_get_simpop ( pop_name,
		 NoInd,
		 NoMar,
		 IndName,
		 MarName,
		 All1,
		 All2,
                 NLoc,
		 map_S,
                 map_str_S,
                 map_str_N,
		 map,
		 mdl,  
		 bg,
		 Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }
  
  /* Determine marker statistics */
  {
    int Srv=0; int* Sretval = &Srv; char*empty="";
    char ** ind_filename = &empty; char ** mar_filename = &empty;
    char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
    gs_marker_stats_01 ( gv ,
			 ind_filename,
			 mar_filename,
			 gen_filename,
			 Sauxfiles,Sauxfiles,Sauxfiles,
			 Sauxfiles,
			 Sretval        );
    if (-2 == *Sretval) { *retval = -2; return; }
  }


  *retval = 1;   return;
}

void st_get_simpop_GV (  char ** pop_name,
			 int*    bg,
			 int*    retval,
			 char** set_name         )
{
  st_get_simpop ( gs_fdta(*set_name),
		  pop_name,
		  bg,
		  retval   );     
}

//+DEBUG++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


extern void ps_get_simpop2 ( char ** pop_name   ,
			    gs_UNSIGNED * NoInd,
			    gs_UNSIGNED * NoMar,
			    char (**IndName) [gs_MAXNME] ,
			    char (**MarName) [gs_MAXNME] ,
			    GS_ALLELE_TYPE *(* All1)    ,
			    GS_ALLELE_TYPE *(* All2)    ,
			    gs_UNSIGNED       *  NLoc ,
			    gs_mappoint_type  ** map_S,
			    gs_mappoint_type *** map  ,
			    gs_UNSIGNED       ** mdl  ,
			    int*    bg,
			    int* Sretval                  );

void st_get_simpop2 (  gs_varset_type *gv   ,
	               char ** pop_name,
		       int*    bg,
		       int*    retval         )
{
  gs_reset (gv);

  gs_UNSIGNED * NoInd = &(gv->NoInd);
  gs_UNSIGNED * NoMar = &(gv->NoMar);

  char   (**IndName) [gs_MAXNME]   = &(gv->IndName); 
  char   (**MarName) [gs_MAXNME]   = &(gv->MarName); 

  GS_ALLELE_TYPE *(* All1)  = &(gv->All1); 
  GS_ALLELE_TYPE *(* All2)  = &(gv->All2); 

  gs_UNSIGNED       *  NLoc  =  &(gv->NLoc);   
  gs_mappoint_type  ** map_S =  &(gv->map_S);  
  gs_mappoint_type *** map   =  &(gv->map);    
  gs_UNSIGNED       ** mdl   =  &(gv->mdl);      

 int Srv=0; int* Sretval = &Srv;
 ps_get_simpop2 ( pop_name,
		 NoInd,
		 NoMar,
		 IndName,
		 MarName,
		 All1,
		 All2,
                 NLoc, 
		 map_S,
		 map,  
		 mdl,  
		 bg,
		 Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }
  
  /* Determine marker statistics */
  {
    int Srv=0; int* Sretval = &Srv; char*empty="";
    char ** ind_filename = &empty; char ** mar_filename = &empty;
    char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
    gs_marker_stats_01 ( gv ,
			 ind_filename,
			 mar_filename,
			 gen_filename,
			 Sauxfiles,Sauxfiles,Sauxfiles,
			 Sauxfiles,
			 Sretval        );
    if (-2 == *Sretval) { *retval = -2; return; }
  }


  *retval = 1;   return;
}

void st_get_simpop2_GV (  char ** pop_name,
			 int*    bg,
			 int*    retval,
			 char** set_name         )
{
  st_get_simpop2 ( gs_fdta(*set_name),
		   pop_name,
		   bg,
		   retval   );     
}

//+DEBUG++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




void il_get_population_size (  gs_varset_type *gv ,
			       int*   NoLines       )
{
  *NoLines = (int) gv->NoInd;
}
				
void il_get_population_size_GV ( int*   NoLines,
				 char** set_name         )
{
  il_get_population_size ( gs_fdta(*set_name),
			   NoLines           ) ;
}


void il_calc_segment_length ( gs_varset_type *gv )
{
  gs_UNSIGNED l,c,m;
  gs_UNSIGNED i1, i2, i3;
  gs_FLOAT    p1, p2, p3;

  gv->dimu = gv->NLoc;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );

  for ( c=0,l=0 ; c < gv->NChrom ; c ++ ) 
    {
      m=0;
      i2 =  gv->chrom[c].first_loc_i + m     ; p2 =  gv->map[i2]->pos;
      i3 =  gv->chrom[c].first_loc_i + m + 1 ; p3 =  gv->map[i3]->pos;
      gv->u[l] = (p3-p2)/2 ; 
      l++;

      for ( m=1 ; m < (-1 + gv->chrom[c].NLocC) ; m++,l++)
	{
	  i1 =  gv->chrom[c].first_loc_i + m - 1 ; p1 =  gv->map[i1]->pos;
	  i2 =  gv->chrom[c].first_loc_i + m     ; p2 =  gv->map[i2]->pos;
	  i3 =  gv->chrom[c].first_loc_i + m + 1 ; p3 =  gv->map[i3]->pos;
	  gv->u[l] = (p2-p1)/2 + (p3-p2)/2; 
	}

      i1 =  gv->chrom[c].first_loc_i + m - 1 ; p1 =  gv->map[i1]->pos;
      i2 =  gv->chrom[c].first_loc_i + m     ; p2 =  gv->map[i2]->pos;
      gv->u[l] = (p2-p1)/2 ; 
      l++;
    }

}

void il_mark_loci  ( gs_varset_type *gv ,
		     int*    chrom      ,
		     double* lower      ,
		     double* upper      ,
		     int*    exclude       )
{

  gs_UNSIGNED l,c,m;
  gs_UNSIGNED i2;
  gs_FLOAT    p2;

  int hit = 0;

  if ( 0 == *chrom ) return ;

  for ( c=0,l=0 ; c < gv->NChrom ; c ++ ) 
    {
      for ( m=0 ; m <  gv->chrom[c].NLocC ; m++,l++)
	{
	  i2 =  gv->chrom[c].first_loc_i + m     ; p2 =  gv->map[i2]->pos;
	  
	  hit = ( ( c == (gs_UNSIGNED) (-1 + *chrom )) &&
		  ( (0 == *upper) || 
		    ( ( *lower <= p2) && ( p2 <= *upper ) ) ) );

	  if (1 == *exclude) hit = !hit;
		  
	  if ( !hit ) gv->u[l] = 0;
	}
    }
}


double il_calc_tot (gs_varset_type *gv)
{
  gs_UNSIGNED l,m;

//  gs_FLOAT    p2;

  int cc, s_c, e_c;
  s_c = 0;
  e_c = (int)gv->NChrom;

  double tot = 0;

  for ( cc = s_c ; cc < e_c ; cc ++ ) 
    {
      for ( m=0 ; m < gv->chrom[cc].NLocC ; m++)
	{
	  l =  gv->chrom[cc].first_loc_i + m ;
//	  p2 = gv->map[l]->pos;
	  tot += gv->u[l];
	}
    }
  tot = tot * 2;
  return (tot);
}

void il_eval_lines  (  gs_varset_type *gv ,
		       int*   allele      ,
		       int*   chrom       ,
		       double* lower      ,
		       double* upper      ,
		       int*   exclude     ,
		       double* abs        ,
		       double* tot        ,
		       double* nseg       ,
		       double* lseg       ,
		       int*   retval        )
{
  gs_UNSIGNED i,l,m,idx;
  gs_UNSIGNED l1,l2,idx1,idx2;

  //  gs_FLOAT    p2;

  int cc, s_c, e_c;

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available");
    *retval = -2; return ;
  }

  if ((*allele < 1) || (*allele > GS_ALLELE_MAX))
    {
      gs_info(gv,-2,"Allele code must be in the GS allele range");
      *retval = -2; return;
    }

  gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  s_c = 0;
  e_c = (int)gv->NChrom;

  if ( ( *chrom <0 ) || ( *chrom > (int)gv->NChrom )) {
    gs_info(gv,-2,"Incorrect chromosome number");
    *retval = -2; 
    return; 
    }

  for ( i=0 ; i < gv->NoInd; i ++ ) 
    {
      
      il_calc_segment_length ( gv );

      il_mark_loci  ( gv           ,
		      &(chrom[i])  ,
		      &(lower[i])  ,
		      &(upper[i])  ,
		      exclude  );

      tot[i] = il_calc_tot (gv);

      abs[i] = 0;

      // Allele frequencies
      for ( cc = s_c ; cc < e_c ; cc ++ ) 
	{
	  for ( m=0 ; m < gv->chrom[cc].NLocC ; m++)
	    {
	      l =  gv->chrom[cc].first_loc_i + m ;
//	      p2 = gv->map[l]->pos;
	      idx = l*gv->NoInd + i;
	      abs[i] +=
		(gv->All1[idx] == *allele) * gv->u[l]  +
		(gv->All2[idx] == *allele) * gv->u[l];

	    }
	}

      // Segments
      nseg[i] = 0;
      for ( cc = s_c ; cc < e_c ; cc ++ ) 
	{
	  m = 0;
	  l2 = gv->chrom[cc].first_loc_i + m ;
	  idx2 = l2*gv->NoInd + i;
	  if  (0 != gv->u[l2]) 
	    {
	       if  (*allele == gv->All1[idx2]) nseg[i] ++ ;
	       if  (*allele == gv->All2[idx2]) nseg[i] ++ ;
	    }

	  for ( m=1 ; m < gv->chrom[cc].NLocC ; m++)
	    {
	      l2 = gv->chrom[cc].first_loc_i + m ;
	      l1 = l2-1;

	      idx1 = l1*gv->NoInd + i;
	      idx2 = l2*gv->NoInd + i;

	      if  ( (0 != gv->u[l1]) && ((0 != gv->u[l2])))
		{
		  if  ( (gv->All1[idx1] != gv->All1[idx2])  &&
			(*allele == gv->All1[idx2])             )
		    nseg[i] ++ ;

		  if  ( (gv->All2[idx1] != gv->All2[idx2])  && 
			(*allele == gv->All2[idx2])             )
		    nseg[i] ++ ;
		} 
	      else if ( (0 == gv->u[l1]) && ((0 != gv->u[l2])))
		{
		  if  (*allele == gv->All1[idx2]) nseg[i] ++ ;
		  if  (*allele == gv->All2[idx2]) nseg[i] ++ ;
		}
	    }
	}
      if ( 0 != nseg[i] )
	lseg[i] = abs[i] / nseg[i];

    } // for i = 0 .. ind

  *retval=0;
}
				
void il_eval_lines_GV ( int*   allele  ,
			int*   chrom   ,
			double* lower  ,
			double* upper  ,
			int* exclude   ,
			double* abs    ,
			double* tot    , 
			double* nseg   ,
			double* lseg   , 
			int*   retval  ,
			char** set_name  )
{
  il_eval_lines ( gs_fdta(*set_name),
		  allele            ,
		  chrom             ,
		  lower             ,
		  upper             ,
		  exclude           ,
		  abs               ,
		  tot               ,
		  nseg              ,
		  lseg              ,
		  retval              ) ;
}

void il_eval_library (  gs_varset_type *gv ,
			double* cov    ,
			double* dep    ,
			double* seglib  ,
			int   * allele ,
			int*    retval       )
{
  gs_UNSIGNED i,j,m,cc,l1,l2,idx,idx1,idx2,cnt;
  double tot;

  if ((*allele < 1) || (*allele > GS_ALLELE_MAX))
    {
      gs_info(gv,-2,"Allele code must be in the GS allele range");
      *retval = -2; return;
    }

  gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  il_calc_segment_length ( gv );

  /* Genome coverage */

  tot = 0;
  *cov = 0;
  for (i=0; i < gv->NoMar; i++)
    {
      tot += gv->u[i] ;
      for (j=0; j < gv->NoInd; j++)
	{
	  idx = i*gv->NoInd + j;
	  if ( (*allele == gv->All1[idx] ) ||
	       (*allele == gv->All2[idx] )    )
	    {
	      *cov += gv->u[i];
	      break;
	    }
	}
      }
  *cov = *cov/tot;

  /* Depth of coverage */

  gs_UNSIGNED *depth =  malloc( gv->NoMar * sizeof(gs_UNSIGNED) );
  if ( NULL == depth ) ERR_M;

  for (i=0; i < gv->NoMar; i++)
    {
      depth[i] = 0;

      for (j=0; j < gv->NoInd; j++)
	{
	  idx = i*gv->NoInd + j;
	  if ( (*allele == gv->All1[idx] ) ||
	       (*allele == gv->All2[idx] )    )
	    {
	      depth[i] ++;
	    }
	}
      }

  *dep = 0;
  cnt  = 0;
  for (i=0; i < gv->NoMar; i++)
    {
      if ( depth[i] > 0 )
	{
	  *dep += (double) depth[i];
	  cnt ++;
	}
    }
  *dep = *dep / (double) cnt;
  
  free(depth);


  /* Number of segments in the library */
  
  *seglib = (double)gv->NChrom ;

  for ( cc = 0 ; cc < gv->NChrom ; cc ++ ) 
    for ( m=1 ; m < gv->chrom[cc].NLocC ; m++)
      for (i=0; i < gv->NoInd; i++)  
	{
	  l2 = gv->chrom[cc].first_loc_i + m ;
	  l1 = l2-1;

	  idx1 = l1*gv->NoInd + i;
	  idx2 = l2*gv->NoInd + i;

	  if (  (gv->All1[idx1] != gv->All1[idx2]) ||
		(gv->All2[idx1] != gv->All2[idx2])    )
	    {
	      (*seglib)++ ;
	      break;
	    }
	}
} 

				
void il_eval_library_GV ( double* cov    ,
			  double* dep    ,
			  double* seglib , 
			  int*    allele ,
			  int*    retval , 
			  char**  set_name  )
{
  il_eval_library ( gs_fdta(*set_name),
		    cov    ,
		    dep    ,
		    seglib ,
		    allele ,
		    retval) ;
}


//##############################################################################
//# BEGIN: Ridge regression after transformation                         2026-02
//##############################################################################

void gs_build_GI_01(gs_varset_type* gv, int* retval)
{
    const double small = 1e-5;

    gs_UNSIGNED ii, jj;
    double d;
    double det;

    const gs_UNSIGNED n = gv->NoInd;
    const gs_UNSIGNED m = gv->NoMar;

    /* ensure storage */
    gv->t9 = realloc(gv->t9, (size_t)n * (size_t)n * sizeof(double));
    gv->GI = realloc(gv->GI, (size_t)n * (size_t)n * sizeof(double));
    gv->u0 = realloc(gv->u0, (size_t)n * sizeof(gs_UNSIGNED));
    if ((NULL == gv->t9) || (NULL == gv->GI) || (NULL == gv->u0)) ERR_M;

    /* G = Z Z'  */
    st_dsymm_ZZt(n, m, gv->Z, gv->t9);


    /* Keep unmodified copy for fallback */
    memcpy(gv->GI, gv->t9, (size_t)n * (size_t)n * sizeof(double));

    /* -------------------------------
       FAST PATH: Cholesky (SPD exploit)
       ------------------------------- */
    {
        /* Work in-place on gv->t9 */
        int info = st_dpotrf_rm_sym(gv->t9, (int)n);
        if (info == 0) {
            info = st_dpotri_rm_sym(gv->t9, (int)n);
            if (info == 0) {

		/* Inverse lives in LOWER triangle of A_rm */
		symmetrize_from_lower_rm(gv->t9, (int)n);

                /* store result into GI */
                memcpy(gv->GI, gv->t9, (size_t)n * (size_t)n * sizeof(double));

                *retval = 0;
                return;
            }
            /* If inversion failed, fall through to eigen fallback */
        }
        /* If factorization failed (not SPD / rank-deficient), fall through */
    }

    /* -------------------------------------------------------
       FALLBACK: keep your existing LU heuristic + eigen ridge
       ------------------------------------------------------- */

    /* Try LU-based inverse; keep it as secondary fast attempt */
    memcpy(gv->t9, gv->GI, (size_t)n * (size_t)n * sizeof(double));
    st_dgetrf(gv->t9, n, gv->u0, &d);

    for (ii = 0, det = 1.0; ii < n; ++ii) {
        if (fabs(gv->t9[(size_t)ii * (size_t)n + ii]) < 1e-10) {
            det = 0.0;
            break;
        }
    }

    if (det != 0.0) {
        st_dgetri(gv->t9, n, gv->u0, gv->GI);
        *retval = 0;
        return;
    }

    /* Eigen ridge approximation  */
    gs_info(gv, -1, "ZZ' has not full rank, using approximation)");

    double *w = (double*)malloc((size_t)n * sizeof(double));
    double *z = (double*)malloc((size_t)n * (size_t)n * sizeof(double));
    if (!w || !z) {
        free(z); free(w);
        *retval = -1;
        return;
    }

    memcpy(gv->t9, gv->GI, (size_t)n * (size_t)n * sizeof(double));
    st_dsyevr_v(gv->t9, w, z, n);

    /* t9 = Z * diag(1/max(w,small)) (scale columns) */
    for (jj = 0; jj < n; ++jj) {
        double denom = (w[jj] < small) ? small : w[jj];
        double dinv  = 1.0 / denom;
        for (ii = 0; ii < n; ++ii) {
            gv->t9[(size_t)ii * (size_t)n + jj] =
                z[(size_t)ii * (size_t)n + jj] * dinv;
        }
    }

    /* GI = t9 %*% t(z) */
    st_dgemm("N","T",
             n, n,
             n, n,
             gv->t9, z, gv->GI);

    free(z);
    free(w);

    *retval = 0;
    return;
}

//##############################################################################
//# BEGIN: Ridge regression after transformation                         2015-01
//##############################################################################

void gs_build_GI_02 ( gs_varset_type* gv ,
		      int*   retval        )
{
  // Inverse of the G matrix: G =ZZ'


  double small = 1e-5;

  //double zero  = 1e-10;
  
  gv->t9  = realloc (gv->t9, gv->NoInd * gv->NoInd *sizeof(double) );
  gv->GI  = realloc (gv->GI, gv->NoInd * gv->NoInd *sizeof(double) );
  if (( NULL == gv->t9 )||( NULL == gv->GI )) ERR_M;

  st_dgemm ("N","T",
	    gv->NoInd, gv->NoMar,  // ZZ'
	    gv->NoInd, gv->NoMar,
	    gv->Z, gv->Z, gv->t9);

  double d;
  gs_UNSIGNED ii,jj; double det=0;
  gv->u0  = realloc ( gv->u0, gv->NoInd * sizeof(gs_UNSIGNED) );

  // LU decomposition
  st_dgetrf ( gv->t9, gv->NoInd, gv->u0, &d);

  // determinant is product of diagonal elements
  // full rank: determinant != 0

  for (ii=0,det=1; ii < gv->NoInd; ii++) 
    if  (   fabs (gv->t9 [ii*gv->NoInd + ii]) < 1e-10)
      det = 0.0;

  if ( 0.0 == det ) {
    gs_info(gv,-1, "ZZ' has not full rank, using approximation)");

    double* w = malloc (gv->NoInd * sizeof(double));
    double* z = malloc (gv->NoInd * gv->NoInd * sizeof(double));
    
    double* t6  = malloc (gv->NoInd * gv->NoInd * sizeof(double) );
    double* t7  = malloc (gv->NoInd * gv->NoInd * sizeof(double) );
    double* t8  = malloc (gv->NoInd * gv->NoInd * sizeof(double) );
    double* t9  = malloc (gv->NoInd * gv->NoInd * sizeof(double) );

    st_dgemm ("N","T",
	      gv->NoInd, gv->NoMar,  
	      gv->NoInd, gv->NoMar,
	      gv->Z, gv->Z, t9);
    
    while ( 0.0 == det )
      // See Joriani et al 2003. J Dairy Sci 86:677-679
      // and J Hered 107:755
      {

	memcpy ( gv->t9, t9, gv->NoInd * gv->NoInd * sizeof(double));
	st_dsyevr_v ( gv->t9, w, z, gv->NoInd );

	for (ii=0; ii < gv->NoInd; ii++)    // delta_n
	  for (jj=0; jj < gv->NoInd; jj++)
	    t6[ii*gv->NoInd + jj] = 0;

	for (ii=0; ii < gv->NoInd; ii++)
	  t6[ii*gv->NoInd + ii] = (w[ii] < small) ? small : w[ii];

	st_dgemm ("N","N",                // U_n * delta_n
		  gv->NoInd, gv->NoInd,  
		  gv->NoInd, gv->NoInd,
		  z, t6, t7);

	st_dgemm ("N","T",                // (U_n * delta_n) * t(U_n)
		  gv->NoInd, gv->NoInd,  
		  gv->NoInd, gv->NoInd,
		  t7, z, t8);

	for (ii=0; ii < gv->NoInd; ii++)   
	  for (jj=0; jj < gv->NoInd; jj++)
	    t9[ii*gv->NoInd + jj] = 
	      t9[ii*gv->NoInd + jj] -
	      (t9[ii*gv->NoInd + jj] - t8[ii*gv->NoInd + jj]);

	memcpy ( gv->t9, t9, gv->NoInd * gv->NoInd * sizeof(double));
	st_dgetrf ( gv->t9, gv->NoInd, gv->u0, &d);

	for (ii=0,det=1; ii < gv->NoInd; ii++) 
	  if  (   fabs (gv->t9 [ii*gv->NoInd + ii]) < 1e-10)
	    det = 0.0;

      }

      free(w); free(z);
      free(t6); free(t7); free(t8); free(t9);

    }    

  st_dgetri ( gv->t9, gv->NoInd, gv->u0 ,gv->GI );

  *retval = 0; return ;
}


/* ------------------------------------------------------------------ */
/* .Call wrapper for NEW implementation gs_build_GI_01                  */
/* Symbol: call_gs_build_GI_01_R                                       */
/* ------------------------------------------------------------------ */

/* Helper: check length safely */
static void check_len_rowmajor(SEXP x, R_xlen_t need, const char *name) {
    if (!isReal(x)) error("%s must be a double vector", name);
    if (XLENGTH(x) < need) error("%s too short: need %lld, got %lld",
                                name, (long long)need, (long long)XLENGTH(x));
}

SEXP call_gs_build_GI_01_R(SEXP Z_rm_, SEXP NoInd_, SEXP NoMar_)
{
    if (!isInteger(NoInd_) || LENGTH(NoInd_) != 1) error("NoInd must be integer scalar");
    if (!isInteger(NoMar_) || LENGTH(NoMar_) != 1) error("NoMar must be integer scalar");

    int NoInd = asInteger(NoInd_);
    int NoMar = asInteger(NoMar_);
    if (NoInd <= 0 || NoMar <= 0) error("NoInd/NoMar must be > 0");

    R_xlen_t lenZ = (R_xlen_t)NoInd * (R_xlen_t)NoMar;
    check_len_rowmajor(Z_rm_, lenZ, "Z_rm");

    /* Build a minimal gv that gs_build_GI_01 needs */
    gs_varset_type gv;
    memset(&gv, 0, sizeof(gv));

    gv.NoInd = (gs_UNSIGNED)NoInd;
    gv.NoMar = (gs_UNSIGNED)NoMar;

    /* Allocate and copy Z (row-major) */
    gv.Z = (double*)malloc((size_t)lenZ * sizeof(double));
    if (!gv.Z) error("Out of memory allocating gv.Z");
    memcpy(gv.Z, REAL(Z_rm_), (size_t)lenZ * sizeof(double));

    /* Ensure these start NULL so realloc() inside gs_build_GI_01 works */
    gv.t9 = NULL;
    gv.GI = NULL;
    gv.u0 = NULL;

    int retval = -999;
    gs_build_GI_01(&gv, &retval);
    if (retval != 0) {
        free(gv.Z);
        free(gv.t9);
        free(gv.GI);
        free(gv.u0);
        error("gs_build_GI_01 failed (retval=%d)", retval);
    }

    /* Return GI as row-major vector */
    R_xlen_t lenGI = (R_xlen_t)NoInd * (R_xlen_t)NoInd;
    SEXP out = PROTECT(allocVector(REALSXP, lenGI));
    memcpy(REAL(out), gv.GI, (size_t)lenGI * sizeof(double));

    /* Cleanup */
    free(gv.Z);
    free(gv.t9);
    free(gv.GI);
    free(gv.u0);

    UNPROTECT(1);
    return out;
}

/* ------------------------------------------------------------------ */
/* .Call wrapper for OLD implementation gs_build_GI_02                  */
/* Symbol: call_gs_build_GI_02_R                                       */
/* ------------------------------------------------------------------ */

SEXP call_gs_build_GI_02_R(SEXP Z_rm_, SEXP NoInd_, SEXP NoMar_)
{
    if (!isInteger(NoInd_) || LENGTH(NoInd_) != 1) error("NoInd must be integer scalar");
    if (!isInteger(NoMar_) || LENGTH(NoMar_) != 1) error("NoMar must be integer scalar");

    int NoInd = asInteger(NoInd_);
    int NoMar = asInteger(NoMar_);
    if (NoInd <= 0 || NoMar <= 0) error("NoInd/NoMar must be > 0");

    R_xlen_t lenZ = (R_xlen_t)NoInd * (R_xlen_t)NoMar;
    check_len_rowmajor(Z_rm_, lenZ, "Z_rm");

    gs_varset_type gv;
    memset(&gv, 0, sizeof(gv));

    gv.NoInd = (gs_UNSIGNED)NoInd;
    gv.NoMar = (gs_UNSIGNED)NoMar;

    gv.Z = (double*)malloc((size_t)lenZ * sizeof(double));
    if (!gv.Z) error("Out of memory allocating gv.Z");
    memcpy(gv.Z, REAL(Z_rm_), (size_t)lenZ * sizeof(double));

    gv.t9 = NULL;
    gv.GI = NULL;
    gv.u0 = NULL;

    int retval = -999;
    gs_build_GI_02(&gv, &retval);
    if (retval != 0) {
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




void gs_mmet_coeff_01 ( gs_varset_type* gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{
  // Constant lambda
  // Calculates GI vor the back transformation

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  if ( NULL == gv->y ) {
      gs_info(gv,-2,
        "No performance data available, run 'read.performance.data' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
      gs_info(gv,-2,
        "No Z matrix available, run 'build.Z' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factor defined");
      *retval = -2; return ;
  }

  /* Allocate memory for the A matrix */

  gv->dimA = 1 + gv->NoInd;
  gv->A =  realloc( gv->A, gv->dimA * gv->dimA * sizeof(double) );
  if ( NULL == gv->A ) ERR_M;

  for ( i = 0; i < gv->dimA ; i++) 
    for ( j = 0; j < gv->dimA ; j++) 
	gv->A [ i*gv->dimA + j ] =  0;

  /* Init part 1*/
  
  gv->A [0] = (double) gv->NoInd;

  /* Init part 2 and 3*/

  for ( j = 0; j < gv->NoInd; j++)  
    {
      gv->A [     0*gv->dimA + (1+j) ] = 1;
      gv->A [ (j+1)*gv->dimA + 0     ] = 1;
    }

  /* Init part 4*/

  gs_build_GI_01(gv,Sretval); // Also used in mmet_solve_01 for back transf.
  if (-2 == *Sretval) { *retval = -2; return; }

  for ( i = 0; i < gv->NoInd; i++)
    for ( j = 0; j < gv->NoInd; j++)
	gv->A [ (1+i)*gv->dimA + (1+j) ] = 
                                      gv->GI[ i*gv->NoInd + j] * gv->ladi[0];

  for ( i = 0; i < gv->NoInd; i++)
    gv->A [ (1+i)*gv->dimA + (1+i) ] += 1;


  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }

      fclose(outfp);

    }

  *retval = 0; return ;
}

void gs_mmet_coeff_01_GV ( char** out_filename,
			  int*   auxfiles, 
			  int*   retval,
			  char** set_name                )
{
  gs_mmet_coeff_01 (gs_fdta(*set_name),
		    out_filename,
		    auxfiles, 
		    retval        ) ;
}

void gs_mmet_coeff_02 ( gs_varset_type* gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{
  // Coefficient matrix for gs_lambda_emstep_tc_01. Constant shrinkage
  // factor. For use iterations, re-uses a stored inv(G) matrix, 
  // calculated by gs_build_GI_01

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  if ( NULL == gv->y ) {
      gs_info(gv,-2,
        "No performance data available, run 'read.performance.data' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
      gs_info(gv,-2,
        "No Z matrix available, run 'build.Z' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factor defined");
      *retval = -2; return ;
  }

  /* Allocate memory for the A matrix */

  gv->dimA = 1 + gv->NoInd;
  gv->A =  realloc( gv->A, gv->dimA * gv->dimA * sizeof(double) );
  if ( NULL == gv->A ) ERR_M;
 
  for ( i = 0; i < gv->dimA ; i++) 
    for ( j = 0; j < gv->dimA ; j++) 
	gv->A [ i*gv->dimA + j ] =  0;

  /* Init part 1*/
  
  gv->A [0] = (double) gv->NoInd;

  /* Init part 2 and 3*/

  for ( j = 0; j < gv->NoInd; j++)  
    {
      gv->A [     0*gv->dimA + (1+j) ] = 1;
      gv->A [ (j+1)*gv->dimA + 0     ] = 1;
    }

  /* Init part 4*/
  
  // Build A 

  for ( i = 0; i < gv->NoInd; i++)
    for ( j = 0; j < gv->NoInd; j++)
	gv->A [ (1+i)*gv->dimA + (1+j) ] = 
	                gv->ladi[0] * gv->GI[ i*gv->NoInd + j];

  for ( i = 0; i < gv->NoInd; i++)
    gv->A [ (1+i)*gv->dimA + (1+i) ] += 1;


  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_mmet_coeff_02_GV ( char** out_filename,
			  int*   auxfiles, 
			  int*   retval,
			  char** set_name                )
{
  gs_mmet_coeff_02 (gs_fdta(*set_name),
		    out_filename,
		    auxfiles, 
		    retval        ) ;
}


void gs_mmet_rhs_01 ( gs_varset_type *gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{
  gs_UNSIGNED i,j, dim_b;
  FILE *outfp;

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available, run 'build.Z' first");
    *retval = -2; return ;
  }

  /* Allocate memory for b  (Ax=b) */
  dim_b = 1 + gv->NoInd;
  gv->b =  realloc( gv->b, dim_b * sizeof(double) );
  if ( NULL == gv->b ) ERR_M;

  /* Init part 1*/
  
  gv->b[0] = 0;
  for ( j = 0; j < gv->NoInd; j++) 
    gv->b[0] += gv->y[j];

  /* Init part 2*/
  
  for ( j = 0; j < gv->NoInd; j++)  
    gv->b [ 1+j ] = gv->y[j];

  /* Write b vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( i = 0; i < dim_b ; i++)
	fprintf(outfp,"%f \n",gv->b [i] );
      
      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mmet_rhs_01_GV ( char** out_filename,
	 		int*   auxfiles, 
	 		int*   retval,
	 		char** set_name    )
{
  gs_mmet_rhs_01 ( gs_fdta(*set_name) ,
		  out_filename,
		  auxfiles, 
		  retval        ) ;
}

void gs_mmet_solve_01 ( gs_varset_type *gv ,
		      char** out_filename,
		      int*   auxfiles, 
		      int*   retval        )
{
  gs_UNSIGNED i;
  FILE *outfp;

  int Srv=0; int* Sretval = &Srv;

  if ( NULL == gv->A ) {
      gs_info(gv,-2, "No coefficient matrix available");
      *retval = -2; return ;
  }

  if ( NULL == gv->b )  {
      gs_info(gv,-2, "No right hand side available");
      *retval = -2; return ;
  }

  /* Allocate memory for the solution vector of the transformed mode*/
  gv->a =  realloc( gv->a, (1+gv->NoInd) * sizeof(double) );
  if ( NULL == gv->a ) ERR_M;

  /* Solution vector for the transformed data*/
  st_dgesv ( gv->dimA, gv->dimA, gv->A, gv->a, gv->b, Sretval );
  if (*Sretval==-2) {*retval = -2; return ;}

  /* Allocate memory for the solution vector of marker effects*/
  gv->dimu = 1 + gv->NoMar;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );
  if ( NULL == gv->u ) ERR_M;

  double *tmp1;
  tmp1 = malloc ( gv->NoInd * gv->NoMar * sizeof(double) );
  if ( NULL == tmp1 ) ERR_M;

  st_dgemm ("T","N",
	    gv->NoInd, gv->NoMar,  // Z' GSI
	    gv->NoInd, gv->NoInd,
	    gv->Z, gv->GI, tmp1);
 
  st_dgemm ("N","N",
	    gv->NoMar, gv->NoInd,  // (Z' GSI ) a
	    gv->NoInd,         1,
	    tmp1, 1+gv->a, 1+gv->u);

  gv->u[0] = gv->a[0];

  free(gv->A); gv->A =NULL;
  free(tmp1);

  /* Write solution vector to output file */
  if (1 == *auxfiles) {
   
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"effect\n");
      for ( i = 0; i < gv->dimu ; i++)
	fprintf(outfp,"%s.%i %f \n",
		gv->EffNme[i],
		gv->EffAll[i],
		gv->u [i]);

      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mmet_solve_01_GV ( char** out_filename ,
			  int*   auxfiles     , 
			  int*   retval       ,
			  char** set_name      )
{
  gs_mmet_solve_01 ( gs_fdta(*set_name) ,
		    out_filename,
		    auxfiles, 
	            retval        ) ;
}


void gs_lambda_rrblup_01 ( gs_varset_type* gv , int* maxiter, double*
		 	  precision, double* hsq, char** out_filename,
		 	  int* auxfiles, int* retval );

void gs_lambda_rmla_02 ( gs_varset_type *gv , double* alpha , int*
			   maxiter , double* precision , double* hsq ,
			   char** out_filename , int* auxfiles , int*
			   retval );

void gs_lambda_rmlc_01 ( gs_varset_type *gv , int* maxiter, double*
			  precision, double* hsq, char** out_filename,
			 int* auxfiles, int* retval );

void gs_mme_coeff_01 ( gs_varset_type* gv , char** out_filename, int*
		       auxfiles, int* retval );

void gs_lambda_rmla_01 ( gs_varset_type *gv , double* alpha , int*
			  maxiter , double* precision , double* hsq ,
			  char** out_filename , int* auxfiles , int*
			  retval );

void gs_mme_rhs_01 ( gs_varset_type *gv , char** out_filename, int*
		     auxfiles, int* retval ); 


void gs_mme_solve_01 ( gs_varset_type *gv , char** out_filename, int*
		       auxfiles, int* retval );

void gs_lambda_rml_02 ( gs_varset_type *gv , double* alpha , int*
			  maxiter , double* precision , double* hsq ,
			  char** out_filename , int* auxfiles , int*
			 retval );

void gs_lambda_rmlv_01 ( gs_varset_type *gv , int* maxiter, double*
			  precision, double* hsq, char** out_filename,
			  int* auxfiles, int* retval );

void gs_vc_rrblup_01 ( gs_varset_type* gv , int* maxiter, double*
		       precision, double* hsq, char** out_filename,
		       int* auxfiles, int* retval );

void gs_w_aov_01 ( gs_varset_type *gv , double* alpha , char**
		   out_filename , int* auxfiles , int* retval );

void gs_w_I_01 ( gs_varset_type *gv , char** out_filename , int*
		 auxfiles , int* retval );

void gs_single_marker_aov_01 ( gs_varset_type *gv , double* alpha ,
			       char** out_filename , int* auxfiles ,
			       int* retval );


void gs_esteff_rr_01 (  gs_varset_type *gv ,
		        char**  methode,
                        int*    maxiter,
		        double* precision,
		        double* hsq,
			double* alpha,
		        char**  out_filename,
		        int*    auxfiles, 
		        int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  if (strcmp (*methode, "BLUP")  == 0) 
    {

      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_vc_rrblup_01 ( gv, maxiter, precision, hsq,
			nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_rrblup_01 ( gv, maxiter, precision, hsq,
			   nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    } 

    else if (strcmp (*methode, "RMLA") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_vc_rrblup_01 ( gv, maxiter, precision, hsq,
			nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_lambda_rrblup_01 ( gv, maxiter, precision, hsq,
			   nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_single_marker_aov_01 ( gv, alpha, nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_lambda_rmla_02 ( gv, alpha, maxiter, precision, hsq,
			   nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  else if (strcmp (*methode, "rmlc") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_rmlc_01 ( gv,
			  maxiter,precision,hsq,
		      nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  else if (strcmp (*methode, "rmla") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_rmla_01 ( gv,
			  alpha,
			  maxiter,precision,hsq,
			  nixp,neinp,Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  else if (strcmp (*methode, "rmlv") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_rmlv_01 ( gv,
			  maxiter,precision,hsq,
			  nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  
  else 
    { 
      gs_info(gv,-2,"Unknown estimation method");
      *retval = -2;  return;
     }

  *retval = 0; return ;
}

void gs_esteff_rr_01_GV (  char**  methode,
			   int*    maxiter,
		           double* precision,
			   double* hsq,
			   double* alpha,
		           char**  out_filename,
		           int*    auxfiles, 
		           int*    retval,
		           char** set_name         )
{
  gs_esteff_rr_01 ( gs_fdta(*set_name),
                    methode,
		    maxiter,   
	            precision, 
		    hsq,
		    alpha,
	            out_filename,
	            auxfiles, 
	            retval         ) ;
}

void gs_cross_validation_02 ( gs_varset_type *gv   ,
                              char ** methode       ,
                              int*    n_estimation ,
                              int*    n_runs       ,
                              double* hsq          ,
                              double* alpha        ,
                              int*    maxiter      ,
                              double* precision    ,
                              char**  out_filename_r ,
                              char**  out_filename_u ,
                              int*    auxfiles     ,
                              int*    retval       )
{
  int ntr = 1;
  int prev_info_level;
  int no_info[1] = {-2};
  int CVNoRuns;
  int parallel_error = 0;

  double *CVcor = NULL;
  gs_FLOAT *CVest = NULL;

  size_t no_mar;
  size_t cv_cols;
  size_t cv_count;

  *retval = -2;

  if (gv->NoInd > (gs_UNSIGNED)INT_MAX)
    {
      gs_info(gv,-2,"More individuals than INT_MAX");
      return;
    }

  if ((0 == gv->NoInd) || (0 == gv->NoMar))
    {
      gs_info(gv,-2,"No marker data available");
      return;
    }

  if (NULL == gv->y)
    {
      gs_info(gv,-2,"No performance data available");
      return;
    }

  if ((*n_estimation < 1) ||
      (*n_estimation >= (int)gv->NoInd))
    {
      gs_info(gv,-2,"Invalid size of training set");
      return;
    }

  if (*n_runs < 1)
    {
      gs_info(gv,-2,"Number of cross-validation runs must be positive");
      return;
    }

       if (strcmp (*methode, "BLUP") == 0) {}
  else if (strcmp (*methode, "RMLA") == 0) {}
  else if (strcmp (*methode, "rmlc") == 0) {}
  else if (strcmp (*methode, "rmlv") == 0) {}
  else if (strcmp (*methode, "rmla") == 0) {}
  else
    {
      gs_info(gv,-2,"Unknown estimation method");
      return;
    }

  CVNoRuns = *n_runs;

  no_mar = (size_t)gv->NoMar;
  if (((gs_UNSIGNED)no_mar != gv->NoMar) ||
      (no_mar == ((size_t)-1)))
    {
      gs_info(gv,-2,"Cross-validation result size is too large");
      return;
    }

  cv_cols = no_mar + 1;

  if ((size_t)CVNoRuns > ((size_t)-1) / cv_cols)
    {
      gs_info(gv,-2,"Cross-validation result size is too large");
      return;
    }

  cv_count = (size_t)CVNoRuns * cv_cols;

  if (((size_t)CVNoRuns > ((size_t)-1) / sizeof(*CVcor)) ||
      (cv_count > ((size_t)-1) / sizeof(*CVest)))
    {
      gs_info(gv,-2,"Cross-validation result size is too large");
      return;
    }

  CVcor = malloc((size_t)CVNoRuns * sizeof(*CVcor));
  CVest = malloc(cv_count * sizeof(*CVest));

  if ((NULL == CVcor) || (NULL == CVest))
    {
      free(CVcor);
      free(CVest);
      gs_info(gv,-2,"Memory allocation error");
      return;
    }

  for (int c=0; c < CVNoRuns; c++)
    CVcor[c] = 0;

#ifdef MT
  #pragma omp parallel
  {
    #pragma omp single
    {
      ntr = omp_get_num_threads();
    }
  }
#endif

  if (ntr < 1)
    ntr = 1;

#if defined(use_openblas) && defined(MT)
  {
    int eins = 1;
    openblas_set_num_threads(eins);
    goto_set_num_threads(eins);
    omp_set_num_threads(ntr);
  }
#endif

#if defined(use_mkl) && defined(MT)
  {
    int eins = 1;
    mkl_set_num_threads(eins);
    omp_set_num_threads(ntr);
  }
#endif

  st_get_info_level(&prev_info_level);
  st_set_info_level(no_info);

#ifdef MT
  #pragma omp parallel reduction(|:parallel_error)
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
    #pragma omp for
#endif
    for (int r=0; r < CVNoRuns; r++)
      {
        int success_rep = 0;
        int no_unsucc_rep = 0;

        while ((0 == success_rep) &&
               (no_unsucc_rep < 1000))
          {
            *Sretval = 0;

            gs_build_esvs_01(gv,es,vs,n_estimation,
                             nixp,nixp,nixp,nixp,neinp,Sretval);

            if (-2 == *Sretval)
              {
                parallel_error = 1;
                break;
              }

            *Sretval = 0;

            gs_esteff_rr_01(es,methode,maxiter,precision,hsq,alpha,
                            nixp,neinp,Sretval);

            if (-2 == *Sretval)
              {
                no_unsucc_rep++;
                continue;
              }

            success_rep = 1;
          }

        if (parallel_error)
          continue;

        if (0 == success_rep)
          {
            parallel_error = 1;
            continue;
          }

        if ((NULL == es->u) ||
            ((size_t)es->dimu != cv_cols))
          {
            parallel_error = 1;
            continue;
          }

        for (i=0; i < es->dimu; i++)
          CVest[(size_t)r * cv_cols + (size_t)i] = es->u[i];

        *Sretval = 0;

        gs_estimate_bv_01(es,vs,nixp,neinp,Sretval);

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

        for (i=0; i < vs->NoInd; i++)
          {
            sx  += vs->est[i];
            sxx += vs->est[i] * vs->est[i];
            sy  += vs->y[i];
            syy += vs->y[i] * vs->y[i];
            sxy += vs->est[i] * vs->y[i];
          }

        CVcor[r] =
          (k*sxy - sx*sy) /
          sqrt((k*sxx - sx*sx) * (k*syy - sy*sy));
      }

    gs_reset(es);
    gs_reset(vs);
  } // parallel

#if defined(use_openblas) && defined(MT)
  openblas_set_num_threads(ntr);
  goto_set_num_threads(ntr);
#endif

#if defined(use_mkl) && defined(MT)
  mkl_set_num_threads(ntr);
#endif

  if (parallel_error)
    {
      st_set_info_level(&prev_info_level);
      free(CVcor);
      free(CVest);
      gs_info(gv,-2,"Cross validation failed in parallel calculation");
      return;
    }

  {
    int Srv = 0;
    int *Sretval = &Srv;
    char *meth_ = "BLUP";
    char **meth = &meth_;

    gs_esteff_rr_01(gv,meth,maxiter,precision,hsq,alpha,
                    nixp,neinp,Sretval);

    if (-2 == *Sretval)
      {
        st_set_info_level(&prev_info_level);
        free(CVcor);
        free(CVest);
        *retval = -2;
        return;
      }
  }

  if ((NULL == gv->u) ||
      ((size_t)gv->dimu != cv_cols))
    {
      st_set_info_level(&prev_info_level);
      free(CVcor);
      free(CVest);
      gs_info(gv,-2,"Inconsistent number of effects in cross validation");
      return;
    }

  for (gs_UNSIGNED i=0; i < gv->dimu; i++)
    {
      gs_FLOAT sum = 0;

      for (int r=0; r < CVNoRuns; r++)
        sum += CVest[(size_t)r * cv_cols + (size_t)i];

      gv->u[i] = (double)(sum / (double)CVNoRuns);
    }

  st_set_info_level(&prev_info_level);

  if (1 == *auxfiles)
    {
      FILE *fp;

      fp = fopen(*out_filename_r,"w");
      if (NULL == fp)
        {
          free(CVcor);
          free(CVest);
          gs_info(gv,-2,"File open error");
          return;
        }

      fprintf(fp,"cor \n");
      for (int c=0; c < CVNoRuns; c++)
        fprintf(fp,"%i %f\n",c,CVcor[c]);
      fclose(fp);

      fp = fopen(*out_filename_u,"w");
      if (NULL == fp)
        {
          free(CVcor);
          free(CVest);
          gs_info(gv,-2,"File open error");
          return;
        }

      fprintf(fp,"effect\n");
      for (gs_UNSIGNED i=0; i < gv->dimu; i++)
        fprintf(fp,"%s.%i %f \n",
                gv->EffNme[i],
                gv->EffAll[i],
                gv->u[i]);
      fclose(fp);
    }

  free(CVcor);
  free(CVest);

  *retval = 0;
} 


void gs_cross_validation_02_GV ( char ** scheme       ,
				 int*    n_estimation ,
				 int*    n_runs       , 
				 double* hsq          , 
				 double* alpha        ,
				 int*    maxiter      ,
				 double* precision    ,
				 char**  out_filename_r , 
				 char**  out_filename_u , 
				 int*    auxfiles     , 
				 int*    retval       ,
				 char** set_name      )
 {
   gs_cross_validation_02 ( gs_fdta(*set_name),
			    scheme       ,
			    n_estimation ,
			    n_runs       , 
			    hsq          , 
			    alpha        ,
			    maxiter      ,
			    precision    ,
			    out_filename_r , 
			    out_filename_u , 
			    auxfiles     , 
			    retval           ) ;
 }

//##############################################################################
//# END: Ridge regression after transformation                           2015-01
//##############################################################################

//##############################################################################
//# BEGIN: LD based haplotypes                                           2015-04
//##############################################################################

// Calculation of pairwise LD for all loci on a chromosome

gs_FLOAT gs_ld_rsquare ( gs_varset_type *gv ,
			 gs_UNSIGNED ii, 
			 gs_UNSIGNED jj      )
{
  gs_UNSIGNED i,k,l,idx_all1,idx_all2;

  gs_UNSIGNED NoAllii = gv->NoAll[ii];
  gs_UNSIGNED NoAlljj = gv->NoAll[jj];
  GS_ALLELE_TYPE *Allii = gv->All[ii];
  GS_ALLELE_TYPE *Alljj = gv->All[jj];
  gs_UNSIGNED *AllCntii = gv->AllCnt[ii];
  gs_UNSIGNED *AllCntjj = gv->AllCnt[jj];

  if ((0 < NoAllii) &&
      (GS_ALLELE_MISSING == Allii[NoAllii-1])) NoAllii--;
  if ((0 < NoAlljj) &&
      (GS_ALLELE_MISSING == Alljj[NoAlljj-1])) NoAlljj--;

  gs_UNSIGNED CNT;
  double N = 2 * (double)gv->NoInd;
  double r2 = 0;

  for (k = 0; k < NoAllii; k++)          // All combinations of alleles at loci
    for (l = 0; l < NoAlljj; l++)        // ii and jj  
      {
	CNT = 0;
	for ( i = 0; i < gv->NoInd ; i++) {
	  idx_all1 = ii*gv->NoInd + i ;
	  idx_all2 = jj*gv->NoInd + i ;
	  if ( ( gv->All1[idx_all1] == Allii[k]) &&
	       ( gv->All1[idx_all2] == Alljj[l])     )  CNT++;
	  if ( ( gv->All2[idx_all1] == Allii[k]) &&
	       ( gv->All2[idx_all2] == Alljj[l])     )  CNT++;
	}
	double fx = (double)AllCntii[k]/N;
	double fy = (double)AllCntjj[l]/N;
	double Dxy = (double)CNT / N  - fx*fy;
	double r2xy = Dxy*Dxy / (fx * (1-fx) * fy * (1-fy) );
	r2 += r2xy * fx *fy ;
      }
  return  ( r2 );

}

gs_FLOAT gs_ld_dprime ( gs_varset_type *gv ,
			 gs_UNSIGNED ii, 
			 gs_UNSIGNED jj      )
{
  gs_UNSIGNED i,k,l,idx_all1,idx_all2;

  gs_UNSIGNED NoAllii = gv->NoAll[ii];
  gs_UNSIGNED NoAlljj = gv->NoAll[jj];
  GS_ALLELE_TYPE *Allii = gv->All[ii];
  GS_ALLELE_TYPE *Alljj = gv->All[jj];
  gs_UNSIGNED *AllCntii = gv->AllCnt[ii];
  gs_UNSIGNED *AllCntjj = gv->AllCnt[jj];

  if ((0 < NoAllii) &&
      (GS_ALLELE_MISSING == Allii[NoAllii-1])) NoAllii--;
  if ((0 < NoAlljj) &&
      (GS_ALLELE_MISSING == Alljj[NoAlljj-1])) NoAlljj--;

  gs_UNSIGNED CNT;
  double N = 2 * (double)gv->NoInd;
  double Dp = 0;
  double fx,fy,Dxy,c1,c2,Dmax,Dpxy;

  for (k = 0; k < NoAllii; k++)          // All combinations of alleles at loci
    for (l = 0; l < NoAlljj; l++)        // ii and jj  
      {
	CNT = 0;
	for ( i = 0; i < gv->NoInd ; i++) {
	  idx_all1 = ii*gv->NoInd + i ;
	  idx_all2 = jj*gv->NoInd + i ;
	  if ( ( gv->All1[idx_all1] == Allii[k]) &&
	       ( gv->All1[idx_all2] == Alljj[l])     )  CNT++;
	  if ( ( gv->All2[idx_all1] == Allii[k]) &&
	       ( gv->All2[idx_all2] == Alljj[l])     )  CNT++;
	}
	fx = (double)(AllCntii[k])/N;
	fy = (double)(AllCntjj[l])/N;
	Dxy = (double)CNT / N  - fx*fy;
	if (Dxy<0) { c1 = fx * fy; c2 = (1-fx) * (1-fy); }
	else       { c1 = fx * (1-fy); c2 = (1-fx) * fy; }
	if (c1 < c2) Dmax = c1; else Dmax = c2;
	Dpxy = fabs( Dxy / Dmax ) ;
	Dp += Dpxy * fx *fy ;
      }
  return  ( Dp );

}


/* ------------------------------------------------------------------------- */
/* LD from unphased biallelic genotypes.                                    */
/*                                                                           */
/* The two phase-estimating LD modes use a native two-locus maximum-         */
/* likelihood calculation.  For each marker pair, the nine observable       */
/* diploid genotype classes are counted without using the order of All1 and  */
/* All2 as phase information.  The ambiguous double heterozygotes are split  */
/* between the two possible diplotypes by an expectation-maximization (EM)   */
/* algorithm under Hardy-Weinberg equilibrium.  Several deterministic        */
/* starting values are tried and the converged solution with the greatest   */
/* observed-data likelihood is retained.                                     */
/* ------------------------------------------------------------------------- */

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


static int gs_ld_unphased_dosage(GS_ALLELE_TYPE a1,
                                 GS_ALLELE_TYPE a2,
                                 const GS_ALLELE_TYPE allele[2])
{
  int d = 0;

  if ((GS_ALLELE_MISSING == a1) || (GS_ALLELE_MISSING == a2))
    return -1;

  if ((a1 != allele[0]) && (a1 != allele[1]))
    return -2;
  if ((a2 != allele[0]) && (a2 != allele[1]))
    return -2;

  if (a1 == allele[1]) d++;
  if (a2 == allele[1]) d++;

  return d;
}


static int gs_ld_unphased_table(const gs_varset_type *gv,
                                gs_UNSIGNED ii,
                                gs_UNSIGNED jj,
                                gs_UNSIGNED table[9])
{
  GS_ALLELE_TYPE ai[2];
  GS_ALLELE_TYPE aj[2];

  if (!gs_ld_unphased_marker_alleles(gv, ii, ai) ||
      !gs_ld_unphased_marker_alleles(gv, jj, aj))
    return 0;

  for (int k = 0; k < 9; k++)
    table[k] = 0;

  for (gs_UNSIGNED ind = 0; ind < gv->NoInd; ind++)
    {
      size_t xi = (size_t)ii * (size_t)gv->NoInd + (size_t)ind;
      size_t xj = (size_t)jj * (size_t)gv->NoInd + (size_t)ind;
      int di = gs_ld_unphased_dosage(gv->All1[xi], gv->All2[xi], ai);
      int dj = gs_ld_unphased_dosage(gv->All1[xj], gv->All2[xj], aj);

      /* Missing genotypes are omitted pairwise. */
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
                             double initial_coupling,
                             double h[4],
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

  /* Haplotype copies contributed by genotype classes with known phase. */
  fixed[0] = 2.0 * (double)table[0] + (double)table[1] + (double)table[3];
  fixed[1] = 2.0 * (double)table[2] + (double)table[1] + (double)table[5];
  fixed[2] = 2.0 * (double)table[6] + (double)table[3] + (double)table[7];
  fixed[3] = 2.0 * (double)table[8] + (double)table[5] + (double)table[7];

  if (p < 0.0) p = 0.0;
  if (p > 1.0) p = 1.0;

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
          if (d > delta) delta = d;
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


static int gs_ld_unphased_haplotypes(const gs_UNSIGNED table[9],
                                     double h[4])
{
  /* Starting values describe the fraction of double heterozygotes initially
     assigned to the coupling diplotype.  Multiple starts protect against
     stationary boundary or saddle-point solutions. */
  static const double starts[] =
    {0.0, 0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99, 1.0};
  double best_ll = -HUGE_VAL;
  double best_h[4] = {0.0, 0.0, 0.0, 0.0};
  int found = 0;

  for (size_t s = 0; s < sizeof(starts)/sizeof(starts[0]); s++)
    {
      double candidate[4];
      double ll;

      if (!gs_ld_unphased_em(table, starts[s], candidate, &ll))
        continue;

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

  return 1;
}


static int gs_ld_unphased_measures(const gs_varset_type *gv,
                                   gs_UNSIGNED ii,
                                   gs_UNSIGNED jj,
                                   double *r2,
                                   double *dprime)
{
  gs_UNSIGNED table[9];
  double h[4];
  double pa0, pa1, pb0, pb1;
  double D, denom, dmax;

  *r2 = 0.0;
  *dprime = 0.0;

  if (!gs_ld_unphased_table(gv, ii, jj, table))
    return 0;
  if (!gs_ld_unphased_haplotypes(table, h))
    return 0;

  pa0 = h[0] + h[1];
  pa1 = h[2] + h[3];
  pb0 = h[0] + h[2];
  pb1 = h[1] + h[3];

  denom = pa0 * pa1 * pb0 * pb1;
  if (!(denom > 0.0))
    return 1;

  D = h[0] - pa0 * pb0;
  *r2 = D * D / denom;

  if (D < 0.0)
    dmax = fmin(pa0 * pb0, pa1 * pb1);
  else
    dmax = fmin(pa0 * pb1, pa1 * pb0);

  if (dmax > 0.0)
    *dprime = fabs(D) / dmax;

  /* Protect against harmless floating-point overshoot at the boundaries. */
  if (*r2 < 0.0) *r2 = 0.0;
  if (*r2 > 1.0 && *r2 < 1.0 + 1.0e-12) *r2 = 1.0;
  if (*dprime < 0.0) *dprime = 0.0;
  if (*dprime > 1.0 && *dprime < 1.0 + 1.0e-12) *dprime = 1.0;

  return 1;
}


static gs_FLOAT gs_ld_rsquare_unphased(gs_varset_type *gv,
                                       gs_UNSIGNED ii,
                                       gs_UNSIGNED jj)
{
  double r2, dprime;

  if (!gs_ld_unphased_measures(gv, ii, jj, &r2, &dprime))
    return 0.0;

  return (gs_FLOAT)r2;
}


static gs_FLOAT gs_ld_dprime_unphased(gs_varset_type *gv,
                                      gs_UNSIGNED ii,
                                      gs_UNSIGNED jj)
{
  double r2, dprime;

  if (!gs_ld_unphased_measures(gv, ii, jj, &r2, &dprime))
    return 0.0;

  return (gs_FLOAT)dprime;
}

void gs_calc_ld_01 ( gs_varset_type *gv  ,
		     char** ld_measure   ,
		     char** out_filename ,
		     int* auxfiles       ,
		     int* retval           )

{
  gs_UNSIGNED c, i, j, idx;
  gs_UNSIGNED ii,jj;
  gs_mappoint_type ** p;
  int Srv=0; int* Sretval = &Srv;

  gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }


  gs_FLOAT (*gs_ld) (gs_varset_type *gv, 
		     gs_UNSIGNED i, gs_UNSIGNED j );

  if ( 0 == strcmp(*ld_measure,"r2"))
    {
      gs_ld = gs_ld_rsquare;
    }

  else if ( 0 == strcmp(*ld_measure,"r2-estimate-phases"))
    {
      if (gs_ld_unphased_data_are_biallelic(gv))
        gs_ld = gs_ld_rsquare_unphased;
      else
        {
          gs_info(gv,-2,"Phase-estimating LD requires exactly two biological alleles at every marker");
          *retval = -2; return;
        }
    }

  else if ( 0 == strcmp(*ld_measure,"Dp"))
    {
      gs_ld = gs_ld_dprime;
    }

  else if ( 0 == strcmp(*ld_measure,"Dp-estimate-phases"))
    {
      if (gs_ld_unphased_data_are_biallelic(gv))
        gs_ld = gs_ld_dprime_unphased;
      else
        {
          gs_info(gv,-2,"Phase-estimating LD requires exactly two biological alleles at every marker");
          *retval = -2; return;
        }
    }

  else 
    {
      gs_info(gv,-2,"LD measure not (yet) implemented");
      *retval = -2; return;
    }

  /* Allocate memory for ld measure  */
  for (c = 0; c < gv->NChrom; c ++ ) 
    {
      gv->chrom[c].ld = realloc (gv->chrom[c].ld, 
	gv->chrom[c].NLocC * gv->chrom[c].NLocC * sizeof(gs_FLOAT) );
    }

  for ( c=0; c < gv->NChrom; c++ )
    {
       p =  gv->chrom[c].first_loc ;
#ifdef MT
       #pragma omp parallel for private(i,j,ii,jj,idx)  
#endif
       for  ( i = 0 ; i < ( gv->chrom[c].NLocC - 1 ) ; i ++)
	 for  ( j = i+1 ; j < gv->chrom[c].NLocC ; j ++)
	   {
	     ii = (p[i])->nr;
	     jj = (p[j])->nr;
	     idx = i * gv->chrom[c].NLocC + j;
	     (gv->chrom[c].ld) [idx] = gs_ld (gv,ii,jj);
	   }
    }

  if (1 == *auxfiles)
    {
      gs_UNSIGNED loc_rgt, loc_lft ;

      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( c=0; c < gv->NChrom; c++ )
	{
	  for  ( i = 0 ; i < ( gv->chrom[c].NLocC - 1 ) ; i ++)
	    for  ( j = i+1 ; j < gv->chrom[c].NLocC ; j ++)
	      {

		idx = i * gv->chrom[c].NLocC + j;

		loc_lft = gv->chrom[c].first_loc_i + i ;  // Index on the map
		loc_rgt = gv->chrom[c].first_loc_i + j ;  

		/*
		sprintf( gv->msg,"%lu %lu %lu %s %s %f \n",
			 1+c,
			 1+i,
			 1+j,
			 gv->map[loc_lft]->name,
			 gv->map[loc_rgt]->name,
			 (double) (gv->chrom[c].ld) [idx] 
			 );
		*/

		snprintf(gv->msg, sizeof gv->msg,
			 "%lu %lu %lu %s %s %f \n",
			 1UL + (unsigned long)c,
			 1UL + (unsigned long)i,
			 1UL + (unsigned long)j,
			 gv->map[loc_lft]->name,
			 gv->map[loc_rgt]->name,
			 (double)(gv->chrom[c].ld)[idx]);
		
		fprintf(outfp,"%s",gv->msg);

	      }
	  
	}

      fclose(outfp);
    }

  *retval=0; return;

 }

void gs_calc_ld_01_GV (	char** ld_measure ,
			char** out_filename ,
			int*   auxfiles     ,
			int*   retval       ,
			char** set_name       )
{
  gs_calc_ld_01 ( gs_fdta(*set_name),
		  ld_measure        ,
		  out_filename      ,
		  auxfiles          ,
		  retval              ) ;
}

void gs_set_hblocks_01 ( gs_varset_type* gv , 
			 int* n_hap,
			 int* blk,
			 char** hap_symbol,
			 char** out_filename,
			 int*   auxfiles,
			 int* retval        )
/* Determining haplotype borders */
{

  if ( 0 == gv->NLoc ) 
    {
      snprintf(gv->msg, sizeof(gv->msg),
         "No linkage map has been loaded");
      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }

  if ( (gs_UNSIGNED)(*n_hap) != gv->NLoc ) 
    {
      snprintf(gv->msg, sizeof(gv->msg),
         "Incorrect number of markers: %lu - %lu",
         (gs_UNSIGNED)(*n_hap), gv->NLoc);

      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }


  //gs_UNSIGNED i,l,u,o,n, fib, lib;

  gs_UNSIGNED i,l, fib;

  gv->h_u   = realloc ( gv->h_u , gv->NLoc*sizeof(gs_UNSIGNED) )  ;
  gv->h_o   = realloc ( gv->h_o , gv->NLoc*sizeof(gs_UNSIGNED) )  ;
  gv->h_n   = realloc ( gv->h_n , gv->NLoc*sizeof(gs_UNSIGNED) )  ;
  if ( (NULL == gv->h_u)||(NULL == gv->h_o)||(gv->h_n == 0) ) ERR_M;

  gv->NHap  = 0;                    

  // Make a list of haplotype blocks from the assignment list of
  // markers to blocks

  gv->NHap = 0;
  fib=0;

  for ( i=1; i<gv->NLoc; i++)
    {
      if  (blk[i-1] != blk[i])
	{
	  gv->h_u [gv->NHap] = fib;
	  gv->h_o [gv->NHap] = i-1;
	  gv->h_n [gv->NHap] = i-fib;
	  gv->NHap++;
	  fib=i;
	}
    }
  gv->h_u [gv->NHap] = fib;
  gv->h_o [gv->NHap] = i-1;
  gv->h_n [gv->NHap] = i-fib;
  gv->NHap++;
  
  
  // Determine haploblock map
  {
    int map_rc = gs_build_haplotype_map(gv,hap_symbol);

    if (0 != map_rc)
      {
        gv->NHap = 0;
        GS_FREE(gv->Hap_S);
        GS_FREE(gv->Hap_str_S);
        gv->Hap_str_N = 0;

        if (-1 == map_rc)
          ERR_M;

        gs_info(gv,-2,"Haplotype marker name too long");
        *retval = -2;
        return;
      }
  }

  if (1 == *auxfiles){   
    // Write individual names to .hap file
    // New Version 2017-08-31

    FILE *outfp;
    if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

    for ( l=0; l < gv->NHap; l++)
      {
	char* lne = (char *) malloc ( gs_MAXLNE * sizeof(char) );
	snprintf(lne, gs_MAXLNE, "%i %f %s %s ",
         gv->Hap_S[l].chrom,
         gv->Hap_S[l].pos,
         gv->Hap_S[l].name,
         gv->Hap_S[l].class);

	for ( i = gv->h_u[l] ; i <= gv->h_o[l] ; i++ ) 
	  {
	    if ( gv->mdl[i] < gs_UNSIGNED_MAX )
	      {
		size_t a=strlen(lne);
		size_t b=strlen(gv->MarName[ gv->mdl[i]]);
		lne = realloc (lne, 2+a+b );
		strcat( lne, gv->MarName[ gv->mdl[i] ] );
		strcat( lne, ";" );
	      }
	  }

	fprintf(outfp,"%s\n",lne);
	free (lne);
      }
    
    fclose(outfp);
  }
  
}

void gs_set_hblocks_01_GV ( int*   n_hap,
			    int*   hap_list,
			    char** hap_symbol,
			    char** out_filename,
			    int*   auxfiles,    
			    int*   retval,
			    char** set_name)
{
  gs_set_hblocks_01 ( gs_fdta(*set_name),
		      n_hap,
		      hap_list,
		      hap_symbol,
		      out_filename,
		      auxfiles,
		      retval            ) ;
}


//##############################################################################
//# END:   LD based haplotypes                                           2015-04
//##############################################################################


//##############################################################################
//# Begin:   Genetic distances for factoials                             2015-11
//##############################################################################

void gs_cross_eval_gd_fct_01 ( gs_varset_type *gv  ,
		               char ** dist,
                               int*    split,
			       int*    retval)
{

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  gs_UNSIGNED sp = (gs_UNSIGNED)*split;

  if ( (sp < 1) || ( sp > (gv->NoInd)-1 ) ) {
    gs_info(gv,-2,"Incorrect split value");
    *retval = -2; return ;
  }

  gs_UNSIGNED i,j,k,l;
  
  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  for ( k = 0; k <  gv->NCrs ; k++)  gv->crs[k]->gd = 0;

  // gs_UNSIGNED fact[ sp * (gv->NoInd-sp) ];

  gs_UNSIGNED* fact = malloc( sp * (gv->NoInd-sp) * sizeof(gs_UNSIGNED) );
  
  k = 0;  l=0;
  for ( i=0; i < (gv->NoInd-1) ; i++ )  
    for ( j = i+1; j < gv->NoInd ; j++ )
      {
	if ( ( i < sp ) && (j >= sp ) )
	  {
	    fact[l]=k;
	    l++;
	  }
	k++;
      }

  
  if ( 0 == strcmp("mrd",*dist) ) 
    {
      #ifdef MT
  #pragma omp parallel for 
#endif
  
      for ( gs_UNSIGNED l = 0; l < (sp * (gv->NoInd-sp) ) ; l++) gs_calc_mrd (gv,fact[l]);
    }
  else if ( 0 == strcmp("rd",*dist) ) 
    {
      #ifdef MT
  #pragma omp parallel for 
#endif
  
      for ( gs_UNSIGNED l = 0; l < (sp * (gv->NoInd-sp) ) ; l++) gs_calc_rd (gv,fact[l]);
    }
  else
    {
      #ifdef MT
  #pragma omp parallel for 
#endif
  
      for ( gs_UNSIGNED l = 0; l < (sp * (gv->NoInd-sp) ) ; l++) gs_calc_euc (gv,fact[l]);
    }

  free(fact);
  *retval = 0; return;
}

void gs_cross_eval_gd_fct_01_GV ( char**  dist    ,
			          int*    split   ,
			          int*    retval  ,
			          char**  set_name   )
{
  gs_cross_eval_gd_fct_01 ( gs_fdta(*set_name) ,
			    dist,
			    split,
			    retval              );
}

//##############################################################################
//# End:  Genetic distances for factoials                                2015-11
//##############################################################################

//##############################################################################
//# Begin: Recode reference for two parents                              2017-09
//##############################################################################

void gs_recode_ref_02 ( gs_varset_type* gv , 
			int* r1,
			int* r2,
			int* descending,
			int* retval        )
/* Recoding to a reference indidivudal on individual marker basis */
{
  gs_UNSIGNED  i, j, idx, idx_R1, idx_R2;
  int a,b,c,d,x,y;

  if (( *r1 < 1 ) || ( *r2 < 1 ) ) {
      gs_info(gv,-2, "Invalid reference individual");
      *retval = -2; return ;
  }

  GS_ALLELE_TYPE *Ref1 =
    malloc(gv->NoMar * gv->NoInd * sizeof(*Ref1));
  GS_ALLELE_TYPE *Ref2 =
    malloc(gv->NoMar * gv->NoInd * sizeof(*Ref2));
  if ( (NULL == Ref1)  || (NULL == Ref2) ) ERR_M;
  
  gs_UNSIGNED R1 = (gs_UNSIGNED) *r1 - 1;
  gs_UNSIGNED R2 = (gs_UNSIGNED) *r2 - 1;
      
  // Recode alleles in local memory
  for (i=0; i < gv->NoMar; i++)
    {
      idx_R1 = i*gv->NoInd + R1;
      a = gv->All1[idx_R1];
      b = gv->All2[idx_R1];
      
      idx_R2 = i*gv->NoInd + R2;
      x = gv->All1[idx_R2];
      y = gv->All2[idx_R2];
      
      for (j=0; j <gv->NoInd; j++)
	{
	  idx = i*gv->NoInd + j;
	  c = gv->All1[idx];
	  d = gv->All2[idx];

          if ((GS_ALLELE_MISSING == a) ||
              (GS_ALLELE_MISSING == b) ||
              (GS_ALLELE_MISSING == c) ||
              (GS_ALLELE_MISSING == d))
            {
              Ref1[idx] = GS_ALLELE_MISSING;
              Ref2[idx] = GS_ALLELE_MISSING;
            } 
	  else
	    {
	      if (0==*descending) {

		if (( (c==a)&&(d==b)) ||
		    ( (c==b)&&(d==a))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 1;
		  } 
		else if (
	            ( (c==x)&&(d==y)) ||
		    ( (c==y)&&(d==x))    ) 
		  {
		    Ref1[idx] = 2;
		    Ref2[idx] = 2;
		  } 
		else if (
			 ( (c==a)||(c==b)) ||
			 ( (d==a)||(d==b))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 3;
		  } 
		else if (
			 ( (c==x)||(c==y)) ||
			 ( (d==x)||(d==y))    ) 
		  {
		    Ref1[idx] = 2;
		    Ref2[idx] = 3;
		  } 
		else
		  {
		    Ref1[idx] = 3;
		    Ref2[idx] = 3;
		  }

	      } else {
		if (( (c==a)||(c==b)) &&
		    ( (d==a)||(d==b))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 1;
		  } 
		else if (
			 ( (c==a)||(c==b)) ||
			 ( (d==a)||(d==b))    ) 
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
  
  // Write back the recoded alleles
  free(gv->All1); gv->All1 = Ref1;
  free(gv->All2); gv->All2 = Ref2;

  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv;
  char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return;

}

void gs_recode_ref_02_GV ( int* r1,
			   int* r2,
			   int* descending,
			   int* retval,
			   char**  set_name)
{
  gs_recode_ref_02 ( gs_fdta(*set_name),
		     r1,
		     r2,
		     descending,
		     retval            ) ;

}
//##############################################################################
//# End: Recode reference for two parents                                2017-09
//##############################################################################


//##############################################################################
//# BEGIN:  New st_mixed implementation                                     2026
//##############################################################################


/* Symmetrize full matrix from upper triangle (column-major) */
static void sym_full_from_upper_cm(int n, double *A)
{
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      A[i + (size_t)n*(size_t)j] = A[j + (size_t)n*(size_t)i];
    }
  }
}

/* 
  st_mxd: internal cleanup version (keeps  public interface)

  Assumptions :
  - X is always an intercept-only column (vector of 1's), i.e. p = *cX == 1.
  - Each VV[i] is n×n, FULLY STORED, symmetric (both triangles valid).
  - Matrices are passed as "C-style" pointers but because they are fully symmetric,
    the flat memory layout is identical under row-major and column-major flattening.
    => We can safely call Fortran BLAS/LAPACK treating them as column-major.
*/

void st_mxd(
  int*     n,
  double*  y,
  int*     cX,
  double*  X,
  int*     kk,
  double** Z,
  double** VV,
  double*  b,
  double*  u,
  double*  varcmp,
  int*     maxcyc,
  double*  tol,
  int*     retval
)
{
  (void)Z; (void)b; (void)u;

  if (!retval) return;
  *retval = 0;

  if (!n || !y || !cX || !X || !kk || !VV || !varcmp || !maxcyc || !tol) {
    *retval = -1; return;
  }
  if (*n <= 0 || *cX <= 0 || *kk < 0 || *maxcyc <= 0 || *tol <= 0.0) {
    *retval = -2; return;
  }
  if (*cX > *n) {
    st_info(-1, "Invalid dimensions: cX > n");
    *retval = -3; return;
  }

  /* We assume intercept-only */
  if (*cX != 1) {
    st_info(-1, "st_mxd: assumes cX==1 (intercept only).");
    *retval = -4; return;
  }

  const int n_int = *n;
  const int p_int = 1;
  const int k_int = 1 + *kk;               /* kernels + residual */
  const int rankQK = n_int - p_int;        /* n - 1 */
  if (rankQK <= 0) {
    st_info(-1, "rankQK <= 0 (n - 1 <= 0). Cannot proceed.");
    *retval = -5; return;
  }

  const BLAS_INT N  = (BLAS_INT)n_int;
  const BLAS_INT INC1 = 1;

  const size_t nn = (size_t)n_int * (size_t)n_int;
  const BLAS_INT NN = (BLAS_INT)nn;
  /* ----------------------------
     Allocations
     ---------------------------- */
  double *Sigma = NULL;      /* n×n (cm); used as factor / inv */
  double *W     = NULL;      /* inv(Sigma) (cm; full) */
  double *WQK   = NULL;      /* projection (cm; full) */

  double *ones  = NULL;      /* length n */
  double *w1    = NULL;      /* W*1 */
  double *Wy    = NULL;      /* W*y (sometimes useful) */

  double *q     = NULL;      /* q = WQK*y */
  double *tmpn  = NULL;      /* length n scratch */

  double *xvec  = NULL;      /* length k_int */
  double *sigma = NULL;      /* length k_int */
  double *coef  = NULL;      /* length k_int */
  double *taper = NULL;      /* length maxcyc */

  double **T    = NULL;      /* T[i] n×n (cm); for i<k-1: WQK*V[i]; last: WQK */

  double *A     = NULL;      /* k×k (cm) */
  double *Acpy  = NULL;      /* k×k (cm) for solve */
  double *delta = NULL;      /* length k */
  BLAS_INT *ipiv= NULL;      /* pivot (k) */

  int cycle = 0;
  double llik = 0.0, llik_prev = 0.0, delta_llik = 0.0;

  Sigma = (double*)malloc(nn * sizeof(double));
  W     = (double*)malloc(nn * sizeof(double));
  WQK   = (double*)malloc(nn * sizeof(double));
  ones  = (double*)malloc((size_t)n_int * sizeof(double));
  w1    = (double*)malloc((size_t)n_int * sizeof(double));
  Wy    = (double*)malloc((size_t)n_int * sizeof(double));
  q     = (double*)malloc((size_t)n_int * sizeof(double));
  tmpn  = (double*)malloc((size_t)n_int * sizeof(double));

  xvec  = (double*)calloc((size_t)k_int, sizeof(double));
  sigma = (double*)malloc((size_t)k_int * sizeof(double));
  coef  = (double*)malloc((size_t)k_int * sizeof(double));
  taper = (double*)malloc((size_t)(*maxcyc) * sizeof(double));

  T     = (double**)malloc((size_t)k_int * sizeof(double*));

  A     = (double*)malloc((size_t)k_int * (size_t)k_int * sizeof(double));
  Acpy  = (double*)malloc((size_t)k_int * (size_t)k_int * sizeof(double));
  delta = (double*)malloc((size_t)k_int * sizeof(double));
  ipiv  = (BLAS_INT*)malloc((size_t)k_int * sizeof(BLAS_INT));

  if (!Sigma || !W || !WQK || !ones || !w1 || !Wy || !q || !tmpn ||
      !xvec || !sigma || !coef || !taper || !T || !A || !Acpy || !delta || !ipiv) {
    st_info(-1, "Out of memory in st_mxd");
    *retval = -1010;
    goto cleanup;
  }

  for (int i = 0; i < n_int; ++i) ones[i] = 1.0;

  for (int i = 0; i < k_int; ++i) {
    T[i] = (double*)malloc(nn * sizeof(double));
    if (!T[i]) {
      st_info(-1, "Out of memory allocating T[i]");
      *retval = -1010;
      goto cleanup;
    }
  }

  /* Starting values */
  if (varcmp[0] == 0.0) {
    double mean_y = 0.0, ssq = 0.0;
    for (int i = 0; i < n_int; ++i) { mean_y += y[i]; ssq += y[i]*y[i]; }
    mean_y /= (double)n_int;
    double var_y = 1.0;
    if (n_int > 1) {
      var_y = (ssq - (double)n_int * mean_y * mean_y) / (double)(n_int - 1);
      if (!(var_y > 0.0) || !(var_y == var_y)) var_y = 1.0;
    }
    for (int i = 0; i < k_int; ++i) sigma[i] = coef[i] = var_y / (double)k_int;
  } else {
    for (int i = 0; i < k_int; ++i) sigma[i] = coef[i] = varcmp[i];
  }

  /* Taper */
  if (*maxcyc >= 1) taper[0] = 0.5;
  if (*maxcyc >= 2) taper[1] = 0.7;
  for (int i = 2; i < *maxcyc; ++i) taper[i] = 0.9;


  
  /* ----------------------------
     Iterations
     ---------------------------- */
  for (cycle = 0; cycle < *maxcyc; ++cycle) {

    /* Sigma = sum_{i<k-1} sigma[i]*V[i] + sigma_e*I
       VV[i] is n×n symmetric full; treat as column-major.
    */
    /* Sigma <- sigma[0] * V0 */
    memcpy(Sigma, VV[0], nn * sizeof(double));
    {
      const double a = sigma[0];
      F77_CALL(dscal)(&NN, &a, Sigma, &INC1);
    }

    /* add other kernels */
    for (int i = 1; i < k_int - 1; ++i) {
      const double a = sigma[i];
      /* Sigma += a * VV[i] */
      /* tmp: we can use daxpy directly since arrays are dense */
      F77_CALL(daxpy)(&NN, &a, VV[i], &INC1, Sigma, &INC1);
    }

    /* add residual on diagonal */
    const double se = sigma[k_int - 1];
    for (int i = 0; i < n_int; ++i) {
      Sigma[i + (size_t)n_int*(size_t)i] += se;
    }

    /* Cholesky: Sigma = U'U in-place (upper) */
    {
      char uplo = 'U';
      BLAS_INT info = 0;
      F77_CALL(dpotrf)(&uplo, &N, Sigma, &N, &info FCONE);
      if (info != 0) {
        st_info(-1, "dpotrf failed for Sigma (not SPD)");
        *retval = -6;
        goto cleanup;
      }
    }

    /* logdet(Sigma) from chol diag */
    double logdetSigma = 0.0;
    for (int i = 0; i < n_int; ++i) {
      double di = Sigma[i + (size_t)n_int*(size_t)i];
      if (!(di > 0.0) || !(di == di)) { *retval = -7; goto cleanup; }
      logdetSigma += 2.0 * log(di);
    }

    /* Invert Sigma via dpotri: Sigma currently has U, becomes inv(Sigma) in upper */
    {
      char uplo = 'U';
      BLAS_INT info = 0;
      F77_CALL(dpotri)(&uplo, &N, Sigma, &N, &info FCONE);
      if (info != 0) {
        st_info(-1, "dpotri failed for Sigma inverse");
        *retval = -8;
        goto cleanup;
      }
    }

    /* Copy inv(Sigma) to W and symmetrize full */
    memcpy(W, Sigma, nn * sizeof(double));
    sym_full_from_upper_cm(n_int, W);

    /* w1 = W*1 ; Wy = W*y (both symmetric) */
    {
      char uplo = 'U';
      const double one = 1.0, zero = 0.0;
      /* w1 */
      F77_CALL(dsymv)(&uplo, &N, &one, W, &N, ones, &INC1, &zero, w1, &INC1 FCONE);
      /* Wy */
      F77_CALL(dsymv)(&uplo, &N, &one, W, &N, y,    &INC1, &zero, Wy, &INC1 FCONE);
    }

    /* den = 1' w1 */
    double den = 0.0;
    for (int i = 0; i < n_int; ++i) den += w1[i];
    if (!(den > 0.0) || !(den == den)) {
      st_info(-1, "Invalid den = 1'W1 (projection denominator)");
      *retval = -9;
      goto cleanup;
    }

    /* WQK = W - (w1 w1')/den  (rank-1 update)
       Use dsyr on upper triangle: WQK_upper = W_upper + alpha*w1*w1'
       with alpha = -1/den
    */
    memcpy(WQK, W, nn * sizeof(double));
    {
      char uplo = 'U';
      const double alpha = -1.0 / den;
      F77_CALL(dsyr)(&uplo, &N, &alpha, w1, &INC1, WQK, &N FCONE);
    }
    sym_full_from_upper_cm(n_int, WQK);

    /* rss = y' WQK y */
    {
      char uplo = 'U';
      const double one = 1.0, zero = 0.0;
      F77_CALL(dsymv)(&uplo, &N, &one, WQK, &N, y, &INC1, &zero, tmpn, &INC1 FCONE);
    }
    double rss = F77_CALL(ddot)(&N, y, &INC1, tmpn, &INC1);
    if (!(rss > 0.0) || !(rss == rss)) {
      st_info(-1, "Non-positive or NaN rss");
      *retval = -10;
      goto cleanup;
    }

    /* REML-ish objective with p=1:
       llik = -0.5 * ( logdet(Sigma) + log(den) + rankQK*log(rss/rankQK) )
    */
    llik = -0.5 * ( logdetSigma + log(den) + (double)rankQK * log(rss / (double)rankQK) );

    if (cycle == 0) llik_prev = llik;
    delta_llik = llik - llik_prev;
    llik_prev = llik;

    /* q = WQK * y (reuse tmpn if you want; keep q separate for clarity) */
    {
      char uplo = 'U';
      const double one = 1.0, zero = 0.0;
      F77_CALL(dsymv)(&uplo, &N, &one, WQK, &N, y, &INC1, &zero, q, &INC1 FCONE);
    }

    /* Build T:
       - For i<k-1: T[i] = WQK * V[i]  (dsymm, left)
       - For residual: T[k-1] = WQK
    */
    {
      char side = 'L', uplo = 'U';
      const double one = 1.0, zero = 0.0;
      for (int i = 0; i < k_int - 1; ++i) {
        /* T[i] = WQK * VV[i] */
        F77_CALL(dsymm)(&side, &uplo, &N, &N, &one, WQK, &N, VV[i], &N, &zero, T[i], &N FCONE FCONE);
      }
    }
    memcpy(T[k_int - 1], WQK, nn * sizeof(double));

    /* x[i] = y' T[i] q - tr(T[i])  (since q = WQK*y)
       Compute t = T[i] * q via dgemv (general) for i<k-1; dsymv for residual.
    */
    for (int i = 0; i < k_int; ++i) {
      double tr = 0.0;
      /* trace */
      for (int a = 0; a < n_int; ++a) tr += T[i][a + (size_t)n_int*(size_t)a];

      /* tmpn = T[i] * q */
      if (i == k_int - 1) {
        char uplo = 'U';
        const double one = 1.0, zero = 0.0;
        F77_CALL(dsymv)(&uplo, &N, &one, T[i], &N, q, &INC1, &zero, tmpn, &INC1 FCONE);
      } else {
        char trans = 'N';
        const double one = 1.0, zero = 0.0;
        F77_CALL(dgemv)(&trans, &N, &N, &one, T[i], &N, q, &INC1, &zero, tmpn, &INC1 FCONE);
      }

      double yt = F77_CALL(ddot)(&N, y, &INC1, tmpn, &INC1);
      xvec[i] = yt - tr;
    }

    /* A[i,j] = tr(T[i] T[j]) = sum_{a,b} T_i[a,b] * T_j[b,a]
       Efficient O(n^2) using BLAS ddot with stride:
         tr(Ti*Tj) = sum_{col=b} dot( Ti[:,b],  Tj[b,:] )
       In column-major:
         Ti[:,b] is contiguous starting at Ti + n*b (inc=1)
         Tj[b,:] is strided starting at Tj + b (inc=n)
    */
    for (int j = 0; j < k_int; ++j) {
      for (int i = 0; i < k_int; ++i) {
        double s = 0.0;
        const double *Ti = T[i];
        const double *Tj = T[j];
        for (int bcol = 0; bcol < n_int; ++bcol) {
          const BLAS_INT incy = (BLAS_INT)n_int;
          const double *xcol = Ti + (size_t)n_int*(size_t)bcol;
          const double *yrow = Tj + (size_t)bcol; /* row index bcol, col=0 */
          s += F77_CALL(ddot)(&N, xcol, &INC1, yrow, &incy);
        }
        A[i + (size_t)k_int*(size_t)j] = s; /* column-major k×k */
      }
    }

    /* Solve A * delta = xvec
       Try dgesv (pivoted). If singular, optionally fall back to SVD ginv.
    */
    memcpy(Acpy, A, (size_t)k_int*(size_t)k_int*sizeof(double));
    memcpy(delta, xvec, (size_t)k_int*sizeof(double));

{
  /* Prefer La_INT for LAPACK (Windows/R uses La_INT in Lapack.h) */
  La_INT K    = (La_INT)k_int;
  La_INT NRHS = 1;
  La_INT info = 0;

  /* dgesv overwrites Acpy and delta; this assumes Acpy=A and delta=xvec were set before */
  F77_CALL(dgesv)(&K, &NRHS, Acpy, &K, ipiv, delta, &K, &info);

  if (info != 0) {
    /* Fallback: solve min ||A*delta - x|| via SVD (DGELSS), no wrappers */

    /* Restore inputs because dgesv overwrote them */
    memcpy(Acpy, A, (size_t)k_int * (size_t)k_int * sizeof(double));
    memcpy(delta, xvec, (size_t)k_int * sizeof(double));

    double *S = (double*)malloc((size_t)k_int * sizeof(double));
    if (!S) { *retval = -1010; goto cleanup; }

    La_INT rank = 0;
    double rcond = -1.0;          /* use machine precision */
    La_INT LDA = K, LDB = K;

    /* Workspace query */
    double work_q = 0.0;
    La_INT lwork = -1;
    La_INT info2 = 0;

    F77_CALL(dgelss)(&K, &K, &NRHS,
                     Acpy, &LDA,
                     delta, &LDB,
                     S, &rcond, &rank,
                     &work_q, &lwork,
                     &info2);

    if (info2 != 0 || !(work_q == work_q) || work_q <= 0.0) {
      free(S);
      *retval = -11;
      goto cleanup;
    }

    lwork = (La_INT)work_q;
    if (lwork < 1) lwork = 1;

    double *work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) { free(S); *retval = -1010; goto cleanup; }

    /* Actual solve */
    F77_CALL(dgelss)(&K, &K, &NRHS,
                     Acpy, &LDA,
                     delta, &LDB,
                     S, &rcond, &rank,
                     work, &lwork,
                     &info2);

    free(work);
    free(S);

    if (info2 != 0) {
      *retval = -12;
      goto cleanup;
    }
    /* solution is now in delta[0..k_int-1] */
  }
}
    

    /* Update coef/sigma: coef += taper * delta */
    {
      const double a = taper[cycle];
      for (int i = 0; i < k_int; ++i) coef[i] += a * delta[i];

      /* optional: keep variances nonnegative (recommended) */
      const double eps = 1e-12;
      for (int i = 0; i < k_int; ++i) {
        if (!(coef[i] == coef[i])) coef[i] = eps;
        if (coef[i] < eps) coef[i] = eps;
        sigma[i] = coef[i];
      }
    }

    /* Convergence checks */
    if (cycle > 1 && fabs(delta_llik) < (*tol) * 10.0) break;

    double maxabsx = 0.0;
    for (int i = 0; i < k_int; ++i) {
      double ax = fabs(xvec[i]);
      if (ax > maxabsx) maxabsx = ax;
    }
    if (maxabsx < *tol) break;

    /* Optional progress message */
    /* char msg[128]; snprintf(msg, sizeof(msg), "Iter %d: llik=%g", cycle+1, llik); st_info(2, msg); */
  }

  if (cycle == *maxcyc) {
    st_info(-1, "No convergence in variance component estimation");
    *retval = 1; /* warning */
  }

  for (int i = 0; i < k_int; ++i) varcmp[i] = sigma[i];

cleanup:
  if (T) {
    for (int i = 0; i < k_int; ++i) free(T[i]);
    free(T);
  }

  free(Sigma);
  free(W);
  free(WQK);

  free(ones);
  free(w1);
  free(Wy);

  free(q);
  free(tmpn);

  free(xvec);
  free(sigma);
  free(coef);
  free(taper);

  free(A);
  free(Acpy);
  free(delta);
  free(ipiv);

  return;
}



//##############################################################################
//# END:  New st_mixed implementation                                       2026
//##############################################################################


SEXP st_get_Z( SEXP dataSet)
{
  char* set_name = (char*) CHAR(STRING_ELT(dataSet,0));
  gs_varset_type  *gv = gs_fdta(set_name);

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available");
    return  (R_NilValue);
  }

  if (NULL == gv->Z) 
    {
      int Srv; int* Sretval = &Srv;
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) return (R_NilValue); 
    }
  
  SEXP ZZ;
  R_len_t n, m;

  n=(R_len_t) gv->NoInd;
  m=(R_len_t) gv->NoMar;
   
  PROTECT ( ZZ = allocMatrix(REALSXP, n, m)) ;
	  
  double* Z = REAL(ZZ);

  size_t ii = (size_t)gv->NoInd;
  size_t mm = (size_t)gv->NoMar;
  
  for( size_t i = 0; i < ii; i++ ) {
    for( size_t m = 0; m < mm; m++ ) {
      Z[i + m*ii ] = gv->Z[ i*mm + m ];
    }
  }
  
  UNPROTECT(1);
  return (ZZ);
}

SEXP st_get_V( SEXP dataSet)
{
  char* set_name = (char*) CHAR(STRING_ELT(dataSet,0));
  gs_varset_type  *gv = gs_fdta(set_name);

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available");
    return  (R_NilValue);
  }

  if (NULL == gv->Z) 
    {
      int Srv; int* Sretval = &Srv;
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) return (R_NilValue); 
    }

  if (NULL == gv->V) 
    {
      int Srv; int* Sretval = &Srv;
      gs_build_V_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) return (R_NilValue); 
    }
  
  SEXP VV;
  //  R_len_t n, m;
  R_len_t n;

  n=(R_len_t) gv->NoInd;
   
  PROTECT ( VV = allocMatrix(REALSXP, n, n)) ;
	  
  double* V = REAL(VV);

  size_t ii = (size_t)gv->NoInd;

  memcpy (V, gv->V, ii*ii*sizeof(double));
    
  UNPROTECT(1);
  return (VV);
}

SEXP st_get_y( SEXP dataSet)
{
  char* set_name = (char*) CHAR(STRING_ELT(dataSet,0));
  gs_varset_type  *gv = gs_fdta(set_name);

  if ( NULL == gv->y ) {
    gs_info(gv,-2,"No performance data available");
    return  (R_NilValue);
  }

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available");
    return  (R_NilValue);
  }
  
  SEXP yy;
  R_len_t n =(R_len_t) gv->NoInd;
   
  PROTECT ( yy = allocVector(REALSXP, n)) ;
	  
  double* y = REAL(yy);

  size_t ii = (size_t)gv->NoInd;

  memcpy (y, gv->y, ii*sizeof(double));
    
  UNPROTECT(1);
  return (yy);
}

SEXP st_get_im( SEXP i, SEXP m, SEXP dataSet)
{
  char* set_name = (char*) CHAR(STRING_ELT(dataSet,0));
  gs_varset_type  *gv = gs_fdta(set_name);

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available");
    return  (R_NilValue);
  }
  
  int* ii = INTEGER(i);
  int* mm = INTEGER(m);

  *ii = (int)gv->NoInd;
  *mm = (int)gv->NoMar;
  
  return (R_NilValue);

}




SEXP st_mxd_interface( SEXP yy,      
		       SEXP XX,      
		       SEXP ZZ,      
		       SEXP VV,
		       SEXP mmaxcyc,
		       SEXP pprecision,
		       SEXP rrlst,
		       SEXP rretval)
{

  // n
  R_xlen_t n_list = XLENGTH(yy);
  int n[1] = {n_list}; 

  // y
  double* y = REAL(yy);

  // cX
  int cX[1] = { (int)XLENGTH(XX)/(int)XLENGTH(yy)};

  // X
  double * X =  REAL(XX);

  double **Z = NULL;
  int k[1];           // No of var.comp.s without var(e)

  // Z
  if ( TYPEOF(ZZ) != NILSXP ) {
    R_xlen_t n_list = XLENGTH(ZZ);
    *k = (int)n_list;
    
    Z = malloc( (size_t)n_list * sizeof(double*));
  
    for (R_xlen_t i = 0; i < n_list; i++)
      Z[(int)i] = REAL(VECTOR_ELT(ZZ, i));
  }

  double **V = NULL;
  
  // V
  if ( TYPEOF(VV) != NILSXP ) {
    R_xlen_t n_list = XLENGTH(VV);
    *k = (int)n_list;
    
    V = malloc( (size_t)n_list * sizeof(double*));
  
    for (R_xlen_t i = 0; i < n_list; i++)
      V[(int)i] = REAL(VECTOR_ELT(VV, i));
  }

  int *maxcyc = INTEGER(mmaxcyc);
  double * precision = REAL(pprecision);
  
  // Return values b, u, varcomp

  double *varcmp  = REAL(VECTOR_ELT(rrlst, 0));
  double *b       = REAL(VECTOR_ELT(rrlst, 1));
  double *u       = REAL(VECTOR_ELT(rrlst, 2));

  int *retval;  retval = INTEGER(rretval);
 
   
  st_mxd (n,     // Number of observations
  	  y,     // Observations, dim: n x 1
  	  cX,    // Numer of columns of X, full rank
  	  X,     // Design matrix, dim: n x cX
  	  k,     // Number of variance components
  	  Z,     // Array of the pointers to Z matrices | either Z or
  	  V,      // Array of the pointers to V matrices | V is NULL
  	  b,     // BLUEs of fixed effects, NULL on input: don't calculate
  	  u,      // BLUPs of random effects, NULL on input: don't calculate
  	  varcmp, // Returned variance components
  	  maxcyc,
  	  precision,
  	  retval
  	  );
 	       
  free(V);
  free(Z);
  
  return (R_NilValue);
}


// ##############################################################################
// BEGIN: Testing st_mxd                                                     2026
// ##############################################################################


static void assert_matrix_double_nxn(SEXP M, int n, const char* what) {
  if (!isReal(M) || !isMatrix(M)) {
    error("%s must be a numeric matrix", what);
  }
  SEXP dim = getAttrib(M, R_DimSymbol);
  if (LENGTH(dim) != 2) error("%s must have 2 dimensions", what);
  int nr = INTEGER(dim)[0], nc = INTEGER(dim)[1];
  if (nr != n || nc != n) {
    error("%s must be %dx%d; got %dx%d", what, n, n, nr, nc);
  }
}

SEXP st_mxd_R(SEXP y_, SEXP X_, SEXP VV_, SEXP varcmp_,
              SEXP maxcyc_, SEXP tol_)
{
  if (!isReal(y_) || !isVector(y_)) error("y must be a numeric vector");
  if (!isReal(X_) || !isMatrix(X_)) error("X must be a numeric matrix");
  if (!isNull(varcmp_) && !isReal(varcmp_)) error("varcmp must be NULL or numeric");
  if (!isInteger(maxcyc_) && !isReal(maxcyc_)) error("maxcyc must be numeric/integer");
  if (!isReal(tol_) && !isInteger(tol_)) error("tol must be numeric/integer");

  int n = (int)XLENGTH(y_);

  // X dims
  SEXP Xdim = getAttrib(X_, R_DimSymbol);
  int nX = INTEGER(Xdim)[0];
  int cX = INTEGER(Xdim)[1];
  if (nX != n) error("nrow(X) must equal length(y)");

  // VV must be a list of matrices
  if (!isNewList(VV_)) error("VV must be a list of n×n numeric matrices");
  int kk = (int)XLENGTH(VV_);
  if (kk < 0) error("VV list length must be >= 0");

  // Build double** pointer array for VV
  // Use R_alloc so it is auto-freed
  double **VV_ptrs = (double**)R_alloc((size_t)kk, sizeof(double*));

  for (int i = 0; i < kk; ++i) {
    SEXP Vi = VECTOR_ELT(VV_, i);
    char buf[64];
    snprintf(buf, sizeof(buf), "VV[[%d]]", i + 1);
    assert_matrix_double_nxn(Vi, n, buf);
    VV_ptrs[i] = REAL(Vi);
  }

  // varcmp: length must be kk+1 (includes residual), or NULL to auto-init in C
  SEXP varcmp_out;
  double *varcmp_ptr;

  if (isNull(varcmp_)) {
    varcmp_out = PROTECT(allocVector(REALSXP, (R_xlen_t)(kk + 1)));
    varcmp_ptr = REAL(varcmp_out);
    // signal “no starting values” the way your C expects: varcmp[0]==0
    for (int i = 0; i < kk + 1; ++i) varcmp_ptr[i] = 0.0;
  } else {
    if (XLENGTH(varcmp_) != (R_xlen_t)(kk + 1))
      error("varcmp must have length kk+1 (kk=%d => %d)", kk, kk + 1);
    varcmp_out = PROTECT(duplicate(varcmp_));
    varcmp_ptr = REAL(varcmp_out);
  }

  int maxcyc = asInteger(maxcyc_);
  double tol = asReal(tol_);
  int retval = 0;

  // Call into C (Z,b,u are not used in your current code path; pass NULL)
  st_mxd(&n, REAL(y_), &cX, REAL(X_), &kk,
         NULL, VV_ptrs,
         NULL, NULL, varcmp_ptr,
         &maxcyc, &tol, &retval);

  // Return list(varcmp=..., retval=...)
  SEXP ans = PROTECT(allocVector(VECSXP, 2));
  SEXP nms = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, mkChar("varcmp"));
  SET_STRING_ELT(nms, 1, mkChar("retval"));
  setAttrib(ans, R_NamesSymbol, nms);

  SET_VECTOR_ELT(ans, 0, varcmp_out);
  SEXP ret = PROTECT(ScalarInteger(retval));
  SET_VECTOR_ELT(ans, 1, ret);

  UNPROTECT(4);
  return ans;
}



// ##############################################################################
// END: Testing st_mxd                                                       2026 
// ##############################################################################



void gs_lambda_rrblup_01 ( gs_varset_type* gv ,
		 	  int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{

  (void) maxiter;
  (void) precision;
  (void) hsq;
    
  // Initialize Vector of shrinkage factors lambda
  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if (NULL==gv->ladi) ERR_M;

  double lambda =  gv->sigma[0]/gv->sigma[1];

  for (gs_UNSIGNED i = 0 ; i < gv->NoMar ; i++)  gv->ladi[i] = lambda;

  
  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for (gs_UNSIGNED m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rrblup_01_GV ( int*    maxiter        ,
			     double* precision      ,
			     double* hsq            ,
			     char**  out_filename   ,
			     int*    auxfiles       , 
			     int*    retval         ,
			     char**  set_name         )
{
  gs_lambda_rrblup_01 ( gs_fdta(*set_name),
		       maxiter,
		       precision,
		       hsq,
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}

void gs_vc_rrblup_01 ( gs_varset_type* gv ,
		       int*    maxiter,
		       double* precision,
		       double* hsq,
		       char**  out_filename,
		       int*    auxfiles, 
		       int*    retval         )
{

  char msg[gs_MAXSTR];

  // gs_UNSIGNED m;
  // int i;
  // double qu, maxqu=0;
  //  int Srv; int* Sretval = &Srv;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available");
    *retval = -2; return ;
  }

  if ( NULL == gv->y ) {
    gs_info(gv,-2,"No performance data available");
    *retval = -2; return ;
  }

  
  int  k  = 1 ;
  int  n[1]  = { (int)gv->NoInd };
  int  cX[1] = {1};

  // Prepare y, X and V matrices

  double* y = malloc ( gv->NoInd * sizeof(double) );
  double* X = malloc ( gv->NoInd * sizeof(double) );
  double** V    = malloc ( k * sizeof(double*) );
  V[0] = malloc ( gv->NoInd * gv->NoInd * sizeof(double) );

if ( NULL==y )    {*retval = -2;return;}
if ( NULL==X )    {*retval = -2;return;}
if ( NULL==V )    {*retval = -2;return;}
if ( NULL==V[0] ) {*retval = -2;return;}

  st_dgemm ("N","T",
            gv->NoInd, gv->NoMar,
  	    gv->NoInd, gv->NoMar,
  	    gv->Z, gv->Z, V[0]);
  
  for (int i = 0; i < *n; i++)
      X[i] = 1.0;

  st_dcopy( *n, y, gv->y );

  // Starting values for variance components from hsq
  double varcmp[2] = {0,0};

  double mean_y = 0;  double var_y = 0;
  for (int i=0; i<*n; i++)
    mean_y += y[i] , var_y += (y[i])*(y[i]);
  mean_y = mean_y / *n ;
  var_y = ( var_y - *n * mean_y * mean_y ) / ( *n - 1 );
  
  double vg = *hsq * var_y;
  double ve = var_y - vg;

  varcmp[0] = vg/(double)gv->NoMar; varcmp[1] = ve;

  
  if (0 < *maxiter) // otherwise this is rrwe with hsq
    {
      snprintf(msg, sizeof(msg),
         "Start values from hsq = %4.2f: var(m) = %f, var(e) = %5f",
         *hsq, varcmp[0], varcmp[1]);
      gs_info(gv,2,msg);
      
      st_mxd (n,     // Number of observations
	      y,     // Observations, dim: n x 1
	      cX,    // Numer of columns of X, full rank 
	      X,     // Design matrix, dim: n x cX
	      &k,     // Number of variance components
	      NULL,   // Array of the pointers to Z matrices | either Z or
	      V,      // Array of the pointers to V matrices | V is NULL
	      NULL,   // BLUEs of fixed effects, NULL on input: don't calculate
	      NULL,   // BLUPs of random effects, NULL on input: don't calculate
	      varcmp, // Returned variance components		   
	      maxiter,
	      precision,
	      retval   
	      );
    }

  snprintf(msg, sizeof(msg),
         "Variance components: var(m) = %f, var(e) = %5f",
         varcmp[0], varcmp[1]);
 
  gs_info(gv,2,msg);

  /* Allocate memory for the variace components */
  gv->sigma =  realloc( gv->sigma, (1+gv->NoMar) * sizeof(double) );
  if ( NULL == gv->sigma ) ERR_M;
  
  gv->sigma[0] = varcmp[1];

  for (gs_UNSIGNED i=1; i < (1+gv->NoMar); i++)
    gv->sigma[i] = varcmp[0];

  free(y); free(X);
  free(V[0]);  free(V);


  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"vc \n");
      for (gs_UNSIGNED m = 0; m < 2 ; m++)
      	fprintf(outfp,"%f\n", gv->sigma[m]);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}


void gs_vc_rrblup_01_GV ( int*    maxiter        ,
			  double* precision      ,
			  double* hsq            ,
			  char**  out_filename   ,
			  int*    auxfiles       , 
			  int*    retval         ,
			  char**  set_name         )
{
  gs_vc_rrblup_01 ( gs_fdta(*set_name),
		    maxiter,
		    precision,
		    hsq,
		    out_filename,
		    auxfiles, 
		    retval         ) ;
}



void gs_single_marker_aov_01 ( gs_varset_type *gv  ,
			       double* alpha       ,
			       char** out_filename ,
			       int*   auxfiles     , 
			       int*   retval         )
{
  /* single marker aov with a random model to estimate the
     variance component of one marker.
     0 1 2 design matrix
  */

  (void) alpha;
  
  gs_UNSIGNED i,m,idx_Z;
  FILE *outfp;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2, "No Z matrix available, run 'build.Z' first");
    *retval = -2; return;
  }

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return;
  }

  /* Check the design matrix Z for 0,1,2 */
  int ext = 0; double zz;
  for ( i = 0; i < gv->NoInd ; i++) 
    for ( m = 0; m < gv->NoMar; m++) {
      idx_Z   = i*gv->NoMar + m ;
      zz = gv->Z[idx_Z] ;
      if ( (0!=zz)&&(1!=zz)&&(2!=zz) ) {ext = 1;break;} ;
    }
  
  if ( 1 == ext ) {
    gs_info(gv,-2,
      "Incorrect design matrix for single marker anova");
    *retval = -2; return;
  }

  /* Allocate memory for the results */
  gv->SMp   =  realloc( gv->SMp, gv->NoMar   * sizeof(double) );
  gv->SMF   =  realloc( gv->SMF, gv->NoMar   * sizeof(double) );
  gv->SMsqm =  realloc( gv->SMsqm, gv->NoMar * sizeof(double) );
  gv->SMsqe =  realloc( gv->SMsqe, gv->NoMar * sizeof(double) );
  gv->SMvar =  realloc( gv->SMvar, gv->NoMar * sizeof(double) );
  if ( ( NULL ==  gv->SMp)   || ( NULL ==  gv->SMF ) ||
       ( NULL ==  gv->SMsqm) || ( NULL ==  gv->SMsqe ) ||
       ( NULL ==  gv->SMvar )                             ) ERR_M;

  /* Initialize vector for p-values */
  for ( m = 0; m < gv->NoMar ; m++) gv->SMp [m] = 0;

  double k = (double) (gv->NoInd);
  double sy=0, syy=0,
    SQT=0, SQM=0, SQE=0, mm, kk;
  
  double sy0, sy1, sy2, 
         ny0, ny1, ny2;

  /* Calculate single marker anova */

  for ( i = 0; i < gv->NoInd; i++ )
    {
      sy  += gv->y[i] ;
      syy += gv->y[i] * gv->y[i] ;
    }
  SQT = syy - sy*sy/k;

  for ( m = 0; m < gv->NoMar ; m++) 
    {
      sy0=0; sy1=0; sy2=0;
      ny0=0; ny1=0; ny2=0;

      for ( i = 0; i < gv->NoInd; i++ )
	{
	  idx_Z   = i*gv->NoMar + m ;
	  if (0==gv->Z[idx_Z]) { ny0++; sy0 += gv->y[i]; } else
	  if (2==gv->Z[idx_Z]) { ny2++; sy2 += gv->y[i]; } else
	  if (1==gv->Z[idx_Z]) { ny1++; sy1 += gv->y[i]; } 
	}

      mm  =    ( 0==ny0 ? 0 : 1)
             + ( 0==ny1 ? 0 : 1)
	     + ( 0==ny2 ? 0 : 1 );

      SQE = syy - ( 0==ny0 ? 0 : sy0*sy0/ny0)
	        - ( 0==ny1 ? 0 : sy1*sy1/ny1)
    	        - ( 0==ny2 ? 0 : sy2*sy2/ny2) ;

      SQM = SQT - SQE ; 

      gv->SMsqm [m] = SQM ;
      gv->SMsqe [m] = SQE ;

      if ( (2>mm) || ((k-mm)<1) ) {
	  gv->SMvar[m] = FLT_MIN;
	  gv->SMF[m] = 0 ;
	  gv->SMp[m] = 1 ;
      } 
      else {
	  kk = ( 2/(mm-1) ) * ( k - (ny0*ny0 + ny1*ny1 + ny2*ny2)/k );
	  gv->SMvar [m] = (SQM/(mm-1) - SQE/(k-mm)) / kk;
	  //if (gv->SMvar[m] <= 0 ) gv->SMvar[m] = FLT_MIN;
	  gv->SMF[m] = ( SQM/(mm-1) ) / ( SQT/(k-mm));
	  gv->SMp[m] = pf( gv->SMF[m] ,(int)(mm-1),(int)(k-mm),0,0);
	  //if ( gv->SMp[m] > *alpha )  gv->SMvar[m] = FLT_MIN;
      }
    }

  double vcmin;
  vcmin = fabs ( gv->SMvar[0] ) ;
  for ( m = 1; m < gv->NoMar ; m++)    // handle negative variances
    if ( fabs(gv->SMvar[m] ) < vcmin) vcmin =  fabs(gv->SMvar[m] );
  vcmin = vcmin/1000.0;						

  for ( m = 0; m < gv->NoMar ; m++)   
    if ( gv->SMvar[m]  < 0 ) gv->SMvar[m] =vcmin;

  
  /* Write solution vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"SQM SQE var F p\n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f %f %f %f %f\n",
		gv->MarName[m],
		gv->SMsqm [m],
		gv->SMsqe [m],
		gv->SMvar [m],
		gv->SMF[m],
		gv->SMp[m]
		);

      fclose(outfp);  
  }

  *retval = 0; return;

}

void gs_single_marker_aov_01_GV ( double* alpha       ,
				  char** out_filename,
				  int*   auxfiles    , 
				  int*   retval      ,
				  char** set_name     )
{
  gs_single_marker_aov_01 ( gs_fdta(*set_name) ,
			    alpha       ,
			    out_filename,
			    auxfiles, 
			    retval         );
}

void gs_lambda_rmla_02 (  gs_varset_type *gv     ,
			   double* alpha          ,
			   int*    maxiter        ,
			   double* precision      ,
			   double* hsq            ,
			   char**  out_filename   ,
			   int*    auxfiles       , 
			   int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;

  (void) alpha;
  (void) maxiter;
  (void) precision;
  (void) hsq ;
  
  //  int Srv; int* Sretval = &Srv;

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factors available that could be adjusted");
      *retval = -2; return ;
  }

  
  long double Svar = 0.0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    Svar += gv->SMvar[m] ;

  Svar = Svar / (long double) gv->NoMar;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] =  ( Svar / gv->SMvar[m]  ) * gv->ladi[m];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rmla_02_GV ( double* alpha,
			     int*    maxiter       ,
			     double* precision     ,
			     double* hsq           ,
			     char**  out_filename  ,
			     int*    auxfiles      , 
			     int*    retval        ,
			     char** set_name        )
{
  gs_lambda_rmla_02 ( gs_fdta(*set_name),
		       alpha             ,
		       maxiter           ,
		       precision         ,
		       hsq               ,
		       out_filename      ,
		       auxfiles          , 
		       retval             ) ;
}

void gs_w_aov_01 (  gs_varset_type *gv     ,
		    double* alpha          ,
		    char**  out_filename   ,
		    int*    auxfiles       , 
		    int*    retval         )

// Caclulates the diagonal of a weighing matrix W
// and stores it in gv->lambda
  
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_aov_01 ( gv, alpha, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }
  
  long double Svar = 0.0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    Svar += gv->SMvar[m] ;

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if ( NULL== gv->ladi )  ERR_M;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] =  (  gv->SMvar[m] / Svar ) * (double) gv->NoMar;
  
  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"W \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_w_aov_01_GV ( double* alpha,
		      char**  out_filename  ,
		      int*    auxfiles      , 
		      int*    retval        ,
		      char** set_name        )
{
  gs_w_aov_01 ( gs_fdta(*set_name),
		alpha             ,
		out_filename      ,
		auxfiles          , 
		retval             ) ;
}

void gs_w_I_01 (  gs_varset_type *gv     ,
		    char**  out_filename   ,
		    int*    auxfiles       , 
		    int*    retval         )

// Caclulates the diagonal of a weighing matrix W
// and stores it in gv->lambda
  
{
  gs_UNSIGNED m;
  FILE *outfp;

  //  int Srv; int* Sretval = &Srv;

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if ( NULL== gv->ladi )  ERR_M;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] =  1; 

  
  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"W \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_w_I_01_GV ( char**  out_filename  ,
		    int*    auxfiles      , 
		    int*    retval        ,
		    char** set_name        )
{
  gs_w_I_01 ( gs_fdta(*set_name),
	      out_filename      ,
	      auxfiles          , 
	      retval             ) ;
}


//##############################################################################
//# End:  NR Algorithm for RR-BLUP                                       2017-11
//##############################################################################


//##############################################################################
//# Begin: EM based genomic prediction                                   2012-18
//##############################################################################

/*
 + The functions in this block provide
 +
 + gs.esteff.lsq  and gs.esteff.lsq  
 + gs.esteff.rmlc     gs.esteff.rmlc 
 + gs.esteff.rmlv     gs.esteff.rmlv 
 + gs.esteff.rmla  
 +
 + using the new linear algebra routines
 +
 + They are not required to run gs.esteff.rr
 +
 */

void gs_lambda_const_02 ( gs_varset_type *gv ,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  gs_UNSIGNED i,m;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available, run 'build.Z' first");
    *retval = -2; return ;
  }

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if ( NULL== gv->ladi )  ERR_M;

  double lambda =  ( 1 / *hsq -1 ) * (double)gv->NoMar; 

  for ( i = 0 ; i < gv->NoMar ; i++)  gv->ladi[i] = lambda;

  /* Write output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;

      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    } 

  *retval = 0; return ;

}

void gs_lambda_const_02_GV ( double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval         ,
			     char** set_name         )
{
  gs_lambda_const_02 ( gs_fdta(*set_name),
		       hsq,
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}



void gs_lambda_aov_01 (  gs_varset_type *gv ,
			 double* hsq,
                	 double* alpha,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_aov_01 ( gv, alpha, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if (NULL==gv->ladi) ERR_M;

  double Svar = 0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    Svar += gv->SMvar[m] ;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    {
      gv->ladi[m] = ( ( 1 / *hsq) - 1 ) * Svar / gv->SMvar[m];
    }

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_aov_01_GV ( double* hsq           ,
			   double* alpha         ,
			   char**  out_filename  ,
			   int*    auxfiles      , 
			   int*    retval        ,
			   char** set_name        )
{
  gs_lambda_aov_01 ( gs_fdta(*set_name),
		     hsq,
		     alpha,
		     out_filename,
		     auxfiles, 
		     retval         ) ;
}


void gs_single_marker_reg_01 ( gs_varset_type *gv  ,
			       char** out_filename ,
			       int*   auxfiles     , 
			       int*   retval         )
{
  gs_UNSIGNED i,m,idx_Z;
  FILE *outfp;


  if ( NULL == gv->Z ) {
    gs_info(gv,-2, "No Z matrix available, run 'build.Z' first");
    *retval = -2; return;
  }

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return;
  }

  /* Allocate memory for the results */
  gv->SMp   =  realloc( gv->SMp, gv->NoMar   * sizeof(double) );
  gv->SMF   =  realloc( gv->SMF, gv->NoMar   * sizeof(double) );
  gv->SMsqm =  realloc( gv->SMsqm, gv->NoMar * sizeof(double) );
  gv->SMsqe =  realloc( gv->SMsqe, gv->NoMar * sizeof(double) );
  if ( ( NULL ==  gv->SMp)   || ( NULL ==  gv->SMF ) ||
       ( NULL ==  gv->SMsqm) || ( NULL ==  gv->SMsqe ) ) ERR_M;

  /* Initialize vector for p-values */
  for ( m = 0; m < gv->NoMar ; m++) gv->SMp [m] = 0;

  double k = (double) (gv->NoInd);
  double sy=0, syy=0, sx=0, sxx=0, sxy=0,
    b,
    SQT=0, SQM=0, SQE=0, F, p;
  
  /* Calculate single marker regression */

  for ( i = 0; i < gv->NoInd; i++ )
    {
      sy  += gv->y[i] ;
      syy += gv->y[i] * gv->y[i] ;
    }
  SQT = syy - sy*sy/k;

  for ( m = 0; m < gv->NoMar ; m++) 
    {
      sx=0; sxx=0; sxy=0;
      for ( i = 0; i < gv->NoInd; i++ )
	{
	  idx_Z   = i*gv->NoMar + m ;
	  sx  += gv->Z[idx_Z];
	  sxx += gv->Z[idx_Z] * gv->Z[idx_Z] ;
	  sxy += gv->Z[idx_Z] * gv->y[i]     ;
	}

      b = ( sxy - sx*sy/k ) / 
          ( sxx - sx*sx/k) ;
    
      SQM = b * (sxy - sx*sy/k);
      SQE = SQT - SQM;
      F   = SQM / (SQE/(k-2)) ;
      p   = pf(F,1,(int)gv->NoInd,0,0);

      gv->SMsqm [m] = SQM ;
      gv->SMsqe [m] = SQE ;
      gv->SMF   [m] = F   ;
      gv->SMp   [m] = p   ; 

    }

  /* Write solution vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"SQM SQE F p\n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f %f %f %f\n",
		gv->MarName[m],
		gv->SMsqm [m],
		gv->SMsqe [m],
		gv->SMF[m],
		gv->SMp[m]
		);

      fclose(outfp);  
  }

  *retval = 0; return;

}

void gs_single_marker_reg_01_GV ( char** out_filename,
				  int*   auxfiles    , 
				  int*   retval      ,
				  char** set_name     )
{
  gs_single_marker_reg_01 ( gs_fdta(*set_name) ,
			    out_filename,
			    auxfiles, 
			    retval         );
}


void gs_lambda_reg_01 (  gs_varset_type *gv ,
			 double* hsq,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_reg_01 ( gv, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if (NULL==gv->ladi) ERR_M;

  double SSQM = 0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    SSQM += gv->SMsqm[m] ;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] = ( ( 1 / *hsq) - 1 ) * SSQM / gv->SMsqm[m];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_reg_01_GV ( double* hsq           ,
			   char**  out_filename  ,
			   int*    auxfiles      , 
			   int*    retval        ,
			   char** set_name        )
{
  gs_lambda_reg_01 ( gs_fdta(*set_name),
		     hsq,
		     out_filename,
		     auxfiles, 
		     retval         ) ;
}


void gs_mme_coeff_01 ( gs_varset_type* gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  if ( NULL == gv->y ) {
      gs_info(gv,-2,
        "No performance data available, run 'read.performance.data' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
      gs_info(gv,-2,
        "No Z matrix available, run 'build.Z' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factor defined");
      *retval = -2; return ;
  }

  /* Allocate memory for the A matrix */

  gv->dimA = 1 + gv->NoMar;
  gv->A =  realloc( gv->A, gv->dimA * gv->dimA * sizeof(double) );
  if ( NULL == gv->A ) ERR_M;

  for ( i = 0; i < gv->dimA ; i++) 
    for ( j = 0; j < gv->dimA ; j++) 
	gv->A [ i*gv->dimA + j ] =  0;

  /* Initialize the 1 and I matrix */

  double *one_vec; 

  one_vec = malloc ( gv->NoInd * sizeof(double) );
  if ( NULL == one_vec ) ERR_M;
  for ( i = 0; i < gv->NoInd ; i++) one_vec[i] = 1;

  /* Init part 1*/
  
  gv->A [0] = (double) gv->NoInd;

  /* Init part 2 and 3*/

  double *tmp1;
  tmp1 = malloc ( gv->NoMar * sizeof(double) );
  if ( NULL == tmp1 ) ERR_M;

    st_dgemm ("T","N",
	    gv->NoInd, 1,
  	     gv->NoInd, gv->NoMar,
  	     one_vec, gv->Z, tmp1);

  for ( j = 0; j < gv->NoMar; j++)  
    {
      gv->A [     0*gv->dimA + (1+j) ] = tmp1[j];
      gv->A [ (j+1)*gv->dimA + 0     ] = tmp1[j];
    }

  /* Init part 4*/

  /* implementation of tmp2 as nr-matrix would save space */
  double *tmp2;
  tmp2 = malloc ( gv->NoMar * gv->NoMar *sizeof(double) );
  if ( NULL == tmp2 ) ERR_M;

    st_dgemm ("T","N",
	    gv->NoInd, gv->NoMar,
  	     gv->NoInd, gv->NoMar,
  	     gv->Z, gv->Z, tmp2);
  
  /* Build A */
  for ( i = 0; i < gv->NoMar; i++)
    for ( j = 0; j < gv->NoMar; j++)
      {
	gv->A [ (1+i)*gv->dimA + (1+j) ] = tmp2[ i*gv->NoMar + j];
      }

  /* Store A without applied lambda */

  if (1) 
    {
      gv->As =  realloc( gv->As, gv->dimA * gv->dimA *sizeof(double) );
      if ( NULL == gv->As ) ERR_M;
      for ( i = 0; i < gv->dimA; i++)
	for ( j = 0; j < gv->dimA; j++)
	  {
	    gv->As [ i*gv->dimA + j ] = gv->A [ i*gv->dimA + j ] ;
	  }
    }

  /* Apply lambda */
  for ( i = 1; i < gv->dimA; i++)
	gv->A [ i*gv->dimA + i ]  += gv->ladi[i-1] ;

  free(one_vec);
  free(tmp1);
  free(tmp2);

  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_mme_coeff_01_GV ( char** out_filename,
			  int*   auxfiles, 
			  int*   retval,
			  char** set_name                )
{
  gs_mme_coeff_01 (gs_fdta(*set_name),
		   out_filename,
		   auxfiles, 
		   retval        ) ;
}

void gs_mme_restcoeff_01 ( gs_varset_type* gv ,
			   char** out_filename,
			   int*   auxfiles, 
			   int*   retval        )
{

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  int Srv; int* Sretval = &Srv;

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factor defined");
      *retval = -2; return ;
  }

  if ( NULL == gv->As ) 
    {
      /* No stored value availabe, Calculate */
      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }
  else
    {
      /* Restore A */
      gv->dimA = 1 + gv->NoMar;
      gv->A =  realloc( gv->A, gv->dimA * gv->dimA * sizeof(double) );
      if ( NULL == gv->A ) ERR_M;

      memcpy ( gv->A , gv->As, (size_t)gv->dimA * 
                               (size_t)gv->dimA * sizeof(double));
      
      /* Apply lambda */
      for ( i = 1; i < gv->dimA; i++)
	gv->A [ i*gv->dimA + i ]  += gv->ladi[i-1] ;
    }

  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_mme_restcoeff_01_GV ( char** out_filename,
			 int*   auxfiles, 
			 int*   retval        )
{
  gs_mme_restcoeff_01 (gv,
		       out_filename,
		       auxfiles, 
		       retval        ) ;
}




void gs_mme_rhs_01 ( gs_varset_type *gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{
  gs_UNSIGNED i,j, dim_b;
  FILE *outfp;

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available, run 'build.Z' first");
    *retval = -2; return ;
  }

  /* Allocate memory for the b and initialize */
  dim_b = 1 + gv->NoMar;
  gv->b =  realloc( gv->b, dim_b * sizeof(double) );
  if ( NULL == gv->b ) ERR_M;
  for ( i = 0; i < dim_b ; i++) gv->b [i] =  0;

  double *one_vec;
  one_vec = malloc ( gv->NoInd * sizeof(double) );
  if ( NULL == one_vec ) ERR_M;
  for ( i = 0; i < gv->NoInd ; i++) one_vec[i] = 1;

  /* Init part 1*/
  
  double *tmp ;
  tmp = malloc ( 1 * sizeof(double) );
  if ( NULL == tmp ) ERR_M;

    st_dgemm ("T","N",
	    gv->NoInd, 1,
  	     gv->NoInd, 1,
  	     one_vec, gv->y, tmp);

  gv->b [0] = tmp[0];

  /* Init part 2*/

  double *tmp2 ;
  tmp2 = malloc ( gv->NoMar *sizeof(double) );
  if ( NULL == tmp2 ) ERR_M;

    st_dgemm ("T","N",
	    gv->NoInd, gv->NoMar,
  	     gv->NoInd, 1,
  	     gv->Z, gv->y, tmp2);

  for ( j = 0; j < gv->NoMar; j++)  
    gv->b [ 1+j ] = tmp2[j];

  /* Store b */
  if (1) 
    {
      gv->bs =  realloc( gv->bs, dim_b * sizeof(double) );
      if ( NULL == gv->bs ) ERR_M;
      for ( j = 0; j < dim_b; j++)  
	gv->bs[j] = gv->b[j];
    }

  free(one_vec);
  free(tmp);
  free(tmp2);

  /* Write b vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( i = 0; i < dim_b ; i++)
	fprintf(outfp,"%f \n",gv->b [i] );
      
      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mme_rhs_01_GV ( char** out_filename,
			int*   auxfiles, 
			int*   retval,
			char** set_name    )
{
  gs_mme_rhs_01 ( gs_fdta(*set_name) ,
		  out_filename,
		  auxfiles, 
		  retval        ) ;
}

void gs_mme_restrhs_01 ( gs_varset_type *gv ,
			 char** out_filename,
			 int*   auxfiles, 
			 int*   retval        )
{
  gs_UNSIGNED i, dim_b ;
  FILE *outfp;

  int Srv; int* Sretval = &Srv;

  dim_b = 1 + gv->NoMar;

  if ( NULL == gv->bs ) 
    {
      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }
  else
    {
      gv->b =  realloc( gv->b, dim_b * sizeof(double) );
      if ( NULL == gv->b ) ERR_M;

      memcpy ( gv->b, gv->bs, (size_t) dim_b );
    }

  /* Write b vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( i = 0; i < dim_b ; i++)
	fprintf(outfp,"%f \n",gv->b [i] );
      
      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mme_restrhs_01_GV ( char** out_filename,
			    int*   auxfiles, 
			    int*   retval        )
{
  gs_mme_restrhs_01 ( gv ,
		     out_filename,
		     auxfiles, 
		     retval        ) ;
}


void gs_mme_solve_01 ( gs_varset_type *gv ,
		      char** out_filename,
		      int*   auxfiles, 
		      int*   retval        )
{
  gs_UNSIGNED i;
  FILE *outfp;

  int Srv=0; int* Sretval = &Srv;

  if ( NULL == gv->A ) {
      gs_info(gv,-2, "No coefficient matrix available");
      *retval = -2; return ;
  }

  if ( NULL == gv->b )  {
      gs_info(gv,-2, "No right hand side available");
      *retval = -2; return ;
  }

  /* Allocate memory for the solution vector */
  gv->dimu = 1 + gv->NoMar;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );
  if ( NULL == gv->u ) ERR_M;

  st_dgesv ( gv->dimA, gv->dimA, gv->A, gv->u, gv->b, Sretval );

  free(gv->A); 
  gv->A =NULL;

  if (*Sretval==-2) {*retval = -2; return ;}

  /* Write solution vector to output file */
  if (1 == *auxfiles) {
   
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"effect\n");
      for ( i = 0; i < gv->dimu ; i++)
	fprintf(outfp,"%s.%i %f \n",
		gv->EffNme[i],
		gv->EffAll[i],
		gv->u [i]);

      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mme_solve_01_GV ( char** out_filename ,
			  int*   auxfiles     , 
			  int*   retval       ,
			  char** set_name      )
{
  gs_mme_solve_01 ( gs_fdta(*set_name) ,
		    out_filename,
		    auxfiles, 
	            retval        ) ;
}

void gs_esteff_rrwr_01 (  gs_varset_type *gv ,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_reg_01 ( gv, hsq, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rrwr_01_GV ( double* hsq,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval,
			    char** set_name         )
{
  gs_esteff_rrwr_01 ( gs_fdta(*set_name),
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_esteff_rrwa_01 (  gs_varset_type *gv ,
			  double* hsq,
			  double* alpha,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_aov_01 ( gv, hsq, alpha, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rrwa_01_GV ( double* hsq,
                            double* alpha,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval,
			    char** set_name         )
{
  gs_esteff_rrwa_01 ( gs_fdta(*set_name),
		      hsq,
		      alpha,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}


void gs_esteff_rrwe_01 (  gs_varset_type *gv ,
			  double* lambda, 
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_const_01 ( gv, lambda,nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rrwe_01_GV ( double* lambda   , 
			    char**  out_filename,
			    int*    auxfiles , 
			    int*    retval   ,
			    char** set_name   )
{
  gs_esteff_rrwe_01 ( gs_fdta(*set_name),
		      lambda            ,
		      out_filename      ,
		      auxfiles          , 
		      retval             ) ;
}

void gs_esteff_rrwe_02 (  gs_varset_type *gv ,
		       double* hsq,
		       char**  out_filename,
		       int*    auxfiles, 
		       int*    retval         )
{
  int Srv; int* Sretval = &Srv;
	   
  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_const_02 ( gv, hsq, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rrwe_02_GV ( double* hsq,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval,
			    char** set_name         )
{
  gs_esteff_rrwe_02 ( gs_fdta(*set_name),
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_lambda_emstep_01 ( gs_varset_type *gv     ,
			   int*    constvar       ,
			   char**  out_filename   ,
			   int*    auxfiles       , 
			   int*    retval         )
{
  gs_UNSIGNED i,j,m;

  if (  NULL== gv->ladi ) {
    gs_info(gv,-2,"No lambda defined");
    *retval = -2; return ;
  }

  int Srv; int* Sretval = &Srv;

  gs_mme_restcoeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_restrhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  /* --- */

  /* Allocate memory for the solution vector */
  gv->dimu = 1 + gv->NoMar;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );
  if ( NULL == gv->u ) ERR_M;

  /* Initialize the solution vector */
  for ( i = 0; i < gv->dimu ; i++) gv->u [i] =  gv->b[i];

  /* Trace of the inverse coefficient matrix */
  gv->t0 = realloc(gv->t0, gv->dimA * sizeof(double) );

  { // LU decomposition 
    double d;

    gv->u0  = realloc ( gv->u0, gv->dimA * sizeof(gs_UNSIGNED) );

    /* LU decomposition */
     st_dgetrf ( gv->A, gv->dimA, gv->u0, &d);

    /* Solutions for u */
     st_dgetrs ( gv->A, gv->dimA, gv->u0, gv->u );

    /* Inverse */
     gv->t4 = realloc (  gv->t4, gv->dimA*gv->dimA * sizeof(double) );

    st_dgetri ( gv->A, gv->dimA, gv->u0 ,gv->t4 );

    /* Diagonal elements of the inverse coefficient matrix */
     for ( j=0 ; j < gv->dimA ; j++ ) 
       gv->t0[j] = gv->t4[ gv->dimA *j + j];
   
  }


  /* Allocate memory for the variace components */
  gv->sigma =  realloc( gv->sigma, gv->dimA * sizeof(double) );
  if ( NULL == gv->sigma ) ERR_M;

  gv->t1 = realloc (gv->t1, sizeof(double));
  gv->t2 = realloc (gv->t2, ( (gv->NoInd > gv->NoMar) ?
		      gv->NoInd : gv->NoMar    ) 
		    * sizeof(double)            );
  gv->t3 = realloc (gv->t3, ( (gv->NoInd > gv->NoMar) ?
		           gv->NoInd : gv->NoMar    ) 
			 * sizeof(double)            );

  for ( i=0 ; i<gv->NoInd ; i++) gv->t2[i] = 1;

    st_dgemm ("T","N",
	    gv->NoInd,1 ,         // y'y
	      gv->NoInd,1 ,	     
	      gv->y       ,
	      gv->y       ,
	      gv->sigma ) ;

    st_dgemm ("T","N",
	    gv->NoInd,1 ,         // 1'y
	      gv->NoInd,1 ,	     
	      gv->t2          ,
	      gv->y       ,
	      gv->t1          ) ;

    st_dgemm ("T","N",
	    gv->NoInd, gv->NoMar, // Z'y	  	     
	      gv->NoInd, 1 ,
	      gv->Z        ,
	      gv->y        ,
	      gv->t2           ) ;

    st_dgemm ("T","N",
	    gv->NoMar, 1 ,	     
	      gv->NoMar, 1 ,
	      1+(gv->u)    ,
	      gv->t2           ,
	      gv->t3           ) ;
  
  if (0 == *constvar) /* Variable variance of the markers     */
    { /* Residual is held constant, and needs to be estimated */
      /* in a previous run with constant variances            */
       gv->sigma[0] = gv->sigma_e ;
       for ( i=1; i < gv->dimA; i++)
	{
	  gv->sigma[i] = 
	     gv->u[i] * gv->u[i] + gv->sigma[0] * gv->t0[i]  ;
	}
    }
  else /* Constant variance across all markers              */
    {  /* Residual error is estimated in each iteration     */
      gv->sigma[0] -= gv->u[0] *gv->t1[0];
      gv->sigma[0] -= gv->t3[0];                     
      gv->sigma[0]  = gv->sigma[0] / (double)gv->NoInd;
      gv->sigma_e   = gv->sigma[0];

        st_dgemm ("T","N",
	    gv->NoMar,1, gv->NoMar,1, 1+(gv->u), 1+(gv->u), gv->t2 ) ;

      double tr = 0; 
      for ( i=1; i < gv->dimA; i++ ) tr +=gv->t0 [ i ];

      for (i=1; i <gv->dimA; i++)
	gv->sigma[i] = (gv->t2[0] + gv->sigma[0] * tr ) / (double) gv->NoMar;

    }


  for (i=1; i <gv->dimA; i++)
    gv->ladi[i-1] = gv->sigma[0] / gv->sigma[i];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;

      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda var \n");
      fprintf(outfp,"err %f %f \n", 0.0, gv->sigma[0] );
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f %f\n",
		gv->MarName[m],
		gv->ladi[m] ,
		gv->sigma[1+m]
		);

      fclose(outfp);  
    } 

  *retval = 0; return ;

}

void gs_lambda_emstep_01_GV ( int*    constvar      ,
			      char**  out_filename  ,
			      int*    auxfiles      , 
			      int*    retval        ,
			      char** set_name       )

{
  gs_lambda_emstep_01 (  gs_fdta(*set_name),
			 constvar,
			 out_filename,
			 auxfiles, 
			 retval         ) ;
}

void gs_lambda_rml_01 (  gs_varset_type* gv ,
			 int*    maxiter,
			 double* precision,
			 int*    constvar,
			 int*    init_lambda,
			 double* hsq,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval         )
{
  gs_UNSIGNED m;
  int i;
  double qu, maxqu=0;
  int Srv; int* Sretval = &Srv;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available");
    *retval = -2; return ;
  }

  double *vc_save = malloc ((1+gv->NoMar)*sizeof(double) );
  if (NULL==vc_save) ERR_M;

  for ( m = 0 ; m < 1+gv->NoMar ; m++) vc_save[m] = 0;

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if (NULL==gv->ladi) ERR_M;


  if ( 1 == *init_lambda )
    {
      // Initialisation with 1 
      // for ( m = 0 ; m < gv->NoMar ; m++)  
      // gv->ladi[m] = 1 ;
       
      /* Initialization with constant lambda */
      double lambda =  ( 1 / *hsq -1 ) * (double)gv->NoMar; 
      for ( m = 0 ; m < gv->NoMar ; m++)  
	gv->ladi[m] = lambda;
    }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ; 
  if (-2 == *Sretval) { *retval = -2; return; }  

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }
      
  for (i=0; i < *maxiter; i ++) 
    {
      if (i>0)  for ( m = 0 ; m < 1+gv->NoMar ; m++) 
		  vc_save[m] = gv->sigma[m];

      gs_lambda_emstep_01 ( gv,constvar,nixp,neinp,Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      maxqu = 0;
      for ( m = 0 ; m < 1+gv->NoMar ; m++) {
	qu = fabs ( (vc_save[m] - gv->sigma[m]) / gv->sigma[m] );
	if (qu>maxqu) maxqu=qu;
      }
      if (maxqu < *precision) break; 
  }
  
  snprintf(gv->msg, sizeof(gv->msg),
         "Number of iterations: %i, precision: %f", i, maxqu);

  gs_info(gv,1,gv->msg);
  if (1 == *constvar) {
    snprintf(gv->msg, sizeof(gv->msg),
	     "Variance components: var(m) = %f, var(e) = %5f",
	     gv->sigma[1], gv->sigma[0]);

    gs_info(gv,1,gv->msg);
  }

  free(vc_save);

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rmlc_01 (  gs_varset_type *gv ,
			  int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int _cc = 1; int* constvar    = &_cc;
  int _il = 1; int* init_lambda = &_il;

  gs_lambda_rml_01  ( gv,
		      maxiter,
		      precision,
		      constvar,
		      init_lambda,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_lambda_rmlv_01 (  gs_varset_type *gv ,
			  int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{

  int Srv =0; int* Sretval = &Srv;

  int _cc = 1; int* constvar = &_cc;
  int _il = 1; int* init_lambda = &_il;

  /* First run with constant variances to estimate residual */
  /* variancs                                              */
  _cc = 1 ;
  _il = 1;
  gs_lambda_rml_01  ( gv,
		      maxiter,
		      precision,
		      constvar,
		      init_lambda,
		      hsq,
		      nixp,
		      neinp, 
		      Sretval         ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  /* Second run with variable variances to estimate marker */
  /* variances                                             */
  _cc = 0; 
  _il = 1;
  gs_lambda_rml_01  ( gv,
		      maxiter,
		      precision,
		      constvar,
		      init_lambda,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}


void gs_lambda_rmlv_01_GV ( int*    maxiter,
			    double* precision,
			    double* hsq,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval         ,
			    char** set_name         )
{
  gs_lambda_rmlv_01 ( gs_fdta(*set_name),
		      maxiter,
		      precision,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_lambda_rmlc_01_GV ( int*    maxiter,
			    double* precision,
			    double* hsq,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval         ,
			    char** set_name         )
{
  gs_lambda_rmlc_01 ( gs_fdta(*set_name),
		      maxiter,
		      precision,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_lambda_rmlr_01 (  gs_varset_type *gv     ,
			  int*    maxiter        ,
			  double* precision      ,
			  double* hsq            ,
			  char**  out_filename   ,
			  int*    auxfiles       , 
			  int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_reg_01 ( gv, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlc_01 ( gv,
		      maxiter,precision,hsq,
		      out_filename,auxfiles,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  long double SSQM = 0.0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    SSQM += gv->SMsqm[m] ;

  SSQM = SSQM / (long double) gv->NoMar;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] =  ( SSQM / gv->SMsqm[m]  ) * gv->ladi[m];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rmlr_01_GV ( int*    maxiter       ,
			    double* precision     ,
			    double* hsq            ,
			    char**  out_filename  ,
			    int*    auxfiles      , 
			    int*    retval        ,
			    char** set_name        )
{
  gs_lambda_rmlr_01 ( gs_fdta(*set_name),
		      maxiter,
		      precision,
		      hsq,
		      out_filename      ,
		      auxfiles          , 
		      retval            ) ;
}

void gs_lambda_rmla_01 (  gs_varset_type *gv     ,
			  double* alpha          ,
			  int*    maxiter        ,
			  double* precision      ,
			  double* hsq            ,
			  char**  out_filename   ,
			  int*    auxfiles       , 
			  int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_aov_01 ( gv, alpha, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlc_01 ( gv,
		      maxiter,precision,hsq,
		      out_filename,auxfiles,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  long double Svar = 0.0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    Svar += gv->SMvar[m] ;

  Svar = Svar / (long double) gv->NoMar;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] =  ( Svar / gv->SMvar[m]  ) * gv->ladi[m];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rmla_01_GV ( double* alpha,
                            int*    maxiter       ,
			    double* precision     ,
			    double* hsq           ,
			    char**  out_filename  ,
			    int*    auxfiles      , 
			    int*    retval        ,
			    char** set_name        )
{
  gs_lambda_rmla_01 ( gs_fdta(*set_name),
		      alpha,
		      maxiter,
		      precision,
		      hsq,
		      out_filename      ,
		      auxfiles          , 
		      retval            ) ;
}


void gs_esteff_rmlc_01 ( gs_varset_type*, int*, double*,double*,
			 char**, int*, int*);

void gs_esteff_lsq_01 (  gs_varset_type *gv ,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  if ( (1+gv->NoMar) > gv->NoInd ) {
    gs_info(gv,-2,"p > n");  
    *retval = -2; return; }

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  double lambda = 0;
  double * lambdaptr = &lambda;

  gs_lambda_const_01(gv, lambdaptr, nixp, neinp, Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_lsq_01_GV ( char**  out_filename,
			   int*    auxfiles, 
			   int*    retval,
			   char** set_name         )
{
  gs_esteff_lsq_01 ( gs_fdta(*set_name),
		     out_filename,
		     auxfiles, 
		     retval         ) ;
}


void gs_lm_ftest_01 ( gs_varset_type *gv ,
		      char** out_filename,
		      int*   auxfiles, 
		      int*   retval        )
/* F-test for full rank model */
{
  int Srv; int* Sretval = &Srv;
  gs_UNSIGNED i,j;
  FILE *outfp;

  /* Allocate memory for the p-values */
  double *Q_val  =  malloc( gv->dimu * sizeof(double) );
  double *F_val  =  malloc( gv->dimu * sizeof(double) );
  gv->p = realloc (gv->p, gv->dimu * sizeof(double) );
  if ( NULL==Q_val ||  NULL==F_val  ||  NULL==gv->p ) ERR_M;

  /* Allocate memory for the contrast vector */
  double *k  =  malloc( gv->dimu * sizeof(double) );
  if ( NULL == k ) ERR_M;

  /* A is X'X */
  gs_mme_restcoeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  /* (X'X)⁻1  is the inverse of A */
  double * ixpx = malloc( gv->dimu * gv->dimu * sizeof(double) );
  if ( NULL == ixpx ) ERR_M;

  /* Invert A */
  /*
  {
    double **a, **y, *col, d;
    gs_UNSIGNED *indx;

    indx=ivector(1,gv->dimA);
    col = vector(1,gv->dimA);
    a = convert_matrix( &(gv->A[0]), 1, gv->dimA, 1, gv->dimA);
    y = convert_matrix( &(ixpx[0]),  1, gv->dimA, 1, gv->dimA);

    ludcmp(a,gv->dimA,indx,&d);

    for ( j=1; j <= gv->dimA; j++ ){
      for ( i=1; i <= gv->dimA; i++) col[i] = 0.0;
      col[j] = 1.0 ;
      lubksb(a,gv->dimA,indx,col);
      for (i=1; i <= gv->dimA; i++) y[i][j] = col[i];
    }

    free_ivector(indx,1);
    free_vector(col,1);
    free_convert_matrix(a,1);
    free_convert_matrix(y,1);
  }
  */

  st_dgetri_inv ( gv->A, ixpx, gv->dimA);

  
  free(gv->A);
  gv->A =NULL;

  /* b ( from Au=b) is X'y */
  gs_mme_restrhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_UNSIGNED  dfe;
  double kb, kxxk, yy, bxy, sse;
  double * kixtx = malloc ( gv->dimA *sizeof(double) );
  
    st_dgemm ("T","N",
	    gv->NoInd,1,gv->NoInd,1,gv->y,gv->y,&yy);
    st_dgemm ("T","N",
	    gv->dimu,1,gv->dimu,1,gv->u,gv->b,&bxy);
  sse = yy - bxy;
  dfe = gv->NoInd - (gv->NoMar+1);

  snprintf(gv->msg, sizeof(gv->msg),
	   "SSE = %9f, DFE = %lu", sse, dfe);
  gs_info(gv,1,gv->msg);

  /* Calculate p-values */
  for ( i = 0; i < gv->dimu ; i++) 
    {
      for ( j = 0; j < gv->dimu ; j++)  k[j]=0;
      k[i]=1;
      
        st_dgemm ("T","N",
	     gv->dimA,1, gv->dimA,1,        k, gv->u,  &kb );
        st_dgemm ("T","N",
	     gv->dimA,1, gv->dimA,gv->dimA, k, ixpx,   kixtx);
        st_dgemm ("T","N",
	     gv->dimA,1, gv->dimA,1,        kixtx, k,  &kxxk);
      Q_val[i] = kb * (1/kxxk) * kb ;
      F_val[i] = Q_val[i] / (sse/(double)dfe);
      gv->p[i] = st_fprob(F_val[i],1,dfe);
    }

  /* Write solution vector and p-values to output file */
  if (1 == *auxfiles) {

    char sig[5];
   
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"effect ssh f pvalue sig\n");
      for ( i = 0; i < gv->dimu ; i++)
	{
	  if (gv->p[i]<0.001) strcpy(sig,"***"); 
	  else if (gv->p[i]<0.01) strcpy(sig,"** "); 
	  else if (gv->p[i]<0.05) strcpy(sig,"*  ");
	  else if (gv->p[i]<0.10) strcpy(sig,".  ");
	  else  strcpy(sig,"");

	fprintf(outfp,"%s.%i %f %f %f %f \"%s\"\n",
		gv->EffNme[i],
		gv->EffAll[i],
		gv->u[i],
		Q_val[i],
		F_val[i],
		gv->p[i],
		sig     );
	}

      fclose(outfp);  
  }

  free(Q_val); 
  free(F_val); 
  free(k); 
  free(ixpx); 
  free(kixtx);

  *retval = 0; return ;
}

void gs_esteff_rmlc_01 (  gs_varset_type *gv ,
                          int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlc_01 ( gv,
		      maxiter,precision,hsq,
		      nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rmlc_01_GV (  int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_esteff_rmlc_01 ( gs_fdta(*set_name),
                      maxiter,   
		      precision, 
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_esteff_rmlv_01 (  gs_varset_type *gv ,
                          int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlv_01 ( gv,
		      maxiter,precision,hsq,
		      nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rmlv_01_GV (  int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_esteff_rmlv_01 ( gs_fdta(*set_name),
                      maxiter,   
		      precision,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}


void gs_permtest_01 ( gs_varset_type *gv ,
                      char *  scheme       ,
                      int*    nperm,
                      int*    maxiter,
                      double* precision,
                      double* hsq,
                      char** out_filename,
                      int*   auxfiles,
                      int*   retval        )
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
      gs_info(gv,-2,"More individuals than INT_MAX");
      return;
    }

  if ((0 == gv->NoInd) || (0 == gv->NoMar) ||
      (NULL == gv->All1) || (NULL == gv->All2))
    {
      gs_info(gv,-2,"No marker data available");
      return;
    }

  if (*nperm < 1)
    {
      gs_info(gv,-2,"Number of permutations must be positive");
      return;
    }

  if ((0 != strcmp(scheme,"lsq")) &&
      (0 != strcmp(scheme,"rmlc")) &&
      (0 != strcmp(scheme,"rmlv")))
    {
      gs_info(gv,-2,"Unknown permutation-test scheme");
      return;
    }

  if ((0 == gv->dimu) || (NULL == gv->u))
    {
      gs_info(gv,-2,"No estimated effects available");
      return;
    }

  dimu = (size_t)gv->dimu;
  if ((gs_UNSIGNED)dimu != gv->dimu)
    {
      gs_info(gv,-2,"Permutation result size is too large");
      return;
    }

  if ((size_t)*nperm > ((size_t)-1) / dimu)
    {
      gs_info(gv,-2,"Permutation result size is too large");
      return;
    }

  eff_count = (size_t)*nperm * dimu;

  if (eff_count > ((size_t)-1) / sizeof(*eff_st))
    {
      gs_info(gv,-2,"Permutation result size is too large");
      return;
    }

  eff_st = calloc(eff_count,sizeof(*eff_st));
  if (NULL == eff_st)
    {
      gs_info(gv,-2,"Memory allocation error");
      return;
    }

#ifdef MT
  #pragma omp parallel
  {
    #pragma omp single
    {
      ntr = omp_get_num_threads();
    }
  }
#endif

  if (ntr < 1)
    ntr = 1;

  if ((size_t)ntr > ((size_t)-1) / sizeof(*thread_es))
    {
      free(eff_st);
      gs_info(gv,-2,"Memory allocation error");
      return;
    }

  thread_es = malloc((size_t)ntr * sizeof(*thread_es));
  if (NULL == thread_es)
    {
      free(eff_st);
      gs_info(gv,-2,"Memory allocation error");
      return;
    }

  /* Serialize because data set management uses global variables. */
  gs_get_info_level(sh,&sh_info_level);
  gs_set_info_level(sh,imz);

  for (int t=0; t < ntr; t++)
    {
      int Srv = 0;
      int *Sretval = &Srv;
      char esname[24];

      snprintf(esname,sizeof(esname),"dta_t%03i",t);
      thread_es[t] = gs_fdta(esname);

      if (NULL == thread_es[t])
        {
          gs_set_info_level(sh,&sh_info_level);
          free(thread_es);
          free(eff_st);
          gs_info(gv,-2,"Memory allocation error");
          return;
        }

      thread_es[t]->info_level = -2;

      gs_copy_marker_data_01(thread_es[t],gv,Sretval);
      if (-2 == *Sretval)
        {
          gs_set_info_level(sh,&sh_info_level);
          free(thread_es);
          free(eff_st);
          return;
        }
    }

  gs_set_info_level(sh,&sh_info_level);

#if defined(use_openblas) && defined(MT)
  {
    int eins = 1;
    openblas_set_num_threads(eins);
    goto_set_num_threads(eins);
    omp_set_num_threads(ntr);
  }
#endif

#if defined(use_mkl) && defined(MT)
  {
    int eins = 1;
    mkl_set_num_threads(eins);
    omp_set_num_threads(ntr);
  }
#endif

#ifdef MT
  #pragma omp parallel reduction(|:parallel_error)
#endif
  {
    gs_UNSIGNED i,j,idx1,idx2;
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
    #pragma omp for schedule(dynamic)
#endif
    for (int r=0; r < *nperm; r++)
      {
        gs_UNSIGNED no_unsucc_rep = 0;
        int success_rep = 0;

        if (local_error)
          continue;

        while ((0 == success_rep) &&
               (no_unsucc_rep < 1000))
          {
            success_rep = 1;

            for (i=1; i < gv->dimu; i++)
              {
                /* Random permutation of the individuals. */
                p[0] = 0;
                for (int ii=1; ii < (int)gv->NoInd; ii++)
                  {
                    int jj = st_rand_index(ii + 1);
                    p[ii] = p[jj];
                    p[jj] = ii;
                  }

                /* Randomize the alleles of marker i-1. */
                for (j=0; j < gv->NoInd; j++)
                  {
                    idx1 = (i-1) * gv->NoInd + j;
                    idx2 = (i-1) * gv->NoInd + p[j];
                    es->All1[idx1] = gv->All1[idx2];
                    es->All2[idx1] = gv->All2[idx2];
                  }

                *Sretval = 0;

                if (0 == strcmp(scheme,"lsq"))
                  gs_esteff_lsq_01(es,nixp,neinp,Sretval);
                else if (0 == strcmp(scheme,"rmlc"))
                  gs_esteff_rmlc_01(es,maxiter,precision,hsq,
                                     nixp,neinp,Sretval);
                else
                  gs_esteff_rmlv_01(es,maxiter,precision,hsq,
                                     nixp,neinp,Sretval);

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
  } // parallel

#if defined(use_openblas) && defined(MT)
  openblas_set_num_threads(ntr);
  goto_set_num_threads(ntr);
#endif

#if defined(use_mkl) && defined(MT)
  mkl_set_num_threads(ntr);
#endif

  free(thread_es);
  thread_es = NULL;

  if (parallel_error)
    {
      free(eff_st);
      gs_info(gv,-2,"Permutation test failed in parallel calculation");
      return;
    }

#ifdef MT
  #pragma omp parallel for
#endif
  for (gs_UNSIGNED i=0; i < gv->dimu; i++)
    qsort(&(eff_st[(size_t)i * (size_t)*nperm]),
          (size_t)*nperm,sizeof(double),gs_dbl_cmp);

  if (dimu > ((size_t)-1) / sizeof(double))
    {
      free(eff_st);
      gs_info(gv,-2,"Permutation result size is too large");
      return;
    }

  p_seq = malloc(dimu * sizeof(*p_seq));
  p_leq = malloc(dimu * sizeof(*p_leq));

  if ((NULL == p_seq) || (NULL == p_leq))
    {
      free(p_seq);
      free(p_leq);
      free(eff_st);
      gs_info(gv,-2,"Memory allocation error");
      return;
    }

  {
    double *tmp_p =
      realloc(gv->p,dimu * sizeof(*gv->p));

    if (NULL == tmp_p)
      {
        free(p_seq);
        free(p_leq);
        free(eff_st);
        gs_info(gv,-2,"Memory allocation error");
        return;
      }

    gv->p = tmp_p;
  }

#ifdef MT
  #pragma omp parallel for
#endif
  for (gs_UNSIGNED i=0; i < gv->dimu; i++)
    {
      gs_UNSIGNED r;

      for (r=0; r < (gs_UNSIGNED)*nperm; r++)
        if (gv->u[i] <
            eff_st[(size_t)i * (size_t)*nperm + (size_t)r]) break;
      p_seq[i] = (double)r;

      for (r=0; r < (gs_UNSIGNED)*nperm; r++)
        if (gv->u[i] >
            eff_st[(size_t)i * (size_t)*nperm +
                   (size_t)(*nperm-1) - (size_t)r]) break;
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

      outfp = fopen(*out_filename,"w");
      if (NULL == outfp)
        {
          free(p_seq);
          free(p_leq);
          free(eff_st);
          gs_info(gv,-2,"File open error");
          return;
        }

      fprintf(outfp,"effect  p.seq p.leq pvalue sig\n");
      for (gs_UNSIGNED i=0; i < gv->dimu; i++)
        {
          if (gv->p[i]<0.001) strcpy(sig,"***");
          else if (gv->p[i]<0.01) strcpy(sig,"** ");
          else if (gv->p[i]<0.05) strcpy(sig,"*  ");
          else if (gv->p[i]<0.10) strcpy(sig,".  ");
          else strcpy(sig,"");

          fprintf(outfp,"%s.%i %f %f %f %f \"%s\"\n",
                  gv->EffNme[i],
                  gv->EffAll[i],
                  gv->u[i],
                  p_seq[i],
                  p_leq[i],
                  gv->p[i],
                  sig);
        }

      fclose(outfp);
    }

  free(p_seq);
  free(p_leq);
  free(eff_st);

  *retval = 0;
}


void gs_testeff_lsq_01 (  gs_varset_type *gv ,
			  int*    nperm,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_esteff_lsq_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  if ( 0 == *nperm)
    {
      gs_lm_ftest_01 ( gv, out_filename, auxfiles, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }
  else
    {
      char scheme[]="lsq";
      int mi=0; int* maxiter = &mi;
      double pr=0; double* precision = &pr;
      double _hsq = 0.9; double* hsq=&_hsq;
      gs_permtest_01 ( gv,scheme,nperm,maxiter,precision,hsq,
		       out_filename,auxfiles,Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  *retval = 0; return ;
}

void gs_testeff_lsq_01_GV ( int*    nperm,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval,
			    char** set_name         )
{
  gs_testeff_lsq_01 ( gs_fdta(*set_name),
		      nperm,
	              out_filename,
		      auxfiles, 
		      retval         ) ;
}


void gs_testeff_rmlc_01 (  gs_varset_type *gv ,
			   int*    nperm,
			   int*    maxiter,
			   double* precision,
			   double* hsq,
			   char**  out_filename,
			   int*    auxfiles, 
			   int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_esteff_rmlc_01 ( gv, maxiter, precision,hsq, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  char scheme[]="rmlc";
  gs_permtest_01 (gv,scheme,nperm,maxiter,precision,hsq,out_filename,
		  auxfiles,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}


void gs_testeff_rmlc_01_GV ( int*    nperm, 
			     int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_testeff_rmlc_01 ( gs_fdta(*set_name),
		       nperm,
		       maxiter,   
		       precision, 
		       hsq,
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}

void gs_testeff_rmlv_01 (  gs_varset_type *gv ,
			   int*    nperm,
			   int*    maxiter,
			   double* precision,
			   double* hsq,
			   char**  out_filename,
			   int*    auxfiles, 
			   int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_esteff_rmlv_01 ( gv, maxiter, precision, hsq, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  char scheme[]="rmlv";
  gs_permtest_01 (gv,scheme,nperm,maxiter,precision,hsq,
		  out_filename,auxfiles,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}


void gs_testeff_rmlv_01_GV ( int*    nperm, 
			     int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_testeff_rmlv_01 ( gs_fdta(*set_name),
		       nperm,
		       maxiter,   
		       precision, 
		       hsq,
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}


void gs_esteff_rmlr_01 (  gs_varset_type *gv ,
                          int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlr_01 ( gv,
		      maxiter,precision,hsq,
		      nixp,neinp,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rmlr_01_GV (  int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_esteff_rmlr_01 ( gs_fdta(*set_name),
                      maxiter,   
		      precision, 
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_esteff_rmla_01 (  gs_varset_type *gv ,
			  double* alpha,
                          int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmla_01 ( gv,
		      alpha,
		      maxiter,precision,hsq,
		      nixp,neinp,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rmla_01_GV (  double* alpha,
			     int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_esteff_rmla_01 ( gs_fdta(*set_name),
		      alpha,
                      maxiter,   
		      precision, 
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_lambda_emstep_tc_01 ( gs_varset_type *gv     ,
			      char**  out_filename   ,
			      int*    auxfiles       , 
			      int*    retval         )
{

  gs_UNSIGNED i,j,m;

  if (  NULL== gv->ladi ) {
    gs_info(gv,-2,"No lambda defined");
    *retval = -2; return ;
  }

  int Srv; int* Sretval = &Srv;

  gs_mmet_coeff_02 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mmet_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  /* Allocate memory for the solution vector */
  gv->a =  realloc( gv->a, gv->dimA * sizeof(double) );
  if ( NULL == gv->a ) ERR_M;

  /* Initialize the solution vector */
  for ( i = 0; i < gv->dimA ; i++) gv->a [i] =  gv->b[i];

  { // LU decomposition 
    double d;

    gv->u0  = realloc ( gv->u0, gv->dimA * sizeof(gs_UNSIGNED) );

    /* LU decomposition */
     st_dgetrf ( gv->A, gv->dimA, gv->u0, &d);

    /* Solutions for a */
     st_dgetrs ( gv->A, gv->dimA, gv->u0, gv->a );

    /* Inverse */
     gv->t4 = realloc (  gv->t4, gv->dimA*gv->dimA * sizeof(double) );

    st_dgetri ( gv->A, gv->dimA, gv->u0 ,gv->t4 );

  }

  /* Allocate memory for the variance components */
  gv->sigma =  realloc( gv->sigma, gv->NoMar * sizeof(double) );
  if ( NULL == gv->sigma ) ERR_M;

  gv->t1 = realloc (gv->t1, sizeof(double));

  gv->t2 = realloc (gv->t2, 
		    ( (gv->NoInd > gv->NoMar) ?  gv->NoInd : gv->NoMar    ) 
		    * sizeof(double)            );

  gv->t3 = realloc (gv->t3, 
		    ( (gv->NoInd > gv->NoMar) ?  gv->NoInd : gv->NoMar    ) 
		    * sizeof(double)            );

  for ( i=0 ; i<gv->NoInd ; i++) gv->t2[i] = 1;

    st_dgemm ("T","N",
	    gv->NoInd,1 ,           // y'y
	         gv->NoInd,1 ,	     
	         gv->y       ,
	         gv->y       ,
	         gv->sigma ) ;

    st_dgemm ("T","N",
	    gv->NoInd,1 ,           // 1'y
	         gv->NoInd,1 ,	     
	         gv->t2          ,
	         gv->y       ,
	         gv->t1          ) ;  

    st_dgemm ("T","N",
	    gv->NoInd, 1 ,	         // a'y
	         gv->NoInd, 1 ,
	         1+(gv->a)    ,
	         gv->y        ,
	         gv->t3           ) ;
  
  // Residual error is estimated 

  gv->sigma[0] -= gv->a[0] *gv->t1[0];  // y'y - b0 1' y
  gv->sigma[0] -= gv->t3[0];            //     - a'Z'y         
  gv->sigma[0]  = gv->sigma[0] / (double)gv->NoInd;
  gv->sigma_e   = gv->sigma[0];

    st_dgemm ("T","N",
	    gv->NoInd, 1,           // a'GI	     
	         gv->NoInd, gv->NoInd ,
	         1+(gv->a)            ,
	         gv->GI               ,
	         gv->t3                ) ;

  st_dgemm ("N","T",
                 1,         gv->NoInd ,  // a' GI a
	         1,         gv->NoInd , 
	         gv->t3               ,
	         1+(gv->a)            ,
	         gv->t2                ) ;

  /* Diagonal elements of the inverse coefficient matrix */

  gv->t5 = realloc ( gv->t5, gv->NoInd * gv->NoInd * sizeof(double) );
  for ( j=1 ; j < gv->dimA ; j++ ) 
    memcpy ( & gv->t5 [ (j-1) * gv->NoInd ] ,
	     & gv->t4 [  j    * gv->dimA   + 1],
	     (gv->NoInd) * sizeof(double) );

  gv->t6 = realloc ( gv->t6, gv->NoInd * gv->NoInd * sizeof(double) );
  st_dgemm ("N","N",
	      gv->NoInd, gv->NoInd ,          
	      gv->NoInd, gv->NoInd ,
	      gv->GI               ,
	      gv->t5               ,
	      gv->t6                ) ;

  double tr = 0; 
  for ( i=0; i < gv->NoInd; i++ ) tr +=gv->t6 [  gv->NoInd *i + i ];

  gv->sigma[1] = (gv->t2[0] + gv->sigma[0] * tr ) / (double) gv->NoInd;

  // Update the shrinkage factor
  for (i=0; i <gv->NoMar; i++)
    gv->ladi[i] = gv->sigma[0] / gv->sigma[1] ;

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;

      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda var \n");
      fprintf(outfp,"err %f %f \n", 0.0, gv->sigma[0] );
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f %f\n",
		gv->MarName[m],
		gv->ladi[m] ,
		gv->sigma[1+m]
		);

      fclose(outfp);  
    } 

  *retval = 0; return ;

}

void gs_lambda_emstep_tc_01_GV ( char**  out_filename  ,
			         int*    auxfiles      , 
			         int*    retval        ,
			         char** set_name       )

{
  gs_lambda_emstep_tc_01 (  gs_fdta(*set_name),
			    out_filename,
			    auxfiles, 
			    retval         ) ;
}

void gs_lambda_rmlct_01 ( gs_varset_type* gv ,
		 	  int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{

  gs_UNSIGNED m;
  int i;
  double qu, maxqu=0;
  int Srv; int* Sretval = &Srv;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available");
    *retval = -2; return ;
  }

  // Save variance components to determine convergence
  double vc_save[2] = {0,0};

  // Initialize Vector of shrinkage factors lambda
  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if (NULL==gv->ladi) ERR_M;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    gv->ladi[m] = ( 1 / *hsq - 1 ) * (double)gv->NoInd;

  // Calculate  inv(ZZ')
  gs_build_GI_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }
      
  for (i=0; i < *maxiter; i ++) 
    {

      if (i>0) memcpy ( vc_save, gv->sigma,  2 * sizeof(double) ); 

      gs_lambda_emstep_tc_01 ( gv,nixp,neinp,Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      maxqu = 0;
      for ( m = 0 ; m < 2 ; m++) {
	qu = fabs ( (vc_save[m] - gv->sigma[m]) / gv->sigma[m] );
	if (qu>maxqu) maxqu=qu;
      }
      if (maxqu < *precision) break; 
  }
  
  snprintf(gv->msg, sizeof(gv->msg),
         "Number of EM iterations: %i, precision: %f", i, maxqu);
  gs_info(gv,1,gv->msg);
  snprintf(gv->msg, sizeof(gv->msg),
         "Variance components: var(m) = %f, var(e) = %5f",
         gv->sigma[1], gv->sigma[0]);
  gs_info(gv,1,gv->msg);

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}


void gs_lambda_rmlct_01_GV ( int*    maxiter        ,
			     double* precision      ,
			     double* hsq            ,
			     char**  out_filename   ,
			     int*    auxfiles       , 
			     int*    retval         ,
			     char**  set_name         )
{
  gs_lambda_rmlct_01 ( gs_fdta(*set_name),
		       maxiter,
		       precision,
		       hsq,
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}

void gs_lambda_rmlat_01 (  gs_varset_type *gv     ,
			   double* alpha          ,
			   int*    maxiter        ,
			   double* precision      ,
			   double* hsq            ,
			   char**  out_filename   ,
			   int*    auxfiles       , 
			   int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_aov_01 ( gv, alpha, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlct_01 ( gv,
		       maxiter,precision,hsq,
		      out_filename,auxfiles,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  long double Svar = 0.0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    Svar += gv->SMvar[m] ;

  Svar = Svar / (long double) gv->NoMar;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] =  ( Svar / gv->SMvar[m]  ) * gv->ladi[m];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rmlat_01_GV ( double* alpha,
			     int*    maxiter       ,
			     double* precision     ,
			     double* hsq           ,
			     char**  out_filename  ,
			     int*    auxfiles      , 
			     int*    retval        ,
			     char** set_name        )
{
  gs_lambda_rmlat_01 ( gs_fdta(*set_name),
		       alpha             ,
		       maxiter           ,
		       precision         ,
		       hsq               ,
		       out_filename      ,
		       auxfiles          , 
		       retval             ) ;
}


//##############################################################################
//# End:   EM based genomic prediction                                   2012-18
//##############################################################################

//##############################################################################
//# Begin: C interfact between genomic selection and simulations            2019
//##############################################################################

static int st_sim_marker_data_valid(gs_varset_type *gv)
{
  if ((0 == gv->NoInd) ||
      (0 == gv->NoMar) ||
      (NULL == gv->All1) ||
      (NULL == gv->All2) ||
      (NULL == gv->MarName))
    {
      gs_info(gv,-2,"No marker data available");
      return 0;
    }

  if (!gs_marker_genotypes_valid(gv))
    {
      gs_info(gv,-2,"Invalid marker genotype representation");
      return 0;
    }

  return 1;
}


static int st_sim_data_valid(gs_varset_type *gv)
{
  gs_UNSIGNED *storage;

  if (!st_sim_marker_data_valid(gv))
    return 0;

  storage = gs_map_storage_index(gv,1,1);

  if (NULL == storage)
    return 0;

  free(storage);

  return 1;
}


void st_set_sim_gp_provide_data ( int*     n,           // Number of chromosomes
				  double** chrlen,      // Chromosome lengths
			          int*     retval,      // Return value
			          char**   set_name     // Name of the data set
			         )
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

  /* Determine chromosome range without assuming map order. */
  for (lne=0; lne < gv->NLoc; lne++)
    {
      chrom = gv->map[lne]->chrom;

      if ((chrom < 1) ||
          !isfinite(gv->map[lne]->pos) ||
          (gv->map[lne]->pos < 0))
        {
          gs_info(gv,-2,
                  "Invalid chromosome or position in linkage map");
          return;
        }

      if (chrom > maxchrom)
        maxchrom = chrom;
    }

  if (maxchrom < 1)
    {
      gs_info(gv,-2,"No valid chromosomes in linkage map");
      return;
    }

  if ((size_t)maxchrom >
      ((size_t)-1) / sizeof(double))
    {
      gs_info(gv,-2,"Memory allocation error");
      return;
    }

  *chrlen = malloc((size_t)maxchrom * sizeof(double));
  maxpos = calloc((size_t)maxchrom, sizeof(double));
  seen = calloc((size_t)maxchrom, sizeof(unsigned char));

  if ((NULL == *chrlen) ||
      (NULL == maxpos) ||
      (NULL == seen))
    {
      gs_info(gv,-2,"Memory allocation error");
      goto fail;
    }

  for (lne=0; lne < gv->NLoc; lne++)
    {
      chrom = gv->map[lne]->chrom - 1;
      seen[chrom] = 1;

      if (gv->map[lne]->pos > maxpos[chrom])
        maxpos[chrom] = gv->map[lne]->pos;
    }

  for (chrom=0; chrom < maxchrom; chrom++)
    {
      if (!seen[chrom])
        {
          gs_info(gv,-2,
                  "Chromosome numbering in linkage map is not contiguous");
          goto fail;
        }

      (*chrlen)[chrom] =
        0.01 + maxpos[chrom] / 100.0;

      if (!isfinite((*chrlen)[chrom]) ||
          ((*chrlen)[chrom] <= 0))
        {
          gs_info(gv,-2,
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


void st_set_sim_mp_provide_data (  gs_UNSIGNED       *  NLoc, // Number of loci
                                   gs_mappoint_type *** map , // Map
                                   int*     retval,      // Return value
 			           char**   set_name     // Name of the data set
			         )
{

  gs_varset_type* gv = gs_fdta(*set_name);   

  if (NULL == gv->map)
    {
      snprintf(gv->msg, sizeof(gv->msg), "No map data available");
      gs_info(gv,-2,gv->msg);
      *retval = -2;
      return;
    }

  {
    gs_UNSIGNED *storage = gs_map_storage_index(gv,0,0);

    if (NULL == storage)
      {
        *retval = -2;
        return;
      }
    free(storage);
  }

  *NLoc = gv->NLoc;
  *map  = gv->map;

  *retval = 0;   return ;
}

static int st_sim_effect_data_state(gs_varset_type *gv)
{
  if ((0 == gv->NoEff) ||
      (0 == gv->dimu) ||
      (NULL == gv->u))
    return 0;

  if ((0 == gv->NoMar) ||
      (NULL == gv->MarName) ||
      !gs_restrict_marker_effects_valid(gv))
    return -1;

  return 1;
}


void st_set_sim_ef_provide_data (  gs_UNSIGNED       *  NoEff,
                                   char             (** EffNme)[gs_MAXNME],
                                   int              ** EffAll,
                                   int              ** CmpAll,
                                   double           ** u,
                                   int               * retval,
                                   char             ** set_name )
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
      gs_info(gv,-2,"Invalid marker effect data");
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


void st_set_sim_pp_provide_data  (  gs_UNSIGNED       *  NoInd,
                                    gs_UNSIGNED       *  NoMar,
                                    char             (** MarName)[gs_MAXNME],
                                    GS_ALLELE_TYPE   ** All1,
                                    GS_ALLELE_TYPE   ** All2,
                                    int               * retval,
                                    char             ** set_name )
{
  gs_varset_type *gv = gs_fdta(*set_name);

  *retval = -2;
  *NoInd = 0;
  *NoMar = 0;
  *MarName = NULL;
  *All1 = NULL;
  *All2 = NULL;

  if (!st_sim_marker_data_valid(gv))
    return;

  *NoInd = gv->NoInd;
  *NoMar = gv->NoMar;
  *MarName = gv->MarName;
  *All1 = gv->All1;
  *All2 = gv->All2;

  *retval = 0;
  return;
}



//##############################################################################
//# End: C interface between genomic selection and simulations              2019
//##############################################################################

//##############################################################################
//# Begin: Write performance data to a marker data set,                     2022
//##############################################################################


void gs_set_performance_data_01 ( gs_varset_type *gv ,
				  int*     n,  
				  double*  y,
         			  int*   retval         )
{

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  if ( *n != (int)gv->NoInd )  {
    gs_info(gv,-2,"Incorrect numer of observations");
    *retval = -2; return ;
  }

  // Allocate memory
  gv->y = realloc( gv->y, gv->NoInd * sizeof(double) );
  if ( NULL== gv->y ) ERR_M;

  // Copy performance data
  st_dcopy( gv->NoInd, gv->y, y );
    
  *retval = 0; return ;

}

void gs_set_performance_data_01_GV (  int*    n,  
				      double* y,
				      int*    retval,
				      char**  set_name         )
{
  gs_set_performance_data_01 ( gs_fdta(*set_name),
			       n,
			       y,
			       retval) ;
}

//##############################################################################
//# End: Write performance data to a marker data set,                     2022
//##############################################################################



//##############################################################################
//# BEGIN: Extend st_mixed for full matrix based genomic precistion         2026
//##############################################################################

static int is_missing(double v) { return !(v == v); } /* NA_real_ is NaN */


/* Copy Z from R column-major (nG×m) into row-major Zrm (nG×m). */
static void Z_cm_to_rm(const double *Zcm, int nG, int m, double *Zrm) {
  for (int i = 0; i < nG; ++i) {
    for (int j = 0; j < m; ++j) {
      Zrm[(size_t)i*(size_t)m + (size_t)j] = Zcm[(size_t)i + (size_t)nG*(size_t)j];
    }
  }
}

/* ============================================================
   K building and coding handling
   ============================================================ */

typedef enum {
  K_STD = 0,
  K_VANRADEN = 1,
  K_SOMMER = 2
} KMethod;

typedef enum { CODING_012=0, CODING_M101=1, CODING_001=2, CODING_OTHER=3 } ZCoding;
typedef enum { TRANS_NONE=0, TRANS_ADD1=1, TRANS_MUL2=2 } ZTransform;

static int is_in_set_012(double v)  { return (v==0.0 || v==1.0 || v==2.0); }
static int is_in_set_m101(double v) { return (v==-1.0 || v==0.0 || v==1.0); }
static int is_in_set_001(double v)  { return (v==0.0 || v==0.5 || v==1.0); }

static const char* coding_name(ZCoding c) {
  switch (c) {
    case CODING_012:  return "0/1/2";
    case CODING_M101: return "-1/0/1";
    case CODING_001:  return "0/0.5/1";
    default:          return "OTHER";
  }
}

static const char* transform_name(ZTransform t) {
  switch (t) {
    case TRANS_NONE: return "none";
    case TRANS_ADD1: return "+1";
    case TRANS_MUL2: return "*2";
    default:         return "unknown";
  }
}

/* Strict compatibility check: all finite values in Zrm belong to given coding */
static int Z_is_compatible_with(const double *Zrm, int nG, int m, ZCoding coding) {
  const size_t nm = (size_t)nG * (size_t)m;
  for (size_t t = 0; t < nm; ++t) {
    double v = Zrm[t];
    if (!(v == v)) continue;
    int ok = 0;
    if (coding == CODING_012)  ok = is_in_set_012(v);
    if (coding == CODING_M101) ok = is_in_set_m101(v);
    if (coding == CODING_001)  ok = is_in_set_001(v);
    if (!ok) return 0;
  }
  return 1;
}

/* detect but ERROR on ambiguity (e.g. only {0,1}) */
static int detect_coding_strict_or_ambiguous(const double *Zrm,
					     int nG, int m, ZCoding *out) {
  int ok012  = Z_is_compatible_with(Zrm, nG, m, CODING_012);
  int okm101 = Z_is_compatible_with(Zrm, nG, m, CODING_M101);
  int ok001  = Z_is_compatible_with(Zrm, nG, m, CODING_001);
  int n_ok = ok012 + okm101 + ok001;

  if (n_ok == 0) return -2001;   /* not recognized */
  if (n_ok > 1)  return -2004;   /* ambiguous */

  if (ok012) *out = CODING_012;
  else if (okm101) *out = CODING_M101;
  else *out = CODING_001;
  return 0;
}

/* parse coding string */
static int parse_coding(SEXP coding_, int *is_detect, ZCoding *out) {
  if (!isString(coding_) || XLENGTH(coding_) != 1) return -2005;
  const char *s = CHAR(STRING_ELT(coding_, 0));
  if (strcmp(s, "detect") == 0) { *is_detect = 1; return 0; }
  *is_detect = 0;
  if (strcmp(s, "0 1 2") == 0) { *out = CODING_012; return 0; }
  if (strcmp(s, "-1 0 1") == 0) { *out = CODING_M101; return 0; }
  if (strcmp(s, "0 0.5 1") == 0) { *out = CODING_001; return 0; }
  return -2005;
}


/* Build K; coding controls how Z is interpreted (strict), then transform to dosage if needed */
static int build_K_from_Z(const double *Zrm, int nG, int m,
                          KMethod method,
                          SEXP coding_,
                          double *Krm_out,
                          double *out_z_mean,
                          double *out_z_sd,
                          ZCoding *out_coding_used,
                          ZTransform *out_transform)
{
  const size_t nm = (size_t)nG * (size_t)m;

  int is_detect = 0;
  ZCoding coding_forced = CODING_OTHER;
  int prc = parse_coding(coding_, &is_detect, &coding_forced);
  if (prc != 0) return -2005;

  ZCoding coding_used = CODING_OTHER;
  if (is_detect) {
    int drc = detect_coding_strict_or_ambiguous(Zrm, nG, m, &coding_used);
    if (drc != 0) return drc;
  } else {
    /* strict enforce chosen coding */
    if (!Z_is_compatible_with(Zrm, nG, m, coding_forced)) return -2006;
    coding_used = coding_forced;
  }
  if (out_coding_used) *out_coding_used = coding_used;

  /* Copy Z to work buffer */
  double *Zwork = (double*)malloc(nm * sizeof(double));
  if (!Zwork) return -1010;
  memcpy(Zwork, Zrm, nm * sizeof(double));


  if (method == K_STD) { // For later ...

    /*    
    st_dsymm_ZZt((gs_UNSIGNED)nG, (gs_UNSIGNED)m, Zwork, Krm_out);

    double tr = 0.0;
    for (int i = 0; i < nG; ++i)
      tr += Krm_out[i*nG + i];     
    
    double scale = tr / (double)nG;
    
    for (size_t t = 0; t < nG*nG; ++t)
      Krm_out[t] /= scale;
  
    
    free(Zwork);
    return 0;
  */
  }


  /* Transform to dosage scale 0/1/2 (only finite entries) */
  ZTransform tr = TRANS_NONE;

  if (coding_used == CODING_M101) {
    /* -1/0/1 -> 0/1/2 */
    tr = TRANS_ADD1;
    for (size_t t = 0; t < nm; ++t) if (Zwork[t] == Zwork[t]) Zwork[t] += 1.0;
  } else if (coding_used == CODING_001) {
    /* 0/0.5/1 -> 0/1/2 */
    tr = TRANS_MUL2;
    for (size_t t = 0; t < nm; ++t) if (Zwork[t] == Zwork[t]) Zwork[t] *= 2.0;
  } else {
    /* CODING_012 already dosage */
    tr = TRANS_NONE;
  }
  if (out_transform) *out_transform = tr;

  /* summary after transform (now always on dosage scale for these methods) */
  {
    double sum = 0.0, sumsq = 0.0;
    size_t cnt = 0;
    for (size_t t = 0; t < nm; ++t) {
      double v = Zwork[t];
      if (!(v == v)) continue;
      sum += v; sumsq += v*v; cnt++;
    }
    double mean = (cnt > 0) ? (sum / (double)cnt) : NAN;
    double var  = (cnt > 1) ? fmax(0.0, (sumsq - (double)cnt*mean*mean) / (double)(cnt - 1)) : NAN;
    double sd   = (cnt > 1) ? sqrt(var) : NAN;
    if (out_z_mean) *out_z_mean = mean;
    if (out_z_sd)   *out_z_sd   = sd;
  }

  /* STRICT check: after transform, finite entries must be 0/1/2 */
  for (size_t t = 0; t < nm; ++t) {
    double v = Zwork[t];
    if (!(v == v)) continue;
    if (!is_in_set_012(v)) { free(Zwork); return -2002; }
  }

  /* Compute allele frequencies p_j = mean(dosage)/2 (ignore NA) */
  double *p = (double*)malloc((size_t)m * sizeof(double));
  if (!p) { free(Zwork); return -1010; }

  for (int j = 0; j < m; ++j) {
    double s = 0.0;
    int c = 0;
    for (int i = 0; i < nG; ++i) {
      double v = Zwork[(size_t)i*(size_t)m + (size_t)j];
      if (!(v == v)) continue;
      s += v; c++;
    }
    /* If a column is all NA, set p=0.5 (neutral) */
    p[j] = (c > 0) ? ((s / (double)c) / 2.0) : 0.5;
  }

  /* denom = 2 * sum p(1-p) */
  double denom = 0.0;
  for (int j = 0; j < m; ++j) denom += p[j] * (1.0 - p[j]);
  denom = 2.0 * denom;

  if (!(denom > 0.0) || !(denom == denom)) {
    free(p); free(Zwork);
    return -2003;
  }

  /* Center markers: Z <- Z - 2p */
  for (int j = 0; j < m; ++j) {
    const double c = 2.0 * p[j];
    for (int i = 0; i < nG; ++i) {
      double *cell = &Zwork[(size_t)i*(size_t)m + (size_t)j];
      if (*cell == *cell) *cell -= c;
    }
  }

  /* K = (Z Z') / denom */
  /*
  st_dgemm("N","T",
           (gs_UNSIGNED)nG, (gs_UNSIGNED)m,
           (gs_UNSIGNED)nG, (gs_UNSIGNED)m,
           Zwork, Zwork, Krm_out);
  */

  st_dsymm_ZZt((gs_UNSIGNED)nG, (gs_UNSIGNED)m, Zwork, Krm_out);


  const double invd = 1.0 / denom;
  for (size_t t = 0; t < (size_t)nG*(size_t)nG; ++t) Krm_out[t] *= invd;

  free(p);
  free(Zwork);
  return 0;
}



static int str_eq(SEXP a, const char *s) {
  return (TYPEOF(a) == STRSXP && XLENGTH(a) == 1 && strcmp(CHAR(STRING_ELT(a,0)), s) == 0);
}


/* ============================================================
   Fast + fallback cores for intercept-only, single K + I model.
   - st_mixed_Xis1_oneZ: FAST reduced-space REML + BLUP
   - st_mixed_Xis1_twoZ: FALLBACK via st_mxd + LU solves
   Wrapper st_mixed_R handles SEXP only.
   ============================================================ */

/* One-variable debug switch:
   set to 1 to force fallback path always. */

static int st_mixed_force_fallback = 0;  

/* --- small helpers (pure C) --- */

static void set_err(char *errbuf, size_t errlen, const char *msg)
{
  if (!errbuf || errlen == 0) return;
  snprintf(errbuf, errlen, "%s", msg ? msg : "error");
}


static void K_rm_times_vec(int nG, const double *Krm, const double *v, double *out)
{
  const char trans = 'N';
  const int n = nG;
  const double alpha = 1.0;
  const double beta  = 0.0;
  const int inc = 1;

  /* Interpreting row-major as column-major gives K^T; if K is symmetric, OK.
     If K is not symmetric and you want the original row-major K, use trans='T'. */
  F77_CALL(dgemv)(&trans, &n, &n,
                  &alpha,
                  (double*)Krm, &n,
                  (double*)v, &inc,
                  &beta,
                  out, &inc FCONE);
}

/* Build active index list: cnt[g] > 0 */
static int build_active_idx(int nG, const double *cnt, int *idx)
{
  int p = 0;
  for (int g = 0; g < nG; ++g) {
    const double c = cnt[g];
    /* treat only truly-observed genotypes as active */
    if (c >= 1.0) idx[p++] = g;
  }
  return p;
}

/* Build S = I + tau * sqrtD * K * sqrtD on active set, column-major p x p.
   Fills ONLY upper triangle (i<=j). */
static void build_S_sqrtDKsqrtD_cm(int nG,
                                  const double *Krm,   /* row-major nG x nG */
                                  const double *cnt,   /* length nG */
                                  const int *idx,      /* length p */
                                  int p,
                                  double tau,
                                  double *Scm)         /* col-major p x p */
{
  for (int j = 0; j < p; ++j) {
    const int gj = idx[j];
    const double sdj = sqrt(cnt[gj]);
    for (int i = 0; i <= j; ++i) {
      const int gi = idx[i];
      const double sdi = sqrt(cnt[gi]);
      const double kij = Krm[(size_t)gi*(size_t)nG + (size_t)gj]; /* K[gi,gj] */
      const double val = (i == j ? 1.0 : 0.0) + tau * sdi * kij * sdj;
      Scm[(size_t)i + (size_t)p*(size_t)j] = val; /* (i,j) in col-major */
    }
  }
}


/* Solve (I + tau*D*K) x = v.
   Workspaces:
     idx: int[nG] (temporary inside; you may externalize if you want speed)
     Scm: double[nG*nG] used as p*p (col-major) - upper triangle filled
     work: double[nG] used as length p RHS
   Output: out[nG]
*/
/* PATCHED: robust against cnt==0 / tiny / NaN on the "active" set */
/* PATCHED: safe when work==out; robust to cnt<=0 slipping into idx */
static int solve_I_tauDK_inv_times_v(int nG,
                                     const double *Krm,   /* row-major */
                                     const double *cnt,   /* length nG */
                                     double tau,
                                     double *Scm,         /* workspace (>= nG*nG) */
                                     double *work,        /* workspace (>= nG) */
                                     double *out,         /* length nG */
                                     const double *v)     /* length nG */
{
  if (nG <= 0) return 0;
  if (tau == 0.0) { memcpy(out, v, (size_t)nG*sizeof(double)); return 0; }

  int rc = 0;
  int *idx = (int*)malloc((size_t)nG * sizeof(int));
  if (!idx) return -1010;

  const int p = build_active_idx(nG, cnt, idx);

  /* If no active genotypes, (I + tau D K)^{-1} = I on empty set => x=v */
  if (p == 0) {
    memcpy(out, v, (size_t)nG*sizeof(double));
    free(idx);
    return 0;
  }

  /* Build S (upper triangle) for active block */
  build_S_sqrtDKsqrtD_cm(nG, Krm, cnt, idx, p, tau, Scm);

  /* IMPORTANT:
     work is treated as a length-p RHS buffer indexed by i=0..p-1,
     whereas out is indexed by genotype gi (0..nG-1).
     If work==out, writing work[i] would corrupt out[0..p-1] for unrelated genotypes.
     Therefore: use a separate RHS buffer when work and out alias.
  */
  double *rhs = work;
  double *rhs_tmp = NULL;
  if (rhs == out) {
    rhs_tmp = (double*)malloc((size_t)p * sizeof(double));
    if (!rhs_tmp) { rc = -1010; goto cleanup; }
    rhs = rhs_tmp;
  }

  /* RHS = D^{-1/2} v on active set (guard cnt) */
  for (int i = 0; i < p; ++i) {
    const int gi = idx[i];
    const double c = cnt[gi];
    if (!(c > 0.0)) rhs[i] = 0.0;
    else rhs[i] = v[gi] / sqrt(c);
  }

  /* Cholesky solve: S * rhs = rhs */
  {
    char uplo='U';
    int nn=p, lda=p, nrhs_n=1, ldb=p, info=0;
    F77_CALL(dpotrf)(&uplo, &nn, Scm, &lda, &info FCONE);
    if (info != 0) { rc = info; goto cleanup; }
    F77_CALL(dpotrs)(&uplo, &nn, &nrhs_n, Scm, &lda, rhs, &ldb, &info FCONE);
    if (info != 0) { rc = info; goto cleanup; }
  }

  /* Start from identity for all inactive: x = v */
  memcpy(out, v, (size_t)nG*sizeof(double));

  /* x_active = D^{1/2} * rhs (guard cnt) */
  for (int i = 0; i < p; ++i) {
    const int gi = idx[i];
    const double c = cnt[gi];
    out[gi] = (c > 0.0) ? (sqrt(c) * rhs[i]) : 0.0;
  }

cleanup:
  if (rhs_tmp) free(rhs_tmp);
  free(idx);
  return rc;
}


/* ============================================================
   FAST REML evaluation (pure C99, no macros)
   Model: y = 1*beta0 + g + e, g~N(0,sigg*K), e~N(0,sige*I)
   Repeated observations per genotype (incidence Z).
   Uses reduced system with A(tau) = I + tau*D*K, tau = sigg/sige.
   ============================================================ */



typedef struct {
  int nobs;
  int nG;
  const double *Krm;     /* row-major nG x nG */
  const double *cnt;     /* length nG */
  const double *sumy;    /* length nG */
  double sumy_total;
  double yTy;

  /* workspace (owned outside) */
  double *A;    /* nG x nG col-major */
  double *b;    /* nG */
  double *Kt;   /* nG */
  double *s;    /* nG */
} FastREMLWork;

/* Evaluate neg-REML objective (up to constant) at logtau.
   On success returns 0 and sets outputs and *out_val.
   On failure returns nonzero (out_val set huge). */


static int solve_I_tauDK_inv_times_v_cholfact(int nG,
                                              const double *cnt,
                                              const int *idx,
                                              int p,
                                              const double *S_cholU, /* p x p col-major, factor U in upper */
                                              double *rhs_p,         /* length p (work) */
                                              double *t_out,         /* length nG */
                                              const double *v)       /* length nG */
{
  /* default for inactive: x=v */
  memcpy(t_out, v, (size_t)nG * sizeof(double));

  /* rhs = D^{-1/2} v on active */
  for (int i = 0; i < p; ++i) {
    const int gi = idx[i];
    rhs_p[i] = v[gi] / sqrt(cnt[gi]);
  }

  /* Solve S * rhs = rhs using the already computed Cholesky factor.
     NOTE: dpotrs overwrites rhs_p, and does not modify S_cholU. */
  {
    char uplo = 'U';
    int nn = p, lda = p, nrhs = 1, ldb = p, info = 0;
    /* dpotrs may overwrite rhs_p only */
    F77_CALL(dpotrs)(&uplo, &nn, &nrhs, (double*)S_cholU, &lda, rhs_p, &ldb, &info FCONE);
    if (info != 0) return info;
  }

  /* map back: t_active = D^{1/2} rhs */
  for (int i = 0; i < p; ++i) {
    const int gi = idx[i];
    t_out[gi] = sqrt(cnt[gi]) * rhs_p[i];
  }

  return 0;
}

static int eval_reml_oneZ(double logtau,
                          FastREMLWork *w,
                          double *out_beta0,
                          double *out_sigg,
                          double *out_sige,
                          double *out_val)
{
  const int nG   = w->nG;
  const int nobs = w->nobs;
  const double tau = exp(logtau);

  double *Scm  = w->A;  /* workspace: we use only p*p col-major */
  double *rhs  = w->b;  /* workspace: length nG, we use first p as RHS */
  double *Kt   = w->Kt; /* length nG */
  double *s    = w->s;  /* length nG */

  /* active set (same as your working code: heap alloc here) */
  int *idx = (int*)malloc((size_t)nG * sizeof(int));
  if (!idx) { *out_val = 1e300; return -1010; }

  const int p = build_active_idx(nG, w->cnt, idx);
  if (p == 0) {
    free(idx);
    *out_val = 1e300;
    return -97;
  }

  /* Build S(tau) and factor once */
  build_S_sqrtDKsqrtD_cm(nG, w->Krm, w->cnt, idx, p, tau, Scm);

  {
    char uplo = 'U';
    int nn = p, lda = p, info = 0;
    F77_CALL(dpotrf)(&uplo, &nn, Scm, &lda, &info FCONE);
    if (info != 0) {
      free(idx);
      *out_val = 1e300;
      return -7;
    }
  }

  /* logdet(S) from Cholesky diag (same factor) */
  double logdetS = 0.0;
  for (int i = 0; i < p; ++i) {
    const double di = Scm[(size_t)i + (size_t)p*(size_t)i];
    if (!(di > 0.0) || !(di == di)) {
      free(idx);
      *out_val = 1e300;
      return -8;
    }
    logdetS += 2.0 * log(di);
  }

  /* ---- t_sumy = (I+tau DK)^{-1} sumy ---- */
  {
    int info = solve_I_tauDK_inv_times_v_cholfact(nG, w->cnt, idx, p, Scm, rhs, Kt, w->sumy);
    if (info != 0) { free(idx); *out_val = 1e300; return -1; }
  }

  K_rm_times_vec(nG, w->Krm, Kt, s); /* s = K * t_sumy */
  double cnt_K_tsumy = 0.0;
  for (int g = 0; g < nG; ++g) cnt_K_tsumy += w->cnt[g] * s[g];
  const double oneBy = w->sumy_total - tau * cnt_K_tsumy;

  /* ---- t_cnt = (I+tau DK)^{-1} cnt ---- */
  {
    int info = solve_I_tauDK_inv_times_v_cholfact(nG, w->cnt, idx, p, Scm, rhs, Kt, w->cnt);
    if (info != 0) { free(idx); *out_val = 1e300; return -2; }
  }

  K_rm_times_vec(nG, w->Krm, Kt, s); /* s = K * t_cnt */
  double cnt_K_t1 = 0.0;
  for (int g = 0; g < nG; ++g) cnt_K_t1 += w->cnt[g] * s[g];
  const double oneB1 = (double)nobs - tau * cnt_K_t1;

  if (!(oneB1 > 0.0) || !(oneB1 == oneB1)) {
    free(idx);
    *out_val = 1e300;
    return -3;
  }

  const double beta0 = oneBy / oneB1;

  /* rr */
  const double rr = w->yTy - 2.0*beta0*w->sumy_total + beta0*beta0*(double)nobs;
  if (!(rr >= 0.0) || !(rr == rr)) {
    free(idx);
    *out_val = 1e300;
    return -4;
  }

  /* s = sumy - beta0*cnt */
  for (int g = 0; g < nG; ++g) s[g] = w->sumy[g] - beta0 * w->cnt[g];

  /* ---- t_s = (I+tau DK)^{-1} s ---- */
  {
    int info = solve_I_tauDK_inv_times_v_cholfact(nG, w->cnt, idx, p, Scm, rhs, Kt, s);
    if (info != 0) { free(idx); *out_val = 1e300; return -5; }
  }

  /* s' K t_s */
  K_rm_times_vec(nG, w->Krm, Kt, w->b); /* reuse w->b as full-length buffer now */
  double s_K_ts = 0.0;
  for (int g = 0; g < nG; ++g) s_K_ts += s[g] * w->b[g];

  const double q  = rr - tau * s_K_ts;
  const double df = (double)(nobs - 1);
  if (!(q > 0.0) || !(q == q) || !(df > 0.0)) {
    free(idx);
    *out_val = 1e300;
    return -6;
  }

  const double sige = q / df;
  const double sigg = tau * sige;

  *out_val   = 0.5 * ( df * log(sige) + logdetS + log(oneB1) + df );
  *out_beta0 = beta0;
  *out_sigg  = sigg;
  *out_sige  = sige;

  free(idx);
  return 0;
}


/* ============================================================
   FAST PATH: X is intercept-only, one K + I
   Returns 0 on success.
   ============================================================ */
static int st_mixed_Xis1_oneZ(int nobs, int nG,
                              const double *y_o,
                              const int *gid_o,
                              const double *Krm,
                              double tol,
                              double *out_beta0,
                              double *out_sigg,
                              double *out_sige,
                              double *out_ghat,
                              int *out_retval,
                              char *errbuf, size_t errlen)
{

  st_info(1,"Fast path");

  int rc = 0;

  double *cnt = NULL, *sumy = NULL;
  double sumy_total = 0.0, yTy = 0.0;

  double *A = NULL, *b = NULL, *Kt = NULL, *s = NULL, *tmp = NULL;

  double best_logtau = 0.0, best_val = 1e300;
  double best_beta0 = 0.0, best_sigg = 0.0, best_sige = 0.0;

  if (out_retval) *out_retval = 0;

  if (nobs <= 1 || nG <= 0) { set_err(errbuf, errlen, "fast: invalid dims"); return -1; }

  cnt  = (double*)calloc((size_t)nG, sizeof(double));
  sumy = (double*)calloc((size_t)nG, sizeof(double));
  if (!cnt || !sumy) { set_err(errbuf, errlen, "fast: OOM cnt/sumy"); rc = -1010; goto cleanup; }

  for (int i = 0; i < nobs; ++i) {
    int g = gid_o[i];
    const double yi = y_o[i];
    if (g < 0 || g >= nG) { set_err(errbuf, errlen, "fast: gid out of range"); rc = -2; goto cleanup; }
    cnt[g] += 1.0;
    sumy[g] += yi;
    sumy_total += yi;
    yTy += yi * yi;
  }

  A   = (double*)malloc((size_t)nG*(size_t)nG*sizeof(double));
  b   = (double*)malloc((size_t)nG*sizeof(double));
  Kt  = (double*)malloc((size_t)nG*sizeof(double));
  s   = (double*)malloc((size_t)nG*sizeof(double));
  tmp = (double*)malloc((size_t)nG*sizeof(double));
  if (!A || !b || !Kt || !s || !tmp) { set_err(errbuf, errlen, "fast: OOM work"); rc = -1010; goto cleanup; }

  FastREMLWork w;
  w.nobs = nobs; w.nG = nG;
  w.Krm = Krm;
  w.cnt = cnt; w.sumy = sumy;
  w.sumy_total = sumy_total; w.yTy = yTy;
  w.A = A; w.b = b; w.Kt = Kt; w.s = s;

  for (int k = 0; k <= 60; ++k) {
    double logtau = -20.0 + (40.0 * (double)k / 60.0);
    double b0=0.0, sg=0.0, se=0.0, val=0.0;
    int erc = eval_reml_oneZ(logtau, &w, &b0, &sg, &se, &val);
    if (erc == 0 && val < best_val) {
      best_val = val;
      best_logtau = logtau;
      best_beta0 = b0;
      best_sigg = sg;
      best_sige = se;
    }
  }
  if (!(best_val < 1e299)) { set_err(errbuf, errlen, "fast: REML grid failed"); rc = -3; goto cleanup; }

  {
    const double gr = 0.6180339887498949;
    double a = best_logtau - 2.0;
    double bnd = best_logtau + 2.0;
    double c = bnd - gr*(bnd - a);
    double d = a + gr*(bnd - a);

    for (int it = 0; it < 80; ++it) {
      double b0c=0.0, sgc=0.0, sec=0.0, fc=1e300;
      double b0d=0.0, sgd=0.0, sed=0.0, fd=1e300;

      (void)eval_reml_oneZ(c, &w, &b0c, &sgc, &sec, &fc);
      (void)eval_reml_oneZ(d, &w, &b0d, &sgd, &sed, &fd);

      if (fc < fd) {
        bnd = d;
        d = c;
        c = bnd - gr*(bnd - a);
        if (fc < best_val) { best_val = fc; best_logtau = c; best_beta0 = b0c; best_sigg = sgc; best_sige = sec; }
      } else {
        a = c;
        c = d;
        d = a + gr*(bnd - a);
        if (fd < best_val) { best_val = fd; best_logtau = d; best_beta0 = b0d; best_sigg = sgd; best_sige = sed; }
      }

      if (fabs(bnd - a) < tol) break;
    }
  }

  /* BLUP ghat:
     ghat = tau * K * ( s - tau D K t ), t=(I+tau K D)^{-1} s
     with stable solve via S = I + tau D^(1/2) K D^(1/2)
  */
 {
    const double tau = exp(best_logtau);

    /* s = Z'(y - beta0 1) */
    for (int g = 0; g < nG; ++g) s[g] = sumy[g] - best_beta0 * cnt[g];

    /* x = (I + tau D K)^{-1} s  -> store in b */
    if (solve_I_tauDK_inv_times_v(nG, Krm, cnt, tau, A, b, b, s) != 0) {
      set_err(errbuf, errlen, "fast: BLUP solve (I+tau DK) failed");
      rc = -4; goto cleanup;
    }
    /* now b = x */
    
    /* ghat = tau * K * x */
    K_rm_times_vec(nG, Krm, b, out_ghat);
    for (int g = 0; g < nG; ++g) out_ghat[g] *= tau;
 }

  *out_beta0 = best_beta0;
  *out_sigg  = best_sigg;
  *out_sige  = best_sige;

cleanup:
  free(tmp);
  free(s);
  free(Kt);
  free(b);
  free(A);
  free(sumy);
  free(cnt);
  return rc;
}

/* ============================================================
   FALLBACK PATH: st_mxd + LU solve 
   ============================================================ */
static int st_mixed_Xis1_twoZ(int nobs, int nG,
                              const double *y_o,
                              const int *gid_o,
                              const double *Krm,
                              int maxcyc,
                              double tol,
                              double varcmp_inout[2],
                              double *out_beta0,
                              double *out_ghat,   /* length nG */
                              int *out_retval,
                              char *errbuf, size_t errlen)
{
  int rc = 0;

  double *Voo = NULL;
  double *Xo  = NULL;
  double *Sigma = NULL;
  double *Sigma_work = NULL;
  gs_UNSIGNED *piv = NULL;
  double *ones = NULL;
  double *Sy   = NULL;
  double *S1   = NULL;
  double *u    = NULL;

#define FAIL_LOCAL(msg) do { set_err(errbuf, errlen, msg); rc = -1; goto cleanup; } while(0)

  if (out_retval) *out_retval = 0;

  Voo = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  Xo  = (double*)malloc((size_t)nobs*sizeof(double));
  if (!Voo || !Xo) FAIL_LOCAL("fallback: OOM Voo/Xo");

  for (int i = 0; i < nobs; ++i) Xo[i] = 1.0;

  for (int i = 0; i < nobs; ++i) {
    int gi = gid_o[i];
    for (int j = 0; j < nobs; ++j) {
      int gj = gid_o[j];
      Voo[(size_t)i*(size_t)nobs + (size_t)j] = Krm[(size_t)gi*(size_t)nG + (size_t)gj];
    }
  }

  double *VV_arr[1]; VV_arr[0] = Voo;
  int nobs_i = nobs, cX_i = 1, kk_i = 1;
  int retval_mxd = 0;

  st_mxd(&nobs_i, (double*)y_o, &cX_i, Xo, &kk_i,
         NULL, VV_arr, NULL, NULL, varcmp_inout,
         &maxcyc, &tol, &retval_mxd);

  if (out_retval) *out_retval = retval_mxd;

  const double sigg = varcmp_inout[0];
  const double sige = varcmp_inout[1];

  Sigma = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  Sigma_work = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  if (!Sigma || !Sigma_work) FAIL_LOCAL("fallback: OOM Sigma");

  for (int t = 0; t < nobs*nobs; ++t) Sigma[t] = sigg * Voo[t];
  for (int i = 0; i < nobs; ++i) Sigma[(size_t)i*(size_t)nobs + (size_t)i] += sige;

  memcpy(Sigma_work, Sigma, (size_t)nobs*(size_t)nobs*sizeof(double));

  piv = (gs_UNSIGNED*)malloc((size_t)nobs*sizeof(gs_UNSIGNED));
  ones = (double*)malloc((size_t)nobs*sizeof(double));
  Sy   = (double*)malloc((size_t)nobs*sizeof(double));
  S1   = (double*)malloc((size_t)nobs*sizeof(double));
  if (!piv || !ones || !Sy || !S1) FAIL_LOCAL("fallback: OOM piv/ones/Sy/S1");

  for (int i = 0; i < nobs; ++i) ones[i] = 1.0;

  {
    double d = 0.0;
    st_dgetrf(Sigma_work, (gs_UNSIGNED)nobs, piv, &d);
    if ((int)d != 0) FAIL_LOCAL("fallback: Sigma LU failed");
  }

  memcpy(Sy, y_o, (size_t)nobs*sizeof(double));
  st_dgetrs(Sigma_work, (gs_UNSIGNED)nobs, piv, Sy);

  memcpy(S1, ones, (size_t)nobs*sizeof(double));
  st_dgetrs(Sigma_work, (gs_UNSIGNED)nobs, piv, S1);

  double num = 0.0, den = 0.0;
  for (int i = 0; i < nobs; ++i) { num += Sy[i]; den += S1[i]; }
  if (!(den != 0.0) || !(den == den)) FAIL_LOCAL("fallback: beta0 denom invalid");

  *out_beta0 = num / den;

  u = (double*)malloc((size_t)nobs*sizeof(double));
  if (!u) FAIL_LOCAL("fallback: OOM u");
  for (int i = 0; i < nobs; ++i) u[i] = y_o[i] - *out_beta0;
  st_dgetrs(Sigma_work, (gs_UNSIGNED)nobs, piv, u);

  for (int g = 0; g < nG; ++g) out_ghat[g] = 0.0;
  for (int i = 0; i < nobs; ++i) {
    int gi = gid_o[i];
    double ui = u[i];
    for (int g = 0; g < nG; ++g) {
      out_ghat[g] += sigg * Krm[(size_t)g*(size_t)nG + (size_t)gi] * ui;
    }
  }

cleanup:
  free(u);
  free(S1);
  free(Sy);
  free(ones);
  free(piv);
  free(Sigma_work);
  free(Sigma);
  free(Xo);
  free(Voo);
  return rc;

#undef FAIL_LOCAL
}


/* ============================================================
   Wrapper: st_mixed_R (SEXP only here)
   ============================================================ */


/* -------- shared core: fixed VC single-kernel -------- */
static int st_mixed_single_core_fixedVC_impl(
  int nobs, int nG,
  const double *y_o,
  const int *gid_o,
  const double *Krm,         /* row-major nG*nG */
  double sigg, double sige,
  double ridge,              /* can be 0 */
  double *out_beta0,
  double *out_ghat,          /* length nG */
  double *out_yhat,          /* length nG */
  char *errbuf, size_t errlen
) {
  if (nobs <= 0 || nG <= 0) {
    snprintf(errbuf, errlen, "single_core_fixedVC: invalid dims");
    return -1;
  }
  if (!(sigg >= 0.0) || !(sige > 0.0) || !(sigg == sigg) || !(sige == sige)) {
    snprintf(errbuf, errlen, "single_core_fixedVC: invalid varcmp");
    return -2;
  }

  double *cnt  = (double*)calloc((size_t)nG, sizeof(double));
  double *sumy = (double*)calloc((size_t)nG, sizeof(double));
  if (!cnt || !sumy) {
    free(cnt); free(sumy);
    snprintf(errbuf, errlen, "single_core_fixedVC: OOM cnt/sumy");
    return -1010;
  }

  /* collapse to genotype sums */
  for (int i = 0; i < nobs; ++i) {
    int g = gid_o[i];
    if (g < 0 || g >= nG) {
      free(cnt); free(sumy);
      snprintf(errbuf, errlen, "single_core_fixedVC: gid out of range");
      return -3;
    }
    cnt[g]  += 1.0;
    sumy[g] += y_o[i];
  }

  int *idx = (int*)malloc((size_t)nG * sizeof(int));
  if (!idx) {
    free(cnt); free(sumy);
    snprintf(errbuf, errlen, "single_core_fixedVC: OOM idx");
    return -1010;
  }
  const int p = build_active_idx(nG, cnt, idx);
  if (p <= 0) {
    /* no active genotypes? shouldn't happen */
    for (int g = 0; g < nG; ++g) out_ghat[g] = 0.0;
    *out_beta0 = NA_REAL;
    for (int g = 0; g < nG; ++g) out_yhat[g] = NA_REAL;
    free(idx); free(cnt); free(sumy);
    snprintf(errbuf, errlen, "single_core_fixedVC: no active genotypes");
    return -4;
  }

  /* tau = sigg/sige. If sigg==0 => ghat=0 and beta0 = mean(y) */
  if (sigg == 0.0) {
    double num = 0.0, den = 0.0;
    for (int ii = 0; ii < p; ++ii) {
      int g = idx[ii];
      num += sumy[g];
      den += cnt[g];
    }
    double beta0 = (den > 0.0) ? (num / den) : NA_REAL;
    *out_beta0 = beta0;
    for (int g = 0; g < nG; ++g) out_ghat[g] = 0.0;
    for (int g = 0; g < nG; ++g) out_yhat[g] = beta0;
    free(idx); free(cnt); free(sumy);
    return 0;
  }

  const double tau = sigg / sige;

  /* Workspaces sized by nG (your solver uses active subset internally) */
  double *Scm  = (double*)malloc((size_t)p * (size_t)p * sizeof(double));
  double *work = (double*)malloc((size_t)p * sizeof(double));
  double *tmpG = (double*)malloc((size_t)nG * sizeof(double));
  double *Wy   = (double*)malloc((size_t)nG * sizeof(double));
  double *W1   = (double*)malloc((size_t)nG * sizeof(double));
  double *r    = (double*)malloc((size_t)nG * sizeof(double));
  double *alpha= (double*)malloc((size_t)nG * sizeof(double));
  if (!Scm || !work || !tmpG || !Wy || !W1 || !r || !alpha) {
    free(Scm); free(work); free(tmpG); free(Wy); free(W1); free(r); free(alpha);
    free(idx); free(cnt); free(sumy);
    snprintf(errbuf, errlen, "single_core_fixedVC: OOM work");
    return -1010;
  }

  /* Compute Wy = (I + tau D K)^{-1} sumy  and W1 = (I + tau D K)^{-1} cnt */
  {
    int rc;

    rc = solve_I_tauDK_inv_times_v(nG, Krm, cnt, tau, Scm, work, Wy, sumy);
    if (rc != 0) {
      snprintf(errbuf, errlen, "single_core_fixedVC: solve Wy failed (rc=%d)", rc);
      goto fail;
    }
    rc = solve_I_tauDK_inv_times_v(nG, Krm, cnt, tau, Scm, work, W1, cnt);
    if (rc != 0) {
      snprintf(errbuf, errlen, "single_core_fixedVC: solve W1 failed (rc=%d)", rc);
      goto fail;
    }
  }

  /* beta0 = (1' Wy)/(1' W1) but only active contribute (inactive are zeros anyway) */
  double num = 0.0, den = 0.0;
  for (int ii = 0; ii < p; ++ii) {
    int g = idx[ii];
    num += Wy[g];
    den += W1[g];
  }
  if (!(den > 0.0) || !(den == den)) {
    snprintf(errbuf, errlen, "single_core_fixedVC: invalid beta0 denom");
    goto fail;
  }
  const double beta0 = num / den;
  *out_beta0 = beta0;

  /* r = sumy - beta0*cnt */
  for (int g = 0; g < nG; ++g) r[g] = sumy[g] - beta0 * cnt[g];

  /* alpha = (I + tau D K)^{-1} r */
  {
    int rc = solve_I_tauDK_inv_times_v(nG, Krm, cnt, tau, Scm, work, alpha, r);
    if (rc != 0) {
      snprintf(errbuf, errlen, "single_core_fixedVC: solve alpha failed (rc=%d)", rc);
      goto fail;
    }
  }

  /* ghat = tau * K * alpha  (same as your fast path BLUP form) */
  K_rm_times_vec(nG, Krm, alpha, out_ghat);
  for (int g = 0; g < nG; ++g) out_ghat[g] *= tau;

  /* yhat per genotype */
  for (int g = 0; g < nG; ++g) out_yhat[g] = beta0 + out_ghat[g];

  free(Scm); free(work); free(tmpG); free(Wy); free(W1); free(r); free(alpha);
  free(idx); free(cnt); free(sumy);
  return 0;

fail:
  free(Scm); free(work); free(tmpG); free(Wy); free(W1); free(r); free(alpha);
  free(idx); free(cnt); free(sumy);
  return -5;
}


static SEXP st_mixed_single_impl(SEXP y_, SEXP id_, SEXP Z_, SEXP Z_ids_,
                                SEXP method_, SEXP coding_, SEXP varcmp_,
                                SEXP maxcyc_, SEXP tol_)
{
  int nprot = 0;
  char errbuf[512];
  errbuf[0] = '\0';

  /* malloc'd pointers (init NULL for safe cleanup) */
  double *y_kept = NULL;
  int    *gid_kept = NULL;
  double *Zrm = NULL;
  double *Krm = NULL;
  int    *obs = NULL;
  double *y_o = NULL;
  int    *gid_o = NULL;
  double *ghat_c = NULL;

  /* SEXPs */
  SEXP idc = R_NilValue;
  SEXP miss_names = R_NilValue;
  SEXP missing_ids = R_NilValue;

  SEXP beta0 = R_NilValue;
  SEXP varcmp_out = R_NilValue;
  SEXP ghat = R_NilValue;
  SEXP yhat = R_NilValue;
  SEXP diag = R_NilValue;
  SEXP dn = R_NilValue;
  SEXP ans = R_NilValue;
  SEXP an = R_NilValue;

  int ny = 0, nG = 0, m = 0;
  int nkeep = 0, nobs = 0, nmis = 0;

  double beta0_val = 0.0;
  double sigg = NAN, sige = NAN;
  double z_mean = 0.0, z_sd = 0.0;
  ZCoding coding_used = CODING_OTHER;
  ZTransform tr = TRANS_NONE;

  int retval_mxd = 0;

  /* -------------------- argument checks -------------------- */
  if (!isReal(y_) || !isVector(y_)) {
    st_info(-2, "st_mixed_single_impl: y must be numeric vector");
    return R_NilValue;
  }
  if (!isVector(id_)) {
    st_info(-2, "st_mixed_single_impl: id must be a vector (character/factor)");
    return R_NilValue;
  }
  if (!isReal(Z_) || !isMatrix(Z_)) {
    st_info(-2, "st_mixed_single_impl: Z must be a numeric matrix");
    return R_NilValue;
  }
  if (!isString(Z_ids_) || XLENGTH(Z_ids_) < 1) {
    st_info(-2, "st_mixed_single_impl: Z_ids must be a character vector");
    return R_NilValue;
  }
  if (!isString(method_) || XLENGTH(method_) != 1) {
    st_info(-2, "st_mixed_single_impl: method must be a single string");
    return R_NilValue;
  }
  if (!isString(coding_) || XLENGTH(coding_) != 1) {
    st_info(-2, "st_mixed_single_impl: coding must be a single string");
    return R_NilValue;
  }

  ny = (int)XLENGTH(y_);
  {
    SEXP Zdim = getAttrib(Z_, R_DimSymbol);
    nG = INTEGER(Zdim)[0];
    m  = INTEGER(Zdim)[1];
  }

  if ((int)XLENGTH(Z_ids_) != nG) {
    st_info(-2, "st_mixed_single_impl: Z_ids length must equal nrow(Z)");
    return R_NilValue;
  }

  /* enforce unique Z_ids (O(n^2), as in your current code) */
  for (int i = 0; i < nG; ++i) {
    const char *si = CHAR(STRING_ELT(Z_ids_, i));
    for (int j = i + 1; j < nG; ++j) {
      if (strcmp(si, CHAR(STRING_ELT(Z_ids_, j))) == 0) {
        st_info(-2, "st_mixed_single_impl: Z_ids must be unique");
        return R_NilValue;
      }
    }
  }

  if (TYPEOF(id_) == STRSXP) idc = id_;
  else { idc = PROTECT(coerceVector(id_, STRSXP)); nprot++; }

  /* map id -> gid (index in Z_ids_) and keep only those found */
  int *gid  = (int*)R_alloc((size_t)ny, sizeof(int));
  int *keep = (int*)R_alloc((size_t)ny, sizeof(int));
  miss_names = PROTECT(allocVector(STRSXP, ny)); nprot++;
  int nmiss = 0;
  nkeep = 0;

  for (int i = 0; i < ny; ++i) {
    const char *idname = CHAR(STRING_ELT(idc, i));
    int found = -1;
    for (int j = 0; j < nG; ++j) {
      if (strcmp(idname, CHAR(STRING_ELT(Z_ids_, j))) == 0) { found = j; break; }
    }
    if (found >= 0) {
      gid[i] = found;
      keep[nkeep++] = i;
    } else {
      gid[i] = -1;
      SET_STRING_ELT(miss_names, nmiss++, STRING_ELT(idc, i));
    }
  }

  /* unique missing ids */
  int *miss_keep = (int*)R_alloc((size_t)nmiss, sizeof(int));
  int nmiss_u = 0;
  for (int i = 0; i < nmiss; ++i) {
    const char *si = CHAR(STRING_ELT(miss_names, i));
    int dup = 0;
    for (int j = 0; j < nmiss_u; ++j) {
      if (strcmp(si, CHAR(STRING_ELT(miss_names, miss_keep[j]))) == 0) { dup = 1; break; }
    }
    if (!dup) miss_keep[nmiss_u++] = i;
  }

  missing_ids = PROTECT(allocVector(STRSXP, nmiss_u)); nprot++;
  for (int i = 0; i < nmiss_u; ++i) {
    SET_STRING_ELT(missing_ids, i, STRING_ELT(miss_names, miss_keep[i]));
  }

  /* keep y/id only for matched ids */
  y_kept   = (double*)malloc((size_t)nkeep * sizeof(double));
  gid_kept = (int*)malloc((size_t)nkeep * sizeof(int));
  if (!y_kept || !gid_kept) {
    st_info(-2, "st_mixed_single_impl: Out of memory (y_kept/gid_kept)");
    goto cleanup_err;
  }
  for (int i = 0; i < nkeep; ++i) {
    int idx = keep[i];
    y_kept[i]   = REAL(y_)[idx];
    gid_kept[i] = gid[idx];
  }

  /* method -> enum */
  KMethod km = K_STD;
  if (str_eq(method_, "std")) km = K_STD;
  else if (str_eq(method_, "vanraden")) km = K_VANRADEN;
  else if (str_eq(method_, "sommer")) km = K_SOMMER;
  else {
    st_info(-2, "st_mixed_single_impl: method must be one of: 'std','vanraden','sommer'");
    goto cleanup_err;
  }

  /* allocate Zrm + Krm */
  Zrm = (double*)malloc((size_t)nG * (size_t)m  * sizeof(double));
  Krm = (double*)malloc((size_t)nG * (size_t)nG * sizeof(double));
  if (!Zrm || !Krm) {
    st_info(-2, "st_mixed_single_impl: Out of memory (Zrm/Krm)");
    goto cleanup_err;
  }

  /* Z: column-major (R) -> row-major (C) */
  Z_cm_to_rm(REAL(Z_), nG, m, Zrm);

  /* build K from Z */
  {
    int rcK = build_K_from_Z(Zrm, nG, m, km, coding_, Krm,
                             &z_mean, &z_sd, &coding_used, &tr);
    if (rcK == -2004) {
      st_info(-2, "st_mixed_single_impl: coding='detect' is ambiguous for this Z; set coding explicitly.");
      goto cleanup_err;
    }
    if (rcK != 0) {
      snprintf(errbuf, sizeof(errbuf), "st_mixed_single_impl: build_K_from_Z failed (code %d).", rcK);
      st_info(-2, errbuf);
      goto cleanup_err;
    }
  }

  /* build observed subset */
  obs = (int*)malloc((size_t)nkeep * sizeof(int));
  if (!obs) {
    st_info(-2, "st_mixed_single_impl: Out of memory (obs)");
    goto cleanup_err;
  }

  nobs = 0;
  for (int i = 0; i < nkeep; ++i) if (!is_missing(y_kept[i])) obs[nobs++] = i;
  nmis = nkeep - nobs;
  if (nobs <= 0) {
    st_info(-2, "st_mixed_single_impl: No non-missing y values after matching ids to Z.");
    goto cleanup_err;
  }

  y_o = (double*)malloc((size_t)nobs * sizeof(double));
  gid_o = (int*)malloc((size_t)nobs * sizeof(int));
  if (!y_o || !gid_o) {
    st_info(-2, "st_mixed_single_impl: Out of memory (y_o/gid_o)");
    goto cleanup_err;
  }

  for (int i = 0; i < nobs; ++i) {
    y_o[i]   = y_kept[ obs[i] ];
    gid_o[i] = gid_kept[ obs[i] ];
  }

  ghat_c = (double*)malloc((size_t)nG * sizeof(double));
  if (!ghat_c) {
    st_info(-2, "st_mixed_single_impl: Out of memory (ghat_c)");
    goto cleanup_err;
  }

  /* ============================================================
     FIXED-VC PATH (no duplication): if varcmp_ is provided,
     skip VC estimation and call the shared core.
     ============================================================ */
  if (!isNull(varcmp_)) {
    if (!isReal(varcmp_) || XLENGTH(varcmp_) != 2) {
      st_info(-2, "st_mixed_single_impl: varcmp must be NULL or numeric length 2");
      goto cleanup_err;
    }
    double sigg_in = REAL(varcmp_)[0];
    double sige_in = REAL(varcmp_)[1];

    if (!(sigg_in >= 0.0) || !(sige_in > 0.0) || !(sigg_in == sigg_in) || !(sige_in == sige_in)) {
      st_info(-2, "st_mixed_single_impl: varcmp invalid (need sigg>=0 and sige>0)");
      goto cleanup_err;
    }

    /* core returns beta0 and ghat (and can also compute yhat) */
    double *yhat_tmp = (double*)malloc((size_t)nG * sizeof(double));
    if (!yhat_tmp) {
      st_info(-2, "st_mixed_single_impl: Out of memory (yhat_tmp)");
      goto cleanup_err;
    }

    int rc_core = st_mixed_single_core_fixedVC_impl(
      nobs, nG,
      y_o, gid_o,
      Krm,
      sigg_in, sige_in,
      0.0,              /* ridge (single path has none) */
      &beta0_val,
      ghat_c,
      yhat_tmp,
      errbuf, sizeof(errbuf)
    );

    free(yhat_tmp);

    if (rc_core != 0) {
      if (errbuf[0] == '\0') snprintf(errbuf, sizeof(errbuf), "st_mixed_single_impl: fixed-VC core failed (rc=%d)", rc_core);
      st_info(-2, errbuf);
      goto cleanup_err;
    }

    sigg = sigg_in;
    sige = sige_in;

    goto outputs;
  }

  /* ============================================================
     ESTIMATE VC PATH: use existing fast/fallback (no duplication)
     ============================================================ */
  if (!st_mixed_force_fallback) {
    int rc_fast = st_mixed_Xis1_oneZ(
      nobs, nG,
      y_o, gid_o,
      Krm,
      asReal(tol_),
      &beta0_val, &sigg, &sige,
      ghat_c,
      NULL, errbuf, sizeof(errbuf)
    );

    if (rc_fast == 0 &&
        (sigg == sigg) && (sige == sige) &&
        (sigg > 0.0) && (sige > 0.0)) {
      /* ok */
    } else {
      /* fall through */
      st_mixed_force_fallback = 0;
      sigg = NAN; sige = NAN;
    }
  }

  if (st_mixed_force_fallback ||
      !(sigg > 0.0) || !(sige > 0.0) ||
      !(sigg == sigg) || !(sige == sige))
  {
    double varcmp2[2];
    if (isNull(varcmp_)) { varcmp2[0] = 0.0; varcmp2[1] = 0.0; }
    else {
      if (!isReal(varcmp_) || XLENGTH(varcmp_) != 2) {
        st_info(-2, "st_mixed_single_impl: varcmp must be NULL or numeric length 2");
        goto cleanup_err;
      }
      varcmp2[0] = REAL(varcmp_)[0];
      varcmp2[1] = REAL(varcmp_)[1];
    }

    int rc_fb = st_mixed_Xis1_twoZ(
      nobs, nG,
      y_o, gid_o,
      Krm,
      asInteger(maxcyc_),
      asReal(tol_),
      varcmp2,
      &beta0_val,
      ghat_c,
      &retval_mxd,
      errbuf, sizeof(errbuf)
    );

    if (rc_fb != 0) {
      if (errbuf[0] == '\0') snprintf(errbuf, sizeof(errbuf), "st_mixed_single_impl: fallback failed (rc=%d)", rc_fb);
      st_info(-2, errbuf);
      goto cleanup_err;
    }

    sigg = varcmp2[0];
    sige = varcmp2[1];
  }

outputs:
  /* outputs */
  beta0 = PROTECT(allocVector(REALSXP, 1)); nprot++;
  varcmp_out = PROTECT(allocVector(REALSXP, 2)); nprot++;
  ghat = PROTECT(allocVector(REALSXP, nG)); nprot++;
  yhat = PROTECT(allocVector(REALSXP, nG)); nprot++;

  REAL(beta0)[0] = beta0_val;
  REAL(varcmp_out)[0] = sigg;
  REAL(varcmp_out)[1] = sige;

  for (int g = 0; g < nG; ++g) REAL(ghat)[g] = ghat_c[g];
  for (int g = 0; g < nG; ++g) REAL(yhat)[g] = beta0_val + REAL(ghat)[g];

  diag = PROTECT(allocVector(VECSXP, 9)); nprot++;
  dn   = PROTECT(allocVector(STRSXP, 9)); nprot++;
  SET_STRING_ELT(dn, 0, mkChar("n_total"));
  SET_STRING_ELT(dn, 1, mkChar("n_used"));
  SET_STRING_ELT(dn, 2, mkChar("n_train"));
  SET_STRING_ELT(dn, 3, mkChar("n_pred"));
  SET_STRING_ELT(dn, 4, mkChar("nG"));
  SET_STRING_ELT(dn, 5, mkChar("m"));
  SET_STRING_ELT(dn, 6, mkChar("Z_mean_sd"));
  SET_STRING_ELT(dn, 7, mkChar("Z_coding_used"));
  SET_STRING_ELT(dn, 8, mkChar("Z_transform_applied"));
  setAttrib(diag, R_NamesSymbol, dn);

  SET_VECTOR_ELT(diag, 0, ScalarInteger(ny));
  SET_VECTOR_ELT(diag, 1, ScalarInteger(nkeep));
  SET_VECTOR_ELT(diag, 2, ScalarInteger(nobs));
  SET_VECTOR_ELT(diag, 3, ScalarInteger(nmis));
  SET_VECTOR_ELT(diag, 4, ScalarInteger(nG));
  SET_VECTOR_ELT(diag, 5, ScalarInteger(m));
  {
    SEXP ms = PROTECT(allocVector(REALSXP, 2)); nprot++;
    REAL(ms)[0] = z_mean;
    REAL(ms)[1] = z_sd;
    SET_VECTOR_ELT(diag, 6, ms);
  }
  SET_VECTOR_ELT(diag, 7, mkString(coding_name(coding_used)));
  SET_VECTOR_ELT(diag, 8, mkString(transform_name(tr)));

  ans = PROTECT(allocVector(VECSXP, 7)); nprot++;
  an  = PROTECT(allocVector(STRSXP, 7)); nprot++;
  SET_STRING_ELT(an, 0, mkChar("beta0"));
  SET_STRING_ELT(an, 1, mkChar("varcmp"));
  SET_STRING_ELT(an, 2, mkChar("ghat"));
  SET_STRING_ELT(an, 3, mkChar("yhat"));
  SET_STRING_ELT(an, 4, mkChar("missing_ids"));
  SET_STRING_ELT(an, 5, mkChar("diag"));
  SET_STRING_ELT(an, 6, mkChar("retval"));
  setAttrib(ans, R_NamesSymbol, an);

  SET_VECTOR_ELT(ans, 0, beta0);
  SET_VECTOR_ELT(ans, 1, varcmp_out);
  SET_VECTOR_ELT(ans, 2, ghat);
  SET_VECTOR_ELT(ans, 3, yhat);
  SET_VECTOR_ELT(ans, 4, missing_ids);
  SET_VECTOR_ELT(ans, 5, diag);
  SET_VECTOR_ELT(ans, 6, ScalarInteger(retval_mxd));

  goto cleanup;

cleanup_err:
  ans = R_NilValue;
  goto cleanup;

cleanup:
  free(ghat_c);
  free(gid_o);
  free(y_o);
  free(obs);
  free(Krm);
  free(Zrm);
  free(gid_kept);
  free(y_kept);

  UNPROTECT(nprot);
  return ans;
}


//##############################################################################
//# Hybrids
//##############################################################################


/* linear search id in STRSXP */
static int find_id_index(SEXP ids, const char *key) {
  int n = (int)XLENGTH(ids);
  for (int i = 0; i < n; ++i) {
    if (strcmp(key, CHAR(STRING_ELT(ids, i))) == 0) return i;
  }
  return -1;
}

/* convert row-major A_rm (n×n) to column-major A_cm (n×n) */
static void rm_to_cm_sym(int n, const double *A_rm, double *A_cm) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      A_cm[i + (size_t)n * (size_t)j] = A_rm[(size_t)i * (size_t)n + (size_t)j];
    }
  }
}

static SEXP st_mixed_hybrid_impl(SEXP y_, SEXP hybrid_id_,
                                SEXP pool1_, SEXP pool2_,
                                SEXP Z1_, SEXP Z1_ids_,
                                SEXP Z2_, SEXP Z2_ids_,
                                SEXP include_sca_,
                                SEXP method_, SEXP coding_, SEXP varcmp_,
                                SEXP maxcyc_, SEXP tol_, SEXP ridge_)
{
  int nprot = 0;
  char errbuf[512]; errbuf[0] = '\0';

  /* malloc pointers init NULL for safe cleanup */
  int *train_pos = NULL;

  double *Z1rm=NULL, *Z2rm=NULL;
  double *Kp1_rm=NULL, *Kp2_rm=NULL;
  double *Kp1_cm=NULL, *Kp2_cm=NULL;

  /* training kernels for st_mxd (row-major) */
  double *K1_tr_rm=NULL, *K2_tr_rm=NULL, *K3_tr_rm=NULL;
  double **VV=NULL;

  double *ytr=NULL;
  double *X=NULL;
  double *varcmp=NULL;

  /* prediction: train×use cross-covariances (column-major, columns are use hybrids) */
  double *K1_truse_cm=NULL, *K2_truse_cm=NULL, *K3_truse_cm=NULL;

  /* Sigma training (column-major) */
  double *K1_tr_cm=NULL, *K2_tr_cm=NULL, *K3_tr_cm=NULL;
  double *Sigma=NULL;

  /* RHS for solves */
  double *B=NULL;      /* nobs×2: col0=ytr, col1=ones */
  double *r=NULL;      /* residual */
  double *alpha=NULL;  /* Sigma^{-1} r */

  /* outputs */
  SEXP yhat_obs = R_NilValue;
  SEXP beta0_out = R_NilValue;
  SEXP varcmp_out = R_NilValue;
  SEXP diag = R_NilValue, dn = R_NilValue;
  SEXP ans = R_NilValue, an = R_NilValue;

  int retval = 0;

#define FAIL_GOTO(MSG_LITERAL) do { \
  st_info(-2, (MSG_LITERAL)); \
  goto cleanup_err; \
} while(0)

#define FAILF_GOTO(FMT, ...) do { \
  int _n = snprintf(errbuf, sizeof(errbuf), (FMT), __VA_ARGS__); \
  if (_n < 0) snprintf(errbuf, sizeof(errbuf), "%s", "error"); \
  else if (_n >= (int)sizeof(errbuf)) { \
    size_t L = sizeof(errbuf); \
    if (L >= 4) { errbuf[L-4]='.'; errbuf[L-3]='.'; errbuf[L-2]='.'; errbuf[L-1]='\0'; } \
    else { errbuf[L-1]='\0'; } \
  } \
  st_info(-2, errbuf); \
  goto cleanup_err; \
} while(0)

  /* -------------------- argument checks -------------------- */
  if (!isReal(y_) || !isVector(y_)) FAIL_GOTO("y must be numeric vector");
  const int ny = (int)XLENGTH(y_);

  if (!isString(hybrid_id_) || (int)XLENGTH(hybrid_id_) != ny) FAIL_GOTO("hybrid_id must be character length(y)");
  if (!isString(pool1_)     || (int)XLENGTH(pool1_)     != ny) FAIL_GOTO("pool1 must be character length(y)");
  if (!isString(pool2_)     || (int)XLENGTH(pool2_)     != ny) FAIL_GOTO("pool2 must be character length(y)");

  if (!isReal(Z1_) || !isMatrix(Z1_)) FAIL_GOTO("Z1 must be numeric matrix");
  if (!isReal(Z2_) || !isMatrix(Z2_)) FAIL_GOTO("Z2 must be numeric matrix");

  if (!isString(Z1_ids_) || XLENGTH(Z1_ids_) < 1) FAIL_GOTO("Z1_ids must be character");
  if (!isString(Z2_ids_) || XLENGTH(Z2_ids_) < 1) FAIL_GOTO("Z2_ids must be character");

  if (!isString(method_) || XLENGTH(method_) != 1) FAIL_GOTO("method must be single string");
  if (!isString(coding_) || XLENGTH(coding_) != 1) FAIL_GOTO("coding must be single string");

  const int include_sca = asLogical(include_sca_) ? 1 : 0;
  const double ridge = asReal(ridge_);
  if (!(ridge >= 0.0) || !(ridge == ridge)) FAIL_GOTO("ridge must be >= 0");

  /* duplicated hybrids (O(n^2) as requested earlier; OK at ~1k) */
  for (int i = 0; i < ny; ++i) {
    const char *hi = CHAR(STRING_ELT(hybrid_id_, i));
    for (int j = i+1; j < ny; ++j) {
      if (strcmp(hi, CHAR(STRING_ELT(hybrid_id_, j))) == 0) {
        FAILF_GOTO("Duplicated hybrids detected (hybrid_id must be unique), e.g. '%s'", hi);
      }
    }
  }

  /* dims */
  SEXP d1 = getAttrib(Z1_, R_DimSymbol);
  SEXP d2 = getAttrib(Z2_, R_DimSymbol);
  const int nP1 = INTEGER(d1)[0], m1 = INTEGER(d1)[1];
  const int nP2 = INTEGER(d2)[0], m2 = INTEGER(d2)[1];
  if ((int)XLENGTH(Z1_ids_) != nP1) FAIL_GOTO("length(Z1_ids) must equal nrow(Z1)");
  if ((int)XLENGTH(Z2_ids_) != nP2) FAIL_GOTO("length(Z2_ids) must equal nrow(Z2)");

  /* method -> enum */
  KMethod km = K_STD;
  if (str_eq(method_, "std")) km = K_STD;
  else if (str_eq(method_, "vanraden")) km = K_VANRADEN;
  else if (str_eq(method_, "sommer")) km = K_SOMMER;
  else FAIL_GOTO("method must be one of: 'std','vanraden','sommer'");

  /* map each obs to pool indices; keep those found in both */
  int *p1_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *p2_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *use    = (int*)R_alloc((size_t)ny, sizeof(int));
  int nuse = 0;

  for (int i = 0; i < ny; ++i) {
    const char *a = CHAR(STRING_ELT(pool1_, i));
    const char *b = CHAR(STRING_ELT(pool2_, i));
    int ia = find_id_index(Z1_ids_, a);
    int ib = find_id_index(Z2_ids_, b);
    p1_idx[i] = ia;
    p2_idx[i] = ib;
    if (ia >= 0 && ib >= 0) use[nuse++] = i;
  }
  if (nuse <= 0) FAIL_GOTO("No rows matched to both pool1 and pool2 IDs.");

  /* train positions within use-space */
  train_pos = (int*)malloc((size_t)nuse * sizeof(int));
  if (!train_pos) FAIL_GOTO("Out of memory");
  int nobs = 0, npred = 0;
  for (int t = 0; t < nuse; ++t) {
    int i = use[t];
    if (!is_missing(REAL(y_)[i])) train_pos[nobs++] = t;
    else npred++;
  }
  if (nobs <= 0) FAIL_GOTO("No non-missing y values after matching parents.");

  /* (5) Precompute use-space parent indices and train-space parent indices */
  int *p1u = (int*)R_alloc((size_t)nuse, sizeof(int));
  int *p2u = (int*)R_alloc((size_t)nuse, sizeof(int));
  for (int t = 0; t < nuse; ++t) {
    int obs = use[t];
    p1u[t] = p1_idx[obs];
    p2u[t] = p2_idx[obs];
  }
  int *p1t = (int*)R_alloc((size_t)nobs, sizeof(int));
  int *p2t = (int*)R_alloc((size_t)nobs, sizeof(int));
  for (int a = 0; a < nobs; ++a) {
    int tpos = train_pos[a];
    p1t[a] = p1u[tpos];
    p2t[a] = p2u[tpos];
  }

  /* ---- build parental kernels with existing builder: row-major -> column-major ---- */
  Z1rm   = (double*)malloc((size_t)nP1*(size_t)m1*sizeof(double));
  Z2rm   = (double*)malloc((size_t)nP2*(size_t)m2*sizeof(double));
  Kp1_rm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  Kp2_rm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  Kp1_cm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  Kp2_cm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  if (!Z1rm || !Z2rm || !Kp1_rm || !Kp2_rm || !Kp1_cm || !Kp2_cm) FAIL_GOTO("Out of memory");

  Z_cm_to_rm(REAL(Z1_), nP1, m1, Z1rm);
  Z_cm_to_rm(REAL(Z2_), nP2, m2, Z2rm);

  double zmean1=0, zsd1=0, zmean2=0, zsd2=0;
  ZCoding cu1 = CODING_OTHER, cu2 = CODING_OTHER;
  ZTransform tr1 = TRANS_NONE, tr2 = TRANS_NONE;

  int rc1 = build_K_from_Z(Z1rm, nP1, m1, km, coding_, Kp1_rm, &zmean1, &zsd1, &cu1, &tr1);
  if (rc1 == -2004) FAIL_GOTO("coding='detect' ambiguous for Z1; set coding explicitly.");
  if (rc1 != 0) FAILF_GOTO("build_K_from_Z(Z1) failed (code %d).", rc1);

  int rc2 = build_K_from_Z(Z2rm, nP2, m2, km, coding_, Kp2_rm, &zmean2, &zsd2, &cu2, &tr2);
  if (rc2 == -2004) FAIL_GOTO("coding='detect' ambiguous for Z2; set coding explicitly.");
  if (rc2 != 0) FAILF_GOTO("build_K_from_Z(Z2) failed (code %d).", rc2);

  rm_to_cm_sym(nP1, Kp1_rm, Kp1_cm);
  rm_to_cm_sym(nP2, Kp2_rm, Kp2_cm);

  /* ---- (3) Build only needed blocks (no full K_use) ---- */

  const int kk = include_sca ? 3 : 2;

  /* training blocks (row-major) for st_mxd */
  K1_tr_rm = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  K2_tr_rm = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  if (include_sca) K3_tr_rm = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  if (!K1_tr_rm || !K2_tr_rm || (include_sca && !K3_tr_rm)) FAIL_GOTO("Out of memory");

  /* fill training kernels directly from parental GRMs */
  for (int i = 0; i < nobs; ++i) {
    const int ia = p1t[i];
    const int ja = p2t[i];
    for (int j = 0; j < nobs; ++j) {
      const int ib = p1t[j];
      const int jb = p2t[j];
      double v1 = Kp1_cm[ ia + (size_t)nP1*(size_t)ib ];
      double v2 = Kp2_cm[ ja + (size_t)nP2*(size_t)jb ];
      K1_tr_rm[(size_t)i*(size_t)nobs + (size_t)j] = v1; /* row-major */
      K2_tr_rm[(size_t)i*(size_t)nobs + (size_t)j] = v2;
      if (include_sca) K3_tr_rm[(size_t)i*(size_t)nobs + (size_t)j] = v1 * v2;
    }
  }

  /* y_train */
  ytr = (double*)malloc((size_t)nobs*sizeof(double));
  if (!ytr) FAIL_GOTO("Out of memory");
  for (int a = 0; a < nobs; ++a) ytr[a] = REAL(y_)[ use[ train_pos[a] ] ];

  /* X intercept */
  int cX = 1;
  X = (double*)malloc((size_t)nobs*sizeof(double));
  if (!X) FAIL_GOTO("Out of memory");
  for (int i = 0; i < nobs; ++i) X[i] = 1.0;

  /* varcmp length = kk+1 (residual last) */
  varcmp = (double*)calloc((size_t)(kk+1), sizeof(double));
  if (!varcmp) FAIL_GOTO("Out of memory");
  if (!isNull(varcmp_)) {
    if (!isReal(varcmp_) || (int)XLENGTH(varcmp_) != (kk+1))
      FAIL_GOTO("varcmp must be NULL or numeric length (kk+1)");
    for (int i = 0; i < kk+1; ++i) varcmp[i] = REAL(varcmp_)[i];
  }

  /* VV pointer array for st_mxd */
  VV = (double**)malloc((size_t)kk*sizeof(double*));
  if (!VV) FAIL_GOTO("Out of memory");
  VV[0] = K1_tr_rm;
  VV[1] = K2_tr_rm;
  if (include_sca) VV[2] = K3_tr_rm;

  int maxcyc = asInteger(maxcyc_);
  double tol = asReal(tol_);
  st_mxd(&nobs, ytr, &cX, X, (int*)&kk, NULL, VV, NULL, NULL, varcmp, &maxcyc, &tol, &retval);

  const double s1 = varcmp[0];
  const double s2 = varcmp[1];
  const double s3 = include_sca ? varcmp[2] : 0.0;
  const double se = varcmp[kk];

  /* ---- (2) Convert training blocks once to column-major for BLAS/LAPACK ---- */
  K1_tr_cm = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  K2_tr_cm = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  if (include_sca) K3_tr_cm = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  if (!K1_tr_cm || !K2_tr_cm || (include_sca && !K3_tr_cm)) FAIL_GOTO("Out of memory");

  rm_to_cm_sym(nobs, K1_tr_rm, K1_tr_cm);
  rm_to_cm_sym(nobs, K2_tr_rm, K2_tr_cm);
  if (include_sca) rm_to_cm_sym(nobs, K3_tr_rm, K3_tr_cm);

  /* ---- Build Sigma_tr in column-major ---- */
  Sigma = (double*)malloc((size_t)nobs*(size_t)nobs*sizeof(double));
  if (!Sigma) FAIL_GOTO("Out of memory");

  const size_t nnobs = (size_t)nobs*(size_t)nobs;
  for (size_t i = 0; i < nnobs; ++i) {
    double v = s1*K1_tr_cm[i] + s2*K2_tr_cm[i];
    if (include_sca) v += s3*K3_tr_cm[i];
    Sigma[i] = v;
  }
  for (int i = 0; i < nobs; ++i) Sigma[i + (size_t)nobs*(size_t)i] += se + ridge;

  /* ---- (1) Cholesky factorize Sigma and use dpotrs solves (NO inverse) ---- */
  {
    char uplo = 'U';
    int info = 0;
    F77_CALL(dpotrf)(&uplo, &nobs, Sigma, &nobs, &info FCONE);
    if (info != 0) FAILF_GOTO("Sigma dpotrf failed (info=%d).", info);
  }

  /* Solve Sigma * [Wy W1] = [y 1] with one dpotrs call (nrhs=2) */
  B = (double*)malloc((size_t)nobs * (size_t)2 * sizeof(double));
  r = (double*)malloc((size_t)nobs * sizeof(double));
  alpha = (double*)malloc((size_t)nobs * sizeof(double));
  if (!B || !r || !alpha) FAIL_GOTO("Out of memory");

  for (int i = 0; i < nobs; ++i) {
    B[i + (size_t)nobs*(size_t)0] = ytr[i];
    B[i + (size_t)nobs*(size_t)1] = 1.0;
  }

  {
    char uplo = 'U';
    int nrhs = 2;
    int info = 0;
    F77_CALL(dpotrs)(&uplo, &nobs, &nrhs, Sigma, &nobs, B, &nobs, &info FCONE);
    if (info != 0) FAILF_GOTO("Sigma dpotrs failed (info=%d).", info);
  }

  /* Wy and W1 are now in B[:,0], B[:,1] */
  double num = 0.0, den = 0.0;
  for (int i = 0; i < nobs; ++i) {
    num += B[i + (size_t)nobs*(size_t)0];
    den += B[i + (size_t)nobs*(size_t)1];
  }
  if (!(den > 0.0) || !(den == den)) FAILF_GOTO("Invalid denominator for beta0 (den=%g).", den);
  const double beta0 = num / den;

  for (int i = 0; i < nobs; ++i) r[i] = ytr[i] - beta0;

  /* Solve Sigma * alpha = r */
  {
    char uplo = 'U';
    int nrhs = 1;
    int info = 0;
    /* alpha <- r, then overwritten by solution */
    memcpy(alpha, r, (size_t)nobs*sizeof(double));
    F77_CALL(dpotrs)(&uplo, &nobs, &nrhs, Sigma, &nobs, alpha, &nobs, &info FCONE);
    if (info != 0) FAILF_GOTO("Sigma dpotrs(alpha) failed (info=%d).", info);
  }

  /* ---- (3) Build train×use cross-covariances in column-major ----
     K*_truse_cm has dim (nobs x nuse), column t is cov(train, use[t]) contiguous.
  */
  K1_truse_cm = (double*)malloc((size_t)nobs*(size_t)nuse*sizeof(double));
  K2_truse_cm = (double*)malloc((size_t)nobs*(size_t)nuse*sizeof(double));
  if (include_sca) K3_truse_cm = (double*)malloc((size_t)nobs*(size_t)nuse*sizeof(double));
  if (!K1_truse_cm || !K2_truse_cm || (include_sca && !K3_truse_cm)) FAIL_GOTO("Out of memory");

  for (int t = 0; t < nuse; ++t) {
    const int ib = p1u[t];
    const int jb = p2u[t];
    for (int a = 0; a < nobs; ++a) {
      const int ia = p1t[a];
      const int ja = p2t[a];
      double v1 = Kp1_cm[ ia + (size_t)nP1*(size_t)ib ];
      double v2 = Kp2_cm[ ja + (size_t)nP2*(size_t)jb ];
      K1_truse_cm[a + (size_t)nobs*(size_t)t] = v1;
      K2_truse_cm[a + (size_t)nobs*(size_t)t] = v2;
      if (include_sca) K3_truse_cm[a + (size_t)nobs*(size_t)t] = v1 * v2;
    }
  }

  /* ---- (4) Predict yhat_obs using ddot on contiguous columns ---- */
  yhat_obs = PROTECT(allocVector(REALSXP, ny)); nprot++;
  for (int i = 0; i < ny; ++i) REAL(yhat_obs)[i] = NA_REAL;

  {
    const BLAS_INT N = (BLAS_INT)nobs;
    const BLAS_INT INC = 1;

    for (int t = 0; t < nuse; ++t) {
      const double *c1 = K1_truse_cm + (size_t)nobs*(size_t)t;
      const double *c2 = K2_truse_cm + (size_t)nobs*(size_t)t;

      double pred = beta0;
      pred += s1 * F77_CALL(ddot)(&N, c1, &INC, alpha, &INC);
      pred += s2 * F77_CALL(ddot)(&N, c2, &INC, alpha, &INC);
      if (include_sca) {
        const double *c3 = K3_truse_cm + (size_t)nobs*(size_t)t;
        pred += s3 * F77_CALL(ddot)(&N, c3, &INC, alpha, &INC);
      }
      REAL(yhat_obs)[ use[t] ] = pred;
    }
  }

  /* outputs */
  beta0_out = PROTECT(allocVector(REALSXP, 1)); nprot++;
  REAL(beta0_out)[0] = beta0;

  varcmp_out = PROTECT(allocVector(REALSXP, kk+1)); nprot++;
  for (int i = 0; i < kk+1; ++i) REAL(varcmp_out)[i] = varcmp[i];

  diag = PROTECT(allocVector(VECSXP, 6)); nprot++;
  dn   = PROTECT(allocVector(STRSXP, 6)); nprot++;
  SET_STRING_ELT(dn, 0, mkChar("n_total"));
  SET_STRING_ELT(dn, 1, mkChar("n_matched"));
  SET_STRING_ELT(dn, 2, mkChar("n_train"));
  SET_STRING_ELT(dn, 3, mkChar("n_pred"));
  SET_STRING_ELT(dn, 4, mkChar("kk"));
  SET_STRING_ELT(dn, 5, mkChar("include_sca"));
  setAttrib(diag, R_NamesSymbol, dn);
  SET_VECTOR_ELT(diag, 0, ScalarInteger(ny));
  SET_VECTOR_ELT(diag, 1, ScalarInteger(nuse));
  SET_VECTOR_ELT(diag, 2, ScalarInteger(nobs));
  SET_VECTOR_ELT(diag, 3, ScalarInteger(npred));
  SET_VECTOR_ELT(diag, 4, ScalarInteger(kk));
  SET_VECTOR_ELT(diag, 5, ScalarLogical(include_sca));

  ans = PROTECT(allocVector(VECSXP, 6)); nprot++;
  an  = PROTECT(allocVector(STRSXP, 6)); nprot++;
  SET_STRING_ELT(an, 0, mkChar("beta0"));
  SET_STRING_ELT(an, 1, mkChar("varcmp"));
  SET_STRING_ELT(an, 2, mkChar("yhat_obs"));
  SET_STRING_ELT(an, 3, mkChar("diag"));
  SET_STRING_ELT(an, 4, mkChar("retval"));
  SET_STRING_ELT(an, 5, mkChar("ridge"));
  setAttrib(ans, R_NamesSymbol, an);

  SET_VECTOR_ELT(ans, 0, beta0_out);
  SET_VECTOR_ELT(ans, 1, varcmp_out);
  SET_VECTOR_ELT(ans, 2, yhat_obs);
  SET_VECTOR_ELT(ans, 3, diag);
  SET_VECTOR_ELT(ans, 4, ScalarInteger(retval));
  SET_VECTOR_ELT(ans, 5, ScalarReal(ridge));

  goto cleanup;

cleanup_err:
  ans = R_NilValue;
  goto cleanup;

cleanup:
  free(train_pos);

  free(Z1rm); free(Z2rm);
  free(Kp1_rm); free(Kp2_rm);
  free(Kp1_cm); free(Kp2_cm);

  free(K1_tr_rm); free(K2_tr_rm); free(K3_tr_rm);
  free(VV);

  free(ytr);
  free(X);
  free(varcmp);

  free(K1_tr_cm); free(K2_tr_cm); free(K3_tr_cm);
  free(Sigma);

  free(K1_truse_cm); free(K2_truse_cm); free(K3_truse_cm);

  free(B);
  free(r);
  free(alpha);

  UNPROTECT(nprot);

#undef FAIL_GOTO
#undef FAILF_GOTO
  return ans;
}


/* small helpers */
static SEXP getListElement(SEXP list, const char *str)
{
  SEXP names = getAttrib(list, R_NamesSymbol);
  if (names == R_NilValue) return R_NilValue;

  for (R_xlen_t i = 0; i < XLENGTH(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

//
#include <Rmath.h>  /* GetRNGstate/PutRNGstate/unif_rand */

/* ---------------- model_type parsing ---------------- */

typedef enum { MODEL_KERNEL=0, MODEL_MARKER=1 } ModelType;

static int parse_model_type(SEXP model_type_, ModelType *out) {
  *out = MODEL_KERNEL;
  if (model_type_ == R_NilValue) return 0;
  if (!isString(model_type_) || XLENGTH(model_type_) != 1) return -3001;
  const char *s = CHAR(STRING_ELT(model_type_, 0));
  if (strcmp(s, "kernel") == 0) { *out = MODEL_KERNEL; return 0; }
  if (strcmp(s, "marker") == 0) { *out = MODEL_MARKER; return 0; }
  return -3002;
}

/* ---------------- bootstrap parsing ---------------- */
typedef struct {
  int nObs;       /* 0 means: use full training set */
  int nMarkers;   /* 0 means: use full marker set per pool */
  int nRep;
  int enabled;
} BootSpec;



/* ============================================================
   Sigma operator for training hybrids:
   Sigma = s1*K1 + s2*K2 + s3*K3 + se*I
   where K1_ij = Kp1[a_i,a_j], K2_ij = Kp2[b_i,b_j], K3_ij = K1_ij*K2_ij
   ============================================================ */

typedef struct {
  int nobs;
  int nP1, nP2;
  const int *a;             /* length nobs */
  const int *b;             /* length nobs */
  const double *Kp1_cm;     /* col-major nP1 x nP1 */
  const double *Kp2_cm;     /* col-major nP2 x nP2 */
  int include_sca;
  double s1, s2, s3, se;

  /* work buffers length nobs */
  double *tmp1;
  double *tmp2;
  double *tmp3;             /* only if include_sca */
} SigmaOp;

/* compute out = K1 * v (nobs) without forming K1 */
static void matvec_K1(const SigmaOp *op, const double *v, double *out)
{
  /* aggregate by parent a: sA[u] = sum_{j:a_j=u} v_j */
  double *sA = (double*)calloc((size_t)op->nP1, sizeof(double));
  if (!sA) Rf_error("OOM in matvec_K1");
  for (int j=0;j<op->nobs;++j) sA[ op->a[j] ] += v[j];

  /* tA = Kp1 * sA */
  double *tA = (double*)malloc((size_t)op->nP1*sizeof(double));
  if (!tA) Rf_error("OOM in matvec_K1");
  {
    char uplo='U';
    BLAS_INT N=(BLAS_INT)op->nP1, INC=1;
    const double one=1.0, zero=0.0;
    F77_CALL(dsymv)(&uplo, &N, &one, (double*)op->Kp1_cm, &N,
                    sA, &INC, &zero, tA, &INC FCONE);
  }

  for (int i=0;i<op->nobs;++i) out[i] = tA[ op->a[i] ];

  free(tA);
  free(sA);
}

/* out = K2 * v */
static void matvec_K2(const SigmaOp *op, const double *v, double *out)
{
  double *sB = (double*)calloc((size_t)op->nP2, sizeof(double));
  if (!sB) Rf_error("OOM in matvec_K2");
  for (int j=0;j<op->nobs;++j) sB[ op->b[j] ] += v[j];

  double *tB = (double*)malloc((size_t)op->nP2*sizeof(double));
  if (!tB) Rf_error("OOM in matvec_K2");
  {
    char uplo='U';
    BLAS_INT N=(BLAS_INT)op->nP2, INC=1;
    const double one=1.0, zero=0.0;
    F77_CALL(dsymv)(&uplo, &N, &one, (double*)op->Kp2_cm, &N,
                    sB, &INC, &zero, tB, &INC FCONE);
  }

  for (int i=0;i<op->nobs;++i) out[i] = tB[ op->b[i] ];

  free(tB);
  free(sB);
}

/* out = K3 * v, with K3 = K1 ? K2 (Hadamard) */
static void matvec_K3(const SigmaOp *op, const double *v, double *out)
{
  const int nobs = op->nobs;
  const int nP1 = op->nP1, nP2 = op->nP2;

  /* A_v[a,b] = sum_{j: (a_j,b_j)=(a,b)} v_j  (col-major nP1 x nP2) */
  double *Av = (double*)calloc((size_t)nP1*(size_t)nP2, sizeof(double));
  if (!Av) Rf_error("OOM in matvec_K3 Av");

  for (int j=0;j<nobs;++j) {
    int a = op->a[j];
    int b = op->b[j];
    Av[a + (size_t)nP1*(size_t)b] += v[j];
  }

  /* W = Kp1 * Av  (nP1 x nP2) */
  double *W = (double*)malloc((size_t)nP1*(size_t)nP2*sizeof(double));
  if (!W) Rf_error("OOM in matvec_K3 W");
  {
    char transN='N';
    BLAS_INT M=(BLAS_INT)nP1, N=(BLAS_INT)nP2, K=(BLAS_INT)nP1;
    BLAS_INT lda=M, ldb=M, ldc=M;
    const double one=1.0, zero=0.0;
    F77_CALL(dgemm)(&transN,&transN,&M,&N,&K,&one,
                    (double*)op->Kp1_cm,&lda,
                    Av,&ldb,
                    &zero,
                    W,&ldc FCONE FCONE);
  }

  /* T = W * Kp2^T  (nP1 x nP2) */
  double *T = (double*)malloc((size_t)nP1*(size_t)nP2*sizeof(double));
  if (!T) Rf_error("OOM in matvec_K3 T");
  {
    char transN='N', transT='T';
    BLAS_INT M=(BLAS_INT)nP1, N=(BLAS_INT)nP2, K=(BLAS_INT)nP2;
    BLAS_INT lda=(BLAS_INT)nP1, ldb=(BLAS_INT)nP2, ldc=(BLAS_INT)nP1;
    const double one=1.0, zero=0.0;
    F77_CALL(dgemm)(&transN,&transT,&M,&N,&K,&one,
                    W,&lda,
                    (double*)op->Kp2_cm,&ldb,
                    &zero,
                    T,&ldc FCONE FCONE);
  }

  for (int i=0;i<nobs;++i) {
    int a = op->a[i], b = op->b[i];
    out[i] = T[a + (size_t)nP1*(size_t)b];
  }

  free(T);
  free(W);
  free(Av);
}

static void Sigma_times_vec(const SigmaOp *op, const double *v, double *out)
{
  const int nobs = op->nobs;

  /* out = s1*K1*v + s2*K2*v + s3*K3*v + se*v */
  matvec_K1(op, v, op->tmp1);
  matvec_K2(op, v, op->tmp2);

  for (int i=0;i<nobs;++i) out[i] = op->s1*op->tmp1[i] + op->s2*op->tmp2[i];

  if (op->include_sca) {
    matvec_K3(op, v, op->tmp3);
    for (int i=0;i<nobs;++i) out[i] += op->s3*op->tmp3[i];
  }

  for (int i=0;i<nobs;++i) out[i] += op->se * v[i];
}

/* Conjugate Gradient solve: Sigma x = b */
static int cg_solve_sigma(const SigmaOp *op, const double *b, double *x, int maxit, double tol)
{
  const int n = op->nobs;
  double *r = (double*)malloc((size_t)n*sizeof(double));
  double *p = (double*)malloc((size_t)n*sizeof(double));
  double *Ap= (double*)malloc((size_t)n*sizeof(double));
  if (!r || !p || !Ap) { free(r); free(p); free(Ap); return -1010; }

  /* r = b - A x  (x assumed initialized by caller) */
  Sigma_times_vec(op, x, Ap);
  for (int i=0; i<n; ++i) {
    r[i] = b[i] - Ap[i];
    p[i] = r[i];
  }

  double rr = 0.0;
  for (int i=0;i<n;++i) rr += r[i]*r[i];
  if (!(rr == rr)) { free(r); free(p); free(Ap); return -1; }

  const double rr0 = rr;
  const double thresh = (tol > 0.0 ? tol*tol*rr0 : 0.0);

  for (int it=0; it<maxit; ++it) {
    Sigma_times_vec(op, p, Ap);

    double pAp = 0.0;
    for (int i=0;i<n;++i) pAp += p[i]*Ap[i];
    if (!(pAp > 0.0) || !(pAp == pAp)) { free(r); free(p); free(Ap); return -2; }

    const double alpha = rr / pAp;
    for (int i=0;i<n;++i) x[i] += alpha * p[i];
    for (int i=0;i<n;++i) r[i] -= alpha * Ap[i];

    double rr_new = 0.0;
    for (int i=0;i<n;++i) rr_new += r[i]*r[i];

    if (rr_new <= thresh) { free(r); free(p); free(Ap); return 0; }

    const double beta = rr_new / rr;
    for (int i=0;i<n;++i) p[i] = r[i] + beta * p[i];
    rr = rr_new;
  }

  free(r); free(p); free(Ap);
  return 1; /* no convergence */
}

static int hybrid_marker_fullpcg_core(
  int ny,
  int nuse,
  int nobs,
  const int *use,          /* length nuse: original-row indices of matched hybrids */
  const int *p1_idx,       /* length ny: pool1 parent index for each original row */
  const int *p2_idx,       /* length ny: pool2 parent index for each original row */
  const int *a_tr,         /* length nobs: pool1 parent index for training rows */
  const int *b_tr,         /* length nobs: pool2 parent index for training rows */
  const double *ytr,       /* length nobs */
  int nP1,
  int nP2,
  const double *Kp1_cm,    /* col-major nP1 x nP1 */
  const double *Kp2_cm,    /* col-major nP2 x nP2 */
  int include_sca,
  double s1, double s2, double s3, double se,
  double tol_user,
  double *out_beta0,
  double *out_yhat_obs     /* length ny; caller allocates */
){
  SigmaOp op;
  op.nobs = nobs;
  op.nP1 = nP1;
  op.nP2 = nP2;
  op.a = a_tr;
  op.b = b_tr;
  op.Kp1_cm = Kp1_cm;
  op.Kp2_cm = Kp2_cm;
  op.include_sca = include_sca;
  op.s1 = s1;
  op.s2 = s2;
  op.s3 = include_sca ? s3 : 0.0;
  op.se = se;

  op.tmp1 = (double*)malloc((size_t)nobs*sizeof(double));
  op.tmp2 = (double*)malloc((size_t)nobs*sizeof(double));
  op.tmp3 = include_sca ? (double*)malloc((size_t)nobs*sizeof(double)) : NULL;
  if (!op.tmp1 || !op.tmp2 || (include_sca && !op.tmp3)) {
    free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
    return -1010;
  }

  /* strict internal CG settings */
  double tol_cg = 1e-10;
  if ((tol_user > 0.0) && (tol_user == tol_user) && (tol_user < tol_cg)) {
    tol_cg = tol_user;
  }
  int maxit = 4000;

  double *Wy   = (double*)calloc((size_t)nobs, sizeof(double));
  double *W1   = (double*)calloc((size_t)nobs, sizeof(double));
  double *ones = (double*)malloc((size_t)nobs*sizeof(double));
  double *r    = (double*)malloc((size_t)nobs*sizeof(double));
  double *alpha= (double*)calloc((size_t)nobs, sizeof(double));

  double *sA = NULL, *sB = NULL, *tA = NULL, *tB = NULL, *Tsca = NULL;

  if (!Wy || !W1 || !ones || !r || !alpha) {
    free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
    free(Wy); free(W1); free(ones); free(r); free(alpha);
    return -1010;
  }

  for (int i=0;i<nobs;++i) ones[i] = 1.0;

  {
    int rc = cg_solve_sigma(&op, ytr, Wy, maxit, tol_cg);
    if (rc != 0) {
      free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
      free(Wy); free(W1); free(ones); free(r); free(alpha);
      return -2001;
    }
  }

  {
    int rc = cg_solve_sigma(&op, ones, W1, maxit, tol_cg);
    if (rc != 0) {
      free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
      free(Wy); free(W1); free(ones); free(r); free(alpha);
      return -2002;
    }
  }

  double num = 0.0, den = 0.0;
  for (int i=0;i<nobs;++i) {
    num += Wy[i];
    den += W1[i];
  }
  if (!(den > 0.0) || !(den == den)) {
    free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
    free(Wy); free(W1); free(ones); free(r); free(alpha);
    return -2003;
  }

  const double beta0 = num / den;
  *out_beta0 = beta0;

  for (int i=0;i<nobs;++i) r[i] = ytr[i] - beta0;

  {
    int rc = cg_solve_sigma(&op, r, alpha, maxit, tol_cg);
    if (rc != 0) {
      free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
      free(Wy); free(W1); free(ones); free(r); free(alpha);
      return -2004;
    }
  }

  /* grouped parent sums */
  sA = (double*)calloc((size_t)nP1, sizeof(double));
  sB = (double*)calloc((size_t)nP2, sizeof(double));
  tA = (double*)malloc((size_t)nP1*sizeof(double));
  tB = (double*)malloc((size_t)nP2*sizeof(double));

  if (!sA || !sB || !tA || !tB) {
    free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
    free(Wy); free(W1); free(ones); free(r); free(alpha);
    free(sA); free(sB); free(tA); free(tB);
    return -1010;
  }

  for (int j=0;j<nobs;++j) {
    sA[a_tr[j]] += alpha[j];
    sB[b_tr[j]] += alpha[j];
  }

  {
    char uplo='U';
    BLAS_INT INC=1;
    const double one=1.0, zero=0.0;
    BLAS_INT N1=(BLAS_INT)nP1, N2=(BLAS_INT)nP2;
    F77_CALL(dsymv)(&uplo, &N1, &one, (double*)Kp1_cm, &N1, sA, &INC, &zero, tA, &INC FCONE);
    F77_CALL(dsymv)(&uplo, &N2, &one, (double*)Kp2_cm, &N2, sB, &INC, &zero, tB, &INC FCONE);
  }

  if (include_sca) {
    const size_t R = (size_t)nP1*(size_t)nP2;
    double *Aalpha = (double*)calloc(R, sizeof(double)); /* col-major nP1 x nP2 */
    double *W = (double*)malloc(R * sizeof(double));
    Tsca = (double*)malloc(R * sizeof(double));
    if (!Aalpha || !W || !Tsca) {
      free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
      free(Wy); free(W1); free(ones); free(r); free(alpha);
      free(sA); free(sB); free(tA); free(tB);
      free(Aalpha); free(W); free(Tsca);
      return -1010;
    }

    for (int j=0;j<nobs;++j) {
      Aalpha[a_tr[j] + (size_t)nP1*(size_t)b_tr[j]] += alpha[j];
    }

    {
      char transN='N';
      BLAS_INT M=(BLAS_INT)nP1, N=(BLAS_INT)nP2, K=(BLAS_INT)nP1;
      BLAS_INT lda=M, ldb=M, ldc=M;
      const double one=1.0, zero=0.0;
      F77_CALL(dgemm)(&transN,&transN,&M,&N,&K,&one,
                      (double*)Kp1_cm,&lda,
                      Aalpha,&ldb,
                      &zero,
                      W,&ldc FCONE FCONE);
    }

    {
      char transN='N', transT='T';
      BLAS_INT M=(BLAS_INT)nP1, N=(BLAS_INT)nP2, K=(BLAS_INT)nP2;
      BLAS_INT lda=(BLAS_INT)nP1, ldb=(BLAS_INT)nP2, ldc=(BLAS_INT)nP1;
      const double one=1.0, zero=0.0;
      F77_CALL(dgemm)(&transN,&transT,&M,&N,&K,&one,
                      W,&lda,
                      (double*)Kp2_cm,&ldb,
                      &zero,
                      Tsca,&ldc FCONE FCONE);
    }

    free(Aalpha);
    free(W);
  }

  /* final predictions for all matched hybrids */
  for (int i=0;i<ny;++i) out_yhat_obs[i] = NA_REAL;

  for (int t=0; t<nuse; ++t) {
    int obs = use[t];
    int a0 = p1_idx[obs];
    int b0 = p2_idx[obs];
    double pred = beta0;
    pred += s1 * tA[a0];
    pred += s2 * tB[b0];
    if (include_sca) pred += s3 * Tsca[a0 + (size_t)nP1*(size_t)b0];
    out_yhat_obs[obs] = pred;
  }

  free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
  free(Wy); free(W1); free(ones); free(r); free(alpha);
  free(sA); free(sB); free(tA); free(tB);
  if (Tsca) free(Tsca);

  return 0;
}

static int parse_bootstrap(SEXP bootstrap_, BootSpec *bs)
{
  bs->enabled = 0;
  bs->nObs = 0;
  bs->nMarkers = 0;
  bs->nRep = 1;

  if (bootstrap_ == R_NilValue) return 0;
  if (TYPEOF(bootstrap_) != VECSXP) return -3101;

  SEXP nObs_ = getListElement(bootstrap_, "nObs");
  SEXP nMarkers_ = getListElement(bootstrap_, "nMarkers");
  SEXP nRep_ = getListElement(bootstrap_, "nRep");

  if (nRep_ != R_NilValue) {
    int nr = asInteger(nRep_);
    if (nr == NA_INTEGER || nr < 1) return -3104;
    bs->nRep = nr;
  }

  if (nObs_ != R_NilValue) {
    int no = asInteger(nObs_);
    if (no == NA_INTEGER || no < 1) return -3102;
    bs->nObs = no;
  } else {
    bs->nObs = 0;
  }

  if (nMarkers_ != R_NilValue) {
    int nm = asInteger(nMarkers_);
    if (nm == NA_INTEGER || nm < 1) return -3103;
    bs->nMarkers = nm;
  } else {
    bs->nMarkers = 0;
  }

  bs->enabled = 1;
  return 0;
}

static int sample_int_without_replacement(int n, int k, int *out)
{
  if (n <= 0 || k < 0) return -1;

  if (k >= n) {
    for (int i = 0; i < n; ++i) out[i] = i;
    return n;
  }

  int *pool = (int*)malloc((size_t)n * sizeof(int));
  if (!pool) return -1010;

  for (int i = 0; i < n; ++i) pool[i] = i;

  for (int i = 0; i < k; ++i) {
    int j = i + (int)floor(unif_rand() * (double)(n - i));
    if (j < i) j = i;
    if (j >= n) j = n - 1;
    int tmp = pool[i];
    pool[i] = pool[j];
    pool[j] = tmp;
    out[i] = pool[i];
  }

  free(pool);
  return k;
}

static int bootstrap_estimate_varcmp_hybrid(
  int nobs_full,
  const double *ytr_full,
  const int *a_full,
  const int *b_full,
  int nP1, int m1,
  int nP2, int m2,
  const double *Z1_cm,
  const double *Z2_cm,
  KMethod km,
  SEXP coding_,
  int include_sca,   /* must be 1 */
  int maxcyc,
  double tol,
  double ridge,
  const BootSpec *bs,
  double *varcmp_out,  /* length 4: s1,s2,s3,se (se excludes ridge) */
  int *out_n_success
){
  if (!bs || !bs->enabled) return -3200;
  if (!include_sca) return -3201;
  if (bs->nRep <= 0) return -3202;

  const int use_full_obs = (bs->nObs <= 0 || bs->nObs >= nobs_full);
  const int use_full_m1  = (bs->nMarkers <= 0 || bs->nMarkers >= m1);
  const int use_full_m2  = (bs->nMarkers <= 0 || bs->nMarkers >= m2);

  const int nObs_use = use_full_obs ? nobs_full : bs->nObs;
  const int nM1_use  = use_full_m1  ? m1        : bs->nMarkers;
  const int nM2_use  = use_full_m2  ? m2        : bs->nMarkers;

  double acc[4] = {0,0,0,0};
  int n_ok = 0;

  int *sampObs = (int*)malloc((size_t)nObs_use * sizeof(int));
  int *sampM1  = (int*)malloc((size_t)nM1_use  * sizeof(int));
  int *sampM2  = (int*)malloc((size_t)nM2_use  * sizeof(int));

  double *y_bs = (double*)malloc((size_t)nObs_use * sizeof(double));
  int *a_bs = (int*)malloc((size_t)nObs_use * sizeof(int));
  int *b_bs = (int*)malloc((size_t)nObs_use * sizeof(int));

  double *Z1rm_sub = (double*)malloc((size_t)nP1*(size_t)nM1_use*sizeof(double));
  double *Z2rm_sub = (double*)malloc((size_t)nP2*(size_t)nM2_use*sizeof(double));

  double *Kp1_rm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_rm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  double *Kp1_cm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_cm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));

  double *K1_tr_rm = (double*)malloc((size_t)nObs_use*(size_t)nObs_use*sizeof(double));
  double *K2_tr_rm = (double*)malloc((size_t)nObs_use*(size_t)nObs_use*sizeof(double));
  double *K3_tr_rm = (double*)malloc((size_t)nObs_use*(size_t)nObs_use*sizeof(double));
  double *X = (double*)malloc((size_t)nObs_use*sizeof(double));

  if (!sampObs||!sampM1||!sampM2||!y_bs||!a_bs||!b_bs||!Z1rm_sub||!Z2rm_sub||
      !Kp1_rm||!Kp2_rm||!Kp1_cm||!Kp2_cm||!K1_tr_rm||!K2_tr_rm||!K3_tr_rm||!X) {
    free(sampObs); free(sampM1); free(sampM2);
    free(y_bs); free(a_bs); free(b_bs);
    free(Z1rm_sub); free(Z2rm_sub);
    free(Kp1_rm); free(Kp2_rm); free(Kp1_cm); free(Kp2_cm);
    free(K1_tr_rm); free(K2_tr_rm); free(K3_tr_rm);
    free(X);
    return -1010;
  }
  for (int i=0;i<nObs_use;++i) X[i]=1.0;

  double *VV_arr[3]; VV_arr[0]=K1_tr_rm; VV_arr[1]=K2_tr_rm; VV_arr[2]=K3_tr_rm;
  int cX_i = 1;
  int kk_i = 3;

  GetRNGstate();

  for (int rep=0; rep<bs->nRep; ++rep) {
    if (use_full_obs) {
      for (int i = 0; i < nObs_use; ++i) sampObs[i] = i;
    } else {
      int rc = sample_int_without_replacement(nobs_full, nObs_use, sampObs);
      if (rc < 0) continue;
    }

    if (use_full_m1) {
      for (int j = 0; j < nM1_use; ++j) sampM1[j] = j;
    } else {
      int rc = sample_int_without_replacement(m1, nM1_use, sampM1);
      if (rc < 0) continue;
    }

    if (use_full_m2) {
      for (int j = 0; j < nM2_use; ++j) sampM2[j] = j;
    } else {
      int rc = sample_int_without_replacement(m2, nM2_use, sampM2);
      if (rc < 0) continue;
    }

    for (int i=0;i<nObs_use;++i) {
      int src = sampObs[i];
      y_bs[i] = ytr_full[src];
      a_bs[i] = a_full[src];
      b_bs[i] = b_full[src];
    }

    for (int i=0;i<nP1;++i) {
      for (int jj=0;jj<nM1_use;++jj) {
        int col = sampM1[jj];
        Z1rm_sub[(size_t)i*(size_t)nM1_use + (size_t)jj] =
          Z1_cm[(size_t)i + (size_t)nP1*(size_t)col];
      }
    }

    for (int i=0;i<nP2;++i) {
      for (int jj=0;jj<nM2_use;++jj) {
        int col = sampM2[jj];
        Z2rm_sub[(size_t)i*(size_t)nM2_use + (size_t)jj] =
          Z2_cm[(size_t)i + (size_t)nP2*(size_t)col];
      }
    }

    double zmean1=0,zsd1=0,zmean2=0,zsd2=0;
    ZCoding cu1=CODING_OTHER, cu2=CODING_OTHER;
    ZTransform tr1=TRANS_NONE, tr2=TRANS_NONE;

    int rc1 = build_K_from_Z(Z1rm_sub, nP1, nM1_use, km, coding_, Kp1_rm,
                             &zmean1, &zsd1, &cu1, &tr1);
    int rc2 = build_K_from_Z(Z2rm_sub, nP2, nM2_use, km, coding_, Kp2_rm,
                             &zmean2, &zsd2, &cu2, &tr2);
    if (rc1 != 0 || rc2 != 0) continue;

    rm_to_cm_sym(nP1, Kp1_rm, Kp1_cm);
    rm_to_cm_sym(nP2, Kp2_rm, Kp2_cm);

    for (int i=0;i<nObs_use;++i) {
      int ia = a_bs[i];
      int ja = b_bs[i];
      for (int j=0;j<nObs_use;++j) {
        int ib = a_bs[j];
        int jb = b_bs[j];
        double v1 = Kp1_cm[ ia + (size_t)nP1*(size_t)ib ];
        double v2 = Kp2_cm[ ja + (size_t)nP2*(size_t)jb ];
        K1_tr_rm[(size_t)i*(size_t)nObs_use + (size_t)j] = v1;
        K2_tr_rm[(size_t)i*(size_t)nObs_use + (size_t)j] = v2;
        K3_tr_rm[(size_t)i*(size_t)nObs_use + (size_t)j] = v1 * v2;
      }
    }

    double varcmp[4] = {0.0,0.0,0.0,0.0};
    int n_i = nObs_use, retval = 0;
    int maxcyc_i = maxcyc;
    double tol_i = tol;

    st_mxd(&n_i, y_bs, &cX_i, X, &kk_i,
           NULL, VV_arr, NULL, NULL, varcmp,
           &maxcyc_i, &tol_i, &retval);

    if (retval != 0) continue;

    int ok = 1;
    for (int t=0;t<4;++t) {
      double v = varcmp[t];
      if (!(v==v) || v < 0.0) { ok = 0; break; }
    }
    if (!ok) continue;
    if (!(varcmp[3] > 0.0)) continue;

    acc[0] += varcmp[0];
    acc[1] += varcmp[1];
    acc[2] += varcmp[2];
    acc[3] += varcmp[3];
    n_ok++;
  }

  PutRNGstate();

  if (out_n_success) *out_n_success = n_ok;

  if (n_ok <= 0) {
    free(sampObs); free(sampM1); free(sampM2);
    free(y_bs); free(a_bs); free(b_bs);
    free(Z1rm_sub); free(Z2rm_sub);
    free(Kp1_rm); free(Kp2_rm); free(Kp1_cm); free(Kp2_cm);
    free(K1_tr_rm); free(K2_tr_rm); free(K3_tr_rm);
    free(X);
    return -3210;
  }

  varcmp_out[0] = acc[0] / (double)n_ok;
  varcmp_out[1] = acc[1] / (double)n_ok;
  varcmp_out[2] = acc[2] / (double)n_ok;
  varcmp_out[3] = acc[3] / (double)n_ok;

  free(sampObs); free(sampM1); free(sampM2);
  free(y_bs); free(a_bs); free(b_bs);
  free(Z1rm_sub); free(Z2rm_sub);
  free(Kp1_rm); free(Kp2_rm); free(Kp1_cm); free(Kp2_cm);
  free(K1_tr_rm); free(K2_tr_rm); free(K3_tr_rm);
  free(X);

  return 0;
}

static SEXP st_mixed_hybrid_marker_impl(
  SEXP y_, SEXP hybrid_id_,
  SEXP pool1_, SEXP pool2_,
  SEXP Z1_, SEXP Z1_ids_,
  SEXP Z2_, SEXP Z2_ids_,
  SEXP include_sca_,
  SEXP method_, SEXP coding_, SEXP varcmp_,
  SEXP maxcyc_, SEXP tol_, SEXP ridge_,
  SEXP bootstrap_)
{
  int nprot = 0;

  if (!isReal(y_) || !isVector(y_)) Rf_error("marker: y must be numeric vector");
  const int ny = (int)XLENGTH(y_);

  if (!isString(hybrid_id_) || (int)XLENGTH(hybrid_id_) != ny) Rf_error("marker: hybrid_id must be character length(y)");
  if (!isString(pool1_) || (int)XLENGTH(pool1_) != ny) Rf_error("marker: pool1 must be character length(y)");
  if (!isString(pool2_) || (int)XLENGTH(pool2_) != ny) Rf_error("marker: pool2 must be character length(y)");

  if (!isReal(Z1_) || !isMatrix(Z1_)) Rf_error("marker: Z1 must be numeric matrix");
  if (!isReal(Z2_) || !isMatrix(Z2_)) Rf_error("marker: Z2 must be numeric matrix");
  if (!isString(Z1_ids_) || XLENGTH(Z1_ids_) < 1) Rf_error("marker: Z1_ids must be character");
  if (!isString(Z2_ids_) || XLENGTH(Z2_ids_) < 1) Rf_error("marker: Z2_ids must be character");

  if (!isString(method_) || XLENGTH(method_) != 1) Rf_error("marker: method must be single string");
  if (!isString(coding_) || XLENGTH(coding_) != 1) Rf_error("marker: coding must be single string");

  const int include_sca = asLogical(include_sca_) ? 1 : 0;
  const double ridge = asReal(ridge_);
  if (!(ridge >= 0.0) || !(ridge == ridge)) Rf_error("marker: ridge must be >= 0");

  SEXP d1 = getAttrib(Z1_, R_DimSymbol);
  SEXP d2 = getAttrib(Z2_, R_DimSymbol);
  const int nP1 = INTEGER(d1)[0], m1 = INTEGER(d1)[1];
  const int nP2 = INTEGER(d2)[0], m2 = INTEGER(d2)[1];
  if ((int)XLENGTH(Z1_ids_) != nP1) Rf_error("marker: length(Z1_ids) must equal nrow(Z1)");
  if ((int)XLENGTH(Z2_ids_) != nP2) Rf_error("marker: length(Z2_ids) must equal nrow(Z2)");

  KMethod km = K_STD;
  if (str_eq(method_, "std")) km = K_STD;
  else if (str_eq(method_, "vanraden")) km = K_VANRADEN;
  else if (str_eq(method_, "sommer")) km = K_SOMMER;
  else Rf_error("marker: method must be one of: 'std','vanraden','sommer'");

  const int kk = include_sca ? 3 : 2;

  int *p1_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *p2_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *use    = (int*)R_alloc((size_t)ny, sizeof(int));
  int nuse = 0;

  for (int i = 0; i < ny; ++i) {
    const char *a = CHAR(STRING_ELT(pool1_, i));
    const char *b = CHAR(STRING_ELT(pool2_, i));
    int ia = find_id_index(Z1_ids_, a);
    int ib = find_id_index(Z2_ids_, b);
    p1_idx[i] = ia;
    p2_idx[i] = ib;
    if (ia >= 0 && ib >= 0) use[nuse++] = i;
  }
  if (nuse <= 0) Rf_error("marker: No rows matched to both pool1 and pool2 IDs.");

  int *train_pos = (int*)malloc((size_t)nuse * sizeof(int));
  if (!train_pos) Rf_error("marker: OOM train_pos");
  int nobs = 0, npred = 0;
  for (int t = 0; t < nuse; ++t) {
    int i = use[t];
    if (!is_missing(REAL(y_)[i])) train_pos[nobs++] = t;
    else npred++;
  }
  if (nobs <= 0) {
    free(train_pos);
    Rf_error("marker: No non-missing y values after matching parents.");
  }

  int *a_tr = (int*)malloc((size_t)nobs * sizeof(int));
  int *b_tr = (int*)malloc((size_t)nobs * sizeof(int));
  double *ytr = (double*)malloc((size_t)nobs * sizeof(double));
  if (!a_tr || !b_tr || !ytr) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: OOM a_tr/b_tr/ytr");
  }
  for (int i = 0; i < nobs; ++i) {
    int obs = use[ train_pos[i] ];
    a_tr[i] = p1_idx[obs];
    b_tr[i] = p2_idx[obs];
    ytr[i]  = REAL(y_)[obs];
  }

  /* variance components: provided or bootstrap */
  double s1=NAN, s2=NAN, s3=0.0, se0=NAN;
  int vc_source_boot = 0;
  int boot_success = 0;

  if (!isNull(varcmp_)) {
    if (!isReal(varcmp_) || (int)XLENGTH(varcmp_) != (kk + 1)) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp must be numeric length (kk+1)=%d", kk+1);
    }
    s1  = REAL(varcmp_)[0];
    s2  = REAL(varcmp_)[1];
    s3  = include_sca ? REAL(varcmp_)[2] : 0.0;
    se0 = REAL(varcmp_)[kk];
  } else {
    if (!include_sca) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp NULL; bootstrap VC only enabled when include_sca=TRUE.");
    }

    BootSpec bs;
    int pr = parse_bootstrap(bootstrap_, &bs);
    if (pr != 0 || !bs.enabled) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp NULL and bootstrap not set correctly.");
    }

    double vc_est[4] = {0,0,0,0};
    int rc = bootstrap_estimate_varcmp_hybrid(
      nobs, ytr, a_tr, b_tr,
      nP1, m1, nP2, m2,
      REAL(Z1_), REAL(Z2_),
      km, coding_, 1,
      asInteger(maxcyc_), asReal(tol_), ridge,
      &bs,
      vc_est, &boot_success
    );
    if (rc != 0) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: bootstrap VC failed (rc=%d).", rc);
    }

    s1 = vc_est[0];
    s2 = vc_est[1];
    s3 = vc_est[2];
    se0= vc_est[3];
    vc_source_boot = 1;
  }

  if (!(s1 >= 0.0) || !(s2 >= 0.0) || !(s3 >= 0.0) || !(se0 > 0.0) ||
      !(s1==s1) || !(s2==s2) || !(s3==s3) || !(se0==se0)) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: invalid variance components.");
  }

  const double se = se0 + ridge;
  if (!(se > 0.0) || !(se==se)) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: residual variance + ridge must be > 0.");
  }

  /* parental kernels from full markers */
  double *Z1rm   = (double*)malloc((size_t)nP1*(size_t)m1*sizeof(double));
  double *Z2rm   = (double*)malloc((size_t)nP2*(size_t)m2*sizeof(double));
  double *Kp1_rm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_rm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  double *Kp1_cm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_cm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  if (!Z1rm||!Z2rm||!Kp1_rm||!Kp2_rm||!Kp1_cm||!Kp2_cm) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    free(Z1rm); free(Z2rm); free(Kp1_rm); free(Kp2_rm); free(Kp1_cm); free(Kp2_cm);
    Rf_error("marker: OOM kernels");
  }

  Z_cm_to_rm(REAL(Z1_), nP1, m1, Z1rm);
  Z_cm_to_rm(REAL(Z2_), nP2, m2, Z2rm);

  double zmean1=0,zsd1=0,zmean2=0,zsd2=0;
  ZCoding cu1=CODING_OTHER, cu2=CODING_OTHER;
  ZTransform tr1=TRANS_NONE, tr2=TRANS_NONE;

  int rc1 = build_K_from_Z(Z1rm, nP1, m1, km, coding_, Kp1_rm, &zmean1, &zsd1, &cu1, &tr1);
  if (rc1 != 0) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    free(Z1rm); free(Z2rm); free(Kp1_rm); free(Kp2_rm); free(Kp1_cm); free(Kp2_cm);
    Rf_error("marker: build_K_from_Z(Z1) failed (code %d).", rc1);
  }

  int rc2 = build_K_from_Z(Z2rm, nP2, m2, km, coding_, Kp2_rm, &zmean2, &zsd2, &cu2, &tr2);
  if (rc2 != 0) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    free(Z1rm); free(Z2rm); free(Kp1_rm); free(Kp2_rm); free(Kp1_cm); free(Kp2_cm);
    Rf_error("marker: build_K_from_Z(Z2) failed (code %d).", rc2);
  }

  rm_to_cm_sym(nP1, Kp1_rm, Kp1_cm);
  rm_to_cm_sym(nP2, Kp2_rm, Kp2_cm);

  /* full-data implicit PCG prediction */
  double beta0 = NA_REAL;
  SEXP yhat_obs = PROTECT(allocVector(REALSXP, ny)); nprot++;

  int rc_core = hybrid_marker_fullpcg_core(
    ny, nuse, nobs,
    use, p1_idx, p2_idx,
    a_tr, b_tr, ytr,
    nP1, nP2,
    Kp1_cm, Kp2_cm,
    include_sca,
    s1, s2, s3, se,
    asReal(tol_),
    &beta0,
    REAL(yhat_obs)
  );
  if (rc_core != 0) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    free(Z1rm); free(Z2rm); free(Kp1_rm); free(Kp2_rm); free(Kp1_cm); free(Kp2_cm);
    Rf_error("marker: full-data PCG core failed (rc=%d).", rc_core);
  }

  SEXP beta0_out = PROTECT(allocVector(REALSXP,1)); nprot++;
  REAL(beta0_out)[0] = beta0;

  SEXP varcmp_out = PROTECT(allocVector(REALSXP, kk+1)); nprot++;
  REAL(varcmp_out)[0] = s1;
  REAL(varcmp_out)[1] = s2;
  if (include_sca) REAL(varcmp_out)[2] = s3;
  REAL(varcmp_out)[kk] = se0;

  int ndiag = 12;
  SEXP diag = PROTECT(allocVector(VECSXP, ndiag)); nprot++;
  SEXP dn   = PROTECT(allocVector(STRSXP, ndiag)); nprot++;
  SET_STRING_ELT(dn, 0, mkChar("n_total"));
  SET_STRING_ELT(dn, 1, mkChar("n_matched"));
  SET_STRING_ELT(dn, 2, mkChar("n_train"));
  SET_STRING_ELT(dn, 3, mkChar("n_pred"));
  SET_STRING_ELT(dn, 4, mkChar("kk"));
  SET_STRING_ELT(dn, 5, mkChar("include_sca"));
  SET_STRING_ELT(dn, 6, mkChar("model_type"));
  SET_STRING_ELT(dn, 7, mkChar("vc_source"));
  SET_STRING_ELT(dn, 8, mkChar("boot_n_success"));
  SET_STRING_ELT(dn, 9, mkChar("boot_nRep"));
  SET_STRING_ELT(dn,10, mkChar("boot_nObs"));
  SET_STRING_ELT(dn,11, mkChar("boot_nMarkers"));
  setAttrib(diag, R_NamesSymbol, dn);

  SET_VECTOR_ELT(diag, 0, ScalarInteger(ny));
  SET_VECTOR_ELT(diag, 1, ScalarInteger(nuse));
  SET_VECTOR_ELT(diag, 2, ScalarInteger(nobs));
  SET_VECTOR_ELT(diag, 3, ScalarInteger(npred));
  SET_VECTOR_ELT(diag, 4, ScalarInteger(kk));
  SET_VECTOR_ELT(diag, 5, ScalarLogical(include_sca));
  SET_VECTOR_ELT(diag, 6, mkString("marker"));
  SET_VECTOR_ELT(diag, 7, mkString(vc_source_boot ? "bootstrap" : "provided"));
  SET_VECTOR_ELT(diag, 8, ScalarInteger(boot_success));

  if (vc_source_boot) {
    BootSpec bs; (void)parse_bootstrap(bootstrap_, &bs);
    SET_VECTOR_ELT(diag, 9, ScalarInteger(bs.nRep));
    SET_VECTOR_ELT(diag,10, ScalarInteger(bs.nObs > 0 ? bs.nObs : nobs));
    SET_VECTOR_ELT(diag,11, ScalarInteger(bs.nMarkers > 0 ? bs.nMarkers : ((m1 < m2) ? m1 : m2)));
  } else {
    SET_VECTOR_ELT(diag, 9, ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(diag,10, ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(diag,11, ScalarInteger(NA_INTEGER));
  }

  SEXP ans = PROTECT(allocVector(VECSXP, 6)); nprot++;
  SEXP an  = PROTECT(allocVector(STRSXP, 6)); nprot++;
  SET_STRING_ELT(an, 0, mkChar("beta0"));
  SET_STRING_ELT(an, 1, mkChar("varcmp"));
  SET_STRING_ELT(an, 2, mkChar("yhat_obs"));
  SET_STRING_ELT(an, 3, mkChar("diag"));
  SET_STRING_ELT(an, 4, mkChar("retval"));
  SET_STRING_ELT(an, 5, mkChar("ridge"));
  setAttrib(ans, R_NamesSymbol, an);

  SET_VECTOR_ELT(ans, 0, beta0_out);
  SET_VECTOR_ELT(ans, 1, varcmp_out);
  SET_VECTOR_ELT(ans, 2, yhat_obs);
  SET_VECTOR_ELT(ans, 3, diag);
  SET_VECTOR_ELT(ans, 4, ScalarInteger(0));
  SET_VECTOR_ELT(ans, 5, ScalarReal(ridge));

  free(train_pos);
  free(a_tr); free(b_tr); free(ytr);
  free(Z1rm); free(Z2rm);
  free(Kp1_rm); free(Kp2_rm);
  free(Kp1_cm); free(Kp2_cm);

  UNPROTECT(nprot);
  return ans;
}


# if 0

static SEXP st_mixed_hybrid_marker_impl(
  SEXP y_, SEXP hybrid_id_,
  SEXP pool1_, SEXP pool2_,
  SEXP Z1_, SEXP Z1_ids_,
  SEXP Z2_, SEXP Z2_ids_,
  SEXP include_sca_,
  SEXP method_, SEXP coding_, SEXP varcmp_,
  SEXP maxcyc_, SEXP tol_, SEXP ridge_,
  SEXP bootstrap_)
{
  int nprot = 0;

  if (!isReal(y_) || !isVector(y_)) Rf_error("marker: y must be numeric vector");
  const int ny = (int)XLENGTH(y_);

  if (!isString(hybrid_id_) || (int)XLENGTH(hybrid_id_) != ny) Rf_error("marker: hybrid_id must be character length(y)");
  if (!isString(pool1_) || (int)XLENGTH(pool1_) != ny) Rf_error("marker: pool1 must be character length(y)");
  if (!isString(pool2_) || (int)XLENGTH(pool2_) != ny) Rf_error("marker: pool2 must be character length(y)");

  if (!isReal(Z1_) || !isMatrix(Z1_)) Rf_error("marker: Z1 must be numeric matrix");
  if (!isReal(Z2_) || !isMatrix(Z2_)) Rf_error("marker: Z2 must be numeric matrix");
  if (!isString(Z1_ids_) || XLENGTH(Z1_ids_) < 1) Rf_error("marker: Z1_ids must be character");
  if (!isString(Z2_ids_) || XLENGTH(Z2_ids_) < 1) Rf_error("marker: Z2_ids must be character");

  if (!isString(method_) || XLENGTH(method_) != 1) Rf_error("marker: method must be single string");
  if (!isString(coding_) || XLENGTH(coding_) != 1) Rf_error("marker: coding must be single string");

  const int include_sca = asLogical(include_sca_) ? 1 : 0;
  const double ridge = asReal(ridge_);
  if (!(ridge >= 0.0) || !(ridge == ridge)) Rf_error("marker: ridge must be >= 0");

  SEXP d1 = getAttrib(Z1_, R_DimSymbol);
  SEXP d2 = getAttrib(Z2_, R_DimSymbol);
  const int nP1 = INTEGER(d1)[0], m1 = INTEGER(d1)[1];
  const int nP2 = INTEGER(d2)[0], m2 = INTEGER(d2)[1];
  if ((int)XLENGTH(Z1_ids_) != nP1) Rf_error("marker: length(Z1_ids) must equal nrow(Z1)");
  if ((int)XLENGTH(Z2_ids_) != nP2) Rf_error("marker: length(Z2_ids) must equal nrow(Z2)");

  KMethod km = K_STD;
  if (str_eq(method_, "std")) km = K_STD;
  else if (str_eq(method_, "vanraden")) km = K_VANRADEN;
  else if (str_eq(method_, "sommer")) km = K_SOMMER;
  else Rf_error("marker: method must be one of: 'std','vanraden','sommer'");

  const int kk = include_sca ? 3 : 2;

  /* match rows to parents */
  int *p1_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *p2_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *use    = (int*)R_alloc((size_t)ny, sizeof(int));
  int nuse = 0;

  for (int i = 0; i < ny; ++i) {
    const char *a = CHAR(STRING_ELT(pool1_, i));
    const char *b = CHAR(STRING_ELT(pool2_, i));
    int ia = find_id_index(Z1_ids_, a);
    int ib = find_id_index(Z2_ids_, b);
    p1_idx[i] = ia;
    p2_idx[i] = ib;
    if (ia >= 0 && ib >= 0) use[nuse++] = i;
  }
  if (nuse <= 0) Rf_error("marker: No rows matched to both pool1 and pool2 IDs.");

  int *train_pos = (int*)malloc((size_t)nuse * sizeof(int));
  if (!train_pos) Rf_error("marker: OOM train_pos");

  int nobs = 0, npred = 0;
  for (int t = 0; t < nuse; ++t) {
    int i = use[t];
    if (!is_missing(REAL(y_)[i])) train_pos[nobs++] = t;
    else npred++;
  }
  if (nobs <= 0) {
    free(train_pos);
    Rf_error("marker: No non-missing y values after matching parents.");
  }

  int *a_tr = (int*)malloc((size_t)nobs * sizeof(int));
  int *b_tr = (int*)malloc((size_t)nobs * sizeof(int));
  double *ytr = (double*)malloc((size_t)nobs * sizeof(double));
  if (!a_tr || !b_tr || !ytr) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: OOM a_tr/b_tr/ytr");
  }

  for (int i = 0; i < nobs; ++i) {
    int obs = use[ train_pos[i] ];
    a_tr[i] = p1_idx[obs];
    b_tr[i] = p2_idx[obs];
    ytr[i]  = REAL(y_)[obs];
  }

  /* variance components: provided or bootstrap-based */
  double s1=NAN, s2=NAN, s3=0.0, se0=NAN;
  int vc_source_boot = 0;
  int boot_success = 0;

  if (!isNull(varcmp_)) {
    if (!isReal(varcmp_) || (int)XLENGTH(varcmp_) != (kk + 1)) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp must be numeric length (kk+1)=%d", kk+1);
    }
    s1  = REAL(varcmp_)[0];
    s2  = REAL(varcmp_)[1];
    s3  = include_sca ? REAL(varcmp_)[2] : 0.0;
    se0 = REAL(varcmp_)[kk];
  } else {
    if (!include_sca) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp NULL; bootstrap VC only enabled when include_sca=TRUE.");
    }

    BootSpec bs;
    int pr = parse_bootstrap(bootstrap_, &bs);
    if (pr != 0 || !bs.enabled) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp NULL and bootstrap not set correctly.");
    }

    double vc_est[4] = {0,0,0,0};
    int rc = bootstrap_estimate_varcmp_hybrid(
      nobs, ytr, a_tr, b_tr,
      nP1, m1, nP2, m2,
      REAL(Z1_), REAL(Z2_),
      km, coding_, 1,
      asInteger(maxcyc_), asReal(tol_), ridge,
      &bs,
      vc_est, &boot_success
    );
    if (rc != 0) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: bootstrap VC failed (rc=%d).", rc);
    }

    s1 = vc_est[0];
    s2 = vc_est[1];
    s3 = vc_est[2];
    se0= vc_est[3];
    vc_source_boot = 1;
  }

  if (!(s1 >= 0.0) || !(s2 >= 0.0) || !(s3 >= 0.0) || !(se0 > 0.0) ||
      !(s1==s1) || !(s2==s2) || !(s3==s3) || !(se0==se0)) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: invalid variance components.");
  }

  const double se = se0 + ridge;

  /* parental kernels from full markers */
  double *Z1rm = (double*)malloc((size_t)nP1*(size_t)m1*sizeof(double));
  double *Z2rm = (double*)malloc((size_t)nP2*(size_t)m2*sizeof(double));
  double *Kp1_rm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_rm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  double *Kp1_cm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_cm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));

  if (!Z1rm||!Z2rm||!Kp1_rm||!Kp2_rm||!Kp1_cm||!Kp2_cm) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    free(Z1rm); free(Z2rm); free(Kp1_rm); free(Kp2_rm); free(Kp1_cm); free(Kp2_cm);
    Rf_error("marker: OOM kernels");
  }

  Z_cm_to_rm(REAL(Z1_), nP1, m1, Z1rm);
  Z_cm_to_rm(REAL(Z2_), nP2, m2, Z2rm);

  double zmean1=0,zsd1=0,zmean2=0,zsd2=0;
  ZCoding cu1=CODING_OTHER, cu2=CODING_OTHER;
  ZTransform tr1=TRANS_NONE, tr2=TRANS_NONE;

  int rc1 = build_K_from_Z(Z1rm, nP1, m1, km, coding_, Kp1_rm, &zmean1, &zsd1, &cu1, &tr1);
  if (rc1 != 0) Rf_error("marker: build_K_from_Z(Z1) failed (code %d).", rc1);

  int rc2 = build_K_from_Z(Z2rm, nP2, m2, km, coding_, Kp2_rm, &zmean2, &zsd2, &cu2, &tr2);
  if (rc2 != 0) Rf_error("marker: build_K_from_Z(Z2) failed (code %d).", rc2);

  rm_to_cm_sym(nP1, Kp1_rm, Kp1_cm);
  rm_to_cm_sym(nP2, Kp2_rm, Kp2_cm);

  /* full-data implicit Vv + PCG prediction */
  double beta0 = NA_REAL;
  SEXP yhat_obs = PROTECT(allocVector(REALSXP, ny)); nprot++;

  int rc_core = hybrid_marker_fullpcg_core(
    ny, nuse, nobs,
    use, p1_idx, p2_idx,
    a_tr, b_tr, ytr,
    nP1, nP2,
    Kp1_cm, Kp2_cm,
    include_sca,
    s1, s2, s3, se,
    asReal(tol_),
    &beta0,
    REAL(yhat_obs)
  );
  if (rc_core != 0) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    free(Z1rm); free(Z2rm); free(Kp1_rm); free(Kp2_rm); free(Kp1_cm); free(Kp2_cm);
    Rf_error("marker: full-data PCG core failed (rc=%d).", rc_core);
  }

  SEXP beta0_out = PROTECT(allocVector(REALSXP,1)); nprot++;
  REAL(beta0_out)[0] = beta0;

  SEXP varcmp_out = PROTECT(allocVector(REALSXP, kk+1)); nprot++;
  REAL(varcmp_out)[0] = s1;
  REAL(varcmp_out)[1] = s2;
  if (include_sca) REAL(varcmp_out)[2] = s3;
  REAL(varcmp_out)[kk] = se0;

  int ndiag = 12;
  SEXP diag = PROTECT(allocVector(VECSXP, ndiag)); nprot++;
  SEXP dn   = PROTECT(allocVector(STRSXP, ndiag)); nprot++;
  SET_STRING_ELT(dn, 0, mkChar("n_total"));
  SET_STRING_ELT(dn, 1, mkChar("n_matched"));
  SET_STRING_ELT(dn, 2, mkChar("n_train"));
  SET_STRING_ELT(dn, 3, mkChar("n_pred"));
  SET_STRING_ELT(dn, 4, mkChar("kk"));
  SET_STRING_ELT(dn, 5, mkChar("include_sca"));
  SET_STRING_ELT(dn, 6, mkChar("model_type"));
  SET_STRING_ELT(dn, 7, mkChar("vc_source"));
  SET_STRING_ELT(dn, 8, mkChar("boot_n_success"));
  SET_STRING_ELT(dn, 9, mkChar("boot_nRep"));
  SET_STRING_ELT(dn,10, mkChar("boot_nObs"));
  SET_STRING_ELT(dn,11, mkChar("boot_nMarkers"));
  setAttrib(diag, R_NamesSymbol, dn);

  SET_VECTOR_ELT(diag, 0, ScalarInteger(ny));
  SET_VECTOR_ELT(diag, 1, ScalarInteger(nuse));
  SET_VECTOR_ELT(diag, 2, ScalarInteger(nobs));
  SET_VECTOR_ELT(diag, 3, ScalarInteger(npred));
  SET_VECTOR_ELT(diag, 4, ScalarInteger(kk));
  SET_VECTOR_ELT(diag, 5, ScalarLogical(include_sca));
  SET_VECTOR_ELT(diag, 6, mkString("marker"));
  SET_VECTOR_ELT(diag, 7, mkString(vc_source_boot ? "bootstrap" : "provided"));
  SET_VECTOR_ELT(diag, 8, ScalarInteger(boot_success));

  if (vc_source_boot) {
    BootSpec bs; (void)parse_bootstrap(bootstrap_, &bs);
    SET_VECTOR_ELT(diag, 9, ScalarInteger(bs.nRep));
    SET_VECTOR_ELT(diag,10, ScalarInteger(bs.nObs));
    SET_VECTOR_ELT(diag,11, ScalarInteger(bs.nMarkers));
  } else {
    SET_VECTOR_ELT(diag, 9, ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(diag,10, ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(diag,11, ScalarInteger(NA_INTEGER));
  }

  SEXP ans = PROTECT(allocVector(VECSXP, 6)); nprot++;
  SEXP an  = PROTECT(allocVector(STRSXP, 6)); nprot++;
  SET_STRING_ELT(an, 0, mkChar("beta0"));
  SET_STRING_ELT(an, 1, mkChar("varcmp"));
  SET_STRING_ELT(an, 2, mkChar("yhat_obs"));
  SET_STRING_ELT(an, 3, mkChar("diag"));
  SET_STRING_ELT(an, 4, mkChar("retval"));
  SET_STRING_ELT(an, 5, mkChar("ridge"));
  setAttrib(ans, R_NamesSymbol, an);

  SET_VECTOR_ELT(ans, 0, beta0_out);
  SET_VECTOR_ELT(ans, 1, varcmp_out);
  SET_VECTOR_ELT(ans, 2, yhat_obs);
  SET_VECTOR_ELT(ans, 3, diag);
  SET_VECTOR_ELT(ans, 4, ScalarInteger(0));
  SET_VECTOR_ELT(ans, 5, ScalarReal(ridge));

  free(train_pos); free(a_tr); free(b_tr); free(ytr);
  free(Z1rm); free(Z2rm);
  free(Kp1_rm); free(Kp2_rm);
  free(Kp1_cm); free(Kp2_cm);

  UNPROTECT(nprot);
  return ans;
}



static SEXP st_mixed_hybrid_marker_impl(
  SEXP y_, SEXP hybrid_id_,
  SEXP pool1_, SEXP pool2_,
  SEXP Z1_, SEXP Z1_ids_,
  SEXP Z2_, SEXP Z2_ids_,
  SEXP include_sca_,
  SEXP method_, SEXP coding_, SEXP varcmp_,
  SEXP maxcyc_, SEXP tol_, SEXP ridge_,
  SEXP bootstrap_)
{
  int nprot = 0;

  if (!isReal(y_) || !isVector(y_)) Rf_error("marker: y must be numeric vector");  
  const int ny = (int)XLENGTH(y_);

  if (!isString(hybrid_id_) || (int)XLENGTH(hybrid_id_) != ny) Rf_error("marker: hybrid_id must be character length(y)");
  if (!isString(pool1_) || (int)XLENGTH(pool1_) != ny) Rf_error("marker: pool1 must be character length(y)");
  if (!isString(pool2_) || (int)XLENGTH(pool2_) != ny) Rf_error("marker: pool2 must be character length(y)");

  if (!isReal(Z1_) || !isMatrix(Z1_)) Rf_error("marker: Z1 must be numeric matrix");
  if (!isReal(Z2_) || !isMatrix(Z2_)) Rf_error("marker: Z2 must be numeric matrix");
  if (!isString(Z1_ids_) || XLENGTH(Z1_ids_) < 1) Rf_error("marker: Z1_ids must be character");
  if (!isString(Z2_ids_) || XLENGTH(Z2_ids_) < 1) Rf_error("marker: Z2_ids must be character");

  if (!isString(method_) || XLENGTH(method_) != 1) Rf_error("marker: method must be single string");
  if (!isString(coding_) || XLENGTH(coding_) != 1) Rf_error("marker: coding must be single string");

  const int include_sca = asLogical(include_sca_) ? 1 : 0;
  const double ridge = asReal(ridge_);
  if (!(ridge >= 0.0) || !(ridge == ridge)) Rf_error("marker: ridge must be >= 0");

  SEXP d1 = getAttrib(Z1_, R_DimSymbol);
  SEXP d2 = getAttrib(Z2_, R_DimSymbol);
  const int nP1 = INTEGER(d1)[0], m1 = INTEGER(d1)[1];
  const int nP2 = INTEGER(d2)[0], m2 = INTEGER(d2)[1];
  if ((int)XLENGTH(Z1_ids_) != nP1) Rf_error("marker: length(Z1_ids) must equal nrow(Z1)");
  if ((int)XLENGTH(Z2_ids_) != nP2) Rf_error("marker: length(Z2_ids) must equal nrow(Z2)");

  KMethod km = K_STD;
  if (str_eq(method_, "std")) km = K_STD;
  else if (str_eq(method_, "vanraden")) km = K_VANRADEN;
  else if (str_eq(method_, "sommer")) km = K_SOMMER;
  else Rf_error("marker: method must be one of: 'std','vanraden','sommer'");

  const int kk = include_sca ? 3 : 2;

  /* match to parents */
  int *p1_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *p2_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *use    = (int*)R_alloc((size_t)ny, sizeof(int));
  int nuse = 0;

  for (int i = 0; i < ny; ++i) {
    const char *a = CHAR(STRING_ELT(pool1_, i));
    const char *b = CHAR(STRING_ELT(pool2_, i));
    int ia = find_id_index(Z1_ids_, a);
    int ib = find_id_index(Z2_ids_, b);
    p1_idx[i] = ia;
    p2_idx[i] = ib;
    if (ia >= 0 && ib >= 0) use[nuse++] = i;
  }
  if (nuse <= 0) Rf_error("marker: No rows matched to both pool1 and pool2 IDs.");

  int *train_pos = (int*)malloc((size_t)nuse * sizeof(int));
  if (!train_pos) Rf_error("marker: OOM train_pos");
  int nobs = 0, npred = 0;
  for (int t = 0; t < nuse; ++t) {
    int i = use[t];
    if (!is_missing(REAL(y_)[i])) train_pos[nobs++] = t;
    else npred++;
  }
  if (nobs <= 0) { free(train_pos); Rf_error("marker: No non-missing y values after matching parents."); }

  int *a_tr = (int*)malloc((size_t)nobs * sizeof(int));
  int *b_tr = (int*)malloc((size_t)nobs * sizeof(int));
  double *ytr = (double*)malloc((size_t)nobs * sizeof(double));
  if (!a_tr || !b_tr || !ytr) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: OOM a_tr/b_tr/ytr");
  }
  for (int i = 0; i < nobs; ++i) {
    int obs = use[ train_pos[i] ];
    a_tr[i] = p1_idx[obs];
    b_tr[i] = p2_idx[obs];
    ytr[i]  = REAL(y_)[obs];
  }

  /* variance components: provided or bootstrap */
  double s1=NAN, s2=NAN, s3=0.0, se0=NAN;
  int vc_source_boot = 0;
  int boot_success = 0;

  if (!isNull(varcmp_)) {
    if (!isReal(varcmp_) || (int)XLENGTH(varcmp_) != (kk + 1)) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp must be numeric length (kk+1)=%d", kk+1);
    }
    s1  = REAL(varcmp_)[0];
    s2  = REAL(varcmp_)[1];
    s3  = include_sca ? REAL(varcmp_)[2] : 0.0;
    se0 = REAL(varcmp_)[kk];
  } else {
    if (!include_sca) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp NULL; bootstrap VC only enabled when include_sca=TRUE.");
    }
    BootSpec bs;
    int pr = parse_bootstrap(bootstrap_, &bs);
    if (pr != 0 || !bs.enabled) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp NULL and bootstrap not set (need nObs and nMarkers).");
    }
    if (bs.nObs > nobs) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: bootstrap$nObs (%d) must be <= n_train (%d).", bs.nObs, nobs);
    }
    if (bs.nMarkers > m1 || bs.nMarkers > m2) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: bootstrap$nMarkers must be <= ncol(Z1) and <= ncol(Z2).");
    }

    double vc_est[4] = {0,0,0,0};
    int rc = bootstrap_estimate_varcmp_hybrid(
      nobs, ytr, a_tr, b_tr,
      nP1, m1, nP2, m2,
      REAL(Z1_), REAL(Z2_),
      km, coding_, 1,
      asInteger(maxcyc_), asReal(tol_), ridge,
      &bs,
      vc_est, &boot_success
    );
    if (rc != 0) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: bootstrap VC failed (rc=%d).", rc);
    }
    s1 = vc_est[0];
    s2 = vc_est[1];
    s3 = vc_est[2];
    se0= vc_est[3];
    vc_source_boot = 1;
  }

  if (!(s1 >= 0.0) || !(s2 >= 0.0) || !(s3 >= 0.0) || !(se0 > 0.0) ||
      !(s1==s1) || !(s2==s2) || !(s3==s3) || !(se0==se0)) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: invalid variance components (need s1,s2,s3>=0 and se>0).");
  }

  const double se = se0 + ridge;
  if (!(se > 0.0) || !(se==se)) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: residual variance + ridge must be > 0.");
  }

  /* parental kernels from FULL markers */
  double *Z1rm = (double*)malloc((size_t)nP1*(size_t)m1*sizeof(double));
  double *Z2rm = (double*)malloc((size_t)nP2*(size_t)m2*sizeof(double));
  double *Kp1_rm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_rm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  double *Kp1_cm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_cm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  if (!Z1rm||!Z2rm||!Kp1_rm||!Kp2_rm||!Kp1_cm||!Kp2_cm) Rf_error("marker: OOM kernels");

  Z_cm_to_rm(REAL(Z1_), nP1, m1, Z1rm);
  Z_cm_to_rm(REAL(Z2_), nP2, m2, Z2rm);

  double zmean1=0,zsd1=0,zmean2=0,zsd2=0;
  ZCoding cu1=CODING_OTHER, cu2=CODING_OTHER;
  ZTransform tr1=TRANS_NONE, tr2=TRANS_NONE;

  int rc1 = build_K_from_Z(Z1rm, nP1, m1, km, coding_, Kp1_rm, &zmean1, &zsd1, &cu1, &tr1);
  if (rc1 != 0) Rf_error("marker: build_K_from_Z(Z1) failed (code %d).", rc1);
  int rc2 = build_K_from_Z(Z2rm, nP2, m2, km, coding_, Kp2_rm, &zmean2, &zsd2, &cu2, &tr2);
  if (rc2 != 0) Rf_error("marker: build_K_from_Z(Z2) failed (code %d).", rc2);

  rm_to_cm_sym(nP1, Kp1_rm, Kp1_cm);
  rm_to_cm_sym(nP2, Kp2_rm, Kp2_cm);

  /* Sigma operator */
  SigmaOp op;
  op.nobs = nobs; op.nP1 = nP1; op.nP2 = nP2;
  op.a = a_tr; op.b = b_tr;
  op.Kp1_cm = Kp1_cm; op.Kp2_cm = Kp2_cm;
  op.include_sca = include_sca;
  op.s1 = s1; op.s2 = s2; op.s3 = include_sca ? s3 : 0.0; op.se = se;
  op.tmp1 = (double*)malloc((size_t)nobs*sizeof(double));
  op.tmp2 = (double*)malloc((size_t)nobs*sizeof(double));
  op.tmp3 = include_sca ? (double*)malloc((size_t)nobs*sizeof(double)) : NULL;
  if (!op.tmp1 || !op.tmp2 || (include_sca && !op.tmp3)) Rf_error("marker: OOM SigmaOp buffers");

  double tol_cg = asReal(tol_);
  if (!(tol_cg > 0.0) || !(tol_cg==tol_cg)) tol_cg = 1e-6;
  int maxit = 400;

  /* solve Wy and W1 */
  double *Wy = (double*)calloc((size_t)nobs, sizeof(double));
  double *W1 = (double*)calloc((size_t)nobs, sizeof(double));
  double *ones = (double*)malloc((size_t)nobs*sizeof(double));
  if (!Wy || !W1 || !ones) Rf_error("marker: OOM Wy/W1/ones");
  for (int i=0;i<nobs;++i) ones[i]=1.0;

  {
    int rc = cg_solve_sigma(&op, ytr, Wy, maxit, tol_cg);
    if (rc != 0) Rf_error("marker: CG failed for Wy (rc=%d)", rc);
  }
  {
    int rc = cg_solve_sigma(&op, ones, W1, maxit, tol_cg);
    if (rc != 0) Rf_error("marker: CG failed for W1 (rc=%d)", rc);
  }

  double num=0.0, den=0.0;
  for (int i=0;i<nobs;++i) { num += Wy[i]; den += W1[i]; }
  if (!(den > 0.0) || !(den==den)) Rf_error("marker: invalid beta0 denominator");
  const double beta0 = num / den;

  double *r = (double*)malloc((size_t)nobs*sizeof(double));
  double *alpha = (double*)calloc((size_t)nobs, sizeof(double));
  if (!r || !alpha) Rf_error("marker: OOM r/alpha");
  for (int i=0;i<nobs;++i) r[i] = ytr[i] - beta0;
  {
    int rc = cg_solve_sigma(&op, r, alpha, maxit, tol_cg);
    if (rc != 0) Rf_error("marker: CG failed for alpha (rc=%d)", rc);
  }

  /* predict for all matched hybrids */
  SEXP yhat_obs = PROTECT(allocVector(REALSXP, ny)); nprot++;
  for (int i=0;i<ny;++i) REAL(yhat_obs)[i] = NA_REAL;

  /* precompute tA = Kp1 * sum_a(alpha), tB = Kp2 * sum_b(alpha) */
  double *sA = (double*)calloc((size_t)nP1, sizeof(double));
  double *sB = (double*)calloc((size_t)nP2, sizeof(double));
  if (!sA || !sB) Rf_error("marker: OOM sA/sB");
  for (int j=0;j<nobs;++j) { sA[a_tr[j]] += alpha[j]; sB[b_tr[j]] += alpha[j]; }

  double *tA = (double*)malloc((size_t)nP1*sizeof(double));
  double *tB = (double*)malloc((size_t)nP2*sizeof(double));
  if (!tA || !tB) Rf_error("marker: OOM tA/tB");

  {
    char uplo='U';
    BLAS_INT INC=1;
    const double one=1.0, zero=0.0;
    BLAS_INT N1=(BLAS_INT)nP1, N2=(BLAS_INT)nP2;
    F77_CALL(dsymv)(&uplo, &N1, &one, Kp1_cm, &N1, sA, &INC, &zero, tA, &INC FCONE);
    F77_CALL(dsymv)(&uplo, &N2, &one, Kp2_cm, &N2, sB, &INC, &zero, tB, &INC FCONE);
  }

  /* SCA precompute Tsca = Kp1 * Aalpha * Kp2^T */
  double *Tsca = NULL;
  if (include_sca) {
    const size_t R = (size_t)nP1*(size_t)nP2;
    double *Aalpha = (double*)calloc(R, sizeof(double)); /* col-major nP1 x nP2 */
    double *W = (double*)malloc(R * sizeof(double));
    Tsca = (double*)malloc(R * sizeof(double));
    if (!Aalpha || !W || !Tsca) Rf_error("marker: OOM SCA matrices");

    for (int j=0;j<nobs;++j) {
      Aalpha[a_tr[j] + (size_t)nP1*(size_t)b_tr[j]] += alpha[j];
    }

    {
      char transN='N';
      BLAS_INT M=(BLAS_INT)nP1, N=(BLAS_INT)nP2, K=(BLAS_INT)nP1;
      BLAS_INT lda=M, ldb=M, ldc=M;
      const double one=1.0, zero=0.0;
      F77_CALL(dgemm)(&transN,&transN,&M,&N,&K,&one, Kp1_cm,&lda, Aalpha,&ldb, &zero, W,&ldc FCONE FCONE);
    }
    {
      char transN='N', transT='T';
      BLAS_INT M=(BLAS_INT)nP1, N=(BLAS_INT)nP2, K=(BLAS_INT)nP2;
      BLAS_INT lda=(BLAS_INT)nP1, ldb=(BLAS_INT)nP2, ldc=(BLAS_INT)nP1;
      const double one=1.0, zero=0.0;
      F77_CALL(dgemm)(&transN,&transT,&M,&N,&K,&one, W,&lda, Kp2_cm,&ldb, &zero, Tsca,&ldc FCONE FCONE);
    }

    free(Aalpha);
    free(W);
  }

  for (int t=0; t<nuse; ++t) {
    int obs = use[t];
    int a0 = p1_idx[obs];
    int b0 = p2_idx[obs];
    double pred = beta0;
    pred += s1 * tA[a0];
    pred += s2 * tB[b0];
    if (include_sca) pred += s3 * Tsca[a0 + (size_t)nP1*(size_t)b0];
    REAL(yhat_obs)[obs] = pred;
  }

  /* outputs */
  SEXP beta0_out = PROTECT(allocVector(REALSXP,1)); nprot++;
  REAL(beta0_out)[0] = beta0;

  SEXP varcmp_out = PROTECT(allocVector(REALSXP, kk+1)); nprot++;
  REAL(varcmp_out)[0] = s1;
  REAL(varcmp_out)[1] = s2;
  if (include_sca) REAL(varcmp_out)[2] = s3;
  REAL(varcmp_out)[kk] = se0;

  int ndiag = 12;
  SEXP diag = PROTECT(allocVector(VECSXP, ndiag)); nprot++;
  SEXP dn   = PROTECT(allocVector(STRSXP, ndiag)); nprot++;
  SET_STRING_ELT(dn, 0, mkChar("n_total"));
  SET_STRING_ELT(dn, 1, mkChar("n_matched"));
  SET_STRING_ELT(dn, 2, mkChar("n_train"));
  SET_STRING_ELT(dn, 3, mkChar("n_pred"));
  SET_STRING_ELT(dn, 4, mkChar("kk"));
  SET_STRING_ELT(dn, 5, mkChar("include_sca"));
  SET_STRING_ELT(dn, 6, mkChar("model_type"));
  SET_STRING_ELT(dn, 7, mkChar("vc_source"));
  SET_STRING_ELT(dn, 8, mkChar("boot_n_success"));
  SET_STRING_ELT(dn, 9, mkChar("boot_nRep"));
  SET_STRING_ELT(dn,10, mkChar("boot_nObs"));
  SET_STRING_ELT(dn,11, mkChar("boot_nMarkers"));
  setAttrib(diag, R_NamesSymbol, dn);

  SET_VECTOR_ELT(diag, 0, ScalarInteger(ny));
  SET_VECTOR_ELT(diag, 1, ScalarInteger(nuse));
  SET_VECTOR_ELT(diag, 2, ScalarInteger(nobs));
  SET_VECTOR_ELT(diag, 3, ScalarInteger(npred));
  SET_VECTOR_ELT(diag, 4, ScalarInteger(kk));
  SET_VECTOR_ELT(diag, 5, ScalarLogical(include_sca));
  SET_VECTOR_ELT(diag, 6, mkString("marker"));
  SET_VECTOR_ELT(diag, 7, mkString(vc_source_boot ? "bootstrap" : "provided"));
  SET_VECTOR_ELT(diag, 8, ScalarInteger(boot_success));
  if (vc_source_boot) {
    BootSpec bs; (void)parse_bootstrap(bootstrap_, &bs);
    SET_VECTOR_ELT(diag, 9, ScalarInteger(bs.nRep));
    SET_VECTOR_ELT(diag,10, ScalarInteger(bs.nObs));
    SET_VECTOR_ELT(diag,11, ScalarInteger(bs.nMarkers));
  } else {
    SET_VECTOR_ELT(diag, 9, ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(diag,10, ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(diag,11, ScalarInteger(NA_INTEGER));
  }

  /* -------- DEBUG payload -------- */
  SEXP debug = PROTECT(allocVector(VECSXP, 9)); nprot++;
  SEXP dbg_n = PROTECT(allocVector(STRSXP, 9)); nprot++;
  SET_STRING_ELT(dbg_n, 0, mkChar("beta0"));
  SET_STRING_ELT(dbg_n, 1, mkChar("a_tr"));
  SET_STRING_ELT(dbg_n, 2, mkChar("b_tr"));
  SET_STRING_ELT(dbg_n, 3, mkChar("alpha"));
  SET_STRING_ELT(dbg_n, 4, mkChar("sA"));
  SET_STRING_ELT(dbg_n, 5, mkChar("sB"));
  SET_STRING_ELT(dbg_n, 6, mkChar("tA"));
  SET_STRING_ELT(dbg_n, 7, mkChar("tB"));
  SET_STRING_ELT(dbg_n, 8, mkChar("p_idx_all"));
  setAttrib(debug, R_NamesSymbol, dbg_n);

  SEXP dbg_beta0 = PROTECT(allocVector(REALSXP, 1)); nprot++;
  REAL(dbg_beta0)[0] = beta0;
  SET_VECTOR_ELT(debug, 0, dbg_beta0);

  SEXP dbg_a_tr = PROTECT(allocVector(INTSXP, nobs)); nprot++;
  SEXP dbg_b_tr = PROTECT(allocVector(INTSXP, nobs)); nprot++;
  for (int i=0;i<nobs;++i) {
    INTEGER(dbg_a_tr)[i] = a_tr[i] + 1; /* back to R 1-based */
    INTEGER(dbg_b_tr)[i] = b_tr[i] + 1;
  }
  SET_VECTOR_ELT(debug, 1, dbg_a_tr);
  SET_VECTOR_ELT(debug, 2, dbg_b_tr);

  SEXP dbg_alpha = PROTECT(allocVector(REALSXP, nobs)); nprot++;
  for (int i=0;i<nobs;++i) REAL(dbg_alpha)[i] = alpha[i];
  SET_VECTOR_ELT(debug, 3, dbg_alpha);

  SEXP dbg_sA = PROTECT(allocVector(REALSXP, nP1)); nprot++;
  SEXP dbg_sB = PROTECT(allocVector(REALSXP, nP2)); nprot++;
  for (int i=0;i<nP1;++i) REAL(dbg_sA)[i] = sA[i];
  for (int i=0;i<nP2;++i) REAL(dbg_sB)[i] = sB[i];
  SET_VECTOR_ELT(debug, 4, dbg_sA);
  SET_VECTOR_ELT(debug, 5, dbg_sB);

  SEXP dbg_tA = PROTECT(allocVector(REALSXP, nP1)); nprot++;
  SEXP dbg_tB = PROTECT(allocVector(REALSXP, nP2)); nprot++;
  for (int i=0;i<nP1;++i) REAL(dbg_tA)[i] = tA[i];
  for (int i=0;i<nP2;++i) REAL(dbg_tB)[i] = tB[i];
  SET_VECTOR_ELT(debug, 6, dbg_tA);
  SET_VECTOR_ELT(debug, 7, dbg_tB);

  SEXP dbg_pidx = PROTECT(allocVector(VECSXP, 2)); nprot++;
  SEXP dbg_pidx_n = PROTECT(allocVector(STRSXP, 2)); nprot++;
  SET_STRING_ELT(dbg_pidx_n, 0, mkChar("p1_idx_all"));
  SET_STRING_ELT(dbg_pidx_n, 1, mkChar("p2_idx_all"));
  setAttrib(dbg_pidx, R_NamesSymbol, dbg_pidx_n);

  SEXP dbg_p1_all = PROTECT(allocVector(INTSXP, ny)); nprot++;
  SEXP dbg_p2_all = PROTECT(allocVector(INTSXP, ny)); nprot++;
  for (int i=0;i<ny;++i) {
    INTEGER(dbg_p1_all)[i] = (p1_idx[i] >= 0) ? (p1_idx[i] + 1) : NA_INTEGER;
    INTEGER(dbg_p2_all)[i] = (p2_idx[i] >= 0) ? (p2_idx[i] + 1) : NA_INTEGER;
  }
  SET_VECTOR_ELT(dbg_pidx, 0, dbg_p1_all);
  SET_VECTOR_ELT(dbg_pidx, 1, dbg_p2_all);
  SET_VECTOR_ELT(debug, 8, dbg_pidx);

  SEXP ans = PROTECT(allocVector(VECSXP, 7)); nprot++;
  SEXP an  = PROTECT(allocVector(STRSXP, 7)); nprot++;
  SET_STRING_ELT(an, 0, mkChar("beta0"));
  SET_STRING_ELT(an, 1, mkChar("varcmp"));
  SET_STRING_ELT(an, 2, mkChar("yhat_obs"));
  SET_STRING_ELT(an, 3, mkChar("diag"));
  SET_STRING_ELT(an, 4, mkChar("retval"));
  SET_STRING_ELT(an, 5, mkChar("ridge"));
  SET_STRING_ELT(an, 6, mkChar("debug"));
  setAttrib(ans, R_NamesSymbol, an);

  SET_VECTOR_ELT(ans, 0, beta0_out);
  SET_VECTOR_ELT(ans, 1, varcmp_out);
  SET_VECTOR_ELT(ans, 2, yhat_obs);
  SET_VECTOR_ELT(ans, 3, diag);
  SET_VECTOR_ELT(ans, 4, ScalarInteger(0));
  SET_VECTOR_ELT(ans, 5, ScalarReal(ridge));
  SET_VECTOR_ELT(ans, 6, debug);

  free(train_pos);
  free(a_tr); free(b_tr); free(ytr);
  free(Z1rm); free(Z2rm);
  free(Kp1_rm); free(Kp2_rm);
  free(Kp1_cm); free(Kp2_cm);
  free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);
  free(Wy); free(W1); free(ones);
  free(r); free(alpha);
  free(sA); free(sB);
  free(tA); free(tB);
  if (Tsca) free(Tsca);

  UNPROTECT(nprot);
  return ans;
}


static SEXP st_mixed_hybrid_marker_impl(
  SEXP y_, SEXP hybrid_id_,
  SEXP pool1_, SEXP pool2_,
  SEXP Z1_, SEXP Z1_ids_,
  SEXP Z2_, SEXP Z2_ids_,
  SEXP include_sca_,
  SEXP method_, SEXP coding_, SEXP varcmp_,
  SEXP maxcyc_, SEXP tol_, SEXP ridge_,
  SEXP bootstrap_)
{
  int nprot = 0;

  if (!isReal(y_) || !isVector(y_)) Rf_error("marker: y must be numeric vector");
  const int ny = (int)XLENGTH(y_);

  if (!isString(hybrid_id_) || (int)XLENGTH(hybrid_id_) != ny) Rf_error("marker: hybrid_id must be character length(y)");
  if (!isString(pool1_) || (int)XLENGTH(pool1_) != ny) Rf_error("marker: pool1 must be character length(y)");
  if (!isString(pool2_) || (int)XLENGTH(pool2_) != ny) Rf_error("marker: pool2 must be character length(y)");

  if (!isReal(Z1_) || !isMatrix(Z1_)) Rf_error("marker: Z1 must be numeric matrix");
  if (!isReal(Z2_) || !isMatrix(Z2_)) Rf_error("marker: Z2 must be numeric matrix");
  if (!isString(Z1_ids_) || XLENGTH(Z1_ids_) < 1) Rf_error("marker: Z1_ids must be character");
  if (!isString(Z2_ids_) || XLENGTH(Z2_ids_) < 1) Rf_error("marker: Z2_ids must be character");

  if (!isString(method_) || XLENGTH(method_) != 1) Rf_error("marker: method must be single string");
  if (!isString(coding_) || XLENGTH(coding_) != 1) Rf_error("marker: coding must be single string");

  const int include_sca = asLogical(include_sca_) ? 1 : 0;
  const double ridge = asReal(ridge_);
  if (!(ridge >= 0.0) || !(ridge == ridge)) Rf_error("marker: ridge must be >= 0");

  SEXP d1 = getAttrib(Z1_, R_DimSymbol);
  SEXP d2 = getAttrib(Z2_, R_DimSymbol);
  const int nP1 = INTEGER(d1)[0], m1 = INTEGER(d1)[1];
  const int nP2 = INTEGER(d2)[0], m2 = INTEGER(d2)[1];
  if ((int)XLENGTH(Z1_ids_) != nP1) Rf_error("marker: length(Z1_ids) must equal nrow(Z1)");
  if ((int)XLENGTH(Z2_ids_) != nP2) Rf_error("marker: length(Z2_ids) must equal nrow(Z2)");

  KMethod km = K_STD;
  if (str_eq(method_, "std")) km = K_STD;
  else if (str_eq(method_, "vanraden")) km = K_VANRADEN;
  else if (str_eq(method_, "sommer")) km = K_SOMMER;
  else Rf_error("marker: method must be one of: 'std','vanraden','sommer'");

  const int kk = include_sca ? 3 : 2;

  /* match to parents */
  int *p1_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *p2_idx = (int*)R_alloc((size_t)ny, sizeof(int));
  int *use    = (int*)R_alloc((size_t)ny, sizeof(int));
  int nuse = 0;

  for (int i = 0; i < ny; ++i) {
    const char *a = CHAR(STRING_ELT(pool1_, i));
    const char *b = CHAR(STRING_ELT(pool2_, i));
    int ia = find_id_index(Z1_ids_, a);
    int ib = find_id_index(Z2_ids_, b);
    p1_idx[i] = ia;
    p2_idx[i] = ib;
    if (ia >= 0 && ib >= 0) use[nuse++] = i;
  }
  if (nuse <= 0) Rf_error("marker: No rows matched to both pool1 and pool2 IDs.");

  int *train_pos = (int*)malloc((size_t)nuse * sizeof(int));
  if (!train_pos) Rf_error("marker: OOM train_pos");
  int nobs = 0, npred = 0;
  for (int t = 0; t < nuse; ++t) {
    int i = use[t];
    if (!is_missing(REAL(y_)[i])) train_pos[nobs++] = t;
    else npred++;
  }
  if (nobs <= 0) { free(train_pos); Rf_error("marker: No non-missing y values after matching parents."); }

  int *a_tr = (int*)malloc((size_t)nobs * sizeof(int));
  int *b_tr = (int*)malloc((size_t)nobs * sizeof(int));
  double *ytr = (double*)malloc((size_t)nobs * sizeof(double));
  if (!a_tr || !b_tr || !ytr) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: OOM a_tr/b_tr/ytr");
  }
  for (int i = 0; i < nobs; ++i) {
    int obs = use[ train_pos[i] ];
    a_tr[i] = p1_idx[obs];
    b_tr[i] = p2_idx[obs];
    ytr[i]  = REAL(y_)[obs];
  }

  /* variance components: provided or bootstrap */
  double s1=NAN, s2=NAN, s3=0.0, se0=NAN;
  int vc_source_boot = 0;
  int boot_success = 0;

  if (!isNull(varcmp_)) {
    if (!isReal(varcmp_) || (int)XLENGTH(varcmp_) != (kk + 1)) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp must be numeric length (kk+1)=%d", kk+1);
    }
    s1  = REAL(varcmp_)[0];
    s2  = REAL(varcmp_)[1];
    s3  = include_sca ? REAL(varcmp_)[2] : 0.0;
    se0 = REAL(varcmp_)[kk];
  } else {
    if (!include_sca) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp NULL; bootstrap VC only enabled when include_sca=TRUE.");
    }
    BootSpec bs;
    int pr = parse_bootstrap(bootstrap_, &bs);
    if (pr != 0 || !bs.enabled) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: varcmp NULL and bootstrap not set (need nObs and nMarkers).");
    }
    if (bs.nObs > nobs) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: bootstrap$nObs (%d) must be <= n_train (%d).", bs.nObs, nobs);
    }
    if (bs.nMarkers > m1 || bs.nMarkers > m2) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: bootstrap$nMarkers must be <= ncol(Z1) and <= ncol(Z2).");
    }

    double vc_est[4] = {0,0,0,0};
    int rc = bootstrap_estimate_varcmp_hybrid(
      nobs, ytr, a_tr, b_tr,
      nP1, m1, nP2, m2,
      REAL(Z1_), REAL(Z2_),
      km, coding_, 1,
      asInteger(maxcyc_), asReal(tol_), ridge,
      &bs,
      vc_est, &boot_success
    );
    if (rc != 0) {
      free(train_pos); free(a_tr); free(b_tr); free(ytr);
      Rf_error("marker: bootstrap VC failed (rc=%d).", rc);
    }
    s1 = vc_est[0];
    s2 = vc_est[1];
    s3 = vc_est[2];
    se0= vc_est[3];
    vc_source_boot = 1;
  }

  if (!(s1 >= 0.0) || !(s2 >= 0.0) || !(s3 >= 0.0) || !(se0 > 0.0) ||
      !(s1==s1) || !(s2==s2) || !(s3==s3) || !(se0==se0)) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: invalid variance components (need s1,s2,s3>=0 and se>0).");
  }

  const double se = se0 + ridge;
  if (!(se > 0.0) || !(se==se)) {
    free(train_pos); free(a_tr); free(b_tr); free(ytr);
    Rf_error("marker: residual variance + ridge must be > 0.");
  }

  /* parental kernels from FULL markers */
  double *Z1rm = (double*)malloc((size_t)nP1*(size_t)m1*sizeof(double));
  double *Z2rm = (double*)malloc((size_t)nP2*(size_t)m2*sizeof(double));
  double *Kp1_rm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_rm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  double *Kp1_cm = (double*)malloc((size_t)nP1*(size_t)nP1*sizeof(double));
  double *Kp2_cm = (double*)malloc((size_t)nP2*(size_t)nP2*sizeof(double));
  if (!Z1rm||!Z2rm||!Kp1_rm||!Kp2_rm||!Kp1_cm||!Kp2_cm) Rf_error("marker: OOM kernels");

  Z_cm_to_rm(REAL(Z1_), nP1, m1, Z1rm);
  Z_cm_to_rm(REAL(Z2_), nP2, m2, Z2rm);

  double zmean1=0,zsd1=0,zmean2=0,zsd2=0;
  ZCoding cu1=CODING_OTHER, cu2=CODING_OTHER;
  ZTransform tr1=TRANS_NONE, tr2=TRANS_NONE;

  int rc1 = build_K_from_Z(Z1rm, nP1, m1, km, coding_, Kp1_rm, &zmean1, &zsd1, &cu1, &tr1);
  if (rc1 != 0) Rf_error("marker: build_K_from_Z(Z1) failed (code %d).", rc1);
  int rc2 = build_K_from_Z(Z2rm, nP2, m2, km, coding_, Kp2_rm, &zmean2, &zsd2, &cu2, &tr2);
  if (rc2 != 0) Rf_error("marker: build_K_from_Z(Z2) failed (code %d).", rc2);

  rm_to_cm_sym(nP1, Kp1_rm, Kp1_cm);
  rm_to_cm_sym(nP2, Kp2_rm, Kp2_cm);

  /* Sigma operator */
  SigmaOp op;
  op.nobs = nobs; op.nP1 = nP1; op.nP2 = nP2;
  op.a = a_tr; op.b = b_tr;
  op.Kp1_cm = Kp1_cm; op.Kp2_cm = Kp2_cm;
  op.include_sca = include_sca;
  op.s1 = s1; op.s2 = s2; op.s3 = include_sca ? s3 : 0.0; op.se = se;
  op.tmp1 = (double*)malloc((size_t)nobs*sizeof(double));
  op.tmp2 = (double*)malloc((size_t)nobs*sizeof(double));
  op.tmp3 = include_sca ? (double*)malloc((size_t)nobs*sizeof(double)) : NULL;
  if (!op.tmp1 || !op.tmp2 || (include_sca && !op.tmp3)) Rf_error("marker: OOM SigmaOp buffers");

  double tol_cg = asReal(tol_);
  if (!(tol_cg > 0.0) || !(tol_cg==tol_cg)) tol_cg = 1e-6;
  int maxit = 400;

  /* solve Wy and W1 */
  double *Wy = (double*)calloc((size_t)nobs, sizeof(double));
  double *W1 = (double*)calloc((size_t)nobs, sizeof(double));
  double *ones = (double*)malloc((size_t)nobs*sizeof(double));
  if (!Wy || !W1 || !ones) Rf_error("marker: OOM Wy/W1/ones");
  for (int i=0;i<nobs;++i) ones[i]=1.0;

  if (cg_solve_sigma(&op, ytr, Wy, maxit, tol_cg) != 0) Rf_error("marker: CG failed for Wy");
  if (cg_solve_sigma(&op, ones, W1, maxit, tol_cg) != 0) Rf_error("marker: CG failed for W1");

  double num=0.0, den=0.0;
  for (int i=0;i<nobs;++i) { num += Wy[i]; den += W1[i]; }
  if (!(den > 0.0) || !(den==den)) Rf_error("marker: invalid beta0 denominator");
  const double beta0 = num / den;

  double *r = (double*)malloc((size_t)nobs*sizeof(double));
  double *alpha = (double*)calloc((size_t)nobs, sizeof(double));
  if (!r || !alpha) Rf_error("marker: OOM r/alpha");
  for (int i=0;i<nobs;++i) r[i] = ytr[i] - beta0;
  if (cg_solve_sigma(&op, r, alpha, maxit, tol_cg) != 0) Rf_error("marker: CG failed for alpha");

  /* predict for all matched hybrids */
  SEXP yhat_obs = PROTECT(allocVector(REALSXP, ny)); nprot++;
  for (int i=0;i<ny;++i) REAL(yhat_obs)[i] = NA_REAL;

  /* precompute tA = Kp1 * sum_a(alpha), tB = Kp2 * sum_b(alpha) */
  double *sA = (double*)calloc((size_t)nP1, sizeof(double));
  double *sB = (double*)calloc((size_t)nP2, sizeof(double));
  if (!sA || !sB) Rf_error("marker: OOM sA/sB");
  for (int j=0;j<nobs;++j) { sA[a_tr[j]] += alpha[j]; sB[b_tr[j]] += alpha[j]; }

  double *tA = (double*)malloc((size_t)nP1*sizeof(double));
  double *tB = (double*)malloc((size_t)nP2*sizeof(double));
  if (!tA || !tB) Rf_error("marker: OOM tA/tB");

  {
    char uplo='U';
    BLAS_INT INC=1;
    const double one=1.0, zero=0.0;
    BLAS_INT N1=(BLAS_INT)nP1, N2=(BLAS_INT)nP2;
    F77_CALL(dsymv)(&uplo, &N1, &one, Kp1_cm, &N1, sA, &INC, &zero, tA, &INC FCONE);
    F77_CALL(dsymv)(&uplo, &N2, &one, Kp2_cm, &N2, sB, &INC, &zero, tB, &INC FCONE);
  }

  /* SCA precompute Tsca = Kp1 * Aalpha * Kp2^T */
  double *Tsca = NULL;
  if (include_sca) {
    const size_t R = (size_t)nP1*(size_t)nP2;
    double *Aalpha = (double*)calloc(R, sizeof(double)); /* col-major nP1 x nP2 */
    double *W = (double*)malloc(R * sizeof(double));
    Tsca = (double*)malloc(R * sizeof(double));
    if (!Aalpha || !W || !Tsca) Rf_error("marker: OOM SCA matrices");

    for (int j=0;j<nobs;++j) {
      Aalpha[a_tr[j] + (size_t)nP1*(size_t)b_tr[j]] += alpha[j];
    }

    {
      char transN='N';
      BLAS_INT M=(BLAS_INT)nP1, N=(BLAS_INT)nP2, K=(BLAS_INT)nP1;
      BLAS_INT lda=M, ldb=M, ldc=M;
      const double one=1.0, zero=0.0;
      F77_CALL(dgemm)(&transN,&transN,&M,&N,&K,&one, Kp1_cm,&lda, Aalpha,&ldb, &zero, W,&ldc FCONE FCONE);
    }
    {
      char transN='N', transT='T';
      BLAS_INT M=(BLAS_INT)nP1, N=(BLAS_INT)nP2, K=(BLAS_INT)nP2;
      BLAS_INT lda=(BLAS_INT)nP1, ldb=(BLAS_INT)nP2, ldc=(BLAS_INT)nP1;
      const double one=1.0, zero=0.0;
      F77_CALL(dgemm)(&transN,&transT,&M,&N,&K,&one, W,&lda, Kp2_cm,&ldb, &zero, Tsca,&ldc FCONE FCONE);
    }

    free(Aalpha);
    free(W);
  }

  for (int t=0; t<nuse; ++t) {
    int obs = use[t];
    int a0 = p1_idx[obs];
    int b0 = p2_idx[obs];
    double pred = beta0;
    pred += s1 * tA[a0];
    pred += s2 * tB[b0];
    if (include_sca) pred += s3 * Tsca[a0 + (size_t)nP1*(size_t)b0];
    REAL(yhat_obs)[obs] = pred;
  }

  /* outputs */
  SEXP beta0_out = PROTECT(allocVector(REALSXP,1)); nprot++;
  REAL(beta0_out)[0] = beta0;

  SEXP varcmp_out = PROTECT(allocVector(REALSXP, kk+1)); nprot++;
  REAL(varcmp_out)[0] = s1;
  REAL(varcmp_out)[1] = s2;
  if (include_sca) REAL(varcmp_out)[2] = s3;
  REAL(varcmp_out)[kk] = se0;

  int ndiag = 12;
  SEXP diag = PROTECT(allocVector(VECSXP, ndiag)); nprot++;
  SEXP dn   = PROTECT(allocVector(STRSXP, ndiag)); nprot++;
  SET_STRING_ELT(dn, 0, mkChar("n_total"));
  SET_STRING_ELT(dn, 1, mkChar("n_matched"));
  SET_STRING_ELT(dn, 2, mkChar("n_train"));
  SET_STRING_ELT(dn, 3, mkChar("n_pred"));
  SET_STRING_ELT(dn, 4, mkChar("kk"));
  SET_STRING_ELT(dn, 5, mkChar("include_sca"));
  SET_STRING_ELT(dn, 6, mkChar("model_type"));
  SET_STRING_ELT(dn, 7, mkChar("vc_source"));
  SET_STRING_ELT(dn, 8, mkChar("boot_n_success"));
  SET_STRING_ELT(dn, 9, mkChar("boot_nRep"));
  SET_STRING_ELT(dn,10, mkChar("boot_nObs"));
  SET_STRING_ELT(dn,11, mkChar("boot_nMarkers"));
  setAttrib(diag, R_NamesSymbol, dn);

  SET_VECTOR_ELT(diag, 0, ScalarInteger(ny));
  SET_VECTOR_ELT(diag, 1, ScalarInteger(nuse));
  SET_VECTOR_ELT(diag, 2, ScalarInteger(nobs));
  SET_VECTOR_ELT(diag, 3, ScalarInteger(npred));
  SET_VECTOR_ELT(diag, 4, ScalarInteger(kk));
  SET_VECTOR_ELT(diag, 5, ScalarLogical(include_sca));
  SET_VECTOR_ELT(diag, 6, mkString("marker"));
  SET_VECTOR_ELT(diag, 7, mkString(vc_source_boot ? "bootstrap" : "provided"));
  SET_VECTOR_ELT(diag, 8, ScalarInteger(boot_success));
  if (vc_source_boot) {
    BootSpec bs; (void)parse_bootstrap(bootstrap_, &bs);
    SET_VECTOR_ELT(diag, 9, ScalarInteger(bs.nRep));
    SET_VECTOR_ELT(diag,10, ScalarInteger(bs.nObs));
    SET_VECTOR_ELT(diag,11, ScalarInteger(bs.nMarkers));
  } else {
    SET_VECTOR_ELT(diag, 9, ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(diag,10, ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(diag,11, ScalarInteger(NA_INTEGER));
  }

  SEXP ans = PROTECT(allocVector(VECSXP, 6)); nprot++;
  SEXP an  = PROTECT(allocVector(STRSXP, 6)); nprot++;
  SET_STRING_ELT(an, 0, mkChar("beta0"));
  SET_STRING_ELT(an, 1, mkChar("varcmp"));
  SET_STRING_ELT(an, 2, mkChar("yhat_obs"));
  SET_STRING_ELT(an, 3, mkChar("diag"));
  SET_STRING_ELT(an, 4, mkChar("retval"));
  SET_STRING_ELT(an, 5, mkChar("ridge"));
  setAttrib(ans, R_NamesSymbol, an);

  SET_VECTOR_ELT(ans, 0, beta0_out);
  SET_VECTOR_ELT(ans, 1, varcmp_out);
  SET_VECTOR_ELT(ans, 2, yhat_obs);
  SET_VECTOR_ELT(ans, 3, diag);
  SET_VECTOR_ELT(ans, 4, ScalarInteger(0));
  SET_VECTOR_ELT(ans, 5, ScalarReal(ridge));

  /* cleanup */
  free(train_pos);
  free(a_tr); free(b_tr); free(ytr);

  free(Z1rm); free(Z2rm);
  free(Kp1_rm); free(Kp2_rm);
  free(Kp1_cm); free(Kp2_cm);

  free(op.tmp1); free(op.tmp2); if (op.tmp3) free(op.tmp3);

  free(Wy); free(W1); free(ones);
  free(r); free(alpha);

  free(sA); free(sB);
  free(tA); free(tB);
  if (Tsca) free(Tsca);

  UNPROTECT(nprot);
  return ans;
}

#endif

static SEXP st_mixed_hybrid_dispatch_impl(
  SEXP y_, SEXP hybrid_id_,
  SEXP pool1_, SEXP pool2_,
  SEXP Z1_, SEXP Z1_ids_,
  SEXP Z2_, SEXP Z2_ids_,
  SEXP include_sca_,
  SEXP method_, SEXP coding_, SEXP varcmp_,
  SEXP maxcyc_, SEXP tol_, SEXP ridge_,
  SEXP model_type_, SEXP bootstrap_)
{
  
  ModelType mt = MODEL_KERNEL;
  int pr = parse_model_type(model_type_, &mt);
  if (pr != 0) Rf_error("hybrid: invalid model_type (must be 'kernel' or 'marker')");

  if (mt == MODEL_MARKER) {
    return st_mixed_hybrid_marker_impl(
      y_, hybrid_id_,
      pool1_, pool2_,
      Z1_, Z1_ids_,
      Z2_, Z2_ids_,
      include_sca_,
      method_, coding_, varcmp_,
      maxcyc_, tol_, ridge_,
      bootstrap_
    );
  }

  /* IMPORTANT: proven kernel path stays exactly as-is */
  return st_mixed_hybrid_impl(
    y_, hybrid_id_,
    pool1_, pool2_,
    Z1_, Z1_ids_,
    Z2_, Z2_ids_,
    include_sca_,
    method_, coding_, varcmp_,
    maxcyc_, tol_, ridge_
  );
}

//

SEXP st_mixed(SEXP args)
{
  
  if (TYPEOF(args) != VECSXP) {
    Rf_error("st_mixed: expected a list of arguments");
  }

  /* calcVC (currently not branching; kept for compatibility) */
  int calcVC = 1;
  (void)calcVC;

  {
    SEXP cvc = getListElement(args, "calcVC");
    if (cvc != R_NilValue) {
      int v = asLogical(cvc);
      if (v != NA_LOGICAL) calcVC = v ? 1 : 0;
    }
  }

  /* ---- IMMEDIATE DISPATCH ---- */
  SEXP mode = getListElement(args, "mode");
  if (mode != R_NilValue && TYPEOF(mode) == STRSXP && XLENGTH(mode) == 1) {

    const char *m = CHAR(STRING_ELT(mode, 0));

    if (strcmp(m, "single") == 0) {
      return st_mixed_single_impl(
        getListElement(args, "y"),
        getListElement(args, "id"),
        getListElement(args, "Z"),
        getListElement(args, "Z_ids"),
        getListElement(args, "method"),
        getListElement(args, "coding"),
        getListElement(args, "varcmp"),
        getListElement(args, "maxcyc"),
        getListElement(args, "tol")
      );
    }

    if (strcmp(m, "hybrid") == 0) {
      /* IMPORTANT: call DISPATCHER, not st_mixed_hybrid_impl */
      return st_mixed_hybrid_dispatch_impl(
        getListElement(args, "y"),
        getListElement(args, "hybrid_id"),
        getListElement(args, "pool1"),
        getListElement(args, "pool2"),
        getListElement(args, "Z1"),
        getListElement(args, "Z1_ids"),
        getListElement(args, "Z2"),
        getListElement(args, "Z2_ids"),
        getListElement(args, "include_sca"),
        getListElement(args, "method"),
        getListElement(args, "coding"),
        getListElement(args, "varcmp"),
        getListElement(args, "maxcyc"),
        getListElement(args, "tol"),
        getListElement(args, "ridge"),
        getListElement(args, "model_type"),  /* NEW */
        getListElement(args, "bootstrap")    /* NEW */
      );
    }

    Rf_error("st_mixed: unknown mode='%s'", m);
  }

  /* Fallback if 'mode' not present: detect hybrid by pool1/pool2 */
  if (getListElement(args, "pool1") != R_NilValue ||
      getListElement(args, "pool2") != R_NilValue) {

    return st_mixed_hybrid_dispatch_impl(
      getListElement(args, "y"),
      getListElement(args, "hybrid_id"),
      getListElement(args, "pool1"),
      getListElement(args, "pool2"),
      getListElement(args, "Z1"),
      getListElement(args, "Z1_ids"),
      getListElement(args, "Z2"),
      getListElement(args, "Z2_ids"),
      getListElement(args, "include_sca"),
      getListElement(args, "method"),
      getListElement(args, "coding"),
      getListElement(args, "varcmp"),
      getListElement(args, "maxcyc"),
      getListElement(args, "tol"),
      getListElement(args, "ridge"),
      getListElement(args, "model_type"),   /* NEW */
      getListElement(args, "bootstrap")     /* NEW */
    );
  }

  /* default: single-kernel */
  return st_mixed_single_impl(
    getListElement(args, "y"),
    getListElement(args, "id"),
    getListElement(args, "Z"),
    getListElement(args, "Z_ids"),
    getListElement(args, "method"),
    getListElement(args, "coding"),
    getListElement(args, "varcmp"),
    getListElement(args, "maxcyc"),
    getListElement(args, "tol")
  );
}
//##############################################################################
//# END: Extend st_mixed for full matrix based genomic precistion           2026
//##############################################################################

#include "SelectionTools01a.h"
