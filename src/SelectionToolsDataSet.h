#ifndef SELECTIONTOOLS_DATA_SET_H
#define SELECTIONTOOLS_DATA_SET_H

#include <stddef.h>
#include <limits.h>
#include "SelectionToolsAllele.h"

#ifndef gs_MAXSTR
#define gs_MAXSTR 4096
#endif
#ifndef gs_MAXNME
#define gs_MAXNME 256
#endif

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
  char (*IndName)[gs_MAXNME];
  char (*MarName)[gs_MAXNME];
  GS_ALLELE_TYPE *All1;
  GS_ALLELE_TYPE *All2;
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
  gs_cross_type *crs_S;
  gs_cross_type **crs;
  gs_UNSIGNED NChrom;
  gs_chrom_type *chrom;
  double *t0,*t1,*t2,*t3,*t4,*t5,*t6,*t7,*t8,*t9;
  gs_UNSIGNED *u0;
} gs_varset_type;

#endif
