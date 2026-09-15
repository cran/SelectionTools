/*
void gs_marker_stats_01_df ( gs_varset_type *gv ,
                             int*   retval        )  {
  
  // Messages on number of markers and individuals
  snprintf(gv->msg, sizeof(gv->msg),
         "No. of individuals: %lu, no. of markers: %lu",
         gv->NoInd, gv->NoMar);
  gs_info(gv,0,gv->msg);
  
  gs_UNSIGNED i,k,m,idx,idx_all;
  int already_defined;
  
  // Determine alleles
  gv->NoAll = realloc ( gv->NoAll, gv->NoMar* sizeof(gs_UNSIGNED) );
  gv->All   = realloc ( gv->All,   gv->NoMar* sizeof(int[gs_MAXALL]));
  if ( (NULL == gv->NoAll)  || (NULL == gv->All) ) ERR_M;
  
  for ( m = 0; m < gv->NoMar; m++)
  {
    gv->NoAll[m] = 0;
    
    for ( i = 0; i < gv->NoInd ; i++)
    {
      idx_all = m*gv->NoInd + i ;
      
      already_defined = 0;
      for (k = 0; k < gv->NoAll[m]; k++){
        if ( gv->All1[idx_all] == gv->All[m][k] ) already_defined = 1;
        if (!already_defined){
          gv->All [m] [ gv->NoAll[m] ] = gv->All1[idx_all] ;
          gv->NoAll[m]++;
          if (gv->NoAll[m] == gs_MAXALL) {
            gs_info(gv,-2,"To many alleles");
            *retval = -2; return ;
          }
        }
      }
      
      
      already_defined = 0;
      for (k = 0; k < gv->NoAll[m]; k++){
        if ( gv->All2[idx_all] == gv->All[m][k] ) already_defined = 1;
        if (!already_defined){
          gv->All [m] [ gv->NoAll[m] ] = gv->All2[idx_all] ;
          gv->NoAll[m]++;
          if (gv->NoAll[m] == gs_MAXALL) {
            gs_info(gv,-2,"To many alleles per marker");
            *retval = -2; return;
          }
        }
      }
    }
    qsort( gv->All[m], gv->NoAll[m], sizeof(int), gs_intcmp);
  }
  
  // Determine allele frequencies
  gv->AllCnt = realloc ( gv->AllCnt, gv->NoMar*sizeof(gs_UNSIGNED[gs_MAXALL]));
  if (NULL == gv->AllCnt)  ERR_M;
  
  for ( m = 0; m < gv->NoMar; m++){
    for (k = 0; k < gv->NoAll[m]; k++)
    {
      gv->AllCnt [m][k] = 0;
      for ( i = 0; i < gv->NoInd ; i++) {
        idx_all = m*gv->NoInd + i ;
        if ( gv->All1[idx_all] == gv->All[m][k]) gv->AllCnt[m][k]++;
        if ( gv->All2[idx_all] == gv->All[m][k]) gv->AllCnt[m][k]++;
      }
    }
  }
  // Marker statistics: missing and pic
  gv->Pic = realloc ( gv->Pic, gv->NoMar* sizeof(double) );
  gv->Mis = realloc ( gv->Mis, gv->NoMar* sizeof(double) );
  if ( (NULL == gv->Pic)  || (NULL == gv->Mis) ) ERR_M;
  
  double observed;
  for ( m = 0; m < gv->NoMar; m++)
  {
    
    // percentage missing per marker
    if ( -1 == gv->All[m][0])
      gv->Mis[m] = ( (double) gv->AllCnt[m][0] ) / 2 / (double) gv->NoInd;
    else
      gv->Mis[m] = 0;
    
    // expected heterozygosity
    gv->Pic[m] = 1;
    observed  = (double)( 2 * gv->NoInd ) ;
    if ( -1 == gv->All[m][0])
      observed  -= gv->AllCnt[m][0];
    else
      gv->Pic[m] -= pow ( (double)gv->AllCnt[m][0]/observed , 2);
    
    if ( observed > 0) {
      for (k = 1; k < gv->NoAll[m]; k++) gv->Pic[m] -= pow ( (double)gv->AllCnt[m][k]/observed , 2);
    }
    else {
      gv->Pic[m] = 0;
    }
  }
  
  // Individual statistics: missing
  gv->IndMis = realloc ( gv->IndMis, gv->NoInd*sizeof(double) );
  if (NULL == gv->IndMis)  ERR_M;
  
  gs_UNSIGNED count_miss ;
  for (i=0; i<gv->NoInd; i++  )
  {
    count_miss = 0;
    for (m=0; m < gv->NoMar; m++)
    {
      idx = m*gv->NoInd + i;
      if (-1 == gv->All1[idx]) count_miss++;
    }
    gv->IndMis[i] = (double)count_miss / (double)gv->NoMar;
  }
  
  *retval = 0; return ;
  
}

*/

void gs_marker_stats_01_df ( gs_varset_type *gv ,
                             int*   retval        )
{
  char *empty = "";
  char **empty_p = &empty;
  int zero = 0;

  gs_marker_stats_01(gv,
                     empty_p,
                     empty_p,
                     empty_p,
                     &zero,
                     &zero,
                     &zero,
                     &zero,
                     retval);
}


SEXP gs_read_marker_data_01_DF(SEXP markermatrix, SEXP dataset)
{
  gs_varset_type *dst;
  gs_varset_type tmp;
  SEXP col_ind, col_mar, col_all;
  R_xlen_t nrows;
  char **ind = NULL;
  char **mar = NULL;
  char **names = NULL;
  unsigned char *geno_state = NULL;
  int failed = 0;
  int rv = -2;

  if ((TYPEOF(dataset) != STRSXP) || (XLENGTH(dataset) < 1) ||
      (NA_STRING == STRING_ELT(dataset,0)))
    return R_NilValue;

  dst = gs_fdta((char *)CHAR(STRING_ELT(dataset,0)));
  gs_tmp_init_from(&tmp,dst);

  if ((TYPEOF(markermatrix) != VECSXP) || (XLENGTH(markermatrix) < 3))
    {
      gs_info(dst,-2,"Marker dataframe must contain at least three columns");
      return R_NilValue;
    }

  col_ind = VECTOR_ELT(markermatrix,0);
  col_mar = VECTOR_ELT(markermatrix,1);
  col_all = VECTOR_ELT(markermatrix,2);

  if ((TYPEOF(col_ind) != STRSXP) ||
      (TYPEOF(col_mar) != STRSXP) ||
      (TYPEOF(col_all) != STRSXP))
    {
      gs_info(dst,-2,"Marker dataframe columns must be character");
      return R_NilValue;
    }

  nrows = XLENGTH(col_ind);
  if ((nrows < 1) ||
      (XLENGTH(col_mar) != nrows) ||
      (XLENGTH(col_all) != nrows) ||
      ((uintmax_t)nrows > (uintmax_t)gs_UNSIGNED_MAX))
    {
      gs_info(dst,-2,"Invalid marker dataframe dimensions");
      return R_NilValue;
    }

  ind = gs_malloc_array((size_t)nrows,sizeof(*ind));
  mar = gs_malloc_array((size_t)nrows,sizeof(*mar));
  names = gs_malloc_array((size_t)nrows,sizeof(*names));
  if ((NULL == ind) || (NULL == mar) || (NULL == names))
    {
      gs_info(dst,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }

  for (R_xlen_t k = 0; k < nrows; k++)
    {
      SEXP si = STRING_ELT(col_ind,k);
      SEXP sm = STRING_ELT(col_mar,k);
      size_t li,lm;

      if ((NA_STRING == si) || (NA_STRING == sm))
        {
          gs_info(dst,-2,"Missing individual or marker name in dataframe");
          failed = 1;
          goto cleanup;
        }

      ind[k] = (char *)CHAR(si);
      mar[k] = (char *)CHAR(sm);
      li = strlen(ind[k]);
      lm = strlen(mar[k]);

      if ((0 == li) || (0 == lm) ||
          (li >= gs_MAXNME) || (lm >= gs_MAXNME))
        {
          gs_info(dst,-2,"Invalid or too long individual/marker name");
          failed = 1;
          goto cleanup;
        }

      if ((int)li > tmp.MaxIL) tmp.MaxIL = (int)li;
      if ((int)lm > tmp.MaxML) tmp.MaxML = (int)lm;
    }

  for (R_xlen_t k = 0; k < nrows; k++) names[k] = ind[k];
  qsort(names,(size_t)nrows,sizeof(*names),gs_string_cmp);
  tmp.NoInd = 1;
  for (R_xlen_t k = 1; k < nrows; k++)
    if (0 != strcmp(names[k-1],names[k])) tmp.NoInd++;

  tmp.IndName = gs_malloc_array((size_t)tmp.NoInd,sizeof(*tmp.IndName));
  if (NULL == tmp.IndName)
    {
      gs_info(dst,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  {
    gs_UNSIGNED i = 0;
    strcpy(tmp.IndName[i++],names[0]);
    for (R_xlen_t k = 1; k < nrows; k++)
      if (0 != strcmp(names[k-1],names[k]))
        strcpy(tmp.IndName[i++],names[k]);
  }

  for (R_xlen_t k = 0; k < nrows; k++) names[k] = mar[k];
  qsort(names,(size_t)nrows,sizeof(*names),gs_string_cmp);
  tmp.NoMar = 1;
  for (R_xlen_t k = 1; k < nrows; k++)
    if (0 != strcmp(names[k-1],names[k])) tmp.NoMar++;

  tmp.MarName = gs_malloc_array((size_t)tmp.NoMar,sizeof(*tmp.MarName));
  if (NULL == tmp.MarName)
    {
      gs_info(dst,-2,"Memory allocation error");
      failed = 1;
      goto cleanup;
    }
  {
    gs_UNSIGNED m = 0;
    strcpy(tmp.MarName[m++],names[0]);
    for (R_xlen_t k = 1; k < nrows; k++)
      if (0 != strcmp(names[k-1],names[k]))
        strcpy(tmp.MarName[m++],names[k]);
  }

  if (0 != gs_alloc_genotype_matrix(&tmp,tmp.NoMar,tmp.NoInd))
    {
      gs_info(dst,-2,"Memory allocation error");
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
        gs_info(dst,-2,"Memory allocation error");
        failed = 1;
        goto cleanup;
      }
    geno_state = calloc(n_genotypes,sizeof(*geno_state));
    if (NULL == geno_state)
      {
        gs_info(dst,-2,"Memory allocation error");
        failed = 1;
        goto cleanup;
      }
  }

  for (R_xlen_t k = 0; k < nrows; k++)
    {
      gs_UNSIGNED m = gs_find_sorted_fixed_name(tmp.MarName,tmp.NoMar,mar[k]);
      gs_UNSIGNED i = gs_find_sorted_fixed_name(tmp.IndName,tmp.NoInd,ind[k]);
      GS_ALLELE_TYPE a,b;
      size_t idx;
      int pair = 1;
      SEXP sa = STRING_ELT(col_all,k);

      if ((gs_UNSIGNED_MAX == m) || (gs_UNSIGNED_MAX == i))
        {
          gs_info(dst,-2,"Internal marker-data name lookup error");
          failed = 1;
          goto cleanup;
        }

      if (NA_STRING == sa)
        {
          a = GS_ALLELE_MISSING;
          b = GS_ALLELE_MISSING;
        }
      else
        {
          const char *allele_text = CHAR(sa);
          size_t allele_len = strlen(allele_text);

          if (0 != gs_determine_alleles(&a,&b,allele_text))
            {
              gs_info(dst,-2,"Invalid allele in marker dataframe");
              failed = 1;
              goto cleanup;
            }

          pair = ((NULL != strchr(allele_text,'/')) ||
                  (!gs_allele_missing_span(allele_text,allele_len) &&
                   (2 == allele_len)));
        }

      idx = (size_t)m * (size_t)tmp.NoInd + (size_t)i;
      if (0 != gs_aggregate_genotype(&(tmp.All1[idx]),
                                     &(tmp.All2[idx]),
                                     &(geno_state[idx]),
                                     a,b,pair))
        {
          gs_info(dst,-2,"Conflicting genotype records in marker dataframe");
          failed = 1;
          goto cleanup;
        }
    }

  gs_marker_stats_01_df(&tmp,&rv);
  if (-2 == rv)
    {
      failed = 1;
      goto cleanup;
    }

  gs_tmp_commit(dst,&tmp);

cleanup:
  free(ind);
  free(mar);
  free(names);
  free(geno_state);

  if (failed)
    gs_reset(&tmp);

  return R_NilValue;
}

