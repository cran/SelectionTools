
SEXP define_map_df(SEXP map)
{
  int i;
  int counter = 0;
  
  if (cMP != 0) {
    info(-2, "First remove previously loaded linkage map\n");
    return R_NilValue;
  }
  
  SEXP chromS = VECTOR_ELT(map, 0);
  SEXP posS   = VECTOR_ELT(map, 1);
  SEXP nameS  = VECTOR_ELT(map, 2);
  SEXP classS = VECTOR_ELT(map, 3);
  
  int n = LENGTH(chromS);
  
  for(i=0;i<n;i++){
    int    chrom = INTEGER(chromS)[i];  
    double pos   = REAL(posS)[i];
    
    PopDescr[ chrom - 1 ].NrMappoints++;
    cMP++;
    
    // If this is the first marker for this chromosome
    // Create a array of 1 TMappoint
    if ( 1 == PopDescr[ chrom - 1 ].NrMappoints )
      PopDescr[ chrom - 1 ].Mappoint = ( TMappoint* ) yCalloc( 1, TMappoint );
    else
      // Else, reallocate memory depending on the size of NrMappoints, probably always one more
      PopDescr[ chrom - 1 ].Mappoint = ( TMappoint* ) yRealloc 
    ( PopDescr[ chrom - 1 ].Mappoint, PopDescr[ chrom - 1 ].NrMappoints, TMappoint );
    
    {// Create a pointer to the current Mapppoint
      TMappoint*tmpPoint = &PopDescr[ chrom - 1 ].Mappoint[ PopDescr[ chrom - 1 ].NrMappoints - 1 ];
      ( *tmpPoint ).chrom = chrom; // Set Chromosome
      ( *tmpPoint ).Pos = pos; // Set position
      ( *tmpPoint ).name = ( char* ) yCalloc( strlen( CHAR(STRING_ELT(nameS, i)) ) + 1, char ); // Allocate memory size of lociname
      ( *tmpPoint ).class = ( char* ) yCalloc( strlen( CHAR(STRING_ELT(classS, i)) ) + 1, char ); // Allocate memory size of classname
      ( *tmpPoint ).counter = counter++; // At an index to the marker
      strcpy( ( *tmpPoint ).name, CHAR(STRING_ELT(nameS, i)) ); // Add lociname to mappoint
      strcpy( ( *tmpPoint ).class, CHAR(STRING_ELT(classS, i)) ); // add classname to mappoint
    }
  }
  
  // Marker times homologues?
  cMPhom = NoHoms * cMP;
  // SEXP must return something, this is like NULL in R
  return R_NilValue;
}

/*
SEXP define_effects_df(SEXP df, SEXP nameS)
{
  // Alt
  int i, j, k, tmpi, maxLoc = 0;
  float tmpf; 
  TMappoint*tmpMP = NULL;
  TEffects *Effect;
  
  // Neu: Spalten aus DF, nrows, Effektname aus R
  SEXP col_marker, col_allele, col_value;
  R_xlen_t nrows;
  const char *effname = CHAR(STRING_ELT(nameS, 0));
  
  // Alt: nimmt jetzt effname aus R
  if (NULL != find_Effect((char*)effname)) {
    info(-2, "This effect is already loaded (same name)\n");
    return R_NilValue;
  }
  
  // Alt
  if (0 == NrEffects++)
    Effects = (TEffects*) yCalloc(1, TEffects);
  else
    Effects = (TEffects*) yRealloc(Effects, NrEffects, TEffects);
  
  Effect = &Effects[NrEffects - 1];
  
  // Alt
  Effect->GValue_avg = 0;
  Effect->weight = 1;
  Effect->NrEff = 0;
  Effect->Pos = NrEffects;
  Effect->Eff = NULL;
  Effect->MPset = NULL;
  
  // Alt: Name kommt jetzt aus R
  strcpy(Effect->Name, effname);
  
  // Neu: statt einem File haben wir jetzt die Spalten aus dem DF
  col_marker = VECTOR_ELT(df, 0); 
  col_allele = VECTOR_ELT(df, 1); 
  col_value  = VECTOR_ELT(df, 2); 
  nrows = XLENGTH(col_marker);
  
  // Iterieren über alle Zeilen im DF
  for (i = 0; i < (int)nrows; i++)
  {
    const char *markerName;
    
    // Locus
    markerName = CHAR(STRING_ELT(col_marker, i));
    // Allel
    tmpi = INTEGER(col_allele)[i];  
    // Effektwert
    tmpf = (float) REAL(col_value)[i];
    // Finde den Mappoint zu diesem Locus
    tmpMP = find_mappoint((char*)markerName);
    
    if (tmpMP == NULL) {
      info(-1, "Marker %s not found.\nrow %d\n", markerName, i + 1);
      return R_NilValue;
    }
    
    // Alt: Array erweitern
    if (0 == (Effect->NrEff)++)
      Effect->Eff = (TEff*) yCalloc(1, TEff);
    else
      Effect->Eff = (TEff*) yRealloc(Effect->Eff, Effect->NrEff, TEff);
    
    // Alt: Neuen Effekt mit null initialisieren
    // Neu: auch .Loc auf NULL initialisieren
    Effect->Eff[Effect->NrEff - 1].Eff    = 0;
    Effect->Eff[Effect->NrEff - 1].NrLoci = 0;
    Effect->Eff[Effect->NrEff - 1].Loc    = NULL;
    
    // Alles alt ab hier
    Insert_ALL_Loc(&Effect->Eff[Effect->NrEff - 1], tmpMP, tmpi);
    
    Effect->Eff[Effect->NrEff - 1].Eff = tmpf;
    
    if (maxLoc < Effect->Eff[Effect->NrEff - 1].NrLoci)
      maxLoc = Effect->Eff[Effect->NrEff - 1].NrLoci;
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
          Effect->MPset = (TMappoint**) yCalloc(1, TMappoint*);
        else
          Effect->MPset = (TMappoint**) yRealloc(Effect->MPset, Effect->MPs, TMappoint*);
        
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
                (Effect->Eff[j].Loc[i].All[k] != Effect->Eff[j].Loc[i].All[n]))
            {
              Effect->NrEff++;
              Effect->Eff = (TEff*) yRealloc(Effect->Eff, Effect->NrEff, TEff);
              
              Effect->Eff[Effect->NrEff - 1].Eff    = Effect->Eff[j].Eff;
              Effect->Eff[Effect->NrEff - 1].NrLoci = Effect->Eff[j].NrLoci;
              Effect->Eff[Effect->NrEff - 1].Loc    = (TGen*) yCalloc(Effect->Eff[j].NrLoci, TGen);
              
              for (l = 0; l < Effect->Eff[j].NrLoci; l++)
              {
                Effect->Eff[Effect->NrEff - 1].Loc[l].MP    = Effect->Eff[j].Loc[l].MP;
                Effect->Eff[Effect->NrEff - 1].Loc[l].NrAll = Effect->Eff[j].Loc[l].NrAll;
                Effect->Eff[Effect->NrEff - 1].Loc[l].All   = (int*) yCalloc(NoHoms, int);
                
                for (m = 0; m < NoHoms; m++)
                {
                  if ((l == i) && (m == k))
                    Effect->Eff[Effect->NrEff - 1].Loc[l].All[m] = Effect->Eff[j].Loc[l].All[n];
                  else if ((l == i) && (m == n))
                    Effect->Eff[Effect->NrEff - 1].Loc[l].All[m] = Effect->Eff[j].Loc[l].All[k];
                  else
                    Effect->Eff[Effect->NrEff - 1].Loc[l].All[m] = Effect->Eff[j].Loc[l].All[m];
                }
              }
            }
          }
        }
      }
    }
  }
  
  // Must return something
  return R_NilValue;
}
*/


SEXP define_effects_df(SEXP df, SEXP nameS)
{
  int i, hom;
  TEffects Effect;
  TEffects *tmpEffects;
  SEXP col_marker, col_allele, col_value;
  R_xlen_t nrows;
  const char *effname;

  effname = CHAR(STRING_ELT(nameS, 0));

  if (NULL != find_Effect((char*)effname)) {
    info(-2, "This effect is already loaded (same name)\n");
    return R_NilValue;
  }

  if (strlen(effname) >= sizeof(Effect.Name)) {
    info(-2, "Effect name too long\n");
    return R_NilValue;
  }

  col_marker = VECTOR_ELT(df, 0);
  col_allele = VECTOR_ELT(df, 1);
  col_value  = VECTOR_ELT(df, 2);
  nrows = XLENGTH(col_marker);

  if (nrows > INT_MAX) {
    info(-2, "Too many effect rows\n");
    return R_NilValue;
  }

  if ((0 < nrows) &&
      ((size_t)nrows > SIZE_MAX / sizeof(*Effect.Eff))) {
    info(-2, "Memory allocation failed in define_effects_df\n");
    return R_NilValue;
  }

  if (NoHoms < 1) {
    info(-2, "Invalid number of homologues\n");
    return R_NilValue;
  }

  if ((size_t)NoHoms > SIZE_MAX / sizeof(int)) {
    info(-2, "Memory allocation failed in define_effects_df\n");
    return R_NilValue;
  }

  memset(&Effect, 0, sizeof(Effect));
  Effect.GValue_avg = 0;
  Effect.weight = 1;
  strcpy(Effect.Name, effname);

  if (0 < nrows) {
    Effect.Eff = malloc((size_t)nrows * sizeof(*Effect.Eff));

    if (NULL == Effect.Eff) {
      info(-2, "Memory allocation failed in define_effects_df\n");
      return R_NilValue;
    }
  }

  for (i = 0; i < (int)nrows; i++) {
    const char *markerName;
    TMappoint *tmpMP;
    TEff *newEff;
    int allele;

    markerName = CHAR(STRING_ELT(col_marker, i));
    allele = INTEGER(col_allele)[i];
    tmpMP = find_mappoint((char*)markerName);

    if (NULL == tmpMP) {
      info(-1, "Marker %s not found.\nrow %d\n", markerName, i + 1);
      goto fail;
    }

    newEff = &(Effect.Eff[Effect.NrEff]);
    newEff->Eff = (float)REAL(col_value)[i];
    newEff->NrLoci = 1;
    newEff->Loc = NULL;
    Effect.NrEff++;

    newEff->Loc = malloc(sizeof(*newEff->Loc));

    if (NULL == newEff->Loc) {
      info(-2, "Memory allocation failed in define_effects_df\n");
      goto fail;
    }

    newEff->Loc[0].MP = tmpMP;
    newEff->Loc[0].NrAll = 1;
    newEff->Loc[0].All = NULL;

    newEff->Loc[0].All =
      malloc((size_t)NoHoms * sizeof(*newEff->Loc[0].All));

    if (NULL == newEff->Loc[0].All) {
      info(-2, "Memory allocation failed in define_effects_df\n");
      goto fail;
    }

    for (hom = 0; hom < NoHoms; hom++)
      newEff->Loc[0].All[hom] = -1;

    newEff->Loc[0].All[0] = allele;
  }

  if (!st_sim_ef_finalize(&Effect)) {
    info(-2, "Memory allocation failed in define_effects_df\n");
    goto fail;
  }

  if (INT_MAX == NrEffects) {
    info(-2, "Too many effects\n");
    goto fail;
  }

  if ((size_t)(NrEffects + 1) > SIZE_MAX / sizeof(*Effects)) {
    info(-2, "Memory allocation failed in define_effects_df\n");
    goto fail;
  }

  tmpEffects = realloc(
    Effects,
    (size_t)(NrEffects + 1) * sizeof(*Effects)
  );

  if (NULL == tmpEffects) {
    info(-2, "Memory allocation failed in define_effects_df\n");
    goto fail;
  }

  Effects = tmpEffects;
  Effect.Pos = NrEffects + 1;
  Effects[NrEffects] = Effect;
  NrEffects++;

  st_clear_all_population_gvalue();

  return R_NilValue;

fail:

  st_sim_ef_free(&Effect);
  return R_NilValue;
}
