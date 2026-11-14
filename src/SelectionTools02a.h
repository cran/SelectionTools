
SEXP define_map_df(SEXP map)
{
    if (cMP != 0)
    {
        info(-2, "First remove previously loaded linkage map\n");
        return R_NilValue;
    }
    if (!validate_map_genome())
        error("A valid genome is required before defining a linkage map");
    if (TYPEOF(map) != VECSXP || XLENGTH(map) < 4)
        error("A linkage map requires chromosome, position, name and class columns");
    SEXP chromosome = VECTOR_ELT(map, 0), position = VECTOR_ELT(map, 1),
         name = VECTOR_ELT(map, 2), class = VECTOR_ELT(map, 3);
    R_xlen_t length = XLENGTH(chromosome);
    if (TYPEOF(chromosome) != INTSXP || TYPEOF(position) != REALSXP ||
        TYPEOF(name) != STRSXP || TYPEOF(class) != STRSXP ||
        XLENGTH(position) != length || XLENGTH(name) != length ||
        XLENGTH(class) != length || length > INT_MAX)
        error("Invalid linkage-map column types or lengths");
    if (!length)
        return R_NilValue;
    if ((size_t)length > SIZE_MAX / sizeof(TMappoint) ||
        length > INT_MAX / NoHoms)
        error("The linkage map is too large");
    st_map_extension *stage = calloc(1, sizeof(*stage));
    if (!stage)
        error("Cannot allocate linkage-map stage");
    st_work_adopt(stage, st_map_extension_free);
    stage->n = (int)length;
    stage->points = calloc((size_t)length, sizeof(*stage->points));
    if (!stage->points)
        error("Cannot allocate linkage-map points");
    for (int i = 0; i < stage->n; i++)
    {
        int chrom = INTEGER(chromosome)[i];
        double pos = REAL(position)[i];
        if (chrom == NA_INTEGER || chrom < 1 || chrom > NoChroms ||
            !R_FINITE(pos) || pos < 0 || pos > PopDescr[chrom - 1].ChromLen ||
            STRING_ELT(name, i) == NA_STRING || STRING_ELT(class, i) == NA_STRING)
            error("Invalid linkage-map chromosome, position, name or class");
        TMappoint *point = stage->points + i;
        point->chrom = chrom;
        point->Pos = pos;
        point->counter = i;
        point->name = st_text_dup(CHAR(STRING_ELT(name, i)));
        point->class = st_text_dup(CHAR(STRING_ELT(class, i)));
        if (!point->name || !point->class)
            error("Cannot allocate linkage-map labels");
    }
    qsort(stage->points, (size_t)stage->n, sizeof(*stage->points), mccmp);
    if (!install_mappoints(stage->points, stage->n))
        error("Cannot install linkage map");
    st_work_detach(stage);
    st_map_extension_free(stage);
    return R_NilValue;
}

static void st_define_effect_stage_free(void *pointer)
{
    st_sim_ef_free(pointer);
    free(pointer);
}

static int st_define_effects_native(const char *effname,
                                    const char **marker_names,
                                    const int *alleles, const double *values,
                                    R_xlen_t nrows)
{
    int i, hom;
    TEffects *Effect = NULL;

    if (NULL != find_Effect((char *)effname))
    {
        info(-2, "This effect is already loaded (same name)\n");
        return -1;
    }

    if (strlen(effname) >= sizeof(Effect->Name))
    {
        info(-2, "Effect name too long\n");
        return -1;
    }

    if (nrows > INT_MAX)
    {
        info(-2, "Too many effect rows\n");
        return -1;
    }

    if ((0 < nrows) && ((size_t)nrows > SIZE_MAX / sizeof(*Effect->Eff)))
    {
        info(-2, "Memory allocation failed in define_effects_df\n");
        return -1;
    }

    if (NoHoms < 1)
    {
        info(-2, "Invalid number of homologues\n");
        return -1;
    }

    if ((size_t)NoHoms > SIZE_MAX / sizeof(int))
    {
        info(-2, "Memory allocation failed in define_effects_df\n");
        return -1;
    }

    for (R_xlen_t row = 0; row < nrows; row++)
        if (!marker_names[row] || !isfinite(values[row]) || fabs(values[row]) > FLT_MAX)
        {
            info(-2, "Effect names and finite float values are required\n");
            return -1;
        }
    Effect = calloc(1, sizeof(*Effect));
    if (!Effect)
    {
        info(-2, "Cannot allocate effect-definition stage\n");
        return -1;
    }
    st_work_adopt(Effect, st_define_effect_stage_free);
    Effect->GValue_avg = 0;
    Effect->weight = 1;
    strcpy(Effect->Name, effname);

    if (0 < nrows)
    {
        Effect->Eff = malloc((size_t)nrows * sizeof(*Effect->Eff));

        if (NULL == Effect->Eff)
        {
            info(-2, "Memory allocation failed in define_effects_df\n");
            goto fail;
        }
    }

    for (i = 0; i < (int)nrows; i++)
    {
        const char *markerName;
        TMappoint *tmpMP;
        TEff *newEff;
        int allele;

        markerName = marker_names[i];
        allele = alleles[i];
        tmpMP = find_mappoint((char *)markerName);

        if (NULL == tmpMP)
        {
            info(-1, "Marker %s not found.\nrow %d\n", markerName, i + 1);
            goto fail;
        }

        newEff = &(Effect->Eff[Effect->NrEff]);
        newEff->Eff = (float)values[i];
        newEff->NrLoci = 1;
        newEff->Loc = NULL;
        Effect->NrEff++;

        newEff->Loc = malloc(sizeof(*newEff->Loc));

        if (NULL == newEff->Loc)
        {
            info(-2, "Memory allocation failed in define_effects_df\n");
            goto fail;
        }

        newEff->Loc[0].MP = tmpMP;
        newEff->Loc[0].NrAll = 1;
        newEff->Loc[0].All = NULL;

        newEff->Loc[0].All =
            malloc((size_t)NoHoms * sizeof(*newEff->Loc[0].All));

        if (NULL == newEff->Loc[0].All)
        {
            info(-2, "Memory allocation failed in define_effects_df\n");
            goto fail;
        }

        for (hom = 0; hom < NoHoms; hom++)
            newEff->Loc[0].All[hom] = -1;

        newEff->Loc[0].All[0] = allele;
    }

    if (!st_sim_ef_finalize(Effect))
    {
        info(-2, "Memory allocation failed in define_effects_df\n");
        goto fail;
    }

    if (!st_effect_commit(Effect))
    {
        info(-2, "Cannot install effect definition\n");
        goto fail;
    }
    st_work_detach(Effect);
    st_define_effect_stage_free(Effect);

    st_clear_all_population_gvalue();

    return 0;

fail:

    st_work_detach(Effect);
    st_define_effect_stage_free(Effect);
    return -1;
}

SEXP define_effects_df(SEXP df, SEXP nameS)
{
    if (TYPEOF(df) != VECSXP || LENGTH(df) < 3 || TYPEOF(nameS) != STRSXP ||
        LENGTH(nameS) != 1 || STRING_ELT(nameS, 0) == NA_STRING)
        error("Invalid effect definition");
    SEXP marker = VECTOR_ELT(df, 0), allele = VECTOR_ELT(df, 1),
         value = VECTOR_ELT(df, 2);
    R_xlen_t n = XLENGTH(marker);
    if (TYPEOF(marker) != STRSXP || TYPEOF(allele) != INTSXP ||
        TYPEOF(value) != REALSXP || XLENGTH(allele) != n ||
        XLENGTH(value) != n || n > INT_MAX)
        error("Invalid effect definition columns");
    const char **names = st_work_alloc(n ? n : 1, sizeof(*names));
    for (R_xlen_t i = 0; i < n; i++)
    {
        if (STRING_ELT(marker, i) == NA_STRING)
            error("Effect locus names must not be missing");
        names[i] = CHAR(STRING_ELT(marker, i));
    }
    st_define_effects_native(CHAR(STRING_ELT(nameS, 0)), names,
                             n ? INTEGER(allele) : NULL,
                             n ? REAL(value) : NULL, n);
    return R_NilValue;
}
