
SEXP define_map_df(SEXP map)
{
    int i;
    int counter = 0;

    if (cMP != 0)
    {
        info(-2, "First remove previously loaded linkage map\n");
        return R_NilValue;
    }

    SEXP chromS = VECTOR_ELT(map, 0);
    SEXP posS = VECTOR_ELT(map, 1);
    SEXP nameS = VECTOR_ELT(map, 2);
    SEXP classS = VECTOR_ELT(map, 3);

    int n = LENGTH(chromS);

    for (i = 0; i < n; i++)
    {
        int chrom = INTEGER(chromS)[i];
        double pos = REAL(posS)[i];

        PopDescr[chrom - 1].NrMappoints++;
        cMP++;

        if (1 == PopDescr[chrom - 1].NrMappoints)
            PopDescr[chrom - 1].Mappoint = (TMappoint *)yCalloc(1, TMappoint);
        else

            PopDescr[chrom - 1].Mappoint = (TMappoint *)yRealloc(
                PopDescr[chrom - 1].Mappoint, PopDescr[chrom - 1].NrMappoints,
                TMappoint);

        {
            TMappoint *tmpPoint =
                &PopDescr[chrom - 1]
                     .Mappoint[PopDescr[chrom - 1].NrMappoints - 1];
            (*tmpPoint).chrom = chrom;
            (*tmpPoint).Pos = pos;
            (*tmpPoint).name =
                (char *)yCalloc(strlen(CHAR(STRING_ELT(nameS, i))) + 1, char);
            (*tmpPoint).class =
                (char *)yCalloc(strlen(CHAR(STRING_ELT(classS, i))) + 1, char);
            (*tmpPoint).counter = counter++;
            strcpy((*tmpPoint).name, CHAR(STRING_ELT(nameS, i)));
            strcpy((*tmpPoint).class, CHAR(STRING_ELT(classS, i)));
        }
    }

    cMPhom = NoHoms * cMP;

    return R_NilValue;
}

static int st_define_effects_native(const char *effname,
                                    const char **marker_names,
                                    const int *alleles, const double *values,
                                    R_xlen_t nrows)
{
    int i, hom;
    TEffects Effect;
    TEffects *tmpEffects;

    if (NULL != find_Effect((char *)effname))
    {
        info(-2, "This effect is already loaded (same name)\n");
        return -1;
    }

    if (strlen(effname) >= sizeof(Effect.Name))
    {
        info(-2, "Effect name too long\n");
        return -1;
    }

    if (nrows > INT_MAX)
    {
        info(-2, "Too many effect rows\n");
        return -1;
    }

    if ((0 < nrows) && ((size_t)nrows > SIZE_MAX / sizeof(*Effect.Eff)))
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

    memset(&Effect, 0, sizeof(Effect));
    Effect.GValue_avg = 0;
    Effect.weight = 1;
    strcpy(Effect.Name, effname);

    if (0 < nrows)
    {
        Effect.Eff = malloc((size_t)nrows * sizeof(*Effect.Eff));

        if (NULL == Effect.Eff)
        {
            info(-2, "Memory allocation failed in define_effects_df\n");
            return -1;
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

        newEff = &(Effect.Eff[Effect.NrEff]);
        newEff->Eff = (float)values[i];
        newEff->NrLoci = 1;
        newEff->Loc = NULL;
        Effect.NrEff++;

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

    if (!st_sim_ef_finalize(&Effect))
    {
        info(-2, "Memory allocation failed in define_effects_df\n");
        goto fail;
    }

    if (INT_MAX == NrEffects)
    {
        info(-2, "Too many effects\n");
        goto fail;
    }

    if ((size_t)(NrEffects + 1) > SIZE_MAX / sizeof(*Effects))
    {
        info(-2, "Memory allocation failed in define_effects_df\n");
        goto fail;
    }

    tmpEffects = realloc(Effects, (size_t)(NrEffects + 1) * sizeof(*Effects));

    if (NULL == tmpEffects)
    {
        info(-2, "Memory allocation failed in define_effects_df\n");
        goto fail;
    }

    Effects = tmpEffects;
    Effect.Pos = NrEffects + 1;
    Effects[NrEffects] = Effect;
    NrEffects++;

    st_clear_all_population_gvalue();

    return 0;

fail:

    st_sim_ef_free(&Effect);
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
    st_define_effects_native(CHAR(STRING_ELT(nameS, 0)), names, INTEGER(allele),
                             REAL(value), n);
    return R_NilValue;
}
