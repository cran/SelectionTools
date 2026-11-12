

void gs_marker_stats_01_df(gs_varset_type *gv, int *retval)
{
    char *empty = "";
    char **empty_p = &empty;
    int zero = 0;

    gs_marker_stats_01(gv, empty_p, empty_p, empty_p, &zero, &zero, &zero,
                       &zero, retval);
}

SEXP gs_read_marker_data_DF_list(SEXP markermatrix, SEXP dataset)
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
        (NA_STRING == STRING_ELT(dataset, 0)))
        return R_NilValue;

    dst = gs_fdta((char *)translateCharUTF8(STRING_ELT(dataset, 0)));
    gs_tmp_init_from(&tmp, dst);

    if ((TYPEOF(markermatrix) != VECSXP) || (XLENGTH(markermatrix) < 3))
    {
        gs_info(dst, -2,
                "Marker dataframe must contain at least three columns");
        return R_NilValue;
    }

    col_ind = VECTOR_ELT(markermatrix, 0);
    col_mar = VECTOR_ELT(markermatrix, 1);
    col_all = VECTOR_ELT(markermatrix, 2);

    if ((TYPEOF(col_ind) != STRSXP) || (TYPEOF(col_mar) != STRSXP) ||
        (TYPEOF(col_all) != STRSXP))
    {
        gs_info(dst, -2, "Marker dataframe columns must be character");
        return R_NilValue;
    }

    nrows = XLENGTH(col_ind);
    if ((nrows < 1) || (XLENGTH(col_mar) != nrows) ||
        (XLENGTH(col_all) != nrows) ||
        ((uintmax_t)nrows > (uintmax_t)gs_UNSIGNED_MAX))
    {
        gs_info(dst, -2, "Invalid marker dataframe dimensions");
        return R_NilValue;
    }

    ind = gs_malloc_array((size_t)nrows, sizeof(*ind));
    mar = gs_malloc_array((size_t)nrows, sizeof(*mar));
    names = gs_malloc_array((size_t)nrows, sizeof(*names));
    if ((NULL == ind) || (NULL == mar) || (NULL == names))
    {
        gs_info(dst, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (R_xlen_t k = 0; k < nrows; k++)
    {
        SEXP si = STRING_ELT(col_ind, k);
        SEXP sm = STRING_ELT(col_mar, k);
        size_t li, lm;

        if ((NA_STRING == si) || (NA_STRING == sm))
        {
            gs_info(dst, -2, "Missing individual or marker name in dataframe");
            failed = 1;
            goto cleanup;
        }

        ind[k] = (char *)translateCharUTF8(si);
        mar[k] = (char *)translateCharUTF8(sm);
        li = strlen(ind[k]);
        lm = strlen(mar[k]);

        if ((0 == li) || (0 == lm) || (li >= INT_MAX) || (lm >= gs_MAXNME))
        {
            gs_info(dst, -2, "Invalid or too long individual/marker name");
            failed = 1;
            goto cleanup;
        }

        if ((int)li > tmp.MaxIL)
            tmp.MaxIL = (int)li;
        if ((int)lm > tmp.MaxML)
            tmp.MaxML = (int)lm;
    }

    for (R_xlen_t k = 0; k < nrows; k++)
        names[k] = ind[k];
    qsort(names, (size_t)nrows, sizeof(*names), gs_string_cmp);
    tmp.NoInd = 1;
    for (R_xlen_t k = 1; k < nrows; k++)
        if (0 != strcmp(names[k - 1], names[k]))
            tmp.NoInd++;

    tmp.IndName = calloc((size_t)tmp.NoInd, sizeof(*tmp.IndName));
    if (NULL == tmp.IndName)
    {
        gs_info(dst, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    {
        gs_UNSIGNED i = 0;
        if (st_identity_set(&tmp, i++, names[0], NULL))
        {
            failed = 1;
            goto cleanup;
        }
        for (R_xlen_t k = 1; k < nrows; k++)
            if (0 != strcmp(names[k - 1], names[k]))
                if (st_identity_set(&tmp, i++, names[k], NULL))
                {
                    failed = 1;
                    goto cleanup;
                }
    }

    for (R_xlen_t k = 0; k < nrows; k++)
        names[k] = mar[k];
    qsort(names, (size_t)nrows, sizeof(*names), gs_string_cmp);
    tmp.NoMar = 1;
    for (R_xlen_t k = 1; k < nrows; k++)
        if (0 != strcmp(names[k - 1], names[k]))
            tmp.NoMar++;

    tmp.MarName = gs_malloc_array((size_t)tmp.NoMar, sizeof(*tmp.MarName));
    if (NULL == tmp.MarName)
    {
        gs_info(dst, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    {
        gs_UNSIGNED m = 0;
        strcpy(tmp.MarName[m++], names[0]);
        for (R_xlen_t k = 1; k < nrows; k++)
            if (0 != strcmp(names[k - 1], names[k]))
                strcpy(tmp.MarName[m++], names[k]);
    }

    if (0 != gs_alloc_genotype_matrix(&tmp, tmp.NoMar, tmp.NoInd))
    {
        gs_info(dst, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }
    gs_set_genotypes_missing(&tmp);

    {
        size_t n_genotypes;
        if (0 !=
            gs_size_mul((size_t)tmp.NoMar, (size_t)tmp.NoInd, &n_genotypes))
        {
            gs_info(dst, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }
        geno_state = calloc(n_genotypes, sizeof(*geno_state));
        if (NULL == geno_state)
        {
            gs_info(dst, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }
    }

    for (R_xlen_t k = 0; k < nrows; k++)
    {
        gs_UNSIGNED m =
            gs_find_sorted_fixed_name(tmp.MarName, tmp.NoMar, mar[k]);
        gs_UNSIGNED i = gs_find_sorted_ind_name(tmp.IndName, tmp.NoInd, ind[k]);
        GS_ALLELE_TYPE a, b;
        size_t idx;
        int pair = 1;
        SEXP sa = STRING_ELT(col_all, k);

        if ((gs_UNSIGNED_MAX == m) || (gs_UNSIGNED_MAX == i))
        {
            gs_info(dst, -2, "Internal marker-data name lookup error");
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
            const char *allele_text = translateCharUTF8(sa);
            size_t allele_len = strlen(allele_text);

            if (0 != gs_determine_alleles(&a, &b, allele_text))
            {
                gs_info(dst, -2, "Invalid allele in marker dataframe");
                failed = 1;
                goto cleanup;
            }

            pair = ((NULL != strchr(allele_text, '/')) ||
                    (!gs_allele_missing_span(allele_text, allele_len) &&
                     (2 == allele_len)));
        }

        idx = (size_t)m * (size_t)tmp.NoInd + (size_t)i;
        if (0 != gs_aggregate_genotype(&(tmp.allele_raw1[idx]),
                                       &(tmp.allele_raw2[idx]),
                                       &(geno_state[idx]), a, b, pair))
        {
            gs_info(dst, -2,
                    "Conflicting genotype records in marker dataframe");
            failed = 1;
            goto cleanup;
        }
    }

    gs_marker_stats_01_df(&tmp, &rv);
    if (-2 == rv)
    {
        failed = 1;
        goto cleanup;
    }

    if (gs_tmp_commit(dst, &tmp))
    {
        failed = 1;
        goto cleanup;
    }

cleanup:
    free(ind);
    free(mar);
    free(names);
    free(geno_state);

    if (failed)
        gs_reset(&tmp);

    return R_NilValue;
}

static int gs_df_matrix_cell(SEXP data, R_xlen_t k, int *a, int *b, int *kind)
{
    *a = *b = GS_ALLELE_MISSING;
    *kind = 0;
    if (TYPEOF(data) == STRSXP)
    {
        SEXP cell = STRING_ELT(data, k);
        if (cell == NA_STRING)
            return 0;
        const char *v = translateCharUTF8(cell);
        if (!strcmp(v, "--") || !strcmp(v, "-"))
            return 0;
        size_t n = strlen(v);
        if (n == 1 && v[0] >= '0' && v[0] <= '2')
        {
            *a = v[0] == '2' ? 2 : 1;
            *b = v[0] == '0' ? 1 : 2;
            *kind = 1;
            return 0;
        }
        if (n == 3 && (v[0] == '1' || v[0] == '2') &&
            (v[2] == '1' || v[2] == '2') && (v[1] == '/' || v[1] == '|'))
        {
            *a = v[0] - '0';
            *b = v[2] - '0';
            *kind = v[1] == '|' ? 2 : 1;
            return 0;
        }
        return -1;
    }
    double d =
        TYPEOF(data) == INTSXP
            ? (INTEGER(data)[k] == NA_INTEGER ? NA_REAL : INTEGER(data)[k])
            : REAL(data)[k];
    if (ISNAN(d))
        return 0;
    if (d != 0 && d != 1 && d != 2)
        return -1;
    *a = d == 2 ? 2 : 1;
    *b = d == 0 ? 1 : 2;
    *kind = 1;
    return 0;
}

SEXP gs_read_marker_data_DF_matrix(SEXP data, SEXP dataset)
{
    SEXP dims = getAttrib(data, R_DimSymbol);
    SEXP dn = getAttrib(data, R_DimNamesSymbol);
    if ((TYPEOF(dataset) != STRSXP) || XLENGTH(dataset) != 1 ||
        STRING_ELT(dataset, 0) == NA_STRING ||
        !strlen(translateCharUTF8(STRING_ELT(dataset, 0))))
        Rf_error("data.set must be one non-empty character string");
    if ((TYPEOF(data) != INTSXP && TYPEOF(data) != REALSXP &&
         TYPEOF(data) != STRSXP) ||
        TYPEOF(dims) != INTSXP || LENGTH(dims) != 2 || INTEGER(dims)[0] < 1 ||
        INTEGER(dims)[1] < 1)
        Rf_error("expected a non-empty numeric or character genotype matrix");
    if (TYPEOF(dn) != VECSXP || LENGTH(dn) != 2)
        Rf_error("individual row names and marker column names are required");

    if (strlen(translateCharUTF8(STRING_ELT(dataset, 0))) >= gs_MAXSTR)
        Rf_error("data.set must contain fewer than %d bytes", gs_MAXSTR);

    int ni = INTEGER(dims)[0], nm = INTEGER(dims)[1];
    if ((uintmax_t)ni > (uintmax_t)gs_UNSIGNED_MAX ||
        (uintmax_t)nm > (uintmax_t)gs_UNSIGNED_MAX)
        Rf_error("matrix dimensions exceed supported limits");

    for (int axis = 0; axis < 2; axis++)
    {
        SEXP names = VECTOR_ELT(dn, axis);
        int n = axis == 0 ? ni : nm;
        if (TYPEOF(names) != STRSXP || XLENGTH(names) != n)
            Rf_error(
                "individual row names and marker column names are required");
        char **sorted = (char **)st_work_alloc(n, sizeof(char *));
        for (int k = 0; k < n; k++)
        {
            if (STRING_ELT(names, k) == NA_STRING)
                Rf_error("names must not be missing");
            sorted[k] = (char *)translateCharUTF8(STRING_ELT(names, k));
            size_t len = strlen(sorted[k]);
            if (!len || len >= (axis == 0 ? INT_MAX : gs_MAXNME))
                Rf_error("names must contain 1 to %d bytes", gs_MAXNME - 1);
            for (size_t j = 0; j < len; j++)
                if (isspace((unsigned char)sorted[k][j]))
                    Rf_error("names must not contain whitespace");
        }
        qsort(sorted, n, sizeof(char *), gs_string_cmp);
        for (int k = 1; k < n; k++)
            if (!strcmp(sorted[k - 1], sorted[k]))
                Rf_error("individual and marker names must each be unique");
    }
    int saw_phased = 0, saw_unphased = 0;
    for (R_xlen_t k = 0; k < XLENGTH(data); k++)
    {
        int a, b, kind;
        if (gs_df_matrix_cell(data, k, &a, &b, &kind))
            Rf_error(
                "invalid genotype at row %lld, column %lld: expected 0/1/2 dosage, allele pair using 1 and 2 with / or |, or NA/--/-",
                (long long)(k % ni + 1), (long long)(k / ni + 1));
        if (kind == 2)
            saw_phased = 1;
        if (kind == 1)
            saw_unphased = 1;
    }

    gs_varset_type *dst =
        gs_fdta((char *)translateCharUTF8(STRING_ELT(dataset, 0)));
    gs_varset_type tmp;
    gs_tmp_init_from(&tmp, dst);
    tmp.phased = saw_phased && !saw_unphased;
    tmp.NoInd = (gs_UNSIGNED)ni;
    tmp.NoMar = (gs_UNSIGNED)nm;
    tmp.IndName = calloc(ni, sizeof(*tmp.IndName));
    tmp.MarName = gs_malloc_array(nm, sizeof(*tmp.MarName));
    if (!tmp.IndName || !tmp.MarName || gs_alloc_genotype_matrix(&tmp, nm, ni))
    {
        gs_reset(&tmp);
        Rf_error("could not allocate marker data");
    }
    for (int i = 0; i < ni; i++)
    {
        const char *name = translateCharUTF8(STRING_ELT(VECTOR_ELT(dn, 0), i));
        if (st_identity_set(&tmp, i, name, NULL))
        {
            gs_reset(&tmp);
            Rf_error("Cannot allocate individual name");
        }
        if ((int)strlen(name) > tmp.MaxIL)
            tmp.MaxIL = (int)strlen(name);
    }
    for (int m = 0; m < nm; m++)
    {
        const char *name = translateCharUTF8(STRING_ELT(VECTOR_ELT(dn, 1), m));
        strcpy(tmp.MarName[m], name);
        if ((int)strlen(name) > tmp.MaxML)
            tmp.MaxML = (int)strlen(name);
        for (int i = 0; i < ni; i++)
        {
            R_xlen_t k = (R_xlen_t)m * ni + i;
            int a, b, kind;

            gs_df_matrix_cell(data, k, &a, &b, &kind);
            tmp.allele_raw1[k] = a;
            tmp.allele_raw2[k] = b;
        }
    }
    int rv = -2;
    gs_marker_stats_01_df(&tmp, &rv);
    if (rv == -2)
    {
        gs_reset(&tmp);
        Rf_error("could not calculate marker statistics");
    }
    if (gs_tmp_commit(dst, &tmp))
    {
        gs_reset(&tmp);
        Rf_error(
            "Cannot prepare allele storage; existing data set was preserved");
    }
    return R_NilValue;
}

static int gs_prepare_marker_map_df(gs_varset_type *gv, gs_mappoint_type *map_s,
                                    gs_UNSIGNED nloc, gs_varset_type *tmp)
{
    gs_mappoint_type **map = NULL;
    gs_name_index_type *map_names = NULL, *marker_names = NULL;
    gs_UNSIGNED *match = NULL;
    gs_UNSIGNED no_matching = 0;
    int failed = 0;
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
    {
        map[lne] = &(map_s[lne]);
    }

    qsort(map, nloc, sizeof(*map), gs_mappoint_cmp);

    for (gs_UNSIGNED first = 0; first < nloc;)
    {
        gs_UNSIGNED last = first + 1;
        double base = map[first]->pos;
        int chrom = map[first]->chrom;

        while ((last < nloc) && (map[last]->chrom == chrom) &&
               (map[last]->pos == base))
        {
            last++;
        }

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
                {
                    step = gap / (double)k;
                }
            }

            for (gs_UNSIGNED j = 1; j < k; j++)
            {
                double pos = base + step * (double)j;
                double prev = map[first + j - 1]->pos;

                if (!(pos > prev))
                {
                    pos = nextafter(prev, INFINITY);
                }

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
            {
                a++;
            }
            else
            {
                b++;
            }
        }
    }

    if (0 == no_matching)
    {
        gs_info(gv, -2, "No linkage-map markers match marker data");
        failed = 1;
        goto cleanup;
    }

    tmp->NoInd = gv->NoInd;
    tmp->NoMar = no_matching;
    tmp->phased = gv->phased;
    tmp->NLoc = no_matching;
    tmp->MaxIL = gv->MaxIL;
    tmp->MaxML = gv->MaxML;

    for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
    {
        if (gs_UNSIGNED_MAX != match[lne])
        {
            if ((0 !=
                 gs_pool_add_string_size(&(tmp->map_str_N), map[lne]->name)) ||
                (0 !=
                 gs_pool_add_string_size(&(tmp->map_str_N), map[lne]->class)))
            {
                gs_info(gv, -2, "Linkage map string storage too large");
                failed = 1;
                goto cleanup;
            }
        }
    }
    tmp->IndName = calloc((size_t)tmp->NoInd, sizeof(*tmp->IndName));
    tmp->MarName = gs_malloc_array((size_t)tmp->NoMar, sizeof(*tmp->MarName));
    tmp->map_S = gs_malloc_array((size_t)tmp->NLoc, sizeof(*tmp->map_S));
    tmp->map_str_S = gs_malloc_array(tmp->map_str_N, sizeof(*tmp->map_str_S));
    tmp->map = gs_malloc_array((size_t)tmp->NLoc, sizeof(*tmp->map));
    tmp->mdl = gs_malloc_array((size_t)tmp->NLoc, sizeof(*tmp->mdl));

    if ((NULL == tmp->IndName) || (NULL == tmp->MarName) ||
        (NULL == tmp->map_S) || (NULL == tmp->map_str_S) ||
        (NULL == tmp->map) || (NULL == tmp->mdl) ||
        (0 != gs_alloc_genotype_matrix(tmp, tmp->NoMar, tmp->NoInd)))
    {
        gs_info(gv, -2, "Memory allocation error");
        failed = 1;
        goto cleanup;
    }

    for (gs_UNSIGNED i = 0; i < tmp->NoInd; i++)
    {
        if (st_identity_copy(tmp, i, gv, i))
        {
            failed = 1;
            goto cleanup;
        }
    }

    {
        gs_UNSIGNED m_new = 0;
        size_t map_str_offset = 0;

        for (gs_UNSIGNED lne = 0; lne < nloc; lne++)
        {
            if (gs_UNSIGNED_MAX != match[lne])
            {
                gs_UNSIGNED m_old = match[lne];

                strcpy(tmp->MarName[m_new], gv->MarName[m_old]);
                st_allele_copy(gv, 0, (size_t)m_old * gv->NoInd,
                               (size_t)tmp->NoInd,
                               &(tmp->allele_raw1[(size_t)m_new * tmp->NoInd]));
                st_allele_copy(gv, 1, (size_t)m_old * gv->NoInd,
                               (size_t)tmp->NoInd,
                               &(tmp->allele_raw2[(size_t)m_new * tmp->NoInd]));

                if (0 != gs_copy_map_point_to_pool(
                             &(tmp->map_S[m_new]), map[lne], m_new,
                             tmp->map_str_S, tmp->map_str_N, &map_str_offset))
                {
                    gs_info(gv, -2, "Invalid linkage map string storage");
                    failed = 1;
                    goto cleanup;
                }
                tmp->map[m_new] = &(tmp->map_S[m_new]);
                tmp->mdl[m_new] = m_new;
                m_new++;
            }
        }

        if (map_str_offset != tmp->map_str_N)
        {
            gs_info(gv, -2, "Invalid linkage map string storage");
            failed = 1;
            goto cleanup;
        }
    }

    if (NULL != gv->y)
    {
        tmp->y = gs_malloc_array((size_t)tmp->NoInd, sizeof(*tmp->y));
        if (NULL == tmp->y)
        {
            gs_info(gv, -2, "Memory allocation error");
            failed = 1;
            goto cleanup;
        }
        memcpy(tmp->y, gv->y, (size_t)tmp->NoInd * sizeof(*tmp->y));
    }

    {
        int srv = 0;
        char *empty = "";
        char **empty_p = &empty;
        int zero = 0;
        int *zero_p = &zero;

        gs_marker_stats_01(tmp, empty_p, empty_p, empty_p, zero_p, zero_p,
                           zero_p, zero_p, &srv);
        if (-2 == srv)
        {
            failed = 1;
            goto cleanup;
        }
    }

cleanup:
    free(map);
    free(map_names);
    free(marker_names);
    free(match);
    return failed ? -1 : 0;
}

SEXP st_read_map_df(SEXP data, SEXP stretch, SEXP dataset)
{
    if (TYPEOF(dataset) != STRSXP || XLENGTH(dataset) != 1 ||
        STRING_ELT(dataset, 0) == NA_STRING ||
        !strlen(CHAR(STRING_ELT(dataset, 0))) ||
        strlen(CHAR(STRING_ELT(dataset, 0))) >= gs_MAXSTR)
    {
        Rf_error("invalid data.set name");
    }
    if (TYPEOF(stretch) != REALSXP || XLENGTH(stretch) != 1 ||
        !R_FINITE(REAL(stretch)[0]) || REAL(stretch)[0] <= 0)
    {
        Rf_error("m.stretch must be finite and positive");
    }
    if (TYPEOF(data) != VECSXP || XLENGTH(data) != 3)
    {
        Rf_error("map must have three MCP columns");
    }
    SEXP names = VECTOR_ELT(data, 0), chrom = VECTOR_ELT(data, 1),
         pos = VECTOR_ELT(data, 2);
    R_xlen_t n = XLENGTH(names);
    if (TYPEOF(names) != STRSXP || TYPEOF(chrom) != REALSXP ||
        TYPEOF(pos) != REALSXP || n < 1 || XLENGTH(chrom) != n ||
        XLENGTH(pos) != n || (uintmax_t)n > gs_UNSIGNED_MAX)
    {
        Rf_error("invalid map columns or dimensions");
    }
    for (R_xlen_t i = 0; i < n; i++)
    {
        if (STRING_ELT(names, i) == NA_STRING)
        {
            Rf_error("marker names must not be missing");
        }
        const char *name = CHAR(STRING_ELT(names, i));
        size_t len = strlen(name);
        if (!len || len >= gs_MAXSTR)
        {
            Rf_error("invalid marker name length");
        }
        for (size_t j = 0; j < len; j++)
        {
            if (isspace((unsigned char)name[j]))
            {
                Rf_error("marker names must not contain whitespace");
            }
        }
        double c = REAL(chrom)[i], p = REAL(pos)[i] * REAL(stretch)[0];
        if (!R_FINITE(c) || c < 1 || c > INT_MAX || c != floor(c))
        {
            Rf_error("chromosomes must be positive integers");
        }
        if (!R_FINITE(p) || p < 0)
        {
            Rf_error("positions must be finite and non-negative");
        }
    }

    const char *name = CHAR(STRING_ELT(dataset, 0));
    gs_varset_type *dst = !strcmp(name, "default") ? gv : NULL;
    for (gs_UNSIGNED i = 0; !dst && i < gs_ndta; i++)
    {
        if (!strcmp(name, gs_data[i]->name))
        {
            dst = gs_data[i];
        }
    }
    if (!dst || !dst->NoInd || !dst->NoMar || !dst->IndName || !dst->MarName ||
        !st_allele_present(dst))
    {
        Rf_error("no marker data available in data.set");
    }

    gs_mappoint_type *points = gs_malloc_array((size_t)n, sizeof(*points));
    if (!points)
    {
        Rf_error("could not allocate map");
    }
    for (R_xlen_t i = 0; i < n; i++)
    {
        points[i].name = (char *)CHAR(STRING_ELT(names, i));
        points[i].class = "m";
        points[i].chrom = (int)REAL(chrom)[i];
        points[i].pos = REAL(pos)[i] * REAL(stretch)[0];
        points[i].nr = gs_UNSIGNED_MAX;
    }
    gs_varset_type tmp;
    gs_tmp_init_from(&tmp, dst);
    int status = gs_prepare_marker_map_df(dst, points, (gs_UNSIGNED)n, &tmp);
    free(points);
    if (status)
    {
        gs_reset(&tmp);
        Rf_error(
            "map loading failed; existing data set was preserved (see preceding message)");
    }
    if (gs_tmp_commit(dst, &tmp))
    {
        gs_reset(&tmp);
        Rf_error(
            "Cannot prepare allele storage; existing data set was preserved");
    }
    return R_NilValue;
}
