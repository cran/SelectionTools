

#include <stdlib.h>
#include <string.h>

#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>

void st_daxpy(gs_UNSIGNED n, double a, double *x, double *y)

{
    if (n == 0 || a == 0.0)
        return;

    for (gs_UNSIGNED i = 0; i < n; ++i)
        y[i] += a * x[i];
}

static void rm_to_cm(const double *A_rm, int m, int n, int lda_rm, double *A_cm,
                     int ldA_cm)
{
    for (int i = 0; i < m; ++i)
    {
        const double *row = A_rm + (size_t)i * (size_t)lda_rm;
        for (int j = 0; j < n; ++j)
        {
            A_cm[i + (size_t)j * (size_t)ldA_cm] = row[j];
        }
    }
}

static void cm_to_rm(const double *A_cm, int m, int n, int ldA_cm, double *A_rm,
                     int lda_rm)
{
    for (int i = 0; i < m; ++i)
    {
        double *row = A_rm + (size_t)i * (size_t)lda_rm;
        for (int j = 0; j < n; ++j)
        {
            row[j] = A_cm[i + (size_t)j * (size_t)ldA_cm];
        }
    }
}

void st_dgemm(char T1[], char T2[], gs_UNSIGNED nra, gs_UNSIGNED nca,
              gs_UNSIGNED nrb, gs_UNSIGNED ncb, double *A, double *B, double *C)
{

    const int A_rows = (int)nra, A_cols = (int)nca;
    const int B_rows = (int)nrb, B_cols = (int)ncb;

    const int transA = (!strcmp(T1, "N")) ? 0 : 1;
    const int transB = (!strcmp(T2, "N")) ? 0 : 1;

    const int M = transA ? A_cols : A_rows;
    const int K = transA ? A_rows : A_cols;
    const int N = transB ? B_rows : B_cols;

    const int Kb = transB ? B_cols : B_rows;
    if (K != Kb)
        return;

    const BLAS_INT ldA = (BLAS_INT)A_cols;
    const BLAS_INT ldB = (BLAS_INT)B_cols;
    const BLAS_INT ldC = (BLAS_INT)N;

    char ta = (transB == 0) ? 'N' : 'T';
    char tb = (transA == 0) ? 'N' : 'T';

    const BLAS_INT m_f = (BLAS_INT)N;
    const BLAS_INT n_f = (BLAS_INT)M;
    const BLAS_INT k_f = (BLAS_INT)K;

    double alpha = 1.0;
    double beta = 0.0;

    F77_CALL(dgemm)(&ta, &tb, &m_f, &n_f, &k_f, &alpha, B, &ldB, A, &ldA, &beta,
                    C, &ldC FCONE FCONE);
}

void st_dsymm_ZZt(gs_UNSIGNED n, gs_UNSIGNED m, const double *Z_rm,
                  double *G_rm)
{

    const int N = (int)n;
    const int K = (int)m;

    char uplo = 'U';
    char trans = 'T';

    const BLAS_INT n_f = (BLAS_INT)N;
    const BLAS_INT k_f = (BLAS_INT)K;

    const BLAS_INT lda_f = (BLAS_INT)K;

    const BLAS_INT ldc_f = (BLAS_INT)N;

    double alpha = 1.0;
    double beta = 0.0;

    F77_CALL(dsyrk)(&uplo, &trans, &n_f, &k_f, &alpha, (double *)Z_rm, &lda_f,
                    &beta, G_rm, &ldc_f FCONE FCONE);

    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            G_rm[(size_t)i * (size_t)N + (size_t)j] =
                G_rm[(size_t)j * (size_t)N + (size_t)i];
}

void st_dcopy(gs_UNSIGNED n, double *dest, double *src)
{
    memcpy(dest, src, (size_t)n * sizeof(double));
}

void st_dscal(gs_UNSIGNED n, double a, double *X)
{
    for (gs_UNSIGNED i = 0; i < n; ++i)
        X[i] *= a;
}

static int dgetrf_rm(int n, double *A_rm, int lda_rm, int *ipiv)
{
    int info = 0;
    int ldA = (n > 0) ? n : 1;

    double *A_cm = (double *)malloc((size_t)ldA * (size_t)n * sizeof(double));
    if (!A_cm)
        return -1010;

    rm_to_cm(A_rm, n, n, lda_rm, A_cm, ldA);

    F77_CALL(dgetrf)(&n, &n, A_cm, &ldA, ipiv, &info);

    cm_to_rm(A_cm, n, n, ldA, A_rm, lda_rm);

    free(A_cm);
    return info;
}

void st_dgetrf(double *a, gs_UNSIGNED nn, gs_UNSIGNED *indx, double *d)

{
    int n = (int)nn;
    int lda_rm = (int)nn;

    int *ipiv = (int *)malloc((size_t)n * sizeof(int));
    if (!ipiv)
    {
        gs_info(gv, -1, "Out of memory in st_dgetrf");
        *d = (double)(-1010);
        return;
    }

    int info = dgetrf_rm(n, a, lda_rm, ipiv);

    for (gs_UNSIGNED i = 0; i < nn; ++i)
        indx[i] = (gs_UNSIGNED)ipiv[i];

    *d = (double)info;

    if (info > 0)
    {
        gs_info(gv, -1, "Singular coefficient matrix in dgetrf");
    }

    free(ipiv);
}

static int dgetrs_rm(int n, const double *A_rm, int lda_rm, const int *ipiv,
                     double *b_rm)
{
    int info = 0;
    int nrhs = 1;

    int ldA = (n > 0) ? n : 1;
    int ldB = (n > 0) ? n : 1;

    double *A_cm = (double *)malloc((size_t)ldA * (size_t)n * sizeof(double));
    double *B_cm = (double *)malloc((size_t)ldB * sizeof(double));

    if (!A_cm || !B_cm)
    {
        free(A_cm);
        free(B_cm);
        return -1010;
    }

    rm_to_cm(A_rm, n, n, lda_rm, A_cm, ldA);

    for (int i = 0; i < n; ++i)
        B_cm[i] = b_rm[i];

    char trans = 'N';
    F77_CALL(dgetrs)(&trans, &n, &nrhs, A_cm, &ldA, (int *)ipiv, B_cm, &ldB,
                     &info FCONE);

    for (int i = 0; i < n; ++i)
        b_rm[i] = B_cm[i];

    free(B_cm);
    free(A_cm);

    return info;
}

void st_dgetrs(double *a, gs_UNSIGNED nn, gs_UNSIGNED *indx, double *b)
{
    int n = (int)nn;
    int lda_rm = (int)nn;

    int *ipiv = (int *)malloc((size_t)n * sizeof(int));
    if (!ipiv)
    {
        gs_info(gv, -1, "Out of memory in st_dgetrs");
        return;
    }

    for (int i = 0; i < n; ++i)
        ipiv[i] = (int)indx[i];

    int info = dgetrs_rm(n, a, lda_rm, ipiv, b);

    if (info != 0)
        gs_info(gv, -1, "Error in dgetrs");

    free(ipiv);
}

static int dgetri_rm(int n, double *A_rm, int lda_rm, const int *ipiv)
{
    int info = 0;
    int ldA = (n > 0) ? n : 1;

    double *A_cm = (double *)malloc((size_t)ldA * (size_t)n * sizeof(double));
    if (!A_cm)
        return -1010;

    rm_to_cm(A_rm, n, n, lda_rm, A_cm, ldA);

    int lwork = -1;
    double work_query = 0.0;

    F77_CALL(dgetri)(&n, A_cm, &ldA, (int *)ipiv, &work_query, &lwork, &info);
    if (info != 0)
    {
        free(A_cm);
        return info;
    }

    lwork = (int)work_query;
    if (lwork < 1)
        lwork = 1;

    double *work = (double *)malloc((size_t)lwork * sizeof(double));
    if (!work)
    {
        free(A_cm);
        return -1010;
    }

    F77_CALL(dgetri)(&n, A_cm, &ldA, (int *)ipiv, work, &lwork, &info);

    cm_to_rm(A_cm, n, n, ldA, A_rm, lda_rm);

    free(work);
    free(A_cm);

    return info;
}

void st_dgetri(double *a, gs_UNSIGNED nn, gs_UNSIGNED *indx, double *inv)
{
    int n = (int)nn;
    int lda_rm = (int)nn;

    int *ipiv = (int *)malloc((size_t)n * sizeof(int));
    if (!ipiv)
    {
        gs_info(gv, -1, "Out of memory in st_dgetri");
        return;
    }

    for (int i = 0; i < n; ++i)
        ipiv[i] = (int)indx[i];

    int info = dgetri_rm(n, a, lda_rm, ipiv);

    if (info != 0)
    {
        gs_info(gv, -1, "Error in dgetri");
        free(ipiv);
        return;
    }

    memcpy(inv, a, (size_t)n * (size_t)n * sizeof(double));

    free(ipiv);
}

void st_dgetri_inv(double *a, double *inva, gs_UNSIGNED n)
{
    double d;
    gs_UNSIGNED *indx = malloc(n * sizeof(gs_UNSIGNED));

    st_dgetrf(a, n, indx, &d);
    st_dgetri(a, n, indx, inva);

    free(indx);
}

void st_inv(double *a, double *inva, int *nn)
{
    gs_UNSIGNED n = (gs_UNSIGNED)*nn;

    st_dgetri_inv(a, inva, n);
}

static int dgesv_rm(int n, double *A_rm, int lda_rm, double *x_rm, int *ipiv)
{
    int info = 0;
    int nrhs = 1;

    int ldA = (n > 0) ? n : 1;
    int ldB = (n > 0) ? n : 1;

    double *A_cm = (double *)malloc((size_t)ldA * (size_t)n * sizeof(double));
    double *B_cm =
        (double *)malloc((size_t)ldB * (size_t)nrhs * sizeof(double));

    if (!A_cm || !B_cm)
    {
        free(A_cm);
        free(B_cm);
        return -1010;
    }

    rm_to_cm(A_rm, n, n, lda_rm, A_cm, ldA);

    for (int i = 0; i < n; ++i)
        B_cm[i] = x_rm[i];

    F77_CALL(dgesv)(&n, &nrhs, A_cm, &ldA, ipiv, B_cm, &ldB, &info);

    for (int i = 0; i < n; ++i)
        x_rm[i] = B_cm[i];

    free(B_cm);
    free(A_cm);

    return info;
}

void st_dgesv(gs_UNSIGNED nra, gs_UNSIGNED nca, double *A, double *x, double *b,
              int *retval)
{
    (void)nca;

    int n = (int)nra;
    int lda_rm = (int)nra;

    int *ipiv = (int *)malloc((size_t)n * sizeof(int));
    if (!ipiv)
    {
        gs_info(gv, -1, "Out of memory in st_dgesv");
        *retval = -2;
        return;
    }

    for (int i = 0; i < n; ++i)
        x[i] = b[i];

    int info = dgesv_rm(n, A, lda_rm, x, ipiv);

    free(ipiv);

    if (info > 0)
    {
        gs_info(gv, -1, "Singular matrix in the MME (dgesv)");
        *retval = -2;
        return;
    }

    *retval = 0;
}

static int dsyevr_vals_rm_m(int n, double *A_rm, int lda_rm, double *w,
                            int *m_out)
{
    int info = 0;

    int ldA = (n > 0) ? n : 1;
    double *A_cm = (double *)malloc((size_t)ldA * (size_t)n * sizeof(double));
    if (!A_cm)
        return -1010;

    rm_to_cm(A_rm, n, n, lda_rm, A_cm, ldA);

    char jobz = 'N';
    char range = 'A';
    char uplo = 'L';

    double vl = 0.0, vu = 0.0;
    int il = 0, iu = 0;
    double abstol = -1.0;

    int m = 0;
    int ldZ = 1;
    double z_dummy = 0.0;
    int isuppz_dummy[2];

    int lwork = -1, liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;

    F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, A_cm, &ldA, &vl, &vu, &il, &iu,
                     &abstol, &m, w, &z_dummy, &ldZ, isuppz_dummy, &work_query,
                     &lwork, &iwork_query, &liwork, &info FCONE FCONE FCONE);

    if (info != 0)
    {
        free(A_cm);
        if (m_out)
            *m_out = m;
        return info;
    }

    lwork = (int)work_query;
    liwork = (int)iwork_query;
    if (lwork < 1)
        lwork = 1;
    if (liwork < 1)
        liwork = 1;

    double *work = (double *)malloc((size_t)lwork * sizeof(double));
    int *iwork = (int *)malloc((size_t)liwork * sizeof(int));
    if (!work || !iwork)
    {
        free(work);
        free(iwork);
        free(A_cm);
        if (m_out)
            *m_out = m;
        return -1010;
    }

    F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, A_cm, &ldA, &vl, &vu, &il, &iu,
                     &abstol, &m, w, &z_dummy, &ldZ, isuppz_dummy, work, &lwork,
                     iwork, &liwork, &info FCONE FCONE FCONE);

    free(iwork);
    free(work);
    free(A_cm);

    if (m_out)
        *m_out = m;
    return info;
}

void st_dsyevr(double *a, double *w, gs_UNSIGNED N)
{
    int n = (int)N;
    int lda_rm = (int)N;

    int m = 0;
    int info = dsyevr_vals_rm_m(n, a, lda_rm, w, &m);

    if ((gs_UNSIGNED)m != N)
    {
        gs_info(gv, -2, "Matrix is not positive definite");
    }

    if (info > 0)
    {
        gs_info(gv, -2, "Computation of eigenvalues failed");
    }
}

void st_eigen(double *a, double *w, int *NN)
{
    gs_UNSIGNED N = (gs_UNSIGNED)(*NN);
    st_dsyevr(a, w, N);
}

static int dsyevr_vecs_rm(int n, double *A_rm, int lda_rm, double *w,
                          double *Z_rm, int ldZ_rm, int *m_out)
{
    int info = 0;

    int ldA = (n > 0) ? n : 1;
    double *A_cm = (double *)malloc((size_t)ldA * (size_t)n * sizeof(double));
    if (!A_cm)
        return -1010;

    rm_to_cm(A_rm, n, n, lda_rm, A_cm, ldA);

    char jobz = 'V';
    char range = 'A';
    char uplo = 'L';

    double vl = 0.0, vu = 0.0;
    int il = 0, iu = 0;
    double abstol = -1.0;

    int m = 0;

    int ldZ_cm = (n > 0) ? n : 1;
    double *Z_cm =
        (double *)malloc((size_t)ldZ_cm * (size_t)n * sizeof(double));
    int *isuppz = (int *)malloc((size_t)2 * (size_t)n * sizeof(int));
    if (!Z_cm || !isuppz)
    {
        free(isuppz);
        free(Z_cm);
        free(A_cm);
        return -1010;
    }

    int lwork = -1, liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;

    F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, A_cm, &ldA, &vl, &vu, &il, &iu,
                     &abstol, &m, w, Z_cm, &ldZ_cm, isuppz, &work_query, &lwork,
                     &iwork_query, &liwork, &info FCONE FCONE FCONE);

    if (info != 0)
    {
        if (m_out)
            *m_out = m;
        free(isuppz);
        free(Z_cm);
        free(A_cm);
        return info;
    }

    lwork = (int)work_query;
    liwork = (int)iwork_query;
    if (lwork < 1)
        lwork = 1;
    if (liwork < 1)
        liwork = 1;

    double *work = (double *)malloc((size_t)lwork * sizeof(double));
    int *iwork = (int *)malloc((size_t)liwork * sizeof(int));
    if (!work || !iwork)
    {
        free(iwork);
        free(work);
        free(isuppz);
        free(Z_cm);
        free(A_cm);
        if (m_out)
            *m_out = m;
        return -1010;
    }

    F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, A_cm, &ldA, &vl, &vu, &il, &iu,
                     &abstol, &m, w, Z_cm, &ldZ_cm, isuppz, work, &lwork, iwork,
                     &liwork, &info FCONE FCONE FCONE);

    if (info == 0)
    {
        cm_to_rm(Z_cm, n, n, ldZ_cm, Z_rm, ldZ_rm);
    }

    if (m_out)
        *m_out = m;

    free(iwork);
    free(work);
    free(isuppz);
    free(Z_cm);
    free(A_cm);

    return info;
}

void st_dsyevr_v(double *a, double *w, double *z, gs_UNSIGNED N)
{
    int n = (int)N;
    int lda_rm = (int)N;
    int ldz_rm = (int)N;

    int m = 0;
    int info = dsyevr_vecs_rm(n, a, lda_rm, w, z, ldz_rm, &m);

    if (info > 0)
    {
        gs_info(gv, -1, "Error in st_dsyevr_v");
    }
}

void st_dsyevr_v_R(double *a, double *w, double *z, int *nn)
{
    if (!a || !w || !z || !nn)
        return;
    st_dsyevr_v(a, w, z, (gs_UNSIGNED)(*nn));
}

static int dgesvd_rm_S(int m, int n, double *A_rm, int lda_rm, double *Svals,
                       double *U_rm, int ldu_rm, double *VT_rm, int ldvt_rm)
{
    int info = 0;

    if (m >= n && lda_rm == n && ldu_rm == n && ldvt_rm == n)
    {
        char jobu = 'S';
        char jobvt = 'S';

        int mm = n;
        int nn = m;

        int ldA = (n > 0) ? n : 1;
        int ldU = (n > 0) ? n : 1;
        int ldVT = (n > 0) ? n : 1;

        int lwork = -1;
        double work_query = 0.0;

        F77_CALL(dgesvd)(&jobu, &jobvt, &mm, &nn, A_rm, &ldA, Svals, VT_rm,
                         &ldU, U_rm, &ldVT, &work_query, &lwork,
                         &info FCONE FCONE);

        if (info != 0)
            return info;

        lwork = (int)work_query;
        if (lwork < 1)
            lwork = 1;

        double *work = (double *)malloc((size_t)lwork * sizeof(double));
        if (!work)
            return -1010;

        F77_CALL(dgesvd)(&jobu, &jobvt, &mm, &nn, A_rm, &ldA, Svals, VT_rm,
                         &ldU, U_rm, &ldVT, work, &lwork, &info FCONE FCONE);

        free(work);
        return info;
    }

    int ldA = (m > 0) ? m : 1;
    int ldU = (m > 0) ? m : 1;
    int ldVT = (n > 0) ? n : 1;

    double *A_cm = (double *)malloc((size_t)ldA * (size_t)n * sizeof(double));
    double *U_cm = (double *)malloc((size_t)ldU * (size_t)n * sizeof(double));
    double *VT_cm = (double *)malloc((size_t)ldVT * (size_t)n * sizeof(double));
    if (!A_cm || !U_cm || !VT_cm)
    {
        free(VT_cm);
        free(U_cm);
        free(A_cm);
        return -1010;
    }

    rm_to_cm(A_rm, m, n, lda_rm, A_cm, ldA);

    char jobu = 'S';
    char jobvt = 'S';

    int lwork = -1;
    double work_query = 0.0;

    F77_CALL(dgesvd)(&jobu, &jobvt, &m, &n, A_cm, &ldA, Svals, U_cm, &ldU,
                     VT_cm, &ldVT, &work_query, &lwork, &info FCONE FCONE);

    if (info != 0)
    {
        free(VT_cm);
        free(U_cm);
        free(A_cm);
        return info;
    }

    lwork = (int)work_query;
    if (lwork < 1)
        lwork = 1;

    double *work = (double *)malloc((size_t)lwork * sizeof(double));
    if (!work)
    {
        free(VT_cm);
        free(U_cm);
        free(A_cm);
        return -1010;
    }

    F77_CALL(dgesvd)(&jobu, &jobvt, &m, &n, A_cm, &ldA, Svals, U_cm, &ldU,
                     VT_cm, &ldVT, work, &lwork, &info FCONE FCONE);

    if (info == 0)
    {

        cm_to_rm(U_cm, m, n, ldU, U_rm, ldu_rm);

        cm_to_rm(VT_cm, n, n, ldVT, VT_rm, ldvt_rm);
    }

    free(work);
    free(VT_cm);
    free(U_cm);
    free(A_cm);

    return info;
}

void st_dgesvd_ginv(double *a, double *inva, gs_UNSIGNED nrow, gs_UNSIGNED ncol)
{
    int m = (int)nrow;
    int n = (int)ncol;

    int lda_rm = n;
    int ldu_rm = n;
    int ldvt_rm = n;

    double *u = (double *)malloc((size_t)m * (size_t)n * sizeof(double));
    double *s = (double *)malloc((size_t)n * sizeof(double));
    double *vt = (double *)malloc((size_t)n * (size_t)n * sizeof(double));

    double *tmp = (double *)malloc((size_t)m * (size_t)n * sizeof(double));

    if (!u || !s || !vt || !tmp)
    {
        free(tmp);
        free(vt);
        free(s);
        free(u);
        gs_info(gv, -1, "Out of memory in st_dgesvd_ginv");
        return;
    }

    int info = dgesvd_rm_S(m, n, a, lda_rm, s, u, ldu_rm, vt, ldvt_rm);
    if (info != 0)
    {
        gs_info(gv, -1, "Error in dgesvd (st_dgesvd_ginv)");
        free(tmp);
        free(vt);
        free(s);
        free(u);
        return;
    }

    const double tol = 1.0e-9;
    for (int i = 0; i < n; ++i)
    {
        double di = (s[i] > tol) ? (1.0 / s[i]) : 0.0;
        double *row = vt + (size_t)i * (size_t)n;
        for (int j = 0; j < n; ++j)
            row[j] *= di;
    }

    st_dgemm("N", "N", (gs_UNSIGNED)m, (gs_UNSIGNED)n, (gs_UNSIGNED)n,
             (gs_UNSIGNED)n, u, vt, tmp);

    for (int i = 0; i < n; ++i)
    {
        double *out_row = inva + (size_t)i * (size_t)m;
        for (int j = 0; j < m; ++j)
        {
            out_row[j] = tmp[(size_t)j * (size_t)n + (size_t)i];
        }
    }

    free(tmp);
    free(vt);
    free(s);
    free(u);
}

static int dgesdd_rm_S(int m, int n, double *A_rm, int lda_rm, double *Svals,
                       double *U_rm, int ldu_rm, double *VT_rm, int ldvt_rm)
{
    int info = 0;

    if (m >= n && lda_rm == n && ldu_rm == n && ldvt_rm == n)
    {
        char jobz = 'S';

        int mm = n;
        int nn = m;
        int ldA = (n > 0) ? n : 1;

        int ldU = (n > 0) ? n : 1;
        int ldVT = (n > 0) ? n : 1;

        int *iwork = (int *)malloc((size_t)8 * (size_t)n * sizeof(int));
        if (!iwork)
            return -1010;

        int lwork = -1;
        double work_query = 0.0;

        F77_CALL(dgesdd)(&jobz, &mm, &nn, A_rm, &ldA, Svals, VT_rm, &ldU, U_rm,
                         &ldVT, &work_query, &lwork, iwork, &info FCONE);

        if (info != 0)
        {
            free(iwork);
            return info;
        }

        lwork = (int)work_query;
        if (lwork < 1)
            lwork = 1;

        double *work = (double *)malloc((size_t)lwork * sizeof(double));
        if (!work)
        {
            free(iwork);
            return -1010;
        }

        F77_CALL(dgesdd)(&jobz, &mm, &nn, A_rm, &ldA, Svals, VT_rm, &ldU, U_rm,
                         &ldVT, work, &lwork, iwork, &info FCONE);

        free(work);
        free(iwork);
        return info;
    }

    return -6;
}

void st_dgesdd_ginv(double *a, double *inva, gs_UNSIGNED nrow, gs_UNSIGNED ncol)
{
    int m = (int)nrow;
    int n = (int)ncol;

    int lda_rm = n;
    int ldu_rm = n;
    int ldvt_rm = n;

    double *u = (double *)malloc((size_t)m * (size_t)n * sizeof(double));
    double *s = (double *)malloc((size_t)n * sizeof(double));
    double *vt = (double *)malloc((size_t)n * (size_t)n * sizeof(double));

    if (!u || !s || !vt)
    {
        free(vt);
        free(s);
        free(u);
        gs_info(gv, -1, "Out of memory in st_dgesdd_ginv");
        return;
    }

    int info = dgesdd_rm_S(m, n, a, lda_rm, s, u, ldu_rm, vt, ldvt_rm);
    if (info != 0)
    {
        gs_info(gv, -1, "Error in dgesdd (st_dgesdd_ginv)");
        free(vt);
        free(s);
        free(u);
        return;
    }

    const double tol = 1.0e-9;
    for (int i = 0; i < n; ++i)
    {
        double di = (s[i] > tol) ? (1.0 / s[i]) : 0.0;
        double *row = vt + (size_t)i * (size_t)n;
        for (int j = 0; j < n; ++j)
            row[j] *= di;
    }

    st_dgemm("T", "T", (gs_UNSIGNED)n, (gs_UNSIGNED)n, (gs_UNSIGNED)m,
             (gs_UNSIGNED)n, vt, u, inva);

    free(vt);
    free(s);
    free(u);
}

void st_ginv2(double *a, double *inva, int *nrow, int *ncol)
{
    gs_UNSIGNED nr = (gs_UNSIGNED)*nrow;
    gs_UNSIGNED nc = (gs_UNSIGNED)*ncol;

    st_dgesvd_ginv(a, inva, nr, nc);
}

void st_ginv(double *a, double *inva, int *nrow, int *ncol)
{
    gs_UNSIGNED nr = (gs_UNSIGNED)*nrow;
    gs_UNSIGNED nc = (gs_UNSIGNED)*ncol;

    st_dgesdd_ginv(a, inva, nr, nc);
}

static int st_dpotrf_rm_sym(double *A_rm, int n)
{
    int info = 0;
    char uplo = 'U';
    int lda = (n > 0) ? n : 1;
    F77_CALL(dpotrf)(&uplo, &n, A_rm, &lda, &info FCONE);
    return info;
}

static int st_dpotri_rm_sym(double *A_rm, int n)
{
    int info = 0;
    char uplo = 'U';
    int lda = (n > 0) ? n : 1;
    F77_CALL(dpotri)(&uplo, &n, A_rm, &lda, &info FCONE);
    return info;
}

static void symmetrize_from_lower_rm(double *A_rm, int n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {

            A_rm[(size_t)i * (size_t)n + (size_t)j] =
                A_rm[(size_t)j * (size_t)n + (size_t)i];
        }
    }
}

SEXP call_st_dgemm_R(SEXP T1_, SEXP T2_, SEXP A_, SEXP B_, SEXP nra_, SEXP nca_,
                     SEXP nrb_, SEXP ncb_)
{
    const char *T1 = CHAR(STRING_ELT(T1_, 0));
    const char *T2 = CHAR(STRING_ELT(T2_, 0));

    unsigned int nra = (unsigned int)asInteger(nra_);
    unsigned int nca = (unsigned int)asInteger(nca_);
    unsigned int nrb = (unsigned int)asInteger(nrb_);
    unsigned int ncb = (unsigned int)asInteger(ncb_);

    int transA = (strcmp(T1, "N") == 0) ? 0 : 1;
    int transB = (strcmp(T2, "N") == 0) ? 0 : 1;
    unsigned int m = transA ? nca : nra;
    unsigned int kA = transA ? nra : nca;
    unsigned int kB = transB ? ncb : nrb;
    unsigned int n = transB ? nrb : ncb;

    if (kA != kB)
        error("nonconformable: kA=%u, kB=%u", kA, kB);

    if (!isReal(A_) || !isReal(B_))
        error("A and B must be double");
    if ((unsigned int)XLENGTH(A_) < nra * nca)
        error("A too short");
    if ((unsigned int)XLENGTH(B_) < nrb * ncb)
        error("B too short");

    SEXP C_ = PROTECT(allocVector(REALSXP, (R_xlen_t)m * (R_xlen_t)n));

    char t1buf[2] = {T1[0], '\0'};
    char t2buf[2] = {T2[0], '\0'};

    st_dgemm(t1buf, t2buf, nra, nca, nrb, ncb, REAL(A_), REAL(B_), REAL(C_));

    UNPROTECT(1);
    return C_;
}

void st_dgetrf_R(double *a, int *nn, int *indx, double *d)
{
    if (!a || !nn || !indx || !d)
        return;

    gs_UNSIGNED n = (gs_UNSIGNED)(*nn);

    gs_UNSIGNED *piv = (gs_UNSIGNED *)malloc((size_t)n * sizeof(gs_UNSIGNED));
    if (!piv)
        return;

    st_dgetrf(a, n, piv, d);

    for (gs_UNSIGNED i = 0; i < n; ++i)
        indx[i] = (int)piv[i];

    free(piv);
}

void st_dgetrs_R(double *a, int *nn, int *indx, double *b)
{
    if (!a || !nn || !indx || !b)
        return;

    gs_UNSIGNED n = (gs_UNSIGNED)(*nn);

    gs_UNSIGNED *piv = (gs_UNSIGNED *)malloc((size_t)n * sizeof(gs_UNSIGNED));
    if (!piv)
        return;

    for (gs_UNSIGNED i = 0; i < n; ++i)
        piv[i] = (gs_UNSIGNED)indx[i];

    st_dgetrs(a, n, piv, b);

    free(piv);
}

void st_dgetri_R(double *a, int *nn, int *indx, double *inv)
{
    if (!a || !nn || !indx || !inv)
        return;

    gs_UNSIGNED n = (gs_UNSIGNED)(*nn);

    gs_UNSIGNED *piv = (gs_UNSIGNED *)malloc((size_t)n * sizeof(gs_UNSIGNED));
    if (!piv)
        return;

    for (gs_UNSIGNED i = 0; i < n; ++i)
        piv[i] = (gs_UNSIGNED)indx[i];

    st_dgetri(a, n, piv, inv);

    free(piv);
}

void st_dgetri_inv_R(double *a, double *inva, int *n)
{
    if (!a || !inva || !n)
        return;
    st_dgetri_inv(a, inva, (gs_UNSIGNED)(*n));
}

void st_inv_R(double *a, double *inva, int *nn)
{
    if (!a || !inva || !nn)
        return;
    st_inv(a, inva, nn);
}

void st_dgesvd_ginv_R(double *a, double *inva, int *nrow, int *ncol)
{
    if (!a || !inva || !nrow || !ncol)
        return;
    st_dgesvd_ginv(a, inva, (gs_UNSIGNED)(*nrow), (gs_UNSIGNED)(*ncol));
}
