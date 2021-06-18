/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : linalg.hpp                                    |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Linear algebra using templatized BLAS and     |
 |  `---'  |                Lapack interface for optimized BLAS library   |
 +========================================================================+
 */

#pragma once
#define CASTOR_LINALG_HPP
#include <castor/matrix.hpp>
#ifdef __CLAPACK_H
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

namespace castor
{

//==========================================================================//
//                               CBLAS                                      //
//==========================================================================//
//==========================================================================
// [tgemm]
/// Matrix-matrix operations C = alpha*A*B + beta*C.
///
/// TGEMM performs one of the matrix-matrix operations using BLAS interface
///    C = alpha*A*B + beta*C,
/// alpha and beta are scalars of any type, and A, B and C are matrices
/// of same type X, with A an m by k matrix, B a k by n matrix
/// and C an m by n matrix.
template<typename R, typename S>
void tgemm(R alpha, matrix<float>const& A, matrix<float>const& B, S beta, matrix<float>& C)
{
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)size(C,1), (int)size(C,2), (int)size(A,2),
                (float)alpha, &A(0), (int)size(A,2), &B(0), (int)size(B,2), (float)beta, &C(0), (int) size(C,2));
}
template<typename R, typename S>
void tgemm(R alpha, matrix<double>const& A, matrix<double>const& B, S beta, matrix<double>& C)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)size(C,1), (int)size(C,2), (int)size(A,2),
                (double)alpha, &A(0), (int)size(A,2), &B(0), (int)size(B,2), (double)beta, &C(0), (int) size(C,2));
}
template<typename R, typename S>
void tgemm(R alpha, matrix<std::complex<float>>const& A,
    matrix<std::complex<float>>const& B, S beta, matrix<std::complex<float>>& C)
{
    std::complex<float> alpha_c = alpha, beta_c = beta;
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)size(C,1), (int)size(C,2), (int)size(A,2),
                &alpha_c, &A(0), (int)size(A,2), &B(0), (int)size(B,2), &beta_c, &C(0), (int) size(C,2));
}
template<typename R, typename S>
void tgemm(R alpha, matrix<std::complex<double>>const& A,
    matrix<std::complex<double>>const& B, S beta, matrix<std::complex<double>>& C)
{
    std::complex<double> alpha_z = alpha, beta_z = beta;
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)size(C,1), (int)size(C,2), (int)size(A,2),
                &alpha_z, &A(0), (int)size(A,2), &B(0), (int)size(B,2), &beta_z, &C(0), (int) size(C,2));
}


//==========================================================================//
//                            LAPACK INTERFACE                              //
//==========================================================================//
//==========================================================================
// [clpk]
extern "C"
{
typedef struct{float r,i;} __CLPK_complex;
}
class clpk : public __CLPK_complex
{
public:
    clpk() {r=0; i=0;};
    clpk(float v, float w=0) {r=v; i=w;};
    clpk(std::complex<float> v) {r=v.real(); i=v.imag();};
    operator std::complex<float>() const {return std::complex<float>(r,i);};
};
std::ostream& operator<<(std::ostream& flux, clpk const& c)
{
    flux << "(" << c.r << "," << c.i << ")"; return flux;
}

//==========================================================================
// [zlpk]
extern "C"
{
typedef struct{double r,i;} __CLPK_doublecomplex;
}
class zlpk : public __CLPK_doublecomplex
{
public:
    zlpk() {r=0; i=0;};
    zlpk(double v, double w=0) {r=v; i=w;};
    zlpk(std::complex<double> v) {r=v.real(); i=v.imag();};
    operator std::complex<double>() const {return std::complex<double>(r,i);};
};
std::ostream& operator<<(std::ostream& flux, zlpk const& c)
{
    flux << "(" << c.r << "," << c.i << ")"; return flux;
}

//==========================================================================
// [extern]
extern "C"
{
// S/D/C/ZGESV
void sgesv_(int *N, int *NRHS, float *A, int *LDA, int *IPIV, float *B, int *LDB, int *INFO);
void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
void cgesv_(int* N, int* NRHS, clpk* A, int* LDA, int* IPIV, clpk* B, int* LDB, int* INFO);
void zgesv_(int* N, int* NRHS, zlpk* A, int* LDA, int* IPIV, zlpk* B, int* LDB, int* INFO);

// S/D/C/ZGESVD
void sgesvd_(char* JOBU, char* JOBVT, int *M, int *N, float *A, int *LDA, float *S, float *U, int *LDU, float *VT, int *LDVT, float *WORK, int *LWORK, int *INFO);
void dgesvd_(char* JOBU, char* JOBVT, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *INFO);
void cgesvd_(char* JOBU, char* JOBVT, int *M, int *N, clpk *A, int *LDA, float *S, clpk *U, int *LDU, clpk *VT, int *LDVT, clpk *WORK, int *LWORK, float *RWORK, int *INFO);
void zgesvd_(char* JOBU, char* JOBVT, int *M, int *N, zlpk *A, int *LDA, double *S, zlpk *U, int *LDU, zlpk *VT, int *LDVT, zlpk *WORK, int *LWORK, double *RWORK, int *INFO);

// S/D/C/ZGESDD
void sgesdd_(char* JOBZ, int *M, int *N, float *A, int *LDA, float *S, float *U, int *LDU, float *VT, int *LDVT, float *WORK, int *LWORK, int *IWORK, int *INFO);
void dgesdd_(char* JOBZ, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *IWORK, int *INFO);
void cgesdd_(char* JOBZ, int *M, int *N, clpk *A, int *LDA, float *S, clpk *U, int *LDU, clpk *VT, int *LDVT, clpk *WORK, int *LWORK, float *RWORK, int *IWORK, int *INFO);
void zgesdd_(char* JOBZ, int *M, int *N, zlpk *A, int *LDA, double *S, zlpk *U, int *LDU, zlpk *VT, int *LDVT, zlpk *WORK, int *LWORK, double *RWORK, int *IWORK, int *INFO);

// S/D/C/ZGEQRF
void sgeqrf_(int *M, int *N, float *A, int *LDA, float *TAU, float *WORK, int *LWORK, int *INFO);
void dgeqrf_(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void cgeqrf_(int *M, int *N, clpk *A, int *LDA, clpk *TAU, clpk *WORK, int *LWORK, int *INFO);
void zgeqrf_(int *M, int *N, zlpk *A, int *LDA, zlpk *TAU, zlpk *WORK, int *LWORK, int *INFO);

// S/D/C/ZGETRF
void sgetrf_(int* M, int* N, float* A, int* LDA, int* IPIV, int* INFO);
void dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);
void cgetrf_(int* M, int* N, clpk* A, int* LDA, int* IPIV, int* INFO);
void zgetrf_(int* M, int* N, zlpk* A, int* LDA, int* IPIV, int* INFO);

// S/D/C/ZGELS
void sgels_(char *TRANS, int *M, int *N, int *NRHS, float *A, int *LDA, float *B, int *LDB, float *WORK, int *LWORK, int *INFO);
void dgels_(char *TRANS, int *M, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *WORK, int *LWORK, int *INFO);
void cgels_(char *TRANS, int *M, int *N, int *NRHS, clpk *A, int *LDA, clpk *B, int *LDB, clpk *WORK, int *LWORK, int *INFO);
void zgels_(char *TRANS, int *M, int *N, int *NRHS, zlpk *A, int *LDA, zlpk *B, int *LDB, zlpk *WORK, int *LWORK, int *INFO);

// S/DORGQR
void sorgqr_(int *M, int *N, int *K, float *A, int *LDA, float *TAU, float *WORK, int *LWORK, int *INFO);
void dorgqr_(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);

// // C/ZUNGQR
void cungqr_(int *M, int *N, int *K, clpk *A, int *LDA, clpk *TAU, clpk *WORK, int *LWORK, int *INFO);
void zungqr_(int *M, int *N, int *K, zlpk *A, int *LDA, zlpk *TAU, zlpk *WORK, int *LWORK, int *INFO);

// C/ZGEEV
void cgeev_(char *JOBVL, char *JOBVR, int *N, clpk *A, int *LDA, clpk *W, clpk *VL, int *LDVL, clpk *VR, int *LDVR, clpk *WORK, int *LWORK, float *RWORK, int *INFO);
void zgeev_(char *JOBVL, char *JOBVR, int *N, zlpk *A, int *LDA, zlpk *W, zlpk *VL, int *LDVL, zlpk *VR, int *LDVR, zlpk *WORK, int *LWORK, double *RWORK, int *INFO);

// S/D/C/ZGGEV
void cggev_(char *JOBVL, char *JOBVR, int *N, clpk *A, int *LDA, clpk *B, int *LDB, clpk *ALPHA, clpk *BETA, clpk *VL, int *LDVL, clpk *VR, int *LDVR, clpk *WORK, int *LWORK, float *RWORK, int *INFO);
void zggev_(char *JOBVL, char *JOBVR, int *N, zlpk *A, int *LDA, zlpk *B, int *LDB, zlpk *ALPHA, zlpk *BETA, zlpk *VL, int *LDVL, zlpk *VR, int *LDVR, zlpk *WORK, int *LWORK, double *RWORK, int *INFO);
}

//==========================================================================
// [mat2lpk]
template<typename T, typename S>
std::vector<T> mat2lpk(matrix<S>const& A, std::size_t L)
{
    if (L<numel(A) || L%size(A,2)!=0)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    std::size_t m=size(A,1), n=size(A,2), M=L/n;
    std::vector<T> V(L,(T)0);
    for (std::size_t j=0; j<n; ++j)
    {
        for (std::size_t i=0; i<m; ++i)
        {
            V[j*M+i] = A(i,j);
        }
    }
    return V;
}

//==========================================================================
// [lpk2mat]
template<typename T, typename S>
matrix<T> lpk2mat(std::vector<S>const& V, std::size_t m, std::size_t n)
{
    if (V.size()<m*n || V.size()%n!=0)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    matrix<T> A(m,n,(T)0);
    std::size_t M = V.size()/n;
    for (std::size_t i=0; i<m; ++i)
    {
        for (std::size_t j=0; j<n; ++j)
        {
            A(i,j) = V[j*M+i];
        }
    }
    return A;
}


//==========================================================================//
//                          LAPACK FUNCTIONS                                //
//==========================================================================//
//========================================================================
// [tgeev]
/// Eigenvalues and the left and/or right eigenvectors.
///
/// TGEEV computes for an N-by-N nonsymmetric matrix A of type T, the
/// eigenvalues and, optionally, the left and/or right eigenvectors.
///
/// The left eigenvectors of A are the same as the right eigenvectors of A^t.
/// If u(j) and v(j) are the left and right eigenvectors, respectively,
/// corresponding to the eigenvalue lambda(j):
///    (u(j)^t)*A = lambda(j)*(u(j)^t),
///    A*v(j) = lambda(j) * v(j).
///
/// The computed eigenvectors are normalized to have Euclidean
/// norm equal to 1 and largest component real.
void xgeev(char& jobl, char& jobr, int& n, std::vector<clpk>& A, std::vector<clpk>& E, std::vector<clpk>& V,
           clpk& wkopt, int& lwork, int& info)
{
    std::vector<float> rwork(2*n);
    cgeev_(&jobl, &jobr, &n, &A[0], &n, &E[0], &V[0], &n, &V[0], &n, &wkopt, &lwork, &rwork[0], &info);
    lwork = (int) wkopt.r;
    std::vector<clpk> work(lwork);
    cgeev_(&jobl, &jobr, &n, &A[0], &n, &E[0], &V[0], &n, &V[0], &n, &work[0], &lwork, &rwork[0], &info);
}
void xgeev(char& jobl, char& jobr, int& n, std::vector<zlpk>& A, std::vector<zlpk>& E, std::vector<zlpk>& V,
           zlpk& wkopt, int& lwork, int& info)
{
    std::vector<double> rwork(2*n);
    zgeev_(&jobl, &jobr, &n, &A[0], &n, &E[0], &V[0], &n, &V[0], &n, &wkopt, &lwork, &rwork[0], &info);
    lwork = (int) wkopt.r;
    std::vector<zlpk> work(lwork);
    zgeev_(&jobl, &jobr, &n, &A[0], &n, &E[0], &V[0], &n, &V[0], &n, &work[0], &lwork, &rwork[0], &info);
}
template<typename T>
void tgeev(std::string typ, int& n, std::vector<T>& A, std::vector<T>& E, std::vector<T>& V)
{
    if (A.size()!=n*n)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix must be square.");
    }
    char jobl = 'N';
    char jobr = 'N';
    if (typ=="none") {}
    else if (typ=="left") {jobl='V'; V.resize(n*n);}
    else if (typ=="right") {jobr='V'; V.resize(n*n);}
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unknown eig type, only use \"none\", \"left\" or \"right\" for 2nd argument.");
    }
    int lwork = -1;
    int info;
    T wkopt;
    xgeev(jobl, jobr, n, A, E, V, wkopt, lwork, info);
    if (info<0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Matrix argument(s) had illegal value(s).");
    }
    else if (info>0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,
                "The QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed.");
    }
}

// [tggev]
/// Generalized eigenvalues and the left and/or right eigenvectors.
///
/// TGEEV computes for an N-by-N nonsymmetric matrix A of type T, the
/// eigenvalues and, optionally, the left and/or right eigenvectors.
void xggev(char &jobl, char &jobr, int &n, std::vector<clpk> &A, std::vector<clpk> &B, std::vector<clpk> &alpha, std::vector<clpk> &beta, std::vector<clpk> &V, clpk &wkopt, int &lwork, int &info)
{
    std::vector<float> rwork(8*n);
    cggev_(&jobl, &jobr, &n, &A[0], &n, &B[0], &n, &alpha[0], &beta[0], &V[0], &n, &V[0], &n, &wkopt, &lwork, &rwork[0], &info); // fetch workspace size
    lwork = static_cast<int>(wkopt.r);
    std::vector<clpk> work(lwork);
    cggev_(&jobl, &jobr, &n, &A[0], &n, &B[0], &n, &alpha[0], &beta[0], &V[0], &n, &V[0], &n, &work[0], &lwork, &rwork[0], &info); // perform computation
}
void xggev(char &jobl, char &jobr, int &n, std::vector<zlpk> &A, std::vector<zlpk> &B, std::vector<zlpk> &alpha, std::vector<zlpk> &beta, std::vector<zlpk> &V, zlpk &wkopt, int &lwork, int &info)
{
    std::vector<double> rwork(8*n);
    zggev_(&jobl, &jobr, &n, &A[0], &n, &B[0], &n, &alpha[0], &beta[0], &V[0], &n, &V[0], &n, &wkopt, &lwork, &rwork[0], &info); // fetch workspace size
    lwork = static_cast<int>(wkopt.r);
    std::vector<zlpk> work(lwork);
    zggev_(&jobl, &jobr, &n, &A[0], &n, &B[0], &n, &alpha[0], &beta[0], &V[0], &n, &V[0], &n, &work[0], &lwork, &rwork[0], &info); // perform computation
}
template<typename T>
void tggev(std::string typ, int &n, std::vector<T> &A, std::vector<T> &B, std::vector<T> &E, std::vector<T> &V)
{
    if(A.size() != n*n)
    {
        error(__FILE__, __LINE__, __FUNCTION__, "Matrix A must be square.");
    }
    if(B.size() != n*n)
    {
        error(__FILE__, __LINE__, __FUNCTION__, "Matrix B must have the same size as A.");
    }
    char jobl = 'N';
    char jobr = 'N';
    if(typ == "none") {}
    else if(typ == "left") {jobl = 'V'; V.resize(n*n);}
    else if(typ == "right"){jobr = 'V'; V.resize(n*n);}
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__, "Unknown eig type, only use \"none\", \"left\" or \"right\" for 2nd argument.");
    }
    int lwork = -1;
    int info  = 0;
    T wkopt;
    std::vector<T> alpha(n), beta(n);
    xggev(jobl, jobr, n, A, B, alpha, beta, V, wkopt, lwork, info);
    // compute the eigenvalues
    for(std::size_t ii=0; ii<n; ++ii)
    {
        auto numr = (alpha[ii].r * beta[ii].r  + alpha[ii].i * beta[ii].i);
        auto numi = (beta[ii].r  * alpha[ii].i - alpha[ii].r * beta[ii].i);
        auto den  = beta[ii].r   * beta[ii].r  + beta[ii].i  * beta[ii].i;
        E[ii].r  = numr/den;
        E[ii].i  = numi/den;
    }
    if(info < 0)
    {
        warning(__FILE__,__LINE__,__FUNCTION__,"Matrix argument(s) had illegal value(s).");
    }else if(info > 0)
    {
        warning(__FILE__,__LINE__,__FUNCTION__,"The algorithm failed to compute the eigenvalues.");
    }
}

//========================================================================
// [tgels]
/// Solve overdetermined or underdetermined linear systems.
///
/// TGELS solves overdetermined or underdetermined linear systems involving
/// an M-by-N matrix A of type T, or its transpose, using a QR or LQ
/// factorization of A. It is assumed that A has full rank.
///
/// If m >= n, it find the least squares solution of an overdetermined system,
/// solving the least squares problem:
///    minimize || B - A*X ||.
///
/// If m < n, it find the minimum norm solution of an underdetermined system
///    A*X = B.
///
/// Several right hand side vectors b and solution vectors x can be handled
/// in a single call; they are stored as the columns of the M-by-NRHS right
/// hand side matrix B and the N-by-NRHS solution matrix X.
void xgels(char& t, int& m, int& n, int& nrhs, std::vector<float>& A,
           std::vector<float>& B, int& l, float& wkopt, int& lwork, int& info)
{
    sgels_(&t, &m, &n, &nrhs, &A[0], &m, &B[0], &l, &wkopt, &lwork, &info);
    lwork = (int) wkopt;
    std::vector<float> work(lwork);
    sgels_(&t, &m, &n, &nrhs, &A[0], &m, &B[0], &l, &work[0], &lwork, &info);
}
void xgels(char& t, int& m, int& n, int& nrhs, std::vector<double>& A,
           std::vector<double>& B, int& l, double& wkopt, int& lwork, int& info)
{
    dgels_(&t, &m, &n, &nrhs, &A[0], &m, &B[0], &l, &wkopt, &lwork, &info);
    lwork = (int) wkopt;
    std::vector<double> work(lwork);
    dgels_(&t, &m, &n, &nrhs, &A[0], &m, &B[0], &l, &work[0], &lwork, &info);
}
void xgels(char& t, int& m, int& n, int& nrhs, std::vector<clpk>& A,
           std::vector<clpk>& B, int& l, clpk& wkopt, int& lwork, int& info)
{
    cgels_(&t, &m, &n, &nrhs, &A[0], &m, &B[0], &l, &wkopt, &lwork, &info);
    lwork = (int) wkopt.r;
    std::vector<clpk> work(lwork);
    cgels_(&t, &m, &n, &nrhs, &A[0], &m, &B[0], &l, &work[0], &lwork, &info);
}
void xgels(char& t, int& m, int& n, int& nrhs, std::vector<zlpk>& A,
           std::vector<zlpk>& B, int& l, zlpk& wkopt, int& lwork, int& info)
{
    zgels_(&t, &m, &n, &nrhs, &A[0], &m, &B[0], &l, &wkopt, &lwork, &info);
    lwork = (int) wkopt.r;
    std::vector<zlpk> work(lwork);
    zgels_(&t, &m, &n, &nrhs, &A[0], &m, &B[0], &l, &work[0], &lwork, &info);
}
template<typename T>
void tgels(int& m, int& n, int& nrhs, std::vector<T>& A, std::vector<T>& B)
{
    if ((A.size()!=m*n) || (B.size()!=std::max(m,n)*nrhs))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    char t = 'N';
    int l  = std::max<int>(m,n);
    T wkopt;
    int lwork = -1;
    int info;
    xgels(t, m, n, nrhs, A, B, l, wkopt, lwork, info);
    if (info<0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Matrix argument(s) had illegal value(s).");
    }
    else if (info>0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,
                "Matrix is not full rank, the least squares solution could not be computed.");
    }
}

//========================================================================
// [tgeqrf]
/// Computes a QR factorization of a M-by-N matrix A of type T (A=Q*R).
void xgeqrf(int&m, int& n, int& l, std::vector<float>&A, std::vector<float>& tau,
            float& wkopt, int& lwork, int& info)
{
    sgeqrf_(&m, &n, &A[0], &m, &tau[0], &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    std::vector<float> work(lwork);
    sgeqrf_(&m, &n, &A[0], &m, &tau[0], &work[0], &lwork, &info);
}
void xorgqr(int&m, int& l, std::vector<float>&A, std::vector<float>& tau,
            float& wkopt, int& lwork, int& info)
{
    sorgqr_(&m, &l, &l, &A[0], &m, &tau[0], &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    std::vector<float> work(lwork);
    sorgqr_(&m, &l, &l, &A[0], &m, &tau[0], &work[0], &lwork, &info);
}
void xgeqrf(int&m, int& n, int& l, std::vector<double>&A, std::vector<double>& tau,
            double& wkopt, int& lwork, int& info)
{
    dgeqrf_(&m, &n, &A[0], &m, &tau[0], &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    std::vector<double> work(lwork);
    dgeqrf_(&m, &n, &A[0], &m, &tau[0], &work[0], &lwork, &info);
}
void xorgqr(int&m, int& l, std::vector<double>&A, std::vector<double>& tau,
            double& wkopt, int& lwork, int& info)
{
    dorgqr_(&m, &l, &l, &A[0], &m, &tau[0], &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    std::vector<double> work(lwork);
    dorgqr_(&m, &l, &l, &A[0], &m, &tau[0], &work[0], &lwork, &info);
}
void xgeqrf(int&m, int& n, int& l, std::vector<clpk>&A, std::vector<clpk>& tau,
            clpk& wkopt, int& lwork, int& info)
{
    cgeqrf_(&m, &n, &A[0], &m, &tau[0], &wkopt, &lwork, &info);
    lwork = (int)wkopt.r;
    std::vector<clpk> work(lwork);
    cgeqrf_(&m, &n, &A[0], &m, &tau[0], &work[0], &lwork, &info);
}
void xorgqr(int&m, int& l, std::vector<clpk>&A, std::vector<clpk>& tau,
            clpk& wkopt, int& lwork, int& info)
{
    cungqr_(&m, &l, &l, &A[0], &m, &tau[0], &wkopt, &lwork, &info);
    lwork = (int)wkopt.r;
    std::vector<clpk> work(lwork);
    cungqr_(&m, &l, &l, &A[0], &m, &tau[0], &work[0], &lwork, &info);
}
void xgeqrf(int&m, int& n, int& l, std::vector<zlpk>&A, std::vector<zlpk>& tau,
            zlpk& wkopt, int& lwork, int& info)
{
    zgeqrf_(&m, &n, &A[0], &m, &tau[0], &wkopt, &lwork, &info);
    lwork = (int)wkopt.r;
    std::vector<zlpk> work(lwork);
    zgeqrf_(&m, &n, &A[0], &m, &tau[0], &work[0], &lwork, &info);
}
void xorgqr(int&m, int& l, std::vector<zlpk>&A, std::vector<zlpk>& tau,
            zlpk& wkopt, int& lwork, int& info)
{
    zungqr_(&m, &l, &l, &A[0], &m, &tau[0], &wkopt, &lwork, &info);
    lwork = (int)wkopt.r;
    std::vector<zlpk> work(lwork);
    zungqr_(&m, &l, &l, &A[0], &m, &tau[0], &work[0], &lwork, &info);
}
template<typename T>
void tgeqrf(int& m, int& n, std::vector<T>& A, std::vector<T>& R)
{
    int l = std::min<int>(m,n);
    if ((A.size()!=m*n) || (R.size()!=l*n))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    std::vector<T> tau(l);
    T wkopt;
    int lwork = -1;
    int info;
    xgeqrf(m, n, l, A, tau, wkopt, lwork, info);
    if (info<0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Matrix argument(s) had illegal value(s).");
    }
    for (int i=0; i<l; ++i)
    {
        for (int j=i; j<n; ++j)
        {
            R[j*l+i] = A[j*m+i];
        }
    }
    A.resize(m*l);
    lwork = -1;
    xorgqr(m, l, A, tau, wkopt, lwork, info);
    if (info<0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Matrix argument(s) had illegal value(s).");
    }
}

//========================================================================
// [tgesdd]
/// Singular value decomposition of a rectangular matrix A.
///
/// TGESDD computes the singular value decomposition (SVD) of a m-by-n
/// matrix A of type T, optionally computing the left and/or right singular
/// vectors. If singular vectors are desired, it uses a divide and conquer
/// algorithm.
///
/// The SVD is written as :
///    A = U*SIGMA*VT
/// where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
/// diagonal elements, U is an m-by-m orthogonal matrix and VT (V transposed)
/// is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
/// are the singular values of A; they are real and non-negative, and are
/// returned in descending order. The first min(m, n) columns of U and V are
/// the left and right singular vectors of A.
///
/// Note that the routine returns VT, not V.
void xgesdd(char& job, int&m, int& n, int& l, std::vector<float>& A,
            std::vector<float>& S, std::vector<float>& U, std::vector<float>& V,
            float& wkopt, int& lwork, std::vector<int>& iwork, int& info)
{
    sgesdd_(&job, &m, &n, &A[0], &m, &S[0], &U[0], &m, &V[0], &l, &wkopt, &lwork, &iwork[0], &info);
    lwork = (int)wkopt;
    std::vector<float> work(lwork);
    sgesdd_(&job, &m, &n, &A[0], &m, &S[0], &U[0], &m, &V[0], &l, &work[0], &lwork, &iwork[0], &info);
}
void xgesdd(char& job, int&m, int& n, int& l, std::vector<double>& A,
            std::vector<double>& S, std::vector<double>& U, std::vector<double>& V,
            double& wkopt, int& lwork, std::vector<int>& iwork, int& info)
{
    dgesdd_(&job, &m, &n, &A[0], &m, &S[0], &U[0], &m, &V[0], &l, &wkopt, &lwork, &iwork[0], &info);
    lwork = (int)wkopt;
    std::vector<double> work(lwork);
    dgesdd_(&job, &m, &n, &A[0], &m, &S[0], &U[0], &m, &V[0], &l, &work[0], &lwork, &iwork[0], &info);
}
void xgesdd(char& job, int&m, int& n, int& l, std::vector<clpk>& A,
            std::vector<float>& S, std::vector<clpk>& U, std::vector<clpk>& V,
            clpk& wkopt, int& lwork, std::vector<int>& iwork, int& info)
{
    std::vector<float> rwork;
    int mn=std::min<int>(m,n), mx=std::max<int>(m,n);
    if (job=='N') {rwork.resize(7*mn);}
    else {rwork.resize(std::max(5*mn*mn+5*mn,2*mx*mn+2*mn*mn+mn));}
    cgesdd_(&job, &m, &n, &A[0], &m, &S[0], &U[0], &m, &V[0], &l, &wkopt, &lwork, &rwork[0], &iwork[0], &info);
    lwork = (int)wkopt.r;
    std::vector<clpk> work(lwork);
    cgesdd_(&job, &m, &n, &A[0], &m, &S[0], &U[0], &m, &V[0], &l, &work[0], &lwork, &rwork[0], &iwork[0], &info);
}
void xgesdd(char& job, int&m, int& n, int& l, std::vector<zlpk>& A,
            std::vector<double>& S, std::vector<zlpk>& U, std::vector<zlpk>& V,
            zlpk& wkopt, int& lwork, std::vector<int>& iwork, int& info)
{
    std::vector<double> rwork;
    int mn=std::min<int>(m,n), mx=std::max<int>(m,n);
    if (job=='N') {rwork.resize(7*mn);}
    else {rwork.resize(std::max(5*mn*mn+5*mn,2*mx*mn+2*mn*mn+mn));}
    zgesdd_(&job, &m, &n, &A[0], &m, &S[0], &U[0], &m, &V[0], &l, &wkopt, &lwork, &rwork[0], &iwork[0], &info);
    lwork = (int)wkopt.r;
    std::vector<zlpk> work(lwork);
    zgesdd_(&job, &m, &n, &A[0], &m, &S[0], &U[0], &m, &V[0], &l, &work[0], &lwork, &rwork[0], &iwork[0], &info);
}
template<typename T, typename T2>
void tgesdd(std::string typ, int&m, int& n, std::vector<T>& A,
            std::vector<T2>& S, std::vector<T>& U, std::vector<T>& V)
{
    int l = std::min<int>(m,n);
    if (A.size()!=m*n || S.size()!=l)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimension must agree.");
    }
    char job;
    if (typ=="none") {job='N';}
    else if (typ=="vect")
    {
        job='S';
        U.resize(m*l);
        V.resize(l*n);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unknown svd type, only use \"none\" or \"vect\" for 2nd argument.");
    }
    int lwork = -1;
    int info;
    T wkopt;
    std::vector<int> iwork(8*l);
    xgesdd(job, m, n, l, A, S, U, V, wkopt, lwork, iwork, info);
    if (info<0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Matrix argument(s) had illegal value(s).");
    }
    else if (info>0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,
                "xBDSDC did not converge, updating process failed.");
    }
}

//========================================================================
// [tgesv]
/// Solve linear system equations  A*X=B.
///
/// TGESV computes the solution to a system of linear equations
///    A*X = B,
/// where A is an N-by-N matrix and X and B are N-by-NRHS matrices,
/// both of type T.
///
/// The LU decomposition with partial pivoting and row interchanges
/// is used to factor A as:
///    A = P * L * U,
/// where P is a permutation matrix, L is unit lower triangular,
/// and U is upper triangular. The factored form of A is then used
/// to solve the system of equations A * X = B.
void xgesv(int& n, int& nrhs, std::vector<float>& A, std::vector<float>& B,
           std::vector<int>& ipiv, int& info)
{
    sgesv_(&n, &nrhs, &A[0], &n, &ipiv[0], &B[0], &n, &info);
}
void xgesv(int& n, int& nrhs, std::vector<double>& A, std::vector<double>& B,
           std::vector<int>& ipiv, int& info)
{
    dgesv_(&n, &nrhs, &A[0], &n, &ipiv[0], &B[0], &n, &info);
}
void xgesv(int& n, int& nrhs, std::vector<clpk>& A, std::vector<clpk>& B,
           std::vector<int>& ipiv, int& info)
{
    cgesv_(&n, &nrhs, &A[0], &n, &ipiv[0], &B[0], &n, &info);
}
void xgesv(int& n, int& nrhs, std::vector<zlpk>& A, std::vector<zlpk>& B,
           std::vector<int>& ipiv, int& info)
{
    zgesv_(&n, &nrhs, &A[0], &n, &ipiv[0], &B[0], &n, &info);
}
template<typename T>
void tgesv(int& n, int& nrhs, std::vector<T>& A, std::vector<T>& B)
{
    if (A.size()!=n*n)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix must be square.");
    }
    else if (B.size()!=n*nrhs)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    std::vector<int> ipiv(n);
    int info;
    xgesv(n, nrhs, A, B, ipiv, info);
    if (info<0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Matrix argument(s) had illegal value(s).");
    }
    else if (info>0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Matrix is singular to working precision.");
    }
}

//========================================================================
// [tgetrf]
/// LU factorization using partial pivoting.
///
/// TGETRF computes an LU factorization of a general M-by-N matrix using
/// partial pivoting with row interchanges. The factorization has the form:
///    P * A = L * U
/// where P is a row-pivot matrix, L is lower triangular with unit diagonal
/// elements (lower trapezoidal if m > n), and U is upper triangular (upper
/// trapezoidal if m < n).
///
/// This is the right-looking Level 3 BLAS version of the algorithm.
void xgetrf(int& m, int& n, std::vector<float>& A, std::vector<int>& ipiv, int& info)
{
    sgetrf_(&m, &n, &A[0], &m, &ipiv[0], &info);
}
void xgetrf(int& m, int& n, std::vector<double>& A, std::vector<int>& ipiv, int& info)
{
    dgetrf_(&m, &n, &A[0], &m, &ipiv[0], &info);
}
void xgetrf(int& m, int& n, std::vector<clpk>& A, std::vector<int>& ipiv, int& info)
{
    cgetrf_(&m, &n, &A[0], &m, &ipiv[0], &info);
}
void xgetrf(int& m, int& n, std::vector<zlpk>& A, std::vector<int>& ipiv, int& info)
{
    zgetrf_(&m, &n, &A[0], &m, &ipiv[0], &info);
}
template<typename T, typename S>
void tgetrf(int& m, int& n, std::vector<T>& A, std::vector<T>& P, matrix<S>& U)
{
    int l = std::min(m,n);
    if (A.size()!=m*n || P.size()!=m*m)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    std::vector<int> ipiv(l), idx(m);
    int info, tmp;
    xgetrf(m, n, A, ipiv, info);
    if (info<0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Matrix argument(s) had illegal value(s).");
    }
    else if (info>0)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"The factorization has been completed, but the factor U is exactly singular.");
    }
    U = zeros<S>(l,n);
    for (int j=0; j<n; ++j)
    {
        for (int i=0; i<=std::min(j,l-1); ++i)
        {
            U(i,j) = A[j*m+i];
            (i==j)? (A[j*m+i]=1) : (A[j*m+i]=0);
        }
    }
    A.resize(m*l);
    for (int i=0; i<m; ++i) {idx[i]=i;}
    for (int i=0; i<l; ++i)
    {
        tmp = idx[i];
        idx[i] = idx[ipiv[i]-1];
        idx[ipiv[i]-1] = tmp;
    }
    for (int i=0; i<m; ++i)
    {
        P[idx[i]*m+i] = 1;
    }
}


//==========================================================================//
//                       MATLAB LIKE FUNCTION                               //
//==========================================================================//
//==========================================================================
// [aca]
/// ACA compression with total or partial pivoting.
///
/// [A,Bt] = aca(M,TOL,RMAX) perform Adaptive Cross Approximation on
/// rectangular m-by-n matrix M, representing M as a low-rank product :
///    M = A * B^t,
/// where A and Bt are low-rank matrices respectively m-by-rank(M)
/// and rank(M)-by-n.
/// Matrix M can be evaluated on the fly by a function FCT(i,j)=M(i,j),
/// associated to row I and column J indices.
/// Default accuracy is TOL=1e-6 and maximum rank is RMAX=1e6.
///
/// Only works for double and std::complex<double> type.
///
/// \code{.cpp}
///    matrix<> M = mtimes(rand(5,3),rand(3,6));
///    matrix<> A, Bt;
///    std::tie(A,Bt) = aca(M,1e-3);
///    disp(M-mtimes(A,Bt));
/// \endcode
///
// \see svd, rank.
auto aca(matrix<std::size_t> I, matrix<std::size_t> J,
         std::function<matrix<double>(matrix<std::size_t>,matrix<std::size_t>)>const& fct,
         double tol=1e-6, std::size_t rmax=1e6, bool acaplus=true);
auto aca(matrix<std::size_t> I, matrix<std::size_t> J,
         std::function<matrix<std::complex<double>>(matrix<std::size_t>,matrix<std::size_t>)>const& fct,
         double tol=1e-6, std::size_t rmax=1e6, bool acaplus=true);
template<typename T>
auto aca(matrix<T>const& M, double tol=1e-6, std::size_t rmax=1e6, bool acaplus=true)
{
    matrix<std::size_t> I=range(0,size(M,1)), J=range(0,size(M,2));
    std::function<matrix<T>(matrix<std::size_t>,matrix<std::size_t>)> fct;
    fct = [&M](matrix<std::size_t> I, matrix<std::size_t> J)
    {
        return eval(M(I,J));
    };
    return aca(I,J,fct,tol,rmax,acaplus);
}
template<typename T>
auto aca(matrix<T>const& A, matrix<T>const& B, double tol=1e-6, std::size_t rmax=1e6, bool acaplus=true)
{
    matrix<std::size_t> I=range(0,size(A,1)), J=range(0,size(B,2));
    std::function<matrix<T>(matrix<std::size_t>,matrix<std::size_t>)> fct;
    fct = [&A,&B](matrix<std::size_t> I, matrix<std::size_t> J)
    {
        matrix<T> C(numel(I),numel(J),0);
        for (std::size_t i=0; i<numel(I); ++i)
        {
            for (std::size_t j=0; j<numel(J); ++j)
            {
                for (std::size_t k=0; k<size(A,2); ++k)
                {
                    C(i,j) += A(I(i),k)*B(k,J(j));
                }
            }
        }
        return C;
    };
    return aca(I,J,fct,tol,rmax,acaplus);
}

//==========================================================================
// [eig]
/// Eigenvalues and eigenvectors.
///
/// E = eig(A) produces a row vector E containing the eigenvalues of
/// a square matrix A.
///
/// [E,U] = eig(A,"left") produces a row vector E containing the eigenvalues of
/// a square matrix A, associated to the left eigenvectors U as:
/// U^t*A = diag(E)*U^t,
/// where U^t is the conjugate transpose.
///
/// [E,V] = eig(A,"right") produces a row vector E containing the eigenvalues of
/// a square matrix A, associated to the right eigenvector V as:
/// A*V = V*diag(E).
///
/// [E,V] = eig(A,"none") is equivalent to E = eig(A).
///
/// E = eig(A,B) produces a row vector E containing the generalized eigenvalues
/// of square matrices matrices A and B
///
/// [E,U] = eig(A,B,"left") produces a row vector and a matrix U containing the
/// generalized eigenvalues and the corresponding left eigenvectors such that
/// U^t*A = diag(E)*U^t*B
/// where U^t is the conjugate transpose.
///
/// [E,V] = eig(A,B,"right") produces a row vector and a matrix V containing the
/// generalized eigenvalues and the corresponding right eigenvectors such that
/// A*V = B*V*diag(E)
///
/// [E,V] = eig(A,B,"none") is equivalent to E = eig(A,B)
///
/// \code{.cpp}
///    matrix<> A = 1-eye(4);
///    matrix<std::complex<double>> E, V;
///    std::tie(E,V) = eig(A,"right");
///    disp(E);
///    disp(mtimes(A,V)-mtimes(V,diag(E)));
/// \endcode
///
// \see qr, svd.
auto eig(matrix<float>const& A, std::string typ)
{
    int n = (int)size(A,1);
    std::vector<clpk> Alpk = mat2lpk<clpk>(A,numel(A));
    std::vector<clpk> Elpk(n), Vlpk;
    tgeev(typ, n, Alpk, Elpk, Vlpk);
    matrix<std::complex<float>> E=Elpk, V;
    if (Vlpk.size()>0) {V = lpk2mat<std::complex<float>>(Vlpk,n,n);}
    return std::make_tuple(E,V);
}
auto eig(matrix<double>const& A, std::string typ)
{
    int n = (int)size(A,1);
    std::vector<zlpk> Alpk = mat2lpk<zlpk>(A,numel(A));
    std::vector<zlpk> Elpk(n), Vlpk;
    tgeev(typ, n, Alpk, Elpk, Vlpk);
    matrix<std::complex<double>> E=Elpk, V;
    if (Vlpk.size()>0) {V = lpk2mat<std::complex<double>>(Vlpk,n,n);}
    return std::make_tuple(E,V);
}
auto eig(matrix<std::complex<float>>const& A, std::string typ)
{
    int n = (int)size(A,1);
    std::vector<clpk> Alpk = mat2lpk<clpk>(A,numel(A));
    std::vector<clpk> Elpk(n), Vlpk;
    tgeev(typ, n, Alpk, Elpk, Vlpk);
    matrix<std::complex<float>> E=Elpk, V;
    if (Vlpk.size()>0) {V = lpk2mat<std::complex<float>>(Vlpk,n,n);}
    return std::make_tuple(E,V);
}
auto eig(matrix<std::complex<double>>const& A, std::string typ)
{
    int n = (int)size(A,1);
    std::vector<zlpk> Alpk = mat2lpk<zlpk>(A,numel(A));
    std::vector<zlpk> Elpk(n), Vlpk;
    tgeev(typ, n, Alpk, Elpk, Vlpk);
    matrix<std::complex<double>> E=Elpk, V;
    if (Vlpk.size()>0) {V = lpk2mat<std::complex<double>>(Vlpk,n,n);}
    return std::make_tuple(E,V);
}
auto eig(matrix<float> const &A, matrix<float> const &B, std::string typ)
{
    int n = (int)size(A,1);
    std::vector<clpk> Alpk = mat2lpk<clpk>(A,numel(A));
    std::vector<clpk> Blpk = mat2lpk<clpk>(B,numel(B));
    std::vector<clpk> Elpk(n), Vlpk;
    tggev(typ, n, Alpk, Blpk, Elpk, Vlpk);
    matrix<std::complex<float>> E=Elpk, V;
    if (Vlpk.size()>0) {V = lpk2mat<std::complex<float>>(Vlpk,n,n);}
    return std::make_tuple(E,V);
}
auto eig(matrix<double> const &A, matrix<double> const &B, std::string typ)
{
    int n = (int)size(A,1);
    std::vector<zlpk> Alpk = mat2lpk<zlpk>(A,numel(A));
    std::vector<zlpk> Blpk = mat2lpk<zlpk>(B,numel(B));
    std::vector<zlpk> Elpk(n), Vlpk;
    tggev(typ, n, Alpk, Blpk, Elpk, Vlpk);
    matrix<std::complex<double>> E=Elpk, V;
    if (Vlpk.size()>0) {V = lpk2mat<std::complex<double>>(Vlpk,n,n);}
    return std::make_tuple(E,V);
}
auto eig(matrix<std::complex<float>> const &A, matrix<std::complex<float>> const &B, std::string typ)
{
    int n = (int)size(A,1);
    std::vector<clpk> Alpk = mat2lpk<clpk>(A,numel(A));
    std::vector<clpk> Blpk = mat2lpk<clpk>(B,numel(B));
    std::vector<clpk> Elpk(n), Vlpk;
    tggev(typ, n, Alpk, Blpk, Elpk, Vlpk);
    matrix<std::complex<float>> E=Elpk, V;
    if (Vlpk.size()>0) {V = lpk2mat<std::complex<float>>(Vlpk,n,n);}
    return std::make_tuple(E,V);
}
auto eig(matrix<std::complex<double>> const &A, matrix<std::complex<double>> const &B, std::string typ)
{
    int n = (int)size(A,1);
    std::vector<zlpk> Alpk = mat2lpk<zlpk>(A,numel(A));
    std::vector<zlpk> Blpk = mat2lpk<zlpk>(B,numel(B));
    std::vector<zlpk> Elpk(n), Vlpk;
    tggev(typ, n, Alpk, Blpk, Elpk, Vlpk);
    matrix<std::complex<double>> E=Elpk, V;
    if (Vlpk.size()>0) {V = lpk2mat<std::complex<double>>(Vlpk,n,n);}
    return std::make_tuple(E,V);
}
template<typename T>
matrix<T> eig(matrix<T> const & A)
{
    matrix<T> E, V;
    std::tie(E,V) = eig(A,"none");
    return E;
}
template<typename T>
matrix<T> eig(matrix<T> const &A, matrix<T> const &B)
{
    matrix<T> E, V;
    std::tie(E,V) = eig(A,B,"none");
    return E;
}

//==========================================================================
// [inv]
/// Matrix inverse.
///
/// inv(X) is the inverse of the square matrix X.
/// A warning message is printed if X is badly scaled or nearly singular.
///
/// \code{.cpp}
///    matrix<> A = rand(3,3);
///    matrix<> X = inv(A);
///    disp(mtimes(X,A));
/// \endcode
///
// \see linsolve, pinv, gmres.
template<typename T>
matrix<T> inv(matrix<T>const& A)
{
    if (size(A,1)!=size(A,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix must be square.");
    }
    matrix<T> B = eye<T>(size(A));
    return linsolve(A,B);
}

//==========================================================================
// [linsolve]
/// linsolve Solve linear system A*X=B.
///
/// X = linsolve(A,B) solves the linear system A*X=B using LU factorization
/// with partial pivoting when A is square, and QR factorization with
/// column pivoting otherwise.
/// linsolve warns if A is ill conditioned (for square matrices) or rank
/// deficient (for rectangular matrices).
///
/// \code{.cpp}
///    matrix<> A = rand(3,3);
///    matrix<> B = eye(3,3);
///    matrix<> X = linsolve(A,B);
///    disp(mtimes(X,A));
/// \endcode
///
// \see inv, pinv, gmres.
matrix<float> linsolve(matrix<float>const& A, matrix<float>const& B)
{
    int m=(int)size(A,1), n=(int)size(A,2), nrhs=(int)size(B,2);
    std::vector<float> Alpk = mat2lpk<float>(A,numel(A));
    std::vector<float> Blpk = mat2lpk<float>(B,std::max(m,n)*nrhs);
    if (m==n) {tgesv(m,nrhs,Alpk,Blpk);}
    else {tgels(m,n,nrhs,Alpk,Blpk);}
    return lpk2mat<float>(Blpk,size(A,2),size(B,2));
}
matrix<double> linsolve(matrix<double>const& A, matrix<double>const& B)
{
    int m=(int)size(A,1), n=(int)size(A,2), nrhs=(int)size(B,2);
    std::vector<double> Alpk = mat2lpk<double>(A,numel(A));
    std::vector<double> Blpk = mat2lpk<double>(B,std::max(m,n)*nrhs);
    if (m==n) {tgesv(m,nrhs,Alpk,Blpk);}
    else {tgels(m,n,nrhs,Alpk,Blpk);}
    return lpk2mat<double>(Blpk,size(A,2),size(B,2));
}
matrix<std::complex<float>> linsolve(matrix<std::complex<float>>const& A, matrix<std::complex<float>>const& B)
{
    int m=(int)size(A,1), n=(int)size(A,2), nrhs=(int)size(B,2);
    std::vector<clpk> Alpk = mat2lpk<clpk>(A,numel(A));
    std::vector<clpk> Blpk = mat2lpk<clpk>(B,std::max(m,n)*nrhs);
    if (m==n) {tgesv(m,nrhs,Alpk,Blpk);}
    else {tgels(m,n,nrhs,Alpk,Blpk);}
    return lpk2mat<std::complex<float>>(Blpk,size(A,2),size(B,2));
}
matrix<std::complex<double>> linsolve(matrix<std::complex<double>>const& A, matrix<std::complex<double>>const& B)
{
    int m=(int)size(A,1), n=(int)size(A,2), nrhs=(int)size(B,2);
    std::vector<zlpk> Alpk = mat2lpk<zlpk>(A,numel(A));
    std::vector<zlpk> Blpk = mat2lpk<zlpk>(B,std::max(m,n)*nrhs);
    if (m==n) {tgesv(m,nrhs,Alpk,Blpk);}
    else {tgels(m,n,nrhs,Alpk,Blpk);}
    return lpk2mat<std::complex<double>>(Blpk,size(A,2),size(B,2));
}

//==========================================================================
// [lu]
/// LU factorization.
///
/// [L,U] = lu(A) returns an upper triangular matrix in U and a permuted
/// lower triangular matrix in L, such that A = L*U. The input matrix A can
/// be rectangular.
///
/// \code{.cpp}
///    matrix<> A = rand(3,4);
///    matrix<> L, U;
///    std::tie(L,U) = lu(A);
///    disp(mtimes(L,U)-A);
/// \endcode
///
// \see qr, linsolve, inv.
auto lu(matrix<float>const& A)
{
    int m=(int)size(A,1), n=(int)size(A,2), l = std::min(m,n);
    std::vector<float> Alpk = mat2lpk<float>(A,numel(A)), P(m*m,0);
    matrix<float> L, U;
    tgetrf(m,n,Alpk,P,U);
    tgesv(m,l,P,Alpk);
    L = lpk2mat<float>(Alpk,m,l);
    return std::make_tuple(L,U);
}
auto lu(matrix<double>const& A)
{
    int m=(int)size(A,1), n=(int)size(A,2), l = std::min(m,n);
    std::vector<double> Alpk = mat2lpk<double>(A,numel(A)), P(m*m,0);
    matrix<double> L, U;
    tgetrf(m,n,Alpk,P,U);
    tgesv(m,l,P,Alpk);
    L = lpk2mat<double>(Alpk,m,l);
    return std::make_tuple(L,U);
}
auto lu(matrix<std::complex<float>>const& A)
{
    int m=(int)size(A,1), n=(int)size(A,2), l = std::min(m,n);
    std::vector<clpk> Alpk = mat2lpk<clpk>(A,numel(A)), P(m*m,0);
    matrix<std::complex<float>> L, U;
    tgetrf(m,n,Alpk,P,U);
    tgesv(m,l,P,Alpk);
    L = lpk2mat<std::complex<float>>(Alpk,m,l);
    return std::make_tuple(L,U);
}
auto lu(matrix<std::complex<double>>const& A)
{
    int m=(int)size(A,1), n=(int)size(A,2), l = std::min(m,n);
    std::vector<zlpk> Alpk = mat2lpk<zlpk>(A,numel(A)), P(m*m,0);
    matrix<std::complex<double>> L, U;
    tgetrf(m,n,Alpk,P,U);
    tgesv(m,l,P,Alpk);
    L = lpk2mat<std::complex<double>>(Alpk,m,l);
    return std::make_tuple(L,U);
}

//==========================================================================
// [pinv]
/// Pseudoinverse.
///
/// X = pinv(A) produces a matrix X of the same dimensions
/// as A' so that A*X*A = A, X*A*X = X and A*X and X*A
/// are Hermitian.
/// The computation is based on QR factorization of A.
///
/// \code{.cpp}
///    matrix<> A = rand(4,3);
///    matrix<> X = pinv(A);
///    disp(mtimes(X,A));
/// \endcode
///
// \see inv, linsolve.
template<typename T>
matrix<T> pinv(matrix<T>const& A)
{
    matrix<T> B = eye<T>(size(A,1));
    return linsolve(A,B);
}

//==========================================================================
// [qr]
/// Orthogonal-triangular decomposition.
///
/// [Q,R] = qr(A), where A is m-by-n, produces an l-by-n upper triangular
/// matrix R and an m-by-l unitary matrix Q so that A = Q*R.
///
/// \code{.cpp}
///    matrix<> A = rand(4,3);
///    matrix<> Q, R;
///    std::tie(Q,R) = qr(A);
///    disp(mtimes(Q,R)-A);
/// \endcode
///
// \see eig, svd, lu.
auto qr(matrix<float>const& A)
{
    int m = (int)size(A,1);
    int n = (int)size(A,2);
    int l = std::min<int>(m,n);
    std::vector<float> Alpk = mat2lpk<float>(A,numel(A));
    std::vector<float> Rlpk(l*n,0.);
    tgeqrf(m,n,Alpk,Rlpk);
    matrix<float> Q = lpk2mat<float>(Alpk,m,l);
    matrix<float> R = lpk2mat<float>(Rlpk,l,n);
    return std::make_tuple(Q,R);
}
auto qr(matrix<double>const& A)
{
    int m = (int)size(A,1);
    int n = (int)size(A,2);
    int l = std::min<int>(m,n);
    std::vector<double> Alpk = mat2lpk<double>(A,numel(A));
    std::vector<double> Rlpk(l*n,0.);
    tgeqrf(m,n,Alpk,Rlpk);
    matrix<double> Q = lpk2mat<double>(Alpk,m,l);
    matrix<double> R = lpk2mat<double>(Rlpk,l,n);
    return std::make_tuple(Q,R);
}
auto qr(matrix<std::complex<float>>const& A)
{
    int m = (int)size(A,1);
    int n = (int)size(A,2);
    int l = std::min<int>(m,n);
    std::vector<clpk> Alpk = mat2lpk<clpk>(A,numel(A));
    std::vector<clpk> Rlpk(l*n,0.);
    tgeqrf(m,n,Alpk,Rlpk);
    matrix<std::complex<float>> Q = lpk2mat<std::complex<float>>(Alpk,m,l);
    matrix<std::complex<float>> R = lpk2mat<std::complex<float>>(Rlpk,l,n);
    return std::make_tuple(Q,R);
}
auto qr(matrix<std::complex<double>>const& A)
{
    int m = (int)size(A,1);
    int n = (int)size(A,2);
    int l = std::min<int>(m,n);
    std::vector<zlpk> Alpk = mat2lpk<zlpk>(A,numel(A));
    std::vector<zlpk> Rlpk(l*n,0.);
    tgeqrf(m,n,Alpk,Rlpk);
    matrix<std::complex<double>> Q = lpk2mat<std::complex<double>>(Alpk,m,l);
    matrix<std::complex<double>> R = lpk2mat<std::complex<double>>(Rlpk,l,n);
    return std::make_tuple(Q,R);
}

//========================================================================
// [qrsvd]
/// Low-rank recompression.
///
/// [U,V] = qrsvd(A,B,TOL) provides a low-rank recompression using
/// QR factorization and SVD.
///
/// \code{.cpp}
///    matrix<> A = rand(5,3), B = rand(3,6), U, V;
///    std::tie(U,V) = qrsvd(A,B,1e-3);
///    disp(mtimes(A,B)-mtimes(U,V));
/// \endcode
///
// \see svd, qr.
template<typename T>
auto qrsvd(matrix<T>const& A, matrix<T>const& B, float tol)
{
    matrix<T> U, V, Qa, Ra, Qb, Rb, Uab, Vab;
    matrix<decltype(std::abs(A(0)))> Sab;
    std::tie(Qa,Ra) = qr(A);
    std::tie(Qb,Rb) = qr(transpose(B));
    std::tie(Sab,Uab,Vab) = svd(mtimes(Ra,transpose(Rb)),"vect");
    std::size_t n = sum<std::size_t>(Sab>=tol*Sab(0));
    Sab = diag(eval(Sab(range(0,n))));
    Uab = eval(Uab(row(Uab),range(0,n)));
    Vab = eval(Vab(range(0,n),col(Vab)));
    U   = mtimes(Qa,mtimes(Uab,Sab));
    V   = mtimes(Vab,transpose(Qb));
    return std::make_tuple(U,V);
}

//========================================================================
// [rank]
/// Matrix rank.
///
/// r = rank(A) provides an estimate of the number of linearly independent
/// rows or columns of a matrix A.
///
/// r = rank(A,TOL) is the number of singular values of A that are larger
/// than TOL. By default, TOL = max(size(A)) * eps(norm(A)).
///
/// \code{.cpp}
///    matrix<> A = mtimes(rand(10,3),rand(3,20));
///    std::size_t r = rank(A);
///    disp(r);
/// \endcode
///
// \see svd, aca.
template<typename T>
std::size_t rank(matrix<T>const& A, float tol=0)
{
    using S = decltype(std::abs(A(0)));
    matrix<S> s = svd(A);
    if (tol==0) {tol = std::max(size(A,1),size(A,2)) * M_EPS(S);}
    else {tol *= s(0);}
    return sum<std::size_t>(s>tol);
}

//==========================================================================
// [svd]
/// Singular value decomposition.
///
/// S = svd(A) returns a row vector S containing the singular values of A,
/// with nonnegative elements sorted in decreasing order.
///
/// [S,U,Vt] = svd(A,"vect") returns singulars values and unitary matrices
/// U and Vt (not V) so that :
///    A = U*diag(S)*V^t
///
/// [U,S,Vt] = svd(A,"none") is equivalent to S = svd(A).
///
/// \code{.cpp}
///    matrix<> A = rand(3,3);
///    matrix<> U, S, Vt;
///    std::tie(S,U,Vt) = svd(A,"vect");
///    disp(mtimes(mtimes(U,diag(S)),Vt)-A);
/// \endcode
///
// \see rank, eig, qr, qrsvd, aca.
auto svd(matrix<float>const& A, std::string typ)
{
    int m = (int)size(A,1);
    int n = (int)size(A,2);
    int l = std::min<int>(m,n);
    std::vector<float> Alpk = mat2lpk<float>(A,numel(A));
    std::vector<float> Slpk(l), Ulpk, Vlpk;
    tgesdd(typ,m,n,Alpk,Slpk,Ulpk,Vlpk);
    matrix<float> S=Slpk, U, Vt;
    if (Ulpk.size()>0) {U = lpk2mat<float>(Ulpk,m,l);}
    if (Vlpk.size()>0) {Vt = lpk2mat<float>(Vlpk,l,n);}
    return std::make_tuple(S,U,Vt);
}
auto svd(matrix<double>const& A, std::string typ)
{
    int m = (int)size(A,1);
    int n = (int)size(A,2);
    int l = std::min<int>(m,n);
    std::vector<double> Alpk = mat2lpk<double>(A,numel(A));
    std::vector<double> Slpk(l), Ulpk, Vlpk;
    tgesdd(typ,m,n,Alpk,Slpk,Ulpk,Vlpk);
    matrix<double> S=Slpk, U, Vt;
    if (Ulpk.size()>0) {U = lpk2mat<double>(Ulpk,m,l);}
    if (Vlpk.size()>0) {Vt = lpk2mat<double>(Vlpk,l,n);}
    return std::make_tuple(S,U,Vt);
}
auto svd(matrix<std::complex<float>>const& A, std::string typ)
{
    int m = (int)size(A,1);
    int n = (int)size(A,2);
    int l = std::min<int>(m,n);
    std::vector<clpk> Alpk = mat2lpk<clpk>(A,numel(A));
    std::vector<float> Slpk(l);
    std::vector<clpk> Ulpk, Vlpk;
    tgesdd(typ,m,n,Alpk,Slpk,Ulpk,Vlpk);
    matrix<float> S = Slpk;
    matrix<std::complex<float>> U, Vt;
    if (Ulpk.size()>0) {U = lpk2mat<std::complex<float>>(Ulpk,m,l);}
    if (Vlpk.size()>0) {Vt = lpk2mat<std::complex<float>>(Vlpk,l,n);}
    return std::make_tuple(S,U,Vt);
}
auto svd(matrix<std::complex<double>>const& A, std::string typ)
{
    int m = (int)size(A,1);
    int n = (int)size(A,2);
    int l = std::min<int>(m,n);
    std::vector<zlpk> Alpk = mat2lpk<zlpk>(A,numel(A));
    std::vector<double> Slpk(l);
    std::vector<zlpk> Ulpk, Vlpk;
    tgesdd(typ,m,n,Alpk,Slpk,Ulpk,Vlpk);
    matrix<double> S = Slpk;
    matrix<std::complex<double>> U, Vt;
    if (Ulpk.size()>0) {U = lpk2mat<std::complex<double>>(Ulpk,m,l);}
    if (Vlpk.size()>0) {Vt = lpk2mat<std::complex<double>>(Vlpk,l,n);}
    return std::make_tuple(S,U,Vt);
}
template<typename T>
auto svd(matrix<T>const& A)
{
    using S = decltype(std::abs(A(0)));
    matrix<S> s;
    matrix<T> u, vt;
    std::tie(s,u,vt) = svd(A,"none");
    return s;
}


//==========================================================================//
//                         DETAILLED ROUTINES                               //
//==========================================================================//
//========================================================================
auto aca(matrix<std::size_t> I, matrix<std::size_t> J,
         std::function<matrix<double>(matrix<std::size_t>,matrix<std::size_t>)>const& fct,
         double tol, std::size_t rmax, bool acaplus)
{
    // Declarations
    auto row = [&fct,&I,&J](std::size_t i) {return fct(I(i),J);};
    auto col = [&fct,&I,&J](std::size_t j) {return fct(I,J(j));};
    std::size_t Nr=numel(I), Nc=numel(J), r, c, n, k;
    std::vector<std::size_t> Ir(Nr), Ic(Nc);
    double alpha, beta, aa, ab, bb, res, err, tmp;
    std::vector<double> a(Nr), b(Nc), anp1(Nr), bnp1(Nc), u(10), v(10);
    bool flag = true, maxsearch = true;
    matrix<double> A, Bt;
    
    // Initialize indices for row and columns pivots
    for(std::size_t l=0; l<Nr; ++l) {Ir[l] = l;}
    for(std::size_t l=0; l<Nc; ++l) {Ic[l] = l;}
    
    // First row (first matrix row)
    r = 0;
    b = row(r).val(); //for(std::size_t l=0; l<Nc; ++l) {b[l] = M(r,l);}
    Ir.erase(Ir.begin()+r);
    
    // Row pivot
    c    = 0;
    beta = 0;
    for(std::size_t l=0; l<Nc; ++l)
    {
        if(std::abs(b[l])>std::abs(beta)) {beta=b[l]; c=l;}
    }
    
    // First column
    a = (col(c)/beta).val(); //for(std::size_t l=0; l<Nr; ++l) {a[l] = M(l,c)/beta;}
    Ic.erase(Ic.begin()+c);
    
    // Frobenius residue to compute error
    aa=0; bb=0;
    for(std::size_t l=0; l<Nr; ++l) {aa += a[l]*a[l];}
    for(std::size_t l=0; l<Nc; ++l) {bb += b[l]*b[l];}
    res = aa*bb;
    err = std::sqrt(aa)*std::sqrt(bb)/std::sqrt(res);
    
    // Initialization for recurisive loop
    anp1 = a;
    bnp1 = b;
    n    = 1;
    
    // Iterative construction using recursive frobenius error
    while (err > tol)
    {
        // Compression failed
        if (n*(Nr+Nc) > Nr*Nc || n>=rmax)
        {
            flag = false;
            break;
        }
        
        // Resize
        u.resize(n);
        v.resize(n);
        
        // Column pivot with previous last column
        k = 0;
        if (maxsearch) // Look for max
        {
            alpha = 0;
            for(std::size_t l=0; l<Ir.size(); ++l)
            {
                tmp = anp1[Ir[l]];
                if(std::abs(tmp)>std::abs(alpha)) {alpha=tmp; k=l;}
            }
        }
        else // Look for min
        {
            alpha = 1e308;
            for(std::size_t l=0; l<Ir.size(); ++l)
            {
                tmp = anp1[Ir[l]];
                if(std::abs(tmp)<std::abs(alpha)) {alpha=tmp; k=l;}
            }
            maxsearch = true;
        }
        r = Ir[k];
        Ir.erase(Ir.begin()+k);
        
        // New row
        for(std::size_t l=0; l<n; ++l)  {u[l] = a[l*Nr+r];}
        bnp1 = row(r).val();//for(std::size_t l=0; l<Nc; ++l) {bnp1[l] = M(r,l);}
        cblas_dgemv(CblasRowMajor, CblasTrans, (int) n, (int) Nc,
                    -1., &b[0], (int) Nc, &u[0], 1, 1., &bnp1[0], 1);
        
        // Row pivot
        k    = 0;
        beta = 0;
        for(std::size_t l=0; l<Ic.size(); ++l)
        {
            tmp = bnp1[Ic[l]];
            if(std::abs(tmp)>std::abs(beta)) {beta=tmp; k=l;}
        }
        c = Ic[k];
        Ic.erase(Ic.begin()+k);
        if (beta==0) {break;}
        
        // New column
        for(std::size_t l=0; l<n; ++l)  {v[l]    = b[l*Nc+c];}
        anp1 = col(c).val();//for(std::size_t l=0; l<Nr; ++l) {anp1[l] = M(l,c);}
        cblas_dgemv(CblasRowMajor, CblasTrans, (int) n, (int) Nr,
                    -1., &a[0], (int) Nr, &v[0], 1, 1., &anp1[0], 1);
        for(std::size_t l=0; l<Nr; ++l) {anp1[l] /= beta;}
        
        // Dot product
        aa=0; bb=0;
        for(std::size_t l=0; l<Nr; ++l) {aa += anp1[l]*anp1[l];}
        for(std::size_t l=0; l<Nc; ++l) {bb += bnp1[l]*bnp1[l];}
        
        // Matrix product
        ab=0;
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (int) n, (int) Nr,
                    1., &a[0], (int) Nr, &anp1[0], 1, 0., &u[0], 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (int) n, (int) Nc,
                    1., &b[0], (int) Nc, &bnp1[0], 1, 0., &v[0], 1);
        for(std::size_t k=0; k<n; ++k) {ab += u[k]*v[k];}
        
        // Relative error (block Frobenius)
        res += 2*ab + aa*bb;
        err = std::sqrt(aa)*std::sqrt(bb)/std::sqrt(res);
        
        // Update matrix
        a.resize(a.size()+Nr);
        b.resize(b.size()+Nc);
        for(std::size_t l=0; l<Nr; ++l) {a[n*Nr+l] = anp1[l];}
        for(std::size_t l=0; l<Nc; ++l) {b[n*Nc+l] = bnp1[l];}
        
        // Incrementation
        ++n;
        
        // ACA plus activation
        if (err<=tol && acaplus)
        {
            err       = 1;
            maxsearch = false;
            acaplus   = false;
        }
    }
    
    // Matrix format (A,Bt)
    A  = a;
    A.reshape(n,Nr);
    A  = transpose(A);
    Bt = b;
    Bt.reshape(n,Nc);
    return std::make_tuple(A,Bt,flag);
}

//========================================================================
auto aca(matrix<std::size_t> I, matrix<std::size_t> J,
         std::function<matrix<std::complex<double>>(matrix<std::size_t>,matrix<std::size_t>)>const& fct,
         double tol, std::size_t rmax, bool acaplus)
{
    // Declarations
    auto row = [&fct,&I,&J](std::size_t i) {return fct(I(i),J);};
    auto col = [&fct,&I,&J](std::size_t j) {return fct(I,J(j));};
    std::size_t Nr=numel(I), Nc=numel(J), r, c, n, k;
    std::vector<std::size_t> Ir(Nr), Ic(Nc);
    std::complex<double> alpha, beta, aa, ab, bb, res, err, tmp;
    std::vector<std::complex<double>> a(Nr), b(Nc), anp1(Nr), bnp1(Nc), u(10), v(10);
    bool flag = true, maxsearch = true;
    matrix<std::complex<double>> A, Bt;
    
    // cblas_zgemv constantes
    std::complex<double> cb_un(1,0), cb_mun(-1,0), cb_z(0,0);
    
    // Initialize indices for row and columns pivots
    for(std::size_t l=0; l<Nr; ++l) {Ir[l] = l;}
    for(std::size_t l=0; l<Nc; ++l) {Ic[l] = l;}
    
    // First row (first matrix row)
    r = 0;
    b = row(r).val(); //for(std::size_t l=0; l<Nc; ++l) {b[l] = M(r,l);}
    Ir.erase(Ir.begin()+r);
    
    // Row pivot
    c    = 0;
    beta = 0;
    for(std::size_t l=0; l<Nc; ++l)
    {
        if(std::abs(b[l])>std::abs(beta)) {beta=b[l]; c=l;}
    }
    
    // First column
    a = (col(c)/beta).val(); //for(std::size_t l=0; l<Nr; ++l) {a[l] = M(l,c)/beta;}
    Ic.erase(Ic.begin()+c);
    
    // Frobenius residue to compute error
    aa=0; bb=0;
    for(std::size_t l=0; l<Nr; ++l) {aa += conj(a[l])*a[l];}
    for(std::size_t l=0; l<Nc; ++l) {bb += conj(b[l])*b[l];}
    res = aa*bb;
    err = sqrt(aa)*sqrt(bb)/sqrt(res);
    
    // Initialization for recurisive loop
    anp1 = a;
    bnp1 = b;
    n    = 1;
    
    // Iterative construction using recursive frobenius error
    while (err.real() > tol)
    {
        // Compression failed
        if (n*(Nr+Nc) > Nr*Nc || n>=rmax)
        {
            flag = false;
            break;
        }

        // Resize
        u.resize(n);
        v.resize(n);
        
        // Column pivot with previous last column
        k = 0;
        if (maxsearch) // Look for max
        {
            alpha = 0;
            for(std::size_t l=0; l<Ir.size(); ++l)
            {
                tmp = anp1[Ir[l]];
                if(std::abs(tmp)>std::abs(alpha)) {alpha=tmp; k=l;}
            }
        }
        else // Look for min
        {
            alpha = 1e308;
            for(std::size_t l=0; l<Ir.size(); ++l)
            {
                tmp = anp1[Ir[l]];
                if(std::abs(tmp)<std::abs(alpha)) {alpha=tmp; k=l;}
            }
            maxsearch = true;
        }
        r = Ir[k];
        Ir.erase(Ir.begin()+k);
        
        // New row
        for(std::size_t l=0; l<n; ++l)  {u[l]    = a[l*Nr+r];}
        bnp1 = row(r).val(); //for(std::size_t l=0; l<Nc; ++l) {bnp1[l] = M(r,l);}
        cblas_zgemv(CblasRowMajor, CblasTrans, (int) n, (int) Nc,
                    &cb_mun, &b[0], (int) Nc, &u[0], 1, &cb_un, &bnp1[0], 1);
        
        // Row pivot
        k    = 0;
        beta = 0;
        for(std::size_t l=0; l<Ic.size(); ++l)
        {
            tmp = bnp1[Ic[l]];
            if(std::abs(tmp)>std::abs(beta)) {beta=tmp; k=l;}
        }
        c = Ic[k];
        Ic.erase(Ic.begin()+k);
        if (std::abs(beta)==0) {break;}
        
        // New column
        for(std::size_t l=0; l<n; ++l)  {v[l]    = b[l*Nc+c];}
        anp1 = col(c).val(); //for(std::size_t l=0; l<Nr; ++l) {anp1[l] = M(l,c);}
        cblas_zgemv(CblasRowMajor, CblasTrans, (int) n, (int) Nr,
                    &cb_mun, &a[0], (int) Nr, &v[0], 1, &cb_un, &anp1[0], 1);
        for(std::size_t l=0; l<Nr; ++l) {anp1[l] /= beta;}
        
        // Dot product
        aa=0; bb=0;
        for(std::size_t l=0; l<Nr; ++l) {aa += conj(anp1[l])*anp1[l];}
        for(std::size_t l=0; l<Nc; ++l) {bb += conj(bnp1[l])*bnp1[l];}
        
        // Matrix product
        ab=0;
        cblas_zgemv(CblasRowMajor, CblasNoTrans, (int) n, (int) Nr,
                    &cb_mun, &a[0], (int) Nr, &anp1[0], 1, &cb_z, &u[0], 1);
        cblas_zgemv(CblasRowMajor, CblasNoTrans, (int) n, (int) Nc,
                    &cb_mun, &b[0], (int) Nc, &bnp1[0], 1, &cb_z, &v[0], 1);
        for(std::size_t k=0; k<n; ++k) {ab += u[k]*v[k];}
        
        // Relative error (block Frobenius)
        res += 2*ab.real() + aa*bb;
        err = sqrt(aa)*sqrt(bb)/sqrt(res);
        
        // Update matrix
        a.resize(a.size()+Nr);
        b.resize(b.size()+Nc);
        for(std::size_t l=0; l<Nr; ++l) {a[n*Nr+l] = anp1[l];}
        for(std::size_t l=0; l<Nc; ++l) {b[n*Nc+l] = bnp1[l];}
        
        // Incrementation
        ++n;
        
        // ACA plus activation
        if (err.real()<=tol && acaplus)
        {
            err       = 1;
            maxsearch = false;
            acaplus   = false;
        }
    }
    
    // Matrix format (A,Bt)
    A  = a;
    A.reshape(n,Nr);
    A  = transpose(A);
    Bt = b;
    Bt.reshape(n,Nc);
    return std::make_tuple(A,Bt,flag);
}

}
