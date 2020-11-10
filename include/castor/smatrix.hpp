/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : smatrix.hpp                                   |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Sparse matrix operations & algebra with       |
 |  `---'  |                matlab-like standard functions and procedures |
 +========================================================================+
 */

#pragma once
#define CASTOR_SMATRIX_HPP
#include <castor/matrix.hpp>

namespace castor
{

//==========================================================================//
//                            ENVIRONMENT DATA                              //
//==========================================================================//
template<typename T>
static const T M_ZERO=0;

template<typename T>
class sview;

template<typename T>
class scview;


//==========================================================================//
//                             MASTER CLASS                                 //
//==========================================================================//
// [smatrix]
/// Array with element of type T (double by default), stored in sparse format,
/// e.g. only non-zeros values are keeped in memory.
///
/// Values are stored in std::vector<T>, associated to linear indices in
/// std::vector<std::size_t> of the same dimensions of the values.
///
/// \code{.cpp}
///    smatrix<logical> s = true;
///    smatrix<int>     X = {1, 2, 3, 4};
///    smatrix<double>  A = {{1, 2, 3, 4},
///                         {5, 6, 7, 8}};
///    auto             B = A;
/// \endcode
///
// \see sview, scview.
template<typename T=double>
class smatrix
{
public:
    // CONSTRUCTORS
    smatrix();                                              // DEFAULT
    smatrix(T v);                                           // SCALAR
    smatrix(std::size_t m, std::size_t n);                  // DIMENSIONS
    smatrix(std::size_t m, std::size_t n, std::size_t nnz); // DIMENSION AND NNZ
    template<typename S>
    smatrix(std::size_t m, std::size_t n,                   // LINEAR INDEX & VALUES
            std::vector<std::size_t>const& L, std::vector<S>const& V);
    template<typename S>
    smatrix(std::size_t m, std::size_t n,                   // BI-LINEAR INDEX & VALUES
            std::vector<std::size_t>const& I, std::vector<std::size_t>const& J, std::vector<S>const& V);
    template<typename S>
    smatrix(matrix<S>const& A);                             // MATRIX
    template<typename S>
    smatrix(smatrix<S>const& As);                           // SPARSE MATRIX (CAST)
    
    // OPERATORS
    template<typename S>
    explicit    operator S() const;                             // CONVERSION TO PRIMITIVE
    smatrix<T>& operator+=(smatrix<T>const& A);                 // A += B;
    smatrix<T>& operator-=(smatrix<T>const& A);                 // A -= B;
    smatrix<T>& operator*=(smatrix<T>const& A);                 // A *= B;
    smatrix<T>& operator/=(smatrix<T>const& A);                 // A /= B;
    const T&    operator()(std::size_t l) const;                // A(l)
    T&          operator()(std::size_t l);                          // A(l)
    const scview<T> operator()(matrix<std::size_t>const& L) const;  // A(L)
    sview<T>        operator()(matrix<std::size_t>const& L);        // A(L)
    const T&    operator()(std::size_t i, std::size_t j) const; // A(i,j)
    T&          operator()(std::size_t i, std::size_t j);       // A(i,j)
    const scview<T> operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J) const; // A(I,J)
    sview<T>        operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J); // A(I,J)
    
    // TOOLS
    void                           check();
    void                           clear();
    void                           reshape(std::size_t m, std::size_t n);
    std::vector<std::size_t>const& ind() const;
    std::size_t const&             ind(std::size_t k) const;
    std::size_t&                   ind(std::size_t k);
    std::size_t                    nnz() const;
    std::size_t                    size(int dim=0) const;
    std::vector<T>const&           val() const;
    T const&                       val(std::size_t k) const;
    T&                             val(std::size_t k);
    
private:
    std::size_t               m_row; // Number of rows
    std::size_t               m_col; // Number of columns
    std::vector<std::size_t>  m_ind; // Linear indices
    std::vector<T>            m_val; // Values
};


//==========================================================================//
//                           COMMON HEADERS                                 //
//==========================================================================//
template<typename T>
std::size_t nnz(smatrix<T>const& As);

template<typename T>
std::size_t numel(smatrix<T>const& As);

template<typename T>
std::size_t size(smatrix<T>const& As, int dim);


//==========================================================================//
//                                VIEW CLASS                                //
//==========================================================================//
// [sview]
///
template<typename T>
class sview
{
public:
    sview(smatrix<T>& As, matrix<std::size_t>const& L):
    m_typ{1}, m_idx{L}, m_mat{As} {};
    sview(smatrix<T>& As, matrix<std::size_t>const& I, matrix<std::size_t>const& J):
    m_typ{2}, m_idx{I}, m_jdx{J}, m_mat{As} {};
    template<typename S>
    smatrix<T>& operator=(smatrix<S>const& As) const
    {
        if (m_typ==1) {set(m_mat,m_idx,smatrix<T>(As));}
        else if (m_typ==2) {set(m_mat,m_idx,m_jdx,smatrix<T>(As));}
        return m_mat;
    };
    smatrix<T> eval() const
    {
        smatrix<T> As;
        if (m_typ==1) {As = get(m_mat,m_idx);}
        else if (m_typ==2) {As = get(m_mat,m_idx,m_jdx);}
        return As;
    };
private:
    int                 m_typ;
    matrix<std::size_t> m_idx;
    matrix<std::size_t> m_jdx;
    smatrix<T>&         m_mat;
};
template<typename T>
inline std::ostream& operator<<(std::ostream& flux, sview<T>const& As) {disp(As); return flux;}
template<typename T>
inline smatrix<T> eval(sview<T>const& Av) {return Av.eval();}
template<typename T>
inline void disp(sview<T>const& Av)
{
    error(__FILE__, __LINE__, __FUNCTION__,"Use 'eval()' function to build the right hand side matrix from sparse view. For example: \nAs = eval(Bs(L)) \nAs = eval(Bs(I,J))");
}

template<typename T>
class scview
{
public:
    scview(smatrix<T>const& As, matrix<std::size_t>const& L):
    m_typ{1}, m_idx{L}, m_mat{As} {};
    scview(smatrix<T>const& As, matrix<std::size_t>const& I, matrix<std::size_t>const& J):
    m_typ{2}, m_idx{I}, m_jdx{J}, m_mat{As} {};
    smatrix<T> eval() const
    {
        smatrix<T> As;
        if (m_typ==1) {As = get(m_mat,m_idx);}
        else if (m_typ==2) {As = get(m_mat,m_idx,m_jdx);}
        return As;
    };
private:
    int                 m_typ;
    matrix<std::size_t> m_idx;
    matrix<std::size_t> m_jdx;
    smatrix<T>const&    m_mat;
};
template<typename T>
inline std::ostream& operator<<(std::ostream& flux, scview<T>const& As) {disp(As); return flux;}
template<typename T>
inline smatrix<T> eval(scview<T>const& Av) {return Av.eval();}
template<typename T>
inline void disp(scview<T>const& Av)
{
    error(__FILE__, __LINE__, __FUNCTION__,"Use 'eval()' function to build the right hand side matrix from sparse view. For example: \nAs = eval(Bs(L)) \nAs = eval(Bs(I,J))");
}


//=========================================================================//
//                             CONSTRUCTORS                                //
//=========================================================================//

/// @name Constructors

/// Builds an empty sparse matrix
/// \code{.cpp}
///    smatrix<> As;
///    disp(As);
/// \endcode
template<typename T>
smatrix<T>::smatrix():
m_row(0), m_col(0) {};

/// Builds a sparse matrix from a scalar
/// \code{.cpp}
///    smatrix<> As = 1;
///    disp(As);
/// \endcode
template<typename T>
smatrix<T>::smatrix(T v):
m_row(1), m_col(1), m_ind(1,0), m_val(1,v) {check();}

/// Builds a sparse matrix from dimensions
/// \code{.cpp}
///    smatrix<> As(3,4);
///    disp(As);
/// \endcode
template<typename T>
smatrix<T>::smatrix(std::size_t m, std::size_t n):
m_row(m), m_col(n) {};

/// Builds a sparse matrix from dimensions and prepare non zeros values
/// \code{.cpp}
///    smatrix<> As(3,4,7);
///    disp(As);
/// \endcode
template<typename T>
smatrix<T>::smatrix(std::size_t m, std::size_t n, std::size_t nnz)
{
    m_row = m;
    m_col = n;
    m_ind = std::vector<std::size_t>(nnz,0);
    m_val = std::vector<T>(nnz,0);
}

/// Builds sparse matrix using linear indexing and values
/// \code{.cpp}
///    smatrix<> As(2,3,{1,0,3,2,5,4},std::vector<double>({0.,0.,1.,2.,0.,3.}));
///    disp(Bs);
/// \endcode
template<typename T>
template<typename S>
smatrix<T>::smatrix(std::size_t m, std::size_t n,
                    std::vector<std::size_t>const& L, std::vector<S>const& V)
{
    if (L.size()!=V.size())
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    m_row = m;
    m_col = n;
    if (L.size()>0)
    {
        matrix<std::size_t> Ia, Ib, A=L;
        std::tie(Ia,Ib) = argunique(A);
        m_ind = eval(A(Ia)).val();
        if (m_ind[m_ind.size()-1]>=m_row*m_col)
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
        }
        m_val = std::vector<T>(numel(Ia),0);
        for (std::size_t l=0; l<numel(A); ++l)
        {
            m_val[Ib(l)] += V[l];
        }
        check();
    }
}

/// Builds sparse matrix using bi-linear indexing and values
/// \code{.cpp}
///    smatrix<> As(2,3,{0,0,1,1},{0,1,0,1},std::vector<double>({1,0,0,2}));
///    disp(As);
/// \endcode
template<typename T>
template<typename S>
smatrix<T>::smatrix(std::size_t m, std::size_t n,
                    std::vector<std::size_t>const& I, std::vector<std::size_t>const& J, std::vector<S>const& V)
{
    if (I.size()!=J.size() || I.size()!=V.size())
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    std::vector<std::size_t> L(I.size(),0);
    for (std::size_t l=0; l<L.size(); ++l)
    {
        if (I[l]>=m || J[l]>=n)
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
        }
        else
        {
            L[l] = I[l]*n+J[l];
        }
    }
    (*this) = smatrix(m,n,L,V);
}

/// Builds a sparse matrix from dense matrix.
/// \code{.cpp}
///    smatrix<> As = M_PI*eye(3,4);
///    disp(As);
/// \endcode
template<typename T>
template<typename S>
smatrix<T>::smatrix(matrix<S>const& A)
{
    m_row = A.size(1);
    m_col = A.size(2);
    for (std::size_t l=0; l<numel(A); ++l)
    {
        m_ind.push_back(l);
        m_val.push_back((T)A(l));
    }
    check();
}

/// Builds a sparse matrix from another sparse matrix.
/// \code{.cpp}
///    smatrix<int> As = eye(3,4);
///    smatrix<> Bs(As);
///    disp(Bs);
/// \endcode
template<typename T>
template<typename S>
smatrix<T>::smatrix(smatrix<S>const& As)
{
    m_row = As.size(1);
    m_col = As.size(2);
    m_ind = std::vector<std::size_t>(As.nnz());
    m_val = std::vector<T>(As.nnz());
    for (std::size_t k=0; k<As.nnz(); ++k)
    {
        m_ind[k] = As.ind(k);
        m_val[k] = As.val(k);
    }
    check();
}


//==========================================================================//
//                         INTERNAL OPERATORS                               //
//==========================================================================//

/// @name Internal operators

/// Converts a (1x1) to fundamental type variable.
/// \code{.cpp}
///    smatrix<int> A(3);
///    a = double(A);
///    disp(a);
/// \endcode
template<typename T>
template<typename S>
smatrix<T>::operator S() const
{
    if (m_row*m_col!=1)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Can't convert implicitly matrix to value because does not have a single element.");
    }
    return m_val[0];
}

/// Performs addition assignment.\n
/// If A is fundamental type variable, add value of A to each element of the sparse matrix
/// If A is a sparse matrix, add each element of A to each element of the sparse matrix,
/// in this case matrix dimensions must agree.
/// \code{.cpp}
///    smatrix<float> As = speye(3,4);
///    As += 1;
///    disp(As);
///    As += spones(3,4);
///    disp(As);
/// \endcode
template<typename T>
smatrix<T>& smatrix<T>::operator+=(smatrix<T>const& As)
{
    if (numel(As)==1)
    {
        for (std::size_t l=0; l<m_row*m_col; ++l)
        {
            m_ind.push_back(l);
            m_val.push_back(As.val(0));
        }
    }
    else if (m_row==As.size(1) && m_col==As.size(2))
    {
        for (std::size_t l=0; l<As.nnz(); ++l)
        {
            m_ind.push_back(As.ind(l));
            m_val.push_back(As.val(l));
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    (*this) = smatrix(m_row,m_col,m_ind,m_val);
    check();
    return (*this);
}

/// Performs substraction assignment.\n
/// If As is fundamental type variable, substract value of As to each element of the sparse matrix
/// If As is a sparse matrix, substract each element of As to each element of the sparse matrix,
/// in this case matrix dimensions must agree.
/// \code{.cpp}
///    smatrix<float> As = speye(3,4);
///    As -= 1;
///    disp(As);
///    As -= spones(3,4);
///    disp(As);
/// \endcode
template<typename T>
smatrix<T>& smatrix<T>::operator-=(smatrix<T>const& As)
{
    if (numel(As)==1)
    {
        for (std::size_t l=0; l<m_row*m_col; ++l)
        {
            m_ind.push_back(l);
            m_val.push_back(-As.val(0));
        }
    }
    else if (m_row==As.size(1) && m_col==As.size(2))
    {
        for (std::size_t l=0; l<As.nnz(); ++l)
        {
            m_ind.push_back(As.ind(l));
            m_val.push_back(-As.val(l));
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    (*this) = smatrix(m_row,m_col,m_ind,m_val);
    check();
    return (*this);
}

/// Performs multiplication assignment.\n
/// If As is fundamental type variable, multiply value of As to each element of the sparse matrix
/// If As is a sparse matrix, multiply each element of As to each element of the sparse matrix,
/// in this case matrix dimensions must agree.
/// \code{.cpp}
///    smatrix<float> As = speye(3,4);
///    As *= 2;
///    disp(As);
/// \endcode
template<typename T>
smatrix<T>& smatrix<T>::operator*=(smatrix<T>const& As)
{
    if (numel(As)==1)
    {
        for (std::size_t k=0; k<nnz(); ++k) {m_val[k] *= As.val(0);}
    }
    else if (m_row==As.size(1) && m_col==As.size(2))
    {
        for (std::size_t k=0; k<nnz(); ++k) {m_val[k] *= As(m_ind[k]);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    check();
    return (*this);
}

/// Performs division assignment.\n
/// If As is fundamental type variable, divide value of As to each element of the sparse matrix
/// If As is a sparse matrix, divide each element of As to each element of the sparse matrix,
/// in this case matrix dimensions must agree.
/// \code{.cpp}
///    smatrix<float> As = speye(3,4);
///    As /= 2;
///    disp(As);
/// \endcode
template<typename T>
smatrix<T>& smatrix<T>::operator/=(smatrix<T>const& As)
{
    if (numel(As)==1)
    {
        for (std::size_t k=0; k<nnz(); ++k) {m_val[k] /= As.val(0);}
    }
    else if (m_row==As.size(1) && m_col==As.size(2))
    {
        (*this) = sparse(full(*this)/full(As));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    check();
    return (*this);
}

/// Returns an element of the sparse matrix by giving the linear indexing of
/// element (const accessor).\n
/// Remark : sparse matrix values are stored in matrix<T> object, associated
/// to matrix<std:size_t> both for row and column indexing.
///
/// \code{.cpp}
///    smatrix<float> As = speye(3,4);
///    float coef = A(11);
///    disp(coef);
/// \endcode
///
/// WARNING : Sparse indexing is relatively slow, please prefer builders or
/// constructors instead.
template<typename T>
const T& smatrix<T>::operator()(std::size_t l) const
{
    if (l>=m_row*m_col)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
    }
    std::size_t n=nnz();
    if (n==0 || l < m_ind[0] || l > m_ind[n-1])
    {
        return M_ZERO<T>;
    }
    else
    {
        std::size_t a=1, b=n, c;
        while (b-a>1)
        {
            c = (a+b)/2;
            (l<m_ind[c-1])? (b=c) : (a=c);
        }
        if (m_ind[a-1]==l) {return m_val[a-1];}
        else if (m_ind[a]==l) {return m_val[a];}
        else {return M_ZERO<T>;}
    }
}

/// Returns an element of the sparse matrix by giving the linear indexing of
/// element (non const accessor).\n
/// Remark : sparse matrix values are stored in matrix<T> object, associated
/// to matrix<std:size_t> both for row and column indexing.
///
/// \code{.cpp}
///    smatrix<float> As = speye(3,4);
///    // access to last element
///    As(11) = 10;
///    disp(As(11));
/// \endcode
///
/// WARNING : Sparse indexing is relatively slow, please prefer builders or
/// constructors instead.
template<typename T>
T& smatrix<T>::operator()(std::size_t l)
{
    if (l>=m_row*m_col)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
    }
    std::size_t n=nnz(), k;
    if (n==0)
    {
        m_ind.push_back(l);
        m_val.push_back(0);
        k = 0;
    }
    else if (l < m_ind[0])
    {
        m_ind.insert(m_ind.begin(),l);
        m_val.insert(m_val.begin(),0);
        k = 0;
    }
    else if (l > m_ind[n-1])
    {
        m_ind.push_back(l);
        m_val.push_back(0);
        k = n;
    }
    else
    {
        std::size_t a=1, b=n, c;
        while (b-a>1)
        {
            c = (a+b)/2;
            (l<m_ind[c-1])? (b=c) : (a=c);
        }
        if (m_ind[a-1]==l) {k=a-1;}
        else if (m_ind[a]==l) {k=a;}
        else
        {
            m_ind.insert(m_ind.begin()+a,l);
            m_val.insert(m_val.begin()+a,0);
            k=a;
        }
    }
    return m_val[k];
}

/// Returns a reference to a submatrix corresponding to linear indexing L (const accessor).
/// \code{.cpp}
///    const smatrix<float> As = speye(3,4);
///    // access to element equal to 1
///    auto Bs = As({0,5,10});
///    disp(eval(B));
/// \endcode
template<typename T>
const scview<T> smatrix<T>::operator()(matrix<std::size_t>const& L) const
{
    return scview<T>(*this,L);
}

/// Returns a reference to a submatrix corresponding to linear indexing L (non-const accessor).
/// \code{.cpp}
///    smatrix<float> As = speye(3,4);
///    // access to element equal to 1
///    auto Bs = As({0,5,10});
///    disp(eval(Bs));
/// \endcode
template<typename T>
sview<T> smatrix<T>::operator()(matrix<std::size_t>const& L)
{
    return sview<T>(*this,L);
}

/// Returns an element of the sparse matrix by giving the bi-linear indexing of
/// element (const accessor).
///
/// \code{.cpp}
///    const smatrix<int> As = speye(3,4);
///    int coef = As(2,3);
///    disp(coef);
/// \endcode
///
/// WARNING : Sparse indexing is relatively slow, please prefer builders or
/// constructors instead.
template<typename T>
const T& smatrix<T>::operator()(std::size_t i, std::size_t j) const
{
    if (i>=m_row || j>=m_col)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Indices exceed array bounds.");
    }
    return (*this)(i*m_col+j);
}

/// Returns an element of the sparse matrix by giving the bi-linear indexing
/// of element (non-const accessor).
///
/// \code{.cpp}
///    smatrix<int> As = speye(3,4);
///    As(2,3) = 5;
///    disp(As(2,3));
/// \endcode
///
/// WARNING : Sparse indexing is relatively slow, please prefer builders or
/// constructors instead.
template<typename T>
T& smatrix<T>::operator()(std::size_t i, std::size_t j)
{
    if (i>=m_row || j>=m_col)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Indices exceed array bounds.");
    }
    return (*this)(i*m_col+j);
}

/// Returns a reference to a submatrix corresponding to bilinear indexing
/// I for rows and J for columns (const accessor).
/// \code{.cpp}
///    smatrix<float> As = speye(3,4);
///    auto Bs = As({0,2}, {0,2});
///    disp(eval(Bs));
/// \endcode
template<typename T>
const scview<T> smatrix<T>::operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J) const
{
    return scview<T>(*this,I,J);
}

/// Return a reference to a submatrix corresponding to bilinear indexing
/// I for rows and J for columns (non-const accessor).
/// \code{.cpp}
///    smatrix<float> As = speye(3,4);
///    auto Bs = As({0,2}, {0,2});
///    disp(eval(Bs));
/// \endcode
template<typename T>
sview<T> smatrix<T>::operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J)
{
    return sview<T>(*this,I,J);
}


//==========================================================================//
//                           INTERNAL TOOLS                                 //
//==========================================================================//
//==========================================================================
// [matrix.check]
/// check(As) removes all zeros values and associated indices and verify if
/// linear indices are well sorted in ascending order.
///
/// \code{.cpp}
///    smatrix<> Xs(range(0,12),eye(1,12),3,4);
///    disp(Xs);
///    check(Xs);
///    disp(Xs);
/// \endcode
template<typename T>
void smatrix<T>::check()
{
    std::size_t k=0, nz=nnz();
    T zero = 0;
    while (k<nz)
    {
        if (m_val[k]==zero)
        {
            m_ind.erase(m_ind.begin()+k);
            m_val.erase(m_val.begin()+k);
            --nz;
        }
        else {++k;}
    }
    if (!std::is_sorted(m_ind.begin(),m_ind.end()))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Sparse matrix indices are not sorted.");
    }
    if (nz>0)
    {
        if (m_ind[nz-1]>=m_row*m_col)
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
        }
    }
}

//==========================================================================
// [matrix.clear]
/// Removes all values of sparse matrix A and fix size to 0x0.\n
/// Clear should be used to free memory without deleting matrix object.
/// 
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    disp(As);
///    As.clear;
///    disp(As);
/// \endcode
template<typename T>
void smatrix<T>::clear()
{
    m_row = 0;
    m_col = 0;
    m_ind.clear();
    m_val.clear();
}

//==========================================================================
// [smatrix.ind]
/// Linear index of the k-th stored value
/// \code{.cpp}
///    smatrix<> As = {{1,2,0},{0,5,6}};
///    disp(As.ind(0));
/// \endcode
template<typename T>
std::size_t const& smatrix<T>::ind(std::size_t k) const
{
#ifdef DEBUG
    if (k>=nnz())
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
    }
#endif
    return m_ind[k];
}
template<typename T>
std::size_t& smatrix<T>::ind(std::size_t k)
{
#ifdef DEBUG
    if (k>=nnz())
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
    }
#endif
    return m_ind[k];
}
template<typename T>
std::vector<std::size_t>const& smatrix<T>::ind() const
{
    return m_ind;
}

//==========================================================================
// [smatrix.nnz]
/// Non zeros values stored
/// \code{.cpp}
///    smatrix<> As = {{1,2,0},{0,5,6}};
///    disp(As.nnz());
/// \endcode
template<typename T>
std::size_t smatrix<T>::nnz() const
{
    return m_val.size();
}

//==========================================================================
// [smatrix.reshape]
/// Modifies the shape of the sparse matrix without changing its values.
/// \code{.cpp}
///    smatrix<> As = speye(3,4);
///    As.reshape(4,3);
///    disp(As);
/// \endcode
template<typename T>
void smatrix<T>::reshape(std::size_t m, std::size_t n)
{
    if (m_row*m_col!=m*n)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"To reshape the number of elements must not change.");
    }
    m_row = m;
    m_col = n;
}

//==========================================================================
// [smatrix.size]
/// If dim=0 (default value) returns the number elements of the sparse matrix
/// else returns the lengths of the specified dimensions dim.
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    disp(As.size(0));
///    disp(As.size(1));
///    disp(As.size(2));
/// \endcode
template<typename T>
std::size_t smatrix<T>::size(int dim) const
{
    std::size_t s=0;
    if (dim==0) {s = m_row*m_col;}
    else if (dim==1) {s = m_row;}
    else if (dim==2) {s = m_col;}
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"When dimension is set, dim=0 count for all elements, dim=1 count for rows and dim=2 for column.");
    }
    return s;
}

//==========================================================================
// [smatrix.val]
/// Value of the k-th stored data.
/// \code{.cpp}
///    smatrix<> As = {{1,2,0},{0,5,6}};
///    disp(As.val(0));
/// \endcode
template<typename T>
T const& smatrix<T>::val(std::size_t k) const
{
#ifdef DEBUG
    if (k>=nnz())
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
    }
#endif
    return m_val[k];
}
template<typename T>
T& smatrix<T>::val(std::size_t k)
{
#ifdef DEBUG
    if (k>=nnz())
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
    }
#endif
    return m_val[k];
}
template<typename T>
std::vector<T>const& smatrix<T>::val() const
{
    return m_val;
}


//==========================================================================//
//                         EXTERNAL OPERATORS                               //
//==========================================================================//
//==========================================================================
// [<<]
/// Flux operator.
///
/// Insert sparse matrix As inside a stream flux respecting standard C++ conventions.
///
/// \code{.cpp}
///    smatrix<> As = speye(3,4);
///    std::cout << As << std::endl;
/// \endcode
///
// \see disp.
template<typename T>
std::ostream& operator<<(std::ostream& flux, smatrix<T>const& As)
{
    disp(As,0,flux,size(As,1));
    return flux;
}

//==========================================================================
// [+]
/// Plus.
///
/// As+Bs adds sparse matrices As and Bs. As and Bs must have compatible sizes.
/// In the simplest cases, they can be the same size or one can be a scalar. Two
/// inputs have compatible sizes if, for every dimension, the dimension
/// sizes of the inputs are the same.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> Bs = {{6,5,4},{3,2,1}};
///    smatrix<> Cs = As+Bs;
///    smatrix<> Ds = As+1;
///    disp(Cs);
///    disp(Ds);
/// \endcode
///
// \see operator-(), operator*(), operator/(), operator+=().
template<typename R, typename S>
auto operator+(smatrix<R>const& As, smatrix<S>const& Bs)
{
    using T = decltype(As(0)+Bs(0));
    smatrix<T> Cs;
    if (numel(As)==1)
    {
        Cs = Bs;
        Cs += As;
    }
    else
    {
        Cs = As;
        Cs += Bs;
    }
    return Cs;
}
template<typename R, typename S>
auto operator+(matrix<R>const& A, smatrix<S>const& Bs)
{
    if (size(A,1)!=size(Bs,1) || size(A,2)!=size(Bs,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(A(0)+Bs(0));
    matrix<T> C = A;
    for (std::size_t k=0; k<nnz(Bs); ++k) {C(Bs.ind(k)) += Bs.val(k);}
    return C;
}
template<typename R, typename S>
auto operator+(smatrix<R>const& As, matrix<S>const& B)
{
    if (size(As,1)!=size(B,1) || size(As,2)!=size(B,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(As(0)+B(0));
    matrix<T> C = B;
    for (std::size_t k=0; k<nnz(As); ++k) {C(As.ind(k)) += As.val(k);}
    return C;
}
template<typename R, typename S>
auto operator+(R a, smatrix<S>const& Bs)
{
    using T = decltype(a+Bs(0));
    matrix<T> C(size(Bs,1),size(Bs,2),a);
    for (std::size_t k=0; k<nnz(Bs); ++k) {C(Bs.ind(k)) += Bs.val(k);}
    return C;
}
template<typename R, typename S>
auto operator+(smatrix<R>const& As, S b)
{
    using T = decltype(As(0)+b);
    matrix<T> C(size(As,1),size(As,2),b);
    for (std::size_t k=0; k<nnz(As); ++k) {C(As.ind(k)) += As.val(k);}
    return C;
}

//==========================================================================
// [-]
/// Unary or binary minus.
///
/// -As negates the elements of As.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> Bs = -As;
///    disp(Bs);
/// \endcode
///
/// As-Bs substracts sparse matrices As and Bs. As and Bs must have compatible sizes.
/// In the simplest cases, they can be the same size or one can be a scalar. Two
/// inputs have compatible sizes if, for every dimension, the dimension
/// sizes of the inputs are the same.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> Bs = {{6,5,4},{3,2,1}};
///    smatrix<> Cs = As-Bs;
///    smatrix<> Ds = As-1;
///    disp(Cs);
///    disp(Ds);
/// \endcode
///
// \see operator+(), operator*(), operator/(), operator-=().
template<typename T>
auto operator-(smatrix<T>const& As)
{
    smatrix<T> Bs = As;
    Bs *= -1;
    return Bs;
}
template<typename R, typename S>
auto operator-(smatrix<R>const& As, smatrix<S>const& Bs)
{
    using T = decltype(As(0)-Bs(0));
    smatrix<T> Cs;
    if (numel(As)==1)
    {
        Cs = -Bs;
        Cs += As;
    }
    else
    {
        Cs = As;
        Cs -= Bs;
    }
    return Cs;
}
template<typename R, typename S>
auto operator-(matrix<R>const& A, smatrix<S>const& Bs)
{
    if (size(A,1)!=size(Bs,1) || size(A,2)!=size(Bs,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(A(0)-Bs(0));
    matrix<T> C = A;
    for (std::size_t k=0; k<nnz(Bs); ++k) {C(Bs.ind(k)) -= Bs.val(k);}
    return C;
}
template<typename R, typename S>
auto operator-(smatrix<R>const& As, matrix<S>const& B)
{
    if (size(As,1)!=size(B,1) || size(As,2)!=size(B,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(As(0)-B(0));
    matrix<T> C = -B;
    for (std::size_t k=0; k<nnz(As); ++k) {C(As.ind(k)) += As.val(k);}
    return C;
}
template<typename R, typename S>
auto operator-(R a, smatrix<S>const& Bs)
{
    using T = decltype(a-Bs(0));
    matrix<T> C(size(Bs,1),size(Bs,2),a);
    for (std::size_t k=0; k<nnz(Bs); ++k) {C(Bs.ind(k)) -= Bs.val(k);}
    return C;
}
template<typename R, typename S>
auto operator-(smatrix<R>const& As, S b)
{
    using T = decltype(As(0)-b);
    matrix<T> C(size(As,1),size(As,2),-b);
    for (std::size_t k=0; k<nnz(As); ++k) {C(As.ind(k)) += As.val(k);}
    return C;
}

//==========================================================================
// [*]
/// Multiply.
///
/// As*Bs multiply sparse matrices As and Bs. As and Bs must have compatible sizes.
/// In the simplest cases, they can be the same size or one can be a scalar. Two
/// inputs have compatible sizes if, for every dimension, the dimension
/// sizes of the inputs are the same.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> Bs = {{6,5,4},{3,2,1}};
///    smatrix<> Cs = As*Bs;
///    smatrix<> Ds = As*2;
///    disp(Cs);
///    disp(Ds);
/// \endcode
///
// \see operator+(), operator-(), operator/(), operator*=().
template<typename R, typename S>
auto operator*(smatrix<R>const& As, smatrix<S>const& Bs)
{
    using T = decltype(As(0)*Bs(0));
    smatrix<T> Cs;
    if (numel(As)==1)
    {
        Cs = Bs;
        Cs *= As;
    }
    else
    {
        Cs = As;
        Cs *= Bs;
    }
    return Cs;
}
template<typename R, typename S>
auto operator*(matrix<R>const& A, smatrix<S>const& Bs)
{
    if (size(A,1)!=size(Bs,1) || size(A,2)!=size(Bs,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(A(0)*Bs(0));
    smatrix<T> Cs = Bs;
    for (std::size_t k=0; k<nnz(Bs); ++k) {Cs.val(k) *= A(Cs.ind(k));}
    return Cs;
}
template<typename R, typename S>
auto operator*(smatrix<R>const& As, matrix<S>const& B)
{
    if (size(As,1)!=size(B,1) || size(As,2)!=size(B,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(As(0)*B(0));
    smatrix<T> Cs = As;
    for (std::size_t k=0; k<nnz(As); ++k) {Cs.val(k) *= B(Cs.ind(k));}
    return Cs;
}
template<typename R, typename S>
inline auto operator*(R As, smatrix<S>const& Bs) {return smatrix<R>(As)*Bs;}
template<typename R, typename S>
inline auto operator*(smatrix<R>const& As, S Bs) {return As*smatrix<S>(Bs);}

//==========================================================================
// [/]
/// Divide.
///
/// As/Bs divide sparse matrices As and Bs. As and Bs must have compatible sizes.
/// In the simplest cases, they can be the same size or one can be a scalar. Two
/// inputs have compatible sizes if, for every dimension, the dimension
/// sizes of the inputs are the same.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> Bs = {{6,5,4},{3,2,1}};
///    smatrix<> Cs = As/Bs;
///    smatrix<> Ds = As/2;
///    disp(Cs);
///    disp(Ds);
/// \endcode
///
// \see operator+(), operator-(), operator*(), operator/=().
template<typename R, typename S>
auto operator/(smatrix<R>const& As, smatrix<S>const& Bs)
{
    using T = decltype(As(0)+Bs(0));
    smatrix<T> Cs;
    if (numel(As)==1)
    {
        Cs = sparse(matrix<T>(size(Bs,1),size(Bs,2),As(0))/full(Bs));
    }
    else
    {
        Cs = As;
        Cs /= Bs;
    }
    return Cs;
}
template<typename R, typename S>
inline auto operator/(matrix<R>const& As, smatrix<S>const& Bs) {return As/full(Bs);}
template<typename R, typename S>
inline auto operator/(smatrix<R>const& As, matrix<S>const& Bs) {return full(As)/Bs;}
template<typename R, typename S>
inline auto operator/(R As, smatrix<S>const& Bs) {return matrix<R>(size(Bs,1),size(Bs,2),As)/Bs;}
template<typename R, typename S>
inline auto operator/(smatrix<R>const& As, S Bs) {return As/smatrix<S>(Bs);}


//==========================================================================//
//                           EXTERNAL TOOLS                                 //
//==========================================================================//
//==========================================================================
// [all]
/// Linear indexing of all elements.
///
/// all(As) return row matrix with all indices of As in linear indexing.
///
/// \code{.cpp}
///    smatrix<> A = eye(3,4);
///    matrix<std::size_t> L = all(A);
///    smatrix<> B = eval(As(L));
///    disp(B);
/// \endcode
///
// \see row, col, get, set, sview, scview.
template<typename T>
matrix<std::size_t> all(smatrix<T>const& As)
{
    return range(0,numel(As));
}

//==========================================================================
// [check]
/// Check sparse matrix data.
///
/// check(As) removes all zeros values and associated indices and verify if
/// linear indices are well sorted in ascending order.
///
/// \code{.cpp}
///    smatrix<> Xs(range(0,12),eye(1,12),3,4);
///    disp(Xs);
///    check(Xs);
///    disp(Xs);
/// \endcode
///
// \see matrix.
template<typename T>
void check(smatrix<T>& As)
{
    As.check();
}

//==========================================================================
// [clear]
/// Clear sparse matrix to free memory.
///
/// clear(As) removes all values and indices of sparse matrix As and fix
/// size to 0x0. clear should be used to free memory without deleting
/// sparse  matrix object.
///
/// \code{.cpp}
///    smatrix<> X = {{1,2,3},{4,5,6}};
///    disp(Xs);
///    clear(Xs);
///    disp(Xs);
/// \endcode
///
// \see matrix.
template<typename T>
void clear(smatrix<T>& As)
{
    As.clear();
}

//==========================================================================
// [col]
/// Bi-linear indexing of column elements.
///
/// col(As) return row matrix with column indices of As in bi-linear indexing.
///
/// \code{.cpp}
///    smatrix<> A = eye(3,4);
///    matrix<std::size_t> J = col(A);
///    smatrix<> B = eval(A(0,J));
///    disp(B);
/// \endcode
///
// \see row, all, get, set, sview, scview.
template<typename T>
matrix<std::size_t> col(smatrix<T>const& As)
{
    return range(0,size(As,2));
}

//==========================================================================
// [disp]
/// Display sparse matrix.
///
/// \code{.cpp}
///    disp(spones(3,4));
///    disp(spones(3,4),2);
///    disp(spones(3,4),2,std::cout);
///    disp(spones(3,4),2,std::cout,12);
/// \endcode
///
// \see operator<<.
template<typename T>
void disp(smatrix<T>const& As, int info=2, std::ostream& flux=std::cout, std::size_t r=3)
{
    // Preparation and informations
    std::size_t m=size(As,1), n=size(As,2), nz=nnz(As);
    double memBytes = nz*sizeof(T);
    if (info>=2)
    {
        flux << "Sparse matrix " << m << "x" << n;
        flux << " of type '" << typeid(T).name() << "'";
        flux << " with " << nz << " elements (";
        if (memBytes<1e3) {flux << memBytes << " B):";}
        else if (memBytes<1e6) {flux << memBytes/1e3 << " KB):";}
        else {flux << memBytes/1e6 << " MB):";}
        flux << std::endl;
    }
    // Empty matrix
    if (nnz(As)==0)
    {
        flux << "-empty-";
    }
    else
    {
        if (nz<=2*r)
        {
            for (std::size_t l=0; l<nz; ++l)
            {
                flux << "(" << As.ind(l)/n << "," << As.ind(l)%n << ")";
                flux << "  " << As.val(l);
                if (l<nz-1) {flux << std::endl;}
            }
        }
        else
        {
            for (std::size_t l=0; l<r; ++l)
            {
                flux << "(" << As.ind(l)/n << "," << As.ind(l)%n << ")";
                flux << "  " << As.val(l) << std::endl;
            }
            flux << "..." << std::endl;
            for (std::size_t l=nz-r; l<nz; ++l)
            {
                flux << "(" << As.ind(l)/n << "," << As.ind(l)%n << ")";
                flux << "  " << As.val(l);
                if (l<nz-1) {flux << std::endl;}
            }
        }
    }
    if (info>=1) {flux << std::endl;}
}

//==========================================================================
// [find]
/// Find bi-linear indices of nonzero elements.
///
/// (I,J,V) = find(As) returns the bi-linear indices corresponding to the
/// nonzero entries of the array As. As may be a logical expression.
///
/// Use index(As) or values(As) to find respectively linear index or values.
///
/// \code{.cpp}
///    smatrix<> As = eye(3,4);
///    smatrix<std::size_t> Ls = find(A);
///    disp(Ls);
/// \endcode
///
// \see index, values.
template<typename T>
auto find(smatrix<T>const& As)
{
    matrix<std::size_t> I, J;
    std::tie(I,J) = ind2sub(size(As),index(As));
    return std::make_tuple(I,J,values(As));
}

//==========================================================================
// [full]
/// Dense conversion.
///
/// A = full(As) convert sparse matrix to dense matrix.
///
/// \code{.cpp}
///    smatrix<> As = speye(3,4);
///    matrix<>  A  = full(As);
///    disp(A);
/// \endcode
///
// \see sparse.
template<typename T>
matrix<T> full(smatrix<T>const& As)
{
    matrix<T> A = zeros(size(As));
    for (std::size_t k=0; k<nnz(As); ++k)
    {
        A(As.ind(k)) = As.val(k);
    }
    return A;
}
template<typename T>
inline matrix<T>const& full(matrix<T>const& A) {return A;}

//==========================================================================
// [get]
/// Get sub-sparse matrix.
///
/// get(As,L) return sparse matrix As(L) with elements taken from As
/// corresponding to linear indexing L.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> Bs = get(As,{0,2,4});
///    disp(Bs);
/// \endcode
///
/// get(As,I,J) return sparse matrix As(I,J) with elements taken from As
/// corresponding to bilinear indexing I for rows and J for columns.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> Bs = get(As,{0,1},{1,2});
///    disp(Bs);
/// \endcode
///
// \see set, all, row, col, sview, scview.
template<typename T>
smatrix<T> get(smatrix<T>const& As, matrix<std::size_t>const& L)
{
    std::vector<std::size_t> ind;
    std::vector<T> val;
    T tmp, zero=0;
    for (std::size_t l=0; l<numel(L); ++l)
    {
        tmp = As(L(l));
        if (tmp!=zero)
        {
            ind.push_back(l);
            val.push_back(tmp);
        }
    }
    return smatrix<T>(size(L,1),size(L,2),ind,val);
}
template<typename T>
smatrix<T> get(smatrix<T>const& As, matrix<std::size_t>const& I, matrix<std::size_t>const& J)
{
    if (std::min(size(I,1),size(I,2))!=1 || std::min(size(J,1),size(J,2))!=1)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Array indices must be vectors.");
    }
    std::vector<std::size_t> ind;
    std::vector<T> val;
    T tmp, zero=0;
    if (numel(I)*numel(J)*std::log2((float)nnz(As)) < numel(I)+numel(J)+nnz(As))
    {
        for (std::size_t i=0; i<numel(I); ++i)
        {
            for (std::size_t j=0; j<numel(J); ++j)
            {
                tmp = As(I(i),J(j));
                if (tmp!=zero)
                {
                    ind.push_back(i*numel(J)+j);
                    val.push_back(tmp);
                }
            }
        }
    }
    else
    {
        std::vector<std::vector<std::size_t>> Ib(size(As,1),std::vector<std::size_t>());
        for (std::size_t i=0; i<numel(I); ++i)
        {
            Ib[I(i)].push_back(i);
        }
        std::vector<std::vector<std::size_t>> Jb(size(As,2),std::vector<std::size_t>());
        for (std::size_t j=0; j<numel(J); ++j)
        {
            Jb[J(j)].push_back(j);
        }
        std::size_t i, j, n=size(As,2), nz=nnz(As), nn=numel(J);
        for (std::size_t k=0; k<nz; ++k)
        {
            i = As.ind(k)/n;
            j = As.ind(k)%n;
            for (std::size_t ii=0; ii<Ib[i].size(); ++ii)
            {
                for (std::size_t jj=0; jj<Jb[j].size(); ++jj)
                {
                    ind.push_back(Ib[i][ii]*nn+Jb[j][jj]);
                    val.push_back(As.val(k));
                }
            }
        }
    }
    return smatrix<T>(numel(I),numel(J),ind,val);
}

//==========================================================================
// [gmres]
/// Generalized Minimum Residual Method.
///
/// X = gmres(As,B,TOL,MAXIT,AsM1,X0) attempts to solve the system of linear
/// equations As*X=B for X. The N-by-N coefficient hmatrix Ah must be square
/// and the right hand side B a N-by-NRHS matrix (MGCR implementation).
/// Default values are :
///   TOL   = 1e-6;
///   MAXIT = 10;
///   AsM1  = Id;
///   X0    = 0;
///
/// \code{.cpp}
///    smatrix<> As = rand(3,3);
///    matrix<>   B = eye(3,3);
///    matrix<>   X = gmres(As,B);
///    disp(mtimes(X,As));
/// \endcode
///
// \see mtimes, tgemm.
template<typename T>
matrix<T> gmres(smatrix<T>const& As, matrix<T>const& B,
                double tol = 1e-6, std::size_t maxit = 10,
                std::function<matrix<T>(matrix<T>const&)>const& Am1 = std::function<matrix<T>(matrix<T>const&)>(),
                matrix<T>const& X0 = matrix<T>())
{
    std::function<matrix<T>(matrix<T>const&)> Afct;
    Afct = [&As](matrix<T>const& X) {return mtimes(As,X);};
    return gmres(Afct,B,tol,maxit,Am1,X0);
}
template<typename T>
matrix<T> gmres(smatrix<T>const& As, matrix<T>const& B, double tol, std::size_t maxit,
                smatrix<T>const& Asm1, matrix<T>const& X0 = matrix<T>())
{
    std::function<matrix<T>(matrix<T>const&)> Afct, Am1fct;
    Afct   = [&As](matrix<T>const& X) {return mtimes(As,X);};
    Am1fct = [&Asm1](matrix<T>const& X) {return mtimes(Asm1,X);};
    return gmres(Afct,B,tol,maxit,Am1fct,X0);
}

//==========================================================================
// [index]
/// Copy indices in a matrix.
///
/// I = index(As) returns matrix with all linear indices.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    matrix<std::size_t> val = index(As);
///    disp(val);
/// \endcode
///
// \see values, find.
template<typename T>
matrix<std::size_t> index(smatrix<T>const& As)
{
    return matrix<std::size_t>(As.ind());
}

//==========================================================================
// [length]
/// Length of sparse vector.
///
/// length(Vs) returns the length of sparse vector V. For sparse matrix A,
/// it is equivalent to max(size(A)).
///
/// \code{.cpp}
///    smatrix<> As = spones(1,3);
///    smatrix<> Bs = spones(4,3);
///    disp(length(As));
///    disp(length(Bs));
/// \endcode
///
// \see numel, size.
template<typename T>
std::size_t length(smatrix<T>const& As)
{
    return std::max(size(As,1),size(As,2));
}

//==========================================================================
// [mtimes]
/// Matrix multiply.
///
/// mtimes(As,B) is the matrix product of As (sparse) and B (dense),
/// the number of columns of As must equal the number of rows of B.
///
/// mtimes(A,Bs) is the matrix product of A (dense) and Bs (sparse),
/// the number of columns of A must equal the number of rows of Bs.
///
/// mtimes(As,Bs) is the matrix product of A (sparse) and Bs (sparse),
/// the number of columns of As must equal the number of rows of Bs.
/// This is a naive implementation.
///
/// \code{.cpp}
///    smatrix<> As = speye(3);
///    matrix<>  B  = ones(3,1);
///    matrix<>  C  = mtimes(As,B);
///    disp(C);
/// \endcode
///
// \see tgemm.
template<typename R, typename S>
auto mtimes(smatrix<R>const& As, matrix<S>const& B)
{
    if (size(As,2)!=size(B,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(As(0)*B(0));
    matrix<T> C(size(As,1),size(B,2));
    std::size_t i, k;
    for (std::size_t l=0; l<nnz(As); ++l)
    {
        i = As.ind(l)/size(As,2);
        k = As.ind(l)%size(As,2);
        for (std::size_t j=0; j<size(B,2); ++j)
        {
            C(i,j) += As.val(l)*B(k,j);
        }
    }
    return C;
}
template<typename R, typename S>
auto mtimes(matrix<R>const& A, smatrix<S>const& Bs)
{
    if (size(A,2)!=size(Bs,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(A(0)*Bs(0));
    matrix<T> C(size(A,1),size(Bs,2));
    std::size_t k, j;
    for (std::size_t l=0; l<nnz(Bs); ++l)
    {
        k = Bs.ind(l)/size(Bs,2);
        j = Bs.ind(l)%size(Bs,2);
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            C(i,j) += A(i,k)*Bs.val(l);
        }
    }
    return C;
}
template<typename R, typename S>
auto mtimes(smatrix<R>const& As, smatrix<S>const& Bs)
{
    if (size(As,2)!=size(Bs,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(As(0)*Bs(0));
    smatrix<T> Cs(size(As,1),size(Bs,2));
    tgemm((T)1,As,Bs,(T)0,Cs);
    return Cs;
}

//==========================================================================
// [nnz]
/// Number of nonzero matrix elements.
///
/// nz = nnz(As) is the number of nonzero elements in As.
///
/// \code{.cpp}
///    smatrix<>   As = speye(3,4);
///    std::size_t nz = nnz(As);
///    disp(nz);
/// \endcode
///
// \see find, size.
template<typename T>
std::size_t nnz(smatrix<T>const& As)
{
    return As.nnz();
}

//==========================================================================
// [numel]
/// Number of elements in a sparse array (including zeros).
///
/// N = numel(As) returns the number of elements N of sparse array As
/// (including zeros), equivalent to prod(size(As)).
///
/// \code{.cpp}
///    smatrix<>   As = speye(3,4);
///    std::size_t n  = numel(As);
///    disp(n);
/// \endcode
///
// \see size, nnz.
template<typename T>
std::size_t numel(smatrix<T>const& As)
{
    return As.size();
}

//==========================================================================
// [reshape]
/// Reshape sparse matrix.
///
/// reshape(As,M,N) returns the M-by-N sparse matrix whose elements are taken
/// from A. An error results if A does not have M*N elements.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> Bs = reshape(A,3,2);
///    disp(Bs);
/// \endcode
///
// \see transpose.
template<typename T>
smatrix<T> reshape(smatrix<T>const& As, std::size_t m, std::size_t n)
{
    smatrix<T> Bs = As;
    Bs.reshape(m,n);
    return Bs;
}

//==========================================================================
// [row]
/// Bi-linear indexing of row elements.
///
/// row(As) return row matrix with row indices of As in bi-linear indexing.
///
/// \code{.cpp}
///    smatrix<> A = eye(3,4);
///    matrix<std::size_t> I = row(A);
///    smatrix<> B = eval(A(I,0));
///    disp(B);
/// \endcode
///
// \see all, col, get, set, sview, scview.
template<typename T>
matrix<std::size_t> row(smatrix<T>const& As)
{
    return range(0,size(As,1));
}

//==========================================================================
// [set]
/// Set sub-sparse matrix.
///
/// set(As,L,Bs) modify As with elements taken from Bs corresponding to
/// linear indexing L, as As(L) = Bs.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    set(As,{0,2,4},-spones(1,3));
///    disp(As);
/// \endcode
///
/// set(As,I,J,Bs) modify As with elements taken from Bs corresponding to
/// bilinear indexing I for rows and J for columns, as As(I,J) = Bs.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    set(As,{0,1},{1,2},speye(2));
///    disp(As);
/// \endcode
///
// \see get, all, row, col, sview, scview.
template<typename T, typename U>
void set(smatrix<T>& As, matrix<std::size_t>const& L, smatrix<U>const& Bs)
{
    if (size(L,1)!=size(Bs,1) || size(L,2)!=size(Bs,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions for indices and values must agree.");
    }
    U tmp;
    for (std::size_t l=0; l<numel(L); ++l)
    {
        tmp = Bs(l);
        if (tmp!=(U)0.)
        {
            As(L(l)) = tmp;
        }
        else if (As(L(l))!=(T)0.)
        {
            As(L(l)) = 0.;
        }
    }
    check(As);
}
template<typename T, typename U>
void set(smatrix<T>& As, matrix<std::size_t>const& I, matrix<std::size_t>const& J, smatrix<U>const& Bs)
{
    if (std::min(size(I,1),size(I,2))!=1 || std::min(size(J,1),size(J,2))!=1)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Subscript indices must be vectors.");
    }
    if (numel(I)!=size(Bs,1) || numel(J)!=size(Bs,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions for indices and values must agree.");
    }
    if (numel(I)*numel(J) < numel(I)+numel(J)+nnz(As)+nnz(Bs))
    {
        for (std::size_t i=0; i<numel(I); ++i)
        {
            for (std::size_t j=0; j<numel(J); ++j)
            {
                if (Bs(i,j)!=(U)0.)
                {
                    As(I(i),J(j)) = Bs(i,j);
                }
                else if (As(I(i),J(j))!=(T)0.)
                {
                    As(I(i),J(j)) = 0.;
                }
            }
        }
    }
    else
    {
        matrix<logical> Il(size(As,1),1); Il(I) = true;
        matrix<logical> Jl(1,size(As,2)); Jl(J) = true;
        std::size_t i, j;
        for (std::size_t k=0; k<nnz(As); ++k)
        {
            i = As.ind(k)/size(As,2);
            j = As.ind(k)%size(As,2);
            if (Il(i) && Jl(j)) {As.val(k)=(T)0.;}
        }
        for (std::size_t k=0; k<nnz(Bs); ++k)
        {
            i = Bs.ind(k)/size(Bs,2);
            j = Bs.ind(k)%size(Bs,2);
            As(I(i),J(j)) = Bs.val(k);
        }
    }
    check(As);
}

//==========================================================================
// [size]
/// Size of sparse matrix.
///
/// S = size(As) for m-by-n sparse matrix As returns the two-element vector [m,n]
/// containing the number of rows and columns in the matrix.
///
/// S = size(As,dim) returns the lengths of the specified dimensions dim.
///
/// \code{.cpp}
///    smatrix<> As = speye(3,4);
///    disp(size(As));
///    disp(size(As,1));
///    disp(size(As,2));
/// \endcode
///
// \see numel, nnz.
template<typename T>
matrix<std::size_t> size(smatrix<T>const& As)
{
    return {As.size(1),As.size(2)};
}
template<typename T>
std::size_t size(smatrix<T>const& As, int dim)
{
    return As.size(dim);
}

//==========================================================================
// [sparse]
/// Sparse conversion.
///
/// As = sparse(A) convert dense to sparse matrix.
///
/// S = sparse(L,V,m,n) uses vectors L and V to generate an m-by-n sparse
/// matrix S, such that S(L(k)) = V(k). Vectors L and V are all the
/// same length. Any elements of V that have duplicate values of L are
/// added together.
///
/// S = sparse(I,J,V,m,n) uses vectors I, J, and V to generate an m-by-n sparse
/// matrix S, such that S(I(k),J(k)) = V(k). Vectors I, J, and V are all the
/// same length. Any elements of V that have duplicate values of I and J are
/// added together.
///
/// S = sparse(I,J,V) use m = max(I)+1 and n = max(J)+1.
///
/// \code{.cpp}
///    matrix<>  A  = eye(3,4);
///    smatrix<> As = sparse(A);
///    disp(As);
/// \endcode
///
// \see full.
template<typename T>
smatrix<T> sparse(matrix<T>const& A)
{
    return smatrix<T>(A);
}
template<typename T>
smatrix<T> sparse(matrix<std::size_t>const& L, matrix<T>const& V, std::size_t m, std::size_t n)
{
    return smatrix<T>(m,n,L.val(),V.val());
}
template<typename T>
smatrix<T> sparse(matrix<std::size_t>const& I, matrix<std::size_t>const& J, matrix<T>const& V)
{
    return smatrix<T>(max(I)+1,max(J)+1,I.val(),J.val(),V.val());
}
template<typename T>
smatrix<T> sparse(matrix<std::size_t>const& I, matrix<std::size_t>const& J, matrix<T>const& V,
		  std::size_t m, std::size_t n)
{
    return smatrix<T>(m,n,I.val(),J.val(),V.val());
}
template<typename T>
inline smatrix<T>const& sparse(smatrix<T>const& As) {return As;}

//==========================================================================
// [speye]
/// Identity sparse matrix.
///
/// speye(N) is the N-by-N identity sparse matrix.
///
/// speye(M,N) or speye({M,N}) is an M-by-N sparse matrix with 1's on the
/// diagonal and zeros elsewhere.
///
/// \code{.cpp}
///    smatrix<> As = speye(2,3);
///    disp(As);
/// \endcode
///
// \see spzeros, spones, sprand, eye.
template<typename T=double>
smatrix<T> speye(std::size_t m, long n=-1)
{
    if (n==-1) {n=m;}
    std::vector<std::size_t> I(std::min<std::size_t>(m,n));
    for (std::size_t l=0; l<I.size(); ++l) {I[l]=l;}
    std::vector<T> V(I.size(),1);
    return smatrix<T>(m,n,I,I,V);
}
template<typename T=double>
smatrix<T> speye(matrix<std::size_t>const& S)
{
    return speye<T>(S(0),S(1));
}

//==========================================================================
// [spones]
/// Ones sparse matrix.
///
/// spones(N) is an N-by-N sparse matrix of ones.
///
/// spones(M,N) or spones({M,N}) is an M-by-N sparse matrix of ones.
///
/// \code{.cpp}
///    matrix<> As = spones(2,3);
///    disp(As);
/// \endcode
///
// \see spzeros, speye, sprand, ones.
template<typename T=double>
smatrix<T> spones(std::size_t m, long n=-1)
{
    return ones<T>(m,n);
}
template<typename T=double>
smatrix<T> spones(matrix<std::size_t>const& S)
{
    return spones<T>(S(0),S(1));
}

//==========================================================================
// [sprand]
/// Sparse matrix of uniformly distributed pseudorandom numbers.
///
/// sprand(N) returns an N-by-N sparse matrix containing pseudorandom values
/// drawn from the standard uniform distribution on the open interval (0,1).
///
/// sprand(M,N) or sprand({M,N}) returns an M-by-N sparse matrix.
///
/// sprand(M,N,true) use standard time to initialize seed.
///
/// \code{.cpp}
///    smatrix<> As = sprand(2,3);
///    smatrix<> Bs = sprand(2,3,true);
///    disp(As);
///    disp(Bs);
/// \endcode
///
// \see spzeros, speye, spones, rand
template<typename T=double>
smatrix<T> sprand(std::size_t m, long n=-1, bool seed=false)
{
    return rand<T>(m,n,seed);
}
template<typename T=double>
smatrix<T> sprand(matrix<std::size_t>const& S, bool seed=false)
{
    return sprand<T>(S(0),S(1),seed);
}

//==========================================================================
// [spzeros]
/// Zeros sparse matrix.
///
/// spzeros(N) is an N-by-N sparse matrix of zeros.
///
/// spzeros(M,N) or spzeros({M,N}) is an M-by-N sparse matrix of zeros.
///
/// \code{.cpp}
///    matrix<> As = spzeros(2,3);
///    disp(As);
/// \endcode
///
// \see spones, speye, sprand, zeros.
template<typename T=double>
smatrix<T> spzeros(std::size_t m, long n=-1)
{
    if (n==-1) {n=m;}
    return smatrix<T>(m,n);
}
template<typename T=double>
smatrix<T> spzeros(matrix<std::size_t>const& S, bool seed=false)
{
    return spzeros<T>(S(0),S(1));
}

//==========================================================================
// [tgemm]
/// In-place sparse matrix product.
///
/// tgemm(alpha,As,Bs,beta,Cs) performs the in-place sparse matrix-matrix operations
///    C = alpha*As*Bs + beta*Cs,
/// where alpha, beta are scalars and As, Bs, Cs are matrices with compatible size.
///
/// NOTE : This is a naive implementation, you can use external sparse solver to
/// get better performance.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> Bs = speye(3,3);
///    smatrix<> Cs = spones(2,3);
///    tgemm(1,As,Bs,1,Cs);
///    disp(Cs);
/// \endcode
///
// \see mtimes.
template<typename P, typename Q, typename R, typename S, typename T>
void tgemm(P alpha, smatrix<Q>const& As, smatrix<R>const& Bs, S beta, smatrix<T>& Cs)
{
    if (size(As,2)!=size(Bs,1) || size(As,1)!=size(Cs,1) || size(Bs,2)!=size(Cs,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    T a=alpha, b=beta;
    for (std::size_t lc=0; lc<nnz(Cs); ++lc)
    {
        Cs.val(lc) *= b;
    }
    std::size_t n=size(Cs,2), p=size(As,2), i, j, k, lb;
    std::vector<std::vector<std::size_t>> Lb(p,std::vector<std::size_t>());
    for (std::size_t lb=0; lb<nnz(Bs); ++lb)
    {
        k = Bs.ind(lb)/n;
        Lb[k].push_back(lb);
    }
    for (std::size_t la=0; la<nnz(As); ++la)
    {
        i = As.ind(la)/p;
        k = As.ind(la)%p;
        for (std::size_t l=0; l<Lb[k].size(); ++l)
        {
            lb = Lb[k][l];
            j  = Bs.ind(lb)%n;
            Cs(i,j) += a*As.val(la)*Bs.val(lb);
        }
    }
    check(Cs);
}

//==========================================================================
// [transpose]
/// Non-conjugate transpose.
///
/// transpose(A) is called for the syntax A^t.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    smatrix<> At = transpose(As);
///    disp(At);
/// \endcode
///
// \see reshape.
template<typename T>
smatrix<T> transpose(smatrix<T>const& As)
{
    std::size_t i, j, m=size(As,1), n=size(As,2);
    smatrix<T> Ast(n,m,nnz(As));
    for (std::size_t k=0; k<nnz(As); ++k)
    {
        i          = As.ind(k)/n;
        j          = As.ind(k)%n;
        Ast.ind(k) = j*m+i;
        Ast.val(k) = As.val(k);
    }
    return Ast;
}

//==========================================================================
// [values]
/// Copy non-zeros values in a matrix.
///
/// V = values(As) returns matrix with all non zeros values of sparse matrix As.
///
/// \code{.cpp}
///    smatrix<> As = {{1,2,3},{4,5,6}};
///    matrix<double> val = values(As);
///    disp(val);
/// \endcode
///
// \see index, find.
template<typename T>
matrix<T> values(smatrix<T>const& As)
{
    return matrix<T>(As.val());
}

}
