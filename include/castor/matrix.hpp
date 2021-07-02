/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : matrix.hpp                                    |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Dense matrix operations and algebra with      |
 |  `---'  |                matlab-like standard functions and procedures |
 +========================================================================+
 */

#pragma once

#define CASTOR_MATRIX_HPP

#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <math.h>
#include <string>
#include <tuple>
#include <typeinfo>
#include <vector>

/// \file matrix.hpp

namespace castor 
{

//==========================================================================//
//                            ENVIRONMENT DATA                              //
//==========================================================================//
#define M_1I  std::complex<double>(0,1)
#define M_1If std::complex<float>(0,1)

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288  
#endif

#define M_EPS(T) std::numeric_limits<T>::epsilon()

static std::vector<std::string> documentationFiles =
{
    "./matrix.hpp"
};

static auto ticTimer = std::chrono::high_resolution_clock::now();

class logical
{
public:
    logical():m_val(){};
    logical(bool v):m_val(v){};
    operator bool() const {return m_val;};
    template<typename T>
    operator std::complex<T>() const {return m_val;}
    bool* operator&() {return &m_val;}
    const bool* operator&() const {return &m_val;}
private:
    bool m_val;
};

template<typename T>
class view;

template<typename T>
class cview;


//==========================================================================//
//                             MASTER CLASS                                 //
//==========================================================================//
// [matrix]
/// Array with element of type T (double by default).
///
/// Matrix values are stored in std::vector<T> object, with row-major layout.
/// C++ numbering [0,n( is used for linear indexing A(l) and bi-linear A(i,j).
/// Dimensions and indexing are of type 'std::size_t', and boolean are given
/// in 'logical' type (which is 'std::uint8_t'), instead of standard 'bool'.
///
/// \code{.cpp}
///    matrix<logical> s = true;
///    matrix<int>     X = {1, 2, 3, 4};
///    matrix<double>  A = {{1, 2, 3, 4},
///                         {5, 6, 7, 8}};
///    auto            B = A;
/// \endcode
///
// \see view, cview.
template<typename T=double>
class matrix
{
public:
    // CONSTRUCTORS
    matrix();                                                      // DEFAULT
    matrix(T v);                                                   // SCALAR
    matrix(std::initializer_list<T>const& v);                      // 1D INITIALIZATION LIST
    matrix(std::initializer_list<std::vector<T>>const& v);         // 2D INITIALIZATION LIST
    matrix(std::size_t m, std::size_t n, T v=0);                   // FILL
    template<typename S>
    matrix(std::vector<S>const& v);                                // STD::VECTOR
    template<typename S=double>
    matrix(std::size_t m, std::size_t n, std::vector<S>const& v);  // FILL STD::VECTOR
    template<typename S=double>                                    // FILL STD::VECTOR(STD::VECTOR)
    matrix(std::size_t m, std::size_t n, std::vector<std::vector<S>>const& v);
    template<typename S>
    matrix(matrix<S>const& A);                                     // CAST
    
    // OPERATORS
    template<typename S>
    explicit       operator S() const;                                 // CONVERSION TO PRIMITIVE
    matrix<T>&     operator+=(matrix<T>const& A);                      // A += B;
    matrix<T>&     operator-=(matrix<T>const& A);                      // A -= B;
    matrix<T>&     operator*=(matrix<T>const& A);                      // A *= B;
    matrix<T>&     operator/=(matrix<T>const& A);                      // A /= B;
    const T&       operator()(std::size_t l) const;                    // A(l)
    T&             operator()(std::size_t l);                          // A(l)
    const cview<T> operator()(matrix<std::size_t>const& L) const;      // A(L)
    view<T>        operator()(matrix<std::size_t>const& L);            // A(L)
    const T&       operator()(std::size_t i, std::size_t j) const;     // A(i,j)
    T&             operator()(std::size_t i, std::size_t j);           // A(i,j)
    const cview<T> operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J) const; // A(I,J)
    view<T>        operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J); // A(I,J)
    
    // TOOLS
    void                 clear();
    void                 reshape(std::size_t m, std::size_t n);
    void                 resize(std::size_t m, std::size_t n, T v=(T)NAN);
    std::size_t          size(int dim=0) const;
    std::vector<T>const& val() const;
    
private:
    // ATTRIBUTS
    std::size_t    m_row;  // Number of rows
    std::size_t    m_col;  // Number of columns
    std::vector<T> m_val;  // Values
};


//==========================================================================//
//                           COMMON HEADERS                                 //
//==========================================================================//
template<typename T>
matrix<std::size_t> col(matrix<T>const& A);

template<typename T>
void disp(T const& A, int info=1, std::ostream& flux=std::cout);

template<typename T>
void disp(T const* A, int info=1, std::ostream& flux=std::cout);

template<typename T>
void disp(matrix<T>const& A, int info=2, std::ostream& flux=std::cout, std::size_t m=3, std::size_t n=3);

template<typename T>
inline void disp(matrix<T>const& A, std::size_t m, std::size_t n) {disp(A,2,std::cout,m,n);};

inline void error(std::string file, int line, std::string function, std::string comment);

inline void help(std::string name="help", std::vector<std::string> filename=documentationFiles);

template<typename T>
std::size_t numel(matrix<T>const& A);

inline matrix<std::size_t> range(std::size_t j, std::size_t k);

template<typename T>
matrix<std::size_t> row(matrix<T>const& A);

template<typename T>
std::size_t size(matrix<T>const& A, int dim);

inline void tic();

inline double toc(bool disp=true);

inline void warning(std::string file, int line, std::string function, std::string comment);


//==========================================================================//
//                                VIEW CLASS                                //
//==========================================================================//
// [view] 
/// 
template<typename T>
class view
{
public:
    // CONSTRUCTORS
    view(matrix<T>& A, matrix<std::size_t>const& L):
    m_typ{1}, m_idx{L}, m_mat{A} {};
    view(matrix<T>& A, matrix<std::size_t>const& I, matrix<std::size_t>const& J):
    m_typ{2}, m_idx{I}, m_jdx{J}, m_mat{A} {};
    
    // OPERATORS
    matrix<T>& operator=(T v) const
    {
        if (m_typ==1) {set(m_mat,m_idx,v);}
        else if (m_typ==2) {set(m_mat,m_idx,m_jdx,v);}
        return m_mat;
    }
    matrix<T>& operator=(std::initializer_list<T>const& v)
    {
        if (m_typ==1) {set(m_mat,m_idx,matrix<T>(v));}
        else if (m_typ==2) {set(m_mat,m_idx,m_jdx,matrix<T>(v));}
        return m_mat;
    }
    matrix<T>& operator=(std::initializer_list<std::vector<T>>const& v)
    {
        if (m_typ==1) {set(m_mat,m_idx,matrix<T>(v));}
        else if (m_typ==2) {set(m_mat,m_idx,m_jdx,matrix<T>(v));}
        return m_mat;
    }
    template<typename S>
    matrix<T>& operator=(matrix<S>const& A) const
    {
        if (m_typ==1) {set(m_mat,m_idx,A);}
        else if (m_typ==2) {set(m_mat,m_idx,m_jdx,A);}
        return m_mat;
    };
    
    // TOOLS
    matrix<T> eval() const
    {
        matrix<T> A;
        if (m_typ==1) {A = get(m_mat,m_idx);}
        else if (m_typ==2) {A = get(m_mat,m_idx,m_jdx);}
        return A;
    };
    
    // ATTRIBUTS
private:
    int                 m_typ;
    matrix<std::size_t> m_idx;
    matrix<std::size_t> m_jdx;
    matrix<T>&          m_mat;
};
template<typename T>
inline std::ostream& operator<<(std::ostream& flux, view<T>const& A) {disp(A); return flux;}
template<typename T>
inline matrix<T> eval(view<T>const& Av) {return Av.eval();}
template<typename T>
inline void disp(view<T>const& Av)
{
    error(__FILE__, __LINE__, __FUNCTION__,"Use 'eval()' function to build the right hand side matrix from view. For example: \nA = eval(B(L)) \nA = eval(B(I,J))");
}

template<typename T>
class cview
{
public:
    cview(matrix<T>const& A, matrix<std::size_t>const& L):
    m_typ{1}, m_idx{L}, m_mat{A} {};
    cview(matrix<T>const& A, matrix<std::size_t>const& I, matrix<std::size_t>const& J):
    m_typ{2}, m_idx{I}, m_jdx{J}, m_mat{A} {};
    matrix<T> eval() const
    {
        matrix<T> A;
        if (m_typ==1) {A = get(m_mat,m_idx);}
        else if (m_typ==2) {A = get(m_mat,m_idx,m_jdx);}
        return A;
    };
private:
    int                 m_typ;
    matrix<std::size_t> m_idx;
    matrix<std::size_t> m_jdx;
    matrix<T>const&     m_mat;
};
template<typename T>
inline std::ostream& operator<<(std::ostream& flux, cview<T>const& A) {disp(A); return flux;}
template<typename T>
inline matrix<T> eval(cview<T>const& Av) {return Av.eval();}
template<typename T>
inline void disp(cview<T>const& Av)
{
    error(__FILE__, __LINE__, __FUNCTION__,"Use 'eval()' function to build the right hand side matrix from view. For example: \nA = eval(B(L)) \nA = eval(B(I,J))");
}

//=========================================================================//
//                             CONSTRUCTORS                                //
//=========================================================================//

/// @name Constructors

//=========================================================================//
// [matrix.constructors]
/// Default constructor. 
///
/// Builds an empty matrix.
/// \code{.cpp}
///    matrix<logical> A;
///    disp(A);
/// \endcode
template<typename T>
matrix<T>::matrix():
m_row{0}, m_col{0}, m_val{std::vector<T>()} {}

/// Builds a (1x1) matrix from fundamental type constant or variable.
/// \code{.cpp}
///    matrix<int> A(M_PI);
///    disp(A);
/// \endcode
template<typename T>
matrix<T>::matrix(T v):
m_row{1}, m_col{1}, m_val{std::vector<T>(1,v)} {}

/// Builds a row matrix from initializer list.
/// \code{.cpp}
///    matrix<float> A({0,1,2,M_PI});
///    disp(A);
/// \endcode
template<typename T>
matrix<T>::matrix(std::initializer_list<T>const& v):
m_row{1}, m_col{v.size()}, m_val{v} {}

/// Builds a matrix from nested initializer list.
/// \code{.cpp}
///    matrix<double> A({{0,1,2,M_PI},{4,5,6,7},{8,9,10,11}});
///    disp(A);
/// \endcode
template<typename T>
matrix<T>::matrix(std::initializer_list<std::vector<T>>const& v)
{
    std::vector<std::vector<T>> v2(v);
    m_row = v2.size(),
    m_col = v2[0].size();
    m_val.resize(m_row*m_col);
    for (std::size_t i=0; i<m_row; ++i)
    {
        for (std::size_t j=0; j<m_col; ++j) {m_val[i*m_col+j] = v2[i][j];}
    }
}

/// Builds a (nxm) matrix whose all coefficients are equal to value of v.
/// \code{.cpp}
///    matrix<std::complex<double>> A(2,3,std::complex<double>(1,M_PI));
///    disp(A);
/// \endcode
template<typename T>
matrix<T>::matrix(std::size_t m, std::size_t n, T v):
m_row{m}, m_col{n}, m_val{std::vector<T>(m*n,v)} {}

/// Builds a row matrix from vector of standard library.
/// \code{.cpp}
///    matrix<double> A(std::vector<int>({0, 1, 2, 3}));
///    disp(A);
/// \endcode
template<typename T>
template<typename S>
matrix<T>::matrix(std::vector<S>const& v):
m_row{1}, m_col{v.size()}
{
    m_val.resize(v.size());
    for (std::size_t l=0; l<v.size(); ++l) {m_val[l] = v[l];}
}

/// Builds a (nxm) matrix from vector of standard library or initializer list.
/// \code{.cpp}
///    matrix<> A(2,2,{0, 1, 2, 3});
///    disp(A);
/// \endcode
template<typename T>
template<typename S>
matrix<T>::matrix(std::size_t m, std::size_t n, std::vector<S>const& v):
m_row{m}, m_col{n}
{
    m_val.resize(m*n);
    for (std::size_t l=0; l<m*n; ++l) {m_val[l] = v[l];}
}

/// Builds a (nxm) matrix from nested vector standard library or nested initializer list.
/// \code{.cpp}
///    matrix<> A(2,2,{{0, 1},{2,3}});
///    disp(A);
/// \endcode
template<typename T>
template<typename S>
matrix<T>::matrix(std::size_t m, std::size_t n, std::vector<std::vector<S>>const& v)
{
    m_row = v.size();
    m_col = v[0].size();
    m_val.resize(m_row*m_col);
    for (std::size_t i=0; i<m_row; ++i)
    {
        for (std::size_t j=0; j<m_col; ++j) {m_val[i*m_col+j] = v[i][j];}
    }
}

/// Builds a matrix from another matrix.
/// \code{.cpp}
///    matrix<> A(2,2,{{0, 1},{2,3}});
///    matrix<> B(A);
///    disp(A);
/// \endcode
template<typename T>
template<typename S>
matrix<T>::matrix(matrix<S>const& A):
m_row{A.size(1)}, m_col{A.size(2)}
{
    m_val.resize(A.size());
    for (std::size_t l=0; l<A.size(); ++l) {m_val[l] = A(l);}
}


//==========================================================================//
//                         INTERNAL OPERATORS                               //
//==========================================================================//

/// @name Internal operators

/// Converts a (1x1) to fundamental type variable.
/// \code{.cpp}
///    matrix<int> A(3);
///    a = double(A);
///    disp(a);
/// \endcode
template<typename T>
template<typename S>
matrix<T>::operator S() const
{
    if (m_row*m_col!=1)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Can't convert implicitly matrix to value because does not have a single element.");
    }
    return m_val[0];
}

/// Performs addition assignment.\n
/// If A is fundamental type variable, add value of A to each element of the matrix
/// If A is a matrix, add each element of A to each element of the matrix, 
/// in this case matrix dimensions must agree.
/// \code{.cpp}
///    matrix<float> V({0,1,2,3}); 
///    V += 1;
///    disp(V);
///    V += {1,1,1,1};
///    disp(V);
/// \endcode
template<typename T>
matrix<T>& matrix<T>::operator+=(matrix<T>const& A)
{
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<m_row*m_col; ++l) {m_val[l] += A(0);}
    }
    else if (m_row==A.size(1) && m_col==A.size(2))
    {
        for (std::size_t l=0; l<m_row*m_col; ++l) {m_val[l] += A(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return (*this);
}

/// Performs substraction assignment.\n
/// If A is fundamental type variable, substract value of A to each element of the matrix.\n
/// If A is a matrix, substract each element of A to each element of the matrix, 
/// in this case matrix dimensions must agree.\n
/// \code{.cpp}
///    matrix<float> V({0,1,2,3}); 
///    V -= 1;
///    disp(V);
///    V -= matrix<int>(1,4,1);
///    disp(V);
/// \endcode
template<typename T>
matrix<T>& matrix<T>::operator-=(matrix<T>const& A)
{
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<m_row*m_col; ++l) {m_val[l] -= A(0);}
    }
    else if (m_row==A.size(1) && m_col==A.size(2))
    {
        for (std::size_t l=0; l<m_row*m_col; ++l) {m_val[l] -= A(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return (*this);
}

/// Performs multiplication assignment.\n
/// If A is fundamental type variable, multiply value of A to each element of the matrix.\n
/// If A is a matrix, multiply each element of A to each element of the matrix, 
/// in this case matrix dimensions must agree.\n
/// \code{.cpp}
///    matrix<float> V({0,1,2,3}); 
///    V *= 2;
///    disp(V);
///    V *= matrix<int>(1,4,2);
///    disp(V);
/// \endcode
template<typename T>
matrix<T>& matrix<T>::operator*=(matrix<T>const& A)
{
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<m_row*m_col; ++l) {m_val[l] *= A(0);}
    }
    else if (m_row==A.size(1) && m_col==A.size(2))
    {
        for (std::size_t l=0; l<m_row*m_col; ++l) {m_val[l] *= A(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return (*this);
}

/// Performs division assignment.\n
/// If A is fundamental type variable, divide value of A to each element of the matrix.\n
/// If A is a matrix, divide each element of A to each element of the matrix, 
/// in this case matrix dimensions must agree.\n
/// \code{.cpp}
///    matrix<float> V({0,1,2,3}); 
///    V /= 2;
///    disp(V);
///    V /= matrix<int>(1,4,2);
///    disp(V);
/// \endcode
template<typename T>
matrix<T>& matrix<T>::operator/=(matrix<T>const& A)
{
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<m_row*m_col; ++l) {m_val[l] /= A(0);}
    }
    else if (m_row==A.size(1) && m_col==A.size(2))
    {
        for (std::size_t l=0; l<m_row*m_col; ++l) {m_val[l] /= A(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return (*this);
}

/// Returns an element of the matrix by giving the linear indexing of element (const accessor).\n
/// Remark : matrix values are stored in std::vector<T> object, with row-major layout.
/// \code{.cpp}
///    const matrix<float> A = eye(3,4); 
///    // access to last element
///    float coef = A(11);
///    disp(coef);
/// \endcode
template<typename T>
const T& matrix<T>::operator()(std::size_t l) const
{
#ifdef DEBUG
    if (l>=m_row*m_col)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
    }
#endif
    return m_val[l];
}

/// Returns an element of the matrix by giving the linear indexing of element (non-const accessor).\n
/// Remark : matrix values are stored in std::vector<T> object, with row-major layout.
/// \code{.cpp}
///    matrix<float> A = eye(3,4); 
///    // access to last element
///    A(11) = 10;
///    disp(A(11));
/// \endcode
template<typename T>
T& matrix<T>::operator()(std::size_t l)
{
#ifdef DEBUG
    if (l>=m_row*m_col)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Index exceeds array bounds.");
    }
#endif
    return m_val[l];
}

/// Returns a reference to a submatrix corresponding to linear indexing L (const accessor).
/// \code{.cpp}
///    const matrix<float> A = eye(3,4); 
///    // access to element equal to 1
///    auto B = A({0,5,10});
///    disp(eval(B));
/// \endcode
template<typename T>
const cview<T> matrix<T>::operator()(matrix<std::size_t>const& L) const
{
    return cview<T>(*this,L);
}

/// Returns a reference to a submatrix corresponding to linear indexing L (non-const accessor).
/// \code{.cpp}
///    matrix<float> A = eye(3,4); 
///    // access to element equal to 1
///    auto B = A({0,5,10});
///    disp(eval(B));
/// \endcode
template<typename T>
view<T> matrix<T>::operator()(matrix<std::size_t>const& L)
{
    return view<T>(*this,L);
}

/// Returns an element of the matrix by giving the bi-linear indexing of element (const accessor).
/// \code{.cpp}
///    const matrix<int> A = eye(3,4); 
///    int coef = A(2,3);
///    disp(coef);
/// \endcode
template<typename T>
const T&  matrix<T>::operator()(std::size_t i, std::size_t j) const
{
#ifdef DEBUG
    if (i>=m_row || j>=m_col)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Indices exceed array bounds.");
    }
#endif
    return m_val[i*m_col+j];
}

/// Returns an element of the matrix by giving the bi-linear indexing of element (non-const accessor).
/// \code{.cpp}
///    matrix<int> A = eye(3,4); 
///    A(2,3) = 5;
///    disp(A(2,3));
/// \endcode
template<typename T>
T& matrix<T>::operator()(std::size_t i, std::size_t j)
{
#ifdef DEBUG
    if (i>=m_row || j>=m_col)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Indices exceed array bounds.");
    }
#endif
    return m_val[i*m_col+j];
}

/// Returns a reference to a submatrix corresponding to bilinear indexing 
/// I for rows and J for columns (const accessor).
/// \code{.cpp}
///    matrix<float> A = eye(3,4); 
///    auto B = A({0,2}, {0,2});
///    disp(eval(B));
/// \endcode
template<typename T>
const cview<T> matrix<T>::operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J) const
{
    return cview<T>(*this,I,J);
}

/// Return a reference to a submatrix corresponding to bilinear indexing 
/// I for rows and J for columns (non-const accessor).
/// \code{.cpp}
///    matrix<float> A = eye(3,4); 
///    auto B = A({0,2}, {0,2});
///    disp(eval(B));
/// \endcode
template<typename T>
view<T> matrix<T>::operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J)
{
    return view<T>(*this,I,J);
}


//==========================================================================//
//                           INTERNAL TOOLS                                 //
//==========================================================================//

/// @name Internal tools

//==========================================================================
// [matrix.clear]
/// Removes all values of matrix A and fix size to 0x0.\n
/// Clear should be used to free memory without deleting matrix object.
/// \code{.cpp}
///    matrix<> X = {{1,2,3},{4,5,6}};
///    disp(X);
///    X.clear;
///    disp(X);
/// \endcode
template<typename T>
void matrix<T>::clear()
{
    m_row = 0;
    m_col = 0;
    m_val.clear();
}

//==========================================================================
// [matrix.reshape]
/// Modifies the shape of the matrix without changing its values.
/// \code{.cpp}
///    matrix<> A = eye(3,2);
///    A.reshape(2,3);
///    disp(A);
/// \endcode
template<typename T>
void matrix<T>::reshape(std::size_t m, std::size_t n)
{
    if (m_row*m_col!=m*n)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"To reshape the number of elements must not change.");
    }
    m_row = m;
    m_col = n;
}

//==========================================================================
// [matrix.resize]
/// Resizes and replaces the values of the matrix, eventually filled with the
/// values v specified (by default v is Nan).
/// \code{.cpp}
///    matrix<> A = eye(3,2);
///    A.resize(2,3,0);
///    disp(A);
/// \endcode
template<typename T>
void matrix<T>::resize(std::size_t m, std::size_t n, T v)
{
    if (m_row==0 || m_col==0)  {m_val.resize(m*n,v);}
    else if (m_row==1 && m==1) {m_val.resize(n,v);}
    else if (m_col==1 && n==1) {m_val.resize(m,v);}
    else if (m_col==n)         {m_val.resize(m*n,v);}
    else
    {
        std::vector<T> val(m*n,v);
        for (std::size_t i=0; i<std::min(m_row,m); ++i)
        {
            for (std::size_t j=0; j<std::min(m_col,n); ++j)
            {
                val[i*n+j] = m_val[i*m_col+j];
            }
        }
        m_val = val;
    }
    m_row = m;
    m_col = n;
}

//==========================================================================
// [matrix.size]
/// If dim=0 (default value) returns the number elements of the matrix
/// else returns the lengths of the specified dimensions dim.
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    disp(A.size(0));
///    disp(A.size(1));
///    disp(A.size(2));
/// \endcode
template<typename T>
std::size_t matrix<T>::size(int dim) const
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
// [matrix.val]
/// Returns std::vector containing matrix values.
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    std::vector val = A.val();
///    disp(matrix<>(val));
/// \endcode
template<typename T>
std::vector<T>const& matrix<T>::val() const
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
/// Insert matrix A inside a stream flux respecting standard C++ conventions.
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    std::cout << A << std::endl;
/// \endcode
///
// \see disp.
template<typename T>
std::ostream& operator<<(std::ostream& flux, matrix<T>const& A)
{
    disp(A,0,flux,size(A,1),size(A,2));
    return flux;
}

//==========================================================================
// [!]
/// Logical not.
///
/// !A performs a logical not of input array A, and returns an array
/// containing elements set to either logical 1 (TRUE) or logical 0 (FALSE).
/// An element of the output array is set to 1 if A contains a zero value
/// element at that same array location.  Otherwise, that element is set to
/// 0.
///
/// \code{.cpp}
///    matrix<logical> A = {{0,0,0},{1,1,1}};
///    matrix<logical> B = !A;
///    disp(B);
/// \endcode
///
// \see operator&&(), operator ||().
template<typename T>
matrix<logical> operator!(matrix<T>const& A)
{
    matrix<logical> B(size(A,1),size(A,2));
    for (std::size_t l=0; l<numel(B); ++l) {B(l) = !A(l);}
    return B;
}

//==========================================================================
// [&&]
/// Logical and.
///
/// A&&B performs a logical and of arrays A and B and returns an array
/// containing elements set to either logical 1 (TRUE) or logical 0
/// (FALSE). An element of the output array is set to 1 if both input
/// arrays contain a non-zero element at that same array location.
/// Otherwise, that element is set to 0. A and B must have compatible
/// sizes. Two inputs have compatible sizes if, for every dimension, the
/// dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<logical> A = {{0,0,0},{1,1,1}};
///    matrix<logical> B = {{0,1,0},{1,0,1}};
///    matrix<logical> C = A&&B;
///    matrix<logical> D = A&&1;
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator||(), operator!().
template<typename R, typename S>
matrix<logical> operator&&(matrix<R>const& A, matrix<S>const& B)
{
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<logical> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)&&B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)&&B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)&&B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline matrix<logical> operator&&(R A, matrix<S>const& B) {return matrix<R>(A)&&B;}
template<typename R, typename S>
inline matrix<logical> operator&&(matrix<R>const& A, S B) {return A&&matrix<S>(B);}

//==========================================================================
// [||]
/// Logical or.
///
/// A||B performs a logical or of arrays A and B and returns an array
/// containing elements set to either logical 1 (TRUE) or logical 0
/// (FALSE). An element of the output array is set to 1 if either input
/// array contains a non-zero element at that same array location.
/// Otherwise, that element is set to 0. A and B must have compatible
/// sizes. Two inputs have compatible sizes if, for every dimension, the
/// dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<logical> A = {{0,0,0},{1,1,1}};
///    matrix<logical> B = {{0,1,0},{1,0,1}};
///    matrix<logical> C = A||B;
///    matrix<logical> D = A||1;
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator&&(), operator!().
template<typename R, typename S>
matrix<logical> operator||(matrix<R>const& A, matrix<S>const& B)
{
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<logical> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)||B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)||B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)||B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline matrix<logical> operator||(R A, matrix<S>const& B) {return matrix<R>(A)||B;}
template<typename R, typename S>
inline matrix<logical> operator||(matrix<R>const& A, S B) {return A||matrix<S>(B);}

//==========================================================================
// [==]
/// Equal.
///
/// A==B does element by element comparisons between A and B and returns
/// an array with elements set to logical 1 (TRUE) where the relation is
/// true and elements set to logical 0 (FALSE) where it is not. A and B
/// must have compatible sizes. In the simplest cases, they can be the same
/// size or one can be a scalar. Two inputs have compatible sizes if, for
/// every dimension, the dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<>        A = {{1,2,3},{4,5,6}};
///    matrix<>        B = {{0,2,0},{4,0,6}};
///    matrix<logical> C = (A==B);
///    matrix<logical> D = (A==1);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator!=(), operator<=(), operator>=().
template<typename R, typename S>
matrix<logical> operator==(matrix<R>const& A, matrix<S>const& B)
{
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<logical> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)==B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)==B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)==B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline matrix<logical> operator==(R A, matrix<S>const& B) {return matrix<R>(A)==B;}
template<typename R, typename S>
inline matrix<logical> operator==(matrix<R>const& A, S B) {return A==matrix<S>(B);}

//==========================================================================
// [!=]
/// Not equal.
///
/// A!=B does element by element comparisons between A and B and returns
/// an array with elements set to logical 1 (TRUE) where the relation is
/// true and elements set to logical 0 (FALSE) where it is not. A and B
/// must have compatible sizes. In the simplest cases, they can be the same
/// size or one can be a scalar. Two inputs have compatible sizes if, for
/// every dimension, the dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<>        A = {{1,2,3},{4,5,6}};
///    matrix<>        B = {{0,2,0},{4,0,6}};
///    matrix<logical> C = (A!=B);
///    matrix<logical> D = (A!=1);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator==(), operator<=(), operator>=().
template<typename R, typename S>
matrix<logical> operator!=(matrix<R>const& A, matrix<S>const& B)
{
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<logical> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)!=B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)!=B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)!=B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline matrix<logical> operator!=(R A, matrix<S>const& B) {return matrix<R>(A)!=B;}
template<typename R, typename S>
inline matrix<logical> operator!=(matrix<R>const& A, S B) {return A!=matrix<S>(B);}

//==========================================================================
// [<=]
/// Less than or equal.
///
/// A<=B does element by element comparisons between A and B and returns
/// an array with elements set to logical 1 (TRUE) where the relation is
/// true and elements set to logical 0 (FALSE) where it is not. A and B
/// must have compatible sizes. In the simplest cases, they can be the same
/// size or one can be a scalar. Two inputs have compatible sizes if, for
/// every dimension, the dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<>        A = {{1,2,3},{4,5,6}};
///    matrix<>        B = {{0,2,0},{4,0,6}};
///    matrix<logical> C = (A<=B);
///    matrix<logical> D = (A<=3);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator==(), operator>=(), operator<().
template<typename R, typename S>
matrix<logical> operator<=(matrix<R>const& A, matrix<S>const& B)
{
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<logical> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)<=B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)<=B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)<=B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline matrix<logical> operator<=(R A, matrix<S>const& B) {return matrix<R>(A)<=B;}
template<typename R, typename S>
inline matrix<logical> operator<=(matrix<R>const& A, S B) {return A<=matrix<S>(B);}

//==========================================================================
// [<]
/// Less than.
///
/// A<B does element by element comparisons between A and B and returns
/// an array with elements set to logical 1 (TRUE) where the relation is
/// true and elements set to logical 0 (FALSE) where it is not. A and B
/// must have compatible sizes. In the simplest cases, they can be the same
/// size or one can be a scalar. Two inputs have compatible sizes if, for
/// every dimension, the dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<>        A = {{1,2,3},{4,5,6}};
///    matrix<>        B = {{0,2,0},{4,0,6}};
///    matrix<logical> C = (A<B);
///    matrix<logical> D = (A<3);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator<=(), operator>(), operator>=().
template<typename R, typename S>
matrix<logical> operator<(matrix<R>const& A, matrix<S>const& B)
{
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<logical> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)<B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)<B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)<B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline matrix<logical> operator<(R A, matrix<S>const& B) {return matrix<R>(A)<B;}
template<typename R, typename S>
inline matrix<logical> operator<(matrix<R>const& A, S B) {return A<matrix<S>(B);}

//==========================================================================
// [>=]
/// Greater than or equal.
///
/// A>=B does element by element comparisons between A and B and returns
/// an array with elements set to logical 1 (TRUE) where the relation is
/// true and elements set to logical 0 (FALSE) where it is not. A and B
/// must have compatible sizes. In the simplest cases, they can be the same
/// size or one can be a scalar. Two inputs have compatible sizes if, for
/// every dimension, the dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<>        A = {{1,2,3},{4,5,6}};
///    matrix<>        B = {{0,2,0},{4,0,6}};
///    matrix<logical> C = (A>=B);
///    matrix<logical> D = (A>=3);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator>(), operator<=(), operator<().
template<typename R, typename S>
matrix<logical> operator>=(matrix<R>const& A, matrix<S>const& B)
{
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<logical> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)>=B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)>=B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)>=B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline matrix<logical> operator>=(R A, matrix<S>const& B) {return matrix<R>(A)>=B;}
template<typename R, typename S>
inline matrix<logical> operator>=(matrix<R>const& A, S B) {return A>=matrix<S>(B);}

//==========================================================================
// [>]
/// Greater.
///
/// A>B does element by element comparisons between A and B and returns
/// an array with elements set to logical 1 (TRUE) where the relation is
/// true and elements set to logical 0 (FALSE) where it is not. A and B
/// must have compatible sizes. In the simplest cases, they can be the same
/// size or one can be a scalar. Two inputs have compatible sizes if, for
/// every dimension, the dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<>        A = {{1,2,3},{4,5,6}};
///    matrix<>        B = {{0,2,0},{4,0,6}};
///    matrix<logical> C = (A>B);
///    matrix<logical> D = (A>3);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator>=(), operator<(), operator<=().
template<typename R, typename S>
matrix<logical> operator>(matrix<R>const& A, matrix<S>const& B)
{
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<logical> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)>B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)>B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)>B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline matrix<logical> operator>(R A, matrix<S>const& B) {return matrix<R>(A)>B;}
template<typename R, typename S>
inline matrix<logical> operator>(matrix<R>const& A, S B) {return A>matrix<S>(B);}

//==========================================================================
// [+]
/// Plus.
///
/// A+B adds matrices A and B. A and B must have compatible sizes. In the
/// simplest cases, they can be the same size or one can be a scalar. Two
/// inputs have compatible sizes if, for every dimension, the dimension
/// sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = {{6,5,4},{3,2,1}};
///    matrix<> C = A+B;
///    matrix<> D = A+1;
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator-(), operator*(), operator/().
template<typename R, typename S>
auto operator+(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)+B(0));
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<T> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)+B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)+B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)+B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline auto operator+(R A, matrix<S>const& B) {return matrix<R>(A)+B;}
template<typename R, typename S>
inline auto operator+(matrix<R>const& A, S B) {return A+matrix<S>(B);}

//==========================================================================
// [-]
/// Unary or binary minus.
///
/// -A negates the elements of A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = -A;
///    disp(B);
/// \endcode
///
/// A-B subtracts matrix B from A. A and B must have compatible sizes. In
/// the simplest cases, they can be the same size or one can be a scalar.
/// Two inputs have compatible sizes if, for every dimension, the dimension
/// sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = {{6,5,4},{3,2,1}};
///    matrix<> C = A-B;
///    matrix<> D = A-1;
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator+(), operator*(), operator/().
template<typename S>
auto operator-(matrix<S>const& A)
{
    using T = decltype(-A(0));
    matrix<T> B(size(A,1),size(A,2));
    for (std::size_t l=0; l<numel(B); ++l) {B(l) = -A(l);}
    return B;
}
template<typename R, typename S>
auto operator-(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)-B(0));
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<T> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)-B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)-B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)-B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline auto operator-(R A, matrix<S>const& B) {return matrix<R>(A)-B;}
template<typename R, typename S>
inline auto operator-(matrix<R>const& A, S B) {return A-matrix<S>(B);}

//==========================================================================
// [*]
/// Array multiply.
///
/// A*B denotes element-by-element multiplication. A and B must have
/// compatible sizes. In the simplest cases, they can be the same size or
/// one can be a scalar. Two inputs have compatible sizes if, for every
/// dimension, the dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = {{6,5,4},{3,2,1}};
///    matrix<> C = A*B;
///    matrix<> D = A*2;
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator+(), operator-(), operator/().
template<typename R, typename S>
auto operator*(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)*B(0));
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<T> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)*B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)*B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)*B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline auto operator*(R A, matrix<S>const& B) {return matrix<R>(A)*B;}
template<typename R, typename S>
inline auto operator*(matrix<R>const& A, S B) {return A*matrix<S>(B);}

//==========================================================================
// [/]
/// Right array divide.
///
/// A/B denotes element-by-element division. A and B must have compatible
/// sizes. In the simplest cases, they can be the same size or one can be a
/// scalar. Two inputs have compatible sizes if, for every dimension, the
/// dimension sizes of the inputs are the same.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = {{6,5,4},{3,2,1}};
///    matrix<> C = A/B;
///    matrix<> D = A/2;
///    disp(C);
///    disp(D);
/// \endcode
///
// \see operator+(), operator-(), operator*().
template<typename R, typename S>
auto operator/(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)/B(0));
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<T> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(0)/B(l);}
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)/B(0);}
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l) {C(l) = A(l)/B(l);}
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline auto operator/(R A, matrix<S>const& B) {return matrix<R>(A)/B;}
template<typename R, typename S>
inline auto operator/(matrix<R>const& A, S B) {return A/matrix<S>(B);}


//==========================================================================//
//                           EXTERNAL TOOLS                                 //
//==========================================================================//
//==========================================================================
// [abs]
/// Absolute value.
///
/// abs(X) is the absolute value of the elements of X. When X is complex,
/// abs(X) is the complex modulus (magnitude) of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {{-1,-2,-3},{4,5,6}};
///    matrix<> Y = abs(X);
///    disp(Y);
/// \endcode
///
// \see angle, sign.
template<typename S>
auto abs(matrix<S>const& X)
{
    using T = decltype(std::abs(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::abs(X(l));}
    return Y;
}

//==========================================================================
// [acos]
/// Inverse cosine, result in radians.
///
/// acos(X) is the arccosine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1};
///    matrix<> Y = acos(X);
///    disp(Y);
/// \endcode
///
// \see acosd, cos.
template<typename S>
auto acos(matrix<S>const& X)
{
    using T = decltype(std::acos(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::acos(X(l));}
    return Y;
}

//==========================================================================
// [acosd]
/// Inverse cosine, result in degrees.
///
/// acosd(X) is the arccosine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1};
///    matrix<> Y = acosd(X);
///    disp(Y);
/// \endcode
///
// \see acos, cosd.
template<typename S>
auto acosd(matrix<S>const& X)
{
    using T = decltype(std::acos(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    T d = 180/M_PI;
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = d*std::acos(X(l));}
    return Y;
}

//==========================================================================
// [acosh]
/// Inverse hyperbolic cosine.
///
/// acosh(X) is the inverse hyperbolic cosine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {{1,2,3},{4,5,6}};
///    matrix<> Y = acosh(X);
///    disp(Y);
/// \endcode
///
// \see cosh.
template<typename S>
auto acosh(matrix<S>const& X)
{
    using T = decltype(std::acosh(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::acosh(X(l));}
    return Y;
}

//==========================================================================
// [all]
/// Linear indexing of all elements.
///
/// all(A) return row matrix with all indices of A in linear indexing.
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    matrix<std::size_t> L = all(A);
///    matrix<> B = eval(A(L));
///    disp(B);
/// \endcode
///
// \see row, col, get, set, view, cview.
template<typename T>
matrix<std::size_t> all(matrix<T>const& A)
{
    return range(0,numel(A));
}

//==========================================================================
// [angle]
/// Phase angle.
///
/// angle(X) returns the phase angles, in radians, of a matrix with
/// complex elements.
///
/// \code{.cpp}
///    matrix<> A = {{1,0,-1},{1,0,-1}};
///    matrix<> B = {{0,0,0},{1,1,1}};
///    auto     X = A + M_1I*B;
///    matrix<> Y = angle(X);
///    disp(Y);
/// \endcode
///
// \see abs.
template<typename S>
auto angle(matrix<S>const& X)
{
    using T = decltype(std::arg(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::arg(X(l));}
    return Y;
}

//==========================================================================
// [argintersect]
/// Index of intersection.
///
/// [IA,IB] = intersect(A,B) returns index vectors IA and IB such that
/// C = A(IA) and C = B(IB). If there are repeated common values in A or B
/// then the index of the first occurrence of each repeated value is returned.
///
/// \code{.cpp}
///    matrix<> A = {1,2,3};
///    matrix<> B = {2,3,4};
///    matrix<std::size_t> IA, IB;
///    std::tie(IA,IB) = argintersect(A,B);
///    disp(eval(A(IA)));
///    disp(eval(B(IB)));
/// \endcode
///
// \see intersect, argsetdiff, argunique.
template<typename R, typename S>
auto argintersect(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)+B(0));
    matrix<T> C = intersect(A,B);
    matrix<std::size_t> Ia(1,numel(C)), Ib(1,numel(C));
    if (numel(C)>0)
    {
        matrix<std::size_t> Is = argsort(A);
        std::size_t k=0;
        for (std::size_t l=0; l<numel(A); ++l)
        {
            if (A(Is(l))==C(k))
            {
                Ia(k)=Is(l);
                ++k;
                if (k==numel(C)) {break;}
            }
        }
        Is = argsort(B);
        k  = 0;
        for (std::size_t l=0; l<numel(B); ++l)
        {
            if (B(Is(l))==C(k))
            {
                Ib(k) = Is(l);
                ++k;
                if (k==numel(C)) {break;}
            }
        }
    }
    return std::make_tuple(Ia,Ib);
}

//==========================================================================
// [argmax]
/// Index of max.
///
/// I = argmax(A) returns the index corresponding to the maximum value.
/// 
/// \code{.cpp}
///    matrix<>    A = {{1,2,3},{4,5,6}};
///    std::size_t I = argmax(A);
///    disp(I);
/// \endcode
///
/// I = argmax(A,DIM) returns a vector of indices corresponding to the
/// maximum values into operating dimension.
///
/// \code{.cpp}
///    matrix<>            A = {{1,2,3},{4,5,6}};
///    matrix<std::size_t> I = argmax(A,1);
///    matrix<std::size_t> J = argmax(A,2);
///    disp(I);
///    disp(J);
/// \endcode
///
// \see max, maximum, argmin, argsort.
template<typename T>
matrix<std::size_t> argmax(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<std::size_t> B;
    if (dim==0)
    {
        B.resize(1,1,0);
        for (std::size_t l=1; l<numel(A); ++l) {if (A(l)>A(B(0))) {B(0)=l;}}
    }
    else if (dim==1)
    {
        B.resize(1,size(A,2));
        for (std::size_t j=0; j<size(A,2); ++j)
        {
            B(j) = 0;
            for (std::size_t i=1; i<size(A,1); ++i) {if (A(i,j)>A(B(j),j)) {B(j)=i;}}
        }
    }
    else if (dim==2)
    {
        B.resize(size(A,1),1);
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            B(i) = 0;
            for (std::size_t j=1; j<size(A,2); ++j) {if (A(i,j)>A(i,B(i))) {B(i)=j;}}
        }
    }
    return B;
}
template<typename T>
inline std::size_t argmax(matrix<T>const& A) {return (std::size_t)argmax(A,0);}

//==========================================================================
// [argmin]
/// Index of min.
///
/// I = argmin(A) returns the index corresponding to the minimum value.
///
/// \code{.cpp}
///    matrix<>    A = {{1,2,3},{4,5,6}};
///    std::size_t I = argmin(A);
///    disp(I);
/// \endcode
///
/// I = argmin(A,DIM) returns a vector of indices corresponding to the
/// minimum values into operating dimension.
///
/// \code{.cpp}
///    matrix<>            A = {{1,2,3},{4,5,6}};
///    matrix<std::size_t> I = argmin(A,1);
///    matrix<std::size_t> J = argmin(A,2);
///    disp(I);
///    disp(J);
/// \endcode
///
// \see min, minimum, argmax, argsort.
template<typename T>
matrix<std::size_t> argmin(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<std::size_t> B;
    if (dim==0)
    {
        B.resize(1,1,0);
        for (std::size_t l=1; l<numel(A); ++l) {if (A(l)<A(B(0))) {B(0)=l;}}
    }
    else if (dim==1)
    {
        B.resize(1,size(A,2));
        for (std::size_t j=0; j<size(A,2); ++j)
        {
            B(j) = 0;
            for (std::size_t i=1; i<size(A,1); ++i) {if (A(i,j)<A(B(j),j)) {B(j)=i;}}
        }
    }
    else if (dim==2)
    {
        B.resize(size(A,1),1);
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            B(i) = 0;
            for (std::size_t j=1; j<size(A,2); ++j) {if (A(i,j)<A(i,B(i))) {B(i)=j;}}
        }
    }
    return B;
}
template<typename T>
inline std::size_t argmin(matrix<T>const& A) {return (std::size_t)argmin(A,0);}

//==========================================================================
// [argsetdiff]
/// Index of difference.
///
/// IA = argsetdiff(A,B) returns index vectors IA such that C = A(IA).
/// If there are repeated values in A that are not in B, then the index
/// of the first occurrence of each repeated value is returned
///
/// \code{.cpp}
///    matrix<> A = {1,2,3,0};
///    matrix<> B = {2,3,4};
///    matrix<std::size_t> IA, IB;
///    IA = argsetdiff(A,B);
///    disp(eval(A(IA)));
/// \endcode
///
// \see setdiff, argintersect, argunique.
template<typename R, typename S>
matrix<std::size_t> argsetdiff(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)+B(0));
    matrix<T> C = setdiff(A,B);
    matrix<std::size_t> Ia(1,numel(C));
    if (numel(C)>0)
    {
        matrix<std::size_t> Is = argsort(A);
        std::size_t k=0;
        for (std::size_t l=0; l<numel(A); ++l)
        {
            if (A(Is(l))==C(k))
            {
                Ia(k)=Is(l);
                ++k;
                if (k==numel(C)) {break;}
            }
        }
    }
    return Ia;
}

//==========================================================================
// [argsort]
/// Index of sort.
///
/// I = argsort(A) returns a sort index I which specifies how the elements
/// of A were rearranged to obtain the sorted output B = sort(A).
///
/// \code{.cpp}
///    matrix<>            A = {{6,5,4},{3,2,1}};
///    matrix<std::size_t> I = argsort(A);
///    disp(I);
/// \endcode
///
/// I = argsort(A,DIM) returns a sort index I which specifies how
/// the elements of A were rearranged to obtain the sorted output
/// B = sort(A,DIM) into operating dimension.
///
/// \code{.cpp}
///    matrix<>            A = {{6,5,4},{3,2,1}};
///    matrix<std::size_t> I = argsort(A,1);
///    matrix<std::size_t> J = argsort(A,2);
///    disp(I);
///    disp(J);
/// \endcode
///
/// The sort ordering is stable. Namely, when more than one element has the
/// same value, the order of the equal elements is preserved in the sorted
/// output B and the indices I relating to equal elements are ascending.
///
// \see sort, argmin, argmax.
template<typename T>
matrix<std::size_t> argsort(matrix<T>const& A, int dim=0)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    std::size_t m=size(A,1), n=size(A,2);
    matrix<std::size_t> B(m,n);
    std::vector<std::size_t> v;
    if (dim==0)
    {
        v.resize(m*n);
        for (std::size_t l=0; l<v.size(); ++l) {v[l] = l;}
        std::stable_sort(v.begin(), v.end(), [&A](std::size_t p, std::size_t q) {return A(p) < A(q);});
        for (std::size_t l=0; l<v.size(); ++l) {B(l) = v[l];}
    }
    else if (dim==1)
    {
        v.resize(m);
        std::vector<T> w(m);
        for (std::size_t j=0; j<n; ++j)
        {
            for (std::size_t i=0; i<m; ++i)
            {
                v[i] = i;
                w[i] = A(i,j);
            }
            std::stable_sort(v.begin(), v.end(), [&w](std::size_t p, std::size_t q) {return w[p] < w[q];});
            for (std::size_t i=0; i<m; ++i) {B(i,j) = v[i];}
        }
    }
    else if (dim==2)
    {
        v.resize(n);
        std::vector<T> w(n);
        for (std::size_t i=0; i<m; ++i)
        {
            for (std::size_t j=0; j<n; ++j)
            {
                v[j] = j;
                w[j] = A(i,j);
            }
            std::stable_sort(v.begin(), v.end(), [&w](std::size_t p, std::size_t q) {return w[p] < w[q];});
            for (std::size_t j=0; j<n; ++j) {B(i,j)=v[j];}
        }
    }
    return B;
}

//==========================================================================
// [argunique]
/// Index of unique.
/// 
/// [IA,IB] = unique(A) returns index vectors IA and IB such that
/// B = A(IA) and A = B(IB).
///
/// \code{.cpp}
///    matrix<> A = {1,2,1,3,2,3};
///    matrix<> B = unique(A);
///    matrix<std::size_t> IA, IB;
///    std::tie(IA,IB) = argunique(A);
///    disp(eval(A(IA)));
///    disp(eval(B(IB)));
/// \endcode
///
// \see unique, argsetdiff, argintersect.
template<typename T>
auto argunique(matrix<T>const& A)
{
    matrix<T> B = unique(A);
    matrix<std::size_t> Is = argsort(A);
    matrix<std::size_t> Ia(1,numel(B));
    matrix<std::size_t> Ib(size(A,1),size(A,2));
    std::size_t k=0;
    Ia(0) = Is(0);
    for (std::size_t l=0; l<numel(A); ++l)
    {
        if (A(Is(l))>B(k)) {++k; Ia(k)=Is(l);}
        else if (A(Is(l))<B(k))
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
        }
        Ib(Is(l)) = k;
    }
    return std::make_tuple(Ia,Ib);
}

//==========================================================================
// [asin]
/// Inverse sine, result in radians.
///
/// asin(X) is the arcsine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1};
///    matrix<> Y = asin(X);
///    disp(Y);
/// \endcode
///
// \see sin, asind.
template<typename S>
auto asin(matrix<S>const& X)
{
    using T = decltype(std::asin(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::asin(X(l));}
    return Y;
}

//==========================================================================
// [asind]
/// Inverse sine, result in degrees.
///
/// asind(X) is the arcsine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1};
///    matrix<> Y = asind(X);
///    disp(Y);
/// \endcode
///
// \see sind, asin.
template<typename S>
auto asind(matrix<S>const& X)
{
    using T = decltype(std::asin(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    T d = 180/M_PI;
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = d*std::asin(X(l));}
    return Y;
}

//==========================================================================
// [asinh]
/// Inverse hyperbolic sine.
///
/// asinh(X) is the inverse hyperbolic sine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {{1,2,3},{4,5,6}};
///    matrix<> Y = asinh(X);
///    disp(Y);
/// \endcode
///
// \see sinh.
template<typename S>
auto asinh(matrix<S>const& X)
{
    using T = decltype(std::asinh(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::asinh(X(l));}
    return Y;
}

//==========================================================================
// [atan]
/// Inverse tangent, result in radians.
///
/// atan(X) is the arctangent of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1};
///    matrix<> Y = atan(X);
///    disp(Y);
/// \endcode
///
// \see tan, atand.
template<typename S>
auto atan(matrix<S>const& X)
{
    using T = decltype(std::atan(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::atan(X(l));}
    return Y;
}

//==========================================================================
// [atand]
/// Inverse tangent, result in degrees.
///
/// atand(X) is the arctangent of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1};
///    matrix<> Y = atand(X);
///    disp(Y);
/// \endcode
///
// \see tand, atan.
template<typename S>
auto atand(matrix<S>const& X)
{
    using T = decltype(std::asin(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    T d = 180/M_PI;
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = d*std::atan(X(l));}
    return Y;
}

//==========================================================================
// [atanh]
/// Inverse hyperbolic tangent.
///
/// atanh(X) is the inverse hyperbolic tangent of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {{1,2,3},{4,5,6}};
///    matrix<> Y = atanh(X);
///    disp(Y);
/// \endcode
///
// \see tanh.
template<typename S>
auto atanh(matrix<S>const& X)
{
    using T = decltype(std::atanh(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::atanh(X(l));}
    return Y;
}

//==========================================================================
// [cart2pol]
/// Transform Cartesian to polar coordinates.
///
/// (THE,RHO) = cart2pol(X,Y) transforms corresponding elements of data
/// stored in cartesian coordinates (X,Y) to polar coordinates (THE,RHO).
///
/// The arrays X and Y must be the same size and angle THE is returned
/// in radians.
///
/// <b>Remark:</b> This function supports scalar arguments.
///
/// \code{.cpp}
///    matrix<> X = {1,-1,-1,1};
///    matrix<> Y = {1,1,-1,-1};
///    matrix<> THE, RHO;
///    std::tie(THE,RHO) = cart2pol(X,Y);
///    disp(rad2deg(THE));
///    disp(RHO);
/// \endcode
///
// \see pol2cart, cart2sph, sph2cart.
template<typename R, typename S>
auto cart2pol(matrix<R>const& X, matrix<S>const& Y)
{
    if (size(X,1)!=size(Y,1) || size(X,2)!=size(Y,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(X(0)+Y(0));
    matrix<T> THE(size(X,1),size(X,2)), RHO(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(X); ++l)
    {
        RHO(l) = std::sqrt(std::pow(X(l),2)+std::pow(Y(l),2));
        THE(l) = std::atan2(Y(l),X(l));
    }
    return std::make_tuple(THE,RHO);
}
template<typename R, typename S>
inline auto cart2pol(R const &X, S const &Y)
{
    using T = decltype(X + Y);
    auto rho     = std::sqrt(X*X + Y*Y);
    auto theta   = std::atan2(Y,X);
    return std::make_tuple(theta,rho);
}

//==========================================================================
// [cart2sph]
/// Transform Cartesian to spherical coordinates.
///
/// (THE,PHI,RHO) = cart2sph(X,Y,Z) transforms corresponding elements
/// of data stored in Cartesian coordinates (X,Y,Z) to spherical
/// coordinates (THE,PHI,RHO). THE is the counterclockwise angle in
/// the xy plane measured from the positive x axis and PHI is the
/// elevation angle from the xy plane.
///
/// The arrays X,Y and Z must be the same size and angles (THE,PHI) are
/// returned in radians.
///
/// <b>Remark:</b> This function supports scalar arguments.
///
/// \code{.cpp}
///    matrix<> X = {1,-1,-1,1};
///    matrix<> Y = {1,1,-1,-1};
///    matrix<> Z = {1,1,-1,-1};
///    matrix<> THE, PHI, RHO;
///    std::tie(THE,PHI,RHO) = cart2sph(X,Y,Z);
///    disp(rad2deg(THE));
///    disp(rad2deg(PHI));
///    disp(RHO);
/// \endcode
///
// \see sph2cart, cart2pol, pol2cart.
template<typename Q, typename R, typename S>
auto cart2sph(matrix<Q>const& X, matrix<R>const& Y, matrix<S>const& Z)
{
    if (size(X,1)!=size(Y,1) || size(X,2)!=size(Y,2) ||
        size(X,1)!=size(Z,1) || size(X,2)!=size(Z,2) )
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(X(0)+Y(0)+Z(0));
    matrix<T> THE(size(X,1),size(X,2)), PHI(size(X,1),size(X,2)), RHO(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(X); ++l)
    {
        RHO(l) = std::sqrt(std::pow(X(l),2)+std::pow(Y(l),2)+std::pow(Z(l),2));
        PHI(l) = std::asin(Z(l)/RHO(l));
        THE(l) = std::atan2(Y(l),X(l));
    }
    return std::make_tuple(THE,PHI,RHO);
}
template<typename Q, typename R, typename S>
inline auto cart2sph(Q const &X, R const &Y, S const &Z)
{
    using T = decltype(X+Y+Z);
    T rho = std::sqrt(X*X + Y*Y + Z*Z);
    T phi = std::asin(Z/rho);
    T the = std::atan2(Y,X);
    return std::make_tuple(the,phi,rho);
}

//==========================================================================
// [cast]
/// Cast an array to a different data type.
///
/// B = cast(A,T()) or B = cast<T>(A) casts A to newtype T. A must be
/// convertible to type T.
///
/// Classical type T are logical, int, float, double (default), etc.
///
/// \code{.cpp}
///    matrix<int>   A = {1,2,3,4};
///    matrix<long>  B = cast(A,long());
///    matrix<float> C = cast<float>(A);
///    auto          D = cast(A);
///    disp(A);
///    disp(B);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see matrix.
template<typename T=double, typename S>
inline matrix<T> cast(matrix<S>const& A, T b=0) {return matrix<T>(A);}

//==========================================================================
// [cat]
/// Concatenate array.
///
/// cat(DIM,A,B) concatenates the arrays A and B along the dimension DIM.
/// cat(1,A,B) is the same as vertcat(A,B).
/// cat(2,A,B) is the same as horzcat(A,B).
///
/// \code{.cpp}
///    matrix<> A = {1,2,3};
///    matrix<> B = {4,5,6};
///    matrix<> C = cat(1,A,B);
///    matrix<> D = cat(2,A,B);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see vertcat, horzcat.
template<typename R, typename S>
auto cat(int dim, matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)+B(0));
    matrix<T> C;
    if (dim==1 && size(A,2)==size(B,2))
    {
        std::size_t m1=size(A,1), m2=size(B,1), n=size(A,2);
        C.resize(m1+m2,n);
        for (std::size_t i=0; i<m1; ++i)
        {
            for (std::size_t j=0; j<n; ++j) {C(i,j) = A(i,j);}
        }
        for (std::size_t i=0; i<m2; ++i)
        {
            for (std::size_t j=0; j<n; ++j) {C(m1+i,j) = B(i,j);}
        }
    }
    else if (dim==2 && size(A,1)==size(B,1))
    {
        std::size_t m=size(A,1), n1=size(A,2), n2=size(B,2);
        C.resize(m,n1+n2);
        for (std::size_t i=0; i<m; ++i)
        {
            for (std::size_t j=0; j<n1; ++j) {C(i,j) = A(i,j);}
            for (std::size_t j=0; j<n2; ++j) {C(i,n1+j) = B(i,j);}
        }
    }
    else if (dim<1 || dim>2)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 1 for rows or 2 for columns.");
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions of arrays being concatenated are not consistent.");
    }
    return C;
}
template<typename R, typename S>
inline auto cat(int dim, R A, matrix<S>const& B) {return cat(dim,matrix<R>(A),B);}
template<typename R, typename S>
inline auto cat(int dim, matrix<R>const& A, S B) {return cat(dim,A,matrix<S>(B));}
template<typename T=double, typename R, typename S>
inline auto cat(int dim, R const& A, S const& B) {return cat(dim,matrix<T>(A),matrix<T>(B));}

//==========================================================================
// [ceil]
/// Round towards plus infinity.
///
/// ceil(X) rounds the elements of X to the nearest integers
/// towards infinity.
///
/// \code{.cpp}
///    matrix<> X = {{1,2,3},{4,5,6}};
///    matrix<> Y = ceil(X+0.5);
///    disp(Y);
/// \endcode
///
// \see floor, round.
template<typename S>
auto ceil(matrix<S>const& X)
{
    using T = decltype(std::ceil(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::ceil(X(l));}
    return Y;
}

//==========================================================================
// [clear]
/// Clear matrix to free memory.
///
/// clear(A) removes all values of matrix A and fix size to 0x0.
/// clear should be used to free memory without deleting matrix object.
///
/// \code{.cpp}
///    matrix<> X = {{1,2,3},{4,5,6}};
///    disp(X);
///    clear(X);
///    disp(X);
/// \endcode
///
// \see matrix, resize.
template<typename T>
void clear(matrix<T>& A)
{
    A.clear();
}

//==========================================================================
// [col]
/// Bi-linear indexing of column elements.
///
/// col(A) return row matrix with column indices of A in bi-linear indexing.
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    matrix<std::size_t> J = col(A);
///    matrix<> B = eval(A(0,J));
///    disp(B);
/// \endcode
///
// \see row, all, get, set, view, cview.
template<typename T>
matrix<std::size_t> col(matrix<T>const& A)
{
    return range(0,size(A,2));
}

//==========================================================================
// [colon]
/// Colon operator (:).
///
/// colon(J,K) is the same as [J, J+1, ..., J+m], where J+m <= K. In the
/// case where both J and K are integers, this is simply [J, J+1, ..., K].
/// This syntax returns an empty matrix if J > K.
///
/// \code{.cpp}
///    matrix<int> I = colon(1,4);
///    disp(I);
///    matrix<> X = colon(0.5,4);
///    disp(X);
/// \endcode
///
/// colon(J,I,K) is the same as [J, J+I, ..., J+m*I], where J+m*I <= K.
/// This syntax returns an empty matrix when I == 0, I > 0 and J > K, or
/// I < 0 and J < K.
///
/// \code{.cpp}
///    matrix<int> I = colon(1,2,9);
///    disp(I);
///    matrix<> X = colon(0.5,0.5,2.1);
///    disp(X);
/// \endcode
///
// \see linspace.
template<typename U, typename V, typename W>
auto colon(U j, V i, W k)
{
    using T = decltype(j+i+k);
    matrix<T> A;
    if ( (i==0) || ((j>k) && (i>0)) || ((j<k) && (i<0)) ) {}
    else
    {
        std::size_t n = std::floor((k-j)/i)+1;
        A.resize(1,n);
        for (std::size_t l=0; l<n; ++l) {A(l) = j + (T)l*i;};
    }
    return A;
}
template<typename U, typename V>
auto colon(U j, V k)
{
    return colon(j,1,k);
}

//==========================================================================
// [conj]
/// Complex conjugate.
///
/// conj(X) is the complex conjugate of X.
/// For a complex matrix X, conj(X) = real(X) - M_1i*imag(X).
///
/// \code{.cpp}
///    matrix<> A = {1,0,-1};
///    matrix<> B = {-1,0,1};
///    auto     X = A + M_1I*B;
///    auto     Y = conj(X);
///    disp(Y);
/// \endcode
///
// \see real, imag.
template<typename S>
auto conj(matrix<S>const& X)
{
    using T = decltype(std::conj(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::conj(X(l));}
    return Y;
}
inline matrix<float> conj(matrix<float>const& A) {return A;}
inline matrix<double> conj(matrix<double>const& A) {return A;}

//==========================================================================
// [conv]
/// Convolution product and polynomial multiplication.
///
/// C = conv(A,B) convolves array A and B. The resulting array is vector
/// of length NUMEL(A)+NUMEL(B)-1. If A and B are vectors of polynomial
/// coefficients, convolving them is equivalent to multiplying the two
/// polynomials.
///
/// C = conv(A,B,DIM) convolve along dimension DIM if A and B have
/// compatible size.
///
/// \code{.cpp}
///    matrix<float> A = eye(3,4);
///    matrix<float> B = rand(3,10);
///    disp(conv(A,B,2));
/// \endcode
///
// \see fftconv.
template<typename T>
matrix<T> conv(matrix<T>const& A, matrix<T>const& B, int dim=0)
{
    matrix<T> C;
    if (dim==0)
    {
        C = matrix<T>(1,numel(A)+numel(B)-1);
        for (std::size_t i=0; i<numel(A); ++i)
        {
            for (std::size_t j=0; j<numel(B); ++j)
            {
                C(i+j) += A(i)*B(j);
            }
        }
    }
    else if (dim==1)
    {
        if (size(A,2)!=size(B,2))
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
        }
        C = matrix<T>(size(A,1)+size(B,1)-1,size(A,2));
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            for (std::size_t j=0; j<size(B,1); ++j)
            {
                for (std::size_t k=0; k<size(A,2); ++k)
                {
                    
                    C(i+j,k) += A(i,k)*B(j,k);
                }
            }
        }
    }
    else if (dim==2)
    {
        if (size(A,1)!=size(B,1))
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
        }
        C = matrix<T>(size(A,1),size(A,2)+size(B,2)-1);
        for (std::size_t k=0; k<size(A,1); ++k)
        {
            for (std::size_t i=0; i<size(A,2); ++i)
            {
                for (std::size_t j=0; j<size(B,2); ++j)
                {
                    C(k,i+j) += A(k,i)*B(k,j);
                }
            }
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    return C;
}

//==========================================================================
// [cos]
/// Cosine of argument in radians.
///
/// cos(X) is the cosine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {0,M_PI/2,M_PI,-M_PI/2};
///    matrix<> Y = cos(X);
///    disp(Y);
/// \endcode
///
// \see acos, cosd.
template<typename S>
auto cos(matrix<S>const& X)
{
    using T = decltype(std::cos(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::cos(X(l));}
    return Y;
}

//==========================================================================
// [cosd]
/// Cosine of argument in degrees.
///
/// cosd(X) is the cosine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {0,90,180,-90};
///    matrix<> Y = cosd(X);
///    disp(Y);
/// \endcode
///
// \see acosd, cos.
template<typename S>
auto cosd(matrix<S>const& X)
{
    using T = decltype(std::cos(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    T r = M_PI/180;
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::cos(r*X(l));}
    return Y;
}

//==========================================================================
// [cosh]
/// Hyperbolic cosine.
///
/// cosh(X) is the hyperbolic cosine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {{1,2,3},{4,5,6}};
///    matrix<> Y = cosh(X);
///    disp(Y);
/// \endcode
///
// \see acosh.
template<typename S>
auto cosh(matrix<S>const& X)
{
    using T = decltype(std::cosh(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::cosh(X(l));}
    return Y;
}

//==========================================================================
// [cross]
/// Matrix cross product.
///
/// C = cross(A,B) returns the cross product of the matrix A and B, C=AxB.
/// A and B must have same number of row with 3 columns.
///
/// \code{.cpp}
///    matrix<> A = {{1,0,0},{0,1,0},{0,0,1}};
///    matrix<> B = {{0,1,0},{0,0,1},{1,0,0}};
///    matrix<> C = cross(A,B);
///    disp(C);
/// \endcode
///
// \see dot.
template<typename R, typename S>
auto cross(matrix<R>const& A, matrix<S>const& B)
{
    if(size(A,2)!=3 || size(B,2)!= 3)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Number of columns should be 3.");
    }
    if(size(A,1)!=size(B,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix rows must agree.");
    }
    using T = decltype(A(0)*B(0));
    matrix<T> C(size(A,1),3);
    for(std::size_t l=0; l<size(C,1); ++l)
    {
        C(l,0) = A(l,1)*B(l,2) - A(l,2)*B(l,1);
        C(l,1) = A(l,2)*B(l,0) - A(l,0)*B(l,2);
        C(l,2) = A(l,0)*B(l,1) - A(l,1)*B(l,0);
    }
    return C;
}

//==========================================================================
// [cumprod]
/// Cumulative product of elements.
///
/// B = cumprod(A) computes the cumulative product of all elements of A.
/// B is the same size as A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = cumprod(A);
///    disp(B);
/// \endcode
///
/// B = cumprod(A,DIM) computes the cumulative product along the dimension DIM.
/// if DIM==0, B is the cumulative product.
/// if DIM==1, B is a vector containing the cumulative product from each column.
/// if DIM==2, B is a vector containing the cumulative product from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = cumprod(A,1);
///    matrix<> C = cumprod(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see cumsum, prod.
template<typename T>
matrix<T> cumprod(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<T> B = A;
    if (dim==0)
    {
        for (std::size_t l=1; l<numel(A); ++l) {B(l) *= B(l-1);}
    }
    else if (dim==1)
    {
        for (std::size_t i=1; i<size(A,1); ++i)
        {
            for (std::size_t j=0; j<size(A,2); ++j) {B(i,j) *= B(i-1,j);}
        }
    }
    else if (dim==2)
    {
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            for (std::size_t j=1; j<size(A,2); ++j) {B(i,j) *= B(i,j-1);}
        }
    }
    return B;
}
template<typename T>
inline matrix<T> cumprod(matrix<T>const& A) {return cumprod(A,0);}

//==========================================================================
// [cumsum]
/// Cumulative sum of elements.
///
/// B = cumsum(A) computes the cumulative sum of all elements of A.
/// B is the same size as A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = cumsum(A);
///    disp(B);
/// \endcode
///
/// B = cumsum(A,DIM) computes the cumulative sum along the dimension DIM.
/// if DIM==0, B is the cumulative sum.
/// if DIM==1, B is a vector containing the cumulative sum from each column.
/// if DIM==2, B is a vector containing the cumulative sum from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = cumsum(A,1);
///    matrix<> C = cumsum(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see cumprod, sum.
template<typename T>
matrix<T> cumsum(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<T> B = A;
    if (dim==0)
    {
        for (std::size_t l=1; l<numel(A); ++l) {B(l) += B(l-1);}
    }
    else if (dim==1)
    {
        for (std::size_t i=1; i<size(A,1); ++i)
        {
            for (std::size_t j=0; j<size(A,2); ++j) {B(i,j) += B(i-1,j);}
        }
    }
    else if (dim==2)
    {
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            for (std::size_t j=1; j<size(A,2); ++j) {B(i,j) += B(i,j-1);}
        }
    }
    return B;
}
template<typename T>
inline matrix<T> cumsum(matrix<T>const& A) {return cumsum(A,0);}

//==========================================================================
// [deg2rad]
/// Convert angles from degrees to radians.
///
/// deg2rad(X) converts angle units from degrees to radians for each
/// element of X.
///
/// \code{.cpp}
///    matrix<> X = {0,90,180,-90};
///    matrix<> Y = deg2rad(X);
///    disp(Y);
/// \endcode
///
// \see rad2deg.
template<typename T>
matrix<T> deg2rad(matrix<T>const& X)
{
    matrix<T> Y = X;
    double c = M_PI/180;
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) *= c;}
    return Y;
}
template<typename T>
auto deg2rad(T x)
{
    return x*M_PI/180;
}

//==========================================================================
// [dft]
/// Discrete Fourier Transform of an array.
///
/// Y = dft(X) computes the discrete Fourier transform (DFT) of all the
/// N values of array X using a naive algorithm :
///
/// Y(k) = sum_{n=0}^(N-1) e^{-2i\pi*k * n/N} X(n).
///
/// Y = dft(X,n) returns the n-point DFT. If no value is specified, Y is the
/// dft of all values of X, in linear indexing.
///
/// Y = dft(X,n,dim) returns the Fourier transform along the dimension dim.
/// For example, if X is a matrix, then dft(X,n,2) returns the n-point Fourier
/// transform of each row.
///
/// Y = dft(X,n,dim,flag) specify is inverse fft is computed.
///
/// \code{.cpp}
///    matrix<float> X = eye(1,4);
///    disp(fft(X));
/// \endcode
///
// \see idft, fft.
template<typename S>
auto dft(matrix<S>const& X, std::size_t n=0, int dim=0, bool isinverse=false)
{
    using T = std::complex<decltype(std::abs(X(0)))>;
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    if (n==0) {n = X.size(dim);}
    matrix<T> Y;
    if (dim==0)
    {
        Y = matrix<T>(1,n);
        T arg(0,-2*M_PI/n);
        if (isinverse) {arg *= -1;}
        for (std::size_t k=0; k<numel(Y); ++k)
        {
            for (std::size_t l=0; l<std::min(n,numel(X)); ++l)
            {
                Y(k) += std::exp(arg*(T)(k*l)) * X(l);
            }
            if (isinverse) {Y(k) /= n;}
        }
    }
    else if (dim==1)
    {
        Y = matrix<T>(n,size(X,2));
        T arg(0,-2*M_PI/n);
        if (isinverse) {arg *= -1;}
        for (std::size_t i=0; i<size(Y,1); ++i)
        {
            for (std::size_t j=0; j<size(Y,2); ++j)
            {
                for (std::size_t l=0; l<std::min(n,size(X,1)); ++l)
                {
                    Y(i,j) += std::exp(arg*(T)(i*l)) * X(l,j);
                }
                if (isinverse) {Y(i,j) /= n;}
            }
        }
    }
    else if (dim==2)
    {
        Y = matrix<T>(size(X,1),n);
        T arg(0,-2*M_PI/n);
        if (isinverse) {arg *= -1;}
        for (std::size_t i=0; i<size(Y,1); ++i)
        {
            for (std::size_t j=0; j<size(Y,2); ++j)
            {
                for (std::size_t l=0; l<std::min(n,size(X,2)); ++l)
                {
                    Y(i,j) += std::exp(arg*(T)(j*l)) * X(i,l);
                }
                if (isinverse) {Y(i,j) /= n;}
            }
        }
    }
    return Y;
}

//==========================================================================
// [diag]
/// Diagonal matrices and diagonals of a matrix.
///
/// A = diag(V) is the same as A = diag(V,0) and puts V on the main diagonal.
///
/// \code{.cpp}
///    matrix<> V = {0,1,2,3};
///    matrix<> A = diag(V);
///    disp(A);
/// \endcode
///
/// diag(V,K) when V is a vector with N components is a square matrix
/// of order N+ABS(K) with the elements of V on the K-th diagonal.
/// K = 0 is the main diagonal, K > 0 is above the main diagonal and K < 0
/// is below the main diagonal.
///
/// \code{.cpp}
///    matrix<> V  = {0,1,2,3};
///    matrix<> Al = diag(V,-1);
///    matrix<> Au = diag(V,2);
///    disp(Al);
///    disp(Au);
/// \endcode
///
/// diag(A) is the main diagonal of A. diag(diag(A)) is a diagonal matrix.
///
/// \code{.cpp}
///    matrix<> A = reshape(colon(1,9),3,3);
///    matrix<> V = diag(A);
///    disp(V);
/// \endcode
///
/// diag(A,K) when A is a matrix is a row vector formed from the elements
/// of the K-th diagonal of A.
///
/// \code{.cpp}
///    matrix<> A  = reshape(colon(1,9),3,3);
///    matrix<> Vl = diag(A,-1);
///    matrix<> Vu = diag(A,2);
///    disp(Vl);
///    disp(Vu);
/// \endcode
///
// \see eye.
template<typename T>
matrix<T> diag(matrix<T> const& A, long k=0)
{
    if (std::min(size(A,1),size(A,2))==1)
    {
        std::size_t n = numel(A)+std::abs(k);
        matrix<T> B(n,n,0);
        if (k<0) {for (std::size_t l=0; l<numel(A) ;++l) {B(l-k,l) = A(l);}}
        else {for (std::size_t l=0; l<numel(A) ;++l) {B(l,l+k) = A(l);}}
        return B;
    }
    else
    {
        long r=size(A,1), c=size(A,2), n;
        if (k<0) {n = std::min(r+k,c);}
        else {n = std::min(c-k,r);}
        n = std::max<long>(n,0);
        matrix<T> B(std::min<long>(1,n),n,0);
        if (k<0) {for (std::size_t l=0; l<n ;++l) {B(l) = A(l-k,l);}}
        else {for (std::size_t l=0; l<n ;++l) {B(l) = A(l,l+k);}}
        return B;
    }
}

//==========================================================================
// [diff]
/// Difference of elements.
///
/// B = diff(A) computes the difference of all elements of A.
/// B is a row vector containing numel(A)-1 elements.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = diff(A);
///    disp(B);
/// \endcode
///
/// B = cumsum(A,DIM) computes the cumulative sum along the dimension DIM.
/// if DIM==0, B is the difference of all elements.
/// if DIM==1, B is a vector containing the difference from each column.
/// if DIM==2, B is a vector containing the difference from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = diff(A,1);
///    matrix<> C = diff(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see sum, prod.
template<typename T>
matrix<T> diff(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<T> B;
    if (dim==0)
    {
        B.resize(1,numel(A)-1);
        for (std::size_t l=0; l<numel(B); ++l) {B(l) = A(l+1)-A(l);}
    }
    else if (dim==1)
    {
        B.resize(size(A,1)-1,size(A,2));
        for (std::size_t i=0; i<size(B,1); ++i)
        {
            for (std::size_t j=0; j<size(B,2); ++j) {B(i,j) = A(i+1,j)-A(i,j);}
        }
    }
    else if (dim==2)
    {
        B.resize(size(A,1),size(A,2)-1);
        for (std::size_t i=0; i<size(B,1); ++i)
        {
            for (std::size_t j=0; j<size(B,2); ++j) {B(i,j) = A(i,j+1)-A(i,j);}
        }

    }
    return B;
}
template<typename T>
inline matrix<T> diff(matrix<T>const& A) {return diff(A,0);}

//==========================================================================
// [disp]
/// Display object or array.
///
/// disp(A,INFO) displays object A printing additional description information:\n
/// INFO=0 : No informations, C++ convention, \n
/// INFO=1 : No informations and add end of line, \n
/// INFO=2 : Add additional informations with end of line.
///
/// \code{.cpp}
///    std::string txt = "hello world!";
///    std::cout << std::endl;
///    disp(txt,0);
///    disp(txt,1);
///    disp(txt,2);
/// \endcode
///
/// disp(A,INFO,FLUX) displays object or array A in the stream FLUX.
/// FLUX=std::cout by default.
///
/// \code{.cpp}
///    std::string txt = "hello world!";
///    disp(txt,2,std::cout);
/// \endcode
///
/// disp(A,INFO,FLUX,M,N) only for arrays, specify number of rows M
/// and column N printed at the beggining and the end of the array.
/// M=4 and N=4 by default.
///
/// \code{.cpp}
///    matrix<int> A = eye<int>(8,8);
///    disp(A);
///    disp(A,2,std::cout,8,8);
/// \endcode
///
// \see help, error, warning.
template<typename T>
void disp(T const& A, int info, std::ostream& flux)
{
    if (info<=0)
    {
        flux << A;
    }
    else if (info==1)
    {
        flux << A << std::endl;
    }
    else if (info>=2)
    {
        flux << "Object of type '" << typeid(T).name() << "':" <<std::endl;
        flux << A << std::endl;
    }
}
template<typename T>
void disp(T const* A, int info, std::ostream& flux)
{
    if (info<=0)
    {
        flux << A;
    }
    else if (info==1)
    {
        flux << A << std::endl;
    }
    else if (info>=2)
    {
        flux << "Object of type '" << typeid(T).name() << "':" <<std::endl;
        flux << A << std::endl;
    }
}
template<typename T>
void disp(matrix<T>const& A, int info, std::ostream& flux, std::size_t r, std::size_t c)
{
    // Preparation and informations
    std::size_t m=size(A,1), n=size(A,2);
    double memBytes = m*n*sizeof(T);
    if (info>=2)
    {
        flux << "Matrix " << m << "x" << n;
        flux << " of type '" << typeid(T).name() << "' (";
        if (memBytes<1e3) {flux << memBytes << " B):";}
        else if (memBytes<1e6) {flux << memBytes/1e3 << " KB):";}
        else {flux << memBytes/1e6 << " MB):";}
        flux << std::endl;
    }
    
    // Empty matrix
    if (numel(A)==0)
    {
        flux << "-empty-";
    }
    else
    {
        // Sub-matrix
        matrix<std::size_t> idx, jdx;
        if (m<=2*r) {idx = row(A);}
        else {idx = cat(2,range(0,r),range(m-r,m));}
        if (n<=2*c) {jdx = col(A);}
        else {jdx = cat(2,range(0,c),range(n-c,n));}
        matrix<T> B = eval(A(idx,jdx));
        
        // Adapt precision for integer
        int prc = 4;
        bool isfixed = true;
        if (std::is_integral<T>::value)
        {
            prc = 0;
            for (std::size_t l=0; l<numel(B); ++l)
            {
                prc = std::max<int>(std::log10((double)std::abs<long>(B(l))),prc);
            }
        }
        else
        {
            for (std::size_t l=0; l<numel(B); ++l)
            {
                if ((std::abs<double>(B(l))>1e4 || std::abs<double>(B(l))<1e-3)
                    && std::abs<double>(B(l))!=0)
                {
                    isfixed = false;
                    break;
                }
            }
        }
        
        // Visu
        for (std::size_t i=0; i<size(B,1); ++i)
        {
            if (m>2*r && i==r) {flux << "..." << std::endl;}
            for (std::size_t j=0; j<size(B,2); ++j)
            {
                if (n>2*c && j==c) { flux << "...  ";}
                // Logical
                if (std::is_same<T,logical>::value)
                {
                    flux << B(i,j) << "  ";
                }
                // Integral (int, long, size_t, etc.)
                else if (std::is_integral<T>::value)
                {
                    flux << std::setw(prc+1) << B(i,j) << "  ";
                }
                // Complex
                else if (std::is_same<T,std::complex<float>>::value || std::is_same<T,std::complex<double>>::value)
                {
                    if (isfixed)
                    {
                        flux << std::fixed << std::setprecision(prc+1) << std::setw(2*(prc+7)+3) << B(i,j) << "  ";
                    }
                    else
                    {
                        flux << std::scientific << std::setprecision(prc) << std::setw(2*(prc+7)+3) << B(i,j) << "  ";
                    }
                }
                // Zero
                else if (std::abs<double>(B(i,j))==0)
                {
                    flux << std::defaultfloat << std::setprecision(prc) << std::setw(prc+7) << B(i,j) << "  ";
                }
                // Fixed
                else if (isfixed)
                {
                    flux << std::fixed << std::setprecision(prc+1) << std::setw(prc+7) << B(i,j) << "  ";
                }
                // Scientific
                else if (!isfixed)
                {
                    flux << std::scientific << std::setprecision(prc) << std::setw(prc+7) << B(i,j) << "  ";
                }
            }
            if (i<size(B,1)-1) {flux << std::endl;}
        }
    }
    
    // End
    flux << std::defaultfloat << std::setprecision(6);
    if (info>=1) {flux << std::endl;}
}

//==========================================================================
// [dot]
/// Matrix dot product.
///
/// C = dot(A,B) returns the dot product of the matrix A and B, C = A.B.
/// A and B must have same number of rows with 3 columns.
///
/// \code{.cpp}
///    matrix<> A = {{1,0,0},{0,1,0},{0,0,1}};
///    matrix<> B = {{0,1,0},{0,0,1},{1,0,0}};
///    matrix<> C = dot(A,B);
///    disp(C);
/// \endcode
///
// \see cross.
template<typename R, typename S>
auto dot(matrix<R>const& A, matrix<S>const& B)
{
    if(size(A,1)!=size(B,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix rows must agree.");
    }
    if(size(A,2)!=3 || size(B,2)!=3)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Number of columns should be 3.");
    }
    using T = decltype(A(0)*B(0));
    matrix<T> C(size(A,1),1);
    for(std::size_t l=0; l<numel(C); ++l)
    {
        C(l) = A(l,0)*B(l,0) + A(l,1)*B(l,1) + A(l,2)*B(l,2);
    }
    return C;
}

//==========================================================================
// [error]
/// Display message in command window and abort execution.
///
/// error(__FILE__, __LINE__, __FUNCTION__,"message") displays a descriptive
/// message when the currently-running program encounters an error condition
/// and exit the program with code 1.
///
/// \code{.cpp}
///    error(__FILE__, __LINE__, __FUNCTION__,"This is an error message.");
/// \endcode
///
// \see warning, disp, help.
inline void error(std::string file, int line, std::string function, std::string comment)
{
    std::cout << std::endl;
    std::cout << "Error in " << file << " at line " << line << " with function '";
    std::cout << function << "':" <<std::endl;
    std::cout << comment << std::endl << std::endl;
    exit(1);
}

//==========================================================================
// [exp]
/// Exponential.
///
/// exp(X) is the exponential of the elements of X, e to the X.
/// For complex Z = X + M_1i*Y, exp(Z) = exp(X) * (cos(Y) + M_1i*sin(Y)).
///
/// \code{.cpp}
///    matrix<> A = {1,1,1};
///    matrix<> B = {-M_PI,0,M_PI};
///    auto     X = A + M_1I*B;
///    auto     Y = exp(X);
///    disp(Y);
/// \endcode
///
// \see log, log10.
template<typename S>
auto exp(matrix<S>const& X)
{
    using T = decltype(std::exp(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::exp(X(l));}
    return Y;
}

//==========================================================================
// [eye]
/// Identity matrix.
///
/// eye(N) is the N-by-N identity matrix.
///
/// eye(M,N) or eye({M,N}) is an M-by-N matrix with 1's on the diagonal
/// and zeros elsewhere.
///
/// \code{.cpp}
///    matrix<> A = eye(2);
///    matrix<> B = eye(2,3);
///    matrix<> C = eye(size(B));
///    disp(A);
///    disp(B);
///    disp(C);
/// \endcode
///
// \see ones, zeros, rand.
template<typename T=double>
matrix<T> eye(std::size_t m, long n=-1)
{
    if (n==-1) {n=m;}
    matrix<T> A = matrix<T>(m,n,0);
    for (std::size_t l=0; l<std::min<std::size_t>(m,n); ++l) {A(l,l) = 1;}
    return A;
}
template<typename T=double>
matrix<T> eye(matrix<std::size_t>const& S)
{
    return eye<T>(S(0),S(1));
}

//==========================================================================
// [find]
/// Find linear indices of nonzero elements.
///
/// L = find(X) returns the linear indices corresponding to the nonzero
/// entries of the array X.  X may be a logical expression.
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    matrix<std::size_t> L = find(A);
///    disp(L);
/// \endcode
///
/// Use ind2sub(size(X),L) to calculate bi-linear indices I, J from
/// the linear indices L.
///
/// \code{.cpp}
///    matrix<>            A = eye(3,4);
///    matrix<std::size_t> I, J;
///    std::tie(I,J) = ind2sub(size(A),find(A));
///    disp(I);
///    disp(J);
/// \endcode
///
// \see nnz, ind2sub.
template<typename T>
matrix<std::size_t> find(matrix<T>const& A)
{
    matrix<std::size_t> L(1,nnz(A));
    std::size_t k=0;
    T z=0;
    for(std::size_t l=0; l<numel(A); ++l) {if (A(l)!=z) {L(k)=l; ++k;}}
    return L;
}

//==========================================================================
// [floor]
/// Round towards minus infinity.
///
/// floor(X) rounds the elements of X to the nearest integers
/// towards minus infinity.
///
/// \code{.cpp}
///    matrix<> X = {{1,2,3},{4,5,6}};
///    matrix<> Y = floor(X+0.5);
///    disp(Y);
/// \endcode
///
// \see ceil, round.
template<typename S>
auto floor(matrix<S>const& X)
{
    using T = decltype(std::floor(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::floor(X(l));}
    return Y;
}

//==========================================================================
// [get]
/// Get sub-matrix.
///
/// get(A,L) return matrix A(L) with elements taken from A corresponding to
/// linear indexing L.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = get(A,{0,2,4});
///    disp(B);
/// \endcode
///
/// get(A,I,J) return matrix A(I,J) with elements taken from A corresponding to
/// bilinear indexing I for rows and J for columns.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = get(A,{0,1},{1,2});
///    disp(B);
/// \endcode
///
// \see set, all, row, col, view, cview.
template<typename T>
matrix<T> get(matrix<T>const& A, matrix<std::size_t>const& L)
{
    matrix<T> B(size(L,1),size(L,2));
    for (std::size_t l=0; l<numel(L); ++l)
    {
        B(l) = A(L(l));
    }
    return B;
}
template<typename T>
matrix<T> get(matrix<T>const& A, matrix<std::size_t>const& I, matrix<std::size_t>const& J)
{
    if (std::min(size(I,1),size(I,2))!=1 || std::min(size(J,1),size(J,2))!=1)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Array indices must be vectors.");
    }
    matrix<T> B(numel(I),numel(J));
    for (std::size_t i=0; i<numel(I); ++i)
    {
        for (std::size_t j=0; j<numel(J); ++j)
        {
            B(i,j) = A(I(i),J(j));
        }
    }
    return B;
}

//==========================================================================
// [gmres]
/// Generalized Minimum Residual Method.
///
/// X = gmres(A,B) attempts to solve the system of linear equations A*X = B
/// for X. The N-by-N coefficient matrix A must be square and the right
/// hand side B a N-by-NRHS matrix (MGCR implementation).
///
/// X = gmres(AFUN,B) accepts a lambda function AFUN instead of the matrix
/// A. AFUN(X) accepts a vector input X and returns the matrix-vector
/// product A*X. In all of the following syntaxes, you can replace A by
/// AFUN.
///
/// X = gmres(A,B,TOL) specifies the tolerance of the method.
/// Default is 1e-6.
///
/// X = gmres(A,B,TOL,MAXIT) specifies the maximum number of outer iterations.
/// Default is 10;
///
/// X = gmres(A,B,TOL,MAXIT,AM1) use matrix AM1~=inv(A) as preconditioner
/// and effectively solve the system AM1*A*X = AM1*B for X.
/// Default is AM1 = ID;
/// GMRES accepts a lambda function AM1FUN(X) = AM1*X instead of the explicit
/// matrix AM1.
///
/// X = gmres(A,B,TOL,MAXIT,AM1,X0) specifies the first initial guess.
/// Default is zeros matrix.
///
/// \code{.cpp}
///    matrix<> A = rand(3,3);
///    matrix<> B = eye(3,3);
///    matrix<> X = gmres(A,B);
///    disp(mtimes(X,A));
/// \endcode
///
// \see mtimes, tgemm.
template<typename T>
matrix<T> gmres(std::function<matrix<T>(matrix<T>const&)>const& A, matrix<T>const& B,
                double tol = 1e-6, std::size_t maxit = 10,
                std::function<matrix<T>(matrix<T>const&)>const& Am1 = std::function<matrix<T>(matrix<T>const&)>(),
                matrix<T>const& X0 = matrix<T>(), int info=1)
{
    // Initialize
    if (info>0)
    {
        disp("Start GMRES using MGCR implementation (Multiple Generalized Conjugate Residual).");
    }
    std::size_t nrhs=size(B,2), kmax, iter = 0;
    matrix<std::size_t> I = row(B);
    matrix<T> xk, rk, nrk0, nrk, alphak(1,nrhs), p, Ap, n2Ap, Ar, sigma, tmp;
    T err = 1;
    
    // First step
    if (isempty(X0)) {xk = matrix<T>(size(B,1),size(B,2));}
    else if (size(X0,1)==size(B,1) && size(X0,2)==size(B,2)) {xk = X0;}
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    rk = A(xk)-B;
    if (Am1!=nullptr) {rk = Am1(rk);}
    nrk0    = sqrt(sum(conj(rk)*rk,1));
    p       = -eval(rk(row(rk),0));
    Ap      = A(p);
    if (Am1!=nullptr) {Ap = Am1(Ap);}
    n2Ap    = mtimes(transpose(conj(Ap)),Ap);

    // Loop
    while (std::abs(err)>tol && iter<maxit)
    {
        // Clock
        if (info>1)
        {
            tic();
        }
        
        // Update all RHS
        alphak = - mtimes(transpose(conj(eval(Ap(I,iter)))),rk) / n2Ap(iter);
        tgemm(1,eval(p(I,iter)),alphak,1,xk);
        tgemm(1,eval(Ap(I,iter)),alphak,1,rk);
        
        // Error
        nrk  = sqrt(sum(conj(rk)*rk,1))/nrk0;
        kmax = argmax(abs(nrk));
        err  = nrk(kmax);
        
        // Weights
        Ar = A(eval(rk(I,kmax)));
        if (Am1!=nullptr) {Ar = Am1(Ar);}
        sigma = mtimes(transpose(conj(Ap)),Ar) / transpose(n2Ap);
        
        // Incrementation
        tmp  = -eval(rk(I,kmax)) + mtimes(p,sigma);
        p    = horzcat(p,tmp);
        tmp  = -Ar + mtimes(Ap,sigma);
        Ap   = horzcat(Ap,tmp);
        n2Ap = horzcat(n2Ap,mtimes(transpose(conj(tmp)),tmp));
        ++iter;
        
        // Infos
        if (info>1)
        {
            std::cout << " + Iteration " << iter << " in " << toc(0) <<
            " seconds with relative residual " << err << "." << std::endl;
        }
    }

    // Infos
    if (std::abs(err)<=tol)
    {
        if (info>0)
        {
            std::cout << "GMRES converged at iteration " << iter <<
            " to a solution with relative residual " << err << std::endl;
            std::cout << "and relative error " << max(abs(A(xk)-B))/max(abs(B)) << "." << std::endl;
        }
    }
    else
    {
        std::cout << "GMRES stopped at iteration " << iter <<
        " without converging to the desired tolerance " << tol << std::endl;
        std::cout << "because the maximum number of iterations was reached. The iterate returned has" << std::endl;
        std::cout << "relative residual " << err;
        std::cout << " and relative error " << max(abs(A(xk)-B))/max(abs(B)) << "." << std::endl;
    }
    return xk;
}
template<typename T>
matrix<T> gmres(std::function<matrix<T>(matrix<T>const&)>const& A, matrix<T>const& B,
                double tol = 1e-6, std::size_t maxit = 10,
                matrix<T>const& Am1 = matrix<T>(), matrix<T>const& X0 = matrix<T>(), int info=1)
{
    std::function<matrix<T>(matrix<T>const&)> Am1fct;
    if (!isempty(Am1)) {Am1fct = [&Am1](matrix<T>const& X) {return mtimes(Am1,X);};}
    return gmres(A,B,tol,maxit,Am1fct,X0,info);
}
template<typename T>
matrix<T> gmres(matrix<T>const& A, matrix<T>const& B,
                double tol = 1e-6, std::size_t maxit = 10,
                std::function<matrix<T>(matrix<T>const&)>const& Am1 = std::function<matrix<T>(matrix<T>const&)>(),
                matrix<T>const& X0 = matrix<T>(), int info=1)
{
    std::function<matrix<T>(matrix<T>const&)> Afct;
    Afct = [&A](matrix<T>const& X) {return mtimes(A,X);};
    return gmres(Afct,B,tol,maxit,Am1,X0,info);
}
template<typename T>
matrix<T> gmres(matrix<T>const& A, matrix<T>const& B, double tol, std::size_t maxit,
                matrix<T>const& Am1, matrix<T>const& X0 = matrix<T>(), int info=1)
{
    std::function<matrix<T>(matrix<T>const&)> Afct, Am1fct;
    Afct = [&A](matrix<T>const& X) {return mtimes(A,X);};
    if (!isempty(Am1)) {Am1fct = [&Am1](matrix<T>const& X) {return mtimes(Am1,X);};}
    return gmres(Afct,B,tol,maxit,Am1fct,X0,info);
}

//==========================================================================
// [help]
/// Display help text in command window.
///
/// help("name") displays the help for the functionality specified by name,
/// such as a function, operator symbol, method, class, etc.
/// "name" as to be written in lowercase.
/// Source file(s) as to be given in static variable 'documentationFiles'.
///
/// \code{.cpp}
///    documentationFiles = {"path1/file1.hpp", "path2/file2.hpp", ...};
///    help("help");
/// \endcode
///
/// help("name",{"filename1.hpp","filename2.hpp"}) specify source file(s)
/// which stand for the functionality specified by name.
///
/// \code{.cpp}
///    help("help",{"path1/file1.hpp", "path2/file2.hpp", ...});
/// \endcode
///
/// The documentation offered here is inspired by that of the Matlab software:
/// https://fr.mathworks.com/help/matlab/
///
// \see disp, error, warning.
inline void help(std::string name, std::vector<std::string> filename)
{
    std::cout << "============================ DOCUMENTATION ============================" << std::endl;
    std::vector<std::ifstream> file(filename.size());
    for (int i=0; i<file.size(); ++i) {file[i] = std::ifstream(filename[i]);}
    std::string line, keyword="// ["+name+"]";
    bool print = false;
    for (int i=0; i<file.size(); ++i)
    {
        if (!print)
        {
            if (file[i].is_open())
            {
                while (getline(file[i],line))
                {
                    if (print && line.find("//",0)==std::string::npos) {break;}
                    if (print && line.find("\\endcode")==std::string::npos)
                    {
                        if (line.find("\\code{.cpp}")!=std::string::npos)
                        {
                            line.replace(line.find("\\code{.cpp}"),11,"Example(s):");
                        }
                        if (line.find("\\see")!=std::string::npos)
                        {
                            line.replace(line.find("\\see "),5,"/See also: \n");
                        }
                        if (line.find("\\n")!=std::string::npos)
                        {
                            line.replace(line.find("\\n"),2,"");
                        }
                        std::cout << line.substr(std::min<std::size_t>(4,line.size()),line.size()) << std::endl;
                    }
                    if (line.find(keyword,0)!=std::string::npos)
                    {
                        print = true;
                        std::cout << "Help on \"" + name + "\":\n";
                    }
                }
                file[i].close();
            }
            else
            {
                std::cout << "Source file '" << filename[i] << "' not found." << std::endl;
            }
        }
    }
    if (!print)
    {
        std::cout << "Function \"" << name << "\" not found in file(s):" << std::endl;
        for (int i=0; i<file.size(); ++i) {std::cout << "  " << filename[i] << std::endl;}
    }
    std::cout << "=======================================================================" << std::endl;
}

//==========================================================================
// [horzcat]
/// Horizontal concatenation.
///
/// horzcat(A,B) is the horizontal concatenation of matrices A and B.
/// A and B must have the same number of rows.
/// horzcat(A,B) is the same as cat(2,A,B).
///
/// \code{.cpp}
///    matrix<> A = {1,2,3};
///    matrix<> B = {4,5,6};
///    matrix<> C = horzcat(A,B);
///    disp(C);
/// \endcode
///
// \see vertcat, cat.
template<typename R, typename S>
inline auto horzcat(R const& A, S const& B) {return cat(2,A,B);}

//==========================================================================
// [idft]
/// Inverse Discrete Fourier Transform of an array.
///
/// Y = idft(X) computes the inverse discrete Fourier transform (IDFT) of
/// all the N values of array X using a naive algorithm :
///
/// Y(n) = 1/N sum_{k=0}^(N-1) e^{2i\pi*n * k/N} X(k).
///
/// Y = idft(X,n) returns the n-point IDFT. If no value is specified, Y is the
/// idft of all values of X, in linear indexing.
///
/// Y = idft(X,n,dim) returns the Inverse Fourier transform along the dimension dim.
/// For example, if X is a matrix, then idft(X,n,2) returns the n-point Fourier
/// transform of each row.
///
/// \code{.cpp}
///    matrix<float> X = eye(1,4);
///    disp(idft(X));
/// \endcode
///
// \see dft, ifft.
template<typename S>
inline auto idft(matrix<S>const& X, std::size_t n=0, int dim=0)
{
    return dft(X, n, dim, true);
}

//==========================================================================
// [idx2sph]
/// Discretization index to sherical coordinates (in radian).
///
/// (AZM,ELV) = idx2sph(IDX) gives azimut AZM and elevation ELV coordinates
/// of corresponding indices IDX for a 1-degree sphere.
///
/// \code{.cpp}
///    matrix<std::size_t> I = range(0,360*181);
///    matrix<> AZM, ELV;
///    std::tie(AZM,ELV) = idx2sph(I);
///    matrix<> J = sph2idx(AZM,ELV);
///    disp(max(J-I));
/// \endcode
///
/// (AZM,ELV) = idx2sph(IDX,M) uses M values both in azimut and elevation.
///
/// (AZM,ELV) = idx2sph(IDX,M,N) convert indices matrix IDX, corresponding
/// to M-by-N spherical discretization, to azimut and elevation angles.
///
// \see sph2idx, sphere.
template<typename T=double>
auto idx2sph(matrix<std::size_t>const& idx, std::size_t m=0, std::size_t n=0)
{
    if (m==0 && n==0) {m=181; n=360;}
    else if (n==0) {n=m;}
    matrix<T> azm(size(idx,1),size(idx,2));
    matrix<T> elv(size(idx,1),size(idx,2));
    matrix<T> rho(size(idx,1),size(idx,2),1);
    T dazm = 2*M_PI/n, delv = M_PI/(m-1), pitwo = M_PI/2;
    for (std::size_t l=0; l<numel(idx); ++l)
    {
        azm(l) = (idx(l)%n) * dazm;
        elv(l) = (idx(l)/n) * delv - pitwo;
    }
    std::tie(azm,elv,rho) = sph2cart(azm,elv,rho);
    std::tie(azm,elv,rho) = cart2sph(azm,elv,rho);
    return std::make_tuple(azm,elv);
}

//==========================================================================
// [imag]
/// Complex imaginary part.
///
/// imag(X) is the imaginary part of X.
///
/// \code{.cpp}
///    matrix<> A = {1,0,-1};
///    matrix<> B = {-1,0,1};
///    auto     X = A + M_1I*B;
///    auto     Y = imag(X);
///    disp(Y);
/// \endcode
///
// \see real, conj, angle, abs.
template<typename S>
auto imag(matrix<S>const& X)
{
    using T = decltype(std::imag(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::imag(X(l));}
    return Y;
}
inline matrix<float> imag(matrix<float>const& A) {return matrix<float>(size(A,1),size(A,2),0.);}
inline matrix<double> imag(matrix<double>const& A) {return matrix<double>(size(A,1),size(A,2),0.);}

//==========================================================================
// [ind2sub]
/// Multiple subscripts from linear index.
///
/// ind2sub is used to determine the equivalent subscript values
/// corresponding to a given single index into an array.
///
/// (I,J) = ind2sub(S,L) returns the arrays I and J containing the
/// equivalent row and column subscripts corresponding to the index
/// matrix L for a matrix of size S.
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    matrix<std::size_t> I, J;
///    std::tie(I,J) = ind2sub(size(A),find(A));
///    disp(I);
///    disp(J);
/// \endcode
///
// \see sub2ind, find.
template<typename T=std::size_t>
std::tuple<matrix<std::size_t>,matrix<std::size_t>> ind2sub(matrix<std::size_t> S, matrix<T>const& L)
{
    if (numel(S)!=2)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Size vector must have 2 elements.");
    }
    matrix<std::size_t> I(size(L,1),size(L,2)), J(size(L,1),size(L,2));
    for (std::size_t l=0; l<numel(L); ++l)
    {
        I(l) = L(l)/S(1);
        J(l) = L(l)%S(1);
    }
    return std::make_tuple(I,J);
}

//==========================================================================
// [intersect]
/// Set intersection of two arrays.
/// C = intersect(A,B) for matrix A and B returns a vector with the values
/// common to the two matrix, with no repetitions. C will be sorted.
///
/// \code{.cpp}
///    matrix<> A = {1,2,3};
///    matrix<> B = {2,3,4};
///    matrix<> C = intersect(A,B);
///    disp(C);
/// \endcode
///
// \see argintersect, setdiff, unique, union2.
template<typename R, typename S>
auto intersect(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)+B(0));
    std::vector<T> a=A.val(), b=B.val(), c;
    std::sort(a.begin(),a.end());
    std::sort(b.begin(),b.end());
    std::set_intersection(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(c));
    auto last = std::unique(c.begin(),c.end());
    c.erase(last,c.end());
    return matrix<T>(c);
}

//==========================================================================
// [isempty]
/// True for empty array.
///
/// isempty(X) returns 'true' if X is an empty array and 'false' otherwise.
/// An empty array has no elements, that is prod(size(X))==0.
///
/// \code{.cpp}
///    bool a = isempty(zeros(0,0));
///    bool b = isempty(zeros(2,3));
///    disp(a);
///    disp(b);
/// \endcode
///
// \see isequal, isvector.
template<typename T>
inline bool isempty(matrix<T>const& X) {return numel(X)==0;}

//==========================================================================
// [isequal]
/// True if arrays are numerically equal.
///
/// isequal(A,B) returns 'true' if arrays A and B are the same
/// size and contain the same values, and 'false' otherwise.
///
/// \code{.cpp}
///    bool a = isequal(ones(3,4),ones(3,4));
///    bool b = isequal(ones(3,4),zeros(3,4));
///    disp(a);
///    disp(b);
/// \endcode
///
// \see isempty, isvector.
template<typename R, typename S>
bool isequal(matrix<R>const& A, matrix<S>const& B)
{
    if (size(A,1)!=size(B,1) || size(A,2)!=size(B,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    bool b = true;
    for (std::size_t l=0; l<numel(A); ++l)
    {
        if (A(l)!=B(l))
        {
            b = false;
            break;
        }
    }
    return b;
}

//==========================================================================
// [isfinite]
/// True for finite element.
///
/// isfinite(A) returns an array B that contains 1's where the elements of A
/// are finite and 0's where they are not.
///
/// \code{.cpp}
///    matrix<> A = eye(2,3);
///    A(0) = NAN;
///    A(1) = INFINITY;
///    A(2) = exp(800);
///    disp(isfinite(A);
/// \endcode
///
// \see isnan, isinf.
template<typename T>
matrix<logical> isfinite(matrix<T>const& A)
{
    matrix<logical> B(size(A,1),size(A,2));
    for (std::size_t l=0; l<numel(A); ++l)
    {
        B(l) = std::isfinite(std::abs(A(l)));
    }
    return B;
}

//==========================================================================
// [isinf]
/// True for infinite element.
///
/// isinf(A) returns an array B that contains 1's where the elements of A
/// are infinite and 0's where they are not.
///
/// \code{.cpp}
///    matrix<> A = eye(2,3);
///    A(0) = INFINITY;
///    disp(isinf(A);
/// \endcode
///
// \see isnan, isfinite.
template<typename T>
matrix<logical> isinf(matrix<T>const& A)
{
    matrix<logical> B(size(A,1),size(A,2));
    for (std::size_t l=0; l<numel(A); ++l)
    {
        B(l) = std::isinf(std::abs(A(l)));
    }
    return B;
}

//==========================================================================
// [isnan]
/// True for Not-A-Number.
///
/// isnan(A) returns an array B that contains 1's where the elements of A
/// are Not-A-Number's and 0's where they are not.
///
/// \code{.cpp}
///    matrix<> A = eye(2,3);
///    A(0) = NAN;
///    disp(isnan(A);
/// \endcode
///
// \see isinf, isfinite.
template<typename T>
matrix<logical> isnan(matrix<T>const& A)
{
    matrix<logical> B(size(A,1),size(A,2));
    for (std::size_t l=0; l<numel(A); ++l)
    {
        B(l) = std::isnan(std::abs(A(l)));
    }
    return B;
}

//==========================================================================
// [isvector]
/// True for monodimensional matrix.
///
/// isvector(X) returns 'true' if X is a monodimensional array and 'false'
/// otherwise. A monodimensional array has at least one dimension equal to 1.
///
/// \code{.cpp}
///    bool a = isvector(ones(1,4));
///    bool b = isvector(ones(3,4));
///    disp(a);
///    disp(b);
/// \endcode
///
// \see isequal, isempty.
template<typename T>
inline bool isvector(matrix<T>const& X) {return std::min(size(X,1),size(X,2))==1;}

//==========================================================================
// [kron]
/// Kronecker tensor product.
///
/// kron(A,B) is the Kronecker tensor product of A and B.
/// The result is a large matrix formed by taking all possible
/// products between the elements of A and those of B. For
/// example, if A is 2 by 3, then kron(A,B) is
/// [ A(1,1)*B  A(1,2)*B  A(1,3)*B
///   A(2,1)*B  A(2,2)*B  A(2,3)*B ]
///
/// \code{.cpp}
///    matrix<> A = eye(2,3);
///    matrix<> B = ones(3,2);
///    matrix<> C = kron(A,B);
///    disp(C);
/// \endcode
///
// \see mtimes.
template<typename R, typename S>
auto kron(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)*B(0));
    std::size_t ma=size(A,1), mb=size(B,1), na=size(A,2), nb=size(B,2);
    matrix<T> C(ma*mb,na*nb);
    for (std::size_t ia=0; ia<ma; ++ia)
    {
        for (std::size_t ja=0; ja<na; ++ja)
        {
            for (std::size_t ib=0; ib<mb; ++ib)
            {
                for (std::size_t jb=0; jb<nb; ++jb)
                {
                    C(ia*mb+ib,ja*nb+jb) = A(ia,ja)*B(ib,jb);
                }
            }
        }
    }
    return C;
}
template<typename R, typename S>
inline auto kron(R A, matrix<S>const& B) {return kron(matrix<R>(A),B);}
template<typename R, typename S>
inline auto kron(matrix<R>const& A, S B) {return kron(A,matrix<S>(B));}

//==========================================================================
// [length]
/// Length of vector.
///
/// length(V) returns the length of vector V. For matrix A, it is equivalent
/// to max(size(A)).
///
/// \code{.cpp}
///    matrix<> A = ones(1,3);
///    matrix<> B = ones(4,3);
///    disp(length(A));
///    disp(length(B));
/// \endcode
///
// \see numel, size.
template<typename T>
std::size_t length(matrix<T>const& A)
{
    return std::max(size(A,1),size(A,2));
}

//==========================================================================
// [linspace]
/// Linearly spaced vector.
///
/// linspace(X1, X2) generates a row vector of 100 linearly equally spaced
/// points between X1 and X2.
///
/// linspace(X1, X2, N) generates N points. For N=1, linspace returns X2.
///
/// \code{.cpp}
///    matrix<> A = linspace(1,2);
///    matrix<> B = linspace(1,2,5);
///    disp(A);
///    disp(B);
/// \endcode
///
// \see logspace, colon.
template<typename T=double, typename U, typename V>
matrix<T> linspace(U x1, V x2, std::size_t n=100)
{
    if (n==0) {return matrix<T>();}
    else if (n==1) {return matrix<T>(x2);}
    else
    {
        matrix<T> A(1,n);
        T stp = (x2-x1)/((T)n-1);
        for (std::size_t l=0; l<n; ++l) {A(l) = x1 + l*stp;}
        return A;
    }
}

//==========================================================================
// [log]
/// Natural logarithm.
///
/// log(X) is the natural logarithm of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1,std::exp(1)};
///    auto     Y = log(X);
///    disp(Y);
/// \endcode
///
// \see log2, log10, exp.
template<typename S>
auto log(matrix<S>const& X)
{
    using T = decltype(std::log(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::log(X(l));}
    return Y;
}

//==========================================================================
// [log2]
/// Common (base 2) logarithm.
///
/// log2(X) is the base 2 logarithm of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1,2,4};
///    auto     Y = log2(X);
///    disp(Y);
/// \endcode
///
// \see log, log10, exp.
template<typename S>
auto log2(matrix<S>const& X)
{
    using T = decltype(std::log2(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::log2(X(l));}
    return Y;
}

//==========================================================================
// [log10]
/// Common (base 10) logarithm.
///
/// log10(X) is the base 10 logarithm of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1,10,100};
///    auto     Y = log10(X);
///    disp(Y);
/// \endcode
///
// \see log, log2, exp.
template<typename S>
auto log10(matrix<S>const& X)
{
    using T = decltype(std::log10(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::log10(X(l));}
    return Y;
}

//==========================================================================
// [logspace]
/// Logarithmically spaced vector.
///
/// logspace(X1, X2) generates a row vector of 50 logarithmically
/// equally spaced points between decades 10^X1 and 10^X2.
///
/// logspace(X1, X2, N) generates N points. For N=1, logspace returns 10^X2.
///
/// \code{.cpp}
///    matrix<> A = logspace(0,3);
///    matrix<> B = logspace(0,3,4);
///    disp(A);
///    disp(B);
/// \endcode
///
// \see linspace, colon, log.
template<typename T=double, typename U, typename V>
matrix<T> logspace(U x1, V x2, std::size_t n=50)
{
    matrix<T> A = linspace<T>(x1,x2,n);
    for (std::size_t l=0; l<numel(A); ++l) {A(l) = std::pow(10,A(l));}
    return A;
}

//==========================================================================
// [max]
/// Maximum elements of an array.
///
/// M = max(A) is the largest element in the array A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    double   m = max(A);
///    disp(m);
/// \endcode
///
/// M = max(A,DIM) operates along the dimension DIM:
/// if DIM==0, M is the maximum element of the array.
/// if DIM==1, M is a vector containing the maximum element from each column.
/// if DIM==2, M is a vector containing the maximum element from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = max(A,1);
///    matrix<> C = max(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see argmax, maximum, min, sort.
template<typename T>
matrix<T> max(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<T> B;
    if (dim==0)
    {
        B = A(0);
        for (std::size_t l=1; l<numel(A); ++l)  {if(A(l)>B(0)) {B(0)=A(l);}}
    }
    else if (dim==1)
    {
        B.resize(1,size(A,2));
        for (std::size_t j=0; j<size(A,2); ++j)
        {
            B(j) = A(0,j);
            for (std::size_t i=1; i<size(A,1); ++i) {if(A(i,j)>B(j)) {B(j)=A(i,j);}}
        }
    }
    else if (dim==2)
    {
        B.resize(size(A,1),1);
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            B(i) = A(i,0);
            for (std::size_t j=1; j<size(A,2); ++j) {if(A(i,j)>B(i)) {B(i)=A(i,j);}}
        }
    }
    return B;
}
template<typename T>
inline T max(matrix<T>const& A) {return (T)max(A,0);}

//==========================================================================
// [maximum]
/// Maximum elements between two arrays.
///
/// C = maximum(A,B) returns an array with the largest elements taken from A
/// or B. A and B must have compatible sizes. In the simplest cases, they can
/// be the same size or one can be a scalar. Two inputs have compatible
/// sizes if, for every dimension, the dimension sizes of the inputs are
/// the same.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = {{6,5,4},{3,2,1}};
///    matrix<> C = maximum(A,B);
///    matrix<> D = maximum(A,3);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see minimum, max, min, sort.
template<typename R, typename S>
auto maximum(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)+B(0));
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<T> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l)
        {
            if (A(0)>B(l)) {C(l) = A(0);}
            else {C(l) = B(l);}
        }
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l)
        {
            if (A(l)>B(0)) {C(l) = A(l);}
            else {C(l) = B(0);}
        }
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l)
        {
            if (A(l)>B(l)) {C(l) = A(l);}
            else {C(l) = B(l);}
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline auto maximum(R A, matrix<S>const& B) {return maximum(matrix<R>(A),B);}
template<typename R, typename S>
inline auto maximum(matrix<R>const& A, S B) {return maximum(A,matrix<S>(B));}

//==========================================================================
// [mean]
/// Average or mean value.
///
/// M = mean(A) is the mean value of all the elements in the array A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    double   m = mean(A);
///    disp(m);
/// \endcode
///
/// M = mean(A,DIM) takes the mean along the dimension DIM of A.
/// if DIM==0, M is the mean value of the array.
/// if DIM==1, M is a vector containing the mean value from each column.
/// if DIM==2, M is a vector containing the mean value from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = mean(A,1);
///    matrix<> C = mean(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see max, min, median, variance, stddev.
template<typename T>
matrix<T> mean(matrix<T>const& A, int dim)
{
    matrix<T> B = sum(A,dim);
    if (dim==0) {B /= (T)numel(A);}
    else {B /= (T)size(A,dim);}
    return B;
}
template<typename T>
inline T mean(matrix<T>const& A) {return (T)mean(A,0);}

//==========================================================================
// [median]
/// Median value.
///
/// M = median(A) is the median value of the elements in A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    double   m = median(A);
///    disp(m);
/// \endcode
///
/// M = median(A,DIM) takes the median along the dimension DIM of A.
/// if DIM==0, M is the mediane value of the array.
/// if DIM==1, M is a vector containing the mediane value from each column.
/// if DIM==2, M is a vector containing the mediane value from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = median(A,1);
///    matrix<> C = median(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see max, min, mean, variance, stddev.
template<typename T>
matrix<T> median(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<T> As = sort(A,dim), Av;
    if (dim==0)
    {
        Av = eval( As({(numel(A)-1)/2,numel(A)/2}) );
    }
    else if (dim==1)
    {
        Av = eval( As({(size(A,1)-1)/2,size(A,1)/2},col(A)) );
    }
    else if (dim==2)
    {
        Av = eval( As(row(A),{(size(A,2)-1)/2,size(A,2)/2}) );
    }
    return mean(Av,dim);
}
template<typename T>
inline T median(matrix<T>const& A) {return (T)median(A,0);}

//==========================================================================
// [meshgrid]
/// Cartesian rectangular grid in 2-D.
///
/// [X,Y] = meshgrid(x,y) returns 2-D grid coordinates based on the
/// coordinates contained in vectors x and y. X is a matrix where each row
/// is a copy of x, and Y is a matrix where each column is a copy of y. The
/// grid represented by the coordinates X and Y has length(y) rows and
/// length(x) columns.
///
/// [X,Y] = meshgrid(x) is the same as [X,Y] = meshgrid(x,x), returning
/// square grid coordinates with grid size length(x)-by-length(x).
///
/// \code{.cpp}
///    matrix<> X,Y;
///    std::tie(X,Y) = meshgrid(linspace(-1,1,10),linspace(0,2,5));
///    disp(X);
///    disp(Y);
/// \endcode
///
// \see kron.
template<typename T>
auto meshgrid(matrix<T>const& x, matrix<T>const& y={})
{
    std::size_t m, n=numel(x);
    if (numel(y)>0) {m=numel(y);}
    else {m=numel(x);}
    matrix<T> X(m,n), Y(m,n);
    for (std::size_t i=0; i<m; ++i)
    {
        for (std::size_t j=0; j<n; ++j)
        {
            X(i,j) = x(j);
            if (numel(y)>0) {Y(i,j) = y(i);}
            else {Y(i,j) = x(i);}
        }
    }
    return std::make_tuple(X,Y);
}

//==========================================================================
// [min]
/// Minimum elements of an array.
///
/// M = min(A) is the smallest element in the array A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    double   m = min(A);
///    disp(m);
/// \endcode
///
/// M = min(A,DIM) operates along the dimension DIM:
/// if DIM==0, M is the minimum element.
/// if DIM==1, M is a vector containing the minimum element from each column.
/// if DIM==2, M is a vector containing the minimum element from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = min(A,1);
///    matrix<> C = min(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see argmin, minimum, max, sort.
template<typename T>
matrix<T> min(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<T> B;
    if (dim==0)
    {
        B = A(0);
        for (std::size_t l=1; l<numel(A); ++l) {if(A(l)<B(0)) {B(0)=A(l);}}
    }
    else if (dim==1)
    {
        B.resize(1,size(A,2));
        for (std::size_t j=0; j<size(A,2); ++j)
        {
            B(j) = A(0,j);
            for (std::size_t i=1; i<size(A,1); ++i) {if(A(i,j)<B(j)) {B(j)=A(i,j);}}
        }
    }
    else if (dim==2)
    {
        B.resize(size(A,1),1);
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            B(i) = A(i,0);
            for (std::size_t j=1; j<size(A,2); ++j) {if(A(i,j)<B(i)) {B(i)=A(i,j);}}
        }
    }
    return B;
}
template<typename T>
T min(matrix<T>const& A) {return (T)min(A,0);}

//==========================================================================
// [minimum]
/// Minimum elements between two arrays.
///
/// C = minimum(A,B) returns an array with the smallest elements taken from A
/// or B. A and B must have compatible sizes. In the simplest cases, they
/// can be the same size or one can be a scalar. Two inputs have compatible
/// sizes if, for every dimension, the dimension sizes of the inputs are
/// the same.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = {{6,5,4},{3,2,1}};
///    matrix<> C = minimum(A,B);
///    matrix<> D = minimum(A,3);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see maximum, min, max, sort.
template<typename R, typename S>
auto minimum(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)+B(0));
    std::size_t m = std::max(size(A,1),size(B,1))*(numel(A)>0)*(numel(B)>0);
    std::size_t n = std::max(size(A,2),size(B,2))*(numel(A)>0)*(numel(B)>0);
    matrix<T> C(m,n);
    if (numel(A)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l)
        {
            if (A(0)<B(l)) {C(l) = A(0);}
            else {C(l) = B(l);}
        }
    }
    else if (numel(B)==1)
    {
        for (std::size_t l=0; l<numel(C); ++l)
        {
            if (A(l)<B(0)) {C(l) = A(l);}
            else {C(l) = B(0);}
        }
    }
    else if (size(A,1)==size(B,1) && size(A,2)==size(B,2))
    {
        for (std::size_t l=0; l<numel(C); ++l)
        {
            if (A(l)<B(l)) {C(l) = A(l);}
            else {C(l) = B(l);}
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    return C;
}
template<typename R, typename S>
inline auto minimum(R A, matrix<S>const& B) {return minimum(matrix<R>(A),B);}
template<typename R, typename S>
inline auto minimum(matrix<R>const& A, S B) {return minimum(A,matrix<S>(B));}

//==========================================================================
// [mtimes]
/// Matrix multiply.
///
/// mtimes(A,B) is the matrix product of A and B, the number of columns
/// of A must equal the number of rows of B.
///
/// \code{.cpp}
///    matrix<> A = ones(3,1);
///    matrix<> B = ones(1,2);
///    matrix<> C = mtimes(A,B);
///    disp(C);
/// \endcode
///
// \see tgemm, kron.
template<typename R, typename S>
auto mtimes(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)*B(0));
    matrix<T> C(size(A,1),size(B,2));
    if (std::is_same<R,T>::value && std::is_same<S,T>::value)
    {
        tgemm((T)1,A,B,(T)0,C);
    }
    else if (std::is_same<R,T>::value)
    {
        tgemm((T)1,A,cast<T>(B),(T)0,C);
    }
    else if (std::is_same<S,T>::value)
    {
        tgemm((T)1,cast<T>(A),B,(T)0,C);
    }
    else
    {
        tgemm((T)1,cast<T>(A),cast<T>(B),(T)0,C);
    }
    return C;
}
template<typename R, typename S>
inline auto mtimes(R A, matrix<S>const& B) {return mtimes(matrix<R>(A),B);}
template<typename R, typename S>
inline auto mtimes(matrix<R>const& A, S B) {return mtimes(A,matrix<S>(B));}

//==========================================================================
// [nnz]
/// Number of nonzero matrix elements.
///
/// nz = nnz(A) is the number of nonzero elements in A.
///
/// \code{.cpp}
///    matrix<>    A  = eye(3,4);
///    std::size_t nz = nnz(A);
///    disp(nz);
/// \endcode
///
// \see find, size.
template<typename T>
std::size_t nnz(matrix<T>const& A)
{
    std::size_t nz = 0;
    for (std::size_t l=0; l<numel(A); ++l)
    {
        if (A(l)!=0) {++nz;}
    }
    return nz;
}

//==========================================================================
// [norm]
/// Vectorial norm applied to all values of a matrix.
/// Vectorial norms are not matrix norm (e.g. Frobenius, SVD, etc.), but
/// norms applied to all element of a matrix, using linear indexing.
///
/// norm(X) is the euclidian norm.
///
/// \code{.cpp}
///    matrix<> A   = eye(3,4);
///    double   nrm = norm(A);
///    disp(nrm);
/// \endcode
///
/// norm(X,TYP) is the TYP-norm of X :
/// TYP=1 returns the 1-norm of X,
/// TYP=2 returns the 2-norm of X,
/// TYP="inf" returns the infinite norm of X.
///
/// \code{.cpp}
///    matrix<> A    = eye(3,4);
///    double   nrm1 = norm(A,"1");
///    double   nrm2 = norm(A,"2");
///    double   nrmI = norm(A,"inf");
///    disp(nrm1);
///    disp(nrm2);
///    disp(nrmI);
/// \endcode
///
/// norm(X,TYP,DIM) returns the TYP-norm of X along the dimension DIM.
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    matrix<> R = norm(A,"inf",1);
///    matrix<> C = norm(A,"inf",2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see max, sum.
template<typename S>
auto norm(matrix<S>const& A, std::string typ, int dim)
{
    using T = decltype(std::abs(A(0)));
    matrix<T> nrm;
    if (typ=="1") {nrm = sum(abs(A),dim);}
    else if (typ=="2") {nrm = sqrt(sum(pow(abs(A),2),dim));}
    else if (typ=="inf" || typ=="INF" || typ=="Inf") {nrm = max(abs(A),dim);}
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"The only norms available are 1, 2 and inf.");
    }
    return nrm;
}
template<typename S>
auto norm(matrix<S>const& A, std::string typ="2")
{
    using T = decltype(std::abs(A(0)));
    return (T)norm(A,typ,0);
}

//==========================================================================
// [numel]
/// Number of elements in a array.
///
/// N = numel(A) returns the number of elements, N, in array A, equivalent
/// to prod(size(A)).
///
/// \code{.cpp}
///    matrix<>    A = ones(3,4);
///    std::size_t n = numel(A);
///    disp(n);
/// \endcode
///
// \see length, size.
template<typename T>
std::size_t numel(matrix<T>const& A)
{
    return A.size();
}

//==========================================================================
// [ones]
/// Ones matrix.
///
/// ones(N) is an N-by-N matrix of ones.
///
/// ones(M,N) or ones({M,N}) is an M-by-N matrix of ones.
///
/// \code{.cpp}
///    matrix<> A = ones(2);
///    matrix<> B = ones(2,3);
///    matrix<> C = ones(size(B));
///    disp(A);
///    disp(B);
///    disp(C);
/// \endcode
///
// \see zeros, eye, rand.
template<typename T=double>
matrix<T> ones(std::size_t m, long n=-1)
{
    if (n==-1) {n=m;}
    return matrix<T>(m,n,1);
}
template<typename T=double>
matrix<T> ones(matrix<std::size_t>const& S)
{
    return matrix<T>(S(0),S(1),1);
}

//==========================================================================
// [pol2cart]
/// Transform polar to Cartesian coordinates.
///
/// (X,Y) = pol2cart(THE,RHO) transforms corresponding elements of data
/// stored in polar coordinates (THE,RHO) to Cartesian coordinates (X,Y).
/// The arrays THE and RHO must the same size and angle THE must be in radians.
///
/// <b>Remark:</b> This function supports scalar arguments.
///
/// \code{.cpp}
///    matrix<> THE = {0,M_PI/2,M_PI,-M_PI/2};
///    matrix<> RHO = {1,1,1,1};
///    matrix<> X, Y;
///    std::tie(X,Y) = pol2cart(THE,RHO);
///    disp(X);
///    disp(Y);
/// \endcode
///
// \see cart2pol, cart2sph, sph2cart.
template<typename R, typename S>
auto pol2cart(matrix<R>const& THE, matrix<S>const& RHO)
{
    if (size(THE,1)!=size(RHO,1) || size(THE,2)!=size(RHO,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(THE(0)+RHO(0));
    matrix<T> X(size(THE,1),size(THE,2)), Y(size(THE,1),size(THE,2));
    for (std::size_t l=0; l<numel(X); ++l)
    {
        X(l) = RHO(l) * std::cos(THE(l));
        Y(l) = RHO(l) * std::sin(THE(l));
    }
    return std::make_tuple(X,Y);
}
template<typename R, typename S>
inline auto pol2cart(R const &THE, S const &RHO)
{
    using T = decltype(THE + RHO);
    T X = RHO*std::cos(THE);
    T Y = RHO*std::sin(THE);
    return std::make_tuple(X,Y);
}

//==========================================================================
// [pow]
/// Power of argument.
///
/// pow(A,B) denotes element-by-element powers. A and B must have
/// compatible sizes. In the simplest cases, they can be the same size or
/// one can be a scalar. Two inputs have compatible sizes if, for every
/// dimension, the dimension sizes of the inputs are either the same or one
/// of them is 1.
///
/// \code{.cpp}
///    matrix<> A = {-2,-1,0,1,2};
///    matrix<> B = {-2,-1,0,1,2};
///    matrix<> C = pow(A,B);
///    matrix<> D = pow(A,2);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see exp, log.
template<typename R, typename S>
auto pow(matrix<R>const& X, matrix<S>const& Y)
{
    if (size(X,1)!=size(Y,1) || size(X,2)!=size(Y,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"matrix dimensions must agree.");
    }
    using T = decltype(std::pow(X(0),Y(0)));
    matrix<T> Z(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Z); ++l) {Z(l) = std::pow(X(l),Y(l));}
    return Z;
}
template<typename R, typename S>
auto pow(matrix<R>const& X, S y)
{
    using T = decltype(std::pow(X(0),y));
    matrix<T> Z(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Z); ++l) {Z(l) = std::pow(X(l),y);}
    return Z;
}
template<typename R, typename S>
auto pow(R x, matrix<S>const& Y)
{
    using T = decltype(std::pow(x,Y(0)));
    matrix<T> Z(size(Y,1),size(Y,2));
    for (std::size_t l=0; l<numel(Z); ++l) {Z(l) = std::pow(x,Y(l));}
    return Z;
}

//==========================================================================
// [prod]
/// Product of elements.
///
/// P = prod(A) is the product of all the elements of the array A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    double   p = prod(A);
///    disp(p);
/// \endcode
///
/// prod(A,DIM) operates along the dimension DIM:
/// if DIM==0, M is the product of all elements.
/// if DIM==1, M is a vector containing the prpduct from each column.
/// if DIM==2, M is a vector containing the product from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = prod(A,1);
///    matrix<> C = prod(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see sum, diff.
template<typename T>
matrix<T> prod(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<T> B;
    if (dim==0)
    {
        B.resize(1,1,1);
        for (std::size_t l=0; l<numel(A); ++l) {B(0) *= A(l);}
    }
    else if (dim==1)
    {
        B.resize(1,size(A,2),1);
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            for (std::size_t j=0; j<size(A,2); ++j) {B(j) *= A(i,j);}
        }
    }
    else if (dim==2)
    {
        B.resize(size(A,1),1,1);
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            for (std::size_t j=0; j<size(A,2); ++j) {B(i) *= A(i,j);}
        }
    }
    return B;
}
template<typename T>
inline T prod(matrix<T>const& A) {return (T)prod(A,0);}

//==========================================================================
// [rad2deg]
/// Convert angles from radians to degrees.
///
/// rad2deg(X) converts angle units from radians to degrees for each
/// element of X.
///
/// \code{.cpp}
///    matrix<> X = {0,M_PI/2,M_PI,-M_PI/2};
///    matrix<> Y = rad2deg(X);
///    disp(Y);
/// \endcode
///
// \see deg2rad.
template<typename T>
matrix<T> rad2deg(matrix<T>const& X)
{
    matrix<T> Y = X;
    double c = 180/M_PI;
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) *= c;}
    return Y;
}
template<typename T>
auto rad2deg(T x)
{
    return x*180/M_PI;
}

//==========================================================================
// [rand]
/// matrix of uniformly distributed pseudorandom numbers.
///
/// rand(N) returns an N-by-N matrix containing pseudorandom values drawn
/// from the standard uniform distribution on the open interval (0,1).
///
/// rand(M,N) or rand({M,N}) returns an M-by-N matrix.
///
/// rand(M,N,true) use standard time to initialize seed.
///
/// \code{.cpp}
///    matrix<> A = rand(2);
///    matrix<> B = rand(2,3);
///    matrix<> C = rand(size(B));
///    matrix<> D = rand(2,3,true);
///    disp(A);
///    disp(B);
///    disp(C);
///    disp(D);
/// \endcode
///
// \see zeros, eye, ones.
template<typename T=double>
matrix<T> rand(std::size_t m, long n=-1, bool seed=false)
{
    if (n==-1) {n=m;}
    matrix<T> A = matrix<T>(m,n);
    T tmp = RAND_MAX;
    if (seed) {std::srand((int)std::time(0));}
    for (std::size_t l=0; l<m*n; ++l) {A(l) = std::rand()/tmp;}
    return A;
}
template<typename T=double>
matrix<T> rand(matrix<std::size_t>const& S, bool seed=false)
{
    return rand<T>(S(0),S(1),seed);
}

//==========================================================================
// [range]
/// Indices for matrix views.
///
/// range(J,K) is simply [J, J+1, ..., K[, respecting c++ numbering.
/// This syntax returns an empty matrix if J >= K.
///
/// \code{.cpp}
///    matrix<> A = range(0,5);
///    disp(A);
/// \endcode
///
// \see colon, get, set, view.
inline matrix<std::size_t> range(std::size_t j, std::size_t k)
{
    return colon(j,1,k-1);
}

//==========================================================================
// [readbin]
/// Read matrix stored in *.bin file.
///
/// readbin("path","filename") read binary file with matrix stored
/// in the following order :
/// m n M(0,0) M(0,1) ... M(i,j) ... M(m-1,n-2) M(m-1,n-1)
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    writebin("./","matrix.bin",A);
///    matrix<> B = readbin("./","matrix.bin");
///    disp(B);
/// \endcode
///
// \see writebin, readtxt.
template<typename T=double>
matrix<T> readbin(std::string path, std::string filename)
{
    std::ifstream file(path+filename, std::ios::in | std::ios::binary);
    matrix<T> A;
    std::size_t m, n;
    bool t;
    T tmp;
    if (file)
    {
        file.read((char*)&m,sizeof(std::size_t));
        file.read((char*)&n,sizeof(std::size_t));
        A.resize(m,n);
        for (std::size_t l=0; l<m*n; ++l)
        {
            file.read((char*)&A(l),sizeof(T));
        }
        t = file.eof();
        file.read((char*)&tmp,sizeof(T));
        if (t || !file.eof())
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Error reading data in "+path+filename+" file.");
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"File "+path+filename+" not found.");
    }
    file.close();
    return A;
}

//==========================================================================
// [readtxt]
/// Read matrix stored in *.txt file.
///
/// readtxt("path","filename") read ascii text file with matrix stored
/// in the following format :
/// m       n
/// M(0,0)    M(0,1)    ...  M(0,n-2)    M(0,n-1)
///  ...                                  ...
/// M(m-1,0)  M(m-1,1)  ...  M(m-1,n-2)  M(m-1,n-1)
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    writetxt("./","matrix.txt",A);
///    matrix<> B = readtxt("./","matrix.txt");
///    disp(B);
/// \endcode
///
// \see writetxt, readbin.
template<typename T=double>
matrix<T> readtxt(std::string path, std::string filename)
{
    std::ifstream file(path+filename);
    matrix<T> A;
    std::size_t m, n;
    bool t;
    T tmp;
    if(file)
    {
        file >> m >> n;
        A.resize(m,n);
        for (std::size_t i=0; i<m; ++i)
        {
            for (std::size_t j=0; j<n; ++j)
            {
                file >> A(i,j);
            }
        }
        t = file.eof();
        file >> tmp;
        if (t || !file.eof())
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Error reading data in "+path+filename+" file.");
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"File "+path+filename+" not found.");
    }
    file.close();
    return A;
}

//==========================================================================
// [real]
/// Complex real part.
///
/// real(X) is the real part of X.
///
/// \code{.cpp}
///    matrix<> A = {1,0,-1};
///    matrix<> B = {-1,0,1};
///    auto     X = A + M_1I*B;
///    auto     Y = real(X);
///    disp(Y);
/// \endcode
///
// \see imag, conj, angle, abs.
template<typename S>
auto real(matrix<S>const& X)
{
    using T = decltype(std::real(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::real(X(l));}
    return Y;
}
inline matrix<float> real(matrix<float>const& A) {return A;}
inline matrix<double> real(matrix<double>const& A) {return A;}

//==========================================================================
// [reshape]
/// Reshape array.
///
/// reshape(A,M,N) returns the M-by-N matrix whose elements are taken from A.
/// An error results if A does not have M*N elements.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = reshape(A,3,2);
///    disp(B);
/// \endcode
///
// \see resize, transpose.
template<typename T>
matrix<T> reshape(matrix<T>const& A, std::size_t m, std::size_t n)
{
    matrix<T> B = A;
    B.reshape(m,n);
    return B;
}

//==========================================================================
// [resize]
/// Resize array.
///
/// resize(A,M,N,V) resize matrix A with M row and N columns, eventually
/// filled with value V (default is NAN). Original values of A are
/// conserved.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = resize(A,3,4);
///    disp(B);
/// \endcode
///
// \see reshape.
template<typename T>
matrix<T> resize(matrix<T>const& A, std::size_t m, std::size_t n, T v=(T)NAN)
{
    matrix<T> B = A;
    B.resize(m,n,v);
    return B;
}

//==========================================================================
// [round]
/// Rounds towards nearest decimal or integer.
///
/// round(X) rounds each element of X to the nearest integer.
///
/// \code{.cpp}
///    matrix<> X = {{1,2,3},{4,5,6}};
///    matrix<> Y = round(X+X/10);
///    disp(Y);
/// \endcode
///
// \see floor, ceil.
template<typename S>
auto round(matrix<S>const& X)
{
    using T = decltype(std::round(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::round(X(l));}
    return Y;
}

//==========================================================================
// [row]
/// Bi-linear indexing of row elements.
///
/// row(A) return row matrix with row indices of A in bi-linear indexing.
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    matrix<std::size_t> I = row(A);
///    matrix<> B = eval(A(I,0));
///    disp(B);
/// \endcode
///
// \see all, col, get, set, view, cview.
template<typename T>
matrix<std::size_t> row(matrix<T>const& A)
{
    return range(0,size(A,1));
}

//==========================================================================
// [set]
/// Set sub-matrix.
///
/// set(A,L,B) modify A with elements taken from B corresponding to
/// linear indexing L, as A(L) = B.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    set(A,{0,2,4},-ones(1,3));
///    disp(A);
/// \endcode
///
/// set(A,I,J,B) modify A with elements taken from B corresponding to
/// bilinear indexing I for rows and J for columns, as A(I,J) = B.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    set(A,{0,1},{1,2},eye(2));
///    disp(A);
/// \endcode
///
// \see get, all, row, col, view, cview.
template<typename T, typename U>
void set(matrix<T>& A, matrix<std::size_t>const& L, U b)
{
    for (std::size_t l=0; l<numel(L); ++l)
    {
        A(L(l)) = b;
    }
}
template<typename T, typename U>
void set(matrix<T>& A, matrix<std::size_t>const& L, matrix<U>const& B)
{
    if (size(L,1)!=size(B,1) || size(L,2)!=size(B,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions for indices and values must agree.");
    }
    for (std::size_t l=0; l<numel(L); ++l)
    {
        A(L(l)) = B(l);
    }
}
template<typename T, typename U>
void set(matrix<T>& A, matrix<std::size_t>const& I, matrix<std::size_t>const& J, U b)
{
    if (std::min(size(I,1),size(I,2))!=1 || std::min(size(J,1),size(J,2))!=1)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Subscript indices must be vectors.");
    }
    for (std::size_t i=0; i<numel(I); ++i)
    {
        for (std::size_t j=0; j<numel(J); ++j)
        {
            A(I(i),J(j)) = b;
        }
    }
}
template<typename T, typename U>
void set(matrix<T>& A, matrix<std::size_t>const& I, matrix<std::size_t>const& J, matrix<U>const& B)
{
    if (std::min(size(I,1),size(I,2))!=1 || std::min(size(J,1),size(J,2))!=1)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Subscript indices must be vectors.");
    }
    if (numel(I)!=size(B,1) || numel(J)!=size(B,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions for indices and values must agree.");
    }
    for (std::size_t i=0; i<numel(I); ++i)
    {
        for (std::size_t j=0; j<numel(J); ++j)
        {
            A(I(i),J(j)) = B(i,j);
        }
    }
}

//==========================================================================
// [setdiff]
/// Set difference of two arrays.
/// C = setdiff(A,B) for matrix A and B returns a vector with the values
/// in A that are not in B, with no repetitions. C will be sorted.
///
/// \code{.cpp}
///    matrix<> A = {1,2,3};
///    matrix<> B = {2,3,4};
///    matrix<> C = setdiff(A,B);
///    disp(C);
/// \endcode
///
// \see argsetdiff, intersect, unique, union2.
template<typename R, typename S>
auto setdiff(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)+B(0));
    std::vector<T> a=A.val(), b=B.val(), c;
    std::sort(a.begin(),a.end());
    std::sort(b.begin(),b.end());
    std::set_difference(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(c));
    auto last = std::unique(c.begin(),c.end());
    c.erase(last,c.end());
    return matrix<T>(c);
}

//==========================================================================
// [sign]
/// Signum function.
///
/// For each element of X, sign(X) returns 1 if the element
/// is greater than zero, 0 if it equals zero and -1 if it is
/// less than zero.  For the nonzero elements of complex X,
/// sign(X) = X ./ ABS(X).
///
/// \code{.cpp}
///    matrix<> X = {{-1,-2,0},{3,4,0}};
///    matrix<> Y = sign(X);
///    disp(Y);
/// \endcode
///
// \see abs.
template<typename S>
auto sign(matrix<S>const& X)
{
    using T = decltype(X(0)/std::abs(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l)
    {
        if (std::abs(X(l))==0) {Y(l) = 0;}
        else {Y(l) = X(l)/std::abs(X(l));}
    }
    return Y;
}

//==========================================================================
// [sin]
/// Sine of argument in radians.
///
/// sin(X) is the sine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {0,M_PI/2,M_PI,-M_PI/2};
///    matrix<> Y = sin(X);
///    disp(Y);
/// \endcode
///
// \see asin, sind.
template<typename S>
auto sin(matrix<S>const& X)
{
    using T = decltype(std::sin(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::sin(X(l));}
    return Y;
}

//==========================================================================
// [sind]
/// Sine of argument in degrees.
///
/// sind(X) is the sine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {0,90,180,-90};
///    matrix<> Y = sind(X);
///    disp(Y);
/// \endcode
///
// \see asind, sin.
template<typename S>
auto sind(matrix<S>const& X)
{
    using T = decltype(std::sin(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    T r = M_PI/180;
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::sin(r*X(l));}
    return Y;
}

//==========================================================================
// [sinh]
/// Hyperbolic sine.
///
/// sinh(X) is the hyperbolic sine of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {0,1,2,3};
///    matrix<> Y = sinh(X);
///    disp(Y);
/// \endcode
///
// \see asinh.
template<typename S>
auto sinh(matrix<S>const& X)
{
    using T = decltype(std::sinh(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::sinh(X(l));}
    return Y;
}

//==========================================================================
// [size]
/// Size of array.
///
/// S = size(A) for m-by-n matrix A returns the two-element vector [m,n]
/// containing the number of rows and columns in the matrix.
///
/// S = size(A,dim) returns the lengths of the specified dimensions dim.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    disp(size(A));
///    disp(size(A,1));
///    disp(size(A,2));
/// \endcode
///
// \see length, numel.
template<typename T>
matrix<std::size_t> size(matrix<T>const& A)
{
    return {A.size(1),A.size(2)};
}
template<typename T>
std::size_t size(matrix<T>const& A, int dim)
{
    return A.size(dim);
}

//==========================================================================
// [sort]
/// Sort in ascending order.
///
/// B = sort(A) sorts the elements of A in ascending order. The sorted
/// output B has the same type and size as A.
///
/// \code{.cpp}
///    matrix<> A = {{6,5,4},{3,2,1}};
///    matrix<> B = sort(A);
///    disp(B);
/// \endcode
///
/// B = sort(A,DIM) also specifies a dimension DIM to sort along.
/// if DIM==0, M is the sorted element of the array.
/// if DIM==1, M is a vector containing the sorted elements from each column.
/// if DIM==2, M is a vector containing the sorted elements from each row.
///
/// \code{.cpp}
///    matrix<> A = {{6,5,4},{3,2,1}};
///    matrix<> R = sort(A,1);
///    matrix<> C = sort(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see argsort, min, max.
template<typename T>
matrix<T> sort(matrix<T> const& A, int dim=0)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    std::size_t m=size(A,1), n=size(A,2);
    matrix<T> B(m,n);
    std::vector<T> v;
    if (dim==0)
    {
        v.resize(m*n);
        for (std::size_t l=0; l<v.size(); ++l) {v[l]=A(l);}
        std::sort(v.begin(),v.end());
        for (std::size_t l=0; l<v.size(); ++l) {B(l)=v[l];}
    }
    else if (dim==1)
    {
        v.resize(m);
        for (std::size_t j=0; j<n; ++j)
        {
            for (std::size_t i=0; i<m; ++i) {v[i]=A(i,j);}
            std::sort(v.begin(),v.end());
            for (std::size_t i=0; i<m; ++i) {B(i,j)=v[i];}
        }
    }
    else if (dim==2)
    {
        v.resize(n);
        for (std::size_t i=0; i<m; ++i)
        {
            for (std::size_t j=0; j<n; ++j) {v[j]=A(i,j);}
            std::sort(v.begin(),v.end());
            for (std::size_t j=0; j<n; ++j) {B(i,j)=v[j];}
        }
    }
    return B;
}

//==========================================================================
// [sph2cart]
/// Transform spherical to Cartesian coordinates.
///
/// (X,Y,Z) = sph2cart(THE,PHI,RHO) transforms corresponding elements
/// of data stored in spherical coordinates (THE,PHI,RHO) to Cartesian
/// coordinates X,Y,Z. The arrays THE, PHI, and RHO must be the same size.
/// THE and PHI must be in radians.
///
/// THE is the counterclockwise angle in the xy plane measured from the
/// positive x axis. PHI is the elevation angle from the xy plane.
///
/// <b>Remark:</b> This function supports scalar arguments.
///
/// \code{.cpp}
///    matrix<> THE = {0,M_PI/2,M_PI,-M_PI/2};
///    matrix<> PHI = {0,0,0,0};
///    matrix<> RHO = {1,1,1,1};
///    matrix<> X, Y, Z;
///    std::tie(X,Y,Z) = sph2cart(THE,PHI,RHO);
///    disp(X);
///    disp(Y);
///    disp(Z);
/// \endcode
///
// \see cart2sph, cart2pol, pol2cart.
template<typename Q, typename R, typename S>
auto sph2cart(matrix<Q>const& THE, matrix<R>const& PHI, matrix<S>const& RHO)
{
    if (size(THE,1)!=size(PHI,1) || size(THE,2)!=size(PHI,2) ||
        size(THE,1)!=size(RHO,1) || size(THE,2)!=size(RHO,2) )
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(THE(0)+PHI(0)+RHO(0));
    matrix<T> X(size(THE,1),size(THE,2)), Y(size(THE,1),size(THE,2)), Z(size(THE,1),size(THE,2));
    for (std::size_t l=0; l<numel(X); ++l)
    {
        X(l) = RHO(l) * std::cos(THE(l)) * std::cos(PHI(l));
        Y(l) = RHO(l) * std::sin(THE(l)) * std::cos(PHI(l));
        Z(l) = RHO(l) * std::sin(PHI(l));
    }
    return std::make_tuple(X,Y,Z);
}
template<typename Q, typename R, typename S>
inline auto sph2cart(Q const &THE, R const &PHI, S const &RHO)
{
    using T = decltype(THE+PHI+RHO);
    T X = RHO*std::cos(THE)*std::cos(PHI);
    T Y = RHO*std::sin(THE)*std::cos(PHI);
    T Z = RHO*std::sin(PHI);
    return std::make_tuple(X,Y,Z);
}

//==========================================================================
// [sph2idx]
/// Sherical coordinate (in radian) to discretization index.
///
/// IDX = sph2idx(AZM,ELV) give indices of a 1-degree spherical discretization
/// cooresponding to azimut AZM and elevation ELV
///
/// \code{.cpp}
///    matrix<std::size_t> I = range(0,360*181);
///    matrix<> AZM, ELV;
///    std::tie(AZM,ELV) = idx2sph(I);
///    matrix<> J = sph2idx(AZM,ELV);
///    disp(max(J-I));
/// \endcode
///
/// IDX = sph2idx(AZM,ELV,M) use M azimut and elevation.
///
/// IDX = sph2idx(AZM,ELV,M,N) convert azimuth AZM and elevation ELV angles
/// to indices IDX for spherical coordinates of M-by-N size.
///
// \see idx2sph, sphere.
template<typename R, typename S>
matrix<std::size_t> sph2idx(matrix<R>const& azm, matrix<S>const& elv, std::size_t m=0, std::size_t n=0)
{
    if (m==0 && n==0) {m=181; n=360;}
    else if (n==0) {n=m;}
    if (size(azm,1)!=size(elv,1) || size(azm,2)!=size(elv,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    using T = decltype(azm(0)+elv(0));
    matrix<std::size_t> idx(size(azm,1),size(azm,2));
    matrix<T> t,p,r=matrix<T>(size(azm,1),size(azm,2),1);
    T dazm = 2*M_PI/n, delv = M_PI/(m-1), pitwo = M_PI/2, twopi = 2*M_PI;
    long i, j;
    std::tie(t,p,r) = sph2cart(azm,elv,r);
    std::tie(t,p,r) = cart2sph(t,p,r);
    for (std::size_t l=0; l<numel(t); ++l)
    {
        i = std::round(std::fmod(t(l)+twopi,twopi) / dazm);
        j = std::round((p(l)+pitwo) / delv);
        if (i==n) {i=0;}
        idx(l) = j*n + i;
    }
    return idx;
}

//==========================================================================
// [sphere]
/// Generate sphere.
///
/// [X,Y,Z] = sphere() generate Carthesian coordinates of a uniform 1-degree
/// discretisation sphere.
///
/// [X,Y,Z] = sphere(M) generate sphere with M azimuts and elevations.
///
/// [X,Y,Z] = sphere(M,N) generates carthesian coordinates (X,Y,Z) of a unit
/// sphere with M elevations and N azimuts.
///
/// \code{.cpp}
///    matrix<> X,Y,Z;
///    std::tie(X,Y,Z) = sphere(3,4);
///    disp(X);
///    disp(Y);
///    disp(Z);
/// \endcode
///
// \see idx2sph, sph2idx.
template<typename T=double>
auto sphere(std::size_t m=0, std::size_t n=0)
{
    if (m==0 && n==0) {m=181; n=360;}
    else if (n==0) {n=m;}
    matrix<T> azm = linspace<T>(0.,2*M_PI*(1-1./n),n);
    matrix<T> elv = linspace<T>(-M_PI/2,M_PI/2,m);
    azm = mtimes(ones<T>(m,1),azm);
    elv = mtimes(transpose(elv),ones<T>(1,n));
    return sph2cart<T>(azm,elv,ones<T>(m,n));
}

//==========================================================================
// [sphere2]
/// Generate sphere using fibonacci rules.
///
/// [X,Y,Z] = sphere2(N) generate a quasi-uniform unit sphere with
/// N vertices.
///
/// [X,Y,Z] = sphere2(N,rho) generate N-sphere of radius rho.
///
/// \code{.cpp}
///    matrix<> X,Y,Z;
///    std::tie(X,Y,Z) = sphere2(100);
///    disp(X);
///    disp(Y);
///    disp(Z);
/// \endcode
///
// \see sphere.
template<typename T=double>
auto sphere2(std::size_t n, T rho=1)
{
    matrix<T> azm(1,n), elv(1,n);
    T gold = 0.5*(1.+std::sqrt(5));
    for (std::size_t l=0; l<n; ++l)
    {
        azm(l) = fmod(2.*M_PI/gold*l+M_PI,2.*M_PI)-M_PI;
        elv(l) = std::asin(-1.+2./(n-1.)*l);
    }
    return sph2cart<T>(azm,elv,rho*ones<T>(1,n));
}

//==========================================================================
// [sqrt]
/// Square root.
///
/// sqrt(X) is the square root of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {-1,0,1,2,4};
///    auto     Y = sqrt(X);
///    disp(Y);
/// \endcode
///
// \see pow.
template<typename S>
auto sqrt(matrix<S>const& X)
{
    using T = decltype(std::sqrt(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::sqrt(X(l));}
    return Y;
}

//==========================================================================
// [stddev]
/// Standard deviation of elements.
///
/// S = stddev(A) is the standard deviation of all elements of the array A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    double   m = stddev(A);
///    disp(m);
/// \endcode
///
/// S = stddev(A,DIM) is the standard deviation along the dimension DIM.
/// if DIM==0, M is the standard deviation of the array.
/// if DIM==1, M is a vector containing the standard deviation from each column.
/// if DIM==2, M is a vector containing the standard deviation from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = stddev(A,1);
///    matrix<> C = stddev(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see max, min, mean, median, variance.
template<typename T>
matrix<T> stddev(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    return sqrt(variance(A,dim));
}
template<typename T>
inline T stddev(matrix<T>const& A) {return (T)stddev(A,0);}

//==========================================================================
// [sub2ind]
/// Linear index from multiple subscripts.
///
/// sub2ind is used to determine the equivalent single index
/// corresponding to a given set of subscript values.
///
/// L = sub2ind(S,I,J) returns the linear index equivalent to the
/// row and column subscripts in the arrays I and J for a matrix of
/// size S.
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    matrix<std::size_t> I, J, L;
///    std::tie(I,J) = ind2sub(size(A),find(A));
///    L = sub2ind(size(A),I,J);
///    disp(L);
/// \endcode
///
// \see ind2sub, find.
template<typename T=std::size_t>
matrix<std::size_t> sub2ind(matrix<std::size_t> S, matrix<T>const& I, matrix<T>const& J)
{
    if (numel(S)!=2)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Size vector must have 2 elements.");
    }
    if (size(I,1)!=size(J,1)||size(I,2)!=size(J,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"The subscript vectors must be of the same size.");
    }
    matrix<std::size_t> L(size(I,1),size(I,2));
    for (std::size_t l=0; l<numel(L); ++l) {L(l) = I(l)*S(1)+J(l);}
    return L;
}

//==========================================================================
// [sum]
/// Sum of elements.
///
/// S = sum(A) is the sum of all the elements of the array A.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    double   s = sum(A);
///    disp(s);
/// \endcode
///
/// S = sum(A,DIM) sums along the dimension DIM.
/// if DIM==0, M is the sum of all element of the array.
/// if DIM==1, M is a vector containing the sum from each column.
/// if DIM==2, M is a vector containing the sum from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = sum(A,1);
///    matrix<> C = sum(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see prod, diff, cumsum.
template<typename T>
matrix<T> sum(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    matrix<T> B;
    if (dim==0)
    {
        B.resize(1,1,0);
        for (std::size_t l=0; l<numel(A); ++l) {B(0) += A(l);}
    }
    else if (dim==1)
    {
        B.resize(1,size(A,2),0);
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            for (std::size_t j=0; j<size(A,2); ++j) {B(j) += A(i,j);}
        }
    }
    else if (dim==2)
    {
        B.resize(size(A,1),1,0);
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            for (std::size_t j=0; j<size(A,2); ++j) {B(i) += A(i,j);}
        }
    }
    return B;
}
template<typename T>
inline T sum(matrix<T>const& A) {return (T)sum(A,0);}

//==========================================================================
// [tan]
/// Tangent of argument in radians.
///
/// tan(X) is the tangent of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {0,M_PI/2,M_PI,-M_PI/2};
///    matrix<> Y = tan(X);
///    disp(Y);
/// \endcode
///
// \see atan, tand.
template<typename S>
auto tan(matrix<S>const& X)
{
    using T = decltype(std::tan(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::tan(X(l));}
    return Y;
}

//==========================================================================
// [tand]
/// Tangent of argument in degrees.
///
/// tand(X) is the tangent of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {0,90,180,-90};
///    matrix<> Y = tand(X);
///    disp(Y);
/// \endcode
///
// \see atand, tan.
template<typename S>
auto tand(matrix<S>const& X)
{
    using T = decltype(std::tan(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    T r = M_PI/180;
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::tan(r*X(l));}
    return Y;
}

//==========================================================================
// [tanh]
/// Hyperbolic tangent.
///
/// tanh(X) is the hyperbolic tangent of the elements of X.
///
/// \code{.cpp}
///    matrix<> X = {0,1,2,3};
///    matrix<> Y = tanh(X);
///    disp(Y);
/// \endcode
///
// \see atanh.
template<typename S>
auto tanh(matrix<S>const& X)
{
    using T = decltype(std::tanh(X(0)));
    matrix<T> Y(size(X,1),size(X,2));
    for (std::size_t l=0; l<numel(Y); ++l) {Y(l) = std::tanh(X(l));}
    return Y;
}

//==========================================================================
// [tic]
/// Start a stopwatch timer.
///
/// tic() and toc() functions work together to measure elapsed time.
///
/// tic(), by itself, saves the current time that toc() uses later to
/// measure the time elapsed between the two.
/// Use toc(0) to avoid printing.
///
/// \code{.cpp}
///    tic();
///    while (toc(0)<0.314159) {}
///    toc();
/// \endcode
///
// \see toc.
void tic()
{
    ticTimer = std::chrono::high_resolution_clock::now();
}

//==========================================================================
// [tgemm]
/// In-place matrix product.
///
/// tgemm(alpha,A,B,beta,C) performs the in-place matrix-matrix operations
///    C = alpha*A*B + beta*C,
/// where alpha, beta are scalars and A, B, C are matrices with compatible size.
///
/// NOTE : This is a naive implementation, you can use the cblas interface
/// proposed in "linalg.hpp" to get better performance.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> B = eye(3,3);
///    matrix<> C = zeros(2,3);
///    tgemm(1,A,B,1,C);
///    disp(C);
/// \endcode
///
// \see mtimes, xgemm.
template<typename P, typename Q, typename R, typename S, typename T>
void tgemm(P alpha, matrix<Q>const& A, matrix<R>const& B, S beta, matrix<T>& C)
{
    if (size(A,2)!=size(B,1) || size(A,1)!=size(C,1) || size(B,2)!=size(C,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    T a=(T)alpha, b=(T)beta;
    for (std::size_t i=0; i<size(C,1); ++i)
    {
        for (std::size_t j=0; j<size(C,2); ++j)
        {
            C(i,j) = b * ((T)C(i,j));
            for (std::size_t k=0; k<A.size(2); ++k)
            {
                C(i,j) += a * ((T)A(i,k)) * ((T)B(k,j));
            }
        }
    }
}

//==========================================================================
// [transpose]
/// Non-conjugate transpose.
///
/// transpose(A) is called for the syntax A^t.
///
/// \code{.cpp}
///    matrix<> A  = {{1,2,3},{4,5,6}};
///    matrix<> At = transpose(A);
///    disp(At);
/// \endcode
///
// \see reshape.
template<typename T>
matrix<T> transpose(matrix<T>const& A)
{
    matrix<T> B(size(A,2),size(A,1));
    for (std::size_t i=0; i<size(B,1); ++i)
    {
        for (std::size_t j=0; j<size(B,2); ++j)
        {
            B(i,j) = A(j,i);
        }
    }
    return B;
}

//==========================================================================
// [toc]
/// Read the stopwatch timer.
///
/// tic() and toc() functions work together to measure elapsed time.
///
/// toc(), by itself, displays the elapsed time, in seconds, since
/// the most recent execution of the tic() command.
/// t = toc(); saves the elapsed time in t as a double scalar.
/// Use toc(0) to avoid printing.
///
/// \code{.cpp}
///    tic();
///    while (toc(0)<0.314159) {}
///    toc();
/// \endcode
///
// \see tic.
double toc(bool disp)
{
    const auto tocTimer = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time_span = tocTimer - ticTimer;
    if (disp) {std::cout << "Elapsed time is " << time_span.count() << " seconds." << std::endl;}
    return time_span.count();
}

//==========================================================================
// [union2]
/// Set union of two arrays.
/// C = union2(A,B) for matrix A and B, returns a vector with the combined
/// values of the two matrix, with no repetitions. C will be sorted.
///
/// \code{.cpp}
///    matrix<> A = {1,2,3};
///    matrix<> B = {2,3,4};
///    matrix<> C = union2(A,B);
///    disp(C);
/// \endcode
///
// \see intersect, setdiff, unique.
template<typename R, typename S>
auto union2(matrix<R>const& A, matrix<S>const& B)
{
    using T = decltype(A(0)+B(0));
    std::vector<T> a=A.val(), b=B.val(), c;
    std::sort(a.begin(),a.end());
    std::sort(b.begin(),b.end());
    std::set_union(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(c));
    auto last = std::unique(c.begin(),c.end());
    c.erase(last,c.end());
    return matrix<T>(c);
}

//==========================================================================
// [unique]
/// Set unique of an array.
/// B = unique(A) for the matrix A returns a vector with the same values
/// as in A, but with no repetitions. B will be sorted.
///
/// \code{.cpp}
///    matrix<> A = {1,2,1,3,2,3};
///    matrix<> B = unique(A);
///    disp(B);
/// \endcode
///
// \see argunique, intersect, setdiff, union2.
template<typename T>
matrix<T> unique(matrix<T>const& A)
{
    std::vector<T> b = A.val();
    std::sort(b.begin(), b.end());
    auto last = std::unique(b.begin(),b.end());
    b.erase(last, b.end());
    return matrix<T>(b);
}

//==========================================================================
// [values]
/// Copy all values in a matrix.
///
/// V = values(A) returns matrix containing a copy of all values.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> V = values(A);
///    disp(V);
/// \endcode
template<typename T>
matrix<T> values(matrix<T>const& A)
{
    return matrix<T>(A.val());
}

//==========================================================================
// [variance]
/// Variance of elements.
///
/// S = variance(X) is the variance of all elements of the array X.
/// If N is the sample size, variance normalizes by N and produces
/// the second moment of the sample about its mean.
/// 
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    double   v = variance(A);
///    disp(v);
/// \endcode
///
/// S = variance(X,DIM) is the variance along the dimension DIM:
/// if DIM==0, M is the variance of the array.
/// if DIM==1, M is a vector containing the variance from each column.
/// if DIM==2, M is a vector containing the variance from each row.
///
/// \code{.cpp}
///    matrix<> A = {{1,2,3},{4,5,6}};
///    matrix<> R = variance(A,1);
///    matrix<> C = variance(A,2);
///    disp(R);
///    disp(C);
/// \endcode
///
// \see max, min, median, mean, stddev.
template<typename T>
matrix<T> variance(matrix<T>const& A, int dim)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    return 1./size(A,dim)*sum(pow(A,2),dim) - pow(mean(A,dim),2);
}
template<typename T>
inline T variance(matrix<T>const& A) {return (T)variance(A,0);}

//==========================================================================
// [writebin]
/// Write matrix in *.bin file.
///
/// writebin("path","filename") write binary file with matrix stored
/// in the following order :
/// m n M(0,0) M(0,1) ... M(i,j) ... M(m-1,n-2) M(m-1,n-1)
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    writebin("./","matrix.bin",A);
///    matrix<> B = readbin("./","matrix.bin");
///    disp(B);
/// \endcode
///
// \see readbin, writetxt.
template<typename T>
void writebin(std::string path, std::string filename, matrix<T>const& A)
{
    std::ofstream file(path+filename, std::ios::out | std::ios::binary);
    if (file)
    {
        std::size_t m=size(A,1), n=size(A,2);
        file.write((char*)&m,sizeof(std::size_t));
        file.write((char*)&n,sizeof(std::size_t));
        for (std::size_t l=0; l<numel(A); ++l)
        {
            file.write((char*)&A(l),sizeof(T));
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"File "+path+filename+" not found.");
    }
    file.close();
}

//==========================================================================
// [writetxt]
/// Write matrix in *.txt file.
///
/// writetxt("path","filename") write ascii text file with matrix stored
/// in the following format :
/// m       n
/// M(0,0)    M(0,1)    ...  M(0,n-2)    M(0,n-1)
/// ...                                  ...
/// M(m-1,0)  M(m-1,1)  ...  M(m-1,n-2)  M(m-1,n-1)
///
/// \code{.cpp}
///    matrix<> A = eye(3,4);
///    writetxt("./","matrix.txt",A);
///    matrix<> B = readtxt("./","matrix.txt");
///    disp(B);
/// \endcode
///
// \see readtxt, writebin.
template<typename T>
void writetxt(std::string path, std::string filename, matrix<T>const& A)
{
    std::ofstream file(path+filename);
    if (file)
    {
        file << size(A,1) << "  " << size(A,2) << std::endl;
        for (std::size_t i=0; i<size(A,1); ++i)
        {
            for (std::size_t j=0; j<size(A,2); ++j)
            {
                file << A(i,j) << "  ";
            }
            file << std::endl;
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"File "+path+filename+" not found.");
    }
    file.close();
}

//==========================================================================
// [vertcat]
/// Vertical concatenation.
///
/// vertcat(A,B) is the vertical concatenation of array A and B.
/// A and B must have the same number of columns.
/// vertcat(A,B) is the same as cat(1,A,B).
///
/// \code{.cpp}
///    matrix<> A = {1,2,3};
///    matrix<> B = {4,5,6};
///    matrix<> C = vertcat(A,B);
///    disp(C);
/// \endcode
///
// \see horzcat, cat.
template<typename R, typename S>
inline auto vertcat(R const& A, S const& B) {return cat(1,A,B);}

//==========================================================================
// [warning]
/// Display warning message in command window.
///
/// warning(__FILE__, __LINE__, __FUNCTION__,"message") displays a descriptive
/// message when the currently-running program encounters a warning condition.
///
/// \code{.cpp}
///    warning(__FILE__, __LINE__, __FUNCTION__,"This is a warning message.");
/// \endcode
///
// \see error, disp, help.
inline void warning(std::string file, int line, std::string function, std::string comment)
{
    std::cout << std::endl;
    std::cout << "WARNING! In " << file << " at line " << line << " with function '";
    std::cout << function << "':" <<std::endl;
    std::cout << comment << std::endl << std::endl;
}

//==========================================================================
// [zeros]
/// Zeros matrix.
///
/// zeros(N) is an N-by-N matrix of zeros.
///
/// zeros(M,N) or zeros({M,N}) is an M-by-N matrix of zeros.
///
/// \code{.cpp}
///    matrix<> A = zeros(2);
///    matrix<> B = zeros(2,3);
///    matrix<> C = zeros(size(B));
///    disp(A);
///    disp(B);
///    disp(C);
/// \endcode
///
// \see ones, eye, rand.
template<typename T=double>
matrix<T> zeros(std::size_t m, long n=-1)
{
    if (n==-1) {n=m;}
    return matrix<T>(m,n,0);
}
template<typename T=double>
matrix<T> zeros(matrix<std::size_t>const& S)
{
    return matrix<T>(S(0),S(1),0);
}

}
