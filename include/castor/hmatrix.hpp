/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : hmatrix.hpp                                   |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Hierarchical matrix operations & algebra with |
 |  `---'  |                matlab-like standard functions and procedures |
 +========================================================================+
 */

#pragma once
#define CASTOR_HMATRIX_HPP
#include <castor/matrix.hpp>
#include <castor/linalg.hpp>
#include <castor/smatrix.hpp>

namespace castor
{

//==========================================================================//
//                         BINARY TREE CLASS                                //
//==========================================================================//
// [bintree]
///
template<typename T>
class bintree
{
public:
    // CONSTRUCTOR
    bintree(){};
    bintree(matrix<T>const& X, std::size_t n=100);
    
    // FUNCTIONS
    matrix<T> leaf() const;
    
    // ACCESSORS
    bool const&               isleaf()   const {return m_isleaf;};
    matrix<T>const&           crd()      const {return m_crd;};
    matrix<T>const&           ctr()      const {return m_ctr;};
    matrix<T>const&           edg()      const {return m_edg;};
    matrix<std::size_t>const& ind(int i) const {return m_ind[i];};
    bintree<T>const&          sub(int i) const {return m_sub[i];};
    
private:
    bool                    m_isleaf; // TRUE IF LEAF
    matrix<T>               m_crd;    // POINTS COORDINATES
    matrix<T>               m_min;    // MINIMUM
    matrix<T>               m_max;    // MAXIMUM
    matrix<T>               m_ctr;    // CENTER
    matrix<T>               m_edg;    // EDGE LENGTH
    matrix<std::size_t>     m_ind[2]; // CHILDREN INDICES
    std::vector<bintree<T>> m_sub;    // CHILDREN (RECURSION)
};

//==========================================================================
//
///
template<typename T>
bintree<T>::bintree(matrix<T>const& X, std::size_t n)
{
    // Initialize
    std::size_t N = size(X,1);
    
    // Box definition
    m_isleaf = (N<n);
    m_crd = X;
    m_min = min(X,1);
    m_max = max(X,1);
    m_ctr = (m_min+m_max)/2.;
    m_edg = m_max-m_min;
    
    // Children recursion
    if (!m_isleaf)
    {
        std::size_t         dim  = argmax(m_edg);
        matrix<std::size_t> iloc = argsort(eval(X(range(0,N),dim)));
        m_ind[0] = eval(iloc(range(0,N/2)));
        m_ind[1] = eval(iloc(range(N/2,N)));
        m_sub.push_back(bintree(eval(X(m_ind[0],col(X))),n));
        m_sub.push_back(bintree(eval(X(m_ind[1],col(X))),n));
    }
}

//==========================================================================
//
///
template<typename T>
matrix<T> bintree<T>::leaf() const
{
    matrix<T> v(1,size(m_crd,1));
    if (m_sub.size()==2)
    {
        v(m_ind[0]) = m_sub[0].leaf();
        v(m_ind[1]) = m_sub[1].leaf();
    }
    else
    {
        v = rand(1)*ones(size(v));
    }
    return v;
}


//==========================================================================//
//                           H-MATRIX CLASS                                 //
//==========================================================================//
// [hmatrix]
///
template<typename T>
class hmatrix
{
public:
    // CONSTRUCTOR
    hmatrix(){};
    hmatrix(std::size_t m, std::size_t n, double tol, T v=0);
    template<typename S>
    hmatrix(bintree<S>const& X, bintree<S>const& Y, double tol,
            matrix<std::size_t> I={}, matrix<std::size_t> J={});
    template<typename S>
    hmatrix(bintree<S>const& X, bintree<S>const& Y, double tol, smatrix<T>const& Ms);
    template<typename S>
    hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol);
    template<typename S>
    hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol, matrix<T>const& M);
    template<typename S>
    hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol, smatrix<T>const& Ms);
    template<typename S>
    hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol,
            std::function<matrix<T>(matrix<std::size_t>,matrix<std::size_t>)>const& fct);
    template<typename S>
    hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol,
            std::function<matrix<T>(matrix<std::size_t>,matrix<std::size_t>)>const& fct,
            matrix<T>const& V, matrix<T>& MV);
    
    // FUNCTIONS
    void check();
    void full2lowrank();
    void fusion();
    void hfull(matrix<T>& M, matrix<std::size_t> I={}, matrix<std::size_t> J={}) const;
    void hinv();
    void hllowsolve(hmatrix<T>const& Lh);
    void hllowsolve(matrix<T>& B, matrix<std::size_t> I={}) const;
    void hlmtimes(matrix<T>const& B, matrix<T>& C, matrix<std::size_t> I={}, matrix<std::size_t> J={}) const;
    void hlu(hmatrix<T>& U);
    void hlupsolve(matrix<T>& B, matrix<std::size_t> I={}) const;
    void hrmtimes(matrix<T>const& A, matrix<T>& C, matrix<std::size_t> I={}, matrix<std::size_t> J={}) const;
    void hrupsolve(hmatrix<T>const& Uh);
    void hrupsolve(matrix<T>& B, matrix<std::size_t> J={}) const;
    void htranspose();
    void leafptr(std::vector<hmatrix<T>*>& ptr,
                 std::vector<matrix<std::size_t>>& idx, std::vector<matrix<std::size_t>>& jdx,
                 matrix<std::size_t> I={}, matrix<std::size_t> J={});
    void recompress(double tol);
    void spy(matrix<T>& M, matrix<std::size_t> I={}, matrix<std::size_t> J={}) const;
    void stat(matrix<std::size_t>& info) const;
    void tgeabhm(T alpha, matrix<T>const& A, matrix<T>const& B, T beta, matrix<std::size_t> I={}, matrix<std::size_t> J={});
    void tgehmhm(T alpha, hmatrix<T>const& Ah, hmatrix<T>const& Bh, T beta);
    
    // OPERATORS
    hmatrix<T>& operator+=(T b);
    hmatrix<T>& operator*=(T b);
    hmatrix<T>& operator+=(hmatrix<T>const& Bh);
    
    // ACCESSORS
    hmatrix<T>const&          getChildren(int i) const {return m_chd[i];};
    matrix<std::size_t>const& getColumn(int i)   const {return m_col[i];};
    matrix<T>const&           getData(int i)     const {return m_dat[i];};
    matrix<std::size_t>const& getRow(int i)      const {return m_row[i];};
    std::size_t               getSize(int i)     const {return m_siz(i);};
    double                    getTolerance()     const {return m_tol;};
    int                       getType()          const {return m_typ;};
    std::vector<hmatrix<T>>&  setChildren()      {return m_chd;};
    hmatrix<T>&               setChildren(int i) {return m_chd[i];};
    matrix<std::size_t>&      setColumn(int i)   {return m_col[i];};
    matrix<T>&                setData(int i)     {return m_dat[i];};
    matrix<std::size_t>&      setRow(int i)      {return m_row[i];};
    matrix<std::size_t>&      setSize()          {return m_siz;};
    double&                   setTolerance()     {return m_tol;};
    int&                      setType()          {return m_typ;};
    
private:
    int                      m_typ=-1;    // LEAF TYPE (-1:UNDEF; 0:RECURSIVE; 1:COMPRESSED; 2:FULL)
    matrix<std::size_t>      m_siz={0,0}; // LEAF SIZE
    double                   m_tol;       // ACCURACY OF COMPRESSION
    matrix<T>                m_dat[2];    // MATRIX DATA (M or A*B)
    matrix<std::size_t>      m_row[4];    // CHILDREN ROWS INDICES (M11,M12,M21,M22)
    matrix<std::size_t>      m_col[4];    // CHILDREN COLUMNS INDICES
    std::vector<hmatrix<T>>  m_chd;       // CHILDREN (RECURSION)
};

//==========================================================================//
//                              CONSTRUCTORS                                //
//==========================================================================//
//==========================================================================
// [.hmatrix]
///
template<typename T>
hmatrix<T>::hmatrix(std::size_t m, std::size_t n, double tol, T v)
{
    m_typ    = 1;
    m_siz    = {m,n};
    m_tol    = tol;
    m_dat[0] = matrix<T>(m,1,v);
    m_dat[1] = matrix<T>(1,n,1.);
}

//==========================================================================
// [.hmatrix]
///
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(bintree<S>const& X, bintree<S>const& Y, double tol,
                    matrix<std::size_t> I, matrix<std::size_t> J)
{
    // Global Indices
    if (isempty(I)) {I = range(0,size(X.crd(),1));}
    if (isempty(J)) {J = range(0,size(Y.crd(),1));}
    m_siz = {numel(I),numel(J)};
    m_tol = tol;
    
    // Far test
    matrix<S> XYctr = abs(X.ctr()-Y.ctr());
    matrix<S> XYedg = X.edg()+Y.edg();
    bool isfar = max( (XYctr>=0.75*XYedg) && (XYctr>=0.75*mean(XYedg)) );
    
    // Far leaf
    if (isfar)
    {
        m_typ = 1;
    }
    // Full leaf
    else if (X.isleaf() || Y.isleaf())
    {
        m_typ = 2;
    }
    // Recursion
    else
    {
        m_typ = 0;
        m_chd.resize(4,hmatrix<T>());
        for (int h=0; h<4; ++h)
        {
            m_row[h] = X.ind(h/2);
            m_col[h] = Y.ind(h%2);
            m_chd[h] = hmatrix<T>(X.sub(h/2),Y.sub(h%2),tol,eval(I(m_row[h])),eval(J(m_col[h])));
        }
    }
}

//==========================================================================
// [.hmatrix]
///
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(bintree<S>const& X, bintree<S>const& Y, double tol, smatrix<T>const& Ms)
{
    // Global data
    m_siz = size(Ms);
    m_tol = tol;
    
    // Empty leaf
    if (nnz(Ms)==0)
    {
        m_typ = 1;
        m_dat[0] = zeros<T>(m_siz(0),1);
        m_dat[1] = zeros<T>(1,m_siz(1));
    }
    // Full leaf
    else if (X.isleaf() || Y.isleaf())
    {
        m_typ = 2;
        m_dat[0] = full(Ms);
    }
    // Recursion
    else
    {
        m_typ = 0;
        m_chd.resize(4,hmatrix<T>());
        for (int h=0; h<4; ++h)
        {
            m_row[h] = X.ind(h/2);
            m_col[h] = Y.ind(h%2);
            m_chd[h] = hmatrix<T>(X.sub(h/2),Y.sub(h%2),tol,eval(Ms(m_row[h],m_col[h])));
        }
    }
}

//==========================================================================
// [.hmatrix]
///
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol)
{
    bintree<S> Xtree(X);
    bintree<S> Ytree(Y);
    (*this) = hmatrix<T>(Xtree,Ytree, tol);
}

//==========================================================================
// [.hmatrix]
///
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol, matrix<T>const& M)
{
    auto fct = [&M](matrix<std::size_t> Ix,matrix<std::size_t> Iy)
    {
        return eval(M(Ix,Iy));
    };
    (*this) = hmatrix(X,Y,tol,fct);
}

//==========================================================================
// [.hmatrix]
///
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol, smatrix<T>const& Ms)
{
    bintree<S> Xtree(X);
    bintree<S> Ytree(Y);
    (*this) = hmatrix<T>(Xtree,Ytree, tol, Ms);
}

//==========================================================================
// [.hmatrix]
///
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol,
                    std::function<matrix<T>(matrix<std::size_t>,matrix<std::size_t>)>const& fct)
{
    // Build tree and block interactions
    (*this) = hmatrix<T>(X,Y,tol);
    
    // Get leaves references
    std::vector<hmatrix<T>*> ptr;
    std::vector<matrix<std::size_t>> idx;
    std::vector<matrix<std::size_t>> jdx;
    leafptr(ptr,idx,jdx);
    bool flag;
    
    // Build leaf data with partial pivoting
#pragma omp parallel for
    for (std::size_t l=0; l<ptr.size(); ++l)
    {
        if (ptr[l]->m_siz(0)!=numel(idx[l]) || ptr[l]->m_siz(1)!=numel(jdx[l]))
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
        }
        if (ptr[l]->m_typ==1 && tol>0)
        {
            std::tie(ptr[l]->m_dat[0],ptr[l]->m_dat[1],flag) = aca(idx[l],jdx[l],fct,tol);
            if (!flag)
            {
                error(__FILE__, __LINE__, __FUNCTION__,"ACA compression failed.");
            }
        }
        else if (ptr[l]->m_typ==2 || tol<=0)
        {
            ptr[l]->m_dat[0] = fct(idx[l],jdx[l]);
            ptr[l]->m_typ = 2;
            ptr[l]->full2lowrank();
        }
        else
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
        }
    }
    
    // Check and fusion
    check();
}

//==========================================================================
// [.hmatrix]
///
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol,
                    std::function<matrix<T>(matrix<std::size_t>,matrix<std::size_t>)>const& fct,
                    matrix<T>const& V, matrix<T>& MV)
{
    if (size(V,1)!=size(Y,1) || size(MV,1)!=size(X,1) || size(V,2)!=size(MV,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    
    // Build tree and block interactions
    (*this) = hmatrix<T>(X,Y,tol);
    
    // Get leaves references
    std::vector<hmatrix<T>*> ptr;
    std::vector<matrix<std::size_t>> idx;
    std::vector<matrix<std::size_t>> jdx;
    leafptr(ptr,idx,jdx);
    
    // Temporary data
    matrix<std::size_t> C = col(V);
    matrix<T> A, B, M;
    bool flag;
    
    // Build leaf data with partial pivoting
#pragma omp parallel for
    for (std::size_t l=0; l<ptr.size(); ++l)
    {
        if (ptr[l]->m_siz(0)!=numel(idx[l]) || ptr[l]->m_siz(1)!=numel(jdx[l]))
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
        }
        if (ptr[l]->m_typ==1 && tol>0)
        {
            std::tie(A,B,flag) = aca(idx[l],jdx[l],fct,tol);
            if (!flag)
            {
                error(__FILE__, __LINE__, __FUNCTION__,"ACA compression failed.");
            }
            MV(idx[l],C) = eval(MV(idx[l],C)) + mtimes(A,mtimes(B,eval(V(jdx[l],C))));
        }
        else if (ptr[l]->m_typ==2 || tol<=0)
        {
            M = fct(idx[l],jdx[l]);
            MV(idx[l],C) = eval(MV(idx[l],C)) + mtimes(M,eval(V(jdx[l],C)));
        }
        else
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
        }
    }
}


//==========================================================================//
//                               FUNCTIONS                                  //
//==========================================================================//
//==========================================================================
// [.check]
///
template<typename T>
void hmatrix<T>::check()
{
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].check();
        }
        fusion();
    }
    else if (m_typ==1)
    {
        if (size(m_dat[0],1)!=m_siz(0) || size(m_dat[1],2)!=m_siz(1))
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Compressed leaf is empty.");
        }
        if (size(m_dat[0],2)!=size(m_dat[1],1))
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Compressed leaf has uncompatible rank.");
        }
    }
    else if (m_typ==2)
    {
        if (size(m_dat[0],1)!=m_siz(0) || size(m_dat[0],2)!=m_siz(1))
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Full leaf is empty.");
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.full2lowrank]
///
template<typename T>
void hmatrix<T>::full2lowrank()
{
    if (m_typ==2 && m_tol>0)
    {
        std::size_t rk = rank(m_dat[0]);
        if (rk<min(size(m_dat[0]))-1)
        {
            matrix<T> A, B;
            bool flag;
            std::tie(A,B,flag) = aca(m_dat[0],m_tol,min(size(m_dat[0]))/4);
            if (flag)
            {
                m_dat[0] = A;
                m_dat[1] = B;
                m_typ    = 1;
            }
        }
    }
}

//==========================================================================
// [.fusion]
///
template<typename T>
void hmatrix<T>::fusion()
{
    if (m_typ==0)
    {
        // Leaves
        matrix<int> leaf(1,4);
        for (int h=0; h<4; ++h)
        {
            leaf(h) = m_chd[h].m_typ;
        }
        
        // Low-rank fusion
        if (nnz(leaf==1)==4)
        {
            std::size_t l=0, n=0;
            matrix<std::size_t> r(1,4);
            for (int h=0; h<4; ++h)
            {
                r(h) = size(m_chd[h].m_dat[0],2);
                n += numel(m_chd[h].m_dat[0]) + numel(m_chd[h].m_dat[1]);
            }
            matrix<T> A(m_siz(0),sum(r)), B(sum(r),m_siz(1));
            for (int h=0; h<4; ++h)
            {
                A(m_row[h],range(l,l+r(h))) = m_chd[h].m_dat[0];
                B(range(l,l+r(h)),m_col[h]) = m_chd[h].m_dat[1];
                l += r(h);
            }
            bool flag;
            std::size_t rmax = (n/(m_siz(0)+m_siz(1)))+1;
            std::tie(A,B,flag) = aca(A,B,m_tol,rmax);
            if (flag)
            {
                m_typ    = 1;
                m_dat[0] = A;
                m_dat[1] = B;
                for (int h=0; h<4; ++h)
                {
                    m_row[h].resize(0,0);
                    m_col[h].resize(0,0);
                }
                m_chd.clear();
            }
        }
    }
}

//==========================================================================
// [.hfull]
///
template<typename T>
void hmatrix<T>::hfull(matrix<T>& M, matrix<std::size_t> I, matrix<std::size_t> J) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].hfull(M,eval(I(m_row[h])),eval(J(m_col[h])));
        }
    }
    else if (m_typ==1) {M(I,J) = mtimes(m_dat[0],m_dat[1]);}
    else if (m_typ==2) {M(I,J) = m_dat[0];}
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.hinv]
///
template<typename T>
void hmatrix<T>::hinv()
{
    if (m_siz(0)!=m_siz(1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    if (m_typ==0)
    {
        // Am1 -> M11
        m_chd[0].hinv();

        // Usefull data
        hmatrix<T> Am1B = mtimes(m_chd[0],m_chd[1]);
        hmatrix<T> CAm1 = mtimes(m_chd[2],m_chd[0]);

        // S = D - C * Am1 * B -> M22
        m_chd[3].tgehmhm((T)-1,m_chd[2],Am1B,(T)1);

        // Sm1 -> M22
        m_chd[3].hinv();

        // - Am1 * B * Sm1 -> M12
        m_chd[1].tgehmhm((T)-1,Am1B,m_chd[3],(T)0);

        // Am1 + Am1 * B * Sm1 * C * Am1 -> M11
        m_chd[0].tgehmhm((T)-1,m_chd[1],CAm1,(T)1);

        // - Sm1 * C * Am1 -> M21
        m_chd[2].tgehmhm((T)-1,m_chd[3],CAm1,(T)0);

        // Leaf fusion
        fusion();
    }
    else if (m_typ==2)
    {
        m_dat[0] = inv(m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.hllowsolve]
///
template<typename T>
void hmatrix<T>::hllowsolve(hmatrix<T>const& Lh)
{
    if (m_siz(0)!=Lh.getSize(0))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    // hmatrix \ hmatrix -> hmatrix
    if (Lh.getType()==0 && m_typ==0)
    {
        // X11 -> L11 \ B11
        m_chd[0].hllowsolve(Lh.getChildren(0));

        // X12 -> L11 \ B12
        m_chd[1].hllowsolve(Lh.getChildren(0));

        // X21 -> L22 \ (B21 - L21*X11)
        m_chd[2].tgehmhm((T)-1,Lh.getChildren(2),m_chd[0],(T)1);
        m_chd[2].hllowsolve(Lh.getChildren(3));

        // X22 -> L22 \ (B22 - L21 * X12)
        m_chd[3].tgehmhm((T)-1,Lh.getChildren(2),m_chd[1],(T)1);
        m_chd[3].hllowsolve(Lh.getChildren(3));

        // Fusion
        fusion();
    }
    // hmatrix \ compr or full -> compr or full
    else if (Lh.getType()==0 && m_typ>0)
    {
        Lh.hllowsolve(m_dat[0]);
    }
    // full \ compr or full -> compr or full
    else if (Lh.getType()==2 && m_typ>0)
    {
        m_dat[0] = linsolve(Lh.getData(0),m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.hllowsolve]
///
template<typename T>
void hmatrix<T>::hllowsolve(matrix<T>& B, matrix<std::size_t> I) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(I))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    matrix<std::size_t> J = col(B);
    if (m_typ==0)
    {
        // X1 -> L11 \ B1
        matrix<std::size_t> I1 = eval(I(m_row[0]));
        m_chd[0].hllowsolve(B,I1);

        // X2 -> L22 \ (B2 - L21*X1)
        matrix<std::size_t> I2 = eval(I(m_row[2]));
        B(I2,J) = eval(B(I2,J)) - mtimes(m_chd[2],eval(B(I1,J)));
        m_chd[3].hllowsolve(B,I2);
    }
    else if (m_typ==2)
    {
        B(I,J) = linsolve(m_dat[0],eval(B(I,J)));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.hlupsolve]
///
template<typename T>
void hmatrix<T>::hlupsolve(matrix<T>& B, matrix<std::size_t> I) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(I))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    matrix<std::size_t> J = col(B);
    if (m_typ==0)
    {
        // X2 -> U22 \ B2
        matrix<std::size_t> I2 = eval(I(m_row[2]));
        m_chd[3].hlupsolve(B,I2);

        // X1 -> U11 \ (B1 - U12*X2)
        matrix<std::size_t> I1 = eval(I(m_row[0]));
        B(I1,J) = eval(B(I1,J)) - mtimes(m_chd[1],eval(B(I2,J)));
        m_chd[0].hlupsolve(B,I1);
    }
    else if (m_typ==2)
    {
        B(I,J) = linsolve(m_dat[0],eval(B(I,J)));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.hlmtimes]
///
template<typename T>
void hmatrix<T>::hlmtimes(matrix<T>const& B, matrix<T>& C, matrix<std::size_t> I, matrix<std::size_t> J) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].hlmtimes(B,C,eval(I(m_row[h])),eval(J(m_col[h])));
        }
    }
    else if (m_typ==1)
    {
        matrix<std::size_t> c = col(B);
        C(I,c) = eval(C(I,c)) + mtimes(m_dat[0],mtimes(m_dat[1],eval(B(J,c))));
    }
    else if (m_typ==2)
    {
        matrix<std::size_t> c = col(B);
        C(I,c) = eval(C(I,c)) + mtimes(m_dat[0],eval(B(J,c)));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.hlu]
///
template<typename T>
void hmatrix<T>::hlu(hmatrix<T>& Uh)
{
    if (m_siz(0)!=m_siz(1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    Uh.setTolerance() = m_tol;
    Uh.setSize()      = m_siz;
    if (m_typ==0)
    {
        // Initialize Uh
        Uh.setType() = 0;
        Uh.setChildren().resize(4,hmatrix<T>());
        for (int h=0; h<4; ++h)
        {
            Uh.setRow(h)    = m_row[h];
            Uh.setColumn(h) = m_col[h];
        }
        
        // [L11,U11] -> M11
        m_chd[0].hlu(Uh.setChildren(0));

        // U12 -> L11 \ M12   &&   L12 -> 0
        Uh.setChildren(1) = std::move(m_chd[1]);
        Uh.setChildren(1).hllowsolve(m_chd[0]);
        m_chd[1] = hmatrix<T>(numel(m_row[1]),numel(m_col[1]),m_tol);
        
        // L21 -> M21 / U11   &&   U21 -> 0
        m_chd[2].hrupsolve(Uh.getChildren(0));
        Uh.setChildren(2) = hmatrix<T>(numel(m_row[2]),numel(m_col[2]),m_tol);

        // M22 -> M22 - L21*U12
        m_chd[3].tgehmhm((T)-1,m_chd[2],Uh.getChildren(1),(T)1);

        // [L22,U22] -> M22
        m_chd[3].hlu(Uh.setChildren(3));

        // Fusion
        fusion();
        Uh.fusion();
    }
    else if (m_typ==2)
    {
        Uh.setType() = 2;
        std::tie(m_dat[0],Uh.setData(0)) = lu(m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.hrmtimes]
///
template<typename T>
void hmatrix<T>::hrmtimes(matrix<T>const& A, matrix<T>& C, matrix<std::size_t> I, matrix<std::size_t> J) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].hrmtimes(A,C,eval(I(m_row[h])),eval(J(m_col[h])));
        }
    }
    else if (m_typ==1)
    {
        matrix<std::size_t> r = row(A);
        C(r,J) = eval(C(r,J)) + mtimes(mtimes(eval(A(r,I)),m_dat[0]),m_dat[1]);
    }
    else if (m_typ==2)
    {
        matrix<std::size_t> r = row(A);
        C(r,J) = eval(C(r,J)) + mtimes(eval(A(r,I)),m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.hrupsolve]
///
template<typename T>
void hmatrix<T>::hrupsolve(hmatrix<T>const& Uh)
{
    if (m_siz(1)!=Uh.getSize(1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    // hmatrix \ hmatrix -> hmatrix
    if (Uh.getType()==0 && m_typ==0)
    {
        // X11 -> B11 / U11
        m_chd[0].hrupsolve(Uh.getChildren(0));

        // X21 -> B21 / U11
        m_chd[2].hrupsolve(Uh.getChildren(0));
        
        // X12 -> (B12 - X11*U12) / U22
        m_chd[1].tgehmhm((T)-1,m_chd[0],Uh.getChildren(1),(T)1);
        m_chd[1].hrupsolve(Uh.getChildren(3));

        // X22 -> (B22 - X21*U12) / U22
        m_chd[3].tgehmhm((T)-1,m_chd[2],Uh.getChildren(1),(T)1);
        m_chd[3].hrupsolve(Uh.getChildren(3));

        // Fusion
        fusion();
    }
    // compr / hmatrix -> compr
    else if (Uh.getType()==0 && m_typ==1)
    {
        Uh.hrupsolve(m_dat[1]);
    }
    // full / hmatrix -> full
    else if (Uh.getType()==0 && m_typ==2)
    {
        Uh.hrupsolve(m_dat[0]);
    }
    // compr / full -> compr
    else if (Uh.getType()==2 && m_typ==1)
    {
        m_dat[1] = linsolve(transpose(Uh.getData(0)),transpose(m_dat[1]));
        m_dat[1] = transpose(m_dat[1]);
    }
    else if (Uh.getType()==2 && m_typ==2)
    {
        m_dat[0] = linsolve(transpose(Uh.getData(0)),transpose(m_dat[0]));
        m_dat[0] = transpose(m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.hrupsolve]
///
template<typename T>
void hmatrix<T>::hrupsolve(matrix<T>& B, matrix<std::size_t> J) const
{
    if (isempty(J)) {J = range(0,m_siz(0));}
    if (m_siz(0)!=numel(J) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    matrix<std::size_t> I = row(B);
    if (m_typ==0)
    {
        // X1 -> B1 / U11
        matrix<std::size_t> J1 = eval(J(m_col[0]));
        m_chd[0].hrupsolve(B,J1);

        // X2 -> (B2 - X1*U12) / U22
        matrix<std::size_t> J2 = eval(J(m_col[1]));
        B(I,J2) = eval(B(I,J2)) - mtimes(eval(B(I,J1)),m_chd[1]);
        m_chd[3].hrupsolve(B,J2);
    }
    else if (m_typ==2)
    {
        B(I,J) = transpose(linsolve(transpose(m_dat[0]),transpose(eval(B(I,J)))));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.htranspose]
///
template<typename T>
void hmatrix<T>::htranspose()
{
    m_siz = {m_siz(1),m_siz(0)};
    if (m_typ==0)
    {
        auto Mh  = std::move(m_chd[1]);
        m_chd[1] = std::move(m_chd[2]);
        m_chd[2] = std::move(Mh);
        for (int h=0; h<4; ++h)
        {
            m_chd[h].htranspose();
        }
        matrix<std::size_t> row[4] = {m_row[0],m_row[1],m_row[2],m_row[3]};
        m_row[0] = m_col[0];  m_row[1] = m_col[2];
        m_row[2] = m_col[1];  m_row[3] = m_col[3];
        m_col[0] = row[0];  m_col[1] = row[2];
        m_col[2] = row[1];  m_col[3] = row[3];
    }
    else  if (m_typ==1)
    {
        matrix<T> At = transpose(m_dat[0]);
        m_dat[0] = transpose(m_dat[1]);
        m_dat[1] = At;
    }
    else  if (m_typ==2)
    {
        m_dat[0] = transpose(m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.leafptr]
///
template<typename T>
void hmatrix<T>::leafptr(std::vector<hmatrix<T>*>& ptr,
                         std::vector<matrix<std::size_t>>& idx, std::vector<matrix<std::size_t>>& jdx,
                         matrix<std::size_t> I, matrix<std::size_t> J)
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].leafptr(ptr,idx,jdx,eval(I(m_row[h])),eval(J(m_col[h])));
        }
    }
    else
    {
        ptr.push_back(this);
        idx.push_back(I);
        jdx.push_back(J);
    }
}

//==========================================================================
// [.recompress]
///
template<typename T>
void hmatrix<T>::recompress(double tol)
{
    m_tol = std::max(tol,m_tol);
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].recompress(tol);
        }
        fusion();
    }
    else if (m_typ==1)
    {
        std::tie(m_dat[0],m_dat[1]) = qrsvd(m_dat[0],m_dat[1],m_tol);
    }
    else if (m_typ==2)
    {
        full2lowrank();
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.spy]
///
template<typename T>
void hmatrix<T>::spy(matrix<T>& M, matrix<std::size_t> I, matrix<std::size_t> J) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        matrix<std::size_t> I1=range(0,numel(I)/2), I2=range(numel(I)/2,numel(I));
        matrix<std::size_t> J1=range(0,numel(J)/2), J2=range(numel(J)/2,numel(J));
        m_chd[0].spy(M,eval(I(I1)),eval(J(J1)));
        m_chd[1].spy(M,eval(I(I1)),eval(J(J2)));
        m_chd[2].spy(M,eval(I(I2)),eval(J(J1)));
        m_chd[3].spy(M,eval(I(I2)),eval(J(J2)));
    }
    else
    {
        M(I,J) = m_typ;
    }
}

//==========================================================================
// [.stat]
///
template<typename T>
void hmatrix<T>::stat(matrix<std::size_t>& info) const
{
    if (m_typ==-1) {return;}
    info(m_typ) += 1;
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].stat(info);
        }
    }
    else if (m_typ==1) {info(3) += numel(m_dat[0])+numel(m_dat[1]);}
    else if (m_typ==2) {info(3) += numel(m_dat[0]);}
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.tgeabhm]
///
template<typename T>
void hmatrix<T>::tgeabhm(T alpha, matrix<T>const& A, matrix<T>const& B, T beta, matrix<std::size_t> I, matrix<std::size_t> J)
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].tgeabhm(alpha,A,B,beta,eval(I(m_row[h])),eval(J(m_col[h])));
        }
        fusion();
    }
    else if (m_typ==1)
    {
        matrix<std::size_t> c = col(A);
        m_dat[0] *= beta;
        m_dat[0] = horzcat(m_dat[0],alpha*eval(A(I,c)));
        m_dat[1] = vertcat(m_dat[1],eval(B(c,J)));
        std::tie(m_dat[0],m_dat[1]) = qrsvd(m_dat[0],m_dat[1],m_tol);
    }
    else if (m_typ==2)
    {
        matrix<std::size_t> c = col(A);
        m_dat[0] *= beta;
        m_dat[0] += mtimes(alpha*eval(A(I,c)),eval(B(c,J)));
        full2lowrank();
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [.tgehmhm]
///
template<typename T>
void hmatrix<T>::tgehmhm(T alpha, hmatrix<T>const& Ah, hmatrix<T>const& Bh, T beta)
{
    if (m_siz(0)!=Ah.getSize(0) || m_siz(1)!=Bh.getSize(1) || Ah.getSize(1)!=Bh.getSize(0))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    // -------------------------------------------------------------- Ch is hmatrix
    // alpha * hmatrix * hmatrix + beta * hmatrix -> hmatrix
    if (Ah.getType()==0 && Bh.getType()==0 && m_typ==0)
    {
//        std::cout << " H1 ";
        // alpha * (A11*B11 + A12*B21) + beta C11 -> C11
        m_chd[0].tgehmhm(alpha,Ah.getChildren(0),Bh.getChildren(0),beta);
        m_chd[0].tgehmhm(alpha,Ah.getChildren(1),Bh.getChildren(2),(T)1);
        
        // alpha * (A11*B12 + A12*B22) + beta C12 -> C12
        m_chd[1].tgehmhm(alpha,Ah.getChildren(0),Bh.getChildren(1),beta);
        m_chd[1].tgehmhm(alpha,Ah.getChildren(1),Bh.getChildren(3),(T)1);
        
        // alpha * (A21*B11 + A22*B21) + beta C21 -> C21
        m_chd[2].tgehmhm(alpha,Ah.getChildren(2),Bh.getChildren(0),beta);
        m_chd[2].tgehmhm(alpha,Ah.getChildren(3),Bh.getChildren(2),(T)1);
        
        // alpha * (A21*B12 + A22*B22) + beta C22 -> C223
        m_chd[3].tgehmhm(alpha,Ah.getChildren(2),Bh.getChildren(1),beta);
        m_chd[3].tgehmhm(alpha,Ah.getChildren(3),Bh.getChildren(3),(T)1);

        // Fusion
        fusion();
    }
    // alpha * compr * hmatrix + beta * hmatrix -> hmatrix
    else if (Ah.getType()==1 && Bh.getType()==0 && m_typ==0)
    {
//        std::cout << " H2 ";
        (*this).tgeabhm(alpha,Ah.getData(0),mtimes(Ah.getData(1),Bh),beta);
    }
    // alpha * hmatrix * compr + beta * hmatrix -> hmatrix
    else if (Ah.getType()==0 && Bh.getType()==1 && m_typ==0)
    {
//        std::cout << " H3 ";
        (*this).tgeabhm(alpha,mtimes(Ah,Bh.getData(0)),Bh.getData(1),beta);
    }
    // alpha * compr * compr + beta * hmatrix -> hmatrix
    else if (Ah.getType()==1 && Bh.getType()==1 && m_typ==0)
    {
//        std::cout << " H4 ";
        (*this).tgeabhm(alpha,mtimes(Ah.getData(0),mtimes(Ah.getData(1),Bh.getData(0))),Bh.getData(1),beta);
    }
    
    // -------------------------------------------------------------- Ch is compr
    // alpha * hmatrix * hmatrix + beta * compr -> compr
    else if (Ah.getType()==0 && Bh.getType()==0 && m_typ==1)
    {
//        std::cout << " C1 ";
        m_typ = 0;
        m_chd.resize(4,hmatrix<T>());
        for (int h=0; h<4; ++h)
        {
            m_row[h] = Ah.getRow(h);
            m_col[h] = Bh.getColumn(h);
            m_chd[h] = hmatrix(numel(m_row[h]),numel(m_col[h]),m_tol);
            m_chd[h].setData(0) = eval(m_dat[0](m_row[h],col(m_dat[0])));
            m_chd[h].setData(1) = eval(m_dat[1](row(m_dat[1]),m_col[h]));
        }
        m_dat[0].clear();
        m_dat[1].clear();
        (*this).tgehmhm(alpha,Ah,Bh,beta);
//        (*this).tgeabhm(alpha,full(Ah),full(Bh),beta);
    }
    // alpha * compr * hmatrix + beta * compr -> compr
    else if (Ah.getType()==1 && Bh.getType()==0 && m_typ==1)
    {
//        std::cout << " C2 ";
        (*this).tgeabhm(alpha,Ah.getData(0),mtimes(Ah.getData(1),Bh),beta);
    }
    // alpha * hmatrix * compr + beta * compr -> compr
    else if (Ah.getType()==0 && Bh.getType()==1 && m_typ==1)
    {
//        std::cout << " C3 ";
        (*this).tgeabhm(alpha,mtimes(Ah,Bh.getData(0)),Bh.getData(1),beta);
    }
    // alpha * compr * compr + beta * compr -> compr
    else if (Ah.getType()==1 && Bh.getType()==1 && m_typ==1)
    {
//        std::cout << " C4 ";
        (*this).tgeabhm(alpha,mtimes(Ah.getData(0),mtimes(Ah.getData(1),Bh.getData(0))),Bh.getData(1),beta);
    }
    // alpha * full * compr + beta * compr -> compr
    else if (Ah.getType()==2 && Bh.getType()==1 && m_typ==1)
    {
//        std::cout << " C5 ";
        (*this).tgeabhm(alpha,mtimes(Ah.getData(0),Bh.getData(0)),Bh.getData(1),beta);
    }
    // alpha * compr * full + beta * compr -> compr
    else if (Ah.getType()==1 && Bh.getType()==2 && m_typ==1)
    {
//        std::cout << " C6 ";
        (*this).tgeabhm(alpha,Ah.getData(0),mtimes(Ah.getData(1),Bh.getData(0)),beta);
    }
    // alpha * full * full + beta * compr -> full or compr
    else if (Ah.getType()==2 && Bh.getType()==2 && m_typ==1)
    {
//        std::cout << " C7 ";
        m_typ = 2;
        m_dat[0] = mtimes(m_dat[0],m_dat[1]);
        tgemm(alpha,Ah.getData(0),Bh.getData(0),beta,m_dat[0]);
        m_dat[1].clear();
        full2lowrank();
    }
    
    // -------------------------------------------------------------- Ch is full
    // alpha * compr * compr + beta * full -> full
    else if (Ah.getType()==1 && Bh.getType()==1 && m_typ==2)
    {
//        std::cout << " F1 ";
        matrix<T> tmp = mtimes(Ah.getData(0),mtimes(Ah.getData(1),Bh.getData(0)));
        tgemm(alpha,tmp,Bh.getData(1),beta,m_dat[0]);
        full2lowrank();
    }
    // alpha * full * compr + beta * full -> full
    else if (Ah.getType()==2 && Bh.getType()==1 && m_typ==2)
    {
//        std::cout << " F2 ";
        matrix<T> tmp = mtimes(Ah.getData(0),Bh.getData(0));
        tgemm(alpha,tmp,Bh.getData(1),beta,m_dat[0]);
        full2lowrank();
    }
    // alpha * compr * full + beta * full -> full
    else if (Ah.getType()==1 && Bh.getType()==2 && m_typ==2)
    {
//        std::cout << " F3 ";
        matrix<T> tmp = mtimes(Ah.getData(1),Bh.getData(0));
        tgemm(alpha,Ah.getData(0),tmp,beta,m_dat[0]);
        full2lowrank();
    }
    // alpha * full * full + beta * full -> full
    else if (Ah.getType()==2 && Bh.getType()==2 && m_typ==2)
    {
//        std::cout << " F4 ";
        tgemm(alpha,Ah.getData(0),Bh.getData(0),beta,m_dat[0]);
    }
    // Unavailable cases
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}


//==========================================================================//
//                               OPERATORS                                  //
//==========================================================================//
//==========================================================================
// [operator+=]
///
template<typename T>
hmatrix<T>& hmatrix<T>::operator+=(T b)
{
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h] += b;
        }
    }
    else if (m_typ==1)
    {
        m_dat[0] = horzcat(m_dat[0],b*ones<T>(size(m_dat[0],1),1));
        m_dat[1] = vertcat(m_dat[1],ones<T>(1,size(m_dat[1],2)));
    }
    else  if (m_typ==2)
    {
        m_dat[0] += b;
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
    return (*this);
}

//==========================================================================
// [operator*=]
///
template<typename T>
hmatrix<T>& hmatrix<T>::operator*=(T b)
{
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h] *= b;
        }
    }
    else  if (m_typ==1)
    {
        m_dat[0] *= b;
    }
    else  if (m_typ==2)
    {
        m_dat[0] *= b;
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
    return (*this);
}

//==========================================================================
// [operator+=]
///
template<typename T>
hmatrix<T>& hmatrix<T>::operator+=(hmatrix<T>const& Bh)
{
    if (m_siz(0)!=Bh.getSize(0) || m_siz(1)!=Bh.getSize(1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    // hmatrix + hmatrix -> hmatrix
    if (m_typ==0 && Bh.getType()==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h] += Bh.getChildren(h);
        }
        fusion();
    }
    // hmatrix + compr -> hmatrix
    else if (m_typ==0 && Bh.getType()==1)
    {
        tgeabhm((T)1,Bh.getData(0),Bh.getData(1),(T)1);
    }
    // compr + --- -> ---
    else if (m_typ==1)
    {
        matrix<T> A = m_dat[0];
        matrix<T> B = m_dat[1];
        (*this) = Bh;
        tgeabhm((T)1,A,B,(T)1);
    }
    // full + compr -> full or compr
    else if (m_typ==2 && Bh.getType()==1)
    {
        tgeabhm((T)1,Bh.getData(0),Bh.getData(1),(T)1);
    }
    // full + full -> full
    else if (m_typ==2 && Bh.getType()==2)
    {
        m_dat[0] += Bh.getData(0);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
    return (*this);
}

//==========================================================================
// [operator+]
///
template<typename T>
inline hmatrix<T> operator+(hmatrix<T>const& Ah, T b)
{
    hmatrix<T> Ch = Ah;
    Ch += b;
    return Ch;
}
template<typename T>
inline hmatrix<T> operator+(T a, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Bh;
    Ch += a;
    return Ch;
}
template<typename T>
inline hmatrix<T> operator+(hmatrix<T>const& Ah, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Ah;
    Ch += Bh;
    return Ch;
}

//==========================================================================
// [operator-]
///
template<typename T>
inline hmatrix<T> operator-(hmatrix<T>const& Ah)
{
    hmatrix<T> Bh = Ah;
    Bh *= ((T)-1);
    return Bh;
}
template<typename T>
inline hmatrix<T> operator-(hmatrix<T>const& Ah, T b)
{
    hmatrix<T> Ch = Ah;
    Ch += (-b);
    return Ch;
}
template<typename T>
inline hmatrix<T> operator-(T a, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Bh;
    Ch *= ((T)-1);
    Ch += a;
    return Ch;
}
template<typename T>
inline hmatrix<T> operator-(hmatrix<T>const& Ah, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Bh;
    Ch *= -1;
    Ch += Ah;
    return Ch;
}

//==========================================================================
// [operator*]
///
template<typename T>
inline hmatrix<T> operator*(hmatrix<T>const& Ah, T b)
{
    hmatrix<T> Ch = Ah;
    Ch *= b;
    return Ch;
}
template<typename T>
inline hmatrix<T> operator*(T a, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Bh;
    Ch *= a;
    return Ch;
}

//==========================================================================
// [operator/]
///
template<typename T>
inline hmatrix<T> operator/(hmatrix<T>const& Ah, T b)
{
    hmatrix<T> Ch = Ah;
    Ch *= ((T)1)/b;
    return Ch;
}

//==========================================================================
// [operator<<]
///
template<typename T>
std::ostream& operator<<(std::ostream& flux, hmatrix<T> const& Ah)
{
    matrix<std::size_t> info(1,4);
    Ah.stat(info);
    double memBytes = info(3)*sizeof(T);
    flux << "H-matrix " << size(Ah,1) << "x" << size(Ah,2);
    flux << " of type '" << typeid(T).name() << "' (";
    if (memBytes<1e3) {flux << memBytes << " B)";}
    else if (memBytes<1e6) {flux << memBytes/1e3 << " KB)";}
    else {flux << memBytes/1e6 << " MB)";}
    flux << std::endl;
    flux << "  Leaf : " << info(0) << "h " << info(1) << "c " << info(2) << "f " << std::endl;
    flux << "  Rate : " << 100.*info(3)/(size(Ah,1)*size(Ah,2)) << " %" << std::endl;
    flux << "  tol  : " << Ah.getTolerance();
    return flux;
}


//==========================================================================//
//                       MATLAB-LIKE FUNCTION                               //
//==========================================================================//
//==========================================================================
// [full]
/// Hierarchical to dense conversion.
template<typename T>
inline matrix<T> full(hmatrix<T>const& Ah)
{
    matrix<T> A(size(Ah,1),size(Ah,2));
    Ah.hfull(A);
    return A;
}

//==========================================================================
// [gmres]
/// Iterative solver.
template<typename T>
matrix<T> gmres(hmatrix<T>const& Ah, matrix<T>const& B,
                double tol = 1e-6, std::size_t maxit = 10,
                std::function<matrix<T>(matrix<T>const&)>const& Am1 = std::function<matrix<T>(matrix<T>const&)>(),
                matrix<T>const& X0 = matrix<T>())
{
    std::function<matrix<T>(matrix<T>const&)> Afct;
    Afct = [&Ah](matrix<T>const& X) {return mtimes(Ah,X);};
    return gmres(Afct,B,tol,maxit,Am1,X0);
}
template<typename T>
matrix<T> gmres(hmatrix<T>const& Ah, matrix<T>const& B, double tol, std::size_t maxit,
                hmatrix<T>const& Ahm1, matrix<T>const& X0 = matrix<T>())
{
    std::function<matrix<T>(matrix<T>const&)> Afct, Am1fct;
    Afct   = [&Ah](matrix<T>const& X) {return mtimes(Ah,X);};
    Am1fct = [&Ahm1](matrix<T>const& X) {return mtimes(Ahm1,X);};
    return gmres(Afct,B,tol,maxit,Am1fct,X0);
}
template<typename T>
matrix<T> gmres(hmatrix<T>const& Ah, matrix<T>const& B, double tol, std::size_t maxit,
                hmatrix<T>const& Lh, hmatrix<T>const& Uh, matrix<T>const& X0 = matrix<T>())
{
    std::function<matrix<T>(matrix<T>const&)> Afct, Am1fct;
    Afct   = [&Ah](matrix<T>const& X) {return mtimes(Ah,X);};
    Am1fct = [&Lh,&Uh](matrix<T>const& B)
    {
        matrix<T> X = B;
        Lh.hllowsolve(X);
        Uh.hlupsolve(X);
        return X;
    };
    return gmres(Afct,B,tol,maxit,Am1fct,X0);
}

//==========================================================================
// [inv]
/// Hierarchical matrix inverse.
template<typename T>
inline hmatrix<T> inv(hmatrix<T>const& Ah)
{
    hmatrix<T> Ahm1 = Ah;
    Ahm1.hinv();
    return Ahm1;
}

//==========================================================================
// [linsolve]
/// Solve linear system Ah*X=B.
///
/// mtimes(Ah,B) is the hierarchical matrix product of A and B, where Ah is
/// H-Matrix and B dense matrix. The number of columns of Ah must equal
/// to the number of rows of B. Use Hierarchical LU factorization.
template<typename T>
inline matrix<T> linsolve(hmatrix<T>const& Ah, matrix<T>const& B)
{
    hmatrix<T> Lh, Uh;
    std::tie(Lh,Uh) = lu(Ah);
    matrix<T> X = B;
    Lh.hllowsolve(X);
    Uh.hlupsolve(X);
    return X;
}

//==========================================================================
// [lu]
/// Hierarchical lower/upper factorization.
template<typename T>
inline auto lu(hmatrix<T>const& Ah, double tol=0)
{
    hmatrix<T> Lh=Ah, Uh;
    if (tol>0) {Lh.recompress(tol);}
    Lh.hlu(Uh);
    return std::make_tuple(Lh,Uh);
}

//==========================================================================
// [mtimes]
/// Hierarchical-matrix multiply.
///
/// mtimes(A,B) is the hierarchical matrix product of A and B, where A and
/// B are H-Matrix and/or dense matrix. The number of columns of A must equal
/// the number of rows of B.
template<typename T>
inline matrix<T> mtimes(hmatrix<T>const& Ah, matrix<T>const& B)
{
    matrix<T> C(size(Ah,1),size(B,2));
    Ah.hlmtimes(B,C);
    return C;
}
template<typename T>
inline matrix<T> mtimes(matrix<T>const& A, hmatrix<T>const& Bh)
{
    matrix<T> C(size(A,1),size(Bh,2));
    Bh.hrmtimes(A,C);
    return C;
}
template<typename T>
inline hmatrix<T> mtimes(hmatrix<T>const& Ah, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch(size(Ah,1),size(Bh,2),std::max(Ah.getTolerance(),Bh.getTolerance()));
    tgemm((T)1,Ah,Bh,(T)1,Ch);
    return Ch;
}

//==========================================================================
// [recompress]
/// Leaves recompression with new tolerance.
template<typename T>
inline hmatrix<T> recompress(hmatrix<T>const& Ah, double tol)
{
    hmatrix<T> Bh = Ah;
    Bh.recompress(tol);
    return Bh;
}

//==========================================================================
// [size]
/// Size of H-matrix.
///
/// S = size(Ah) for m-by-n H-matrix A returns the two-element vector [m,n]
/// containing the number of rows and columns in the matrix.
///
/// S = size(Ah,dim) returns the lengths of the specified dimensions dim.
///
/// \see length, numel.
template<typename T>
inline matrix<std::size_t> size(hmatrix<T>const& Ah)
{
    return {size(Ah,1),size(Ah,2)};
}
template<typename T>
inline std::size_t size(hmatrix<T>const& Ah, int dim)
{
    return Ah.getSize(dim-1);
}

//==========================================================================
// [spy]
/// Spy hierarchical structure.
template<typename T>
inline matrix<T> spy(hmatrix<T>const& Ah)
{
    if (size(Ah,1)*size(Ah,2)>1e7)
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"H-matrix is too large to be spied on.");
        return 0;
    }
    else
    {
        matrix<T> A(size(Ah,1),size(Ah,2));
        Ah.spy(A);
        return A;
    }
}

//==========================================================================
// [tgeabm]
/// In-place low-rank matrix product with H-matrix.
///
/// tgeabm(alpha,A,B,beta,Ch) performs the in-place matrix-hmatrix operations
///    C = alpha*A*B + beta*Ch,
/// where alpha, beta are scalars and A, B are matrices with compatible size
/// and Ch is H-Matrix.
template<typename T>
inline void tgeabm(T alpha, matrix<T>const& A, matrix<T>const& B, T beta, hmatrix<T>& Ch)
{
    Ch.tgeabhm(alpha,A,B,beta);
}

//==========================================================================
// [tgemm]
/// In-place H-matrix product.
///
/// tgemm(alpha,Ah,Bh,beta,Ch) performs the in-place hmatrix-hmatrix operations
///    C = alpha*Ah*Bh + beta*Ch,
/// where alpha, beta are scalars and Ah, Bh, Ch are hmatrices with compatible size.
template<typename T>
inline void tgemm(T alpha, hmatrix<T>const& Ah, hmatrix<T>const& Bh, T beta, hmatrix<T>& Ch)
{
    Ch.tgehmhm(alpha,Ah,Bh,beta);
}

//==========================================================================
// [transpose]
/// Non-conjugate transpose.
///
/// transpose(Ah) is called for the syntax A^t.
template<typename T>
inline hmatrix<T> transpose(hmatrix<T>const& Ah)
{
    hmatrix<T> Aht = Ah;
    Aht.htranspose();
    return Aht;
}

}





//==========================================================================
//        using T2 = decltype(std::abs(m_dat[0](0)));
//        matrix<T2> S;
//        matrix<T> U, Vt;
//        std::size_t rk;
//        std::tie(S,U,Vt) = svd(m_dat[0],"vect");
//        rk = sum(S>=S(0)*m_tol);
//        if (S(numel(S)-1)<max(size(m_dat[0]))*M_EPS(T2) && rk<min(size(m_dat[0]))/4)
//        {
//            matrix<std::size_t> R = range(0,rk);
//            m_dat[0] = mtimes( eval(U(row(U),R)) , diag(eval(S(R))) );
//            m_dat[1] = eval(Vt(R,col(Vt)));
//            m_typ    = 1;
//        }
