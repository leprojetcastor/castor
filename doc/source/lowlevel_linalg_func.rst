
Low-level interface
+++++++++++++++++++

The user will find here the low-level linear algebra interface. Most of the functions here are wrappers around the CBLAS or LAPACK API and the user should be careful when using them, see :ref:`label-blaslapack-issue` to understand why.

.. _label-lpk2mat:

lpk2mat
-------
.. doxygenfunction:: lpk2mat(std::vector<S>const &V, std::size_t m, std::size_t n)
   :project: castor

See :ref:`label-mat2lpk`.


.. _label-mat2lpk:

mat2lpk
-------
.. doxygenfunction:: mat2lpk(matrix<S> const &A, std::size_t L)
   :project: castor

See :ref:`label-lpk2mat`.


.. _label-tgeev:

tgeev
-----
.. doxygenfunction:: tgeev(std::string typ, int &n, std::vector<T> &A, std::vector<T> &E, std::vector<T> &V)
   :project: castor
.. doxygenfunction:: xgeev(char &jobl, char &jobr, int &n, std::vector<zlpk> &A, std::vector<zlpk> &E, std::vector<zlpk> &V, zlpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgeev(char &jobl, char &jobr, int &n, std::vector<clpk> &A, std::vector<clpk> &E, std::vector<clpk> &V, clpk &wkopt, int &lwork, int &info)
   :project: castor


.. _label_tgels:

tgels
-----
.. doxygenfunction:: tgels(int &m, int &n, int &nrhs, std::vector<T> &A, std::vector<T> &B)
   :project: castor
.. doxygenfunction:: xgels(char &t, int &m, int &n, int &nrhs, std::vector<zlpk> &A, std::vector<zlpk> &B, int &l, zlpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgels(char &t, int &m, int &n, int &nrhs, std::vector<double> &A, std::vector<double> &B, int &l, double &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgels(char &t, int &m, int &n, int &nrhs, std::vector<clpk> &A, std::vector<clpk> &B, int &l, clpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgels(char &t, int &m, int &n, int &nrhs, std::vector<float> &A, std::vector<float> &B, int &l, float &wkopt, int &lwork, int &info)
   :project: castor


.. _label-tgemm-blas:

tgemm
-----
.. doxygenfunction:: tgemm(R alpha, matrix<double> const &A, matrix<double> const &B, S beta, matrix<double> &C)
   :project: castor
.. doxygenfunction:: tgemm(R alpha, matrix<std::complex<float>> const &A, matrix<std::complex<float>> const &B, S beta, matrix<std::complex<float>> &C)
   :project: castor
.. doxygenfunction:: tgemm(R alpha, matrix<std::complex<double>> const &A, matrix<std::complex<double>> const &B, S beta, matrix<std::complex<double>> &C)
   :project: castor
.. doxygenfunction:: tgemm(R alpha, matrix<float> const &A, matrix<float> const &B, S beta, matrix<float> &C)
   :project: castor

See :ref:`label-tgemm-naive`.


.. _label-tgeqrf:

tgeqrf
------
.. doxygenfunction:: tgeqrf(int &m, int &n, std::vector<T> &A, std::vector<T> &R)
   :project: castor
.. doxygenfunction:: xgeqrf(int &m, int &n, int &l, std::vector<zlpk> &A, std::vector<zlpk> &tau, zlpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgeqrf(int &m, int &n, int &l, std::vector<double> &A, std::vector<double> &tau, double &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgeqrf(int &m, int &n, int &l, std::vector<clpk> &A, std::vector<clpk> &tau, clpk &wkopt, int &lwork, int &info)
   :project: castor
.. doxygenfunction:: xgeqrf(int &m, int &n, int &l, std::vector<float> &A, std::vector<float> &tau, float &wkopt, int &lwork, int &info)
   :project: castor


.. _label-tgesdd:

tgesdd
------
.. doxygenfunction:: tgesdd(std::string typ, int &m, int &n, std::vector<T> &A, std::vector<T2> &S, std::vector<T> &U, std::vector<T> &V)
   :project: castor
.. doxygenfunction:: xgesdd(char &job, int &m, int &n, int &l, std::vector<zlpk> &A, std::vector<double> &S, std::vector<zlpk> &U, std::vector<zlpk> &V, zlpk &wkopt, int &lwork, std::vector<int> &iwork, int &info)
   :project: castor
.. doxygenfunction:: xgesdd(char &job, int &m, int &n, int &l, std::vector<double> &A, std::vector<double> &S, std::vector<double> &U, std::vector<double> &V, double &wkopt, int &lwork, std::vector<int> &iwork, int &info)
   :project: castor
.. doxygenfunction:: xgesdd(char &job, int &m, int &n, int &l, std::vector<clpk> &A, std::vector<float> &S, std::vector<clpk> &U, std::vector<clpk> &V, clpk &wkopt, int &lwork, std::vector<int> &iwork, int &info)
   :project: castor
.. doxygenfunction:: xgesdd(char &job, int &m, int &n, int &l, std::vector<float> &A, std::vector<float> &S, std::vector<float> &U, std::vector<float> &V, float &wkopt, int &lwork, std::vector<int> &iwork, int &info)
   :project: castor


.. _label-tgesv:

tgesv
-----
.. doxygenfunction:: tgesv(int &n, int &nrhs, std::vector<T> &A, std::vector<T> &B)
   :project: castor
.. doxygenfunction:: xgesv(int &n, int &nrhs, std::vector<zlpk> &A, std::vector<zlpk> &B, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgesv(int &n, int &nrhs, std::vector<double> &A, std::vector<double> &B, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgesv(int &n, int &nrhs, std::vector<clpk> &A, std::vector<clpk> &B, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgesv(int &n, int &nrhs, std::vector<float> &A, std::vector<float> &B, std::vector<int> &ipiv, int &info)
   :project: castor


.. _label-tgetrf:

tgetrf
------
.. doxygenfunction:: tgetrf(int &m, int &n, std::vector<T> &A, std::vector<T> &P, matrix<S> &U)
   :project: castor
.. doxygenfunction:: xgetrf(int &m, int &n, std::vector<zlpk> &A, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgetrf(int &m, int &n, std::vector<double> &A, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgetrf(int &m, int &n, std::vector<clpk> &A, std::vector<int> &ipiv, int &info)
   :project: castor
.. doxygenfunction:: xgetrf(int &m, int &n, std::vector<float> &A, std::vector<int> &ipiv, int &info)
   :project: castor
