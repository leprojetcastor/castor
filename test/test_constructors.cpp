#include <catch.hpp>
#include <castor/matrix.hpp>

using namespace castor;

TEST_CASE( "Constructor : matrix(T v) ", "[constructor]" ) 
{
    matrix<int> A(1);
    REQUIRE( A.size() == 1 );
    REQUIRE( A.size(1) == 1 );
    REQUIRE( A.size(2) == 1 );
    REQUIRE( A(0) == 1 );
}

TEST_CASE( "Constructor : matrix(std::initializer_list<T>const& v)", "[constructor]" ) 
{
    matrix<float> A({0,1,2,M_PI});
    REQUIRE( A.size() == 4 );
    REQUIRE( A.size(1) == 1 );
    REQUIRE( A.size(2) == 4 );
    REQUIRE( A(0) == 0 );
    REQUIRE( A(1) == 1 );
    REQUIRE( A(2) == 2 );
    REQUIRE( A(3) == Approx(M_PI) );
}

TEST_CASE( "Constructor : matrix(std::initializer_list<std::vector<T>>const& v)", "[constructor]" )
{
    matrix<double> A({{0,1,2,M_PI},{4,5,6,7},{8,9,10,11}});
    REQUIRE( A.size() == 12 );
    REQUIRE( A.size(1) == 3 );
    REQUIRE( A.size(2) == 4 );
    REQUIRE( A(0) == 0 );
    REQUIRE( A(1) == 1 );
    REQUIRE( A(2) == 2 );
    REQUIRE( A(3) == M_PI );
    REQUIRE( A(4) == 4 );
    REQUIRE( A(5) == 5 );
    REQUIRE( A(6) == 6 );
    REQUIRE( A(7) == 7);
    REQUIRE( A(8) == 8 );
    REQUIRE( A(9) == 9 );
    REQUIRE( A(10) == 10 );
    REQUIRE( A(11) == 11);
}

TEST_CASE( "Constructor : matrix(std::size_t m, std::size_t n, T v=0)", "[constructor]")
{
    matrix<std::complex<double>> A(2,3,std::complex<double>(1,M_PI));
    REQUIRE( A.size() == 6 );
    REQUIRE( A.size(1) == 2 );
    REQUIRE( A.size(2) == 3 );
    for (std::size_t i=0; i<A.size(); ++i)
    {
        REQUIRE( A(i) == 1. + M_PI*M_1I );
    }
}

TEST_CASE( "Constructor : matrix(std::vector<S>const& v)", "[constructor]")
{
    matrix<double> A(std::vector<int>({0, 1, 2, 3}));
    REQUIRE( A.size() == 4 );
    REQUIRE( A.size(1) == 1 );
    REQUIRE( A.size(2) == 4 );
    REQUIRE( A(0) == 0 );
    REQUIRE( A(1) == 1 );
    REQUIRE( A(2) == 2 );
    REQUIRE( A(3) == 3 );
}

TEST_CASE( "Constructor : matrix(std::size_t m, std::size_t n, std::vector<S>const& v)", "[constructor]")
{
    matrix<double> A(2,2,{0, 1, 2, 3});
    REQUIRE( A.size() == 4 );
    REQUIRE( A.size(1) == 2 );
    REQUIRE( A.size(2) == 2 );
    REQUIRE( A(0) == 0 );
    REQUIRE( A(1) == 1 );
    REQUIRE( A(2) == 2 );
    REQUIRE( A(3) == 3 );
}

TEST_CASE( "Constructor : matrix(std::size_t m, std::size_t n, std::vector<std::vector<S>>const& v)", "[constructor]")
{
    matrix<> A(2,2,{{0, 1},{2,3}});
    REQUIRE( A.size() == 4 );
    REQUIRE( A.size(1) == 2 );
    REQUIRE( A.size(2) == 2 );
    REQUIRE( A(0) == 0 );
    REQUIRE( A(1) == 1 );
    REQUIRE( A(2) == 2 );
    REQUIRE( A(3) == 3 );
}

TEST_CASE( "Constructor : matrix(matrix<S>const& A)", "[constructor]")
{
    matrix<> B(2,2,{{0, 1},{2,3}});
    matrix<> A(B);
    REQUIRE( A.size() == 4 );
    REQUIRE( A.size(1) == 2 );
    REQUIRE( A.size(2) == 2 );
    REQUIRE( A(0) == 0 );
    REQUIRE( A(1) == 1 );
    REQUIRE( A(2) == 2 );
    REQUIRE( A(3) == 3 );
}
