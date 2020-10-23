#include <catch.hpp>
#include <castor/matrix.hpp>

using namespace castor;

TEST_CASE( "Internal operator : matrix<T>& operator+=(matrix<T>const& A)", "[Internal operator]" ) 
{
    matrix<int> V = {0,1,2,3};
    V += matrix<int>(1,4,1);
    REQUIRE( V(0) == 1 );
    REQUIRE( V(1) == 2 );
    REQUIRE( V(2) == 3 );
    REQUIRE( V(3) == 4 );
}

TEST_CASE( "Internal operator : matrix<T>& operator-=(matrix<T>const& A)", "[Internal operator]" ) 
{
    matrix<int> V = {0,1,2,3};
    V -= matrix<int>(1,4,1);
    REQUIRE( V(0) == -1 );
    REQUIRE( V(1) == 0 );
    REQUIRE( V(2) == 1 );
    REQUIRE( V(3) == 2 );
}

TEST_CASE( "Internal operator : matrix<T>& operator*=(matrix<T>const& A)", "[Internal operator]" ) 
{
    matrix<int> V = {0,1,2,3};
    V *= matrix<int>(1,4,2);
    REQUIRE( V(0) == 0 );
    REQUIRE( V(1) == 2 );
    REQUIRE( V(2) == 4 );
    REQUIRE( V(3) == 6 );
}

TEST_CASE( "Internal operator : matrix<T>& operator/=(matrix<T>const& A)", "[Internal operator]" ) 
{
    matrix<int> V = {0,2,4,6};
    V /= matrix<int>(1,4,2);
    REQUIRE( V(0) == 0 );
    REQUIRE( V(1) == 1 );
    REQUIRE( V(2) == 2 );
    REQUIRE( V(3) == 3 );
}

TEST_CASE( "Internal operator : std::size_t size(int dim=0) const", "[Internal operator]" ) 
{
    matrix<int> A({{0,1,2,3},{4,5,6,7},{8,9,10,11}});
    REQUIRE( A.size() == 12 );
    REQUIRE( A.size(1) == 3 );
    REQUIRE( A.size(2) == 4 );
}
