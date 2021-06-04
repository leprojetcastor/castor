/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : overview_matrix.cpp                           |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Overview of all fonctionalities of matrix     |
 |  `---'  |                framework                                     |
 +========================================================================+
 */

#include <castor/matrix.hpp>

using namespace castor;

//==============================================================================
int main (int argc, char* argv[])
{
    //================================================== TOOLS
    std::cout << "+=====================+" << std::endl;
    std::cout << "|       TOOLS         |" << std::endl;
    std::cout << "+=====================+" << std::endl;
    documentationFiles =
    {
        "/usr/local/include/castor/matrix.hpp" // Full name of your file matrix.hpp
    };
    help("help");

    tic();
    while (toc(0)<0.314159) {}
    toc();

    disp("hello world");
    disp(M_PI);

    disp(eye<logical>(10,10));
    disp(10*eye<int>(10,10));
    disp(100*eye(10,10));
    disp(100*eye(10,10) + M_1I*eye(10,10));

    warning(__FILE__, __LINE__, __FUNCTION__,"This is just a funny warning, enjoy!.");
    
    // =============================================================== LOGICAL
    std::cout << "+===============+" << std::endl;
    std::cout << "|    LOGICAL    |" << std::endl;
    std::cout << "+===============+" << std::endl;

    logical yes=true, no=0.;
    disp(yes,2);
    disp(no,2);
    disp(yes+no,2);
    disp(yes*no,2);
    disp(yes&&no,2);
    disp(yes||no,2);
    disp(yes+1,2);
    yes = yes + 1;
    disp(yes,2);

    //================================================== CONSTRUCTORS
    std::cout << "+=====================+" << std::endl;
    std::cout << "|    CONSTRUCTORS     |" << std::endl;
    std::cout << "+=====================+" << std::endl;

    matrix<logical> A1;
    disp(A1);
    matrix<int> A2(3);
    disp(A2);
    matrix<float> A3({0,1,2,M_PI});
    disp(A3);
    matrix<double> A4({{0,1,2,M_PI},{4,5,6,7},{8,9,10,11}});
    disp(A4);
    matrix<std::complex<float>> A5(3,4,M_PI);
    disp(A5);
    matrix<double> A6(std::vector<int>({0, 1, 2, 3}));
    disp(A6);
    matrix<float> A7(2,2,{0, 1, 2, 3});
    disp(A7);
    matrix<float> A8(2,2,{{0, 1},{2,3}});
    disp(A8);
    matrix<int> A9(A8);
    disp(A9);

    matrix<int> s = M_PI;
    int u = (double)s;
    disp(u);
    matrix<int> X = {0,1,2,3};
    disp(X);
    matrix<logical> Y = X;
    disp(Y);
    matrix<> A = {{0,1,2,M_PI},
        {4,5,6,7},
        {8,9,10,11}};
    disp(A);
    matrix<> B = A;
    disp(B);
    matrix<> V = X;
    disp(V);
    std::cout << V << std::endl;
    
    // =============================================================== INTERNAL OPERATORS
    std::cout << "+====================+" << std::endl;
    std::cout << "| INTERNAL OPERATORS |" << std::endl;
    std::cout << "+====================+" << std::endl;
    
    std::cout << "matrix : " << V << std::endl;
    V += 1;
    std::cout << "+= 1   : " << V << std::endl;
    V += {1,1,1,1};
    std::cout << "+= {}  : " << V << std::endl;
    V += std::vector<int>(4,1);
    std::cout << "+= v   : " << V << std::endl;
    V += matrix<double>(1,4,1);
    std::cout << "+= A   : " << V << std::endl;
    V += matrix<int>(1,4,1);
    std::cout << "+= A<> : " << V << std::endl;
    
    V = {0,1,2,3};
    std::cout << "matrix : " << V << std::endl;
    V -= 1;
    std::cout << "-= 1   : " << V << std::endl;
    V -= matrix<int>(1,4,1);
    std::cout << "-= A<> : " << V << std::endl;
    
    V = {0,1,2,3};
    std::cout << "matrix : " << V << std::endl;
    V *= 2;
    std::cout << "*= 2   : " << V << std::endl;
    V *= matrix<int>(1,4,2);
    std::cout << "*= A<> : " << V << std::endl;
    
    V = {0,1,2,3};
    std::cout << "matrix : " << V << std::endl;
    V /= 2;
    std::cout << "/= 2   : " << V << std::endl;
    V /= matrix<int>(1,4,2);
    std::cout << "/= A<> : " << V << std::endl;
    
    V = {0,1,2,3};
    std::cout << "matrix    : " << V << std::endl;
    std::cout << "V(0)      : " << V(0) << std::endl;
    V(0) = -1;
    std::cout << "V(0)=-1   : " << V << std::endl;
    std::cout << "V(0,0)    : " << V(0,0) << std::endl;
    V(0,0) = -2;
    std::cout << "V(0,0)=-2 : " << V << std::endl;

    // =============================================================== INTERNAL OPERATORS
    std::cout << "+====================+" << std::endl;
    std::cout << "|   INTERNAL TOOLS   |" << std::endl;
    std::cout << "+====================+" << std::endl;

    disp(A.size());
    disp(A.size(1));
    disp(A.size(2));

    A.reshape(4,3);
    disp(A);
    
    A.resize(4,5);
    disp(A);
    A.clear();
    disp(A);
    A = B;
    disp(matrix<int>(A.val()));
    
    // =============================================================== VIEW
    std::cout << "+================+" << std::endl;
    std::cout << "|      VIEW      |" << std::endl;
    std::cout << "+================+" << std::endl;
    
    // disp(A({1,3,5}));
    // std::cout << (A({1,3,5}));
    
    matrix<float> C = eval(A({1,3,5}));
    disp(C);
    C = eval(A(range(0,4)));
    disp(C);
    C = -eval(A(range(0,4)));
    disp(C);
    
    C = A;
    V = {-1,-1,-1,-1};
    C(range(0,4)) = 0;
    disp(C);
    C(range(0,4)) = {1,1,1,1};
    disp(C);
    C(range(0,4)) = matrix<int>({2,2,2,2});
    disp(C);
    C(range(0,4)) = eval(C(range(4,8)));
    disp(C);
    
    C = eval(A({1,2},col(A)));
    disp(C);
    C = eval(A(range(0,2),col(A)));
    disp(C);
    C = -eval(A(range(0,2),col(A)));
    disp(C);
    
    C = A;
    C(0,col(C)) = 0;
    disp(C);
    C({0,1},col(C)) = {{1,1,1,1},{1,1,1,1}};
    disp(C);
    C({0,1},col(C)) = matrix<int>({{2,2,2,2},{1,1,1,1}});
    disp(C);
    C(range(0,2),col(C)) = eval(C({1,0},col(C)));
    disp(C);
    
    matrix<>const M = eye(3,4);
    disp(eval(M(all(M))));
    disp(eval(M(row(M),col(M))));

    // =============================================================== EXTERNAL OPERATORS
    std::cout << "+====================+" << std::endl;
    std::cout << "| EXTERNAL OPERATORS |" << std::endl;
    std::cout << "+====================+" << std::endl;

    V = {0,1,2,3};
    std::cout << "matrix : " << V << std::endl;
    std::cout << "A+1    : " << (V + 1) << std::endl;
    std::cout << "1+A    : " << (1 + V) << std::endl;
    std::cout << "A+v    : " << (V + sum(V)) << std::endl;
    std::cout << "v+A    : " << (sum(V) + V) << std::endl;
    std::cout << "A+B<>  : " << (V + matrix<int>(1,4,1)) << std::endl;
    std::cout << "A<>+B  : " << (matrix<int>(1,4,1) + V) << std::endl;
    std::cout << "[]+1   : " << (matrix<>() + 1) << std::endl;
    std::cout << "1+[]   : " << (1 + matrix<>()) << std::endl;

    std::cout << "matrix : " << V << std::endl;
    std::cout << "-A     : " << (-V) << std::endl;
    std::cout << "A-1    : " << (V - 1) << std::endl;
    std::cout << "1-A    : " << (1 - V) << std::endl;
    std::cout << "A-B<>  : " << (V - matrix<int>(1,4,1)) << std::endl;
    std::cout << "A<>-B  : " << (matrix<int>(1,4,1) - V) << std::endl;

    std::cout << "matrix : " << V << std::endl;
    std::cout << "A*2    : " << (V * 2) << std::endl;
    std::cout << "2*A    : " << (2 * V) << std::endl;
    std::cout << "A*B<>  : " << (V * matrix<int>(1,4,2)) << std::endl;
    std::cout << "A<>*B  : " << (matrix<int>(1,4,2) * V) << std::endl;

    std::cout << "matrix : " << V << std::endl;
    std::cout << "A/2    : " << (V / 2) << std::endl;
    std::cout << "2/A    : " << (2 / V) << std::endl;
    std::cout << "A/B<>  : " << (V / matrix<int>(1,4,2)) << std::endl;
    std::cout << "A<>/B  : " << (matrix<int>(1,4,2) / V) << std::endl;

    // =============================================================== COMPARE
    std::cout << "+===============+" << std::endl;
    std::cout << "|    COMPARE    |" << std::endl;
    std::cout << "+===============+" << std::endl;

    matrix<int> I = {1,1,0,0};
    matrix<logical> J = {0,1,0,1};
    std::cout << "original: " << I << std::endl;
    std::cout << "          " << J << std::endl;
    std::cout << "!A: " << (!I) << std::endl;
    std::cout << "A&&B: " << (I&&J) << std::endl;
    std::cout << "A&&v: " << (I&&0) << std::endl;
    std::cout << "v&&B: " << (0&&J) << std::endl;
    std::cout << "A||B: " << (I||J) << std::endl;
    std::cout << "A||v: " << (I||1) << std::endl;
    std::cout << "v||B: " << (1||J) << std::endl;
    std::cout << "A==B: " << (I==J) << std::endl;
    std::cout << "A==v: " << (I==1) << std::endl;
    std::cout << "v==A: " << (1==I) << std::endl;
    std::cout << "A!=B: " << (I!=J) << std::endl;
    std::cout << "A!=v: " << (I!=1) << std::endl;
    std::cout << "v!=A: " << (1!=I) << std::endl;
    std::cout << "A<=B: " << (I<=J) << std::endl;
    std::cout << "A<=v: " << (I<=0) << std::endl;
    std::cout << "v<=A: " << (1<=I) << std::endl;
    std::cout << "A<B: " << (I<J) << std::endl;
    std::cout << "A<v: " << (I<1) << std::endl;
    std::cout << "v<A: " << (0<I) << std::endl;
    std::cout << "A>=B: " << (I>=J) << std::endl;
    std::cout << "A>=v: " << (I>=1) << std::endl;
    std::cout << "v>=A: " << (0>=I) << std::endl;
    std::cout << "A>B: " << (I>J) << std::endl;
    std::cout << "A>v: " << (I>0) << std::endl;
    std::cout << "v>A: " << (1>I) << std::endl;

    //================================================== DIMENSIONS
    std::cout << "+=====================+" << std::endl;
    std::cout << "|     DIMENSIONS      |" << std::endl;
    std::cout << "+=====================+" << std::endl;

    std::cout << length(s) << std::endl;
    std::cout << length(X) << std::endl;
    std::cout << length(A) << std::endl;
    std::cout << nnz(A) << std::endl;
    std::cout << numel(s) << std::endl;
    std::cout << numel(X) << std::endl;
    std::cout << numel(A) << std::endl;
    std::cout << size(A,0) << std::endl;
    std::cout << size(A,1) << std::endl;
    std::cout << size(A,2) << std::endl;
    disp(size(s));
    disp(size(X));
    disp(size(A));

    //================================================== BUILDERS
    std::cout << "+================+" << std::endl;
    std::cout << "|    BUILDERS    |" << std::endl;
    std::cout << "+================+" << std::endl;

    disp(colon(0,7));
    disp(colon(0,1.5,7));
    disp(diag(X));
    disp(diag(X,2));
    disp(diag(X,-2));
    disp(diag(A));
    disp(diag(A,2));
    disp(diag(A,-2));
    disp(eye<logical>(size(A)));
    disp(eye(4,3));
    disp(linspace(1,0,5));
    disp(logspace(0,3,4));
    disp(ones<logical>(size(A)));
    disp(ones(3,4));
    disp(rand<float>(size(A)));
    disp(rand(3,4));
    disp(zeros<logical>(size(A)));
    disp(zeros(3,4));

    //================================================== I/O
    std::cout << "+=====================+" << std::endl;
    std::cout << "|     INPUT/OUTPUT    |" << std::endl;
    std::cout << "+=====================+" << std::endl;

    matrix<> in = eye(3,4);
    in(2,0) = -M_PI;
    disp(in);
    writetxt("./","matrix_IOtest.txt",in);
    auto txt = readtxt("./","matrix_IOtest.txt");
    disp(txt);
    writebin("./","matrix_IOtest.bin",in);
    auto bin = readbin("./","matrix_IOtest.bin");
    disp(bin);

    //================================================== MANIPULATIONS
    std::cout << "+=====================+" << std::endl;
    std::cout << "|    MANIPULATIONS    |" << std::endl;
    std::cout << "+=====================+" << std::endl;

    disp(cast<int>(A));
    disp(cast(A,int()));
    disp(cast(A,logical()));
    disp(cat(1,0,0));
    disp(cat(2,0,0));
    disp(cat(1,0,ones(3,1)));
    disp(cat(2,ones(1,4),0));
    disp(cat(1,cast<int>(A),eye<float>(size(A))));
    disp(cat(2,A,eye(size(A))));
    clear(C);
    disp(C);
    disp(find(eye(3,4)));
    disp(all(A));
    disp(row(A));
    disp(col(A));
    disp(get(A,all(A)));
    disp(get(A,2*eye<std::size_t>(3,4)));
    disp(get(A,row(A),0));
    disp(get(A,0,col(A)));
    disp(get(A,range(0,3),range(0,4)));
    disp(horzcat(A,eye<int>(3,4)));
    matrix<std::size_t> idx,jdx;
    std::tie(idx,jdx) = ind2sub({2,3},range(0,6));
    disp(vertcat(idx,jdx));
    disp(isempty(matrix<>()));
    disp(isempty(V));
    disp(isequal(V,V));
    disp(isequal(V,V+1));
    disp(isvector(V));
    disp(isvector(ones(3)));
    disp(reshape(A,4,3));
    disp(resize(A,4,4));
    B = A;
    set(B,all(B),M_PI);
    disp(B);
    set(B,all(B),rand(1,12));
    disp(B);
    set(B,row(B),0,M_PI);
    disp(B);
    set(B,row(B),0,ones(3,1));
    disp(B);
    set(B,0,col(B),M_PI);
    disp(B);
    set(B,0,col(B),ones(1,4));
    disp(B);
    disp(sub2ind({2,3},{0,1,2},{0,1,2}));
    disp(transpose(A));
    disp(vertcat(A,eye<int>(3,4)));
    disp(values(A));
    B = eye(2,3);
    B(0) = NAN;
    B(1) = INFINITY;
    B(2) = std::exp(800);
    disp(B);
    disp(isnan(B));
    disp(isinf(B));
    disp(isfinite(B));

    //================================================== ALGORITHMS
    std::cout << "+=====================+" << std::endl;
    std::cout << "|     ALGORITHMS      |" << std::endl;
    std::cout << "+=====================+" << std::endl;

    disp(argmax(A));
    disp(argmax(A,1));
    disp(argmax(A,2));
    disp(argmin(A));
    disp(argmin(A,1));
    disp(argmin(A,2));
    disp(argsort(eye(3,4)));
    disp(argsort(eye(3,4),1));
    disp(argsort(eye(3,4),2));
    disp(max(A));
    disp(max(A,1));
    disp(max(A,2));
    disp(maximum(cast(A,int()),A+1));
    disp(maximum(cast(A,int()),5.5));
    disp(maximum(5,cast(A,int())));
    disp(mean(A));
    disp(mean(A,1));
    disp(mean(A,2));
    disp(median(A));
    disp(median(A,1));
    disp(median(A,2));
    disp(min(A));
    disp(min(A,1));
    disp(min(A,2));
    disp(minimum(A-1,A));
    disp(minimum(A,5));
    disp(minimum(5,A));
    disp(norm(A));
    disp(norm(A,"inf"));
    disp(norm(A,"inf",2));
    disp(prod(A));
    disp(prod(A,1));
    disp(prod(A,2));
    disp(cumprod(A));
    disp(cumprod(A,1));
    disp(cumprod(A,2));
    disp(stddev(A));
    disp(stddev(A,1));
    disp(stddev(A,2));
    disp(sum(A));
    disp(sum(A,1));
    disp(sum(A,2));
    disp(sum<int>(cast(A,logical())));
    disp(cumsum(A));
    disp(cumsum(A,1));
    disp(cumsum(A,2));
    disp(diff(A));
    disp(diff(A,1));
    disp(diff(A,2));
    disp(variance(eye(2)));
    disp(variance(eye(2),1));
    disp(variance(eye(2),2));
    disp(sort(eye(4,3)));
    disp(sort(eye(4,3),1));
    disp(sort(eye(3,4),2));
    disp(cross(eye(4,3),ones(4,3)));
    disp(dot(eye(4,3),ones(4,3)));
    disp(unique(eye(3)));
    disp(intersect(eye(3),zeros(4)));
    disp(setdiff(eye(3),zeros(4)));
    disp(union2(eye(3),zeros(4)));
    disp(gmres(eye(3),ones(3,2)));
    disp(conv(ones(3,4),eye(3),2));
    disp(dft(eye(1,3)));
    disp(idft(ones<std::complex<double>>(1,3)));
    
    // =============================================================== FUNCTIONS
    std::cout << "+================================+" << std::endl;
    std::cout << "|     MATHEMATICAL FUNCTIONS     |" << std::endl;
    std::cout << "+================================+" << std::endl;

    V = {-M_PI,-2,-1,0,1,2,M_PI};
    std::cout << "original : " << V << std::endl;
    std::cout << "abs      : " << abs(V) << std::endl;
    std::cout << "acos     : " << acos(V) << std::endl;
    std::cout << "acosd    : " << acosd(V) << std::endl;
    std::cout << "acosh    : " << acosh(V) << std::endl;
    std::cout << "asin     : " << asin(V) << std::endl;
    std::cout << "asind    : " << asind(V) << std::endl;
    std::cout << "asinh    : " << asinh(V) << std::endl;
    std::cout << "atan     : " << atan(V) << std::endl;
    std::cout << "atand    : " << asind(V) << std::endl;
    std::cout << "atanh    : " << atanh(V) << std::endl;
    std::cout << "ceil     : " << ceil(V) << std::endl;
    std::cout << "cos      : " << cos(V) << std::endl;
    std::cout << "cosd     : " << cosd(V*180/M_PI) << std::endl;
    std::cout << "cosh     : " << cosh(V) << std::endl;
    std::cout << "deg2rad  : " << deg2rad(colon(0,90,360)) << std::endl;
    std::cout << "exp      : " << exp(V) << std::endl;
    std::cout << "floor    : " << floor(V) << std::endl;
    std::cout << "log      : " << log(V) << std::endl;
    std::cout << "log2     : " << log2(V) << std::endl;
    std::cout << "log10    : " << log10(V) << std::endl;
    std::cout << "pow      : " << pow(V,V) << std::endl;
    std::cout << "pow      : " << pow(V,2) << std::endl;
    std::cout << "pow      : " << pow(2,V) << std::endl;
    std::cout << "rad2deg  : " << rad2deg(V) << std::endl;
    std::cout << "round    : " << round(V) << std::endl;
    std::cout << "sign     : " << sign(V) << std::endl;
    std::cout << "sin      : " << sin(V) << std::endl;
    std::cout << "sind     : " << sind(V*180/M_PI) << std::endl;
    std::cout << "sinh     : " << sinh(V) << std::endl;
    std::cout << "sqrt     : " << sqrt(V) << std::endl;
    std::cout << "tan      : " << tan(V) << std::endl;
    std::cout << "tand     : " << tand(V*180/M_PI) << std::endl;
    std::cout << "tanh     : " << tanh(V) << std::endl;
    disp(abs<logical>(V));

    // =============================================================== COMPLEX
    std::cout << "+================+" << std::endl;
    std::cout << "|    COMPLEX     |" << std::endl;
    std::cout << "+================+" << std::endl;

    disp(M_1I);

    matrix<> Re = {1,0,-1,0};
    matrix<> Im = {0,1,0,-1};
    auto     Ac = Re + M_1I*Im;

    disp(Re);
    disp(Im);
    disp(Ac);
    disp(abs(Ac));
    disp(rad2deg(angle(Ac)));
    disp(conj(Ac));
    disp(exp(Ac*M_PI));
    disp(imag(Ac));
    disp(real(Ac));
    disp(sqrt(Ac));
    
    // =============================================================== ADVANCED
    std::cout << "+================+" << std::endl;
    std::cout << "|    GEOMETRY    |" << std::endl;
    std::cout << "+================+" << std::endl;

    matrix<> x = {1,-1,-1,1,1,-1,-1,1}, y={1,1,-1,-1,1,1,-1,-1}, z={1,1,1,1,-1,-1,-1,-1};
    matrix<> th, ph, r;
    double s_x=-1.0, s_y=1.0, s_z=1.0;
    double s_th=0.0, s_ph = 0.0, s_r = 0.0;

    // POLAR <-> CARTHESIAN
    disp(cat(1,x,y));
    std::tie(th,r) = cart2pol(x,y);
    disp(cat(1,rad2deg(th),r));
    std::tie(x,y) = pol2cart(th,r);
    disp(cat(1,x,y));
    disp(matrix<>({s_x,s_y}));
    std::tie(s_th,s_r) = cart2pol(s_x,s_y);
    disp(matrix<>({rad2deg(s_th),s_r}));
    std::tie(s_x,s_y) = pol2cart(s_th,s_r);
    disp(matrix<>({s_x,s_y}));

    // SPHERICAL <-> CARTHESIAN
    disp(cat(1,cat(1,x,y),z));
    std::tie(th,ph,r) = cart2sph(x,y,z);
    disp(cat(1,cat(1,rad2deg(th),rad2deg(ph)),r));
    std::tie(x,y,z) = sph2cart(th,ph,r);
    disp(cat(1,cat(1,x,y),z));
    disp(matrix<>({s_x,s_y,s_z}));
    std::tie(s_th,s_ph,s_r) = cart2sph(s_x,s_y,s_z);
    disp(matrix<>({rad2deg(s_th),rad2deg(s_ph),s_r}));
    std::tie(s_x,s_y,s_z) = sph2cart(s_th,s_ph,s_r);
    disp(matrix<>({s_x,s_y,s_z}));

    // SPHERE GENERATOR
    std::tie(x,y,z) = sphere(6);
    std::tie(th,ph,r) = cart2sph(x,y,z);
    disp(rad2deg(th));
    disp(rad2deg(ph));

    // SPHERICAL <-> DISCRETIZATION
    matrix<> azm, elv;
    idx = sph2idx(th,ph,6);
    std::tie(azm,elv) = idx2sph(idx,6);
    disp(idx);
    disp(rad2deg(azm));
    disp(rad2deg(elv));

    // VALIDATE BIJECTION
    std::tie(x,y,z)   = sphere(1e2,2e2);
    std::tie(th,ph,r) = cart2sph(x,y,z);
    idx               = sph2idx(th,ph,1e2,2e2);
    std::tie(azm,elv) = idx2sph(idx,1e2,2e2);
    matrix<std::size_t> ind = sph2idx(azm,elv,1e2,2e2);
    std::tie(azm,elv) = idx2sph(ind,1e2,2e2);
    disp( norm(azm-th,"inf") );
    disp( norm(elv-ph,"inf") );
    disp( max(ind-idx) );
    
    // FIBONACCI SPHERE GENERATOR
    disp("Sphere using fibonacci rules :");
    std::tie(x,y,z) = sphere2(6);
    disp(x);
    disp(y);
    disp(z);
    
    // MESH GRID
    disp("Mesh grid (2D) :");
    std::tie(x,y) = meshgrid(linspace(-1,1,10),linspace(0,2,5));
    disp(x);
    disp(y);
 
    // =============================================================== ADVANCED
    std::cout << "+================+" << std::endl;
    std::cout << "|    ADVANCED    |" << std::endl;
    std::cout << "+================+" << std::endl;
    
    int m=3, n=4;
    matrix<float> D = ones<float>(m,n);
    
    tgemm(M_PI,ones<int>(m,n),eye<logical>(n,n),2,D);
    disp(D);
    disp(kron(1,ones<float>(3,2)));
    disp(kron(eye<int>(2,3),1));
    disp(kron(eye<int>(2,3),ones<float>(3,2)));
    disp(mtimes(ones<int>(m,1),colon<float>(0,n)));
    disp(mtimes(eye(m,m),mtimes(ones(m,1),colon(0,n))));
    
    // =============================================================== PERFOS
    std::cout << "+================+" << std::endl;
    std::cout << "|     PERFOS     |" << std::endl;
    std::cout << "+================+" << std::endl;
    
    std::cout << "==> DGEMM :" << std::endl;
    matrix<float> AA(21,16);
    matrix<float> BB(16,512);
    matrix<float> CC(21,512);
    for (int i=0; i<10; i++)
    {
        AA = rand<float>(21,16);
        BB = rand<float>(16,512);
        tic();
        tgemm(1,AA,BB,0,CC);
        toc();
    }
    
    std::cout << "==> VIEWS :" << std::endl;
    tic();
    AA = rand<float>(1e3,1e3);
    toc();
    tic();
    BB = eval(AA(row(AA),0));
    toc();
    tic();
    BB = eval(AA(0,col(AA)));
    toc();
    tic();
    BB = eval(AA(row(AA),col(AA)));
    toc();
    
    
    // End of file
    std::cout << "done !" << std::endl;
    return 0;
}
