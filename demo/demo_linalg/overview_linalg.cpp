/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : overview_linalg.cpp                           |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Linear algebra using templatized BLAS and     |
 |  `---'  |                Lapack interface for fast library             |
 +========================================================================+
 */

#include <castor/linalg.hpp>
#include <castor/matrix.hpp>

using namespace castor;

int main()
{
    //===============================================================
    std::cout << "+=====================+" << std::endl;
    std::cout << "|      INITIALIZE     |" << std::endl;
    std::cout << "+=====================+" << std::endl;
    
    // Documentation
    documentationFiles =
    {
        "/usr/local/include/castor/matrix.hpp",
        "/usr/local/include/castor/linalg.hpp"
    };
    help("svd");
    
    // Dimensions
    std::size_t m=100, n=200, k=30;
    
    //Random number generation
    time_t t;
    srand((unsigned) time(&t));
    
    //===============================================================
    std::cout << "+===================+" << std::endl;
    std::cout << "| CONVERSION LAPACK |" << std::endl;
    std::cout << "+===================+" << std::endl;
    
    // Conversion MATRIX
    matrix<double> A = rand(m,n);
    std::vector<float> V = mat2lpk<float>(A,numel(A));
    matrix<double> B = lpk2mat<double>(V,size(A,1),size(A,2));
    disp( norm(A-B,"inf") );
    
    // Conversion COMPLEX FLOAT
    clpk cl = std::complex<float>(1,1);
    cl.r = 2;
    cl.i = 3;
    disp(cl);
    std::complex<float> cc = cl;
    disp(cc);
    
    // Conversion COMPLEX DOUBLE
    zlpk zl = std::complex<double>(1,1);
    zl.r = 4;
    zl.i = 5;
    disp(zl);
    std::complex<double> zc = zl;
    disp(zc);
    
    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|   MULTIPLY    |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // FLOAT
    matrix<float> Af = rand<float>(1,m);
    matrix<float> Bf = rand<float>(m,1);
    float Cf = 0;
    for (std::size_t l=0; l<m; ++l) {Cf+=Af(l)*Bf(l);}
    disp( norm(mtimes(Af,Bf)-Cf,"inf") );
    
    // DOUBLE
    matrix<> Ad = rand(1,m);
    matrix<> Bd = rand(m,1);
    double Cd = 0;
    for (std::size_t l=0; l<m; ++l) {Cd+=Ad(l)*Bd(l);}
    disp( norm(mtimes(Ad,Bd)-Cd,"inf") );
    
    // COMPLEX FLOAT
    auto Ac = rand<float>(1,m) + M_1If*rand<float>(1,m);
    auto Bc = rand<float>(m,1) + M_1If*rand<float>(m,1);
    std::complex<float> Cc = 0;
    for (std::size_t l=0; l<m; ++l) {Cc+=Ac(l)*Bc(l);}
    disp( norm(mtimes(Ac,Bc)-Cc,"inf") );
    
    // COMPLEX DOUBLE
    auto Az = rand(1,m) + M_1I*rand(1,m);
    auto Bz = rand(m,1) + M_1I*rand(m,1);
    std::complex<double> Cz = 0;
    for (std::size_t l=0; l<m; ++l) {Cz+=Az(l)*Bz(l);}
    disp( norm(mtimes(Az,Bz)-Cz,"inf") );
    
    //===============================================================
    std::cout << "+==============+" << std::endl;
    std::cout << "|   LINSOLVE   |" << std::endl;
    std::cout << "+==============+" << std::endl;
    
    // FLOAT
    Af = rand<float>(m);
    Bf = rand<float>(m,2);
    matrix<float> Xf = linsolve(Af,Bf);
    disp( norm(mtimes(Af,Xf)-Bf,"inf") );
    
    Af = cat(1,Af,Af);
    Bf = cat(1,Bf,Bf);
    Xf = linsolve(Af,Bf);
    disp( norm(mtimes(Af,Xf)-Bf,"inf") );
    
    Af = rand<float>(m,n);
    Bf = rand<float>(m,2);
    Xf = linsolve(Af,Bf);
    disp( norm(mtimes(Af,Xf)-Bf,"inf") );
    
    // DOUBLE
    Ad = rand(m);
    Bd = rand(m,2);
    matrix<> Xd = linsolve(Ad,Bd);
    disp( norm(mtimes(Ad,Xd)-Bd,"inf") );
    
    Ad = cat(1,Ad,Ad);
    Bd = cat(1,Bd,Bd);
    Xd = linsolve(Ad,Bd);
    disp( norm(mtimes(Ad,Xd)-Bd,"inf") );
    
    Ad = rand(m,n);
    Bd = rand(m,2);
    Xd = linsolve(Ad,Bd);
    disp( norm(mtimes(Ad,Xd)-Bd,"inf") );
    
    // COMPLEX FLOAT
    Ac = rand<float>(m) + M_1If*rand<float>(m);
    Bc = rand<float>(m,2) + M_1If*rand<float>(m,2);
    auto Xc = linsolve(Ac,Bc);
    disp( norm(mtimes(Ac,Xc)-Bc,"inf") );
    
    Ac = cat(1,Ac,Ac);
    Bc = cat(1,Bc,Bc);
    Xc = linsolve(Ac,Bc);
    disp( norm(mtimes(Ac,Xc)-Bc,"inf") );
    
    Ac = rand<float>(m,n) + M_1If*rand<float>(m,n);
    Bc = rand<float>(m,2) + M_1If*rand<float>(m,2);
    Xc = linsolve(Ac,Bc);
    disp( norm(mtimes(Ac,Xc)-Bc,"inf") );
    
    // COMPLEX DOUBLE
    Az = rand(m) + M_1I*rand(m);
    Bz = rand(m,2) + M_1I*rand(m,2);
    auto Xz = linsolve(Az,Bz);
    disp( norm(mtimes(Az,Xz)-Bz,"inf") );
    
    Az = cat(1,Az,Az);
    Bz = cat(1,Bz,Bz);
    Xz = linsolve(Az,Bz);
    disp( norm(mtimes(Az,Xz)-Bz,"inf") );
    
    Az = rand(m,n) + M_1I*rand(m,n);
    Bz = rand(m,2) + M_1I*rand(m,2);
    Xz = linsolve(Az,Bz);
    disp( norm(mtimes(Az,Xz)-Bz,"inf") );

    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|      INV      |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // FLOAT
    Af = rand<float>(m);
    Xf = inv(Af);
    disp( norm(mtimes(Af,Xf)-eye<float>(m),"inf") );
    
    // DOUBLE
    Ad = rand(m);
    Xd = inv(Ad);
    disp( norm(mtimes(Ad,Xd)-eye(m),"inf") );
    
    // COMPLEX FLOAT
    Ac = rand<float>(m) + M_1If*rand<float>(m);
    Xc = inv(Ac);
    disp( norm(mtimes(Ac,Xc)-eye<float>(m),"inf") );
    
    // COMPLEX DOUBLE
    Az = rand(m) + M_1I*rand(m);
    Xz = inv(Az);
    disp( norm(mtimes(Az,Xz)-eye(m),"inf") );
    
    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|     PINV      |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // FLOAT
    Af = rand<float>(m,n);
    Xf = pinv(Af);
    disp(norm(mtimes(mtimes(Xf,Af),Xf)-Xf,"inf"));
    disp(norm(mtimes(mtimes(Af,Xf),Af)-Af,"inf"));
    
    Af = rand<float>(n,m);
    Xf = pinv(Af);
    disp(norm(mtimes(mtimes(Xf,Af),Xf)-Xf,"inf"));
    disp(norm(mtimes(mtimes(Af,Xf),Af)-Af,"inf"));
    
    // DOUBLE
    Ad = rand(m,n);
    Xd = pinv(Ad);
    disp(norm(mtimes(mtimes(Xd,Ad),Xd)-Xd,"inf"));
    disp(norm(mtimes(mtimes(Ad,Xd),Ad)-Ad,"inf"));
    
    Ad = rand(n,m);
    Xd = pinv(Ad);
    disp(norm(mtimes(mtimes(Xd,Ad),Xd)-Xd,"inf"));
    disp(norm(mtimes(mtimes(Ad,Xd),Ad)-Ad,"inf"));
    
    // COMPLEX FLOAT
    Ac = rand<float>(m,n) + M_1If*rand<float>(m,n);
    Xc = pinv(Ac);
    disp(norm(mtimes(mtimes(Xc,Ac),Xc)-Xc,"inf"));
    disp(norm(mtimes(mtimes(Ac,Xc),Ac)-Ac,"inf"));
    
    Ac = rand<float>(n,m) + M_1If*rand<float>(n,m);
    Xc = pinv(Ac);
    disp(norm(mtimes(mtimes(Xc,Ac),Xc)-Xc,"inf"));
    disp(norm(mtimes(mtimes(Ac,Xc),Ac)-Ac,"inf"));
    
    // COMPLEX DOUBLE
    Az = rand(m,n) + M_1I*rand(m,n);
    Xz = pinv(Az);
    disp(norm(mtimes(mtimes(Xz,Az),Xz)-Xz,"inf"));
    disp(norm(mtimes(mtimes(Az,Xz),Az)-Az,"inf"));
    
    Az = rand(n,m) + M_1I*rand(n,m);
    Xz = pinv(Az);
    disp(norm(mtimes(mtimes(Xz,Az),Xz)-Xz,"inf"));
    disp(norm(mtimes(mtimes(Az,Xz),Az)-Az,"inf"));
    
    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|       LU      |" << std::endl;
    std::cout << "+===============+" << std::endl;
    
    // FLOAT
    matrix<float> Lf, Uf;
    Af = rand<float>(m,n);
    std::tie(Lf,Uf) = lu(Af);
    disp( norm(mtimes(Lf,Uf)-Af,"inf") );
    
    Af = rand<float>(n,m);
    std::tie(Lf,Uf) = lu(Af);
    disp( norm(mtimes(Lf,Uf)-Af,"inf") );
    
    // DOUBLE
    matrix<> Ld, Ud;
    Ad = rand(m,n);
    std::tie(Ld,Ud) = lu(Ad);
    disp( norm(mtimes(Ld,Ud)-Ad,"inf") );
    
    Ad = rand(n,m);
    std::tie(Ld,Ud) = lu(Ad);
    disp( norm(mtimes(Ld,Ud)-Ad,"inf") );
    
    // COMPLEX FLOAT
    matrix<std::complex<float>> Lc, Uc;
    Ac = rand<float>(m,n) + M_1If*rand<float>(m,n);
    std::tie(Lc,Uc) = lu(Ac);
    disp( norm(mtimes(Lc,Uc)-Ac,"inf") );
    
    Ac = rand<float>(n,m) + M_1If*rand<float>(n,m);
    std::tie(Lc,Uc) = lu(Ac);
    disp( norm(mtimes(Lc,Uc)-Ac,"inf") );
    
    // COMPLEX DOUBLE
    matrix<std::complex<double>> Lz, Uz;
    Az = rand(m,n) + M_1I*rand(m,n);
    std::tie(Lz,Uz) = lu(Az);
    disp( norm(mtimes(Lz,Uz)-Az,"inf") );

    Az = rand(n,m) + M_1I*rand(n,m);
    std::tie(Lz,Uz) = lu(Az);
    disp( norm(mtimes(Lz,Uz)-Az,"inf") );
    
    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|      EIG      |" << std::endl;
    std::cout << "+===============+" << std::endl;

    // FLOAT
    Af = 1-eye<float>(m);
    matrix<float> Ef, Vf;
    std::tie(Ef,Uf) = eig(Af,"left");
    Uf = transpose(Uf);
    disp(norm(mtimes(Uf,Af)-mtimes(diag(Ef),Uf),"inf"));
    std::tie(Ef,Vf) = eig(Af,"right");
    disp(norm(mtimes(Af,Vf)-mtimes(Vf,diag(Ef)),"inf"));

    // DOUBLE
    Ad = 1-eye(m);
    matrix<> Ed, Vd;
    std::tie(Ed,Ud) = eig(Ad,"left");
    Ud = transpose(Ud);
    disp(norm(mtimes(Ud,Ad)-mtimes(diag(Ed),Ud),"inf"));
    std::tie(Ed,Vd) = eig(Ad,"right");
    disp(norm(mtimes(Ad,Vd)-mtimes(Vd,diag(Ed)),"inf"));

    // COMPLEX FLOAT
    Ac = rand<float>(m);
    matrix<std::complex<float>> Ec, Vc;
    std::tie(Ec,Uc) = eig(Ac,"left");
    Uc = conj(transpose(Uc));
    disp(norm(mtimes(Uc,Ac)-mtimes(diag(Ec),Uc),"inf"));
    std::tie(Ec,Vc) = eig(Ac,"right");
    disp(norm(mtimes(Ac,Vc)-mtimes(Vc,diag(Ec)),"inf"));

    // COMPLEX DOUBLE
    Az = rand(m);
    matrix<std::complex<double>> Ez, Vz;
    std::tie(Ez,Uz) = eig(Az,"left");
    Uz = conj(transpose(Uz));
    disp(norm(mtimes(Uz,Az)-mtimes(diag(Ez),Uz),"inf"));
    std::tie(Ez,Vz) = eig(Az,"right");
    disp(norm(mtimes(Az,Vz)-mtimes(Vz,diag(Ez)),"inf"));

    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|      QR       |" << std::endl;
    std::cout << "+===============+" << std::endl;

    // FLOAT
    Af = rand<float>(m,n);
    matrix<float> Qf, Rf;
    std::tie(Qf,Rf) = qr(Af);
    disp(norm(Af-mtimes(Qf,Rf),"inf"));

    Af = rand<float>(n,m);
    std::tie(Qf,Rf) = qr(Af);
    disp(norm(Af-mtimes(Qf,Rf),"inf"));

    // DOUBLE
    Ad = rand(m,n);
    matrix<double> Qd, Rd;
    std::tie(Qd,Rd) = qr(Ad);
    disp(norm(Ad-mtimes(Qd,Rd),"inf"));

    Ad = rand(n,m);
    std::tie(Qd,Rd) = qr(Ad);
    disp(norm(Ad-mtimes(Qd,Rd),"inf"));

    // COMPLEX FLOAT
    Ac = rand<float>(m,n)*M_1If;
    matrix<std::complex<float>> Qc, Rc;
    std::tie(Qc,Rc) = qr(Ac);
    disp(norm(Ac-mtimes(Qc,Rc),"inf"));

    Ac = rand<float>(n,m)*M_1If;
    std::tie(Qc,Rc) = qr(Ac);
    disp(norm(Ac-mtimes(Qc,Rc),"inf"));

    // COMPLEX DOUBLE
    Az = rand(m,n)*M_1I;
    matrix<std::complex<double>> Qz, Rz;
    std::tie(Qz,Rz) = qr(Az);
    disp(norm(Az-mtimes(Qz,Rz),"inf"));

    Az = rand(n,m)*M_1I;
    std::tie(Qz,Rz) = qr(Az);
    disp(norm(Az-mtimes(Qz,Rz),"inf"));

    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|      SVD      |" << std::endl;
    std::cout << "+===============+" << std::endl;

    // FLOAT
    Af = rand<float>(m,n);
    matrix<float> Sf;
    std::tie(Sf,Uf,Vf) = svd(Af,"vect");
    disp(norm(Af-mtimes(Uf,mtimes(diag(Sf),Vf)),"inf"));

    Af = rand<float>(n,m);
    std::tie(Sf,Uf,Vf) = svd(Af,"vect");
    disp(norm(Af-mtimes(Uf,mtimes(diag(Sf),Vf)),"inf"));

    // DOUBLE
    Ad = rand(m,n);
    matrix<> Sd;
    std::tie(Sd,Ud,Vd) = svd(Ad,"vect");
    disp(norm(Ad-mtimes(Ud,mtimes(diag(Sd),Vd)),"inf"));

    Ad = rand(n,m);
    std::tie(Sd,Ud,Vd) = svd(Ad,"vect");
    disp(norm(Ad-mtimes(Ud,mtimes(diag(Sd),Vd)),"inf"));

    // COMPLEX FLOAT
    Ac = rand<float>(m,n)*M_1If;
    std::tie(Sf,Uc,Vc) = svd(Ac,"vect");
    disp(norm(Ac-mtimes(Uc,mtimes(diag(Sf),Vc)),"inf"));

    Ac = rand<float>(n,m)*M_1If;
    std::tie(Sf,Uc,Vc) = svd(Ac,"vect");
    disp(norm(Ac-mtimes(Uc,mtimes(diag(Sf),Vc)),"inf"));

    // COMPLEX DOUBLE
    Az = rand(m,n)*M_1I;
    std::tie(Sd,Uz,Vz) = svd(Az,"vect");
    disp(norm(Az-mtimes(Uz,mtimes(diag(Sd),Vz)),"inf"));

    Az = rand(n,m)*M_1I;
    std::tie(Sd,Uz,Vz) = svd(Az,"vect");
    disp(norm(Az-mtimes(Uz,mtimes(diag(Sd),Vz)),"inf"));

    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|      RANK     |" << std::endl;
    std::cout << "+===============+" << std::endl;

    // FLOAT
    Af = mtimes(rand<float>(m,k),rand<float>(k,n));
    disp(rank(Af));

    // DOUBLE
    Ad = mtimes(rand(m,k),rand(k,n));
    disp(rank(Ad));

    // COMPLEX FLOAT
    Ac = mtimes(rand<float>(m,k),rand<float>(k,n))*M_1If;
    disp(rank(Ac));

    // COMPLEX DOUBLE
    Az = mtimes(rand(m,k),rand(k,n))*M_1I;
    disp(rank(Az));

    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|      ACA      |" << std::endl;
    std::cout << "+===============+" << std::endl;

    // DOUBLE
    bool flag;
    matrix<> Md = mtimes(rand(m,k),rand(k,n));
    std::tie(Ad,Bd,flag) = aca(Md,1e-3);
    disp(size(Ad));
    disp(norm(Md-mtimes(Ad,Bd),"inf"));

    // DOUBLE
    matrix<std::size_t> I=range(0,m/2), J=range(0,n/2);
    std::function<matrix<double>(matrix<std::size_t>,matrix<std::size_t>)> Fd;
    Fd = [&Md](matrix<std::size_t> I, matrix<std::size_t> J)
    {
        return eval(Md(I,J));
    };
    std::tie(Ad,Bd,flag) = aca(I,J,Fd,1e-3);
    disp(size(Ad));
    disp(norm(eval(Md(I,J))-mtimes(Ad,Bd),"inf"));

    // COMPLEX DOUBLE
    matrix<std::complex<double>> Mz = mtimes(rand(m,k),rand(k,n))*M_1I;
    std::tie(Az,Bz,flag) = aca(Mz,1e-3);
    disp(size(Az));
    disp(norm(Mz-mtimes(Az,Bz),"inf"));
        
    // COMPLEX DOUBLE
    std::function<matrix<std::complex<double>>(matrix<std::size_t>,matrix<std::size_t>)> Fz;
    Fz = [&Mz](matrix<std::size_t> I, matrix<std::size_t> J)
    {
        return eval(Mz(I,J));
    };
    std::tie(Az,Bz,flag) = aca(I,J,Fz,1e-3);
    disp(size(Az));
    disp(norm(eval(Mz(I,J))-mtimes(Az,Bz),"inf"));
    
    //===============================================================
    std::cout << "+===============+" << std::endl;
    std::cout << "|     QRSVD     |" << std::endl;
    std::cout << "+===============+" << std::endl;

    // FLOAT
    Af = rand(m,k); Bf = rand(k,n);
    std::tie(Uf,Vf) = qrsvd(Af,Bf,1e-3);
    disp(norm(mtimes(Af,Bf)-mtimes(Uf,Vf),"inf"));

    // DOUBLE
    Ad = rand(m,k); Bd = rand(k,n);
    std::tie(Ud,Vd) = qrsvd(Ad,Bd,1e-3);
    disp(norm(mtimes(Ad,Bd)-mtimes(Ud,Vd),"inf"));

    // COMPLEX FLOAT
    Ac = rand<float>(m,k); Bc = rand<float>(k,n);
    std::tie(Uc,Vc) = qrsvd(Ac,Bc,1e-3);
    disp(norm(mtimes(Ac,Bc)-mtimes(Uc,Vc),"inf"));

    // COMPLEX DOUBLE
    Az = rand(m,k); Bz = rand(k,n);
    std::tie(Uz,Vz) = qrsvd(Az,Bz,1e-3);
    disp(norm(mtimes(Az,Bz)-mtimes(Uz,Vz),"inf"));
    
    // End of file
    std::cout << "done !" << std::endl;
    return 0;
}
