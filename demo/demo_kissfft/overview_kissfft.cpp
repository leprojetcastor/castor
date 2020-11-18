/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : overview_kissfft.cpp                          |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 25.11.2020                                    |
 |  / 0 \  |   LAST MODIF : 25.11.2020                                    |
 | ( === ) |   SYNOPSIS   : Interface of kissfft library for complex and  |
 |  `---'  |                real fourier transform.                       |
 +========================================================================+
 */

#include <castor/matrix.hpp>
#include <castor/kissfft.hpp>

using namespace castor;

int main(int argc, const char * argv[])
{
    disp("+===============+");
    disp("|     CKISS     |");
    disp("+===============+");
    
    ckiss a(3,4);
    ckiss b(3,-4);
    a *= b;
    disp(a);
    
    disp("+===============+");
    disp("|      FFT      |");
    disp("+===============+");
    
    matrix<std::complex<float>> A = rand(3,4) + M_1I*rand(3,4), B;
    
    // FFT for all values
    disp("FFT for all values :");
    B = ifft(fft(A));
    disp( norm(eval(A(range(0,12)))-B,"inf") );
    B = ifft(fft(A,6),6);
    disp( norm(eval(A(range(0,6)))-B,"inf") );
    B = ifft(fft(A,24));
    disp( norm(eval(A(range(0,12)))-eval(B(range(0,12))),"inf") );

    // FFT along first dimension
    disp("FFT WITH along 1st dimension :");
    B = ifft(fft(A,0,1),0,1);
    disp( norm(A-B,"inf") );
    B = ifft(fft(A,2,1),2,1);
    disp( norm(eval(A({0,1},col(A)))-B,"inf") );
    B = ifft(fft(A,6,1),6,1);
    disp( norm(A-eval(B(row(A),col(A))),"inf") );
    
    // FFT along second dimension
    disp("FFT WITH along 2nd dimension:");
    B = ifft(fft(A,0,2),0,2);
    disp( norm(A-B,"inf") );
    B = ifft(fft(A,2,2),2,2);
    disp( norm(eval(A(row(A),{0,1}))-B,"inf") );
    B = ifft(fft(A,6,2),6,2);
    disp( norm(A-eval(B(row(A),col(A))),"inf") );
    
    disp("+===============+");
    disp("|      DFT      |");
    disp("+===============+");
    
    matrix<std::complex<float>> ref, sol;
    
    // FFT versus DFT for all values
    disp("DFT vs FFT for all values :");
    A = rand(30,30);
    tic();
    ref = dft(A);
    toc();
    tic();
    sol = fft(A);
    toc();
    disp( norm(ref-sol,"inf") / norm(ref,"inf") );
    
    // FFT versus DFT for first dimension
    disp("DFT vs FFT along 1st dimension :");
    ref = dft(A,0,1);
    sol = fft(A,0,1);
    disp( norm(ref-sol,"inf") / norm(ref,"inf") );

    // FFT versus DFT for second dimension
    disp("DFT vs FFT along 2nd dimension :");
    ref = dft(A,0,2);
    sol = fft(A,0,2);
    disp( norm(ref-sol,"inf") / norm(ref,"inf") );

    // DFT AND IDFT
    disp("Discrete fourier transform forward and backward : ");
    ref = A;
    sol = reshape(idft(dft(A)),size(A,1),size(A,2));
    disp( norm(ref-sol,"inf") / norm(ref,"inf") );
    
    disp("+===============+");
    disp("|     PERFO     |");
    disp("+===============+");
    
    float time = 0;
    A = rand(1,1024) + M_1I*rand(1,1024);

    // Direct Fourier Transform
    disp("FFT forward (s) :");
    for (int h=0; h<100; ++h)
    {
        tic();
        B = fft(A);
        time += toc(0);
    }
    disp(time/100);
    
    // Inverse Fourier Transform
    disp("FFT backward (s) :");
    time = 0;
    for (int h=0; h<100; ++h)
    {
        tic();
        B = ifft(A);
        time += toc(0);
    }
    disp(time/100);
    
    
    disp("Done !");
    return 0;
}
