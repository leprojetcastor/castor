#include <iostream>
#include <complex>
#include <functional>

#include "castor/matrix.hpp"
#include "castor/graphics.hpp"

using namespace castor;

struct FractalPars
{
    double xmin{0.}, ymin{0.};
    double zoom{0};
    std::size_t nx{0}, ny{0};
    std::size_t itMax{0};
    double divcrit{2.};
};

matrix<std::size_t> julia(FractalPars const &pars, std::complex<double> c);
matrix<std::size_t> mandelbrot(FractalPars const &pars);
matrix<std::size_t> generalized_mandelbrot(FractalPars const &pars, std::function<std::complex<double>(std::complex<double>, std::complex<double>)> const &_func);
FractalPars getPars(short _fracType); // PRESET PARAMETERS FOR THIS EXAMPLE

int main()
{
    // GENERAL SETTINGS

    std::complex<double> c(0.285,0.01); // ONLY RELEVANT FOR JULIA

    // CREATE PLOT
    std::vector<figure> figs;
    for(short fractal_type: {0,1,2}) // 0: JULIA,  1: MANDELBROT, 2: GENERALIZED MANDELBROT
    {
        matrix<std::size_t> fractal;
        auto pars = getPars(fractal_type);
        tic();
        switch(fractal_type)
        {
            case 0:
                fractal = julia(pars, c);
                break;
            case 1:
                fractal = mandelbrot(pars);
                break;
            case 2:
                fractal = generalized_mandelbrot(pars, [](std::complex<double> _z, std::complex<double> _c) -> std::complex<double>{return std::cos(_z) + 1./_c;});
                break;
            default:
                return EXIT_FAILURE;
        }
        toc();
        // INVERT COLORS
        fractal = pars.itMax - fractal;
        // ADD SOME PLOT
        figs.push_back(figure());
        imagesc<float>(figs.back(),fractal);
    }
    drawnow(figs.back());
    return EXIT_SUCCESS;
}

FractalPars getPars(short fracType)
{
    FractalPars pars;
    double xmax, ymax;
    switch(fracType)
    {
        case 0:
            pars.xmin   = -1.;
            xmax = 1.;
            pars.ymin   = -1.2;
            ymax = 1.2;
            pars.zoom   = 500;
            pars.itMax  = 150;
            pars.nx     = static_cast<std::size_t>((xmax - pars.xmin)*pars.zoom);
            pars.ny     = static_cast<std::size_t>((ymax - pars.ymin)*pars.zoom);
            pars.divcrit = 2.;
            break;
        case 1:
            pars.xmin   = -2.1;
            xmax = 0.6;
            pars.ymin   = -1.2;
            ymax = 1.2;
            pars.zoom   = 500;
            pars.itMax  = 150;
            pars.nx     = static_cast<std::size_t>((xmax - pars.xmin)*pars.zoom);
            pars.ny     = static_cast<std::size_t>((ymax - pars.ymin)*pars.zoom);
            pars.divcrit = 2.;
            break;
        case 2:
            pars.xmin   = -1;
            xmax = 1;
            pars.ymin   = -1;
            ymax = 1;
            pars.zoom   = 500;
            pars.itMax  = 50;
            pars.nx     = static_cast<std::size_t>((xmax - pars.xmin)*pars.zoom);
            pars.ny     = static_cast<std::size_t>((ymax - pars.ymin)*pars.zoom);
            pars.divcrit = 1024*1024;
            break;
        default:
            std::exit(EXIT_FAILURE);
    }
    return pars;
}

matrix<std::size_t> julia(FractalPars const &pars, std::complex<double> c)
{
    std::cout << "+---------------------------+" << std::endl;
    std::cout << "| Computing Julia ensembles |" << std::endl;
    std::cout << "+---------------------------+" << std::endl;
    matrix<std::size_t> frac = zeros(pars.ny,pars.nx);
    //
    for(auto iy = 0; iy<pars.ny; ++iy)
    {
        for(auto ix = 0; ix<pars.nx; ++ix)
        {
            auto iy_ = pars.ny - iy - 1;
            auto z = std::complex<double>(ix/pars.zoom+pars.xmin, iy/pars.zoom + pars.ymin);
            //
            do
            {
                z = z*z + c;
                ++frac(iy_,ix);
            } while (std::abs(z)<pars.divcrit && frac(iy_,ix) < pars.itMax);
        }
    }
    return frac;
}

matrix<std::size_t> mandelbrot(FractalPars const &pars)
{
    std::cout << "+----------------------------------+" << std::endl;
    std::cout << "| Computing the Mandelbrot fractal |" << std::endl;
    std::cout << "+----------------------------------+" << std::endl;
    matrix<std::size_t> frac = zeros(pars.ny,pars.nx);
    //
    for(auto iy = 0; iy<pars.ny; ++iy)
    {
        for(auto ix = 0; ix<pars.nx; ++ix)
        {
            auto iy_ = pars.ny - iy - 1;
            auto z = std::complex<double>(0,0);
            auto c = std::complex<double>(ix/pars.zoom+pars.xmin, iy/pars.zoom + pars.ymin);
            //
            do
            {
                z = z*z + c;
                ++frac(iy_,ix);
            } while (std::abs(z)<pars.divcrit && frac(iy_,ix) < pars.itMax);
        }
    }
    return frac;
}

matrix<std::size_t> generalized_mandelbrot(FractalPars const &pars, std::function<std::complex<double>(std::complex<double>, std::complex<double>)> const &func)
{
    // JUST LIKE MANDELBROT, BUT WE REPLACED z <- z*z + c BY SOME OTHER FUNCTION
    std::cout << "+--------------------------------------------+" << std::endl;
    std::cout << "| Computing a Generalized Mandelbrot fractal |" << std::endl;
    std::cout << "+--------------------------------------------+" << std::endl;
    matrix<std::size_t> frac = zeros(pars.ny,pars.nx);
    //
    for(auto iy = 0; iy<pars.ny; ++iy)
    {
        for(auto ix = 0; ix<pars.nx; ++ix)
        {
            auto iy_ = pars.ny - iy - 1;
            auto z = std::complex<double>(0,0);
            auto c = std::complex<double>(ix/pars.zoom+pars.xmin, iy/pars.zoom + pars.ymin);
            //
            do
            {
                z = func(z,c);
                ++frac(iy_,ix);
            } while (std::abs(z)<pars.divcrit && frac(iy_,ix) < pars.itMax);
        }
    }
    return frac;
}
