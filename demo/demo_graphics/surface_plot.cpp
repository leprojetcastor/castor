/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : surface_plot.cpp                              |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Regular mesh gridding and plot as surface     |
 |  `---'  |                or wireframe                                  |
 +========================================================================+
 */

#include "castor/matrix.hpp"
#include "castor/graphics.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    matrix<> X,Y;
    std::tie(X,Y) = meshgrid(linspace(-M_PI,M_PI,100));
    matrix<> Z = 2*sin(X)/X * sin(Y)/Y;
    // create the figure
    figure fig;
    caxis(fig,{-1,1}); // scaled color in the range [-1,1]
    mesh(fig,X,Y,Z);
    // display
    drawnow(fig);
    return EXIT_SUCCESS;
}
