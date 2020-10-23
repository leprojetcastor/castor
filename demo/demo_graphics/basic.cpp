/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : basic.cpp                                     |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Basic 2D plot with style and legend           |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/matrix.hpp"
#include "castor/graphics.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    matrix<> X = linspace(0,10,100);
    matrix<> Y = cos(X);
    figure fig;
    plot(fig,X,Y,{"r-d"},{"sin(x)"});
    X = linspace(-2,2,30);
    Y = sqrt(X);
    plot(fig,X,Y,{"b"},{"sqrt(x)"});
    // display
    drawnow(fig);
    return EXIT_SUCCESS;
}
