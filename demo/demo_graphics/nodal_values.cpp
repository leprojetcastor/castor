/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : nodal_values.cpp                              |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Nodal values on triangular mesh               |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/matrix.hpp"
#include "castor/graphics.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    // geometric data
    matrix<> X,Y;
    std::tie(X,Y) = meshgrid(linspace(-1,1,10),linspace(-1,1,5));
    X = vertcat(X,X); 
    Y = vertcat(Y,Y);
    matrix<> Z = zeros(size(X));

    // create mesh
    matrix<> vtx;
    matrix<std::size_t> elt;
    std::tie(elt,vtx) = tridelaunay(X,Y,Z);

    // display
    figure fig;
    trimesh(fig,elt,vtx,eval(vtx(row(vtx),0)));
    drawnow(fig);
    return EXIT_SUCCESS;
}
