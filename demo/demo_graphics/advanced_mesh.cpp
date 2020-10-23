/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : advanced_mesh.cpp                             |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Advanced use of 3D meshing and plot           |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/matrix.hpp"
#include "castor/graphics.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    std::size_t nvtx=100;
    // 1. Create the mesh
    matrix<> vtx;
    matrix<std::size_t> elt;
    
    matrix<> X,Y,Z;
    std::tie(X,Y,Z) = sphere2(nvtx); // Fibonacci sphere
    // Delaunay tetrahedrisation
    std::tie(elt,vtx) = tetdelaunay(X,Y,Z);
    // Boundary
    std::tie(elt,vtx) = tetboundary(elt,vtx);

    // 2. Compute the normal vectors and the centers of the triangles
    std::size_t nelt = size(elt,1);
    matrix<> nrm = zeros(nelt,3);
    matrix<> ctr = zeros(nelt,3);
    for(std::size_t ie=0; ie<nelt; ++ie)
    {
        // center
        for(std::size_t i=0; i<3; ++i)
        {
            ctr(ie,i) = (vtx(elt(ie,0),i)+vtx(elt(ie,1),i)+vtx(elt(ie,2),i))/3.;
        }
        // normal vector to triangle {A,B,C}
        matrix<> AB=zeros(1,3), BC=zeros(1,3), nv = zeros(1,3);
        for(std::size_t i=0; i<3; ++i)
        {
            AB(i) = vtx(elt(ie,1),i) - vtx(elt(ie,0),i);
            BC(i) = vtx(elt(ie,2),i) - vtx(elt(ie,1),i);
        }
        // for single vectors, this code is faster than
        // a call to 'cross' (for the cross-product) or 'norm'.
        nv(0) = AB(1)*BC(2) - AB(2)*BC(1);
        nv(1) = AB(2)*BC(0) - AB(0)*BC(2);
        nv(2) = AB(0)*BC(1) - AB(1)*BC(0);
        double l = std::sqrt(nv(0)*nv(0)+nv(1)*nv(1)+nv(2)*nv(2));
        for(std::size_t i=0; i<3; ++i) 
        {
            nrm(ie,i) = nv(i)/(2*l); // arrows of length 0.5
        }
    }

    // 3. Plot everything
    figure fig;
    trimesh(fig,elt,vtx);
    quiver(fig,ctr,nrm);
    drawnow(fig);

    // // 4. save to .ply
    // std::string path="./", name="testfile.ply";
    // triwrite(path,name,elt,vtx);

    // the end
    return EXIT_SUCCESS;
}
