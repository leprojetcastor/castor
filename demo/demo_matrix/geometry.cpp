/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : geometry.cpp                                  |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : ND geometry tools                             |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include <castor/matrix.hpp>

using namespace castor;

//==============================================================================
int main (int argc, char* argv[])
{
    disp("+================+");
    disp("|    GEOMETRY    |");
    disp("+================+");

    matrix<std::size_t> idx;
    matrix<> x = {1,-1,-1,1,1,-1,-1,1}, y={1,1,-1,-1,1,1,-1,-1}, z={1,1,1,1,-1,-1,-1,-1};
    matrix<> x2, y2, z2;
    matrix<> th, ph, r;
    matrix<> azm, elv;

    double s_x=-1.0, s_y=1.0, s_z=1.0, s_x2=0.0, s_y2=0.0, s_z2=0.0;
    double s_th=0.0, s_ph = 0.0, s_r = 0.0;
    
    disp("Carthesian coordinates of the vertices of a 3D cube :");
    disp(x);
    disp(y);
    disp(z);
    
    disp("Convert carthesian coordinates to polar (degrees):");
    std::tie(th,r) = cart2pol(x,y);
    std::tie(s_th,s_r) = cart2pol(s_x,s_y);
    disp(rad2deg(th));
    disp(r);
    disp(rad2deg(s_th));
    disp(s_r);
    
    disp("Convert polar coordinates to carthesian (error) :");
    std::tie(x2,y2) = pol2cart(th,r);
    std::tie(s_x2,s_y2) = pol2cart(s_th,s_r);
    disp(max(cat(2,x-x2,y-y2)));
    disp(std::max(std::abs(s_x-s_x2),std::abs(s_y-s_y2)));
    
    disp("Convert carthesian coordinates to spherical (degrees) :");
    std::tie(th,ph,r) = cart2sph(x,y,z);
    std::tie(s_th,s_ph,s_r) = cart2sph(s_x,s_y,s_z);
    disp(rad2deg(th));
    disp(rad2deg(ph));
    disp(r);
    disp(rad2deg(s_th));
    disp(rad2deg(s_ph));
    disp(rad2deg(s_r));
    
    disp("Convert spherical coordinates to carthesian (error) :");
    std::tie(x2,y2,z2) = sph2cart(th,ph,r);
    std::tie(s_x2,s_y2,s_z2) = sph2cart(s_th,s_ph,s_r);
    disp(max(cat(2,cat(2,x-x2,y-y2),z-z2)));
    disp(std::max({std::abs(s_x-s_x2),std::abs(s_y-s_y2),std::abs(s_z-s_z2)}));
        
    disp("Azimut and elevation coordinates of a sphere :");
    std::tie(x,y,z) = sphere(6);
    std::tie(th,ph,r) = cart2sph(x,y,z);
    disp(rad2deg(th));
    disp(rad2deg(ph));

    disp("Conversion carthesian <-> spherical coordinates (error) :");
    std::tie(x2,y2,z2) = sph2cart(th,ph,r);
    disp(max(cat(2,cat(2,x-x2,y-y2),z-z2)));
    
    disp("Conversion spherical coordinates to linear indice :");
    idx = sph2idx(th,ph,6);
    disp(idx);
    
    disp("Conversion spherical coordinates <-> linear indice (error) :");
    std::tie(azm,elv) = idx2sph(idx,6);
    idx = sph2idx(azm,elv,6);
    std::tie(azm,elv) = idx2sph(idx,6);
    std::tie(x2,y2,z2) = sph2cart(azm,elv,ones(size(th)));
    disp(max(cat(2,cat(2,x-x2,y-y2),z-z2)));
    
    disp("Sphere using fibonacci rules :");
    std::tie(x,y,z) = sphere2(6);
    disp(x);
    disp(y);
    disp(z);
    
    disp("Mesh grid (2D) :");
    std::tie(x,y) = meshgrid(linspace(-1,1,10),linspace(0,2,5));
    disp(x);
    disp(y);
    
    disp("done !");
    return 0;
}
