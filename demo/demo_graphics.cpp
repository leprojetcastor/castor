/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : demo_graphics.cpp                         |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Running test for VTK interface                |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include "castor/matrix.hpp"
#include "castor/smatrix.hpp"
#include "castor/graphics.hpp"

using namespace castor;

int main (int argc, char* argv[])
{
    //===============================================================
    std::cout << "+=====================+" << std::endl;
    std::cout << "|      INITIALIZE     |" << std::endl;
    std::cout << "+=====================+" << std::endl;

    // Documentation
    documentationFiles =
    {
        "/usr/local/include/castor/matrix.hpp",
        "/usr/local/include/castor/graphics.hpp"
    };
    help("plot");
    
    //===============================================================
    std::cout << "+=====================+" << std::endl;
    std::cout << "|         PLOT        |" << std::endl;
    std::cout << "+=====================+" << std::endl;
    
    figure fig;
    matrix<> X = linspace(0,10,100);
    matrix<> Y = cat(1,cos(X),sin(X));
    plot(fig,X,Y,{"k","r-"},{"cos(X)","sin(X)"});
    X = linspace(-2,2,30);
    plot(fig,X,sqrt(X),{"bx-"});
    plot(fig,X,X/2);
    plot(fig,{10},{1});
    xlim(fig,{-4,10});
    ylim(fig,{-2,2});
    
    //===============================================================
    std::cout << "+=====================+" << std::endl;
    std::cout << "|        PLOT3        |" << std::endl;
    std::cout << "+=====================+" << std::endl;

    X = linspace(-1,1,10);
    Y = X;
    matrix<> Z = X;
    figure fig2;
    plot3(fig2,X,Y,Z,"-r");
    plot3(fig2,X,zeros(size(Y)),zeros(size(Z)),"b");

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|        IMAGESC        |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    matrix<> M = linspace(-1,1,100);
    M = reshape(M,5,20);
    figure fig3;
    imagesc(fig3,M);

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|          SPY          |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    smatrix<> Ms = speye(10,15);
    figure fig31;
    spy(fig31,Ms,"b","Ms");

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|         MESH          |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    std::tie(X,Y) = meshgrid(linspace(-M_PI,M_PI,100));
    Z = 2 * sin(X)/X * sin(Y)/Y;

    figure fig4;
    mesh(fig4,Z);

    figure fig5;
    caxis(fig5,{-1,1});
    mesh(fig5,X,Y,Z);
    mesh(fig5,X,Y,-Z);

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|       VERMESH         |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    matrix<std::size_t> elt = range(0,100);
    matrix<> vtx = -1+2*rand(numel(elt),3);
    figure fig6;
    vermesh(fig6,elt,vtx,eval(vtx(row(vtx),0)));

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|       EDGMESH         |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    elt = transpose(range(0,100));
    elt = cat(2,elt,elt+1);
    elt(size(elt,1)-1,1) = 0;
    vtx = -1+2*rand(numel(elt),3);
    figure fig7;
    edgmesh(fig7,elt,vtx,eval(vtx(row(vtx),0)));

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|       TRIMESH         |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    // 2D grid with duplicate vertices
    std::tie(X,Y) = meshgrid(linspace(-1,1,10),linspace(-1,1,5));
    X = vertcat(X,X);
    Y = vertcat(Y,Y);
    Z = zeros(size(X));

    // Delaunay triangulation
    std::tie(elt,vtx) = tridelaunay(X,Y,Z);

    // Triangular mesh with vextex data
    figure fig8;
    trimesh(fig8,elt,vtx,eval(vtx(row(vtx),0)));

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|       TETMESH         |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    // 3D grid
    std::tie(X,Y) = meshgrid(linspace(-1,1,3));
    X = vertcat(X,X);
    Y = vertcat(Y,Y);
    Z = vertcat(zeros(3),ones(3));

    // Delaunay tetrahedrisation
    std::tie(elt,vtx) = tetdelaunay(X,Y,Z);

    // Tetrahedral mesh with vertex data
    figure fig9;
    tetmesh(fig9,elt,vtx,eval(vtx(row(vtx),0)));

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|    BOUNDARY MESH      |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    // Sphere grid
    std::tie(X,Y,Z) = sphere2(100);
    X = horzcat(0,reshape(X,1,numel(X)));
    Y = horzcat(0,reshape(Y,1,numel(Y)));
    Z = horzcat(0,reshape(Z,1,numel(Z)));

    // Delaunay tetrahedrisation
    std::tie(elt,vtx) = tetdelaunay(X,Y,Z);

    // Free boundary
    std::tie(elt,vtx) = tetboundary(elt,vtx);

    // Triangular mesh with vertex data
    figure fig10;
    trimesh(fig10,elt,vtx);

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|     TRIMESH I/O       |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    // Read .ply
    std::string path="./", name="testfile.ply";
    triwrite(path,name,elt,vtx);
    std::tie(elt,vtx) = triread(path,name);

    // Read .vtk
    name="testfile.vtk";
    triwrite(path,name,elt,vtx,rand(size(vtx,1),1));
    std::tie(elt,vtx) = triread(path,name);

    // Display
    figure fig11;
    trimesh(fig11,elt,vtx);

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|        QUIVER         |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    matrix<> dir = vtx;
    figure fig12;
    trimesh(fig12,elt,vtx);
    quiver(fig12,vtx,dir);

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|        DRAWNOW        |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    drawnow(fig);

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|     IMAGE WRITE       |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    std::vector<std::string> ext = {{""},{".png"},{".jpg"},{".ps"},
        {".tiff"},{".bmp"},{".pnm"}};
    for (int i=0; i<ext.size(); ++i)
    {
        writeimg(fig,"testfile"+ext[i]);
    }

    //===============================================================
    std::cout << "+=======================+" << std::endl;
    std::cout << "|     MOVIE WRITE       |" << std::endl;
    std::cout << "+=======================+" << std::endl;

    // Initialize source and movie
    vtkNew<vtkWindowToImageFilter> source;
    vtkNew<vtkOggTheoraWriter> movie;
    movie->SetInputConnection(source->GetOutputPort());
    movie->SetFileName("testfile.avi");
    movie->SetQuality(1); // in [0,2]
    movie->SetRate(25);   // frame per seconds

    // Recording frame by frame
    movie->Start();
    for (int i = 0; i < 50; i++)
    {
        figure fig13;
        matrix<> X = linspace(0,10,100);
        matrix<> Y = cos(X+i/10.);
        plot(fig13,X,Y,{"r-x"});
        source->SetInput(fig13.GetView()->GetRenderWindow());
        source->SetInputBufferTypeToRGB();
        source->ReadFrontBufferOff();
        movie->Write();
    }
    movie->End();
    
    disp("done !");
    return 0;
}
