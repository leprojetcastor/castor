/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : graphics.hpp                                  |
 |    #    |   VERSION    : 1.0.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Graphical representation based on VTK library |
 |  `---'  |                                                              |
 +========================================================================+
 */

#pragma once

#define CASTOR_GRAPHICS_HPP

#include <castor/matrix.hpp>
#include <castor/smatrix.hpp>

#include <vtkArrowSource.h>
#include <vtkAxesActor.h>
#include <vtkAxis.h>
#include <vtkBMPWriter.h>
#include <vtkCellData.h>
#include <vtkChartLegend.h>
#include <vtkChartXY.h>
#include <vtkChartXYZ.h>
#include <vtkCleanPolyData.h>
#include <vtkConeSource.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkImageWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkLine.h>
#include <vtkLookupTable.h>
#include <vtkMapper.h>
#include <vtkNew.h>
#include <vtkOggTheoraWriter.h>
#include <vtkPen.h>
#include <vtkPlaneSource.h>
#include <vtkPlotLine.h>
#include <vtkPlotLine3D.h>
#include <vtkPlotPoints.h>
#include <vtkPlotPoints3D.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkPNGWriter.h>
#include <vtkPNMWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPostScriptWriter.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkTetra.h>
#include <vtkTextProperty.h>
#include <vtkTIFFWriter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangle.h>
#include <vtkVertex.h>
#include <vtkUnstructuredGrid.h>
#include <vtkWindow.h>
#include <vtkWindowToImageFilter.h>


namespace castor
{

//==========================================================================//
//                             MASTER CLASS                                 //
//==========================================================================//
// [figure]
/// Container to display 2- or 3-dimensional plots
///
/// Data is added to the figure by passing it as an argument to the
/// corresponding functions. Figure are fully interactable.
class figure
{
public:
    figure(int w=0, int h=0);
    
    void axes();
    void caxis(matrix<double>const& val);
    void colorbar(vtkMapper* mapper);
    void drawnow();
    void interactor();
    void xlim(matrix<double>const& limits);
    void ylim(matrix<double>const& limits);

    template<typename T>
    void edgmesh(matrix<std::size_t>const& edg, matrix<T>const& vtx, matrix<T>const& val);
    template<typename T>
    void imagesc(matrix<T>const& M);
    template<typename T>
    void mesh(matrix<T>const& X, matrix<T>const& Y, matrix<T>const& Z, std::string const& options);
    template<typename T=double>
    void plot(matrix<T>const& X, matrix<T>const&Y, std::vector<std::string>const& style,
              std::vector<std::string>const& label);
    template<typename T=double>
    void plot3(matrix<T>const& X, matrix<T>const& Y, matrix<T>const& Z, std::string const& style);
    template<typename T>
    void quiver(matrix<T>const& vtx, matrix<T>const& dir, matrix<T>const& val);
    template<typename T>
    void tetmesh(matrix<std::size_t>const& tet, matrix<T>const& vtx, matrix<T>const& val);
    template<typename T>
    void trimesh(matrix<std::size_t>const& tri, matrix<T>const& vtx, matrix<T>const& val);
    template<typename T>
    void vermesh(matrix<std::size_t>const& ver, matrix<T>const& vtx, matrix<T>const& val);
    void writeimg(std::string const& filename) const;
    
    vtkSmartPointer<vtkContextView> GetView() {return m_view;}
    
private:
    vtkSmartPointer<vtkScalarBarActor> m_bar;
    vtkSmartPointer<vtkChartXY>        m_chartXY;
    vtkSmartPointer<vtkChartXYZ>       m_chartXYZ;
    matrix<double>                     m_range;
    vtkSmartPointer<vtkContextView>    m_view;
};

/// @name Constructor

//==========================================================================//
//                           MEMBER FUNCTIONS                               //
//==========================================================================//
//==========================================================================
// [figure.constructors]
/// Creates a figure object of size (\e w, \e h) pixels.
///
/// The code below builds a figure of size (\e w = W/3,\e h = H/3) where
/// W is the width of the screen and H is its height.
/// \code{.cpp}
///     figure fig;
/// \endcode
///
/// Build a figure with size (300,200)
/// \code{.cpp}
///     figure fig(300,200);
/// \endcode
figure::figure(int w, int h)
{
    m_view = vtkSmartPointer<vtkContextView>::New();
    int ws = m_view->GetRenderWindow()->GetScreenSize()[0];
    int hs = m_view->GetRenderWindow()->GetScreenSize()[1];
    if (w==0 && h==0) {w=ws/3; h=hs/3;}
    m_view->GetRenderWindow()->SetSize(w,h);
    m_view->GetRenderWindow()->SetPosition((ws-w)/2,(hs-h)/2);
    m_view->GetRenderer()->SetBackground(1,1,1);
}

//==========================================================================
// [figure.axes]
///
void figure::axes()
{
    if (m_view->GetInteractor()->GetInitialized()==0)
    {
        vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
        m_view->GetRenderer()->AddActor(axes);
        m_view->GetRenderer()->GradientBackgroundOn();
        m_view->GetRenderer()->SetBackground(1,1,1);
        m_view->GetRenderer()->SetBackground2(0.3,0.3,0.3);
    }
}

//==========================================================================
// [figure.caxis]
///
void figure::caxis(matrix<double>const& val)
{
    if (numel(val)!=2)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Range values must be a two element vector.");
    }
    if (numel(m_range)==0)
    {
        m_range = val;
    }
}

//==========================================================================
// [figure.colorbar]
///
void figure::colorbar(vtkMapper* mapper)
{
    // Setup mapper and lookup table
    mapper->SetScalarRange(m_range(0),m_range(1));
    vtkSmartPointer<vtkLookupTable> table = vtkSmartPointer<vtkLookupTable>::New();
    mapper->SetLookupTable(table);
    table->SetHueRange(0.65,0);
    table->Build();
    
    // Build colorbar at first call
    if (m_bar==nullptr)
    {
        vtkSmartPointer<vtkScalarBarActor> bar = vtkSmartPointer<vtkScalarBarActor>::New();
        bar->SetNumberOfLabels(5);
        bar->SetLookupTable(table);
        m_view->GetRenderer()->AddActor2D(bar);
    }
}

//==========================================================================
// [figure.interactor]
///
void figure::interactor()
{
    if (m_view->GetInteractor()->GetInitialized()==0)
    {
        vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        interactor->SetRenderWindow(m_view->GetRenderWindow());
        m_view->SetInteractor(interactor);
        m_view->GetInteractor()->ReInitialize();
        m_view->GetRenderWindow()->Render();
    }
    m_view->GetRenderer()->ResetCamera();
}

//==========================================================================
// [figure.drawnow]
///
void figure::drawnow()
{
    m_view->GetInteractor()->Start();
}

//==========================================================================
// [figure.edgmesh]
///
template<typename T>
void figure::edgmesh(matrix<std::size_t>const& edg, matrix<T>const& vtx, matrix<T>const& val)
{
    // Check input
    if (size(vtx,2)!=3 || size(edg,2)!=2)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Vertices should be a Nvtx-by-3 real matrix and elements a Nelt-by-2 integer matrix.");
    }
    if (numel(val)!=size(vtx,1) && numel(val)!=size(edg,1) && numel(val)!=0)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Values should be empty or have the same number of entry as vertices or elements.");
    }
    
    // Points list
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (std::size_t i=0; i<size(vtx,1); ++i)
    {
        points->InsertNextPoint(vtx(i,0),vtx(i,1),vtx(i,2));
    }
    
    // Edges list
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    for (std::size_t i=0; i<size(edg,1); ++i)
    {
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0,edg(i,0));
        line->GetPointIds()->SetId(1,edg(i,1));
        lines->InsertNextCell(line);
    }
    
    // Create polydata object
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetLines(lines);
    
    // Create mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);
    
    // Create mesh actor
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->GetProperty()->SetLineWidth(2);
    actor->GetProperty()->SetRenderLinesAsTubes(true);
    actor->SetMapper(mapper);
    m_view->GetRenderer()->AddActor(actor);
    
    // If data values
    if (numel(val)>0)
    {
        // Add data values
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfValues(numel(val));
        for (std::size_t l=0; l<numel(val); ++l)
        {
            data->SetValue(l,(double)val(l));
        }
        
        // Configure polydata and mapper
        if (numel(val)==size(vtx,1))
        {
            polydata->GetPointData()->SetScalars(data);
            mapper->SetScalarModeToUsePointData();
        }
        else if (numel(val)==size(edg,1))
        {
            polydata->GetCellData()->SetScalars(data);
            mapper->SetScalarModeToUseCellData();
        }
        mapper->ScalarVisibilityOn();
        mapper->SetColorModeToMapScalars();
        
        // Colorbar
        caxis({min(val),max(val)});
        colorbar(mapper);
    }
    
    // Configuration
    axes();
    interactor();
}

//==========================================================================
// [figure.imagesc]
///
template<typename T>
void figure::imagesc(matrix<T>const& M)
{
    // Create a rectangular grid for matrix representation
    vtkSmartPointer<vtkPlaneSource> grid = vtkSmartPointer<vtkPlaneSource>::New();
    grid->SetResolution((int)size(M,2),(int)size(M,1));
    grid->SetOrigin(0,0,0);
    grid->SetPoint1(1,0,0);
    double scale = size(M,1)/(double)size(M,2);
    grid->SetPoint2(0,scale,0);
    grid->Update();
    
    // Create Double array for cells data
    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
    data->SetNumberOfValues(numel(M));
    std::size_t l;
    for(std::size_t i=0; i<size(M,1); ++i)
    {
        for(std::size_t j=0; j<size(M,2); ++j)
        {
            l = (size(M,1)-1-i)*size(M,2)+j;
            data->SetValue(l,(double)M(i,j));
        }
    }
    
    // Create polydata using grid mesh and array data
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->DeepCopy(grid->GetOutput());
    polydata->GetCellData()->SetScalars(data);
    
    // Create mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);
    mapper->ScalarVisibilityOn();
    mapper->SetScalarModeToUseCellData();
    mapper->SetColorModeToMapScalars();
    
    // Create actor
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    m_view->GetRenderer()->AddActor(actor);
    
    // Configuration
    caxis({min(M),max(M)});
    colorbar(mapper);
    m_view->GetInteractor()->Initialize();
}

//==========================================================================
// [figure.mesh]
///
template<typename T>
void figure::mesh(matrix<T>const& X, matrix<T>const& Y, matrix<T>const& Z, std::string const& options)
{
    // Check input
    if (numel(X)!=numel(Y) || numel(X)!=numel(Z) || size(X,1)<2 || size(X,2)<2)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Input X, Y, and Z must have same size, with at least two elements per rows and columns.");
    }
    
    // Create a rectangular grid for matrix representation
    vtkSmartPointer<vtkPlaneSource> grid = vtkSmartPointer<vtkPlaneSource>::New();
    grid->SetResolution((int)size(X,2)-1,(int)size(X,1)-1);
    grid->Update();
    
    // Create polydata using grid mesh and array data
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->DeepCopy(grid->GetOutput());
    
    // Adjust points mesh to X, Y and Z values
    if (options.find("axis:normalize")!=std::string::npos)
    {
        T mxx = max(X), mnx = min(X);
        T mxy = max(Y), mny = min(Y);
        T mxz = max(Z), mnz = min(Z);
        for (std::size_t l=0; l<numel(X); ++l)
        {
            polydata->GetPoints()->SetPoint(l,-1+2*(X(l)-mnx)/(mxx-mnx),
                                            -1+2*(Y(l)-mny)/(mxy-mny),
                                            -1+2*(Z(l)-mnz)/(mxz-mnz));
        }
    }
    else
    {
        for (std::size_t l=0; l<numel(X); ++l)
        {
            polydata->GetPoints()->SetPoint(l,X(l),Y(l),Z(l));
        }
    }
    
    // Create double array for points data
    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
    data->SetNumberOfValues(numel(Z));
    std::size_t l;
    for(std::size_t i=0; i<size(Z,1); ++i)
    {
        for(std::size_t j=0; j<size(Z,2); ++j)
        {
            l = i*size(Z,2)+j;
            data->SetValue(l,(double)Z(i,j));
        }
    }
    polydata->GetPointData()->SetScalars(data);
    
    // Create mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);
    mapper->ScalarVisibilityOn();
    mapper->SetScalarModeToUsePointData();
    mapper->SetColorModeToMapScalars();
    
    // Create mesh actor
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetLineWidth(2.0);
    actor->GetProperty()->SetRenderLinesAsTubes(true);
    actor->SetMapper(mapper);
    m_view->GetRenderer()->AddActor(actor);
    
    // Configuration
    axes();
    caxis({min(Z),max(Z)});
    colorbar(mapper);
    interactor();
}

//==========================================================================
// [figure.plot]
///
template<typename T>
void figure::plot(matrix<T>const& X, matrix<T>const&Y, std::vector<std::string>const& style,
                  std::vector<std::string>const& label)
{
    // Check input
    if (!isvector(X) || size(Y,2)!=numel(X))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Input X must be a vector of size N and Y=f(X) a M-by-N matrix.");
    }
    
    // Configure chartXY
    if (m_view->GetScene()->GetNumberOfItems()==0)
    {
        m_chartXY = vtkSmartPointer<vtkChartXY>::New();
        m_chartXY->SetShowLegend(true);
        m_chartXY->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(12);
        m_chartXY->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(15);
        m_chartXY->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(12);
        m_chartXY->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(15);
        m_chartXY->GetLegend()->SetLabelSize(15);
        m_view->GetScene()->AddItem(m_chartXY);
    }
    
    // Add multiple plots with stylesheet
    std::size_t m=size(Y,1), n=size(Y,2);
    vtkSmartPointer<vtkTable> table;
    vtkSmartPointer<vtkDoubleArray> arrX, arrY;
    vtkPlot *plot;
    for (std::size_t i=0; i<m; ++i)
    {
        // Configure plot style
        if (i<style.size())
        {
            // Define plot (points are default)
            if (style[i].find("-")!=std::string::npos)
            {
                plot = m_chartXY->AddPlot(vtkChart::LINE);
            }
            else {plot = m_chartXY->AddPlot(vtkChart::POINTS);}
            
            // Set marker
            if (style[i].find("x")!=std::string::npos)
            {
                vtkPlotPoints::SafeDownCast(plot)->SetMarkerStyle(vtkPlotPoints::CROSS);
            }
            else if (style[i].find("o")!=std::string::npos)
            {
                vtkPlotPoints::SafeDownCast(plot)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
            }
            else if (style[i].find("+")!=std::string::npos)
            {
                vtkPlotPoints::SafeDownCast(plot)->SetMarkerStyle(vtkPlotPoints::PLUS);
            }
            else if (style[i].find("d")!=std::string::npos)
            {
                vtkPlotPoints::SafeDownCast(plot)->SetMarkerStyle(vtkPlotPoints::DIAMOND);
            }
            else if (style[i].find("s")!=std::string::npos)
            {
                vtkPlotPoints::SafeDownCast(plot)->SetMarkerStyle(vtkPlotPoints::SQUARE);
            }
            
            // Set color
            if (style[i].find("r")!=std::string::npos) {plot->SetColor(255,0,0,255);}
            else if (style[i].find("g")!=std::string::npos) {plot->SetColor(0,255,0,255);}
            else if (style[i].find("b")!=std::string::npos) {plot->SetColor(0,0,255,255);}
            else if (style[i].find("y")!=std::string::npos) {plot->SetColor(255,255,0,255);}
            else if (style[i].find("c")!=std::string::npos) {plot->SetColor(0,255,255,255);}
            else if (style[i].find("m")!=std::string::npos) {plot->SetColor(255,0,255,255);}
            else if (style[i].find("k")!=std::string::npos) {plot->SetColor(0,0,0,255);}
            else if (style[i].find("w")!=std::string::npos) {plot->SetColor(255,255,255,255);}
        }
        else
        {
            plot = m_chartXY->AddPlot(vtkChart::POINTS);
        }
        
        // Convert input data to table
        arrX = vtkSmartPointer<vtkDoubleArray>::New();
        arrX->SetName("X");
        arrY = vtkSmartPointer<vtkDoubleArray>::New();
        if (i<label.size()) {arrY->SetName(label[i].c_str());}
        else {arrY->SetName("");}
        table = vtkSmartPointer<vtkTable>::New();
        table->AddColumn(arrX);
        table->AddColumn(arrY);
        table->SetNumberOfRows(n);
        for (std::size_t j=0; j<n; ++j)
        {
            table->SetValue(j, 0, X(j));
            table->SetValue(j, 1, Y(i,j));
        }
        
        // Add plot data
        plot->SetInputData(table,0,1);
    }
    
    // Configuration
    m_view->GetInteractor()->Initialize();
}

//==========================================================================
// [figure.plot3]
///
template<typename T>
void figure::plot3(matrix<T>const& X, matrix<T>const& Y, matrix<T>const& Z, std::string const& style)
{
    // Check input
    if (numel(Y)!=numel(X) || numel(Z)!=numel(X))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Input X, Y and Z must be same size.");
    }
    
    // Configure chartXYZ
    if (m_view->GetScene()->GetNumberOfItems()==0)
    {
        m_chartXYZ = vtkSmartPointer<vtkChartXYZ>::New();
        int w = m_view->GetRenderWindow()->GetSize()[0];
        int h = m_view->GetRenderWindow()->GetSize()[1];
        int b = 10;
        m_chartXYZ->SetGeometry(vtkRectf(b, b, w-2*b, h-2*b));
        m_view->GetScene()->AddItem(m_chartXYZ);
    }
    
    // Configure plot style
    vtkSmartPointer<vtkPlot3D> plot;
    if (style.find("-")!=std::string::npos)
    {
        plot = vtkSmartPointer<vtkPlotLine3D>::New();
    }
    else
    {
        plot = vtkSmartPointer<vtkPlotPoints3D>::New();
    }
    if (style.find("r")!=std::string::npos) {plot->GetPen()->SetColor(255,0,0,255);}
    else if (style.find("g")!=std::string::npos) {plot->GetPen()->SetColor(0,255,0,255);}
    else if (style.find("b")!=std::string::npos) {plot->GetPen()->SetColor(0,0,255,255);}
    else if (style.find("y")!=std::string::npos) {plot->GetPen()->SetColor(255,255,0,255);}
    else if (style.find("c")!=std::string::npos) {plot->GetPen()->SetColor(0,255,255,255);}
    else if (style.find("m")!=std::string::npos) {plot->GetPen()->SetColor(255,0,255,255);}
    else if (style.find("k")!=std::string::npos) {plot->GetPen()->SetColor(0,0,0,255);}
    else if (style.find("w")!=std::string::npos) {plot->GetPen()->SetColor(255,255,255,255);}
    
    // Convert input data to table
    vtkNew<vtkTable> table;
    vtkNew<vtkDoubleArray> arrX, arrY, arrZ;
    arrX->SetName("X");
    table->AddColumn(arrX);
    arrY->SetName("Y");
    table->AddColumn(arrY);
    arrZ->SetName("Z");
    table->AddColumn(arrZ);
    std::size_t n=numel(X);
    table->SetNumberOfRows(n);
    for (std::size_t l=0; l<n; ++l)
    {
        table->SetValue(l, 0, X(l));
        table->SetValue(l, 1, Y(l));
        table->SetValue(l, 2, Z(l));
    }
    
    // Add plot data and add to chart
    plot->SetInputData(table);
    m_chartXYZ->AddPlot(plot);
    
    // Configuration
    m_view->GetInteractor()->Initialize();
}

//==========================================================================
// [figure.quiver]
///
template<typename T>
void figure::quiver(matrix<T>const& vtx, matrix<T>const& dir, matrix<T>const& val)
{
    // Check input
    if (size(vtx,1)!=size(dir,1) || size(vtx,2)!=3 || size(dir,2)!=3)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Vertices and directions should be N-by-3 real matrix.");
    }
    if (numel(val)!=size(vtx,1) && numel(val)!=0)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Values should be empty or have the same number of entry as vertices or elements.");
    }
    
    // Distance
    matrix<double> dst = sqrt(sum(pow(dir,2),2));
    matrix<double> ctr = vtx + dir/2;
    
    // Initialize geometric object
    vtkSmartPointer<vtkConeSource> coneSource = vtkSmartPointer<vtkConeSource>::New();
    vtkSmartPointer<vtkPoints>     points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray>  cells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData>   cone;
    vtkSmartPointer<vtkTriangle>   triangle;

    // For each direction
    std::size_t n = 0;
    for (std::size_t i=0; i<size(vtx,1); ++i)
    {
        // Define current cone
        coneSource->SetDirection(dir(i,0),dir(i,1),dir(i,2));
        coneSource->SetHeight(dst(i));
        coneSource->SetRadius(dst(i)/20);
        coneSource->SetCenter(ctr(i,0),ctr(i,1),ctr(i,2));
        coneSource->Update();
        cone = coneSource->GetOutput();
        
        // Points list
        for (std::size_t l=0; l<cone->GetNumberOfPoints(); ++l)
        {
            points->InsertNextPoint(cone->GetPoint(l));
        }
        
        // Cell list
        for (std::size_t l=0; l<cone->GetNumberOfCells(); ++l)
        {
            triangle = vtkSmartPointer<vtkTriangle>::New();
            for(std::size_t j=0; j<3; ++j)
            {
                triangle->GetPointIds()->SetId(j,n+cone->GetCell((int)l)->GetPointId((int)j));
            }
            cells->InsertNextCell(triangle);
        }
        
        // Increment element table
        n = n + cone->GetNumberOfPoints();
    }
    
    // Define polydata
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetPolys(cells);
    
    // Add mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);
    
    // Create mesh actor
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    m_view->GetRenderer()->AddActor(actor);
    
    // Configuration
    axes();
    interactor();
}

//==========================================================================
// [figure.tetmesh]
///
template<typename T>
void figure::tetmesh(matrix<std::size_t>const& tet, matrix<T>const& vtx, matrix<T>const& val)
{
    // Check input
    if (size(vtx,2)!=3 || size(tet,2)!=4)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Vertices should be a Nvtx-by-3 real matrix and elements a Nelt-by-4 integer matrix.");
    }
    if (numel(val)!=size(vtx,1) && numel(val)!=size(tet,1) && numel(val)!=0)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Values should be empty or have the same number of entry as vertices or elements.");
    }
    
    // Points list
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (std::size_t i=0; i<size(vtx,1); ++i)
    {
        points->InsertNextPoint(vtx(i,0),vtx(i,1),vtx(i,2));
    }
    
    // Tetras list
    vtkSmartPointer<vtkCellArray> tetras = vtkSmartPointer<vtkCellArray>::New();
    for (std::size_t i=0; i<size(tet,1); ++i)
    {
        vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
        tetra->GetPointIds()->SetId(0,tet(i,0));
        tetra->GetPointIds()->SetId(1,tet(i,1));
        tetra->GetPointIds()->SetId(2,tet(i,2));
        tetra->GetPointIds()->SetId(3,tet(i,3));
        tetras->InsertNextCell(tetra);
    }
    
    // Create unstructured grid object
    vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ungrid->SetPoints(points);
    ungrid->SetCells(VTK_TETRA,tetras);
    
    // Create mapper
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(ungrid);
    
    // Create mesh actor
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetLineWidth(2.0);
    actor->GetProperty()->SetRenderLinesAsTubes(true);
    actor->SetMapper(mapper);
    m_view->GetRenderer()->AddActor(actor);
    
    // If data values
    if (numel(val)>0)
    {
        // Add data values
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfValues(numel(val));
        for (std::size_t l=0; l<numel(val); ++l)
        {
            data->SetValue(l,(double)val(l));
        }
        
        // Configure grid and mapper
        if (numel(val)==size(vtx,1))
        {
            ungrid->GetPointData()->SetScalars(data);
            mapper->SetScalarModeToUsePointData();
        }
        else if (numel(val)==size(tet,1))
        {
            ungrid->GetCellData()->SetScalars(data);
            mapper->SetScalarModeToUseCellData();
        }
        mapper->ScalarVisibilityOn();
        mapper->SetColorModeToMapScalars();
        mapper->SetScalarRange(min(val),max(val));
        
        // Colorbar
        caxis({min(val),max(val)});
        colorbar(mapper);
    }
    
    // Configuration
    axes();
    interactor();
}

//==========================================================================
// [figure.trimesh]
///
template<typename T>
void figure::trimesh(matrix<std::size_t>const& tri, matrix<T>const& vtx, matrix<T>const& val)
{
    // Check input
    if (size(vtx,2)!=3 || size(tri,2)!=3)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Vertices should be a Nvtx-by-3 real matrix and elements a Nelt-by-3 integer matrix.");
    }
    if (numel(val)!=size(vtx,1) && numel(val)!=size(tri,1) && numel(val)!=0)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Values should be empty or have the same number of entry as vertices or elements.");
    }
    
    // Points list
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (std::size_t i=0; i<size(vtx,1); ++i)
    {
        points->InsertNextPoint(vtx(i,0),vtx(i,1),vtx(i,2));
    }
    
    // Triangles list
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    for (std::size_t i=0; i<size(tri,1); ++i)
    {
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        triangle->GetPointIds()->SetId(0,tri(i,0));
        triangle->GetPointIds()->SetId(1,tri(i,1));
        triangle->GetPointIds()->SetId(2,tri(i,2));
        triangles->InsertNextCell(triangle);
    }
    
    // Create polydata object
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetPolys(triangles);
    
    // Create mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);
    
    // Create mesh actor
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetLineWidth(2.0);
    actor->GetProperty()->SetRenderLinesAsTubes(true);
    actor->SetMapper(mapper);
    m_view->GetRenderer()->AddActor(actor);
    
    // If data values
    if (numel(val)>0)
    {
        // Add data values
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfValues(numel(val));
        for (std::size_t l=0; l<numel(val); ++l)
        {
            data->SetValue(l,(double)val(l));
        }
        
        // Configure polydata and mapper
        if (numel(val)==size(vtx,1))
        {
            polydata->GetPointData()->SetScalars(data);
            mapper->SetScalarModeToUsePointData();
        }
        else if (numel(val)==size(tri,1))
        {
            polydata->GetCellData()->SetScalars(data);
            mapper->SetScalarModeToUseCellData();
        }
        mapper->ScalarVisibilityOn();
        mapper->SetColorModeToMapScalars();
        
        // Colorbar
        caxis({min(val),max(val)});
        colorbar(mapper);
    }
    
    // Configuration
    axes();
    interactor();
}

//==========================================================================
// [figure.vermesh]
///
template<typename T>
void figure::vermesh(matrix<std::size_t>const& ver, matrix<T>const& vtx, matrix<T>const& val)
{
    // Check input
    if (size(vtx,2)!=3 || !isvector(ver))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Vertices should be a Nvtx-by-3 real matrix and elements a Nelt integer vector.");
    }
    if (numel(val)!=size(vtx,1) && numel(val)!=numel(ver) && numel(val)!=0)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Values should be empty or have the same number of entry as vertices or elements.");
    }
    
    // Points list
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (std::size_t i=0; i<size(vtx,1); ++i)
    {
        points->InsertNextPoint(vtx(i,0),vtx(i,1),vtx(i,2));
    }
    
    // Vertices list
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
    for (std::size_t l=0; l<numel(ver); ++l)
    {
        vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
        vertex->GetPointIds()->SetId(0,ver(l));
        vertices->InsertNextCell(vertex);
    }
    
    // Create polydata object
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetVerts(vertices);
    
    // Create mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);
    
    // Create mesh actor
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->GetProperty()->VertexVisibilityOn();
    actor->GetProperty()->SetPointSize(8);
    actor->GetProperty()->SetRenderPointsAsSpheres(true);
    actor->SetMapper(mapper);
    m_view->GetRenderer()->AddActor(actor);
    
    // If data values
    if (numel(val)>0)
    {
        // Add data values
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfValues(numel(val));
        for (std::size_t l=0; l<numel(val); ++l)
        {
            data->SetValue(l,(double)val(l));
        }
        
        // Configure polydata and mapper
        if (numel(val)==size(vtx,1))
        {
            polydata->GetPointData()->SetScalars(data);
            mapper->SetScalarModeToUsePointData();
        }
        else if (numel(val)==numel(ver))
        {
            polydata->GetCellData()->SetScalars(data);
            mapper->SetScalarModeToUseCellData();
        }
        mapper->ScalarVisibilityOn();
        mapper->SetColorModeToMapScalars();
        
        // Colorbar
        caxis({min(val),max(val)});
        colorbar(mapper);
    }
    
    // Configuration
    axes();
    interactor();
}

//==========================================================================
// [figure.writeimg]
///
void figure::writeimg(std::string const& filename) const
{
    if (!filename.empty())
    {
        // Check filename and extension (.png is default)
        std::string fn = filename;
        std::string ext;
        auto found = fn.find_last_of(".");
        if (found == std::string::npos)
        {
            ext = ".png";
            fn += ext;
        }
        else
        {
            ext = filename.substr(found,filename.size());
        }
        std::locale loc;
        std::transform(ext.begin(), ext.end(), ext.begin(),
                       [=](char const& c) { return std::tolower(c, loc); });

        // Select file write from extension
        auto writer = vtkSmartPointer<vtkImageWriter>::New();
        if (ext == ".bmp")
        {
            writer = vtkSmartPointer<vtkBMPWriter>::New();
        }
        else if (ext == ".jpg")
        {
            writer = vtkSmartPointer<vtkJPEGWriter>::New();
        }
        else if (ext == ".pnm")
        {
            writer = vtkSmartPointer<vtkPNMWriter>::New();
        }
        else if (ext == ".ps")
        {
            writer = vtkSmartPointer<vtkPostScriptWriter>::New();
        }
        else if (ext == ".tiff")
        {
            writer = vtkSmartPointer<vtkTIFFWriter>::New();
        }
        else
        {
            writer = vtkSmartPointer<vtkPNGWriter>::New();
        }

        // Image buffer
        vtkNew<vtkWindowToImageFilter> window_to_image_filter;
        window_to_image_filter->SetInput(m_view->GetRenderWindow());
        window_to_image_filter->SetScale(1); // image quality
        window_to_image_filter->SetInputBufferTypeToRGB();
        window_to_image_filter->ReadFrontBufferOff();
        window_to_image_filter->Update();
        
        // Write
        writer->SetFileName(fn.c_str());
        writer->SetInputConnection(window_to_image_filter->GetOutputPort());
        writer->Write();
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"No filename provided.");
    }
}

//==========================================================================
// [figure.xlim]
///
void figure::xlim(matrix<double>const& limit)
{
    if (!isvector(limit) || numel(limit) != 2)
    {
        error(__FILE__, __LINE__, __FUNCTION__, "Limits must be a 2-element vector of increasing numeric values.");
    }
    m_chartXY->GetAxis(1)->SetBehavior(vtkAxis::FIXED);
    m_chartXY->GetAxis(1)->SetRange(limit(0),limit(1));
}

//==========================================================================
// [figure.ylim]
///
void figure::ylim(matrix<double>const& limit)
{
    if (!isvector(limit) || numel(limit) != 2)
    {
        error(__FILE__, __LINE__, __FUNCTION__, "Limits must be a 2-element vector of increasing numeric values.");
    }
    m_chartXY->GetAxis(0)->SetBehavior(vtkAxis::FIXED);
    m_chartXY->GetAxis(0)->SetRange(limit(0),limit(1));
}

//==========================================================================//
//                          MATLAB-LIKE FUNCTIONS                           //
//==========================================================================//
//==========================================================================
// [caxis]
/// Sets the colormap limit for the current axis.
///
/// \e fig is the figure object for which those limits must be set and \e val is
/// a two-element matrix containing the lower- and upper-bound.
///
/// \code{.cpp}
///     matrix<> X,Y;
///     std::tie(X,Y) = meshgrid(linspace(-M_PI,M_PI,100));
///     auto Z = 2*sin(X)/X * sin(Y)/Y;
///
///     figure fig;
///     caxis(fig,{-1,1});
///     mesh(fig,X,Y,Z);
///     drawnow(fig);
/// \endcode
inline void caxis(figure& fig, matrix<double>const& val)
{
    fig.caxis(val);
}

//==========================================================================
// [drawnow]
/// Displays all the figure defined since the beginning of the program or
/// the last call to \e drawnow.
///
/// Calling \e drawnow is \b blocking meaning that it will suspend the
/// execution of the program. Plotting cannot be deported to another thread
/// as it is not permitted by the VTK framework. The execution will resume
/// all the currently displayed figures are closed.
inline void drawnow(figure& fig)
{
    fig.drawnow();
}

//==========================================================================
// [edgmesh]
/// Displays a set of edges colored using vertex values.
///
/// The edges must be given in a mesh-like format where \e vtx is the list
/// of all the nodes which must be displayed and \e edg is the list of the
/// elements. The color on the edges is obtained by setting \e val which
/// must contain a list of vertex-values.
///
// \see vermesh, trimesh, tetmesh.
template<typename T>
inline void edgmesh(figure& fig, matrix<std::size_t>const& edg, matrix<T>const& vtx,
                    matrix<T>const& val={})
{
    fig.edgmesh(edg,vtx,val);
}

//==========================================================================
// [imagesc]
/// Displays the data in a matrix using scaled colors.
template<typename T>
inline void imagesc(figure& fig, matrix<T>const& M)
{
    fig.imagesc(M);
}

//==========================================================================
// [mesh]
/// Displays a surface Z = f(X,Y) using scaled colors.
///
/// \code{.cpp}
///     matrix<> X,Y;
///     std::tie(X,Y) = meshgrid(linspace(-M_PI,M_PI,100));
///     auto Z = 2*sin(X)/X * sin(Y)/Y;
///     figure fig;
///     mesh(fig,X,Y,Z);
///     drawnow(fig);
/// \endcode
template<typename T>
inline void mesh(figure& fig, matrix<T>const& X, matrix<T>const& Y, matrix<T>const& Z, std::string const& options="")
{
    fig.mesh(X,Y,Z,options);
}
template<typename T>
void mesh(figure& fig, matrix<T>const& M, std::string options="")
{
    matrix<T> X, Y;
    X = linspace<T>(1,size(M,2),size(M,2));
    Y = linspace<T>(1,size(M,1),size(M,1));
    std::tie(X,Y) = meshgrid(X,Y);
    options += " axis:normalize ";
    fig.mesh(X,Y,M,options);
}

//==========================================================================
// [plot]
/// Plots a 2D-curve Y = f(X) with customizable options.
///
/// Multiple plots can be given at once by giving Y as a multi-dimensional
/// array. In that case, each line corresponds to a different plot. The
/// customization options must be given as an array of strings (one string
/// per plot). Each string \e may (but it is not mandatory) contain one of the
/// following elements:\n
///  - a line-style. If "-" is specified, the style is \e solid; otherwise
/// the line will be \e dotted.\n
///  - a marker-style which can be any of the following : "x" (cross), "o"
/// (circle), "+" (plus sign), "d" (diamond), "s" (square).\n
///  - a color which can be any of the following : "r" (red), "g" (green),
/// "b" (blue), "y" (yellow), "c" (cyan), "m" (magenta), "k" (black) or
/// "w" (white).
///
/// The same format should be used for the labels.
///
/// \code{.cpp}
///     matrix<> X = linspace(0,10,100);
///     matrix<> Y = cos(X);
///
///     // plot y=sin(x) as a red solid line with cross markers
///     figure fig;
///     plot(fig,X,Y,{"r-x"});
///     drawnow(fig);
/// \endcode
///
// \see plot3, spy.
template<typename T=double>
inline void plot(figure& fig, matrix<T>const& X, matrix<T>const&Y, std::vector<std::string>const& style={""},
                 std::vector<std::string>const& label={""})
{
    fig.plot(X,Y,style,label);
}

//==========================================================================
// [plot3]
/// Plots curves defined by (X,Y,Z) in 3D-space, eventually linked by a
/// styled line.
///
/// On the contrary to \e plot, only one curve can be given at a time.
///
// \see plot
template<typename T=double>
inline void plot3(figure& fig, matrix<T>const& X, matrix<T>const& Y, matrix<T>const& Z, std::string const& style="")
{
    fig.plot3(X,Y,Z,style);
}


//==========================================================================
// [quiver]
/// Plots a set of vectors.
///
/// The vectors are defined by their origin \e vtx and their direction \e dir.
///
/// \code{.cpp}
///     matrix<> vtx({{0.,0.,0.},{0.,1.,0.}});
///     matrix<> dir({{1.,1.,1.},{-0.2,1.,3.}});
///
///     figure fig;
///     quiver(fig,vtx,dir);
///     drawnow(fig);
/// \endcode
///
// \see plot, plot3.
template<typename T>
inline void quiver(figure& fig, matrix<T>const& vtx, matrix<T>const& dir, matrix<T>const& val={})
{
    fig.quiver(vtx,dir,val);
}

//==========================================================================
// [spy]
/// Spy sparse matrix non-zeros values.
///
/// \code{.cpp}
///     smatrix<> Ms = speye(10);
///     figure fig;
///     spy(Ms,"b","Ms");
///     drawnow(fig);
/// \endcode
///
// \see plot.
template<typename T>
inline void spy(figure& fig, smatrix<T>const& Ms, std::string const& style="", std::string const& label="")
{
    matrix<long> X, Y;
    long m=size(Ms,1)-1, n=size(Ms,2)-1;
    X = {0,n,n,0,0};
    Y = {0,0,-m,-m,0};
    fig.plot(X,Y,{"-k"},{""});
    matrix<T> V;
    std::tie(Y,X,V) = find(Ms);
    fig.plot(X,-Y,{style},{label+" (nnz="+std::to_string(nnz(Ms))+")"});
}

//==========================================================================
// [tetboundary]
/// Extracts the boundary of a tetrahedral mesh defined by its elements
/// \e elt and its vertices \e vtx and returns the result as a mesh in the
/// same format.
///
/// \code{.cpp}
///     matrix<> X({0.,1,0,0}), Y({0.,0,1,0}), Z({0,0,0,1});
///     matrix<> vtx, bvtx;
///     matrix<> elt, belt;
///     // Delaunay triangulation
///     std::tie(elt,vtx) = tetdelaunay(X,Y,Z);
///     // Boundary extraction
///     std::tie(belt,bvtx) = tetboundary(elt,vtx);
/// \endcode
///
// \see tetdelaunay, tetmesh, trimesh.
template<typename T>
auto tetboundary(matrix<std::size_t>const& tet, matrix<T>const& vtx)
{
    // Check input
    if (size(vtx,2)!=3 || size(tet,2)!=4)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Vertices should be a Nvtx-by-3 real matrix and elements a Nelt-by-4 integer matrix.");
    }
    
    // Points list
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (std::size_t i=0; i<size(vtx,1); ++i)
    {
        points->InsertNextPoint(vtx(i,0),vtx(i,1),vtx(i,2));
    }
    
    // Tetras list
    vtkSmartPointer<vtkCellArray> tetras = vtkSmartPointer<vtkCellArray>::New();
    for (std::size_t i=0; i<size(tet,1); ++i)
    {
        vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
        tetra->GetPointIds()->SetId(0,tet(i,0));
        tetra->GetPointIds()->SetId(1,tet(i,1));
        tetra->GetPointIds()->SetId(2,tet(i,2));
        tetra->GetPointIds()->SetId(3,tet(i,3));
        tetras->InsertNextCell(tetra);
    }
    
    // Create unstructured grid object
    vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ungrid->SetPoints(points);
    ungrid->SetCells(VTK_TETRA,tetras);
    
    // Surface filtering
    vtkSmartPointer<vtkDataSetSurfaceFilter> filter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    filter->SetInputData(ungrid);
    filter->Update();
    vtkSmartPointer<vtkPolyData> polydata = filter->GetOutput();
    
    // Clean polydata mesh
    vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanPolyData->SetInputData(polydata);
    cleanPolyData->Update();
    polydata = cleanPolyData->GetOutput();
    
    // Extract vertices from polydata
    matrix<T> vtx2(polydata->GetNumberOfPoints(),3);
    for(std::size_t i=0; i<size(vtx2,1); ++i)
    {
        for(std::size_t j=0; j<size(vtx2,2); ++j)
        {
            vtx2(i,j) = polydata->GetPoint(int(i))[j];
        }
    }
    
    // Extract triangles
    matrix<std::size_t> tri(polydata->GetNumberOfCells(),3);
    for(std::size_t i=0; i<size(tri,1); ++i)
    {
        for(std::size_t j=0; j<size(tri,2); ++j)
        {
            tri(i,j) = polydata->GetCell((int)i)->GetPointId((int)j);
        }
    }
    
    // Return surfacic mesh
    return std::make_tuple(tri,vtx2);
}

//==========================================================================
// [tetdelaunay]
/// Creates a Delaunay tetrahedral mesh from a set of nodes (X,Y,Z).
///
/// \code{.cpp}
///     matrix<> X({0.,1,0,0}), Y({0.,0,1,0}), Z({0,0,0,1});
///     matrix<> vtx;
///     matrix<std::size_t> elt;
///     // Delaunay triangulation
///     std::tie(elt,vtx) = tetdelaunay(X,Y,Z);
/// \endcode
///
// \see tetmesh, tetboundary.
template<typename T>
auto tetdelaunay(matrix<T>const& X, matrix<T>const& Y, matrix<T>const& Z)
{
    // Test input compatibility
    if (numel(X)!=numel(Y) || numel(X)!=numel(Z))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Input X, Y and Z must be same size.");
    }
    
    // Define points using input grid in linear indexing
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(std::size_t l=0; l<numel(X); ++l)
    {
        points->InsertNextPoint(X(l),Y(l),Z(l));
    }
    
    // Delaunay meshing using polydata input
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    
    //    // Clean polydata mesh
    //    vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
    //    cleanPolyData->SetInputData(polydata);
    //    cleanPolyData->Update();
    //    polydata = cleanPolyData->GetOutput();
    
    // Tetra meshing
    vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
    delaunay->SetInputData(polydata);
    delaunay->Update();
    vtkSmartPointer<vtkUnstructuredGrid> ungrid = delaunay->GetOutput();
    
    // Extract vertices
    matrix<T> vtx(ungrid->GetNumberOfPoints(),3);
    for(std::size_t i=0; i<size(vtx,1); ++i)
    {
        for(std::size_t j=0; j<size(vtx,2); ++j)
        {
            vtx(i,j) = ungrid->GetPoint(int(i))[j];
        }
    }
    
    // Extract tetras
    matrix<std::size_t> tet(ungrid->GetNumberOfCells(),4);
    for(std::size_t i=0; i<size(tet,1); ++i)
    {
        for(std::size_t j=0; j<size(tet,2); ++j)
        {
            tet(i,j) = ungrid->GetCell((int)i)->GetPointId((int)j);
        }
    }
    
    // Return tables
    return std::make_tuple(tet,vtx);
}

//==========================================================================
// [tetmesh]
/// Plots a tetrahedral mesh colored with vertex values.
///
/// \code{.cpp}
///     matrix<> X({0.,1,0,0}), Y({0.,0,1,0}), Z({0,0,0,1});
///     matrix<> vtx;
///     matrix<std::size_t> elt;
///     // Delaunay triangulation
///     std::tie(elt,vtx) = tetdelaunay(X,Y,Z);
///
///     figure fig;
///     tetmesh(fig,elt,vtx);
///     drawnow(fig);
/// \endcode
///
// \see tetdelaunay, vermesh, edgmesh, trimesh.
template<typename T>
inline void tetmesh(figure& fig, matrix<std::size_t>const& tet, matrix<T>const& vtx, matrix<T>const& val={})
{
    fig.tetmesh(tet,vtx,val);
}

//==========================================================================
// [tridelaunay]
/// Computes a Delaunay triangulation from a set of nodes.
///
// \see trimesh, tetdelaunay.
template<typename T>
auto tridelaunay(matrix<T>const& X, matrix<T>const& Y, matrix<T>const& Z={})
{
    // Test input compatibility
    if (numel(X)!=numel(Y) || (numel(Z)>0 && numel(X)!=numel(Z)))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Input coordinates X, Y and Z must have same number of elements.");
    }
    
    // Define points using input grid in linear indexing
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(std::size_t l=0; l<numel(X); ++l)
    {
        if (numel(Z)>0) {points->InsertNextPoint(X(l),Y(l),Z(l));}
        else {points->InsertNextPoint(X(l),Y(l),0);}
    }
    
    // Delaunay meshing using polydata input/output
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
    delaunay->SetInputData(polydata);
    delaunay->Update();
    polydata = delaunay->GetOutput();
    
    // Clean polydata mesh
    vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanPolyData->SetInputData(polydata);
    cleanPolyData->Update();
    polydata = cleanPolyData->GetOutput();
    
    // Extract vertices
    matrix<T> vtx(polydata->GetNumberOfPoints(),3);
    for(std::size_t i=0; i<size(vtx,1); ++i)
    {
        for(std::size_t j=0; j<size(vtx,2); ++j)
        {
            vtx(i,j) = polydata->GetPoint(int(i))[j];
        }
    }
    
    // Extract triangles
    matrix<std::size_t> tri(polydata->GetNumberOfCells(),3);
    for(std::size_t i=0; i<size(tri,1); ++i)
    {
        for(std::size_t j=0; j<size(tri,2); ++j)
        {
            tri(i,j) = polydata->GetCell((int)i)->GetPointId((int)j);
        }
    }
    
    // Return tables
    return std::make_tuple(tri,vtx);
}

//==========================================================================
// [trimesh]
/// Plots a triangular mesh colored using vertex values.
///
// \see tridelaunay, vermesh, edgmesh, tetmesh.
template<typename T>
inline void trimesh(figure& fig, matrix<std::size_t>const& tri, matrix<T>const& vtx, matrix<T>const& val={})
{
    fig.trimesh(tri,vtx,val);
}

//==========================================================================
// [triread]
/// Reads a triangular mesh file in \e .ply or \e .vtk format.
///
// \see triwrite, trimesh.
template<typename T=double>
auto triread(std::string const& path, std::string const& name)
{
    // Read file
    vtkSmartPointer<vtkPolyData> polydata;
    std::string filename = path+name;
    if (name.find(".ply")!=std::string::npos)
    {
        vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();
        polydata = reader->GetOutput();
    }
    else if (name.find(".vtk")!=std::string::npos)
    {
        vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();
        polydata = reader->GetPolyDataOutput();
    }
    else
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Unsupported format, please use .ply or .vtk.");
    }
    
    // Clean polydata mesh
    vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanPolyData->SetInputData(polydata);
    cleanPolyData->Update();
    polydata = cleanPolyData->GetOutput();
    
    // Extract vertices
    matrix<T> vtx(polydata->GetNumberOfPoints(),3);
    for(std::size_t i=0; i<size(vtx,1); ++i)
    {
        for(std::size_t j=0; j<size(vtx,2); ++j)
        {
            vtx(i,j) = polydata->GetPoint(int(i))[j];
        }
    }
    
    // Extract elements
    matrix<std::size_t> tri(polydata->GetNumberOfCells(),3);
    for(std::size_t i=0; i<size(tri,1); ++i)
    {
        for(std::size_t j=0; j<size(tri,2); ++j)
        {
            tri(i,j) = polydata->GetCell((int)i)->GetPointId((int)j);
        }
    }
    
    // Return tables
    return std::make_tuple(tri,vtx);
}

//==========================================================================
// [triwrite]
/// Writes a triangular mesh file in \e .ply or \e .vtk format.
///
// \see triread, trimesh.
template<typename T>
void triwrite(std::string const& path, std::string const& name,
              matrix<std::size_t>const& tri, matrix<T>const& vtx, matrix<T>const& val={})
{
    // Check input
    if (size(vtx,2)!=3 || size(tri,2)!=3)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Vertices should be a Nvtx-by-3 real matrix and elements a Nelt-by-3 integer matrix.");
    }
    if (numel(val)!=size(vtx,1) && numel(val)!=size(tri,1) && numel(val)!=0)
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Values should be empty or have the same number of entry as vertices or elements.");
    }
    
    // Points list
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (std::size_t i=0; i<size(vtx,1); ++i)
    {
        points->InsertNextPoint(vtx(i,0),vtx(i,1),vtx(i,2));
    }
    
    // Triangles list
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    for (std::size_t i=0; i<size(tri,1); ++i)
    {
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        triangle->GetPointIds()->SetId(0,tri(i,0));
        triangle->GetPointIds()->SetId(1,tri(i,1));
        triangle->GetPointIds()->SetId(2,tri(i,2));
        triangles->InsertNextCell(triangle);
    }
    
    // Create polydata object
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetPolys(triangles);
    
    // If data values
    if (numel(val)>0)
    {
        // Add data values
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfValues(numel(val));
        for (std::size_t l=0; l<numel(val); ++l)
        {
            data->SetValue(l,(double)val(l));
        }
        
        // Configure polydata and mapper
        if (numel(val)==size(vtx,1))
        {
            polydata->GetPointData()->SetScalars(data);
        }
        else if (numel(val)==size(tri,1))
        {
            polydata->GetCellData()->SetScalars(data);
        }
    }
    
    // Write file
    std::string filename = path+name;
    if (name.find(".ply")!=std::string::npos)
    {
        vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetInputData(polydata);
        writer->Write();
    }
    else if (name.find(".vtk")!=std::string::npos)
    {
        vtkSmartPointer<vtkGenericDataObjectWriter> writer = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetInputData(polydata);
        writer->Write();
    }
    else
    {
        warning(__FILE__, __LINE__, __FUNCTION__,"Unsupported format, please use .ply or .vtk.");
    }
}

//==========================================================================
// [vermesh]
/// Plots colored vertices using vertex values.
///
/// \code{.cpp}
///     matrix<std::size_t> elt = range(0,100);
///     matrix<> vtx = -1+2*rand(numel(elt),3);
///     figure fig;
///     vermesh(fig,elt,vtx,eval(vtx(row(vtx),0)));
///     drawnow(fig);
/// \endcode
///
// \see edgmesh, trimesh, tetmesh.
template<typename T>
inline void vermesh(figure& fig, matrix<std::size_t>const& ver, matrix<T>const& vtx, matrix<T>const& val={})
{
    fig.vermesh(ver,vtx,val);
}

template<typename T>
inline void vermesh(figure& fig, matrix<T>const& vtx)
{
    fig.vermesh(range(0,size(vtx,1)),vtx,{});
}

//==========================================================================
// [writeimg]
/// Selects what image writer to use based on the file extenstion and then
/// writes the render window to the file. The following formats are supported:
/// BMP, JPEG, PNM, PNG, PostScript, TIFF.
///
/// If no file extension is specified, PNG is assumed. Function derived from:
/// https://kitware.github.io/vtk-examples/site/Cxx/IO/ImageWriter/
///
/// \code{.cpp}
///     matrix<> X = linspace(0,10,100);
///     matrix<> Y = cos(X);
///
///     figure fig;
///     plot(fig,X,Y,{"r-x"});
///     std::vector<std::string> ext = {{""},{".png"}, {".jpg"}, {".ps"},{".tiff"}, {".bmp"}, {".pnm"}};
///     for (int i=0; i<ext.size(); ++i)
///     {
///         writeimg(fig,"output"+ext[i]);
///     }
///     drawnow(fig);
/// \endcode
///
// \see figure.
inline void writeimg(figure const& fig, std::string const& filename)
{
    fig.writeimg(filename);
    fig.writeimg(filename);
}

//==========================================================================
// [xlim]
/// Set x-axis limits
///
/// xlim(limits) specifies the x-axis limits for the current axes. Specify
/// limits as a two-element vector of the form [xmin xmax], where xmax is a
/// numeric value greater than xmin.
inline void xlim(figure& fig, matrix<double>const& limits)
{
    fig.xlim(limits);
}

//==========================================================================
// [ylim]
/// Set y-axis limits
///
/// ylim(limits) specifies the y-axis limits for the current axes. Specify
/// limits as a two-element vector of the form [ymin ymax], where ymax is a
/// numeric value greater than ymin.
inline void ylim(figure& fig, matrix<double>const& limits)
{
    fig.ylim(limits);
}

// End of namespace
}

