set(CMAKE_CXX_STANDARD 14)

#include(${VTK_USE_FILE})

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES basic.cpp overview_graphics.cpp nodal_values.cpp surface_plot.cpp advanced_mesh.cpp some_fractals.cpp)

add_executable(basic basic.cpp)
target_link_libraries(basic matrix graphics ${VTK_LIBRARIES})

add_executable(overview_graphics overview_graphics.cpp)
target_link_libraries(overview_graphics matrix graphics ${VTK_LIBRARIES})

add_executable(nodal_values nodal_values.cpp)
target_link_libraries(nodal_values matrix graphics ${VTK_LIBRARIES})

add_executable(surface_plot surface_plot.cpp)
target_link_libraries(surface_plot matrix graphics ${VTK_LIBRARIES})

add_executable(advanced_mesh advanced_mesh.cpp)
target_link_libraries(advanced_mesh matrix graphics ${VTK_LIBRARIES})

add_executable(some_fractals some_fractals.cpp)
target_link_libraries(some_fractals matrix graphics ${VTK_LIBRARIES})

vtk_module_autoinit(
  TARGETS overview_graphics
  MODULES ${VTK_LIBRARIES}
)

set_target_properties(basic overview_graphics nodal_values surface_plot advanced_mesh some_fractals PROPERTIES FOLDER "demo")
