source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES demo_matrix.cpp demo_smatrix.cpp demo_hmatrix.cpp 
                                                    demo_kissfft.cpp demo_linalg.cpp demo_graphics.cpp)

add_executable(demo_matrix demo_matrix.cpp)
target_link_libraries(demo_matrix castor)
set_target_properties(demo_matrix PROPERTIES FOLDER "demo")

add_executable(demo_smatrix demo_smatrix.cpp)
target_link_libraries(demo_smatrix castor)
set_target_properties(demo_smatrix PROPERTIES FOLDER "demo")

add_executable(demo_hmatrix demo_hmatrix.cpp)
target_link_libraries(demo_hmatrix castor)
set_target_properties(demo_hmatrix PROPERTIES FOLDER "demo")

add_executable(demo_kissfft demo_kissfft.cpp)
target_link_libraries(demo_kissfft castor)
set_target_properties(demo_kissfft PROPERTIES FOLDER "demo")

if(LAPACK_FOUND)
  add_executable(demo_linalg demo_linalg.cpp)
  target_link_libraries(demo_linalg castor)
  set_target_properties(demo_linalg PROPERTIES FOLDER "demo")
endif()

if(VTK_FOUND)
  add_executable(demo_graphics demo_graphics.cpp)
  target_link_libraries(demo_graphics castor)
  set_target_properties(demo_graphics PROPERTIES FOLDER "demo")
  vtk_module_autoinit(
    TARGETS demo_graphics
    MODULES ${VTK_LIBRARIES}
  )
endif()
