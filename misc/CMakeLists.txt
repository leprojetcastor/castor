source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES perf_hmat.cpp some_fractals.cpp)

add_executable(perf_hmat perf_hmat.cpp)
target_link_libraries(perf_hmat castor)
set_target_properties(perf_hmat PROPERTIES FOLDER "misc")

if(VTK_FOUND)
  add_executable(some_fractals some_fractals.cpp)
  target_link_libraries(some_fractals castor)
  set_target_properties(some_fractals PROPERTIES FOLDER "misc")
  vtk_module_autoinit(TARGETS some_fractals MODULES ${VTK_LIBRARIES})
endif()
