source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES main.cpp)

add_executable(main main.cpp)
target_link_libraries(main castor)
##set_target_properties(main PROPERTIES FOLDER "getting_started")

if(VTK_FOUND)
  vtk_module_autoinit(TARGETS main MODULES ${VTK_LIBRARIES})
endif()
