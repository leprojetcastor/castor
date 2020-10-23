set(CMAKE_CXX_STANDARD 14)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES overview_smatrix.cpp)

add_executable(overview_smatrix overview_smatrix.cpp)
target_link_libraries(overview_smatrix smatrix matrix)

set_target_properties(overview_smatrix PROPERTIES FOLDER "demo")

