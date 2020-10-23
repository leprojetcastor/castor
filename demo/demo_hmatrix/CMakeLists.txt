set(CMAKE_CXX_STANDARD 14)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES overview_hmatrix.cpp)

add_executable(overview_hmatrix overview_hmatrix.cpp)
target_link_libraries(overview_hmatrix linalg hmatrix matrix ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

set_target_properties(overview_hmatrix PROPERTIES FOLDER "demo")

