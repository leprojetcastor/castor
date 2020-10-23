set(CMAKE_CXX_STANDARD 14)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES algorithms.cpp builders.cpp complex.cpp constructors.cpp external.cpp geometry.cpp getting_started.cpp input_output.cpp internal.cpp manipulations.cpp mathematics.cpp overview_matrix.cpp view.cpp)

add_executable(algorithms algorithms.cpp)
target_link_libraries(algorithms matrix)

add_executable(builders builders.cpp)
target_link_libraries(builders matrix)

add_executable(complex complex.cpp)
target_link_libraries(complex matrix)

add_executable(constructors constructors.cpp)
target_link_libraries(constructors matrix)

add_executable(external external.cpp)
target_link_libraries(external matrix)

add_executable(geometry geometry.cpp)
target_link_libraries(geometry matrix)

add_executable(getting_started getting_started.cpp)
target_link_libraries(getting_started matrix)

add_executable(input_output input_output.cpp)
target_link_libraries(input_output matrix)

add_executable(internal internal.cpp)
target_link_libraries(internal matrix)

add_executable(manipulations manipulations.cpp)
target_link_libraries(manipulations matrix)

add_executable(mathematics mathematics.cpp)
target_link_libraries(mathematics matrix)

add_executable(overview_matrix overview_matrix.cpp)
target_link_libraries(overview_matrix matrix)

add_executable(view view.cpp)
target_link_libraries(view matrix)

set_target_properties(algorithms builders complex constructors external geometry getting_started input_output internal manipulations mathematics overview_matrix view PROPERTIES FOLDER "demo")

