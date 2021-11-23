set(CMAKE_CXX_STANDARD 14)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES test_main.cpp test_constructors.cpp test_internal_op.cpp)

add_executable(test test_main.cpp test_constructors.cpp test_internal_op.cpp)
target_link_libraries(test castor catch)
