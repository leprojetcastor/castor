set(CMAKE_CXX_STANDARD 14)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES overview_kissfft.cpp)

add_executable(overview_kissfft overview_kissfft.cpp)
target_link_libraries(overview_kissfft kissfft matrix)

set_target_properties(overview_kissfft PROPERTIES FOLDER "demo")

