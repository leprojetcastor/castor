source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/castor FILES castor/matrix.hpp castor/linalg.hpp castor/graphics.hpp castor/smatrix.hpp castor/hmatrix.hpp)

add_custom_target(include SOURCES castor/matrix.hpp castor/linalg.hpp castor/graphics.hpp castor/smatrix.hpp castor/hmatrix.hpp)
