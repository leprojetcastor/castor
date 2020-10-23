add_subdirectory(demo_matrix)

add_subdirectory(demo_smatrix)

add_subdirectory(demo_hmatrix)

if(LAPACK_FOUND)
    add_subdirectory(demo_linalg)
endif()

if(VTK_FOUND)
     add_subdirectory(demo_graphics)
endif()
