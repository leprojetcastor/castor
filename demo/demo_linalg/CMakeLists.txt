set(CMAKE_CXX_STANDARD 14)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES lapack_conversion.cpp multiply.cpp linsolve.cpp inv_pinv.cpp lu.cpp eig.cpp qr.cpp svd.cpp rank.cpp aca.cpp qrsvd.cpp overview_linalg.cpp)

add_executable(lapack_conversion lapack_conversion.cpp)
target_link_libraries(lapack_conversion matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(multiply multiply.cpp)
target_link_libraries(multiply matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(linsolve linsolve.cpp)
target_link_libraries(linsolve matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(inv_pinv inv_pinv.cpp)
target_link_libraries(inv_pinv matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(lu lu.cpp)
target_link_libraries(lu matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(eig eig.cpp)
target_link_libraries(eig matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(qr qr.cpp)
target_link_libraries(qr matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(svd svd.cpp)
target_link_libraries(svd matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(rank rank.cpp)
target_link_libraries(rank matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(aca aca.cpp)
target_link_libraries(aca matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(qrsvd qrsvd.cpp)
target_link_libraries(qrsvd matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(overview_linalg overview_linalg.cpp)
target_link_libraries(overview_linalg matrix linalg ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

set_target_properties(lapack_conversion multiply linsolve inv_pinv lu eig qr svd rank aca qrsvd overview_linalg PROPERTIES FOLDER "demo")
