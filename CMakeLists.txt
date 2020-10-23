cmake_minimum_required(VERSION 3.1)

project(castor LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_library(matrix INTERFACE)
target_include_directories(matrix INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

add_library(smatrix INTERFACE)
target_include_directories(smatrix INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

add_library(hmatrix INTERFACE)
target_include_directories(hmatrix INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

find_package(LAPACK)
if(LAPACK_FOUND)
  message(STATUS "BLAS_LIBRARIES : " ${BLAS_LIBRARIES})
  add_library(linalg INTERFACE)
  if(BLAS_LIBRARIES MATCHES "Accelerate")
    target_compile_definitions(linalg INTERFACE -D__CLAPACK_H)
  else()
    find_path(CBLAS_INCLUDE_DIR cblas.h) 
  endif()
  target_include_directories(linalg INTERFACE
     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
     $<BUILD_INTERFACE:${CBLAS_INCLUDE_DIR}>
     $<INSTALL_INTERFACE:include>)  
else()
  message(WARNING "BLAS or LAPACK NOT FOUND: linalg library will not be used")
endif()

find_package(VTK 8.2.0 EXACT NO_MODULE)
if(VTK_FOUND)
  message(STATUS "Find VTK version : ${VTK_VERSION}")
  message(STATUS "VTK_INCLUDE_DIRS : ${VTK_INCLUDE_DIRS}")
  message(STATUS "VTK_LIBRARIES    : ${VTK_LIBRARIES}")
  message(STATUS "VTK_USE_FILE     : ${VTK_USE_FILE}")
  add_library(graphics INTERFACE)
  target_include_directories(graphics INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${VTK_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>)
else()
  message(WARNING "VTK NOT FOUND: graphics library will not be used")
endif()

add_subdirectory(include)
add_subdirectory(extern)
add_subdirectory(test)
add_subdirectory(demo)

# install
include(GNUInstallDirs)
install(FILES  ${CMAKE_CURRENT_SOURCE_DIR}/include/castor/matrix.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/castor)
install(FILES  ${CMAKE_CURRENT_SOURCE_DIR}/include/castor/smatrix.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/castor)
if(LAPACK_FOUND)
  install(FILES  ${CMAKE_CURRENT_SOURCE_DIR}/include/castor/linalg.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/castor)
  install(FILES  ${CMAKE_CURRENT_SOURCE_DIR}/include/castor/hmatrix.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/castor)
endif()
if(VTK_FOUND)
  install(FILES  ${CMAKE_CURRENT_SOURCE_DIR}/include/castor/graphics.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/castor)
endif()
