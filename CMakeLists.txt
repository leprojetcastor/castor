cmake_minimum_required(VERSION 3.18)

project(castor LANGUAGES CXX VERSION 1.0.0)

set(CASTOR_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

message(STATUS "Building castor v${${PROJECT_NAME}_VERSION}")

###############################
# dependencies

find_package(LAPACK)
if(LAPACK_FOUND)
  message(STATUS "BLAS and LAPACK found")
else()
  message(STATUS "BLAS and LAPACK not found")
  set(BLAS_LIBRARIES "")
endif()
message(STATUS "BLAS_LIBRARIES : " ${BLAS_LIBRARIES})
message(STATUS "LAPACK_LIBRARIES : " ${LAPACK_LIBRARIES})

find_package(VTK COMPONENTS
  vtkChartsCore
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkInteractionStyle
  vtkRenderingContext2D
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkViewsContext2D

  vtkRenderingAnnotation
  vtkIOImage
  vtkIOPLY
  vtkIOLegacy
  vtkIOOggTheora
  vtkFiltersGeometry
  QUIET
)
if(VTK_FOUND)
  message(STATUS "Find VTK version : ${VTK_VERSION}")
else()
  message(STATUS "VTK not found")
endif()
message(STATUS "VTK_LIBRARIES    : ${VTK_LIBRARIES}")

###############################
# target

add_library(castor INTERFACE)

if(LAPACK_FOUND)
  if(BLAS_LIBRARIES MATCHES "Accelerate")
    target_compile_definitions(castor INTERFACE -D__CLAPACK_H)
  else()
    find_path(CBLAS_INCLUDE_DIR cblas.h)
  endif()
endif()

target_include_directories(castor INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CBLAS_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${VTK_INCLUDE_DIRS}>
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:${CBLAS_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${VTK_INCLUDE_DIR}>)

target_link_libraries(castor INTERFACE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${VTK_LIBRARIES})

target_compile_features(castor INTERFACE cxx_std_14)

add_subdirectory(include)
add_subdirectory(extern)
add_subdirectory(test)
add_subdirectory(demo)
add_subdirectory(misc)
add_subdirectory(getting_started)

###############################
# installation

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS castor
        EXPORT ${PROJECT_NAME}-targets)

# Makes the project importable from the build directory
export(EXPORT ${PROJECT_NAME}-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")

install(DIRECTORY ${CASTOR_INCLUDE_DIRS}/castor
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
####install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/castor/matrix.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/castor)

set(CASTOR_CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for matrixConfig.cmake")

configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION ${CASTOR_CMAKECONFIG_INSTALL_DIR})

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
                                 VERSION ${${PROJECT_NAME}_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CASTOR_CMAKECONFIG_INSTALL_DIR})
install(EXPORT ${PROJECT_NAME}-targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${CASTOR_CMAKECONFIG_INSTALL_DIR})

configure_file(${PROJECT_NAME}.pc.in
               "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
                @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")
