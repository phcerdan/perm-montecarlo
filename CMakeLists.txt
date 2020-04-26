cmake_minimum_required(VERSION 3.8.1)

if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

project(PERM
  VERSION 0.1.0
  LANGUAGES CXX)
message(STATUS "PERM version: ${PERM_VERSION}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(GNUInstallDirs) # Define CMAKE_INSTALL_xxx: LIBDIR, INCLUDEDIR

# Set a default build type if none was specified
set(_default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${_default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${_default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(PERM_BUILD_TESTING "Enable Tests" OFF)
option(PERM_BUILD_ENABLE_VALGRIND "Enable Valgrind as a memchecker for tests (require debug symbols)" OFF)
mark_as_advanced(PERM_BUILD_ENABLE_VALGRIND)
option(PERM_ENABLE_GPU "Enable GPU optimization (CUDA)" OFF)
option(PERM_WRAP_PYTHON "Use pybind11 to create bindings to python." OFF)
if(NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
endif()

set(perm_export_file "${PROJECT_BINARY_DIR}/PERMTargets.cmake")

if(PERM_BUILD_TESTING)
  if(PERM_BUILD_ENABLE_VALGRIND)
    find_program(MEMORYCHECK_COMMAND valgrind)
    set(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full")
    set(MEMORYCHECK_SUPPRESSIONS_FILE "${CMAKE_CURRENT_BINARY_DIR}/valgrind_suppress.txt")
    include(Dart)
  endif()
  enable_testing()
  set(GTEST_LIBRARIES gtest gtest_main gmock)
  include(GoogleTest)
  #############################################################################
  # Fetch GTest
  include(FetchContent)

  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.8.x
  )

  set(CMAKE_POLICY_DEFAULT_CMP0048 NEW) # google test raises warning about it
  FetchContent_GetProperties(googletest)
  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
  endif()
  #############################################################################
  # add_subdirectory(test/fixtures)
endif()

if(PERM_WRAP_PYTHON)
  # Fetch pybind11
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.3
  )
  FetchContent_GetProperties(pybind11)
  if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
  endif()
endif()

include(PERMModuleMacros) # module/test macros

# Code
add_subdirectory(src)

if(PERM_BUILD_TESTING)
  add_subdirectory(test)
endif()

if(PERM_WRAP_PYTHON)
  add_subdirectory(wrap)
endif()

# Install/package
set(install_cmake_dir "${CMAKE_INSTALL_LIBDIR}/cmake/perm")
install (EXPORT PERMTargets
  NAMESPACE PERM::
  DESTINATION ${install_cmake_dir} )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/PERMConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/PERMConfigVersion.cmake
              DESTINATION ${install_cmake_dir} )

include(CMakePackageConfigHelpers)

write_basic_package_version_file(PERMConfigVersion.cmake
  VERSION ${PERM_VERSION}
  COMPATIBILITY SameMajorVersion)

# Build tree
set(PERM_TARGETS_FILE ${perm_export_file})
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/PERMConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/PERMConfig.cmake
  INSTALL_DESTINATION ${install_cmake_dir}
  PATH_VARS PERM_TARGETS_FILE
  NO_CHECK_REQUIRED_COMPONENTS_MACRO # PERM does not provide components
  INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}
  )

# Install tree
set(PERM_TARGETS_FILE ${CMAKE_INSTALL_PREFIX}/${install_cmake_dir}/PERMTargets.cmake)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/PERMConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/PERMConfig.cmake
  INSTALL_DESTINATION ${install_cmake_dir}
  PATH_VARS PERM_TARGETS_FILE
  NO_CHECK_REQUIRED_COMPONENTS_MACRO # PERM does not provide components
  )

