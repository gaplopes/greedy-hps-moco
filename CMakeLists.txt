cmake_minimum_required(VERSION 3.18)

project(greedy-hps-moco LANGUAGES CXX)

option(GHPS_BUILD_APPS "Build the solver executable" ON)
option(GHPS_FETCH_INSTANCES "Download the optional experiment dataset" OFF)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Set the build type to Release by default
if(NOT CMAKE_BUILD_TYPE)
  # set(CMAKE_BUILD_TYPE Debug)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Print the build type
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Warning flags
set(GREEDY_HPS_CXX_WARN_FLAGS "")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  list(APPEND GREEDY_HPS_CXX_WARN_FLAGS
    -Wall
    -Wextra
    -Wpedantic
    # -pedantic-errors
    # -Wold-style-cast
    # -Wundef
    # -Wredundant-decls
    # -Wcast-qual
    # -Wmissing-include-dirs
    # -Wdisabled-optimization
    # -Wconversion
    # -Wdouble-promotion
    # -Wshadow
    # -Wnull-dereference
  )
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    list(APPEND GREEDY_HPS_CXX_WARN_FLAGS
      -Wvector-operation-performance
      -Wduplicated-cond
    )
  endif()
else()
  message(WARNING "Unknown compiler (not GNU, Clang, or MSVC), not setting any warning flags.")
  set(GREEDY_HPS_CXX_WARN_FLAGS "")
endif()

# Find dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

# Include FetchContent module
include(FetchContent)

find_package(OpenMP REQUIRED COMPONENTS CXX)

find_package(mooutils CONFIG QUIET)
# If mooutils not found fallback to github
if (NOT mooutils_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    mooutils
    GIT_REPOSITORY https://github.com/adbjesus/mooutils.git
    GIT_TAG        0ed6dc52b88483be1f80f4417641814de60ae959
  )
  FetchContent_MakeAvailable(mooutils)
  add_library(mooutils::mooutils ALIAS mooutils)
endif()

find_package(GUROBI 12.0 REQUIRED)

# Fetch mobkp-instances
if(GHPS_FETCH_INSTANCES)
  FetchContent_Declare(
    mobkp_instances
    GIT_REPOSITORY https://github.com/gaplopes/mobkp-instances.git
    GIT_TAG        1e053bbcafb7792b7778dce222434038880cfb49
  )
  FetchContent_GetProperties(mobkp_instances)
  if(NOT mobkp_instances_POPULATED)
    FetchContent_Populate(mobkp_instances)
  endif()
  message(STATUS "Experiment instances: ${mobkp_instances_SOURCE_DIR}/instances")
endif()

# Create modern INTERFACE library
add_library(moco-hps INTERFACE)
add_library(moco::ghps ALIAS moco-hps)
set_target_properties(moco-hps PROPERTIES EXPORT_NAME ghps)
target_compile_features(moco-hps INTERFACE cxx_std_20)

# Add include directories to the library
target_include_directories(moco-hps INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Link dependencies to the library
target_link_libraries(moco-hps INTERFACE
  mooutils::mooutils
  OpenMP::OpenMP_CXX
  GUROBI::GUROBI
)

install(TARGETS moco-hps EXPORT moco-hps-targets)
install(DIRECTORY include/ghps DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT moco-hps-targets NAMESPACE moco::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/moco-hps)
configure_package_config_file(cmake/moco-hps-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/moco-hps-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/moco-hps)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/moco-hps-config.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/moco-hps)
install(FILES cmake/modules/FindGUROBI.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/moco-hps/modules)

if(GHPS_BUILD_APPS)
  find_package(CLI11 CONFIG QUIET)
  if(NOT CLI11_FOUND)
    FetchContent_Declare(
      cli11
      GIT_REPOSITORY https://github.com/CLIUtils/CLI11
      GIT_TAG        6c7b07a878ad834957b98d0f9ce1dbe0cb204fc9
    )
    FetchContent_MakeAvailable(cli11)
  endif()

  # Add the executable and specify the source files
  add_executable(greedy-hps-moco
    ${CMAKE_CURRENT_SOURCE_DIR}/apps/main.cpp
  )

  # Link the executable with the interface library
  target_link_libraries(greedy-hps-moco
    PRIVATE
      moco-hps
      CLI11::CLI11
  )
  target_compile_options(greedy-hps-moco PRIVATE ${GREEDY_HPS_CXX_WARN_FLAGS})

  # Set the output name of the executable
  set_target_properties(greedy-hps-moco PROPERTIES OUTPUT_NAME greedy-hps-moco)

  # Install the target
  install(TARGETS greedy-hps-moco)
endif()
