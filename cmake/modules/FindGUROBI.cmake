# Expand SDK directories: find_path does not expand shell wildcards.
file(GLOB _GUROBI_ROOTS /opt/gurobi* /opt/gurobi*/linux64)
find_path(GUROBI_INCLUDE_DIRS
    NAMES gurobi_c++.h
    HINTS "$ENV{GUROBI_HOME}" ${_GUROBI_ROOTS}
    PATH_SUFFIXES include
)

set(GUROBI_VERSION "")
if(EXISTS "${GUROBI_INCLUDE_DIRS}/gurobi_c.h")
    file(READ "${GUROBI_INCLUDE_DIRS}/gurobi_c.h" _GUROBI_HEADER)
    foreach(component MAJOR MINOR TECHNICAL)
        string(REGEX MATCH "#[ \t]*define[ \t]+GRB_VERSION_${component}[ \t]+([0-9]+)" _match "${_GUROBI_HEADER}")
        set(_GUROBI_${component} "${CMAKE_MATCH_1}")
    endforeach()
    if(NOT _GUROBI_MAJOR STREQUAL "" AND NOT _GUROBI_MINOR STREQUAL "" AND NOT _GUROBI_TECHNICAL STREQUAL "")
        set(GUROBI_VERSION "${_GUROBI_MAJOR}.${_GUROBI_MINOR}.${_GUROBI_TECHNICAL}")
    endif()
endif()
get_filename_component(_GUROBI_SDK "${GUROBI_INCLUDE_DIRS}" DIRECTORY)
find_library(GUROBI_LIBRARY
    NAMES gurobi${_GUROBI_MAJOR}${_GUROBI_MINOR} gurobi
    PATHS "${_GUROBI_SDK}/lib"
    NO_DEFAULT_PATH
)

find_library(GUROBI_CXX_LIBRARY
    NAMES gurobi_c++
    PATHS "${_GUROBI_SDK}/lib"
    NO_DEFAULT_PATH
)

# Cached paths survive SDK changes; reject a core library naming another version.
set(GUROBI_LIBRARY_VERSION_MATCH TRUE)
get_filename_component(_GUROBI_CORE "${GUROBI_LIBRARY}" REALPATH)
get_filename_component(_GUROBI_CORE "${_GUROBI_CORE}" NAME)
if(_GUROBI_CORE MATCHES "^(lib)?gurobi([0-9]+)[.]")
    if(NOT CMAKE_MATCH_2 STREQUAL "${_GUROBI_MAJOR}${_GUROBI_MINOR}")
        set(GUROBI_LIBRARY_VERSION_MATCH FALSE)
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI
    REQUIRED_VARS GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_INCLUDE_DIRS
        GUROBI_VERSION GUROBI_LIBRARY_VERSION_MATCH
    VERSION_VAR GUROBI_VERSION
    REASON_FAILURE_MESSAGE "Select one SDK with matching headers and libraries. Clear cached GUROBI paths when changing SDKs.")

if(GUROBI_FOUND AND NOT TARGET GUROBI::GUROBI)
    add_library(GUROBI::GUROBI INTERFACE IMPORTED)
    set_target_properties(GUROBI::GUROBI PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}"
    )
endif()
