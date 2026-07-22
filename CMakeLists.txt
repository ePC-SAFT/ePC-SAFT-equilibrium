cmake_minimum_required(VERSION 3.20)
project(epcsaft_equilibrium LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)

include(FetchContent)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_CXX_EXE OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_EXTRAS_LIB OFF CACHE BOOL "" FORCE)
set(ZLIB OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  highs
  URL https://github.com/ERGO-Code/HiGHS/releases/download/v1.15.1/source-archive.tar.gz
  URL_HASH SHA256=30132da2b21899c62f07c4d6566325a3b869a7da4ce61d8b25f833b38e1562b1
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_GetProperties(highs)
if(NOT highs_POPULATED)
  FetchContent_Populate(highs)
  add_subdirectory("${highs_SOURCE_DIR}" "${highs_BINARY_DIR}" EXCLUDE_FROM_ALL)
endif()

set(NLOPT_CXX ON CACHE BOOL "" FORCE)
set(NLOPT_FORTRAN OFF CACHE BOOL "" FORCE)
set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
set(NLOPT_MATLAB OFF CACHE BOOL "" FORCE)
set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
set(NLOPT_JAVA OFF CACHE BOOL "" FORCE)
set(NLOPT_SWIG OFF CACHE BOOL "" FORCE)
set(NLOPT_LUKSAN OFF CACHE BOOL "" FORCE)
set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  nlopt
  URL https://github.com/stevengj/nlopt/archive/refs/tags/v2.11.0.tar.gz
  URL_HASH SHA256=53e552d83e9294d67db37f0f4a23f15933a9ef698485301a18b98b40004cf0de
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_GetProperties(nlopt)
if(NOT nlopt_POPULATED)
  FetchContent_Populate(nlopt)
  add_subdirectory("${nlopt_SOURCE_DIR}" "${nlopt_BINARY_DIR}" EXCLUDE_FROM_ALL)
endif()

execute_process(
  COMMAND
    "${Python_EXECUTABLE}" -c
    "import importlib.metadata as m, json, sys; d=m.distribution('epcsaft'); u=json.loads(d.read_text('direct_url.json') or '{}'); sys.exit('epcsaft version mismatch') if d.version != '0.1.0.dev0' else None; sys.exit('editable epcsaft provider is forbidden') if u.get('dir_info', {}).get('editable', False) else None; print(d.locate_file('epcsaft/include'))"
  RESULT_VARIABLE epcsaft_header_status
  OUTPUT_VARIABLE epcsaft_include_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT epcsaft_header_status EQUAL 0)
  message(FATAL_ERROR "Build requires the non-editable epcsaft 0.1.0.dev0 provider distribution")
endif()
if(NOT EXISTS "${epcsaft_include_dir}/epcsaft/native_sdk_v1.h")
  message(FATAL_ERROR "Installed epcsaft wheel does not contain native_sdk_v1.h")
endif()

pybind11_add_module(
  _equilibrium MODULE
  cpp/src/module.cpp
  cpp/src/provider.cpp
  cpp/src/held.cpp
  cpp/src/held_stage_i.cpp
  cpp/src/held_stage_ii.cpp
  cpp/src/held_stage_iii.cpp
  cpp/src/held2.cpp
  cpp/src/held2_progress.cpp
  cpp/src/held2_stage_i_direct.cpp
  cpp/src/held2_stage_ii_basin.cpp
  cpp/src/held2_stage_ii_upper.cpp
  cpp/src/held2_stage_iii.cpp
  cpp/src/held2_stage_iii_generic.cpp
  cpp/src/saturation.cpp
)
target_compile_features(_equilibrium PRIVATE cxx_std_17)
target_include_directories(_equilibrium PRIVATE "${epcsaft_include_dir}")
target_link_libraries(_equilibrium PRIVATE PkgConfig::Ipopt highs nlopt)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(_equilibrium PRIVATE -Wall -Wextra -Wpedantic)
endif()

install(TARGETS _equilibrium LIBRARY DESTINATION epcsaft_equilibrium)
