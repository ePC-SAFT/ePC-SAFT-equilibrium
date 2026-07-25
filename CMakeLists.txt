cmake_minimum_required(VERSION 3.20)
project(epcsaft_equilibrium LANGUAGES CXX)

option(
  EPCSAFT_EQUILIBRIUM_BUILD_DIAGNOSTIC
  "Build the development-only native equilibrium diagnostic executable"
  OFF
)
option(
  EPCSAFT_EQUILIBRIUM_ENABLE_SANITIZERS
  "Enable AddressSanitizer and UndefinedBehaviorSanitizer"
  OFF
)

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
    "import importlib.metadata as m, json, sys; d=m.distribution('epcsaft'); u=json.loads(d.read_text('direct_url.json') or '{}'); sys.exit('epcsaft version mismatch') if d.version != '0.1.0.dev0' else None; sys.exit('editable epcsaft provider is forbidden') if u.get('dir_info', {}).get('editable', False) else None; print(d.locate_file('epcsaft/cmake'))"
  RESULT_VARIABLE epcsaft_header_status
  OUTPUT_VARIABLE epcsaft_cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT epcsaft_header_status EQUAL 0)
  message(FATAL_ERROR "Build requires the non-editable epcsaft 0.1.0.dev0 provider distribution")
endif()
find_package(epcsaft CONFIG REQUIRED PATHS "${epcsaft_cmake_dir}" NO_DEFAULT_PATH)

add_library(
  epcsaft_equilibrium_core STATIC
  cpp/src/provider.cpp
  cpp/src/flash.cpp
  cpp/src/result_json.cpp
  cpp/src/held.cpp
  cpp/src/held_stage_i.cpp
  cpp/src/held_stage_ii.cpp
  cpp/src/held_stage_iii.cpp
  cpp/src/held2.cpp
  cpp/src/held2_controller.cpp
  cpp/src/held2_progress.cpp
  cpp/src/held2_step1.cpp
  cpp/src/held2_step2.cpp
  cpp/src/held2_step3.cpp
  cpp/src/held2_step4.cpp
  cpp/src/held2_step5.cpp
  cpp/src/held2_stage_iii.cpp
  cpp/src/held2_stage_iii_generic.cpp
  cpp/src/saturation.cpp
)
set_target_properties(epcsaft_equilibrium_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(epcsaft_equilibrium_core PUBLIC cxx_std_17)
target_include_directories(
  epcsaft_equilibrium_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cpp/src"
)
target_link_libraries(
  epcsaft_equilibrium_core
  PUBLIC epcsaft::native_sdk
  PRIVATE PkgConfig::Ipopt highs nlopt
)

pybind11_add_module(_equilibrium MODULE cpp/src/module.cpp)
target_link_libraries(_equilibrium PRIVATE epcsaft_equilibrium_core)
set(epcsaft_equilibrium_native_targets epcsaft_equilibrium_core _equilibrium)
set(epcsaft_equilibrium_link_targets _equilibrium)
if(EPCSAFT_EQUILIBRIUM_BUILD_DIAGNOSTIC)
  add_executable(
    epcsaft-equilibrium-diagnostic cpp/src/diagnostic_main.cpp
  )
  target_link_libraries(
    epcsaft-equilibrium-diagnostic PRIVATE epcsaft_equilibrium_core
  )
  add_executable(
    epcsaft-equilibrium-manufactured-check
    cpp/src/manufactured_check_main.cpp
    cpp/tests/held2_step_checks.cpp
  )
  target_link_libraries(
    epcsaft-equilibrium-manufactured-check
    PRIVATE epcsaft_equilibrium_core
  )
  list(
    APPEND epcsaft_equilibrium_native_targets
    epcsaft-equilibrium-diagnostic
    epcsaft-equilibrium-manufactured-check
  )
  list(
    APPEND epcsaft_equilibrium_link_targets
    epcsaft-equilibrium-diagnostic
    epcsaft-equilibrium-manufactured-check
  )
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  foreach(target IN LISTS epcsaft_equilibrium_native_targets)
    target_compile_options(
      ${target} PRIVATE -Wall -Wextra -Wpedantic
    )
    if(EPCSAFT_EQUILIBRIUM_ENABLE_SANITIZERS)
      target_compile_options(
        ${target}
        PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer
      )
    endif()
  endforeach()
  if(EPCSAFT_EQUILIBRIUM_ENABLE_SANITIZERS)
    foreach(target IN LISTS epcsaft_equilibrium_link_targets)
      target_link_options(
        ${target} PRIVATE -fsanitize=address,undefined
      )
    endforeach()
  endif()
endif()

install(TARGETS _equilibrium LIBRARY DESTINATION epcsaft_equilibrium)
