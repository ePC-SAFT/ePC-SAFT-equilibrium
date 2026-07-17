cmake_minimum_required(VERSION 3.20)
project(epcsaft_equilibrium LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)

execute_process(
  COMMAND
    "${Python_EXECUTABLE}" -c
    "import importlib.resources; print(importlib.resources.files('epcsaft') / 'include')"
  RESULT_VARIABLE epcsaft_header_status
  OUTPUT_VARIABLE epcsaft_include_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT epcsaft_header_status EQUAL 0)
  message(FATAL_ERROR "Build requires an installed epcsaft provider wheel")
endif()
if(NOT EXISTS "${epcsaft_include_dir}/epcsaft/native_sdk_v1.h")
  message(FATAL_ERROR "Installed epcsaft wheel does not contain native_sdk_v1.h")
endif()

pybind11_add_module(_equilibrium MODULE cpp/src/module.cpp cpp/src/saturation.cpp)
target_compile_features(_equilibrium PRIVATE cxx_std_17)
target_include_directories(_equilibrium PRIVATE "${epcsaft_include_dir}")
target_link_libraries(_equilibrium PRIVATE PkgConfig::Ipopt)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(_equilibrium PRIVATE -Wall -Wextra -Wpedantic)
endif()

install(TARGETS _equilibrium LIBRARY DESTINATION epcsaft_equilibrium)
