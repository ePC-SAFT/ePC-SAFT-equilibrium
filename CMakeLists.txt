cmake_minimum_required(VERSION 3.20)
project(epcsaft_equilibrium LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)

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
  cpp/src/saturation.cpp
  cpp/src/two_phase_flash.cpp
)
target_compile_features(_equilibrium PRIVATE cxx_std_17)
target_include_directories(_equilibrium PRIVATE "${epcsaft_include_dir}")
target_link_libraries(_equilibrium PRIVATE PkgConfig::Ipopt)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(_equilibrium PRIVATE -Wall -Wextra -Wpedantic)
endif()

install(TARGETS _equilibrium LIBRARY DESTINATION epcsaft_equilibrium)
