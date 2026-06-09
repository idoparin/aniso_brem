cmake_minimum_required(VERSION 3.15...3.30)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

set(f2py_module_name "cseqs")
set(fortran_src_files
    "${CMAKE_CURRENT_SOURCE_DIR}/fortran/cseqs.f95"
    "${CMAKE_CURRENT_SOURCE_DIR}/fortran/dilog490.f"
)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}module.c"
         "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers.f"
         "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers2.f90"
  DEPENDS ${fortran_src_files}
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          ${fortran_src_files}
          -m ${f2py_module_name}
          --lower
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

python_add_library(
  ${f2py_module_name} MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers2.f90"
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers.f"
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}module.c"
  ${fortran_src_files}
  WITH_SOABI
)
target_link_libraries(${f2py_module_name} PRIVATE fortranobject)

install(TARGETS ${f2py_module_name} DESTINATION aniso_brem)
