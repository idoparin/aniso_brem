cmake_minimum_required(VERSION 3.17.2...3.30)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

# ---------------------------------------------------------------------------
# Shared cross-section objects. cseqs.f95 is compiled exactly once here and the
# resulting `cseqs` Fortran module (.mod) plus object code are reused by both
# Python extensions below. This avoids compiling cseqs.f95 twice (which races
# on the cseqs.mod file) and keeps a single source of truth for the analytic
# bremsstrahlung cross sections.
# ---------------------------------------------------------------------------
set(CSEQS_MOD_DIR "${CMAKE_CURRENT_BINARY_DIR}/cseqs_mod")
file(MAKE_DIRECTORY "${CSEQS_MOD_DIR}")

add_library(cseqs_obj OBJECT
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/cseqs.f95"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/dilog490.f"
)
set_target_properties(cseqs_obj PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY "${CSEQS_MOD_DIR}"
)
# Expose the generated cseqs.mod to anything that links cseqs_obj.
target_include_directories(cseqs_obj PUBLIC "${CSEQS_MOD_DIR}")

# ---------------------------------------------------------------------------
# Extension 1: cseqs -- the cross-section routines themselves.
# ---------------------------------------------------------------------------
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/cseqsmodule.c"
         "${CMAKE_CURRENT_BINARY_DIR}/cseqs-f2pywrappers.f"
         "${CMAKE_CURRENT_BINARY_DIR}/cseqs-f2pywrappers2.f90"
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/fortran/cseqs.f95"
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/dilog490.f"
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/cseqs.f95"
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/dilog490.f"
          -m cseqs --lower
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

python_add_library(
  cseqs MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/cseqs-f2pywrappers2.f90"
  "${CMAKE_CURRENT_BINARY_DIR}/cseqs-f2pywrappers.f"
  "${CMAKE_CURRENT_BINARY_DIR}/cseqsmodule.c"
  WITH_SOABI
)
# The f2py F90 wrappers do `use cseqs`, so they need cseqs.mod + the objects.
target_link_libraries(cseqs PRIVATE cseqs_obj fortranobject)
install(TARGETS cseqs DESTINATION aniso_brem)

# ---------------------------------------------------------------------------
# Extension 2: _flux -- the bremsstrahlung flux routines (module `flux`).
# f2py only wraps `flux` here because it `use`s `cseqs`; the cseqs procedures
# come from the shared cseqs_obj objects linked below.
# ---------------------------------------------------------------------------
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/_fluxmodule.c"
         "${CMAKE_CURRENT_BINARY_DIR}/_flux-f2pywrappers.f"
         "${CMAKE_CURRENT_BINARY_DIR}/_flux-f2pywrappers2.f90"
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/fortran/cseqs.f95"
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/flux.f95"
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/flux_aniso.f95"
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/dilog490.f"
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/cseqs.f95"
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/flux.f95"
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/flux_aniso.f95"
          "${CMAKE_CURRENT_SOURCE_DIR}/fortran/dilog490.f"
          -m _flux --lower
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

python_add_library(
  _flux MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/_flux-f2pywrappers2.f90"
  "${CMAKE_CURRENT_BINARY_DIR}/_flux-f2pywrappers.f"
  "${CMAKE_CURRENT_BINARY_DIR}/_fluxmodule.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/flux.f95"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/flux_aniso.f95"
  WITH_SOABI
)
set_target_properties(_flux PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_flux_mod"
)
# flux.f95 does `use cseqs`; cseqs_obj (PUBLIC include) provides cseqs.mod and
# the cross-section object code.
target_link_libraries(_flux PRIVATE cseqs_obj fortranobject)
# Ensure cseqs_obj (and thus a complete cseqs.mod) is fully built before any of
# _flux's own Fortran sources are compiled.
add_dependencies(_flux cseqs_obj)
install(TARGETS _flux DESTINATION aniso_brem)
