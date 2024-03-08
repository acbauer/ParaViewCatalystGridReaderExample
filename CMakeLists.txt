cmake_minimum_required(VERSION 3.13)
project(CxxDualMeshExample LANGUAGES C CXX)

include (GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

#------------------------------------------------------------------------------
# since we use C++11 in this example.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Since this example uses MPI, find and link against it.
find_package(MPI COMPONENTS C CXX)
if (NOT MPI_FOUND)
  message(STATUS
    "ERROR: Requires MPI.")
  return ()
endif ()

#------------------------------------------------------------------------------
add_executable(CxxDualMeshExample
  FEDriver.cxx)
target_link_libraries(CxxDualMeshExample
  PRIVATE
    MPI::MPI_C
    MPI::MPI_CXX)

#------------------------------------------------------------------------------
find_package(VTK COMPONENTS
  CommonCore IOXML
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "CxxDualMeshExample: Unable to find the VTK build folder.")
endif()

# Assume that the VTK we link to is the same one that is in ParaView. Then
# find the ParaView Catalyst library path relative to VTK_DIR.
find_path(ParaView_CATALYST_DIR libcatalyst-paraview.so HINTS ${VTK_DIR}/lib64/catalyst/ ${VTK_DIR}/../../../lib64/catalyst/)

find_package(catalyst REQUIRED
  PATHS "${ParaView_DIR}/catalyst")
target_compile_definitions(CxxDualMeshExample
  PRIVATE
  "PARAVIEW_IMPL_DIR=\"${ParaView_CATALYST_DIR}\"")

target_link_libraries(CxxDualMeshExample
  PRIVATE
  catalyst::catalyst VTK::IOXML VTK::CommonCore)

vtk_module_autoinit(
  TARGETS CxxDualMeshExample
  MODULES ${VTK_LIBRARIES}
)


