########################################
# Ata Karakci
#
########################################

# This is CMakeLists.txt

cmake_minimum_required(VERSION 2.8)

PROJECT(mesh)

FIND_PACKAGE(VXL)

IF(VXL_FOUND)
  # include the standard VXL stuff
  SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
  INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
  INCLUDE_DIRECTORIES(${mesh_SOURCE_DIR})

  SUBDIRS(mesh2d)

ENDIF(VXL_FOUND)

