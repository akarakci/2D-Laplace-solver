########################################
# Ata Karakci
#
########################################

# This is mesh2d/CMakeLists.txt

cmake_minimum_required(VERSION 2.8)

SET(mesh2d_sources karamete_mesh_generator.cpp karamete_mesh_generator.h mesh2d.c mesh2d.h mesh_boundary.h mesh_generator.h mesh_vertex.h v_mesh.cpp v_mesh.h vtriangle.h laplace_solver.cpp laplace_solver.h potential_field.h)

ADD_LIBRARY(mesh2d ${mesh2d_sources})
TARGET_LINK_LIBRARIES(mesh2d vnl_algo vnl)

IF( BUILD_TESTING )
  SUBDIRS(tests)
ENDIF( BUILD_TESTING )
