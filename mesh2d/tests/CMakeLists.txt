########################################
# Ata Karakci
#
########################################

# This is mesh2d/tests/CMakeLists.txt

cmake_minimum_required(VERSION 2.8)

ADD_EXECUTABLE( mesh_generator_test
  # Driver source and utilities
  test_driver_mesh_generator.cpp

  # The tests
  test_mesh_generator.cpp)

ADD_EXECUTABLE( laplace_solver_test
  # Driver source and utilities
  test_driver_laplace_solver.cpp

  # The tests
  test_laplace_solver.cpp)

TARGET_LINK_LIBRARIES(laplace_solver_test mesh2d testlib)

ADD_TEST (laplace_solver_test test_laplace_solver)

TARGET_LINK_LIBRARIES(mesh_generator_test mesh2d testlib)

ADD_TEST (mesh_generator_test test_laplace_solver)
