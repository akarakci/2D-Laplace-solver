////////////////////////////////////////
// Ata Karakci
// December 15, 2010
// ENGN2912B
//
// This file contains the implementation of mesh_generator test
////////////////////////////////////////

#include <testlib/testlib_test.h>
#include <mesh2d/karamete_mesh_generator.h>
#include <fstream>
#include <iostream>
using namespace std;

static void test_mesh_generator()
{
	mesh_boundary mb;
	ifstream ist("circle.dat");
	ist >> mb;
	karamete_mesh_generator kmg;
	kmg.set_boundary(mb);
	kmg.generate(0.4, 0.4);
	v_mesh vm = kmg.mesh();
	
	ofstream ost("mesh.wrl");
	vm.write_triangles(ost);
	ost.close();
}

TESTMAIN(test_mesh_generator);




