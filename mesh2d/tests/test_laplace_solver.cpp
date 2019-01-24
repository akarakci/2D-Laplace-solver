////////////////////////////////////////
// Ata Karakci
// December 15, 2010
// ENGN2912B
//
// This file contains the implementation of test_laplace_solver test
////////////////////////////////////////

#include <testlib/testlib_test.h>
#include <mesh2d/laplace_solver.h>
#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

static void test_laplace_solver()
{
	/* boundary geometry */
	mesh_boundary mb;
	ifstream ist("circle.dat");
	ist >> mb;

	/* boundary conditions */
	vector<double> bc;
	for(unsigned i=0; i<4; ++i)
		bc.push_back(10);
	bc[1]=20; bc[2]=20;
	
	laplace_solver ls;

	bool bl = ls.solve(mb, bc, 0.1, 0.1);
	if(bl)
	{
		ofstream omesh("laplace_mesh.wrl");
		ls.write_potential(omesh);
		omesh.close();
		
		ofstream onode("laplace_node.wrl");
		ls.write_nodes(onode);
		onode.close();
	}
	else
		cout << "cannot solve PDE" << endl;

}		

TESTMAIN(test_laplace_solver);




