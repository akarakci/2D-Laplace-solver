////////////////////////////////////////
// Ata Karakci
// December 15, 2010
// ENGN2912B
//
// This file contains the laplace_solver class
////////////////////////////////////////

#ifndef laplace_solver_h_
#define laplace_solver_h_

#include <mesh2d/v_mesh.h>
#include <mesh2d/potential_field.h>

class laplace_solver
{
public:

	laplace_solver(){}
	
	~laplace_solver(){}
	
	void set_boundary_conditions(vector<double>&);
	bool solve(mesh_boundary const&, vector<double>&, double, double);
	double update_phi(unsigned);

	/* vrml functions */
	void face_color(vtriangle vt, double& r, double& g, double& b);
	void vertex_color(unsigned i, double& r, double& g, double& b);
	void write_potential(ofstream&);
	void write_nodes(ofstream&);

	v_mesh vm_;
	potential_field field_;
};



#endif //laplace_solver_h_