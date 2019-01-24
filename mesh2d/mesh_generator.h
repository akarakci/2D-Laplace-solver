////////////////////////////////////////
// Ata Karakci
// December 12, 2010
// ENGN2912B
//
// This file contains the mesh_generator class
////////////////////////////////////////

#ifndef mesh_generator_h_
#define mesh_generator_h_

#include <mesh2d/v_mesh.h>
#include <mesh2d/mesh_boundary.h>
#include <iostream>

class mesh_generator
{
public:

	mesh_generator(){}
	
	~mesh_generator(){}
	
	virtual void set_boundary (mesh_boundary const& mb)=0;
	virtual bool generate(double alpha, double beta)=0;
	v_mesh mesh() {return vm_;}
	
protected:
	
	v_mesh vm_;
};



#endif //mesh_generator_h_