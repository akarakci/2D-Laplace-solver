////////////////////////////////////////
// Ata Karakci
// December 12, 2010
// ENGN2912B
//
// This file contains the karamete_mesh_generator class
////////////////////////////////////////

#ifndef karamete_mesh_generator_h_
#define karamete_mesh_generator_h_

#include <mesh2d/mesh_generator.h>

class karamete_mesh_generator : public mesh_generator
{
public:

	karamete_mesh_generator();
	
	void set_boundary (mesh_boundary const& mb);
	bool generate(double alpha, double beta);
	
private:
	
	void set_neighbors();
};



#endif //karamete_mesh_generator_h_