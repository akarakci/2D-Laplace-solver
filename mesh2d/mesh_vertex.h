////////////////////////////////////////
// Ata Karakci
// December 5, 2010
// ENGN2912B
//
// This file contains the mesh_vertex class
////////////////////////////////////////

#ifndef mesh_vertex_h_
#define mesh_vertex_h_

class mesh_vertex
{
public:

	mesh_vertex() : x_(0), y_(0) {}
	mesh_vertex(double x, double y) : x_(x), y_(y) {}
	~mesh_vertex(){}
	
	bool operator == (mesh_vertex const& mv)
	{
		if(this->x_==mv.x_ && this->y_==mv.y_)
			return true;
		return false;
	}

	double x_;
	double y_;
};



#endif //mesh_vertex_h_