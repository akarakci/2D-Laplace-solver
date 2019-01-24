////////////////////////////////////////
// Ata Karakci
// December 15, 2010
// ENGN2912B
//
// This file contains the v_mesh class
////////////////////////////////////////

#ifndef v_mesh_h_
#define v_mesh_h_

#include <vector>
#include <set>
#include <mesh2d/mesh_boundary.h>
#include <mesh2d/mesh_vertex.h>
#include <mesh2d/vtriangle.h>
#include <fstream>

class v_mesh
{
public:

	v_mesh(){}
	
	~v_mesh(){}
	
	void resize(unsigned n_points, unsigned n_triangles);
	
	void set_boundary(mesh_boundary const& mb){bnd_ = mb;}
	void set_vertex(unsigned i, mesh_vertex const& v);
	void set_triangle(unsigned i, vtriangle const& t);
	
	void add_neighbor(unsigned i, unsigned j){nbrs_[i].insert(j);}
	mesh_vertex& vert(unsigned i){return verts_[i];}
	vtriangle& tri(unsigned i){return tri_[i];}
	unsigned n_verts(){return verts_.size();}
	unsigned n_tri(){return tri_.size();}
	
	unsigned n_neighbors(unsigned i){return nbrs_[i].size();}
	unsigned neighbor(unsigned i, unsigned k);
	
	/* vector of indices of vertices of the boundary */
	vector<unsigned> bnd_verts();
	
	void write_triangles(ofstream&);
	void write_sphere(ofstream&, double, double, double, double, const float,
					  const float, const float, const float);
	void write_verts(ofstream&);
	
private:
	
	vector<mesh_vertex> verts_;
	vector<set<unsigned> > nbrs_;
	vector<vtriangle> tri_;
	mesh_boundary bnd_;
};



#endif //v_mesh_h_