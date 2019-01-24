////////////////////////////////////////
// Ata Karakci
// December 15, 2010
// ENGN2912B
//
// This file contains the v_mesh class
////////////////////////////////////////

#include <mesh2d/v_mesh.h>

	
void v_mesh::resize(unsigned n_points, unsigned n_triangles)
{
	verts_.resize(n_points);
	nbrs_.resize(n_points);
	tri_.resize(n_triangles);
}

void v_mesh::set_vertex(unsigned i, mesh_vertex const& v)
{
	verts_[i].x_ = v.x_; 
	verts_[i].y_ = v.y_;
}

void v_mesh::set_triangle(unsigned i, vtriangle const& t)
{
	for(unsigned j=0; j<3; ++j)
	tri_[i].set_vert(j, t.vert(j));
}
	
unsigned v_mesh::neighbor(unsigned i, unsigned k)
{
	set<unsigned>::iterator sit = nbrs_[i].begin();
	unsigned j = 0;
	for(; sit != nbrs_[i].end() && j != k; ++j)
		++sit;
	if(j==k) return *sit;
	return -1;
}
	
vector<unsigned> v_mesh::bnd_verts()
{
	vector<unsigned> vx, bx;
	for(unsigned i=0; i<verts_.size(); ++i)
	{
		unsigned n=0; /* index of the vertex in bnd_ */
		for(unsigned j=0; j<bnd_.n_boundaries(); ++j)
			for(unsigned k=0; k<bnd_.size(j); ++k)
			{
				if(verts_[i] == bnd_.vert(j, k))
				{
					vx.push_back(i); /* ith vertex of v_mesh */
					bx.push_back(n); /* nth vertex of bnd_ */
				}
				n++;
			}
	}
	vector<unsigned> v(vx.size());
	for(unsigned j=0; j<vx.size(); ++j)
		v[bx[j]] = vx[j];
	
	return v;
}

void v_mesh::write_triangles(ofstream& os)
{
	
	os << "#VRML V2.0 utf8\n\n";
	os << "Shape {\n";
	os << " geometry IndexedFaceSet\n";
	os << "  { \n";
	os << "   coord Coordinate{\n";
	os << "    point [ \n";
	unsigned n = verts_.size();
	for(unsigned i=0; i<n; ++i) 
	{
		mesh_vertex v = verts_[i];
		os << v.x_ << ' ' << v.y_ << ' ' << 0.0 << ",\n";
	}
	os << "    ]} \n";
	os << "  coordIndex [ \n";
	unsigned NNE = tri_.size();
	for(unsigned i=0; i<NNE; ++i)
	{
		vtriangle t = tri_[i];
		os << t.vert(0) << ", "
		   << t.vert(1) << ", "
		   << t.vert(2) << ", " << -1 << '\n';
	}
	os << "  ]\n";
	os << "  colorPerVertex FALSE\n";
	os << "  color Color {\n";
	os << "   color [\n";
	
	for(unsigned i=1; i<=NNE; ++i)
		os << (0.25 + (float)i/NNE) <<' '<< (0.1 + (float)(NNE - i)/NNE) <<' '<< 0.6 <<",\n";

	os << "   ]\n";
	os << "  }\n";
	os << "  colorIndex [ \n";
	for(unsigned i=1; i<=NNE; ++i)
		os << i << ' ';
	os << "  ]\n";
	os << " } \n";         
	os << "} \n";
}
void v_mesh::write_sphere(ofstream& str, double x0, double y0, double z0, double rad,
						  const float r, const float g, const float b, const float t)
{
	str << "Transform {\n"
	<< "translation " << x0 << ' ' << y0 << ' ' << z0 << '\n'
	<< "childeren [\n"
	<< " Shape {\n"
	<< " appearance Appearance{\n"
	<< "   material Material\n"
	<< "    {\n"
	<< "	diffuseColor " << r << ' ' << g << ' ' << b << '\n'
	<< "	transparency " << t << '\n'
	<< "	}\n"
	<< "  }\n"
	<< " geometry Sphere\n"
	<< "  {\n"
	<< "  radius " << rad << '\n'
	<< "	}\n"
	<< "   }\n"
	<< "  ]\n"
	<< " }\n";
}

void v_mesh::write_verts(ofstream& os)
{
	unsigned nv = verts_.size();
	
	for(unsigned i=0; i<nv; ++i)
	{
		os << "#VRML V2.0 utf8\n\n";
		os << "Group {\n"
		   << "childeren [";
		mesh_vertex v = verts_[i];
		float r = (0.1 + 3*(float)i/nv), g = (0.1 + (nv - 3*(float)i)/nv);
		write_sphere(os, v.x_, v.y_, 0.0, 0.05, r, g, 0.4, 0.0);
	}
	os << "  ]\n"
	   << "}\n";
}
	
	
	
	

