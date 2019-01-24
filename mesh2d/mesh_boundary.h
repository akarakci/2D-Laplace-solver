////////////////////////////////////////
// Ata Karakci
// December 5, 2010
// ENGN2912B
//
// This file contains the mesh_boundary class and global I/O stream operators
////////////////////////////////////////

#ifndef mesh_boundary_h_
#define mesh_boundary_h_

#include <mesh2d/mesh_vertex.h>
#include <sstream>
#include <vector>
#include <fstream>

using namespace std;

class mesh_boundary
{
public:

	mesh_boundary(){}
	~mesh_boundary(){}
	
	unsigned n_boundaries() const {return boundaries_.size();}
	unsigned size(unsigned k) const {return (boundaries_[k]).size();}
	mesh_vertex vert(unsigned i, unsigned j) const {return boundaries_[i][j];}
	double density(unsigned i, unsigned j) const {return linear_density_[i][j];}
	
	friend inline ostream& operator << (ofstream&, mesh_boundary const&);
	friend inline istream& operator >> (ifstream&, mesh_boundary&);
	
	vector<vector<mesh_vertex> > boundaries_;
	vector<vector<double> > linear_density_;
};

inline ostream& operator << (ofstream& os, mesh_boundary const& mb)
{
	vector<vector<mesh_vertex> > bd = mb.boundaries_;
	vector<vector<double> > ld = mb.linear_density_;
	
	/* number of boundaries */
	unsigned n = bd.size();
	
	/* number of vertices on each boundary */
	vector<unsigned> v(n);
	for(unsigned i=0; i<n; ++i)
		v[i] = bd[i].size();
	
	if(os)
	{
		os << 0.4 <<" "<< 0.4 << '\n' << n << '\n';
		
		unsigned w=0;
		for(unsigned i=0; i<n; ++i){
			os << w+1 <<" "<< v[i]+w << '\n';
			w += v[i];}
		
		os << w << '\n';
		
		for(unsigned i=0; i<n; ++i)
		{
			for(unsigned j=0; j<v[i]; ++j)
				os << bd[i][j].x_ <<" "<< bd[i][j].y_ <<" "<< ld[i][j] << '\n';
		}
	}
	
	os.close();
	
	return os;
}

inline istream& operator >> (ifstream& is, mesh_boundary& mb)
{
	double x, y, d; /* coordinates of the mesh_vertex and density */
	unsigned v, bd; /* num. of vertices in each boundary and num. of boundaries */

	vector<unsigned> nv; /* storage for v */
	vector<double> ld;   /* storage for d */
	vector<mesh_vertex> vx; /* storage for vt */
	
	while(is)
	{
		string line;
		getline(is, line);
		stringstream sst;
		sst.str(line);
		size_t bindx, vindx;
		
		string temp;
		
		bindx = line.find("Num_boundaries");
		vindx = line.find("Num_vertices");
		
		if(bindx!=string::npos)
			sst >> temp >> bd;
		else if(vindx!=string::npos)
		{
			sst >> temp >> v;
			nv.push_back(v);
		}
		else
		{
			sst >> x >> y >> d;
			mesh_vertex vt(x, y);
			vx.push_back(vt);
			ld.push_back(d);
		}
	}
		
	unsigned k = 0; /* counts the number of elements in vx and ld */
	for (unsigned i=0; i<bd; ++i)
	{
		vector<mesh_vertex> vert; /* temp. storage for vertices of each boundary */
		vector<double> dens;  /* temp. storage for linear densities at each vertex of boundary */
		unsigned j=0;
		while(j<nv[i])
		{
			vert.push_back(vx[k]);
			dens.push_back(ld[k]);
			++j;
			++k;
		}
		(mb.boundaries_).push_back(vert);
		(mb.linear_density_).push_back(dens);
		
	}
	return is;
}

#endif //mesh_boundary_h_