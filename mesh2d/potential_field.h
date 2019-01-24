////////////////////////////////////////
// Ata Karakci
// December 15, 2010
// ENGN2912B
//
// This file contains the potential_field class
////////////////////////////////////////

#ifndef potential_field_h_
#define potential_field_h_

#include <vector>
#include <mesh2d/v_mesh.h>
#include <math.h>

class potential_field
{
public:

	potential_field(){}
	
	~potential_field(){}

	/* sets reciprocal distances */
	void set_rdist(v_mesh vm)
	{
		unsigned n = vm.n_verts();
		field_.resize(n); rdist_.resize(n); rdist_sum_.resize(n);
		
		for(unsigned i=0; i<n; ++i)
		{
			double xi = vm.vert(i).x_, yi = vm.vert(i).y_;
			unsigned m = vm.n_neighbors(i);
			rdist_[i].resize(m);
			for(unsigned j=0; j<m; ++j)
			{
				unsigned nj = vm.neighbor(i, j);
				double xj = vm.vert(nj).x_, yj = vm.vert(nj).y_;
				rdist_[i][j] = 1/sqrt(((xj-xi)*(xj-xi)) + ((yj-yi)*(yj-yi)));
				rdist_sum_[i] += rdist_[i][j];
			}
		}
	}
	
	void set_node(unsigned i, double const& phi)
	{
		field_[i] = phi;
	}

	double node(unsigned index){return field_[index];}
	double rdist(unsigned node, unsigned neigh){return rdist_[node][neigh];}
	double rdsum(unsigned node){return rdist_sum_[node];}
		
private:
	
	vector<double> field_;
	vector<vector<double> > rdist_;
	vector<double> rdist_sum_;
};



#endif //potential_field_h_