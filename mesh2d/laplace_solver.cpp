////////////////////////////////////////
// Ata Karakci
// December 15, 2010
// ENGN2912B
//
// This file contains the laplace_solver class
////////////////////////////////////////

#include <mesh2d/laplace_solver.h>
#include <mesh2d/karamete_mesh_generator.h>
#include <algorithm>
#include <math.h>
#include <vnl/vnl_vector.h>


void laplace_solver::set_boundary_conditions(vector<double>& bc)
{
	/* boundary vertices of v_mesh */
	vector<unsigned> bv = vm_.bnd_verts();
	
	/* set field values of the boundary rest being default 0 */
	for(unsigned i=0; i<bv.size(); ++i)
		field_.set_node(bv[i], bc[i]);
}
	
bool laplace_solver::solve(mesh_boundary const& mb,
						   vector<double>& bc, double alpha, double beta)
{
	/* generate mesh from mb */
	karamete_mesh_generator kmg;
	kmg.set_boundary(mb);
	kmg.generate(alpha, beta);
	vm_ = kmg.mesh();
	
	/* set members of field_ */
	field_.set_rdist(vm_);

	set_boundary_conditions(bc);
	
	/* vector containing the indices of boundary vertices */
	vector<unsigned> bv = vm_.bnd_verts();
	
	/* total number of vertices */
	unsigned nv = vm_.n_verts();
	
	double ftol = 1e-11;
	double maxd;
	do{
		maxd = -vnl_numeric_traits<double>::maxval;
		
		/* vector of potential differences at each node at each iteration */
		vector<double> diff(nv);
		for(unsigned i=0; i<nv; ++i)
		{
			double phi = field_.node(i);
			
			/* if the node is not a boundary vertex update the field value */
			vector<unsigned>::iterator vit = find(bv.begin(), bv.end(), i);
			if(vit == bv.end())
				field_.set_node(i, update_phi(i));
			phi -= field_.node(i);
			diff[i] = fabs(phi);
		}
		
		/* find the maximum change */
		for(unsigned i=0; i<nv; ++i)
			if(diff[i] > maxd) maxd = diff[i];
	}while(ftol < maxd);

	return true;	
}

double laplace_solver::update_phi(unsigned i)
{
	double potential_sum = 0;
	for(unsigned k=0; k<vm_.n_neighbors(i); ++k)
	{
		unsigned nk = vm_.neighbor(i, k);
		double d = field_.rdist(i, k);
		double phi_k = field_.node(nk);
		potential_sum += phi_k*d;
	}
	return potential_sum/field_.rdsum(i);
}

void laplace_solver::face_color(vtriangle vt, double& r, double& g, double& b)
{
	double phi_min = vnl_numeric_traits<double>::maxval;
	double phi_max = -vnl_numeric_traits<double>::maxval;
	for(unsigned i=0; i< vm_.n_verts(); ++i)
	{
		if(field_.node(i)<phi_min) phi_min = field_.node(i);
		if(field_.node(i)>phi_max) phi_max = field_.node(i);
	}
	
	/* find the average of potential values at vertices of the triangle */
	unsigned v0=vt.vert(0), v1=vt.vert(1), v2=vt.vert(2);
	double phi = (field_.node(v0)+field_.node(v1)+field_.node(v2))/3;
	
	/* scale the potential */
	phi -= phi_min;
	phi /= (phi_max - phi_min);
	
	/* linearly assign a color to each potential value (0 = b < g < r = 1) */
	if(phi<=0.25){b=1; g=2*phi; r=0;}
	if(phi>0.25 && phi<=0.75){b=2-4*phi; g=2*phi; r=0;}
	if(phi>0.5 && phi<=0.75){b=0; g=2-2*phi; r=4*phi-2;}
	if(phi>0.75){b=0; g=2-2*phi; r=1;}
 }

void laplace_solver::vertex_color(unsigned vrt, double& r, double& g, double& b)
{
	double phi_min = vnl_numeric_traits<double>::maxval;
	double phi_max = -vnl_numeric_traits<double>::maxval;
	for(unsigned i=0; i< vm_.n_verts(); ++i)
	{
		if(field_.node(i)<phi_min) phi_min = field_.node(i);
		if(field_.node(i)>phi_max) phi_max = field_.node(i);
	}
	
	double phi = field_.node(vrt);
	phi -= phi_min;
	phi /= (phi_max - phi_min);
	
	if(phi<=0.25){b=1; g=2*phi; r=0;}
	if(phi>0.25 && phi<=0.75){b=2-4*phi; g=2*phi; r=0;}
	if(phi>0.5 && phi<=0.75){b=0; g=2-2*phi; r=4*phi-2;}
	if(phi>0.75){b=0; g=2-2*phi; r=1;}
}

void laplace_solver::write_potential(ofstream& os)
{
	os << "#VRML V2.0 utf8\n\n";
	os << "Shape {\n";
	os << " geometry IndexedFaceSet\n";
	os << "  { \n";
	os << "   coord Coordinate{\n";
	os << "    point [ \n";
	unsigned n = vm_.n_verts();
	for(unsigned i=0; i<n; ++i) 
	{
		mesh_vertex v = vm_.vert(i);
		os << v.x_ << ' ' << v.y_ << ' ' << 0.0 << ",\n";
	}
	os << "    ]} \n";
	os << "  coordIndex [ \n";
	unsigned NNE = vm_.n_tri();
	for(unsigned i=0; i<NNE; ++i)
	{
		vtriangle t = vm_.tri(i);
		os << t.vert(0) << ", "
		<< t.vert(1) << ", "
		<< t.vert(2) << ", " << -1 << '\n';
	}
	os << "  ]\n";
	os << "  colorPerVertex FALSE\n";
	os << "  color Color {\n";
	os << "   color [\n";
	
	for(unsigned i=0; i<NNE; ++i)
	{
		/* set face color */
		double r=0., g=0., b=0.;
		face_color(vm_.tri(i), r, g, b);
		os << r <<' '<< g <<' '<< b <<",\n";
	}
	
	os << "   ]\n";
	os << "  }\n";
	os << "  colorIndex [\n";
	for(unsigned i=0; i<NNE; ++i)
		os << i << ' ';
	os << "]\n";
	os << "}\n";         
	os << "}\n";
}

void laplace_solver::write_nodes(ofstream& os)
{
	os << "#VRML V2.0 utf8\n\n";
	os << "Shape {\n";
	os << " geometry PointSet\n";
	os << "  { \n";
	os << "   coord Coordinate{\n";
	os << "    point [ \n";
	unsigned n = vm_.n_verts();
	for(unsigned i=0; i<n; ++i) 
	{
		mesh_vertex v = vm_.vert(i);
		os << v.x_ << ' ' << v.y_ << ' ' << 0.0 << ",\n";
	}
	os << "    ]} \n";
	os << "  color Color {\n";
	os << "      color [\n";
	
	for(unsigned i=0; i<n; ++i)
	{
		/* set vertex color */
		double r=0., g=0., b=0.;
		vertex_color(i, r, g, b);
		os << r <<' '<< g <<' '<< b <<",\n";
	}
	os << "   ]\n";
	os << "  }\n";
	os << "}\n";         
	os << "}\n";
}
