////////////////////////////////////////
// Ata Karakci
// December 12, 2010
// ENGN2912B
//
// This file contains the karamete_mesh_generator class
////////////////////////////////////////


#include <mesh2d/karamete_mesh_generator.h>
#include <iostream>
#include <stdio.h>
#include <mesh2d/mesh2d.h> /* include as a normal c++ header */
using namespace std;

karamete_mesh_generator::karamete_mesh_generator(){}

void karamete_mesh_generator::set_boundary (mesh_boundary const& mb)
{
	int nb = mb.n_boundaries();
	int* boundary_start = new int[nb];
	int* boundary_end = new int[nb];
	unsigned v_index = 0;

	for(unsigned b=0; b<nb; ++b)
	{
		unsigned nv = mb.size(b);
		unsigned end = v_index + nv -1;
		boundary_start[b] = v_index;
		boundary_end[b] = end;
		v_index += nv;
	}

	double* Xb = new double[v_index];
	double* Yb = new double[v_index];
	double* Sb = new double[v_index];
	unsigned j = 0;
	for(unsigned b=0; b<nb; ++b)
		for(unsigned i=0; i<mb.size(b); ++i)
		{
			Xb[j] = mb.vert(b, i).x_;
			Yb[j] = mb.vert(b, i).y_;
			Sb[j] = mb.density(b, i);
			++j;
		}
	
	SetBoundary(nb, boundary_start, boundary_end, v_index, Xb, Yb, Sb);
	vm_.set_boundary(mb);

	delete boundary_start;
	delete boundary_end;
	delete Xb; delete Yb; delete Sb;
}

bool karamete_mesh_generator::generate(double alpha, double beta)
{
	int tri[MAXESIZE], node[MAXPSIZE], nnode[MAXPSIZE], ntri[MAXPSIZE];
	int i, nexact, NNE, NNP;
	
	triangle tcm[MAXESIZE];
	tc = tcm;
	
	Initialize();
	
	do { ++NP; Engine(NP);} while(NP<NBP);
	
	while((MNP=Find_Missing_Edges()) neq 0)
		for(i=NP-MNP+1; i<=NP; ++i) Engine(i);
	
	for(i=1; i<=NE; ++i) tc[i].t[0]=0;
	NE1=0;
	
	Remove_Unwanted_Triangles();
	
	Insert_Nodes(alpha, beta);
	
	if(!Node_Renumber(tri, node, nnode, ntri, &NNE, &NNP))
		return false;
	
	cout << NP-4 <<"nodes and "<< NE-NE1 <<"triangles are generated\n";
		
	/* set the size of v_mesh */
	vm_.resize(NNP, NNE);
	
	/* copy vertices into v_mesh */
	for(i=1; i<=NNP; ++i)
	{
		mesh_vertex v(pc[nnode[i]].x[1], pc[nnode[i]].x[2]);
		vm_.set_vertex(i-1, v);
	}
	
	/* copy the triangles into v_mesh */
	for(i=1; i<=NNE; ++i)
	{
		vtriangle t(node[tc[ntri[i]].n[1]]-1,
					node[tc[ntri[i]].n[2]]-1,
					node[tc[ntri[i]].n[3]]-1);
		vm_.set_triangle(i-1,t);
	}
	/* set the vertex neighborhood */
	this->set_neighbors();
	
	return true;
}

void karamete_mesh_generator::set_neighbors()
{
/* maps to the next and previous vertex on the triangle */
	unsigned forward[3] = {1, 2, 0};
	unsigned backward[3] = {2, 0, 1};

/* iterate through the triangles */
	for (unsigned i=0; i<vm_.n_tri(); ++i)
	{
		vtriangle& t = vm_.tri(i);
		
		for(unsigned j=0; j<3; ++j)
		{			
			vm_.add_neighbor(t.vert(j), t.vert(forward[j]));
			vm_.add_neighbor(t.vert(j), t.vert(backward[j]));
		}
	}
}


	





