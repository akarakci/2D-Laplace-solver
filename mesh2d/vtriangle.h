////////////////////////////////////////
// Ata Karakci
// December 12, 2010
// ENGN2912B
//
// This file contains the vtriangle class
////////////////////////////////////////

#ifndef vtriangle_h_
#define vtriangle_h_

class vtriangle
{
public:

	vtriangle()
	{
		for(unsigned i=0; i<3; ++i) verts_[i]=0;
	}
	vtriangle(unsigned v0, unsigned v1, unsigned v2)
	{
		verts_[0]=v0; verts_[1]=v1; verts_[2]=v2;
	}
	
	~vtriangle(){}
	
	void set_vert(unsigned i, unsigned v){verts_[i]=v;}
	unsigned vert(unsigned i) const {return verts_[i];}
	
private:
	
	unsigned verts_[3];
};



#endif //vtriangle_h_