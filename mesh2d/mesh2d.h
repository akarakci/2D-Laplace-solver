//make extern for C++
#ifdef __cplusplus //since extern doesn't exist for C compiler
extern "C" {
#endif

#ifndef mesh2d_h_
#define mesh2d_h_

/* Defined symbols */
#define eq       ==
//#define or       ||
#define neq      !=
//#define and      &&
#define true     1
#define false    0
#define MAXPSIZE 20001 /* Size of point array */
#define MAXESIZE 50001 /* Size of triangle array */

/* Data structures */
typedef struct point_list {double x[3];int bnd,next;int m;} point;

typedef struct triangle_list 
{
  int n[4];
  int t[4],nt;
  point xc;
  double R;
}
triangle;

typedef struct boundary_index {int n1,n2;} boundary;
typedef struct integer_point {long x[3];} ipoint;

/* global variables */
int NP,NB,NBP,MNP,NE,NE1,nbnd;
double alpha,beta;

/* storage arrays for points and triangles */
point*    pc;
triangle* tc;

/* Function interfaces */
void Load_Triangle(int,int,int,int,int);
double Cross(point,point);
double Dot(point,point);
int Node_In(int,int);
void Add_Boundary_Node(int);
short InCircleTest(point,triangle);
void Initialize(void);
void Suffered(int,int);
void Construct_New_Triangles(int);
int Newplaces(int);
void Write2tecfile(void);
int Engine(int);
void Construct_New_Neighbors1(void);
void Construct_New_Neighbors2(void);
int Find_Missing_Edges(void);
void Set_Next_Fields(void);
int Node_Renumber(int*, int*, int*, int*, int*, int*);
void Initialize(void);
void Intro(int , char **);
int Readfile(const char* file);
void Remove_Unwanted_Triangles(void);
void Insert_Nodes(double, double);
point Vector(point,point);
triangle Compute_Circumcircle(triangle);
void print_point (point*);
void print_triangle(triangle*);
void save_vrml(const char*, int*, int*, int*, int, int);

/* New C function for Laplace Solver */
void SetBoundary(int, int*, int*, int, double*, double*, double*);
	
	
#endif // mesh2d_h_
	
#ifdef __cplusplus
}
#endif
