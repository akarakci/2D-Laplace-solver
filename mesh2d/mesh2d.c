/* 
   Static Encoding of Unconstrainted Delaunay Triangulation 
               Algorithm by Dr. B. Kaan Karamete 

  Mesh Generator:mesh2d.c
   Must be compiled like: cc mesh2d.c -lm -O -o mesh2d

  Input : mesh2d.dat
  Input Generator: model2d.c with modler2d.h
  Input Format: 
          alpha beta <- local and global adaptation parameters.
          NB         <- number of boundaries
          bindx[1].n1 bindx[1].n2 <- Each boundary start and end nodes
            ..        ..
          bindx[NB].n1 bindx[NB].n2
          NBP        <-  Total number of nodes
          pc[1].x[1]  pc[1].x[2]    pc[1].x[0] <- Each node's coords and
                                                  point spacing value. 
            ..          ..             .. 
            ..          ..             ..
          pc[NBP].x[1] pc[NBP].x[2]  pc[NBP].x[0]  

  Output: mesh2d.out

*/
#include "mesh2d.h"
#include <stdio.h>
#include <math.h>
/*
   if maximum number of nodes is exceeded or you have more memory 
   you can increase both MAXPSIZE and/or MAXESIZE 
*/

#define SUFFERED 50

#define sqr(x) (x)*(x)


int NBP,NT,NSE,NewE,nblocks;
int      suf[SUFFERED+1];
triangle ntc[3*SUFFERED+1];
boundary bindx[10];
int near[10],NNEARS;
int nadapt=8; /*default number of adaptation cycles*/
int smoothing=1; /*default flag for coordinate smoothing*/
point max,min,a,b;
int sufdel[SUFFERED];
int extra=0; /* extra nodes to be put in case of missing edges*/ 
char *LastUpdate="  Last Update: 4 Dec 1997";

/* JLM print point */
void print_point (point* p)
{
  double x = p->x[1];
  double y = p->x[2];
  printf("p[%d]( %f %f )",p->bnd, x, y); 
}
/* JLM print triangle */
void print_triangle(triangle* t)
{
int i,n, pci;
point p;
  printf("Triangle:\n");
  for(i = 1; i<=3; ++i){
    pci = t->n[i];
    p = pc[pci];
    print_point(&p);
    printf("\n");
  }
  printf("Neighbors\n");
  for(n = 1; n<=(t->nt); ++n)
    printf("%d ",t->t[n]);
  printf("\n");
}

/* Implementation of Functions */

void Load_Triangle(id,n1,n2,n3,nt)
int id,n1,n2,n3,nt;
{ 
  tc[id].n[1]=n1;
  tc[id].n[2]=n2;
  tc[id].n[3]=n3;
  tc[id].nt=nt;
}
 

point Vector(p,q)
point p,q;
{ int i;
  point pq;
 for(i=1;i<=2;++i) pq.x[i]=q.x[i]-p.x[i];
 return(pq);
}

ipoint Vectori(p,q)
point p,q;
{ int i;
  ipoint pp,qq; 
  ipoint pq;
  int j;
   for(j=1;j<=2;++j)
   {
    pp.x[j]=(long)(a.x[j]*p.x[j]+b.x[j]);
    qq.x[j]=(long)(a.x[j]*q.x[j]+b.x[j]);
   } 
 
 for(i=1;i<=2;++i) pq.x[i]=qq.x[i]-pp.x[i];
 return(pq);
}

double Cross(p1,p2)
point p1,p2;
{
  return((p1.x[1]*p2.x[2]-p1.x[2]*p2.x[1]));
} 

long Crossi(p1,p2)
ipoint p1,p2;
{ 
  return((p1.x[1]*p2.x[2]-p1.x[2]*p2.x[1]));
} 

double Dot(p1,p2)
point p1,p2;
{
 return((p1.x[1]*p2.x[1]+p1.x[2]*p2.x[2]));
}


int Node_In(node,prev)
int node,prev;
{
  double epsnodein;
/*  long epsnodein;*/
  int i,j,tt,k,search,z;
  int n1,n2,n3,n4;
  short notin,ok;
  

  epsnodein=-1e-8; 
/*  epsnodein=0; */
  search=prev;
  do
  { 
   notin=false;
   for(i=1;i<=3;++i)
 if(Cross(Vector(pc[tc[search].n[i]],pc[tc[search].n[(i%3)+1]]),
           Vector(pc[tc[search].n[i]],pc[node]))
    < epsnodein)
      { notin=true;
        n1=tc[search].n[i];n2=tc[search].n[(i%3)+1];
        ok=false; 
        for(j=1;j<=tc[search].nt;++j)
        { 
         tt=tc[search].t[j];
         for(k=1;k<=3;++k)
         { 
           n3=tc[tt].n[k]; n4=tc[tt].n[(k%3)+1];          
           if(((n1 eq n3) && (n2 eq n4)) ||  
              ((n1 eq n4) && (n2 eq n3)))    
            { 
              search=tt;
              ok=true;
              break;
            }
         }
        if(ok) break; 
        } 
       if(ok) break;
      }
    } while(notin);
 
  return(search);
} 
  
short InCircleTest(node,tri)
point node;
triangle tri;
{ double epsradius=-1e-5;
  point point2R ;
  point2R=Vector(node,tri.xc);
  if(Dot(point2R,point2R)<(1.0+epsradius)*tri.R) return(true);
  else return(false);
}

triangle Compute_Circumcircle(tri)
triangle tri;
{
int i;
point pq,pr,pxc,xc;
double Area, f[4];
triangle triupdate;

pq=Vector(pc[tri.n[1]],pc[tri.n[2]]);
pr=Vector(pc[tri.n[1]],pc[tri.n[3]]);
Area=Cross(pq,pr);
for(i=1;i<=3;++i) f[i]=Dot(pc[tri.n[i]],pc[tri.n[i]]);
xc.x[1]=((f[2]-f[1])*pr.x[2]-(f[3]-f[1])*pq.x[2])/(2.0*Area);
xc.x[2]=((f[3]-f[1])*pq.x[1]-(f[2]-f[1])*pr.x[1])/(2.0*Area);
triupdate=tri;
triupdate.xc=xc;
pxc=Vector(pc[tri.n[1]],xc);
triupdate.R=Dot(pxc,pxc);
return(triupdate);
} 

Test_Cavity(node)
int node;
{
  int i,k,n1,n2,j,l;
  short olmadi;
  long Area;
 for(i=1;i<=NSE;++i) sufdel[i]=1;
 for(i=1;i<=NSE;++i)
  for(k=1;k<=3;++k)
   {
     n1=tc[suf[i]].n[k];
     n2=tc[suf[i]].n[(k%3)+1];
     olmadi=false;
    for(j=1;j<=NSE;++j)
    if(i neq j)
    {
      for(l=1;l<=3;++l)
      if(( (tc[suf[j]].n[l] eq n1) &&
           (tc[suf[j]].n[(l%3)+1] eq n2)
          ) ||
         ( (tc[suf[j]].n[l] eq n2) &&
           (tc[suf[j]].n[(l%3)+1] eq n1) ))
       {olmadi=true;break;}
 
     }
      if(!olmadi)
      { ntc[1].n[1]=n1;
        ntc[1].n[2]=n2;
        ntc[1].n[3]=node;
      
        Area=Crossi(Vectori(pc[ntc[1].n[1]],pc[ntc[1].n[2]]),
                    Vectori(pc[ntc[1].n[1]],pc[ntc[1].n[3]]));
        if(Area<0) sufdel[i]=-1;
      }
   }
}
   
/* Suffered finds a list of triangles, suf, whose circumcircle contains 
   node and can be found as transitive closure as neighbors of pivot */
void Suffered(node,pivot)
int node,pivot;
{ int i;
  short olmadi;
  int ind,j;
/*  int update[SUFFERED],nnse;  */
 NSE=1;
 suf[NSE]=pivot;  /*int suf[SUFFERED+1];*/
 ind=0;
while(ind<NSE && NSE<SUFFERED)
{++ind;
 pivot=suf[ind];/*next triangle in suf*/
 for(i=1;i<=tc[pivot].nt;++i)/*search over neighbors of pivot*/
 { olmadi=false; 
 for(j=1;j<=NSE;++j) /*interate over suf if neighbor is in suf, break*/
   if(tc[pivot].t[i] eq suf[j]) {olmadi=true;break;}
  if(!olmadi)
    /* if node is inside the circumcircle of the neighbor then add it to suf*/
    if(InCircleTest(pc[node],tc[tc[pivot].t[i]])) 
     if(tc[tc[pivot].t[i]].t[0] neq 1) /*not one of the bounding box verts*/
     { ++NSE;
       suf[NSE]=tc[pivot].t[i];
      }
  }
}
/*
 Test_Cavity(node);
 nnse=0;
 for(i=1;i<=NSE;++i)
  if(sufdel[i] neq -1) { ++nnse; update[nnse]=suf[i];}
 NSE=nnse;
 for(i=1;i<=NSE;++i) suf[i]=update[i];
*/
 printf("Leaving Suffered \n");
 for(i = 1; i<=NSE; ++i)
   print_triangle(&tc[suf[i]]);
 printf("==End suf triangles== \n");
}

/* Iterate through the set of triangles, suf and find a triangle
   that has an edge not common with another triangle in suf */
void Construct_New_Triangles(Newnode)
int Newnode;
{
 int i,j,k,l;
 short olmadi;
 int n1,n2;
 printf("In Construct_New_Triangles\n");

/*triangle ntc[3*SUFFERED+1]*/
 for(i=1;i<=3*SUFFERED-1;++i) ntc[i].nt=0;
 NewE=0;

 for(i=1;i<=NSE;++i) /* NSE is the size of suf */
  for(k=1;k<=3;++k)
   {
     n1=tc[suf[i]].n[k]; /* vertex k of suf[i] */
     n2=tc[suf[i]].n[(k%3)+1]; /* next ccw vertex k+1 */
     olmadi=false; /* fail = false */
   for(j=1;j<=NSE;++j) /* iterate over suf again */
    if(i neq j) /* skip i=j */
    { 
      for(l=1;l<=3;++l) /* vertex l of suf[j] */
      if(( (tc[suf[j]].n[l] eq n1) &&
           (tc[suf[j]].n[(l%3)+1] eq n2)
          ) ||
         ( (tc[suf[j]].n[l] eq n2) &&
           (tc[suf[j]].n[(l%3)+1] eq n1) ))
       {olmadi=true;
       break;}/* break if an edge of suf[i] is and edge of suf[j]*/
    }    

   /* found a triangle that doesn't share a common edge
      with a triangle in suf*/
	printf("found a triangle without common edge\n");
	print_triangle(&tc[suf[i]]);
      if(!olmadi) { ++NewE;
                    ntc[NewE].n[1]=n1; 
                    ntc[NewE].n[2]=n2; 
                    ntc[NewE].n[3]=Newnode; 
                    ntc[NewE].n[0]=suf[i];
                    ntc[NewE]=Compute_Circumcircle(ntc[NewE]);
                    
      } 
   }
}  

int Newplaces(newel)
int newel;
{
 if(newel<=NSE) return(suf[newel]);
 else return(NE+(newel-NSE));
}
/* construct neighborhood info for new triangles
  for new triangles 
  if any two nodes of the new triangle=the other new triangles
  or
  the neighborhood triangles of its suffered triangle
  stored at ntc[i].n[0].
*/

void Construct_New_Neighbors1(void)
{
 int i,j,k,l,n1,n2,kk;
 short oldu;
 int neybor;

/* check for neighborhood triangles of its suffered triangle*/

 for(i=1;i<=NewE;++i)
  for(k=1;k<=3;++k)
  {
    n1=ntc[i].n[k];
    n2=ntc[i].n[(k%3)+1];
    for (j=1;j<=tc[ntc[i].n[0]].nt;++j)
    { oldu=false;
      neybor=tc[ntc[i].n[0]].t[j];
      for(l=1;l<=3;++l)
     if( ((tc[neybor].n[l] eq n1) && (tc[neybor].n[(l%3)+1] eq n2))
          ||
        ((tc[neybor].n[l] eq n2) && (tc[neybor].n[(l%3)+1] eq n1)))
     { oldu=true; break;}
   if(oldu) break;
  }
  if(oldu)
  {
    ++ntc[i].nt;
    ntc[i].t[ntc[i].nt]=neybor;
    for(kk=1;kk<=tc[neybor].nt;++kk)
     if(tc[neybor].t[kk] eq ntc[i].n[0]) break;
    tc[neybor].t[kk]=Newplaces(i);
    
  }   
 }         

}

void Construct_New_Neighbors2(void)
{
int i,j,k,l,n1,n2;
short oldu;
int neybor;
 oldu=false;

 for(i=1;i<=NewE;++i)
  for(k=1;k<=3;++k)
  {
    n1=ntc[i].n[k];
    n2=ntc[i].n[(k%3)+1];
    for(j=1;j<=NewE;++j)
    if(i neq j)
    { oldu=false;
      for(l=1;l<=3;++l)
      if( ((ntc[j].n[l] eq n1) && (ntc[j].n[(l%3)+1] eq n2))
          ||
          ((ntc[j].n[l] eq n2) && (ntc[j].n[(l%3)+1] eq n1)))
      { oldu=true; break;}
     if (oldu) break;
    }
    if(oldu)
    {
    ++ntc[i].nt;
    ntc[i].t[ntc[i].nt]=Newplaces(j);
    }
  }

 for(i=1;i<=NewE;++i){
   tc[Newplaces(i)]=ntc[i];}
 NE=NE+NewE-NSE;
}

Engine(Newnode)
int Newnode;
{
 static int pivot=1;
 int i;
 /* jlm debug */
 printf("In Engine\n");
 printf("NewNode ");
 print_point(&pc[Newnode]);

 pivot=Node_In(Newnode,pivot);


 /* jlm debug */
 printf("pivot ");
 print_triangle(&tc[pivot]);

 Suffered(Newnode,pivot);

 /* jlm debug */
 printf("pivot after suffered ");
 print_triangle(&tc[pivot]);


 Construct_New_Triangles(Newnode);
 Construct_New_Neighbors1();
 Construct_New_Neighbors2();
}
void Set_Next_Fields(void)
{ int i,j;
 for(i=1;i<=NBP;++i) {pc[i].next=0;pc[i].bnd=0;}
 for(i=1;i<=NB;++i)
 {
  for(j=bindx[i].n1;j<=bindx[i].n2;++j)
  {
   pc[j].next=j+1;
   pc[j].bnd=i;
   }
   pc[bindx[i].n2].next=bindx[i].n1;
 }
 for(i=1;i<=NBP;++i) pc[i].m=pc[i].next;
}

int Readfile(const char* file)
{ 
  FILE *in;
  int i,j;
  
  in=fopen(file,"r");
  if(in eq NULL) return(0);
    
  fscanf(in,"%lf %lf\n",&alpha,&beta);
  fscanf(in,"%d\n",&NB);
  for(i=1;i<=NB;++i)
  {
   fscanf(in,"%d %d\n",&bindx[i].n1,&bindx[i].n2);
   bindx[i].n1=bindx[i].n1+4;
   bindx[i].n2=bindx[i].n2+4;
  }
  fscanf(in,"%d\n",&NT); /* NP olacak */
  NBP=NT+4;
  fscanf(in,"%lf %lf %lf\n",&pc[5].x[1],&pc[5].x[2],&pc[5].x[0]);
   min=pc[5];
   max=pc[5];
  for(i=6;i<=NBP;++i)
  {
   fscanf(in,"%lf %lf %lf\n",&pc[i].x[1],&pc[i].x[2],&pc[i].x[0]);
    for(j=1;j<=2;++j)
    {
     if(pc[i].x[j]<min.x[j]) min.x[j]=pc[i].x[j]; 
     if(pc[i].x[j]>max.x[j]) max.x[j]=pc[i].x[j]; 
    }
   } 
  fclose(in);
  nbnd=bindx[NB].n2;
  Set_Next_Fields();
  return(1);
}

void Initialize()
{
  int i;
  point boxmax,boxmin;

 /* create initial convex hull */

 for(i=1;i<=2;++i)
 {
  pc[1].x[i]=min.x[i]-fabs((max.x[i]-min.x[i])*0.2); 
  pc[3].x[i]=max.x[i]+fabs((max.x[i]-min.x[i])*0.2);
 }
 pc[2].x[1]=pc[3].x[1];
 pc[2].x[2]=pc[1].x[2];
 pc[4].x[1]=pc[1].x[1];
 pc[4].x[2]=pc[3].x[2];

 Load_Triangle(1,1,2,3,1);
 tc[1].t[1]=2;
 Load_Triangle(2,1,3,4,1);
 tc[2].t[1]=1;
 NP=4; 
 tc[1]=Compute_Circumcircle(tc[1]);
 tc[2]=Compute_Circumcircle(tc[2]);
 NE=2;
 /* JLM */
 print_triangle(&tc[1]);
 print_triangle(&tc[2]);
 
/*
 Points 1-4 and Triangles 1-2 are used for creating  the initial template 
*/
 for(i=1;i<=2;++i)
 { 
  boxmin.x[i]=1.0;
  boxmax.x[i]=10000.0;
  a.x[i]=(boxmax.x[i]-boxmin.x[i])/(pc[3].x[i]-pc[1].x[i]);
  b.x[i]=boxmax.x[i]-a.x[i]*pc[3].x[i];
 }
puts("**********************************************************");
printf("\nGenerator in progress starting with %d nodes and\n",NBP-4);
printf(" with the parameters n1=%d  n3=%d\n\n",
                  nadapt,smoothing);
puts("**********************************************************");

}

void Write2tecfile(void)
{
 FILE *t,*f;
 int i,j,trino;

 t=fopen("mesh2d.out","w");
 fprintf(t,"VARIABLES=x y dpi bnd\n");
 fprintf(t,"ZONE N=%d,E=%d,F=FEPOINT,ET=TRIANGLE\n",NP,NE-NE1);
 for(i=1;i<=NP;++i) 
  fprintf(t,"%f %f %f %d\n",
pc[i].x[1],pc[i].x[2],pc[i].x[0],pc[i].bnd);
 for(i=1;i<=NE;++i)
  if(tc[i].t[0] neq 1)
  fprintf(t,"%d %d %d\n",tc[i].n[1],tc[i].n[2],tc[i].n[3]);
 fclose(t);
 printf("Output file mesh2d.out is generated!\n");

 f=fopen("ptc.dat","w");
 fprintf(f,"%d\n",NE);
 for(i=1;i<=NE;++i)
  { fprintf(f,"%d ",i);
     for(j=1;j<=tc[i].nt;++j)
     { trino=tc[i].t[j];
       if(tc[trino].t[0] neq 1)
        fprintf(f,"%d ",trino);
       else fprintf(f,"-1 ");
      }
    fprintf(f,"\n");
  }
 fclose(f);
 printf("Triangles info file ptc.dat is generated!\n");

}  

void Remove_Unwanted_Triangles(void)
{ int i,j;
  int pivot;
  int temp;

 puts("Unwanted triangles are being swept off...");
 for(i=1;i<=NE;++i)
  for(j=1;j<=3;++j)
   if(tc[i].n[j] <= 4) {tc[i].t[0]=1;++NE1;break;}

 for(i=1;i<=NE;++i)
 if (
     pc[tc[i].n[2]].bnd eq pc[tc[i].n[1]].bnd &&
     pc[tc[i].n[3]].bnd eq pc[tc[i].n[1]].bnd && 
     pc[tc[i].n[1]].bnd neq 0)
 { 
   j=1; 
   pivot=tc[i].n[1];
   do
   { 
    if(pc[pivot].next eq tc[i].n[(j%3)+1]) break;
     else if(pc[pivot].next eq tc[i].n[(((j%3)+1)%3)+1])
                        {tc[i].t[0]=1;++NE1;break;}
      else if(j<3) pivot=tc[i].n[++j];
       else pivot=pc[pivot].next;
   }while(true);
 }
} 

int Find_Missing_Edges(void)
{
 int i,j,k,pivot;
 int missing[100];
 int MNP=0;

 for(i=1;i<=NE;++i)
  for(j=1;j<=3;++j)
   for(k=1;k<=3;++k)
    if(pc[tc[i].n[j]].m eq tc[i].n[k])
        pc[tc[i].n[j]].m=0;
 
  for(i=1;i<=NP;++i)
   if(pc[i].m neq 0) missing[++MNP]=i; 

 for(i=1;i<=MNP;++i)
 {
  ++NBP;
  ++NP;
  pc[NP].next=pc[missing[i]].next;
  pc[NP].m=pc[missing[i]].next;
  pc[NP].bnd=pc[missing[i]].bnd;
  pc[NP].x[1]=(pc[missing[i]].x[1]+pc[pc[missing[i]].next].x[1])/2.0;
  pc[NP].x[2]=(pc[missing[i]].x[2]+pc[pc[missing[i]].next].x[2])/2.0;
  pc[NP].x[0]=pc[missing[i]].x[0];
  printf("%lf %lf\n",pc[NP].x[1],pc[NP].x[2]);
  pc[missing[i]].next=NP;
  pc[missing[i]].m=NP;
  pc[missing[i]].bnd=pc[missing[i]].bnd;
  nbnd=bindx[NB].n2+1;
 }
 printf("Detecting and Curing %d Missing Edges...\n",MNP);
 extra=extra+MNP;
return(MNP);
}

void Near_Nodes(tno,nn)
int tno;
int *nn;
{ int i,j,k,l;
   int pivot[5];
  pivot[1]=tno;
  for(i=1;i<=tc[tno].nt;++i) pivot[i+1]=tc[tno].t[i]; 
  k=0;
 for(l=1;l<=tc[tno].nt+1;++l)
  for(i=1;i<=tc[pivot[l]].nt;++i)
    if(tc[tc[pivot[l]].t[i]].n[0] neq 0)
      {++k;near[k]=tc[tc[pivot[l]].t[i]].n[0];}
  *nn=k;
} 

void Insert_Nodes(alpha,beta)
double alpha,beta;
{
 int  i,j,k;
 point p[5];
 double dm[4],sj;
 short  reject;
 short EXCEEDED=false;
 int nold=0;
 int cnt=0;

 while(cnt<nadapt && (NP-nold)>(0.1*nold))
 { ++cnt;
  nold=NP;
  for(i=1;i<=NE;++i) tc[i].n[0]=0;
  for(i=1;i<=NE;++i)
  if(tc[i].t[0] eq 0)
  {
  reject=0;
  /* p[4] is centroid of the 3 vertices of tc[i] */
  for (j=1;j<=3;++j) p[j]=pc[tc[i].n[j]];
  p[4].x[0]=(p[1].x[0]+p[2].x[0]+p[3].x[0])/3.0;
  p[4].x[1]=(p[1].x[1]+p[2].x[1]+p[3].x[1])/3.0;
  p[4].x[2]=(p[1].x[2]+p[2].x[2]+p[3].x[2])/3.0;
  /* dm[j] is distance from the centroid to vertex j */
  for (j=1;j<=3;++j)
  dm[j]=sqrt((p[j].x[1]-p[4].x[1])*(p[j].x[1]-p[4].x[1])
           +(p[j].x[2]-p[4].x[2])*(p[j].x[2]-p[4].x[2]));
  /* if any distance is less that alpha*spacing then reject */
  for (j=1;j<=3;++j)  if (dm[j]<(alpha*p[4].x[0])) { reject=1; break;}

  /* if not rejected then Near_Nodes */
  if(!reject)
  {
  Near_Nodes(i,&NNEARS);
  for(j=1;j<=NNEARS;++j)
  {
   sj=sqrt(sqr(pc[near[j]].x[1]-p[4].x[1])+
        sqr(pc[near[j]].x[2]-p[4].x[2]));
     if (sj<(beta*p[4].x[0])) { reject=1;break;}
   }  
  } 
  if(!reject)
   { 
    ++NP;
    pc[NP]=p[4];
    pc[NP].next=0;
    pc[NP].bnd=0;
    tc[i].n[0]=NP;
    /* jlm -- debug print */
    printf("\nadding point %d",NP);
    print_point(&pc[NP]);
    if(NP eq MAXPSIZE-1) { EXCEEDED=true;break;}
    
   }
  }
 printf("%d Adaptation: %d nodes",cnt,NP-4);
 for(i=nold+1;i<=NP;++i) Engine(i);
 printf(" %d triangles after insert\n",NE-NE1);

 /* JLM -- debug print */
 for(i=NE1+1;i<=NE;++i)
   /*  if(tc[i].t[0] neq 1)*/
    print_triangle(&tc[i]);

 if(EXCEEDED) { printf("MAXPSIZE exceeded! Stopping...\n");break;}
 }
}

Compute_Centroids(tno,xc,yc)
int tno;
double *xc;
double *yc;
{
  int i;
  double sumx,sumy;
 sumx=0.0;sumy=0.0;
 for(i=1;i<=3;++i)
 {
  sumx=sumx+pc[tc[tno].n[i]].x[1];
  sumy=sumy+pc[tc[tno].n[i]].x[2];
 }
 *xc=sumx/3.0;
 *yc=sumy/3.0;
}
  
int Node_Renumber(tri, node, nnode, ntri, NNE_p, NNP_p)
  int* tri;
  int* node;
  int* nnode;
  int* ntri;
  int* NNE_p; int* NNP_p;
{
  /*int tri[MAXESIZE],node[MAXPSIZE],nnode[MAXPSIZE],ntri[MAXESIZE];*/
int NNP, NNE;
int i,j,pivot,cpivot,trino,cnode,ind,ok;
double xc,yc,xsc,ysc,minimum,diff;
int k,m,z,maxz,band[4],maxband,neks;
int piv,sakla,stack[500],nstack,maxinter;
point orta;
/*FILE *t;*/
puts("Node and Triangle Renumbering...");      
NNP=0;
NNE=1;
for(i=1;i<=NE;++i) tri[i]=0;
for(i=1;i<=NP;++i) node[i]=0;
if(true){   
  orta.x[1]=min.x[1]+(max.x[1]-min.x[1])/2.0;
   orta.x[2]=min.x[2]+(max.x[2]-min.x[2])/2.0;}

 minimum=sqr(max.x[1]);  
 for(i=1;i<=NE;++i)
   if(tc[i].t[0] neq 1)
   {
      Compute_Centroids(i,&xc,&yc);
       switch (3)
        { case 1: diff=sqr(xc-max.x[1]);
                  if(diff<minimum) {minimum=diff;pivot=i;}
                  break;
          case 2: diff=sqr(yc-max.x[2]);
                  if(diff<minimum) {minimum=diff;pivot=i;}
                  break;
          case 3:
                  diff=sqr(xc-orta.x[1])+sqr(yc-orta.x[2]);
                  if(diff<minimum) {minimum=diff;pivot=i;}
                  break; 
          default:
                  diff=sqr(xc-orta.x[1])+sqr(yc-orta.x[2]);
                  if(diff<minimum) {minimum=diff;pivot=i;}
                  break; 
        }
     }  
tri[pivot]=NNE;
for(i=1;i<=3;++i) {++NNP;node[tc[pivot].n[i]]=NNP;}
Compute_Centroids(pivot,&xc,&yc);
xsc=xc;ysc=yc;
nstack=0;
stack[++nstack]=pivot;
maxinter=0;
maxband=0;
do
{ ind=0;
  if(nstack>maxinter) maxinter=nstack;
  for(j=1;j<=nstack;++j)
  { pivot=stack[j];
    sakla=ind;
  for(i=1;i<=tc[pivot].nt;++i)
  { trino=tc[pivot].t[i];
    if((tri[trino] eq 0) && (tc[trino].t[0] neq 1))
    {
      Compute_Centroids(trino,&xc,&yc);
      ++ind;
      switch (3)
      { case 1: diff=sqr(xc-xsc);break;
        case 2: diff=sqr(yc-ysc);break;
        case 3: diff=sqr(xc-xsc)+sqr(yc-ysc);break;
        case 4: z=0;
                for(k=1;k<=3;++k)
                  for(m=1;m<=3;++m)
                    if(tc[pivot].n[k]==tc[trino].n[m])
                      band[++z]=tc[pivot].n[k];
                  for(m=1;m<=3;++m)
                  { ok=0;
                   for(k=1;k<=z;++k)
                   if(band[k] eq tc[trino].n[m]) {ok=1; break;}
                    if(!ok) {neks=tc[trino].n[m];break;}
                   }
               if(node[neks] eq 0) neks=NNP+1; else neks=node[neks];
                maxz=0;
                for(k=1;k<=z;++k) 
                  if(abs(node[band[k]]-neks)>maxz)
                      maxz=abs(node[band[k]]-neks);
                 diff=(sqr(xc-xsc)+sqr(yc-ysc))*(maxz);
                 break;
        default: diff=sqr(xc-xsc)+sqr(yc-ysc);break;
       }
      if(ind eq 1) { minimum=diff; cpivot=trino;piv=pivot;}
      else if (minimum>diff) {minimum=diff; cpivot=trino;piv=pivot;}
    }
  }
  if(sakla eq ind) stack[j]=0;
  }   
  for(i=1;i<=nstack;++i)
    if(stack[i] eq 0)
    { 
      for(j=i+1;j<=nstack;++j)
       stack[j-1]=stack[j];
      --nstack;
     }
    
  if(ind>0)
  {
    stack[++nstack]=cpivot;
    Compute_Centroids(cpivot,&xc,&yc);
    xsc=(xsc*NNE+xc)/(NNE+1);
    ysc=(ysc*NNE+yc)/(NNE+1);
    ++NNE;
    tri[cpivot]=NNE;
    for(i=1;i<=3;++i)
    {  
      ok=0;
      for(j=1;j<=3;++j)
      if(tc[cpivot].n[i] eq tc[piv].n[j])
      { ok=1;break;}
     if(!ok)
     {
      cnode=tc[cpivot].n[i];
      break;
     }
    }
    if(node[cnode] eq 0) {++NNP;node[cnode]=NNP;}

     for(k=1;k<=3;++k) 
      if(abs(node[tc[cpivot].n[k]]-node[tc[cpivot].n[(k%3)+1]])>maxband)
        maxband=abs(node[tc[cpivot].n[k]]-node[tc[cpivot].n[(k%3)+1]]);
   }    
  else
  break;
 }while(NNE<(NE-NE1));
printf("Maximum interface elements= %d\n",maxinter);
printf("Maximum band= %d\n",maxband);
if(NP-4 eq NNP) 
 printf("ok.\n"); 
 else {printf("Renumbering could not be implemented properly\n"); 
       return(0);}

for(i=1;i<=NP;++i) 
  nnode[node[i]]=i;
for(i=1;i<=NE;++i)
  ntri[tri[i]]=i;
/* output number of nodes and triangles */
 *NNE_p = NNE; *NNP_p = NNP;
return(1);
}
void save_vrml(fname, node, nnode, ntri, NNE, NNP)
     const char* fname;
     int* node;
     int* nnode;
     int* ntri;
     int NNE, NNP;
{
 int i;
 FILE *t;
 t=fopen(fname,"w"); 


 fprintf(t,"#VRML V2.0 utf8\n\n");
 fprintf(t,"Shape {\n");
 fprintf(t," geometry IndexedFaceSet\n");
 fprintf(t,"  { \n");
 fprintf(t,"   coord Coordinate{\n");
 fprintf(t,"    point [ \n");
 for(i=1;i<=NNP;++i) 
  fprintf(t,"%f %f %f,\n",
  pc[nnode[i]].x[1],pc[nnode[i]].x[2],0.0);
 fprintf(t,"    ]} \n");
 fprintf(t,"  coordIndex [ \n");
 
 for(i=1;i<=NNE;++i)
  fprintf(t,"%d, %d, %d, %d,\n",node[tc[ntri[i]].n[1]]-1,
 node[tc[ntri[i]].n[2]]-1,node[tc[ntri[i]].n[3]]-1, -1);
 fprintf(t,"  ]\n");
 fprintf(t,"  colorPerVertex FALSE\n");
 fprintf(t,"  color Color {\n");
 fprintf(t,"   color [\n");
 for(i=1;i<=NNE;++i){
   
	 fprintf(t, "%f %f %f,\n", 0.25 + (float)i/NNE, 0.1 + (float)(NNE - i)/NNE, 0.6);
 }
 fprintf(t,"   ]\n");
 fprintf(t,"  }\n");
 fprintf(t,"  colorIndex [ \n");
 for(i=1;i<=NNE;++i)
 fprintf(t,"%d  ", i);
 fprintf(t,"  ]\n");
 fprintf(t," } \n");         
 fprintf(t,"} \n"); 

 fclose(t);
 printf("Output file mesh2d.out is generated!\n");
}

void Intro(c,v)
int c;
char *v[];
{
  puts("\nA Delaunay 2D Constrained Unstructured Mesh");
  puts(" Generator developed by B. Kaan Karamete");
  puts(LastUpdate);
  puts(" Usage: mesh2d <n1> <n2> <n3>");
  puts("   n1: number of adaptation cycles<1..10>");
  puts("   n2: centroidal node renumbering style<1..3>");
  puts("   n3: Flag for coordinate Smoothing<0..1>");
  puts(" defaults: n1=8; n2=3; n3=1");    
if(c>2) nadapt=atoi(v[2]);  if(nadapt<0) nadapt=8;
/*if(c>3) bloking=atoi(v[3]); if(!(bloking>=0 and bloking<5)) bloking=3;*/
if(c>4) smoothing=atoi(v[4]); if(smoothing neq 0 && smoothing neq 1) 
 smoothing=1;

}    

/* not used in this project - for reference only */
Smooth()
{
  int ptc[MAXPSIZE][20];
  int i,j,tno,k;
  double xc,yc,xsc,ysc;

puts("Smoothing three times...");      
for(i=1;i<=NP;++i) ptc[i][0]=0;

 for(i=1;i<=NE;++i)
  if(tc[i].t[0] neq 1)
   for(j=1;j<=3;++j)
   if(tc[i].n[j]>bindx[NB].n2+4+extra)
   {++ptc[tc[i].n[j]][0];
    tno=ptc[tc[i].n[j]][0];
    ptc[tc[i].n[j]][tno]=i;
   }

for(k=1;k<=3;++k)
 for(i=1;i<=NP;++i)
 if(i>bindx[NB].n2+4+extra)
 { xsc=0.0;
   ysc=0.0;
  for(j=1;j<=ptc[i][0];++j)
  {
    Compute_Centroids(ptc[i][j],&xc,&yc);
    xsc=xsc+xc;
    ysc=ysc+yc;
  }
  if(ptc[i][0]>0)
  {
   pc[i].x[1]=xsc/ptc[i][0];
   pc[i].x[2]=ysc/ptc[i][0];
  }
 }

}
#if 0
main(int c,char *v[])
{
int i,nexact;
Intro(c,v);
if(!Readfile(v[1])){
    printf("There is no mesh2d.dat file present!!!\a\n"); return(0);}
Initialize();

do { ++NP;Engine(NP); } while(NP<NBP);

while((MNP=Find_Missing_Edges()) neq 0) 
 for(i=NP-MNP+1;i<=NP;++i) Engine(i);

for(i=1;i<=NE;++i) tc[i].t[0]=0;
NE1=0;
 printf("Just before Unwanted triangles\n"); 
 for(i=1;i<=NE;++i)
   print_triangle(&tc[i]);
Remove_Unwanted_Triangles();

 printf("=====Just after Unwanted triangles====\n"); 
#if 0
 for(i=1;i<=NE;++i)
 if (
     pc[tc[i].n[2]].bnd eq pc[tc[i].n[1]].bnd and
     pc[tc[i].n[3]].bnd eq pc[tc[i].n[1]].bnd and 
     pc[tc[i].n[1]].bnd neq 0)
   print_triangle(&tc[i]);
#endif

printf("Bnd. Triangulation: %d nodes %d triangles\n",NP-4,NE-NE1);
Insert_Nodes(alpha,beta);
#if 0
if(smoothing) Smooth();
#endif

if(bloking!=0) 
{
  if(!Node_Renumber()) Write2tecfile();
} else Write2tecfile();  
printf("With %d nodes %d triangles are generated!\n",NP-4,NE-NE1);
nexact=2*((NP-4)-(nbnd-4))+(nbnd-4)-2+2*(NB-1); 
if(nexact neq (NE-NE1)) 
printf("Something wrong! Should be=%d\n\a",nexact);
return 0;

}
#endif

/* New C function for Laplace Solver */

void SetBoundary(n_boundaries, boundary_start, boundary_end, nverts, Xb, Yb, Sb)
	int n_boundaries;
	int nverts;
	int* boundary_start;
	int* boundary_end;
	double* Xb; 
	double* Yb;
	double* Sb;
{	
	int i, j;
	NB = n_boundaries;
	for(i=1; i<=NB; ++i)
	{
		bindx[i].n1 = boundary_start[i-1]+1;
		bindx[i].n2 = boundary_end[i-1]+1;
		bindx[i].n1 = bindx[i].n1+4;
		bindx[i].n2 = bindx[i].n2+4;
	}
	NT = nverts;
	NBP=NT+4;
	
	point pcm[MAXPSIZE];
	pc = pcm;
	
	pc[5].x[1]=Xb[0]; pc[5].x[2]=Yb[0]; pc[5].x[0]=Sb[0];
			
	min=pc[5];
	max=pc[5];
	
	for(i=6; i<=NBP; ++i)
	{
		pc[i].x[1]=Xb[i-5]; pc[i].x[2]=Yb[i-5]; pc[i].x[0]=Sb[i-5];
		for(j=1; j<=2; ++j)
		{
			if(pc[i].x[j]<min.x[j]) min.x[j]=pc[i].x[j];
			if(pc[i].x[j]>max.x[j]) max.x[j]=pc[i].x[j];
		}
	}
	nbnd = bindx[NB].n2;
	Set_Next_Fields();
}
	
	

