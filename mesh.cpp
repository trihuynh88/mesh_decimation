#include "mesh.h"
//#include "MemLeakCheck.h"

Quadric::Quadric()
{
	a2 = ab = ac = ad = b2 = bc = bd = c2 = cd = d2 = area = 0;
	memset(matrix,0,9*sizeof(double));
}

Quadric::Quadric(double a, double b, double c, double v1, double v2, double v3, double fArea)
{
	double d = -a*v1 - b*v2 - c*v3;
	matrix[0][0] = a2 = fArea*a*a;
	matrix[0][1] = matrix[1][0] = ab = fArea*a*b;
	matrix[0][2] = matrix[2][0] = ac = fArea*a*c;
	ad = fArea*a*d;
	matrix[1][1] = b2 = fArea*b*b;
	matrix[1][2] = matrix[2][1] = bc = fArea*b*c;
	bd = fArea*b*d;
	matrix[2][2] = c2 = fArea*c*c;
	cd = fArea*c*d;
	d2 = fArea*d*d;
	area = fArea;
}

void Quadric::sum(Quadric q)
{
	a2+=q.a2;
	ab+=q.ab;
	ac+=q.ac;
	ad+=q.ad;

	b2+=q.b2;
	bc+=q.bc;
	bd+=q.bd;
	
	c2+=q.c2;
	cd+=q.cd;

	d2+=q.d2;

	matrix[0][0] = a2;
	matrix[0][1] = matrix[1][0] = ab;
	matrix[0][2] = matrix[2][0] = ac;	
	matrix[1][1] = b2;
	matrix[1][2] = matrix[2][1] = bc;	
	matrix[2][2] = c2;
}


double Quadric::evaluateCost(double x, double y, double z)
{
	return x*x*a2 + 2*x*y*ab + 2*x*z*ac + 2*x*ad
	+ y*y*b2 + 2*y*z*bc + 2*y*bd
	+ z*z*c2 + 2*z*cd
	+ d2;
}

Vertex::Vertex()
{	
	match = -1;
	neighborEdges.reserve(1000);
	neighborFaces.reserve(1000);
}

Vertex::~Vertex()
{
	neighborEdges.clear();
	neighborFaces.clear();
}
/*
Vertex& Vertex::operator=( const Vertex& v ) {
      x = other.x;
      c = other.c;
      s = other.s;
      return *this;
}
*/

Face::Face()
{
	isValid = true;
}

Edge::Edge()
{
	isValid = true;
}

Edge::Edge(int iV1, int iV2, double fCost)
{
	v1 = iV1;
	v2 = iV2;
	cost = fCost;
	isValid = true;
}