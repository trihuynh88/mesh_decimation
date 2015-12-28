#include "vector"
using namespace std;

#pragma once


typedef struct _MeshSegment 
{ 
    std::vector<double> m_vCoord; 
    std::vector<int> m_vTriangleIndex; 
    std::vector<double> m_vTriangleNormal; 

} MeshSegment;


class Edge
{
public:
	int v1, v2;
	double cost;
	int posInHeap;
	bool isValid;	
	Edge();
	Edge(int iV1, int iV2, double fCost);
};

class Face
{
public:	
	int v1, v2, v3;	
	//normal
	double n1, n2, n3;
	bool isValid;
	Face();
};


class Quadric
{
public:
	double a2, ab, ac, ad;
	double b2, bc, bd;
	double c2, cd;
	double d2;
	double area;
	double matrix[3][3];
	void sum(Quadric q);
	double evaluateCost(double x, double y, double z);
	Quadric();
	Quadric(double n1, double n2, double n3, double v1, double v2, double v3, double area);	
};


class Vertex
{
public: 	
	double x,y,z;
	vector<int> neighborFaces;
	vector<Edge*> neighborEdges;
	Quadric quad;
	int match;
	Vertex();
	~Vertex();
};


