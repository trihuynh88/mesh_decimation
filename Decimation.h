/*
Mesh Decimation version 1.0

  Author: Tri Huynh (tri_huynh@siemens.com)
  Functionality: reduce the number of triangles from the input mesh to a desired target number.
  Acknowledgement: the algorithm is adapted from the article "Surface Simplification Using Quadric Error Metrics" of Michael Garland. & Paul S. Heckbert
					and Michael Garland's thesis on "Quadric-Based Polygonal Surface Simplification"
*/

#include "mesh.h"
#include "heap.h"

class Decimation
{
public:
/*
	typedef struct _MeshSegment 
	{
		std::vector<double> m_vCoord; 
		std::vector<int> m_vTriangleIndex; 
		std::vector<double> m_vTriangleNormal; 

	} MeshSegment;
	*/
	
	void decimate(int nFace);
	void setModel(MeshSegment & mesh);
	MeshSegment getResultModel() { return m_mesh; }
	Decimation();
	~Decimation();

private:
	int m_numFace;
	int m_heapCount;
	MeshSegment m_mesh;
	bool m_isFirstRun;

	vector<Vertex> m_originalVertexList;
	vector<Face> m_originalFaceList;		

	vector<Vertex> m_vertexList;
	vector<Face> m_faceList;	
	vector<Edge*> m_edgeList;
	int* m_countRemovedVertex;

	void initializeQuadrics();
	void createEdge(int v1, int v2);
	void getMiddlePoint(int v1, int v2, double &x, double &y, double &z);
	void initializeHeap();
	void contract(Edge* edge);
	double invert(double A[3][3], double B[3][3]);
	void findOptimumPos(int v1, int v2, double &x, double &y, double &z);			
	void calculateNormal(int f);
	bool isEdgeExist(int v1, int v2);
	void initializeEdgeList();
	void copyOriginalFaceAndVertex();
	void createFinalModel();
	
};