/*
Mesh Decimation version 1.0

  Author: Tri Huynh (hquoctri@gmail.com)
  Functionality: reduce the number of triangles from the input mesh to a desired target number.
  Acknowledgement: the algorithm is adapted from the article "Surface Simplification Using Quadric Error Metrics" of Michael Garland. & Paul S. Heckbert
					and Michael Garland's thesis on "Quadric-Based Polygonal Surface Simplification"
*/

#include "Decimation.h"
//#include "MemLeakCheck.h"
using namespace std;

#define SWAP(a, b, t)   {t = a; a = b; b = t;}

void Decimation::decimate(int nFace)
{
	if (m_isFirstRun)
	{
		initializeQuadrics();
		m_isFirstRun = false;
	}

	copyOriginalFaceAndVertex();

	initializeEdgeList();
	initializeHeap();
	m_numFace = m_faceList.size();
	
	while (m_numFace > nFace)
	{
		Edge* top = popHeapTop(&m_edgeList,m_heapCount);
		if (top == NULL)
			break;
		contract(top);
	}
	createFinalModel();
}

void Decimation::contract(Edge* edge)
{
	
	//replace v2 by v1 in the faces, also mark invalid faces, add new neighbor faces to v1
	for (int i = 0; i<m_vertexList[edge->v2].neighborFaces.size(); i++)
		if (m_faceList[m_vertexList[edge->v2].neighborFaces[i]].isValid)
		{
			if ((edge->v1 == m_faceList[m_vertexList[edge->v2].neighborFaces[i]].v1)
				|| (edge->v1 == m_faceList[m_vertexList[edge->v2].neighborFaces[i]].v2)
				|| (edge->v1 == m_faceList[m_vertexList[edge->v2].neighborFaces[i]].v3))
			{
				m_faceList[m_vertexList[edge->v2].neighborFaces[i]].isValid = false;
				m_numFace--;
			}
			else
			{
				m_vertexList[edge->v1].neighborFaces.push_back(m_vertexList[edge->v2].neighborFaces[i]);

				if (edge->v2 == m_faceList[m_vertexList[edge->v2].neighborFaces[i]].v1)
					m_faceList[m_vertexList[edge->v2].neighborFaces[i]].v1 = edge->v1;
				else
					if (edge->v2 == m_faceList[m_vertexList[edge->v2].neighborFaces[i]].v2)
						m_faceList[m_vertexList[edge->v2].neighborFaces[i]].v2 = edge->v1;
					else
						if (edge->v2 == m_faceList[m_vertexList[edge->v2].neighborFaces[i]].v3)
							m_faceList[m_vertexList[edge->v2].neighborFaces[i]].v3 = edge->v1;
			}
		}

	//calculate new quadric
	m_vertexList[edge->v1].quad.sum(m_vertexList[edge->v2].quad);	
	findOptimumPos(edge->v1,edge->v2,m_vertexList[edge->v1].x,m_vertexList[edge->v1].y,m_vertexList[edge->v1].z);
	//update Edges
	int oldV2 = edge->v2;
	//connect edges with v2 to v1
	for (i = 0; i<m_vertexList[oldV2].neighborEdges.size(); i++)
		if (m_vertexList[oldV2].neighborEdges[i]->isValid)
		{
			if (m_vertexList[oldV2].neighborEdges[i]->v1 == oldV2)
			{			
				m_vertexList[oldV2].neighborEdges[i]->v1 = edge->v1;
				
			}
			else
			{				
				m_vertexList[oldV2].neighborEdges[i]->v2 = edge->v1;
				
			}
			if (m_vertexList[oldV2].neighborEdges[i]->v2 == m_vertexList[oldV2].neighborEdges[i]->v1)
			{
				m_vertexList[oldV2].neighborEdges[i]->isValid = false;				
			}
			if (m_vertexList[oldV2].neighborEdges[i]->isValid)
				m_vertexList[edge->v1].neighborEdges.push_back(m_vertexList[oldV2].neighborEdges[i]);
		}

	//recalculate the cost
	double cost = 0;
	for (i = 0; i<m_vertexList[edge->v1].neighborEdges.size(); i++)
		if (m_vertexList[edge->v1].neighborEdges[i]->isValid)
		{			
				Quadric q;
				q.sum(m_vertexList[m_vertexList[edge->v1].neighborEdges[i]->v1].quad);
				q.sum(m_vertexList[m_vertexList[edge->v1].neighborEdges[i]->v2].quad);
				double x,y,z;
				findOptimumPos(m_vertexList[edge->v1].neighborEdges[i]->v1,m_vertexList[edge->v1].neighborEdges[i]->v2,x,y,z);				
				cost = q.evaluateCost(x,y,z);
				updateHeap(&m_edgeList,m_heapCount,m_vertexList[edge->v1].neighborEdges[i]->posInHeap,-cost);
		}

	m_vertexList[oldV2].match = edge->v1;
}

void Decimation::initializeQuadrics()
{
	
	for (int i = 0; i<m_originalFaceList.size(); i++)
	{
		double a1, a2, a3, b1,b2,b3,c1,c2,c3;
		a1 = m_originalVertexList[m_originalFaceList[i].v2].x - m_originalVertexList[m_originalFaceList[i].v1].x;
		a2 = m_originalVertexList[m_originalFaceList[i].v2].y - m_originalVertexList[m_originalFaceList[i].v1].y;
		a3 = m_originalVertexList[m_originalFaceList[i].v2].z - m_originalVertexList[m_originalFaceList[i].v1].z;

		b1 = m_originalVertexList[m_originalFaceList[i].v3].x - m_originalVertexList[m_originalFaceList[i].v1].x;
		b2 = m_originalVertexList[m_originalFaceList[i].v3].y - m_originalVertexList[m_originalFaceList[i].v1].y;
		b3 = m_originalVertexList[m_originalFaceList[i].v3].z - m_originalVertexList[m_originalFaceList[i].v1].z;

		c1 = a2*b3 - a3*b2;
		c2 = a3*b1 - a1*b3;
		c3 = a1*b2 - a2*b1;

		double area = 0.5*(c1*c1+c2*c2+c3*c3);

		Quadric q(c1,c2,c3,m_originalVertexList[m_originalFaceList[i].v1].x,m_originalVertexList[m_originalFaceList[i].v1].y,m_originalVertexList[m_originalFaceList[i].v1].z,area);
		m_originalVertexList[m_originalFaceList[i].v1].quad.sum(q);
		m_originalVertexList[m_originalFaceList[i].v2].quad.sum(q);
		m_originalVertexList[m_originalFaceList[i].v3].quad.sum(q);
	}	
}


Decimation::Decimation()
{	
}

void Decimation::createEdge(int v1, int v2)
{
	double cost = 0;
	
	Quadric q;
	q.sum(m_vertexList[v1].quad);
	q.sum(m_vertexList[v2].quad);
	
	double x,y,z;


	double inverseMat[3][3];
	if (invert(q.matrix,inverseMat)>(1e-12))
	{
		x = -inverseMat[0][0]*q.ad-inverseMat[0][1]*q.bd-inverseMat[0][2]*q.cd;
		y = -inverseMat[1][0]*q.ad-inverseMat[1][1]*q.bd-inverseMat[1][2]*q.cd;
		z = -inverseMat[2][0]*q.ad-inverseMat[2][1]*q.bd-inverseMat[2][2]*q.cd;
	}
	else
	{
		x = (m_vertexList[v1].x + m_vertexList[v2].x)/2;
		y = (m_vertexList[v1].y + m_vertexList[v2].y)/2;
		z = (m_vertexList[v1].z + m_vertexList[v2].z)/2;
	}	

	cost = q.evaluateCost(x,y,z);

	Edge* edge = new Edge(v1, v2, -1*cost);
	m_vertexList[v1].neighborEdges.push_back(edge);
	m_vertexList[v2].neighborEdges.push_back(edge);
	m_edgeList.push_back(edge);
	edge->posInHeap = m_edgeList.size() - 1;

}

void Decimation::getMiddlePoint(int v1, int v2, double &x, double &y, double &z)
{
	x = (m_vertexList[v1].x + m_vertexList[v2].x)/2;
	y = (m_vertexList[v1].y + m_vertexList[v2].y)/2;
	z = (m_vertexList[v1].z + m_vertexList[v2].z)/2;
}

bool Decimation::isEdgeExist(int v1, int v2)
{
	for (int i = 0 ; i<m_vertexList[v1].neighborEdges.size(); i++)
		if (((m_vertexList[v1].neighborEdges[i]->v1 == v1) && (m_vertexList[v1].neighborEdges[i]->v2 == v2))
			|| ((m_vertexList[v1].neighborEdges[i]->v2 == v1) && (m_vertexList[v1].neighborEdges[i]->v1 == v2)))
		{
			return true;
		}
	return false;
}


void Decimation::initializeHeap()
{
	m_heapCount = m_edgeList.size();
	buildHeap(&m_edgeList,m_heapCount);
}

double Decimation::invert(double A[3][3], double B[3][3])
{
    int N = 3;
    int i, j, k;
    double max, t, det, pivot;
	double AA[3][3];
	memcpy(AA,A,9*sizeof(double));
    /*---------- forward elimination ----------*/

    for (i=0; i<N; i++)                 /* put identity matrix in B */
        for (j=0; j<N; j++)
            B[i][j] = (double)(i==j);

    det = 1.0;
    for (i=0; i<N; i++) {               /* eliminate in column i, below diag */
        max = -1.;
        for (k=i; k<N; k++)             /* find pivot for column i */
            if (abs(AA[k][i]) > max) {
                max = abs(AA[k][i]);
                j = k;
            }
        if (max<=0.) return 0.;         /* if no nonzero pivot, PUNT */
        if (j!=i) {                     /* swap rows i and j */
            for (k=i; k<N; k++)
                SWAP(AA[i][k], AA[j][k], t);
            for (k=0; k<N; k++)
                SWAP(B[i][k], B[j][k], t);
            det = -det;
        }
        pivot = AA[i][i];
        det *= pivot;
        for (k=i+1; k<N; k++)           /* only do elems to right of pivot */
            AA[i][k] /= pivot;
        for (k=0; k<N; k++)
            B[i][k] /= pivot;
        /* we know that A(i, i) will be set to 1, so don't bother to do it */

        for (j=i+1; j<N; j++) {         /* eliminate in rows below i */
            t = AA[j][i];                /* we're gonna zero this guy */
            for (k=i+1; k<N; k++)       /* subtract scaled row i from row j */
                AA[j][k] -= AA[i][k]*t;   /* (ignore k<=i, we know they're 0) */
            for (k=0; k<N; k++)
                B[j][k] -= B[i][k]*t;
        }
    }

    /*---------- backward elimination ----------*/

    for (i=N-1; i>0; i--) {             /* eliminate in column i, above diag */
        for (j=0; j<i; j++) {           /* eliminate in rows above i */
            t = AA[j][i];                /* we're gonna zero this guy */
            for (k=0; k<N; k++)         /* subtract scaled row i from row j */
                B[j][k] -= B[i][k]*t;
        }
    }

    return det;
}

void Decimation::findOptimumPos(int v1, int v2, double &x, double &y, double &z)
{
	Quadric q;
	q.sum(m_vertexList[v1].quad);
	q.sum(m_vertexList[v2].quad);
	double inverseMat[3][3];
	if (invert(q.matrix,inverseMat)>(1e-12))
	{
		x = -inverseMat[0][0]*q.ad-inverseMat[0][1]*q.bd-inverseMat[0][2]*q.cd;
		y = -inverseMat[1][0]*q.ad-inverseMat[1][1]*q.bd-inverseMat[1][2]*q.cd;
		z = -inverseMat[2][0]*q.ad-inverseMat[2][1]*q.bd-inverseMat[2][2]*q.cd;
	}
	else
	{
		getMiddlePoint(v1,v2,x,y,z);		
	}
}

void Decimation::calculateNormal(int f)
{
		double a1, a2, a3, b1,b2,b3;
		a1 = m_vertexList[m_faceList[f].v2].x - m_vertexList[m_faceList[f].v1].x;
		a2 = m_vertexList[m_faceList[f].v2].y - m_vertexList[m_faceList[f].v1].y;
		a3 = m_vertexList[m_faceList[f].v2].z - m_vertexList[m_faceList[f].v1].z;

		b1 = m_vertexList[m_faceList[f].v3].x - m_vertexList[m_faceList[f].v1].x;
		b2 = m_vertexList[m_faceList[f].v3].y - m_vertexList[m_faceList[f].v1].y;
		b3 = m_vertexList[m_faceList[f].v3].z - m_vertexList[m_faceList[f].v1].z;

		m_faceList[f].n1 = a2*b3 - a3*b2;
		m_faceList[f].n2 = a3*b1 - a1*b3;
		m_faceList[f].n3 = a1*b2 - a2*b1;
}

void Decimation::initializeEdgeList()
{
	m_edgeList.clear();	
	for (int i = 0; i<m_faceList.size(); i++)
	{
			if (!isEdgeExist(m_faceList[i].v1,m_faceList[i].v2))
			{
				createEdge(m_faceList[i].v1,m_faceList[i].v2);
			}
			if (!isEdgeExist(m_faceList[i].v2,m_faceList[i].v3))
			{
				createEdge(m_faceList[i].v2,m_faceList[i].v3);
			}
			if (!isEdgeExist(m_faceList[i].v3,m_faceList[i].v1))
			{
				createEdge(m_faceList[i].v3,m_faceList[i].v1);
			}
	}
}

void Decimation::copyOriginalFaceAndVertex()
{
	m_vertexList.resize(m_originalVertexList.size());
	m_faceList.resize(m_originalFaceList.size());
	for (int i = 0; i<m_originalVertexList.size(); i++)
		m_vertexList[i] = m_originalVertexList[i];
	for (i = 0; i<m_originalFaceList.size(); i++)
		m_faceList[i] = m_originalFaceList[i];
}

Decimation::~Decimation()
{
	m_faceList.clear();
	m_vertexList.clear();	

	m_originalFaceList.clear();
	m_originalVertexList.clear();

	for (int i = 0; i<m_edgeList.size(); i++)
		delete m_edgeList[i];
	
	m_edgeList.clear();	
}


void Decimation::setModel(_MeshSegment & mesh)
{	
	m_isFirstRun = true;
	m_originalFaceList.clear();
	m_originalVertexList.clear();
	Vertex v;
	Face f;
	for (int i = 0; i<(mesh.m_vCoord.size()/3); i++)
	{		
		v.x = mesh.m_vCoord[i*3];	
		v.y = mesh.m_vCoord[i*3+1];
		v.z = mesh.m_vCoord[i*3+2];
		m_originalVertexList.push_back(v);
	}

	for (i = 0; i<(mesh.m_vTriangleIndex.size()/3); i++)
	{
		f.v1 = mesh.m_vTriangleIndex[i*3];
		f.v2 = mesh.m_vTriangleIndex[i*3+1];
		f.v3 = mesh.m_vTriangleIndex[i*3+2];
		m_originalFaceList.push_back(f);
		m_originalVertexList[f.v1].neighborFaces.push_back(m_originalFaceList.size()-1);
		m_originalVertexList[f.v2].neighborFaces.push_back(m_originalFaceList.size()-1);
		m_originalVertexList[f.v3].neighborFaces.push_back(m_originalFaceList.size()-1);
	}
}

void Decimation::createFinalModel()
{
	m_mesh.m_vCoord.clear();
	m_mesh.m_vTriangleIndex.clear();
	m_mesh.m_vTriangleNormal.clear();

	int v1, v2, v3;
	double x,y,z;
	m_countRemovedVertex = (int *) malloc( (m_vertexList.size()+1) * sizeof(int) );
	m_countRemovedVertex[0] = 0;
	
	for (int i = 0; i<m_vertexList.size(); i++)
		if (m_vertexList[i].match<0)
		{
			x = m_vertexList[i].x;
			y = m_vertexList[i].y;
			z = m_vertexList[i].z;
			m_mesh.m_vCoord.push_back(m_vertexList[i].x);
			m_mesh.m_vCoord.push_back(m_vertexList[i].y);
			m_mesh.m_vCoord.push_back(m_vertexList[i].z);
			m_countRemovedVertex[i+1] = m_countRemovedVertex[i];
		}
		else
		{
			m_countRemovedVertex[i+1] = m_countRemovedVertex[i] + 1;
		}
	for (i = 0; i<m_faceList.size(); i++)
		if (m_faceList[i].isValid)
		{
			v1 = m_faceList[i].v1 - m_countRemovedVertex[m_faceList[i].v1+1];
			v2 = m_faceList[i].v2 - m_countRemovedVertex[m_faceList[i].v2+1];
			v3 = m_faceList[i].v3 - m_countRemovedVertex[m_faceList[i].v3+1];

			calculateNormal(i);
			m_mesh.m_vTriangleIndex.push_back(v1);
			m_mesh.m_vTriangleIndex.push_back(v2);
			m_mesh.m_vTriangleIndex.push_back(v3);

			m_mesh.m_vTriangleNormal.push_back(m_faceList[i].n1);
			m_mesh.m_vTriangleNormal.push_back(m_faceList[i].n2);
			m_mesh.m_vTriangleNormal.push_back(m_faceList[i].n3);
		}
        
    for (i = 0; i<m_edgeList.size(); i++)
        delete m_edgeList[i]; 

    m_edgeList.clear(); 

    free(m_countRemovedVertex);
}
