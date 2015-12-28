// SiemensDecimation.cpp : Defines the entry point for the console application.
//


//#include "stdafx.h"
#include "heap.h"
#include "mesh.h"
#include "vector"
#include "decimation.h"
#include "time.h"
#include "fstream"
//#include "MemLeakCheck.h"

using namespace std;



int main(int argc, char* argv[])
{
//	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	clock_t t = clock();
	float totalTime = 0;
	Decimation decimation;	
	//int* a = new int[10];

	MeshSegment mesh;
	
	ifstream infile("mesh_withNormal.smf");
	char line[50];
	Vertex vertex;
	Face face;
	
	while (!infile.eof())
	{
			
		if (infile.getline(line,50,'\n').good())
		{	
			switch (line[0])
			{
			case 'v':			
				sscanf(line,"v %lf %lf %lf",&vertex.x,&vertex.y,&vertex.z);

				mesh.m_vCoord.push_back(vertex.x);
				mesh.m_vCoord.push_back(vertex.y);
				mesh.m_vCoord.push_back(vertex.z);
				break;
			case 'f':			
				sscanf(line,"f %d %d %d",&face.v1,&face.v2,&face.v3);
				face.v1-=1;
				face.v2-=1;
				face.v3-=1;
				
				mesh.m_vTriangleIndex.push_back(face.v1);
				mesh.m_vTriangleIndex.push_back(face.v2);
				mesh.m_vTriangleIndex.push_back(face.v3);

				
				break;
			default:
				break;
			}
			
		}
	}

	infile.close();
	

	decimation.setModel(mesh);
	
	printf("Loading file time: %f s\n",(float)(clock()-t)/CLOCKS_PER_SEC);
	totalTime+=(float)(clock()-t)/CLOCKS_PER_SEC;
	t = clock();
		
	decimation.decimate(3000);
	printf("Decimation time: %f s\n",(float)(clock()-t)/CLOCKS_PER_SEC);
	totalTime+=(float)(clock()-t)/CLOCKS_PER_SEC;
	t = clock();
	
	
	mesh = decimation.getResultModel();

	ofstream outfile("mesh3000.smf");
	for (int i = 0; i<mesh.m_vCoord.size()/3; i++)
		outfile << "v " << mesh.m_vCoord[i*3] << " " << mesh.m_vCoord[i*3+1] << " " << mesh.m_vCoord[i*3+2] << "\n";
	for (i = 0; i<mesh.m_vTriangleIndex.size()/3; i++)
		outfile << "f " << mesh.m_vTriangleIndex[i*3] + 1 << " " << mesh.m_vTriangleIndex[i*3+1]+1 << " " << mesh.m_vTriangleIndex[i*3+2]+1 << "\n";
	outfile << "bind n face\n";
	for (i = 0; i<mesh.m_vTriangleNormal.size()/3; i++)
		outfile << "n " << mesh.m_vTriangleNormal[i*3] << " " << mesh.m_vTriangleNormal[i*3+1] << " " << mesh.m_vTriangleNormal[i*3+2] << "\n";
	outfile.close();


	//Test
	printf("Writing file time: %f s\n",(float)(clock()-t)/CLOCKS_PER_SEC);
	totalTime+=(float)(clock()-t)/CLOCKS_PER_SEC;

	printf("Total time is: %f s\n",totalTime);

	
	ofstream outfileNormal("meshNormal3000.txt");
	for (i = 0; i<mesh.m_vTriangleNormal.size()/3; i++)
		outfileNormal << mesh.m_vTriangleNormal[i*3] << " " << mesh.m_vTriangleNormal[i*3+1] << " " << mesh.m_vTriangleNormal[i*3+2] << "\n";
	outfileNormal.close();

	ofstream outfileTriangle("meshTriangle3000.txt");
	for (i = 0; i<mesh.m_vTriangleIndex.size()/3; i++)
		outfileTriangle << mesh.m_vTriangleIndex[i*3] << " " << mesh.m_vTriangleIndex[i*3+1] << " " << mesh.m_vTriangleIndex[i*3+2] << "\n";
	outfileTriangle.close();

	ofstream outfilePoint("meshPoint3000.txt");
	for (i = 0; i<mesh.m_vCoord.size()/3; i++)
		outfilePoint << mesh.m_vCoord[i*3] << " " << mesh.m_vCoord[i*3+1] << " " << mesh.m_vCoord[i*3+2] << "\n";
	outfilePoint.close();

	return 0;
}

