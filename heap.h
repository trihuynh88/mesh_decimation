#include "mesh.h"
#include "vector"
using namespace std;

void siftDown(vector<Edge*>* edges, int root, int bottom);
void updateHeap(vector<Edge*>* edges, int array_size, int pos, double value);
void buildHeap(vector<Edge*>* edges, int array_size);
void upHeap(vector<Edge*>* edges, int array_size, int pos);
Edge* popHeapTop(vector<Edge*>* edges, int &array_size);