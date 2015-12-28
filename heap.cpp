#include "heap.h"
#include "vector"
//#include "MemLeakCheck.h"
using namespace std;

void siftDown(vector<Edge*>* edges, int root, int bottom) {
  int maxChild = root * 2 + 1;
 
  // Find the biggest child
  if(maxChild < bottom) {
    int otherChild = maxChild + 1;
    // Reversed for stability
    maxChild = ((*edges)[otherChild]->cost > (*edges)[maxChild]->cost)?otherChild:maxChild;
  } else {
    // Don't overflow
    if(maxChild > bottom) return;
  }
 
  // If we have the correct ordering, we are done.
  if((*edges)[root]->cost >= (*edges)[maxChild]->cost) return;
 
  // Swap
  Edge* temp = (*edges)[root];
  (*edges)[root] = (*edges)[maxChild];
  (*edges)[maxChild] = temp;
  (*edges)[maxChild]->posInHeap = maxChild;
  (*edges)[root]->posInHeap = root;

  // Tail queue recursion. Will be compiled as a loop with correct compiler switches.
  siftDown(edges, maxChild, bottom);
}

void buildHeap(vector<Edge*>* edges, int array_size)
{
	 for (int i = (array_size / 2); i >= 0; i--) {
		siftDown(edges, i, array_size - 1);
	  }
}

void upHeap(vector<Edge*>* edges, int array_size, int pos) {
	if (pos == 0)
		return;
	int parent = int((pos - 1)/2);
	Edge* temp;
	if ((*edges)[pos]->cost > (*edges)[parent]->cost) 
	{
		temp = (*edges)[pos];
		(*edges)[pos] = (*edges)[parent];
		(*edges)[parent] = temp;		

		(*edges)[pos]->posInHeap = pos;
		(*edges)[parent]->posInHeap = parent;		


		upHeap(edges,array_size,parent);
	}
}

void updateHeap(vector<Edge*>* edges, int array_size, int pos, double value) {
	(*edges)[pos]->cost = value;
	int parent = int((pos - 1)/2);
	Edge* temp;
	if ((*edges)[pos]->cost > (*edges)[parent]->cost) 
	{
		temp = (*edges)[pos];
		(*edges)[pos] = (*edges)[parent];
		(*edges)[parent] = temp;		

		(*edges)[pos]->posInHeap = pos;
		(*edges)[parent]->posInHeap = parent;		

		upHeap(edges,array_size,parent);
	}
	else
	{
		siftDown(edges,pos,array_size - 1);
	}
}


Edge* popHeapTop(vector<Edge*>* edges, int &array_size)
{
	Edge* temp;
	if (array_size == 0)
		return NULL;
	do
	{
		temp = (*edges)[0];
		(*edges)[0] = (*edges)[array_size-1];
		(*edges)[array_size-1] = temp;
		(*edges)[0]->posInHeap = 0;
		(*edges)[array_size-1]->posInHeap = array_size-1;
		array_size--;
		siftDown(edges, 0, array_size-1);
	}
	while ((!temp->isValid)&&(array_size>0));
	if (temp->isValid)
		return temp;
	else
		return NULL;
}