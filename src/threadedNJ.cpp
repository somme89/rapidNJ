/**The neighbour-joining method using two cores*/

#include "stdinclude.h"
#include "threadedNJ.h"
#include <pthread.h>
#include <sched.h>
#include "minFinder.h"
#include "float.h"

using namespace std;

distType** matrix;   	
node** nodes;
distType* separationsums;
distType* separations;
int matrixSize;
int clusterCount;
threadState * state1;
threadState * state2;
pthread_attr_t attrib1;
sched_param s_param;

threadedNJ::threadedNJ(datareader* reader) {
  matrixSize = reader->getSize();
  matrix = reader->getMatrix();
  nodes = reader->getNodes();
  separationsums = new distType[matrixSize];
  separations = new distType[matrixSize];
  clusterCount = matrixSize;                  
}

void* threadedNJ::njThread( void* ptr ) {
  threadState* state = (threadState*) ptr;	
  executefind(state);
  return 0;
}

node* threadedNJ::run() {		
  initialize();
  while(clusterCount > 2){
    findMin(); 		
    mergeMinNodes(); 
    updateMatrix();                       
  }
  // finish by joining the two remaining clusters
  node* node1 = NULL;
  node* node2 = NULL;
  int index1 = -1;
  int index2 = -1;
  // find the last nodes
  for(int i = 0; i < matrixSize; i++){
    node* node = nodes[i];
    if(node != NULL){
      if(index1 == -1){
	node1 = node;
	index1 = i;
      } else {
	node2 = node;
	index2 = i;
	break;
      }            
    }
  }
  double distance = matrix[index1][index2];
  node1->addEdge(node2,distance);
  node2->addEdge(node1,distance);
  return node2;
}

void threadedNJ::initialize(){
   //calculate initial seperation rows
  for(int i = 0; i < matrixSize; i++){
    double sum = 0;
    for(int j = 0; j < matrixSize; j++){
      sum += matrix[i][j];
    }
    separationsums[i] = sum;
    separations[i] = sum / (clusterCount - 2); 
  }  
  // create threads
  state1 = new threadState();  
  state2 = new threadState();  
  state1->rowstart = 0;
  state1->rowEnd = (matrixSize / 2);  
  state2->rowstart = matrixSize / 2;
  state2->rowEnd = matrixSize;
}

void threadedNJ::updateMatrix(){    
  double newSeparationsum = 0;
  double mutualDistance = matrix[min1][min2];
  distType* row1 = matrix[min1];
  distType* row2 = matrix[min2];
  for(int i = 0; i < matrixSize; i++){
    if(i == min1 || i == min2 || nodes[i] == NULL){
      row1[i] = 0;
    } else {
      double val1 = row1[i];
      double val2 = row2[i];
      double dist = (val1 + val2 - mutualDistance) / 2.0;
      newSeparationsum += dist;
      // update the separationsum of cluster i.
      separationsums[i] += (dist - val1 - val2);
      separations[i] = separationsums[i] / (clusterCount -2); 
      row1[i] = dist;
      matrix[i][min1] = dist;
    }
  }
  separationsums[min1] = newSeparationsum;
  separations[min1] = newSeparationsum / (clusterCount - 2);
  separationsums[min2] = 0;
}

void threadedNJ::mergeMinNodes(){        
  // calculate distances
  node* node1 = nodes[min1];
  node* node2 = nodes[min2];
  node* newNode = new node();
  double dist = matrix[min1][min2];
  double sep1 = separations[min1];
  double sep2 = separations[min2];
  double dist1 = (0.5 * dist) + (0.5 * (sep1 - sep2));
  double dist2 = (0.5 * dist) + (0.5 * (sep2 - sep1));
  // update tree  
  newNode->addEdge(node1, dist1);
  node1->addEdge(newNode, dist1);
  newNode->addEdge(node2, dist2);
  node2->addEdge(newNode, dist2);
  // update data  
  nodes[min1] = newNode;
  nodes[min2] = NULL;
  clusterCount--;
}

void threadedNJ::findMin() {
  pthread_t thread1, thread2;
  pthread_create( &thread1, NULL, threadedNJ::njThread, (void*) state1);
  pthread_create( &thread2, NULL, threadedNJ::njThread, (void*) state2);  
  pthread_join(thread1, NULL);
  pthread_join(thread2, NULL);
  if(state1->min < state2->min){
    min1 = state1->min1;
    min2 = state1->min2;
  } else {
    min1 = state2->min1;
    min2 = state2->min2;	
  }
}

void executefind(threadState *state){  
  state->min1 = -1;
  state->min2 = -1;
  double min = DBL_MAX;
  for (int i = state->rowstart; i < state->rowEnd; i++) {
    if(nodes[i] != NULL){
      distType* row = matrix[i];
      double sep1 = separations[i];
      for(int j = 0; j < matrixSize; j++){
	if(nodes[j] != NULL && i != j){
	  double sep2 = separations[j];
	  double val = row[j] - sep1 - sep2;                   
	  if(val < min){
	    // new minimum
	    state->min1 = i;
	    state->min2 = j;
	    min = val;
	  }
	}
      }
    }
  }
  state->min = min;
}

threadedNJ::~threadedNJ()
{
}
