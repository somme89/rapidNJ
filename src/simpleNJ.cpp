/*A simple implementation of the neighbour-joining method*/

#include "stdinclude.h"
#include "simpleNJ.h"
#include "float.h"

using namespace std;

void printMatrix(distType** matrix, int size);
void printArray(distType* a, int size);
int countersim = 0;

simpleNJ::simpleNJ(distMatrixReader* reader, int matrixSize, bool negative_branches, ProgressBar* pb)
{
  simpleNJ::matrixSize = matrixSize;
  simpleNJ::reader = reader;
  simpleNJ::negative_branches = negative_branches;
  simpleNJ::pb = pb;
  matrix = reader->getMatrix();
  separationsums = new distType[matrixSize];
  separations = new distType[matrixSize];
  clusterCount = matrixSize;      
  activeRows = new int[matrixSize*2];
  nextId = matrixSize;
}

polytree* simpleNJ::run(){
  initialize();
  while(clusterCount > 2){
    findMin();  //O(n^2)
    mergeMinNodes();
    updateMatrix();                       
  }
    // finish by joining the two remaining clusters
  int index1 = -1;
  int index2 = -1;
  // find the last nodes
  for(int i = 0; i < matrixSize; i++){
    if(activeRows[i] != -1){
      if(index1 == -1){
	index1 = i;
      } else {
	index2 = i;	
	break;
      }            
    }
  } 
  double distance = matrix[index1][index2];
  mytree->set_serialization_indices(activeRows[index1],activeRows[index2], distance / 2.0);  
  pb->finish();
  return mytree;
}

void simpleNJ::updateMatrix(){    
  distType newSeparationsum = 0;
  distType mutualDistance = matrix[min1][min2];
  distType* row1 = matrix[min1];
  distType* row2 = matrix[min2];
  for(int i = 0; i < matrixSize; i++){
    if(i == min1 || i == min2 || activeRows[i] == -1){
      row1[i] = 0;
    } else {
      distType val1 = row1[i];
      distType val2 = row2[i];
      distType dist = (val1 + val2 - mutualDistance) / ((distType)2.0);
      if(simpleNJ::negative_branches && dist < 0) {
        dist = 0;
      }
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
  activeRows[min2] = -1;
  activeRows[min1] =  nextId++ ;  
}

void simpleNJ::mergeMinNodes() {
  // calculate distances
  double dist = matrix[min1][min2];
  double sep1 = separations[min1];
  double sep2 = separations[min2];
  double dist1 = (0.5 * dist) + (0.5 * (sep1 - sep2));
  double dist2 = (0.5 * dist) + (0.5 * (sep2 - sep1));
  if(simpleNJ::negative_branches) {
    if(dist1 < 0.0) {
      dist2 += dist1;
      dist1 = 0;
    }
    if(dist2 < 0.0) {
      dist1 += dist2;
      if(dist1 < 0) {
        dist1 = 0;
      }
      dist2 = 0;
    }
  }
  // update tree
  mytree->addInternalNode(dist1, dist2, activeRows[min1], activeRows[min2]);
  clusterCount--;
  pb->setProgress((matrixSize - clusterCount) / (double) matrixSize);
}

void simpleNJ::initialize(){
  //calculate initial seperation rows
  mytree = new polytree(matrixSize, reader->getSequenceNames());  
  for(int i = 0; i < matrixSize; i++){
    distType sum = 0;
    for(int j = 0; j < matrixSize; j++){
      sum += matrix[i][j];
    }
    separationsums[i] = sum;
    separations[i] = sum / (clusterCount - 2); 
    activeRows[i] = i;
  }
}

void simpleNJ::findMin() {
  min1 = -1;
  min2 = -1;
  double min = DBL_MAX;        
  for (int i = 0; i < matrixSize; i++) {    
    if(activeRows[i] != -1){
      distType* row = matrix[i];
      double sep1 = separations[i];
      for(int j = 0; j < matrixSize; j++){
	if(activeRows[j] != -1 && i != j){                    
	  double sep2 = separations[j];
	  double val = row[j] - sep1 - sep2;                   

	  if(val < min){
	    // new minimum
	    min1 = i;
	    min2 = j;
	    min = val;
	  }					
	}
      }
    }
  }
  //cout << "JOINING min1=" << min1 << " min2=" << min2 << " global_min=" << min << " clusterCount=" << clusterCount << "\n";
}

simpleNJ::~simpleNJ(void)
{
}
