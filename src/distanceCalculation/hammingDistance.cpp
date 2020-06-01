#include "simpleDistanceCalculator.hpp"
#include "bitDistanceGap.hpp"
#include "bitDistanceProtein.hpp"
#include "DistanceEstimate.hpp"
#include "hammingDistance.hpp"
#include <math.h>

hammingDistance::hammingDistance(bool fastdist, dataloader* loader) {
  hammingDistance::seqCount = loader->getSequenceCount();
  hammingDistance::seqLength = loader->getSequenceLength();  
  hammingDistance::sequenceNames = *loader->getSequenceNames();
  hammingDistance::type = loader->type;
  hammingDistance::loader = loader;
  if(loader->fastdist) {
    if(type == DNA) {
      distEstimator = new bitDistanceGap(loader);
    } else if(type == PROTEIN) {
      distEstimator = new bitDistanceProtein(loader);
    } else {
      cerr << "ERROR: Unknown sequence type \"" << type << "\"" << endl;
    }
  } else {
    distEstimator = new simpleDistanceCalculator(loader);
  }
}

void hammingDistance::computeDistanceMatrix() {
  distanceMatrix = new distType*[seqCount];
  unsigned long long data[3];
  
  for(int i = 0; i < seqCount; i++){
    distanceMatrix[i] = new distType[seqCount];
    for(int j = i+1; j < seqCount; j++){
      distEstimator->computeDistance(i,j, data);
      distanceMatrix[i][j] = (distType)(data[0]+data[1]) / (distType) data[2];
    }
  }
  postProcessDistanceMatrix();
}

void hammingDistance::postProcessDistanceMatrix(){
  for (int i = 0; i < seqCount; i++) {
    for (int j = 0; j < seqCount; j++) {
      if(i > j){
	distanceMatrix[i][j] = distanceMatrix[j][i];
      } else if(j == i){
	distanceMatrix[i][j] = 0;
      } 
    }
  }
}

void hammingDistance::printDistances(){
  for (int i = 0; i < seqCount; i++) {
    cout << sequenceNames[i] << "\t";
    for (int j = 0; j < seqCount; j++) {
      cout << distanceMatrix[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
}

// Test method
/*void hammingDistance::getTsTv() {
  simpleDistanceCalculator* sc = new simpleDistanceCalculator(verbose, loader);
  unsigned int distances[2];
  distanceMatrixTS = new int*[seqCount];
  distanceMatrixTV = new int*[seqCount];
  for(int i = 0; i < seqCount; i++){
    distanceMatrixTS[i] = new int[seqCount];
    distanceMatrixTV[i] = new int[seqCount];
    for(int j = i+1; j < seqCount; j++){
      sc->calculateDistanceDNA_TV_TS(i,j, distances);
      distanceMatrixTS[i][j] = distances[0];
      distanceMatrixTV[i][j] = distances[1];
    }
  }
  for (int i = 0; i < seqCount; i++) {
    cout << sequenceNames[i] << "\t";
    for (int j = 0; j < seqCount; j++) {
      if(i <= j){
	cout << "[" << distanceMatrixTS[i][j] << "," << distanceMatrixTV[i][j] << "]\t";
      } else {
	cout << "[" << distanceMatrixTS[j][i] << "," << distanceMatrixTV[j][i] << "]\t";
      }
    }
    cout << endl;
  }
}*/

distType** hammingDistance::getDistanceMatrix(){
  return distanceMatrix;
}

