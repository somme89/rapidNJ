#ifndef NAIVEHAMMING_H
#define NAIVEHAMMING_H

#include <fstream>
#include "stdinclude.h"
#include "dataloader.hpp"

using namespace std;

class hammingDistance {
  
 public:
  hammingDistance(bool fastdist, dataloader* dataloader);
  void computeDistanceMatrix();
  distType** getDistanceMatrix();  
  ~hammingDistance();
  void printDistances();
  
 private:
  void postProcessDistanceMatrix();

  char** sequences;
  int seqLength;
  int seqCount;
  InputType type;
  int** distanceMatrixTS;
  int** distanceMatrixTV;
  distType** distanceMatrix;
  vector<string> sequenceNames;
  static const int BUFFER_SIZE = 16384;
  void readHeader();
  void parseName(char* input, int seq);
  int compareSequences(int i, int j);
  void createDataStructures();
  void initializeData(string);
  dataloader* loader;
  DistanceEstimate* distEstimator;
};

#endif
