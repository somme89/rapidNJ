#ifndef SIMPLENJ_H
#define SIMPLENJ_H

#include "polytree.h"
#include "distMatrixReader.hpp"

class simpleNJ
{
 public:
  simpleNJ(distMatrixReader* datareader, int matrixSize, bool negative_branches, ProgressBar* pb);
  ~simpleNJ(void);
  polytree* run();
  
 private:
  distType** matrix;
  polytree* mytree;
  int matrixSize;	
  ProgressBar* pb;
  distType* separationsums;
  distType* separations;
  int clusterCount;
  int min1;
  int min2;
  int* activeRows;
  int nextId;
  distMatrixReader* reader;
  bool negative_branches;
  
  void findMin();
  void initialize();
  void mergeMinNodes();
  void updateMatrix();
};
#endif
