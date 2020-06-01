#ifndef TESTNJ_H
#define TESTNJ_H

#include "polytree.h"
#include "distMatrixReader.hpp"

class testNJ
{
 public:
  testNJ(distMatrixReader* datareader, int matrixSize);
  ~testNJ(void);
  void step(int,int,distType);
  void stepNoVal(int idx1, int idx2);
  
 private:
  distType** matrix;
  polytree* mytree;
  int matrixSize;	
  distType* separationsums;
  distType* separations;
  int clusterCount;
  int min1;
  int min2;
  double min;
  int* activeRows;
  int nextId;
  distMatrixReader* reader;
  
  void findMin();
  void initialize();
  void mergeMinNodes();
  void updateMatrix();
};
#endif
