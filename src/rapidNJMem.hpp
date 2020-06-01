#ifndef RAPIDNJMEM_H_
#define RAPIDNJMEM_H_

#include "polytree.h"
#include "distMatrixReader.hpp"
#include <algorithm>
#include "cluster_pair.h"
#include <list>

class rapidNJMem {
	
 public:
  rapidNJMem(distMatrixReader* reader, int matrixSize, int sortedMatrixSize, bool verbose, bool negative_branches, ProgressBar* pb);
  ~rapidNJMem(void);
  polytree* run();
  
 private:
  distType** distanceMatrix;
  polytree* mytree;
  int matrixSize;	
  bool negative_branches;
  ProgressBar* pb;
  distType* separationsums;
  distType* separations;
  int clusterCount;
  int min1;
  int min2;
  cluster_pair **cluster_data;
  cluster_pair *cluster_buffer;
  cluster_pair min_pair;
  int currentId;    
  int *idToIndex;
  int *indexToId;
  distType max_separation;
  double global_min;
  int* garbage_flags;    
  int min_row_cache1;
  int min_row_cache2;
  int *row_lengths;
  distMatrixReader* reader;
  std::list<int> *redundantMap;
  int *redundantCount;  
  int newRowIndex;
  int obsoleteRowIndex;
  int *redundantToMasterMap;
  int sortedMatrixSize;
  int totalCount;
  distType* maxRowSeparations;
  bool verbose;

  void findMin();
  int findRowMin(int i, int rowsize);
  void findRowMinGarbage(int i, int rowsize);
  void initialize();
  void mergeMinNodes();
  void updateData();
  int next_id();
  inline void handleRedundancy(int i, int j);
  inline void buildSortedRow(int i);
  void searchRedundantEntries(int i);

    inline distType getDist(int i, int j){
    if(i >= j){
      return distanceMatrix[i][j];
    } else {
      return distanceMatrix[j][i];
    }
  } 

  inline void setDist(int i, int j, distType value){
    if(i >= j){
      distanceMatrix[i][j] = value;
    } else {
      distanceMatrix[j][i] = value;
    }
  } 
};

#endif
