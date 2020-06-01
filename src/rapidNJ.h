#ifndef RAPIDNJ_H_
#define RAPIDNJ_H_

#include "polytree.h"
#include "distMatrixReader.hpp"
#include <algorithm>
#include "cluster_pair.h"
#include <list>

class rapidNJ {
	
 public:
  rapidNJ(distMatrixReader* datareader, int matrix_size, bool negative_branches, ProgressBar* pb);
  ~rapidNJ(void);
  polytree* run();
  
 private:
  distType** matrix;
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
  int currentId;    
  int *idToIndex;
  int *indexToId;
  distType max_separation;
  distType global_min;
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
  distType* maxRowSeparations;

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
};

#endif
