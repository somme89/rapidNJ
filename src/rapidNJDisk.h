#ifndef RAPIDNJDISK_H_
#define RAPIDNJDISK_H_

#include "polytree.h"
#include "cluster_pair.h"
#include "diskMatrix.h"
#include "rdDataInitialiser.h"
#include <list>
#include <vector>

class rapidNJDisk {
	
 public:
  rapidNJDisk(rdDataInitialiser* datareader, bool verbose, bool negative_branches, ProgressBar* pb);
  ~rapidNJDisk(void);
  polytree* run();

 private:
  diskMatrix* disk_matrix;
  polytree* mytree;
  bool negative_branches;
  ProgressBar* pb;
  int matrixSize;
  int matrixSizeCache;
  int newRowsStartId;
  distType* separationsums;
  distType* separations;
  int clusterCount;
  int min1;
  int min2;
  cluster_pair **cluster_data;
  cluster_pair min_pair;
  cluster_pair *cluster_buffer;
  int currentId;    
  int *idToIndex;
  int *indexToId;
  distType max_separation;
  distType global_min;
  /*update_flags[i]==0 if the i'th row hasn't been cached and an id > 0 otherwise. The id is an increasing number starting from 1*/
  int* update_flags;
  int updateCount;
  int min_row_cache1;
  int min_row_cache2;
  int *sorted_row_lengths;
  distType *rowBuffer1;
  distType *rowBuffer2;
  distType* maxRowSeparations;

  /*Contains all rows of the distance matrix which has been modified since the last cache flush.*/
  distType **column_cache;
  int dataStructureSize;
  bool verbose;
  int *redundantToMasterMap;
  int newRowIndex;
  int obsoleteRowIndex;
  std::list<int> *redundantMap;
  int *redundantCount;
  int totalCount;
  short *activeRowFlags;
  short *possibleRedundantRows;

  void findMin();
  void findRowMin(int i);
  void initialize();
  void mergeMinNodes();
  void updateData();
  int next_id();
  void flushUpdateCache();
  void handleRedundancy(int i, int j);
  void rebuildDataStructures(int newRowPos, int oldRowPos, int newDataStructureSize);
  void searchRedundantEntries(int i);
};

#endif

