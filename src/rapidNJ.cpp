/*An implementation of RapidNJ */

#include "stdinclude.h"
#include "rapidNJ.h"
#include "float.h"
#include <algorithm>

#if TESTSORT
#include "distMatrixReader.hpp"
#include "testNJ.h"
testNJ* test;
#endif

using namespace std;

/* initialize datastructures which doesn't require processing of data */
rapidNJ::rapidNJ(distMatrixReader* reader, int matrixSize, bool negative_branches, ProgressBar* pb) {
  rapidNJ::reader = reader;
  rapidNJ::matrixSize = matrixSize;
  rapidNJ::negative_branches = negative_branches;
  rapidNJ::pb = pb;
  matrix = reader->getMatrix();
  separationsums = new distType[matrixSize];
  separations = new distType[matrixSize];
  row_lengths = new int[matrixSize];
  idToIndex = new int[2*matrixSize];
  indexToId = new int[matrixSize];
  garbage_flags = new int[matrixSize];
  redundantCount = new int[matrixSize];
  redundantMap = new list<int>[matrixSize];
  cluster_data = new cluster_pair*[matrixSize];
  redundantToMasterMap = new int[matrixSize];
  maxRowSeparations = new distType[matrixSize];  
  for(int i = 0; i < matrixSize; i++){
    cluster_data[i] = new cluster_pair[matrixSize];    
  } 
#if TESTSORT
  distMatrixReader* reader2 = new distMatrixReader(false,reader->getFileName(), matrixSize, false);
  reader2->read_data();
  test = new testNJ(reader2, matrixSize);
#endif
}

/* starts the algorithm */
polytree* rapidNJ::run(){  
  initialize();
  while(clusterCount > 2){	  
    findMin();  //O(n^2)		
    mergeMinNodes();
    clusterCount--;
    updateData();  
	//cout << clusterCount << ": " << min1 << " " << min2 << " " << indexToId[1503] << endl;
  }
  // finish by joining the two remaining clusters
  int index1 = -1;
  int index2 = -1;
  // find the last nodes
  for(int i = 0; i < matrixSize; i++){
    if(indexToId[i] != -1){
      if(index1 == -1){
        if(redundantCount[i] == 2){
          index1 = i;
          index2 = redundantMap[i].front();
          indexToId[index2] = index2;
          break;
        } else {
          index1 = i;
        }
      } else {
        index2 = i;
        break;
      }
    }
  }

#if TESTSORT
  exit(0);
#endif
  if(index1 == -1 || index2 == -1) {
    cerr << "ERROR: an error occured while constructing the tree" << endl;
    exit(1);
  }
  distType distance = matrix[index1][index2];
  mytree->set_serialization_indices(indexToId[index1],indexToId[index2], distance);
  pb->finish();

  return mytree;
}

/* Initialize datastructures which require pre-processing */
void rapidNJ::initialize(){
  mytree = new polytree(matrixSize, reader->getSequenceNames());  
  clusterCount = matrixSize;
  currentId = matrixSize;
  min_row_cache1 = 0;
  min_row_cache2 = 0;
  max_separation = FLT_MIN;
  min1 = 0;
  min2 = 0;
  for(int i = 0; i < matrixSize; i++){
    indexToId[i] = i;
    idToIndex[i] = i;
    garbage_flags[i] = 0;
    redundantCount[i] = 1;
    //cluster_data[i] = new cluster_pair[matrixSize];    
    redundantToMasterMap[i] = i;
    distType sum = 0;
    for(int j = 0; j < matrixSize; j++){      
      sum += matrix[i][j];
    }
    separationsums[i] = sum;
    distType separation = sum / (clusterCount -2);
    separations[i] = separation;				
    //Set maxRowSeparation here to avoid getting the rows own separation as max
    maxRowSeparations[i] = max_separation;
    // update the maximum separation
    if(separation > max_separation){
      max_separation = separation;
    }
  }
  // Handle redundancy
  for(int i = matrixSize-1; i >= 0; i--){
    if(indexToId[i] == -1){
      // this row has allready been identified as redundant
      continue;
    }
    for(int j = 0; j < i; j++){      
      if(matrix[i][j] == 0){
        handleRedundancy(i,j);
      }
    }
  }
  // Build the sorted matrix
  for(int i = 0; i < matrixSize; i++){
    if(indexToId[i] != -1){
      buildSortedRow(i);
    }
  }
}

/* makes a sorted list of clusters for a row */
inline void rapidNJ::buildSortedRow(int i){  
  int idx = 0;
  //lower left triangle
  for(int j = 0; j < i; j++){
    if(indexToId[j] != -1){
      cluster_data[i][idx].distance = matrix[i][j];			
      cluster_data[i][idx].id = j;
      idx++;
    }
  }
  // sort the row				
  sort(&cluster_data[i][0],&cluster_data[i][idx]);
  row_lengths[i] = idx;  
}

inline void rapidNJ::handleRedundancy(int i, int j){  
  for(int k = 0; k < matrixSize; k++){
    if(matrix[i][k] != matrix[j][k]){
      // not redundant
      return;
    }
  }
  indexToId[j] = -1;  
  redundantToMasterMap[j] = i;
  redundantMap[i].push_front(j);
  redundantCount[i]++;  
}

void rapidNJ::updateData(){
  distType newSeparationsum = 0;
  //  distType mutualDistance = min_pair.distance;  
  distType mutualDistance = matrix[redundantToMasterMap[min1]][redundantToMasterMap[min2]];
  distType* row1 = matrix[redundantToMasterMap[min1]];
  distType* row2 = matrix[redundantToMasterMap[min2]];  
  distType* newRow = matrix[newRowIndex];
  max_separation = FLT_MIN;
  for(int curId = 0; curId < currentId; curId++){
    int i = idToIndex[curId];    
    if(i == -1 || indexToId[i] == -1){
      continue;
    }
    if(i == newRowIndex || i == obsoleteRowIndex){
      newRow[i] = 0;
    } else {
      distType val1 = row1[i];
      distType val2 = row2[i];
      distType dist = ((val1 + val2 - mutualDistance) / 2.0);
      if(rapidNJ::negative_branches && dist < 0) {
        dist = 0;
      }
      newSeparationsum += dist * redundantCount[i];
      // update the separationsum of cluster i.
      separationsums[i] += (dist - val1 - val2);
      distType separation = separationsums[i] / (clusterCount - 2);
      separations[i] = separation;
      // update the maximum separation
      maxRowSeparations[i] = max_separation;
      if(separation > max_separation){
        max_separation = separation;
      }  
      newRow[i] = dist;
      matrix[i][newRowIndex] = dist;
    }
  }
  // delete obsolete data  
  idToIndex[indexToId[newRowIndex]] = -1;
  idToIndex[indexToId[obsoleteRowIndex]] = -1;
  indexToId[obsoleteRowIndex] = -1;
  separationsums[newRowIndex] = newSeparationsum;
  separationsums[obsoleteRowIndex] = 0;
  row_lengths[obsoleteRowIndex] = 0;
  separations[newRowIndex] = newSeparationsum / (clusterCount - 2);
  maxRowSeparations[newRowIndex] = max_separation; 
  separations[obsoleteRowIndex] = FLT_MIN;
  redundantToMasterMap[newRowIndex] = newRowIndex;
  redundantToMasterMap[obsoleteRowIndex] = obsoleteRowIndex;

  // create new cluster pair
  cluster_pair *new_data = cluster_data[newRowIndex];
  int nextId = next_id();
  indexToId[newRowIndex] = nextId;
  idToIndex[nextId] = newRowIndex;
  int idx = 0;

  for(int i = 0; i < matrixSize; i++) {
    if(i != newRowIndex && indexToId[i] != -1){
      new_data[idx].distance = newRow[i];
      new_data[idx].id = indexToId[i];   
      idx++;
    }
  }  
  row_lengths[newRowIndex] = idx;
  sort(&new_data[0],&new_data[idx]);
  cluster_data[newRowIndex] = new_data;
}


/* merge two clusters */
void rapidNJ::mergeMinNodes(){  
  // calculate distances
  double dist = matrix[redundantToMasterMap[min1]][redundantToMasterMap[min2]];
  double sep1 = separations[redundantToMasterMap[min1]];
  double sep2 = separations[redundantToMasterMap[min2]];
  double dist1 = (0.5 * dist) + (0.5 * (sep1 - sep2));
  double dist2 = (0.5 * dist) + (0.5 * (sep2 - sep1));
  if(rapidNJ::negative_branches) {
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
  // handle redundancy
  newRowIndex = redundantToMasterMap[min1];
  obsoleteRowIndex = redundantToMasterMap[min2];
  if(redundantCount[redundantToMasterMap[min1]] != 1){
    // use one of the redundant rows as new row
    newRowIndex = redundantMap[redundantToMasterMap[min1]].front();
    redundantMap[redundantToMasterMap[min1]].pop_front();
    redundantCount[redundantToMasterMap[min1]]--;
    // this works because redundant rows have not been joined before.	
    indexToId[newRowIndex] = newRowIndex;
  }
  if(redundantCount[redundantToMasterMap[min2]] != 1){
    // remove one of the redundant rows
    obsoleteRowIndex = redundantMap[redundantToMasterMap[min2]].front();
    redundantMap[redundantToMasterMap[min2]].pop_front();
    redundantCount[redundantToMasterMap[min2]]--;
    // this works because redundant rows have not been joined before.	
    indexToId[obsoleteRowIndex] = obsoleteRowIndex;
  }
  mytree->addInternalNode(dist1, dist2, indexToId[newRowIndex], indexToId[obsoleteRowIndex]);
  pb->setProgress((matrixSize - clusterCount) / (double) matrixSize);
#if TESTSORT
  test->step(newRowIndex,obsoleteRowIndex,global_min);
#endif
}

int rapidNJ::next_id(){
  return currentId++;
}

void rapidNJ::findMin() {
  global_min = FLT_MAX;  
  int startidx = 0;
  if(indexToId[min_row_cache1] != -1){
    startidx = min_row_cache1;
  } else {
    startidx = min_row_cache2;
  }
  int i = startidx;
  // start at the next best row found in the previous iteration, and hope the minimum value found in this row
  // is a good approximation for the global minimum. If the row has been joined and deleted, we just take the next row
  do {
    if(indexToId[i] == -1){
      // inactive row						
    } else if(garbage_flags[i] == 1){
      int rowsize = row_lengths[i];
      findRowMinGarbage(i, rowsize);
      garbage_flags[i] = 0;
    } else {		
      int rowsize = row_lengths[i];
      int dead_count = findRowMin(i, rowsize);
      if(dead_count > 30){			
        // mark row for gabagecollection
        garbage_flags[i] = 1;
      }
    }
    i++;		
    i = i % matrixSize;
  } while (i != startidx);
  //cout << "JOINING min1=" << min1 << " min2=" << min2 << " global_min=" << global_min << " clusterCount=" << clusterCount << "\n";
  //cout << global_min << endl;   
}

/* Search for the global minimum in row i */
int rapidNJ::findRowMin(int i, int rowsize) {
  int dead_count = 0;
  distType row_sep = separations[i];
  //distType row_max_sep = row_sep + max_separation;
  distType row_max_sep = row_sep + maxRowSeparations[i];
  int j;
  cluster_pair* row = cluster_data[i];
  //Check for redundant join
  if(redundantCount[i] > 1){
    searchRedundantEntries(i);    
  }
  for(j = 0; j < rowsize; j++){
    cluster_pair pair = row[j];
    // check if pair is active	
    if(idToIndex[pair.id] != -1) {
      // check if we have looked at enough elements      
      if(pair.distance - row_max_sep >= global_min){		
        break;
      }
      // calculate the value we're optimizing over
      int other_cluster_index = idToIndex[pair.id];
      double value = pair.distance - separations[redundantToMasterMap[other_cluster_index]] - row_sep;
	    if(value < global_min){
        // cache last minimum row
        if(i != min1){
          min_row_cache2 = min_row_cache1;
          min_row_cache1 = min1;
        }
        global_min = value;
        min1 = i;
        min2 = other_cluster_index;
      }	
    } else {
      dead_count++;
    }    
  }
  return dead_count;
}

/* search row i for the global minimum and preform garbage collection of dead pairs in the row at the same time by shifting 
* pairs which is alive to the left in the arrays, thereby erasing the dead pairs.
*/
void rapidNJ::findRowMinGarbage(int i, int rowsize) {	
  distType row_sep = separations[i];
  distType row_max_sep = row_sep + maxRowSeparations[i];
  //distType row_max_sep = row_sep + max_separation;
  int j;
  int buf_index = 0;
  cluster_pair* row = cluster_data[i];	
  //Check for redundant join
  if(redundantCount[i] > 1){
    searchRedundantEntries(i);    
  }
  for(j = 0; j < rowsize; j++){
    cluster_pair pair = row[j];
    // check if pair is active
    if(idToIndex[pair.id] != -1) {
      // check if we have looked at enough elements		  	
      if(pair.distance - row_max_sep >= global_min) {							
        break;
      }			
      // calculate the value we're optimizing over
      int other_cluster_index = idToIndex[pair.id];
      double value = pair.distance - separations[redundantToMasterMap[other_cluster_index]] - row_sep;
      if(value < global_min){
        if(i != min1){
          min_row_cache2 = min_row_cache1;
          min_row_cache1 = min1;
        }
        global_min = value;
        min1 = i;
        min2 = other_cluster_index;
      }
      row[buf_index] = pair;
      buf_index++;
    }
  }
  // continue garbage collecting
  for(;j < rowsize; j++){
    cluster_pair pair = row[j];
    if(idToIndex[pair.id] != -1) {
      row[buf_index] = pair;
      buf_index++;
    }
  }
  row_lengths[i] = buf_index;
}


void rapidNJ::searchRedundantEntries(int i){
  double value =  -2 * separations[i];
  if(value < global_min){
    if(i != min1){
      min_row_cache2 = min_row_cache1;
      min_row_cache1 = min1;
    }
    global_min = value;
    min2 = redundantMap[i].front();
    min1 = i;
  }
}

rapidNJ::~rapidNJ(void)
{
  delete[] separationsums;
  delete[] separations;
  delete[] idToIndex;
  delete[] indexToId;
  for(int i = 0; i < matrixSize; i++){
    delete[] matrix[i];
    delete[] cluster_data[i];
  }
  delete[] matrix;
  delete[] cluster_data;
  delete reader;
  delete[] row_lengths;
  delete[] garbage_flags;
  delete[] redundantMap;
  delete[] redundantCount;
  delete[] redundantToMasterMap;
  delete[] maxRowSeparations;
}
