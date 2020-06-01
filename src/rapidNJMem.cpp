/*An implementation of RapidNJ */

#include "stdinclude.h"
#include "rapidNJMem.hpp"
#include "float.h"
#include <algorithm>

#if TESTNJMEM
#include "distMatrixReader.hpp"
#include "testNJ.h"
testNJ* test;
#endif


using namespace std;

/* initialize datastructures which doesn't require processing of data */
rapidNJMem::rapidNJMem(distMatrixReader* reader, int matrixSize, int sortedMatrixSize, bool verbose, bool negative_branches, ProgressBar* pb) {
  rapidNJMem::matrixSize = matrixSize;
  rapidNJMem::negative_branches = negative_branches;
  rapidNJMem::pb = pb;
  distanceMatrix = reader->getMatrix();
  rapidNJMem::reader = reader;
  rapidNJMem::sortedMatrixSize = sortedMatrixSize;
  clusterCount = matrixSize;
  currentId = 0;
  min_row_cache1 = 0;
  min_row_cache2 = 0;
  max_separation = FLT_MIN;       
  min1 = 0;
  min2 = 0;  
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
  cluster_buffer = new cluster_pair[matrixSize];
  maxRowSeparations = new distType[matrixSize];
  rapidNJMem::verbose = verbose;

#if TESTNJMEM
  distMatrixReader* reader2 = new distMatrixReader(false,reader->getFileName(), matrixSize, false);
  reader2->read_data();
  test = new testNJ(reader2, matrixSize);
#endif
}

/* starts the algorithm */
polytree* rapidNJMem::run(){  
  initialize();
  while(clusterCount > 2){
    findMin();
    mergeMinNodes();
    clusterCount--;
    updateData();  
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
#if TESTNJMEM
  exit(0);
#endif
  double distance = getDist(index1,index2);
  mytree->set_serialization_indices(indexToId[index1],indexToId[index2], distance);  
  pb->finish();
  return mytree;  
}

/* Initialize datastructures which require pre-processing */
void rapidNJMem::initialize(){
  mytree = new polytree(matrixSize, reader->getSequenceNames());
  //calculate initial seperation rows
  for(int i = 0; i < matrixSize; i++){
    int nextId = next_id();
    indexToId[i] = nextId;
    idToIndex[nextId] = i;
    garbage_flags[i] = 0;
    // }
    //for(int i = 0; i < matrixSize; i++){
    distType sum = 0;
    redundantCount[i] = 1;
    cluster_data[i] = new cluster_pair[sortedMatrixSize];    
    redundantToMasterMap[i] = i;
    for(int j = 0; j < matrixSize; j++){      
      sum += getDist(i,j);
    }  
    separationsums[i] = sum;
    distType separation = sum / (clusterCount -2);
    separations[i] = separation;				
    maxRowSeparations[i] = max_separation;
    if(separation > max_separation){
      max_separation = separation;
    }    
  }
  // Handle redundancy
  for(int i = matrixSize-1; i > 0; i--){
    if(indexToId[i] == -1){
      // this row has allready been identified as redundant
      continue;
    }
    for(int j = 0; j < i; j++){
      if(getDist(i,j) == 0){
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
inline void rapidNJMem::buildSortedRow(int i){  
  //It's important that the row repesenting the redundant rows is the longest. 
  int rowLength = min(i,sortedMatrixSize); 
  int idx = 0;
  for(int j = 0; j < i; j++){
    if(indexToId[j] != -1){
      cluster_buffer[idx].distance = getDist(i,j);			
      cluster_buffer[idx].id = j;
      idx++;
    }
  }
  // sort the row
  rowLength = min(idx,sortedMatrixSize); 
  sort(&cluster_buffer[0],&cluster_buffer[idx]);
  memcpy(cluster_data[i], &cluster_buffer[0], rowLength * sizeof(cluster_pair));
  row_lengths[i] = rowLength;
}

inline void rapidNJMem::handleRedundancy(int i, int j){  
  for(int k = 0; k < matrixSize; k++){
    if(getDist(i,k) != getDist(j,k)){
      // not redundant
      return;
    }
  }
  indexToId[j] = -1;
  redundantToMasterMap[j] = i;
  redundantMap[i].push_front(j);
  redundantCount[i]++;  
}

void rapidNJMem::updateData(){
  distType newSeparationsum = 0;
  distType mutualDistance = min_pair.distance;  
  int min1_master = redundantToMasterMap[min1];
  int min2_master = redundantToMasterMap[min2];  
  max_separation = FLT_MIN;
  for(int curId = 0; curId < currentId; curId++){
    int i = idToIndex[curId];    
    if(i == -1 || indexToId[i] == -1){
      continue;
    }
    if(i == newRowIndex || i == obsoleteRowIndex){
      setDist(newRowIndex,i,0);
    } else {
      distType val1 = getDist(min1_master,i);
      distType val2 = getDist(min2_master,i);
      distType dist = ((val1 + val2 - mutualDistance) / 2.0f);
      if(rapidNJMem::negative_branches && dist < 0) {
        dist = 0;
      }
      newSeparationsum += dist * redundantCount[i];
      // update the separationsum of cluster i.
      separationsums[i] += (dist - val1 - val2);
      distType separation = separationsums[i] / (clusterCount - 2);
      separations[i] = separation;
      //Set maxRowSeparation here to avoid getting the rows own separation as max
      maxRowSeparations[i] = max_separation;
      // update the maximum separation
      if(separation > max_separation){
        max_separation = separation;
      }
      setDist(newRowIndex, i, dist);
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
  redundantToMasterMap[newRowIndex] = newRowIndex;
  redundantToMasterMap[obsoleteRowIndex] = obsoleteRowIndex;

  // create new cluster pair
  int nextId = next_id();
  indexToId[newRowIndex] = nextId;
  idToIndex[nextId] = newRowIndex;
  int idx = 0;

  for(int i = 0; i < matrixSize; i++) {
    if(i != newRowIndex && indexToId[i] != -1){
      cluster_buffer[idx].distance = getDist(newRowIndex,i);
      cluster_buffer[idx].id = indexToId[i];
      idx++;
    }
  }
  sort(&cluster_buffer[0],&cluster_buffer[idx]);
  row_lengths[newRowIndex] = min(idx,sortedMatrixSize);
  memcpy(cluster_data[newRowIndex], &cluster_buffer[0], row_lengths[newRowIndex] * sizeof(cluster_pair));
  pb->setProgress((matrixSize - clusterCount) / (double) matrixSize);
}


/* merge two clusters */
void rapidNJMem::mergeMinNodes(){  
  // calculate distances
  double dist = min_pair.distance;
  double sep1 = separationsums[redundantToMasterMap[min1]] / (clusterCount - 2);
  double sep2 = separationsums[redundantToMasterMap[min2]] / (clusterCount - 2);
  double dist1 = (0.5 * dist) + (0.5 * (sep1 - sep2));
  double dist2 = (0.5 * dist) + (0.5 * (sep2 - sep1));
  if(rapidNJMem::negative_branches) {
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
  // update tree
  mytree->addInternalNode(dist1, dist2, indexToId[newRowIndex], indexToId[obsoleteRowIndex]);
  pb->setProgress((matrixSize - clusterCount) / (double) matrixSize);
#if TESTNJMEM
  test->step(newRowIndex,obsoleteRowIndex,global_min);  
#endif

}

int rapidNJMem::next_id(){
  return currentId++;
}

void rapidNJMem::findMin() {
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
  //  cout << "JOINING min1=" << min1 << " min2=" << min2 << " global_min=" << global_min << " clusterCount=" << clusterCount << "\n";
  //  cout << global_min << endl;
}

/* Search for the global minimum in row i */
int rapidNJMem::findRowMin(int i, int rowsize) {
  int dead_count = 0;
  distType row_sep = separations[i];
  //  distType row_max_sep = row_sep + max_separation;
  distType row_max_sep = row_sep + maxRowSeparations[i];
  int j;
  cluster_pair* row = cluster_data[i];

  //Check for redundant join
  if(redundantCount[i] > 1){
    //cout << i << " searching redundant " << endl;
    searchRedundantEntries(i);    
  }
  for(j = 0; j < rowsize; j++){
    cluster_pair pair = row[j];
    // check if pair is active
    if(idToIndex[pair.id] != -1) {
      // check if we have looked at enough elements      
      if(pair.distance - row_max_sep >= global_min){	
        return dead_count;
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
        min_pair = pair;
        min1 = i;
        min2 = other_cluster_index;
      }	
    } else {
      dead_count++;
    }    
  }  
  // we need to search the distance matrix now
  for(int k = 0; k < matrixSize; k++){
    if(indexToId[k] != -1 && i != k) {
      double value = getDist(i,k) - separations[k] - row_sep;
      if(value < global_min){
        // cache last minimum row
        if(i != min1) {
          min_row_cache2 = min_row_cache1;
          min_row_cache1 = min1;
        }
        global_min = value;
        min_pair.distance = getDist(i,k);
        min_pair.id = -1;
        min1 = i;
        min2 = k;
      }	
    }
  }
  return dead_count;
}

/* search row i for the global minimum and preform garbage collection of dead pairs in the row at the same time by shifting 
* pairs which is alive to the left in the arrays, thereby erasing the dead pairs.
*/
void rapidNJMem::findRowMinGarbage(int i, int rowsize) {	
  distType row_sep = separations[i];
  //  distType row_max_sep = row_sep + max_separation;
  distType row_max_sep = row_sep + maxRowSeparations[i];
  int j;
  int buf_index = 0;
  cluster_pair* row = cluster_data[i];	
  bool reachedUpperBound = false;
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
        reachedUpperBound = true;
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
        min_pair = pair;
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

  if(!reachedUpperBound){
    //search the distance matrix
    for(int k = 0; k < matrixSize; k++){
      if(indexToId[k] != -1 && i != k) {
        double value = getDist(i,k) - separations[k] - row_sep;
        if(value < global_min){
          // cache last minimum row
          if(i != min1){
            min_row_cache2 = min_row_cache1;
            min_row_cache1 = min1;
          }
          global_min = value;
          min_pair.distance = getDist(i,k);
          min_pair.id = -1;
          min1 = i;
          min2 = k;
        }	
      }
    } 
  }
}

inline void rapidNJMem::searchRedundantEntries(int i){
  double value =  -2 * separations[i];
  if(value < global_min){
    if(i != min1){
      min_row_cache2 = min_row_cache1;
      min_row_cache1 = min1;
    }
    global_min = value;
    min2 = redundantMap[i].front();
    min1 = i;
    min_pair.distance = 0;
    min_pair.id = -1;
  }
}

rapidNJMem::~rapidNJMem(void)
{
  delete[] separationsums;
  delete[] separations;
  delete[] idToIndex;
  delete[] indexToId;
  for(int i = 0; i < matrixSize; i++){
    delete[] distanceMatrix[i];
    delete[] cluster_data[i];
  }
  delete[] distanceMatrix;
  delete[] cluster_data;
  delete reader;
  delete[] row_lengths;
  delete[] garbage_flags;
  delete[] redundantCount;
  delete[] redundantToMasterMap;
  delete[] maxRowSeparations;
  delete[] cluster_buffer;
  delete[] redundantMap;
}
