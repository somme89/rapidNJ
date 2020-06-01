/*An implementation of RapidNJ with disk caching*/

#include "stdinclude.h"
#include "rapidNJDisk.h"
#include "float.h"
#include <algorithm>

#if TESTDISK
#include "distMatrixReader.hpp"
#include "testNJ.h"
testNJ* test;
int* translate;
#endif

using namespace std;

/* initialize datastructures which doesn't require processing of any data */
rapidNJDisk::rapidNJDisk(rdDataInitialiser* reader, bool verbose, bool negative_branches, ProgressBar* pb) {
  rapidNJDisk::verbose = verbose;
  rapidNJDisk::negative_branches = negative_branches;
  rapidNJDisk::pb = pb;
  matrixSize = reader->getSize();
  matrixSizeCache = matrixSize;
  possibleRedundantRows = reader->getPossibleRedundantRows();
  clusterCount = matrixSize;  
  currentId = matrixSize;
  newRowsStartId = matrixSize;
  min_row_cache1 = 0;
  min_row_cache2 = 0;
  max_separation = FLT_MIN;       
  maxRowSeparations = new distType[matrixSize];
  min1 = 0;
  min2 = 0;
  updateCount = 0;  
  dataStructureSize = reader->getDataStructureSize();
  cluster_data = reader->getMatrix();
  mytree = reader->getTree();
  sorted_row_lengths = reader->getSortedRowLengths();
  separationsums = reader->getSeparationSums();
  disk_matrix = reader->getDiskMatrix();
  column_cache = new distType*[dataStructureSize];
  for(int i = 0; i < dataStructureSize; i++){
    column_cache[i] = new distType[matrixSize];
  }  
  separations = new distType[matrixSize];
  idToIndex = new int[2*matrixSize];
  indexToId = new int[matrixSize];
  update_flags = new int[matrixSize];
  rowBuffer1 = new distType[matrixSize];
  rowBuffer2 = new distType[matrixSize];
  cluster_buffer = new cluster_pair[matrixSize];  
  redundantCount = new int[matrixSize];
  redundantMap = new list<int>[matrixSize];
  redundantToMasterMap = new int[matrixSize];
  activeRowFlags = new short[matrixSize];

#if TESTDISK
  translate = new int[matrixSize];      
  distMatrixReader* reader2 = new distMatrixReader(false,reader->getFileName(), matrixSize, false);
  reader2->read_data();
  test = new testNJ(reader2, matrixSize);
#endif
}

/* starts the algorithm */
polytree* rapidNJDisk::run(){  
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
    if(activeRowFlags[i]){
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
  double distance;
  if(update_flags[index1] > update_flags[index2]){
    distance = column_cache[update_flags[index1]-1][index2];
  } else {
    distance = disk_matrix->readEntry(index2,index1);
  }  
#if TESTDISK
  exit(0);
#endif
  mytree->set_serialization_indices(indexToId[index1],indexToId[index2], distance);
  pb->finish();
  return mytree;
}

/* Initialize the datastructure which requires some processing */
void rapidNJDisk::initialize(){	
  for(int i = 0; i < matrixSize; i++){
#if TESTDISK
    translate[i] = i;
#endif
    update_flags[i] = 0;
    indexToId[i] = i;
    activeRowFlags[i] = 1;
    idToIndex[i] = i;
    distType separation = separationsums[i] / (clusterCount -2);
    separations[i] = separation;
    // set the maxRowSeparation here to avoid getting the rows own separation
    maxRowSeparations[i] = max_separation;
    if(separation > max_separation){
      max_separation = separation;      
    }
  }
  // handle redundancy. We iterate through the rows backwards so the longest sorted rows become master rows
  for(int i = matrixSize-1; i >= 0; i--){
    // must be set to one because the values in redundantCount is used as a multiplier in the update-phase
    redundantCount[i] = 1;    
    if(!activeRowFlags[i]){
      // this row has allready been identified as redundant
      continue;
    }
    redundantToMasterMap[i] = i;
    if(!possibleRedundantRows[i]){
      //this row has not been identified as a possible redundant row
      continue;
    }
    disk_matrix->readArray(rowBuffer1,i,matrixSize);    
    for(int j = i-1; j >= 0; j--){
      if(rowBuffer1[j] == 0 && i!=j){
        handleRedundancy(i,j);
      }
    }
  }
  for(int i = 0; i < matrixSize; i++){    
    if(!activeRowFlags[i]){
      continue;
    }
    disk_matrix->readArray(rowBuffer1,i,matrixSize);
    /*    for(int j = 0; j < matrixSize; j++){
      cout << rowBuffer1[j] << " ";
    }
    cout << endl;*/
    int idx = 0;
    for(int j = 0; j < i; j++){
      if(activeRowFlags[j]){
        cluster_buffer[idx].distance = rowBuffer1[j];
        cluster_buffer[idx].id = j;
        idx++;
      }
    }
    sort(&cluster_buffer[0],&cluster_buffer[idx]);
    int sortedRowLength = min(idx, dataStructureSize);
    memcpy(cluster_data[i], &cluster_buffer[0], sortedRowLength * sizeof(cluster_pair));
    sorted_row_lengths[i] = sortedRowLength;
  }
}

void rapidNJDisk::handleRedundancy(int i, int j){
  disk_matrix->readArray(rowBuffer2,j,matrixSize);
  for(int k = matrixSize-1; k >= 0; k--){
    if(rowBuffer1[k] != rowBuffer2[k]){
      // not redundant
      return;
    }
  }
  activeRowFlags[j] = 0;
  redundantToMasterMap[j] = i;
  redundantMap[i].push_front(j);
  redundantCount[i]++;
  //  cout << "REDUNDANT: " << i << " <- " << j << endl;
}

/* updates datastructures after 2 clusters has been joined */
void rapidNJDisk::updateData(){
  distType newSeparationsum = 0;
  distType mutualDistance = min_pair.distance;
  int master1 = redundantToMasterMap[min1];
  int master2 = redundantToMasterMap[min2];
  if(update_flags[master1] != 0){
    memcpy(rowBuffer1,column_cache[update_flags[master1]-1],matrixSize * sizeof(distType));
  } else {
    disk_matrix->readArray(rowBuffer1,master1,matrixSize);
  }
  if(update_flags[master2] != 0){
    memcpy(rowBuffer2,column_cache[update_flags[master2]-1],matrixSize * sizeof(distType));
  } else {
    disk_matrix->readArray(rowBuffer2,master2,matrixSize);
  }
  max_separation = FLT_MIN;

  for(int curId = 0; curId < currentId; curId++){
    int i = idToIndex[curId];    
    if(i == -1){
      continue;
    } 
    distType val1, val2, dist, separation;
    if(i == newRowIndex || i == obsoleteRowIndex || !activeRowFlags[i]){
      rowBuffer1[i] = 0;
    } else {
      if(update_flags[i] > update_flags[master1]){
        // read from cache. update_flags[i] > update_flags[master1] means the row i has been updated after row 'master1' was updated
        val1 = column_cache[update_flags[i]-1][master1];
      } else {
        val1 = rowBuffer1[i];
      }
      if(update_flags[i] > update_flags[master2]){
        val2 = column_cache[update_flags[i]-1][master2];
      } else {
        val2 = rowBuffer2[i];
      }
      dist = (val1 + val2 - mutualDistance) / 2.0f;
      if(rapidNJDisk::negative_branches && dist < 0) {
        dist = 0;
      }
      newSeparationsum += dist * redundantCount[i];
      // update the separationsum of cluster i.
      separationsums[i] += (dist - val1 - val2);
      separation = separationsums[i] / (clusterCount - 2);
      separations[i] = separation;
      //Set maxRowSeparation here to avoid getting the rows own separation as max
      maxRowSeparations[i] = max_separation;
      // update the maximum separation
      if(separation > max_separation){
        max_separation = separation;
      } 
      //new row and column entries      
      rowBuffer1[i] = dist;
    }
  }
  disk_matrix->writeArray(rowBuffer1,newRowIndex,matrixSize);
  memcpy(column_cache[updateCount], rowBuffer1, matrixSize * sizeof(distType));
  updateCount++;
  update_flags[newRowIndex] = updateCount;
  separationsums[newRowIndex] = newSeparationsum;
  separations[newRowIndex] = newSeparationsum / (clusterCount - 2);  
  maxRowSeparations[newRowIndex] = max_separation; 

  // delete obsolete data
  delete[] cluster_data[obsoleteRowIndex];
  cluster_data[obsoleteRowIndex] = NULL;
  idToIndex[indexToId[newRowIndex]] = -1;
  idToIndex[indexToId[obsoleteRowIndex]] = -1;
  indexToId[obsoleteRowIndex] = -1;
  activeRowFlags[obsoleteRowIndex] = 0;
  separationsums[obsoleteRowIndex] = 0;
  sorted_row_lengths[obsoleteRowIndex] = 0;
  redundantToMasterMap[newRowIndex] = newRowIndex;
  redundantToMasterMap[obsoleteRowIndex] = obsoleteRowIndex;


  // create new cluster data row
  int nextId = next_id();
  indexToId[newRowIndex] = nextId;
  idToIndex[nextId] = newRowIndex;
  int idx = 0;

  for(int i = 0; i < matrixSize; i++) {
    if(i != newRowIndex && activeRowFlags[i]){
      cluster_buffer[idx].distance = rowBuffer1[i];			
      cluster_buffer[idx].id = indexToId[i];
      idx++;
    }
  } 
  sort(cluster_buffer,&cluster_buffer[idx]);  
  int clusterRowSize = min(idx,dataStructureSize);
  memcpy(cluster_data[newRowIndex],cluster_buffer,clusterRowSize*sizeof(cluster_pair));
  sorted_row_lengths[newRowIndex] = clusterRowSize;
  if(updateCount == dataStructureSize){
    //flush column cache to HDD
    flushUpdateCache();
  }
  pb->setProgress((newRowsStartId - clusterCount) / (double) newRowsStartId);    
}


void rapidNJDisk::flushUpdateCache(){    
  if(verbose){
    cerr << "flushing update cache...\n";
  }
  // Space have been free'd. Increase the size of datastructures if there's enough free space.
  int deletedEntriesSize = (matrixSizeCache - clusterCount) * dataStructureSize * sizeof(cluster_pair);
  int newDataStructureSize = dataStructureSize + (deletedEntriesSize / ((sizeof(distType) + sizeof(cluster_pair)) * clusterCount));
  newDataStructureSize = min(clusterCount,newDataStructureSize);
  if(!(newDataStructureSize >= dataStructureSize * (DATA_STRUCTURE_SIZE_TRESSHOLD / 100.0 + 1))){
    // do not increase the size of data structures
    newDataStructureSize = -1;
  } 

  int currentRow = 0;
  for(int i = 0; i < matrixSize; i++){
    if(activeRowFlags[i]){
      //update distance matrix
      if(update_flags[i] != 0){
        memcpy(rowBuffer1,column_cache[update_flags[i]-1], matrixSize * sizeof(distType));
      } else {	
        disk_matrix->readArray(rowBuffer1,i,matrixSize);
      }
      rebuildDataStructures(currentRow, i, newDataStructureSize);
      currentRow++;
    } else if(!activeRowFlags[i] && redundantToMasterMap[i] != i){
      // just leave this row as it is. Redundant rows must be preserved in case they are needed for joining
      disk_matrix->updateRowIndex(currentRow, clusterCount);
      if(currentRow != i){
        cluster_data[currentRow] = cluster_data[i];
        cluster_data[i] = NULL;
      }
      if(newDataStructureSize != -1){
        //resize the sorted row for future use
        delete[] cluster_data[currentRow];
        cluster_data[currentRow] = new cluster_pair[newDataStructureSize];
      }
      sorted_row_lengths[currentRow] = sorted_row_lengths[i];
      currentRow++;
    }
  }
  //update various datastructures to match the new matrix
  int currentIdx = 0;
  for(int i = 0; i < matrixSize; i++) {
    update_flags[i] = 0;
    if(activeRowFlags[i] || redundantToMasterMap[i] != i){
      // row i is either an active row or a redundant row
      separationsums[currentIdx] = separationsums[i];
      separations[currentIdx] = separations[i];      
      redundantMap[currentIdx] = redundantMap[i];
      idToIndex[indexToId[i]] = currentIdx;
      indexToId[currentIdx] = indexToId[i];
      maxRowSeparations[currentIdx] = maxRowSeparations[i];

      activeRowFlags[currentIdx] = activeRowFlags[i];
#if TESTDISK
      translate[currentIdx] = translate[i];
#endif
      if(currentIdx != i){
        activeRowFlags[i] = 0;
      }
      //make sure all the redundant clusters who has this row as master, points to the correct row.
      for(list<int>::iterator it = redundantMap[currentIdx].begin(); it !=  redundantMap[currentIdx].end(); it++){
        redundantToMasterMap[*it] = currentIdx;
      }
      if(redundantToMasterMap[i] != i) {
        //a redundant row. Update the redundant map of the master row
        for(list<int>::iterator it = redundantMap[redundantToMasterMap[i]].begin(); it !=  redundantMap[redundantToMasterMap[i]].end(); it++){
          // locate the correct variable
          if(*it == i){
            *it = currentIdx;
          }
        }
      } else {	
        redundantToMasterMap[currentIdx] = currentIdx;
      }
      redundantCount[currentIdx] = redundantCount[i];            
      currentIdx++;
    }
  }

  if(newDataStructureSize != -1){
    // delete the old cache    
    for(int i = 0; i < dataStructureSize; i++){
      delete[] column_cache[i];
    }
    delete[] column_cache;

    //create a new column cache
    column_cache = new distType*[newDataStructureSize];
    for(int i = 0; i < newDataStructureSize; i++){
      column_cache[i] = new distType[clusterCount];
    }

    // set the new dataStructureSize variable
    dataStructureSize = newDataStructureSize;
    matrixSizeCache = matrixSize;
  }

  matrixSize = clusterCount;
  updateCount = 0;
  min_row_cache1 = 0;
  min_row_cache2 = 0;
  min1 = 0;
  min2 = 0;  
  if(verbose) {
    cerr << "finished flushing. New matrix size: " << clusterCount << ". New dataStructureSize: " << dataStructureSize << "\n";
  }
}

void rapidNJDisk::rebuildDataStructures(int newRowPos, int oldRowPos, int newDataStructureSize){
  int currentIdx = 0;
  int clusterBufferIdx = 0;
  //  int clusterBufferIdx = 0;
  for(int j = 0; j < matrixSize; j++){
    if(activeRowFlags[j] || redundantToMasterMap[j] != j){
      // copy this entry to the new distance matrix
      if(update_flags[j] > update_flags[oldRowPos]){
        rowBuffer1[currentIdx] = column_cache[update_flags[j]-1][oldRowPos];	
      } else {
        rowBuffer1[currentIdx] = rowBuffer1[j];	
      }
      if(newDataStructureSize != -1 && j != oldRowPos && redundantToMasterMap[j] == j){	
        cluster_buffer[clusterBufferIdx].distance = rowBuffer1[currentIdx];
        cluster_buffer[clusterBufferIdx].id = indexToId[j];	
        clusterBufferIdx++;
      }
      currentIdx++;
    }
  }
  disk_matrix->writeArrayNewSize(rowBuffer1,newRowPos,clusterCount);
  if(newRowPos != oldRowPos){
    cluster_data[newRowPos] = cluster_data[oldRowPos];
    cluster_data[oldRowPos] = NULL;
    sorted_row_lengths[newRowPos] = sorted_row_lengths[oldRowPos];
  } 

  if(newDataStructureSize != -1){
    //resize the array
    delete[] cluster_data[newRowPos];
    cluster_data[newRowPos] = new cluster_pair[newDataStructureSize];
    if(newRowsStartId <= indexToId[oldRowPos]){
      //A joined row. Sort the entire row
      sort(cluster_buffer,&cluster_buffer[clusterBufferIdx]);      
      int clusterRowSize = min(clusterBufferIdx,newDataStructureSize);
      memcpy(cluster_data[newRowPos],cluster_buffer,clusterRowSize*sizeof(cluster_pair));
      sorted_row_lengths[newRowPos] = clusterRowSize;
    } else {
      //An original row.
      //Try to limit the size of this row. indexToId[oldRowPos] = the original row index
      int deletedRows = indexToId[oldRowPos] - newRowPos;
      int upperSizeBound = indexToId[oldRowPos] - deletedRows;
      int sortSize = min(upperSizeBound, clusterBufferIdx);
      sort(cluster_buffer,&cluster_buffer[sortSize]);      
      int clusterRowSize = min(newDataStructureSize, sortSize);
      clusterRowSize = min(clusterBufferIdx, clusterRowSize);
      memcpy(cluster_data[newRowPos],cluster_buffer,clusterRowSize*sizeof(cluster_pair));
      sorted_row_lengths[newRowPos] = clusterRowSize;
    } 
  }
}

/* merge two clusters */
void rapidNJDisk::mergeMinNodes(){  
  // calculate distances
  double dist = min_pair.distance;
  double sep1 = separationsums[redundantToMasterMap[min1]] / (clusterCount - 2);
  double sep2 = separationsums[redundantToMasterMap[min2]] / (clusterCount - 2);
  double dist1 = (0.5 * dist) + (0.5 * (sep1 - sep2));
  double dist2 = (0.5 * dist) + (0.5 * (sep2 - sep1));
  if(rapidNJDisk::negative_branches) {
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
    // activate this row for consistency 
    activeRowFlags[newRowIndex] = 1;
  }
  if(redundantCount[redundantToMasterMap[min2]] != 1){
    // remove one of the redundant rows
    obsoleteRowIndex = redundantMap[redundantToMasterMap[min2]].front();
    redundantMap[redundantToMasterMap[min2]].pop_front();
    redundantCount[redundantToMasterMap[min2]]--;
    // activate this row for consistency 
    activeRowFlags[obsoleteRowIndex] = 1;
  }
  // update tree
  mytree->addInternalNode(dist1, dist2, indexToId[newRowIndex], indexToId[obsoleteRowIndex]);
#if TESTDISK
  test->step(translate[newRowIndex],translate[obsoleteRowIndex],global_min);  
#endif
  //  cout << "RESOLVED: " << newRowIndex << " " << obsoleteRowIndex << endl;
}

int rapidNJDisk::next_id(){
  return currentId++;
}

void rapidNJDisk::findMin() {  
  global_min = FLT_MAX;  
  int startidx = 0;
  if(activeRowFlags[min_row_cache1]){
    startidx = min_row_cache1;
  } else {
    startidx = min_row_cache2;
  }  
  int i = startidx;
  // start at the next best row found in the previous iteration, and hope the minimum value found in this row
  // is a good approximation for the global minimum. If the row has been joined and deleted, we just take the next row
  do {    
    if(activeRowFlags[i]){
      findRowMin(i);
    }
    i++;		
    i = i % matrixSize;
  } while (i != startidx);  

  //  cout << "JOINING min1=" << min1 << " min2=" << min2 << " global_min=" << global_min <<  " clusterCount=" << clusterCount << "\n";
  //  cout << global_min << endl;
}

/* Search for the global minimum in row i */
void rapidNJDisk::findRowMin(int rowIdx) {
  int rowsize = sorted_row_lengths[rowIdx];  
  distType row_sep = separations[rowIdx];
  distType row_max_sep = row_sep + maxRowSeparations[rowIdx];
  cluster_pair* row = cluster_data[rowIdx];
  cluster_pair pair;
  //Check for redundant join
  if(redundantCount[rowIdx] > 1){
    searchRedundantEntries(rowIdx);    
  }
  for(int i = 0; i < rowsize; i++){
    pair = row[i];
    // check if pair is active
    if(idToIndex[pair.id] != -1) {
      // check if we have looked at enough elements
      if(pair.distance - row_max_sep >= global_min){
        return;
      }
      // calculate the value we're optimizing over
      int other_cluster_index = idToIndex[pair.id];
      distType value = pair.distance - separations[redundantToMasterMap[other_cluster_index]] - row_sep;
      if(value < global_min){
        // cache last minimum row
        if(rowIdx != min1){
          min_row_cache2 = min_row_cache1;
          min_row_cache1 = min1;
        }
        global_min = value;
        min_pair = pair;
        min1 = rowIdx;
        min2 = other_cluster_index;
      }
    } 
  }
  if(rowsize < dataStructureSize){
    return;
  }
  // now we have to search the entire row
  if(update_flags[rowIdx] != 0){
    //load the row from cache
    memcpy(rowBuffer1,column_cache[update_flags[rowIdx]-1],matrixSize * sizeof(distType));
  } else {
    // load the row from HDD
    disk_matrix->readArray(rowBuffer1,rowIdx,matrixSize);
  }
  for(int i = 0; i < matrixSize; i++){    
    if(i != rowIdx && activeRowFlags[i]) {
      distType value;
      distType distance;
      if(update_flags[i] > update_flags[rowIdx]){
        distance = column_cache[update_flags[i]-1][rowIdx];
      } else {
        distance = rowBuffer1[i];
      }      
      value = distance - row_sep - separations[i];
      if(value < global_min){
        // cache last minimum row
        if(rowIdx != min1){
          min_row_cache2 = min_row_cache1;
          min_row_cache1 = min1;
        }
        global_min = value;
        min_pair.distance = distance;
        min_pair.id = indexToId[i];
        min1 = rowIdx;
        min2 = i;
      }

    } 
  }  
  return;
}

void rapidNJDisk::searchRedundantEntries(int i){
  distType value =  -2 * separations[i];
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

rapidNJDisk::~rapidNJDisk(void){  
  delete[] separations;
  delete[] idToIndex;
  delete[] indexToId;
  delete disk_matrix;
  delete[] activeRowFlags;
  delete[] maxRowSeparations;
  delete[] update_flags;
  delete[] rowBuffer1;
  delete[] rowBuffer2;
  delete[] redundantCount;
  delete[] redundantMap;
  delete[] redundantToMasterMap;
  delete[] cluster_buffer;  
  for(int i = 0; i < dataStructureSize; i++){
    delete[] column_cache[i];
  }
  delete[] column_cache;
}

//CHECK MAXROWSEPARATIONS INITIALIZE
/*
cout << "------------------------------------------------ \n";
for(int i = 0; i < currentId; i++){
if(idToIndex[i] != -1){
int rowIdx = idToIndex[i];
if(!activeRowFlags[rowIdx]){
continue;
}
double _max = 0; 
for(int j = 0; j < sorted_row_lengths[rowIdx]; j++){
int colIdx = idToIndex[cluster_data[rowIdx][j].id];
if(colIdx == -1 || !activeRowFlags[colIdx]) {
continue;
}	
double sep = separations[colIdx];
if(rowIdx == 38){
cout << colIdx << ": " << sep << endl;
}
if(sep > _max){
_max = sep;
}	  
}
if(redundantCount[i] != 1){
if(separations[rowIdx] > _max){
_max = separations[rowIdx];
}	  

}
if(sorted_row_lengths[rowIdx] > 0 && _max != maxRowSeparations[rowIdx]){
cout << "INITIALIZE " << rowIdx << " " <<  _max << "!=" << maxRowSeparations[rowIdx] << ": " << sorted_row_lengths[rowIdx] << endl;
printf("TEST:  %.20f %.20f \n",_max,maxRowSeparations[rowIdx]);
exit(0);
}
}
}
*/

//CHECK MAXROWSEPARATIONS UPDATEDATA
/*
for(int i = 0; i < currentId; i++){
if(idToIndex[i] != -1){
int rowIdx = idToIndex[i];
if(!activeRowFlags[rowIdx]){
continue;
}
double _max = 0; 
for(int j = 0; j < sorted_row_lengths[rowIdx]; j++){
int colIdx = idToIndex[cluster_data[rowIdx][j].id];
if(colIdx == -1 || !activeRowFlags[colIdx]) {
continue;
}	
double sep = separations[colIdx];
if(sep > _max){
_max = sep;
}	  
}      
if(sorted_row_lengths[rowIdx] > 0 && _max != maxRowSeparations[rowIdx]){
cout << rowIdx << " " <<  _max << "!=" << maxRowSeparations[rowIdx] << ": " << sorted_row_lengths[rowIdx] << endl;
printf("TEST:  %.20f %.20f \n",_max,maxRowSeparations[rowIdx]);
exit(0);
}
}
}
*/
