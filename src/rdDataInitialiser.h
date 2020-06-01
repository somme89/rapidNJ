#ifndef RDDATAINITIALISER_H
#define RDDATAINITIALISER_H
#include "polytree.h"
#include <string>
#include "diskMatrix.h"
#include "cluster_pair.h"

class rdDataInitialiser {
  
 public:
  rdDataInitialiser(bool verbose, int sortedMatrixSize, string cacheDir, string filename);
  rdDataInitialiser(bool verbose, int sortedMatrixSize, string cacheDir, int matrixSize);
  ~rdDataInitialiser(void);
  bool read_data();
  void initializeFromExistingMatrix(vector<string>* sequenceNames, diskMatrix* dm);
  cluster_pair** getMatrix();
  polytree* getTree();
  int getSize();
  int* getRowLengths();
  int* getSortedRowLengths();
  distType* getSeparationSums();
  diskMatrix* getDiskMatrix();
  int getDataStructureSize();  
  short* getPossibleRedundantRows();
  std::string getFileName();
  
 private:
  cluster_pair** matrix;
  polytree* mytree;  
  distType* rowBuffer;
  int *row_lengths;
  int *sorted_row_lengths;
  distType *separation_sums;
  short *possibleRedundantRows;
  diskMatrix *dm;
  int matrixSize;
  int curRow;
  int curCol;
  bool verbose;
  distType curSeparationSum;
  int sortedMatrixSize;
  distType memSize;
  std::string cacheDir;
  std::string filename;
  
  void createDatastructures();
  void parseData(char* input);	
  void processRow(bool writeToDisk);  
};

#endif
