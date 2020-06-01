#ifndef DISTMATRIXREADER_H
#define DISTMATRIXREADER_H
#include "polytree.h"
#include "stdinclude.h"
#include <fstream>
#include <vector>
#include "diskMatrix.h"

using namespace std;

class distMatrixReader {
  
 public:
   /*Constructer for reading a distance matrix from disk*/
   distMatrixReader(bool verbose, std::string filename, int matrixSize, bool halfMatrix);
  
  /*Constructor for distance matrices computed from sequence data*/
  distMatrixReader(bool verbose, int matrixSize, bool halfMatrix, vector<string>* sequenceNames, distType** matrix);
  ~distMatrixReader(void);
  void read_data(diskMatrix* dm);
  distType** getMatrix();
  vector<string>* getSequenceNames();
  string getFileName();
  void initializeData();  

  void printMatrix(){
    for(int i = 0; i < matrixSize; i++){
      int rowSize = matrixSize;
      if(halfMatrix){
        rowSize = i+1;
      }
      for(int j = 0; j < rowSize; j++){
        cout << matrix[i][j] << "\t";
      }
      cout << endl;
    }
  }  


 private:
  distType** matrix;
  //  node** nodes;
  int matrixSize;
  int curRow;
  int curCol;
  bool verbose;
  string filename;
  static const int bufsize = 65536;
  char* buf;
  bool halfMatrix;
  vector<string>* sequenceNames;
  void createDatastructures();
  inline void parseData(char* input);	
};

#endif
