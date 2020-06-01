#include "stdinclude.h"
#include <fstream>
#include "rdDataInitialiser.h"
#include <algorithm>
#include "cluster_pair.h"

using namespace std;

rdDataInitialiser::rdDataInitialiser(bool verbose, int sortedMatrixSize, string cacheDir, string filename){
  rdDataInitialiser::verbose = verbose;
  matrixSize = -1;
  curRow = 0;
  curCol = 0;    
  curSeparationSum = 0;
  rdDataInitialiser::sortedMatrixSize = sortedMatrixSize;
  rdDataInitialiser::cacheDir = cacheDir;
  rdDataInitialiser::filename = filename;
  rdDataInitialiser::dm = NULL;
}

rdDataInitialiser::rdDataInitialiser(bool verbose, int sortedMatrixSize, string cacheDir, int matrixSize){
  rdDataInitialiser::verbose = verbose;
  curRow = 0;
  curCol = 0;    
  curSeparationSum = 0;
  rdDataInitialiser::sortedMatrixSize = sortedMatrixSize;
  rdDataInitialiser::cacheDir = cacheDir;
  rdDataInitialiser::filename = "";
  rdDataInitialiser::dm = NULL;
  rdDataInitialiser::matrixSize = matrixSize;
}

void rdDataInitialiser::initializeFromExistingMatrix(vector<string>* sequenceNames, diskMatrix* dm) {
  rdDataInitialiser::dm = dm;
  createDatastructures();
  for(int i = 0; i < matrixSize; i++) {
    dm->readArray(rowBuffer, i, matrixSize);    
    for(int j = 0; j < matrixSize; j++){
      curSeparationSum += rowBuffer[j];
    }
    processRow(false);
    curSeparationSum = 0;
    mytree->addLeaf(sequenceNames->at(i));    
    curRow++;
  }  
  delete[] rowBuffer;
}

/* Reads data from a file. The data must be formated as a symmetrical phylip matix.*/
bool rdDataInitialiser::read_data(){
  int bufsize = 65536;
  char * buf = new char[bufsize];
  char charInput[256];
  ifstream is;
  // open the file
  is.open(filename.data(),ifstream::in);  
  if(!is.is_open()){
    return false;
  }
  int j = 0;
  // main loop. Keep reading until no more data or enough data has been read
  while(is.good()){
    is.read(buf,bufsize);
    int i = 0;
    // read the size of the matrix
    if(matrixSize == -1){
      while(i < is.gcount()) {				
        char c = buf[i];
        if(c == 10){
          charInput[j] = '\0';
          i++;		
          break;
        } else if(c >= 47 && c <= 57){
          charInput[j] = c;
          j++;
        }
        i++;
      }
      matrixSize = atoi(charInput);
      j = 0;
      createDatastructures();
    }
    // build matrix		
    bool trailingSpace = false;
    while(i < is.gcount()) {
      char c = buf[i];			
      if(c == 32 || c == 9){
        if(j == 0){
          // ignore. leading space
        } else if(curCol == matrixSize){
          trailingSpace = true;
        } else {					
          // parse the data found
          charInput[j] = '\0';
          parseData(charInput);
          j = 0;
          curCol++;
          if(curCol > matrixSize){
            //cout << curCol << " " << matrixSize << endl;
            cerr << "Matrix doesn't match the stated size \n";
            return false;
          }
        }
      } else if(c > 32 && c < 127){
        // found data
        if(trailingSpace){
          cerr << "Matrix doesn't match the stated size \n";
          return false;
        }
        charInput[j] = c;
        j++;
      } else if(c == 10){
        // new line				
        if(curCol == matrixSize){
          charInput[j] = '\0';	  
          parseData(charInput);
          j = 0;
          trailingSpace = false;
          processRow(true);
          curSeparationSum = 0;
          curCol = 0;	  
          curRow++;	  
        } else {
          cerr << "Matrix doesn't match the stated size \n";
          return false;
        }
        if(curRow == matrixSize){
          // ignore the rest of the data if any.
          return true;
        }
      } 
      i++;
    }    
  }
  if(j > 0){
    // parse last piece of data, in case of missing newline on last line
    parseData(charInput);
    processRow(true);
  }
  is.close();
  return true;
}

/* Parse a piece of data as either text or a distType depending on the position of the data in the matrix */
inline void rdDataInitialiser::parseData(char* input){
  if(curCol > 0){
    // parse as distType
    rowBuffer[curCol-1] = (distType) atof(input);
    //cout << "row " << rowBuffer[curCol-1] << endl;
    curSeparationSum += rowBuffer[curCol-1];
  } else {
    // parse as string
    mytree->addLeaf(input);
  }
}

/* Create matrix datastructure and tree datastructure*/
void rdDataInitialiser::createDatastructures(){  
  matrix = new cluster_pair*[matrixSize];
  for(int k = 0; k < matrixSize; k++){
    matrix[k] = new cluster_pair[sortedMatrixSize];
  }
  mytree = new polytree(matrixSize, NULL);
  if(dm == NULL) {
    dm = new diskMatrix(cacheDir, matrixSize);
  }

  rowBuffer = new distType[matrixSize];
  sorted_row_lengths = new int[matrixSize];
  separation_sums = new distType[matrixSize];
  possibleRedundantRows = new short[matrixSize];
}

void rdDataInitialiser::processRow(bool writeToDisk){
  // write the distance matrix row to the HDD
  if(writeToDisk) {
    dm->writeArray(rowBuffer,curRow,matrixSize);
  }
  //search for possible redundant rows
  possibleRedundantRows[curRow] = 0;
  for(int i = 0; i < matrixSize; i++){
    if(rowBuffer[i] == 0 && i != curRow){
      possibleRedundantRows[curRow] = 1;
      break;
    }
  }
  separation_sums[curRow] = curSeparationSum;
}

int rdDataInitialiser::getSize(){
  return matrixSize;
}

cluster_pair** rdDataInitialiser::getMatrix(){
  return matrix;
}

polytree* rdDataInitialiser::getTree(){
  return mytree;
}

int* rdDataInitialiser::getSortedRowLengths(){
  return sorted_row_lengths;
}

distType* rdDataInitialiser::getSeparationSums(){
  return separation_sums;
}

diskMatrix* rdDataInitialiser::getDiskMatrix(){
  return dm;
}

int rdDataInitialiser::getDataStructureSize(){
  return sortedMatrixSize;
}

short* rdDataInitialiser::getPossibleRedundantRows(){
  return possibleRedundantRows;
}

string rdDataInitialiser::getFileName(){
  return filename;
}

rdDataInitialiser::~rdDataInitialiser(void)
{  
  for(int i = 0; i < matrixSize; i++){
    delete[] matrix[i];
  }
  delete[] matrix;
  delete[] possibleRedundantRows;
  delete[] separation_sums;
  delete[] sorted_row_lengths;
}
