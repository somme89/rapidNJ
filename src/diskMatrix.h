#ifndef DISKMATRIX_H
#define DISKMATRIX_H

#include "stdinclude.h"
#include <string>
#include <fstream>

class diskMatrix {
  
 public:
  diskMatrix(const std::string data_dir, int size);
  ~diskMatrix(void);
  void writeArray(distType *values, int row, int size);
  void writeEntry(distType value, int row, int col);
  string getTempFile(string directory);
  distType readEntry(int row, int col);
  void readArray(distType *values, int row, int size);
  void writeArrayNewSize(distType *values, int newRow, int newSize);
  void setNewSize(int newSize);
  void updateRowIndex(int newRow, int newSize);
  void flush();
  void finishWrite();
  
 private:
  int size;
  const static int MAX_FILE_SIZE = 1073741824; //bytes (1GB)
  std::string data_dir;
  int *rowToFileId;
  int *rowToFileIndex;
  std::fstream* fstreams;
  string* file_names;
  int rowsPerFile;
  int numberOfFiles;  
  pthread_t thread;
};

#endif

