#include "stdinclude.h"
#include "diskMatrix.h"
#include <assert.h>
#include <sched.h>

using namespace std;

/* A distance matrix stored in external memory. The matrix is stored in files of size 1 GB.
*/

diskMatrix::diskMatrix(std::string data_dir, int size) {
  diskMatrix::size = size;
  diskMatrix::data_dir = data_dir;
  rowToFileId = new int[size];
  rowToFileIndex = new int[size];
  rowsPerFile = MAX_FILE_SIZE / (size * ((int)sizeof(distType)));
  if (rowsPerFile == 0) rowsPerFile = 1;
  numberOfFiles = size / rowsPerFile;
  if(size % rowsPerFile != 0.0){
    numberOfFiles++;
  }

  fstreams = new fstream[numberOfFiles];
  file_names = new string[numberOfFiles];
  int rowSize = size * sizeof(distType);
  int fileId = 0;
  int fileRowCount = 0;
  int totalRowCount = 0;
  
  //cout << "data_dir " << data_dir << endl;
  for(int i = 0; i < numberOfFiles; i++){
    string fname = getTempFile(data_dir);
    fstreams[fileId].open(fname.data(),fstream::out | fstream::in | fstream::binary | fstream::trunc);
    file_names[fileId] = fname;
    if(!fstreams[fileId].is_open()){      
      cerr << "Could not open file " << fname << "\n";
      cerr << "Check if the cache directory is accessible\n";
      exit(1);
    }	
    for(int j = 0; j < rowsPerFile; j++){
      rowToFileIndex[totalRowCount] = fileRowCount * rowSize;
      rowToFileId[totalRowCount] = fileId;
      if(totalRowCount == size-1){
        return;
      }
      fileRowCount++;
      totalRowCount++;
    }
    fileRowCount = 0;
    fileId++;
  }  
}

string diskMatrix::getTempFile(string directory) {
  string fname;
  if(directory == "") {
    if(getenv("TMPDIR") != NULL) {
      fname = getenv("TMPDIR");
    } else if (getenv("TMP") != NULL) {
      fname = getenv("TMP");
    } else if (getenv("TEMP") != NULL) {
      fname = getenv("TEMP");
    } else {
      fname = P_tmpdir;
    }
  } else {
    fname = directory;
  }
  char* buf = NULL;
#ifdef __WINDOWS__
  buf = new char[L_tmpnam];
  tmpnam(buf);
  fname += buf;
#else
  string s = fname + "/XXXXXX";
  buf = new char [s.size()+1];
  strcpy (buf, s.c_str());
  int ret = mkstemp(buf);
  if(ret != -1) {
    close(ret);  
  }
  fname = buf;
#endif
  delete[] buf;
  return fname;
}


void diskMatrix::writeEntry(distType value, int row, int col){
  int fileId = rowToFileId[row];  
  fstreams[fileId].seekp(rowToFileIndex[row] + col * sizeof(distType));
  fstreams[fileId].write((char*)&value, sizeof(distType));
  fstreams[fileId].flush();  
}

void diskMatrix::writeArray(distType *values, int row, int size){
  int fileId = rowToFileId[row];
  fstreams[fileId].seekp(rowToFileIndex[row]);
  fstreams[fileId].write((char*)values, sizeof(distType) * size);  
  fstreams[fileId].flush();  
}

distType diskMatrix::readEntry(int row, int col){
  char char_buf[sizeof(distType)];
  distType distType_buf[1];
  int fileId = rowToFileId[row];
  fstreams[fileId].seekg(rowToFileIndex[row] + col * sizeof(distType));
  fstreams[fileId].read(char_buf, sizeof(distType));
  memcpy(distType_buf,char_buf,sizeof(distType));
  return distType_buf[0];
}

void diskMatrix::readArray(distType *values, int row, int size){  
  int fileId = rowToFileId[row];
  fstreams[fileId].seekg(rowToFileIndex[row]);
  fstreams[fileId].read((char*)values, sizeof(distType) * size);
}

void diskMatrix::writeArrayNewSize(distType *values, int newRow, int newSize){  
  int fileId = rowToFileId[newRow];  
  int newIndex = (newRow % rowsPerFile) * newSize * sizeof(distType);
  fstreams[fileId].seekp(newIndex);  
  rowToFileIndex[newRow] = newIndex;
  fstreams[fileId].write((char*)values, sizeof(distType) * newSize);  
  fstreams[fileId].flush();
}

void diskMatrix::setNewSize(int newSize){
  size = newSize;
}

void diskMatrix::updateRowIndex(int newRow, int newSize){  
  rowToFileIndex[newRow] = (newRow % rowsPerFile) * newSize * sizeof(distType);;
}

void diskMatrix::flush(){
  for(int i = 0; i < numberOfFiles; i++){
    fstreams[i].flush();
  }
}

diskMatrix::~diskMatrix(){  
  for(int i = 0; i < numberOfFiles; i++){
    fstreams[i].close();
    remove(file_names[i].c_str());
  }
  delete[] file_names;
  delete[] fstreams;
  delete[] rowToFileId;
  delete[] rowToFileIndex;
}
