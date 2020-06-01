#include "distMatrixReader.hpp"

distMatrixReader::distMatrixReader(bool verbose, string filename, int matrixSize, bool halfMatrix) {	
  distMatrixReader::matrixSize = matrixSize;
  curRow = 0;
  curCol = 0;
  distMatrixReader::verbose = verbose;
  distMatrixReader::filename = filename;  
  buf = new char[bufsize];
  distMatrixReader::halfMatrix = halfMatrix;
  distMatrixReader::sequenceNames = NULL;
}

distMatrixReader::distMatrixReader(bool verbose, int matrixSize, bool halfMatrix, vector<string>* sequenceNames, distType** matrix) {	
  distMatrixReader::matrixSize = matrixSize;
  curRow = 0;
  curCol = 0;  
  distMatrixReader::verbose = verbose;
  buf = new char[bufsize];
  distMatrixReader::halfMatrix = halfMatrix;
  distMatrixReader::matrix = matrix;
  distMatrixReader::sequenceNames = sequenceNames;
}

/*Used to initialise data structures when an alignment has been used as input */
void distMatrixReader::initializeData() {
  if(halfMatrix){
    for(int k = 0; k < matrixSize; k++){
      distType* newRow = new distType[k+1];
      memcpy(newRow,matrix[k], sizeof(distType)*(k+1));
      delete[] matrix[k];
      matrix[k] = newRow;      
    }   
  }
}

/* Reads data from a file. The data must be formated as a phylip matrix.*/
void distMatrixReader::read_data(diskMatrix* dm){
  char charInput[256];
  ifstream is;
  // open the file
  is.open(filename.data(),ifstream::in);
  if(!is.is_open()){
    cerr << "Could not read file: " << filename << endl;
    exit(1);
  }
  createDatastructures();

  int j = 0;
  int i = 0;
  bool firstLine = true;

  while(is.good()) {
    is.read(buf,bufsize);
    i = 0;
    if(firstLine){
      // scan over the first line containing the size  
      while(firstLine && i < is.gcount()) {			
        char c = buf[i];
        i++;
        if(c == 10){
          firstLine = false;
          break;
        }
      }
    }
    // fill matrix
    while(i < is.gcount()) {
      char c = buf[i];	
      if(c == 32 || c == 9){
        if(j == 0){
          // ignore. leading space
        } else {
          // parse the data found
          charInput[j] = '\0';
          parseData(charInput);
          j = 0;
          curCol++;	  
        }
      } else if(c > 32 && c < 127){
        // found data
        charInput[j] = c;
        j++;
      } else if(c == 10){
        // new line	
        charInput[j] = '\0';
        if(j != 0){
          // if j = 0 then we have read some trailing whitespaces.
          parseData(charInput);
          curCol++;
        }
        if(curCol < matrixSize+1) {
          cout << "Row " << curRow << " has fewer columns than the stated size of " << matrixSize << endl;
          exit(1);
        }
        j = 0;
        curCol = 0;
        curRow++;	  
        if(curRow > matrixSize){
          // ignore the rest of the data, we have read enough
          is.close();						   
          return;
        }
      } 
      i++;
    }    
  }
  if(j > 0){
    // parse last piece of data, in case of missing newline on last line
    parseData(charInput);
  }
  is.close();						   
  return;
}

/* Parse a piece of data as either text or a distType depending on the position of the data in the matrix */
inline void distMatrixReader::parseData(char* input){
  if(curCol > matrixSize){
    cout << "Row " << curRow << " has more columns than the stated size of " << matrixSize << endl;
    exit(1);
  }
  if(curCol > 0){    
    // parse as distType    
    if(halfMatrix){
      if(curRow >= curCol-1){
        matrix[curRow][curCol-1] = atof(input);
      }
    } else {
      matrix[curRow][curCol-1] = atof(input);
    }
  } else {
    // parse as string
    sequenceNames->push_back(input);
  }
}

/* Create matrix datastructure and tree datastructure*/
void distMatrixReader::createDatastructures(){
  if(halfMatrix){
    matrix = new distType*[matrixSize];
    for(int k = 0; k < matrixSize; k++){
      matrix[k] = new distType[k+1];
    }
  } else {
    matrix = new distType*[matrixSize];
    for(int k = 0; k < matrixSize; k++){
      matrix[k] = new distType[matrixSize];
    }  
  }
  if(sequenceNames == NULL) {
    sequenceNames = new vector<string>();
  }
}

distType** distMatrixReader::getMatrix(){
  return matrix;
}

vector<string>* distMatrixReader::getSequenceNames() {
  return sequenceNames;
}

string distMatrixReader::getFileName(){
  return distMatrixReader::filename;
}

distMatrixReader::~distMatrixReader(void)
{
  delete[] buf;  
}

