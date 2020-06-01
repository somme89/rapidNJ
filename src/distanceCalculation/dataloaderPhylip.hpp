#ifndef DATALOADER_PHILIP_H
#define DATALOADER_PHILIP_H

#include "../stdinclude.h"
#include "dataloader.hpp"
#include <string>
#include <fstream>
#include "dnaBitString.hpp"

class dataloaderPhylip : public dataloader {

public:
  dataloaderPhylip();
  ~dataloaderPhylip();

  void load(string filename);
  unsigned int** getBitStrings();
  unsigned int** getGapFilters();
  unsigned int getSequenceCount();
  unsigned int getSequenceLength();
  unsigned int getBitStringsCount();
  vector<string>* getSequenceNames();
  vector<char*>* getSequences();
  void setSequences(vector<char*>* val);

  vector<std::string>* sequenceNames;
  int seqLength;
  int seqCount;
  unsigned int** bitStrings;
  unsigned int** gapFilters;
  char** buffers;
  int* bufferSizes;
  int bitStringsCount;
  vector<char*> sequences;
  int* sequencesIdx;  
  static const int BUFFER_SIZE = 16384;
  
private:
  int* dnaBitStringIdx;
  bool parsedNames;
  void readHeader(std::ifstream &is);
  void parseName(char* input, int seq);
  inline void parseBuffers();
  inline void parseBuffersWithGaps();
  inline void parseBuffersProtein();
  void createBuffers();
  void createBitStrings();  
  void createDataStructuresNaive();
  void loadFastDist(string filename);
  void loadNormal(string filename);
  int counter;
};

#endif
