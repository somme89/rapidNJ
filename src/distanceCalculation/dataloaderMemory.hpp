#ifndef DATALOADER_MEMORY_H
#define DATALOADER_MEMORY_H

#include "stdinclude.h"
#include "dataloader.hpp"

/*Wrapper class used for bootstrapping where the dataset is kept in memory */
class dataloaderMemory : public dataloader{

public:

  dataloaderMemory(unsigned int sequenceLength, unsigned int sequenceCount, vector<string>* sequenceNames, InputType type){
    dataloaderMemory::sequenceLength = sequenceLength;
    dataloaderMemory::sequenceCount = sequenceCount;
    dataloaderMemory::sequenceNames = sequenceNames;
    dataloader::type = type;    
  }

  void initialize(vector<char*>* sequences){
    dataloaderMemory::sequences = sequences;
  }

  void initialize(unsigned int** bitStrings, unsigned int** gapFilters, unsigned int bitStringsCount, int paddingLength){
    dataloaderMemory::bitStrings = bitStrings;
    dataloaderMemory::gapFilters = gapFilters;
    dataloaderMemory::bitStringsCount = bitStringsCount;
    dataloaderMemory::paddingLength = paddingLength;
  }

  unsigned int** getBitStrings(){
    return bitStrings;
  }

  unsigned int** getGapFilters(){
    return gapFilters;
  }

  unsigned int getSequenceCount(){
    return sequenceCount;
  }

  unsigned int getSequenceLength(){
    return sequenceLength;
  }
  
  unsigned int getBitStringsCount(){
    return bitStringsCount;
  }
  
  vector<string>* getSequenceNames(){
    return sequenceNames;
  }

  vector<char*>* getSequences(){
    return sequences;
  }

  void setSequences(vector<char*>* val){
    sequences = val;
  }

  ~dataloaderMemory(){}

private:
  unsigned int** bitStrings;
  unsigned int** gapFilters;
  vector<char*>* sequences;
  vector<string>* sequenceNames;
  unsigned int sequenceLength;
  unsigned int sequenceCount;
  int bitStringsCount;
  int paddingLength;
};

#endif

