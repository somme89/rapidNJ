#ifndef DATALOADER_H
#define DATALOADER_H

#include <stdinclude.h>
#include <fstream>
#ifdef ENABLEGPU
#include "cutil_inline.h"
#include "gpu/constants.hpp"
#endif

using namespace std;

class dataloader {

public:
  dataloader(void);
  virtual ~dataloader();  
  void sample_sequences();  
  virtual void load(string filename){};
  virtual unsigned int** getBitStrings() = 0;
  virtual unsigned int** getGapFilters() = 0;
  virtual unsigned int getSequenceCount() = 0;
  virtual unsigned int getSequenceLength() = 0;
  virtual unsigned int getBitStringsCount() = 0;
  virtual vector<string>* getSequenceNames() = 0;
  virtual vector<char*>* getSequences() = 0;
  virtual void setSequences(vector<char*>*) = 0;
  unsigned int* bitStringsGPU;
  unsigned int* gapFiltersGPU; 
  vector<char*>* sequences;
  vector<unsigned int*>* bitStrings;
  vector<unsigned int*>* gapFilters;
  /*Used to store the original data when performing bootstrapping*/
  vector<char*>* sequences_original;
  vector<unsigned int*>* bitStrings_original;
  vector<unsigned int*>* gapFilters_original;

  void initGPUData(){
#ifdef ENABLEGPU
    unsigned int elmPrInt = 4;
    if(type == DNA){
      elmPrInt = 16;
    }
    bsStride = getSequenceLength() / elmPrInt;
    if(getSequenceLength() % elmPrInt != 0) {
      bsStride++;
    }
    bsStride = (bsStride + threads_pr_block-1) & ~(threads_pr_block-1);
    unsigned int memSize = bsStride * getSequenceCount() * sizeof(unsigned int);
    cudaHostAlloc((void**)&bitStringsGPU,memSize,cudaHostAllocDefault);
    cudaHostAlloc((void**)&gapFiltersGPU,memSize,cudaHostAllocDefault);
    for(unsigned int i = 0; i < bsStride * getSequenceCount(); i++){
      bitStringsGPU[i] = 0;
      gapFiltersGPU[i] = 0;
    }
    unsigned int** bitStrings = getBitStrings();
    unsigned int** gapFilters = getGapFilters();
    unsigned int dataSize = getBitStringsCount() * 4;
    dataSize = min(dataSize,bsStride);
    for(unsigned int i = 0; i < getSequenceCount(); i++) {
      memcpy(&bitStringsGPU[i*bsStride], bitStrings[i], dataSize*sizeof(unsigned int));
      if(type == DNA){	
	memcpy(&gapFiltersGPU[i*bsStride], gapFilters[i], dataSize*sizeof(unsigned int));
      }
    }
#else
    cout << "GPU distance computation is not available" << endl;
    exit(1);
#endif
  };

  void readData(string filename, bool fastdist, bool verbose, InputType type, bool gpu) {
    dataloader::verbose = verbose;     
    dataloader::fastdist = fastdist;
    
    dataloader::type = type;
    dataloader::initialised = false;
    load(filename);
    if(gpu){
      initGPUData();
    }
  };

  InputType type;  
  unsigned int bsStride;  
  bool fastdist;

protected:
  static const int BUFFER_SIZE = 16384;
  bool verbose;
  unsigned int sequenceLength;
  unsigned int sequenceCount;
  unsigned int bitStringsCount;
  int paddingLength;
  bool initialised;
  unsigned int BLOCK_SIZE; 
  
private:
  void copy_original_data();
  void sample_fast_dist();
  void sample_char();
};
#endif
