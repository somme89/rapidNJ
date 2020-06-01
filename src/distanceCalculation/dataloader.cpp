#include "dataloader.hpp"
#include "bitStringUtils.hpp"

dataloader::dataloader(void){
  type = UNKNOWN;
  bitStrings = NULL;
  gapFilters = NULL;
  sequences_original = NULL;
  bitStrings_original = NULL;
  gapFilters_original = NULL;
  BLOCK_SIZE = sizeof(unsigned int) * 8;
  //srand (time(NULL));
  srand(0); 
}

void dataloader::sample_fast_dist() {  
  for(unsigned int i = 0; i < sequenceLength; i++){    
    //select random column    
    unsigned int col = rand() % sequenceLength;
//    unsigned int col = i; //TEST
    unsigned int bits_per_site = 2;
    unsigned int mask = 3;
    if(type == PROTEIN) {
      bits_per_site = 8;
      mask = 255;
    }    
    //cout << "Column: " << col << " -> " << i << endl;    
    for(unsigned int row = 0; row < sequenceCount; row++){      
      unsigned int offset_target = i % (BLOCK_SIZE / bits_per_site);
      unsigned int offset_source = col % (BLOCK_SIZE / bits_per_site);
      unsigned int bitStringIdx_target = i / (BLOCK_SIZE / bits_per_site);
      unsigned int bitStringIdx_source = col / (BLOCK_SIZE / bits_per_site);
      unsigned int* bitString = bitStrings->at(row);
      unsigned int* gapFilter = NULL;
      if(gapFilters != NULL) {
        gapFilter = gapFilters->at(row);
      }
      //cout << "ROW: " << row << " bitstringindex: " << bitStringIdx_source << " " << bitStringIdx_target << endl;
      if(offset_target == 0) {
        bitString[bitStringIdx_target] = 0;
        if(gapFilter != NULL) {
          gapFilter[bitStringIdx_target] = 0;
        }
      }
      unsigned int temp;
      temp = bitStrings_original->at(row)[bitStringIdx_source] >> (offset_source * bits_per_site);      
      temp &= mask;
      bitString[bitStringIdx_target] |= temp << (offset_target * bits_per_site);
      if(gapFilter != NULL) {
        temp = gapFilters_original->at(row)[bitStringIdx_source] >> (offset_source * bits_per_site);      
        temp &= mask;
        gapFilter[bitStringIdx_target] |= temp << (offset_target * bits_per_site);      
      }      
    }
  }
  //printProteinSequence(((v4ui*)(bitStrings_original->at(0)))[0]);
  //printProteinSequence(((v4ui*)bitStrings->at(0))[0]);
}

void dataloader::sample_char() {
  for(unsigned int i = 0; i < sequenceLength; i++){
    //select random column    
    unsigned int col = rand() % sequenceLength;
    //unsigned int col = i;
    for(unsigned int row = 0; row < sequenceCount; row++){
      sequences->at(row)[i] = sequences_original->at(row)[col];
    }
  }
}

void dataloader::copy_original_data() {
  if(fastdist) {    
    bitStrings_original = bitStrings;
    bitStrings = new vector<unsigned int*>;
    if(gapFilters != NULL) {
      gapFilters_original = gapFilters;
      gapFilters = new vector<unsigned int*>;
    }
    for(unsigned int i = 0; i < sequenceCount; i++) {      
      unsigned int* bitString = (unsigned int*) _mm_malloc(bitStringsCount*4*sizeof(unsigned int), 16);
      for(unsigned int j = 0; j < bitStringsCount * 4; j++) {
        bitString[j] = bitStrings_original->at(i)[j];
      }
      bitStrings->push_back(bitString);
      if(gapFilters != NULL) {
        unsigned int* gapFilter = (unsigned int*) _mm_malloc(bitStringsCount*4*sizeof(unsigned int), 16);
        for(unsigned int j = 0; j < bitStringsCount * 4; j++) {
          gapFilter[j] = gapFilters_original->at(i)[j];
        }
        gapFilters->push_back(gapFilter);
      }      
    }
  } else {
    sequences_original = sequences;
    sequences = new vector<char*>;
    for(unsigned int i = 0; i < sequenceCount; i++) {
      char* data = new char[sequenceLength];
      sequences->push_back(data);
    }
  }
}

void dataloader::sample_sequences(){
  if((fastdist && bitStrings_original == NULL) || (!fastdist && sequences_original == NULL)) {
    copy_original_data();
  }
  if(fastdist) {
    sample_fast_dist();
  } else {
    sample_char();
  }
}

dataloader::~dataloader(){
  if(sequences_original != NULL) {
    for(unsigned int i = 0; i < sequences_original->size(); i++) {
      delete[] sequences_original->at(i);
    }
    delete sequences_original;
  }
  if(bitStrings_original != NULL) {
    for(unsigned int i = 0; i < bitStrings_original->size(); i++) {
      _mm_free(bitStrings_original->at(i));
    }
    delete bitStrings_original;
  }
  if(gapFilters_original != NULL) {
    for (unsigned int i = 0; i < gapFilters_original->size(); i++) {
      _mm_free(gapFilters_original->at(i));
    }
    delete gapFilters_original;
  }
}
