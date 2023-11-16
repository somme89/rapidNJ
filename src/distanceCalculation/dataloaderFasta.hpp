#ifndef DATALOADER_FASTA_H
#define DATALOADER_FASTA_H

#include "stdinclude.h"
#include "dataloader.hpp"
#include "bitStringUtils.hpp"

class dataloaderFasta : public dataloader{

public:

  dataloaderFasta(){
    sequences = NULL;
    bitStrings = NULL;
    gapFilters = NULL;
    buf = new char[BUFFER_SIZE];
    firstSequence = true;
    sequenceLength = 0;
    sequenceCount = 0;
    is = new ifstream();
    sequenceNames = new vector<string>;
  }

  void load(string filename) {
    if(fastdist){
      bitStrings = new vector<unsigned int*>;
      parseData(filename);      
    } else {
      sequences = new vector<char*>;
      parseData(filename);
    }
    if(sequenceCount <= 1) {
      cerr << "Found " << sequenceCount << " sequences in the input file. At least 2 sequences are needed to build a tree." << endl;
      exit(1);
    }
    if(verbose) {
      cerr << "Input type determined as ";
      if(type == UNKNOWN){
        cerr << "UNKNOWN. Please supply the type as parameter." << endl;
        exit(1);
      } else if(type == DNA){
        cerr << "DNA." << endl;
      } else if(type == PROTEIN){
        cerr << "PROTEIN." << endl;
      } else {
        cerr << "UNKNOWN. Please supply the type as parameter." << endl;
        exit(1);
      }
      cerr << "Number of sequences: " << sequenceCount << endl;
      cerr << "Sequence length: " << sequenceLength << endl;
    } 
  }

  void discoverInputType(){
    if(type != UNKNOWN) return;
    type = DNA;
    for(unsigned int i = 0; i < sequenceLength; i++){
      char c = charBuffer.at(i); 
      if(!(c < 65   ||
        c == 'A' || 
        c == 'a' ||
        c == 'C' ||
        c == 'c' ||
        c == 'G' ||
        c == 'g' ||
        c == 'T' ||
        c == 't' ||
        c == 'U' ||
        c == 'u' ||
        c == 'N' ||
        c == 'n')) {
        type = PROTEIN;
        return;
      }
    }
  }

  void parseData(string filename) {    
    // open the file
    is->open(filename.data(),ifstream::in);
    if (!is->is_open()) {
      cerr << "ERROR: Could not read file: " << filename << "\n";
      exit(1);
    }
    fillBuffer();
    while(buf != NULL) {
      curChar = getNextChar();
      if(curChar == ';') {
        parseComment();
      } else if(curChar == '>') {
        storeSequence(); //Ignored if there is no sequence data
        curChar = getNextChar();
        parseSequenceName();
      } else {
        parseSequenceLine();
      }
    }	
    storeSequence();
  }

  inline void fillBuffer(){    
    if(!is->good()){
      delete[] buf;
      buf = NULL;
      is->close();
      return;
    }    
    is->read(buf,BUFFER_SIZE);
    if (is->gcount() == 0){
      delete[] buf;
      buf = NULL;
      is->close();
      return;
    }
    bufPos = 0;    
  }

  inline char getNextChar() {
    if(bufPos == is->gcount()){
      fillBuffer();
    }
    if(buf != NULL) {
      return buf[bufPos++];
    }
    return 10;
  }

  //Finds the next new line
  void parseComment(){
    while (curChar != 10){
      curChar = getNextChar(); 
    }
  }  

  inline void encodeProteinSequence(unsigned int* bitString, vector<char>& data){
    for(unsigned int i = 0; i < sequenceLength+paddingLength; i++) {
      int offset = i % 4;
      unsigned int bitStringIdx = i / 4;
      if(offset == 0){
        bitString[bitStringIdx] = 0;
      }
      char c = resolveChar(data[i]);
      bitString[bitStringIdx] += (c << (offset*8));
    }
  }

  inline void encodeDNASequence(unsigned int* bitString, unsigned int* gapFilter, vector<char>& data){    
    for(unsigned int i = 0; i < sequenceLength+paddingLength; i++) {
      int offset = i % (BLOCK_SIZE / 2);
      unsigned int bitStringIdx = i / (BLOCK_SIZE / 2);
      if(offset == 0){
        bitString[bitStringIdx] = 0;
        gapFilter[bitStringIdx] = 0;
      }
      char c = data[i];
      switch (c) {
        //Nothing is done for gaps and ambigious nucleotides
      case 'A': 
      case 'a':
        //A has value 00 so nothing is added to the bitStrings
        gapFilter[bitStringIdx] += (Gbin << (offset*2));
        break;
      case 'C':
      case 'c':
        bitString[bitStringIdx] += (Cbin<<(offset*2));
        gapFilter[bitStringIdx] += (Gbin << (offset*2));
        break;      
      case 'G': 
      case 'g':
        bitString[bitStringIdx] += (Gbin<<(offset*2));
        gapFilter[bitStringIdx] += (Gbin << (offset*2));
        break;
      case 'T': 
      case 't':
        bitString[bitStringIdx] += (Tbin<<(offset*2));
        gapFilter[bitStringIdx] += (Gbin << (offset*2));
        break;      
      default:
        break;
      }
    }
  }

  inline char resolveChar(char c){
    //ignore positions with ambigious/unknown nucleotides
    if(type == DNA){
      if (!(c == 'a' || c == 'A' || c == 'c' || c == 'C' || c == 'g' || c == 'G' || c == 't' || c == 'T' || c == 'u' || c == 'U')){
        return '-';
      } else {
        return c;
      }
    } else {
      if(c == '-' || c == '.' || c == 'X' || c == 'x' || c == 'z' || c == 'Z'  || c == 'b' || c == 'B' || c == 'J' || c == 'j' || c == '?') {
        return '-';
      } else {
        return c;
      }
    }
  }

  void parseSequenceLine(){
    while (curChar != 10){	
      if(curChar > 32){
        charBuffer.push_back(curChar);
      }
      curChar = getNextChar();
    }     
  }

  void parseSequenceName(){    
    string name = "";
    while(curChar > 31){
      name.push_back(curChar);
      curChar = getNextChar();
    }
    sequenceNames->push_back(name);
  }

  void storeSequenceFastDist(){
    if(firstSequence){
      if(type == DNA){
        bitStringsCount = sequenceLength / 64 + 6;
        paddingLength = bitStringsCount*64 - sequenceLength;
        gapFilters = new vector<unsigned int*>;
      } else {
        bitStringsCount = sequenceLength / 16 + 8;
        paddingLength = bitStringsCount*16 - sequenceLength;
      }
    }
    // insert paddings
    for(int i = 0; i < paddingLength; i++){
      charBuffer.push_back('-');
    }
    // Encode the data in char vectors.      
    // unsigned int* bitString = new unsigned int[bitStringsCount*4+3];
    // bitString = (unsigned int*) ((((unsigned int) bitString) + 12) & ~15);
    unsigned int* bitString = (unsigned int*) _mm_malloc(bitStringsCount*4*sizeof(unsigned int), 16);
    if(type == DNA){
      unsigned int* gapFilter;
      gapFilter = (unsigned int*) _mm_malloc(bitStringsCount*4*sizeof(unsigned int), 16);
      encodeDNASequence(bitString, gapFilter, charBuffer);  
      gapFilters->push_back(gapFilter);
    } else {
      encodeProteinSequence(bitString, charBuffer);
    }

    bitStrings->push_back(bitString);
    firstSequence = false;    
    sequenceCount++;
    charBuffer.clear();
  }

  void storeSequence(){
    if(charBuffer.size() != 0) {

      if(firstSequence) {
        sequenceLength = charBuffer.size();
        discoverInputType();	
      } else if(charBuffer.size() != sequenceLength){
	cerr << "ERROR: Sequence " << sequenceCount << " has length " << charBuffer.size() << " while the preceding sequences in the alignment which have length " << sequenceLength << "." << endl;	
        exit(1);
      }
      if(fastdist) {
        storeSequenceFastDist();
        return;
      }
      firstSequence = false;
      char* data = new char[sequenceLength];
      for(unsigned int i = 0; i < sequenceLength; i++){	
        data[i] = resolveChar(charBuffer[i]);
      }
      sequences->push_back(data);
      charBuffer.clear();
      sequenceCount++;
    }
  }

  unsigned int** getBitStrings(){
    return &(*bitStrings)[0];
  }

  unsigned int** getGapFilters(){
    return &(*gapFilters)[0];
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

  ~dataloaderFasta(){
    if(sequenceNames != NULL) {
      delete sequenceNames;
    }
    delete is;
    if(sequences != NULL) {
      for(unsigned int i = 0; i < sequences->size(); i++) {
        delete[] sequences->at(i);
      }
      delete sequences;
    }
    if(bitStrings != NULL) {
      for(unsigned int i = 0; i < bitStrings->size(); i++) {
        _mm_free(bitStrings->at(i));
      }
      delete bitStrings;
    }
    if(gapFilters != NULL) {
      for (unsigned int i = 0; i < gapFilters->size(); i++) {
        _mm_free(gapFilters->at(i));
      }
      delete gapFilters;
    }
  }

private:  
  int bufPos;
  char* buf;  
  vector<string>* sequenceNames;
  bool firstSequence;  
  char curChar;  
  ifstream* is;  
  vector<char> charBuffer;  
};

#endif
