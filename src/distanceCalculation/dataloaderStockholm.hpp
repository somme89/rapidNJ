#ifndef DATALOADER_STOCKHOLM_H
#define DATALOADER_STOCKHOLM_H

#include "stdinclude.h"
#include "dataloader.hpp"
#include "bitStringUtils.hpp"

class dataloaderStockholm : public dataloader{

public:
  
  dataloaderStockholm() {
    buf = new char[BUFFER_SIZE];
    firstSequence = true;
    sequenceLength = 0;
    sequenceCount = 0;
    tempBuffer = "";
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
    if(verbose){
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

  void discoverInputType(char* data){
    if(type != UNKNOWN) return;
    type = DNA;
    for(unsigned int i = 0; i < sequenceLength; i++){
      if(!(data[i] < 65   ||
        data[i] == 'A' || 
        data[i] == 'a' ||
        data[i] == 'C' ||
        data[i] == 'c' ||
        data[i] == 'G' ||
        data[i] == 'g' ||
        data[i] == 'T' ||
        data[i] == 't' ||
        data[i] == 'U' ||
        data[i] == 'u' ||
        data[i] == 'N' ||
        data[i] == 'n'))
      {
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
    curChar = getNextChar();

    while(true) {
      while(curChar < 32){
        curChar = getNextChar();
      }
      if(curChar == 35) {
        // # markup line.
        parseMarkupLine();
      } else if(curChar == 47) {
        // '/' might be the end of the alignment.
        curChar = getNextChar();
        if(curChar == 47){
          //finished
          is->close();
          return;
        } else {
          tempBuffer = "/";
          parseSequenceLine();
        }
      } else {
        parseSequenceLine();
      }
      curChar = getNextChar();
    }
  }

  inline void fillBuffer(){    
    if(is->good()){
      is->read(buf,BUFFER_SIZE);
      bufPos = 0;
    } else {
      is->close();
      cerr << "Stream ended unexpectedly. Did not find '//' marking the end of the alignment" << endl;
      exit(1);
    }
  }

  inline char getNextChar() {
    if(bufPos == is->gcount()){
      fillBuffer();
    }
    return buf[bufPos++];
  }

  //Finds the next new line
  void parseMarkupLine(){
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
      int offset = i % 16;
      unsigned int bitStringIdx = i / 16;
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

  void parseSequenceLineFastDist(){
    // find the length and type of the sequence
    parseSequenceName();
    while (curChar != 10){
      if(curChar > 32){
        charBuffer.push_back(curChar);	 
      }
      curChar = getNextChar();
    }
    if(firstSequence){
      sequenceLength = charBuffer.size();
      discoverInputType(&charBuffer[0]);
      if(type == DNA){
        bitStringsCount = sequenceLength / 64 + 6;
        paddingLength = bitStringsCount*64 - sequenceLength;
        gapFilters = new vector<unsigned int*>;
      } else {	
        bitStringsCount = sequenceLength / 16 + 8;
        paddingLength = bitStringsCount*16 - sequenceLength;
      }
    }
    if(charBuffer.size() != sequenceLength) {
      cerr << "ERROR: Sequence " << sequenceCount << " has length " << charBuffer.size() << " while the preceding sequences in the alignment has length " << sequenceLength << "." << endl;
      exit(1);
    }
    // insert paddings
    for(int i = 0; i < paddingLength; i++){
      charBuffer.push_back('-');
    }

    // Encode the data in char vectors.      
    // unsigned int* bitString = new unsigned int[bitStringsCount*4+3];
    // bitString = (unsigned int*) ((((unsigned int) bitString) + 12) & ~15);
    unsigned int* bitString;
    bitString = (unsigned int*) _mm_malloc(bitStringsCount*4*sizeof(unsigned int), 16);
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
    if(fastdist){
      parseSequenceLineFastDist();
      return;
    }
    if(firstSequence){
      // find the length of the sequence
      parseSequenceName();      
      while (curChar != 10){	
        if(curChar > 32) {
          charBuffer.push_back(curChar);	 
        }
        curChar = getNextChar();
      }
      // We have the length. Use char arrays to store the char data in.
      sequenceLength = charBuffer.size();
      char* data = new char[sequenceLength];
      sequences->push_back(data);

	  for(unsigned int i = 0; i < sequenceLength; i++){	
        data[i] =  resolveChar(charBuffer[i]);
      }
      charBuffer.clear();
      discoverInputType(sequences->at(0));      
      firstSequence = false;
    } else {
      parseSequenceName();      
      unsigned int counter = 0;
      char* data = new char[sequenceLength];
      sequences->push_back(data);
      while (curChar != 10){        
        if(curChar > 32){
          data[counter] = resolveChar(curChar);
          counter++;
        }
        curChar = getNextChar();
      }
      if(counter != sequenceLength){
        cerr << "ERROR: Sequence " << sequenceCount << " has length " << counter << " while the preceding sequences in the alignment has length " << sequenceLength << "." << endl;
        exit(1);
      }
    }
    sequenceCount++;
  }

  void parseSequenceName(){    
    string name = "";
    if(tempBuffer != ""){
      name.append(tempBuffer);
      tempBuffer = "";
    }
    while(curChar > 32){
      name.push_back(curChar);
      curChar = getNextChar();
    }
    sequenceNames->push_back(name);
    while(curChar < 32){
      curChar = getNextChar();
      if(curChar == 10) {
        cerr << "ERROR: alignment with name \"" << sequenceNames->back() << "\" has zero length" << endl;
        exit(1);
      }
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

  ~dataloaderStockholm(){
    delete sequenceNames;
  }

private:  
  int bufPos;
  char* buf;  
  vector<string>* sequenceNames;
  bool firstSequence;  
  char curChar;  
  string tempBuffer;
  ifstream* is;  
  vector<char> charBuffer;
  
};

#endif
