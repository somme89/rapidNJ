#include "dataloaderPhylip.hpp"
#include "bitStringUtils.hpp"

using namespace std;

dataloaderPhylip::dataloaderPhylip() {
  dataloaderPhylip::parsedNames = false;
  type = UNKNOWN;
  sequenceNames = new vector<string>;
}

void dataloaderPhylip::load(string filename) {
  if(fastdist){
    loadFastDist(filename);
  } else {
    loadNormal(filename);
  }
}

// Parse the data for fast dist. DNA alignments only!
void dataloaderPhylip::loadFastDist(string filename) {
  char * buf = new char[BUFFER_SIZE];
  char name_buffer[11];
  int curSeq = 0, curCol = 0, curLine = 0;
  // open the file
  ifstream is;
  is.open(filename.data(),ifstream::in);
  if (!is.is_open()) {
    cerr << "ERROR: Could not read file: " << filename << "\n";
    exit(1);
  }
  readHeader(is);  
  createBuffers();
  while (is.good()) {
    is.read(buf,BUFFER_SIZE);
    int i = 0;
    while (i < is.gcount()) {
      char c = buf[i];
      if (c == 10) {
        // new line
        if (curCol > 0) {
          // we have read some data so data on the next line belongs to the next sequence.
          curSeq++;
          curSeq = curSeq % seqCount;
          curCol = 0;
        }
        curLine++;
      } else if (curCol < 10) {
        if (!parsedNames) {
          //part of the sequence name
          name_buffer[curCol] = c;
        }
        curCol++;
      } else if (curCol == 10) {
        //end of name
        if (!(c == 32 || c == 9)) {
          cerr << "ERROR: Sequence name in line: " << curLine << " was not of length 10" << endl;
          exit(1);
        }
        if (!parsedNames) {
          name_buffer[10] = '\0';
          parseName(name_buffer, curSeq);
        }
        curCol++;
      } else {
        // part of sequence
        if (c > 32) {
          buffers[curSeq][bufferSizes[curSeq]] = c;
          bufferSizes[curSeq]++;
          curCol++;
        }
      }
      i++;
    }
    parseBuffers();
  }
  //Pad the buffers
  for (int seq = 0; seq < seqCount; seq++) {
    if (dnaBitStringIdx[seq]+bufferSizes[seq] > seqLength) {
      cerr << "ERROR: sequence with name " << (*sequenceNames)[seq] << " has length >= " << (dnaBitStringIdx[seq]+bufferSizes[seq]) << " which is longer than the stated sequence length of " << seqLength << endl;
      exit(1);
    }
    if (bufferSizes[seq] > 0) {
      for (int k = bufferSizes[seq]; k < 64; k++) {
        if(type == PROTEIN){
          buffers[seq][k] = '-';
        } else {
          buffers[seq][k] = 'A';
        }
      }
      bufferSizes[seq] = 64;
    }
  }
  parseBuffers();  

  if(type == DNA){
    //pad dnaBitStrings
    for (int seq = 0; seq < seqCount; seq++) {
      while (dnaBitStringIdx[seq] != bitStringsCount) {
        bitStrings[seq][dnaBitStringIdx[seq]*4] = 0;
        bitStrings[seq][dnaBitStringIdx[seq]*4+1] = 0;
        bitStrings[seq][dnaBitStringIdx[seq]*4+2] = 0;
        bitStrings[seq][dnaBitStringIdx[seq]*4+3] = 0;
        gapFilters[seq][dnaBitStringIdx[seq]*4] = 0;
        gapFilters[seq][dnaBitStringIdx[seq]*4+1] = 0;
        gapFilters[seq][dnaBitStringIdx[seq]*4+2] = 0;
        gapFilters[seq][dnaBitStringIdx[seq]*4+3] = 0;
        dnaBitStringIdx[seq]++;	
      }
    }  
  } else {
    //pad protein BitStrings
    for (int seq = 0; seq < seqCount; seq++) {
      while (dnaBitStringIdx[seq] != bitStringsCount) {
        bitStrings[seq][dnaBitStringIdx[seq]*4] = 0x2D2D2D2D;
        bitStrings[seq][dnaBitStringIdx[seq]*4+1] = 0x2D2D2D2D;
        bitStrings[seq][dnaBitStringIdx[seq]*4+2] = 0x2D2D2D2D;
        bitStrings[seq][dnaBitStringIdx[seq]*4+3] = 0x2D2D2D2D;
        dnaBitStringIdx[seq]++;
      }
    }  
  }
  is.close();
}

void dataloaderPhylip::readHeader(ifstream& is) {
  bool foundSeqCount = false;
  char c;
  int i = 0;
  char charInput[256];
  while (is.good()) {
    is.get(c);
    if (c == 10 || i >= 256) {
      if (!foundSeqCount) {
        cerr << "ERROR: Number of sequences was not stated" << endl;
        exit(1);
      }
      if (i == 0) {
        cerr << "ERROR: Sequence length was not stated" << endl;
        exit(1);
      }
      charInput[i] = '\0';
      dataloaderPhylip::seqLength = atoi(charInput);
      break;
    } else if ((c == 32 || c == 9) && i != 0) {
      if (foundSeqCount) {
        charInput[i] = '\0';
        i++;
      } else {
        foundSeqCount = true;
        charInput[i] = '\0';
        dataloaderPhylip::seqCount = atoi(charInput);
        i = 0;
      }
    } else if (c >= 48 && c <= 57) {
      charInput[i] = c;
      i++;
    }
  }
  if (verbose) {
    cerr << "Sequence count: " << seqCount <<"\n";
    cerr << "Sequence length: " << seqLength <<"\n";
  }
}

void dataloaderPhylip::parseName(char* input, int seq) {
  sequenceNames->push_back(string(input));
  if (seq == seqCount) {
    parsedNames = true;
  }
}

void dataloaderPhylip::parseBuffers() {
  if(type == UNKNOWN) {
    //Discover inputType
    type = DNA;
    for(int i = 0; i < bufferSizes[0]; i++){
      char c = buffers[0][i];
      if(!(c == 'A' || 
        c == 'a' ||
        c == 'C' ||
        c == 'c' ||
        c == 'G' ||
        c == 'g' ||
        c == 'T' ||
        c == 't' ||
        c == '.' ||
        c == 'N' ||
        c == 'n' ||
        c == '?' ||
        c == '-')){
          type = PROTEIN;
          break;
      }
    }    
  }
  if(!initialised){
    if(verbose){
      cerr << "Input type determined as ";
      if(type == DNA){
        cerr << "DNA " << endl;
      } else {
        cerr << "protein" << endl;
      }
    }
    if(type == DNA){
      bitStringsCount = seqLength / 64 + 6;
    } else {
      bitStringsCount = seqLength / 16 + 8;
    }
    createBitStrings();    
    initialised = true;
  }
  //Determine how to parse the buffers
  if(type == PROTEIN){
    parseBuffersProtein();
    return;
  }  
  parseBuffersWithGaps();
  return;
  

  //parse the buffers as DNA without gaps
  /*for (int seqNumber = 0; seqNumber < seqCount; seqNumber++) {
    int remainder = bufferSizes[seqNumber] % 64;
    int maxNewStrings = bufferSizes[seqNumber] / 64;
    for (int i = 0; i < maxNewStrings; i++) {
      char* seq = &buffers[seqNumber][i*64];
      unsigned int* ptr = &bitStrings[seqNumber][dnaBitStringIdx[seqNumber]*4];
      for (int idx = 0; idx < 4; idx++) {
        ptr[idx] = 0;
        for (int j = 0; j < 16; j++) {
          switch (seq[j + (idx*16)]) {
            //A has value 00 so nothing is added
          case 'C': {
            ptr[idx] += (Cbin<<(j*2));
            break;
                    }
          case 'G': {
            ptr[idx] += (Gbin<<(j*2));
            break;
                    }
          case 'T': {
            ptr[idx] += (Tbin<<(j*2));
            break;
                    } 
          default:
            break;
          }
        }
      }
      dnaBitStringIdx[seqNumber]++;
    }
    bufferSizes[seqNumber] = remainder;
    for (int i = maxNewStrings*64, j = 0; j < remainder; j++, i++) {
      //copy the remaining part of the buffer into the start of the buffer
      buffers[seqNumber][j] = buffers[seqNumber][i];
    }
  }*/
}

void dataloaderPhylip::parseBuffersWithGaps() {
  for (int seqNumber = 0; seqNumber < seqCount; seqNumber++) {
    int remainder = bufferSizes[seqNumber] % 64;
    int maxNewStrings = bufferSizes[seqNumber] / 64;
    for (int i = 0; i < maxNewStrings; i++) {
      char* seq = &buffers[seqNumber][i*64];
      unsigned int* bitStringPtr = &bitStrings[seqNumber][dnaBitStringIdx[seqNumber]*4];
      unsigned int* gapFilterPtr = &gapFilters[seqNumber][dnaBitStringIdx[seqNumber]*4];
      for (int idx = 0; idx < 4; idx++) {
        bitStringPtr[idx] = 0;
        gapFilterPtr[idx] = 0;
        for (int j = 0; j < 16; j++) {
          switch (seq[j + (idx*16)]) {
            //Nothing is done for gaps and 'N'
          case 'A': 
          case 'a':
            //A has value 00 so nothing is added to the bitStrings
            gapFilterPtr[idx] += (Gbin << (j*2));
            break;	  
          case 'C':
          case 'c':
            bitStringPtr[idx] += (Cbin<<(j*2));
            gapFilterPtr[idx] += (Gbin << (j*2));
            break;

          case 'G':
          case 'g':
            bitStringPtr[idx] += (Gbin<<(j*2));
            gapFilterPtr[idx] += (Gbin << (j*2));
            break;          
          case 'T': 
          case 't':
            bitStringPtr[idx] += (Tbin<<(j*2));
            gapFilterPtr[idx] += (Gbin << (j*2));
            break;          	  
          default:
            break;
          }
        }
      }
      dnaBitStringIdx[seqNumber]++;
    }
    bufferSizes[seqNumber] = remainder;
    for (int i = maxNewStrings*64, j = 0; j < remainder; j++, i++) {
      //copy the remaining part of the buffer into the start of the buffer
      buffers[seqNumber][j] = buffers[seqNumber][i];
    }
  }
}

void dataloaderPhylip::parseBuffersProtein() {
  for (int seqNumber = 0; seqNumber < seqCount; seqNumber++) {
    int remainder = bufferSizes[seqNumber] % 16;
    int maxNewStrings = bufferSizes[seqNumber] / 16;
    for (int i = 0; i < maxNewStrings; i++) {
      char* seq = &buffers[seqNumber][i*16];
      unsigned int* ptr = &bitStrings[seqNumber][dnaBitStringIdx[seqNumber]*4];
      for (int idx = 0; idx < 4; idx++) {
        ptr[idx] = 0;
        for (int j = 0; j < 4; j++) {
          char c = seq[j + (idx*4)];	  
          if(c == '-' || c == '.' || c == 'X' || c == 'x' || c == 'z' || c == 'Z'  || c == 'b' || c == 'B' || c == 'J' || c == 'j' || c == '?' ) {
            ptr[idx] += ('-' << (j*8));
          } else {
            ptr[idx] += (c << (j*8));
          }
        }
      }
      dnaBitStringIdx[seqNumber]++;
    }
    bufferSizes[seqNumber] = remainder;
    for (int i = maxNewStrings*16, j = 0; j < remainder; j++, i++) {
      //copy the remaining part of the buffer into the start of the buffer
      buffers[seqNumber][j] = buffers[seqNumber][i];
    }
  }
}

void dataloaderPhylip::createBuffers() {
  buffers = new char*[seqCount];
  bufferSizes = new int[seqCount];
  for (int i = 0; i < seqCount; i++) {
    buffers[i] = new char[BUFFER_SIZE];
    bufferSizes[i] = 0;
  }
}

void dataloaderPhylip::createBitStrings(){
  dnaBitStringIdx = new int[seqCount];
  bitStrings = (unsigned int**) _mm_malloc(sizeof(unsigned int*) * seqCount, 16);
  for (int i = 0; i < seqCount; i++) {
    bitStrings[i] = (unsigned int*) _mm_malloc(sizeof(unsigned int) * bitStringsCount*4 , 16);
    dnaBitStringIdx[i] = 0;
  }
  
  gapFilters = (unsigned int**) _mm_malloc(sizeof(unsigned int*) * seqCount, 16);
  for (int i = 0; i < seqCount; i++) {
    gapFilters[i] = (unsigned int*) _mm_malloc(sizeof(unsigned int) * bitStringsCount*4 , 16);
  }

}

void dataloaderPhylip::loadNormal(string filename) {
  char * buf = new char[BUFFER_SIZE];
  char name_buffer[11];
  int curSeq = 0, curCol = 0, curLine = 0;
  // open the file
  ifstream is;
  is.open(filename.data(),ifstream::in);
  if (!is.is_open()) {
    cerr << "ERROR: Could not read file: " << filename << "\n";
    exit(1);
  }
  readHeader(is);
  createDataStructuresNaive();
  while (is.good()) {
    is.read(buf,BUFFER_SIZE);
    int i = 0;
    while (i < is.gcount()) {
      char c = buf[i];
      if (c == 10) {
        // new line
        if (curCol > 0) {
          // we have read some data so data on the next line belongs to the next sequence.
          curSeq++;
          curSeq = curSeq % seqCount;
          curCol = 0;
        }
        curLine++;
      } else if (curCol < 10) {
        if (!parsedNames) {
          //part of the sequence name
          name_buffer[curCol] = c;
        }
        curCol++;
      } else if (curCol == 10) {
        //end of name
        if (!(c == 32 || c == 9)) {
          cerr << "ERROR: Sequence name in line: " << curLine << " was not of length 10" << endl;
          exit(1);
        }
        if (!parsedNames) {
          name_buffer[10] = '\0';
          parseName(name_buffer, curSeq);
        }
        curCol++;
      } else {
        // part of sequence
        if (c > 32) {
          sequences[curSeq][sequencesIdx[curSeq]] = c;
          sequencesIdx[curSeq]++;
          curCol++;
        }
      }
      i++;
    }
  }
  is.close();
}

void dataloaderPhylip::createDataStructuresNaive() {
  sequenceNames = new vector<string>;
  sequences = *(new vector<char*>());
  sequences.resize(seqCount);
  sequencesIdx = new int[seqCount];
  for (int i = 0; i < seqCount; i++) {
    sequences[i] = new char[seqLength+1];
    sequencesIdx[i] = 0;
  }
}

unsigned int** dataloaderPhylip::getBitStrings(){
  return bitStrings;
}

unsigned int** dataloaderPhylip::getGapFilters(){
  return gapFilters;
}

unsigned int dataloaderPhylip::getSequenceCount(){
  return seqCount;
}

unsigned int dataloaderPhylip::getSequenceLength(){
  return seqLength;
}

unsigned int dataloaderPhylip::getBitStringsCount(){
  return bitStringsCount;
}

vector<string>* dataloaderPhylip::getSequenceNames(){
  return sequenceNames;
}

vector<char*>* dataloaderPhylip::getSequences(){
  return &sequences;
}

void dataloaderPhylip::setSequences(vector<char*>* val){
  sequences = *val;
}

dataloaderPhylip::~dataloaderPhylip() {

}
