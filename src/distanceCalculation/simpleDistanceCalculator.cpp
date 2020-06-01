#include "simpleDistanceCalculator.hpp"

simpleDistanceCalculator::simpleDistanceCalculator(dataloader* loader) : DistanceEstimate(loader) {
  simpleDistanceCalculator::sequences = *loader->getSequences();
  simpleDistanceCalculator::seqLength = loader->getSequenceLength();
}

void simpleDistanceCalculator::computeDistance(int i, int j, unsigned long long* data) {
  if(type == DNA) {
    calculateTsTv(i, j, data);
  } else if(type == PROTEIN){
    calculateHammingDistance(i, j, data);
  } else {
    cerr << "ERROR: simpleDistanceCalculator - unkown sequence type \"" << type << "\"" << endl; 
  }
}


/*
distType simpleDistanceCalculator::calculateDistance(int i, int j) {
distType distance = 0;
distType nucleotideCount = 0;
for (int k = 0; k < seqLength; k++) {
if(sequences[i][k] != '.' && sequences[i][k] != '-' && sequences[j][k] != '.' && sequences[j][k] != '-'){
// not a gap
nucleotideCount++;
if(sequences[i][k] != sequences[j][k]){
distance++;
}
}
}
if(nucleotideCount == 0){    
return -1;
}
//cout << distance << "/" << nucleotideCount << endl;
return distance / nucleotideCount;
}
*/

void simpleDistanceCalculator::calculateHammingDistance(int i, int j, unsigned long long* data) {
  unsigned int distance = 0;
  unsigned int ungappedLength = 0;
  for (int k = 0; k < seqLength; k++) {
    if(sequences[i][k] != '.' && sequences[i][k] != '-' && sequences[j][k] != '.' && sequences[j][k] != '-'){
      // not a gap
      if(sequences[i][k] != sequences[j][k]){
        distance++;        
      }
      ungappedLength++;
    }
  }
  data[0] = distance;
  data[1] = 0;
  data[2] = ungappedLength;
}

// Calculates the distance between two sequences. Transversion and Transitions are calculated
void simpleDistanceCalculator::calculateTsTv(int i, int j, unsigned long long* data) {
  unsigned int ts = 0;
  unsigned int tv = 0;
  unsigned int length = 0;
  for (int k = 0; k < seqLength; k++) {
    if( sequences[i][k] == '-' || sequences[j][k] == '-'){
      continue;
    }
    if(sequences[i][k] != sequences[j][k]){    
      if (sequences[i][k] == 'A') {
        if (sequences[j][k] == 'G') {
          ts++;
        } else {
          tv++;
        }
      } else if (sequences[i][k] == 'C') {
        if (sequences[j][k] == 'T') {
          ts++;
        } else {
          tv++;
        }
      } else if (sequences[i][k] == 'G') {
        if (sequences[j][k] == 'A') {
          ts++;
        } else {
          tv++;	  
        }
      } else if (sequences[i][k] == 'T') {
        if (sequences[j][k] == 'C') {
          ts++;
        } else {
          tv++;
        }
      }
    }
    length++;
  }
  data[0] = ts;
  data[1] = tv;
  data[2] = length;
}

simpleDistanceCalculator::~simpleDistanceCalculator(){}
