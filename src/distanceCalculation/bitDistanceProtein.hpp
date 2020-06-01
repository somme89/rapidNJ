#ifndef bitDistanceProtein_HPP_INCLUDED
#define bitDistanceProtein_HPP_INCLUDED

#include "stdinclude.h"
#include "dataloader.hpp"
#include "DistanceEstimate.hpp"
#include <math.h>
#include "bitStringUtils.hpp"

class bitDistanceProtein : public DistanceEstimate{

public:
  bitDistanceProtein(dataloader* loader);
  ~bitDistanceProtein();
  void computeDistance(int i, int j, unsigned long long* data);
  void printDistances();
  long** distanceMatrixTS;
  long** distanceMatrixTV;

private:
  v4ui** bitStrings;
  int seqCount;
  int seqLength;
  int bitStringsSize;
  int bitStringsCount;
  static v4ui mask1of8;
  static v4ui mask8;
  static v4ui mask16;
  static v4ui gapMask;
  std::vector<string> sequenceNames;
  int seq1, seq2, idx;
  int numItr8Bit, numItr16Bit, numItr32Bit, numItrTop, numItr8Bit_cache, numItr16Bit_cache, numItr32Bit_cache, numItrTop_cache;
  int paddedSeqLength;

  v4ui increaseBlockSize16(v4ui);
  v4ui increaseBlockSize32(v4ui);
  void sum8_8BitVectors(v4ui&,v4ui&);
  void combineInto16BitVectors(v4ui&,v4ui&);
  void combineInto32BitVectors(v4ui&,v4ui&);
  void sum32Bit(v4ui&,v4ui&);
};

#endif // bitDistanceProtein_HPP_INCLUDED
