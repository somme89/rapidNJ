#ifndef bitDistanceGap_HPP_INCLUDED
#define bitDistanceGap_HPP_INCLUDED

#include "stdinclude.h"
#include "dataloader.hpp"
#include "DistanceEstimate.hpp"


class bitDistanceGap : public DistanceEstimate {

public:
  bitDistanceGap(dataloader* loader);
  ~bitDistanceGap();
  void computeDistance(int seq1, int seq2, unsigned long long* retVal);
  void printDistances();
  long** distanceMatrixTS;
  long** distanceMatrixTV;

private:
  v4ui** bitStrings;
  v4ui** gapFilters;
  int seqCount;
  int seqLength;
  int bitStringsSize;
  int bitStringsCount;
  static v4ui mask1;
  static v4ui negate_mask1;
  static v4ui mask2;
  static v4ui mask4;
  static v4ui mask8;
  static v4ui mask16;
  static v4ui mask32;
  static v4ui mask64;
  std::vector<string> sequenceNames;
  int seq1, seq2, idx1, idx2;
  int numItr8Bit, numItr16Bit, numItr32Bit, numItrTop, numItr8Bit_cache, numItr16Bit_cache, numItr32Bit_cache, numItrTop_cache;

  v4ui increaseBlockSize4(v4ui);
  v4ui increaseBlockSize8(v4ui);
  v4ui increaseBlockSize16(v4ui);
  v4ui increaseBlockSize32(v4ui);
  v4ui increaseBlockSize64(v4ui);
  void initialiseDataStructures();
  void combine2And4Bit(v4ui&, v4ui&, v4ui&);
  void combine8Bit(v4ui&, v4ui&, v4ui&);
  void combine16Bit(v4ui&, v4ui&, v4ui&);
  void combine32Bit(v4ui&, v4ui&, v4ui&);
  int printSums(v4ui, int);
};

#endif // bitDistanceGap_HPP_INCLUDED
