#include "bitDistanceProtein.hpp"
#include <emmintrin.h>

v4ui bitDistanceProtein::mask1of8 = _mm_set_epi32(0x01010101,0x01010101,0x01010101,0x01010101);
v4ui bitDistanceProtein::mask8 = _mm_set_epi32(0x00ff00ff,0x00ff00ff,0x00ff00ff,0x00ff00ff);
v4ui bitDistanceProtein::mask16 = _mm_set_epi32(0x0000ffff,0x0000ffff,0x0000ffff,0x0000ffff);
v4ui bitDistanceProtein::gapMask = _mm_set_epi32(0x2d2d2d2d,0x2d2d2d2d,0x2d2d2d2d,0x2d2d2d2d);

bitDistanceProtein::bitDistanceProtein(dataloader* loader) : DistanceEstimate(loader) {
  bitDistanceProtein::bitStrings = (v4ui**) loader->getBitStrings();
  bitDistanceProtein::seqCount = loader->getSequenceCount();
  bitDistanceProtein::seqLength = loader->getSequenceLength();
  bitDistanceProtein::bitStringsCount = loader->getBitStringsCount();
  bitDistanceProtein::sequenceNames = *loader->getSequenceNames();
  bitDistanceProtein::paddedSeqLength = (seqLength + 127) & ~127;
  numItr8Bit_cache = bitStringsCount / 8; // bitstrings are padded to a multiple of 8
  numItr16Bit_cache = numItr8Bit_cache / 31;
  if(numItr8Bit_cache % 31 != 0){
    numItr16Bit_cache++;
  }
  numItr32Bit_cache = numItr16Bit_cache / 128;
  if(numItr16Bit_cache % 128 != 0){
    numItr32Bit_cache++;
  }
  numItrTop_cache = numItr32Bit_cache / 32768;
  if(numItr32Bit_cache % 32768 != 0){
    numItrTop_cache++;
  }
}

void bitDistanceProtein::computeDistance(int i, int j, unsigned long long* data) {
  numItr8Bit = numItr8Bit_cache;
  numItr16Bit = numItr16Bit_cache;
  numItr32Bit = numItr32Bit_cache;
  numItrTop = numItrTop_cache;
  long totalDist = 0;
  long totalGap = 0;
  bitDistanceProtein::seq1 = i;
  bitDistanceProtein::seq2 = j;
  idx = 0;
  for (int k = 0; k < numItrTop; k++) {
    v4ui_vector sumDist;
    v4ui_vector sumGap;
    sumDist.v = _mm_setzero_si128();
    sumGap.v = _mm_setzero_si128();

    sum32Bit(sumDist.v, sumGap.v);
    totalDist += sumDist.d[0];
    totalDist += sumDist.d[1];
    totalDist += sumDist.d[2];
    totalDist += sumDist.d[3];
    totalGap += sumGap.d[0];
    totalGap += sumGap.d[1];
    totalGap += sumGap.d[2];
    totalGap += sumGap.d[3];
  }
  data[0] = totalDist;
  data[1] = 0;
  data[2] = paddedSeqLength - totalGap;
  //cout << seq1 << " " << seq2 << ": " << totalDist << " " <<  paddedSeqLength - totalGap << endl;
}

/** Increases block size from 8 bits to 16 bits */
inline v4ui bitDistanceProtein::increaseBlockSize16(v4ui v) {
    return _mm_add_epi32(_mm_and_si128(v,mask8), _mm_and_si128(_mm_srli_epi32(v,8),mask8));
}

/** Increases block size from 16 bits to 32 bits */
inline v4ui bitDistanceProtein::increaseBlockSize32(v4ui v) {
  return _mm_add_epi32(_mm_and_si128(v,mask16), _mm_and_si128(_mm_srli_epi32(v,16),mask16));
}

inline void bitDistanceProtein::sum8_8BitVectors(v4ui& sumDist, v4ui& sumGap) {
  v4ui r, g;  
  //unrolled by the compiler
  for(int i = 0; i < 8; i++){    
    r = _mm_cmpeq_epi8(bitStrings[seq1][idx],bitStrings[seq2][idx]);
    r = _mm_andnot_si128(r,mask1of8);
    g = _mm_and_si128(_mm_or_si128(_mm_cmpeq_epi8(bitStrings[seq1][idx],gapMask),_mm_cmpeq_epi8(bitStrings[seq2][idx],gapMask)),mask1of8);  
    sumGap = _mm_add_epi32(g,sumGap);
    sumDist = _mm_add_epi32(_mm_andnot_si128(g,r),sumDist);  
    idx++;
  }
}

/*sumDist is a 16 bit vector. This methods adds less than 512 to each element of the vector*/
void bitDistanceProtein::combineInto16BitVectors(v4ui& sumDist, v4ui& sumGap) {
  v4ui sum_dist_temp = _mm_setzero_si128();
  v4ui sum_gap_temp = _mm_setzero_si128();
  //TODO inline the sum function
  int iterations = min(numItr8Bit, 31);
  for (int i = 0; i < iterations; i++) {
    sum8_8BitVectors(sum_dist_temp, sum_gap_temp);    
  }
  numItr8Bit -= iterations;
  sumDist = _mm_add_epi32(increaseBlockSize16(sum_dist_temp), sumDist);
  sumGap = _mm_add_epi32(increaseBlockSize16(sum_gap_temp), sumGap);
}

/* sumDist is a 32 bit vector. This methods adds less than 2*65536 to each element of the vector */
void bitDistanceProtein::combineInto32BitVectors(v4ui& sumDist, v4ui& sumGap) {
  int iterations = min(numItr16Bit, 128);
  v4ui sum_dist_temp = _mm_setzero_si128();
  v4ui sum_gap_temp = _mm_setzero_si128();
  for (int i = 0; i < iterations; i++) {
    combineInto16BitVectors(sum_dist_temp, sum_gap_temp);
  }
  sumDist = _mm_add_epi32(increaseBlockSize32(sum_dist_temp), sumDist);
  sumGap = _mm_add_epi32(increaseBlockSize32(sum_gap_temp), sumGap);
  numItr16Bit -= iterations;
}

void bitDistanceProtein::sum32Bit(v4ui& sumDist, v4ui& sumGap) {
  int iterations = min(numItr32Bit, 32768);
  for (int i = 0; i < iterations; i++) {
    combineInto32BitVectors(sumDist, sumGap);
  }
  numItr32Bit -= iterations;
}

bitDistanceProtein::~bitDistanceProtein(){
  
}

