#include "bitDistanceGap.hpp"
#include <math.h>
#include <pthread.h>
#include "bitStringUtils.hpp"

v4ui bitDistanceGap::mask1 = _mm_set_epi32(0x55555555,0x55555555,0x55555555,0x55555555);
v4ui bitDistanceGap::mask2 = _mm_set_epi32(0x33333333,0x33333333,0x33333333,0x33333333);
v4ui bitDistanceGap::mask4 = _mm_set_epi32(0x0f0f0f0f,0x0f0f0f0f,0x0f0f0f0f,0x0f0f0f0f);
v4ui bitDistanceGap::mask8 = _mm_set_epi32(0x00ff00ff,0x00ff00ff,0x00ff00ff,0x00ff00ff);
v4ui bitDistanceGap::mask16 = _mm_set_epi32(0x0000ffff,0x0000ffff,0x0000ffff,0x0000ffff);

bitDistanceGap::bitDistanceGap(dataloader* loader) : DistanceEstimate(loader) {
  bitDistanceGap::bitStrings = (v4ui**) loader->getBitStrings();
  bitDistanceGap::gapFilters = (v4ui**) loader->getGapFilters();
  bitDistanceGap::seqCount = loader->getSequenceCount();
  bitDistanceGap::seqLength = loader->getSequenceLength();
  bitDistanceGap::bitStringsCount = loader->getBitStringsCount();
  bitDistanceGap::sequenceNames = *loader->getSequenceNames();
  numItr8Bit_cache = bitStringsCount / 6; // at least one since bitstrings are padded to a multiple of 6
  numItr16Bit_cache = numItr8Bit_cache / 8;
  if(numItr8Bit_cache % 8 != 0){
    numItr16Bit_cache++;
  }
  numItr32Bit_cache = numItr16Bit_cache / 128;
  if(numItr16Bit_cache % 128 != 0){
    numItr32Bit_cache++;
  }
  numItrTop_cache = numItr32Bit_cache / 32768 + 1;
  if(numItr32Bit_cache % 32768 != 0){
    numItrTop_cache++;
  }
  /*    pthread_t tid;
  sched_param param;
  int priority;
  int policy;
  int ret;
  // scheduling parameters of target thread
  ret = pthread_getschedparam (pthread_self(), &policy, &param);
  // sched_priority contains the priority of the thread 
  priority = param.sched_priority;
  cout << pthread_self() << ": priority=" << priority << " policy=" << policy << endl; 
  */
}

int bitDistanceGap::printSums(v4ui value, int bits) {
  v4ui_vector v;
  v.v = value;
  int sum = 0;
  for (int i = 0; i < 4; i++) {
    for (unsigned int j = 0; j < (unsigned) (32/bits); j++) {
      unsigned int val = (unsigned int) v.d[i] >> (j*bits) & ((unsigned int) pow((float)2,bits) - 1);
      sum += val;
      cout << val << " + ";
    }
  }
  cout << " = " << sum << endl;
  return sum;
}

void bitDistanceGap::computeDistance(int seq1, int seq2, unsigned long long* retVal) {
  numItr8Bit = numItr8Bit_cache;
  numItr16Bit = numItr16Bit_cache;
  numItr32Bit = numItr32Bit_cache;
  numItrTop = numItrTop_cache;
  long totalTS = 0;
  long totalTV = 0;  
  long totalResidues = 0;
  bitDistanceGap::seq1 = seq1;
  bitDistanceGap::seq2 = seq2;
  idx1 = 0;
  idx2 = 0;
  for (int k = 0; k < numItrTop; k++) {
    v4ui_vector sum_ts; 
    v4ui_vector sum_tv;
    v4ui_vector sum_gf;
    sum_ts.v = _mm_setzero_si128();
    sum_tv.v = _mm_setzero_si128();
    sum_gf.v = _mm_setzero_si128();

    combine32Bit(sum_ts.v,sum_tv.v,sum_gf.v);
    totalTS += sum_ts.d[0];
    totalTS += sum_ts.d[1];
    totalTS += sum_ts.d[2];
    totalTS += sum_ts.d[3];
    totalTV += sum_tv.d[0];
    totalTV += sum_tv.d[1];
    totalTV += sum_tv.d[2];
    totalTV += sum_tv.d[3];
    totalResidues += sum_gf.d[0];
    totalResidues += sum_gf.d[1];
    totalResidues += sum_gf.d[2];
    totalResidues += sum_gf.d[3];
  }

  retVal[0] = totalTS;
  retVal[1] = totalTV;
  retVal[2] = totalResidues;
}

void bitDistanceGap::printDistances() {
  for (int i = 0; i < seqCount; i++) {
    cout << sequenceNames[i] << "\t";
    for (int j = 0; j < seqCount; j++) {
      if (i <= j) {
        cout << "[" << distanceMatrixTS[i][j] << "," << distanceMatrixTV[i][j] << "]\t";
      } else {
        cout << "[" << distanceMatrixTS[j][i] << "," << distanceMatrixTV[j][i] << "]\t";
      }
    }
    cout << endl;
  }
}

/** Increases block size from 2 bits to 4 bits */
inline v4ui bitDistanceGap::increaseBlockSize4(v4ui v) {
  return _mm_add_epi32(_mm_and_si128(v,mask2), _mm_and_si128(_mm_srli_epi32(v,2),mask2));
}

/** Increases block size from 4 bits to 8 bits */
inline v4ui bitDistanceGap::increaseBlockSize8(v4ui v) {
  return _mm_add_epi32(_mm_and_si128(v,mask4), _mm_and_si128(_mm_srli_epi32(v,4),mask4));
}

/** Increases block size from 8 bits to 16 bits */
inline v4ui bitDistanceGap::increaseBlockSize16(v4ui v) {
  return _mm_add_epi32(_mm_and_si128(v,mask8), _mm_and_si128(_mm_srli_epi32(v,8),mask8));
}

/** Increases block size from 16 bits to 32 bits */
inline v4ui bitDistanceGap::increaseBlockSize32(v4ui v) {
  return _mm_add_epi32(_mm_and_si128(v,mask16), _mm_and_si128(_mm_srli_epi32(v,16),mask16));
}

void bitDistanceGap::combine2And4Bit(v4ui& sum_ts, v4ui& sum_tv, v4ui& sum_gf) {
  v4ui r0, r1, r2, tv0, tv1, tv2, ts0, ts1, ts2, gf0, gf1, gf2, sum_tv_temp1, sum_tv_temp2, sum_ts_temp1, sum_ts_temp2, sum_gf_temp1, sum_gf_temp2;

  // Load 3x128 bit strings and gap filters.
  r0 = _mm_xor_si128(bitStrings[seq1][idx1],bitStrings[seq2][idx2]);
  gf0 = _mm_and_si128(gapFilters[seq1][idx1++],gapFilters[seq2][idx2++]);
  r1 = _mm_xor_si128(bitStrings[seq1][idx1],bitStrings[seq2][idx2]);
  gf1 = _mm_and_si128(gapFilters[seq1][idx1++],gapFilters[seq2][idx2++]);
  r2 = _mm_xor_si128(bitStrings[seq1][idx1],bitStrings[seq2][idx2]);
  gf2 = _mm_and_si128(gapFilters[seq1][idx1++],gapFilters[seq2][idx2++]);  

  // Count transversions
  tv0 = _mm_and_si128(_mm_srli_epi32(r0, 1),mask1);
  tv1 = _mm_and_si128(_mm_srli_epi32(r1, 1),mask1);
  tv2 = _mm_and_si128(_mm_srli_epi32(r2, 1),mask1);

  // Count transitions
  ts0 = _mm_andnot_si128(tv0,_mm_and_si128(r0,mask1));
  ts1 = _mm_andnot_si128(tv1,_mm_and_si128(r1,mask1));
  ts2 = _mm_andnot_si128(tv2,_mm_and_si128(r2,mask1));

  // Handle gaps
  tv0 = _mm_and_si128(tv0,gf0);
  tv1 = _mm_and_si128(tv1,gf1);
  tv2 = _mm_and_si128(tv2,gf2);
  ts0 = _mm_and_si128(ts0,gf0);
  ts1 = _mm_and_si128(ts1,gf1);
  ts2 = _mm_and_si128(ts2,gf2);

  //combine three blocks:
  sum_tv_temp1 = _mm_add_epi32(_mm_add_epi32(tv0, tv1),tv2);
  sum_ts_temp1 = _mm_add_epi32(_mm_add_epi32(ts0, ts1),ts2);
  sum_gf_temp1 = _mm_add_epi32(_mm_add_epi32(gf0, gf1),gf2);

  // Load 3x128 bit strings and gap filters.
  r0 = _mm_xor_si128(bitStrings[seq1][idx1],bitStrings[seq2][idx2]);
  gf0 = _mm_and_si128(gapFilters[seq1][idx1++],gapFilters[seq2][idx2++]);
  r1 = _mm_xor_si128(bitStrings[seq1][idx1],bitStrings[seq2][idx2]);
  gf1 = _mm_and_si128(gapFilters[seq1][idx1++],gapFilters[seq2][idx2++]);
  r2 = _mm_xor_si128(bitStrings[seq1][idx1],bitStrings[seq2][idx2]);
  gf2 = _mm_and_si128(gapFilters[seq1][idx1++],gapFilters[seq2][idx2++]);

  // Count transversions
  tv0 = _mm_and_si128(_mm_srli_epi32(r0, 1),mask1);
  tv1 = _mm_and_si128(_mm_srli_epi32(r1, 1),mask1);
  tv2 = _mm_and_si128(_mm_srli_epi32(r2, 1),mask1);

  // Count transitions
  ts0 = _mm_andnot_si128(tv0,_mm_and_si128(r0,mask1));
  ts1 = _mm_andnot_si128(tv1,_mm_and_si128(r1,mask1));
  ts2 = _mm_andnot_si128(tv2,_mm_and_si128(r2,mask1));

  // Handle gaps
  tv0 = _mm_and_si128(tv0,gf0);
  tv1 = _mm_and_si128(tv1,gf1);
  tv2 = _mm_and_si128(tv2,gf2);
  ts0 = _mm_and_si128(ts0,gf0);
  ts1 = _mm_and_si128(ts1,gf1);
  ts2 = _mm_and_si128(ts2,gf2);

  //combine three blocks:
  sum_tv_temp2 = _mm_add_epi32(_mm_add_epi32(tv0, tv1),tv2);
  sum_ts_temp2 = _mm_add_epi32(_mm_add_epi32(ts0, ts1),ts2);
  sum_gf_temp2 = _mm_add_epi32(_mm_add_epi32(gf0, gf1),gf2);  

  //Increase blocksize to 4 bits and combine two blocks:
  sum_ts = _mm_add_epi32(increaseBlockSize4(sum_ts_temp1), increaseBlockSize4(sum_ts_temp2));
  sum_tv = _mm_add_epi32(increaseBlockSize4(sum_tv_temp1), increaseBlockSize4(sum_tv_temp2));
  sum_gf = _mm_add_epi32(increaseBlockSize4(sum_gf_temp1), increaseBlockSize4(sum_gf_temp2));
}

void bitDistanceGap::combine8Bit(v4ui& sum_ts, v4ui& sum_tv, v4ui& sum_gf) {
  v4ui sum_tv_temp;
  v4ui sum_ts_temp;
  v4ui sum_gf_temp;
  if (numItr8Bit < 8) {
    for (int i = 0; i < numItr8Bit; i++) {
      combine2And4Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
      sum_ts = _mm_add_epi32(sum_ts, increaseBlockSize8(sum_ts_temp));
      sum_tv = _mm_add_epi32(sum_tv, increaseBlockSize8(sum_tv_temp));
      sum_gf = _mm_add_epi32(sum_gf, increaseBlockSize8(sum_gf_temp));
    }
    numItr8Bit = 0;
  } else {    
    //loop 1
    combine2And4Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
    sum_ts = _mm_add_epi32(sum_ts, increaseBlockSize8(sum_ts_temp));
    sum_tv = _mm_add_epi32(sum_tv, increaseBlockSize8(sum_tv_temp));
    sum_gf = _mm_add_epi32(sum_gf, increaseBlockSize8(sum_gf_temp));

    //loop 2
    combine2And4Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);    
    sum_ts = _mm_add_epi32(sum_ts, increaseBlockSize8(sum_ts_temp));
    sum_tv = _mm_add_epi32(sum_tv, increaseBlockSize8(sum_tv_temp));
    sum_gf = _mm_add_epi32(sum_gf, increaseBlockSize8(sum_gf_temp));

    //loop 3
    combine2And4Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
    sum_ts = _mm_add_epi32(sum_ts, increaseBlockSize8(sum_ts_temp));
    sum_tv = _mm_add_epi32(sum_tv, increaseBlockSize8(sum_tv_temp));
    sum_gf = _mm_add_epi32(sum_gf, increaseBlockSize8(sum_gf_temp));

    //loop 4
    combine2And4Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
    sum_ts = _mm_add_epi32(sum_ts, increaseBlockSize8(sum_ts_temp));
    sum_tv = _mm_add_epi32(sum_tv, increaseBlockSize8(sum_tv_temp));
    sum_gf = _mm_add_epi32(sum_gf, increaseBlockSize8(sum_gf_temp));

    //loop 5
    combine2And4Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
    sum_ts = _mm_add_epi32(sum_ts, increaseBlockSize8(sum_ts_temp));
    sum_tv = _mm_add_epi32(sum_tv, increaseBlockSize8(sum_tv_temp));
    sum_gf = _mm_add_epi32(sum_gf, increaseBlockSize8(sum_gf_temp));

    //loop 6
    combine2And4Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
    sum_ts = _mm_add_epi32(sum_ts, increaseBlockSize8(sum_ts_temp));
    sum_tv = _mm_add_epi32(sum_tv, increaseBlockSize8(sum_tv_temp));
    sum_gf = _mm_add_epi32(sum_gf, increaseBlockSize8(sum_gf_temp));

    //loop 7
    combine2And4Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
    sum_ts = _mm_add_epi32(sum_ts, increaseBlockSize8(sum_ts_temp));
    sum_tv = _mm_add_epi32(sum_tv, increaseBlockSize8(sum_tv_temp));
    sum_gf = _mm_add_epi32(sum_gf, increaseBlockSize8(sum_gf_temp));

    //loop 8
    combine2And4Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
    sum_ts = _mm_add_epi32(sum_ts, increaseBlockSize8(sum_ts_temp));
    sum_tv = _mm_add_epi32(sum_tv, increaseBlockSize8(sum_tv_temp));
    sum_gf = _mm_add_epi32(sum_gf, increaseBlockSize8(sum_gf_temp));

    numItr8Bit -= 8;
  }
}

void bitDistanceGap::combine16Bit(v4ui& sum_ts, v4ui& sum_tv, v4ui& sum_gf) {
  if (numItr16Bit < 128) {
    for (int i = 0; i < numItr16Bit; i++) {
      v4ui sum_tv_temp = _mm_setzero_si128();
      v4ui sum_ts_temp = _mm_setzero_si128();
      v4ui sum_gf_temp = _mm_setzero_si128();

      combine8Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
      sum_ts_temp = increaseBlockSize16(sum_ts_temp);
      sum_tv_temp = increaseBlockSize16(sum_tv_temp);
      sum_gf_temp = increaseBlockSize16(sum_gf_temp);
      sum_ts = _mm_add_epi32(sum_ts, sum_ts_temp);
      sum_tv = _mm_add_epi32(sum_tv, sum_tv_temp);
      sum_gf = _mm_add_epi32(sum_gf, sum_gf_temp);
    }
    numItr16Bit = 0;
  } else {
    for (int i = 0; i < 128; i++) {
      v4ui sum_tv_temp = _mm_setzero_si128();
      v4ui sum_ts_temp = _mm_setzero_si128();
      v4ui sum_gf_temp = _mm_setzero_si128();

      combine8Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
      sum_ts_temp = increaseBlockSize16(sum_ts_temp);
      sum_tv_temp = increaseBlockSize16(sum_tv_temp);
      sum_gf_temp = increaseBlockSize16(sum_gf_temp);
      sum_ts = _mm_add_epi32(sum_ts, sum_ts_temp);
      sum_tv = _mm_add_epi32(sum_tv, sum_tv_temp);
      sum_gf = _mm_add_epi32(sum_gf, sum_gf_temp);
    }
    numItr16Bit -= 128;
  }
}

void bitDistanceGap::combine32Bit(v4ui& sum_ts, v4ui& sum_tv, v4ui& sum_gf) {
  if (numItr32Bit < 32768) {
    for (int i = 0; i < numItr32Bit; i++) {
      v4ui sum_tv_temp = _mm_setzero_si128();
      v4ui sum_ts_temp = _mm_setzero_si128();
      v4ui sum_gf_temp = _mm_setzero_si128();

      combine16Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
      sum_ts_temp = increaseBlockSize32(sum_ts_temp);
      sum_tv_temp = increaseBlockSize32(sum_tv_temp);
      sum_gf_temp = increaseBlockSize32(sum_gf_temp);
      sum_ts = _mm_add_epi32(sum_ts, sum_ts_temp);
      sum_tv = _mm_add_epi32(sum_tv, sum_tv_temp);
      sum_gf = _mm_add_epi32(sum_gf, sum_gf_temp);
    }
    numItr32Bit = 0;
  } else {
    for (int i = 0; i < 32768; i++) {
      v4ui sum_tv_temp = _mm_setzero_si128();
      v4ui sum_ts_temp = _mm_setzero_si128();
      v4ui sum_gf_temp = _mm_setzero_si128();

      combine16Bit(sum_ts_temp, sum_tv_temp, sum_gf_temp);
      sum_ts_temp = increaseBlockSize32(sum_ts_temp);
      sum_tv_temp = increaseBlockSize32(sum_tv_temp);
      sum_gf_temp = increaseBlockSize32(sum_gf_temp);
      sum_ts = _mm_add_epi32(sum_ts, sum_ts_temp);
      sum_tv = _mm_add_epi32(sum_tv, sum_tv_temp);
      sum_gf = _mm_add_epi32(sum_gf, sum_gf_temp);
    }
    numItr32Bit -= 32768;
  }
}

bitDistanceGap::~bitDistanceGap(){
}

