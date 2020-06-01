#ifndef JCDISTANCE_HPP_INCLUDED
#define JCDISTANCE_HPP_INCLUDED

#include "../stdinclude.h"
#include "dataloader.hpp"
#include "bitDistanceGap.hpp"
#include "dnaBitString.hpp"
#include "DistanceEstimate.hpp"
#include "diskMatrix.h"

class JCdistance{

public:
  JCdistance(bool verbose, bool fastdist, dataloader* loader, diskMatrix* dm);
  //  void calculateDistances();  
  static void* distJCThread(void* ptr);
  distType** getDistanceMatrix();
  void computeDistanceMatrix(int numThreads);

private:
  DistanceEstimate* getDistanceEstimateInstance(dataloader* loader);
  void postProcessDistanceMatrix();
  void computeDistanceMatrixMT(int);
  bool verbose;
  bool fastdist;
  int seqCount;
  int seqLength;
  dataloader* loader;
  vector<std::string> sequenceNames;
  long** distanceMatrixTS;
  long** distanceMatrixTV;
  distType** jcDistMatrix;
  distType maxDistance;
  diskMatrix* dm;
};


struct threadStateJC{
  int seqIdx;
  unsigned int seqCount;
  unsigned int seqLength;
  dataloader* loader;
  distType** jcDistMatrix;
  bool verbose;
  int numThreads;
  distType maxDistance;  
  diskMatrix* dm;
  volatile int availableBuffers;
  DistanceEstimate* estimator;
  pthread_mutex_t mutex;
};
#endif
