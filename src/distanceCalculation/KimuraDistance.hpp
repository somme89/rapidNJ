#ifndef KIMURADISTANCE_HPP
#define KIMURADISTANCE_HPP

#include "stdinclude.h"
#include "dataloader.hpp"
#include "DistanceEstimate.hpp"
#include "diskMatrix.h"

class KimuraDistance {

public:
  KimuraDistance(bool verbose, bool fastdist, dataloader* loader, diskMatrix* dm);
  KimuraDistance(bool verbose, bool fastdist, dataloader* loader, distType** _distMatrix, diskMatrix* dm);
  ~KimuraDistance(void);
  static void* distThread(void* ptr);
  distType** getDistanceMatrix();
  void computeDistances(int);
  void computeDistancesGPU();

private:
  DistanceEstimate* getDistanceEstimateInstance(dataloader* loader);
  void postProcessDistanceMatrix();
  void computeDistanceMatrixMT(int);
  void computeDistancesDNAGPU();
  void computeDistancesProteinGPU();
  void computeDistancesDNAGPU2();

  bool verbose;
  bool fastdist;
  bool popcnt;
  int seqCount;
  dataloader* loader;
  vector<string> sequenceNames;
  long** distanceMatrixTS;
  long** distanceMatrixTV;
  distType** distMatrix;
  distType maxDistance;
  bool gpuInitialised;
  diskMatrix* dm;
};

struct threadStateKimura{
  int seqIdx;
  unsigned int seqCount;
  dataloader* loader;
  distType** distMatrix;
  bool verbose;
  int numThreads;
  distType maxDistance;
  DistanceEstimate* distanceCalculator;
  volatile int availableBuffers;
  diskMatrix* dm;
  pthread_mutex_t mutex;
};
#endif
